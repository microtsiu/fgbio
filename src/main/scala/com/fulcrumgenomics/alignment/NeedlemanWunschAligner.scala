/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.alignment

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Mode.{Global, Glocal, Local}
import com.fulcrumgenomics.alignment.NeedlemanWunschAligner._
import enumeratum.EnumEntry
import htsjdk.samtools.CigarOperator

import scala.collection.{immutable, mutable}
import scala.math.max

/** Trait that entires in Mode will extend. */
sealed trait Mode extends EnumEntry

/** Enum to represent alignment modes supported by the NeedlenameWunschAligner */
object Mode extends FgBioEnum[Mode] {
  /** Alignment mode for global pairwise alignment. */
  case object Global extends Mode
  /** Alignment mode for global alignment of query sequence to local region of target sequence. */
  case object Glocal extends Mode
  /** Alignment mode for local pairwise alignment. */
  case object Local extends Mode

  override def values: immutable.IndexedSeq[Mode] = findValues
}


object NeedlemanWunschAligner {
  /** Generates a simple scoring function using the match and mismatch scores. */
  def simpleScoringFunction(matchScore: Int, mismatchScore: Int): (Byte, Byte) => Int = {
    (lhs: Byte, rhs: Byte) => if (lhs == rhs) matchScore else mismatchScore
  }

  /** Creates a NW aligner with fixed match and mismatch scores. */
  def apply(matchScore: Int, mismatchScore: Int, gapOpen: Int, gapExtend: Int, mode: Mode = Global): NeedlemanWunschAligner = {
    new NeedlemanWunschAligner(scoringFunction=simpleScoringFunction(matchScore, mismatchScore), gapOpen, gapExtend, mode=mode)
  }

  /** Directions within the trace back matrix. */
  type Direction = Byte
  val Left: Direction     = 1.toByte
  val Up  : Direction     = 2.toByte
  val Diagonal: Direction = 3.toByte
  val Done: Direction     = 4.toByte

  // NB: the order of LeftAndDiagonal and UpAndDiagonal matters when breaking ties!
  val AllDirections: Seq[Direction]   = Seq(Diagonal, Left, Up)
  val LeftAndUp: Seq[Direction]       = Seq(Left, Up)
  val LeftAndDiagonal: Seq[Direction] = Seq(Diagonal, Left)
  val UpAndDiagonal: Seq[Direction]   = Seq(Diagonal, Up)

  /** The minimum score allowed to start an alignment.  This prevents underflow. */
  val MinStartScore: Int = Int.MinValue / 2
}


object AlignmentMatrix {
  def apply(direction: Direction, queryLength: Int, targetLength: Int): AlignmentMatrix = {
    AlignmentMatrix(direction=direction, scoring=Matrix[Int](queryLength+1, targetLength+1), trace=Matrix[Direction](queryLength+1, targetLength+1))
  }
}

/** A single alignment matrix for a given [[Direction]] storing both the scoring and traceback matrices produce by the aligner. */
case class AlignmentMatrix(direction: Direction, scoring: Matrix[Int], trace: Matrix[Direction])


/**
  * Implementation of the Needleman-Wunsch algorithm for global alignment of two sequences with support for an
  * affine gap penalty.
  *
  * @param scoringFunction a function to score the alignment of a pair of bases
  * @param gapOpen the gap opening penalty
  * @param gapExtend the gap extension penalty
  * @param useEqualsAndX if true use the = and X cigar operators for matches and mismatches,
  *                      else use the M operator for both.
  * @param mode alignment mode to use when generating alignments
  */
class NeedlemanWunschAligner(val scoringFunction: (Byte,Byte) => Int,
                             val gapOpen: Int,
                             val gapExtend: Int,
                             useEqualsAndX: Boolean = true,
                             val mode: Mode = Global) {

  private val (matchOp, mismatchOp) = if (useEqualsAndX) (CigarOperator.EQ, CigarOperator.X) else (CigarOperator.M, CigarOperator.M)

  /** Convenience method that starts from Strings instead of Array[Byte]s. */
  def align(query: String, target: String): Alignment = align(query.getBytes, target.getBytes)

  /**
    * Align two sequences with the current scoring system and mode.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return an [[Alignment]] object describing the optimal global alignment of the two sequences
    */
  def align(query: Array[Byte], target: Array[Byte]): Alignment = {
    val matrices = buildMatrices(query, target)
    generateAlignment(query, target, matrices)
  }

  /**
    * Constructs both the scoring and traceback matrices.
    *
    * Matrix is constructed such that:
    *   - Matrix as nRow * nColumn (i.e. coordinates are (row, column))
    *   - Query sequence is down the side (each row represents a query base)
    *   - Target sequence is along the top (each column represents a target base)
    *   - i is iterating over rows/query bases
    *   - j is iterating over columns/target bases
    *   - Up   == Consumes query base  but not target base, i.e. Insertion vs. the target
    *   - Left == Consumes target base but not query base, i.e. Deletion vs. the target
    *
    * @param query the query sequence
    * @param target the target sequence
    * @return a tuple of the scoring matrix and the traceback matrix
    */
  protected def buildMatrices(query: Array[Byte], target: Array[Byte]): Map[Direction, AlignmentMatrix] = {
    val matrices = Seq(Diagonal, Left, Up).map { direction =>
      direction -> AlignmentMatrix(direction=direction, queryLength=query.length, targetLength=target.length)
    }.toMap

    // Top left corner - allow all to be zero score but then must be careful to initialize Up and Left to have a gap open
    AllDirections.foreach { direction =>
      matrices(direction).scoring(0, 0) = 0
      matrices(direction).trace(0, 0)   = Done
    }

    // Down the left hand side - for global/glocal we have to keep going up, for local we're done
    mode match {
      case Global | Glocal =>
        forloop(from=1, until=query.length+1) { i =>
          LeftAndDiagonal.foreach { direction =>
            matrices(direction).scoring(i, 0) = MinStartScore
            matrices(direction).trace(i, 0)   = Done
          }
          matrices(Up).scoring(i, 0) = matrices(Up).scoring(i-1, 0) + gapExtend + (if (i == 1) gapOpen else 0)
          matrices(Up).trace(i, 0)   = if (i == 1) Diagonal else Up
        }
      case Local =>
        forloop(from=1, until=query.length+1) { i =>
          forloop(from=1, until=query.length+1) { i =>
            LeftAndUp.foreach { direction =>
              matrices(direction).scoring(i, 0) = MinStartScore
              matrices(direction).trace(i, 0)   = Done
            }
            matrices(Diagonal).scoring(i, 0) = 0
            matrices(Diagonal).trace(i, 0)   = Done
          }
        }
    }

    // Along the top - for global we have to keep going left, for glocal/local we're done
    mode match {
      case Global =>
        forloop(from=1, until=target.length+1) { j =>
          UpAndDiagonal.foreach { direction =>
            matrices(direction).scoring(0, j) = MinStartScore
            matrices(direction).trace(0 , j)   = Done
          }
          matrices(Left).scoring(0, j) = matrices(Left).scoring(0, j-1) + gapExtend + (if (j == 1) gapOpen else 0)
          matrices(Left).trace(0, j)   = if (j == 1) Diagonal else Left
        }
      case Glocal | Local =>
        forloop(from=1, until=target.length+1) { j =>
          LeftAndUp.foreach { direction =>
            matrices(direction).scoring(0, j) = MinStartScore
            matrices(direction).trace(0, j)   = Done
          }
          matrices(Diagonal).scoring(0, j) = 0
          matrices(Diagonal).trace(0, j)   = Done
        }
    }

    // The interior of the matrix
    forloop(from=1, until=query.length+1) { i =>
      forloop(from=1, until=target.length+1) { j =>
        // Diagonal matrix can come from the previous diagonal, up, or left
        val maxDiagDirection = AllDirections.maxBy(direction => matrices(direction).scoring(i-1, j-1))
        matrices(Diagonal).scoring(i, j) = matrices(maxDiagDirection).scoring(i-1, j-1) + scoringFunction(query(i-1), target(j-1))
        matrices(Diagonal).trace(i, j)   = maxDiagDirection
        assert(maxDiagDirection == Diagonal || matrices(maxDiagDirection).scoring(i-1, j-1) > matrices(Diagonal).scoring(i-1, j-1))
        // For local mode, we can start anywhere.  We use < 0 to prefer longer alignments
        if (mode == Local && matrices(Diagonal).scoring(i, j) < 0) {
          matrices(Diagonal).scoring(i, j) = 0
          matrices(Diagonal).trace(i, j)   = Done
        }

        // Up matrix can come from diagonal or up
        val maxUpDirection = UpAndDiagonal.maxBy { direction =>
          val gapPenalty = gapExtend + (if (direction == Up) 0 else gapOpen)
          matrices(direction).scoring(i-1, j) + gapPenalty
        }
        matrices(Up).scoring(i, j) = matrices(maxUpDirection).scoring(i-1, j) + gapExtend + (if (maxUpDirection == Up) 0 else gapOpen)
        matrices(Up).trace(i, j)   = maxUpDirection

        // Left matrix can come from diagonal or left
        val maxLeftDirection = LeftAndDiagonal.maxBy { direction =>
          val gapPenalty = gapExtend + (if (direction == Left) 0 else gapOpen)
          matrices(direction).scoring(i, j-1) + gapPenalty
        }
        matrices(Left).scoring(i, j) = matrices(maxLeftDirection).scoring(i, j-1) + gapExtend + (if (maxLeftDirection == Left) 0 else gapOpen)
        matrices(Left).trace(i, j)   = maxLeftDirection
      }
    }

    matrices
  }

  /**
    * Given the scoring and trace back matrices, construct the [[Alignment]] object for the (or one of the)
    * optimal alignments.
    *
    * @param query the query sequence
    * @param target the target sequence
    * @param matrices the scoring and trace back matrices for the [[Left]], [[Up]], and [[Diagonal]] directions.
    * @return an [[Alignment]] object representing the alignment
    */
  protected def generateAlignment(query: Array[Byte], target: Array[Byte], matrices: Map[Direction, AlignmentMatrix]): Alignment = {
    var currOperator: CigarOperator = null
    var currLength: Int = 0
    val elems = new mutable.ArrayBuffer[CigarElem]()

    // For the target sequence, start at the end for Global, or the highest scoring cell in the last row for Glocal, or
    // the highest scoring cell anywhere in the matrix for Local
    var (curI, curJ, curD) = this.mode match {
      case Global =>
        (query.length, target.length, AllDirections.maxBy(d => matrices(d).scoring(query.length, target.length)))
      case Glocal =>
        var (maxScore, maxJ, maxD) = (MinStartScore, -1, Done)
        forloop(from=1, until=target.length+1) { j =>
          val direction = AllDirections.maxBy(d => matrices(d).scoring(query.length, j))
          val score     = matrices(direction).scoring(query.length, j)
          if (score > maxScore) {
            maxScore = score
            maxJ = j
            maxD = direction
          }
        }
        (query.length, maxJ, maxD)
      case Local =>
        var (maxScore, maxI, maxJ, maxD) = (MinStartScore, -1, -1, Done)
        forloop(from=1, until=query.length+1) { i =>
          forloop(from=1, until=target.length+1) { j =>
            val direction = AllDirections.maxBy(d => matrices(d).scoring(i, j))
            val score     = matrices(direction).scoring(i, j)
            if (score > maxScore) {
              maxScore = score
              maxI = i
              maxJ = j
              maxD = direction
            }
          }
        }
        (maxI, maxJ, maxD)
      }

    // The score is always the score from the starting cell
    val score = matrices(curD).scoring(curI, curJ)

    // For global we have to reach the origin, for glocal we just have to reach the top row, and for local we can stop
    // anywhere.  Fortunately, we have initialized the appropriate cells to "Done"
    val done = () => matrices(curD).trace(curI,curJ) == Done
    while (!done()) {
      val nextD = matrices(curD).trace(curI,curJ)
      val op    = curD match {
        case Up       =>
          curI -= 1
          CigarOperator.INSERTION
        case Left     =>
          curJ -= 1
          CigarOperator.DELETION
        case Diagonal =>
          curI -= 1
          curJ -= 1
          if (query(curI) == target(curJ)) this.matchOp else this.mismatchOp
      }
      curD = nextD

      val extending = op == currOperator
      if (extending) {
        currLength += 1
      }
      else {
        if (currLength > 0) elems += CigarElem(currOperator, currLength)
        currOperator = op
        currLength   = 1
      }

      if (done()) elems += CigarElem(currOperator, currLength)
    }

    Alignment(query=query, target=target, queryStart=curI+1, targetStart=curJ+1, cigar=Cigar(elems.reverse), score=score)
  }
}
