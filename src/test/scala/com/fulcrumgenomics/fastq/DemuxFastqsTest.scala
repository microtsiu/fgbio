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
 *
 */

package com.fulcrumgenomics.fastq

import java.nio.file.Files

import com.fulcrumgenomics.FgBioDef.{DirPath, FilePath, PathToFastq}
import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.fastq.FastqDemultiplexer.{DemuxRecord, DemuxResult}
import com.fulcrumgenomics.illumina.{Sample, SampleSheet}
import com.fulcrumgenomics.testing.{ErrorLogLevel, UnitSpec}
import com.fulcrumgenomics.util.{Io, Metric, ReadStructure, SampleBarcodeMetric}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.sopt.cmdline.ValidationException
import org.scalatest.OptionValues

import scala.collection.mutable.ListBuffer
import scala.util.Try

class DemuxFastqsTest extends UnitSpec with OptionValues with ErrorLogLevel {
  import DemuxFastqs._

  private val sampleBarcode1 = "AAAAAAAAGATTACAGA"
  private val sampleBarcode2 = "CCCCCCCCGATTACAGA"
  private val sampleBarcode3 = "GGGGGGGGGATTACAGA"
  private val sampleBarcode4 = "GGGGGGTTGATTACAGA"

  private val sampleSheetPath: FilePath = {
    val path = makeTempFile("SampleSheet", ".csv")
    val lines = s"""
      |[Header],,,,,,,,,
      |IEMFileVersion,4,,,,,,,,
      |Investigator Name,Joe,,,,,,,,
      |Experiment Name,EXPID,,,,,,,,
      |Date,01/01/2000,,,,,,,,
      |Workflow,GenerateFASTQ,,,,,,,,
      |Application,FASTQ Only,,,,,,,,
      |Assay,Assay Name,,,,,,,,
      |Description,The Description,,,,,,,,
      |Chemistry,Amplicon,,,,,,,,
      |,,,,,,,,,
      |[Reads],,,,,,,,,
      |151,,,,,,,,,
      |151,,,,,,,,,
      |,,,,,,,,,
      |[Settings],,,,,,,,,
      |ReverseComplement,0,,,,,,,,
      |Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
      |AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
      |,,,,,,,,,
      |[Data],,,,,,,,,
      |Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Barcode,Sample_Project,Description
      |20000101-EXPID-1,Sample_Name_1,,,,,$sampleBarcode1,Sample_Project_1,Description_1
      |20000101-EXPID-2,Sample_Name_2,,,,,$sampleBarcode2,Sample_Project_2,Description_2
      |20000101-EXPID-3,Sample_Name_3,,,,,$sampleBarcode3,Sample_Project_3,Description_3
      |20000101-EXPID-4,Sample_Name_4,,,,,$sampleBarcode4,Sample_Project_4,Description_4
    """.stripMargin.trim
    Io.writeLines(path=path, lines=lines.split("\n"))
    path
  }

  private val metadataPath: FilePath = {
    val path = makeTempFile("metadata", ".csv")
    val lines = s"""
                   |Sample_ID,Sample_Name,Sample_Plate,Sample_Well,index,index2,Sample_Barcode,Sample_Project,Description
                   |20000101-EXPID-1,Sample_Name_1,,,,,$sampleBarcode1,Sample_Project_1,Description_1
                   |20000101-EXPID-2,Sample_Name_2,,,,,$sampleBarcode2,Sample_Project_2,Description_2
                   |20000101-EXPID-3,Sample_Name_3,,,,,$sampleBarcode3,Sample_Project_3,Description_3
                   |20000101-EXPID-4,Sample_Name_4,,,,,$sampleBarcode4,Sample_Project_4,Description_4
    """.stripMargin.trim
    Io.writeLines(path=path, lines=lines.split("\n"))
    path
  }

  /** Helper method to create a the sample infos from the sample sheet. */
  private def toSampleInfos(structures: Seq[ReadStructure]): Seq[SampleInfo] = {
    val samples: Seq[Sample] = {
      val sampleSheet = SampleSheet(this.sampleSheetPath).map(s => withCustomSampleBarcode(s, "Sample_Barcode"))
      sampleSheet.toSeq :+ DemuxFastqs.unmatchedSample(sampleSheet.size, structures)
    }
    samples.map { sample => SampleInfo(sample=sample, isUnmatched=sample.sampleName==UnmatchedSampleId) }
  }

  /** Helper method to create a [[FastqDemultiplexer]] */
  private def dx(structures: Seq[ReadStructure], mm: Int = 2, md: Int = 1, mn: Int = 1): FastqDemultiplexer = {
    new FastqDemultiplexer(sampleInfos=toSampleInfos(structures), readStructures=structures,
      maxMismatches=mm, minMismatchDelta=md, maxNoCalls=mn)
  }

  private def fq(name: String, bases: String, quals: Option[String]=None, comment: Option[String] = None, readNumber: Option[Int]=None): FastqRecord =
    FastqRecord(name=name, bases=bases, quals=quals.getOrElse("I"*bases.length), comment=comment, readNumber=readNumber)

  private def verifyFragUnmatchedSample(demuxRecord: DemuxResult): Unit = {
    demuxRecord.records.length shouldBe 1
    val record = demuxRecord.records.headOption.value

    record.name shouldBe "frag"
    record.bases shouldBe "A"*100
    record.quals shouldBe "I"*100
    record.pairedEnd shouldBe false
    demuxRecord.sampleInfo.isUnmatched shouldBe true
    record.molecularBarcode shouldBe 'empty
  }

  private def outputDir(): DirPath = {
    val dir = Files.createTempDirectory("DemuxFromInlineSampleBarcodeTest")
    dir.toFile.deleteOnExit()
    dir
  }

  "FastqDemultiplexer" should "demux fragment reads" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B100T")))
    val fastqRecord = fq(name="frag", bases=sampleBarcode1 + "A"*100)

    val demuxRecord = demuxer.demultiplex(fastqRecord)
    demuxRecord.records.length shouldBe 1
    val record = demuxRecord.records.headOption.value

    record.name shouldBe "frag"
    record.bases shouldBe "A"*100
    record.quals shouldBe "I"*100
    record.pairedEnd shouldBe false
    record.molecularBarcode shouldBe 'empty
  }

  it should "demux paired reads with in-line sample barcodes" in {
    val demuxer = dx(structures=Seq(ReadStructure("8B100T"),ReadStructure("9B100T")))
    val fq1 = fq(name="pair", bases="AAAAAAAA" + "A"*100)
    val fq2 = fq(name="pair", bases="GATTACAGA" + "T"*100)

    val demuxRecord = demuxer.demultiplex(fq1, fq2)
    demuxRecord.records.length shouldBe 2
    val r1 = demuxRecord.records.headOption.value
    val r2 = demuxRecord.records.lastOption.value

    r1.name shouldBe "pair"
    r1.bases shouldBe "A"*100
    r1.quals shouldBe "I"*100
    r1.pairedEnd shouldBe true
    r1.readNumber shouldBe 1
    r1.molecularBarcode shouldBe 'empty

    r2.name shouldBe "pair"
    r2.bases shouldBe "T"*100
    r2.quals shouldBe "I"*100
    r2.pairedEnd shouldBe true
    r2.readNumber shouldBe 2
    r2.molecularBarcode shouldBe 'empty
  }

  it should "demux dual-indexed paired end reads" in {
    val demuxer = dx(structures=Seq(ReadStructure("8B"), ReadStructure("100T"), ReadStructure("100T"), ReadStructure("9B")))
    val fq1 = fq(name="pair", bases="AAAAAAAA")
    val fq2 = fq(name="pair", bases="A"*100)
    val fq3 = fq(name="pair", bases="T"*100)
    val fq4 = fq(name="pair", bases="GATTACAGA")

    val demuxRecord = demuxer.demultiplex(fq1, fq2, fq3, fq4)
    demuxRecord.records.length shouldBe 2
    val r1 = demuxRecord.records.headOption.value
    val r2 = demuxRecord.records.lastOption.value

    r1.name shouldBe "pair"
    r1.bases shouldBe "A"*100
    r1.quals shouldBe "I"*100
    r1.pairedEnd shouldBe true
    r1.pairedEnd shouldBe true
    r1.readNumber shouldBe 1
    r1.molecularBarcode shouldBe 'empty

    r2.name shouldBe "pair"
    r2.bases shouldBe "T"*100
    r2.quals shouldBe "I"*100
    r2.pairedEnd shouldBe true
    r2.readNumber shouldBe 2
    r2.molecularBarcode shouldBe 'empty
  }

  it should "demux a very weird set of reads" in {
    val demuxer = dx(structures=Seq(ReadStructure("4B4M8S"), ReadStructure("4B100T"), ReadStructure("100S3B"), ReadStructure("6B1S1M1T")))
    val fq1 = fq(name="pair", bases="AAAACCCCGGGGTTTT")
    val fq2 = fq(name="pair", bases="A"*104)
    val fq3 = fq(name="pair", bases="T"*100 + "GAT")
    val fq4 = fq(name="pair", bases="TACAGAAAT")

    val demuxRecord = demuxer.demultiplex(fq1, fq2, fq3, fq4)
    demuxRecord.records.length shouldBe 2
    val r1 = demuxRecord.records.headOption.value
    val r2 = demuxRecord.records.lastOption.value

    r1.name shouldBe "pair"
    r1.bases shouldBe "A"*100
    r1.quals shouldBe "I"*100
    r1.pairedEnd shouldBe true
    r1.readNumber shouldBe 1
    r1.molecularBarcode.value shouldBe "CCCC-A"
    r1.sampleBarcode.value shouldBe "AAAA-AAAA-GAT-TACAGA"

    r2.name shouldBe "pair"
    r2.bases shouldBe "T"
    r2.quals shouldBe "I"
    r2.pairedEnd shouldBe true
    r2.readNumber shouldBe 2
    r2.molecularBarcode.value shouldBe  "CCCC-A"
    r2.sampleBarcode.value shouldBe "AAAA-AAAA-GAT-TACAGA"
  }

  it should "demux a read structure with multiple independent template segments in one read" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B20T20S20T20S20T"))) // 3 distinct templates segments!
    val fastqRecord = fq(name="frag", bases=sampleBarcode1 + "A"*20 + "C"*20+ "A"*20 + "C"*20+ "A"*20) // template should be A*60

    val demuxRecord = demuxer.demultiplex(fastqRecord)
    demuxRecord.records.length shouldBe 1
    val record = demuxRecord.records.headOption.value

    record.name shouldBe "frag"
    record.bases shouldBe "A"*60
    record.quals shouldBe "I"*60
    record.pairedEnd shouldBe false
    record.molecularBarcode shouldBe 'empty
  }

  it should "fail if zero or more than two read structures have template bases" in {
    an[Exception] should be thrownBy dx(structures=Seq.empty)
    an[Exception] should be thrownBy dx(structures = Seq(ReadStructure("100S")))
    an[Exception] should be thrownBy dx(structures = Seq(ReadStructure("100S"), ReadStructure("100S"), ReadStructure("100S")))
  }

  it should "fail if not enough or too many fastq records are passed to demultiplex()" in {
    // too many
    {
      val demuxer = dx(structures = Seq(ReadStructure("8B92T")))
      val f1      = fq(name="frag", bases="A" * 100)
      val f2      = fq(name="frag", bases="A" * 100)
      val f3      = fq(name="frag", bases="A" * 100)
      an[Exception] should be thrownBy demuxer.demultiplex(f1, f2)
      an[Exception] should be thrownBy demuxer.demultiplex(f1, f2, f3)
    }

    // not enough
    {
      val demuxer = dx(structures = Seq(ReadStructure("100T"), ReadStructure("100T"), ReadStructure("10B")))
      val f1      = fq(name="frag", bases="A" * 100)
      val f2      = fq(name="frag", bases="A" * 100)
      an[Exception] should be thrownBy demuxer.demultiplex(f1, f2)
      an[Exception] should be thrownBy demuxer.demultiplex(f2)
      an[Exception] should be thrownBy demuxer.demultiplex()
    }
  }

  it should "set molecular barcodes" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B5M100T")))
    val fastqRecord = fq(name="frag", bases=sampleBarcode1 + "NNNNN" + "A"*100)

    val demuxRecord = demuxer.demultiplex(fastqRecord)
    demuxRecord.records.length shouldBe 1
    val record = demuxRecord.records.headOption.value

    record.name shouldBe "frag"
    record.bases shouldBe "A"*100
    record.quals shouldBe "I"*100
    record.pairedEnd shouldBe false
    record.molecularBarcode.value shouldBe "NNNNN"
  }

  it should "assign to the 'unmatched' sample if there are too many mismatches" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B100T")), mm=0) // no mismatches allowed
    val fastqRecord = fq(name="frag", bases="AAAAAAAAGATTACAGT" + "A"*100) // last sample barcode base mismatches
    verifyFragUnmatchedSample(demuxer.demultiplex(fastqRecord))
  }

  it should "assign to the 'unmatched' sample if the read matched two sample barcodes within the mismatch delta" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B100T")), mm=10, md=3) // so many mismatches allowed!
    val fastqRecord = fq(name="frag", bases=sampleBarcode4 + "A"*100) // matches the 4th barcode perfectly and the 3rd barcode with two mismatches
    verifyFragUnmatchedSample(demuxer.demultiplex(fastqRecord))
  }

  it should "assign to the 'unmatched' sample if the read's sample barcode has too many Ns" in {
    val demuxer     = dx(structures=Seq(ReadStructure("17B100T")), mm=10, md=1, mn=0) // so many mismatches allowed!
    val fastqRecord = fq(name="frag", bases="GGGGGGTTGATTACAGN" + "A"*100) // one N
    verifyFragUnmatchedSample(demuxer.demultiplex(fastqRecord))
  }

  "FastqDemultiplexer.countMismatches" should "find no mismatches" in {
    FastqDemultiplexer.countMismatches("GATTACA".getBytes, "GATTACA".getBytes) shouldBe 0
  }

  it should "find two mismatches" in {
    FastqDemultiplexer.countMismatches("GATTACA".getBytes, "GACCACA".getBytes) shouldBe 2
  }

  it should "not count no calls" in {
    FastqDemultiplexer.countMismatches("GATTACA".getBytes, "GANNACA".getBytes) shouldBe 0
  }

  it should "find compare two sequences that have all mismatches" in {
    FastqDemultiplexer.countMismatches("GATTACA".getBytes, "CTAATGT".getBytes) shouldBe 7
  }

  private def throwableMessageShouldInclude(msg: String)(r: => Unit): Unit = {
    val result = Try(r)
    result.isFailure shouldBe true
    val failure   = result.failed
    val throwable = failure.get
    val message = throwable match {
      case validationException: ValidationException => validationException.messages.mkString("\n")
      case thr => thr.getMessage
    }
    message should not be null
    message should include (msg)
  }

  // A file containing valid FASTQ records.
  private val fastqPath = {
    val fastqs = new ListBuffer[FastqRecord]()
    fastqs += fq(name="frag1", bases=sampleBarcode1 + "A"*100) // matches the first sample -> first sample
    fastqs += fq(name="frag2", bases="AAAAAAAAGATTACAGT" + "A"*100) // matches the first sample, one mismatch -> first sample
    fastqs += fq(name="frag3", bases="AAAAAAAAGATTACTTT" + "A"*100) // matches the first sample, three mismatches -> unmatched
    fastqs += fq(name="frag4", bases=sampleBarcode4 + "A"*100) // matches the 4th barcode perfectly and the 3rd barcode with two mismatches, delta too small -> unmatched
    fastqs += fq(name="frag5", bases="AAAAAAAAGANNNNNNN" + "A"*100) // matches the first sample, too many Ns -> unmatched

    val path = makeTempFile("test", ".fastq")
    Io.writeLines(path, fastqs.map(_.toString))
    path
  }
  
  "DemuxFastqs"  should "fail if the number of fastqs does not equal the number of read structures" in {
    val structures = Seq(ReadStructure("8B100T"), ReadStructure("100T"))
    val fastq      = PathUtil.pathTo("/path/to/nowhere", ".fastq")
    val output     = PathUtil.pathTo("/path/to/nowhere", "output")
    val metrics    = PathUtil.pathTo("/path/to/nowhere", "metrics")
    throwableMessageShouldInclude("same number of read structures should be given as FASTQs") {
      new DemuxFastqs(inputs=Seq(fastq), output=output, metadata=sampleSheetPath,
        readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3)
    }
  }

  it should "fail if no sample barcodes are found in the read structures" in {
    val structures = Seq(ReadStructure("100T"))
    val fastq      = PathUtil.pathTo("/path/to/nowhere", ".fastq")
    val output     = PathUtil.pathTo("/path/to/nowhere", "output")
    val metrics    = PathUtil.pathTo("/path/to/nowhere", "metrics")
    throwableMessageShouldInclude("No sample barcodes found in read structures") {
      new DemuxFastqs(inputs=Seq(fastq), output=output, metadata=sampleSheetPath,
        readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3)
    }
  }

  it should "fail if there are not 1-2 read structures with template bases" in {
    // 0 read structures with template bases
    {
      val structures = Seq(ReadStructure("8B"))
      val fastq      = PathUtil.pathTo("/path/to/nowhere", ".fastq")
      val output     = PathUtil.pathTo("/path/to/nowhere", "output")
      val metrics    = PathUtil.pathTo("/path/to/nowhere", "metrics")
      throwableMessageShouldInclude("with template bases but expected 1 or 2.") {
        new DemuxFastqs(inputs=Seq(fastq), output=output, metadata=sampleSheetPath,
          readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3)
      }
    }

    // 3 read structures with template bases
    {
      val structures = Seq(ReadStructure("8B92T"), ReadStructure("100T"), ReadStructure("100T"))
      val fastq      = PathUtil.pathTo("/path/to/nowhere", ".fastq")
      val output     = PathUtil.pathTo("/path/to/nowhere", "output")
      val metrics    = PathUtil.pathTo("/path/to/nowhere", "metrics")
      throwableMessageShouldInclude("with template bases but expected 1 or 2.") {
        new DemuxFastqs(inputs=Seq(fastq, fastq, fastq), output=output, metadata=sampleSheetPath,
          readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3)
      }
    }
  }

  it should "fail if there are a different number of sample barcode bases in the read structure compare to the sample sheet" in {
    val output: DirPath = outputDir()
    val metrics = makeTempFile("metrics", ".txt")
    val structures = Seq(ReadStructure("18B100T"))

    val demuxFastqs = new DemuxFastqs(inputs=Seq(fastqPath), output=output, metadata=sampleSheetPath,
      readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3)
    throwableMessageShouldInclude("The number of sample barcodes bases did not match") {
      demuxFastqs.execute()
    }
  }

  it should "fail if there are missing sample barcodes in the sample sheet" in {
    val output: DirPath = outputDir()
    val metrics = makeTempFile("metrics", ".txt")
    val structures = Seq(ReadStructure("17B100T"))

    // Add a sample without a sample barcode
    val sampleSheetLines: Seq[String] = Io.readLines(sampleSheetPath).toSeq
    val buggySampleSheet = makeTempFile("SampleSheet", ".csv")
    Io.writeLines(buggySampleSheet, sampleSheetLines ++ Seq("20000101-EXPID-5,Sample_Name_5,,,,,,Sample_Project_5,Description_5"))

    val demuxFastqs = new DemuxFastqs(inputs=Seq(fastqPath), output=output, metadata=buggySampleSheet,
      readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3)
    throwableMessageShouldInclude("Sample barcode not found in column") {
      demuxFastqs.execute()
    }
  }

  Seq(1, 2).foreach { threads =>
    Seq(true, false).foreach { useSampleSheet =>
      OutputType.values.foreach { outputType =>
        it should s"run end-to-end with $threads threads using a ${if (useSampleSheet) "sample sheet" else "metadata CSV"} $outputType output" in {

          val output: DirPath = outputDir()

          val metrics = makeTempFile("metrics", ".txt")
          val structures = Seq(ReadStructure("17B100T"))

          val metadata = if (useSampleSheet) sampleSheetPath else metadataPath

          new DemuxFastqs(inputs=Seq(fastqPath), output=output, metadata=metadata,
            readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3,
            threads=threads, outputType=Some(outputType)).execute()

          val sampleInfos = toSampleInfos(structures)
          val samples: Seq[Sample] = sampleInfos.map(_.sample)

          // Check the BAMs
          sampleInfos.foreach { sampleInfo =>
            val sample = sampleInfo.sample
            val prefix = toSampleOutputPrefix(sample, sampleInfo.isUnmatched, false, output, UnmatchedSampleId)

            def checkOutput(names: Seq[String], sampleBarcodes: Seq[String]): Unit = {
              if (sample.sampleOrdinal == 1) {
                names.length shouldBe 2
                names should contain theSameElementsInOrderAs Seq("frag1", "frag2")
                sampleBarcodes should contain theSameElementsInOrderAs Seq("AAAAAAAAGATTACAGA", "AAAAAAAAGATTACAGT")
              }
              else if (sample.sampleId == UnmatchedSampleId) {
                names.length shouldBe 3
                names should contain theSameElementsInOrderAs Seq("frag3", "frag4", "frag5")
                sampleBarcodes should contain theSameElementsInOrderAs Seq("AAAAAAAAGATTACTTT", "GGGGGGTTGATTACAGA", "AAAAAAAAGANNNNNNN") // NB: raw not assigned
              }
              else {
                names shouldBe 'empty
                sampleBarcodes.isEmpty shouldBe true
              }
            }

            if (outputType.producesFastq) {
              val fastq = PathUtil.pathTo(prefix + ".fastq.gz")
              val records = FastqSource(fastq).toSeq
              checkOutput(records.map(_.name.replaceAll(":.*", "")), records.map(_.name.replaceAll(".*:", "")))
            }
            if (outputType.producesBam) {
              val bam = PathUtil.pathTo(prefix + ".bam")
              val records = readBamRecs(bam)
              checkOutput(records.map(_.name), records.map(_[String]("BC")))
            }
          }

          // Check the metrics
          val sampleBarcodMetrics = Metric.read[SampleBarcodeMetric](metrics)

          sampleBarcodMetrics.length shouldBe samples.length

          sampleBarcodMetrics.zip(samples).foreach { case (metric, sample) =>
            metric.barcode_name shouldBe sample.sampleName
            metric.barcode      shouldBe sample.sampleBarcodeString
            metric.library_name shouldBe sample.libraryId
            if (sample.sampleOrdinal == 1) {
              metric.templates                   shouldBe 2
              metric.pf_templates                shouldBe 2
              metric.perfect_matches         shouldBe 1
              metric.one_mismatch_matches    shouldBe 1
              metric.pf_perfect_matches      shouldBe 1
              metric.pf_one_mismatch_matches shouldBe 1
            }
            else if (sample.sampleId == UnmatchedSampleId) {
              metric.templates                   shouldBe 3
              metric.pf_templates                shouldBe 3
              metric.perfect_matches         shouldBe 0
              metric.one_mismatch_matches    shouldBe 0
              metric.pf_perfect_matches      shouldBe 0
              metric.pf_one_mismatch_matches shouldBe 0
            }
            else {
              metric.templates                   shouldBe 0
              metric.pf_templates                shouldBe 0
              metric.perfect_matches         shouldBe 0
              metric.one_mismatch_matches    shouldBe 0
              metric.pf_perfect_matches      shouldBe 0
              metric.pf_one_mismatch_matches shouldBe 0
            }
          }
        }
      }
    }
  }

  Seq(true, false).foreach { illuminaStandards =>
    // A file containing illumina-style read names
    val illuminaReadNamesFastqPath = {
      val fastqs = new ListBuffer[FastqRecord]()
      fastqs += fq(name="Instrument:RunID:FlowCellID:Lane:Tile:X:1", comment=Some("ReadNum:FilterFlag:0:SampleNumber"), bases=sampleBarcode1 + "A"*100) // matches the first sample -> first sample
      fastqs += fq(name="Instrument:RunID:FlowCellID:Lane:Tile:X:2", comment=Some("ReadNum:FilterFlag:0:SampleNumber"), bases="AAAAAAAAGATTACAGT" + "A"*100) // matches the first sample, one mismatch -> first sample
      fastqs += fq(name="Instrument:RunID:FlowCellID:Lane:Tile:X:3", comment=Some("ReadNum:FilterFlag:0:SampleNumber"), bases="AAAAAAAAGATTACTTT" + "A"*100) // matches the first sample, three mismatches -> unmatched
      fastqs += fq(name="Instrument:RunID:FlowCellID:Lane:Tile:X:4", comment=Some("ReadNum:FilterFlag:0:SampleNumber"), bases=sampleBarcode4 + "A"*100) // matches the 4th barcode perfectly and the 3rd barcode with two mismatches, delta too small -> unmatched
      fastqs += fq(name="Instrument:RunID:FlowCellID:Lane:Tile:X:5", comment=Some("ReadNum:FilterFlag:0:SampleNumber"), bases="AAAAAAAAGANNNNNNN" + "A"*100) // matches the first sample, too many Ns -> unmatched

      val path = makeTempFile("test", ".fastq")
      Io.writeLines(path, fastqs.map(_.toString))
      path
    }

    it should s"${if (illuminaStandards) "not " else ""}append /1 and /2 to read pairs when outputting FASTQs" in {
      val output     = outputDir()
      val metrics    = makeTempFile("metrics", ".txt")
      val structures = Seq(ReadStructure("17B100T"), ReadStructure("117T"))
      new DemuxFastqs(inputs=Seq(illuminaReadNamesFastqPath, illuminaReadNamesFastqPath), output=output, metadata=sampleSheetPath,
        readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3,
        outputType=Some(OutputType.Fastq), illuminaStandards=illuminaStandards).execute()

      val sampleInfos = toSampleInfos(structures)

      // Check the BAMs
      sampleInfos.foreach { sampleInfo =>
        val sample     = sampleInfo.sample
        val prefix     = toSampleOutputPrefix(sample, isUnmatched=sampleInfo.isUnmatched, illuminaStandards=illuminaStandards, output, UnmatchedSampleId)
        val extensions = FastqRecordWriter.extensions(pairedEnd=true, illuminaStandards=illuminaStandards)

        def toHeader(rec: FastqRecord): String = {
          if (illuminaStandards) rec.header.endsWith(" ReadNum:FilterFlag:0:SampleNumber") shouldBe true
          else if (rec.readNumber.value == 1) rec.header.endsWith("/1") shouldBe true
          else  rec.header.endsWith("/2") shouldBe true
          rec.header.split(" ").head.replaceAll(".*:", "").replaceAll("/[12]$", "")
        }

        val headersR1 = {
          val fastq = PathUtil.pathTo(prefix + extensions.head)
          FastqSource(fastq).toSeq.map(toHeader).toList
        }
        val headersR2 = {
          val fastq = PathUtil.pathTo(prefix + extensions.last)
          FastqSource(fastq).toSeq.map(toHeader).toList
        }

        headersR1 should contain theSameElementsInOrderAs headersR2

        if (sample.sampleOrdinal == 1) {
          headersR1.length shouldBe 2
          if (illuminaStandards) headersR1 should contain theSameElementsInOrderAs Seq("1", "2")
          else headersR1 should contain theSameElementsInOrderAs Seq("AAAAAAAAGATTACAGA", "AAAAAAAAGATTACAGT")
        }
        else if (sample.sampleId == UnmatchedSampleId) {
          headersR1.length shouldBe 3
          if (illuminaStandards) headersR1 should contain theSameElementsInOrderAs Seq("3", "4", "5")
          else headersR1 should contain theSameElementsInOrderAs Seq("AAAAAAAAGATTACTTT", sampleBarcode4, "AAAAAAAAGANNNNNNN")
        }
        else {
          headersR1 shouldBe 'empty
        }
      }
    }
  }

  it should "demultiplex fragment reads with standard qualities" in {
    val output     = outputDir()
    val metrics    = makeTempFile("metrics", ".txt")

    val structures=Seq(ReadStructure("4B+T"))
    val fastqRecord = fq(name="frag", bases="TCGTGACAGGAATCAAATGAAAACACTTGGT", quals=Some("1>1>11AFF?FFGGGGGGCBDGBGCFGH11B"))
    val illuminaReadNamesFastqPath = {
      val fastqs = Seq(fastqRecord)
      val path = makeTempFile("test", ".fastq")
      Io.writeLines(path, fastqs.map(_.toString))
      path
    }

    val metadata = {
      val data = "Sample_Barcode,Sample_ID,Sample_Name,Library_ID,Description\nTCGT,foo,S1,,"
      val path = makeTempFile("test", ".fastq")
      Io.writeLines(path, data.split("\n"))
      path
    }

    new DemuxFastqs(inputs=Seq(illuminaReadNamesFastqPath), output=output, metadata=metadata,
      readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3,
      outputType=Some(OutputType.Bam)).execute()

    val record = readBamRecs(output.resolve("foo-S1-TCGT.bam")).head

    record.name shouldBe "frag"
    record.basesString shouldBe fastqRecord.bases.drop(4)
    record.qualsString shouldBe fastqRecord.quals.drop(4)
  }

  it should "create the correct output prefix" in {
    val structures = Seq(ReadStructure("17B100T"))
    val sampleInfos = toSampleInfos(structures)

    val matchedSampleInfo   = sampleInfos.find(!_.isUnmatched).get
    val unmatchedSampleInfo = sampleInfos.find(_.isUnmatched).get
    val output = outputDir()

    {
      val prefix = toSampleOutputPrefix(matchedSampleInfo.sample, matchedSampleInfo.isUnmatched, false, output, UnmatchedSampleId)
      prefix.toString shouldBe output.resolve("20000101-EXPID-1-Sample_Name_1-AAAAAAAAGATTACAGA").toString
    }

    {
      val prefix = toSampleOutputPrefix(unmatchedSampleInfo.sample, unmatchedSampleInfo.isUnmatched, false, output, UnmatchedSampleId)
      prefix.toString shouldBe output.resolve("unmatched").toString
    }

    {
      val prefix = toSampleOutputPrefix(matchedSampleInfo.sample, matchedSampleInfo.isUnmatched, true, output, UnmatchedSampleId)
      prefix.toString shouldBe output.resolve("Sample_Name_1_S1_L001").toString
    }

    {
      val prefix = toSampleOutputPrefix(unmatchedSampleInfo.sample, unmatchedSampleInfo.isUnmatched, true, output, UnmatchedSampleId)
      prefix.toString shouldBe output.resolve("unmatched_S4_L001").toString
    }
  }

  it should "fail if one FASTQ has fewer records than the other" in {
    val fastq1 = makeTempFile("test", ".R1.fastq")
    val fastq2 = makeTempFile("test", ".R2.fastq")

    def toFq(i: Int): FastqRecord = fq(name=s"frag$i", bases=sampleBarcode1 + "A"*100) // matches the first sample -> first sample

    Io.writeLines(fastq1, Seq(1,2).map(toFq(_).toString))
    Io.writeLines(fastq2, Seq(1).map(toFq(_).toString))

    val output = outputDir()

    val metrics = makeTempFile("metrics", ".txt")

    throwableMessageShouldInclude("out of sync") {
      new DemuxFastqs(inputs=Seq(fastq1, fastq2), output=output, metadata=sampleSheetPath,
        readStructures=Seq(ReadStructure("17B100T"), ReadStructure("100T")), metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3
      ).execute()
    }
  }


  private def write(rec: FastqRecord): PathToFastq = {
    val path = makeTempFile("fastq", ".fastq")
    val writer = FastqWriter(path)
    writer.write(rec)
    writer.close()
    path
  }

  it should "demux dual-indexed paired end reads" in {
    val fq1 = write(fq(name="frag", bases="AAAAAAAA"))
    val fq2 = write(fq(name="frag", bases="A"*100))
    val fq3 = write(fq(name="frag", bases="T"*100))
    val fq4 = write(fq(name="frag", bases="GATTACAGA"))
    val structures = Seq(ReadStructure("8B"), ReadStructure("100T"), ReadStructure("100T"), ReadStructure("9B"))

    val output: DirPath = outputDir()
    val metrics = makeTempFile("metrics", ".txt")
    
    new DemuxFastqs(inputs=Seq(fq1, fq2, fq3, fq4), output=output, metadata=sampleSheetPath,
      readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3).execute()

    val sampleBarcodMetrics = Metric.read[SampleBarcodeMetric](metrics)
    val sampleInfos = toSampleInfos(structures)
    val samples: Seq[Sample] = sampleInfos.map(_.sample)

    sampleBarcodMetrics.zip(samples).foreach { case (metric, sample) =>
      metric.barcode_name shouldBe sample.sampleName
      metric.barcode      shouldBe sample.sampleBarcodeString
      metric.library_name shouldBe sample.libraryId
      if (sample.sampleOrdinal == 1) {
        metric.templates                   shouldBe 1
        metric.pf_templates                shouldBe 1
        metric.perfect_matches         shouldBe 1
        metric.one_mismatch_matches    shouldBe 0
        metric.pf_perfect_matches      shouldBe 1
        metric.pf_one_mismatch_matches shouldBe 0
      }
      else if (sample.sampleId == UnmatchedSampleId) {
        metric.templates                   shouldBe 0
        metric.pf_templates                shouldBe 0
        metric.perfect_matches         shouldBe 0
        metric.one_mismatch_matches    shouldBe 0
        metric.pf_perfect_matches      shouldBe 0
        metric.pf_one_mismatch_matches shouldBe 0
      }
      else {
        metric.templates                   shouldBe 0
        metric.pf_templates                shouldBe 0
        metric.perfect_matches         shouldBe 0
        metric.one_mismatch_matches    shouldBe 0
        metric.pf_perfect_matches      shouldBe 0
        metric.pf_one_mismatch_matches shouldBe 0
      }
    }
  }

  it should "include all bases if the include-all-bases-in-fastqs option is used" in {
    val fq1 = write(fq(name="frag", bases="AAAAAAAA"))
    val fq2 = write(fq(name="frag", bases="A"*100))
    val fq3 = write(fq(name="frag", bases="T"*100))
    val fq4 = write(fq(name="frag", bases="GATTACAGA"))
    val structures = Seq(ReadStructure("8B"), ReadStructure("5M10S5M10S70T"), ReadStructure("10S90T"), ReadStructure("9B"))

    val output: DirPath = outputDir()
    val metrics = makeTempFile("metrics", ".txt")

    new DemuxFastqs(inputs=Seq(fq1, fq2, fq3, fq4), output=output, metadata=sampleSheetPath,
      readStructures=structures, metrics=Some(metrics), maxMismatches=2, minMismatchDelta=3,
      outputType=Some(OutputType.BamAndFastq), includeAllBasesInFastqs=true).execute()

    val sampleInfos         = toSampleInfos(structures)
    val matchedSampleInfo   = sampleInfos.find(!_.isUnmatched).get

    matchedSampleInfo.sample.sampleOrdinal shouldBe 1

    val prefix = toSampleOutputPrefix(matchedSampleInfo.sample, matchedSampleInfo.isUnmatched, false, output, UnmatchedSampleId)
    prefix.toString shouldBe output.resolve("20000101-EXPID-1-Sample_Name_1-AAAAAAAAGATTACAGA").toString

    val extensions = FastqRecordWriter.extensions(pairedEnd=true, illuminaStandards=false)
    def fastqRecord(which: Int): FastqRecord = {
      val path = PathUtil.pathTo(prefix + extensions(which-1))
      FastqSource(path).toSeq.headOption.value
    }

    val fastqR1    = fastqRecord(1)
    val fastqR2    = fastqRecord(2)
    val records    = readBamRecs(output.resolve(prefix + ".bam"))
    val recR1      = records.find(_.firstOfPair).value
    val recR2      = records.find(_.secondOfPair).value

    fastqR1.bases shouldBe "A"*100
    fastqR2.bases shouldBe "T"*100

    recR1.basesString shouldBe "A"*70
    recR2.basesString shouldBe "T"*90
  }

  "FastqRecordWriter.readName" should "set the read name based if it should follow Illumina standards" in {
    val output = outputDir()
    val baseRec = DemuxRecord(name="name", bases="", quals="", molecularBarcode=Some("MB"), sampleBarcode=Some("SB"), readNumber=1, pairedEnd=false, comment=None)

    {
      val rec    = baseRec.copy(pairedEnd=false, comment=None)
      val writer = new FastqRecordWriter(output.resolve("prefix"), pairedEnd=false, illuminaStandards=false)
      writer.readName(rec) shouldBe "name:SB:MB"
      writer.readName(rec.copy(molecularBarcode=None)) shouldBe "name:SB"
      writer.readName(rec.copy(sampleBarcode=None)) shouldBe "name:MB"
      writer.readName(rec.copy(molecularBarcode=None, sampleBarcode=None)) shouldBe "name"
    }

    {
      val rec    = baseRec.copy(pairedEnd=true, name="Instrument:RunID:FlowCellID:Lane:Tile:X:Y", comment=Some("ReadNum:FilterFlag:0:SampleNumber"))
      val writer = new FastqRecordWriter(output.resolve("prefix"), pairedEnd=true, illuminaStandards=true)
      writer.readName(rec) shouldBe "Instrument:RunID:FlowCellID:Lane:Tile:X:Y ReadNum:FilterFlag:0:SampleNumber"
    }
  }

  it should "fail if there was no comment in the given record when following Illumina standards" in {
    val output = outputDir()
    val baseRec = DemuxRecord(name="name", bases="", quals="", molecularBarcode=Some("MB"), sampleBarcode=Some("SB"), readNumber=1, pairedEnd=false, comment=None)

    {
      val rec    = baseRec.copy(pairedEnd=true, name="Instrument:RunID:FlowCellID:Lane:Tile:X:Y", comment=None)
      val writer = new FastqRecordWriter(output.resolve("prefix"), pairedEnd=true, illuminaStandards=true)
      an[Exception] should be thrownBy writer.readName(rec)
    }
  }

  it should "fail if the read name or comment does not follow Illumina standards" in {
    val output = outputDir()
    val baseRec = DemuxRecord(name="name", bases="", quals="", molecularBarcode=Some("MB"), sampleBarcode=Some("SB"), readNumber=1, pairedEnd=false, comment=None)

    {
      val rec    = baseRec.copy(pairedEnd=true, name="RunID:FlowCellID:Lane:Tile:X:Y", comment=Some("ReadNum:FilterFlag:0:SampleNumber"))
      val writer = new FastqRecordWriter(output.resolve("prefix"), pairedEnd=true, illuminaStandards=true)
      an[Exception] should be thrownBy writer.readName(rec)
    }

    {
      val rec    = baseRec.copy(pairedEnd=true, name="Instrument:RunID:FlowCellID:Lane:Tile:X:Y:Z", comment=Some("ReadNum:FilterFlag:0:SampleNumber"))
      val writer = new FastqRecordWriter(output.resolve("prefix"), pairedEnd=true, illuminaStandards=true)
      an[Exception] should be thrownBy writer.readName(rec)
    }

    {
      val rec    = baseRec.copy(pairedEnd=true, name="Instrument:RunID:FlowCellID:Lane:Tile:X:Y", comment=Some("FilterFlag:0:SampleNumber"))
      val writer = new FastqRecordWriter(output.resolve("prefix"), pairedEnd=true, illuminaStandards=true)
      an[Exception] should be thrownBy writer.readName(rec)
    }

    {
      val rec    = baseRec.copy(pairedEnd=true, name="Instrument:RunID:FlowCellID:Lane:Tile:X:Y", comment=Some("ReadNum:FilterFlag:0:SampleNumber:Z"))
      val writer = new FastqRecordWriter(output.resolve("prefix"), pairedEnd=true, illuminaStandards=true)
      an[Exception] should be thrownBy writer.readName(rec)
    }
  }
}
