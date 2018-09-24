#!/usr/bin/env python

import argparse
import os
import sys
from pyfaidx import Fasta
import subprocess
import bisect
import random


class variant():
    """
    General class to store a somatic variant
    """

    def __init__(self, chrom, position, slices, alt=None, otherChrom=None, otherPos=None,
                 copyNum=2, phase=None):
        self.chrom = chrom
        self.position = position
        self.slices = slices
        self.alt = alt
        self.otherChrom = otherChrom
        self.otherPos = otherPos
        self.copyNum = copyNum
        self.phase = phase

def loadStructVar(inFile, coverage, sliceSize, ploidy):

    if inFile is None:
        # No SV file was specified
        return []

    # So the highest coverage region corresponding to each variant is identified first, iterate over the dictionary
    # from highest to lowest coverage
    hiToLowKeys = {}
    for chrom, coverageInfo in coverage.items():
        hiToLowKeys[chrom] = tuple(sorted(coverageInfo.keys(), reverse=True))

    # For assigning a phase to variants later
    cnStates = set(range(1, ploidy + 1))

    svs = []
    # Read in Structural variants
    with open(inFile) as f:
        for line in f:
            if line.startswith("#"):  # i.e. header/comment line. Ignore it
                continue
            # The format of the input line should be:
            # "br1Chrom br1Pos  br2Chrom    br2Pos  Cellularity"
            cols = line.split()
            try:
                chrom1 = cols[0]
                pos1 = int(cols[1])
                chrom2 = cols[2]
                pos2 = int(cols[3])
                cellularity = float(cols[4])
                # Sanity check
                if cellularity > 1:
                    raise ValueError(
                        "In \'%s\', line \'%s\', cellularity cannot be greater than 1" % (inFile, line.strip("\n").strip("\r")))
            except IndexError:
                raise MalformedInputException("We were unable to process line \'%s\' in \'%s\', as it is missing one or more columns" % (line.strip("\n").strip("\r"), inFile))

            # Create two breakpoint objects for this event, since the different breakpoints
            # will have different coverages
            slice1, phase1 = assignSliceAndPhase(chrom1, pos1, coverage, hiToLowKeys, ploidy, sliceSize, cnStates, cellularity)
            slice2, phase2 = assignSliceAndPhase(chrom2, pos2, coverage, hiToLowKeys, ploidy, sliceSize, cnStates, cellularity)

            # Take whichever state is the most constricted (i.e. has the highest copy number state)
            if len(phase1) > len(phase2):
                # Create an SV variant for each breakpoint
                var1 = variant(chrom1, pos1, slice1, otherChrom=chrom2, otherPos=pos2, phase=phase1)
                svs.append(var1)
                var2 = variant(chrom2, pos2, slice1, otherChrom=chrom1, otherPos=pos1, phase=phase1)
                svs.append(var2)
            else:
                # Create an SV variant for each breakpoint
                var1 = variant(chrom1, pos1, slice2, otherChrom=chrom2, otherPos=pos2, phase=phase2)
                svs.append(var1)
                var2 = variant(chrom2, pos2, slice2, otherChrom=chrom1, otherPos=pos1, phase=phase2)
                svs.append(var2)

    return tuple(svs)


def loadCopyNumberVar(inFile, coverage, sliceSize, ploidy=2, verbose=False):
    """
    Read CNVs, and calculate the appropriate slice levels for the cnvs
    :param inFile:
    :param coverage:
    :param sliceSize:
    :param ploidy:
    :return:
    """

    if inFile is None:  # i.e. there are no CNVs
        return {}

    cnvs = {}
    skipVar = 0
    totalVar = 0

    # So the highest coverage region corresponding to each variant is identified first, iterate over the dictionary
    # from highest to lowest coverage
    hiToLowKeys = {}
    for chrom, coverageInfo in coverage.items():
        hiToLowKeys[chrom] = tuple(sorted(coverageInfo.keys(), reverse=True))

    # Read in CNVs
    cnStates = set(range(1, ploidy + 1))
    with open(inFile) as f:
        for line in f:
            if line.startswith("#"):  # i.e. header/comment line. Ignore it
                continue
            # The file format should be "chrom  start   end CopyNum Cellularity
            try:
                cols = line.split("\t")
                chrom = cols[0]
                start = int(cols[1]) - 1  # Offset 1 for 0-based indexing, as used here
                end = int(cols[2]) - 1
                cn = int(cols[3])
                cellularity = float(cols[4])
                # Sanity check
                if cellularity > 1:
                    raise ValueError(
                        "In \'%s\', line \'%s\', cellularity cannot be greater than 1" % (inFile, line.strip("\n").strip("\r")))
            except IndexError:
                raise MalformedInputException(
                    "We were unable to process \'%s\' in \'%s\', as it is missing one or more columns" % (
                    line.strip("\n").strip("\r"), inFile))

            # Sanity check that the end of the segment is after the start of the segment
            if start >= end:
                raise MalformedInputException("When processing \'%s\' in \'%s\', the CN start position is greater than the end position" % (line.strip("\n").strip("\r"), inFile))
            if cn == ploidy:  # i.e. there is no copy number change. Ignore this event
                skipVar += 1
                continue

            # Find the highest depth slice which overlaps this event
            for covLevel in hiToLowKeys[chrom]:
                regions = coverage[chrom][covLevel]
                startPoint = bisect.bisect(regions, start)
                endPoint = bisect.bisect(regions, end)
                if startPoint != endPoint or startPoint % 2 == 1:  # An odd number means that it falls within a region. Different bisect points means it spans different regions
                    eventCov = covLevel
                    break  # Since we are going from highest to lowest coverage, if a slice level has been found, it is the max

            # How many reads are needed to support the clonality of this event
            readNum = int(eventCov * cellularity)
            # Which slices need to support this event to give it the desired effect?
            coverages = set(range(sliceSize, eventCov + sliceSize, sliceSize))
            slices = set(random.sample(coverages, round(readNum/sliceSize)))

            # Which "phases" will this event effect?
            alt = cn - ploidy
            if alt < 0:  # Deletion
                phase = set(random.sample(cnStates, alt * -1))
            else:  # Amplification. Assume this always influences the same copy
                phase = set(random.sample(cnStates, 1))
            # Create an variant object for this event
            cnv = variant(chrom, start, slices, otherPos=end, copyNum=cn, phase=phase, alt=alt)

            # Save this event
            if chrom not in cnvs:
                cnvs[chrom] = {}
            for p in phase:
                if p not in cnvs[chrom]:
                    cnvs[chrom][p] = {}
                for slice in slices:
                    if slice not in cnvs[chrom][p]:
                        cnvs[chrom][p][slice] = []
                    cnvs[chrom][p][slice].append(cnv)
            totalVar += 1

    if verbose:
        sys.stderr.write("Sucessfully loaded %s CNVs. %s copy-neutral segments were ignored" % (totalVar, skipVar) + os.linesep)
    return cnvs


def assignSliceAndPhase(chrom, pos, coverage, coverageLevels, ploidy, sliceSize, cnStates, vaf):
    """
    Document later
    :return:
    """
    # Obtain the coverage for this chromosome
    try:
        chromCov = coverage[chrom]
    except KeyError:  # i.e. this contig has no coverage. Ignore this variant
        return {}, {}

    varTotalCov = None
    # Find the coverage for this variant
    for covLevel in coverageLevels[chrom]:
        regions = coverage[chrom][covLevel]
        insertPoint = bisect.bisect(regions, pos)
        if insertPoint % 2 == 1:  # An odd number means that it falls within a region
            varTotalCov = covLevel
            break  # Since we are going from the deepest to shallowest coverage, we don't need to search anymore

    if varTotalCov is None:  # i.e. this variant has no coverage. Ignore it
        return {}, {}

    # How many reads need to support this variant to reach the desired VAF?
    varReadCount = int(vaf * varTotalCov) * ploidy
    # Calculate which slices should support this variant
    # We randomly select slices to minimize biases in the reads which support a variant
    coverages = set(range(sliceSize, varTotalCov + sliceSize, sliceSize))

    copyNum = 0
    findingSlices = True
    # Determine which slice intervals should support this variant
    while findingSlices:
        try:
            copyNum += 1
            slices = set(random.sample(coverages, int(varReadCount / sliceSize / copyNum)))  # This may throw ValueError
            findingSlices = False
        except ValueError:
            # If we end up in here, it means that, given the current copy number number state
            # assigned to this variant, it is impossible to meet the VAF specified by the user
            # In this case, we need to assign a higher copy number state to this variant, and try again

            # If the current copy number state matches the ploidy of this region, then it is impossible to
            # introduce this variant
            if copyNum == ploidy:
                raise ValueError("Unable to process SNV \'%s\', as it is impossible to meet the given VAF specified" % (
                    line.strip("\n").strip("\r")))
            pass

    # Which copy of the genome is this event in?
    phase = set(random.sample(cnStates, copyNum))

    return slices, phase


def loadSSMs(inFile, coverage, sliceSize, ploidy, verbose=False):
    """
    Read SNVs, and determine the appropriate slice level for the snv
    :param inFile: A string containing a filepath to a tab-deliniated file containing SNVs
    :param coverage: A dictionary containing a coverage profile for the synthetic genome.
    :param sliceSize: An int listing the coverage for any given slice
    :return: 
    
    The format of "coverage" is {chrom: {coverageLevel: [start, end, start2, end2, ...]...}...}
    """

    if inFile is None:  # i.e. the user did not specify a SNV file
        return {}

    skipVar = 0
    totalVar = 0
    snvs = {}  # Will store {chromosome: {phase: {sliceLevel: [var1, var2, var3,...]}}}

    # Used to assign a variant a specific "phase" later
    cnStates = set(range(1, ploidy+1))

    # So the highest coverage region corresponding to each variant is identified first, iterate over the dictionary
    # from highest to lowest coverage
    hiToLowKeys = {}
    for chrom, coverageInfo in coverage.items():
        hiToLowKeys[chrom] = tuple(sorted(coverageInfo.keys(), reverse=True))

    # Read in variants
    with open(inFile) as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip header lines
            # Assume the file format is "chrom  position    altAllele   VAF"
            try:
                cols = line.split("\t")
                chrom = cols[0]
                pos = int(cols[1]) - 1  # Subtract 1 to account for 0-based indexing
                altAllele = cols[2]
                vaf = float(cols[3])

                # Sanity check
                if vaf > 1:
                    raise ValueError("In \'%s\', line \'%s\', VAF cannot be greater than 1" % (inFile, line.strip("\n").strip("\r")))
            except IndexError:  # i.e. we are missing one or more columns
                raise MalformedInputException(
                    "We were unable to process \'%s\' in \'%s\', as it is missing one or more columns" % (
                    line.strip("\n").strip("\r"), inFile))

            slices, phase = assignSliceAndPhase(chrom, pos, coverage, hiToLowKeys, ploidy, sliceSize, cnStates, vaf)
            if slices is None:  # i.e. this variant has no coverage
                skipVar += 1
                continue

            # Create this variant
            var = variant(chrom, pos, slices, alt=altAllele, phase=phase)

            # Save this variant
            if chrom not in snvs:
                snvs[chrom] = {}
            for p in phase:
                if p not in snvs[chrom]:
                    snvs[chrom][p] = {}
                for sliceLevel in slices:
                    if sliceLevel not in snvs[chrom][p]:
                        snvs[chrom][p][sliceLevel] = []
                    snvs[chrom][p][sliceLevel].append(var)
            totalVar += 1

    # To speedup downstream processing, convert the variable lists to tuples
    for chrom in snvs.keys():
        for phase in snvs[chrom].keys():
            for sliceLevel in snvs[chrom][phase].keys():
                snvs[chrom][phase][sliceLevel] = tuple(snvs[chrom][phase][sliceLevel])

    if verbose:
        sys.stderr.write("Sucessfully loaded %s SNVs. %s SNVs were ignored" % (totalVar, skipVar) + os.linesep)
    return snvs


def loadSlices(inFile):
    """
    Load slice coverage profile
    
    Create a dictionary, listing all regions with higher than x coverage

    :param inFile: A string containing a filepath to the slice coverage file 
    :return: 
    """

    sliceRegions = {}
    coverageProfile = {}
    sliceIntervals = []

    with open(inFile) as f:
        for line in f:

            if line.startswith("#"):
                continue  # Skip header lines
            chrom, start, end, sliceCov = line.split("\t")

            sliceCov = int(sliceCov)

            if chrom not in coverageProfile:
                coverageProfile[chrom] = {}
            if sliceCov not in coverageProfile[chrom]:
                coverageProfile[chrom][sliceCov] = []
            if sliceCov not in sliceRegions:
                sliceRegions[sliceCov] = {}
                # For calculating coverage at a given variant later

            if chrom not in sliceRegions[sliceCov]:
                sliceRegions[sliceCov][chrom] = []

            # Add the coordinates of this slice
            sliceRegions[sliceCov][chrom].append((int(start), int(end)))

            # So we can determine the maximum coverage withing a given region later (i.e when adding mutations),
            # also save this region and coverage in a data structure which can be searched easily (and quickly) later
            coverageProfile[chrom][sliceCov].extend((int(start), int(end)))

            # Store the slice coverage so we can calculate the distance between slices
            sliceIntervals.append(sliceCov)

    # Calculate the slice coverage
    # Assuming the input file is unmodified, the mimumum coverage level should be the slice interval
    sliceCov = min(sliceIntervals)

    # As a sanity check, confirm that all coverage levels are a multiple of this
    for sliceLevel in sliceIntervals:
        if sliceLevel % sliceCov != 0:
            raise AttributeError("The input file \'%s\' appears to use slice levels of \'%s\', yet \'%s\' is not a multiple of the slice level" % (inFile, sliceCov, sliceLevel))

    return sliceRegions, sliceCov, coverageProfile


def isValidFile(file, parser):
    """
    Checks to determine if the input file exists, and throws an error if it does not
    """

    if not os.path.exists(file):
        raise parser.error("Unable to locate \'%s\': No such file or directory" % file)
    else:
        return file


def generateFasta(refGenome, sliceSegments, snvs, cnvs, svs, outName, sliceDepth, readLength, copyNum):
    """
    Generates a FASTA entry for each region inside region, incorperating the specified variants
    """

    offset = 60

    with open(outName, "a") as o:

        regions = sliceSegments[sliceDepth]

        # Generate a reference sequence for each region
        for chrom, coordinates in regions.items():
            for start, end in coordinates:

                # Store deletions which overlap this segment for later processing
                deletions = []
                # Obtain the FASTA sequence
                # Since ART doesn't generate exactly linear coverage across the specified region, add some
                # flanking bases
                start -= offset
                end += offset
                seq = list(refGenome[chrom][start:end].seq)

                # Process snvs
                try:
                    for snv in snvs[chrom][copyNum][sliceDepth]:
                        if start <= snv.position < end:  # This variant is covered by this reference
                            snvPos = snv.position - start
                            seq[snvPos] = snv.alt
                except KeyError:  #i.e. there are no variants on this chromosome, at this slice level, or at this phase
                    pass

                # Process copy number events
                removeRegion = False
                try:
                    for cnv in cnvs[chrom][copyNum][sliceDepth]:
                        # Determine if this event overlaps this region
                        if cnv.position < end and cnv.otherPos > start:  # cnv.position is start of the event, and vise-versa
                            # Is this an amplification or deletion?
                            if cnv.alt > 0:  # Amplification
                                # Create additional copy/copies of the reference sequence in the region which overlaps
                                # this amplification
                                if cnv.position > start:
                                    eventStart = cnv.position - start  # This is relative to the sequence already prepared for this region
                                else:
                                    eventStart = 0
                                if cnv.otherPos > end:
                                    eventEnd = end - start  # Relative to the sequence already parsed for this region.
                                else:
                                    eventEnd = cnv.otherPos - start
                                eventSeq = seq[eventStart:eventEnd]

                                # Create extra copies
                                for j in range(1, cnv.alt + 1):

                                    # Create a name for this event
                                    eventName = ">" + chrom + "-" + str(eventStart) + "-" + str(eventEnd) + "-" + str(sliceDepth) + "-" + str(copyNum) + "-AMP" + str(j)

                                    # Write out this event
                                    o.write(eventName + os.linesep)
                                    o.write("".join(eventSeq) + os.linesep)
                            else:  # i.e. this is a deletion

                                # Does this deletion completely encompass this region?
                                if cnv.position < start and cnv.otherPos > end:
                                    # If so, we can simply not make a reference for this region, thus creating a deletion
                                    removeRegion = True
                                    break
                                else:
                                    # Store this deletion for later
                                    eventStart = cnv.position - start
                                    if eventStart < 0:
                                        eventStart = 0
                                    eventEnd = cnv.otherPos - start
                                    if eventEnd > end - start:
                                        eventEnd = end - start
                                    deletions.append((eventStart, eventEnd))
                except KeyError:  # i.e. no events on this chromosome, phase, or slice level
                    pass

                if removeRegion:  # i.e. a deletion spans this entire region. Don't generate it
                    continue

                splitRegions = []
                if not deletions:
                    # No deletions overlap this region. Use the unaltered region
                    splitRegions.append((0, end - start))
                else:
                    # We have to split the reference into multiple sections

                    # First, sort the deleted regions
                    deletions.sort(key=lambda x:x[0])
                    # Next, remove any deletions which overlap one another
                    lastDelEnd = 0
                    i = 0
                    try:
                        while True:
                            currentDel = deletions[i]
                            if currentDel[0] < lastDelEnd:
                                # This deletion overlaps the previous deletion. Ignore it
                                # in theory, this could be handeled a lot better
                                # i.e. look at if a larger deletion overlaps many smaller deletions etc
                                # But this should work
                                deletions.pop(i)
                            else:
                                lastDelEnd = currentDel[1]
                                i += 1
                    except IndexError:
                        # We have processed all deletions
                        pass

                    # Finally, create the subdivided regions for this region
                    x = 0
                    assert len(deletions) != 0
                    for delStart, delEnd in deletions:
                        if delStart > x:
                            splitRegions.append((x, delStart))
                        x = delEnd
                    if delEnd < end - start:
                        splitRegions.append((x, end - start))

                for rStart, rEnd in splitRegions:
                    rSeq = seq[rStart: rEnd]
                    # Restore the original coordinates of this slice
                    rStart += start
                    rEnd += end

                    # Generate a unique name for this entry
                    fastaName = ">" + chrom + "-" + str(rStart) + "-" + str(rEnd) + "-" + str(sliceDepth) + "-" + str(
                        copyNum)

                    # Process SVs
                    for sv in svs:
                        if sv.chrom != chrom:
                            continue  # This structural variant falls on a different chromosome
                        # Are we even processing this SV in this slice?
                        elif sliceDepth not in sv.slices:
                            continue  # Don't process this event now
                        elif copyNum not in sv.phase:
                            continue
                        elif rStart <= sv.position <= rEnd:  # This breakpoint falls within the region
                            otherSeq = refGenome[sv.otherChrom][sv.otherPos - readLength: sv.otherPos + readLength].seq
                            otherBreakPoint = readLength
                            # Create both breakpoints
                            svSeq1 = "".join(seq[:sv.position - start]) + otherSeq[otherBreakPoint:]
                            svSeq2 = otherSeq[otherBreakPoint:] + "".join(seq[:sv.position - start])
                            # We can add both sequences to the FASTA, since each sequence only has half of this region
                            svName1 = fastaName + "-BR1-" + str(sv.otherChrom) + "-" + str(sv.otherPos)
                            svName2 = fastaName + "-BR2-" + str(sv.otherChrom) + "-" + str(sv.otherPos)
                            o.write(svName1 + os.linesep)
                            o.write(svSeq1 + os.linesep)
                            o.write(svName2 + os.linesep)
                            o.write(svSeq2 + os.linesep)

                            # To avoid double-counting the breakpoint sequences, add a deletion
                            delStart = sv.position - readLength - start
                            if delStart < 0:
                                delStart = 0
                            delEnd = sv.position + readLength - start
                            if delEnd > end - start:
                                delEnd = end - start
                            deletions.append((delStart, delEnd))

                    # Generate one FASTA for each region.
                    o.write(fastaName + os.linesep)
                    o.write("".join(rSeq) + os.linesep)


def runArt(artCom):
    # Run ART
    artStdErr = []
    artStdOut = []
    artRun = subprocess.Popen(artCom, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    for oLine in artRun.stdout:
        oLine = oLine.decode("utf-8")
        artStdOut.append(oLine)
    for eLine in artRun.stderr:
        eLine = eLine.decode("utf-8")
        artStdErr.append(eLine)

    artRun.stdout.close()
    artRun.stderr.close()
    artRun.wait()

    if artRun.returncode != 0:
        sys.stderr.write(os.linesep.join(artStdErr))
        sys.stdout.write(os.linesep.join(artStdOut))


def createOfftargetRef(refGenome, outFile, cnvs, ploidy, minSlice):

    # Create a modified reference with all FASTA entries merged into one
    firstHeader = True
    with open(outFile, "w") as o:
        # Create a sequence for each copy of the genome
        o.write(">OfftargetReads" + os.linesep)
        for copyNum in range(1, ploidy + 1):
            for chrom in refGenome.faidx.index.keys():
                # Since copy number alterations also affect off-target reads, and some tools use off-target reads to identify
                # CNAs, we should simulate those as well

                # First, find all copy number events on this chromosome and ploidy
                # To handle sub-clonal events, I will sample the lowest slice level to (somewhat poorly) simulate the
                # likelihood that off-target reads for a region are obtained from the cells which contain such an event
                try:
                    cnaEvents = cnvs[chrom][copyNum][minSlice]
                except KeyError:  # There are no copy number events to worry about
                    cnaEvents = []

                cnaEvents.sort(key=lambda x: x.position)

                # Process events
                refPos = 0  # We have processed up to this position so far
                refLength = refGenome.faidx.index[chrom].rlen
                for event in cnaEvents:
                    cnChange = event.alt

                    # Sanity check that this event does not stretch past the end of the chromosome
                    if event.otherPos > refLength:
                        eventEnd = refLength
                    else:
                        eventEnd = event.otherPos

                    if cnChange < 0:  # i.e. deletion

                        # If this event was ecliped by a previous event, don't include it
                        if refPos > event.otherPos:
                            continue

                        # Don't include the reference sequence corresponding to this event
                        o.write(refGenome[chrom][refPos:event.position].seq)
                        o.write(os.linesep)
                        refPos = eventEnd
                    else:  # Amplification
                        # Write out the sequence up to this point, unless that has already been done
                        if refPos < eventEnd:
                            o.write(refGenome[chrom][refPos:eventEnd].seq)
                            o.write(os.linesep)
                            refPos = eventEnd
                        # Duplicate this sequence
                        # In essence, we are simulating a tandem duplication
                        for i in range(0, cnChange):
                            o.write(refGenome[chrom][event.position:eventEnd].seq)
                            o.write(os.linesep)

                # Finally, process the rest of the sequence
                o.write(refGenome[chrom][refPos:refLength].seq)
                o.write(os.linesep)


def createDup(fastq1, dupRate, fastq2=None):
    """
    Create PCR and optical duplicates

    Read through each FASTQ file. Create copies of a given read (or read pair) depending upon the duplicate rate

    Note that duplicates will have the exact same sequence and quality scores as the original reads

    :param fastq1:
    :param dupRate:
    :param fastq2:
    :return:
    """

    # Create tmp output files for duplicate FASTQs
    dupFastq1 = fastq1 + "dup"
    if fastq2 is not None:
        dupFastq2 = fastq2 + "dup"

    with open(dupFastq1, "w") as o1, open(dupFastq2 if fastq2 is not None else os.devnull, "w") as o2:
        with open(fastq1) as f1, open(fastq2) if fastq2 is not None else None as f2:
            while True:
                # Read a FASTQ entry
                name1 = f1.readline()
                if name1 == "":  # i.e. we have reached the end of the file
                    break
                seq1 = f1.readline()
                strand1 = f1.readline()
                qual1 = f1.readline()

                if fastq2 is not None:  # i.e. paired-end
                    # Read a FASTQ entry
                    name2 = f2.readline()
                    seq2 = f2.readline()
                    strand2 = f2.readline()
                    qual2 = f2.readline()

                i = 0
                # Should we create a duplicate for this read?
                while random.random() < dupRate:
                    i += 1
                    # Change the name of the duplicate to avoid read pairs which share the same name
                    oName1 = name1.replace("@", "@Dup" + str(i) + "-")
                    o1.write(oName1)
                    o1.write(seq1)
                    o1.write(strand1)
                    o1.write(qual1)

                    if fastq2 is not None:
                        oName2 = name2.replace("@", "@Dup" + str(i) + "-")
                        o2.write(oName2)
                        o2.write(seq2)
                        o2.write(strand2)
                        o2.write(qual2)


    # Finally, add the duplicate reads to the input FASTQ files
    with open(fastq1, "a") as f1, open(dupFastq1) as r1:
        for line in r1:
            f1.write(line)

    os.remove(dupFastq1)
    if fastq2 is not None:
        with open(fastq2, "a") as f2, open(dupFastq2) as r2:
            for line in r2:
                f2.write(line)
        os.remove(dupFastq2)

def getArgs():
    """
    Processes command line arguments
    """

    parser = argparse.ArgumentParser(description="Generates synthetic FASTQ files with the following characteristics and features")
    inputArgs = parser.add_argument_group("Input Arguments")
    inputArgs.add_argument("-i", "--input", required=True, metavar="TSV",
                           type=lambda x: isValidFile(x, parser),
                           help="Input slice coverage file")
    inputArgs.add_argument("-r", "--reference", required=True, metavar="FASTA",
                           type=lambda x: isValidFile(x, parser), help="Reference genome, in FASTA format")

    mutArgs = parser.add_argument_group("Germline and somatic mutations")
    mutArgs.add_argument("--snv", required=False, type=lambda x: isValidFile(x, parser),
                         help="A tab-deliniated file listing the positions of single nucleotide variants")
    mutArgs.add_argument("--cnv", required=False, type=lambda x: isValidFile(x, parser),
                         help="A tab-deliniated file listing the position of all copy number variants")
    mutArgs.add_argument("--str", required=False, type=lambda x: isValidFile(x, parser),
                         help="A tab-deliniated file listing the positions of all structural variants")

    outputArgs = parser.add_argument_group("Output arguments")
    outputArgs.add_argument("-o", "--outdir", required=True, help="Output directory")
    outputArgs.add_argument("--samples", default=1, type=int, help="Number of samples to generate [default: %(default)s]")
    outputArgs.add_argument("-p", "--prefix", default="LBS",
                            help="Prefix to prepend to each output file [default: %(default)s]")

    miscArgs = parser.add_argument_group("Tool paths and other miscellaneous arguments")
    miscArgs.add_argument("--tmpdir", default="." + os.sep, help="Path for temporary files")
    miscArgs.add_argument("--keep_tmp", action="store_true",
                          help="Don't delete intermediate files")
    miscArgs.add_argument("--no_paired", action="store_true",
                          help="Don't generate paired-end data")
    miscArgs.add_argument("--art_path", required=True, type=lambda x: isValidFile(x, parser),
                          help="Path to \'art_illumina\', the read simulator")
    miscArgs.add_argument("--max_coverage", default=-1, type=int,
                          help="Maximum output fold coverage. Set to a negative number to disable [default: %(default)s]")
    miscArgs.add_argument("--capture_efficiency", metavar="FLOAT", default=0.90, type=float, help="Overall capture efficiency [default: %(default)s]")
    miscArgs.add_argument("--duplicate_rate", metavar="FLOAT", default=0.10, type=float, help="Read duplicate rate [default: %(default)s]")
    miscArgs.add_argument("--ploidy", default=2, type=int, help="Baseline genome ploidy [default: %(default)s]")
    miscArgs.add_argument("--verbose", action="store_true", help="Print debugging information")
    simulationArgs = parser.add_argument_group("Parameters passed to ART")
    simulationArgs.add_argument("--frag_length", default=166, type=int,
                                help="Median insert size [default: %(default)s]")
    simulationArgs.add_argument("--frag_length_stdev", default=20, type=int,
                                help="Insert size standard deviation [default: %(default)s]")
    simulationArgs.add_argument("--read_length", default=135, type=int,
                                help="Read Length [default: %(default)s]")
    simulationArgs.add_argument("-q1", "--quality1", metavar="TXT", required=True, type=lambda x: isValidFile(x, parser),
                                help="Quality profile for read1. Generated using \'art_profiler_illumina\'")
    simulationArgs.add_argument("-q2", "--quality2", metavar="TXT", type=lambda x: isValidFile(x, parser),
                                help="Quality profile for read2. Generated using \'art_profiler_illumina\'")

    args = parser.parse_args()

    # Sanity check args
    if args.capture_efficiency <= 0 or args.capture_efficiency > 1:
        raise parser.error("--capture_efficiency must be between 0 and 1")
    # Check paired-end args
    if not args.no_paired and not args.quality2:
        raise parser.error("-q2/--quality2 must be specified for paired-end simulation")

    # Standardize path
    args.outdir = os.path.abspath(args.outdir)
    args.tmpdir = os.path.abspath(args.tmpdir)

    if args.tmpdir is None:
        args.tmpdir = args.outdir
    return args


def main(args=None):
    if args is None:
        args = getArgs()

    # Prepare FASTA
    refGenome = Fasta(args.reference, one_based_attributes=False)

    # Store tmp FASTQ files
    fastq1 = {}
    fastq2 = {}
    for i in range(1, args.samples+1):
        fastq1[i] = []
        fastq2[i] = []

    # Load coverage slices
    sliceCov, sliceSize, coverageProfile = loadSlices(args.input)

    # Load mutations
    sys.stderr.write("Loading specified variants" + os.linesep)
    snvs = loadSSMs(args.snv, coverageProfile, sliceSize, args.ploidy, args.verbose)
    cnvs = loadCopyNumberVar(args.cnv, coverageProfile, sliceSize, args.ploidy, args.verbose)
    svs = loadStructVar(args.str, coverageProfile, sliceSize, args.ploidy)

    # Since we need to add a flank to each slice to ensure that sufficient coverage is generated in the regions
    # we care about (i.e. the actual target regions), we need to down-scale the target coverage so it is more
    # representative of the input profile
    effectiveCoverage = sliceSize * ( 1 - 60.0 / float(args.read_length))

    # Since we are generating multiple copies of each reference region (representing haplotypes), we also need to reduce
    # coverage to account for ploidy
    effectiveCoverage = effectiveCoverage / args.ploidy

    # ART does not generate uniform coverage across the specified region. This inconsistency will propagate the more
    # copies of each region which are generated (i.e. higher ploidy)
    # To compensate, we need to increase coverage slightly
    # Use a scaling factor of 9%, based upon some ad-hoc examinations
    effectiveCoverage = round(effectiveCoverage * 1.09)

    currentSlice = sliceSize

    sys.stderr.write("Generating synthetic reads" + os.linesep)
    while args.max_coverage < 0 or currentSlice <= args.max_coverage:

        if args.verbose:
            sys.stderr.write("Processing slice %s" % currentSlice  + os.linesep)
        # Parse the FASTA sequences which overlap this slice, and add variants
        tmpFasta = args.tmpdir + os.sep + str(currentSlice) + ".fa"
        if os.path.exists(tmpFasta):
            os.remove(tmpFasta)
        try:
            for currentCopy in range(1, args.ploidy + 1):
                generateFasta(refGenome, sliceCov, snvs, cnvs, svs, tmpFasta, currentSlice,
                              args.read_length, currentCopy)
        except KeyError as e:  # i.e. there are no regions which are covered by a slice this deep.
            break  # We have finished simulating on-target reads

        # Now, run ART to simulate reads for this region
        for sampleNum in range(1, args.samples + 1):
            fastqPrefixName = args.tmpdir + os.sep + str(currentSlice) + "_" + str(sampleNum) + "_"
            artCom = [args.art_path, "-1", args.quality1, "--fcov", str(effectiveCoverage), "--in",
                      tmpFasta, "-na", "--len", str(args.read_length), "-o", fastqPrefixName]
            if not args.no_paired:  # i.e. paired end simulation
                artCom.extend(("--paired", "--mflen", str(args.frag_length), "--sdev", str(args.frag_length_stdev),
                               "-2", args.quality2))
            if args.verbose:
                sys.stderr.write("Generating reads for sample %s" % sampleNum + os.linesep)
                sys.stderr.write("Running ART command " + " ".join(artCom) + os.linesep)
            runArt(artCom)

            # track the output files for merging later
            fastq1[sampleNum].append(fastqPrefixName + "1.fq")
            fastq2[sampleNum].append(fastqPrefixName + "2.fq")
        currentSlice += sliceSize

        # Remove the tmp FASTA file
        if not args.keep_tmp:
            os.remove(tmpFasta)

    sys.stderr.write("Merging temporary files" + os.linesep)
    # Finally, merge the FASTQ files
    outFastq1 = {}
    outFastq2 = {}
    for sampleNum in range(1, args.samples + 1):
        outFastq1[sampleNum] = args.outdir + os.sep + args.prefix + "." + str(sampleNum) + ".1.fastq"
        outFastq2[sampleNum] = args.outdir + os.sep + args.prefix + "." + str(sampleNum) + ".2.fastq"
        lines = 0
        with open(outFastq1[sampleNum], "w") as o:
            for f1 in fastq1[sampleNum]:
                with open(f1) as f:
                    for line in f:
                        lines += 1
                        o.write(line)

        # Delete tmp files
        if not args.keep_tmp:
            for f1 in fastq1[sampleNum]:
                os.remove(f1)

        if not args.no_paired:
            with open(outFastq2[sampleNum], "w") as o:
                for f2 in fastq2[sampleNum]:
                    with open(f2) as f:
                        for line in f:
                            o.write(line)

            if not args.keep_tmp:
                for f2 in fastq2[sampleNum]:
                    os.remove(f2)

    # Now, generate off-target reads
    sys.stderr.write("Generating off-target reads" + os.linesep)
    # Note: We are using the counts from the last set of FASTQs specified
    onTargetCount = lines/4

    if not args.no_paired:
        onTargetCount = onTargetCount * 2

    offTargetCount = int(float(onTargetCount) / float(args.capture_efficiency) - float(onTargetCount))

    if args.verbose:
        sys.stderr.write("Generated %s on-target reads" % onTargetCount + os.linesep)
        sys.stderr.write("Generating %s off-target reads" % offTargetCount  + os.linesep)

    tmpRef = args.tmpdir + os.path.sep + os.path.basename(args.reference) + ".tmpOfftarget.fa"
    createOfftargetRef(refGenome, tmpRef, cnvs, args.ploidy, sliceSize)

    # Actually generate off-target reads
    for sampleNum in range(1, args.samples + 1):
        fastqPrefixName = args.tmpdir + os.sep + "offtarget." + str(sampleNum) + "."

        artCom = [args.art_path, "-1", args.quality1, "-c", str(offTargetCount), "--in",
                  tmpRef, "-na", "--len", str(args.read_length), "-o", fastqPrefixName]
        if not args.no_paired:  # i.e. paired end simulation
            artCom.extend(("--paired", "--mflen", str(args.frag_length), "--sdev", str(args.frag_length_stdev),
                           "-2", args.quality2))
        runArt(artCom)

        # Finally, merge the on-target and off-target FASTQ files
        offFastq1 = fastqPrefixName + "1.fq"
        with open(outFastq1[sampleNum], "a") as o:
            with open(offFastq1) as f:
                for line in f:
                    lines += 1
                    o.write(line)

        # Remove the off-target FASTQ file
        if not args.keep_tmp:
            os.remove(offFastq1)

        if not args.no_paired:
            offFastq2 = fastqPrefixName + "2.fq"
            with open(outFastq2[sampleNum], "a") as o:
                with open(offFastq2) as f:
                    for line in f:
                        o.write(line)

            if not args.keep_tmp:
                os.remove(offFastq2)

    # Create PCR duplicates
    # TODO: Merge with above section to improve performance
    sys.stderr.write("Creating PCR Duplicates" + os.linesep)
    for sampleNum in range(1, args.samples + 1):
        createDup(outFastq1[sampleNum], args.duplicate_rate, outFastq2[sampleNum] if not args.no_paired else None)


# Error classes
class MalformedInputException(Exception):
    pass

if __name__ == "__main__":
    main()

