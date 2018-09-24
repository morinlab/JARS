#!/usr/bin/env python

import argparse
import os
import pysam
import sys


class segment():

    def __init__(self, chrom, start):
        self.chrom = chrom
        self.start = start
        self.end = None


def isValidFile(file, parser):
    """
    Checks to determine if the specified file path is valid
    
    :param file: A string containing a filepath 
    :param parser: An argparse.ArgumentParser() object
    :return: file
    :raises: argparse.ArgumentParser.error() if the file does not exist
    """

    if not os.path.exists(file):
        raise parser.error("Unable to locate \'%s\': No such file or directory" % file)
    else:
        return file


def addSlice(finishedSlices, currentSlice, sliceLevel):
    """
    
    :param finishedSlices: 
    :param currentSlice: 
    :return: 
    """

    # Store this slice
    if sliceLevel not in finishedSlices:
        finishedSlices[sliceLevel] = [currentSlice]
    else:
        # Check to see if there is another fragment at this slice depth which overlaps this fragment
        previousRegion = finishedSlices[sliceLevel][-1]
        if previousRegion.end > currentSlice.start and previousRegion.chrom == currentSlice.chrom:
            # These slices overlap
            # Don't create a new region. Instead, extend the previous region
            previousRegion.end = currentSlice.end
        else:
            finishedSlices[sliceLevel].append(currentSlice)


def getArgs():
    """
    Process command line arguments

    :return: An arguments namespace
    """

    parser = argparse.ArgumentParser(description="Generates a slice profile, listing regions which are above a specified coverage threshold")
    parser.add_argument("-i", "--input", metavar="BAM", required=True, nargs="+", type=lambda x: isValidFile(x, parser), help="Input BAM/CRAM file(s)")
    parser.add_argument("-s", "--slice", metavar="INT", type=int, default=20, help="Slice depth [default: %(default)s]")
    parser.add_argument("-o", "--output", metavar="TSV", required=True, help="Output tab-deliniated slice profile")
    parser.add_argument("--minimum_slice_length", metavar="INT", type=int, default=180, help="Minimum genomic length of a slice [default: %(default)s]")
    parser.add_argument("--max_depth", metavar="INT", type=int, default=20000, help="Maximum coverage ceiling [default: %(default)s]")
    parser.add_argument("--flank", metavar="INT", default=60, help="Add this many bases onto each region")
    parser.add_argument("--verbose", action="store_true", help="Print status messages")

    args = parser.parse_args()

    return args


def main(args=None):
    if args is None:
        args = getArgs()

    # If there are more than one input file, we need to merge them and downsample them
    if len(args.input) > 1:
        if args.verbose:
            sys.stderr.write("Merging and downsampling BAM files" + os.linesep)
        tmpOutMerged = args.output + "tmpMergedBAM.bam"
        mergeArgs = [tmpOutMerged]
        mergeArgs.extend(args.input)
        pysam.merge(*mergeArgs)
        downsampleRate = 1.0 / float(len(args.input))
        inSeq = args.output + "tmpDownsampled.bam"
        pysam.view("-b", "-s", downsampleRate, "-o", inSeq, tmpOutMerged)
        os.remove(tmpOutMerged)
    else:
        inSeq = args.input[0]

    # Check to make sure the input BAM is indexed
    if not os.path.exists(inSeq + ".bai"):
        if args.verbose:
            sys.stderr.write("Indexing input BAM file" + os.linesep)
        pysam.index(inSeq)

    inFile = pysam.AlignmentFile(inSeq)

    openSlices = {}
    finishedSlices = {}
    previousSlice = 0
    lastPos = 0
    previousChrom = None
    posProcessed = 0

    for pileupCol in inFile.pileup(max_depth=args.max_depth):

        currentSlice = int(pileupCol.nsegments / args.slice)

        if previousChrom != pileupCol.reference_name:
            # We have changed chromosomes
            if previousChrom is not None:
                # Finish processing all open slices
                for i, finishedSlice in openSlices.items():
                    finishedSlice.end = lastPos + args.flank

                    addSlice(finishedSlices, finishedSlice, i)

                openSlices = {}
            previousChrom = pileupCol.reference_name
            previousSlice = 0

        # Create new slices
        for i in range(previousSlice, currentSlice):
            if i not in openSlices:
                openSlices[i] = segment(pileupCol.reference_name, pileupCol.reference_pos + 1 - args.flank)

        # We have reduced in slice number. Process the previous slices
        for i in range(currentSlice, previousSlice):
            finishedSlice = openSlices[i]
            finishedSlice.end = lastPos
            del openSlices[i]

            addSlice(finishedSlices, finishedSlice, i)

        lastPos = pileupCol.reference_pos + 1
        previousSlice = currentSlice

        # Status messages
        posProcessed += 1
        if posProcessed % 1000000 == 0:
            if args.verbose:
                sys.stderr.write("Positions Processed: %s. Current Position: %s." % (posProcessed,
                                pileupCol.reference_name + ":" + str(pileupCol.reference_pos + 1)) + os.linesep)

    # Now that we have finished processing the BAM file, process any open slices
    for i, finishedSlice in openSlices.items():
        finishedSlice.end = lastPos  + args.flank
        addSlice(finishedSlices, finishedSlice, i)

    # Write out slices
    if args.verbose:
        sys.stderr.write("Writing slices" + os.linesep)
    addLeft = True
    with open(args.output, "w") as o:
        for sliceLevel, regions in finishedSlices.items():
            sliceLevel += 1  # To offset 0-based indexing
            for region in regions:
                # Add bases to regions that are smaller than the specified length
                sliceDiff = region.end - region.start - args.minimum_slice_length
                if sliceDiff < 0:
                    if addLeft:
                        region.start += sliceDiff
                    else:
                        region.end -= sliceDiff
                    addLeft = not addLeft
                o.write("\t".join([region.chrom, str(region.start),
                                   str(region.end), str(sliceLevel * args.slice) + os.linesep]))

        # Finally, if a temp BAM file was created when we merged multiple input files for processing,
        # we need to delete it
        if len(args.input) > 1:
            os.remove(inSeq)

if __name__ == "__main__":
    main()