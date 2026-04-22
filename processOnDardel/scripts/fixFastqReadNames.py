#
# Python script which takes in args to pass on to function to fix read names in fastq files
#
# Outline of the fix: 
#   expected read name structure:
#       - if the read contains a UMI: readID_UMI someString:cellBarcode
#         e.g.: 'A00621:811:HTTHFDRX2:1:2101:1253:1063_CGGGGGTT 1:N:0:CCATCACTAT+TCCAATGCAA'
#       - if the read does not contain a UMI: readID someString:cellBarcode
#         e.g.: 'A00621:811:HTTHFDRX2:1:2101:1253:1063 1:N:0:CCATCACTAT+TCCAATGCAA'
#
#   see https://samtools.github.io/hts-specs/SAMtags.pdf for the choice
#   to use OX and CR to referer to the uncorrected UMI and cell barcode respectively
#
#   how the new read name should look: 
#       - if the read contains a UMI: readID_OX:UMI_CR:index1-index2
#         e.g.: 'A00621:811:HTTHFDRX2:1:2101:1253:1063_OX:CGGGGGTT_CR:CCATCACTAT-TCCAATGCAA'
#       - if the read does not contain a UMI: readID_CR:index1-index2
#         e.g.: 'A00621:811:HTTHFDRX2:1:2101:1253:1063_CR:CCATCACTAT-TCCAATGCAA'
#
#   if the read ID contains "_" at length(ID)-8 then it has a UMI (UMI is always 8 bases long)
#       this way of checking is to circumvent the unlikely case where "_" is present in the rest of the ID


## packages
import dnaio
import argparse

## function to fix read names

def fixReadNames(inR1, inR2, outR1, outR2, nThreads=1):

    with dnaio.open(file1=inR1, file2=inR2, mode="r", open_threads=nThreads) as inFastq, \
        dnaio.open(file1=outR1, file2=outR2, mode="w", open_threads=nThreads) as outFastq:

        ## read counter
        readCounter = 0

        for r1, r2 in inFastq:
            
            countUnderscore = r1.id.count("_")

            ## case 1: read name contains UMI
            if countUnderscore == 1:

                # sanity check: "_" should be at position len(id)-8
                if r1.id[-9] != "_" or r2.id[-9] != "_":
                    raise ValueError(f"The '_' was found at an unexpected position for read {r1.id}")

                # sanity check: r1 and r2 should have the same UMI sequence
                if r1.id[-8:] != r2.id[-8:]: 
                    raise ValueError(f"R1 and R2 have different UMIs for read {r1.id}")
            
                # fix read names
                r1Mod = dnaio.SequenceRecord(
                    name = r1.id[:-9] + "_OX:" + r1.id[-8:] + "_CR:" + r1.comment.replace("+", "-").split(":")[-1], 
                    sequence = r1.sequence, 
                    qualities = r1.qualities
                )
                r2Mod = dnaio.SequenceRecord(
                    name = r2.id[:-9] + "_OX:" + r2.id[-8:] + "_CR:" + r2.comment.replace("+", "-").split(":")[-1], 
                    sequence = r2.sequence, 
                    qualities = r2.qualities
                )
            
            ## case 2: read name contains no UMI
            elif countUnderscore == 0: 

                # fix read names
                r1Mod = dnaio.SequenceRecord(
                    name = r1.id + "_CR:" + r1.comment.replace("+", "-").split(":")[-1], 
                    sequence = r1.sequence, 
                    qualities = r1.qualities
                )
                r2Mod = dnaio.SequenceRecord(
                    name = r2.id + "_CR:" + r2.comment.replace("+", "-").split(":")[-1], 
                    sequence = r2.sequence, 
                    qualities = r2.qualities
                )

            ## case 3: sanity check for cases where the "_" count is greater than 1
            else:
                raise ValueError(f"More than one underscore found. Stopping at read {r1.id}.")
            
            ## update read counter
            readCounter += 1

            ## write to output fastq files
            outFastq.write(r1Mod, r2Mod)

        ## print out the number of processed reads
        print(f"Processed {readCounter} reads.")


## read in passed args and fix read names

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process given R1 and R2 fastq files")
    parser.add_argument("--inFastq1", required=True, help="Input R1 FASTQ file")
    parser.add_argument("--inFastq2", required=True, help="Input R2 FASTQ file")
    parser.add_argument("--outFastq1", required=True, help="Output R1 FASTQ file")
    parser.add_argument("--outFastq2", required=True, help="Output R2 FASTQ file")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")

    args = parser.parse_args()

    fixReadNames(args.inFastq1, args.inFastq2, args.outFastq1, args.outFastq2, nThreads=args.threads)

