#!/usr/bin/env python3
#
# Either install dependencies manually (see list below) or run the script
# with 'uv run cellbarcodes.py'
#
# /// script
# dependencies = [
#   "tinyalign",
#   "pysam",
# ]
# ///
"""
Error-correct cell barcodes
"""

from itertools import product
import random
import sys
from collections import Counter
from typing import Dict, Tuple, Iterator
from tinyalign import hamming_distance
from argparse import ArgumentParser
from pysam import AlignmentFile, AlignedSegment

# TODO should become CR and CB
RAW_CELL_BARCODE_TAG = "CB"
CORRECTED_CELL_BARCODE_TAG = "cc"


def main():
    #    logging.basicConfig(format="%(message)s", level="INFO")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        "-o",
        help="Output BAM file. For demultiplexing, add {barcode} to the output file name to turn it into a filename template",
    )
    # parser.add_argument("--length", "-l", default=8, type=int, help="Barcode length")
    parser.add_argument(
        "--mismatches", "-k", default=1, type=int, help="No. of allowed mismatches per cell barcode"
    )
    parser.add_argument("bam")
    args = parser.parse_args()

    # Count how often cell barcodes occur
    counts1 = Counter()
    counts2 = Counter()
    with AlignmentFile(args.bam) as bam:
        for record in bam:
            # What is in the CB tag was concatenated from two index sequences
            cb1, cb2 = extract_cell_barcodes(record.get_tag(RAW_CELL_BARCODE_TAG))
            counts1[cb1] += 1
            counts2[cb2] += 1

    # TODO hardcoded numbers
    barcodes1 = [seq for (seq, count) in counts1.most_common()][:24]
    barcodes2 = [seq for (seq, count) in counts2.most_common()][:16]
    print("Detected cell barcodes:")
    print(*barcodes1)
    print(*barcodes2)

    corrector = CellbarcodeCorrector(barcodes1, barcodes2, args.mismatches)

    if not args.output:
        print("No output file given, exiting")
        return

    with AlignmentFile(args.bam) as inbam:
        with AlignmentFile(args.output, mode="wb", template=inbam) as outbam:
            for record in inbam:
                corrector.correct_one_read(record)
                outbam.write(record)

    print("Error statistics")
    print("Errors     count")
    error_counts = corrector.error_counts
    for errors, count in error_counts.most_common():
        print(f"{errors:6} {count:9}")


class CellbarcodeCorrector:
    def __init__(self, barcodes1, barcodes2, mismatches: int):
        self.barcodes1 = {f"{i+1}": seq for i, seq in enumerate(sorted(barcodes1))}
        self.barcodes2 = {f"{i+1}": seq for i, seq in enumerate(sorted(barcodes2))}

        self.index1 = make_index(self.barcodes1, mismatches)
        self.index2 = make_index(self.barcodes2, mismatches)

        self.error_counts = Counter()

    def correct_one_read(self, record: AlignedSegment):
        """
        Write an error-corrected cell barcode tag to the record

        The record is modified in place.
        """
        cb1, cb2 = extract_cell_barcodes(record.get_tag(RAW_CELL_BARCODE_TAG))
        if cb1 in self.index1 and cb2 in self.index2:
            name1, errors1 = self.index1[cb1]
            name2, errors2 = self.index2[cb2]
            corrected = self.barcodes1[name1] + self.barcodes2[name2]
            self.error_counts[errors1 + errors2] += 1
            record.set_tag(CORRECTED_CELL_BARCODE_TAG, corrected)


def extract_cell_barcodes(cb: str) -> tuple[str, str]:
    """Extract the two cell barcodes from the CB tag"""
    # TODO currently hardcoded
    return cb[:10], cb[10:]


def hamming_sphere(s: str, k: int) -> Iterator[str]:
    """
    Yield all strings t for which the hamming distance between s and t is exactly k,
    assuming the alphabet is A, C, G, T.
    """
    assert k >= 0
    if k == 0:
        yield s
        return
    n = len(s)

    # i is the first position that is varied
    for i in range(n - k + 1):
        prefix = s[:i]
        c = s[i]
        suffix = s[i + 1 :]
        for ch in "ACGT":
            if ch == c:
                continue
            for t in hamming_sphere(suffix, k - 1):
                y = prefix + ch + t
                assert len(y) == n
                yield y


def hamming_environment(s: str, k: int):
    """
    Find all strings t for which the hamming distance between s and t is at most k,
    assuming the alphabet is A, C, G, T.

    Yield tuples (t, e), where e is the hamming distance between s and t.
    """
    n = len(s)
    for e in range(k + 1):
        for t in hamming_sphere(s, e):
            yield t, e


def make_index(
    sequences: Dict[str, str], max_errors: int
) -> Dict[str, Tuple[str, int]]:
    """
    Pre-compute all possible strings that are at most max_errors (Hamming
    distance) away from the input sequences.
    """
    index: Dict[str, Tuple[str, int]] = dict()
    has_warned = False
    for name, sequence in sequences.items():
        for t, errors in hamming_environment(sequence, max_errors):
            if t in index:
                other_name, other_errors = index[t]
                if errors > other_errors:
                    continue
                if other_errors == errors and not has_warned:
                    print(
                        f"Sequences '{name}' ({sequence}) and '{other_name}' "
                        f"({sequences[other_name]}) are very similar. At {max_errors} allowed "
                        f"error(s), it is not possible to assign the sequence '{t}' uniquely to one "
                        f"of them because the number of errors is {errors} compared to both.",
                        file=sys.stderr,
                    )
                    # has_warned = True
            else:
                index[t] = (name, errors)
    print(f"Built an index containing {len(index)} strings.", file=sys.stderr)
    return index


if __name__ == "__main__":
    main()
