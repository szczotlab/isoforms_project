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
from pysam import AlignmentFile


def main():
    #    logging.basicConfig(format="%(message)s", level="INFO")
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output", "-o", help="Output BAM file")
    # parser.add_argument("--length", "-l", default=8, type=int, help="Barcode length")
    # parser.add_argument("--mismatches", "-k", default=2, type=int, help="No. of mismatches")
    parser.add_argument("bam")
    args = parser.parse_args()

    # Count how often cell barcodes occur
    counts1 = Counter()
    counts2 = Counter()
    with AlignmentFile(args.bam) as bam:
        for record in bam:
            cell_barcode = record.get_tag("CB")
            # What is in the CB tag was concatenated from two index sequences
            cb1, cb2 = cell_barcode[:10], cell_barcode[10:]
            counts1[cb1] += 1
            counts2[cb2] += 1

    barcodes1 = [seq for (seq, count) in counts1.most_common()][:24]
    barcodes2 = [seq for (seq, count) in counts2.most_common()][:16]
    print("Detected cell barcodes:")
    print(*barcodes1)
    print(*barcodes2)

    barcodes1 = {f"{i+1}": seq for i, seq in enumerate(sorted(barcodes1))}
    barcodes2 = {f"{i+1}": seq for i, seq in enumerate(sorted(barcodes2))}

    index1 = make_index(barcodes1, 2)
    index2 = make_index(barcodes2, 2)

    if not args.output:
        print("No output file given, exiting")
        return

    error_counts = Counter()
    with AlignmentFile(args.bam) as inbam:
        with AlignmentFile(args.output, mode="wb", template=inbam) as outbam:
            for record in inbam:

                cell_barcode = record.get_tag("CB")
                cb1, cb2 = cell_barcode[:10], cell_barcode[10:]

                if cb1 in index1 and cb2 in index2:
                    name1, errors1 = index1[cb1]
                    name2, errors2 = index2[cb2]
                    corrected = barcodes1[name1] + barcodes2[name2]
                    error_counts[errors1 + errors2] += 1
                    record.set_tag("cc", corrected)

                outbam.write(record)

    print("Error statistics")
    print("Errors     count")
    for (errors, count) in error_counts.most_common():
        print(f"{errors:6} {count:9}")


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
