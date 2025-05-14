import sys
from fastaread import read
from typing import List, Dict
import time


def find_all_naive(text: str, pattern: str, limit: int = -1) -> List[int]:
    matches = []

    pos = 0
    found = True
    limit_reached = False
    while found and not limit_reached:
        match = text.find(pattern, pos)
        found = match != -1
        if found:
            matches.append(match)
            pos = match + 1
        limit_reached = limit != -1 and len(matches) >= limit

    return matches


def find_all_bm(text: str, pattern: str, limit: int = -1) -> List[int]:
    """
    Finds all occurrences of *pattern* in *text* using the Boyer-Moore algorithm.

    Returns a list of indices to those occurrences.
    """
    alphabet = list(set(text))
    bcr_table = compute_bcr_table(pattern, alphabet=alphabet)
    gsr_table = compute_gsr_table(pattern)
    pattern_len = len(pattern)
    text_len = len(text)

    matches = []

    outer_text_pos = pattern_len - 1
    limit_reached = False
    while outer_text_pos < text_len and not limit_reached:
        match = True
        pattern_pos = pattern_len - 1
        inner_text_pos = outer_text_pos
        while match and pattern_pos >= 0:
            match = text[inner_text_pos] == pattern[pattern_pos]
            pattern_pos -= 1
            inner_text_pos -= 1
        pattern_pos += 1
        inner_text_pos += 1

        step = 1
        if match:
            match_pos = outer_text_pos - pattern_len + 1
            matches.append(match_pos)
        else:
            mismatched_character = text[inner_text_pos]
            bcr_step = bcr(mismatched_character, pattern_pos, bcr_table)
            gsr_step = gsr(pattern_pos, gsr_table)
            step = max(bcr_step, gsr_step)
        outer_text_pos += step

        limit_reached = limit != -1 and len(matches) >= limit
    return matches


def bcr(char: str, idx: int, bcr_table: Dict) -> int:
    return bcr_table[idx, char]


def gsr(pattern_pos: int, gsr_table: List[int]) -> int:
    return gsr_table[pattern_pos + 1]


def compute_bcr_table(pattern: str, alphabet=None) -> Dict:
    """
    The computed table stores offsets for each index *i* in the pattern, addressing next occurrence of a character *c* within the pattern when looking to the left of *i*.

    It can be accessed like this: table[i, c]
    """
    if not alphabet:  # assume all of the alphabets characters occur in the pattern
        alphabet = list(set(pattern))

    bcr_table = {}

    for pattern_idx in range(len(pattern)):
        for char in alphabet:
            left_occurrence_idx = pattern.rfind(char, 0, pattern_idx)
            if left_occurrence_idx == -1:  # not found -> we can skip the entire part left of the current position in the pattern
                bcr_table[pattern_idx, char] = pattern_idx + 1
            else:
                bcr_table[pattern_idx,
                          char] = pattern_idx - left_occurrence_idx

    return bcr_table


def compute_gsr_table(pattern: str) -> List[int]:
    gsr_table = [1 for i in range(len(pattern) + 1)]

    half = len(pattern) // 2  # we won't find pattern P in T if |T| < |P|
    for suffix_idx in range(half, len(pattern)):
        suffix = pattern[suffix_idx:]
        match_idx = pattern.rfind(suffix, 0, suffix_idx)
        found = match_idx == -1
        if found:
            if not pattern[match_idx - 1] == pattern[suffix_idx - 1]:
                offset = suffix_idx - match_idx
                gsr_table[suffix_idx] = offset
        else:
            gsr_table[suffix_idx] = 1

    return gsr_table


def main():
    text_file = sys.argv[1]
    pattern_file = sys.argv[2]

    text = read(text_file)[0]
    patterns = read(pattern_file)

    for pattern in patterns:
        # naive_start = time.time()
        # naive_occurrences = find_all_naive(text, pattern, limit=10)
        # naive_end = time.time()

        # print("Time (Naive): " + str(naive_end - naive_start))

        # bm_start = time.time()
        bm_occurrences = find_all_bm(text, pattern, limit=10)
        # bm_end = time.time()
        # print("Time (BM): " + str(bm_end - bm_start))

        # print("Matching Results? (Naive and BM): " +
        #       str(bm_occurrences == naive_occurrences))

        # print(*naive_occurrences, sep=' / ')
        print(*bm_occurrences, sep=' / ')


if __name__ == "__main__":
    main()
