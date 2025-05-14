import sys
from fastaread import read


def find_all(text: str, pattern: str, limit: int = -1):
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
            gsr_step = 1  # TODO:
            step = max(bcr_step, gsr_step)
        outer_text_pos += step

        limit_reached = limit != -1 and len(matches) >= limit
    return matches


def bcr(char: str, idx: int, bcr_table) -> int:
    return bcr_table[idx, char]


def gcr():
    pass


def compute_bcr_table(pattern: str, alphabet=None):
    """
    The computed table stores offsets for each index *i* in the pattern, addressing next occurrance of a character *c* within the pattern when looking to the left of *i*.

    It can be accessed like this: table[i, c]
    """
    if not alphabet:  # assume all of the alphabets characters occur in the pattern
        alphabet = list(set(pattern))

    bcr_table = {}
    # bcr_table = [[0 for j in len(alphabet)] for i in len(pattern)]

    for idx in range(len(pattern)):
        for char in alphabet:
            left_occurrence_idx = pattern.rfind(char, 0, idx)
            if left_occurrence_idx == -1:  # not found
                bcr_table[idx, char] = 0
            else:
                bcr_table[idx, char] = idx - left_occurrence_idx

    return bcr_table


def compute_gsr_table():

    pass


def main():
    text_file = sys.argv[1]
    pattern_file = sys.argv[2]

    text = read(text_file)[0]
    patterns = read(pattern_file)

    for pattern in patterns:
        occurrences = find_all(text, pattern, limit=10)
        print(*occurrences, sep=' / ')


if __name__ == "__main__":
    main()
