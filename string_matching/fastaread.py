import sys
import time


def read(path: str):
    sequences = []
    with open(path) as f:
        data = f.read()
        raw_sequences = data.split(">")
        raw_sequences = raw_sequences[
            1:]  # remove (typically empty) string before first ">"

    for seq in raw_sequences:
        start_idx = seq.find("\n")
        seq = seq[start_idx + 1:]  # skip description line
        seq = seq.replace("\n", "")
        sequences.append(seq)

    return sequences


def main():
    file_paths = sys.argv[1:]
    sequences = []

    for path in file_paths:
        sequences.extend(read(path))

    for seq in sequences:
        print(len(seq))
        # print(seq)


if __name__ == "__main__":
    main()
