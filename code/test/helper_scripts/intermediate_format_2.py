import sys


def main():
    precision = 10
    with sys.stdin as reader:
        for line in reader:
            line = line.strip()
            tokens = line.split('\t')
            # limit float sizes to 10 characters
            for j in range(len(tokens)-2, len(tokens)):
                float_tokens = tokens[j].split(',')
                for i in range(len(float_tokens)):
                    try:
                        float(float_tokens[i])
                    except ValueError:
                        continue
                    if len(float_tokens[i]) > precision:
                        float_tokens[i] = float_tokens[i][:precision]
                tokens[j] = ','.join(float_tokens)
            # check if alt ids are equal, sorting is messed up from py2 to 3
            id_toks = tokens[-2].split(',')
            if len(id_toks) > 1 and id_toks[0] == id_toks[1]:
                tokens[-3] = ','.join(sorted(tokens[-3].split(',')))
            line = '\t'.join(tokens)
            print(line)


if __name__ == "__main__":
    main()
