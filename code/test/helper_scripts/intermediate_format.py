import sys


def main():
    precision = 10
    with sys.stdin as reader:
        for line in reader:
            line = line.strip()
            tokens = line.split('\t')
            if '=' in tokens[-1]:
                eq_tokens = tokens[-1].split(' ')
                for i in range(len(eq_tokens)):
                    try:
                        float(eq_tokens[i])
                    except ValueError:
                        continue
                    if len(eq_tokens[i]) > precision:
                        eq_tokens[i] = eq_tokens[i][:precision]
                tokens[-1] = ' '.join(eq_tokens)
                line = '\t'.join(tokens)
            print(line)


if __name__ == "__main__":
    main()
