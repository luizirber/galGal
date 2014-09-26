#! /usr/bin/env python


def append_barcode(inputfile, outputfile):
    with open(outputfile, 'w') as output:
        with open(inputfile, 'r') as fd:
            for line in fd:
                line = line.strip()
                if 'Barcode' in line:
                    current_barcode = line[1:].replace(' ', '_').replace(':', '')
                elif line.startswith('>'):
                    output.write("".join([line, "_", current_barcode, "\n"]))
                else:
                    output.write(line + "\n")


if __name__ == "__main__":
    import sys
    append_barcode(sys.argv[1], sys.argv[2])
