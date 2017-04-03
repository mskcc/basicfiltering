#!/usr/bin/python
"""postprocess"""

import argparse
import ruamel.yaml
import os

def read(filename):
    """return file contents"""

    with open(filename, 'r') as file_in:
        return file_in.read()


def write(filename, cwl):
    """write to file"""

    with open(filename, 'w') as file_out:
        file_out.write(cwl)


def main():
    """main function"""

    parser = argparse.ArgumentParser(description='postprocess')

    parser.add_argument(
        '-f',
        action="store",
        dest="filename_cwl",
        help='Name of the cwl file',
        required=True
    )

    params = parser.parse_args()
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cwl = ruamel.yaml.load(read(params.filename_cwl),
                           ruamel.yaml.RoundTripLoader)

    script_path = os.path.join(dir_path,'filter_vardict.py')
    cwl['baseCommand'] = ['python',script_path]
    cwl['inputs']['inputVcf']['type'] = ['string', 'File']
    cwl['inputs']['hotspotVcf']['type'] = ['null', 'string', 'File']

    write(params.filename_cwl, ruamel.yaml.dump(
        cwl, Dumper=ruamel.yaml.RoundTripDumper))


if __name__ == "__main__":

    main()