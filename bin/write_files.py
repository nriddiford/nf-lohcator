import fnmatch
import os
import sys
import re
from collections import defaultdict
from optparse import OptionParser
import pandas as pd
import glob
import json

# import fnmatch


def find_normal(options):
    config = pd.read_csv(options.config, delimiter="\t", index_col=False, header=0, na_filter=False, names=['tumour', 'normal'])

    d = dict(zip(config['tumour'].tolist(),config['normal'].tolist()))

    return d


def get_reads(options, id):
    file_path = os.path.join(options.dir)
    reads = glob.glob(file_path + id + '*')

    d = {}

    if reads:
        f_match = re.compile(".*(R1|forward).*" + options.ext)
        r_match = re.compile(".*(R2|reverse).*" + options.ext)
        # f_read = glob.glob(t_reads + '*R1*' + options.ext)
        r1 = list(filter(f_match.match, reads))[0] # Read Note
        r2 = list(filter(r_match.match, reads))[0] # Read Note

        d =  {'id': id, 'r1': r1, 'r2': r2}

    return d


def write_plan(options):
    print("Looking for '*%s' files in %s" % (options.ext, options.dir))
    print("Writing sample plan to '%s'" % (options.out_file))

    tumour_mapping = find_normal(options)

    d =  {}
    files = {}

    for tumour_id in tumour_mapping.keys():
        d['tumour'] = get_reads(options, tumour_id)
        if d['tumour']:
            d['normal'] = get_reads(options, tumour_mapping[tumour_id])
            files[tumour_id] = { 'tumour': d['tumour'], 'normal': d['normal'] }

    # print(json.dumps(files, indent=4, sort_keys=True))

    with open(options.out_file, 'w') as out:
        out.write('tumour_id,tr1,tr2,normal_id,nr1,nr2\n')
        for t_id in sorted(files):
            t_id = files[t_id]['tumour']['id']
            t_r1 = files[t_id]['tumour']['r1']
            t_r2 = files[t_id]['tumour']['r2']
            n_id = files[t_id]['normal']['id']
            n_r1 = files[t_id]['normal']['r1']
            n_r2 = files[t_id]['normal']['r2']
            print("-> %s is paired with %s" % (t_id, n_id))
            l = ','.join(map(str,[t_id, t_r1, t_r2, n_id, n_r1, n_r2]))
            out.write(l + '\n')

    return True


def get_args():
    parser = OptionParser()

    parser.add_option("--out_file", "-o", dest="out_file", action="store", help="File to write config to [Default 'my_sample_plan.csv']")
    parser.add_option("--config", "-c", dest="config", action="store", help="Mapping for tumour/normal samples [Default 'data/samples.tsv']")
    parser.add_option("--directory", "-d", dest="dir", action="store", help="Directory to look for files in [Default data/]")
    parser.add_option("--extension", "-e", dest="ext", action="store", help="Extension to search for [Default 'fastq.gz']")

    parser.set_defaults(dir='data/', ext='fastq.gz', config='data/samples2.txt', out_file='my_sample_plan.csv')
    options, args = parser.parse_args()

    return options, args


def main():
    options, args = get_args()

    try:
        write_plan(options)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return


if __name__ == "__main__":
    sys.exit(main())
