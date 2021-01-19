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
    config = pd.read_csv(options.config, delimiter="\t", index_col=False, na_filter=False, names=['sample', 'assay'])

    samples = config['sample'].tolist()
    it = iter(samples)

    s = {}
    for x in it:
        s[x] = next(it)

    return s

def get_reads(options, id, assay, d):
    file_path = os.path.join(options.dir)
    reads = glob.glob(file_path + id + '*')

    if reads:
        f_match = re.compile(".*(R1|forward).*" + options.ext)
        r_match = re.compile(".*(R2|reverse).*" + options.ext)
        # f_read = glob.glob(t_reads + '*R1*' + options.ext)
        r1 = list(filter(f_match.match, reads))[0] # Read Note
        r2 = list(filter(r_match.match, reads))[0] # Read Note

        print("Assay: %s, R1: %s, R2: %s" % (assay, r1, r2))

        # d.push({"tumour_id": id})

        return d

def write_plan(options):
    print("Looking for '*%s' files in %s" % (options.ext, options.dir))
    print("Writing sample plan to '%s'" % (options.out_file))

    tumour_mapping = find_normal(options)

    file_mapping = [
        {
            'tumour_id': [],
            'normal_id': []
        }
    ]
            # 'tumour_id': [],
            # 'normal_id': []

    for tumour_id in tumour_mapping.keys():

        # for file in glob.glob('data/HUM-1*'):
        # file_mapping[glob.glob('data/HUM-1*')] = glob.glob('data/HUM-3*')
        # file_mapping[tumour_id] = glob.glob('data/HUM-1*')
        # t_reads = glob.glob('data/' + tumour_id + '*')

        file_mapping = get_reads(options, tumour_id, 'tumour', file_mapping)

        # if t_reads:
        #     r1 = re.compile(".*R1.*" + options.ext)
        #     r2 = re.compile(".*R2.*" + options.ext)
        #     # f_read = glob.glob(t_reads + '*R1*' + options.ext)
        #     tr1 = list(filter(r1.match, t_reads))[0] # Read Note
        #     tr2 = list(filter(r2.match, t_reads))[0] # Read Note
        #
        #     file_mapping[tumour_id] = [tr1, tr2]

    # with open(options.out_file, 'w') as out:
    #     out.write('tumour_id,tr1,tr2,normal_id,nr1,nr2\n')
    #     r1 = ''
    #     for file in sorted(os.listdir(options.dir)):
    #         if file.endswith(options.ext):
    #             sample = str(file.split('_')[0])
    #             if len(sample.split('.')) > 1:
    #                 sample = str(file.split('.')[0])
    #             condition = 'normal'
    #             if sample in tumour_mapping.keys():
    #                 condition = 'tumour'
    #             if fnmatch.fnmatch(file, '*.R1.*'):
    #                 r1 = os.path.join(options.dir, file)
    #                 continue
    #             r2 = os.path.join(options.dir, file)
    #             print("Sample %s is %s: [%s, %s]" % (sample, condition, r1, r2))
    #             l = ','.join(map(str,[sample, r1, r2]))
    #             if condition == 'tumour':
    #                 out.write(l + ',')
    #             else:
    #                 out.write(l + '\n')

    print(json.dumps(file_mapping, indent=4, sort_keys=True))

    return True


def get_args():
    parser = OptionParser()

    parser.add_option("--out_file", "-o", dest="out_file", action="store", help="File to write config to [Default 'my_sample_plan.csv']")
    parser.add_option("--config", "-c", dest="config", action="store", help="Mapping for tumour/normal samples [Default 'data/samples.tsv']")
    parser.add_option("--directory", "-d", dest="dir", action="store", help="Directory to look for files in [Default data/]")
    parser.add_option("--extension", "-e", dest="ext", action="store", help="Extension to search for [Default 'fastq.gz']")

    parser.set_defaults(dir='data/', ext='fastq.gz', config='data/samples.tsv', out_file='my_sample_plan.csv')
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
