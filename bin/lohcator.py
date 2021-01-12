import sys, os
import fnmatch
from optparse import OptionParser
import pandas as pd
import ntpath
from collections import defaultdict
import json
import vcf

pd.set_option('display.float_format', lambda x: '%.3f' % x)

def find_normal(options):
    sample = ntpath.basename(options.freebayes_file).split("_")[0]
    config = pd.read_csv(options.config, delimiter="\t", index_col=False, na_filter=False, names=['sample', 'assay'])

    samples = config['sample'].tolist()
    it = iter(samples)

    s = {}
    for x in it:
        s[x] = next(it)

    if(s[sample]):
        print("Tumour: %s" % sample)
        print("Normal: %s" % s[sample])
        return sample, s[sample]
    else:
        print("Cant find corresponding normal sample for %s" % sample)


def parse_freebayes(options):

    tumour, normal = find_normal(options)

    mode = 'r'
    if options.freebayes_file.endswith('.gz'):
        mode = 'rb'

    vcf_reader = vcf.Reader(open(options.freebayes_file, mode))

    for record in vcf_reader.fetch(options.chromosome):
        if(record.genotype(tumour)['DP'] and record.genotype(tumour)['DP'] > 20 and record.genotype(normal)['DP'] and record.genotype(normal)['DP'] > 20 and record.genotype(tumour)['GQ'] > 1 ):
            if (record.genotype(tumour)['GT'] == '0/0' and record.genotype(normal)['GT'] == '0/0') or (record.genotype(tumour)['GT'] == '1/1' and record.genotype(normal)['GT'] == '1/1'):
                continue
            if 'snp' not in record.INFO['TYPE'] or len(record.INFO['TYPE']) > 1:
                continue

            try:
                status = record.INFO['VT']
            except KeyError:
                status = 'germline'
                pass

            # Need to skip complex calls
            # Should also skip non Germline...
            accepted_genotypes = ['0/0', '0/1', '1/0', '1/1']

            if record.genotype(tumour)['GT'] not in accepted_genotypes or record.genotype(normal)['GT'] not in accepted_genotypes:
                continue
            taf = round((record.genotype(tumour)['AO'] / (record.genotype(tumour)['AO'] + record.genotype(tumour)['RO'])),2)
            naf = round((record.genotype(normal)['AO'] / (record.genotype(normal)['AO'] + record.genotype(normal)['RO'])),2)

            af_diff = abs(naf - taf)

            t_a_count = record.genotype(tumour)['AO']
            n_a_count = record.genotype(normal)['AO']

            if(isinstance(t_a_count, list)):
                t_a_count = t_a_count[0] # Some records contain two values ... (?)
            if(isinstance(n_a_count, list)):
                n_a_count = n_a_count[0] # Some records contain two values ... (?)

            t_freq = round((t_a_count/record.genotype(tumour)['DP'])*100, 2)
            n_freq = round((n_a_count/record.genotype(normal)['DP'])*100, 2)

            difference = abs(t_freq - n_freq)

            # if difference > 20:
                # print("LOH", record.INFO, record.genotype(tumour))


def parse_varscan(options):

    sample = ntpath.basename(options.varscan_file).split(".")[0]

    bed_file = sample + '_LOH_regions.bed'

    df = pd.read_csv(options.varscan_file, delimiter="\t", dtype={"chrom": str, "position": int, "normal_var_freq": str, "tumor_var_freq": str, "somatic_status": str, "somatic_p_value": str})
    df = df.sort_values(['chrom', 'position'])

    chroms = ['2L', '2R', '3L', '3R', 'X', 'Y', '4']
    loh = defaultdict(lambda: defaultdict(dict))
    start = False
    start_chain = False
    loh_count = defaultdict(int)
    last_informative_snp = defaultdict(lambda: defaultdict(int))

    min_count = 10
    min_size = 30000
    upper_threshold = 100 - options.loh_threshold
    lower_threshold = options.loh_threshold


    if(options.lenient):
        print("--lenient set to False. Will look for smaller LOH regions")
        bed_file = sample + '_LOH_regions_lenient.bed'
        min_count = 4
        min_size = 10000

    for row in df.itertuples(index=True, name='Pandas'):
        chrom, pos, n_freq, t_freq, snv_type, p_val = getattr(row, "chrom"), getattr(row, "position"), getattr(row, "normal_var_freq"), getattr(row, "tumor_var_freq"), getattr(row, "somatic_status"), float(getattr(row, "somatic_p_value"))

        if chrom not in chroms: continue

        t_freq = float(t_freq.rstrip("%"))
        n_freq = float(n_freq.rstrip("%"))
        t_freq += 0.001
        n_freq += 0.001

        if snv_type == 'LOH':
            start_chain = True
            chain = True
        elif snv_type == 'Germline' and t_freq <= lower_threshold or t_freq >= upper_threshold:
            chain = True
            start_chain = False
        elif snv_type == 'Germline' and p_val > 0.5:
            chain = True
            start_chain = False
        elif snv_type == 'Germline' and abs(n_freq-t_freq) >= 10:
            chain = True
            start_chain = False
        elif snv_type == 'Somatic':
            chain = True
            start_chain = False
        else:
            chain = False
            start = False
            # This ^ is wrong

        if chrom not in loh:
            start = pos

        if chain and snv_type == 'LOH' and loh_count[start] > 5:
            last_informative_snp[chrom] = pos

        if chain and start:
            loh[chrom].setdefault(start, []).append(pos)
            start_chain = False
            if snv_type == 'LOH':
                loh_count[start] += 1
        elif start_chain:
            start = pos

    loh_run = defaultdict(lambda: defaultdict(dict))

    with open(bed_file, 'w') as bed_out:
        for c in sorted(loh):
            for p in sorted(loh[c]):
                start = int(p)
                end = int(max(loh[c][p]))
                if p in loh_count:
                    if loh_count[p] >= min_count or loh_count[p] >= (min_count/2) and (end - start > min_size):
                    # if loh_count[p] >= 10 or loh_count[p] >= 5 and (end - start > 30000):

                        loh_run[c][p] = loh[c][p]
                        bed_out.write('%s\t%s\t%s\n' %(c, p, end))

    # print(json.dumps(last_informative_snp[options.chromosome], indent=4, sort_keys=True))

    if(loh_run[options.chromosome]):
        max_start = max(loh_run[options.chromosome])
        max_end   = max(loh_run[options.chromosome][max_start])

        print("Last LOH snp on %s: %s" % (options.chromosome, last_informative_snp[options.chromosome]))

        breakpoint_window_start = last_informative_snp[options.chromosome]
        breakpoint_window_end   = max_end + options.window

        print("Breakpoint window on %s (+/- %s): %s:%s-%s" % (options.chromosome, options.window, options.chromosome, breakpoint_window_start, breakpoint_window_end))

        if options.write_breakpoint:
            breakpoint = sample + '_breakpoint_region.bed'
            print("Writing breakpoint region to %s" % breakpoint)
            with open(breakpoint, 'w') as breakpoint_out:
                breakpoint_out.write('%s\t%s\t%s\n' %(options.chromosome, breakpoint_window_start, breakpoint_window_end))

    # print(json.dumps(loh_run, indent=4, sort_keys=True))
    return True


def main():
    parser = OptionParser()

    parser.add_option("-v", "--varscan_file", dest="varscan_file", help="Varscan native file")
    parser.add_option("-f", "--freebayes_file", dest="freebayes_file", help="Freebayes VCF file")
    parser.add_option("-w", "--window", dest="window", action="store", help="Print window at breakpoint on 2L")
    parser.add_option("-c", "--chromosome", dest="chromosome", action="store", help="Chromosome to look for LOH on [Default = 2L]")
    parser.add_option("--lenient", dest="lenient", action="store_true", help="Make more, smaller, lower confidence LOH region calls")
    parser.add_option("--loh_threshold", dest="loh_threshold", action="store", type="int", help="Set the threshold at which to call LOH [Default = 25%]")

    parser.add_option("--write_breakpoint", dest="write_breakpoint", action="store_true", help="Write window at breakpoint on 2L as bed file")

    parser.add_option("--config", dest="config", action="store", help="mapping for tumour/normal samples")


    parser.set_defaults(config='/Users/Nick_curie/Desktop/script_test/alleleFreqs/data/samples.tsv',
                        chromosome='2L',
                        loh_threshold=25,
                        window=5000)

    options, args = parser.parse_args()

    if options.varscan_file is None and options.freebayes_file is None:
        parser.print_help()
        print
    else:
        try:
            if options.freebayes_file:
                parse_freebayes(options)
            else:
                parse_varscan(options)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())
