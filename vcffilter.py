#!/usr/bin/env python
from __future__ import division

# Version Check
###
import sys
pyversion = sys.version_info[0:2]
if pyversion < (2, 7):
    print '%s requires Python 2.7 (you have %s)' % (sys.argv[0],
                                                    '.'.join(pyversion))
    exit(1)

# Imports
###
import itertools
import datetime
import os
import argparse
import gzip
import bz2

from functools import wraps


# Arguement parsing
###
parser = argparse.ArgumentParser(description='Filter variants in a VCF file')
parser.add_argument('-f', '--file', required=True, metavar='vcffile',
                    help='VCF file for processing')
parser.add_argument('-o', '--out', required=False, dest='outfile', metavar='outfile',
                    help='File for output', default=os.devnull)
parser.add_argument('-r', '--region', required=False, nargs=3,
                    help='Constrain to region')
parser.add_argument('--rs', action='store_true', dest='rsonly',
                    help='Only return variants with rs numbers')
parser.add_argument('--snv', action='store_true', dest='snvonly',
                    help='Only return SNVs')
parser.add_argument('-q', '--qual', required=False, type=float, dest='minqual',
                    help='Minimum quality score to include')
parser.add_argument('-g', '--geno', required=False, default=None, dest='min_call_rate',
                    help="Minimum genotype call rate", type=float)
parser.add_argument('--no-qc', required=False, action='store_false', dest='qcfilter',
                    help="Do not filter on FILTER column == PASS")
parser.add_argument('--info_filter', dest='ifilters', nargs=3, action='append',
                    help='Filter on info string')
parser.add_argument('--model', dest='model', nargs=1, action='store', required=False,
                    help='Filter variants consistent with mendelian inheritance',
                    choices=['dom', 'rec'], default=None)
parser.add_argument('--quiet', action='store_true',
                    help='Suppress some output')
args = parser.parse_args()

# Functions
###


def make_numeric(value):
    ''' Makes a value numeric if possible '''
    try:
        return float(value)
    except ValueError:
        return value


def make_info_condition_function(condition):
    # on is the field we're filtering on
    # op is the operation we're performing on the field
    # value is the threshold/condition
    on, op, value = condition
    value = make_numeric(value)

    def err_on_nonnumeric(v):
        try:
            float(v)
        except ValueError:
            raise ValueError(
                'Nonnumeric value (%s) for operator %s' % (value, op))

    def false_on_keyerror(f):
        @wraps(f)
        def wrapper(*args, **kwds):
            try:
                return f(*args, **kwds)
            except KeyError:
                return False
        return wrapper
    if op == 'is' and value == 'set':
        return lambda x: on in x['INFO']
    elif op == 'not' and value == 'set':
        return lambda x: on not in x['INFO']
    elif op == 'gt':
        err_on_nonnumeric(value)

        @false_on_keyerror
        def f(x):
            return x['INFO'][on] > value
        return f
    elif op == 'gte':
        err_on_nonnumeric(value)

        @false_on_keyerror
        def f(x):
            return x['INFO'][on] >= value
        return f
    elif op == 'lt':
        err_on_nonnumeric(value)

        @false_on_keyerror
        def f(x):
            return x['INFO'][on] < value
        return f
    elif op == 'lte':
        err_on_nonnumeric(value)

        @false_on_keyerror
        def f(x):
            return x['INFO'][on] <= value
        return f
    elif op == 'eq':
        @false_on_keyerror
        def f(x):
            return x['INFO'][on] == value
        return f
    elif op == 'neq':
        @false_on_keyerror
        def f(x):
            return x['INFO'][on] != value
        return f
    elif op == 'contains':
        @false_on_keyerror
        def f(x):
            return value in x['INFO'][on]
        return f
    elif op == 'ncontains':
        @false_on_keyerror
        def f(x):
            return value not in x['INFO'][on]
        return f
    else:
        raise ValueError('Unknown operation for filter')


def parse_info_conditions(conditionlist):
    ''' Makes closures for the conditions '''
    return [make_info_condition_function(condition)
            for condition in conditionlist]


def parse_info(inf):
    ''' Parses a VCF info string into a key-value dictionary '''
    inf = inf.split(';')
    info_dict = {}
    for s in inf:
        s = s.split('=')
        if len(s) == 1:
            s = s[0]
            info_dict[s] = s
        else:
            key, value = s
            value = make_numeric(value)
            info_dict[key] = value
    return info_dict


def parse_genotype(g):
    if '/' in g:
        g = g.split('/')
    elif '|' in g:
        g = g.split('|')
    else:
        raise ValueError('Bad genotype: %s' % g)
    if g == ['.', '.']:
        return None
    else:
        return int(g[0]), int(g[1])


def get_genotypes_from_record(record):
    vcfcols = {'CHROM', 'POS', 'ID', 'REF', 'ALT',
               'QUAL', 'FILTER', 'INFO', 'FORMAT'}

    genotype_cols = {k: record[k].split(':') for k in record if k not in vcfcols}
    format = record['FORMAT'].split(':')

    genotypes = [dict(zip(format, genotype_cols[k]))['GT']
                 for k in genotype_cols]
    genotypes = [parse_genotype(g) for g in genotypes]
    return genotypes


def meets_conditions(inf, conditions):
    ''' Returns a list of whether variants meet the conditions for filtering '''
    return [condition(inf) for condition in conditions]


def ibs(g1, g2):
    """
    Returns the number of alleles identical by state between two genotypes
    Arguements: Two tuples
    Returns: an integer
    """
    g1, g2 = sorted(g1), sorted(g2)
    if g1 == g2:
        return 2
    g2s = set(g2)
    if g1[0] in g2s or g1[1] in g2s:
        return 1
    return 0


def consistent_dominant(record, strong=True):
    genotypes = [g for g in get_genotypes_from_record(record) if g]
    # Consistency with dominance requires two conditions
    # 1) Everyone must have (at least) one alternate call.
    # Since these are represented in the VCF files as numbers
    # greater than 0, we can just sum them.
    if not all(sum(g) > 0 for g in genotypes):
        return False
    # Everybody must have the same genotype
    if strong and not all(ibs(a, b) == 2 for a, b
                          in itertools.combinations(genotypes, 2)):
        return False
    return True


def consistent_recessive(record, strong=True, altcallsonly=True):
    genotypes = [g for g in get_genotypes_from_record(record) if g]
    # Consistency with recessive inheritance requires three conditions:
    # 1) All genotypes must be homozygous
    if not all(g[0] == g[1] for g in genotypes):
        return False
    # All alleles must be alt calls
    if altcallsonly and not all(0 not in g for g in genotypes):
        return False
    # Everybody should have the same genotype
    if strong and not all(ibs(a, b) > 0 for a, b
                          in itertools.combinations(genotypes, 2)):
        return False
    return True


def call_rate(record):
    ''' Returns the percent of nonmissing genotypes for a record '''
    genotypes = get_genotypes_from_record(record)
    return sum(1 for g in genotypes if g) / len(genotypes)


def smartopen(filename, mode='r'):
    """
    Seamlessly open gzipped and bzipped2 files. Use like regular open
    """
    if filename.endswith('.gz'):
        return gzip.GzipFile(filename, mode, compresslevel=5)
    elif filename.endswith('.bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename, mode)

# Program logic
###

conditions = []
condition_desc = []

if args.region:
    chr, start, stop = args.region
    try:
        start = int(start) if start != '-1' else float('-inf')
        stop = int(stop) if stop != '-1' else float('inf')
    except ValueError:
        print 'Error: bounds not numeric!'
        exit(1)

    def regioncheck(record):
        return record['CHROM'] == chr and (start <= int(record['POS']) <= stop)

    conditions.append(regioncheck)
    condition_desc.append('chr%s:%s-%s' % (chr, start, stop))

if args.minqual:
    conditions.append(lambda x: x['QUAL'] > args.minqual)
    condition_desc.append('Qual > %s' % args.minqual)

if args.qcfilter:
    conditions.append(lambda x: x['FILTER'] == 'PASS')
    condition_desc.append('FILTER=PASS')

if args.min_call_rate:
    conditions.append(lambda r: call_rate(r) >= args.min_call_rate)
    condition_desc.append('Call rate >= %s' % args.min_call_rate)

if args.ifilters:
    operators = {'gt': '>', 'gte': '>=', 'lt': '<', 'lte': '<=', 'eq': '=', 'neq': '!=',
                 'contains': 'contains', 'ncontains': 'does not contain',
                 'is': 'is', 'not': 'not'}

    def ifilter_describer(ifilter):
        ifilter[1] = operators[ifilter[1]]
        return ' '.join(ifilter)

    conditions.extend(parse_info_conditions(args.ifilters))
    condition_desc.extend([ifilter_describer(x) for x in args.ifilters])

if args.model:
    if args.model[0] == 'dom':
        conditions.append(consistent_dominant)
        condition_desc.append('Dominant')
    if args.model[0] == 'rec':
        conditions.append(consistent_recessive)
        condition_desc.append('Recessive')

print '%s filters in place' % len(conditions)
variants_passing_filters = [0] * len(conditions)
variants_passing_sequential = [0] * len(conditions)

with smartopen(args.file) as vcf, smartopen(args.outfile, 'w') as outfile:

    def outwrite(string):
        if args.outfile == os.devnull:
            pass
        else:
            outfile.write(string)
    header = None
    # Make sure this is a VCF file. The first line of a VCS file
    # has to have VCS in it
    firstline = vcf.readline()
    if 'VCF' not in firstline:
        print 'ERROR: File given is not a VCF file'
        exit(1)
    outwrite(firstline)
    # Skip the metadata
    for record in vcf:
        outwrite(record)
        if record.startswith('##'):
            l = record[2:]
            key, value = l.strip().split('=', 1)
            if key == 'fileDate':
                d = datetime.date.today()
                value = d.strftime('%Y%m%d')
                outwrite('##%s=%s\n' % (key, value))
        elif record.startswith('#'):
            # This is the header line, after this the data starts
            header = record[1:].strip().split()
            break
        else:
            print 'ERROR: poorly formed VCF skipped header line'
            exit(1)
    # Now start working with the actual data
    for variant_count, line in enumerate(vcf):
        l = line.strip().split()
        record = dict(zip(header, l))

        if args.rsonly and not record['ID'].startswith('rs'):
            continue

        if args.snvonly and not (len(record['ALT']) == len(record['REF']) == 1):
            continue


        record['QUAL'] = float(record['QUAL'])
        record['INFO'] = parse_info(record['INFO'])
        filters_passed = meets_conditions(record, conditions)
        # This might look a little odd because filters_passed is a list
        # of True/False values. Bools in python are just special integers and
        # you can do arithmetic with them:
        #    True == 1, False == 0. So True + 1 == (1 + 1) == 2.
        # Does it pass the filter (without regard for order)?
        for i, x in enumerate(filters_passed):
            variants_passing_filters[i] += x
        # How does it pass the filters (sequentially)?
        for i, x in enumerate(itertools.takewhile(lambda x: x, filters_passed)):
            variants_passing_sequential[i] += x
        if all(filters_passed):
            if not args.quiet:
                print 'Variant passed: %s %s %s' % (record['CHROM'],
                                                    record['POS'],
                                                    record['ID'])
            outwrite(line)
print
print 'Tested %d varaints' % variant_count
print '\t'.join(['Filter', 'Filter Description', 'Variants passing', 'Variants passing sequentially'])
for i, v in enumerate(itertools.izip(condition_desc, variants_passing_filters, variants_passing_sequential)):
    print '\t'.join([str(i+1)] + [str(x) for x in v])
