#!/usr/bin/env python

### Version Check
###
import sys
pyversion = sys.version_info[0:2]
if pyversion < (2, 7):
    print '%s requires Python 2.7 (you have %s)' % (sys.argv[0],
                                                    '.'.join(pyversion))
    exit(1)

### Imports
###
import itertools
import datetime
import os
from functools import wraps
import argparse

### Arguement parsing
###
parser = argparse.ArgumentParser(description='Filter variants in a VCF file')
parser.add_argument('-f','--file', required = True, metavar='vcffile',
                     help='VCF file for processing')
parser.add_argument('-o','--out', required = False, dest='outfile', metavar='outfile',
                     help='File for output', default=os.devnull)
parser.add_argument('-r','--region', required = False, nargs=3, help = 'Constrain to region') 
parser.add_argument('-q','--qual', required = False, type=float, dest='minqual',
                    help = 'Minimum quality score to include')
parser.add_argument('--no-qc', required = False, action='store_false', dest='qcfilter',
                    help="Do not filter on FILTER column == PASS")
parser.add_argument('--info_filter', dest='ifilters', nargs=3, action='append',
                     help='Filter on info string') 
args = parser.parse_args()

### Functions
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
            raise ValueError('Nonnumeric value (%s) for operator %s' % (value,op))
    def false_on_keyerror(f):
        @wraps(f)
        def wrapper(*args, **kwds):
            try:
                return f(*args,**kwds)
            except KeyError:
                return False
        return wrapper
    if op == 'gt':
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


def meets_conditions(inf, conditions):
    ''' Returns a list of whether variants meet the conditions for filtering '''
    return [condition(inf) for condition in conditions]

### Program logic
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
    condition_desc.append('chr%s:%s-%s' % (chr,start,stop))
if args.minqual:
    conditions.append(lambda x: x['QUAL'] > args.minqual)
    condition_desc.append('Qual > %s' % args.minqual)
if args.qcfilter:
    conditions.append(lambda x: x['FILTER'] == 'PASS')
    condition_desc.append('FILTER=PASS')
if args.ifilters:
    operators = {'gt':'>', 'gte':'>=', 'lt':'<', 'lte':'<=', 'eq': '=', 'neq': '!=',
                 'contains': 'contains', 'ncontains': 'does not contain'}
    def ifilter_describer(ifilter):
        ifilter[1] = operators[ifilter[1]]
        return ' '.join(ifilter)
    conditions.extend(parse_info_conditions(args.ifilters))
    condition_desc.extend([ifilter_describer(x) for x in args.ifilters])
print '%s filters in place' % len(conditions)
variants_passing_filters = [0] * len(conditions)
variants_passing_sequential = [0] * len(conditions)

with open(args.file) as vcf, open(args.outfile,'w') as outfile:

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
            l=record[2:]
            key, value = l.strip().split('=',1)
            if key == 'fileDate':
                d = datetime.date.today()
                value = d.strftime('%Y%m%d')
                outwrite('##%s=%s\n' % (key,value))
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
        record = dict(zip(header,l)[0:9])
        record['QUAL'] = float(record['QUAL'])
        record['INFO'] = parse_info(record['INFO'])
        filters_passed = meets_conditions(record, conditions)
        # This might look a little odd because filters_passed is a list
        # of True/False values. Bools in python are just special integers and
        # you can do arithmetic with them: 
        #    True == 1, False == 0. So True + 1 == (1 + 1) == 2.
        # Does it pass the filter (without regard for order)?
        for i,x in enumerate(filters_passed):
            variants_passing_filters[i] += x
        # How does it pass the filters (sequentially)?
        for i,x in enumerate(itertools.takewhile(lambda x: x, filters_passed)):
            variants_passing_sequential[i] += x 
        if all(filters_passed):
            print 'Variant passed: %s %s %s' % (record['CHROM'],
                                                   record['POS'],
                                                   record['ID'])
            outwrite(line)
print
print '\t'.join(['Filter','Filter Description','Variants passing','Variants passing sequentially'])
for i,v in enumerate(itertools.izip(condition_desc, variants_passing_filters, variants_passing_sequential)):
    print '\t'.join([str(i+1)] + [str(x) for x in v])
