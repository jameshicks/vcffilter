#!/usr/bin/env python

### Imports
###
import argparse
import itertools
import os # For os.devnull

### Arguement parsing
###
parser = argparse.ArgumentParser(description='Filter variants in a VCF file')
parser.add_argument('-f','--file', required = True, metavar='vcffile',
                     help='VCF file for processing')
parser.add_argument('-o','--out', required = False, metavar='outfile',
                     help='File for output')
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
    if op == 'gt':
        err_on_nonnumeric(value)
        return lambda x: x['INFO'][on] > value
    elif op == 'gte':
        err_on_nonnumeric(value)
        return lambda x: x['INFO'][on] >= value
    elif op == 'lt':
        err_on_nonnumeric(value)
        return lambda x: x['INFO'][on] < value
    elif op == 'lte':
        err_on_nonnumeric(value)
        return lambda x: x['INFO'][on] <= value
    elif op == 'eq':
        return lambda x: x['INFO'][on] == value
    elif op == 'neq':
        return lambda x: x['INFO'][on] != value
    elif op == 'contains':
        return lambda x: value in x['INFO'][on] 
    elif op == 'ncontains':
        return lambda x: value not in x['INFO'][on]
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
# FIXME add qual filter option
if True:
    conditions.append(lambda x: x['FILTER'] == 'PASS')
if args.ifilters:
    conditions.extend(parse_info_conditions(args.ifilters))

variants_passing_filters = [0] * len(conditions)
variants_passing_sequential = [0] * len(conditions)

with open(args.file) as vcf:

    def outwrite(string):
        pass
    header = None
    # Make sure this is a VCF file. The first line of a VQS file
    # has to have VQS in it
    if 'VCF' not in vcf.readline():
        print 'ERROR: File given is not a VCF file'
        exit(1)
    # Skip the metadata
    for record in vcf:
        outwrite(record)
        if record.startswith('##'):
            continue
        elif record.startswith('#'):
            header = record[1:].strip().split()
        else:
            # The first line without a # is the start of the data
            break
    # Now start working with the actual data
    for variant_count, line in enumerate(vcf):
        l = line.strip().split()
        record = dict(zip(header,l)[0:9])
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
            outwrite(line.strip())
