# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:01:02 2015

@author: James C. Pino
"""
import re
import sys
#name = sys.argv[1]
name = 'Tests/gene_expr_func.bngl'
File = open(name,'r')
data = File.read().splitlines()
File.close()
global parameters_all
parameters_all = []


class Converter():
    def __init__(self):
        self.parameters = []
        self.observables = []
        self.functions = []
        self.rules = []


def extract_params():
    cont = False
    params = []
    parameters = ''
    for line in data:
        line = line.strip()
        if line.startswith('begin parameters'):
            cont = True
        elif line.startswith('end parameters'):
            cont = False
        elif cont:
            if line.startswith('#'):
                continue
            name, value = line.split(' ',1)
            if len(value.split('#')) > 1:
                value = value.split('#')[0]
            name = name.strip()
            # print "Parameter('%s' , %s )" % (name,value)
            parameters_all.append(name)
            parameters += "Parameter('%s', %s)\n" % (name, value)
            params.append(name)
    return parameters

def extract_molecule_types():
    cont = False
    species = ''
    spec = []
    for line in data:
        if line.startswith('begin molecule types'):
            cont = True
        elif line.startswith('end molecule types'):
            cont = False
        elif cont:
            # print line
            name, tail = line.split('(')
            name = name.strip()
            tail = tail.strip()
            spec.append(name)
            # print name, tail
            sites = []
            states = {}
            for each in tail.split(','):
                if each != '':
                    sites.append(each)
                if len(each.split('~')) > 1:
                    sites.pop()
                    sites.append(each.split('~')[0])
                    for i in range(1, len(each.split('~'))):
                        states[each.split('~')[0]] = each.split('~')[i]
                        # states[each.split('~')[0]] = each.split('~')
            if len(sites) == 0:
                # print "Monomer('%s')" % (name)
                species += "Monomer('%s')\n" % (name)
            elif len(states) > 0:
                # print "Monomer('%s' , %s , %s)" % (name,sites,states)
                species += "Monomer('%s', %s, %s)\n" % (name, sites, states)
            else:
                # print "Monomer('%s' , %s )" % (name,sites)
                species += "Monomer('%s', %s)\n" % (name, sites)
    return species , spec


def convert_specie(specie):
    specie = specie.replace('0 ', 'None ')
    specie = specie.replace('~', '=')
    specie = specie.replace('!', '=')
    specie = specie.replace('.', '%')
    return specie

def extract_rules():
    # gather rules
    cont = False
    rules = ''
    count = 0
    for line in data:
        line = line.strip()
        if line.startswith('#') or line == '':
            continue
        if line.startswith('begin reaction rules'):
            cont = True
        elif line.startswith('end reaction rules'):
            cont = False
        elif cont:
            try:
                name, rule = line.split(':')
            except:
                rule = line
                name = "rule_%s" % count
                count += 1
            name = name.strip()
            rule = rule.strip()
            if len(rule.split('<->')) == 2:
                reactants, products = rule.split('<->')
            else:
                reactants, products = rule.split('->')
            reactants = convert_specie(reactants)
            products = convert_specie(products)

            rate = products.split(' ')[-1]
            products = products.rstrip(rate)
            reactants = reactants.split('@')[0]
            products = products.split('@')[0]
            reactants = reactants.replace('Source()', 'None')
            products = products.replace('Sink()', 'None')
            # print 'Rule("%s" , %s >> %s, %s)' % (name,reactants,products,rate)
            rules += 'Rule("%s", %s >> %s, %s)\n' % (name, reactants, products, rate.strip('()'))
    return rules

def extract_obs():
    # gather     observables
    cont = False
    observables = ''
    observales_dict = {}
    for line in data:
        line = line.strip()
        if line.startswith('begin observables'):
            cont = True
        elif line.startswith('end observables'):
            cont = False
        elif cont:
            if line.startswith('#'):
                continue
            obs = line.split('#', 1)[0]
            obs = obs.strip()
            obs = obs.lstrip('Species')
            obs = obs.lstrip('Molecules')
            obs = obs.strip()
            obs = convert_specie(obs)
            name, obs = obs.split()
            observales_dict[name] = "obs_%s" % name
            # print 'Observable("obs_%s" , %s )' % (name,obs)
            observables += 'Observable("obs_%s", %s)\n' % (name, obs)
    return observables , observales_dict

def extract_seed_species():
    cont = False
    initial_conditions = ''
    for line in data:
        line = line.strip()
        if line.startswith('begin seed species'):
            cont = True
        elif line.startswith('end seed species'):
            cont = False
        elif cont:
            print line
            if line.startswith('#'):
                continue
            if len(line.split(':')) == 2:
                compartment, name = line.split(':')
            else:
                name = line
            # print compartment,name
            if len(name.split(':')) == 2:
                name, comment = name.split('#', 1)
            name = name.rstrip(' ')
            name = name.lstrip('$')
            name, value = name.split()
            name = convert_specie(name)
            # print 'Initial(%s ,Parameter("%s_0",%s) )' % (name,name.rstrip('()'),value)
            initial_conditions += 'Initial(%s, Parameter("%s_0",%s))\n' % (name, name.rstrip('()'), value)
    return initial_conditions


def extract_functions(observales_dict,spec):
    cont = False
    functions = ''
    for line in data:
        line = line.strip()
        if line.startswith('begin functions'):
            cont = True
        elif line.startswith('end functions'):
            cont = False
        elif cont:
            if line.startswith('#'):
                continue
            print line
            name, rate = line[1:-1].split('=')
            rate = rate.strip()
            name = name.strip().strip("()")
            # print 'Expression("%s",%s )' %(name,rate)
            for i in observales_dict.keys():
                print i,rate
                rate = rate.replace(" " + i, " " + observales_dict[i])
            functions += 'Expression("%s", %s )\n' % (name, rate)
    return functions


def convert():
    parameters = extract_params()
    species,spec = extract_molecule_types()
    print spec
    rules = extract_rules()
    observables, observales_dict = extract_obs()
    print observales_dict
    initial_conditions = extract_seed_species()
    functions = extract_functions(observales_dict,spec)

    return parameters+species+observables+initial_conditions+functions+rules



params = convert()
#print params

#print func
with open('gene_expr_func.py','w') as output:
    output.write('from pysb import *\nModel()\n')
    output.write(params)

from gene_expr_func import model
from pysb.bng import generate_network
generate_network(model,cleanup=False)
