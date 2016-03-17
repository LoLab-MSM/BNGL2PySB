# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:01:02 2015

@author: James C. Pino
"""
import re
import sys
#name = sys.argv[1]
name = 'Tests/p53_mdm2.bngl'
File = open(name,'r')
data = File.readlines()
File.close()
global parameters_all
parameters_all = []
def convert():
    cont = False
    params = []
    parameters = ''
    for line in data:
        if line.startswith('begin parameters'):
            cont = True
        elif line.startswith('end parameters'):        
            cont = False
        elif cont == True:
            name,value = line[:-1].split(' ')
            name = name.strip()
            #print "Parameter('%s' , %s )" % (name,value)
            parameters_all.append(name)
            parameters+="Parameter('%s', %s)\n" % (name,value)
            params.append(name)

    cont = False
    species = ''
    spec = []
    for line in data:
        if line.startswith('begin molecule types'):
            cont = True
        elif line.startswith('end molecule types'):        
            cont = False
        elif cont == True:
            #print line
            name , tail =   line[1:-2].split('(')
            name = name.strip()
            tail = tail.strip()
            spec.append(name)
            #print name, tail
            sites = []
            states = {}
            for each in tail.split(','):
                if each != '':
                    sites.append(each)
                if len(each.split('~')) >1:     
                    sites.pop()
                    sites.append(each.split('~')[0])
                    for i in range(1,len(each.split('~'))):
                        states[each.split('~')[0]] = each.split('~')[i]
                    #states[each.split('~')[0]] = each.split('~')
            if len(sites) == 0:
                #print "Monomer('%s')" % (name)
                species+="Monomer('%s')\n" % (name)
            elif len(states) >0:
                #print "Monomer('%s' , %s , %s)" % (name,sites,states)
                species +="Monomer('%s', %s, %s)\n" % (name,sites,states)
            else:
                #print "Monomer('%s' , %s )" % (name,sites)
                species +="Monomer('%s', %s)\n" % (name,sites)
    # gather rules
    cont = False
    rules = ''
    for line in data:
        if line.startswith('begin reaction rules'):
            cont = True
        elif line.startswith('end reaction rules'):        
            cont = False
        elif cont == True:
            name , rule =   line[1:-1].split(':')
            name = name.strip()
            rule = rule.strip()
            reactants ,products= rule.split('->')
            reactants = reactants.replace('0 ','None ')
            products = products.replace('0 ', 'None ')
            products = products.replace('~','=')
            products = products.replace('!','=')
            products = products.replace('.','%')
            reactants = reactants.replace('~','=')
            reactants = reactants.replace('!','=')
            reactants = reactants.replace('.','%')
            rate = products.split(' ')[-1]
            products = products.rstrip(rate)
            reactants = reactants.split('@')[0]
            products = products.split('@')[0]
            reactants = reactants.replace('Source()', 'None')
            products = products.replace('Sink()', 'None')
            #print 'Rule("%s" , %s >> %s, %s)' % (name,reactants,products,rate)
            rules+= 'Rule("%s", %s >> %s, %s)\n' % (name,reactants,products,rate.strip('()'))
    # gather     observables
    cont = False
    observables = ''
    observales_dict = {}
    for line in data:
        if line.startswith('begin observables'):
            cont = True
        elif line.startswith('end observables'):        
            cont = False
        elif cont == True:
            obs = line[:-1].split('#',1)[0]
            obs = obs.strip()
            obs = obs.lstrip('Species')
            obs = obs.strip()
            name,obs = obs.split(' ')
            observales_dict[name] = "obs_%s"%name
            #print 'Observable("obs_%s" , %s )' % (name,obs)
            observables += 'Observable("obs_%s", %s)\n' % (name,obs)
    print observables
    cont = False
    initial_conditions = ''
    for line in data:
        if line.startswith('begin seed species'):
            cont = True
        elif line.startswith('end seed species'):        
            cont = False
        elif cont == True:
            compartment , name =   line[1:-1].split(':')
            #print compartment,name
            name , comment = name.split('#',1)
            #print name, comment
            name = name.rstrip(' ')
            name = name.lstrip('$')
            name,value = name.split(' ')
            #print 'Initial(%s ,Parameter("%s_0",%s) )' % (name,name.rstrip('()'),value)
            initial_conditions += 'Initial(%s, Parameter("%s_0",%s))\n' % (name,name.rstrip('()'),value)

    cont = False
    functions = ''
    for line in data:
        if line.startswith('begin functions'):
            cont = True
        elif line.startswith('end functions'):        
            cont = False
        elif cont:
            name, rate = line[1:-1].split('=')
            rate = rate.strip()
            name = name.strip().strip("()")
            #print 'Expression("%s",%s )' %(name,rate)
            for i in spec:
                rate = rate.replace(" "+i," "+observales_dict[i])
            functions += 'Expression("%s", %s )\n' % (name, rate)

    return parameters+species+observables+initial_conditions+functions+rules



params = convert()
#print params

#print func
with open('p53_2.py','w') as output:
    output.write('from pysb import *\nModel()\n')
    output.write(params)

from p53_2 import model
from pysb.bng import generate_network
generate_network(model,cleanup=False)
