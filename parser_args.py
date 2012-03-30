#!/usr/bin/env python

"""
A module utilizing argparse which will handle command line input
based on a specification dictionary provided by the app. The dictionary is
indexed by a descriptive keyword and contains a Tuple with the following
organization:

(           name or flags -         Either a name or a list of option strings, e.g. foo or -f, --foo.
            action -                The basic type of action to be taken when this argument 
                                        is encountered at the command line.
            nargs -                 The number of command-line arguments that should be consumed.
            const -                 A constant value required by some action and nargs selections.
            default -               The value produced if the argument is absent from the command line.
            type -                  The type to which the command-line argument should be converted.
            choices -               A container of the allowable values for the argument.
            required -              Whether or not the command-line option may be omitted (optionals only).
            help -                  A brief description of what the argument does.
            metavar -               A name for the argument in usage messages.
            dest -                  The name of the attribute to be added to the object returned by parse_args().

)

"""
import argparse

parser=argparse.ArgumentParser(description="A program to set parameters for parsing BLAST output. Default is 0.")

parser.add_argument('-v', "--verbose",
                    help="verbose output", )
                
parser.add_argument('-i', "--identity", type=int, default=100,
                    help= "Minimum Percent Identity", dest='IDENTITY')
                
parser.add_argument("-l", "--length", type=int, default=0,
                    help= "Minimum Length", dest = 'LENGTH')
                 
parser.add_argument("-e", "--e_value", type=float, default=0.00, 
                    help= "Minimum E-value", dest ='E_VALUE')
                 
parser.add_argument("-b", "--bit_score", type=float, default=0.00, 
                    help= "Minimum bit score", dest='BIT_SCORE'  )

results=parser.parse_args()



'''
Created on Mar 7, 2012

@author: scro934
'''
'''
Created on Mar 8, 2012

@author: scro934
'''
