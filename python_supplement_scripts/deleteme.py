# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 10:53:50 2024

@author: kgrun
"""


import argparse



def set_up_global_options():
    parser = argparse.ArgumentParser(description='Generate sensitivity curves for PTA datasets.')

    parser.add_argument('--chain', type=str, nargs='*', default=None,
                        help='name of the chains')
    return parser.parse_args()



args = set_up_global_options()

print(len(args.chain), args.chain[0])