# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 11:44:30 2023

@author: kgrun
"""

import json, argparse, glob
import numpy as np
import matplotlib.pyplot as plt



def setup_options():
    parser = argparse.ArgumentParser(description='create free spectrum')
    parser.add_argument('--chain', type=str, nargs='*', default=None)
    parser.add_argument('--nbins', type=int, default=None)
    parser.add_argument('--start', type=int, default=None)
    parser.add_argument('--out', type=str, default=None)
    parser.add_argument('--labels', type=str, nargs='*', default=None)
    return parser.parse_args()




def main():
    args = setup_options()
    print('Found {} chains to plot'.format(len(args.chain)))
    for ii,chain in enumerate(args.chain):
        print('... loading chain from {}'.format(chain))
        c = np.loadtxt(chain)
        burn = int(0.3*len(c[:,0]))
        chain_cut = c[burn:,:]
        print('... plotting chain from {}'.format(chain))
        plt.violinplot(chain_cut[:,args.start:args.start+args.nbins])
    
    plt.legend(args.labels, labelcolor=['blue','orange'])
    plt.savefig(args.out, bbox_inches='tight', dpi=400)
    plt.clf()
    print('Done.')
    return None

if __name__ == '__main__':
    main()


