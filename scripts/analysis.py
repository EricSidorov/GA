#! /usr/bin/env python
import roslib; roslib.load_manifest('GA')
import os
import yaml
import sys
from GA import *
import matplotlib.pyplot as plt
def plot_best_avg(pop_dir):
    done = False
    k = 0;
    avg = [];
    best = [];
    ks = [];
    while not done:
        file_path = os.path.join(pop_dir,str(k)+'.yaml')
        try:
            pop = yaml.load(file(file_path))
            best.append(pop.Best())
            avg.append(pop.Average())
            ks.append(k)
            k+=1
        except IOError:
            done = True

    N = len(best[0])
    plt.figure(1)
    for k in xrange(N):
        sub = 100*N + 10 + (k+1)
        average = [a[k] for a in avg]
        bst = [b[k] for b in best]
        plt.subplot(sub)
        plt.plot(ks,bst, label='best')
        plt.plot(ks,average,label='average')
        plt.legend()
    plt.show()
    
def scatter(file_path):
    pop = yaml.load(file(file_path))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for member in pop:
        ax.scatter(member.GetFitness()[0], member.GetFitness()[1])
    plt.show()

if __name__ == '__main__':
    # plot_best_avg(sys.argv[1])
    scatter(sys.argv[1])    



