#! /usr/bin/env python
import random, copy
import numpy as np
import math
import collections

class GAError(Exception):
    """
        For information as to why the exception was raised, see GAError.Message attribute
    """
    def __init__(self,strMessage):
        self.Message = strMessage

class Genom(object):
    #defines a genetic sequence
    def __init__(self,seq,sigma_mut,clamp):  
        if len(clamp) != len(seq):
            raise GAError("clamp and seq must have same length")
        self._clamp = clamp
        self._seq = seq
        self._fitness = []
        if len(sigma_mut) == len(seq):
            self._sigma_mut = sigma_mut
        else:
            raise GAError("sigma_mut must have same length as seq")

    def GetSequence(self):
        return copy.copy(self._seq)

    def GetFitness(self):
        return copy.copy(self._fitness)

    def SetFitness(self,fitness):
        self._fitness = fitness

    def Dominates(self,partner):
        if len(self._fitness) == len(partner._fitness):
            dom = [self._fitness[x] >= partner._fitness[x] for x in xrange(0,len(self._fitness))]
            # print 'a: ', self._fitness, ' b: ',partner._fitness, 'dom: ', dom , '\n'
            return (sum(dom) == (len(self.GetFitness())))
        raise GAError("Dominating(): fitness length mismatch")

    def Crossover(self,partner,p1,p2):
        if isinstance(partner, Genom):
            if len(partner.GetSequence()) == len(self.GetSequence()):
                if (len(self.GetSequence())-1) >= p2 > p1 >= 0:
                    new_seq = self.GetSequence()
                    new_seq[p1:p2]=partner.GetSequence()[p1:p2]
                    return Genom(new_seq,self._sigma_mut,self._clamp)
                else:
                    raise GAError("Crossover: p2 must be bigger than p1 and smaller than the seq. length")
            else: 
                raise GAError("Crossover:Both genes must have the same sequence length")
        else:
            raise GAError("Corssover: both partners must be of calss Genom")

    def Mutate(self):
        new_seq = self.GetSequence()
        for ind,gene in enumerate(new_seq):
            mut_gene = random.gauss(gene,self._sigma_mut[ind])
            new_seq[ind] = sorted([self._clamp[ind][0],self._clamp[1],mut_gene])[1]
        return Genom(new_seq,self._sigma_mut,self._clamp)



class Population(collections.MutableSequence):
    def  __init__(self,*members):
        self._genpool=list()
        self.extend(members)
        self._avg = []
        self._best = []
    def _check(self,v):
        if not isinstance(v,Genom):
            raise GAError("population members must be of class Genom")

    def __len__(self): return len(self._genpool)

    def __getitem__(self, i): return self._genpool[i]

    def __delitem__(self, i): del self._genpool[i]

    def __setitem__(self, i, v):
        self._check(v)
        self._genpool[i] = v

    def insert(self, i, v):
        self._check(v)
        self._genpool.insert(i, v)

    def ParetoDivide(self):

        """ returns a population of the paretto front members and a population of the rest """
        front = copy.deepcopy(self)
        rest = Population()
        k = 0
        for ind_a, mem_a in enumerate(front):
            for ind_b, mem_b in enumerate(front):
                if mem_a.GetSequence() != mem_b.GetSequence():
                    if mem_a.Dominates(mem_b):
                        rest.append(front.pop(ind_b))
        return front,rest

    def SelectNfittest(self,N):
        """ returns a population of the fittest: sampling pareto fronts until is colected N membrs """
        fittest = Population()
        front, rest = self.ParetoDivide()
        fittest.extend(random.sample(front,min(N,len(front))))
        while len(fittest)<N and len(rest):
            front, rest = rest.ParetoDivide()
            # print 'len front: ', len(front), 'len rest: ', len(rest),'\n'
            fittest.extend(random.sample(front,min(N-len(fittest),len(front))))
        return fittest
    def Average(self):
        avg = [0.0 for k in self[0].GetFitness()];
        for k in self: 
            avg = [avg[j] + k.GetFitness()[j] for j in xrange(len(avg))]
        avg = [avg[j]/len(self) for j in xrange(len(avg))]
        return avg
    def Best(self):
        best = [0.0 for k in self[0].GetFitness()];
        for score_key in xrange(len(self[0].GetFitness())):
            srt = sorted(self._genpool,key = lambda gene: gene.GetFitness()[score_key], reverse = True)
            best[score_key] = srt[0].GetFitness()[score_key]
        return copy.copy(best)

class GAtester(object):
    def __init__(self):
        pass
    def Evaluate(self,seq):
        return []


class GA(object):
    def __init__(self,tester):
        if isinstance(tester,GAtester):
            self._tester = tester
        else:
            raise GAError('tester must be of class GAtester')

    def _choose(self,prob,x):
        for k in xrange(len(prob)):
            p1 = sum(prob[0:k])
            p2 = sum(prob[0:k+1])
            if p1<=x<=p2:
                return k

    def _breed(self,population,Nchld):
        fronts = []
        front, rest = population.ParetoDivide()
        sorted_pop = front
        chld = Population()

        k = 1;
        while len(rest):
           front, rest = rest.ParetoDivide()
           sorted_pop.extend(front)
           fronts.extend([k for j in front])
           k +=1

        prob = [(fronts[-1] -j +1) for j in fronts]
        prob = [float(p)/sum(prob) for p in prob]

        while len(chld)<Nchld:            
            rnd1 = random.random()
            rnd2 = random.random()
            id1 = self._choose(prob,rnd1)
            id2 = self._choose(prob,rnd2)
            if id2!=id1:
                m1 = population[id1]
                m2 = population[id2]
                pt = sorted(random.sample(xrange(len(m1.GetSequence())),2))
                chld.append(m1.Crossover(m2,pt[0],pt[1]))

        return chld

    def Step(self,population,Ntop):
        if len(population)-2*Ntop<0:
            raise GAError("len(population)-2*Ntop<0")
        for member in population:
            member.SetFitness(self._tester.Evaluate(member.GetSequence()))
        top = population.SelectNfittest(Ntop)
        chld = self._breed(population,len(population)-2*Ntop)
        mut = Population()
        for k in top:
            mut.append(k.Mutate())
        nxt = Population()

        nxt.extend(top)
        nxt.extend(chld)
        nxt.extend(mut)
        return population, nxt

if __name__ == '__main__':
    #test script
    class test_eval(GAtester):
        def Evaluate(self,seq):
            fit = [math.cos(seq[1])**2*seq[2]**2 + abs(math.sin(seq[4]))*100 +seq[0]**2, seq[3]**2, abs(sum(seq))]
            return fit

    pop = Population()
    for k in xrange(100):
        gen = Genom([random.uniform(0,3.14) for j in xrange(10)],[0.2 for j in xrange(10)],[[0, 2*3.14] for j in xrange(10)])
        pop.append(gen)
    # print 'before'
    test = test_eval()
    ga = GA(test)
    prev,pop =  ga.Step(pop,50)
    print 'best: ', prev.Best(), 'avg :', prev.Average()
    for k in xrange(50):
        # print 'iteration: ', k
        prev,pop =  ga.Step(pop,20)
    print 'best: ', prev.Best(), 'avg :', prev.Average()


    # print 'population fitness: \n'
    # for k in pop:
    #     print k.GetFitness()
    # front, rest = pop.ParetoDivide()
    # print '\n paretto front: \n'
    # for k in front:
    #     print k.GetFitness()
    # print 'the rest: \n'
    # for k in rest:
    #     print k.GetFitness() 
    # print '\n average: ', pop.Average(),'best: ', pop.Best()

    # print '\n select 4: \n'
    # fit = pop.SelectNfittest(4)
    # for k in fit:
    #     print k.GetFitness()



