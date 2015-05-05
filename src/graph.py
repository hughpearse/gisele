#!/bin/python

import sys, getopt, os, pdb
import logging
logging.getLogger("scapy.runtime").setLevel(logging.ERROR)

import random
from array import *

import numpy as np
from numpy import *
from numpy.linalg import *
from numpy.random import *

import binascii

from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA

import scipy as scipy
from scipy import spatial

#mean shift imports
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets.samples_generator import make_blobs

#plot imports
import pylab as pl
from itertools import cycle

#needleman wunsch imports
from math import *
import copy

try:
    from scapy.all import * #Required for Scapy 2.0 and above
except:
    from scapy import * #Scapy 1.0

from scapy.utils import rdpcap
from scapy.layers.inet import IP,UDP,TCP
from scapy.packet import Raw

#Number of packets to sample
packetsToSample = 30

#Baum-Welch Variables
aFirst = np.array([0.25,0.25,0.25,0.25])
aLast = np.array([0.25,0.25,0.25,0.25])

#Hidden Markov Model - Forward Algorithm
#
# source: http://seat.massey.ac.nz/personal/s.r.marsland/Code/15/HMM.py
#
def HMMfwd(a,b,obs):

    nStates = shape(b)[0]
    T = shape(obs)[0]

    alpha = zeros((nStates,T))

    alpha[:,0] = aFirst*b[:,obs[0]]

    for t in range(1,T):
        for s in range(nStates):
            alpha[s,t] = b[s,obs[t]] * sum(alpha[:,t-1] * a[:,s])

    #alpha[q,T] = sum(alpha[:,T-1]*aLast[:,q])
    #print max(alpha[:,T-1])
    return alpha

#Hidden Markov Model - Backward Algorithm
#
# source: http://seat.massey.ac.nz/personal/s.r.marsland/Code/15/HMM.py
#
def HMMbwd(a,b,obs):

    nStates = shape(b)[0]
    T = shape(obs)[0]

    beta = zeros((nStates,T))

    beta[:,T-1] = aLast

    for t in range(T-2,0,-1):
        for s in range(nStates):
            beta[s,t] = b[s,obs[t+1]] * sum(beta[:,t+1] * a[:,s])

    beta[:,0] = b[:,obs[0]] * sum(beta[:,1] * aFirst)
    return beta

#Baum-Welch Algorithm
#
# source: http://seat.massey.ac.nz/personal/s.r.marsland/Code/15/HMM.py
#
def BaumWelch(obs,nStates):

    T = shape(obs)[0]
    a = random.rand(nStates,nStates)
    b = random.rand(nStates,T)
    olda = zeros((nStates,nStates)) 
    oldb = zeros((nStates,T)) 
    maxCount = 50
    tolerance = 1e-5

    count = 0
    while (abs(a-olda)).max() > tolerance and (abs(b-oldb)).max() > tolerance and count < maxCount:
        # E-step

        alpha = HMMfwd(a,b,obs)
        beta = HMMbwd(a,b,obs)
        gamma = zeros((nStates,nStates,T))

        for t in range(T-1):
            for s in range(nStates):
                gamma[:,s,t] = alpha[:,t] * a[:,s] * b[s,obs[t+1]] * beta[s,t+1] / max(alpha[:,T-1])
    
        # M-step
        olda = a.copy()
        oldb = b.copy()

        for i in range(nStates):
            for j in range(nStates):
                a[i,j] = sum(gamma[i,j,:])/sum(sum(gamma[i,:,:]))

        for o in range(max(obs)):
            for j in range(nStates):
                places = (obs==o).nonzero()
                tally = sum(gamma[j,:,:],axis=0)
                b[j,o] = sum(tally[places])/sum(sum(gamma[j,:,:]))
                #print b[j,o], sum(gamma[j,places])/sum(gamma[j,:])
    
        count += 1
    print count
    return a,b

#Smith Waterman sequence alignment
#
# source: https://gist.github.com/Puriney/6286305
#
def sequencealignment(strA, strB):
    totalDifference = 0;
    totalMatch = 0;
    strG = ''.join(strA)
    strR = ''.join(strB)
    matrix = []
    path = []
    row = len(strR) 
    col = len(strG)
    strR = "^"+strR
    strG = "^"+strG
    for i in range(row+1):
        matrix.append([0]*(col+1))
        path.append(["N"]*(col+1))
    
    def print_matrix(matrix):
        print '\t'+('\t'.join(map(str,list(strG))))
        i = 0
        for line in matrix:
            print strR[i]+"\t"+('\t'.join(map(str,line)))
            i +=1
    
    # print_matrix(matrix)
    indelValue = -1
    matchValue = 2
    for i in range(1,row+1):
        for j in range(1,col+1):
            # penalty map
            from_left = matrix[i][j-1] + indelValue
            from_top = matrix[i-1][j] + indelValue
            if strR[i]==strG[j]:
                from_diag = matrix[i-1][j-1] + matchValue
            else:
                from_diag = matrix[i-1][j-1] + indelValue
    
            matrix[i][j]= max(from_left,from_top,from_diag)
            # path map
            if matrix[i][j]==from_left:
                path[i][j]="-"
                totalDifference += indelValue
            elif matrix[i][j]==from_top:
                path[i][j] = "|"
                totalDifference += indelValue
            elif matrix[i][j] == from_diag:
                path[i][j] = "M"
                totalMatch += matchValue
            else:
                pass
    
            if matrix[i][j]<0:
                matrix[i][j]=0
    pass
    #print_matrix(matrix)
    #print
    #print_matrix(path)
    
    iRow = len(matrix)-1
    jCol = len(matrix[0])-1
    
    while iRow>=0:
        maxPnltyVaue = max(matrix[iRow])
        while jCol>=0:
            if matrix[iRow][jCol]==maxPnltyVaue:
                ipair = iRow
                jpair = jCol
                reportR=[]
                reportG=[]
                while (1):
                    if ipair ==0 and jpair==0:
                        break
                    # else:
                    if path[ipair][jpair]=="M":
                        reportR.append(strR[ipair])
                        reportG.append(strG[jpair])
                        ipair -= 1
                        jpair -= 1
                    elif path[ipair][jpair]=="-":
                        # reportR.append(strR[ipair])
                        reportR.append("-")
                        reportG.append(strG[jpair])
                        # ipair -= 1
                        jpair -= 1
                    elif path[ipair][jpair]=="|":
                        reportR.append(strR[ipair])
                        reportG.append("-")
                        ipair -= 1
                    elif path[ipair][jpair]=="N":
                        if ipair > 0:
                            reportR.append(strR[ipair])
                            ipair -=1
                        if jpair > 0:
                            reportG.append(strG[jpair])
                            jpair -= 1
                s1 = "".join(reportR)
                s2 = "".join(reportG)
                
            jCol -=1
        iRow -=1
    
    #Custom code to construct GMS
    s3 = []
    i=0
    while i < max(len(s1), len(s2)):
        try:
            if(s1[i] == s2[i]):
                s3.append(s1[i])
            else:
                s3.append('-')
        except IndexError:
            s3.append('-')
        i+=1
        
    return list("".join(reversed(s3))), abs(totalMatch), abs(totalDifference), list("".join(reversed(s1))), list("".join(reversed(s2)))
                                                                                                                       
#Multi-Dimensional Scaling
#
# Source: http://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#example-manifold-plot-mds-py
#
def mds(d, dimensions = 2):
        """
        Multidimensional Scaling - Given a matrix of interpoint distances,
        find a set of low dimensional points that have similar interpoint
        distances.
        """

        (n,n) = d.shape
        E = (-0.5 * d**2)

        # Use mat to get column and row means to act as column and row means.
        Er = mat(mean(E,1))
        Es = mat(mean(E,0))

        # From Principles of Multivariate Analysis: A User's Perspective (page 107).
        F = array(E - transpose(Er) - Es + mean(E))

        [U, S, V] = svd(F)

        Y = U * sqrt(S)

        return (Y[:,0:dimensions], S)

#Longest Common Substring - Smith Waterman
#
# Source: http://hipersayanx.blogspot.co.uk/2012/08/longest-common-subsequence-in-python.html
#
def lcs(a=[], b=[]):
    ja = -1
    jb = -1
    n = 0
 
    if a == [] or b == []:
        return ja, jb, n
 
    l = len(a) + len(b) - 1
    ia = len(a) - 1
    ib = 0
    s = 1
 
    for k in range(l):
        nCur = 0
 
        for r in range(s):
            if a[ia + r] == b[ib + r]:
                nCur += 1
 
                if nCur > n:
                    ja = ia + r - nCur + 1
                    jb = ib + r - nCur + 1
                    n = nCur
            else:
                nCur = 0
 
        if k < min(len(a), len(b)) - 1:
            ia -= 1
            s += 1
        elif k > l - min(len(a), len(b)) - 1:
            ib += 1
            s -= 1
        elif ia > 0:
            ia -= 1
        else:
            ib += 1
 
    #return ja, jb, n
    return n

def main(argv):
        inputfile = ''
        outputfile = ''
        hookfile = ''
        try:
                opts, args = getopt.getopt(argv,"hi:o:t:",["ifile=","ofile=","transformation="])
        except getopt.GetoptError:
                print sys.argv[0] + ' -i inputfile -o outputfile -t [hook script]'
                print sys.argv[0] + ' -i sample.pcap -o result.pcap -t example'
                sys.exit(2)
        for opt, arg in opts:
                if opt == '-h':
                        print 'test.py -i <inputfile> -o <outputfile>'
                        sys.exit()
                elif opt in ("-i", "--ifile"):
                        inputfile = arg
                elif opt in ("-o", "--ofile"):
                        outputfile = arg
                elif opt in ("-t", "--transformation"):
                        hookfile = arg
        
        #copy packet data to dictionary object
        pktdict={}
        pkts=rdpcap(inputfile)
        i=0
        for pkt in pkts:
            try:
                my_array = []
                if pkt.haslayer(TCP):
                    for d in str(pkt.getlayer(TCP).payload):
                        my_array.append(d)
                if pkt.haslayer(UDP):
                    for d in str(pkt.getlayer(UDP).payload):
                        my_array.append(d)

                #reverse packet for backtrace on needleman wunch
                pktdict[i] = list("".join(reversed(my_array)))
                i=i+1
            except:
                raise

        #Create distance matrix
        dictSize = len(pktdict)
        diffMatrix = zeros((packetsToSample,packetsToSample))
        
        x=0
        while x < packetsToSample:
                y=0
                print ""
                print "Packet " + str(x) + ": " + str(pktdict[x])
                while y < packetsToSample:
                        #calculate common substring length between packets
                        #similarity = lcs(pktdict[x], pktdict[y])
                        gms, similarity, distance, alignedseq1Discard, alignedseq2Discard = sequencealignment(pktdict[x], pktdict[y])
                        #distance = 1 - (similarity + 1)/2
                        print "Packet " + str(x) + " similarity to packet " + str(y) + " = " + str(similarity)
                        print "Packet " + str(x) + " distance from packet " + str(y) + " = " + str(distance)
                        
                        #assign value to symmetrically opposite cells
                        #as Smith-Waterman score follows triangle equality rule
                        diffMatrix[x][y]=distance
                        diffMatrix[y][x]=distance
                        y=y+1
                x=x+1
        
        print " "
        print "Distance Matrix:"
        print diffMatrix
        print ""
        
        #Multi-Dimensional Scaling from distances to XY points
        #
        # Source: http://scikit-learn.org/stable/auto_examples/manifold/plot_mds.html#example-manifold-plot-mds-py
        #
        seed = np.random.RandomState(seed=3)
        mds = manifold.MDS(n_components=2, max_iter=1000, eps=0.8, random_state=seed, dissimilarity="precomputed", n_jobs=1)
        pos = mds.fit(diffMatrix).embedding_
        pos = manifold.MDS(dissimilarity="precomputed").fit_transform(diffMatrix)
        pos *= np.sqrt(100000) / np.sqrt((pos ** 2).sum())
        clf = PCA(n_components=2)
        pos = clf.fit_transform(pos)
        
        #Display distance matrix
        print "Coordinates of plotted packets: "
        for p in pos:
            print p.astype(int)
        
        #Calculate number of clusters
        #
        # Source: http://scikit-learn.org/stable/auto_examples/cluster/plot_mean_shift.html
        #
        print ""
        bandwidth = estimate_bandwidth(pos, quantile=0.2, n_samples=500)
        ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
        ms.fit(pos)
        labels = ms.labels_
        cluster_centers = ms.cluster_centers_
        labels_unique = np.unique(labels)
        n_clusters_ = len(labels_unique)
        print("Estimated number of clusters (k): %d" % n_clusters_)
        
        #Plot on graph
        pl.figure(1)
        pl.clf()
        colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
        for k, col in zip(range(n_clusters_), colors):
            my_members = labels == k
            cluster_center = cluster_centers[k]
            pl.plot(pos[my_members, 0], pos[my_members, 1], col + '.')
            print ""
            print "Cluster: " + str(k)
            #print str(my_members)
            
            #Create GMS for cluster using Needleman-Wunch
            #
            #This section is not part of the clustering code
            #
            clusterPackets = [];
            origionalClusterPackets = [];
            offset = 0;
            #extract packets from each cluster
            for val in my_members:
                if(str(val) == "True"):
                    clusterPackets.append(pktdict[offset]);
                offset += 1;
            origionalClusterPackets = copy.deepcopy(clusterPackets);
            
            print 'Compressing GMS .',
            #compress all GMS pairs to single GMS for the cluster
            while len(clusterPackets) > 1:
                print '.',
                gmsList1 = [];
                #calculate generic message sequence for each pair of messages
                for i in xrange(len(clusterPackets) - 1):
                    current_item, next_item = clusterPackets[i], clusterPackets[i + 1]
                    gms, totalMatch, totalDifference, alignedseq1Discard, alignedseq2Discard = sequencealignment(current_item, next_item)
                    gmsList1.append(gms)
                clusterPackets = copy.deepcopy(gmsList1)
            print ""
            gmspkt = list(reversed(clusterPackets[0]))
            #gmsbin = str("".join(gmspkt))
            
            #compress all substitution characters to a single character
            beforeGmsLen = len(gmspkt)+1
            afterGmsLen = len(gmspkt)
            while(beforeGmsLen > afterGmsLen):
                beforeGmsLen = len(gmspkt)
                for i in xrange(len(gmspkt) - 1, 0, -1):
                    if(gmspkt[i] == "-" and gmspkt[i-1] == "-"):
                        del gmspkt[i]
                afterGmsLen = len(gmspkt)
            print str(gmspkt)
            #print list(reversed(gmspkt))
            print ""
            #enumerate ngrams in variable data
            clusterTokens = [];
            for clusPkt in origionalClusterPackets:
                gmsDiscard, similarityDiscard, distanceDiscard, alignedseq1Keep, alignedseq2Keep = sequencealignment(clusPkt, list(reversed(gmspkt)))
                pktAlignedGMS1 = list(reversed(alignedseq1Keep))
                GMSAlignedData = list(reversed(alignedseq2Keep))
                tmpData = copy.deepcopy(GMSAlignedData)
                packettokens = [];
                packettoken = [];
                gmsoffset = 0
                while gmsoffset < len(GMSAlignedData):
                    if pktAlignedGMS1[gmsoffset]:
                        if(pktAlignedGMS1[gmsoffset] != "-"):
                            tmpData[gmsoffset] = "-"
                    gmsoffset+=1;
                packettokens = copy.deepcopy(tmpData)
                
                splittoken = [];
                splittokens = [];
                gmsoffset = 0
                while gmsoffset < len(packettokens):
                    if packettokens[gmsoffset]:
                        if(packettokens[gmsoffset] != "-"):
                            splittoken.append(copy.deepcopy(packettokens[gmsoffset]))
                        if(gmsoffset+1 < len(packettokens)):
                            if(packettokens[gmsoffset+1] == "-"):
                                if(len(splittoken) > 0):
                                    #TODO:
                                    # having problems with passing by reference
                                    # the beginning of the list vanishes.
                                    splittokens.append(copy.deepcopy(splittoken))
                                    del splittoken[:]
                    gmsoffset+=1
                
                clusterTokens.append(splittokens)
                
                print
                print "GMS: " + str(pktAlignedGMS1)
                print "Data: " + str(GMSAlignedData)
                print "Masked Data: " + str(tmpData)
                print "Tokens: " + str(packettokens)
                print "Split Tokens: " + str(splittokens)
                print
                
            print ""
            
            #infer token data type
            #integer, float, character, string
            fieldtype = "Blob"
            for tokens in clusterTokens:
                for token in tokens:
                    singleToken = ''.join(token)
                    newToken = str(singleToken)
                    if newToken == 'True' or newToken == 'False':
                        fieldtype = "Flag"
                    else:
                        try:
                            int(newToken)
                            fieldtype = "Number"
                        except ValueError:
                            try:
                                float(newToken)
                                fieldtype = "Number"
                            except ValueError:
                                fieldtype = "Blob"
            
            chunkOffset = 0;
            staticFieldBuff = [];
            fieldSwitch = 0;
            fieldLength = 1;
            print '<DataModel name="cluster' + str(k) + '">'
            while chunkOffset < len(gmspkt):
                
                #print conditions
                if gmspkt[chunkOffset]:
                    if(gmspkt[chunkOffset] != "-"):
                        staticFieldBuff.append(gmspkt[chunkOffset])
                    if(chunkOffset+1 < len(gmspkt)):
                        if(gmspkt[chunkOffset] != "-" and gmspkt[chunkOffset+1] == "-"):
                            #print '<Blob valueType="hex" value="' + str(staticFieldBuff).replace("[", "").replace("]", "").replace("'", "").replace("\\x", "") + '" mutable="false"/>'
                            print '<Blob valueType="hex" value="',
                            for c in staticFieldBuff:
                                print binascii.hexlify(c),
                            print '" mutable="false"/>'
                            del staticFieldBuff[:]
                           
                    if(chunkOffset == 0 and gmspkt[chunkOffset] == "-"):
                        fieldSwitch = 1
                    if(chunkOffset-1 > 0):
                        if(gmspkt[chunkOffset] == "-" and gmspkt[chunkOffset-1] != "-"):
                            fieldSwitch = 1
                    
                if(fieldSwitch == 1):
                    print '<' + str(fieldtype) + ' mutable="true"/>'
                    fieldLength=0
                    fieldSwitch = 0
                    
                chunkOffset+=1
                fieldLength+=1
                
            print '</DataModel>'
            
            
            pl.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col, markeredgecolor='k', markersize=8)
        pl.title('Estimated number of clusters: %d' % n_clusters_)
        pl.show()
        
        
        

#pdb.set_trace()
total = len(sys.argv)
main(sys.argv[1:])

