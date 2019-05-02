#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 22:13:41 2019

@author: preranaparthasarathy
"""

#Diwaan

from igraph import Graph, mean
from igraph import summary
from igraph import plot
import igraph
import time
import sys
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
import networkx as nx
import heapq
import operator

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import scipy.cluster.hierarchy as sch
#import louvain



#This function returns a sanitized protein with spaces and digits removed, and case set to lower		
def sanitize(string):
	string=re.sub(r'[^\w]','',string)	
	string=re.sub(r'[\d_]','',string)		
	string=string.lower()
	return string

#returns a list as a string with sep separators.
def listtostring(list, sep):
	string=""
	if len(list)>0:
		string+=list[0]#.decode('utf-8')
	for item in list[1:]:
		string+=sep+item#.decode('utf-8')
	return string

#function to read fasta files
def readFasta(file, keepFile="NA"):
	import re
	sequences=[]
	kp=set()
	with open(file,'r') as input:
		if keepFile!="NA":
			with open(keepFile,'r') as keep:
				for line in keep:
					kp.add(line[:-1])
		#print kp
		for line in input:
			if line[0]!=">":
				if len(kp)>0:
					if sequences[-1][0] in kp:
						#print "keep"
						sequences[-1][1]+=sanitize(line)
					else:
						#print sequences[-1], "not on keep list"
						continue
				else:
					sequences[-1][1]+=sanitize(line)
			else:
				if len(sequences)>0 and sequences[-1][1]=="":
					del sequences[-1]
				sequences.append( [re.sub('[^0-9]', '',line[1:]),""] )
	print(len(sequences), sequences[0])
	return sequences

#Edit Distance
def levenshteinDistance(s1,s2):
		
    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            #scoremod=matrix[(char1,char2)]-4 #Value between 0 and -8
            if char1 == char2:
                newDistances.append(distances[index1])
                #newDistances.append(distances[index1]-scoremod)
            else:
                newDistances.append(1 + min((distances[index1],
                                             distances[index1+1],
                                             newDistances[-1])))
                #newDistances.append(1 + min((distances[index1],
                #                             distances[index1+1],
                #                             newDistances[-1])) - scoremod)
        distances = newDistances
#   print distances[-1]
    return distances[-1]
#END

#bounded edit distance    
def levenshteinDistanceThresh(s1,s2,k):
#    print "LD:", k
#    print s1
#    print s2

    from Bio.SubsMat import MatrixInfo as matlist
    matrix = matlist.blosum62

    if len(s1) > len(s2):
        s1,s2 = s2,s1
    distances = range(len(s1) + 1)
    for index2,char2 in enumerate(s2):
        newDistances = [index2+1]
        for index1,char1 in enumerate(s1):
            if abs(index1-index2) > k:
                newDistances.append(sys.maxsize)
            elif char1 == char2:
                newDistances.append(distances[index1])
            else:
                newDistances.append(1 + min((distances[index1],distances[index1+1],newDistances[-1])))	
        distances = newDistances

        if min(distances)>k:
            return sys.maxsize		
    return distances[-1]

#DiWANN graph construction
def makeGraphEfficiently(species, repeatmatrix, filename, sparse=False):
	
	G = Graph(directed=True)
	G.add_vertices(len(species))
	
	#summary (G)
	
	names=[]
	for r in labelss:
		names.append('; '.join(r))
	#G.vs["label"] = names
	
###
	for row in	 range(len(species)):
		#populate graph
		#print len(repeatmatrix[row])
		if sparse==True:
			import scipy
			m=min(list(scipy.sparse.find(repeatmatrix.tocsc()[:,row])[2]) + list(scipy.sparse.find(repeatmatrix.tocsc()[row,:])[2]))
			edgeTo = getSparseMinED(repeatmatrix, row, m)
			#return
		else:
			m=min(l for l in repeatmatrix[row] if l > 0)
			edgeTo = getMinED(repeatmatrix, row, m)
		#print edgeTo
		for e in edgeTo:
			G.add_edge(e[0], e[1], weight=m)
		
	#summary(G)
	G.simplify(combine_edges=min)
	summary(G)
	G.write(filename+"DWRNgraphEff.gml","gml")
	#summary(G)
	return G

#threshold based network construction
def makeGraph(species, RepeatMatrix, filename):

	G_new = Graph(directed=False)
	G_new.add_vertices(len(species))
	
	G = Graph()
	G.add_vertices(len(species))
	
	names=[]
	for r in labelss:
		names.append('; '.join(r))

	#G_new.vs["label"] = names
	#G.vs["label"] = names
	
	its=len(species)*len(species)
	it=0
	for i in range(len(species)):
		# print RepeatMatrix[i]
		minED = min(l for l in RepeatMatrix[i] if l > 0)
		
		for j in range(len(species)):
			#print "adding"
			if (RepeatMatrix[i][j] == minED):
				G_new.add_edge(i,j,weight=RepeatMatrix[i][j])
				#weights.append(RepeatMatrix[i][j])
			if i>j and RepeatMatrix[i][j]<2000000:
				G.add_edge(i,j,weight=RepeatMatrix[i][j])
			#print "added"
			it+=1
		print(round((it*1.0)/its,3))
	summary(G_new)
	summary(G)
	
	
	G_new.write(filename+"DWN-base.gml","gml")
	G.write(filename+"Thresh.gml","gml")
 
	return G

#finds lowest in a vector
def getMinED(repeatmatrix, row, minVal):
#	minVal=sys.maxint
	
	#print repeatmatrix
	ret=[]
	
	for col in range(len(repeatmatrix[row])):
#		#print "a"
#		if repeatmatrix[row][col] < minVal:
#			minVal = repeatmatrix[row][col]
#			ret = [(row,col)]
#		elif repeatmatrix[row][col] == minVal:
		if repeatmatrix[row][col] == minVal:
			ret.append((row,col))
		if repeatmatrix[row][col] < minVal and repeatmatrix[row][col] < 0:
			print("Warning, min passed was wrong")
			ret = [(row,col)]
			minVal=repeatmatrix[row][col] 
	
	#print minVal
	#print ret
	return ret
	
def getSparseMinED(repeatmatrix, row, minVal):
#	minVal=sys.maxint
	
	#print repeatmatrix
	ret=[]
	
	for col in range(len(repeatmatrix[row])):
#		#print "a"
#		if repeatmatrix[row][col] < minVal:
#			minVal = repeatmatrix[row][col]
#			ret = [(row,col)]
#		elif repeatmatrix[row][col] == minVal:
		if repeatmatrix[row,col] == minVal:
			ret.append((row,col))
		if repeatmatrix[row,col] < minVal and repeatmatrix[row,col] < 0:
			print("Warning, min passed was wrong")
			ret = [(row,col)]
			minVal=repeatmatrix[row,col] 
	
	#print minVal
	#print ret
	return ret

#creats DiWANN similarity matrix
def createMinMatrix(species, useBoundedLD=False):
    #print("started creatMinMatrix")
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2
    matrix = matlist.blosum62
	#print(matrix)
    skippedCells=0	
    repeatmatrix=[]
    its=(len(species)*len(species))/2
    it=len(species)
    #print("its=",its)
    #print("it=",it)
    for row in range(len(species)):
        #Add new row to matrix
        repeatmatrix.append([])
        repeatmatrix[row].extend([sys.maxsize]*len(species))
    #print(repeatmatrix)
    for row in range(len(species)):
        #print "row:", row
        if row>0:
            print("Percent complete:", round((it*1.0)/its,3))
            #calculate ranges
            maxED=[]
            minED=[]
            if row>1:
                rowMin=min(repeatmatrix[row][0:row-1])
                #print repeatmatrix[row][0:row-1]
                #print rowMin
            else:
                rowMin=sys.maxsize
			
            for col in range(row+1,len(species)):
                minED.append(abs(repeatmatrix[0][col]-repeatmatrix[0][row]))
                maxED.append(repeatmatrix[0][col]+repeatmatrix[0][row])
			#then only compare possible mins
			#print row, len(minED)
			
            if len(maxED)>0:
                lowestMax = min(maxED)
            else:
                lowestMax = sys.maxsize
                #print lowestMax
			
            for col in range(row+1,len(species)):
                it+=1
                colMin = min(repeatmatrix[col])
                #print colMin
                #if col - (row+1) < 0:
                #print "ERROR"
                if (minED[col-(row+1)] > lowestMax or minED[col-(row+1)] > rowMin) and (minED[col-(row+1)] > colMin): 
                    repeatmatrix[row][col]=sys.maxsize
                    repeatmatrix[col][row]=sys.maxsize
					#print "skipping ", row, col
                    skippedCells+=1
                else:
                    if useBoundedLD == True:
                        repeatmatrix[row][col]=levenshteinDistanceThresh(species[row][0],species[col][0], max(colMin,rowMin))
					 #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                        if repeatmatrix[row][col] > max(colMin,rowMin):
                            repeatmatrix[row][col]=sys.maxsize
                        #print max(colMin,rowMin)
                    else:
                        repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
                        #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                    repeatmatrix[col][row]=repeatmatrix[row][col]
                    #if repeatmatrix[row][col] == 0:
                        #print "WARNING!!!"
                    if repeatmatrix[row][col] < rowMin:
                        rowMin=repeatmatrix[row][col]

        else:
            for col in range(len(species)):
                repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
                #repeatmatrix[row][col]=pairwise2.align.globaldx(species.repeats[row].sequence.upper(), species.repeats[col].sequence.upper(), matrix, score_only=True, one_alignment_only=True)
                #print repeatmatrix[col][row]
                repeatmatrix[col][row]=repeatmatrix[row][col]
            repeatmatrix[0][0] = sys.maxsize

    #print("Skipped computing", skippedCells, "of", ((len(species)*len(species))-len(species))/2, "cells")
    return repeatmatrix

#Creats similarity matric for threshold-based methods
def createRepeatMatrix(species, thresh=2):
    from Bio.SubsMat import MatrixInfo as matlist
    from Bio import pairwise2
    matrix = matlist.blosum62

    repeatmatrix=[]
    for row in range(len(species)):
        #Add new row to matrix
        repeatmatrix.append([ listtostring(labelss[row],";") ])
        repeatmatrix[row].extend([0]*len(species))
    its=len(species)*len(species)/2
    i=0
    for row in range(len(species)):
        for col in range(row,len(species)):
            #repeatmatrix[row][col]=levenshteinDistance(species[row][0],species[col][0])
            repeatmatrix[row][col]=levenshteinDistanceThresh(species[row][0],species[col][0],thresh)
            #repeatmatrix[row][col]=pairwise2.align.globaldx(species[row][0].upper(), species[col][0].upper(), matrix, score_only=True, one_alignment_only=True)
            repeatmatrix[col][row]=repeatmatrix[row][col]
            
            i+=1
        print("Percent complete:", round((i*1.0)/its,3))

    return repeatmatrix


def MatrixtoCSV(mat,file):
	with open(file, 'wb') as csv:
		csv.write(u'\ufeff'.encode('utf-8'))
		for row in mat:
			for col in row:
				csv.write(str(col)+",")
			csv.write("\n")

def create_dendrogram(df_matrix):
    hc=AgglomerativeClustering(n_clusters=7,affinity='precomputed',linkage='complete')
    y_hc=hc.fit_predict(df_matrix)
    linkage_matrix=sch.linkage(df_matrix,'complete')
    #visualize clusters
    plt.figure(figsize=(35, 15))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    sch.dendrogram(linkage_matrix,leaf_rotation=90.,leaf_font_size=8.,labels=labelss)
    plt.show()
    plt.savefig("Dendrogram01.jpg")
    return y_hc

def makeNetwork(filename,thresh,bthresh,type="", species=None,BL=""):
    graph=""
    if type=="PB":
        print("Prune + bound*****************")
        print(filename)
        start = time.time()
        C = createMinMatrix(species, False)
        C = np.asarray(C)
        end = time.time()
       # y_hc=create_dendrogram(C)
        print("Pruning and Bounding:", (end-start))
        #MatrixtoCSV(C,"DWRN"+filename+".csv")
        graph = makeGraphEfficiently(species, C, filename+"minNet")
        summary(graph)
    elif type=="BFT":
        print("Brute force*****************")
        #thr=len(species[0][0])/3
        thr=thresh
        print(filename)#, "up to", thr
        start = time.time()
        C = createRepeatMatrix(species, thr)
        #C = createRepeatMatrix(species)
        end = time.time()
        print("Brute Force (my implementation):", (end-start))
        #MatrixtoCSV(C,"BF"+filename+".csv")
        graph = makeGraph(species, C, filename+"fullNet")
        graph.simplify(combine_edges=min)
        summary(graph)
    elif type=="Blast":
        print("Approximate distances*****************")
        print(filename, BL)
        import csv
        nodes={}
        graph=Graph()
        sum=0
        with open(BL, "r") as file:
            reader=csv.reader(file,delimiter="\t")
            for line in reader:
                if line[0] not in nodes:
                    nodes[line[0]]=len(nodes)
                    #print(line[0])
                    #print(nodes)
                    graph.add_vertex(name=line[0])
                if line[1] not in nodes:
                    nodes[line[1]]=len(nodes)
                    #print(line[1])
                    graph.add_vertex(name=line[1])
                #if float(line[10])<1.0e-10:
                #if float(line[2])>92.0:
                if float(line[11])>bthresh:
					
                    graph.add_edge(nodes[line[0]],nodes[line[1]],weight=float(line[11]))
                    sum+=1
            #summary(graph)
            graph=graph.simplify(combine_edges=min)
            summary(graph)
            graph.write(filename+"BlastGraph.gml","gml")
    return graph


#centrality
#marginale
def create_centrality_plot_marginale(network):
    am_centrality=network.pagerank()
    central=sorted(range(len(am_centrality)), key=lambda x: am_centrality[x])[-10:]
    central_labels=[]
    for i in central:
        central_labels.append(labelss[i][0])
    amarg_centrality=Graph()
    amarg_centrality=network
    nodes=amarg_centrality.vs.select(central)
    nodes["color"]="yellow"
    #nodes["size"]=30
    amarg_centrality.vs['size']=[i * 1000 for i in am_centrality]
    return amarg_centrality
    
#threshold
def create_centrality_plot(network):
    t_centrality=network.pagerank()
    central_t=sorted(range(len(t_centrality)), key=lambda x: t_centrality[x])[-10:]
    central_labeless_t=[]
    for i in central_t:
        central_labeless_t.append(labelss[i][0])
    th_centrality=Graph()
    th_centrality=network
    nodes=th_centrality.vs.select(central_t)
    nodes["color"]="yellow"
    #th_centrality.vs['size']=[i * 1000 for i in t_centrality]
    return th_centrality
    
style = {}
style["edge_curved"] = False
style["margin"] = 100
style["edge_arrow_size"]=0.5
style["vertex.label.cex"]=0.75
style["vertex.label.family"]="Helvetica"
style["vetrex.label.color"]='black'
style["edge_width"]=0.5
style["edge_color"]="blue"
style["vertex_size"]=6
style["vertex.label.font"]=2

dstyle={}
#dstyle["vertex_size"]=d_centrality*10000
dstyle["edge_curved"] = False
dstyle["margin"] = 100
dstyle["edge_arrow_size"]=0.5
dstyle["vertex.label.cex"]=0.75
dstyle["vertex.label.family"]="Helvetica"
dstyle["vetrex.label.color"]='black'
dstyle["edge_width"]=0.5
dstyle["edge_color"]="blue"
dstyle["vertex.label.font"]=2
dstyle["vertex_size"]=6
#dstyle["vertex_color"]="blue"    

#Program starts here
for choice in range(1,5):
    
    if(choice == 1):
        with open('AM.csv',encoding='UTF-16') as fh:
            test = pd.read_csv(fh,header=0,delimiter="\t")
        fh.close()
        x_labels = test.iloc[:,[0]]
        labelss=x_labels.values.tolist()
        sequence=test.iloc[:,[1]]
        species=[]
        species=sequence.values.tolist()
        rfile="amarginale"
        
        #DiWANN marginale
        diwann_am=Graph()				
        diwann_am=makeNetwork(rfile,thresh="",bthresh="",type="PB",species=species)
        plot(diwann_am,'diwann_plot_marginale.png',**style)
        
        #marginale threshold
        threshold_exact_am=Graph()
        threshold_exact_am=makeNetwork('AM.csv', thresh=2,bthresh="" ,type="BFT", species=species)
        plot(threshold_exact_am,'threshold_exact_plot_marginale.png',**style)
        
        #blast marginale
        blast_graph_am=Graph()
        blast_graph_am=makeNetwork('AM.csv',thresh="",bthresh=44.8,type="Blast",species=species,BL="amarHit.csv")
        plot(blast_graph_am,'blast_graph_plot_marginale.png',**style)

        #degree distribution
        #Diwann
        xs, ys = zip(*[(left, count) for left, _, count in 
        diwann_am.degree_distribution().bins()])
        pylab.bar(xs, ys)
        pylab.xlabel('Degree')
        pylab.ylabel('Frequency')
        pylab.title('DiWANN Degree Distribution(A. Marginale)')
        pylab.show()

        #exact
        xs, ys = zip(*[(left, count) for left, _, count in 
        threshold_exact_am.degree_distribution().bins()])
        pylab.bar(xs, ys)
        pylab.xlabel('Degree')
        pylab.ylabel('Frequency')
        pylab.title('Exact Threshold Degree Distribution(A. Marginale)')
        pylab.show()

        #blast
        xs, ys = zip(*[(left, count) for left, _, count in 
        blast_graph_am.degree_distribution().bins()])
        pylab.bar(xs, ys)
        pylab.xlabel('Degree')
        pylab.ylabel('Frequency')
        pylab.title('Inexact Threshold Degree Distribution (A. Marginale)')
        pylab.show()
        
        #am centrality    
        diwann_am_centrality=Graph()
        diwann_am_centrality=create_centrality_plot_marginale(diwann_am)    
        plot(diwann_am_centrality,'am_diwann_centrality.png',**dstyle)
        
        exact_am_centrality=Graph()
        exact_am_centrality=create_centrality_plot_marginale(threshold_exact_am)    
        plot(exact_am_centrality,'am_exact_centrality.png',**dstyle)
        
        inexact_am_centrality=Graph()
        inexact_am_centrality=create_centrality_plot_marginale(blast_graph_am)    
        plot(inexact_am_centrality,'am_inexact_centrality.png',**dstyle)
        
    if(choice == 2):
        with open('rice_protein_seq.txt',encoding='UTF-16') as fh:
            test_rice = pd.read_csv(fh,header=0,delimiter="\t")
        fh.close()
        x_labels = test_rice.iloc[:,[0]]
        labelss=x_labels.values.tolist()
        sequence=test_rice.iloc[:,[1]]
        species=[]
        species=sequence.values.tolist()
        rfile="rice"
        
        #rice diwann
        diwann_rice=Graph()				
        diwann_rice=makeNetwork(rfile, thresh="", bthresh="", type="PB", species=species)
        plot(diwann_rice,'diwann_plot_rice.png',**style)

        #rice threshold
        threshold_exact_rice=Graph()
        threshold_exact_rice=makeNetwork('rice_protein_seq.txt', thresh=81.2, bthresh="", type="BFT", species=species)
        plot(threshold_exact_rice,'threshold_exact_plot_rice.png',**style)
        
        #rice blast
        blast_graph_rice=Graph()
        blast_graph_rice=makeNetwork('rice_protein_seq.txt', thresh="", bthresh=19.8,type="Blast",species=species,BL="riceSimilarity")
        plot(blast_graph_rice,'blast_graph_plot_rice.png',**style)
        
        #rice centrality
        diwann_rice_centrality=Graph()
        diwann_rice_centrality=create_centrality_plot(diwann_rice)    
        plot(diwann_rice_centrality,'rice_diwann_centrality.png',**dstyle)
        
    if(choice == 3):
        with open('wheat_protein_seq.txt',encoding='UTF-16') as fh:
            test_wheat = pd.read_csv(fh,header=0,delimiter="\t")
        fh.close()
        x_labels = test_wheat.iloc[:,[0]]
        labelss=x_labels.values.tolist()
        sequence=test_wheat.iloc[:,[1]]
        species=[]
        species=sequence.values.tolist()
        rfile="wheat"
        
        #wheat diwann
        diwann_wheat=Graph()				
        diwann_wheat=makeNetwork(rfile, thresh="", bthresh="", type="PB", species=species)
        plot(diwann_wheat,'diwann_plot_wheat.png',**style)

        #wheat threshold
        threshold_exact_whaet=Graph()
        threshold_exact_wheat=makeNetwork('wheat_protein_seq.txt', thresh=75.5, bthresh="",type="BFT", species=species)
        plot(threshold_exact_wheat,'threshold_exact_plot_wheat.png',**style)
        
        #wheat blast
        blast_graph_wheat=Graph()
        blast_graph_wheat=makeNetwork('wheat_protein_seq.txt', thresh="", bthresh=19.8,type="Blast",species=species, BL="wheatSimilarity")
        plot(blast_graph_wheat,'blast_graph_plot_wheat.png',**style)
        
        #wheat centrality
        diwann_wheat_centrality=Graph()
        diwann_wheat_centrality=create_centrality_plot(diwann_wheat)    
        plot(diwann_wheat_centrality,'wheat_diwann_centrality.png',**dstyle)
        
    if(choice==4):
        with open('rice_indica_protein_seq.txt',encoding='UTF-16') as fh:
            test_indica = pd.read_csv(fh,header=0,delimiter="\t")
        fh.close()
        x_labels = test_indica.iloc[:,[0]]
        labelss=x_labels.values.tolist()
        sequence=test_indica.iloc[:,[1]]
        species=[]
        species=sequence.values.tolist()
        rfile="indica"
        
        #DiWANN indica
        diwann_indica=Graph()				
        diwann_indica=makeNetwork(rfile,thresh="",bthresh="",type="PB",species=species)
        plot(diwann_indica,'diwann_plot_indica.png',**style)
        
        diwann_indica_centrality=Graph()
        diwann_indica_centrality=create_centrality_plot(diwann_indica)    
        plot(diwann_indica_centrality,'indica_diwann_centrality.png',**dstyle)
        

