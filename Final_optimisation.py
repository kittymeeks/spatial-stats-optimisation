# -*- coding: utf-8 -*-


import numpy as np
import copy



'''MATRIX/GRAPH/SET MANIPULATION UTILITIES'''

def get_nbs(index,matrix):
    '''returns a list of the neighbours of vertex index from the adjacency matrix'''
    nbs = []
    for i in range(0,len(matrix[index])):
        if matrix[index][i]==1:
            nbs.append(i)
    return nbs


def powerset(seq):
    ''' generator returning the powerset of a set '''
    if len(seq) <= 1:
        yield seq
        yield []
    else:
        for item in powerset(seq[1:]):
            yield [seq[0]]+item
            yield item
            

def get_edges(matrix):
    '''computes an edge list from the given adjacency matrix'''
    edges = []
    for i in range(0,len(matrix)):
        for j in range(i+1,len(matrix[i])):
            if matrix[i][j] == 1:
                edges.append([i,j])
    return edges



def degree(vertex,matrix):
    deg = 0
    for i in range(0,len(matrix[vertex])):
        deg = deg + matrix[vertex][i]
    return deg

def remove_edges(matrix,edges):
    new_matrix = copy.deepcopy(matrix)
    for e in edges:
        new_matrix[e[0]][e[1]]=0
        new_matrix[e[1]][e[0]]=0
    return new_matrix





'''OPTIMISATION FUNCTIONS'''




def cond_sorted_nbs(i,nbs,options,data,avs,n):
    fixed_nbs = [nb for nb in nbs if nb not in options]
    nbr_options = [fixed_nbs + poss_nbs for poss_nbs in powerset(options) if len(fixed_nbs + poss_nbs) > 0]
    nbr_list = [(nbs,new_vertex_val(i,nbs,data,avs,n)) for nbs in nbr_options]
    return sorted(nbr_list,key=lambda x: -x[1])

def cond_sorted_nbs_mult(i,nbs,options,data,avs,n):
    fixed_nbs = [nb for nb in nbs if nb not in options]
    nbr_options = [fixed_nbs + poss_nbs for poss_nbs in powerset(options) if len(fixed_nbs + poss_nbs) > 0]
    nbr_list = [(nbs,vertex_val_multiple(i,nbs,data,avs,n)) for nbs in nbr_options]
    return sorted(nbr_list,key=lambda x: -x[1])



def find_best_cond(nbr_list,with_nbs,without_nbs):
    i = 0
    found = False
    while (i < len(nbr_list)):
        cand_nbs = nbr_list[i][0]
        found = True
        for w in with_nbs:
            if w not in cand_nbs:
                found = False
        for w in without_nbs:
            if w in cand_nbs:
                found = False
        if found:
            return nbr_list[i]
        i = i + 1


def greedy_opt(matrix,values,avsum):
    n = len(matrix)
    newmatrix = matrix
    for vx in range(0, len(matrix)):
    # for every vertex in the matrix
        if len(get_nbs(vx,newmatrix)) > 1:
        # if it has degree greater than one (after updates so far)
            nb_options = [nb for nb in get_nbs(vx,newmatrix) if nb > vx and len(get_nbs(nb,newmatrix)) > 1]
            # nb_options is the set of neighbours from which we can potentially delete - those that have degree greater than one,
            # and which we haven't yet considered
            sorted_nbd_list = cond_sorted_nbs(vx,get_nbs(vx,newmatrix),nb_options,values,avsum,n)
            # sorted_nbd_list is the sorted list of all neighbour options for vx
            nbr_sorted_nbd_lists = []
            for u in nb_options:
                u_options = [w for w in get_nbs(u,newmatrix) if w >= vx and len(get_nbs(w,newmatrix)) > 1]
                # u_options is the set of neighbours for u that we could still change
                nbr_sorted_nbd_lists += [cond_sorted_nbs(u,get_nbs(u,newmatrix),u_options,values,avsum,n)]
                # add the sorted list of neighbour options for u to the nbr_sorted_nbd_lists
            remove = []
            for i in range(0,len(nb_options)):
                nb = nb_options[i]
                vx_gain = find_best_cond(sorted_nbd_list,[nb],[])[1] - find_best_cond(sorted_nbd_list,[],[nb])[1]
                nb_gain = find_best_cond(nbr_sorted_nbd_lists[i],[vx],[])[1] - find_best_cond(nbr_sorted_nbd_lists[i],[],[vx])[1]
                combined_gain = vx_gain + nb_gain
                if (combined_gain < 0):
                    remove.append([[vx,nb],combined_gain])
            if len(remove) < len(get_nbs(vx,newmatrix)):
                newmatrix = remove_edges(newmatrix,[e[0] for e in remove])
            else:
                sorted_remove = sorted(remove,key=lambda x: x[1])
                newmatrix = remove_edges(newmatrix,[e[0] for e in sorted_remove[1:]])
    return newmatrix


            
def greedy_opt_multiple(matrix,data,avsums):
    n = len(matrix)
    newmatrix = matrix
    for vx in range(0, len(matrix)):
        if len(get_nbs(vx,newmatrix)) > 1:
            nb_options = [nb for nb in get_nbs(vx,newmatrix) if nb > vx and len(get_nbs(nb,newmatrix)) > 1]
            sorted_nbd_list = cond_sorted_nbs_mult(vx,get_nbs(vx,newmatrix),nb_options,data,avsums,n)
            nbr_sorted_nbd_lists = []
            for u in nb_options:
                u_options = [w for w in get_nbs(u,newmatrix) if w >= vx and len(get_nbs(w,newmatrix)) > 1]
                nbr_sorted_nbd_lists += [cond_sorted_nbs_mult(u,get_nbs(u,newmatrix),u_options,data,avsums,n)]
            remove = []
            for i in range(0,len(nb_options)):
                nb = nb_options[i]
                vx_gain = find_best_cond(sorted_nbd_list,[nb],[])[1] - find_best_cond(sorted_nbd_list,[],[nb])[1]
                nb_gain = find_best_cond(nbr_sorted_nbd_lists[i],[vx],[])[1] - find_best_cond(nbr_sorted_nbd_lists[i],[],[vx])[1]
                combined_gain = vx_gain + nb_gain
                if (combined_gain < 0):
                    remove.append([[vx,nb],combined_gain])
            if len(remove) < len(get_nbs(vx,newmatrix)):
                newmatrix = remove_edges(newmatrix,[e[0] for e in remove])
            else:
                sorted_remove = sorted(remove,key=lambda x: x[1])
                newmatrix = remove_edges(newmatrix,[e[0] for e in sorted_remove[1:]])
    return newmatrix           
        




def default_avsum(matrix,values):
    return sum([disc(i,matrix,values) for i in range(0,len(matrix))])

            

 
def new_vertex_val(index,nbs,values,avsum,n):
    vx_av = len(nbs)*(values[index] - sum([values[u] for u in nbs])/len(nbs))**2
    ub = np.log(len(nbs))/2 - (n/2)*(np.log(1 + vx_av/max(0.001,avsum - vx_av)))
    return ub


def vertex_val_multiple(index,nbs,data,avsums,n):
    contrib = 0
    deg = len(nbs)
    q = len(data)
    for i in range(0,q):
        vx_av = len(nbs)*(data[i][index] - sum([data[i][u] for u in nbs])/deg)**2
        val = q*np.log(deg)/2 - (n*q/2)*(np.log(1 + vx_av/max(0.001,avsums[i]-vx_av)))
        contrib += val
    return contrib

    


def disc(index,matrix,values):
    deg = degree(index,matrix)
    result = deg*(values[index] - sum([values[j] for j in get_nbs(index,matrix)])/deg)**2
    return result

def matrix_score(matrix,values):
    n = len(matrix)
    pluspart = sum([np.log(degree(i,matrix)) for i in range(0,n)])/2
    penalty = (n/2)*np.log(sum([disc(index,matrix,values) for index in range(0,n)])/n)
    return pluspart - penalty


def matrix_score_mult(matrix,data):
    n = len(matrix)
    q = len(data)
    pluspart = (q/2)*sum([np.log(degree(i,matrix)) for i in range(0,n)])
    penalty = (n*q/2)*np.log(sum([sum([disc(index,matrix,values) for index in range(0,n)]) for values in data])/(n*q))
    return pluspart- penalty
    

def iterative_opt(matrix,values):
    newscore = matrix_score(matrix,values)
    newmatrix = matrix
    oldscore = -100000
    while oldscore < newscore:
        oldmatrix = newmatrix
        oldscore = newscore
        def_avsum = default_avsum(oldmatrix,values)
        newmatrix = greedy_opt(oldmatrix,values,def_avsum)
        newscore = matrix_score(newmatrix,values)
    return oldmatrix


def iterative_opt_multiple(matrix,data):
    newscore = matrix_score_mult(matrix,data)
    newmatrix = matrix
    oldscore = -100000
    while oldscore < newscore:
        oldmatrix = newmatrix
        oldscore = newscore
        def_avsums = [default_avsum(oldmatrix,dat) for dat in data]
        newmatrix = greedy_opt_multiple(oldmatrix,data,def_avsums)
        newscore = matrix_score_mult(newmatrix,data)
    return oldmatrix   


