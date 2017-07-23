#/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx

#==============================================
def cria_rede_barabasi_albert(N=9,m=1):
    """
    Retorna a matriz de adjacência e a rede tipo 
    barabasi-albert de N vertices e parametro m.    
    """
    g = nx.barabasi_albert_graph(N,m)
    # A matriz de adjacencia é calculada por meio da função
    # nx.adjacency_matrix(g). O valor retornado é uma
    # matriz esparsa (forma eficiente de armazenar numeros).
    # Como não estamos preocupados com memória, vamos deixar
    # a matriz numa forma densa convencional.
    return np.asarray((nx.adjacency_matrix(g)).todense(),dtype=bool);
#==============================================
def cria_rede_quadrada(N=9,dim=2):
    """
    Retorna a matriz de adjacência e a rede regular quadrada 
    de dimensão dim, com N vertices, e condições periodicas.
    """                    
    A=np.zeros((N,N),dtype=bool);
    l=np.int( N**(1.0/dim) );
    L=np.array([l**k for k in np.arange(dim+1)])
    k=0;
    while k < N:
        c=np.array([np.mod(k,L[n+1])//L[n] for n in np.arange(dim)])
        j=0;
        while j< dim:
            v   =c.copy();
            v[j]=np.mod(v[j]+1,l);
            k1  =np.sum(L[:-1]*v);
            A[k1,k]=True;
            j+=1
        i=0;
        while i< dim:
            b   =c.copy();
            b[i]=np.mod(b[i]-1+l,l);
            k2  =np.sum(L[:-1]*b);
            A[k2,k]=True;
            i+=1
        k+=1
    return A                    
#==============================================
def densidade_grau(A):
    """
    Retorna uma lista com o grau de cada vértice 
    divido pelo total de vértices.
    """
    return np.sum(A,0)*1.0/len(A)
#==============================================
def inicializa_rede(N):
    """
    Retorna uma lista com valores +1 ou -1 distribuidos
    uniformemente.
    """
    # configuracao de spin
    # se rand > 0.5, produz spin =+1 == pessoa doente
    # se rand < 0.5, produz spin =-1 == pessoa saudavel
    return  np.asarray(2*(np.random.rand(N) > 0.5) -1,dtype=np.int8); 
#==============================================
def metropolis(beta,delta_energia):
    """
    Retorna uma variavel booleana para determinar se 
    deve ou não atualizar uma configuracao dado o valor de
    delta energia  e beta.
    Se delta energia < 0 : aceita
    Se delta energia > 0 : aceita com probabilidade exp(-beta*delta)
    """
    return any((delta_energia < 0, np.random.rand()<np.exp(-beta*delta_energia)))
#==============================================
def atualiza(conf,beta,A,kappa,params):
    """
    Atualiza a configuração conf via metropolis
    
    conf  : array_like. Cada entrada pode ser +1 ou -1
    beta  : float. Inverso da temperatura
    A     : array_like. Matriz de adjacência (NxN)
    kappa : array_like. Graus de cada vértice/N
    params: array_like. [0]=N; [1]=alpha ; [2]=gamma
    """
    N    =params[0];
    alpha=params[1]; # params[1] é o parametro epid (infe)
    gamma=params[2];    # params[2] é o parametro epid (cura)    
    k=0;
    while k < N:
        # energia do k-esimo sitio
        energia =(-(alpha*0.25/N)*np.sum(A[:,k]*conf)-(alpha*kappa[k]-gamma)*0.5)*conf[k] ;
        # variacao da energia caso troque o spin
        delta_energia = -2.0*energia;
        # METROPOLIS
        # se metropolis = True , conf[k]=-conf[k]
        # se metropolis = False, conf[k]=+conf[k]
        conf[k]=(1-2*metropolis(beta,delta_energia))*conf[k]
        k+=1    
    return None
#==============================================
def argumentos(argv):
    # PARAMETROS
    # numero de agentes da comunidade a
    N=int(argv[1])
    
    medidas=100
    alpha  =4*N
    gamma  =16
    beta_min=0.01
    beta_max=1.01
    delta_beta=0.025
    try:
        medidas   =int(argv[2])
        alpha     =float(argv[3])
        gamma     =float(argv[4])
        beta_min  =float(argv[5])
        beta_max  =float(argv[6])
        delta_beta=float(argv[7])
    except:
        pass
    return (N,medidas,alpha,gamma,beta_min,beta_max,delta_beta)
