#/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx
import time
#import matplotlib.pyplot as plt

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
    alpha=params[1];
    gamma=params[2];        
    k=0;
    while k < N:
        # energia do k-esimo sitio
        energia = ((-alpha*0.25/N)*np.sum(A[:,k]*conf)
                   -(alpha*kappa[k]-gamma)*0.5)*conf[k] ;
        # variacao da energia caso troque o spin
        delta_energia = -2.0*energia;
        # METROPOLIS
        # se metropolis = True , conf[k]=-conf[k]
        # se metropolis = False, conf[k]=+conf[k]
        conf[k]=(1-2*metropolis(beta,delta_energia))*conf[k]
        k+=1    
    return None
#==============================================    



#==============================================    
if __name__ == "__main__":

    
    # PARAMETROS
    # numero de agentes da comunidade a
    N    = 100
    # parametros epidemiologicos (recupera ising-2d)
    alpha= N*4; 
    gamma= alpha*(4.0/N);
    params=(N,alpha,gamma)
    # parametros da simulacao de MC
    beta_min=0.01;
    beta_max= 1.01;
    delta_beta=0.05;
    medidas=100;
    num_corr=np.sqrt(N)//2;
    #------------------------------------
    betas = np.arange(beta_min,beta_max,delta_beta)
    mag = np.zeros( (len(betas),2) )
    chi = np.zeros( (len(betas),2) )
    m1=np.zeros(medidas);
    m2=np.zeros(medidas);
    #------------------------------------

    A    =cria_rede_quadrada(N,dim=2)
    kappa=densidade_grau(A);
    conf =inicializa_rede(N);


    i=0
    for beta in betas:
        start = time.time()
        k=0
        while k < medidas:
            #===========================
            # laço para remover autocorrelaçao entre configuracoes
            # por se tratar de uma cadeia de markov, uma atualizacao
            # de metropolis não remove a correlacao entre duas
            # configurações sucessivas. Para evitar carregar essas
            # autocorrelacoes, descarta-se num_corr atualizacoes
            #===========================
            j=0
            while j <num_corr:
                atualiza(conf,beta,A,kappa,params)        
                j+=1
            #==========================
            m1[k]=np.sum(conf)*1.0/N
            k+=1
        m2=m1*m1 ;
        mag[i,0] = np.mean(m1);
        mag[i,1] = np.sqrt(np.var(m1));
        chi[i,0] = np.mean(m2) - mag[i,0]**2 ;
        chi[i,1] = np.sqrt(np.var(m2));
        i+=1

        end = time.time()
        print('beta = %f, mag =%f, elapsed time = %f' % (beta,mag[i-1,0],end-start))


    # plt.figure(1)
    # plt.plot(betas,mag[:,0],'.-');
    # plt.fill_between(betas,mag[:,0]-mag[:,1],mag[:,0]+mag[:,1],facecolor='blue',alpha=0.25)
    # plt.xlabel('beta')
    # plt.ylabel('magnetizacao')
    # # plt.savefig('mag_barabasi_N%d_m%d.eps' % (N,barabasi))
    # plt.figure(2)
    # plt.plot(betas,chi[:,0],'.-',color='blue');
    # plt.fill_between(betas,chi[:,0]-chi[:,1],chi[:,0]+chi[:,1],facecolor='blue',alpha=0.25)
    # plt.xlabel('beta')
    # plt.ylabel('susceptibilidade')
    # # plt.savefig('chi_barabsi_N%d_m%d.eps' % (N,barabasi))

    f=open('mag_N%d.dat'%(N),'w')
    for k in range(len(betas)):
        f.write('%f %f %f \n' %(betas[k],mag[k,0],mag[k,1]));
    f.close()
    f=open('chi_N%d.dat'%(N),'w')
    for k in range(len(betas)):
        f.write('%f %f %f \n' %(betas[k],chi[k,0],chi[k,1]));
    f.close()


    # plt.show()
    
