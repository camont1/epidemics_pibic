#/usr/bin/python3
# -*- coding: utf-8 -*-


import time
from montecarlo import *

#==============================================    
# EXECUCAO COMO SCRIPT
#==============================================    
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import sys
    if len(sys.argv) < 2:
        print(u"""
        Número incorreto de argumentos
        Exemplo:
        
        python3 montecarlo N medidas alpha gamma beta0 beta1 delta
        
        N : número de vertices
        medidas : número de medidas realizadas para cada beta
        alpha   : taxa de infecção (não é alpha/N)
        gamma   : taxa de cura
        beta0   : valor inicial de beta (inv. temperatura)
        beta1   : valor final de beta
        delta   : incremento de beta

        Valores padrão: medidas=100, alpha=4*N, gamma=0,
                        beta0=0.01, beta1=2.01, delta=0.1        
        """)
        quit();

    (N,medidas,alpha,gamma,beta_min,beta_max,delta_beta)=argumentos(sys.argv)
    params=(N,alpha,gamma)
    num_corr=np.sqrt(N)//2;
    #------------------------------------
    betas = np.arange(beta_min,beta_max,delta_beta)
    mag = np.zeros( (len(betas),2) )
    chi = np.zeros( (len(betas),2) )
    m1=np.zeros(medidas);
    m2=np.zeros(medidas);
    #------------------------------------
    # mude aqui a rede aqui
    A    =cria_rede_barabasi_albert(N,3) 
    kappa=densidade_grau(A);
    conf =inicializa_rede(N);

    print(u"""
    ------------------------------------------------------------
    PARAMETROS:
    
    N=%d, medidas=%d, alpha=%f, gamma=%f, betas=[%f:%f]
    ------------------------------------------------------------
    """ % (N,medidas,alpha,gamma,beta_min,beta_max))

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

    f=open('mag_N%d.dat'%(N),'w')
    for k in range(len(betas)):
        f.write('%f %f %f \n' %(betas[k],mag[k,0],mag[k,1]));
    f.close()
    f=open('chi_N%d.dat'%(N),'w')
    for k in range(len(betas)):
        f.write('%f %f %f \n' %(betas[k],chi[k,0],chi[k,1]));
    f.close()

    
    plt.figure(1)
    plt.plot(betas,mag[:,0],'.-');
    # 1 sigmas (68% de conf)
    plt.fill_between(betas,mag[:,0]-mag[:,1],mag[:,0]+mag[:,1],facecolor='blue',alpha=0.25)
    plt.xlabel('beta')
    plt.ylabel('magnetizacao')
    plt.figure(2)
    plt.plot(betas,chi[:,0],'.-',color='blue');
    # 1 sigma (68% de conf)
    plt.fill_between(betas,chi[:,0]-chi[:,1],chi[:,0]+chi[:,1],facecolor='blue',alpha=0.25)
    plt.xlabel('beta')
    plt.ylabel('susceptibilidade')
    
    plt.show()
