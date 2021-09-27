"""
Nome: Alessandro Melo de Oliveira
NºUSP: 10788662

Trabalho 1 - Aerodinâmica Computacional SAA0199 
Resolver a equação do calor pelo método implícito (Crank-Nicolson)
"""

import numpy as np # Biblioteca para cálculos de álgebra
import matplotlib.pyplot as plt # Biblioteca para fazer os plots
from scipy.linalg import solve # Solver para resolver o sistema Ax=b

class EqCalor:
    
    def __init__(self):
        
        self.L = 1 # Comprimento da haste [m]
        self.T0 = 300 # Temperatura no começo da haste [K]
        self.TM = 300 # Temperatura média [K]
        self.TL = 1000 # Temperatura no final da haste [K]
        self.alpha = 0.1 # Coeficiente de difusividade térmica [m^2/s] 
        self.nx = 100 # Intervalos em que a barra foi dividida
        self.nt = 10000 # Intervalos em que o tempo foi dividido
        self.dx = self.L/self.nx # Discretização da haste
        self.dt = (self.dx**2)/(2*self.alpha) # Discretização do tempo de simulacao
        self.C = (self.dx**2)/(self.alpha*self.dt) # Parâmetro auxiliar
        
        self.x = np.linspace(0, self.L, self.nx+1) # Vetor discretizado das posições
        self.t = [self.dt*it for it in np.arange(0, self.nt)] # Vetor discretizado dos tempos
        
        self.Ta = np.zeros((self.nt, self.nx+1)) # Matriz com os resultados de temperatura utilizando a solução analítica
        self.T = np.zeros((self.nt, self.nx+1)) # Matriz com os resultados de temperatura utilizando a aproximação numérica
        self.erro_max = np.zeros((self.nt))
        
        self.A = np.zeros((self.nx+1, self.nx+1)) # Matriz A, usada na solução do problema
        self.B = np.zeros((self.nt, self.nx+1)) # Vetor B, utilizado na solução do problema
        
        # Montando a estrutura da matriz A, que é independente do tempo
        self.A[0,0] = 1
        self.A[self.nx, self.nx] = 1
        for i in np.arange(1, self.nx):
            for j in np.arange(0, self.nx+1):
                if i == 0 and j==0:
                    self.A[i,j] = 1
                elif i == self.nx and j == self.nx:
                    self.A[i,j] = 1
                elif i == j:
                    self.A[i,j] = 2*(self.C+1)
                    
                elif j == i-1 or j == i+1:
                    self.A[i,j] = -1
                else:
                    self.A[i,j] = 0
            
            #Temperatura inicial da barra       
            self.Ta[0,i] = self.TM*np.sin(self.x[i]*np.pi/self.L) + self.x[i]*(self.TL - self.T0)/self.L + self.T0 # temperatura na barra no instante t = 0
            self.T[0, i] = self.TM*np.sin(self.x[i]*np.pi/self.L) + self.x[i]*(self.TL - self.T0)/self.L + self.T0 # temperatura na barra no instante t = 0
        
        #Condições de contorno
        self.T[:, 0] = self.T0
        self.T[:, -1] = self.TL
    
    #Essa função efetivamente implementa o método numérico
    def solucao(self):
        
        self.B[:,0] = self.T0
        self.B[:,-1] = self.TL
        
        for it in np.arange(0, self.nt-1):
            
            self.Ta[it,0] = self.T0
            self.Ta[it,self.nx] = self.TL
            
            for ix in np.arange(1, self.nx):
                
                self.B[it,ix] = self.T[it,ix-1] + 2*(self.C-1)*self.T[it,ix]+self.T[it,ix+1]
                self.Ta[it+1,ix]=self.TM*np.exp(-self.alpha*(np.pi/self.L)**2*self.t[it+1])*np.sin(np.pi/self.L*self.x[ix])+(self.TL-self.T0)/self.L*self.x[ix]+self.T0

            self.T[it+1] = solve(self.A, self.B[it])
            
        self.error = abs((self.T - self.Ta)/(self.Ta)) # Calcula o erro relativo
        self.error_max = np.amax(self.error, axis = 1) # Calculando o erro máximo para cada tempo
        
    
    def plot(self, t_range, save = False):
        fig_plot = plt.figure(figsize=(10,5))
        
        for t in t_range:
            plt.plot(self.x, self.T[t,:], label = "t: {}".format(self.t[t]), lw = 2)
        plt.xlabel("Comprimento da barra [m]", fontsize = 12)
        plt.legend()
        plt.ylabel("Temperatura [K]", fontsize = 12)
        plt.grid()
        
        if save == True:
            plt.savefig("plot_graphs_distribution_t_range.pdf")
        plt.show()
        
    def plot_erro_relativo(self, t_range, save = False):
        fig_plot_erro = plt.figure(figsize=(10,5))
        
        for t in t_range:
            plt.plot(self.x, self.error[t,:], label = "t: {}".format(self.t[t]), lw = 2)
        plt.xlabel("Comprimento da barra [m]", fontsize = 12)
        plt.legend()
        plt.ylabel("Erro relativo", fontsize = 12)
        plt.grid()
        
        if save == True:
            plt.savefig("plot_graphs_error_t_range.pdf")
        plt.show()
        
        
    def plot_erro_maximo(self, save = False):
        fig_plot_erro_maximo = plt.figure(figsize=(10,5))
        plt.plot(self.t, self.error_max, lw = 2)
        plt.xlabel("Intervalo de tempo [s]", fontsize = 12)
        plt.legend()
        plt.ylabel("Erro máximo", fontsize = 12)
        plt.grid()
        
        if save == True:
            plt.savefig("plot_graphs_erro_maximo.pdf")
        plt.show()
