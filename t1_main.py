"""
Nome: Alessandro Melo de Oliveira
NºUSP: 10788662

Trabalho 1 - Aerodinâmica Computacional SAA0199 
Arquivo principal
"""

from t1_10788662 import EqCalor

%%
eq_calor = EqCalor() #Instancia o problema, juntamente com suas variáveis
eq_calor.solucao() # Implementa as soluções do problema (numérica e analítica)
eq_calor.plot(t_range = [0, 250, 500, 1000, 1500, 2000, 4000, 8000, 8500], save = True) # Plota a variação da temperatura ao longo do tempo
eq_calor.plot_erro_relativo(t_range = [0, 250, 500, 1000, 1500, 2000, 4000, 8000, 8500], save = True) # Plota a variação da distribuição do erro ao longo do tempo
eq_calor.plot_erro_maximo(save = True) # Plota o erro máximo em cada instante de tempo da simulação
