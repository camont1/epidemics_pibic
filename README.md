# epidemics
Simulação de Monte Carlo para propagação de doenças endêmicas (modelo SIS), no equilíbrio.

## Contexto

A propagação de doenças descritas pelo modelo epidemiológico SIS geralmente utiliza uma descrição de fenômenos fora do equilíbrio. Argumentamos em (colocar referencia) que também é possível descrevê-la como se fosse um processo no regime de equilíbrio.


## Instalação

- para ter acesso aos métodos dentro do interpretador, digite
```python
import montecarlo
```


## Notas:

- Código fonte: montecarlo.py;
- A rede pode ser quadrada (qualquer dimensão) ou barabasi-albert;
- Saídas (magnetização e susceptibilidade) são escritas em arquivos;

### Variáveis relevantes:

1. $N$ : número de vértices
1. $\alpha$: parâmetro epidemiológico. Taxa de infecção
1. $\gamma$: parâmetro epidemiológico. Taxa de cura


