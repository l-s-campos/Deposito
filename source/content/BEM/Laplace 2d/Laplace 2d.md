- [[#Formulação]]
    - [[#Método indireto para o cálculo da diagonal da matriz @import url('https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.css')HHH﻿]]
    - [[#Exemplo ]]
    - [[#Pontos internos]]
    - [[#Exercícios]]
    - [[#Desafio]]

- Código
    
    [https://github.com/l-s-campos/BEM](https://github.com/l-s-campos/BEM)  
    faça o download e rode de dentro da pasta:
    
    ```Julia
    using Pkg
    
    Pkg.activate(pwd())
    
    Pkg.instantiate()
    ```
    
    @which
    

# Formulação

A equação de Laplace é dada por:

$$\nabla^2T=0,$$

usando o método dos resíduos ponderados e aplicando a segunda identidade de Green:

$$\int_\Omega(v\nabla^2u-u\nabla^2v) d\Omega=\int_\Gamma\biggl(v\frac{\partial u}{\partial n}-u\frac{\partial v}{\partial n}\biggr)ds$$

obtém-se a equação integral de contorno:

$$0.5 T(x_{d},y_{d})=\int_{\Gamma}Tq^*ds-\int_{\Gamma}T^*qds,$$

onde $T^*$ e $q^*$ são as soluções fundamentais:

$$T^* = \frac{-1}{2\pi k}\ln r$$

$$q^*=\frac{1}{2\pi r^2}[(x-x_{d})n_{x}+(y-y_{d})n_{y}].$$

Podemos representar a temperatura e fluxo em um elemento descontinuo com $m$ nós como:$T = N_1 T_1 + N_2 T_2 + N_3 T_3 +\ldots+ N_m T_m$ $q = N_1 q_1 + N_2 q_2 + N_3 q_3 +\ldots+ N_m q_m$  
onde estamos aproximando nossa distribuição de temperatura e fluxo no elemento por uma função polinomial de grau $(m-1)$.

Discretizando em $n$ elementos de contorno descontínuos:  
$0.5T(x_{d},y_{d}) = \sum_{j = 1}^{n_{elem}} \left[ \int_{\Gamma_j} Tq^* d\Gamma \right] - \sum_{j = 1}^{n_{elem}} \left[ \int_{\Gamma_j} T^* qd\Gamma \right]$  

Usando a representação de $T$ e $q$ na equação acima temos:  

$$0.5T(x_{d},y_{d})=\sum_{j=1}^{n_{elem}} \left\{ \int_{\Gamma_j} \left[\begin{array}{rrrrr}N_1&N_2&N_3&\cdots&N_m\end{array}\right]\left[\begin{array}{c}T_1\\T_2\\T_3\\\vdots\\T_m\end{array}\right]_j q^*d\Gamma\right\} \nonumber \\ -  
\sum_{j=1}^{n_{elem}}{\left\{\int_{\Gamma_j} T^* \left[\begin{array}{rrrrr}N_1&N_2&N_3&\cdots&N_m\end{array}\right]\left[\begin{array}{c}q_1\\q_2\\q_3\\\vdots\\q_m\end{array}\right]_j d\Gamma\right\}} \nonumber$$

Repetindo essa equação para cada diferente nó podemos montar um sistema matricial:

$$HT=Gq$$

## Método indireto para o cálculo da diagonal da matriz $H$

A diagonal da matriz $H$ contém uma singularidade forte que dificulta o cálculo numérico desses termos. Uma alternativa não faz a integração de maneira explícita mas usa uma propriedade da matriz $H$ decorrente da modelagem de um corpo sob temperatura constante. Sem perder a generalidade, considere que todos os nós de um corpo encontre-se com a temperatura $T=1$. Neste caso, o fluxo será nulo em todos os nós, ou seja, $q=0$ em todos os nós. Desta forma, a equação matricial é reescrita como $H\{1\}=G\{0\}$.

Daí, os termos da diagonal da matriz $[H]$ pode ser calculado da seguinte forma:

$H_{ii}=-\sum_{j=1}^{N}{H_{ij}}, \text{ com } {\it i} \neq j, \text{ para }{\it i} = 1, 2, ..., {\it N},$

uma vez que todos os termos de fora da diagonal são integrais regulares e já foram previamente calculados.

Um procedimento parecido, usando uma outra distribuição de temperatura, pode ser feito com a diagonal da matriz $G$ mas isso normalmente não é necessário por se tratar de uma singularidade fraca

## Exemplo

A fim de ilustrar como se aplica as condições de contorno e se calcula as variáveis desconhecidas será analisado um problema de condução de calor unidirecional com uma discretização de um elemento por lado.

  

![[image 69.png|image 69.png]]

As equações obtidas podem ser escritas na forma matricial, como:

$$\left(\begin{array}{cccc}  
H_{11} & H_{12} & H_{13} & H_{14}\\H_{21} & H_{22} & H_{23} & H_{24}\\H_{31} & H_{32} & H_{33} & H_{34}\\H_{41} & H_{42} & H_{43} & H_{44}\\  
\end{array}\right)  
\left(\begin{array}{c}  
\bar{T_{1}}\\T_{2}\\\bar{T_{3}}\\T_{4}  
\end{array}\right)  
=  
\left(\begin{array}{cccc}  
G_{11} & G_{12} & G_{13} & G_{14}\\G_{21} & G_{22} & G_{23} & G_{24}\\G_{31} & G_{32} & G_{33} & G_{34}\\G_{41} & G_{42} & G_{43} & G_{44}\\  
\end{array}\right)  
\left(\begin{array}{c}  
q_{1}\\\bar{q_{2}}\\q_{3}\\\bar{q_{4}}  
\end{array}\right)$$

onde $\bar{T}$ e $\bar{q}$ são termos conhecidos.

Separando os termos conhecidos dos desconhecidos:

$$\left(\begin{array}{cccc}-G_{11} & H_{12} & -G_{13} & H_{14}\\-G_{21} & H_{22} & -G_{23} & H_{24}\\-G_{31} & H_{32} & -G_{33} & H_{34}\\-G_{41} & H_{42} & -G_{43} & H_{44}\\\end{array}\right)\left(\begin{array}{c}q_{1}\\T_{2}\\q_{3}\\T_{4}\end{array}\right)=\left(\begin{array}{cccc}-H_{11} & G_{12} & -H_{13} & G_{14}\\-H_{21} & G_{22} & -H_{23} & G_{24}\\-H_{31} & G_{32} & -H_{33} & G_{34}\\-H_{41} & G_{42} & -H_{43} & G_{44}\\\end{array}\right) \left(\begin{array}{c}\bar{T_{1}}\\\bar{q_{2}}\\\bar{T_{3}}\\\bar{q_{4}}\end{array}\right)$$

Assim, pode-se escrever $Ax=b$ e resolvendo o sistema linear calcula-se os valores das variáveis desconhecidas.

## Pontos internos

A equação integral para pontos internos é ligeiramente modificada:

$$T(x_{d},y_{d})=\int_{\Gamma}Tq^*ds-\int_{\Gamma}T^*qds$$

Aqui temos um detalhe importante, uma vez conhecida $T$ e $q$ no contorno a temperatura pode ser calculada em qualquer ponto interno sem ser necessário resolver algum sistema linear.

## Exercícios

[[erros]]

1. Resolva esse problema com elementos lineares, quadráticos e cúbicos e calcule o erro médio, erro máximo e a norma l2 do erro no contorno para diferentes discretizações. Faça um gráfico comparando a convergência dos dos 3 tipos de elementos. A solução analítica é dada por:

$$\begin{aligned}u(\theta)&=\frac{\theta}{\pi}\\q(x)&=-\frac{1}{\pi x}\quad \text{em}\quad y=0\end{aligned}$$

![[image 1 8.png|image 1 8.png]]

  

  

  

  

  

1. Analise o seguinte problema:

Analise uma placa com as seguintes condições de contorno:

- $T(x=0)=0$
- $T(x = 1) = \cos(\pi y)$
- $q(y = 0) = q(y = 1) = -k \frac{\partial u}{\partial n} = 0$

![[image 2 7.png|image 2 7.png]]

A solução analítica para este problema é dada por:

$T^{an}=\sinh(\pi x)\cos(\pi y)/\sinh(\pi)$

Analise o problema usando diferentes números de elementos e faça uma tabela mostrando o erro percentual em um ponto interno de coordenadas $(x,y)=(\sqrt{2}/2,\sqrt{2}/2)$ para os diversos casos analisados.

  

  

1. Faça um mapa de cor da distribuição de temperatura em uma placa com dimensões e condições de contorno mostradas na figura.  
    

![[image 3 7.png|image 3 7.png]]

## Desafio

Baseado nesse [artigo](https://onlinelibrary.wiley.com/doi/epdf/10.1002/fld.1650030504) calcule o coeficiente de sustentação de um perfil NACA usando o BEM.

[gravação](https://youtu.be/_G4yNayAPPE)