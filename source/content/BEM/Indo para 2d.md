## Cálculo do perímetro de figuras planas

Uma vez que a integral ao longo do contorno no $s$ pode ser bastante difícil, ou mesmo impossível, uma estratégia para o cálculo do perímetro é a divisão do contorno em uma soma de pequenos  
pedaços $s_{1},s_{2},...,s_{n}$, ou seja, $s=\sum_{i=1}^{n}s$, onde $n$ é o número de pedaços em que o contorno foi dividido. Uma vez  
que estes pedaços podem ter uma forma qualquer, cada pedaço $s_{i}$ será aproximado por uma forma conhecida. Por simplicidade, esta forma é quase sempre dada por um polinômio (linha reta, parábola, etc).

![Untitled 44.png](../attachments/Untitled%2044.png)

Dessa maneira, cada pedaço $s_{1},s_{2},...,s_{n}$ é aproximado por formas conhecidas $\Gamma_{1},\Gamma_{2},...,\Gamma_{n}$, chamados Elementos de Contorno.

Considere agora que os elementos de contorno $\Gamma_{i}$ sejam parabólicos, ou seja, são descritos por polinômios de 2ª ordem (equação de uma parábola). Desta forma, são necessários 3 pontos de $\Gamma_i$ para se definir a parábola. Estes pontos são dados por $(x_1,y_1)$, $(x_2,y_2)$ e $(x_3,y_3)$, que correspondem respectivamente a $\xi=-1$, $\xi=0$ e $\xi=1$.

![Untitled 1 30.png](../attachments/Untitled%201%2030.png)

Criando uma função parabólica para relacionar $x$ com $\xi$, tem-se: $x=a\xi^{2}+b\xi+c$ sendo que:

$$x(\xi=-1)=x_{1}=a(-1)^{2}+b(-1)+c\\x(\xi=0)=x_{2}=a(0)^{2}+b(0)+c\\x(\xi=+1)=x_{3}=a(+1)^{2}+b(+1)+c$$

resolvendo esse sistema obtem-se $x=N_{1}x_{1}+N_{2}x_{2}+N_{3}x_{3}$ onde $N_{1}=\frac{\xi}{2}(\xi-1)$, $N_{2}=(1-\xi)(1+\xi)$ e $N_{3}=\frac{\xi}{2}(\xi+1)$ são as funções de forma quadráticas contínuas.

![Untitled 2 29.png](../attachments/Untitled%202%2029.png)

A derivada de $x$ em relação a $\xi$ é dada por $\frac{dx}{d\xi}=\frac{d[N_{1}(\xi)]}{d\xi}x_{1}+\frac{d[N_{2}(\xi)]}{d\xi}x_{2}+\frac{d[N_{3}(\xi)]}{d\xi}x_{3}$ onde $\frac{dN_1}{d\xi}=\xi-\frac{1}{2},$$\frac{dN_2}{d\xi}=-2\xi$ e $\frac{dN_3}{d\xi}=\xi+\frac{1}{2}$.

O comprimento $\Gamma$ do perímetro da figura é então dado por:

$\Gamma=\sum_{i=1}^{n}\int_{-1}^{1}\frac{d\Gamma}{d\xi}d\xi=\sum_{i=1}^{n}\int_{-1}^{1}\sqrt{\left(\frac{dx}{d\xi}\right)^{2}+\left(\frac{dy}{d\xi}\right)^{2}}d\xi=\sum_{i=1}^{n}\left[\int_{-1}^{1}J(\xi)d\xi\right]$ em que a função $J(\xi)=\sqrt{\left(\frac{dx}{d\xi}\right)^{2}+\left(\frac{dy}{d\xi}\right)^{2}}$é chamada jacobiano da transformação. Usando quadratura de Gauss com $p$ pontos de Gauss, conclui-se que $\Gamma=\sum_{i=1}^{n}\left[\sum_{k=1}^{p}w_{k}J(\xi_{k})\right]$.

[código para propriedade geométrica](https://1drv.ms/f/s!AmfyGvdmTYongqYkwQK1hIXEROJbEA?e=omRWUR)

```Julia
include("dad.jl")
include("format.jl")
include("calcula_PropGeom.jl")
include("calc_fforma.jl")
using FastGaussQuadrature,GLMakie

PONTOS,SEGMENTOS,MALHA=dad_1(10,1) \#Arquivo de entrada de dados

NOS,ELEM=format_dad(PONTOS,SEGMENTOS,MALHA)# formata os dados (cria as
  # matrizes NOS e ELEM a partir das matrizes PONTOS, SEGMENTOS e MALHA)
display(mostra_geo(SEGMENTOS,PONTOS,ELEM,NOS))

Perimetro,Area,xbarra,ybarra=calcula_PropGeom(ELEM,NOS,4,4); # Calcula as propriedades geom�tricas de figuras plana
println("Perimetro calculado numericamente: $Perimetro")
println("Área calculada numericamente: $Area")
println("Centróide calculado numericamente: ( $xbarra, $ybarra)")
```

![Untitled 3 24.png](../attachments/Untitled%203%2024.png)

## Área

A integral de uma função $f(x,y)$ sobre a área $\Omega$ de uma figura plana é dada pela integral $I=\int_{\Omega}f(x,y)dxdy.$ A integral pode ser escrita, em coordenadas polares,  
como: $I=\int_{\Omega}f(x(\rho,\theta),y(\rho,\theta))\rho d\rho d\theta=\int_{\theta}\int_{0}^{r}f(x(\rho,\theta),y(\rho,\theta))\rho d\rho d\theta$ definindo $F$ como $F(\rho,\theta)=\int_{0}^{r}f(x(\rho,\theta),y(\rho,\theta))\rho d\rho$, pode-se escrever:$I=\int_{\theta}F(\rho,\theta)d\theta$.

como $\cos\alpha=\frac{r\frac{d\theta}{2}}{\frac{d\Gamma}{2}}$

$\cos\alpha=\frac{r\frac{d\theta}{2}}{\frac{d\Gamma}{2}}$$d\theta=\frac{\vec{n}.\vec{r}}{r}d\Gamma$ e finalmente $I=\int_{\Gamma}F\frac{\vec{n}.\vec{r}}{r}d\Gamma$.

Dividindo $\Gamma$ em uma soma de pequenos pedaços do contorno:

$I=\sum_{i=1}^{NE}\int_{\Gamma_{i}}F\frac{\vec{n}\cdot\vec{r}}{r}d\Gamma$.

![Untitled 4 22.png](../attachments/Untitled%204%2022.png)

![Untitled 5 19.png](../attachments/Untitled%205%2019.png)

### Cálculo do vetor normal $\vec{n}$

![Untitled 6 15.png](../attachments/Untitled%206%2015.png)

![Untitled 7 12.png](../attachments/Untitled%207%2012.png)

$\vec{S}$= vetor tangente ao contorno $\Gamma$. $\vec{S}$$=dx\vec{i}+dy\vec{j}$ vetor não unitário

$\vec{s}$= Vetor unitário tangente ao contorno $\Gamma$

$\Gamma$$\vec{s}=\frac{\vec{S}}{|\vec{S|}}=\frac{dx\vec{i}+dy\vec{j}}{\sqrt{dx^{2}+dy^{2}}}$

Dividindo tudo por $d\xi$, tem-se:

$\vec{s}=\frac{\frac{dx}{d\xi}}{\sqrt{(\frac{dx}{d\xi})^{2}+(\frac{dy}{d\xi})^{2}}}\vec{i}+\frac{\frac{dy}{d\xi}}{\sqrt{(\frac{dx}{d\xi})^{2}+(\frac{dy}{d\xi})^{2}}}\vec{j}$

$\vec{s}=\frac{\frac{dx}{d\xi}}{\frac{d\Gamma}{d\xi}}\vec{i}+\frac{\frac{dy}{d\xi}}{\frac{d\Gamma}{d\xi}}\vec{j}=s_{x}\vec{i}+s_{y}\vec{j}$

$\vec{n}=$ vetor unitário normal ao contorno $\Gamma$ apontando para fora do domínio $\Omega$

$\vec{n}=n_{x}\vec{i}+n_{y}\vec{j}$  
Temos duas possibilidades, uma vez que $\vec{s}$ é unitário:

$n_{x}= -s_{y}\quad n_{y}= s_{x}\Rightarrow(-s_{y})s_{x}+(s_{x})s_{y}=0$

$n_{x}= s_{y}\quad n_{y}= -s_{x} \Rightarrow(s_{y})s_{x}+(-s_{x})s_{y}=0$  
O vetor $\vec{s}$ é sempre definido percorrendo $\Gamma$ no sentido anti-horário, se $\Gamma$ for um contorno externo e horário, caso $\Gamma$ seja um contorno interno. Portanto, a segunda possibilidade é escolhida de modo que a normal esteja apontada para fora do domínio.

### Exercícios

![Untitled 8 10.png](../attachments/Untitled%208%2010.png)

$I_{x}=\frac{a^{4}}{96} (9\sqrt{3}-2\pi)$

![Untitled 9 8.png](../attachments/Untitled%209%208.png)

$I_x= 2 \cdot 10^4\pi - \frac{20^2 \cdot \pi}{2} \left( \frac{80}{3 \cdot \pi} \right)^2 +\left(\frac{20^2\pi}{2}\right)\left(15+\frac{80}{3\pi}\right)^2$

Calcule, modificando o programa propgeo, os momentos de inércia em relação ao eixo x para as duas figuras e compare com o resultado analítico.

  

  

- Gravações
    
    [https://youtu.be/TZvHS2JzKOs](https://youtu.be/TZvHS2JzKOs)
    
    [https://youtu.be/08vADm8JguE](https://youtu.be/08vADm8JguE)
    
    [https://youtu.be/O14y7y9ZCLs](https://youtu.be/O14y7y9ZCLs)
    
    [https://youtu.be/YbKmFAwBbwo](https://youtu.be/YbKmFAwBbwo)