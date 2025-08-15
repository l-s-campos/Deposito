- [MECID](#MECID)
    - [Exercícios](#Exercícios)
    - [Gravação](#Gravação)

# MECID

A equação de Poisson é dada por:

$$\nabla^2T=f(x,y),$$

um novo termo aparece na equação integral de contorno:

$$0.5 T(x_{d},y_{d})=\int_{\Gamma}Tq^*ds-\int_{\Gamma}T^*qds+\red{\int_{\Omega}T^*fd\Omega},$$

esse termo precisa de um tratamento extra. Vamos aproximar $T^*f$ por funções de base radial:

$$T^*f=\sum_{j=1}^n\phi(r_j)\alpha_j.$$

$\phi$ é uma função de base radial que será tomada aqui por $r^2 \log(r)$. $\alpha$ são coeficientes desconhecidos que podem ser obtidos aplicando essa equação $n$ vezes e montando um sistema matricial:

$$\begin{bmatrix}\phi(\mathrm X_1-\mathrm X_1)&\phi(\mathrm X_1-\mathrm X_2)&\cdots&\phi(\mathrm X_1-\mathrm X_N)\\\phi(\mathrm X_2-\mathrm X_1)&\phi(\mathrm X_2-\mathrm X_2)&\cdots&\phi(\mathrm X_2-\mathrm X_N)\\\vdots&\vdots&\ddots&\vdots\\\phi(\mathrm X_N-\mathrm X_1)&\phi(\mathrm X_N-\mathrm X_2)&\cdots&\phi(\mathrm X_N-\mathrm X_N)\end{bmatrix}\begin{Bmatrix}\alpha_1\\\alpha_2\\\vdots\\\alpha_N\end{Bmatrix}=\begin{Bmatrix}\mathrm f(\mathrm X_1)T^*(X_d,X_1)\\\mathrm f(\mathrm X_2)T^*(X_d,X_2)\\\vdots\\\mathrm f(\mathrm X_N)T^*(X_d,X_N)\end{Bmatrix}=\begin{Bmatrix}\mathrm T^*(X_d,X_1)&\cdots&0\\0&\mathrm T^*(X_d,X_2)&0\\\vdots&\ddots&\vdots\\0&0&\mathrm T^*(X_d,X_N)\end{Bmatrix}\begin{Bmatrix}\mathrm f(\mathrm X_1)\\\mathrm f(\mathrm X_2)\\\vdots\\\mathrm f(\mathrm X_N)\end{Bmatrix}$$

escrito de maneira sintética:

$$\alpha=F^{-1}Df.$$

Voltando a integral de domínio:

$$\int_{\Omega}T^*fd\Omega=\int_{\Omega}\sum_{j=1}^n\phi(r_j)\alpha_jd\Omega=\sum_{j=1}^n\int_{\Omega}\phi(r_j)d\Omega F^{-1}Df=sF^{-1}Df=s_fDf$$

$$\int_{\Omega}T^*fd\Omega=\begin{Bmatrix}\mathrm s_1&\mathrm s_2&\dots&\mathrm s_n\end{Bmatrix}F^{-1}\begin{Bmatrix}\mathrm T^*(X_d,X_1)&\cdots&0\\0&\mathrm T^*(X_d,X_2)&0\\\vdots&\ddots&\vdots\\0&0&\mathrm T^*(X_d,X_N)\end{Bmatrix}\begin{Bmatrix}\mathrm f(\mathrm X_1)\\\mathrm f(\mathrm X_2)\\\vdots\\\mathrm f(\mathrm X_N)\end{Bmatrix}=\begin{Bmatrix}\mathrm s_{f1}\mathrm T^*(X_d,X_1)&\mathrm s_{f2}\mathrm T^*(X_d,X_2)&\dots&\mathrm s_{fN}\mathrm T^*(X_d,X_N)\end{Bmatrix}\begin{Bmatrix}\mathrm f(\mathrm X_1)\\\mathrm f(\mathrm X_2)\\\vdots\\\mathrm f(\mathrm X_N)\end{Bmatrix}$$

Fazendo isso para os $n$ pontos fontes:

$$\begin{Bmatrix}\mathrm s_{f1}\mathrm T^*(X_1,X_1)&\mathrm s_{f2}\mathrm T^*(X_1,X_2)&\dots&\mathrm s_{fN}\mathrm T^*(X_1,X_N)\\\mathrm s_{f1}\mathrm T^*(X_2,X_1)&\mathrm s_{f2}\mathrm T^*(X_2,X_2)&\dots&\mathrm s_{fN}\mathrm T^*(X_2,X_N)\\\vdots\\\mathrm s_{f1}\mathrm T^*(X_N,X_1)&\mathrm s_{f2}\mathrm T^*(X_N,X_2)&\dots&\mathrm s_{fN}\mathrm T^*(X_N,X_N)\end{Bmatrix}\begin{Bmatrix}\mathrm f(\mathrm X_1)\\\mathrm f(\mathrm X_2)\\\vdots\\\mathrm f(\mathrm X_N)\end{Bmatrix}=Mf$$

a diagonal dessa matriz não é bem definida devido a singularidade da solução fundamental. Ela será calculada de maneira indireta como feito com a matriz $H$.

Esse procedimento tem de ser capaz de calcular essa integral quando a função $f$ for constante e igual a 1.

$$I_{1d}=\int_{\Omega}T^*(X_d,X)d\Omega$$

Igualando as equações:

$$M[1]=[I_1]$$

a diagonal da matriz $M$ tem de ser dada por

$$M_{ii}=I_{1i}-\sum_{j=1}^{N}{M_{ij}}, \text{ com } {\it i} \neq j, \text{ para }{\it i} = 1, 2, ..., {\it N},$$

```Julia
Ht, Gt = BEM.calc_HeGt(dad)
A, b = BEM.aplicaCDC(Ht, Gt, dad)

M = BEM.Monta_M_RIMd(dad, npg)# calc_HeG_potencial linha 310
f=ones(nc(dad)+ni(dad))*10
x = A \ (b+M*f)\#carga distribuída
x = A \ (b+M*u̇)\#difusão
T, q = separa(dad, x) \#format 479
Ti=x[nc(dad)+1:end]
```

## Exercícios

Todos os exercícios tem de ser resolvidos com diferentes discretizações.

1. Determine a superfície de deflexão de uma membrana elástica em forma de um triângulo equilátero com comprimento lateral a = 5,0 m. A membrana está fixa ao longo de sua borda e é submetida a uma carga distribuída uniformemente f = 10 kN/m² e uma tensão S = 1 kN/m. Os eixos coordenados são tomados como mostrado na figura.

A equação que rege esse problema é: $S\nabla^2w=−f$

![image 70.png](../attachments/image%2070.png)

A solução analítica é dada por:

$$w = -\frac{f}{2S} \left[ \frac{1}{2} \left( x^2 + y^2 \right) - \frac{1}{a\sqrt{3}} \left( y^3 - 3x^2 y \right) - \frac{1}{18} a^2 \right]$$

Calcule o erro médio e a norma L2 do erro em 100 pontos internos. Faça um gráfico com a distribuição de deflexão.

1. Considere um cubo unitário inicialmente a temperatura zero e submetido a um aquecimento súbito em uma de suas faces.

![image 1 9.png](../attachments/image%201%209.png)

A face aquecida é elevada e mantida a uma temperatura unitária. Supõe-se que as propriedades do material sejam unitárias. A solução analítica para este problema de exemplo pode ser encontrada como:

$$T(Y,t)=1-\frac{4}{\pi}\sum_{n=0}^{\infty}\frac{(-1)^{n}}{2n+1}\exp\left\{-\frac{(2n+1)^{2}\pi^{2}\kappa t}{4L^{2}}\right\}\cos\frac{(2n+1)\pi Y}{2L}$$

Usando algum método de análise transiente da aula 2, resolva esse problema e compare com a solução analítica.

1. Considere uma elipse que é governada por: $\nabla^{2}u=4-x^{2}$

a solução analítica é dada por:

  

$u=\left[1.6-{\frac{1}{246}}\left(50x^{2}-8y^{2}+33.6\right)\right]\left({\frac{x^{2}}{4}}+y^{2}-1\right)$e $\begin{array}{c}{{q=0.4\left(x^{2}+8y^{2}\right)+\frac{1}{246}\left(-50x^{3}-96x y^{2}+83.2x\right)\frac{x}{2}+}}\\ {{\frac{1}{246}\left(-96x^{2}y+32y^{3}-83.2y\right)y}}\end{array}$calcule os erros do potencial nos pontos internos e do fluxo no contorno.

  

![image 2 8.png](../attachments/image%202%208.png)

extra

Um cilindro oco com temperatura inicial zero é considerado. O raio interno a=1 e o raio externo b = 2. A superfície interna do cilindro é mantida a uma temperatura=1. Nesse caso, são impostas restrições de simetria. As propriedades do material são unitárias. A solução analítica é dada por:

$$T(r,t)=\frac{\mathrm{ln}(b/r)}{\mathrm{ln}(b/a)}T_{1}+\pi\sum_{n=1}^{\infty}\frac{J_{0}(b\alpha_{n})J_{0}(a\alpha_{n})}{J_{0}^{2}(a\alpha_{n})-J_{0}^{2}(b\alpha_{n})}\{J_{0}(r\alpha_{n})Y_{0}(b\alpha_{n})-J_{0}(b\alpha_{n})Y_{0}(r\alpha_{n})\}\mathrm{e}^{-i\alpha_{n}^{2}},$$

onde $T_i$ é a temperatura constante na superfície interna, $J_0$ e $Y_0$ são as funções de Bessel de primeira e segunda espécies, respectivamente, e $\alpha$ é a raiz de:

$$J_{0}(a x)Y_{0}(b x)-J_{0}(b x)Y_{0}(a x)=0.$$

você pode usar a função [fzeros](https://juliamath.github.io/Roots.jl/stable/roots/#Searching-for-all-zeros-in-an-interval) para encontrar $\alpha$

## Gravação

[https://youtu.be/uSREar_ejnM](https://youtu.be/uSREar_ejnM)