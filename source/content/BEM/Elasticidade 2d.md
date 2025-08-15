- [Equações importantes](#Equações%20importantes)
    - [Exercícios](#Exercícios)
        - [Cilindro pressurizado](#Cilindro%20pressurizado)
        - [Placa com furo](#Placa%20com%20furo)
        - [Viga](#Viga)

# Equações importantes

A relação de equilíbrio de tensão pode ser escrita como:

$$\frac{\partial\sigma_{xx}}{\partial x}+\frac{\partial\sigma_{xy}}{\partial y}+f_{x}=0\\\frac{\partial\sigma_{yx}}{\partial x}+\frac{\partial\sigma_{yy}}{\partial y}+f_{y}=0$$

$$\frac{\partial\sigma_{ij}}{\partial x_j}+f_i=0$$

Relação deformação-deslocamento pode ser escrita como:

$$\varepsilon_{xx}=\frac{\partial u_{x}}{\partial x};\quad\varepsilon_{yy}=\frac{\partial u_{y}}{\partial y};\quad\varepsilon_{xy}=\frac{1}{2}\left(\frac{\partial u_{x}}{\partial y}+\frac{\partial u_{y}}{\partial x}\right)$$

$$\varepsilon_{ij}=\frac{1}{2}\biggl(\frac{\partial u_{i}}{\partial x_{j}}+\frac{\partial u_{j}}{\partial x_{i}}\biggr)$$

Relação tensão-deformação pode ser escrita como:

$$\begin{aligned}&\mathbf{\varepsilon}_{xx} =\left(\frac{1}{E}\right)\sigma_{xx}+\left(\frac{-v}{E}\right)\sigma_{yy} \\&\varepsilon_{yy} =\left(\frac{-v}{E}\right)\sigma_{xx}+\left(\frac{1}{E}\right)\sigma_{yy} \\&\varepsilon_{xy} =\frac1{2\mu}\sigma_{xy} \end{aligned}$$

$$\mu=\frac E{2(1+v)}$$

Essas três relações podem ser usadas para chegar na seguinte equação diferencial de deslocamentos:

$$\frac{\partial^2u_x}{\partial x^2}+\frac{\partial^2u_x}{\partial y^2}+\frac1{1-2v}\Bigg(\frac{\partial^2u_x}{\partial x^2}+\frac{\partial^2u_y}{\partial x \partial y}\Bigg)=\frac{-f_x}\mu\\\frac{\partial^2u_y}{\partial x^2}+\frac{\partial^2u_y}{\partial y^2}+\frac1{1-2v}\Bigg(\frac{\partial^2u_y}{\partial y^2}+\frac{\partial^2u_x}{\partial x \partial y}\Bigg)=\frac{-f_y}\mu$$

$$\frac{\partial^2u_i}{\partial x_j\partial x_j}+\left(\frac1{1-2v}\right)\frac{\partial^2u_j}{\partial x_i\partial x_j}=\frac{-f_i}\mu$$

A solução fundamental dessa equação de deslocamento é dada por:

$$U_{ij}(p,Q)=\frac{1}{8\pi\mu(1-v)}\left\{(3-4v)\ln\left[\frac{1}{r(p,Q)}\right]\delta_{ij}+\frac{\partial r(p,Q)}{\partial x_{i}}\frac{\partial r(p,Q)}{\partial x_{j}}\right\}$$

e para a tração:

$$\begin{aligned}  
T_{ij}(p,Q)& =\frac{-1}{4\pi(1-\nu)r(p,Q)}{\left(\frac{\partial r(p,Q)}{\partial n}\right)} \\  
&\times\left[(1-2\nu)\delta_{ij}+2\frac{\partial r(p,Q)}{\partial x_i}\frac{\partial r(p,Q)}{\partial x_j}\right] \\  
&+\frac{1-2v}{4\pi(1-v)r(p,Q)}\biggl[\frac{\partial r(p,Q)}{\partial x_j}n_i-\frac{\partial r(p,Q)}{\partial x_i}n_j\biggr]  
\end{aligned}$$

Usando um procedimento semelhante ao apresentado para o problema potencial obtemos a equação integral de contorno:

$$\begin{aligned}\begin{bmatrix}u_x(p)\\u_y(p)\end{bmatrix}&+\int_\Gamma\begin{bmatrix}T_{xx}(p,Q)&T_{xy}(p,Q)\\T_{yx}(p,Q)&T_{yy}(p,Q)\end{bmatrix}\begin{bmatrix}u_x(Q)\\u_y(Q)\end{bmatrix}\mathrm{d}\Gamma(Q)\\&=\int_\Gamma\begin{bmatrix}U_{xx}(p,Q)&U_{xy}(p,Q)\\U_{yx}(p,Q)&U_{yy}(p,Q)\end{bmatrix}\begin{bmatrix}t_x(Q)\\t_y(Q)\end{bmatrix}\mathrm{d}\Gamma(Q)\end{aligned}$$

$$u_i(p)+\int_\Gamma T_{ij}(p,Q)u_j(Q)\operatorname{d}\Gamma(Q)=\int_\Gamma U_{ij}(p,Q)t_j(Q)\operatorname{d}\Gamma(Q)$$

Essa equação pode ser derivada e junto com as relações tensão deformação obtém-se a equação que pode ser usada para calcular a tensão:

$$\begin{gathered}\sigma_{ij}(p)+\int_{\Gamma}\left\{\frac{2\mu\nu}{1-2\nu}\delta_{ij}\frac{\partial T_{mk}(p,Q)}{\partial x_{m}}\right. \\+\mu\left[\frac{\partial T_{ik}(p,Q)}{\partial x_j}+\frac{\partial T_{jk}(p,Q)}{\partial x_i}\right]\biggr\}u_k(Q) \mathrm{d}\Gamma(Q) \\=\int_\Gamma\left\{\frac{2\mu\nu}{1-2\nu}\delta_{ij}\frac{\partial U_{mk}(p,Q)}{\partial x_m}\right.+\mu\biggl[\frac{\partial U_{ik}(p,Q)}{\partial x_j}+\frac{\partial U_{jk}(p,Q)}{\partial x_i}\biggr]\biggr\}t_k(Q) \mathrm{d}\Gamma(Q)\end{gathered}$$

que pode ser reescrita como:

$$\sigma_{ij}(p)+\int_\Gamma S_{kij}(p,Q)u_k(Q)\mathrm{~d}\Gamma(Q)=\int_\Gamma D_{kij}(p,Q)t_k(Q)\mathrm{~d}\Gamma(Q)$$

onde os tensores $S_{kij}$ e $D_{kij}$ são dados por:

$$\begin{aligned}  
S_{kij}(p,Q)& =\frac\mu{2\pi(1-\nu)}\left(\frac1{r^2}\right)n_i\left[2\nu\frac{\partial r}{\partial x_j}\frac{\partial r}{\partial x_k}+(1-2\nu)\delta_{jk}\right] \\  
&+\frac\mu{2\pi(1-\nu)}\left(\frac1{r^2}\right)n_j\left[2\nu\frac{\partial r}{\partial x_i}\frac{\partial r}{\partial x_k}+(1-2\nu)\delta_{ik}\right] \\  
&+\frac\mu{2\pi(1-v)}\biggl(\frac1{r^2}\biggr)n_k\biggl[2(1-2v)\frac{\partial r}{\partial x_i}\frac{\partial r}{\partial x_j}-(1-4v)\delta_{ij}\biggr] \\  
&+\frac\mu{\pi(1-\nu)}\biggl(\frac1{r^2}\biggr)\biggl(\frac{\partial r}{\partial n}\biggr)\biggl[(1-2\nu)\delta_{ij}\frac{\partial r}{\partial x_k}+\nu\biggl(\delta_{jk}\frac{\partial r}{\partial x_i}+\delta_{ik}\frac{\partial r}{\partial x_j}\biggr) \\  
&-4\frac{\partial r}{\partial x_i}\frac{\partial r}{\partial x_j}\frac{\partial r}{\partial x_k}\biggl]\end{aligned}$$

$$D_{kij}(p,Q)=\frac{1}{4\pi(1-\nu)}\biggl(\frac{1}{r}\biggr)\biggl[(1-2\nu)\biggl(\delta_{jk}\frac{\partial r}{\partial x_{i}}+\delta_{ik}\frac{\partial r}{\partial x_{j}}-\delta_{ij}\frac{\partial r}{\partial x_{k}}\biggr)\\+2\frac{\partial r}{\partial x_{i}}\frac{\partial r}{\partial x_{j}}\frac{\partial r}{\partial x_{k}}\biggr]$$

## Exercícios

### Cilindro pressurizado

Raio interno $R_a = 50 mm$

Raio externo $R_b= 100 mm$

Modulo de elasticidade $E= 200 GPa$

Poisson $\nu=0.32$

Pressão $P=100 N/mm$

  

![image 71.png](../attachments/image%2071.png)

$$\begin{gathered}  
u_r=\frac{(1+v) p r_a^2}{\left(r_b^2-r_a^2\right) E}\left[(1-2 v) r+\frac{r_b^2}{r}\right] \\  
\sigma_r=\frac{p r_a^2}{r_b^2-r_a^2}-\frac{r_a^2 r_b^2 p}{r_b^2-r_a^2} \frac{1}{r^2} \\  
\sigma_h=\frac{p r_a^2}{r_b^2-r_a^2}+\frac{r_a^2 r_b^2 p}{r_b^2-r_a^2} \frac{1}{r^2}  
\end{gathered}$$

### Placa com furo

Raio $R = 50 mm$

Modulo de elasticidade $E= 100 GPa$

Poisson $\nu=0.25$

Pressão $P=1 N/mm$

Esse problema é de estado plano de tensão. Para ser tratado pelas soluções fundamentais de estado plano de deformação as propriedades tem de ser ajustadas:

![image 1 10.png](../attachments/image%201%2010.png)

$$\nu^{'}=\frac{\nu}{1+\nu}\\E^{'}=E\left[1-\frac{\nu^{'2}}{(1+\nu^{'})^{2}}\right]$$

A solução analítica é dada por:

$$\begin{array}{c}{{\sigma_{r}=\frac{\sigma}{2}\left(1-\frac{a^{2}}{r^{2}}\right)+\frac{\sigma}{2}\left(1+\frac{3a^{4}}{r^{2}}-\frac{4a^{2}}{r^{2}}\right)\cos(2\theta)}}\\ {\sigma_{\theta}=\frac{\sigma}{2}\left(1+\frac{a^{2}}{r^{2}}\right)-\frac{\sigma}{2}\left(1+\frac{3a^{4}}{r^{4}}\right)\cos(2\theta)}\\{\tau_{r\theta}=-\frac{\sigma}{2}\left(1-\frac{3a^{4}}{r^{4}}+\frac{4a^{2}}{r^{2}}\right)\sin(2\theta)}\end{array}$$

### Viga

Uma viga com um carregamento parabólico $t_2(y)=-\frac{P}{2I}\left(\frac{D^2}{4}-y^2\right)$ tem solução analítica:

Comprimento $L= 48 mm$

Altura $D= 12 mm$

Modulo de elasticidade $E= 300 GPa$

Poisson $\nu=0.3$

Pressão $P=1000 N$

![image 2 9.png](../attachments/image%202%209.png)

$$u_1(x,y)=-\frac{Py}{6EI}\left[(6L-3x)x+(2+v)\left(y^2-\frac{D^2}{4}\right)\right]$$

$$u_2(x,y)=\frac{P}{6EI}\left[(3v)y^2(L-x)+(4+5v)\frac{D^2x}{4}+(3L-x)x^2\right]$$

$$\sigma_{xx}(x,y)=-\frac{P(L-x)y}{I}\\\tau_{xy}(x,y)=-\frac{P}{2I}\left(\frac{D^{2}}{4}-y^{2}\right)$$

[gravação](https://youtu.be/R6-_ECEQXRk)