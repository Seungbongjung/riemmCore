# riemmCore

## Partial Isotropy Core Shrinkage Estimator (PICSE)

Suppose $\Sigma$ is a $p_1p_2\times p_1p_2$ covariance matrix of $p_1\times p_2$ matrix-variate data with a partial-isotropy core. Namely, 

$$
\Sigma=K^{1/2}((1-\lambda)AA^\top+\lambda I_p)K^{1/2,\top},
$$

where $K^{1/2}$ is either $L_2\otimes L_1$ or $S_2\otimes S_1$, where $L_i$ (resp. $S_i$) is a $p_i\times p_i$ Cholesky factor (resp. strictly positive definite matrix). Also, $\lambda\in(0,1)$ and $A\in\mathbb{R}^{p\times r}$ of full-column rank with $p:=p_1p_2$, $p_1/p_2+p_2/p_1<r\ll p$, and the Kronecker MLE being $I_p$. The separable component of $\Sigma$, representing the most separable part of $\Sigma$, is $K^{1/2}K^{1/2,\top}$, whereas the core component, representing the non-separable part, is $(1-\lambda)AA^\top+\lambda I_p$ in terms of Kronecker-core decomposition (KCD) proposed by Hoff et al. (2023).    

Assume $Y_1,\ldots,Y_n\overset{i.i.d.}{\sim}N_{p_1\times p_2}(0,\Sigma)$. Then we propose a partial isotropy core shrinkage estimator ($\texttt{PICSE}$) to compute the MLE of $\Sigma$ via alternating minimization of the negative log-likelihood in each parameter constituting $\Sigma$. Specifically, we implement Algorithm 6.1 from Sung (2025) and provide a replication code for Section 7 of Sung (2025).

For some parameters, we adopt a second-order Riemannian optimization. For efficient optimization, we use the $\texttt{R}$ code, $\texttt{lsmr.R}$, based on the conjugate-gradient type method for solving linear equations, originally available from $\texttt{R}$ package $\texttt{sTR}$ (https://rdrr.io/cran/stR/src/R/lsmr.R). 

## Reference 

Hoff, McCormack and Zhang (2023), Core Shrinkage Covariance Estimation for Matrix-variate Data J. R. Stat. Soc., B: Stat., 85 (2023), pp. 1659--1679.

Sung (2025), sfdfsdf
