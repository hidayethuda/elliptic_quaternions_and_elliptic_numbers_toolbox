# elliptic_quaternions_and_elliptic_numbers_toolbox
A MATLAB toolbox that performs advanced algebraic calculations in the algebra of elliptic quaternions and elliptic numbers has been developed.

In this project presents a toolbox developed in MATLAB version 0.1 for elliptic numbers and elliptic quaternions. Elliptic numbers and elliptic quaternions are numerical systems for solving mathematical and physical problems. On the other hand, as shown below Figure, elliptic quaternions are a general form of elliptic numbers, reduced biquaternions (commutative quaternions), complex numbers, and real numbers. This generalization offers a more expansive framework, enhancing the versatility and applicability of the developed toolbox.

![class drawio](https://github.com/user-attachments/assets/afece45f-ad78-4926-a08e-1a863120ae2c)



\begin{table}[H]
\centering

\caption{Summary of key functions provided by the elliptic number toolbox}
\renewcommand{\arraystretch}{2.3}
\resizebox{\textwidth}{!}{%
%\begin{tabular}{|l|l|}
\begin{tabular}{|l|p{13cm}|l|}
\hline
\textbf{Function name}       & \textbf{Description} \\
\hline
\texttt{ellipticcomplex2complexnumber($q_{(e)}$)}              & Returns the complex number (or matrix) to which the elliptic number  (or matrix) is isomorphic. \\ \cline{1-2}
\texttt{complexnumber2ellipticcomplex($p$-value,complex)}         & Returns the elliptic number (or matrix) to which the complex number (or matrix) is isomorphic. \\ \cline{1-2}
\texttt{elliptic\_n\_addition($q_{(e)}$,$p_{(e)}$)}          & Returns the sum of two elliptic numbers (or matrices). \\ \cline{1-2}
\texttt{elliptic\_n\_subtraction($q_{(e)}$,$p_{(e)}$)}             & Returns the difference of two elliptic numbers (or matrices). \\ \cline{1-2}
\texttt{elliptic\_n\_product($q_{(e)}$,$p_{(e)}$)}             & Returns the product of two elliptic numbers (or matrices).\\ \cline{1-2}
\texttt{elliptic\_n\_transpose($Q_{(e)}$)}              & Returns the transpose of the elliptic matrix $Q_{(e)}$. \\ \cline{1-2}
\texttt{elliptic\_n\_conjugate($Q_{(e)}$)}          & Returns the conjugate of the elliptic matrix $Q_{(e)}$. \\ \cline{1-2}
\texttt{elliptic\_n\_hermitianconjugate($Q_{(e)}$)}        & Returns the Hermitian conjugate of the elliptic matrix $Q_{(e)}$. \\ \cline{1-2}
\texttt{elliptic\_n\_eig($Q_{(e)}$)}     & [$D_{(e)}$,$V_{(e)}$] = elliptic\_n\_eig($Q_{(e)}$) returns diagonal matrix $D_{(e)}$ of  eigenvalues  and  matrix $V_{(e)}$  whose columns are the corresponding eigenvectors,  so that $Q_{(e)}V_{(e)} = V_{(e)}*D_{(e)}$. \\ \cline{1-2}
\texttt{ elliptic\_n\_svd($Q_{(e)}$)}             & [$U_e$,$S_e$,$V_e$] = elliptic\_n\_svd($Q_{(e)}$) performs a singular value  decomposition of elliptic matrix $Q_{(e)}$, such that $Q_{(e)} = U_{(e)}S_{(e)}V^*_{(e)}$.\\ \cline{1-2}
\texttt{elliptic\_n\_pinv($Q_{(e)}$)}          & Returns the pseudo inverse of the elliptic matrix $Q_{(e)}$.  \\ \cline{1-2}
\texttt{elliptic\_n\_rank($Q_{(e)}$)} & Returns the rank of the elliptic matrix $Q_{(e)}$. \\ \cline{1-2}
\texttt{elliptic\_n\_lss($Q_{(e)},P_{(e)}$)}             &  Returns the least squares solution of the elliptic matrix equation $Q_{(e)}X_{(e)}=P_{(e)}$  and the least square error. \\

\hline
\end{tabular}
\label{class1}
}
\end{table}
