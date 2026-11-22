# Reference answers: open after your independent calculation

These answers correspond to `element-hand-check.md` and `element_demo`.
They are AI-provided reference values, not evidence of independent human review.

## Case A

P0=(0,0), P1=(1,0), P2=(0,1); f=1.
B and B^(-T) are the identity, D=1, A_T=1/2.
G0=(-1,-1), G1=(1,0), G2=(0,1).

\[
E=\begin{pmatrix}1&-1/2&-1/2\\-1/2&1/2&0\\-1/2&0&1/2\end{pmatrix},
\qquad r=(1/6,1/6,1/6)^T.
\]

## Case B

P0=(1,1), P1=(3,1), P2=(2,2); f=1+x+2*y.

\[
B=\begin{pmatrix}2&1\\0&1\end{pmatrix},\quad D=2,\quad A_T=1,
\]

\[
B^{-1}=\begin{pmatrix}1/2&-1/2\\0&1\end{pmatrix},\quad
B^{-T}=\begin{pmatrix}1/2&0\\-1/2&1\end{pmatrix}.
\]

\[
G_0=(-1/2,-1/2)^T,\quad G_1=(1/2,-1/2)^T,\quad G_2=(0,1)^T.
\]

\[
E=\begin{pmatrix}1/2&0&-1/2\\0&1/2&-1/2\\-1/2&-1/2&1\end{pmatrix}.
\]

Since A_T=1, the dot-product table H equals E in this case.

| q | Q_q | F_q | Physical weight abs(D)*w_q |
|---|---|---|---|
| 0 | (3/2,7/6) | 29/6 | 1/3 |
| 1 | (5/2,7/6) | 35/6 | 1/3 |
| 2 | (2,5/3) | 19/3 | 1/3 |

The nine load contributions C_iq=abs(D)*w_q*F_q*Phi_iq are

\[
C=\begin{pmatrix}
29/27&35/108&19/54\\
29/108&35/27&19/54\\
29/108&35/108&38/27
\end{pmatrix}.
\]

Sum each row:

\[
r=(7/4,\ 23/12,\ 2)^T,
\qquad \sum_i r_i=17/3.
\]

Independent affine-source check: V0=4, V1=6, V2=7.
Thus r0=(2*4+6+7)/12=7/4,
r1=(4+2*6+7)/12=23/12,
and r2=(4+6+2*7)/12=2.

Decimal values: r=(1.75, 1.9166666666666667, 2).
The local matrix has zero row sums and is singular, as expected before global
assembly and boundary conditions.
