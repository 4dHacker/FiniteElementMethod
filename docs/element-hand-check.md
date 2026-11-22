# P1 triangle: step-by-step hand verification

This worksheet verifies only the local element for -Laplace(u)=f, with K=1.
It follows the reference triangle and three interior quadrature points in
Peter Bastian's *Minimal Finite Element Formulation*, 4 November 2026.
No global matrix, boundary conditions or linear solve are performed here.

First work through case A. Then calculate case B on paper before opening
`element-reference-results.md` or running the skew demo. Retain exact fractions
until the final comparison. The AI-generated tests are supporting checks, not
a substitute for your independent calculation.

## 1. What is given, computed, or still unknown?

| Symbol | Type / size | Meaning | Status at this stage |
|---|---|---|---|
| T | One physical triangle | The element being checked | Given |
| P_i=(x_i,y_i), i=0,1,2 | Three 2-vectors | Ordered physical vertices | Given |
| K | Scalar | Diffusion coefficient | Fixed to 1 |
| f(x,y) | Scalar function | Source/load in the PDE | Given |
| xi, eta | Two scalars | Reference coordinates, not physical x,y | Independent reference variables |
| phi_hat_i | Three known functions | Reference shape functions | Given below |
| G_hat_i | Three 2-vectors | Reference gradients | Given below |
| q=0,1,2 | Index | Which integration point | Summation index |
| xi_q, eta_q, w_q | Three points and weights | Reference quadrature rule | Given below |
| B | 2 by 2 matrix | Reference-to-physical linear map | Compute |
| D | Signed scalar | det(B), equals twice signed area | Compute |
| A_T | Positive scalar | Physical triangle area | Compute |
| B^(-T) | 2 by 2 matrix | Inverse transpose of B | Compute |
| G_i | Three 2-vectors | Physical gradients of shape functions | Compute |
| Q_q | Three 2-vectors | Physical integration points | Compute |
| F_q=f(Q_q) | Three scalars | Source sampled in the physical triangle | Compute from the given f |
| Phi_iq | 3 by 3 table | Shape i at reference point q | Compute once; same for all triangles |
| E_ij | 3 by 3 matrix | Local stiffness matrix | Main output |
| r_i | Length-3 vector | Local load vector | Main output |
| u_i | Three coefficients | Values of the unknown solution at vertices | Still unknown; NOT solved here |

The approximation on T will eventually be u_h=sum_i u_i*phi_i. E and r are
coefficients for global assembly, not the unknown temperature values.
The local matrix is singular: do not try to solve E*u=r as a standalone PDE.

## 2. Fixed reference data

The reference vertices are (0,0), (1,0), (0,1).

\[
\hat\phi_0=1-\xi-\eta,\quad \hat\phi_1=\xi,\quad \hat\phi_2=\eta.
\]

\[
\hat G_0=(-1,-1)^T,\quad\hat G_1=(1,0)^T,\quad\hat G_2=(0,1)^T.
\]

Use column vectors for gradients. Define the integration data:

| q | xi_q | eta_q | w_q |
|---|---|---|---|
| 0 | 1/6 | 1/6 | 1/6 |
| 1 | 2/3 | 1/6 | 1/6 |
| 2 | 1/6 | 2/3 | 1/6 |

Evaluate all three shapes at all three points:

\[
\Phi=\begin{pmatrix}
2/3&1/6&1/6\\
1/6&2/3&1/6\\
1/6&1/6&2/3
\end{pmatrix}.
\]

Rows are shape index i; columns are quadrature index q. Every column sums to 1.
The weights sum to 1/2, not 1: they already contain the reference area.

## 3. Case A: worked example

Given P0=(0,0), P1=(1,0), P2=(0,1), and f(x,y)=1.

### Step A1: construct B

Place the two edge vectors in the COLUMNS:

\[
B=[P_1-P_0\;P_2-P_0]
=\begin{pmatrix}x_1-x_0&x_2-x_0\\y_1-y_0&y_2-y_0\end{pmatrix}
=\begin{pmatrix}1&0\\0&1\end{pmatrix}.
\]

The map is Q=P0+B*(xi,eta)^T. As a check, the three reference vertices must
map to P0, P1 and P2, respectively.

### Step A2: compute determinant and area

For B=[[a,b],[c,d]], D=a*d-b*c and A_T=abs(D)/2.
Here D=1 and A_T=1/2. Positive D means counter-clockwise vertex order.

### Step A3: invert, then transpose

\[
B^{-1}=\frac1D\begin{pmatrix}d&-b\\-c&a\end{pmatrix},
\qquad
B^{-T}=\frac1D\begin{pmatrix}d&-c\\-b&a\end{pmatrix}.
\]

Here both are the identity. They are generally different matrices; case B
will check that difference.

### Step A4: transform each gradient

Use G_i=B^(-T)*G_hat_i. Multiply a 2 by 2 matrix by a column vector:

\[
G_0=(-1,-1)^T,\qquad G_1=(1,0)^T,\qquad G_2=(0,1)^T.
\]

Check G0+G1+G2=(0,0). This follows from phi0+phi1+phi2=1.

### Step A5: compute the local stiffness matrix

First form all dot products G_i dot G_j:

\[
H=\begin{pmatrix}2&-1&-1\\-1&1&0\\-1&0&1\end{pmatrix}.
\]

Then E=A_T*H:

\[
E=\begin{pmatrix}1&-1/2&-1/2\\-1/2&1/2&0\\-1/2&0&1/2\end{pmatrix}.
\]

For example, E01=(1/2)*[(-1)*1+(-1)*0]=-1/2.
Check all nine entries, symmetry and zero row sums. These properties alone
are not enough to prove the implementation correct; also compare the values.

### Step A6: map the quadrature points and evaluate f

Compute Q_q=P0+B*(xi_q,eta_q)^T. Since B is the identity and P0=0:

| q | Q_q | F_q=f(Q_q) |
|---|---|---|
| 0 | (1/6,1/6) | 1 |
| 1 | (2/3,1/6) | 1 |
| 2 | (1/6,2/3) | 1 |

In general, f must be evaluated at the PHYSICAL point Q_q, while Phi_iq
is evaluated at the REFERENCE point (xi_q,eta_q).

### Step A7: compute the local load vector

For each i, sum over q:

\[
r_i=|D|\sum_{q=0}^2 w_qF_q\Phi_{iq}.
\]

For i=0, calculate the three contributions separately:

\[
r_0=1\left[\frac16\cdot1\cdot\frac23
+\frac16\cdot1\cdot\frac16
+\frac16\cdot1\cdot\frac16\right]
=\frac19+\frac1{36}+\frac1{36}=\frac16.
\]

Repeat for i=1 and i=2. The result is r=(1/6,1/6,1/6).
Independent check for constant f: each entry must equal f*A_T/3, and the
sum of entries must equal f*A_T=1/2.
Do not multiply by A_T again: abs(D) and the reference weights already
produce the physical integration measure.

## 4. Case B: independent worksheet

Given P0=(1,1), P1=(3,1), P2=(2,2), and f(x,y)=1+x+2*y.
Use exactly this vertex order. This triangle is translated, stretched and
sheared; its inverse is not its inverse transpose.

Follow the same seven steps. Record the following before looking at answers:

| Step | Quantity to write down | Your result |
|---|---|---|
| B1 | P1-P0 and P2-P0 | __________ |
| B1 | B, with edge vectors as columns | __________ |
| B2 | D and A_T | __________ |
| B3 | B^(-1), then B^(-T) | __________ |
| B4 | G0, G1, G2, then their sum | __________ |
| B5 | The 3 by 3 dot-product table H | __________ |
| B5 | E=A_T*H | __________ |
| B6 | Q0, Q1, Q2 in physical coordinates | __________ |
| B6 | F0, F1, F2 from 1+x+2*y | __________ |
| B7 | Each of the nine terms abs(D)*w_q*F_q*Phi_iq | __________ |
| B7 | r0, r1, r2 and their sum | __________ |

For example, start B6 by writing

\[
Q_0=P_0+B\begin{pmatrix}1/6\\1/6\end{pmatrix},
\qquad F_0=1+(Q_0)_x+2(Q_0)_y.
\]

For B7 use the same Phi table as in case A. Changing the physical triangle
does not change the reference shape-value table or reference weights.

An independent exact check is available because f is affine. Define vertex
source values V_i=f(P_i), which are NOT the quadrature values F_q. Then

\[
r_i=\frac{A_T}{12}(2V_i+V_j+V_k),\quad\{i,j,k\}=\{0,1,2\},
\qquad \sum_i r_i=\frac{A_T}{3}(V_0+V_1+V_2).
\]

This follows from integral(phi_i^2)=A_T/6 and
integral(phi_i*phi_j)=A_T/12 for i != j. The three-point rule is exact here
because f*phi_i has degree 2. It need not be exact for a nonlinear source.

## 5. Compare paper and program

In the repository root, build and run tests:

```powershell
.\scripts\build-element.ps1
```

This deliberately does not print reference answers. After your calculation:

```powershell
$env:PATH = "C:\msys64\mingw64\bin;$env:PATH"
.\build\element_demo.exe reference
.\build\element_demo.exe skew
```

The demo prints B, signed determinant, area, inverse transpose, gradients,
physical quadrature points, sampled source values, shape values, E and r.
Compare them in that order; the first difference normally identifies the
incorrect step. Use exact fractions or at least 12 significant digits when
comparing at the automated test tolerance of 1e-12*max(1,abs(expected)).
If you round paper calculations earlier, record the rounding and use a
correspondingly looser comparison rather than treating it as a code failure.

For an orientation check, exchange P1 and P2. Area must remain positive;
D changes sign. Rows and columns 1 and 2 of E, and entries 1 and 2 of r,
must exchange places. Compare after that permutation, not entry-for-entry
with the original local numbering. The implementation preserves input order.

## 6. Record the independent review

Keep your dated calculations or a scan. Record the case, your values for B,
B^(-T), E and r, the program output, discrepancies and whether each was
resolved. A useful review statement is: "Independently checked the geometry,
gradient transformation, stiffness and load of cases A and B."

The implementation rejects non-finite input/results, zero area, and numerically
flat triangles. Its relative flatness threshold is 64 times double machine
epsilon after scaling edge components. It also rejects unrepresentable area
or intermediate arithmetic; it does not support arbitrary extreme coordinates.

No full solver or global boundary treatment is implemented by this milestone.
