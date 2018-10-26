function [S] = Helmholtz_Sv222noah(bc,k,xi,Lref,nz,lz)
%
% Solves:  
%        Lref*Lref d/dz{xi*dS/dz }  -  S/k   =  f
%        where S = darcy flux at cell centers
%
%        where Lref = sqrt(kref*xiref)and f = - 1
%

% -- Mesh discretization --
cdz  = nz/lz;
cdz2 = cdz^2;
Lsq=Lref.*Lref;

% here solving for S(2:nz+1) with bc S(1) and S(nz+2)
% equations(matrix rows)numbered 1:nz corresponding to S(2:nz+1)
% -- Generic diagonals and R.H.S. --
xi_n=0.5*(xi(2:nz+1)+xi(3:nz+2));
xi_s=0.5*(xi(1:nz)+xi(2:nz+1));
A_Diag   = - Lsq(2:nz+1).*cdz2.*(xi_n+xi_s)-1./k(2:nz+1);
A_Diag_m = [Lsq(2:nz)'.*cdz2.*xi_s(2:nz); 0];
A_Diag_p = [0; Lsq(1:nz-1)'.*cdz2.*xi_n(1:nz-1)];  % incorrect version was xi_n(2:nz)
                                         % changed direction 9/26/12
b=ones(nz,1);

%  -- apply boundary conditions (dirichlet)
%%bc.S0=k(2);
%%bc.Snz=k(nz+1);
% bc.Snz=0;
% bc.S0=0;
% for dSz/dz=0 at boundary
A_Diag(1)  = A_Diag(1)  + Lsq(1)*cdz2*xi_s(1);
b(1)       = b(1);
%  for Sz=S0 at boundary
% A_Diag(1)  = A_Diag(1)  - Lsq*cdz2*xi_s(1);
% b(1)       = b(1)       - Lsq*cdz2*xi_s(1) * 2*bc.S0;
% for Sz=Snz at boundary
A_Diag(nz) = A_Diag(nz) - Lsq(nz)*cdz2*xi_n(nz);
b(nz)      = b(nz)      - Lsq(nz)*cdz2*xi_n(nz) * 2*bc.Snz;
% for dSz/dz=0 at boundary
% A_Diag(nz) = A_Diag(nz) + Lsq*cdz2*xi_n(nz);
% b(nz)      = b(nz);

% ghost-cell population:
%boundS = @(v) [(2 * bc.S0 - v(2)) v(2:nz + 1) (2 * bc.Snz - v(nz+1))];
% boundS = @(v) [(2 * bc.S0 - v(2)) v(2:nz + 1) v(nz + 1)];
boundS = @(v) [v(2) v(2:nz + 1) (2 * bc.Snz - v(nz+1))];

% -- Assemble sparse operator --
A = spdiags([A_Diag_m A_Diag A_Diag_p],[-1 0 1],nz,nz);
%condest(A)

% -- Solution --
Sdir = [0 (inv(A) * b)' 0];
S    = boundS(Sdir);

% BiCG
% bicg_tol   = 1e-9;
% bicg_maxit = 2048;
% Sbicg = [0 bicg(A,b,bicg_tol,bicg_maxit,[],[],(inv(A) * b))' 0];
% %Sbicg = [0 bicg(A,b,bicg_tol,bicg_maxit)' 0];
% S = boundS(Sbicg);

end




