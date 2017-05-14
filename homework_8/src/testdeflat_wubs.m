n=100; 
m=10; % m should be a divisor of n
P=kron(eye(m),ones(n/m,1)); %space used for Galerkin projection of B
% the vectors are rather smooth, since we are only interested in 
%finding an approximation to the eigenvector corresponding to the
%smallest eigenvector.
A=linbuck(n); 
B=-(A+10*speye(n)); % matrix shifted such that one eigenvalue is negative
%Find approximation to the eigenvector corresponding to smallest 
%eigenvalue. Note that Bsmall has only positive eigenvalues, while B has
%one negative eigenvalue.
Bsmall=P'*B*P; [V,D]=eig(Bsmall); v=P*V(:,1); v=v/norm(v);
proj1=speye(n)-v*v';
Bproj1=proj1*B*proj1; %Galerkin approximation of 
%B on space perpendicular to v
eig(Bproj1) %note that all eigenvalues are positive now
proj2=eye(n)-proj1*(Bproj1\(proj1*B)); 
norm(proj2'*B*proj1) % should be zero by construction
Bproj2=proj2'*B*proj2; % Another Galerkin projection of B
% In fact we will use the following projection/rewrite of the equations
%  [proj1, proj2]'*(B*[proj1, proj2][xproj1,xproj2]'-rhs)=0
eig(Bproj2) %note that there is only one negative eigenvalue
opts.type='ilutp';
opts.droptol=1.0;
[L,U]=ilu(B,opts); %makes incomplete LU, in fact a diagonal matrix
norm(B-L*U,1) % shows how far we are off
B-L*U % """
eig(full(L*U),full(B)) %shows the eigenvalues of preconditioned matrix
rhs=sin(0.01*(1:n)'); % An arbitrary smooth rhs.
rhsproj1=proj1*rhs; %   
rhsproj2=proj2'*rhs; 
%straightforward solution of system
[x,FLAG,RELRES,ITER,resvec]=gmres(B,rhs,20,1e-12,20,L,U);
resvec' % it fails
norm(B*x-rhs) % check result
RELRES,ITER
[xproj1,FLAG,RELRES,ITER,resvec]=gmres(Bproj1,rhsproj1,20,1e-12,20, L,U);
%[xproj1,FLAG,RELRES,ITER,resvec]=pcg(Bproj1,rhsproj1,1e-12,400, L,U);
xproj1=proj1*xproj1; 
resvec' %it converges
RELRES,ITER
%[xproj2,FLAG,RELRES,ITER,resvec]=gmres(Bproj2,rhsproj2,20,1e-12,20);
[xproj2,FLAG,RELRES,ITER,resvec]=pcg(-Bproj2,-rhsproj2,1e-12,400);
xproj2=proj2*xproj2;
resvec' % converges in one step.
RELRES,ITER
x=xproj1+xproj2; %construction of final solution. 
norm(B*x-rhs) % check of result
