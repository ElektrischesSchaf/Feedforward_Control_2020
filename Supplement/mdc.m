% MDC   Decouple a real matrix into stable, unstable and pure nonhyperbolic
%       submatrices.
%       [As, Au, Anh, A_dec, T]= mdc(A, CRITERIA, te) decouples the input matrix
%       A into the stable submatrix As, unstable submatrix Au, and the 
%       nonhyperbolic submatrix Anh. The returned mapping matrix T is
%       given as:
%
%             /  As   0   0 \%      A =  T |  0  Anh  0  | T^{-1}   %             \  0  0   Au  /
%      
%     The returned matrix is empty if the corresponding portion does not 
%     exsit in the input matrix A. The string CRITERIA selects the method
%     to do the decoupling among the following:
%        'c'       Decouple A in the continuouse time sense, i.e. according
%                  to the imaginary axis.
%        'd'       Decouple A in the discrete time sense, i.e. according to
%                  the unit circle.
%     te is the truncation error for the eigenvalues of A, i.e. if the 
%     imaginary portion of the eigenvalue is smaller than te, it is 
%     regarded as a real eigenvalue, the default value of 10^{-12} is used
%     if this parameter is missed.
%     Note A must be real and have distinct eigenvectors.
%     Version : 0.1
%     Date    : July, 2001.
function [As, Au, Anh, A_dec, T]= mdc(A, CRITERIA, te);%......................Input Paramters.....................................%..........................................................................% A        : The matrix to be decoupled.
% CRITERIA : The criteria to do the decoupling:
%    'c'   : Do the decoupling in the continuous time sense, i.e. according
%            to the imaginary axis.
%    'd'   : Do the decoupling in the discrete time sense, i.e.according to 
%            the unit circle.
%......................Output Paramters.................................... 
%..........................................................................
% As   : The stable submatrix.
% Au   : The unstable subsystem.
% T    : The transform matrix from undecoupled coordinate to the 
%          decoupled coordinate.
%...........................................................................
% nargin: Variable provided by MATLAB (number of arguments of input).
if nargin == 2
   te = 1e-12;
end% Get the diagonal matrix with eigenvalues of the input matrix and the
% corresponding transformation matrix.
[Tq, Jq]  = jordan(A);
Index_us    = 0;
Index_s     = 0;
Index_nh    = 0;
dem         = length(A);
% Obtain the index of the stable and unstable eigenvalues and their
% eigenvectors.
for jj=1:dem
 switch(CRITERIA)
 case 'c' 
    Eig_check = real(Jq(jj,jj));
 case 'd'
    Eig_check = abs(Jq(jj,jj))-1;
 end
       % The nonhyperbolic portion.
 if(abs(Eig_check)<=te)
       Index_nh = Index_nh+1;
       % Eig_nh: The nonhyperbolic matrix.
       Eig_nh(Index_nh)=Jq(jj,jj);
       NH_Index(Index_nh)=jj;
 else
      if(Eig_check>te)
      % The unstable portion.
       Index_us=Index_us+1;
       % Eig_us: The unstable matrix.
       Eig_us(Index_us)=Jq(jj,jj);  
       U_Index(Index_us)=jj;
     else
      Index_s=Index_s+1;
      S_Index(Index_s)=jj;
      Eig_s(Index_s)=Jq(jj,jj);
     end
  end
end
% 1. Rearrange the transform matrix, so that the upper left submatrix of 
%    the decoupled matrix will be the stable portion, the middle submatrix
%    will be the nonhyperbolic portion, and the lower right 
%    submatrix will be the unstable portion.
% 2. If the eigenvalue is a complex number, replace its eigenvector with
%    the real part of the eigenvector, and replace the eigenvector of its 
%    conjugate eigenvalue with the imaginary part of that eigenvector.
     jj=1;
     while(jj<=Index_s)
     if(abs(imag(Eig_s(jj)))>te)
         T(:,jj)=real(Tq(:,S_Index(jj)));
         T(:,jj+1) = imag(Tq(:, S_Index(jj)));
	 jj = jj+2;
      else
         T(:,jj)=real(Tq(:,S_Index(jj)));
         jj = jj+1;      end
     end
     kk = Index_s+1;
     jj = 1;
     while(jj<=Index_nh)
     if(abs(imag(Eig_nh(jj)))>te)
         T(:,kk)   = real(Tq(:,NH_Index(jj)));
         T(:,kk+1) = imag(Tq(:, NH_Index(jj)));
	 jj = jj+2;
         kk = kk+2;
      else
         T(:,jj+Index_s) = real(Tq(:,NH_Index(jj)));
         jj = jj+1;
         kk = kk+1;
      end
    end
   kk = Index_s+Index_nh+1;
   jj = 1;   while(jj<=Index_us)
     if(abs(imag(Eig_us(jj)))>te)
         T(:,kk)=real(Tq(:,U_Index(jj)));
         T(:,kk+1) = imag(Tq(:, U_Index(jj+1)));
		jj = jj+2;
         kk = kk+2;
      else
         T(:,kk)=real(Tq(:,U_Index(jj)));
         jj = jj+1;
         kk = kk+1;
      end
    end
% Transform the matrix.
A_dec = inv(T)*A*T;
As    = A_dec(1:Index_s, 1:Index_s);
Anh   = A_dec(Index_s+1:Index_s+Index_nh, Index_s+1:Index_s+Index_nh);
Au    = A_dec(dem-Index_us+1:dem, dem-Index_us+1:dem);







