function aa=cleanup(a,tol)
%function cleanup(a,tol) removes small terms 
%and replace them with a perfect zero.

if nargin~=2
	tol=eps;
end
[n,m]=size(a);
aa=a;
for i=1:n
	for j=1:m
		if(abs(a(i,j)) <= tol)
			aa(i,j)=0;
		else
		end
	end
end
