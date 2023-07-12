function A= call_mat(M)
 [n,m]=size(M);
A=zeros(n*m,1);
for i =1:n
   for j=1:m 
    A(m*(i-1)+j,1)=M(i,j);
   end
end
end