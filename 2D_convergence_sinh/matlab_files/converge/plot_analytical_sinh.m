%plot 
close all;
x=[0:0.01:1];
y=[1:-0.01:0];

N=length(x);

p=zeros(N,N);
u=zeros(N,N);
v=zeros(N,N);

H=zeros(N,N);


for i=1:N
   for j=1:N
       
       X=x(j);
       Y=y(i);
       
      p(i,j)=-sin(X)*sinh(Y)-(cos(1)-1)*(cosh(1)-1);
      u(i,j)=cos(X)*sinh(Y);
      v(i,j)=sin(X)*cosh(Y);
      H(i,j)=X*Y;
   end
    
    
end



figure;
imagesc(p);

figure;
imagesc(u);

figure;
imagesc(v);