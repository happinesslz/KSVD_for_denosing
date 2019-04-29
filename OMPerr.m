function [A]=OMPerr(D,X,errorGoal); 
%=============================================
% Sparse coding of a group of signals based on a given 
% dictionary and specified number of atoms to use. 
% input arguments: D - the dictionary
%                  X - the signals to represent
%                  errorGoal - the maximal allowed representation error for
%                  each siganl.
% output arguments: A - sparse coefficient matrix.
%=============================================
[n,P]=size(X);%n=64  P= 62001=249*249
[n,K]=size(D);%n=64 K=256
E2 = errorGoal^2*n;
maxNumCoef = n/2;%%%%%%32
A = sparse(size(D,2),size(X,2));%参考稀疏矩阵的帮助256*10000
for k=1:1:P,
    a=[];
    x=X(:,k);
    residual=x;
	indx = [];
	a = [];
	currResNorm2 = sum(residual.^2);
	j = 0;

    while currResNorm2>E2 & j < maxNumCoef,
		j = j+1;
        proj=D'*residual;%参考pinv函数的帮助 256*1
        pos=find(abs(proj)==max(abs(proj)));%看看D（256列）中哪一列的值最大
        pos=pos(1);
        indx(j)=pos;%%%index的值为1到256
        %c++的opm优化速度的算法     http://blog.csdn.net/pi9nc/article/details/26593003
        a=pinv(D(:,indx(1:j)))*x;%j*64  *64*1=j*1    
        residual=x-D(:,indx(1:j))*a;
		currResNorm2 = sum(residual.^2);
   end;
   if (length(indx)>0)
       A(indx,k)=a;%%%a是j*1的矩阵,其中j=maxNumCoef
   end
end;
return;
