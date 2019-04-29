%%%%此为伪代码不能运行，是别人论文里的

DCT=zeros(8,16);
Pn=10;
bb=10;
for k=0:1:Pn-1
    V=cos([0:1:bb-1]'*k*pi/16);
    if k>0
        V=V-mean(V);
    end
    DCT(:,k+1)=V/norm(V);
end
DCT=kron(DCT,DCT);