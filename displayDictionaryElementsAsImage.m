function I = displayDictionaryElementsAsImage(D, numRows, numCols,X,Y,sortVarFlag)
% function I = displayDictionaryElementsAsImage(D, numRows, numCols, X,Y)
% displays the dictionary atoms as blocks. For activation, the dictionary D
% should be given, as also the number of rows (numRows) and columns
% (numCols) for the atoms to be displayed. X and Y are the dimensions of
% each atom.


borderSize = 1;
columnScanFlag = 1;
strechEachVecFlag = 1;
showImFlag = 1;

if (length(who('X'))==0)
    X = 8;
    Y = 8;
end
if (length(who('sortVarFlag'))==0)
    a=1%%%%没有显示，说明没进入此循环
    sortVarFlag = 1;
end

numElems = size(D,2);
if (length(who('numRows'))==0)
    numRows = floor(sqrt(numElems));
    numCols = numRows;
end
if (length(who('strechEachVecFlag'))==0) 
    strechEachVecFlag = 0;
end
if (length(who('showImFlag'))==0) 
    showImFlag = 1;
end

%%% sort the elements, if necessary.
%%% construct the image to display (I)
sizeForEachImage = sqrt(size(D,1))+borderSize;
I = zeros(sizeForEachImage*numRows+borderSize,sizeForEachImage*numCols+borderSize,3);
%%% fill all this image in blue
I(:,:,1) = 0;%min(min(D));
I(:,:,2) = 0; %min(min(D));
I(:,:,3) = 1; %max(max(D));

%%% now fill the image squares with the elements (in row scan or column
%%% scan).
if (strechEachVecFlag)
    for counter = 1:size(D,2)
        D(:,counter) = D(:,counter)-min(D(:,counter));
        if (max(D(:,counter)))
            D(:,counter) = D(:,counter)./max(D(:,counter));%%%归一化
        end
    end
end
% 要注意的是var函数所采用公式中，分母不是 ，而是 。这是因为var函数实际上求的并不是方差，而是误差理论中“有限次测量数据的标准偏差的估计值”。
% X = [4 -2 1; 9 5 7]   sum((X(:,2)-mean(X(:,2))).^2)/(length(X(:,2))-1)等价于：var(X,0,1)   
% 注：var(X,W)  % W可以取0或1，取0求样本方差的无偏估计值（除以N-1；对应取1求得的是方差（除以N）， W也可以是向量，但必须与X中的第一个维度数相同，即length(W)= size(X,1)          
%   所以还存在：  var(X ,0 ,dim) % 除以N     dim =1 对每列操作      dim = 2 对每行操作
%                 var(X ,1 ,dim) % 除以N-1   dim =1 对每列操作    dim = 2 对每行操作
if (sortVarFlag)%%%%%%%%%%%%%%%%%按照方差从小打大排序
    vars = var(D);
    [V,indices] = sort(vars');%%从小到大的顺序将字典排序
    indices = fliplr(indices);%%a=[1 2 3;4 5 6;7 8 9]  a=fliplr(a)  左右翻转  如果对于只有一列的向量，翻转后是不变的
    D = [D(:,1:sortVarFlag-1),D(:,indices+sortVarFlag-1)];
    signs = sign(D(1,:));
    signs(find(signs==0)) = 1;
    D = D.*repmat(signs,size(D,1),1);%%%repmat(signs,size(D,1),1)  此矩阵的大小为64*256     D也为64*256
    D = D(:,1:numRows*numCols);
end

counter=1;
for j = 1:numRows
    for i = 1:numCols
%         if (strechEachVecFlag)
%             D(:,counter) = D(:,counter)-min(D(:,counter));
%             D(:,counter) = D(:,counter)./max(D(:,counter));
%         end
%         if (columnScanFlag==1)
%             I(borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,1)=reshape(D(:,counter),8,8);
%             I(borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,2)=reshape(D(:,counter),8,8);
%             I(borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,3)=reshape(D(:,counter),8,8);
%         else
            % Go in Column Scan:
            I(borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,1)=reshape(D(:,counter),X,Y);%%%8*8
            I(borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,2)=reshape(D(:,counter),X,Y);
            I(borderSize+(j-1)*sizeForEachImage+1:j*sizeForEachImage,borderSize+(i-1)*sizeForEachImage+1:i*sizeForEachImage,3)=reshape(D(:,counter),X,Y);
%         end
        counter = counter+1;
    end
end

if (showImFlag) 
    I = I-min(min(min(I)));
    I = I./max(max(max(I)));
    imshow(I,[]);
end
