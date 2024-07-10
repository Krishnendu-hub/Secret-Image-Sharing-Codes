clear all;
close all;
im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\lena.bmp'));
%im=im2double(imread('cameraman.tif'));
%im=im2double(imread('rice.png'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\baboon512.bmp'));
%im=im2double(imread('text.png'));
[R,C]=size(im);
n=4;
BlockSize=32;
[dividedImage,TOTAL_BLOCKS]=divideIntoNonOverlapBlocks(im,BlockSize);
syms x;
k=0;
for j=1:n+2
    p(x)=(x^2-1)^k;
    for i=1:k
        p(x)=diff(p(x)); 
    end
    pk(j)=1/(2^k*factorial(k))*p(x);
    k=k+1;
end

for i=1:BlockSize*BlockSize
         %A(i)=(2*(i-1)-(BlockSize*BlockSize-1))/(BlockSize*BlockSize-1);
         A(i)=i/(BlockSize*BlockSize);
end
 for j=1:n+2
    for i=1:BlockSize*BlockSize
       PolynomialMatrix(i,j)=(subs(pk(j),A(i)));
    end
 end
COUNT=1;
SS1=zeros(R,C);
SS2=zeros(R,C);
SS3=zeros(R,C);
SS4=zeros(R,C);

  U_matrix=(3:n+2);
  V_matrix=perms(U_matrix);
  %I=inv(transpose(PolynomialMatrix)*PolynomialMatrix)*transpose(PolynomialMatrix);
  PolynomialMatrix_1=pinv(PolynomialMatrix);
  %PolynomialMatrix_1*PolynomialMatrix
  %I*PolynomialMatrix;
for p=1:(C/BlockSize)
    for q=1:(R/BlockSize)
        im1=dividedImage(:,:,COUNT);
        im2 = reshape(im1,BlockSize*BlockSize,1);
%         x=double(linsolve(PolynomialMatrix,im2));
%         PolynomialMatrix*x;
         x=PolynomialMatrix_1*im2;
%        double(PolynomialMatrix*PolynomialMatrix_1)
%        x= gausselimination(PolynomialMatrix,im2);
        ran=randi(factorial(n),1);
%        double(PolynomialMatrix(1,1)*x(1)+PolynomialMatrix(1,2)*x(2)+PolynomialMatrix(1,3)*x(3)+PolynomialMatrix(1,4)*x(4));
%         s1=PolynomialMatrix(:,V_matrix(ran,1))*x(V_matrix(ran,1))+(im2-PolynomialMatrix*x)/n;
%         s2=PolynomialMatrix(:,V_matrix(ran,2))*x(V_matrix(ran,2))+(im2-PolynomialMatrix*x)/n;
%         s3=PolynomialMatrix(:,V_matrix(ran,3))*x(V_matrix(ran,3))+(im2-PolynomialMatrix*x)/n;
%         s4=PolynomialMatrix(:,V_matrix(ran,4))*x(V_matrix(ran,4))+(im2-PolynomialMatrix*x)/n;
        if mod(COUNT,n)==1   
        s1=PolynomialMatrix(:,V_matrix(ran,1))*x(3)+ PolynomialMatrix(:,V_matrix(ran,1))*(x(V_matrix(ran,1))-x(3))  + (im2-PolynomialMatrix*x)/n +(x(1)*PolynomialMatrix(:,1));
        s2=PolynomialMatrix(:,V_matrix(ran,2))*x(4)+ PolynomialMatrix(:,V_matrix(ran,2))*(x(V_matrix(ran,2))-x(4))  + (im2-PolynomialMatrix*x)/n;
        s3=PolynomialMatrix(:,V_matrix(ran,3))*x(5)+ PolynomialMatrix(:,V_matrix(ran,3))*(x(V_matrix(ran,3))-x(5))  + (im2-PolynomialMatrix*x)/n;
        s4=PolynomialMatrix(:,V_matrix(ran,4))*x(6)+ PolynomialMatrix(:,V_matrix(ran,4))*(x(V_matrix(ran,4))-x(6))  + (im2-PolynomialMatrix*x)/n +(x(2)*PolynomialMatrix(:,2));
       
       
        elseif mode(COUNT,n)==2
        s1=PolynomialMatrix(:,V_matrix(ran,1))*x(3)+ PolynomialMatrix(:,V_matrix(ran,1))*(x(V_matrix(ran,1))-x(3))  + (im2-PolynomialMatrix*x)/n;
        s2=PolynomialMatrix(:,V_matrix(ran,2))*x(4)+ PolynomialMatrix(:,V_matrix(ran,2))*(x(V_matrix(ran,2))-x(4))  + (im2-PolynomialMatrix*x)/n + (x(1)*PolynomialMatrix(:,1));
        s3=PolynomialMatrix(:,V_matrix(ran,3))*x(5)+ PolynomialMatrix(:,V_matrix(ran,3))*(x(V_matrix(ran,3))-x(5))  + (im2-PolynomialMatrix*x)/n + (x(2)*PolynomialMatrix(:,2));
        s4=PolynomialMatrix(:,V_matrix(ran,4))*x(6)+ PolynomialMatrix(:,V_matrix(ran,4))*(x(V_matrix(ran,4))-x(6))  + (im2-PolynomialMatrix*x)/n;
        
  
        
        elseif mode(COUNT,n)==3
        s1=PolynomialMatrix(:,V_matrix(ran,1))*x(3)+ PolynomialMatrix(:,V_matrix(ran,1))*(x(V_matrix(ran,1))-x(3)) + (im2-PolynomialMatrix*x)/n;
        s2=PolynomialMatrix(:,V_matrix(ran,2))*x(4)+ PolynomialMatrix(:,V_matrix(ran,2))*(x(V_matrix(ran,2))-x(4)) + (im2-PolynomialMatrix*x)/n + (x(2)*PolynomialMatrix(:,2));
        s3=PolynomialMatrix(:,V_matrix(ran,3))*x(5)+ PolynomialMatrix(:,V_matrix(ran,3))*(x(V_matrix(ran,3))-x(5)) + (im2-PolynomialMatrix*x)/n + (x(1)*PolynomialMatrix(:,1));
        s4=PolynomialMatrix(:,V_matrix(ran,4))*x(6)+ PolynomialMatrix(:,V_matrix(ran,4))*(x(V_matrix(ran,4))-x(6)) + (im2-PolynomialMatrix*x)/n;
        
     
        
        else %mode(COUNT,n)==0
        s1=PolynomialMatrix(:,V_matrix(ran,1))*x(3)+ PolynomialMatrix(:,V_matrix(ran,1))*(x(V_matrix(ran,1))-x(3))  + (im2-PolynomialMatrix*x)/n + (x(2)*PolynomialMatrix(:,2));
        s2=PolynomialMatrix(:,V_matrix(ran,2))*x(4)+ PolynomialMatrix(:,V_matrix(ran,2))*(x(V_matrix(ran,2))-x(4))  + (im2-PolynomialMatrix*x)/n;
        s3=PolynomialMatrix(:,V_matrix(ran,3))*x(5)+ PolynomialMatrix(:,V_matrix(ran,3))*(x(V_matrix(ran,3))-x(5))  + (im2-PolynomialMatrix*x)/n;
        s4=PolynomialMatrix(:,V_matrix(ran,4))*x(6)+ PolynomialMatrix(:,V_matrix(ran,4))*(x(V_matrix(ran,4))-x(6))  + (im2-PolynomialMatrix*x)/n + (x(1)*PolynomialMatrix(:,1));
     
        end
        
        ss1 =reshape(s1, BlockSize, BlockSize);
        ss2 =reshape(s2, BlockSize, BlockSize);
        ss3 =reshape(s3, BlockSize, BlockSize);
        ss4 =reshape(s4, BlockSize, BlockSize);        
        
        startrow =BlockSize*(p-1)+1;     
        startcol =BlockSize*(q-1)+1;
        endrow = startrow + BlockSize -1;
        endcol = startcol + BlockSize-1;
        SS1(startrow : endrow, startcol : endcol)=ss1;
        SS2(startrow : endrow, startcol : endcol)=ss2;
        SS3(startrow : endrow, startcol : endcol)=ss3;
        SS4(startrow : endrow, startcol : endcol)=ss4;
        COUNT=COUNT+1;
    end
end
PolynomialMatrix;
SS=SS1+SS2+SS3+SS4;%+SS5;%+SS6+SS7+SS8+SS9+SS10+SS11+SS12+SS13+SS14+SS15+SS16;
figure(1)
imshow(mat2gray(SS1));
figure(2)
imshow(mat2gray(SS2));
figure(3)
imshow(mat2gray(SS3));
figure(4)
imshow(mat2gray(SS4));

figure(15)
imshow(mat2gray(SS1+SS2));
figure(16)
imshow(mat2gray(SS1+SS2+SS3));
figure(17)
imshow(mat2gray(SS));