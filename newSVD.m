clear all;
close all;
clc;
%im=im2double(rgb2gray(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\news_paper512.jpg')));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\TextMask256.png'));
im=(imread('cameraman.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\TextMask256.png'));
%im=im2double(imread('text.png'));
%im=(imread('rice.png'));
%im=im2double(imread('testpat1.png'));
%im=im2double(imread('circles.png'));
%im=im2double(imread('liftingbody.png'));
%im=(imread('mri.tif')); 
%im=(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\male.tiff'));
%im=(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\lena.bmp'));
%im=im2double(imread('Mandrill.bmp'));
%im=im2double(rgb2gray(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\lena1.bmp')));
%im=im2double(rgb2gray(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Airplane1.tiff')));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\baboon512.bmp'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\ruler.tiff'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Clock.tiff'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Couple.tiff'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\boat.512.tiff'));
%im=im2double((imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Corn.tif')));
%im=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Airplane.tif')
figure(1)
imshow(mat2gray(im));
[R,C]=size(im);
m=round(mean(mean(im)));
m1=sqrt(var(double(reshape(im,R*C,1))));

for i=1:R
    for j=1:C
        i1=mod((i+j),R)+1;
        j1=mod((m*i+(m+1)*j),R)+1;
        im1(i1,j1)=im(i,j);
    end
end
%im1=im;
im1=single(im1);
figure(2)
imshow(mat2gray(im1));
numberOfShares=4;
blockSize=4;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ;
%[R,C]=size(im);
rows=R/2^(log2(blockSize));% define how many rows of block
cols=C/2^(log2(blockSize));% define how many colls of block 
sizeI = size(im1);
blocks =mat2tiles(im1, ceil(sizeI(1:2)./[rows cols]));%[rows cols]
[a,b]=size(blocks);
U_matrix=(1:numberOfShares);
V_matrix=perms(U_matrix);
[p,q]=size(V_matrix);
p1=1;
index=0;
ss1=zeros(R,C);
ss2=zeros(R,C); 
ss3=zeros(R,C);
ss4=zeros(R,C);
for i = 1:rows/2
    for j=1:cols/2
    k=i+rows/2;
    l=j+cols/2;
    [u1,s1,v1]=svd(blocks{j,i});
    [u2,s2,v2]=svd(blocks{l,k});
    v1=v1';
    v2=v2';
    r = randi(factorial(numberOfShares),1);

    s1_1=s2(1,1)*u2(:,V_matrix(r,1))*v2(V_matrix(r,2),:)+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,2))*v1(V_matrix(r,2),:)-s2(1,1)*u2(:,V_matrix(r,2))*v2(V_matrix(r,1),:);%+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)-s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)
    s2_1=s2(1,1)*u2(:,V_matrix(r,2))*v2(V_matrix(r,1),:)+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,1))*v1(V_matrix(r,1),:)-s2(1,1)*u2(:,V_matrix(r,1))*v2(V_matrix(r,2),:);%+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)-s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)
    s3_1=s2(1,1)*u2(:,V_matrix(r,3))*v2(V_matrix(r,2),:)+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,3))*v1(V_matrix(r,3),:)-s2(1,1)*u2(:,V_matrix(r,3))*v2(V_matrix(r,4),:);%+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)-s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)
    s4_1=s2(1,1)*u2(:,V_matrix(r,3))*v2(V_matrix(r,4),:)+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,4))*v1(V_matrix(r,4),:)-s2(1,1)*u2(:,V_matrix(r,3))*v2(V_matrix(r,2),:);%+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)-s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)
    
    r = randi(factorial(numberOfShares),1);
    s1_2=s1(1,1)*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,2))*v2(V_matrix(r,2),:)-s1(1,1)*u1(:,V_matrix(r,2))*v1(V_matrix(r,1),:);%+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)-s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)
    s2_2=s1(1,1)*u1(:,V_matrix(r,2))*v1(V_matrix(r,1),:)+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,1))*v2(V_matrix(r,1),:)-s1(1,1)*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:);%+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)-s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)
    s3_2=s1(1,1)*u1(:,V_matrix(r,3))*v1(V_matrix(r,4),:)+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,3))*v2(V_matrix(r,3),:)-s1(1,1)*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:);%+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)-s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)
    s4_2=s1(1,1)*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,4))*v2(V_matrix(r,4),:)-s1(1,1)*u1(:,V_matrix(r,3))*v1(V_matrix(r,4),:);%+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)-s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)
    
    rand=randi(factorial(numberOfShares),1);   
     
        if rand==1
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1; 
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==2
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_1;
        ss2(startrow : endrow, startcol : endcol) = s3_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
         
        elseif rand==3
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;

        elseif rand==4
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==5
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==6
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s3_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==7
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==8
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s3_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==9
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==10
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==11
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==12
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==13
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==14
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        
        elseif rand==15
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==16
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==17
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        
        elseif rand==18
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==19
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==20
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        
        elseif rand==21
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==22
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==23
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s3_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        
        elseif rand==24
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        end
         
        rand=randi(factorial(numberOfShares),1);
        
        if rand==1
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==2
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
          
        elseif rand==3
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;

        elseif rand==4
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==5
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==6
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==7
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==8
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==9
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
        
        elseif rand==10
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==11
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==12
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
        
        elseif rand==13
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==14
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==15
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==16
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==17
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==18
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==19
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
        
        elseif rand==20
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==21
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
        
        elseif rand==22
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==23
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==24
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
         end
    p1=p1+1;
    end 
end 
q1=p1;
for i =rows/2+1:rows
    for j=1:cols/2
    k=i-rows/2;
    l=j+cols/2;
    [u1,s1,v1]=svd(blocks{j,i});
    [u2,s2,v2]=svd(blocks{l,k});
    v1=v1';
    v2=v2';

    r = randi(factorial(numberOfShares),1); 
    s1_1=s2(1,1)*u2(:,V_matrix(r,1))*v1(V_matrix(r,4),:)+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,2))*v1(V_matrix(r,2),:)-s2(1,1)*u2(:,V_matrix(r,2))*v1(V_matrix(r,3),:);
    s2_1=s2(1,1)*u2(:,V_matrix(r,2))*v1(V_matrix(r,3),:)+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,1))*v1(V_matrix(r,1),:)-s2(1,1)*u2(:,V_matrix(r,1))*v1(V_matrix(r,4),:);
    s3_1=s2(1,1)*u2(:,V_matrix(r,3))*v1(V_matrix(r,2),:)+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,3))*v1(V_matrix(r,3),:)-s2(1,1)*u2(:,V_matrix(r,4))*v1(V_matrix(r,1),:);
    s4_1=s2(1,1)*u2(:,V_matrix(r,4))*v1(V_matrix(r,1),:)+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,4))*v1(V_matrix(r,4),:)-s2(1,1)*u2(:,V_matrix(r,3))*v1(V_matrix(r,2),:);
    
    r = randi(factorial(numberOfShares),1);
    
    s1_2=s1(1,1)*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,2))*v2(V_matrix(r,2),:)-s1(1,1)*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:);
    s2_2=s1(1,1)*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,1))*v2(V_matrix(r,1),:)-s1(1,1)*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:);
    s3_2=s1(1,1)*u2(:,V_matrix(r,3))*v2(V_matrix(r,2),:)+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,3))*v2(V_matrix(r,3),:)-s1(1,1)*u2(:,V_matrix(r,4))*v2(V_matrix(r,1),:);
    s4_2=s1(1,1)*u2(:,V_matrix(r,4))*v2(V_matrix(r,1),:)+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,4))*v2(V_matrix(r,4),:)-s1(1,1)*u2(:,V_matrix(r,3))*v2(V_matrix(r,2),:);
    
    rand=randi(factorial(numberOfShares),1);  
   
        if rand==1
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==2
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
         
        elseif rand==3
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;

        elseif rand==4
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s3_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==5
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==6
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==7
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s3_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==8
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s3_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==9
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==10
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==11
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==12
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==13
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==14
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        
        elseif rand==15
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s4_1;
        
        elseif rand==16
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s4_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==17
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        
        elseif rand==18
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==19
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s2_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==20
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s1_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        
        elseif rand==21
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s3_1;
        
        elseif rand==22
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s2_1;
        ss3(startrow : endrow, startcol : endcol) = s3_1;
        ss4(startrow : endrow, startcol : endcol) = s1_1;
        
        elseif rand==23
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_1;
        ss2(startrow : endrow, startcol : endcol) = s3_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
        
        elseif rand==24
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_1;
        ss2(startrow : endrow, startcol : endcol) = s4_1;
        ss3(startrow : endrow, startcol : endcol) = s1_1;
        ss4(startrow : endrow, startcol : endcol) = s2_1;
         end
        
        rand=randi(factorial(numberOfShares),1);
        if rand==1
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==2
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s1_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
         
        elseif rand==3
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;

        elseif rand==4
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==5
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==6
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==7
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==8
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==9
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
        
        elseif rand==10
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==11
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==12
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s2_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
        
        elseif rand==13
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==14
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==15
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s4_2;
        
        elseif rand==16
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s4_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==17
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==18
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==19
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s2_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
        
        elseif rand==20
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s1_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==21
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s3_2;
        
        elseif rand==22
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s2_2;
        ss3(startrow : endrow, startcol : endcol) = s3_2;
        ss4(startrow : endrow, startcol : endcol) = s1_2;
        
        elseif rand==23
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s4_2;
        ss2(startrow : endrow, startcol : endcol) = s3_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        
        elseif rand==24
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1(startrow : endrow, startcol : endcol) = s3_2;
        ss2(startrow : endrow, startcol : endcol) = s4_2;
        ss3(startrow : endrow, startcol : endcol) = s1_2;
        ss4(startrow : endrow, startcol : endcol) = s2_2;
        end
    q1=q1+1;
    end                
end
% ss1=uint8(ss1);
% ss2=uint8(ss2);
% ss3=uint8(ss3);
% ss4=uint8(ss4);
mean1=mean(reshape(ss1,R*C,1));
mean2=mean(reshape(ss1,R*C,1));
mean3=mean(reshape(ss1,R*C,1));
mean4=mean(reshape(ss1,R*C,1));
%ss1(21,5)=100000;
tic
%n2=single(ss1+ss2+ss3)%+ss4
%n1=ceil(mean(mean((uint8(ss1+ss2+ss3)))));
% ceil(mean(mean((ss1+ss2+ss3)*100)))
im2=single(ss1+ss2+ss3+ss4);

toc
%im123=single(ss1+ss2+ss3);
%im2=single(ss1+ss2+ss3);
%n=ceil(mean(uint8(reshape(im2,R*C,1))));
%n123=ceil(mean(uint8(reshape(im123,R*C,1))));
n=round(mean(reshape(im2,R*C,1)));
%n=126;
n1=sqrt(var(double(reshape(im2,R*C,1))));
cor12=corr2(ss1,ss2);
%im4=uint8(ss1+ss2+ss3);
for i=1:R
    for j=1:C
        i2=mod(((n+1)*(i-1)-(j-1)),R)+1;
        j2=mod((-n*(i-1)+(j-1)),R)+1; 
        im3(i2,j2)=im2(i,j);

%         i22=mod(((n1+1)*(i-1)-(j-1)),R)+1;
%         j22=mod((-n1*(i-1)+(j-1)),R)+1; 
%         im5(i22,j22)=im4(i,j);
    end
end
%reconstruction from addition of share (123)
% for i=1:R
%     for j=1:C
%         i2=mod(((n123+1)*(i-1)-(j-1)),R)+1;
%         j2=mod((-n123*(i-1)+(j-1)),R)+1; 
%         im3_123(i2,j2)=im123(i,j);
% 
% %         i22=mod(((n1+1)*(i-1)-(j-1)),R)+1;
% %         j22=mod((-n1*(i-1)+(j-1)),R)+1; 
% %         im5(i22,j22)=im4(i,j);
%     end
% end

% for i=1:R
%     for j=1:C
%         i2=mod(((n123+1)*(i-1)-(j-1)),R)+1;
%         j2=mod((-n123*(i-1)+(j-1)),R)+1; 
%         im3_123(i2,j2)=im123(i,j);
% 
% %         i22=mod(((n1+1)*(i-1)-(j-1)),R)+1;
% %         j22=mod((-n1*(i-1)+(j-1)),R)+1; 
% %         im5(i22,j22)=im4(i,j);
%     end
% end
%n1=sqrt(var(double(reshape(im3,R*C,1))));

% j=entropy(im);
figure(3);
imshow(mat2gray(ss1));
%imshow(imhist(ss1))
imwrite(ss1,'share1.png');
%j1=entropy(ss1)
figure(4);
imshow(mat2gray(ss2));
imwrite(ss2,'share2.png');
%j2=entropy(ss2)
figure(5);
imshow(mat2gray(ss3));
imwrite(ss3,'share3.png');
%j3=entropy(ss3)
figure(6);
imshow(mat2gray(ss4));
imwrite(ss4,'share4.png');
%j4=entropy(ss4)
% figure(7);
% imshow(mat2gray(ss1+ss2+ss3+ss4));
figure(8);
imshow(mat2gray(im3));
% 
% % figure(8);
% % imshow(mat2gray(ss1+ss4));
% % % % j7=entropy(ss1+ss4)
% % figure(9);
% % imshow(mat2gray(ss2+ss3));
% % % % j8=entropy(ss2+ss3)
% % figure(10);
% % imshow(mat2gray(ss2+ss4));
% % % % j9=entropy(ss2+ss4)
% % figure(11);
% % imshow(mat2gray(im5));
% % % % j10=entropy(ss3+ss4)
% figure(12);
% imshow(mat2gray(im3_123));
% entropy(uint8(im3_123))
% figure(13);
% imshow(mat2gray(im3));
% figure(14)
% histogram(double(ss1));
% figure(15)
% histogram(uint8(ss2));
% figure(16)
% histogram(uint8(ss3));
% figure(17)
% histogram(uint8(ss4));
%entropy(uint8(im3))
% j123=entropy(ss1+ss2+ss3); 
% figure(13);
% imshow(mat2gray(ss1+ss2+ss4));
% %j7=entropy(ss1+ss2+ss4);
% figure(14);
% imshow(mat2gray(ss2+ss3+ss4));
% % % j13=entropy(ss2+ss3+ss4)
% im=(im);
% ss=(ss1+ss2+ss3+ss4);
% ss_1=(ss1+ss2+ss3);
% ss_2=(ss1+ss2);
% Min_square_error=immse(im,ss);
% % Min_square_error_1=immse(im,ss_1);
% % Min_square_error_2=immse(im,ss_2);%uint16
% PSNR_Value=psnr(im,ss);
% % PSNR_Value_1=psnr(im,ss_1);
% % PSNR_Value_2=psnr(im,ss_2);
% figure(15);
% imshow(mat2gray((ss)));%ss1+ss2+ss3+ss4
% j1234=entropy(ss1+ss2+ss3+ss4);
% % ss=(ss1+ss2+ss3+ss4);
% % im=(im); 
% 
% % OriginalImage_contrast=max(im(:)) - min(im(:));
% % Recontructed_image_contrast = max(ss(:)) - min(ss(:));
% % 
% % % Recontructed_image_contrast_1 = max(ss_1(:)) - min(ss_1(:));
% % % Recontructed_image_contrast_2 = max(ss_2(:)) - min(ss_2(:));
%  
%Stactural_similarity=ssim(im,ss);
% Stactural_similarity_1=ssim(ss1,single(im));
% Stactural_similarity_2=ssim(ss2,single(im));
% Stactural_similarity_3=ssim(ss3,single(im));
% Stactural_similarity_4=ssim(ss4,single(im));
% % 
%  Correlation_coefficint=corr2(im,ss);
% Correlation_coefficint_1=corr2(im,ss_1);
% Correlation_coefficint_2=corr2(im,ss_2);
% 
% j_1234=entropy(ss1+ss2+ss3+ss4);
% j_123=entropy(ss1+ss2+ss3);
% j_12=entropy(ss1+ss2);
%  k=min(min(im));
%  k1=max(max(im));
%  con=k1-k
%  contrast_secret=(k1-k)/(k1+k)
% 
%  k=min(min(ss));
%  k1=max(max(ss));
%  contrast_reconstructed=(k1-k)/(k1+k)
%im_diff=im2bw(im)-im2bw(ss);
%ss=im2bw(ss);
%im=im2bw(im);
% [hash_original,A]=hash1(im);
% [hash_im,A]=hash1(im);
% [hash_ss,A]=hash1(ss);
% [hash_S1,A1]=hash1(ss1);
% [hash_S2,A2]=hash1(ss2);
% [hash_S3,A3]=hash1(ss3);
% [hash_S4,A4]=hash1(ss4);
% [hash_S12,A1]=hash2(hash_S1,hash_S2);
% [hash_S34,A2]=hash2(hash_S3,hash_S4);
% [hash_S1234,A2]=hash2(hash_S12,hash_S34);

