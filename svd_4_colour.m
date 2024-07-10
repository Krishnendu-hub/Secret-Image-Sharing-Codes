clear all;
close all;
clc;
%im=((imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Barbara256rgb.tif')));
%im=(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\airplane.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Lake.tif'));
%im=(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Zebra.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Haight.tif'));
%im=(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\lena.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\pepper.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Flower.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\starfish.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Miss.tif'));
%im=(imread('lena.bmp'));
%im=(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\House.tif'));
im=im2double(imread('D:\Research\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Mickey.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Tower.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Nanna.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Bear.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Mural.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Penester.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Buddhist.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\Cowboy.tif'));
figure(1)
imshow(mat2gray(im));


R    = single(im(:, :, 1));
G    = single(im(:, :, 2));
B    = single(im(:, :, 3));

[ROW,COL,CHANNEL]=size(im);
m1=(mean(mean(R)));
m2=(mean(mean(G)));
m3=(mean(mean(B)));
m=ceil((m1+m2+m3)/3);
%m1=sqrt(var(double(reshape(im,R*C,1))));

for i=1:ROW
    for j=1:COL
        i1=mod((i+j),ROW)+1;
        j1=mod((m*i+(m+1)*j),ROW)+1;
        im_R(i1,j1)=R(i,j);
        im_G(i1,j1)=G(i,j);
        im_B(i1,j1)=B(i,j);
    end
end
im1=cat(3,im_R,im_G,im_B);
figure(1)
imshow(mat2gray(im1));
%entropy_lena_color=entropy(im);
% figure(8)
% imshow(mat2gray(R));
% figure(9)
% imshow(mat2gray(G));
% figure(10)
% imshow(mat2gray(B));
% k=min(min(R));
% k1=max(max(R));
% contrast_secret_R=(k1-k)/(k1+k);
% 
% k=min(min(G));
% k1=max(max(G));
% contrast_secret_G=(k1-k)/(k1+k);
% 
% k=min(min(B));
% k1=max(max(B));
% contrast_secret_B=(k1-k)/(k1+k);
% 
% (contrast_secret_R+contrast_secret_G+contrast_secret_B)/3;

numberOfShares=4;
blockSize=4;

%im1=single(im1);
rows=ROW/2^(log2(blockSize));% define how many rows of block
cols=COL/2^(log2(blockSize));% define how many colls of block 
sizeI = size(im);
blocks_R = mat2tiles(im_R, ceil(sizeI(1:2)./[rows cols]));%[rows cols]
blocks_G = mat2tiles(im_G, ceil(sizeI(1:2)./[rows cols]));
blocks_B = mat2tiles(im_B, ceil(sizeI(1:2)./[rows cols]));
U_matrix=[1 2 3 4];
V_matrix=perms(U_matrix);
[p,q]=size(V_matrix);
p1=1;
index=0;

for i = 1:rows/2
    for j=1:cols/2
    k=i+rows/2;
    l=j+cols/2;
    [u1_R,s1_R,v1_R]=svd(blocks_R{j,i});
    [u1_G,s1_G,v1_G]=svd(blocks_G{j,i});
    [u1_B,s1_B,v1_B]=svd(blocks_B{j,i});
    
    [u2_R,s2_R,v2_R]=svd(blocks_R{l,k});
    [u2_G,s2_G,v2_G]=svd(blocks_G{l,k});
    [u2_B,s2_B,v2_B]=svd(blocks_B{l,k});
    
    v1_R=v1_R';
    v1_G=v1_G';
    v1_B=v1_B';
    
    v2_R=v2_R';
    v2_G=v2_G';
    v2_B=v2_B';
      
    r = randi(factorial(numberOfShares),1);%factorial(numberOfShares)
    
    s1_1_R=s2_R(1,1)*u2_R(:,V_matrix(r,1))*v1_R(V_matrix(r,3),:)+s1_R(V_matrix(r,2),V_matrix(r,2))*u1_R(:,V_matrix(r,2))*v1_R(V_matrix(r,2),:)-s2_R(1,1)*u2_R(:,V_matrix(r,2))*v1_R(V_matrix(r,4),:);%+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)-s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)
    s2_1_R=s2_R(1,1)*u2_R(:,V_matrix(r,2))*v1_R(V_matrix(r,4),:)+s1_R(V_matrix(r,4),V_matrix(r,4))*u1_R(:,V_matrix(r,4))*v1_R(V_matrix(r,4),:)-s2_R(1,1)*u2_R(:,V_matrix(r,1))*v1_R(V_matrix(r,3),:);%+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)-s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)
    s3_1_R=s2_R(1,1)*u2_R(:,V_matrix(r,4))*v1_R(V_matrix(r,2),:)+s1_R(V_matrix(r,1),V_matrix(r,1))*u1_R(:,V_matrix(r,1))*v1_R(V_matrix(r,1),:)-s2_R(1,1)*u2_R(:,V_matrix(r,3))*v1_R(V_matrix(r,1),:);%+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)-s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)
    s4_1_R=s2_R(1,1)*u2_R(:,V_matrix(r,3))*v1_R(V_matrix(r,1),:)+s1_R(V_matrix(r,3),V_matrix(r,3))*u1_R(:,V_matrix(r,3))*v1_R(V_matrix(r,3),:)-s2_R(1,1)*u2_R(:,V_matrix(r,4))*v1_R(V_matrix(r,2),:);%+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)-s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)
    
    r = randi(factorial(numberOfShares),1);
    s1_2_R=s1_R(1,1)*u1_R(:,V_matrix(r,1))*v2_R(V_matrix(r,3),:)+s2_R(V_matrix(r,2),V_matrix(r,2))*u2_R(:,V_matrix(r,2))*v2_R(V_matrix(r,2),:)-s1_R(1,1)*u1_R(:,V_matrix(r,4))*v2_R(V_matrix(r,2),:);%+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)-s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)
    s2_2_R=s1_R(1,1)*u1_R(:,V_matrix(r,2))*v2_R(V_matrix(r,4),:)+s2_R(V_matrix(r,4),V_matrix(r,4))*u2_R(:,V_matrix(r,4))*v2_R(V_matrix(r,4),:)-s1_R(1,1)*u1_R(:,V_matrix(r,1))*v2_R(V_matrix(r,3),:);%+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)-s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)
    s3_2_R=s1_R(1,1)*u1_R(:,V_matrix(r,3))*v2_R(V_matrix(r,1),:)+s2_R(V_matrix(r,1),V_matrix(r,1))*u2_R(:,V_matrix(r,1))*v2_R(V_matrix(r,1),:)-s1_R(1,1)*u1_R(:,V_matrix(r,2))*v2_R(V_matrix(r,4),:);%+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)-s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)
    s4_2_R=s1_R(1,1)*u1_R(:,V_matrix(r,4))*v2_R(V_matrix(r,2),:)+s2_R(V_matrix(r,3),V_matrix(r,3))*u2_R(:,V_matrix(r,3))*v2_R(V_matrix(r,3),:)-s1_R(1,1)*u1_R(:,V_matrix(r,3))*v2_R(V_matrix(r,1),:);%+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)-s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)
     
    r = randi(factorial(numberOfShares),1);
    s1_1_G=s2_G(1,1)*u2_G(:,V_matrix(r,1))*v1_G(V_matrix(r,3),:)+s1_G(V_matrix(r,2),V_matrix(r,2))*u1_G(:,V_matrix(r,2))*v1_G(V_matrix(r,2),:)-s2_G(1,1)*u2_G(:,V_matrix(r,2))*v1_G(V_matrix(r,4),:);%+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)-s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)
    s2_1_G=s2_G(1,1)*u2_G(:,V_matrix(r,2))*v1_G(V_matrix(r,4),:)+s1_G(V_matrix(r,4),V_matrix(r,4))*u1_G(:,V_matrix(r,4))*v1_G(V_matrix(r,4),:)-s2_G(1,1)*u2_G(:,V_matrix(r,1))*v1_G(V_matrix(r,3),:);%+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)-s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)
    s3_1_G=s2_G(1,1)*u2_G(:,V_matrix(r,4))*v1_G(V_matrix(r,2),:)+s1_G(V_matrix(r,1),V_matrix(r,1))*u1_G(:,V_matrix(r,1))*v1_G(V_matrix(r,1),:)-s2_G(1,1)*u2_G(:,V_matrix(r,3))*v1_G(V_matrix(r,1),:);%+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)-s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)
    s4_1_G=s2_G(1,1)*u2_G(:,V_matrix(r,3))*v1_G(V_matrix(r,1),:)+s1_G(V_matrix(r,3),V_matrix(r,3))*u1_G(:,V_matrix(r,3))*v1_G(V_matrix(r,3),:)-s2_G(1,1)*u2_G(:,V_matrix(r,4))*v1_G(V_matrix(r,2),:);%+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)-s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)
    
    
    r = randi(factorial(numberOfShares),1);
    s1_2_G=s1_G(1,1)*u1_G(:,V_matrix(r,1))*v2_G(V_matrix(r,3),:)+s2_G(V_matrix(r,2),V_matrix(r,2))*u2_G(:,V_matrix(r,2))*v2_G(V_matrix(r,2),:)-s1_G(1,1)*u1_G(:,V_matrix(r,4))*v2_G(V_matrix(r,2),:);%+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)-s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)
    s2_2_G=s1_G(1,1)*u1_G(:,V_matrix(r,2))*v2_G(V_matrix(r,4),:)+s2_G(V_matrix(r,4),V_matrix(r,4))*u2_G(:,V_matrix(r,4))*v2_G(V_matrix(r,4),:)-s1_G(1,1)*u1_G(:,V_matrix(r,1))*v2_G(V_matrix(r,3),:);%+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)-s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)
    s3_2_G=s1_G(1,1)*u1_G(:,V_matrix(r,3))*v2_G(V_matrix(r,1),:)+s2_G(V_matrix(r,1),V_matrix(r,1))*u2_G(:,V_matrix(r,1))*v2_G(V_matrix(r,1),:)-s1_G(1,1)*u1_G(:,V_matrix(r,2))*v2_G(V_matrix(r,4),:);%+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)-s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)
    s4_2_G=s1_G(1,1)*u1_G(:,V_matrix(r,4))*v2_G(V_matrix(r,2),:)+s2_G(V_matrix(r,3),V_matrix(r,3))*u2_G(:,V_matrix(r,3))*v2_G(V_matrix(r,3),:)-s1_G(1,1)*u1_G(:,V_matrix(r,3))*v2_G(V_matrix(r,1),:);%+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)-s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)
     
    r = randi(factorial(numberOfShares),1);
    s1_1_B=s2_B(1,1)*u2_B(:,V_matrix(r,1))*v1_B(V_matrix(r,3),:)+s1_B(V_matrix(r,2),V_matrix(r,2))*u1_B(:,V_matrix(r,2))*v1_B(V_matrix(r,2),:)-s2_B(1,1)*u2_B(:,V_matrix(r,2))*v1_B(V_matrix(r,4),:);%+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)-s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)
    s2_1_B=s2_B(1,1)*u2_B(:,V_matrix(r,2))*v1_B(V_matrix(r,4),:)+s1_B(V_matrix(r,4),V_matrix(r,4))*u1_B(:,V_matrix(r,4))*v1_B(V_matrix(r,4),:)-s2_B(1,1)*u2_B(:,V_matrix(r,1))*v1_B(V_matrix(r,3),:);%+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)-s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)
    s3_1_B=s2_B(1,1)*u2_B(:,V_matrix(r,4))*v1_B(V_matrix(r,2),:)+s1_B(V_matrix(r,1),V_matrix(r,1))*u1_B(:,V_matrix(r,1))*v1_B(V_matrix(r,1),:)-s2_B(1,1)*u2_B(:,V_matrix(r,3))*v1_B(V_matrix(r,1),:);%+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)-s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)
    s4_1_B=s2_B(1,1)*u2_B(:,V_matrix(r,3))*v1_B(V_matrix(r,1),:)+s1_B(V_matrix(r,3),V_matrix(r,3))*u1_B(:,V_matrix(r,3))*v1_B(V_matrix(r,3),:)-s2_B(1,1)*u2_B(:,V_matrix(r,4))*v1_B(V_matrix(r,2),:);%+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)-s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)
    
    
    r = randi(factorial(numberOfShares),1);
    s1_2_B=s1_B(1,1)*u1_B(:,V_matrix(r,1))*v2_B(V_matrix(r,3),:)+s2_B(V_matrix(r,2),V_matrix(r,2))*u2_B(:,V_matrix(r,2))*v2_B(V_matrix(r,2),:)-s1_B(1,1)*u1_B(:,V_matrix(r,4))*v2_B(V_matrix(r,2),:);%+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)-s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)
    s2_2_B=s1_B(1,1)*u1_B(:,V_matrix(r,2))*v2_B(V_matrix(r,4),:)+s2_B(V_matrix(r,4),V_matrix(r,4))*u2_B(:,V_matrix(r,4))*v2_B(V_matrix(r,4),:)-s1_B(1,1)*u1_B(:,V_matrix(r,1))*v2_B(V_matrix(r,3),:);%+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)-s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)
    s3_2_B=s1_B(1,1)*u1_B(:,V_matrix(r,3))*v2_B(V_matrix(r,1),:)+s2_B(V_matrix(r,1),V_matrix(r,1))*u2_B(:,V_matrix(r,1))*v2_B(V_matrix(r,1),:)-s1_B(1,1)*u1_B(:,V_matrix(r,2))*v2_B(V_matrix(r,4),:);%+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)-s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)
    s4_2_B=s1_B(1,1)*u1_B(:,V_matrix(r,4))*v2_B(V_matrix(r,2),:)+s2_B(V_matrix(r,3),V_matrix(r,3))*u2_B(:,V_matrix(r,3))*v2_B(V_matrix(r,3),:)-s1_B(1,1)*u1_B(:,V_matrix(r,3))*v2_B(V_matrix(r,1),:);%+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)-s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)
        
    rand=randi(factorial(numberOfShares),1);   
     
        if rand==1
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s3_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s4_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s3_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s4_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s3_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s4_1_B;
        
        
        elseif rand==2
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
         
        elseif rand==3
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;

        elseif rand==4
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==5
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==6
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==7
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==8
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==9
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==10
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==11
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==12
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==13
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==14
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==15
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==16
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==17
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==18
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==19
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==20
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==21
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==22
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==23
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==24
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        end
        
        rand=randi(factorial(numberOfShares),1);
        
        if rand==1
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==2
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
         
        elseif rand==3
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;

        elseif rand==4
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==5
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        elseif rand==6
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==7
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==8
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==9
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==10
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==11
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==12
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==13
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==14
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==15
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==16
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==17
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==18
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==19
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==20
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        elseif rand==21
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==22
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==23
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==24
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        end
    p1=p1+1;
    end 
end  

q1=p1;
for i =rows/2+1:rows
    for j=1:cols/2
    k=i-rows/2;
    l=j+cols/2;
    [u1_R,s1_R,v1_R]=svd(blocks_R{j,i});
    [u1_G,s1_G,v1_G]=svd(blocks_G{j,i});
    [u1_B,s1_B,v1_B]=svd(blocks_B{j,i});
    
    [u2_R,s2_R,v2_R]=svd(blocks_R{l,k});
    [u2_G,s2_G,v2_G]=svd(blocks_G{l,k});
    [u2_B,s2_B,v2_B]=svd(blocks_B{l,k});
    
    v1_R=v1_R';
    v1_G=v1_G';
    v1_B=v1_B';
    
    v2_R=v2_R';
    v2_G=v2_G';
    v2_B=v2_B';
    
    r = randi(factorial(numberOfShares),1);%factorial(numberOfShares)
    
    s1_1_R=s2_R(1,1)*u2_R(:,V_matrix(r,1))*v1_R(V_matrix(r,3),:)+s1_R(V_matrix(r,2),V_matrix(r,2))*u1_R(:,V_matrix(r,2))*v1_R(V_matrix(r,2),:)-s2_R(1,1)*u2_R(:,V_matrix(r,2))*v1_R(V_matrix(r,4),:);%+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)-s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)
    s2_1_R=s2_R(1,1)*u2_R(:,V_matrix(r,2))*v1_R(V_matrix(r,4),:)+s1_R(V_matrix(r,4),V_matrix(r,4))*u1_R(:,V_matrix(r,4))*v1_R(V_matrix(r,4),:)-s2_R(1,1)*u2_R(:,V_matrix(r,1))*v1_R(V_matrix(r,3),:);%+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)-s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)
    s3_1_R=s2_R(1,1)*u2_R(:,V_matrix(r,4))*v1_R(V_matrix(r,2),:)+s1_R(V_matrix(r,1),V_matrix(r,1))*u1_R(:,V_matrix(r,1))*v1_R(V_matrix(r,1),:)-s2_R(1,1)*u2_R(:,V_matrix(r,3))*v1_R(V_matrix(r,1),:);%+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)-s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)
    s4_1_R=s2_R(1,1)*u2_R(:,V_matrix(r,3))*v1_R(V_matrix(r,1),:)+s1_R(V_matrix(r,3),V_matrix(r,3))*u1_R(:,V_matrix(r,3))*v1_R(V_matrix(r,3),:)-s2_R(1,1)*u2_R(:,V_matrix(r,4))*v1_R(V_matrix(r,2),:);%+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)-s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)
    
    r = randi(factorial(numberOfShares),1);
    s1_2_R=s1_R(1,1)*u1_R(:,V_matrix(r,1))*v2_R(V_matrix(r,3),:)+s2_R(V_matrix(r,2),V_matrix(r,2))*u2_R(:,V_matrix(r,2))*v2_R(V_matrix(r,2),:)-s1_R(1,1)*u1_R(:,V_matrix(r,4))*v2_R(V_matrix(r,2),:);%+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)-s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)
    s2_2_R=s1_R(1,1)*u1_R(:,V_matrix(r,2))*v2_R(V_matrix(r,4),:)+s2_R(V_matrix(r,4),V_matrix(r,4))*u2_R(:,V_matrix(r,4))*v2_R(V_matrix(r,4),:)-s1_R(1,1)*u1_R(:,V_matrix(r,1))*v2_R(V_matrix(r,3),:);%+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)-s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)
    s3_2_R=s1_R(1,1)*u1_R(:,V_matrix(r,3))*v2_R(V_matrix(r,1),:)+s2_R(V_matrix(r,1),V_matrix(r,1))*u2_R(:,V_matrix(r,1))*v2_R(V_matrix(r,1),:)-s1_R(1,1)*u1_R(:,V_matrix(r,2))*v2_R(V_matrix(r,4),:);%+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)-s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)
    s4_2_R=s1_R(1,1)*u1_R(:,V_matrix(r,4))*v2_R(V_matrix(r,2),:)+s2_R(V_matrix(r,3),V_matrix(r,3))*u2_R(:,V_matrix(r,3))*v2_R(V_matrix(r,3),:)-s1_R(1,1)*u1_R(:,V_matrix(r,3))*v2_R(V_matrix(r,1),:);%+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)-s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)
     
    r = randi(factorial(numberOfShares),1);
    s1_1_G=s2_G(1,1)*u2_G(:,V_matrix(r,1))*v1_G(V_matrix(r,3),:)+s1_G(V_matrix(r,2),V_matrix(r,2))*u1_G(:,V_matrix(r,2))*v1_G(V_matrix(r,2),:)-s2_G(1,1)*u2_G(:,V_matrix(r,2))*v1_G(V_matrix(r,4),:);%+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)-s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)
    s2_1_G=s2_G(1,1)*u2_G(:,V_matrix(r,2))*v1_G(V_matrix(r,4),:)+s1_G(V_matrix(r,4),V_matrix(r,4))*u1_G(:,V_matrix(r,4))*v1_G(V_matrix(r,4),:)-s2_G(1,1)*u2_G(:,V_matrix(r,1))*v1_G(V_matrix(r,3),:);%+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)-s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)
    s3_1_G=s2_G(1,1)*u2_G(:,V_matrix(r,4))*v1_G(V_matrix(r,2),:)+s1_G(V_matrix(r,1),V_matrix(r,1))*u1_G(:,V_matrix(r,1))*v1_G(V_matrix(r,1),:)-s2_G(1,1)*u2_G(:,V_matrix(r,3))*v1_G(V_matrix(r,1),:);%+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)-s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)
    s4_1_G=s2_G(1,1)*u2_G(:,V_matrix(r,3))*v1_G(V_matrix(r,1),:)+s1_G(V_matrix(r,3),V_matrix(r,3))*u1_G(:,V_matrix(r,3))*v1_G(V_matrix(r,3),:)-s2_G(1,1)*u2_G(:,V_matrix(r,4))*v1_G(V_matrix(r,2),:);%+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)-s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)
    
    
    r = randi(factorial(numberOfShares),1);
    s1_2_G=s1_G(1,1)*u1_G(:,V_matrix(r,1))*v2_G(V_matrix(r,3),:)+s2_G(V_matrix(r,2),V_matrix(r,2))*u2_G(:,V_matrix(r,2))*v2_G(V_matrix(r,2),:)-s1_G(1,1)*u1_G(:,V_matrix(r,4))*v2_G(V_matrix(r,2),:);%+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)-s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)
    s2_2_G=s1_G(1,1)*u1_G(:,V_matrix(r,2))*v2_G(V_matrix(r,4),:)+s2_G(V_matrix(r,4),V_matrix(r,4))*u2_G(:,V_matrix(r,4))*v2_G(V_matrix(r,4),:)-s1_G(1,1)*u1_G(:,V_matrix(r,1))*v2_G(V_matrix(r,3),:);%+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)-s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)
    s3_2_G=s1_G(1,1)*u1_G(:,V_matrix(r,3))*v2_G(V_matrix(r,1),:)+s2_G(V_matrix(r,1),V_matrix(r,1))*u2_G(:,V_matrix(r,1))*v2_G(V_matrix(r,1),:)-s1_G(1,1)*u1_G(:,V_matrix(r,2))*v2_G(V_matrix(r,4),:);%+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)-s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)
    s4_2_G=s1_G(1,1)*u1_G(:,V_matrix(r,4))*v2_G(V_matrix(r,2),:)+s2_G(V_matrix(r,3),V_matrix(r,3))*u2_G(:,V_matrix(r,3))*v2_G(V_matrix(r,3),:)-s1_G(1,1)*u1_G(:,V_matrix(r,3))*v2_G(V_matrix(r,1),:);%+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)-s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)
     
    r = randi(factorial(numberOfShares),1);
    s1_1_B=s2_B(1,1)*u2_B(:,V_matrix(r,1))*v1_B(V_matrix(r,3),:)+s1_B(V_matrix(r,2),V_matrix(r,2))*u1_B(:,V_matrix(r,2))*v1_B(V_matrix(r,2),:)-s2_B(1,1)*u2_B(:,V_matrix(r,2))*v1_B(V_matrix(r,4),:);%+s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)-s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)
    s2_1_B=s2_B(1,1)*u2_B(:,V_matrix(r,2))*v1_B(V_matrix(r,4),:)+s1_B(V_matrix(r,4),V_matrix(r,4))*u1_B(:,V_matrix(r,4))*v1_B(V_matrix(r,4),:)-s2_B(1,1)*u2_B(:,V_matrix(r,1))*v1_B(V_matrix(r,3),:);%+s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)-s1(V_matrix(r,4),V_matrix(r,4))*u1(:,V_matrix(r,2))*v1(V_matrix(r,3),:)
    s3_1_B=s2_B(1,1)*u2_B(:,V_matrix(r,4))*v1_B(V_matrix(r,2),:)+s1_B(V_matrix(r,1),V_matrix(r,1))*u1_B(:,V_matrix(r,1))*v1_B(V_matrix(r,1),:)-s2_B(1,1)*u2_B(:,V_matrix(r,3))*v1_B(V_matrix(r,1),:);%+s1(V_matrix(r,2),V_matrix(r,2))*u1(:,V_matrix(r,3))*v1(V_matrix(r,1),:)-s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)
    s4_1_B=s2_B(1,1)*u2_B(:,V_matrix(r,3))*v1_B(V_matrix(r,1),:)+s1_B(V_matrix(r,3),V_matrix(r,3))*u1_B(:,V_matrix(r,3))*v1_B(V_matrix(r,3),:)-s2_B(1,1)*u2_B(:,V_matrix(r,4))*v1_B(V_matrix(r,2),:);%+s1(V_matrix(r,1),V_matrix(r,1))*u1(:,V_matrix(r,4))*v1(V_matrix(r,2),:)-s1(V_matrix(r,3),V_matrix(r,3))*u1(:,V_matrix(r,1))*v1(V_matrix(r,4),:)
    
    
    r = randi(factorial(numberOfShares),1);
    s1_2_B=s1_B(1,1)*u1_B(:,V_matrix(r,1))*v2_B(V_matrix(r,3),:)+s2_B(V_matrix(r,2),V_matrix(r,2))*u2_B(:,V_matrix(r,2))*v2_B(V_matrix(r,2),:)-s1_B(1,1)*u1_B(:,V_matrix(r,4))*v2_B(V_matrix(r,2),:);%+s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)-s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)
    s2_2_B=s1_B(1,1)*u1_B(:,V_matrix(r,2))*v2_B(V_matrix(r,4),:)+s2_B(V_matrix(r,4),V_matrix(r,4))*u2_B(:,V_matrix(r,4))*v2_B(V_matrix(r,4),:)-s1_B(1,1)*u1_B(:,V_matrix(r,1))*v2_B(V_matrix(r,3),:);%+s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)-s2(V_matrix(r,4),V_matrix(r,4))*u2(:,V_matrix(r,2))*v2(V_matrix(r,3),:)
    s3_2_B=s1_B(1,1)*u1_B(:,V_matrix(r,3))*v2_B(V_matrix(r,1),:)+s2_B(V_matrix(r,1),V_matrix(r,1))*u2_B(:,V_matrix(r,1))*v2_B(V_matrix(r,1),:)-s1_B(1,1)*u1_B(:,V_matrix(r,2))*v2_B(V_matrix(r,4),:);%+s2(V_matrix(r,2),V_matrix(r,2))*u2(:,V_matrix(r,3))*v2(V_matrix(r,1),:)-s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)
    s4_2_B=s1_B(1,1)*u1_B(:,V_matrix(r,4))*v2_B(V_matrix(r,2),:)+s2_B(V_matrix(r,3),V_matrix(r,3))*u2_B(:,V_matrix(r,3))*v2_B(V_matrix(r,3),:)-s1_B(1,1)*u1_B(:,V_matrix(r,3))*v2_B(V_matrix(r,1),:);%+s2(V_matrix(r,1),V_matrix(r,1))*u2(:,V_matrix(r,4))*v2(V_matrix(r,2),:)-s2(V_matrix(r,3),V_matrix(r,3))*u2(:,V_matrix(r,1))*v2(V_matrix(r,4),:)
        
    rand=randi(factorial(numberOfShares),1);   
     
        if rand==1
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s3_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s4_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s3_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s4_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s3_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s4_1_B;
        
        
        elseif rand==2
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
         
        elseif rand==3
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;

        elseif rand==4
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==5
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==6
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==7
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==8
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==9
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==10
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==11
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==12
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==13
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==14
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==15
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==16
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==17
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==18
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==19
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==20
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==21
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==22
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==23
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        
        elseif rand==24
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_1_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_1_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_1_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_1_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_1_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_1_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_1_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_1_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_1_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_1_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_1_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_1_B;
        end
        
        rand=randi(factorial(numberOfShares),1);
        
        if rand==1
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==2
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
         
        elseif rand==3
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;

        elseif rand==4
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==5
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        elseif rand==6
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==7
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==8
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==9
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==10
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==11
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==12
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==13
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==14
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==15
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==16
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==17
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==18
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==19
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==20
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        elseif rand==21
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==22
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==23
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        
        elseif rand==24
        startrow =blockSize*(l-1)+1;     
        startcol =blockSize*(k-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_2_R;
        ss2_R(startrow : endrow, startcol : endcol) = s2_2_R;
        ss3_R(startrow : endrow, startcol : endcol) = s4_2_R;
        ss4_R(startrow : endrow, startcol : endcol) = s3_2_R;
        
        ss1_G(startrow : endrow, startcol : endcol) = s1_2_G;
        ss2_G(startrow : endrow, startcol : endcol) = s2_2_G;
        ss3_G(startrow : endrow, startcol : endcol) = s4_2_G;
        ss4_G(startrow : endrow, startcol : endcol) = s3_2_G;
        
        ss1_B(startrow : endrow, startcol : endcol) = s1_2_B;
        ss2_B(startrow : endrow, startcol : endcol) = s2_2_B;
        ss3_B(startrow : endrow, startcol : endcol) = s4_2_B;
        ss4_B(startrow : endrow, startcol : endcol) = s3_2_B;
        end
    p1=p1+1;
    
    end
end
im2_R=single(ss1_R+ss2_R+ss3_R+ss4_R);
im2_G=single(ss1_G+ss2_G+ss3_G+ss4_G);
im2_B=single(ss1_B+ss2_B+ss3_B+ss4_B);

im2_12_R=single(ss1_R+ss2_R);
im2_12_G=single(ss1_G+ss2_G);
im2_12_B=single(ss1_B+ss2_B);

n1_12=mean(mean(im2_12_R));
n2_12=mean(mean(im2_12_G));
n3_12=mean(mean(im2_12_B));

n1=mean(mean(im2_R));
n2=mean(mean(im2_G));
n3=mean(mean(im2_B));

n=ceil((n1+n2+n3)/3);
n12=ceil((n1_12+n2_12+n3_12)/3);

SS1=cat(3,ss1_R,ss1_G,ss1_B);
SS2=cat(3,ss2_R,ss2_G,ss2_B);
SS3=cat(3,ss3_R,ss3_G,ss3_B);
SS4=cat(3,ss4_R,ss4_G,ss4_B);

 SSS=single(SS1+SS2+SS3+SS4);

%im4=uint8(ss1+ss2+ss3);
for i=1:ROW
    for j=1:COL
        i2=mod(((n+1)*(i-1)-(j-1)),ROW)+1;
        j2=mod((-n*(i-1)+(j-1)),ROW)+1; 
        im3_R(i2,j2)=im2_R(i,j);
        im3_G(i2,j2)=im2_G(i,j);
        im3_B(i2,j2)=im2_B(i,j);

%         i22=mod(((n1+1)*(i-1)-(j-1)),R)+1;
%         j22=mod((-n1*(i-1)+(j-1)),R)+1; 
%         im5(i22,j22)=im4(i,j);
    end
end

for i=1:ROW
    for j=1:COL
        i2=mod(((n12+1)*(i-1)-(j-1)),ROW)+1;
        j2=mod((-n12*(i-1)+(j-1)),ROW)+1; 
        im3_12_R(i2,j2)=im2_12_R(i,j);
        im3_12_G(i2,j2)=im2_12_G(i,j);
        im3_12_B(i2,j2)=im2_12_B(i,j);

%         i22=mod(((n1+1)*(i-1)-(j-1)),R)+1;
%         j22=mod((-n1*(i-1)+(j-1)),R)+1; 
%         im5(i22,j22)=im4(i,j);
    end
end

im3=cat(3,im3_R,im3_G,im3_B);%reconstructed image
im3_12=cat(3,im3_12_R,im3_12_G,im3_12_B);
 
 Reconstructed_R=SSS(:, :, 1);
 Reconstructed_G= SSS(:, :, 2);
 Reconstructed_B= SSS(:, :, 3);
 
%  k=min(min(Reconstructed_R));
%  k1=max(max(Reconstructed_R));
%  contrast_reconstructed_R=(k1-k)/(k1+k);
%  
%  k=min(min(Reconstructed_G));
%  k1=max(max(Reconstructed_G));
%  contrast_reconstructed_G=(k1-k)/(k1+k);
%SSS=int32(SSS);
mean1=mean(reshape(single(SS1),[],1));
mean2=mean(reshape(single(SS2),[],1));
mean3=mean(reshape(single(SS2),[],1));
mean4=mean(reshape(single(SS4),[],1));
j=entropy(im);
figure(2);
imshow(mat2gray((SS1)));%share 1
%j1=entropy(SS1);
figure(3);
imshow(mat2gray((SS2)));%share 2
%j2=entropy(SS2);
figure(4);
imshow(mat2gray((SS3)));%share 3
%j3=entropy(SS3);
figure(5);
imshow(mat2gray((SS4)));%share 4
%j4=entropy(SS4);

figure(6);
imshow(mat2gray((im3)));%whith all 4 shares
%j1234=entropy(SSS);
figure(7)
imshow(mat2gray((im3_12)));%with shares 1 and 2
% 
figure(8);
imshow(mat2gray((SS1+SS2)));
%j12=entropy(SS1+SS2)
figure(9);
imshow(mat2gray(SS1+SS3));
% % j7=entropy(ss1+ss4)
 figure(10);
 imshow(mat2gray(SS1+SS4));
% % j8=entropy(ss2+ss3)
 figure(11);
 imshow(mat2gray(SS2+SS3));
% % j9=entropy(ss2+ss4)
 figure(12);
 imshow(mat2gray(SS2+SS4));
% % j10=entropy(ss3+ss4)
figure(13);
imshow(mat2gray(SS3+SS4));

figure(14);
imshow(mat2gray((SS1+SS2+SS3)));
%j123=entropy(SS1+SS2+SS3); 
 figure(15);
 imshow(mat2gray(SS1+SS2+SS4));
% % j12=entropy(ss1+ss2+ss4)
figure(16);
imshow(mat2gray(SS2+SS3+SS4));
% % j13=entropy(ss2+ss3+ss4)
% figure(15);
% imshow(mat2gray(ss1+ss2+ss3+ss4));
% 
%ss=ss1+ss2+ss3+ss4;
ss_1=SS1+SS2+SS3;
ss_2=SS1+SS2;
% OriginalImage_contrast=max(im(:)) - min(im(:));
% Recontructed_image_contrast = max(ss(:)) - min(ss(:));
% 
% % Recontructed_image_contrast_1 = max(ss_1(:)) - min(ss_1(:));
% % Recontructed_image_contrast_2 = max(ss_2(:)) - min(ss_2(:));
%Min_square_error=immse(im,SSS);
% Min_square_error_1_R=immse(im(:,:,1),ss_1(:,:,1));
% Min_square_error_1_G=immse(im(:,:,2),ss_1(:,:,2));
% Min_square_error_1_B=immse(im(:,:,3),ss_1(:,:,3));
% Min_square_error_1=(Min_square_error_1_R + Min_square_error_1_G + Min_square_error_1_B)/3;
% Min_square_error_2=immse(im,ss_2);

%PSNR_Value=psnr(im,SSS);
% PSNR_Value_1_R=psnr(im(:,:,1),ss_1(:,:,1));
% PSNR_Value_1_G=psnr(im(:,:,2),ss_1(:,:,2));
% PSNR_Value_1_B=psnr(im(:,:,3),ss_1(:,:,3));
% PSNR_Value_1=(PSNR_Value_1_R+PSNR_Value_1_G+PSNR_Value_1_B)/3;
% PSNR_Value_2=psnr(im,ss_2);
% 
%Stactural_similarity=ssim(im,SSS);
Stactural_similarity_1=ssim(single(SS1),single(im));
Stactural_similarity_2=ssim(single(SS2),single(im));
Stactural_similarity_3=ssim(single(SS3),single(im));
Stactural_similarity_4=ssim(single(SS4),single(im));
%Stactural_similarity_2=ssim(im,ss_2);
% 
% Correlation_coefficint_R=corr2(im(:,:,1),SSS(:,:,1));
% Correlation_coefficint_G=corr2(im(:,:,2),SSS(:,:,2));
% Correlation_coefficint_B=corr2(im(:,:,3),SSS(:,:,3));
% Correlation_coefficint=(Correlation_coefficint_R+Correlation_coefficint_G+Correlation_coefficint_B)/3;
% Correlation_coefficint_1_R=corr2(im(:,:,1),ss_1(:,:,1));
% Correlation_coefficint_1_G=corr2(im(:,:,2),ss_1(:,:,2));
% Correlation_coefficint_1_B=corr2(im(:,:,3),ss_1(:,:,3));
% Correlation_coefficint_1=(Correlation_coefficint_1_R+Correlation_coefficint_1_G+Correlation_coefficint_1_B)/3;
% Correlation_coefficint_2_R=corr2(im(:,:,1),ss_2(:,:,1));
% Correlation_coefficint_2_G=corr2(im(:,:,2),ss_2(:,:,2));
% Correlation_coefficint_2_B=corr2(im(:,:,3),ss_2(:,:,3));
% Correlation_coefficint_2=(Correlation_coefficint_2_R+Correlation_coefficint_2_G+Correlation_coefficint_2_B)/3;
% j_1234=entropy(ss1+ss2+ss3+ss4);
% j_123=entropy(ss1+ss2+ss3);
% j_12=entropy(ss1+ss2);
