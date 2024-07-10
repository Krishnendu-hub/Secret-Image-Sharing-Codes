clear all;
close all;
clc;
%im=imread('peppers.png');
im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\images\lena.tif'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\cameraman.bmp'));
%im=im2double(imread('text.png'));
%im=im2double(imread('rice.png'));
%im=im2double(imread('testpat1.png'));
%im=im2double(imread('circles.png'));
%im=im2double(imread('liftingbody.png'));
%im=im2double(imread('mri.tif'));
%im=im2double(imread('male.tiff'));
%im=im2double(imread('lena.bmp'));
%im=im2double(imread('Mandrill.bmp'));
%im=im2double(rgb2gray(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\lena1.bmp')));
%im=im2double(rgb2gray(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Airplane1.tiff')));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\baboon512.bmp'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\ruler.tiff'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Clock.tiff'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Couple.tiff'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\boat.512.tiff'));
%im=im2double(rgb2gray(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\news_paper512.jpg')));

figure(1)
imshow(mat2gray(im));

R    = im(:, :, 1);
G    = im(:, :, 2);
B    = im(:, :, 3);
numberOfShares=8;
blockSize=8;
[ROW,COL,CHANNEL]=size(im);
rows=ROW/2^(log2(blockSize));% define how many rows of block
cols=COL/2^(log2(blockSize));% define how many colls of block 
sizeI = size(im);
blocks_R = mat2tiles(R, ceil(sizeI(1:2)./[rows cols]));%[rows cols]
blocks_G = mat2tiles(G, ceil(sizeI(1:2)./[rows cols]));
blocks_B = mat2tiles(B, ceil(sizeI(1:2)./[rows cols]));
U_matrix=[1:numberOfShares];
V_matrix=perms(U_matrix);
[p,q]=size(V_matrix);
p1=1;
index=0;
% ss1_R=zeros(ROW,COL);
% ss2_G=zeros(ROW,COL);
% ss3_B=zeros(ROW,COL);
% 
% ss2_R=zeros(ROW,COL);
% ss2_G=zeros(ROW,COL);
% ss2_B=zeros(ROW,COL);
% 
% ss3_R=zeros(ROW,COL);
% ss3_G=zeros(ROW,COL);
% ss3_B=zeros(ROW,COL);
% 
% ss4_R=zeros(ROW,COL);
% ss4_G=zeros(ROW,COL);
% ss4_B=zeros(ROW,COL);

for i = 1:rows
    for j=1:cols
    [u_R,s_R,v_R]=svd(blocks_R{j,i});
    [u_G,s_G,v_G]=svd(blocks_G{j,i});
    [u_B,s_B,v_B]=svd(blocks_B{j,i});
    
    v1_R=v_R';
    v1_G=v_G';
    v1_B=v_B';
    
%     k_R=(s_R(1,1)+s_R(2,2)+s_R(3,3)+s_R(4,4))/numberOfShares;
%     k_G=(s_G(1,1)+s_G(2,2)+s_G(3,3)+s_G(4,4))/numberOfShares;
%     k_B=(s_B(1,1)+s_B(2,2)+s_B(3,3)+s_B(4,4))/numberOfShares;
    
    u_R(:,1)=u_R(:,1)*s_R(1,1);
    u_R(:,2)=u_R(:,2)*s_R(2,2);
    u_R(:,3)=u_R(:,3)*s_R(3,3);
    u_R(:,4)=u_R(:,4)*s_R(4,4);
    u_R(:,5)=u_R(:,5)*s_R(5,5);
    u_R(:,6)=u_R(:,6)*s_R(6,6);
    u_R(:,7)=u_R(:,7)*s_R(7,7);
    u_R(:,8)=u_R(:,8)*s_R(8,8);
    
    u_G(:,1)=u_G(:,1)*s_G(1,1);
    u_G(:,2)=u_G(:,2)*s_G(2,2);
    u_G(:,3)=u_G(:,3)*s_G(3,3);
    u_G(:,4)=u_G(:,4)*s_G(4,4);
    u_G(:,5)=u_G(:,5)*s_G(5,5);
    u_G(:,6)=u_G(:,6)*s_G(6,6);
    u_G(:,7)=u_G(:,7)*s_G(7,7);
    u_G(:,8)=u_G(:,8)*s_G(8,8);
    
    u_B(:,1)=u_B(:,1)*s_B(1,1);
    u_B(:,2)=u_B(:,2)*s_B(2,2);
    u_B(:,3)=u_B(:,3)*s_B(3,3);
    u_B(:,4)=u_B(:,4)*s_B(4,4);
    u_B(:,5)=u_B(:,5)*s_B(5,5);
    u_B(:,6)=u_B(:,6)*s_B(6,6);
    u_B(:,7)=u_B(:,7)*s_B(7,7);
    u_B(:,8)=u_B(:,8)*s_B(8,8);
      
     r = randi(factorial(numberOfShares),1); 
     
     s8_R=u_R(:,V_matrix(r,8))*v1_R(V_matrix(r,6),:)+(u_R(:,V_matrix(r,8))-u_R(:,V_matrix(r,6)))*v1_R(V_matrix(r,8),:);
     s8_G=u_G(:,V_matrix(r,8))*v1_G(V_matrix(r,6),:)+(u_G(:,V_matrix(r,8))-u_G(:,V_matrix(r,6)))*v1_G(V_matrix(r,8),:);
     s8_B=u_B(:,V_matrix(r,8))*v1_B(V_matrix(r,6),:)+(u_B(:,V_matrix(r,8))-u_B(:,V_matrix(r,6)))*v1_B(V_matrix(r,8),:);
     
     s7_R=u_R(:,V_matrix(r,7))*v1_R(V_matrix(r,5),:)+(u_R(:,V_matrix(r,7))-u_R(:,V_matrix(r,5)))*v1_R(V_matrix(r,7),:);
     s7_G=u_G(:,V_matrix(r,7))*v1_G(V_matrix(r,5),:)+(u_G(:,V_matrix(r,7))-u_G(:,V_matrix(r,5)))*v1_G(V_matrix(r,7),:);
     s7_B=u_B(:,V_matrix(r,7))*v1_B(V_matrix(r,5),:)+(u_B(:,V_matrix(r,7))-u_B(:,V_matrix(r,5)))*v1_B(V_matrix(r,7),:);
     
     s6_R=u_R(:,V_matrix(r,6))*v1_R(V_matrix(r,8),:)+(u_R(:,V_matrix(r,6))-u_R(:,V_matrix(r,8)))*v1_R(V_matrix(r,6),:);
     s6_G=u_G(:,V_matrix(r,6))*v1_G(V_matrix(r,8),:)+(u_G(:,V_matrix(r,6))-u_G(:,V_matrix(r,8)))*v1_G(V_matrix(r,6),:);
     s6_B=u_B(:,V_matrix(r,6))*v1_B(V_matrix(r,8),:)+(u_B(:,V_matrix(r,6))-u_B(:,V_matrix(r,8)))*v1_B(V_matrix(r,6),:);
     
     s5_R=u_R(:,V_matrix(r,5))*v1_R(V_matrix(r,7),:)+(u_R(:,V_matrix(r,5))-u_R(:,V_matrix(r,7)))*v1_R(V_matrix(r,5),:);
     s5_G=u_G(:,V_matrix(r,5))*v1_G(V_matrix(r,7),:)+(u_G(:,V_matrix(r,5))-u_G(:,V_matrix(r,7)))*v1_G(V_matrix(r,5),:);
     s5_B=u_B(:,V_matrix(r,5))*v1_B(V_matrix(r,7),:)+(u_B(:,V_matrix(r,5))-u_B(:,V_matrix(r,7)))*v1_B(V_matrix(r,5),:);
     
     s4_R=u_R(:,V_matrix(r,4))*v1_R(V_matrix(r,1),:)+(u_R(:,V_matrix(r,4))-u_R(:,V_matrix(r,1)))*v1_R(V_matrix(r,4),:);
     s4_G=u_G(:,V_matrix(r,4))*v1_G(V_matrix(r,1),:)+(u_G(:,V_matrix(r,4))-u_G(:,V_matrix(r,1)))*v1_G(V_matrix(r,4),:);
     s4_B=u_B(:,V_matrix(r,4))*v1_B(V_matrix(r,1),:)+(u_B(:,V_matrix(r,4))-u_B(:,V_matrix(r,1)))*v1_B(V_matrix(r,4),:);
     
     s3_R=u_R(:,V_matrix(r,3))*v1_R(V_matrix(r,2),:)+(u_R(:,V_matrix(r,3))-u_R(:,V_matrix(r,2)))*v1_R(V_matrix(r,3),:);
     s3_G=u_G(:,V_matrix(r,3))*v1_G(V_matrix(r,2),:)+(u_G(:,V_matrix(r,3))-u_G(:,V_matrix(r,2)))*v1_G(V_matrix(r,3),:);
     s3_B=u_B(:,V_matrix(r,3))*v1_B(V_matrix(r,2),:)+(u_B(:,V_matrix(r,3))-u_B(:,V_matrix(r,2)))*v1_B(V_matrix(r,3),:);
     
     s2_R=u_R(:,V_matrix(r,2))*v1_R(V_matrix(r,3),:)+(u_R(:,V_matrix(r,2))-u_R(:,V_matrix(r,3)))*v1_R(V_matrix(r,2),:);
     s2_G=u_G(:,V_matrix(r,2))*v1_G(V_matrix(r,3),:)+(u_G(:,V_matrix(r,2))-u_G(:,V_matrix(r,3)))*v1_G(V_matrix(r,2),:);
     s2_B=u_B(:,V_matrix(r,2))*v1_B(V_matrix(r,3),:)+(u_B(:,V_matrix(r,2))-u_B(:,V_matrix(r,3)))*v1_B(V_matrix(r,2),:);
     
     s1_R=u_R(:,V_matrix(r,1))*v1_R(V_matrix(r,4),:)+(u_R(:,V_matrix(r,1))-u_R(:,V_matrix(r,4)))*v1_R(V_matrix(r,1),:);
     s1_G=u_G(:,V_matrix(r,1))*v1_G(V_matrix(r,4),:)+(u_G(:,V_matrix(r,1))-u_G(:,V_matrix(r,4)))*v1_G(V_matrix(r,1),:);
     s1_B=u_B(:,V_matrix(r,1))*v1_B(V_matrix(r,4),:)+(u_B(:,V_matrix(r,1))-u_B(:,V_matrix(r,4)))*v1_B(V_matrix(r,1),:);
     
        
         if mod(p1,numberOfShares)==0
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize-1;
        ss1_R(startrow : endrow, startcol : endcol) = s1_R;
        ss1_G(startrow : endrow, startcol : endcol) = s1_G;
        ss1_B(startrow : endrow, startcol : endcol) = s1_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s2_R;
        ss2_G(startrow : endrow, startcol : endcol) = s2_G;
        ss2_B(startrow : endrow, startcol : endcol) = s2_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s3_R;
        ss3_G(startrow : endrow, startcol : endcol) = s3_G;
        ss3_B(startrow : endrow, startcol : endcol) = s3_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s4_R;
        ss4_G(startrow : endrow, startcol : endcol) = s4_G;
        ss4_B(startrow : endrow, startcol : endcol) = s4_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s5_R;
        ss5_G(startrow : endrow, startcol : endcol) = s5_G;
        ss5_B(startrow : endrow, startcol : endcol) = s5_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s6_R;
        ss6_G(startrow : endrow, startcol : endcol) = s6_G;
        ss6_B(startrow : endrow, startcol : endcol) = s6_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s7_R;
        ss7_G(startrow : endrow, startcol : endcol) = s7_G;
        ss7_B(startrow : endrow, startcol : endcol) = s7_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s8_R;
        ss8_G(startrow : endrow, startcol : endcol) = s8_G;
        ss8_B(startrow : endrow, startcol : endcol) = s8_B;
        
        elseif mod(p1,numberOfShares)==1
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize-1;
        ss1_R(startrow : endrow, startcol : endcol) = s2_R;
        ss1_G(startrow : endrow, startcol : endcol) = s2_G;
        ss1_B(startrow : endrow, startcol : endcol) = s2_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s3_R;
        ss2_G(startrow : endrow, startcol : endcol) = s3_G;
        ss2_B(startrow : endrow, startcol : endcol) = s3_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s4_R;
        ss3_G(startrow : endrow, startcol : endcol) = s4_G;
        ss3_B(startrow : endrow, startcol : endcol) = s4_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s5_R;
        ss4_G(startrow : endrow, startcol : endcol) = s5_G;
        ss4_B(startrow : endrow, startcol : endcol) = s5_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s6_R;
        ss5_G(startrow : endrow, startcol : endcol) = s6_G;
        ss5_B(startrow : endrow, startcol : endcol) = s6_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s7_R;
        ss6_G(startrow : endrow, startcol : endcol) = s7_G;
        ss6_B(startrow : endrow, startcol : endcol) = s7_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s8_R;
        ss7_G(startrow : endrow, startcol : endcol) = s8_G;
        ss7_B(startrow : endrow, startcol : endcol) = s8_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s1_R;
        ss8_G(startrow : endrow, startcol : endcol) = s1_G;
        ss8_B(startrow : endrow, startcol : endcol) = s1_B;
        
        elseif mod(p1,numberOfShares)==2
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize-1;
        ss1_R(startrow : endrow, startcol : endcol) = s3_R;
        ss1_G(startrow : endrow, startcol : endcol) = s3_G;
        ss1_B(startrow : endrow, startcol : endcol) = s3_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s4_R;
        ss2_G(startrow : endrow, startcol : endcol) = s4_G;
        ss2_B(startrow : endrow, startcol : endcol) = s4_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s5_R;
        ss3_G(startrow : endrow, startcol : endcol) = s5_G;
        ss3_B(startrow : endrow, startcol : endcol) = s5_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s6_R;
        ss4_G(startrow : endrow, startcol : endcol) = s6_G;
        ss4_B(startrow : endrow, startcol : endcol) = s6_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s7_R;
        ss5_G(startrow : endrow, startcol : endcol) = s7_G;
        ss5_B(startrow : endrow, startcol : endcol) = s7_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s8_R;
        ss6_G(startrow : endrow, startcol : endcol) = s8_G;
        ss6_B(startrow : endrow, startcol : endcol) = s8_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s1_R;
        ss7_G(startrow : endrow, startcol : endcol) = s1_G;
        ss7_B(startrow : endrow, startcol : endcol) = s1_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s2_R;
        ss8_G(startrow : endrow, startcol : endcol) = s2_G;
        ss8_B(startrow : endrow, startcol : endcol) = s2_B;
        
        elseif mod(p1,numberOfShares)==3
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s4_R;
        ss1_G(startrow : endrow, startcol : endcol) = s4_G;
        ss1_B(startrow : endrow, startcol : endcol) = s4_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s5_R;
        ss2_G(startrow : endrow, startcol : endcol) = s5_G;
        ss2_B(startrow : endrow, startcol : endcol) = s5_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s6_R;
        ss3_G(startrow : endrow, startcol : endcol) = s6_G;
        ss3_B(startrow : endrow, startcol : endcol) = s6_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s7_R;
        ss4_G(startrow : endrow, startcol : endcol) = s7_G;
        ss4_B(startrow : endrow, startcol : endcol) = s7_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s8_R;
        ss5_G(startrow : endrow, startcol : endcol) = s8_G;
        ss5_B(startrow : endrow, startcol : endcol) = s8_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s1_R;
        ss6_G(startrow : endrow, startcol : endcol) = s1_G;
        ss6_B(startrow : endrow, startcol : endcol) = s1_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s2_R;
        ss7_G(startrow : endrow, startcol : endcol) = s2_G;
        ss7_B(startrow : endrow, startcol : endcol) = s2_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s3_R;
        ss8_G(startrow : endrow, startcol : endcol) = s3_G;
        ss8_B(startrow : endrow, startcol : endcol) = s3_B;
        
        elseif mod(p1,numberOfShares)==3
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s4_R;
        ss1_G(startrow : endrow, startcol : endcol) = s4_G;
        ss1_B(startrow : endrow, startcol : endcol) = s4_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s5_R;
        ss2_G(startrow : endrow, startcol : endcol) = s5_G;
        ss2_B(startrow : endrow, startcol : endcol) = s5_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s6_R;
        ss3_G(startrow : endrow, startcol : endcol) = s6_G;
        ss3_B(startrow : endrow, startcol : endcol) = s6_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s7_R;
        ss4_G(startrow : endrow, startcol : endcol) = s7_G;
        ss4_B(startrow : endrow, startcol : endcol) = s7_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s8_R;
        ss5_G(startrow : endrow, startcol : endcol) = s8_G;
        ss5_B(startrow : endrow, startcol : endcol) = s8_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s1_R;
        ss6_G(startrow : endrow, startcol : endcol) = s1_G;
        ss6_B(startrow : endrow, startcol : endcol) = s1_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s2_R;
        ss7_G(startrow : endrow, startcol : endcol) = s2_G;
        ss7_B(startrow : endrow, startcol : endcol) = s2_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s3_R;
        ss8_G(startrow : endrow, startcol : endcol) = s3_G;
        ss8_B(startrow : endrow, startcol : endcol) = s3_B;
        
        elseif mod(p1,numberOfShares)==4
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s5_R;
        ss1_G(startrow : endrow, startcol : endcol) = s5_G;
        ss1_B(startrow : endrow, startcol : endcol) = s5_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s6_R;
        ss2_G(startrow : endrow, startcol : endcol) = s6_G;
        ss2_B(startrow : endrow, startcol : endcol) = s6_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s7_R;
        ss3_G(startrow : endrow, startcol : endcol) = s7_G;
        ss3_B(startrow : endrow, startcol : endcol) = s7_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s8_R;
        ss4_G(startrow : endrow, startcol : endcol) = s8_G;
        ss4_B(startrow : endrow, startcol : endcol) = s8_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s1_R;
        ss5_G(startrow : endrow, startcol : endcol) = s1_G;
        ss5_B(startrow : endrow, startcol : endcol) = s1_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s2_R;
        ss6_G(startrow : endrow, startcol : endcol) = s2_G;
        ss6_B(startrow : endrow, startcol : endcol) = s2_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s3_R;
        ss7_G(startrow : endrow, startcol : endcol) = s3_G;
        ss7_B(startrow : endrow, startcol : endcol) = s3_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s4_R;
        ss8_G(startrow : endrow, startcol : endcol) = s4_G;
        ss8_B(startrow : endrow, startcol : endcol) = s4_B;
        
        elseif mod(p1,numberOfShares)==5
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s6_R;
        ss1_G(startrow : endrow, startcol : endcol) = s6_G;
        ss1_B(startrow : endrow, startcol : endcol) = s6_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s7_R;
        ss2_G(startrow : endrow, startcol : endcol) = s7_G;
        ss2_B(startrow : endrow, startcol : endcol) = s7_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s8_R;
        ss3_G(startrow : endrow, startcol : endcol) = s8_G;
        ss3_B(startrow : endrow, startcol : endcol) = s8_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s1_R;
        ss4_G(startrow : endrow, startcol : endcol) = s1_G;
        ss4_B(startrow : endrow, startcol : endcol) = s1_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s2_R;
        ss5_G(startrow : endrow, startcol : endcol) = s2_G;
        ss5_B(startrow : endrow, startcol : endcol) = s2_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s3_R;
        ss6_G(startrow : endrow, startcol : endcol) = s3_G;
        ss6_B(startrow : endrow, startcol : endcol) = s3_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s4_R;
        ss7_G(startrow : endrow, startcol : endcol) = s4_G;
        ss7_B(startrow : endrow, startcol : endcol) = s4_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s5_R;
        ss8_G(startrow : endrow, startcol : endcol) = s5_G;
        ss8_B(startrow : endrow, startcol : endcol) = s5_B;
        
        elseif mod(p1,numberOfShares)==6
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s7_R;
        ss1_G(startrow : endrow, startcol : endcol) = s7_G;
        ss1_B(startrow : endrow, startcol : endcol) = s7_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s8_R;
        ss2_G(startrow : endrow, startcol : endcol) = s8_G;
        ss2_B(startrow : endrow, startcol : endcol) = s8_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s1_R;
        ss3_G(startrow : endrow, startcol : endcol) = s1_G;
        ss3_B(startrow : endrow, startcol : endcol) = s1_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s2_R;
        ss4_G(startrow : endrow, startcol : endcol) = s2_G;
        ss4_B(startrow : endrow, startcol : endcol) = s2_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s3_R;
        ss5_G(startrow : endrow, startcol : endcol) = s3_G;
        ss5_B(startrow : endrow, startcol : endcol) = s3_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s4_R;
        ss6_G(startrow : endrow, startcol : endcol) = s4_G;
        ss6_B(startrow : endrow, startcol : endcol) = s4_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s5_R;
        ss7_G(startrow : endrow, startcol : endcol) = s5_G;
        ss7_B(startrow : endrow, startcol : endcol) = s5_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s6_R;
        ss8_G(startrow : endrow, startcol : endcol) = s6_G;
        ss8_B(startrow : endrow, startcol : endcol) = s6_B;
        
        elseif mod(p1,numberOfShares)==7
        startrow =blockSize*(j-1)+1;     
        startcol =blockSize*(i-1)+1;
        endrow = startrow + blockSize -1;
        endcol = startcol + blockSize -1;
        ss1_R(startrow : endrow, startcol : endcol) = s8_R;
        ss1_G(startrow : endrow, startcol : endcol) = s8_G;
        ss1_B(startrow : endrow, startcol : endcol) = s8_B;
        
        ss2_R(startrow : endrow, startcol : endcol) = s1_R;
        ss2_G(startrow : endrow, startcol : endcol) = s1_G;
        ss2_B(startrow : endrow, startcol : endcol) = s1_B;
        
        ss3_R(startrow : endrow, startcol : endcol) = s2_R;
        ss3_G(startrow : endrow, startcol : endcol) = s2_G;
        ss3_B(startrow : endrow, startcol : endcol) = s2_B;
        
        ss4_R(startrow : endrow, startcol : endcol) = s3_R;
        ss4_G(startrow : endrow, startcol : endcol) = s3_G;
        ss4_B(startrow : endrow, startcol : endcol) = s3_B;
        
        ss5_R(startrow : endrow, startcol : endcol) = s4_R;
        ss5_G(startrow : endrow, startcol : endcol) = s4_G;
        ss5_B(startrow : endrow, startcol : endcol) = s4_B;
        
        ss6_R(startrow : endrow, startcol : endcol) = s5_R;
        ss6_G(startrow : endrow, startcol : endcol) = s5_G;
        ss6_B(startrow : endrow, startcol : endcol) = s5_B;
        
        ss7_R(startrow : endrow, startcol : endcol) = s6_R;
        ss7_G(startrow : endrow, startcol : endcol) = s6_G;
        ss7_B(startrow : endrow, startcol : endcol) = s6_B;
        
        ss8_R(startrow : endrow, startcol : endcol) = s7_R;
        ss8_G(startrow : endrow, startcol : endcol) = s7_G;
        ss8_B(startrow : endrow, startcol : endcol) = s7_B;
        end
    p1=p1+1;
    end
     
end
            SS1=cat(3,ss1_R,ss1_G,ss1_B);
            SS2=cat(3,ss2_R,ss2_G,ss2_B);
            SS3=cat(3,ss3_R,ss3_G,ss3_B);
            SS4=cat(3,ss4_R,ss4_G,ss4_B);
            SS5=cat(3,ss5_R,ss5_G,ss5_B);
            SS6=cat(3,ss6_R,ss6_G,ss6_B);
            SS7=cat(3,ss7_R,ss7_G,ss7_B);
            SS8=cat(3,ss8_R,ss8_G,ss8_B);            
            
%              SS1=cat(3,ss1_R(:,:),ss1_G(:,:),ss1_B(:,:));
%              SS2=cat(3,ss2_R(:,:),ss2_G(:,:),ss2_B(:,:));
%              SS3=cat(3,ss3_R(:,:),ss1_G(:,:),ss3_B(:,:));
%              SS4=cat(3,ss4_R(:,:),ss4_G(:,:),ss4_B(:,:));


 SSS=SS1+SS2+SS3+SS4+SS5+SS6+SS7+SS8;

%j=entropy(im)
figure(2);
imshow(mat2gray(SS1));
%j1=entropy(ss1)
figure(3);
imshow(mat2gray(SS2));
%j2=entropy(ss2)
figure(4);
imshow(mat2gray(SS3));
%j3=entropy(ss3)
figure(5);
imshow(mat2gray(SS4));
figure(6);
imshow(mat2gray(SS5));
figure(7);
imshow(mat2gray(SS6));
figure(8);
imshow(mat2gray(SS7));

figure(9)
imshow(mat2gray(SS8));

figure(10)
imshow(mat2gray(SSS));
% 
% figure(7);
% imshow(mat2gray(ss1+ss3));
% % j6=entropy(ss1+ss3)
% figure(8);
% imshow(mat2gray(ss1+ss4));
% % j7=entropy(ss1+ss4)
% figure(9);
% imshow(mat2gray(ss2+ss3));
% % j8=entropy(ss2+ss3)
% figure(10);
% imshow(mat2gray(ss2+ss4));
% % j9=entropy(ss2+ss4)
% figure(11);
% imshow(mat2gray(ss3+ss4));
% % j10=entropy(ss3+ss4)
% figure(12);
% imshow(mat2gray(ss1+ss2+ss3));
% 
% figure(13);
% imshow(mat2gray(ss1+ss2+ss4));
% % j12=entropy(ss1+ss2+ss4)
% figure(14);
% imshow(mat2gray(ss2+ss3+ss4));
% % j13=entropy(ss2+ss3+ss4)
% figure(15);
% imshow(mat2gray(ss1+ss2+ss3+ss4));
% 
% ss=ss1+ss2+ss3+ss4;
% ss_1=ss1+ss2+ss3;
% ss_2=ss1+ss2;
% OriginalImage_contrast=max(im(:)) - min(im(:));
% Recontructed_image_contrast = max(ss(:)) - min(ss(:));
% 
% % Recontructed_image_contrast_1 = max(ss_1(:)) - min(ss_1(:));
% % Recontructed_image_contrast_2 = max(ss_2(:)) - min(ss_2(:));
% Min_square_error=immse(im,ss);
% Min_square_error_1=immse(im,ss_1);
% Min_square_error_2=immse(im,ss_2);
% 
% PSNR_Value=psnr(im,ss);
% PSNR_Value_1=psnr(im,ss_1);
% PSNR_Value_2=psnr(im,ss_2);
% 
% Stactural_similarity=ssim(im,ss);
% Stactural_similarity_1=ssim(im,ss_1);
% Stactural_similarity_2=ssim(im,ss_2);
% 
% Correlation_coefficint=corr2(im,ss);
% Correlation_coefficint_1=corr2(im,ss_1);
% Correlation_coefficint_2=corr2(im,ss_2);
% 
% j_1234=entropy(ss1+ss2+ss3+ss4);
% j_123=entropy(ss1+ss2+ss3);
% j_12=entropy(ss1+ss2);