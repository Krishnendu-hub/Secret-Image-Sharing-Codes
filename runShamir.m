%Addition based Shamir Sharing (10,10) scheme
clear all
close all
clc

%im=int32(imread('male.tiff'));
%im=int32(rgb2gray(imread('board.tif')));
%im=int32(imread('moon.tif'));
%im=im2double(imread('text.png'));
%im=int32(imread('rice.png'));
%im=int32((imread('mri.tif')));
%im=int32(imread('cameraman.bmp'));
%im=im2double(imread('lena.bmp'));
%im=int32(imread('mandi.tif'));
%im=int32(imread('pout.tif'));

im=int32(imread('cameraman.tif'));
%BIT_REPRESENTATION=bitget(im,8);
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\cameraman.bmp'));
%im=im2double(rgb2gray(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Airplane1.tiff')));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Clock.tiff'));
%im=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Couple.tiff'));
%im=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\baboon512.bmp');
%im=int32(rgb2gray(imread('peppers.png')));
%im=rgb2gray(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\news_paper512.jpg'));
figure(1)
imshow(mat2gray(im));
[row,col]=size(im);
k = 10; % the number of pieces of info which are sufficient for reconstruction
n = 10; % total number of pieces of info

gen=0;
%gen=0;


if(gen==1)
x=(1:n);
for i=1:n
    p=1;
    for j=1:n
		 if(j~=i) 
			 p = p* (0- x(j))/(x(i) - x(j));
             l(i)=p;
         end   
    end   
end
U_matrix=(1:n);
V_matrix=perms(U_matrix);

for i=1:row
   for j=1:col
      s=(im(i,j));
      d = ShamirSharing(s,k,n);
      c = d(1:n);
      %Weight_matrix=[4,-6,4,-1];%[4,-6,4,-1];%[2,-1];%[3,-3,1];%
      ran=randi(factorial(n),1);
      
      s1(i,j)=c(V_matrix(ran,1))*l(V_matrix(ran,1));
      s2(i,j)=c(V_matrix(ran,2))*l(V_matrix(ran,2));
      s3(i,j)=c(V_matrix(ran,3))*l(V_matrix(ran,3));
      s4(i,j)=c(V_matrix(ran,4))*l(V_matrix(ran,4));
      s5(i,j)=c(V_matrix(ran,5))*l(V_matrix(ran,5));
      s6(i,j)=c(V_matrix(ran,6))*l(V_matrix(ran,6));
      s7(i,j)=c(V_matrix(ran,7))*l(V_matrix(ran,7));
      s8(i,j)=c(V_matrix(ran,8))*l(V_matrix(ran,8));
      

      
   end
end
im=int32(im);
% s12=int32(s1+s2);
% s13=im2double(s1+s3);
% s14=im2double(s1+s4);
%s23=im2double(s2+s3);
% s24=im2double(s2+s4);
% s34=im2double(s3+s4);
% s123=int32(s1+s2+s3);
s123456=s1+s2+s3+s4+s5+s6;
s1_7=s1+s2+s3+s4+s5+s6+s7;
% s12_123=s123-s12;
% s=im-s1234;
figure(2)
imshow(mat2gray(s1));
% %imwrite(s1,'s1.png');
figure(3)
imshow(mat2gray(s2));
% %j2=entropy(s2);
 figure(4)
 imshow(mat2gray(s3));
% %j3=entropy(s3);
 figure(5)
 imshow(mat2gray(s4));
% %j4=entropy(s4);
% figure(6)
% imshow(mat2gray(s12));
% %j5=entropy(s12);
% figure(7)
% imshow(mat2gray(s23));
% figure(8)
% imshow(mat2gray(s123));
% %j6=entropy(s123);
% 
figure(9)
imshow(mat2gray(s123456));
% %j7=entropy(s1234);
% 
figure(10)
imshow(mat2gray(s1_7));
% 
% figure(11)
% histogram(s123);

% figure(12)
% imshow(mat2gray((s1_7)));

% figure(9)
% histogram(s12);
% figure(10)
% histogram(s13);
% figure(11)
% histogram(s14);
%histogram(s12_123);
% figure(12)
% histogram(s23);
% figure(13)
% histogram(s24);
% figure(14)
% histogram(s34);
% figure(10)
% histogram(s123);


% figure(11)
% histogram(s12);
% figure(13)
% histogram(s123);
% figure(14)
% imhist(s123-s12);
% figure(15)
% imhist(im-s1234);
%Min_square_error=immse(im,s1234);
% PSNR_Value=psnr(im,s1234);
% Correlation_coefficint=corr2(im,s1234);
% Entropy=entropy(s1234);
end
%figure(7)
% imshow(histogram(im));
% figure(8)
% imshow(histogram(ss));
% 
% Min_square_error=immse(im,ss);
% PSNR_Value=psnr(im,s1234)
% Stactural_similarity=ssim(im,ss);
% Normalized_Cross_Correlation=normxcorr2(im,ss);normxcorr2(template,A)
% MAE=mae(im,ss);
% Entropy=entropy(ss);
% ncc=0;
% for i=1:row
%    for j=1:col
%        ncc=ncc+(im(i,j)*ss(i,j))/(im(i,j))^2;
%    end
% end
% OriginalImage_contrast=max(im(:)) - min(im(:));
% Recontructed_image_contrast = max(ss(:)) - min(ss(:));

if(gen==0)
x=(1:n);
for i=1:n
    p=1;
    for j=1:n
		 if(j~=i) 
			 p = p* (0- x(j))/(x(i) - x(j));
             l(i)=p;
         end   
    end   
end
U_matrix=(1:n);
V_matrix=perms(U_matrix);
for i=1:row
   for j=1:col
      s =double(im(i,j));
      d =ShamirSharing(s,k,n);
      c =double(d(1:n));
      %Weight_matrix=[4,-6,4,-1];%[4,-6,4,-1];%[2,-1];%[3,-3,1];%
      ran=randi(factorial(n),1);
      prime=256;%997;
      s1(i,j)=mod(int64(c(V_matrix(ran,1))*l(V_matrix(ran,1))),prime);
      s2(i,j)=mod(int64(c(V_matrix(ran,2))*l(V_matrix(ran,2))),prime);
      s3(i,j)=mod(int64(c(V_matrix(ran,3))*l(V_matrix(ran,3))),prime);
      s4(i,j)=mod(int64(c(V_matrix(ran,4))*l(V_matrix(ran,4))),prime); 
      s5(i,j)=mod(int64(c(V_matrix(ran,5))*l(V_matrix(ran,5))),prime);
      s6(i,j)=mod(int64(c(V_matrix(ran,6))*l(V_matrix(ran,6))),prime);
      s7(i,j)=mod(int64(c(V_matrix(ran,7))*l(V_matrix(ran,7))),prime);
      s8(i,j)=mod(int64(c(V_matrix(ran,8))*l(V_matrix(ran,8))),prime);
      s9(i,j)=mod(int64(c(V_matrix(ran,9))*l(V_matrix(ran,9))),prime);
      s10(i,j)=mod(int64(c(V_matrix(ran,10))*l(V_matrix(ran,10))),prime);
%       s11(i,j)=mod(int64(c(V_matrix(ran,1))*l(V_matrix(ran,1))),prime);
%       s12(i,j)=mod(int64(c(V_matrix(ran,2))*l(V_matrix(ran,2))),prime);
%       s13(i,j)=mod(int64(c(V_matrix(ran,3))*l(V_matrix(ran,3))),prime);
%       s14(i,j)=mod(int64(c(V_matrix(ran,4))*l(V_matrix(ran,4))),prime); 
%       s15(i,j)=mod(int64(c(V_matrix(ran,5))*l(V_matrix(ran,5))),prime);
%       s16(i,j)=mod(int64(c(V_matrix(ran,6))*l(V_matrix(ran,6))),prime);
%       s17(i,j)=mod(int64(c(V_matrix(ran,7))*l(V_matrix(ran,7))),prime);
%       s18(i,j)=mod(int64(c(V_matrix(ran,8))*l(V_matrix(ran,8))),prime);
%       s19(i,j)=mod(int64(c(V_matrix(ran,9))*l(V_matrix(ran,9))),prime);
%       s20(i,j)=mod(int64(c(V_matrix(ran,10))*l(V_matrix(ran,10))),prime);

   end
end

s12=mod((s1+s2),prime);
% s13=mod((s1+s3),prime);
% s14=mod((s1+s4),prime);
% s23=mod((s2+s3),prime);
% s24=mod((s2+s4),prime);
% s34=mod((s3+s4),prime);
% s12=mod((s1+s2),prime);
s123=mod((s1+s2+s3),prime);
s1234=mod((s1+s2+s3+s4),prime);
tic
s12345=mod((s1+s2+s3+s4+s5),prime);
toc
% s123456=mod((s1+s2+s3+s4+s5+s6),prime);
% tic
% s1_7=mod((s1+s2+s3+s4+s5+s6+s7),prime);
% toc
tic
s12345678=mod((s1+s2+s3+s4+s5+s6+s7+s8),prime);
toc
tic
s1_9=mod((s1+s2+s3+s4+s5+s6+s7+s8+s9),prime);
toc
tic
s1_10=mod((s1+s2+s3+s4+s5+s6+s7+s8+s9+s10),prime);
toc
% tic
% s1_20=mod((s1+s2+s3+s4+s5+s6+s7+s8+s9+s10+s11+s12+s13+s14+s15+s16+s17+s18+s19+s20),prime);
% toc
% s12_123=(s123-s12);
% s=(im-s1234);
 figure(2)
 imshow(mat2gray(int32(s1)));
%j1=entropy(s1);
 figure(3)
 imshow(mat2gray(int32(s2)));
%j2=entropy(s2);
 figure(4)
 imshow(mat2gray(int32(s3)));
%j3=entropy(s3);
figure(5)
imshow(mat2gray(int32(s4)));
%j4=entropy(s4);
figure(6)
imshow(mat2gray(int32(s5)));
figure(7)
imshow(mat2gray(int32(s6)));
figure(8)
imshow(mat2gray(int32(s7)));
figure(9)
imshow(mat2gray(int32(s8)));
% figure(10)
% imshow(mat2gray(s123456));
% figure(11)
% imshow(mat2gray(s1_7));
figure(11)
imshow(mat2gray(int32(s12345678)));
 figure(12)
 imshow(mat2gray(int32(s1_9)));
 figure(13)
 imshow(mat2gray(int32(s1_10)));
% figure(14)
% imshow(mat2gray(int32(s1_20)));
% % figure(8)
% % imshow(mat2gray(C1234));
% % figure(9)
% % imshow(mat2gray(D1234));
% % figure(10)
% % imshow(mat2gray(E1234));
% % figure(11)
% % imshow(mat2gray(F1234));
% % figure(12)
% % imshow(mat2gray(G1234));
% % figure(13)
% % imshow(mat2gray(H1234));
% % figure(14)
% % imshow(mat2gray(ORR7));
% %figure(15)
% %imshow(mat2gray(s1234));
% figure(14)
% histogram(s1);
% % figure(10)
% % histogram(s13);
% % figure(11)
% % histogram(s);
% % figure(12)
% % histogram(s23);
% % figure(13)
% % histogram(s24);
% % figure(14)
% % histogram(s34);
% im=int16(im);
% s1234=int16(s1234);
% s=im-s1234;
% Min_square_error=immse(im,s1234);
% PSNR_Value=psnr(im,s1234);
% figure(7)
% histogram(s123);
% figure(8)
% histogram(s12);
% figure(9)
% histogram(s12-s123);
% figure(10)
% imhist(s);
% figure(11)
% histogram(im);
% figure(12)
% histogram(s1234);
% % Correlation_coefficint=corr2(im,s1234);
% % Entropy=entropy(s1234);
% int16(im)-int16(s1_10)
end