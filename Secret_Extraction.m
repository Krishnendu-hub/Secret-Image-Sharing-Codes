
%load aecref.mat
filename = 'share_embedded_segment.wav';
[Y,Fs] = audioread(filename);
%sound(Y,Fs);
%Y=Y+K;
wavelet='db4';
level=1;
[cAr,cDr] = dwt(Y,wavelet,'mode','symw');
cDr1=(cDr+P);
cDr2=round(cDr1*10^4);
%im=imread('mri.tif');
im=imread('rice.png');
%im=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\conference\share4.png');
%im=imread('text.png');
[R,C]=size(im);
ImageSize=R*C;
for i=1:ImageSize%36560
    cDr2(i);
    bin1=dec2bin(cDr2(i),13);%dec2bin(cDr2(i))
    pixel_value_bin = extractBefore(bin1,9);
    pixel_value(i)=bin2dec(pixel_value_bin);
end
LengthDetails1=length(cDr2);
for j=1:G
    cDr2(ImageSize+j);
    bin1=dec2bin(cDr2(ImageSize+j),13);
    authentic_bits=extractBefore(bin1,5);
    Authentic_Str(j)=dec2hex(bin2dec(authentic_bits));
end
Authentic_Str;
pixel_value=uint8(reshape(pixel_value,[ImageSize,1]));
pixel_value1=reshape(pixel_value,[R,C]);
[New_hash,B]=hash1(pixel_value1);
figure(1)
imshow(mat2gray(pixel_value1));
figure(2)
%im=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\baboon512.bmp');
%im=imread('rice.png');

imshow(mat2gray(im));
im1=reshape(im,[R*C,1]);
%diff=(pixel_value-im1);
count=0;
for i=1:R*C
    if abs(diff(i))~=0
        i;
        count=count+1;
    end
end
if(Authentic_Str==New_hash)
    fprintf('Authentic Share');
else
    fprintf('Non Authentic Share');
end