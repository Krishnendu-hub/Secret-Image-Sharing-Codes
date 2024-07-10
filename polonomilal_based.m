clear all;
close all;

%s="enter the value n for making maximum share from image";
%n=input(s);
%s="enter the vlaue k for construction of image";
%k = input(s);
n=10;
k=10;

img_mat = (imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\misc\Couple.tiff'));

%img_mat = imread('cameraman.tif');
%img_mat=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\baboon512.bmp');
%img_mat = imread('moon.tif');
[rows,cols] = size(img_mat);  %size of image
num_of_img_mat_ele = rows*cols;   %total number of element in matrix

subshare = zeros(rows,cols/2,n);

P = 65537; % prime number
temp1 = zeros(rows,cols/2);

 % Sharing Phase start
count=0;
for i=1:rows
j=1;
while j <= cols
    
    pixel1bin=dec2bin( img_mat(i,j),8);
    pixel2bin=dec2bin( img_mat(i,j+1),8);
    pixelbin=strcat(pixel1bin,pixel2bin);
    newpixel=bin2dec(pixelbin);
    temp1(i,(j+1)/2)=newpixel;

    random = zeros(1,k-1);
    for ii=1:k-1
        random(ii) = randi([0,65535]);
    end

    for snum=1:n
     
        val = mod( newpixel + find_poly_val(snum,k,random) , P );

        if val == P-1 % P-1 = 65536
            % s="yes" 
            count=count+1;
            break
        end
        subshare(i,(j+1)/2,snum) = val;
    end

    if val ~= P-1
        j = j + 2;
    end

end
end

% Sharing Phase has done---------- 

% Reconstruction Phase start:------------
tic
temp2 = zeros(rows,cols/2);
for ri=1:rows
    for rj=1:cols/2
        sum = 0;
        for si=1:k
             sum = sum + get_val_of_interpolation(si,k,subshare(ri,rj,si));
        end
         sum = mod(sum,P);
         if sum < 0  
             sum = sum + P;
         end
         temp2(ri,rj)=sum;
    end
end

res_img=zeros(rows,cols);

for ri=1:rows
    rj=1;
   while rj <= cols
     pixel_temp2 = temp2(ri,(rj+1)/2);
     binstr = dec2bin(pixel_temp2,16);
     pixel1str = extractBefore(binstr,9);
     pixel2str = extractAfter(binstr,8);
     first_pixel = bin2dec(pixel1str);
     second_pixel = bin2dec(pixel2str);

     res_img(ri,rj)=first_pixel;
     res_img(ri,rj+1)=second_pixel;

     rj=rj+2;
    end
end
toc
% Reconstruction Phase has done

% original input image
figure(1)
imshow(mat2gray(img_mat));


% N Share image
figure(2)
for j=1:n
    %subplot(1,n,j);
    figure(j+1)
    imshow(uint16(subshare(:,:,j)))
    imwrite(uint16(subshare(:,:,2)),'share2.png');
   
end
%subplot(1,n,j);
%     figure(2)
%     imshow(uint16(subshare(:,:,1)))
%     imwrite(uint16(subshare(:,:,1)),'share1.png');
% 
%     figure(3)
%     imshow(uint16(subshare(:,:,2)))
%     imwrite(uint16(subshare(:,:,2)),'share2.png');
% 
%     figure(4)
%     imshow(uint16(subshare(:,:,3)))
%     imwrite(uint16(subshare(:,:,3)),'share3.png');
% 
%     figure(5)
%     imshow(uint16(subshare(:,:,4)))
%     imwrite(uint16(subshare(:,:,4)),'share4.png');
   
%sh1=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\conference\share1.png');
% reconstructed image
figure, imshow(mat2gray(res_img))

% function for finding a1(x) + a2(x^2) + ......... + a(k-1)(x^(k-1));
count1=0;
for i=1:rows
    for j=1:cols/2
        if subshare(i,j,1)<255
            count1=count1+1;
        end
    end
end

function [res] = find_poly_val(xval,k,random)
x=xval;
res=0;
for i=1:k-1
    res=res + random(i)*x;
    x=x*xval;
end
end

% function for finding langrange value (secret pixel) using  y1, y2, ......, yn so that we can find ao ;
function [res] = get_val_of_interpolation(pi,k,y)
     u=1;
     d=1;

   for ii=1:k
    if(ii ~= pi)
        u=u*(-ii);
        d=d*(pi-ii);
    end
   end
       res=floor((y*u)/d);
end


