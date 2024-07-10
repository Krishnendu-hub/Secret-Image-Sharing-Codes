clear all;
close all;
im1=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\svd\Share 1st\share4.png');
im2=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\svd\Share 2nd\share4.png');
SSIM=ssim(im1,im2)
MSE=mse(im1,im2)
PSNR=psnr(im1,im2)
correlation=corr2(im1,im2);
%function [ NPCR ] = Cal_NPCR( imge,enc)
[rows,columns]=size(im1);
step=0;
for i=1:rows
    for j=1:columns
        if im1(i,j)~= im2(i,j)
           step=step+1;
        else 
             step=step+0;
        end
    end
end
NPCR =(step/(rows*columns))*100
%end

%function [UACI_value] = UACI( after_change,befor_change )
[row, col]=size(im1);
AB=[];
for i=1:row
    for j=1:col
        AB(i,j)=abs(im1(i,j)-im2(i,j));
    end
 end
    UACI_value = sum(AB(:))/(255*row*col)*100

PQ=[];
for i=1:row
    for j=1:col
        PQ(i,j)=(abs(im1(i,j)-im2(i,j)))^2;
    end
 end
    MSE = sum(PQ(:))/(row*col)
%end
im1=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\svd\Share 1st\share4.png'));
im2=im2double(imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\svd\Share 2nd\share4.png'));
[rho,pval] = corr(reshape(im1,row*col,1), reshape(im2,row*col,1), 'type', 'spearman')
[rho1,pval1] = corr(reshape(im1,row*col,1), reshape(im2,row*col,1), 'type', 'Kendall')