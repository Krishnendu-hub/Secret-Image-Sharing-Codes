clear all;
close all;
%load aecref.mat
%filename = 'violin.wav';
%filename = 'handel.wav';
filename = 'guitartune.wav';
%filename = 'H:\Jrf Book\image processing\Dr.S.M\New folder\Images\conference\pca_ica\pca_ica\audio\source2.wav';

%filename = 'H:\Jrf Book\image processing\Dr.S.M\New folder\Images\conference\pca_ica\pca_ica\audio\source4.wav';
%filename='H:\Jrf Book\image processing\Dr.S.M\New folder\Images\conference\pca_ica\pca_ica\Ho-Gaya-Hai-Tujhko-DDLJ.wav';
%filename='source5.wav';
%filename = 'bluewhale.wav';
[y,Fs] = audioread(filename);
N=length(y);
time=N/Fs;
%Read a grayscale image for embedding with Audio Signal
%im=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\baboon512.bmp');
%im=imread('rice.png');
%im=imread('mri.tif');
im=imread('H:\Jrf Book\image processing\Dr.S.M\New folder\Images\conference\share4.png');
%im=imread('text.png');
[hash_im,A]=hash1(im);
[R,C]=size(im);
im1=reshape(im,[R*C,1]);
%sound(y,Fs);

% Scalling the input Audio Signal
K=abs(max(y));
y1=(y+K);
N = length(y1); % sample lenth
slength = N/Fs; % total time span of audio signal
t = linspace(0, N/Fs, N);
plot(t, y); % pplots the audio
xlabel('Time')
ylabel('Audio Signal')
%sound(y1,Fs);

N = length(y1); % sample lenth
slength = N/Fs; % total time span of audio signal
t = linspace(0, N/Fs, N);
wavelet='db4';
level=1;
[cA,cD] = dwt(y1,wavelet,'mode','symw'); %fix(log2(length(y1)))
P=abs(max(cD));
cD1=(cD+P);
cD2=round(cD1*10^4);

%Make condition for varerous length secret image and audio signal. 
ImageSize=size(im1,1);
LengthDetais=length(cD);
if ImageSize<=LengthDetais

%Embedding secret pixels into the Audio signal
    for i=1:ImageSize%36560
        str3=dec2bin(im1(i),8);
        
        str3_1=extractBefore(str3,5);
        str3_2=extractAfter(str3,4);

        bin=dec2bin(cD2(i),13);
        str2_1=extractBefore(bin,10);
        bin1=dec2bin(cD2(i*2),13);
        str2_2=extractBefore(bin1,10);

        str4_1=strcat(str2_1,str3_1);
        str4_2=strcat(str2_2,str3_2);

        cD3(i)=bin2dec(str4_1);
        cD3(i*2)=bin2dec(str4_2);
    end

    for i=1:ImageSize%36560
        str3=dec2bin(im1(i),8);
        bin=dec2bin(cD2(i),13);
        str2=extractAfter(bin,8);
        str4=strcat(str3,str2);
        cD3(i)=bin2dec(str4);
    end
  
%Embedding Authentication Bits into Audio signal
    G=size(hash_im,2);
    for k=1:G
        cD2(ImageSize*2+k);
        hash_im(k);
        bin=dec2bin(cD2(ImageSize*2+k),13);
        str2=extractAfter(bin,4);
        str3=dec2bin(hex2dec(hash_im(k)),4);
        str4=strcat(str3,str2);
        cD3(ImageSize*2+k)=bin2dec(str4);
    end

    G=size(hash_im,2);
    for k=1:G
        cD2(ImageSize+k);
        hash_im(k);
        bin=dec2bin(cD2(ImageSize+k),13);
        str2=extractAfter(bin,4);
        str3=dec2bin(hex2dec(hash_im(k)),4);
        str4=strcat(str3,str2);
        cD3(ImageSize+k)=bin2dec(str4);
    end
   
    for j=ImageSize+k+1:LengthDetais
        cD3(j)=cD2(j);
    end
else
%Error message    
    fprintf ('Length of details coefficient shoud be greater than image size');
    return;
end

%Reverse back to the original form of Detail Coefficient
cD3=reshape(cD3,[length(cD3),1]);
cD4=cD3/(10^4);
cD5=(cD4-P);
cD5=reshape(cD5,[length(cD5),1]);
%Inverse DWT transform to constract the Audio signal embedded with Secret
%Image
reconstructed_y=idwt(cA,cD5,wavelet);
reconstructed_y=reconstructed_y-K;
sound(reconstructed_y,Fs);%
%Write new embedded signal into .wav file for export
audiowrite('share_embedded_segment.wav',reconstructed_y,Fs);
figure;
N1 = length(reconstructed_y); % sample lenth
slength = N1/Fs; % total time span of audio signal
t = linspace(0, N1/Fs, N1);
plot(t, reconstructed_y); % pplots the audio
xlabel('Time')
ylabel('Audio Signal with Embedded Shadow')
plot(reconstructed_y)%(:,1)
% Secret_Size=R*C*16;             
% Audio_Size=N*13;
% Embedding_Rate = (Secret_Size / Audio_Size )
%Embedding Rate (bits per second) = (Size of Secret Image in bits) / (Duration of Audio Signal in seconds)


