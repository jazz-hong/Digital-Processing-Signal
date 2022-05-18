clc;
close all;
clear all;

%Audio read file
filename = 'CoffinDance-sound2.wav';
[y,Fs] = audioread(filename);
info = audioinfo('CoffinDance-sound2.wav');
y = y(:,1);
figure; %Figure 1
plot(y,'k');
title('Original Audio Signal')
xlabel('Time');
ylabel('Amplitude of Audio Signal');

%Data and Information
N_range = 100;
N = -N_range:N_range;
L = 2;
M = 2;

%==============================BUTTERWORTH================================%
Fc_Up = pi./L; %Set at default Cut off frequency
Fc_Down = pi./M;
N_New = 2; %Order of Filter

%UpSampling
x = zeros(1,L*length(y));
x([1:L:length(x)])=y;

%Butterworth filter for UpSampling
[b,a] = butter(N_New,(Fc_Up./pi)); %0 to Pi,but we need the range from 0 to 1
%So need to divide the pi to remove it
filteredSignal = L*filter(b,a,x);
figure; %Figure 2, X-axis is higher range because upSampling double the
%signals due to L=2 as Gain
plot(filteredSignal,'b');
title('Butterworth Filter-UpSampling');
xlabel('Time');
ylabel('Amplitude of Audio Signal')

%Butterworth Filter for DownSampling
[b,a] = butter(N_New,(Fc_Down./pi)); %Returns the transfer function coefficients
%of an nth-order Butterworth filter with normalized cutoff frequency
%freqz(b,a)
filteredSignal = filter(b,a,filteredSignal);
BW_DownSampling = filteredSignal(1:M:end);
if N_New==2 %Here is just to show there is no delay
   Delay = 1;
   BW_DownSampling = BW_DownSampling(1+Delay:end);
end
figure; %Figure 3
hold on;
plot(y);
plot(BW_DownSampling,'b');
hold off;
title('Butterworth Filter-DownSampling');
xlabel('Time');
ylabel('Amplitude of Audio Signal');

%=========================BUTTERWORTH SNR LOSS============================%
%SNR with different N
%Signal Noise Ratio (SNR) Loss Calculation
if isrow(BW_DownSampling) %Convert to row because SNR needs to have same position 
    BW_DownSampling = BW_DownSampling';  %Transpose the Matrix array
end
xk = rms(y); %Root Mean Square the output
ek = y-[BW_DownSampling;zeros(Delay,1)];
ek = rms(ek);
snr = 20*log(xk/ek)

%SNR based on different N and L=M
N_vary_BW = [2 5 8 10 15];
A_SNR = [];   %Empty matrix for storing data at final stage
for i = N_vary_BW     %Keep looping the N from 2 until reaches 15
    SNR = []; %Matrix storage for data then followed by variable L condition

    for L = 2:1:7  %Keep lopping the L from 2 until reaches 7 (6 data)
    M = L;         %M equals L
    Fc_Up = pi./L; %Set at default Cut off frequency
    Fc_Down = pi./M;
        if i==2
           Delay = 1;
        elseif i==5
           Delay = 2;
        elseif i==8
           Delay = 3;
        elseif i==10
           Delay = 4;
        elseif i==15
           Delay = 5;
        end

    %Repeat process of UpSampling again
    x = zeros(1,L*length(y));
    x([1:L:length(x)])=y;

    %Butterworth filter for UpSampling
    [b,a] = butter(i,(Fc_Up./pi)); %0 to Pi,but we need the range from 0 to 1
    %So need to divide the pi to remove it
    filteredSignal = L*filter(b,a,x);

    %Repeat process of DownSampling again
    [b,a] = butter(i,(Fc_Down./pi)); %Returns the transfer function coefficients
    %of an nth-order Butterworth filter with normalized cutoff frequency
    filteredSignal = filter(b,a,filteredSignal);
    BW_DownSampling = filteredSignal(1:M:end);
    BW_DownSampling = BW_DownSampling(1+Delay:end);

        if isrow(BW_DownSampling) %Firstly convert back them into all row format first
            BW_DownSampling = BW_DownSampling'; %Transpose same like previously
        end

%Copy back the SNR calculation steps again
xk = rms(y); %Root Mean Square the output
ek = y-[BW_DownSampling;zeros(Delay,1)];
ek = rms(ek);
snr = 20*log(xk/ek);
SNR = [SNR snr];
    end

%The final Matrix Storage
A_SNR = [A_SNR;SNR]; %Data stored previusly will throw into here, then RESET
                   %and looping back the process from top again
end

%After done, need to convert the Array into Table format
fprintf('\nTABULAR FORMAT for the Two Variables N and L/M\n');
Table_of_Data = array2table(A_SNR,'RowNames',{'When N=2','When N=5','When N=8','When N=10','When N=15'},'VariableNames',{'L=M=2','L=M=3','L=M=4','L=M=5','L=M=6','L=M=7',})
m = [2 3 4 5 6 7]; %6 data
figure; %Figure 4
mesh(m,N_vary_BW,A_SNR); %mesh( X,Y,Z ) creates a mesh plot with 3D surface that
%has solid edge colors and no face colors
title('SNR Butterworth Filter in Mesh');
