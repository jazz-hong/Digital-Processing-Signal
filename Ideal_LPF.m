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

%===============================Ideal LPF=================================%
%UpSampling
x = zeros(1,L*length(y));
x([1:L:length(x)])=y;
figure; %Figure 2
subplot (2,1,1); %Graph position
plot(y(1:length(y)));                       
title('Input signal');
subplot (2,1,2);
plot(x);                       
title('Upsampling');

%Convolution LPF for UpSamling
omegac = pi/L;
LPF = L.*sin(omegac*N)./(pi*N);
LPF(N==0) = L*omegac/pi;
figure; %Figure 3
stem(N,LPF);
title ('Low Pass Filter using Convolution-UpSampling');
xlabel('Number of order(th)');
ylabel('Amplitude of Audio Signal')
xn = conv(LPF,x);
xn = xn(N_range+1:end); %Shifting back to avoid delay
figure; %Figure 4
plot(x,'r');
hold on;
stem(xn,'k');
title ('Ideal LPF-UpSampling');
xlabel('Time');
ylabel('Amplitude of Audio Signal')

%Anti-Aliasing LPF with Convolution
omegac = pi/M;
LPF = sin(omegac*N)./(pi*N);
LPF(N==0)=omegac/pi;
figure; %Figure 5
stem(N,LPF);
title ('Anti-Aliasing LPF with Convolution before DownSampling');
xlabel('Number of order(th)');
ylabel('Amplitude of Audio Signal')
xn_aa = conv(LPF,xn);
xn_aa = xn_aa(N_range+1:end); %Shifting back to avoid delay
figure; %Figure 6
plot(x,'r');
hold on;
plot(xn_aa,'k');
title ('After Anti-Aliasing LPF preparing for DownSampling');
xlabel('Time');
ylabel('Amplitude of Audio Signal')

%DownSampling
y1 = xn_aa([1:M:length(xn_aa)]);
figure; %Figure 7
plot(y);                       
hold on;
stem(y1,'k')
title('Ideal LPF-DownSampling');
xlabel('Time');
ylabel('Amplitude of Audio Signal')

%=============================IDEAL SNR LOSS===============================%
%Signal Noise Ratio (SNR) Loss Calculation
if isrow(y1) %Convert to row because SNR needs to have same position 
    y1 = y1';  %Transpose the Matrix array
end
xk = rms(y); %Root Mean Square the output
ek = y-y1(1:length(y));
ek = rms(ek);
snr = 20*log(xk/ek)

%SNR based on different N and L=M
N_vary = [50 75 100 150 200];
A_SNR = [];   %Empty matrix for storing data at final stage
for i = N_vary     %Keep looping the N from 50 until reaches 200
    SNR = []; %Matrix storage for data then followed by variable L condition

    for L = 2:1:7  %Keep lopping the L from 2 until reaches 7 (6 data)
    M = L;         %M equals L

    %Repeat process of UpSampling again
    x = zeros(1,L*length(y)); %But need to insert zero in each sample again
    x([1:L:length(x)])=y;
    omegac = pi/L;
    N = -i:i; %Replace the N into "i" because we have changed the variable in SNR
    LPF = L.*sin(omegac*N)./(pi*N);
    LPF(N==0) = L*omegac/pi;
    xn = conv(LPF,x);
    xn = xn(i+1:end); %Shifting back to avoid delay

    %Repeat process of DownSampling again
    omegac = pi/M;
    LPF = sin(omegac*N)./(pi*N);
    LPF(N==0) = omegac/pi;
    xn_aa = conv(LPF,xn);
    xn_aa = xn_aa(i+1:end);
    y1 = xn_aa([1:M:length(xn_aa)]);

        if isrow(y1) %Firstly convert back them into all row format first
            y1 = y1'; %Transpose same like previously
        end

%Copy back the SNR calculation steps again
xk = rms(y); %Root Mean Square the output
ek = y-y1(1:length(y));
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
Table_of_Data = array2table(A_SNR,'RowNames',{'When N=50','When N=75','When N=100','When N=150','When N=200'},'VariableNames',{'L=M=2','L=M=3','L=M=4','L=M=5','L=M=6','L=M=7',})
m = [2 3 4 5 6 7]; %6 data
figure; %Figure 8
mesh(m,N_vary,A_SNR); %mesh( X,Y,Z ) creates a mesh plot with 3D surface that
%has solid edge colors and no face colors
title('SNR Ideal LPF in Mesh');
