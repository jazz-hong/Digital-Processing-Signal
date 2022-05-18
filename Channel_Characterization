function channel_Characterization
    clc;
    clear all;
    close all;

    Input_Data = load('ak_Input.txt','B'); %Load input data
    Output_Data = load('rk_Output.txt','B'); %Load output data
    
    Nh = 2:1:20; %Length of channel
    plot_SNR = [];

    for x = 1:length(Nh)
        %Compute the auto-correlation matrix
        R = calcAutocorrMatrix(Nh(x),Input_Data); %Calc Nh x Nh auto-corr matrix
        p = calcCrossCorr(Nh(x),Input_Data,Output_Data); %Calc Nh x 1 cross-correlation vector
        Channel_Coeff = inv(R)*p;
        sk = conv(Channel_Coeff,Input_Data);

        %SNR calculation
        SNR_sk = rms(sk);
        SNR_ek_otw = sk(1:length(Output_Data))-Output_Data;
        SNR_ek = rms(SNR_ek_otw);
        SNR = 20*log10(SNR_sk/SNR_ek);
        plot_SNR(x) = SNR; %Need same dimension
    end
    %SNR calculation
    SNR_sk = rms(sk);
    SNR_ek_otw = sk(1:length(Output_Data))-Output_Data;
    SNR_ek = rms(SNR_ek_otw);
    SNR = 20*log10(SNR_sk/SNR_ek)
    figure;
    plot(Output_Data,'k-x');
    hold on;
    plot(sk(1+floor(Nh/2):end),'r-x');
    title('Output signals rk and sk');
    xlabel('Time');
    ylabel('Amplitude');
    figure;
    plot(Nh,plot_SNR,'k-o');
    title('SNR vs Nh');
    xlabel('Nh');
    ylabel('Amplitude of SNR');
end

%%--------------------------FUNCTION DEFINITION--------------------------%%
%Function to calculate the cross-correlation vector
function p = calcCrossCorr(Nh,Input_Data,Output_Data)
  Input_Data = Input_Data(:);
  Output_Data = Output_Data(:);
  N = length(Input_Data);
  for i = 1:Nh
    p(i) = sum(Input_Data(1:N-i+1)'*Output_Data(i:N)/(N-i+1));
  end
  p = p(:);
end

function R = calcAutocorrMatrix(Nh,Input_Data)
    Input_Data = Input_Data(:);
    Na = length(Input_Data);
        for i=1:Nh
            r(i) = Input_Data(i:end)'*Input_Data(1:end-i+1)/(Na-i+1);
        end
            R = toeplitz(r);
end
