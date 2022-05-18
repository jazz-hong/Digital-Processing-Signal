function target_Modification
  clc;
  clear all;
  close all;
 
  Input_Data = load('ak_Input.txt','B'); %Load input data
  Output_Data = load('rk_Output.txt','B'); %Load output data
  
  Nw = 12;       % Length of channel
  plot_SNR = [];
  
  %gk are the target coefficients
  gk=[1 2 1 ]'; %Short target
  Ng = length(gk);
  go = [1:0.01:3];
  Range_go = 1:0.01:3;

  for x = 1:0.01:3
      go = x;
      gk = [1 go 1]';
      R=calcAutocorrMatrix(Nw,Output_Data); % Calc Nw x Nw auto-corr matrix of rk   
      P=calcCrossCorrMatrix(Ng,Nw,Input_Data,Output_Data); % Calc cross-correlation vector
      wk = inv(R)*P'*gk;
      zk = conv(Output_Data,wk);
      dk = conv(gk,Input_Data);

      %Alignment
      dk = dk(1+floor(Ng/2):end-floor(Ng/2)); 
      zk=zk(1+11:length(dk));  

      %SNR calculation
      SNR_dk = rms(dk);
      SNR_ek_otw = dk(1:length(zk))-zk;
      SNR_ek = rms(SNR_ek_otw);
      SNR = 20*log10(SNR_dk/SNR_ek);
      plot_SNR = [plot_SNR,SNR];
   end
   figure;
   plot(Range_go, plot_SNR,'k-x'); %Need to get a peak
   title('SNR vs go');
   xlabel('go');
   ylabel('Amplitude of SNR');
end

%%--------------------------FUNCTION DEFINITION--------------------------%%
%Function to calculate the cross-correlation vector
function P=calcCrossCorrMatrix(Ng,Nw,Input_Data,Output_Data)
  Input_Data=Input_Data(:);
  Output_Data=Output_Data(:);
  rowCorr = calcCrossCorrVect(Input_Data,Output_Data,Nw);
  colCorr = calcCrossCorrVect(Output_Data,Input_Data,Ng)';
  P=toeplitz(colCorr,rowCorr);
end

function p = calcCrossCorrVect(Input_Data,Output_Data,Nw)
  N=min(length(Input_Data),length(Output_Data));
  for i=1:Nw;
    p(i)=Input_Data(1:N-i+1)'*Output_Data(i:N)/(N-i+1);
  end
end

% Function to calculate the auto-correlation matrix
function R=calcAutocorrMatrix(Nh,ak)
  ak=ak(:);
  Nc=length(ak);
  for i=1:Nh
    r(i) = ak(i:end)'*ak(1:end-i+1)/(Nc-i+1);
  end
  R=toeplitz(r);
end
