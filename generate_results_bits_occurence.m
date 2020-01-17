%%% Script para gerar os resultados dos 3 algoritmos simulados
%%% Gera os Qtd. de ocorrencia dos bits 0,2,4,6,8 no espectro

%% Loading static.mat
SimData=load('static.mat');
%staticSNR   = SimData.Static.DataSNR;
%staticBPRB  = SimData.Static.DataBPRB; 
%staticAlloc = SimData.Static.DataPDF; 
staticBit   = SimData.Static.DataBit;

%% Loading dynamicMaxVazao.mat
DynamicData=load('dynamicMaxVazao.mat');
%dynamicSNR   = DynamicData.Dynamic.DataSNR;
%dynamicBPRB  = DynamicData.Dynamic.DataBPRB;
%dynamicAlloc = DynamicData.Dynamic.DataPDF;
dynamicBit   = DynamicData.Dynamic.DataBit;

%% Loanding dynamicAWM.mat
DynamicAWM_Data=load('dynamicAWM.mat'); 
%dynamicAWM_SNR   = DynamicAWM_Data.DynamicAWM.DataSNR;
%dynamicAWM_BPRB  = DynamicAWM_Data.DynamicAWM.DataBPRB;
%dynamicAWM_alloc = DynamicAWM_Data.DynamicAWM.DataPDF;
dynamicAWM_bit   = DynamicAWM_Data.DynamicAWM.DataBit;


SNR = 3:3:21; 
bits = [0 2 4 6 8];

for i=1:length(SNR)
    
    vals = [staticBit(i,:).'  dynamicBit(i,:).'  dynamicAWM_bit(i,:).'];
                                                       %x  %y
    figure('Renderer', 'painters', 'Position', [10 550 520 390]);
    
    xlabel('SNR [dB]'); 
    ylabel('Bits/RB'); 
    b = bar(bits,vals, 'hist');
    legend('Estática','Dinâmica Max. Vazão','Dinâmica AWM/MOM');
    str = sprintf("Histograma SNR: " + SNR(i) );
    title(str)
    xlabel('Bits'); 
    ylabel('Qtd. de ocorrência dos Bits'); 
    grid on;
    grid minor;
    
end




