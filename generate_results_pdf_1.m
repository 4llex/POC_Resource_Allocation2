%%% Script para gerar os resultados dos 3 algoritmos simulados
%%% Gera os histogramas com a distribuição dos bits alocados pelos usuários

%% Loading static.mat
SimData=load('static.mat');
%staticSNR   = SimData.Static.DataSNR;
%staticBPRB  = SimData.Static.DataBPRB; 
staticAlloc = SimData.Static.DataPDF; 

%% Loading dynamicMaxVazao.mat
DynamicData=load('dynamicMaxVazao.mat');
%dynamicSNR   = DynamicData.Dynamic.DataSNR;
%dynamicBPRB  = DynamicData.Dynamic.DataBPRB;
dynamicAlloc = DynamicData.Dynamic.DataPDF;

%% Loanding dynamicAWM.mat
DynamicAWM_Data=load('dynamicAWM.mat'); 
%dynamicAWM_SNR   = DynamicAWM_Data.DynamicAWM.DataSNR;
%dynamicAWM_BPRB  = DynamicAWM_Data.DynamicAWM.DataBPRB;
dynamicAWM_alloc = DynamicAWM_Data.DynamicAWM.DataPDF;

SNR = 3:3:21; 

for i=1:length(SNR)
    %% test - 3 histograma coloridos
    figure;
    histogram(staticAlloc(i,:),'BinMethod','integers','FaceColor','r','EdgeAlpha',0,'FaceAlpha',1);
    hold on;
    histogram(dynamicAlloc(i,:),'BinMethod','integers','FaceColor','g','EdgeAlpha',0,'FaceAlpha',0.7);
    histogram(dynamicAWM_alloc(i,:),'BinMethod','integers','FaceColor','b','EdgeAlpha',0,'FaceAlpha',0.7);
    xlabel('Bits/Usuario');
    ylabel('Frequencia');
    str = sprintf("Histograma SNR: " + SNR(i) );
    title(str)
    %title('Histograma');
    
    legend('Alocação Estática','Alocação Dinamica', 'AWM/MOM')
    xlim([0 150]); 
    
    grid on;
    grid minor;

end





