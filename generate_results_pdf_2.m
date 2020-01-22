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

nbins = 40;% 15ok 20ok
for i=1:length(SNR)
    %% test - 3 histograma coloridos
                                                      %x  %y
    figure('Renderer', 'painters', 'Position', [10 550 520 390])
    
    h1 = histogram(staticAlloc(i,:), nbins);
    h1.FaceColor = 'b';
    h1.FaceAlpha = 0.5;
    hold on

    h2 = histogram(dynamicAlloc(i,:), nbins);
    hold on

    h3 = histogram(dynamicAWM_alloc(i,:), nbins);
    h3.FaceAlpha = 0.5;
    h3.FaceColor = 'y';

    legend('Estática','Dinamica Max. Vazão','Dinamica AWM/MOM');
    xlabel('bits');
    ylabel('Frequencia');
    
    str = sprintf("Histograma SNR: " + SNR(i) );
    title(str)
    
    grid on;
    grid minor;
end






