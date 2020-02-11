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

staticAllocBar = zeros(7,150);
dynamicaAllocBar = zeros(7,150);
AwmMomAllocBar = zeros(7,150);

a1 = staticAlloc(6,:);
b1 = randi([1 10],[1,40]);
c1 = randi([0 10],[1,40]);

% vec = zeros(1,40); % Calcular ocorrecina de valores em bins de tamanho 10
% u = 0;
% for k=1:40
%     vec(k) = sum((dynamicAlloc(1,:)>= u)&(dynamicAlloc(1,:)< u+10));
%     u = u+10;
% end


for i=1:7
    fig=figure; set(fig,'visible','off');
    % Static allocation vec
    a = histogram(staticAlloc(i,:));
    a.BinLimits = [0 500];
    a.BinWidth = 10;
    testa = a.Values;
    j=1;
    for ii=1:50 % 40 é o qtd. de bins
        staticAllocBar(i,j) = a.Values(ii); 
        j=j+3;
    end
    hold on
    % Dynamic allocation WF vec
    b = histogram(dynamicAlloc(i,:));
    b.BinLimits = [0 500];
    b.BinWidth = 10; 
    testb = b.Values;
    j=2;
    for ii=1:50 % 40 é o qtd. de bins
        dynamicaAllocBar(i,j) = b.Values(ii);
        j=j+3;
    end
    
    % Dynamic allocation AWM_MOM vec
    c = histogram(dynamicAWM_alloc(i,:));
    c.BinLimits = [0 500];
    c.BinWidth = 10; 
    testc = c.Values;
    j=3;
    for ii=1:50 % 40 é o qtd. de bins
        AwmMomAllocBar(i,j) = c.Values(ii);
        j=j+3;
    end
    
    close(fig);
end

%% Wheberth suggestion %works
for i=1:7
    x = 0:3.34:497.66;
    figure;
    bar(x,staticAllocBar(i,:).',1.0 ,'r');
    hold on
    bar(x,dynamicaAllocBar(i,:).',1.0 ,'FaceColor','blue');
    bar(x,AwmMomAllocBar(i,:).',1.0 ,'y' )
    xlim([-2 500])
    legend('Statica','Dinamica Max. Vazao','AWM/MOM');
    xlabel('QoS')
    ylabel('Ocorrencia')
    title('Histograma');
    grid on;
    grid minor;
end

%% Other solution for the same ideia, utilizando retorno do metodo histogram
% i=7;
% x= 0:10:399;
% vals = [testa.'  testb.'  testc.'];
% figure('Renderer', 'painters', 'Position', [10 550 520 390]);
% b = bar(x, vals, 'hist');
% xlim([-3 400])










