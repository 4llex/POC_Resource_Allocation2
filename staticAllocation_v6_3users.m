%%% Simulação de alocação estatica de usuarios em simbolo OFDM
%%% OFDMA with static allocation

%%
N = 132;                                    %% Number of subcarriers
TargetSer = 1e-3;                           %% SER Alvo
SNR = 3:3:21;
b = zeros(1,N);                             %% Vetor de Bits das portadoras / Numerologia 3
Total_bits = zeros(1,length(SNR));          %% Total de bits em um simbolo
bits_per_rb = zeros(1,length(SNR));         %% qtd media de Bits por subcarrier 
nusers = 3;

%% SNR gap para constelação M-QAM:
Gamma=(1/3)*qfuncinv(TargetSer/4)^2; % Gap to channel capacity M-QAM


%% 
%subPower = 20/1854; % 20 seria a potencia max do sistema de transmissao
% LTE EVA CHANNEL
freq_sample = 23.76e6;     %N*15e3; %30.72e6; sample rate do LTE
EVA_SR3072_Delay           = [0 30 150 310 370 710 1090 1730 2510].*1e-9;
EVA_SR3072_PowerdB_Gain    = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]; %  20*log10(0.39)= -8.1787 => -8.1787 dB -> Voltage-ratio = 0.398107

chan_EVA = rayleighchan((1/(freq_sample)),0,EVA_SR3072_Delay,EVA_SR3072_PowerdB_Gain);        
impulse= [1; zeros(N - 1,1)];  

mask1 = [ones(1, 44), zeros(1, 132-44)];
mask2 = circshift(mask1, 44);
mask3 = circshift(mask2, 44);

H    = ones(nusers,N);

num_itr = 3000;
bits_occur = zeros(length(SNR), 5); % Qtd. de bits possiveis apos quantização (0,2,4,6,8)
alloc = zeros(length(SNR),num_itr*nusers); %Para cada SNR e cada iteração, guarda  se a alocação de cada usuario
                                           %guarda a quantidade total bits alocada pra cada usuario em cada iteração 
for i=1:length(SNR)
    i
    j=0;
    while j<num_itr 
        
        % Gera o canal randomico para cada user
        for user=1:nusers
            h = filter(chan_EVA, impulse)';
            Hf = fft(h,N);
            H(user,:) = Hf;
        end
        
       
        SNRLIN = 10^(SNR(i)/10);
        P  = 20;
        Pu = P/3;
        
        [subPower1,~, subCapacity1 ] = fcn_waterfilling(Pu, P/(SNRLIN*N), Gamma, H(1,:), mask1);
        [subPower2,~, subCapacity2 ] = fcn_waterfilling(Pu, P/(SNRLIN*N), Gamma, H(2,:), mask2);
        [subPower3,~, subCapacity3 ] = fcn_waterfilling(Pu, P/(SNRLIN*N), Gamma, H(3,:), mask3);
    
        b = subCapacity1 + subCapacity2 +  subCapacity3;

        % Quantização
        b(b<2) = 0;
        b((b>2)&(b<4)) = 2;
        b((b>4)&(b<6)) = 4;
        b((b>6)&(b<8)) = 6;
        b(b>8) = 8;
        
        % Dados de alocação por usuario no Simbolo OFDM para a PDF
        alloc(i,j*3+1) = sum(b.*mask1);
        alloc(i,j*3+2) = sum(b.*mask2);
        alloc(i,j*3+3) = sum(b.*mask3);
        
        % Construir estrutura de dados com media de ocorrencias dos bits(0,2,4,6,8) no espectro.
        bit = 0;
        for nbit=1:5
            bits_occur(i,nbit) = sum(b==bit) + bits_occur(i,nbit);
            bit = bit + 2;
        end
 
        
        Total_bits(i) = Total_bits(i) + sum(b);        
        j = j+1;
    end  
    
    
    %Calcula a media de ocerrencia dos bits(0,2,4,6,8)
    for nbit=1:5
        bits_occur(i,nbit) = floor(bits_occur(i,nbit)/num_itr);
    end
    
    Total_bits(i) = Total_bits(i)/num_itr;
    bits_per_rb(i) = (Total_bits(i)/N); 
end

%% Saving Vector in a File
Static.DataSNR = SNR;   
Static.DataBPRB = bits_per_rb;
Static.DataPDF = alloc; 
Static.DataBit = bits_occur;
FileName = strcat('C:\Users\alexrosa\Documents\MATLAB\POC_Resource_Allocation2\static.mat'); 
save(FileName,'Static');

%% Gera graficos de Bits/SNR
figure;
plot(SNR, bits_per_rb, '-.r*');
title('Alocação Estática em sistema de multiplo acesso Ortogonal');
xlabel('SNR [dB]'); 
ylabel('Bits/Subcarrier'); 
grid on;
grid minor;
