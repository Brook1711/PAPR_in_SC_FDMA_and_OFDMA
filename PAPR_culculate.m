% �ýű�ʵ�ֶ�OFDMA��SC-FDMA��˲ʱ����CDFͼ���ƣ��Լ������ȼ���
% ͨ�ò�������
clear all;
Fs = 5e6;                           % ϵͳ�Ĳ����ʣ�ͬʱ����ϵͳ����
Ts = 1/Fs;                          % ϵͳ����ʱ��
num_run = 2e4;                      % ���д���

%������β���
pulse_shaping = 1;                  % �����һ��ȷ��ʹ�����������˲���
filter_type = 'rc';                 % ѡ�����������˲��������࣬����ѡ�������ҹ�����rc����Ҳ����ѡ��������ҡ�rr��
N_overlapping = 4;                  % ���ڲ�������������˲�������Ҫ�ڳ����ǿ������ٲ�ͬ
rolloff_factor = 0.09;      % �����ǹ���ϵ��ѡ��
                                    % Ҫ��ֹ���������ʹ��0.099999999������0.1��
        
%modulate_type = 'QPSK';            % ÿ��֧·�Ļ������Ʒ�ʽ���������ȡ��QPSK����16QAM��
number_of_subcarriers = 256;        % ���ﶨ�������ɿռ������п��õ����ŵ�����������ȡΪ�����е�256������ֵ��
number_of_symbols = 8;              % ����ÿ���û�ʹ�ü������ز�
Q = number_of_subcarriers/number_of_symbols;    % ��Ƶϵ��������ͬʱ֧�ֶ��ٸ��û�
%Q_tilda = Q-1;                     % ����Ϊ DFDMA ��ʹ�õ��ŵ�������. Q_tilda < Q.
subcarrier_mapping_scheme = 'IFDMA';% ������SCFDMA���ŵ�ӳ�䷽ʽ����ѡ��IFDMA����LFDMA��
f_carrier_center = 2e6;

%��ʼ��SC-FDMA�� OFDMA����
papr_OFDMA = zeros(1,num_run);      % ��ʼ��OFDMA��PAPR���
papr_SC_IFDMA = zeros(1,num_run);    % ��ʼ��SC-OFDMA��PAPR���
papr_SC_LFDMA = zeros(1,num_run);

papr_OFDMA_max = zeros(1,num_run);      % ��ʼ��OFDMA��PAPR���
papr_SC_IFDMA_max = zeros(1,num_run);    % ��ʼ��SC-OFDMA��PAPR���
papr_SC_LFDMA_max = zeros(1,num_run);

papr_OFDMA_mean = zeros(1,num_run);      % ��ʼ��OFDMA��PAPR���
papr_SC_IFDMA_mean = zeros(1,num_run);    % ��ʼ��SC-OFDMA��PAPR���
papr_SC_LFDMA_mean = zeros(1,num_run);

Pi_OFDMA = zeros(1,num_run);      % ��ʼ��OFDMA��˲ʱ���ʽ��
Pi_SC_IFDMA = zeros(1,num_run);    % ��ʼ��SC-OFDMA��˲ʱ���ʽ��
Pi_SC_LFDMA = zeros(1,num_run);

    % ����ѡ���ʼ����������˲���
    if filter_type == 'rc'               % �������˲���
        psFilter = rcPulse(Ts, N_overlapping, rolloff_factor);
    elseif filter_type == 'rr'           % ���������˲���
        psFilter = rrcPulse(Ts, N_overlapping, rolloff_factor);
    end
    
    
t_last_end = 0;
for n = 1:num_run
    %����׼��������ݣ�
    tmp = round(rand(number_of_symbols,2));
    tmp = tmp*2 - 1;                        %ת��Ϊ˫���Բ�������
    basebond_QPSK_data = (tmp(:,1) + 1i*tmp(:,2))/sqrt(2);%����QPSK��������
    
    % ����OFDM���ƣ�
    % ���ȳ�ʼ��OFDMA�������
    OFDMA_RF_data = zeros(number_of_subcarriers,1); 
    OFDMA_RF_data(1:number_of_symbols) = basebond_QPSK_data;
    OFDMA_RF_data = ifft(OFDMA_RF_data);
    
    % ���е��ز����ƣ�
    % ���ȳ�ʼ��SCFDMA�������
    SCFDMA_IFDMA_RF_data = zeros(number_of_subcarriers,1); 
    SCFDMA_LFDMA_RF_data = zeros(number_of_subcarriers,1); 
    % ��������֧ͬ·�Ļ���QPSK�ź�ͨ��һ��DFT��FFT��
    X = fft(basebond_QPSK_data);
    % ͨ����ͬӳ�䷽ʽ���ź�ӳ�����ͬ�����ز���λ����
 
    SCFDMA_IFDMA_RF_data(1:Q:number_of_subcarriers) = X;

    SCFDMA_LFDMA_RF_data(1:number_of_symbols) = X;

    SCFDMA_IFDMA_RF_data = ifft(SCFDMA_IFDMA_RF_data);%ͨ��ifft�õ��źŸ�����
    SCFDMA_LFDMA_RF_data = ifft(SCFDMA_LFDMA_RF_data);
    %������������
    if pulse_shaping == 1
        % Up-sample the symbols.
        IFDMA_oversampled(1:N_overlapping:N_overlapping*number_of_subcarriers) = SCFDMA_IFDMA_RF_data;
        LFDMA_oversampled(1:N_overlapping:N_overlapping*number_of_subcarriers) = SCFDMA_LFDMA_RF_data;
        OFDM_oversampled(1:N_overlapping:N_overlapping*number_of_subcarriers) = OFDMA_RF_data;
        % Perform filtering.
        SCFDMA_IFDMA_RF_data = filter(psFilter, 1, IFDMA_oversampled);
        SCFDMA_LFDMA_RF_data = filter(psFilter, 1, LFDMA_oversampled);
        OFDMA_RF_data = filter(psFilter, 1, OFDM_oversampled);
    else
        
    end
    %��һ���õ������źŵĸ����磬���������������ϱ�Ƶ��
    
    t = 0 : Ts/N_overlapping : Ts/N_overlapping*(length(SCFDMA_IFDMA_RF_data)-1);
    t = t + t_last_end;
    t_last_end = t(end);
    SCFDMA_IFDMA_RF_data = SCFDMA_IFDMA_RF_data .* exp(1i*2*pi*f_carrier_center*t);
    SCFDMA_LFDMA_RF_data = SCFDMA_LFDMA_RF_data .* exp(1i*2*pi*f_carrier_center*t);
    OFDMA_RF_data = OFDMA_RF_data .* exp(1i*2*pi*f_carrier_center*t);

    %������һ֡ʱ�̵ķ��ֵ
    papr_SC_IFDMA(n) = 10*log10(max(abs(SCFDMA_IFDMA_RF_data).^2) / mean(abs(SCFDMA_IFDMA_RF_data).^2));
    papr_SC_LFDMA(n) = 10*log10(max(abs(SCFDMA_LFDMA_RF_data).^2) / mean(abs(SCFDMA_LFDMA_RF_data).^2));
    papr_OFDMA(n) = 10*log10(max(abs(OFDMA_RF_data).^2) / mean(abs(OFDMA_RF_data).^2));
    
    papr_SC_IFDMA_max(n) = max(abs(SCFDMA_IFDMA_RF_data).^2);
    papr_SC_LFDMA_max(n) = max(abs(SCFDMA_LFDMA_RF_data).^2);
    papr_OFDMA_max(n) = max(abs(OFDMA_RF_data).^2);
    
    papr_SC_IFDMA_mean(n) = mean(abs(SCFDMA_IFDMA_RF_data).^2);
    papr_SC_LFDMA_mean(n) = mean(abs(SCFDMA_LFDMA_RF_data).^2);
    papr_OFDMA_mean(n) = mean(abs(OFDMA_RF_data).^2);
    %������һʱ�̵�˲ʱ���ʵ�λΪdBm
    index_temp = round(mod(n,length(SCFDMA_IFDMA_RF_data)));
    if index_temp ==0
        index_temp =1;
    end
%     Pi_SC_IFDMA(n) = 10*log10(0.5*abs(SCFDMA_IFDMA_RF_data(index_temp))^2/0.001);
%     Pi_SC_LFDMA(n) = 10*log10(0.5*abs(SCFDMA_LFDMA_RF_data(index_temp))^2/0.001);
%     Pi_OFDMA(n) = 10*log10(0.5*abs(OFDMA_RF_data(index_temp))^2/0.001);
     %������һʱ�̵�˲ʱ���ʵ�λΪW
    Pi_SC_IFDMA(n) = 0.5*abs(SCFDMA_IFDMA_RF_data(index_temp))^2;
    Pi_SC_LFDMA(n) = 0.5*abs(SCFDMA_LFDMA_RF_data(index_temp))^2;
    Pi_OFDMA(n) =8 * 0.5*abs(OFDMA_RF_data(index_temp))^2;

end
%�����ܵķ����
total_papr_SC_IFDMA_max = max(papr_SC_IFDMA_max);
total_papr_SC_IFDMA_mean = mean(papr_SC_IFDMA_mean);
total_papr_SC_IFDMA = 10*log10(total_papr_SC_IFDMA_max / total_papr_SC_IFDMA_mean);

total_papr_SC_LFDMA_max = max(papr_SC_LFDMA_max);
total_papr_SC_LFDMA_mean = mean(papr_SC_LFDMA_mean);
total_papr_SC_LFDMA = 10*log10(total_papr_SC_LFDMA_max / total_papr_SC_LFDMA_mean);

total_papr_OFDMA_max = max(papr_OFDMA_max);
total_papr_OFDMA_mean = mean(papr_OFDMA_mean);
total_papr_OFDMA = 10*log10(total_papr_OFDMA_max / total_papr_OFDMA_mean);

[N,X] = hist(papr_OFDMA, 10);
[N2,X2] = hist(papr_SC_LFDMA, 100);
[N3,X3] = hist(papr_SC_IFDMA, 100);
figure(1);
semilogy(X,1-cumsum(N)/max(cumsum(N)), X2,1-cumsum(N2)/max(cumsum(N2)),'--mo',X3,1-cumsum(N3)/max(cumsum(N3)),'-.r*');
title('PAPR CCDF of SC-FDMA and OFDM');
legend(strcat('OFDMA PAPR is : ', string(total_papr_OFDMA),' dB'),...
    strcat('LFDMA PAPR is : ', string(total_papr_SC_LFDMA),' dB'),...
    strcat('IFDMA PAPR is : ', string(total_papr_SC_IFDMA),' dB'));
xlabel("����ȣ���λ��dB");
ylabel("����Ȼ����ۼƷֲ�����");

figure(2);
plot(X,cumsum(N)/max(cumsum(N)), X2,cumsum(N2)/max(cumsum(N2)),'--mo',X3,cumsum(N3)/max(cumsum(N3)),'-.r*');
legend('OFDMA','LFDMA','IFDMA');
xlabel("����ȣ���λ��dB");
ylabel("������ۼƷֲ�");

[N4,X4] = hist(Pi_OFDMA, 20);
[N5,X5] = hist(Pi_SC_LFDMA, 20);
[N6,X6] = hist(Pi_SC_IFDMA, 20);
figure(3);
plot(X4,cumsum(N4)/max(cumsum(N4)), X5,cumsum(N5)/max(cumsum(N5)),'--mo',X6,cumsum(N6)/max(cumsum(N6)),'-.r*')
title('Instantaneous power CDF of SC-FDMA and OFDM');
legend('OFDMA','LFDMA','IFDMA');
xlabel("˲ʱ���ʣ���λ��W");
ylabel("˲ʱ�����ۻ��ֲ�");
function r = rcPulse(Ts, Nos, alpha)

t1 = [-8*Ts:Ts/Nos:-Ts/Nos];
t2 = [Ts/Nos:Ts/Nos:8*Ts];

r1 = (sin(pi*t1/Ts)./(pi*t1)).*(cos(pi*alpha*t1/Ts)./(1-(4*alpha*t1/(2*Ts)).^2));
r2 = (sin(pi*t2/Ts)./(pi*t2)).*(cos(pi*alpha*t2/Ts)./(1-(4*alpha*t2/(2*Ts)).^2));

r = [r1 1/Ts r2];
end
function r = rrcPulse(Ts, Nos, alpha)

t1 = [-6*Ts:Ts/Nos:-Ts/Nos];
t2 = [Ts/Nos:Ts/Nos:6*Ts];

r1 = (4*alpha/(pi*sqrt(Ts)))*(cos((1+alpha)*pi*t1/Ts)+(Ts./(4*alpha*t1)).*sin((1-alpha)*pi*t1/Ts))./(1-(4*alpha*t1/Ts).^2);
r2 = (4*alpha/(pi*sqrt(Ts)))*(cos((1+alpha)*pi*t2/Ts)+(Ts./(4*alpha*t2)).*sin((1-alpha)*pi*t2/Ts))./(1-(4*alpha*t2/Ts).^2);

r = [r1 (4*alpha/(pi*sqrt(Ts))+(1-alpha)/sqrt(Ts)) r2];
end