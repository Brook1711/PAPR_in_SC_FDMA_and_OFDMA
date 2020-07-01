% 该脚本实现对OFDMA和SC-FDMA的瞬时功率CDF图绘制，以及其峰均比计算
% 通用参数定义
clear all;
Fs = 5e6;                           % 系统的采样率，同时决定系统带宽
Ts = 1/Fs;                          % 系统采样时间
num_run = 2e4;                      % 运行次数

%脉冲成形参数
pulse_shaping = 1;                  % 如果是一则确定使用脉冲整形滤波器
filter_type = 'rc';                 % 选择脉冲整形滤波器的种类，这里选择升余弦滚降‘rc’；也可以选择根升余弦‘rr’
N_overlapping = 4;                  % 由于采用了脉冲成形滤波器，需要在成形是考虑码速不同
rolloff_factor = 0.09;      % 这里是滚将系数选择
                                    % 要防止被零除，请使用0.099999999而不是0.1。
        
%modulate_type = 'QPSK';            % 每个支路的基带调制方式，这里可以取‘QPSK’或‘16QAM’
number_of_subcarriers = 256;        % 这里定义了自由空间内所有可用的子信道个数，这里取为论文中的256（典型值）
number_of_symbols = 8;              % 定义每个用户使用几个子载波
Q = number_of_subcarriers/number_of_symbols;    % 扩频系数，即可同时支持多少个用户
%Q_tilda = Q-1;                     % 定义为 DFDMA 中使用的信道分配间隔. Q_tilda < Q.
subcarrier_mapping_scheme = 'IFDMA';% 这里是SCFDMA的信道映射方式，可选‘IFDMA’或‘LFDMA’
f_carrier_center = 2e6;

%初始化SC-FDMA和 OFDMA参数
papr_OFDMA = zeros(1,num_run);      % 初始化OFDMA的PAPR结果
papr_SC_IFDMA = zeros(1,num_run);    % 初始化SC-OFDMA的PAPR结果
papr_SC_LFDMA = zeros(1,num_run);

papr_OFDMA_max = zeros(1,num_run);      % 初始化OFDMA的PAPR结果
papr_SC_IFDMA_max = zeros(1,num_run);    % 初始化SC-OFDMA的PAPR结果
papr_SC_LFDMA_max = zeros(1,num_run);

papr_OFDMA_mean = zeros(1,num_run);      % 初始化OFDMA的PAPR结果
papr_SC_IFDMA_mean = zeros(1,num_run);    % 初始化SC-OFDMA的PAPR结果
papr_SC_LFDMA_mean = zeros(1,num_run);

Pi_OFDMA = zeros(1,num_run);      % 初始化OFDMA的瞬时功率结果
Pi_SC_IFDMA = zeros(1,num_run);    % 初始化SC-OFDMA的瞬时功率结果
Pi_SC_LFDMA = zeros(1,num_run);

    % 根据选择初始化脉冲成型滤波器
    if filter_type == 'rc'               % 升余弦滤波器
        psFilter = rcPulse(Ts, N_overlapping, rolloff_factor);
    elseif filter_type == 'rr'           % 根升余弦滤波器
        psFilter = rrcPulse(Ts, N_overlapping, rolloff_factor);
    end
    
    
t_last_end = 0;
for n = 1:num_run
    %首先准备随机数据：
    tmp = round(rand(number_of_symbols,2));
    tmp = tmp*2 - 1;                        %转变为双极性不归零码
    basebond_QPSK_data = (tmp(:,1) + 1i*tmp(:,2))/sqrt(2);%经过QPSK基带调制
    
    % 进行OFDM调制：
    % 首先初始化OFDMA结果序列
    OFDMA_RF_data = zeros(number_of_subcarriers,1); 
    OFDMA_RF_data(1:number_of_symbols) = basebond_QPSK_data;
    OFDMA_RF_data = ifft(OFDMA_RF_data);
    
    % 进行单载波调制：
    % 首先初始化SCFDMA结果序列
    SCFDMA_IFDMA_RF_data = zeros(number_of_subcarriers,1); 
    SCFDMA_LFDMA_RF_data = zeros(number_of_subcarriers,1); 
    % 将基带不同支路的基带QPSK信号通过一次DFT（FFT）
    X = fft(basebond_QPSK_data);
    % 通过不同映射方式将信号映射进不同的子载波的位置上
 
    SCFDMA_IFDMA_RF_data(1:Q:number_of_subcarriers) = X;

    SCFDMA_LFDMA_RF_data(1:number_of_symbols) = X;

    SCFDMA_IFDMA_RF_data = ifft(SCFDMA_IFDMA_RF_data);%通过ifft得到信号复包络
    SCFDMA_LFDMA_RF_data = ifft(SCFDMA_LFDMA_RF_data);
    %进行脉冲整形
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
    %这一步得到的是信号的复包络，接下来进行数字上变频：
    
    t = 0 : Ts/N_overlapping : Ts/N_overlapping*(length(SCFDMA_IFDMA_RF_data)-1);
    t = t + t_last_end;
    t_last_end = t(end);
    SCFDMA_IFDMA_RF_data = SCFDMA_IFDMA_RF_data .* exp(1i*2*pi*f_carrier_center*t);
    SCFDMA_LFDMA_RF_data = SCFDMA_LFDMA_RF_data .* exp(1i*2*pi*f_carrier_center*t);
    OFDMA_RF_data = OFDMA_RF_data .* exp(1i*2*pi*f_carrier_center*t);

    %计算这一帧时刻的峰均值
    papr_SC_IFDMA(n) = 10*log10(max(abs(SCFDMA_IFDMA_RF_data).^2) / mean(abs(SCFDMA_IFDMA_RF_data).^2));
    papr_SC_LFDMA(n) = 10*log10(max(abs(SCFDMA_LFDMA_RF_data).^2) / mean(abs(SCFDMA_LFDMA_RF_data).^2));
    papr_OFDMA(n) = 10*log10(max(abs(OFDMA_RF_data).^2) / mean(abs(OFDMA_RF_data).^2));
    
    papr_SC_IFDMA_max(n) = max(abs(SCFDMA_IFDMA_RF_data).^2);
    papr_SC_LFDMA_max(n) = max(abs(SCFDMA_LFDMA_RF_data).^2);
    papr_OFDMA_max(n) = max(abs(OFDMA_RF_data).^2);
    
    papr_SC_IFDMA_mean(n) = mean(abs(SCFDMA_IFDMA_RF_data).^2);
    papr_SC_LFDMA_mean(n) = mean(abs(SCFDMA_LFDMA_RF_data).^2);
    papr_OFDMA_mean(n) = mean(abs(OFDMA_RF_data).^2);
    %计算这一时刻的瞬时功率单位为dBm
    index_temp = round(mod(n,length(SCFDMA_IFDMA_RF_data)));
    if index_temp ==0
        index_temp =1;
    end
%     Pi_SC_IFDMA(n) = 10*log10(0.5*abs(SCFDMA_IFDMA_RF_data(index_temp))^2/0.001);
%     Pi_SC_LFDMA(n) = 10*log10(0.5*abs(SCFDMA_LFDMA_RF_data(index_temp))^2/0.001);
%     Pi_OFDMA(n) = 10*log10(0.5*abs(OFDMA_RF_data(index_temp))^2/0.001);
     %计算这一时刻的瞬时功率单位为W
    Pi_SC_IFDMA(n) = 0.5*abs(SCFDMA_IFDMA_RF_data(index_temp))^2;
    Pi_SC_LFDMA(n) = 0.5*abs(SCFDMA_LFDMA_RF_data(index_temp))^2;
    Pi_OFDMA(n) =8 * 0.5*abs(OFDMA_RF_data(index_temp))^2;

end
%计算总的峰均比
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
xlabel("峰均比，单位：dB");
ylabel("峰均比互补累计分布函数");

figure(2);
plot(X,cumsum(N)/max(cumsum(N)), X2,cumsum(N2)/max(cumsum(N2)),'--mo',X3,cumsum(N3)/max(cumsum(N3)),'-.r*');
legend('OFDMA','LFDMA','IFDMA');
xlabel("峰均比，单位：dB");
ylabel("峰均比累计分布");

[N4,X4] = hist(Pi_OFDMA, 20);
[N5,X5] = hist(Pi_SC_LFDMA, 20);
[N6,X6] = hist(Pi_SC_IFDMA, 20);
figure(3);
plot(X4,cumsum(N4)/max(cumsum(N4)), X5,cumsum(N5)/max(cumsum(N5)),'--mo',X6,cumsum(N6)/max(cumsum(N6)),'-.r*')
title('Instantaneous power CDF of SC-FDMA and OFDM');
legend('OFDMA','LFDMA','IFDMA');
xlabel("瞬时功率，单位：W");
ylabel("瞬时功率累积分布");
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