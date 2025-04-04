# ECG signal analysis

# 1. Loading data and defining parameters

clc; clear; close all;
load('100m.mat')
ECGsignal=(val-1024)/200;

fs=360;
deltat=1/fs;
t=0:0.002777777778:length(ECGsignal-1)/fs;
N=length(t);
deltaf=1/(deltat*N);


# 2. DFT calculation and periodogram

f=-0.5/deltat : deltaf : 0.5/deltat-deltaf;
DFT=0;
for n=1:N
   DFT=DFT+ECGsignal(n)*exp(-1i*2*pi*f*t(n));
end
DFT=DFT*deltat;


P=(1/max(t))*(abs(DFT).^2);

figure()
plot(t,ECGsignal)
grid minor
title('ECG signal Arithmia')
xlabel('Time [s]')
ylabel('ECG Signal [Hz]')

figure()
plot(f,abs(DFT))
grid minor
title('ECG signal Arithmia - Module')
xlabel('Frequency [Hz]')
ylabel('Module ECG Signal')
ylim([0 2])

figure()
plot(f,unwrap(angle(DFT)))
grid minor
title('ECG signal Arithmia - Shift Fase')
xlabel('Frequency [Hz]')
ylabel('Shift-fase ECG Signal')

figure()
plot(f,P)
title('Periodogram');
grid minor
xlim([-1 1])
xlabel('Frequency [Hz]')
ylabel('Power spectral density [mV²/Hz]')
ylim([0 0.08])

# 3. Low-frequency noise filtering

DFT_FILTERED=DFT;
for n=1:length(DFT)
    if le(abs(f(n)),0.5)
        DFT_FILTERED(n)=0;
    end
end

figure()
plot(f, DFT_FILTERED)
title('Low-frequency noise filtering')
xlabel('Frequency [Hz]')
ylabel('Filtered module')


index_60Hz=14401;
index_neg_60Hz=7201;
DFT_FILTERED(index_60Hz)=0;
DFT_FILTERED(index_neg_60Hz)=0;

figure()
plot(f, DFT_FILTERED)
title('Low-frequency noise filtering')
xlabel('Frequency [Hz]')
ylabel('Filtered module')



# 4. Return to the time domain

IDFT=0;
for n=1:length(DFT_FILTERED)
   IDFT=IDFT+DFT_FILTERED(n)*exp(-1i*2*pi*f(n)*t);
end

figure()
plot(f,real(IDFT))
title('Real part')
xlabel('Frequency [Hz]')
ylabel('Real part')
grid minor

figure()
plot(f,imag(IDFT))
title('Imaginary part')
xlabel('Frequency [Hz]')
ylabel('Imaginary part')
grid minor


# 5. Identification of peak instants

figure
plot(t, real(IDFT)/max(real(IDFT)))
title('Normalized real part')
xlabel('Time [s]')
ylabel('Normalized real part')

ECG_FILTERED=real(IDFT)/max(abs(real(IDFT)));
soglia=0.5;
h=1;

for n = 1:length(IDFT)-1
    if ge(ECG_FILTERED(n+1),soglia) & le(ECG_FILTERED(n),soglia)
        indice_ECG_FILTERED_up = find(ECG_FILTERED == (ECG_FILTERED(n+1)));
        t_up = t(indice_ECG_FILTERED_up);
    end
    if le(ECG_FILTERED(n+1),soglia) & ge(ECG_FILTERED(n),soglia)
        indice_ECG_FILTERED_down = find(ECG_FILTERED == (ECG_FILTERED(n)));
        t_down = t(indice_ECG_FILTERED_down);
        t_RR(h) = (t_up+t_down)/2;
        h = h+1;
    end
end
for n=1:length(t_RR)-1
    diff_RR(n)=t_RR(n+1)-t_RR(n);
end


# 6. Calculation of cardiac signal characteristics

mean = mean(diff_RR);
standard_dev = std(diff_RR);
bpm = 60/mean;
figure;

for n = 1:length(diff_RR)-1
    hold on;
    plot(diff_RR(n),diff_RR(n+1),'s');
    grid on; grid minor; title('Scatter plot')
end

