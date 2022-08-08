% %% Dados
% %   HP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 0.5 dB, As = 40 dB, GdB = 0 dB)
% %   IIR - Chebyshev II, FIR - Janela Fixa
% 
%% Projeto Filtro IIR - Chebyshev II

clear all;
close all;
clc;

Ap = 0.5; % Ganho na banda de passagem em dB
As = 40; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 0; % dB

fs = f1;
fp = f2;
ws = 2*pi*fs;
wp = 2*pi*fp;
wa = fa * 2 * pi;
tetha_s = ws/(wa/2);
tetha_p = wp/(wa/2);
lambda_s = 2*tan(tetha_s * pi/2);
lambda_p = 2*tan(tetha_p * pi/2);
Os = lambda_p/lambda_s;
Op = 1;

% Filtro Chebyshev 2
Rp = Ap; Rs = As; 
[n,Wn] = cheb2ord(Op,Os,Rp,Rs,'s');
[b,a] = cheby2(n,Rs,Wn,'s');

% Plot protótipo filtro PB
figure(1)
[h1,w1] = freqs(b,a,logspace(-2,1,1000));
semilogx(w1,20*log10(abs(h1)));grid on; ylim([-60 5]);
title('H(p)');xlabel('rad/s');ylabel('dB');
hold on
grid on
plot([10^-2,Os,Os,10^1],[0,0,-As,-As], '--r')
plot([10^-2,1,1],[-Ap,-Ap,-80], '--r')

% Transformação de frequência Lowpass para Highpass
ap = a; bp = b; 
syms p;
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp(p)), 5))

% transformação de frequência
syms s;
Hs(s) = collect(subs(Hp(p),lambda_p/s));%transformação lowpass -> bandstop
[N, D] = numden(Hs(s));
pretty(vpa(Hs(s), 5))

% Normalizando de acordo com p^n
bs = sym2poly(N);
as = sym2poly(D);
an = as(1);
bsn = bs/an;
asn = as/an;
Hsn(s) = poly2sym(bsn, s)/poly2sym(asn, s);
pretty(vpa(Hsn(s), 5))

% Plot filtro HP
figure(2)
[h, w] = freqs(bsn,asn, linspace(0, 100, 10000));
plot(w/pi, 20*log10(abs(h))); grid on;hold on;ylim([-80 5]);xlim([0 2])
title('H(s)');xlabel('rad/s');ylabel('dB');
% Fazer a mascara em cima do LAMBDA
plot([0,lambda_s/pi,lambda_s/pi,2],[-As,-As,0,0], '--r')
plot([lambda_p/pi,lambda_p/pi,2],[-80,-Ap,-Ap], '--r')

syms z;
aux = 2*((z-1)/(z+1));
Hz(z) = collect(subs(Hs(s), aux));
pretty(vpa(Hz(z),3))

[Nz,Dz] = numden(Hz(z));
bz = sym2poly(Nz);
az = sym2poly(Dz);

an = az(1);
bzn = bz/an;
azn = az/an;

Hzn(z) = poly2sym(bzn,z) / poly2sym(azn,z);
pretty(vpa(Hzn(z),5))

figure(3)
subplot(211)
[hz, wz] = freqz(bzn,azn, linspace(0, pi, 1000));
plot(wz/pi*fa/2, 20*log10(abs(hz))); grid on;hold on;ylim([-60 5])
title_txt = ['H(z) - BP - Filtro IIR - Chebyshev II - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara do filtro projetado
plot([0,fs,fs,2000],[-As,-As,0,0], '--r')
plot([fp,fp,2000],[-60,-Ap,-Ap], '--r')
hold off;

subplot(212)
plot(wz/pi*fa/2, 20*log10(abs(hz))); grid on;hold on;ylim([-5 2]);xlim([1250 2000]);
title_txt = ['H(z) - BP - Filtro IIR - Chebyshev II - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara do filtro projetado
plot([0,fs,fs,2000],[-As,-As,0,0], '--r')
plot([fp,fp,2000],[-80,-Ap,-Ap], '--r')
hold off;

figure(4)
subplot(121)
zplane(bzn, azn);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(bzn, azn);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');

%% Projeto Filtro FIR - Janela Fixa

clear all;
clc;

Ap = 0.5; % Ganho na banda de passagem em dB
As = 40; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 0; % dB

g0 = GdB;

f = [1000 1300]; % frequências em Hz

% substituindo de Hz para ômega
w = f/fa*(2*pi); %w = 2*pi*f
ws = w(2);
wp = w(1);

%projeto original
Dw = ws - wp;
%M = ceil(3.32*pi/Dw); % ordem (3.32 tabela Hamming)
M = ceil(3.11*pi/Dw); % ordem (3.11 tabela Hann)

%Ajuste de ordem
    M= M-4; 

%Ajuste do ganho
    %levar o pico para abaixo de 0
    ganho = 0.05501; %ganho dB mediddo no plot do filtro
    g0 = ganho;

% primeiro ajuste de M (N/2)
    wp1 = 0.4998*pi; ws1 = ws; % valores medidos no gráfico
    Dw1 = ws1 - wp1;
    M2 = ceil(M*Dw1/Dw);
    M = M2; 

wc = sqrt(wp*ws);   % frequência de corte, média das frequências
wc = wc + 0.0440;   % ajuste da frequencia de corte
k = 1:M;

% Highpass
bi = -sin(k*wc)./(k*pi); 
b0 = 1 - ( wc/pi);
b = [flip(bi) b0 bi];

m = -M : M;
%dw = 0.04; % filtro hamming
dw = 0.00; % filtro hann
wk = (0.5+dw)+(0.5-dw)*cos(2*pi*m/(2*M+1)); %serve para filtro hamming ou hann
b = b.*wk*10^(-g0/20);

% Acertar a janela no plot
ws = w(2);
wp = w(1);

figure(5)
subplot(211)
[h, w] = freqz(b,1,linspace(0,pi,10000)); 
plot(w/pi*fa/2,20*log10(abs(h))); grid on;ylim([-80 5]);
title_txt = ['H(z) - BP - FIR - Janela Fixa - N = ' num2str(M*2)];
hold on;
title(title_txt);xlabel('Hz');ylabel('dB');
plot([0 f1 f1 2000],[-As -As 0 0], '--red')
plot([f2,f2,2000],[-80 -Ap,-Ap], '--red')
hold off;

subplot(212)
plot(w/pi*fa/2,20*log10(abs(h))); grid on;ylim([-5 2]);xlim([1250 2000]);
hold on;
title_txt = ['H(z) - BP - FIR - Janela Fixa - N = ' num2str(M*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
plot([0 f1 f1 2000],[-As -As 0 0], '--red')
plot([f2,f2,2000],[-80 -Ap,-Ap], '--red')
hold off;

figure(6)
subplot(121)
zplane(b, 1);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(b, 1);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');