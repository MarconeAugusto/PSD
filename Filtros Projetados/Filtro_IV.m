% %% Dados
% %   BS - (fa = 6000 Hz, f1 = 1200 Hz; f2 = 1250 Hz, f3 = 1300 Hz; f4 = 1400 Hz, Ap = 0.5 dB, As = 60 dB, GdB = 0 dB)
% %   IIR - Chebyshev I, FIR - PM

%% Projeto Filtro IIR - Chebyshev I

clear all;
close all;
clc;

Ap = 0.5; % Ganho na banda de passagem em dB
As = 60; % Atenuação no stopband em dB
fa = 6000; % Hz
f1 = 1200; % Hz
f2 = 1250; % Hz
f3 = 1300; % Hz
f4 = 1400; % Hz
GdB = 0; % dB

% Frequências
fs1 = f2;
fs2 = f3;
fp1 = f1;
fp2 = f4;
% Frequência -> omega
ws1 = 2*pi*fs1;
wp1 = 2*pi*fp1;
wp2 = 2*pi*fp2;
ws2 = 2*pi*fs2;
wa = fa*2*pi;
% Cálculo do tetha
tetha_s1 = ws1/(wa/2);
tetha_p1 = wp1/(wa/2);
tetha_s2 = ws2/(wa/2);
tetha_p2 = wp2/(wa/2);
% Cálculo do lambda
lambda_s1 = 2*tan(tetha_s1 * pi/2);
lambda_s2 = 2*tan(tetha_s2 * pi/2);
lambda_p1 = 2*tan(tetha_p1 * pi/2);
lambda_p2 = 2*tan(tetha_p2 * pi/2);
lambda_0 = sqrt(lambda_p2*lambda_p1);
lambda_s = min(lambda_s1,lambda_s2);
% Lowpass -> Bandstop
B = lambda_p2 - lambda_p1;
Os = abs((B*lambda_s)/(-lambda_s^2+lambda_0^2));
Op = 1;

% Filtro Chebyshev 1
Rp = Ap; Rs = As;
[n,Wn] = cheb1ord(Op,Os,Rp,Rs,'s');
[b,a] = cheby1(n,Rp,Wn,'s');

% Plot protótipo filtro PB
figure(1)
[h1,w1] = freqs(b,a,logspace(-2,1,1000));
semilogx(w1,20*log10(abs(h1)));grid on; ylim([-65 5]);
title('H(p)');hold on;xlabel('rad/s');ylabel('dB');
plot([10^-2,Os,Os,10^2],[0,0,-As,-As], '--r')
plot([10^-2,1,1],[-Ap,-Ap,-80], '--r')

% Transformação de frequência Lowpass para Bandstop
ap = a; bp = b; 
syms p;
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp(p)), 5))

% transformação de frequência
syms s;
Hs(s) = collect(subs(Hp(p),((B*s)/(s^2 + lambda_0^2))));%transformação lowpass/bandstop
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

% Plot filtro BS
figure(2)
[h, w] = freqs(bsn,asn, linspace(0, 8, 10000));
plot(w/pi, 20*log10(abs(h))); grid on;hold on;ylim([-65 5]);xlim([0 2])
title('H(s)');xlabel('rad/s');ylabel('dB');
% Fazer a mascara em cima do LAMBDA
plot([0,lambda_s1/pi,lambda_s1/pi,lambda_s2/pi,lambda_s2/pi,2],-[0,0,As,As,0,0], '--r')
plot([0,lambda_p1/pi,lambda_p1/pi],-[Ap,Ap,80], '--r')
plot([lambda_p2/pi,lambda_p2/pi,2],-[80,Ap,Ap], '--r')
hold off;

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
plot(wz/pi*fa/2, 20*log10(abs(hz))); grid on;hold on;ylim([-65 5]);xlim([0 2e3]);
title_txt = ['BS - Filtro IIR - Chebyshev I - N = ' num2str(n*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara do filtro projetado
plot([0,f1,f1,f4,f4,2000],-[Ap,Ap,80,80,Ap,Ap], '--r')
plot([0,f2,f2,f3,f3,2000],-[0,0,As,As,0,0], '--r')
hold off;

subplot(212)
plot(wz/pi*fa/2, 20*log10(abs(hz))); grid on;hold on;ylim([-5 2]);xlim([998 1302]);
title_txt = ['BP - Filtro IIR - Chebyshev I - N = ' num2str(n*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara do filtro projetado
plot([0,f1,f1,f4,f4,2000],-[Ap,Ap,80,80,Ap,Ap], '--r')
plot([0,f2,f2,f3,f3,2000],-[0,0,As,As,0,0], '--r')
hold off;

figure(4)
subplot(121)
zplane(bzn, azn);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(bzn, azn);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');

%% Projeto Filtro FIR - PM

clear all;
clc;

Ap = 0.5; % Ganho na banda de passagem em dB
As = 60; % Atenuação no stopband em dB
fa = 6000; % Hz
f1 = 1200; % Hz
f2 = 1250; % Hz
f3 = 1300; % Hz
f4 = 1400; % Hz
GdB = 0; % dB

f = [1200 1250 1300 1400]; % frequências em Hz
w = f/fa*(2*pi);
ws1 = w(1)/pi;
wp1 = w(2)/pi;
wp2 = w(3)/pi;
ws2 = w(4)/pi;
mags = [1 0 1];

% [n,fo,ao,w] = firpmord(f,a,dev,fs)
devAs = 10^(-(As-3.5)/20);
devAp = 1-10^(-(Ap/2+0.05)/20);
devs = [devAp devAs devAp];

% calculo da ordem com firpmord
f = f + [0 0 40 0];
[n,f0,a0,w0] = firpmord(f,mags,devs,fa);

G0 = -Ap/2;
% calculo algoritmo PM
h_pm = firpm(n,f0,a0,w0);
h_pm = h_pm*10^(G0/20);

%clear Hw W
figure(4)
subplot(211)
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)))
title_txt = ['BS - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
hold on;
% Máscara
Amin = 0;
plot([0,f1,f1,f4,f4,fa/2],-[Ap,Ap,80,80,Ap,Ap], '--r');ylim([-80 5]);xlim([1000 1600]);
plot([0,f2,f2,f3,f3,fa/2],[Amin,Amin,-As,-As,Amin,Amin], '--m');grid on;
hold off;

subplot(212)
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)));ylim([-65 -55]);xlim([1225 1375]);
title_txt = ['BS - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
hold on;
% Máscara
plot([0,f1,f1,f4,f4,fa/2],-[Ap,Ap,80,80,Ap,Ap], '--r');
plot([0,f2,f2,f3,f3,fa/2],[0,Amin,-As,-As,Amin,0], '--m');grid on;
hold off;


figure(6)
subplot(121)
zplane(h_pm,1);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(h_pm,1);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');


