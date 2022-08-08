% %% Dados
% %   BP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1200 Hz, f3 = 1250 Hz; f4 = 1300 Hz, Ap = 1 dB, As = 20 dB, GdB = 0 dB)
% %   IIR - Eliptico, FIR - PM
% %
%% Projeto Filtro IIR - Elíptico

clear all;
close all;
clc;

% Projeto Filtro IIR - Eliptico

Ap = 1; % Ganho na banda de passagem em dB
As = 20; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1200; % Hz
f3 = 1250; % Hz
f4 = 1300; % Hz
GdB = 0; % dB

% Frequências
fs1 = f1;
fp1 = f2;
fp2 = f3;
fs2 = f4;
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
% Lowpass -> Bandpass
B = lambda_p2 - lambda_p1;
Os = abs((-lambda_s^2+lambda_0^2)/(B*lambda_s));
Op = 1;

% Filtro elíptico
[n,Wn] = ellipord(Op,Os,Ap,As,'s');
[b,a] = ellip(n,Ap,As,Wn,'s');

% Plot protótipo filtro PB
figure(1)
[h1,w1] = freqs(b,a,logspace(-2,2,10000));
semilogx(w1,20*log10(abs(h1)));grid on; ylim([-30 5]);hold on;
title('H(p)');xlabel('rad/s');ylabel('dB');
plot([10^-2,Os,Os,10^2],[0,0,-As,-As], '--r')
plot([10^-2,1,1],[-Ap,-Ap,-80], '--r')
hold off;

% Transformação de frequência Lowpass para Bandpass
ap = a; bp = b; 
syms p;
Np(p) = poly2sym(bp, p);
Dp(p) = poly2sym(ap, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp(p)), 5))

% transformação de frequência
syms s;
Hs(s) = collect(subs(Hp(p), (s^2 + lambda_0^2)/(B*s))); %transformação lowpass/bandpass
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

% Plot filtro BP
figure(2)
[h, w] = freqs(bsn,asn,linspace(0, 100, 10000));
plot(w/pi, 20*log10(abs(h))); grid on;hold on;ylim([-60 5]);xlim([0 2]);
title('H(s)');xlabel('rad/s');ylabel('dB');
% Fazer a mascara em cima do LAMBDA
plot([0,lambda_s1/pi,lambda_s1/pi,lambda_s2/pi,lambda_s2/pi,2],-[As,As,0,0,As,As], '--r')
plot([lambda_p1/pi,lambda_p1/pi,lambda_p2/pi,lambda_p2/pi],-[80,Ap,Ap,80], '--r')
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
plot(wz/pi*fa/2, 20*log10(abs(hz))); grid on;hold on;ylim([-60 5])
title_txt = ['BP - Filtro IIR - ELÍPTICO - N = ' num2str(n*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara do filtro projetado
plot([0,f1,f1,f4,f4,2000],-[As,As,0,0,As,As], '--r')
plot([f2,f2,f3,f3],-[40,Ap,Ap,40], '--r')
hold off;

subplot(212)
plot(wz/pi*fa/2, 20*log10(abs(hz))); grid on;hold on;ylim([-5 2]);xlim([1190 1260]);
title_txt = ['BP - Filtro IIR - ELÍPTICO - N = ' num2str(n*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara do filtro projetado
Amin = 40;
plot([0,f1,f1,f4,f4,1],-[As,As,0,0,As,As], '--r'); 
plot([f2,f2,f3,f3],-[Amin,Ap,Ap,Amin], '--r'); 
hold off;

figure(4)
subplot(121)
zplane(bzn, azn);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(bzn, azn);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
subplot(211)
% valores retirados do Filter Designer
Den = [1 1.356285095214844 2.376022338867188 1.302505493164063 0.922378540039063]; 
Num = [0.097470283508301 0.132939338684082 0.234899520874023 0.132939338684082 0.097470283508301];
N = Num;
D = Den;
plot(wz/pi*fa/2, 20*log10(abs(hz)));
grid on;hold on;
[hzq, wzq] = freqz(N, D, linspace(0, pi, 10000));
plot(wzq/pi*fa/2, 20*log10(abs(hzq)));ylim([-40 5]);
title_txt = ['BP - Filtro IIR - ELÍPTICO - N = ' num2str(n*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara do filtro projetado
plot([0,f1,f1,f4,f4,2000],-[As,As,0,0,As,As], '--r')
plot([f2,f2,f3,f3],-[40,Ap,Ap,40], '--r')
legend('Filtro projetado','Filtro quantizado');
hold off;

subplot(212)
plot(wz/pi*fa/2, 20*log10(abs(hz)));ylim([-5 2]);xlim([1190 1260]);
title_txt = ['BP - Filtro IIR - ELÍPTICO - N = ' num2str(n*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
grid on;hold on;
plot(wzq/pi*fa/2, 20*log10(abs(hzq)));
% Máscara do filtro projetado
plot([0,f1,f1,f4,f4,2000],-[As,As,0,0,As,As], '--r')
plot([f2,f2,f3,f3],-[40,Ap,Ap,40], '--r')
legend('Filtro projetado','Filtro quantizado');
hold off;

figure(6)
subplot(121)
zplane(bzn, azn);title('Diagrama de pólos e zeros projetado');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
zplane(N, D);title('Diagrama de pólos e zeros quantizado');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%% Projeto BP - Filtro FIR - PM

clear all;
clc;

Ap = 1; % Ganho na banda de passagem em dB
As = 20; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1200; % Hz
f3 = 1250; % Hz
f4 = 1300; % Hz
GdB = 0; % dB

f = [1000 1200 1250 1300]; % frequências em Hz
w = f/fa*(2*pi);
ws1 = w(1)/pi;
wp1 = w(2)/pi;
wp2 = w(3)/pi;
ws2 = w(4)/pi;
mags = [0 1 0];

% [n,fo,ao,w] = firpmord(f,a,dev,fs)
devAs = 10^(-(As+0.5)/20);
devAp = 1-10^(-(Ap/2-0.05)/20);
devs = [devAs devAp devAs];

% calculo da ordem com firpmord
f = f + [0 -150 0 0];
[n,f0,a0,w0] = firpmord(f,mags,devs,fa);

G0 = -Ap/2;
% calculo algoritmo PM
h_pm = firpm(n,f0,a0,w0);
h_pm = h_pm*10^(G0/20);

%clear Hw W
Amin = 60;
figure(7)
subplot(211);
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)));ylim([-Amin 5]);
title_txt = ['BP - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
hold on
% Máscara
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[As,As,0,0,As,As], '--r'); ylim([-60 5]);xlim([800 1500]);
plot([wp1,wp1,wp2,wp2]*fa/2,-[Amin,Ap,Ap,Amin], '--r'); grid on;
hold off;

subplot(212);
[Hw,w] = freqz(h_pm,1,10000);
plot(w*fa/2/pi,20*log10(abs(Hw)));ylim([-Amin 5]);
title_txt = ['BP - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
hold on
% Máscara
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[As,As,0,0,As,As], '--r'); ylim([-5 2]);xlim([998 1302]);
plot([wp1,wp1,wp2,wp2]*fa/2,-[Amin,Ap,Ap,Amin], '--r'); grid on;
hold off;

figure(8)
subplot(121)
zplane(h_pm,1);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(h_pm,1);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
subplot(211)
% valores retirados do Filter Designer
Num = [-0.04296875 0.00390625 0.00390625 0.00390625 0.0048828125 -0.0009765625 0.001953125 ...
        0.0126953125 -0.001953125 -0.0107421875 0.0166015625 0.0126953125 -0.01953125 ...
        0.0009765625 0.0263671875 -0.00781250 -0.015625 0.0185546875 0.0087890625 -0.009765625 ...
        0.001953125 0.001953125 0.005859375 0.0146484375 -0.0205078125 -0.0107421875 0.0478515625 ...
        -0.0048828125 -0.0615234375 0.044921875 0.060546875 -0.0810546875 -0.0263671875 ...
        0.1103515625 -0.0205078125 -0.1044921875 0.07421875 0.07421875 -0.1044921875 ...
        -0.0205078125 0.1103515625 -0.0263671875 -0.0810546875 0.060546875 0.044921875 ...
        -0.0615234375 -0.0048828125 0.0478515625 -0.0107421875 -0.0205078125 0.0146484375 ...
        0.005859375 0.001953125 0.001953125 -0.009765625 0.0087890625 0.0185546875 -0.015625 ...
        -0.0078125 0.0263671875 0.0009765625 -0.01953125 0.0126953125 0.0166015625 -0.0107421875 ...
        -0.001953125 0.0126953125 0.001953125 -0.0009765625 0.0048828125 0.00390625 0.00390625 ...
        0.00390625 -0.04296875] ;
N = Num;
plot(w*fa/2/pi,20*log10(abs(Hw))); grid on;hold on;
[hq, wq] = freqz(N, 1, linspace(0, pi, 10000));
plot(wq/pi*fa/2, 20*log10(abs(hq)));
title_txt = ['BP - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara
Amin = 60;
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[As,As,0,0,As,As], '--r'); ylim([-60 5]);xlim([800 1500]);
plot([wp1,wp1,wp2,wp2]*fa/2,-[Amin,Ap,Ap,Amin], '--r');
legend('Filtro projetado','Filtro quantizado');
hold off;

subplot(212)
plot(w*fa/2/pi,20*log10(abs(Hw)));
title_txt = ['BP - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
grid on;hold on;
plot(wq/pi*fa/2, 20*log10(abs(hq)));
% Máscara
plot([0,ws1,ws1,ws2,ws2,1]*fa/2,-[As,As,0,0,As,As], '--r'); ylim([-5 2]);xlim([998 1302]);
plot([wp1,wp1,wp2,wp2]*fa/2,-[Amin,Ap,Ap,Amin], '--r'); 
legend('Filtro projetado','Filtro quantizado');
hold off;

figure(10)
subplot(121)
zplane(h_pm, 1);title('Diagrama de pólos e zeros projetado');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
zplane(N, 1);title('Diagrama de pólos e zeros quantizado');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
