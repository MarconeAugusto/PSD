% %% Dados
% %   BS - (fa = 6000 Hz, f1 = 1200 Hz; f2 = 1250 Hz, f3 = 1300 Hz; f4 = 1400 Hz, Ap = 0.5 dB, As = 60 dB, GdB = 0 dB)
% %   IIR - Chebyshev I, FIR - PM

%% Projeto Filtro IIR - Chebyshev I

%clear all;
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

%f1 = f1+1;

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

f1 = 1200; % Hz

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)
subplot(211)
% valores retirados do Filter Designer
% Den = [1 -2.66499543190002 8.79464620351791 -15.0964799523354 27.5316314101219...
%       -33.7583267688751 42.2783703207970 -38.2580281496048 35.0287520885468 ...
%       -23.1519310474396 15.6125823259354 -7.05648535490036 3.38096177577972 ...
%       -0.836935400962830 0.255793273448944]; 
% Num = [0.526867985725403 -1.54203569889069 5.62231540679932 -10.6001011729240 ...
%       21.2989954352379 -28.6634656190872 39.4931897521019 -39.2119771838188 ...
%       39.4931897521019 -28.6634656190872 21.2989954352379 -10.6001011729240 ...
%       5.62231540679932 -1.54203569889069 0.526867985725403];
N = Num;
D = Den;
plot(wz/pi*fa/2, 20*log10(abs(hz)));
grid on;hold on;
[hzq, wzq] = freqz(N, D, linspace(0, pi, 10000));
plot(wzq/pi*fa/2, 20*log10(abs(hzq)));ylim([-65 5]);xlim([0 2e3]);
title_txt = ['BP - Filtro IIR - Chebyshev I - N = ' num2str(n*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara do filtro projetado
plot([0,f1,f1,f4,f4,2000],-[Ap,Ap,80,80,Ap,Ap], '--r')
plot([0,f2,f2,f3,f3,2000],-[0,0,As,As,0,0], '--r')
legend('Filtro projetado','Filtro quantizado');
hold off;

subplot(212)
plot(wz/pi*fa/2, 20*log10(abs(hz)));ylim([-5 2]);xlim([998 1302]);
title_txt = ['BP - Filtro IIR - Chebyshev I - N = ' num2str(n*2)];
title(title_txt);xlabel('Hz');ylabel('dB');
grid on;hold on;
plot(wzq/pi*fa/2, 20*log10(abs(hzq)));
% Máscara do filtro projetado
plot([0,f1,f1,f4,f4,2000],-[Ap,Ap,80,80,Ap,Ap], '--r')
plot([0,f2,f2,f3,f3,2000],-[0,0,As,As,0,0], '--r')
legend('Filtro projetado','Filtro quantizado');
hold off;

figure(6)
subplot(121)
zplane(bzn, azn);title('Diagrama de pólos e zeros projetado');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
zplane(N, D);title('Diagrama de pólos e zeros quantizado');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Projeto Filtro FIR - PM

clear all;
%close all;
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
figure(7)
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
Num = [0.00939178466796875 0.00168609619140625 -0.00123596191406250 0.00145721435546875...
    -0.00170898437500000 0.00103759765625000 -0.000633239746093750 0.00136566162109375 ...
    -0.00128936767578125 0.000236511230468750 -0.000617980957031250 0.00144195556640625 ...
    -0.000404357910156250 -9.91821289062500e-05 -0.00117492675781250 0.000854492187500000 ...
    0.000473022460937500 0.000564575195312500 -0.00129699707031250 -0.000411987304687500 ...
    0.000267028808593750 0.00151062011718750 -0.000167846679687500 -0.00112152099609375 ...
    -0.00112915039062500 0.00125122070312500 0.00146484375000000 -0.000129699707031250 ...
    -0.00218200683593750 -0.000602722167968750 0.00175476074218750 0.00189208984375000 ...
    -0.00128936767578125 -0.00239562988281250 -4.57763671875000e-05 0.00273132324218750 ...
    0.00110626220703125 -0.00210571289062500 -0.00225067138671875 0.00122833251953125 ...
    0.00261688232421875 0.000122070312500000 -0.00253295898437500 -0.00111389160156250 ...
    0.00170135498046875 0.00179290771484375 -0.000724792480468750 -0.00168609619140625 ...
    -0.000259399414062500 0.00111389160156250 0.000656127929687500 -0.000160217285156250 ...
    -0.000457763671875000 -0.000556945800781250 -0.000427246093750000 0.000762939453125000 ...
    0.00154876708984375 -8.39233398437500e-05 -0.00249481201171875 -0.00131988525390625 ...
    0.00265502929687500 0.00316619873046875 -0.00178527832031250 -0.00478363037109375 ...
    -0.000190734863281250 0.00558471679687500 0.00283050537109375 -0.00505065917968750 ...
    -0.00550079345703125 0.00309753417968750 0.00740051269531250 -3.81469726562500e-05 ...
    -0.00788879394531250 -0.00341033935546875 0.00667572021484375 0.00637054443359375 ...
    -0.00401306152343750 -0.00801849365234375 0.000534057617187500 0.00794219970703125 ...
    0.00276947021484375 -0.00621032714843750 -0.00502014160156250 0.00344848632812500 ...
    0.00561523437500000 -0.000640869140625000 -0.00457763671875000 -0.00120544433593750 ...
    0.00248718261718750 0.00142669677734375 -0.000396728515625000 5.34057617187500e-05 ...
    -0.000541687011718750 -0.00257873535156250 -0.000503540039062500 0.00493621826171875 ...
    0.00366210937500000 -0.00570678710937500 -0.00825500488281250 0.00374603271484375 ...
    0.0128097534179688 0.00134277343750000 -0.0154800415039063 -0.00894927978515625 ...
    0.0145721435546875 0.0174026489257813 -0.00914764404296875 -0.0243988037109375 ...
    -0.000549316406250000 0.0275115966796875 0.0129699707031250 -0.0249481201171875 ...
    -0.0254974365234375 0.0161514282226563 0.0350112915039063 -0.00212097167968750 ...
    -0.0387039184570313 -0.0146636962890625 0.0348587036132813 0.0307235717773438 ...
    -0.0234298706054688 -0.0423660278320313 0.00615692138671875 0.0466384887695313 ...
    0.0136718750000000 -0.0421066284179688 -0.0319900512695313 0.0292434692382813 ...
    0.0448760986328125 -0.0104598999023438 0.922111511230469 -0.0104598999023438 ...
    0.0448760986328125 0.0292434692382813 -0.0319900512695313 -0.0421066284179688 ...
    0.0136718750000000 0.0466384887695313 0.00615692138671875 -0.0423660278320313 ...
    -0.0234298706054688 0.0307235717773438 0.0348587036132813 -0.0146636962890625 ...
    -0.0387039184570313 -0.00212097167968750 0.0350112915039063 0.0161514282226563 ...
    -0.0254974365234375 -0.0249481201171875 0.0129699707031250 0.0275115966796875 ...
    -0.000549316406250000 -0.0243988037109375 -0.00914764404296875 0.0174026489257813 ...
    0.0145721435546875 -0.00894927978515625 -0.0154800415039063 0.00134277343750000 ...
    0.0128097534179688 0.00374603271484375 -0.00825500488281250 -0.00570678710937500 ...
    0.00366210937500000 0.00493621826171875 -0.000503540039062500 -0.00257873535156250 ...
    -0.000541687011718750 5.34057617187500e-05 -0.000396728515625000 0.00142669677734375 ...
    0.00248718261718750 -0.00120544433593750 -0.00457763671875000 -0.000640869140625000 ...
    0.00561523437500000 0.00344848632812500 -0.00502014160156250 -0.00621032714843750 ...
    0.00276947021484375 0.00794219970703125 0.000534057617187500 -0.00801849365234375 ...
    -0.00401306152343750 0.00637054443359375 0.00667572021484375 -0.00341033935546875 ...
    -0.00788879394531250 -3.81469726562500e-05 0.00740051269531250 0.00309753417968750 ...
    -0.00550079345703125 -0.00505065917968750 0.00283050537109375 0.00558471679687500 ...
    -0.000190734863281250 -0.00478363037109375 -0.00178527832031250 0.00316619873046875 ...
    0.00265502929687500 -0.00131988525390625 -0.00249481201171875 -8.39233398437500e-05 ...
    0.00154876708984375 0.000762939453125000 -0.000427246093750000 -0.000556945800781250 ...
    -0.000457763671875000 -0.000160217285156250 0.000656127929687500 0.00111389160156250 ...
    -0.000259399414062500 -0.00168609619140625 -0.000724792480468750 0.00179290771484375 ...
    0.00170135498046875 -0.00111389160156250 -0.00253295898437500 0.000122070312500000 ...
    0.00261688232421875 0.00122833251953125 -0.00225067138671875 -0.00210571289062500 ...
    0.00110626220703125 0.00273132324218750 -4.57763671875000e-05 -0.00239562988281250 ...
    -0.00128936767578125 0.00189208984375000 0.00175476074218750 -0.000602722167968750 ...
    -0.00218200683593750 -0.000129699707031250 0.00146484375000000 0.00125122070312500 ...
    -0.00112915039062500 -0.00112152099609375 -0.000167846679687500 0.00151062011718750 ...
    0.000267028808593750 -0.000411987304687500 -0.00129699707031250 0.000564575195312500 ...
    0.000473022460937500 0.000854492187500000 -0.00117492675781250 -9.91821289062500e-05 ...
    -0.000404357910156250 0.00144195556640625 -0.000617980957031250 0.000236511230468750 ...
    -0.00128936767578125 0.00136566162109375 -0.000633239746093750 0.00103759765625000 ...
    -0.00170898437500000 0.00145721435546875 -0.00123596191406250 0.00168609619140625 ...
    0.00939178466796875] ;
N = Num;
plot(w*fa/2/pi,20*log10(abs(Hw))); grid on;hold on;
[hq, wq] = freqz(N, 1, linspace(0, pi, 10000));
plot(wq/pi*fa/2, 20*log10(abs(hq)));
title_txt = ['BS - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
% Máscara
Amin = 0;
plot([0,f1,f1,f4,f4,fa/2],-[Ap,Ap,80,80,Ap,Ap], '--r');ylim([-80 5]);xlim([1000 1600]);
plot([0,f2,f2,f3,f3,fa/2],[0,Amin,-As,-As,Amin,0], '--m');
legend('Filtro projetado','Filtro quantizado');
hold off;

subplot(212)
plot(w*fa/2/pi,20*log10(abs(Hw)));ylim([-65 -55]);xlim([1225 1375]);
title_txt = ['BS - Filtro FIR - PM - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
grid on;hold on;
plot(wq/pi*fa/2, 20*log10(abs(hq)));
% Máscara
plot([0,f1,f1,f4,f4,fa/2],-[Ap,Ap,80,80,Ap,Ap], '--r');
plot([0,f2,f2,f3,f3,fa/2],[0,Amin,-As,-As,Amin,0], '--m');
legend('Filtro projetado','Filtro quantizado');
hold off;

figure(10)
subplot(121)
zplane(h_pm, 1);title('Diagrama de pólos e zeros projetado');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
zplane(N, 1);title('Diagrama de pólos e zeros quantizado');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


