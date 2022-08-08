% %% Dados
% %   LP - (fa = 4000 Hz, f1 = 1000 Hz; f2 = 1300 Hz, Ap = 2 dB, As = 30 dB, GdB = 5 dB)
% %   IIR - Butterworth, FIR - Janela Ajustável

%% Projeto Filtro IIR - Butterworth

clear all;
close all;
clc;

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 5; % dB

fp = f1;
fs = f2;

thetap = fp/(fa/2);
thetas = fs/(fa/2);

lambdap = 2*tan((thetap*pi)/2);
lambdas = 2*tan((thetas*pi)/2);

Ws = lambdas/lambdap;
Wp = 1;

E = sqrt(10^(0.1*Ap)-1);
n = ceil((log((10^(0.1*As)-1)/E^2))/(2*log(Ws))); % ordem do filtro
k = 1:n;
pk = E^(-1/n)*exp((1j*(2*k+n-1)/(2*n)*pi));

a = real(poly(pk));         % denominador
b = 10^(GdB/20)*(a(end));   % numerador

syms p
Np(p) = poly2sym(b, p);
Dp(p) = poly2sym(a, p);
Hp(p) = Np(p) / Dp(p);
pretty(vpa(collect(Hp), 5)); % collect simplifica ao maximo a funcao

figure(1)
[h, w] = freqs(b,a, logspace(-1, 1, 10000));
semilogx(w, 20*log10(abs(h)));ylim([-60 10]);grid on;hold on;
% Máscara
plot([0.1,Ws,Ws,10],[0,0,-As,-As]+GdB, '--r')
plot([0.1,Wp,Wp,],[-Ap,-Ap,-80]+GdB, '--r')
title('H(p)');xlabel('rad/s');ylabel('dB');

syms s
Hs(s) = collect(subs(Hp(p), s/lambdap));
[Ns, Ds] = numden(Hs(s));
pretty(vpa(Hs(s), 3))

bs = sym2poly(Ns);
as = sym2poly(Ds);

an = as(1);
bsn = bs/an;
asn = as/an;

Hsn(s) = poly2sym(bsn,s) / poly2sym(asn,s);
pretty(vpa(Hsn(s),5))

figure(2)
[hs,ws] = freqs(bsn,asn,linspace(0,8,10000));
plot(ws/pi, 20*log10(abs(hs)));ylim([-60 10]);
title('H(s)');xlabel('rad/s');ylabel('dB');
grid on; hold on;
% Fazer a máscara em cima do LAMBDA
plot([0,lambdas/pi,lambdas/pi,2],[0,0,-As,-As]+GdB, '--r')
plot([0,lambdap/pi,lambdap/pi],[-Ap,-Ap,-80]+GdB, '--r')

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
[hz, wz] = freqz(bzn, azn, linspace(0, pi, 10000));
plot(wz/pi*fa/2, 20*log10(abs(hz)));ylim([-40 10]);
title_txt = ['H(z) - BP - Filtro IIR - Butterworth - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
grid on;hold on;
% Máscara
plot([0.01,fs,fs,2000],[0,0,-As,-As]+GdB, '--r')
plot([0.01,fp,fp,],[-Ap,-Ap,-80]+GdB, '--r')

subplot(212)
plot(wz/pi*fa/2, 20*log10(abs(hz)));ylim([0 8]);xlim([800 1100]);
title_txt = ['H(z) - BP - Filtro IIR - Butterworth - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
grid on;hold on;
% Máscara
plot([0.01,fs,fs,2000],[0,0,-As,-As]+GdB, '--r')
plot([0.01,fp,fp,],[-Ap,-Ap,-80]+GdB, '--r')

figure(4)
subplot(121)
zplane(bzn, azn);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(bzn, azn);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');
%% Projeto Filtro FIR - Janela Ajustável

clear all;
clc;

Ap = 2; % Ganho na banda de passagem em dB
As = 30; % Atenuação no stopband em dB
fa = 4000; % Hz
f1 = 1000; % Hz
f2 = 1300; % Hz
GdB = 5; % dB

fp = f1; fs = f2;

ganho = GdB;

wp = fp/(fa)*(2*pi);
ws = fs/(fa)*(2*pi);
dw = (ws - wp);
wc = sqrt(wp*ws); % frequência de corte
betha = 0.5842*(As-21)^0.4+0.07886*(As-21);
n = ceil((As - 8)/(2.285*dw) +1);

if mod(n,2) == 1 
    %impar
    n = n+1;
end 

% Ajuste do ganho
    GdB = GdB - 0.0379;
% % Primeiro ajuste de M
    wp1 = 0.5318*pi; ws1 = 0.6373*pi;
    Dw1 = ws1-wp1;
    n1 = ceil(n*Dw1/dw);
    n = n1;
% % Segundo ajuste de M
    wp2 = 0.5178*pi; ws2 = 0.6641*pi;
    Dw2 = ws2-wp2;
    n2 = ceil(n*Dw2/dw);
    n = n2;
    wc = wc - 0.015*pi;

Jkaiser = kaiser(n+1,betha);

k = 1:(n/2);
b1 = sin(k*wc)./(k*pi);
b0 = wc/pi; % L'Hospital do b0 acima
b = [flip(b1) b0 b1];
b = b.'; % matriz b transposta
b = b.*Jkaiser*10^(GdB/20)*10^(-0.189/20);

figure(5)
subplot(211)
[h, w] = freqz(b, 1, linspace(0,pi,10000));
plot(w/pi*fa/2, 20*log10(abs(h))); grid on;
ylim([-60 10])
hold on;title_txt = ['BP - Filtro FIR - Janela ajustável Kaiser - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');
fmax = fa/2;
% Máscara
plot([0,fs,fs,fmax],[0,0,-As,-As]+ganho, '--red')
plot([0,fp,fp],[-Ap,-Ap,-140]+ganho, '--red');hold off;

subplot(212)
plot(w/pi*fa/2, 20*log10(abs(h))); grid on;hold on;
title_txt = ['BP - Filtro FIR - Janela ajustável Kaiser - N = ' num2str(n)];
title(title_txt);xlabel('Hz');ylabel('dB');xlim([800 1020]); ylim([-5 10]);
% Máscara
plot([0,fs,fs,fmax],[0,0,-As,-As]+ganho, '--red')
plot([0,fp,fp],[-Ap,-Ap,-140]+ganho, '--red');
hold off;

figure(6)
subplot(121)
zplane(b,1);title('Diagrama de pólos e zeros');xlabel('Parte real');ylabel('Parte imaginária');axis([-2 2 -3 3]);
subplot(122)
grpdelay(b,1);title('Atraso de grupo');
xlabel('Frequência normalizada [x\pi rad/amostra]');
ylabel('Atraso de grupo [amostra]');

