function [fig_h, SF] = Denisenko_Nigth_Nh_from_hv_IRI_fval(f, hv, f_N, h_N, ks_title, m, p)
fH =1.29;
tetta = 25.19;

n = length(h_N);
foF = f_N(n);
hmF = h_N(n);
nf = length(f);

for i=1:n-2
    if f_N(i)<=f_N(i+1)&&f_N(i+1)>=f_N(i+2)
        foE = f_N(i+1);
        hoE = h_N(i+1);
        break
    end
end

fF_min = f(1);

for i=1:n-1
    if f_N(i)<=fF_min&&fF_min<=f_N(i+1)
        nE_max = i+1;
    end
end

for i=1:nE_max-2
    if f_N(i)<=f_N(i+1)&&f_N(i+1)>=f_N(i+2)
        hoE = h_N(i+1);
        foE = f_N(i+1);
        nE_min = i+1;
        break
    end
end

for i=nE_min:nE_max-2
    if f_N(i)>=f_N(i+1)&&f_N(i+1)<=f_N(i+2)
        nE_val = i+1;
        hval = h_N(i+1);
        fval = f_N(i+1);
        break
    end
end

h_NE = h_N(1:nE_val+1);
f_NE = f_N(1:nE_val+1);
h0 = h_N(1);
%===============================================

for i=1:nf
    dhvE(i) = HvE_o_tabl(f(i),fH,tetta,h0,hoE,hval,h_NE,f_NE);
    dhvF(i) = hv(i) - dhvE(i);
end

% plot(f,dhvE,'-b',f,dhvF,'-r');

options = optimset('TolX',1e-4);
%m = 2; p = 1;
par = [2,hmF];

[par,SF] = fminsearch(@SUM_dhv_model,par,options,f,dhvF,hval,fval,foF,p,m,fH,tetta);
a = par(1);
hmF = par(2);

sigmaF = sqrt(SF/(length(f)-2));
[dhv_calc] = HV_theor_model(f,fH,tetta,p,m,a,hval,hmF,fval,foF);
hv_calc = dhv_calc + dhvE;

dh = (hmF-hval)/25;
hF = hval:dh:hmF;
fpF = fN_model(hF,p,m,a,hval,hmF,fval,foF);

fig_h = figure('NumberTitle','off','Name',['F-region: ','m=',num2str(m),'  p=',num2str(p),...
    '  sigF=',num2str(sigmaF),'  foF=',num2str(foF),'  hmF=',num2str(hmF),'  a=',num2str(a)]);
plot (f_N,h_N,'-b',fpF,hF,'-r',f,hv,'-*k',f,hv_calc,'-^r','LineWidth',2);
ylabel('h & h_v, km')
xlabel('f & f_N, MHz')
legend('h(IRI)','h_c_o_r(h_v)','h_v(exper)','h_v(calc)','Location','southeast')
title(ks_title)
grid on

%====================================================
function [h_tab] = COR(hm,h_N)
np = length(h_N);
dh = (hm - h_N(1))/(np - 1);
h_tab = h_N(1):dh:hm;
%====================================================
% function [h_tab] = COR(hm,h_N)
% 
% n = length(h_N);
% h0 = h_N(1);
% hm0 = h_N(n);
% dhm = (hm - hm0);
% dh = dhm*((h_N - h0)./(hm0 - h0));
% h_tab = h_N + dh ;
%====================================================
function [hv_calc] = HV_theor_tabl(f,fH,tetta,hoE,h_N,fN_tab)

for i=1:length(f)
    hv_calc(i) = Hv_o_tabl(f(i),fH,tetta,hoE,h_N,fN_tab);
end
%====================================================
function dhv2 = SUM_dhv_tabl(hm,f,hv,fH,tetta,hoE,h_N,fN_tab)

[h_tab] = COR(hm,h_N);
n = length(f);
s = 0;
for i=1:n
    hvt(i) = Hv_o_tabl(f(i),fH,tetta,hoE,h_tab,fN_tab);
    s = s + ((hv(i)-hvt(i)))^2;
end
dhv2 = s;
%====================================================
function hv = Hv_o_tabl(f,fH,tetta,hoE,h_tab,fN_tab)

tol = 1e-4;
h0 = h_tab(1);
Xpd = 0.99;
fpd = f*sqrt(Xpd);
hpd = H_refl_tabl(fpd,h_tab,fN_tab);
qD1 = quadl(@Mvh_o_tabl,h0,hoE,tol,0,f,fH,tetta,h_tab,fN_tab);
qD2 = quadl(@Mvh_o_tabl,hoE,hpd,tol,0,f,fH,tetta,h_tab,fN_tab);
hvD = h0+qD1+qD2;

hr = H_refl_tabl(f,h_tab,fN_tab);
Hr = (hr - hpd)/(1 - Xpd);
q = quadl(@Mvr_o,0,sqrt(1-Xpd),tol,0,f,fH,tetta);
hvr = 2*Hr*q;

hv = hvD + hvr;
%===========================================================
function Mgh_o = Mvh_o_tabl(h,f,fH,tetta,h_tab,fN_tab)
% ядро гр. пок. прел. для интегрирования по высоте;

tetta = pi*tetta/180;
cost = cos(tetta);
cos2 = cost^2;
sin2 = 1 - cos2;
Y = fH/f;
YL = Y*cost;
fp = fN_tabl(h,h_tab,fN_tab);
X = (fp/f).^2;
Mgh_o = Mv_o(X,cost,cos2,sin2,Y,YL)./sqrt(1-X);
%============================================================
% function Mgh_x = Mvh_x(h,f,fH,tetta,h_tab,fN_tab)
% % ядро гр. пок. прел. для интегрирования по высоте;
% 
% tetta = pi*tetta/180;
% cost = cos(tetta);
% cos2 = cost^2;
% sin2 = 1 - cos2;
% Y = fH/f;
% YL = Y*cost;
% fp = fN_tabl(h,h_tab,fN_tab);
% X = (fp/f).^2;
% Mgh_x = Mv_x(X,cost,cos2,sin2,Y,YL)./sqrt(1-X./(1-Y));
%============================================================
function Mgr_o = Mvr_o(t,f,fH,tetta)
% ядро гр. пок. прел. для линейного слоя;

tetta = pi*tetta/180;
cost = cos(tetta);
cos2 = cost^2;
sin2 = 1 - cos2;
Y = fH/f;
YL = Y*cost;
X = 1 - t.^2;
Mgr_o = Mv_o(X,cost,cos2,sin2,Y,YL);
%============================================================
function M_o = Mv_o(X,cost,cos2,sin2,Y,YL)
% Ядро груп. пок. преломления.

d = (2.*(1 - X)./Y).*cost/sin2;
P = d./(1 + sqrt(1 + d.^2));
P2 = P.^2;
PYL = P.*YL;
A = sin2 + X.*cos2.*((1 - P2)./(1 + P2)).*(X - P2)./(1 + PYL).^2;
B = sqrt((1 + PYL)./(1 - P2.*cos2)./sin2);
M_o = B.*A;
%======================================================
% function M_x = Mv_x(X,cost,cos2,sin2,Y,YL)
% % Ядро груп. пок. преломления.
% 
% d = (2.*(1 - X)./Y).*cost/sin2;
% P = d./(1 + sqrt(1 + d.^2));
% P2 = P.^2;
% A = sqrt((1+cost*P).*(P-YL)./(P+cost));
% B = 1+cos2/sin2*X.*(1-P2).*(1-X.*P2)./((1+P2).*(P-YL).^2);
% M_x = B.*A;
%======================================================
function fN = fN_tabl(h,h_tab,fN_tab)
x = h_tab;
y = fN_tab;
np = length(h);
n_tab = length(x);
for i=1:np
    for j=1:(n_tab-1)
        if h(i)>=x(j)&&h(i)<=x(j+1)
            fN(i)=y(j)+(h(i)-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j));
            break
        end
    end
end
%============================================================
function hr = H_refl_tabl(fN,h_tab,fN_tab)
x = fN_tab;
y = h_tab;
n = length(x);
for j=1:(n-1)
    if fN>=x(j)&&fN<=x(j+1)
        hr=y(j)+(fN-x(j))/(x(j+1)-x(j))*(y(j+1)-y(j));
        break
    end
end
%============================================================
function dhvE = HvE_o_tabl(f,fH,tetta,h0,hoE,hmax,h_NE,f_NE)

tol = 1e-4;
qE1 = quadl(@Mvh_o_tabl,h0,hoE,tol,0,f,fH,tetta,h_NE,f_NE);
qE2 = quadl(@Mvh_o_tabl,hoE,hmax,tol,0,f,fH,tetta,h_NE,f_NE);
dhvE = h0+qE1+qE2; 
%===========================================================
function [hv_calc] = HV_theor_model(f,fH,tetta,p,m,a,h0,hm,fp0,fc)

for i=1:length(f)
    hv_calc(i) = Hv_o_model(f(i),fH,tetta,p,m,a,h0,hm,fp0,fc);
end
%====================================================
function dhv2 = SUM_dhv_model(par,f,hv,h0,fp0,fc,p,m,fH,tetta)

a = par(1);
hm = par(2);
n = length(f);
s = 0;
for i=1:n
    hvt(i) = Hv_o_model(f(i),fH,tetta,p,m,a,h0,hm,fp0,fc);
    s = s + ((hv(i)-hvt(i)))^2;
end
dhv2 = s;
%====================================================
function hv = Hv_o_model(f,fH,tetta,p,m,a,h0,hm,fp0,fc)

tol = 1e-4;
Xpd = 0.99;
fpd = f*sqrt(Xpd);
hpd = H_refl_model(fpd,p,m,a,h0,hm,fp0,fc);
qD = quadl(@Mvh_o_model,h0,hpd,tol,0,f,fH,tetta,p,m,a,h0,hm,fp0,fc);
hvD = qD;

hr = H_refl_model(f,p,m,a,h0,hm,fp0,fc);
Hr = (hr - hpd)/(1 - Xpd);
q = quadl(@Mvr_o,0,sqrt(1-Xpd),tol,0,f,fH,tetta);
hvr = 2*Hr*q;

hv = hvD + hvr;
%===========================================================
function Mgh_o = Mvh_o_model(h,f,fH,tetta,p,m,a,h0,hm,fp0,fc)
% ядро гр. пок. прел. для интегрирования по высоте;

tetta = pi*tetta/180;
cot = cos(tetta);
co2 = cot^2;
si2 = 1 - co2;
Y = fH/f;
YL = Y*cot;
fp = fN_model(h,p,m,a,h0,hm,fp0,fc);
X = (fp/f).^2;
Mgh_o = Mv_o(X,cot,co2,si2,Y,YL)./sqrt(1-X);
%============================================================
function fp = fN_model(h,p,m,a,h0,hm,fp0,fc)

x = (hm-h)./(hm-h0);
y0 = log(fc/fp0);
b = a/y0^(1/m)-1;
y = (a*x).^m./(1+b*x.^p).^m;
fp = fc*exp(-y);
%=======================================
function dfr = DF_refl(h,fpr,p,m,a,h0,hm,fp0,fc)

y0 = log(fc/fp0);
x = (hm-h)/(hm-h0);
b = a/y0^(1/m)-1;
y = (a*x)^m/(1+b*x^p)^m;
fp = fc*exp(-y);
dfr = fp - fpr;
%============================================================
function hr = H_refl_model(fpr,p,m,a,h0,hm,fp0,fc)

options = optimset('TolX',1e-6);
hr = fzero(@DF_refl,[h0,hm],options,fpr,p,m,a,h0,hm,fp0,fc);
%============================================================

function [y_m] = mean_3(y)

n = length(y);
y_m(1)=(5*y(1)+2*y(2)-y(3))/6;
y_m(n)=(5*y(n)+2*y(n-1)-y(n-2))/6;
for i=2:(n-1)
    y_m(i)=(y(i-1)+y(i)+y(i+1))/3;
end
% ====================================
function [y_m] = mean_5(y)

n = length(y);
y_m(1)=(3*y(1)+2*y(2)+y(3)-y(5))/5;
y_m(2)=(4*y(1)+3*y(2)+2*y(3)+y(4))/10;
y_m(n-1)=(4*y(n)+3*y(n-1)+2*y(n-2)+y(n-3))/10;
y_m(n)=(3*y(n)+2*y(n-1)+y(n-2)-y(n-4))/5;
for i=3:(n-2)
    y_m(i)=(y(i-2)+y(i-1)+y(i)+y(i+1)+y(i+2))/5;
end
%===============================================



