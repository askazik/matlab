function Nigth_Nh_from_hv_IRI_parab_a_b_c
close all 
clc
fH =1.29;
tetta = 25.19;
cost = cos(pi*tetta/180);
cos2 = cost^2;
sin2 = 1 - cos2;

load ('2016-09-24_15-00r.mat');
% load ('2016-09-30_15-00r.mat');

ks_title = timesound;

tmp = flipud(trace_F2o);
f = tmp(:,1)';
hv = tmp(:,2)';
% foF = foF2;

tmp = iri_profile;
hu_tab = tmp(:,2)';
fu_tab = tmp(:,1)';
ind = find(fu_tab == max(fu_tab));
hu_tab = hu_tab(1:ind);
fu_tab = fu_tab(1:ind);

n = length(hu_tab);
h_N = hu_tab; % высоты из IRI
f_N = fu_tab; % плазм частоты IRI
foF = f_N(n);
hmF = h_N(n);
nf = length(f);

% f = f(1:nf-1);
% hv = hv(1:nf-1);
% nf = length(f);

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
        nF_min = i+1;
    end
end

for i=1:nF_min-2
    if f_N(i)<=f_N(i+1)&&f_N(i+1)>=f_N(i+2)
        hoE = h_N(i+1);
        foE = f_N(i+1);
        nE_min = i+1;
        break
    end
end

for i=nE_min:nF_min-2
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
    dhvE(i) = HvE_o_tabl(f(i),fH,cost,cos2,sin2,h0,hoE,hval,h_NE,f_NE);
    dhvF(i) = hv(i) - dhvE(i);
end

nE = 3;

fE = f(1:nE);
dhv_E = dhvF(1:nE);
d_hvE = dhvE(1:nE);
fF = f(nE+1:nf);
dhvF = dhvF(nE+1:nf);
dhvEF = dhvE(nE+1:nf);
nF = length(fF);

options = optimset('TolX',1e-4);

pv = 2;
Hv = 60;
[Hval,Sval] = fminsearch(@SUM_dhv_parab,Hv,options,pv,fE,dhv_E,hval,fval,fH,cost,cos2,sin2);
% SUM_dhv_parab(par,pv,f,hv,hval,fval,fH,cost,cos2,sin2)

sig_E = sqrt(Sval/(nE));
%===================================================================
[dhv_parab] = HV_theor_parab(fE,fH,cost,cos2,sin2,pv,Hval,hval,fval);
hvE_calc = d_hvE + dhv_parab;

fpF = fE(nE);
% hpF = hval + Hval*(fpF/fval-1)^(1/pv);
hpF = hval + Hval*((fpF/fval)^2-1)^(1/pv);
% dif = pv*fval*((hpF-hval)/Hval)^(pv-1)/Hval; 
dif = (fval^2/fpF)*(hpF-hval)/Hval/Hval;
dh = (hpF-hval)/50;
hp_E = hval:dh:hpF;
% fp_E = fval*(1+((hp_E-hval)./Hval).^pv);
fp_E = fval*sqrt((1+((hp_E-hval)./Hval).^pv));

dhv_parab = [];
for i=1:nF
    dhv_parab(i) = HvE_o_parab(fF(i),fH,cost,cos2,sin2,pv,Hval,hval,hpF,fval);
    dhv_F(i) = dhvF(i) - dhv_parab(i);
end

m = 2; c = 0;

par = [c,hmF];
[par,SF] = fminsearch(@SUM_dhv_model,par,options,fF,dhv_F,hpF,fpF,foF,dif,m,fH,cost,cos2,sin2);
c = par(1);
hmF = par(2);
sig_F = sqrt(SF/(length(fF)-2));

fip = (log(foF/fpF))^(1/m);
fipr = (log(foF/fpF))^(1/m-1)/m*(dif*(hmF-hpF)/fpF);
b = -(-fip+fipr+c*(fip+fipr))/fipr;
a = (1+b+c)*fip;

[dhv_calc] = HV_theor_model(fF,fH,cost,cos2,sin2,c,m,hpF,hmF,fpF,foF,dif);
hvF_calc = dhv_calc + dhvEF + dhv_parab;
hv_calc = [hvE_calc,hvF_calc];

dh = (hmF-hpF)/50;
hp_F = hpF+dh:dh:hmF;
fp_F = fN_model(hp_F,c,m,hpF,hmF,fpF,foF,dif);
hp = [hp_E,hp_F];
fp = [fp_E,fp_F];

figure('NumberTitle','off','Name',['nE=',num2str(nE),'  p=',num2str(pv),'  Hv=',num2str(Hval),'  m=',num2str(m),...
    '  sigE=',num2str(sig_E),'  sigF=',num2str(sig_F),'  foE=',num2str(foE),'  foF=',num2str(foF),...
    '  hmF=',num2str(hmF),'  a=',num2str(a),'  b=',num2str(b),'  c=',num2str(c)]);
plot (f_N,h_N,'-b',fp,hp,'-r',f,hv,'-*k',f,hv_calc,'-^r','LineWidth',2);
ylabel('h & h_v, km')
xlabel('f & f_N, MHz')
legend('h(IRI)','h_c_o_r(h_v)','h_v(exper)','h_v(calc)',0)
title(ks_title)
grid on

%===========================================================
function Mgh_o = Mvh_o_tabl(h,f,fH,cost,cos2,sin2,h_tab,fN_tab)
% ядро гр. пок. прел. для интегрирования по высоте;

Y = fH/f;
YL = Y*cost;
fp = fN_tabl(h,h_tab,fN_tab);
X = (fp/f).^2;
Mgh_o = Mv_o(X,cost,cos2,sin2,Y,YL)./sqrt(1-X);
%============================================================
function Mgr_o = Mvr_o(t,f,fH,cost,cos2,sin2)
% ядро гр. пок. прел. для линейного слоя;

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
function dhvE = HvE_o_tabl(f,fH,cost,cos2,sin2,h0,hoE,hmax,h_NE,f_NE)

tol = 1e-4;
% qE1 = quadl(@Mvh_o_tabl,h0,hoE,tol,0,f,fH,tetta,h_NE,f_NE);
% qE2 = quadl(@Mvh_o_tabl,hoE,hmax,tol,0,f,fH,tetta,h_NE,f_NE);
dhE = hoE-h0;
hp = h0+dhE/2;
y1 = Mvh_o_tabl(h0,f,fH,cost,cos2,sin2,h_NE,f_NE);
y2 = Mvh_o_tabl(hp,f,fH,cost,cos2,sin2,h_NE,f_NE);
y3 = Mvh_o_tabl(hoE,f,fH,cost,cos2,sin2,h_NE,f_NE);
qE1 = dhE/6*(y1+4*y2+y3);

dhV = hmax-hoE;
hp = hoE+dhE/2;
y1 = y3;
y2 = Mvh_o_tabl(hp,f,fH,cost,cos2,sin2,h_NE,f_NE);
y3 = Mvh_o_tabl(hmax,f,fH,cost,cos2,sin2,h_NE,f_NE);
qE2 = dhV/6*(y1+4*y2+y3);

dhvE = h0+qE1+qE2; 
%===========================================================
function [hv_calc] = HV_theor_model(f,fH,cost,cos2,sin2,c,m,h0,hm,fp0,fc,dif)

for i=1:length(f)
    hv_calc(i) = Hv_o_model(f(i),fH,cost,cos2,sin2,c,m,h0,hm,fp0,fc,dif);
end
%====================================================
function dhv2 = SUM_dhv_model(par,f,hv,h0,fp0,fc,dif,m,fH,cost,cos2,sin2)

c = par(1);
hm = par(2);
n = length(f);
s = 0;
for i=1:n
    hvt(i) = Hv_o_model(f(i),fH,cost,cos2,sin2,c,m,h0,hm,fp0,fc,dif);
    s = s + ((hv(i)-hvt(i)))^2;
end
dhv2 = s;
%====================================================
function hv = Hv_o_model(f,fH,cost,cos2,sin2,c,m,h0,hm,fp0,fc,dif)

tol = 1e-4;
Xpd = 0.99;
fpd = f*sqrt(Xpd);
hpd = H_refl_model(fpd,c,m,h0,hm,fp0,fc,dif);
qD = quadl(@Mvh_o_model,h0,hpd,tol,0,f,fH,cost,cos2,sin2,c,m,h0,hm,fp0,fc,dif);
hvD = qD;

hr = H_refl_model(f,c,m,h0,hm,fp0,fc,dif);
Hr = (hr - hpd)/(1 - Xpd);
tm = sqrt(1-Xpd); dh = tm/6;
q = dh*(Mvr_o(0,f,fH,cost,cos2,sin2)+4*Mvr_o(tm/2,f,fH,cost,cos2,sin2)+Mvr_o(tm,f,fH,cost,cos2,sin2));
hvr = 2*Hr*q;

hv = hvD + hvr;
%===========================================================
function Mgh_o = Mvh_o_model(h,f,fH,cost,cos2,sin2,c,m,h0,hm,fp0,fc,dif)
% ядро гр. пок. прел. для интегрирования по высоте;

Y = fH/f;
YL = Y*cost;
fp = fN_model(h,c,m,h0,hm,fp0,fc,dif);
X = (fp/f).^2; 
Mgh_o = Mv_o(X,cost,cos2,sin2,Y,YL)./sqrt(1-X);
%============================================================
function fp = fN_model(h,c,m,h0,hm,fp0,fc,dif)

x = (hm-h)./(hm-h0);
fip = (log(fc/fp0))^(1/m);
fipr = (log(fc/fp0))^(1/m-1)/m*(dif*(hm-h0)/fp0);
b = -(-fip+fipr+c*(fip+fipr))/fipr;
a = (1+b+c)*fip;
y = (a*x).^m./(1+b*x+c*x.^2).^m;
fp = fc*exp(-y);
%=======================================
function dfr = DF_refl(h,fpr,c,m,h0,hm,fp0,fc,dif)

x = (hm-h)/(hm-h0);
fip = (log(fc/fp0))^(1/m);
fipr = (log(fc/fp0))^(1/m-1)/m*(dif*(hm-h0)/fp0);
b = -(-fip+fipr+c*(fip+fipr))/fipr;
a = (1+b+c)*fip;
y = (a*x)^m/(1+b*x+c*x^2)^m;
fp = fc*exp(-y);
dfr = fp - fpr;
%============================================================
function hr = H_refl_model(fpr,c,m,h0,hm,fp0,fc,dif)

options = optimset('TolX',1e-6);
hr = fzero(@DF_refl,[h0,hm],options,fpr,c,m,h0,hm,fp0,fc,dif);
%============================================================
function dhval = HvE_o_parab(f,fH,cost,cos2,sin2,pv,Hval,hval,hm,fval)

tol = 1e-4;
dhval = quadl(@Mvh_o_parab,hval,hm,tol,0,f,fH,cost,cos2,sin2,pv,Hval,hval,fval);
%===========================================================
function [dhv_parab] = HV_theor_parab(f,fH,cost,cos2,sin2,pv,Hval,hval,fval)

for i=1:length(f)
    dhv_parab(i) = Hv_o_parab(f(i),fH,cost,cos2,sin2,pv,Hval,hval,fval);
end
%====================================================
function dhv2 = SUM_dhv_parab(par,pv,f,hv,hval,fval,fH,cost,cos2,sin2)

Hval = par;
n = length(f);
s = 0;
for i=1:n
    hvt(i) = Hv_o_parab(f(i),fH,cost,cos2,sin2,pv,Hval,hval,fval);
    s = s + ((hv(i)-hvt(i)))^2;
end
dhv2 = s;
%====================================================
function hv = Hv_o_parab(f,fH,cost,cos2,sin2,pv,Hval,hval,fval)

tol = 1e-4;
Xpd = 0.99;
fpd = f*sqrt(Xpd);
% hpd = hval+Hval*(fpd/fval-1)^(1/pv);
hpd = hval+Hval*((fpd/fval)^2-1)^(1/pv);
hvD = quadl(@Mvh_o_parab,hval,hpd,tol,0,f,fH,cost,cos2,sin2,pv,Hval,hval,fval);

% hr = hval+Hval*(f/fval-1)^(1/pv);
hr = hval+Hval*((f/fval)^2-1)^(1/pv);
Hr = (hr - hpd)/(1 - Xpd);
tm = sqrt(1-Xpd); dh = tm/6;
q = dh*(Mvr_o(0,f,fH,cost,cos2,sin2)+4*Mvr_o(tm/2,f,fH,cost,cos2,sin2)+Mvr_o(tm,f,fH,cost,cos2,sin2));
hvr = 2*Hr*q;

hv = hvD + hvr;
%===========================================================
function Mgh_o = Mvh_o_parab(h,f,fH,cost,cos2,sin2,pv,Hval,hval,fval)
% ядро гр. пок. прел. для интегрирования по высоте;

Y = fH/f;
YL = Y*cost;
% fp = fval*(1+((h-hval)./Hval).^pv);
% X = (fp/f).^2;
fp2 = fval^2*(1+((h-hval)./Hval).^pv);
X = fp2./f^2;
Mgh_o = Mv_o(X,cost,cos2,sin2,Y,YL)./sqrt(1-X);
%============================================================




