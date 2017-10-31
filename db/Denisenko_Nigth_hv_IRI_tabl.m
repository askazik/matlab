function [hv_IRI] = Denisenko_Nigth_hv_IRI_tabl(f,h_N,f_N)

fH =1.29;
tetta = 25.19;
cost = cos(pi*tetta/180);
cos2 = cost^2;
sin2 = 1 - cos2;

[hv_IRI] = HV_theor_tabl(f,fH,cost,cos2,sin2,h_N,f_N);

%====================================================
function [hv_calc] = HV_theor_tabl(f,fH,cost,cos2,sin2,h_N,fN_tab)

for i=1:length(f)
    hv_calc(i) = Hv_o_tabl(f(i),fH,cost,cos2,sin2,h_N,fN_tab);
end
%====================================================
function hv = Hv_o_tabl(f,fH,cost,cos2,sin2,h_tab,fN_tab)

tol = 1e-6;
h0 = h_tab(1);
Xpd = 0.99;
fpd = f*sqrt(Xpd);
hpd = H_refl_tabl(fpd,h_tab,fN_tab);

qD = quadl(@Mvh_o_tabl,h0,hpd,tol,0,f,fH,cost,cos2,sin2,h_tab,fN_tab);
hvD = qD;

hr = H_refl_tabl(f,h_tab,fN_tab);
Hr = (hr - hpd)/(1 - Xpd);
tm = sqrt(1-Xpd); dh = tm/6;
q = dh*(Mvr_o(0,f,fH,cost,cos2,sin2)+4*Mvr_o(tm/2,f,fH,cost,cos2,sin2)+Mvr_o(tm,f,fH,cost,cos2,sin2));
hvr = 2*Hr*q;

hv = h0 + hvD + hvr;
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



