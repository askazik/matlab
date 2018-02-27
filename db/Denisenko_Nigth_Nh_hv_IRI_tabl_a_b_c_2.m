function [fig_h, sig_E, hmF] = Denisenko_Nigth_Nh_hv_IRI_tabl_a_b_c_2(f, hv, f_N, h_N, ks_title, key_sig_E)

% Ниже hmE таблица IRI;
% Выше hmE в долине масштаб. таблицы IRI до hpF;
% Выше долины модель fN=fmF*(a*x)^m/(1+b*x+c*x^2)^m, x = (hmF-h)/(hmF-hpF);

fig_h = 0;
fH =1.29;
tetta = 25.19;
cost = cos(pi*tetta/180);
cos2 = cost^2;
sin2 = 1 - cos2;

% tmp = flipud(trace_F2o);
% f = tmp(:,1)';
% hv = tmp(:,2)';

% nE - Число частот и действ. высот для расчета параболич. профиля.
if hv(3)<hv(1)
    nE = 3;
else
    nE = 1;
end

% tmp = iri_profile;
% hu_tab = tmp(:,2)';
% fu_tab = tmp(:,1)';
hu_tab = h_N;
fu_tab = f_N;

ind = find(fu_tab == max(fu_tab));
hu_tab = hu_tab(1:ind);
fu_tab = fu_tab(1:ind);

n = length(hu_tab);
h_N = hu_tab; % высоты из IRI
f_N = fu_tab; % плазм частоты IRI
foF = f_N(n);
hmF = h_N(n);
nf = length(f);

fE = f(1:nE);          % Частоты для расчета параболич. профиля.
fF = f(nE+1:nf);       % Частоты для расчета модельного. профиля.

% ===== Поиск части таблицы ниже высоты отражения частоты fF_min = fF(1);
fF_min = fF(1);
for i=1:n-1
    if f_N(i)<=fF_min&&fF_min<=f_N(i+1)
        nF_min = i+1;
        break
    end
end

% ===== Поиск foE и hoE по таблице.
for i=1:nF_min-2
    if f_N(i)<=f_N(i+1)&&f_N(i+1)>=f_N(i+2)
        hoE = h_N(i+1);
        foE = f_N(i+1);
        nE_min = i+1;
        break
    end
end

% ===== Поиск hval и fval по таблице.
for i=nE_min:nF_min-2
    if f_N(i)>=f_N(i+1)&&f_N(i+1)<=f_N(i+2)
        hval = h_N(i+1);
        fval = f_N(i+1);
        break
    end
end

h_NE = h_N(1:nE_min+1); % Часть таблицы высот до hmE;
f_NE = f_N(1:nE_min+1); % Часть таблицы пл. частот до hmE;
h0 = h_N(1);

h_NV = h_N(nE_min:nF_min); % Часть таблицы высот выше hmE;
f_NV = f_N(nE_min:nF_min); % Часть таблицы пл. частот выше hmE;


% ===== Расчет вклада в действ. высоты dhvE ниже hval. 
% ===== Определение вклада в действ. высоты dhvF выше hval.
for i=1:nf
    dhvE(i) = HvE_o_tabl(f(i),fH,cost,cos2,sin2,h0,hoE,h_NE,f_NE);
    dhvF(i) = hv(i) - dhvE(i);
end

dhv_E = dhvF(1:nE);    % Вклад в действ. высоты выше hmE для масштабирования. 
d_hvE = dhvE(1:nE);    % Вклад в действ. высоты ниже hmE для масштабирования.
dhvF = dhvF(nE+1:nf);  % Вклад в действ. высоты выше hmE для модельного профиля.
dhvEF = dhvE(nE+1:nf); % Вклад в действ. высоты ниже hmE для модельного профиля.

% ====== Расчет долины масштабированием hmE<h<hpF.

options = optimset('TolX',1e-6);

hm0 = h_NV(length(h_NV));
[hmV,S] = fminsearch(@SUM_dhv_tabl,hm0,options,fE,dhv_E,fH,cost,cos2,sin2,h_NV,f_NV);
[hV_cor] = COR(hmV,h_NV);
sig_E = sqrt(S/(length(fE)));
[hv_V_calc] = HV_theor_tabl(fE,fH,cost,cos2,sin2,hV_cor,f_NV)+d_hvE;

if(key_sig_E==1)
    % Расчет вклада долины в действ. высоты.
    fpF = fE(nE);
    hpF = H_refl_tabl(fpF,hV_cor,f_NV); % Высота сшивки долины и модели.
    for i=1:length(fF)
        dhvV(i) = HvV_o_tabl(fF(i),fH,cost,cos2,sin2,hoE,hpF,hV_cor,f_NV);
        dhv_F(i) = dhvF(i) - dhvV(i);
    end
    % =========================================================================

    hvV_calc = dhvEF + dhvV;   % Расчетные значения действ. высот до hpF.

    hr = H_refl_tabl(fF(1),hV_cor,f_NV);
    dif = (fF(1)-fpF)/(hr-hpF);

    % Расчет параметров модели.
    options = optimset('TolX',1e-4);
    m = 2; c = 0.75;
    par = [c,hmF];
    [par,SF] = fminsearch(@SUM_dhv_model,par,options,fF,dhv_F,hpF,fpF,foF,dif,m,fH,cost,cos2,sin2);
    c = par(1);
    hmF = par(2);
    sig_F = sqrt(SF/(length(fF)));

    fip = (log(foF/fpF))^(1/m);
    fipr = (log(foF/fpF))^(1/m-1)/m*(dif*(hmF-hpF)/fpF);
    b = -(-fip+fipr+c*(fip+fipr))/fipr;
    a = (1+b+c)*fip;

    [dhv_calc] = HV_theor_model(fF,fH,cost,cos2,sin2,c,m,hpF,hmF,fpF,foF,dif);
    hvF_calc = dhv_calc + hvV_calc; % Расчетные значения действ. высот для модельного профиля.
    hv_calc = [hv_V_calc,hvF_calc];
    fd = [fE,fF];

    dh = (hpF-hoE)/10;
    hp_V = hoE:dh:hpF;
    fp_V = fN_tabl(hp_V,hV_cor,f_NV);

    dh = (hmF-hpF)/20;
    hp_F = hpF+dh:dh:hmF;
    fp_F = fN_model(hp_F,c,m,hpF,hmF,fpF,foF,dif);
    hp = [hp_V,hp_F];
    fp = [fp_V,fp_F];

    %[hv_IRI] = HV_theor_tabl(f,fH,cost,cos2,sin2,h_N,f_N)+h_N(1);

    fig_h = figure('NumberTitle','off','Name',['VAL_2  ','nE=',num2str(nE),'  sE=',num2str(sig_E),...
        '  sF=',num2str(sig_F),'  foE=',num2str(foE),'  fv=',num2str(fval),'  fpF=',num2str(fpF),...
        '  foF=',num2str(foF),'  hmE=',num2str(hoE),'  hv=',num2str(hval),'  hpF=',num2str(hpF),...
        '  hmF=',num2str(hmF),'  a=',num2str(a),'  b=',num2str(b),'  c=',num2str(c)]);
    
    plot (f_N,h_N,'-b','LineWidth',2,'DisplayName','h(IRI)');
    hold on
    plot (fp,hp,'-r','LineWidth',2,'DisplayName','h_c_o_r(h_v)');
    plot (f,hv,'-*k','LineWidth',2,'DisplayName','h_v(exper)');
    plot (fd,hv_calc,'-^r','LineWidth',2,'DisplayName','h_v(calc)');
    %plot (f,hv_IRI,'-ob','LineWidth',2);
    
    ylabel('h & h_v, km')
    xlabel('f & f_N, MHz')
    title(ks_title)
    grid on

end % if key_sig_F

%===========================================================
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
% dhm = (hm - hm0)/(n-1);
% dh = dhm*((h_N - h0)./(hm0 - h0));
% h_tab = h_N + dh ;
%====================================================
function [hv_calc] = HV_theor_tabl(f,fH,cost,cos2,sin2,h_N,fN_tab)

for i=1:length(f)
    hv_calc(i) = Hv_o_tabl(f(i),fH,cost,cos2,sin2,h_N,fN_tab);
end
%====================================================
function dhv2 = SUM_dhv_tabl(hm,f,hv,fH,cost,cos2,sin2,h_N,fN_tab)

[h_tab] = COR(hm,h_N);
n = length(f);
s = 0;
for i=1:n
    hvt(i) = Hv_o_tabl(f(i),fH,cost,cos2,sin2,h_tab,fN_tab);
    s = s + ((hv(i)-hvt(i)))^2;
end
dhv2 = s;
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

hv = hvD + hvr;
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
function dhvE = HvE_o_tabl(f,fH,cost,cos2,sin2,h0,hoE,h_NE,f_NE)

tol = 1e-4;
qE1 = quadl(@Mvh_o_tabl,h0,hoE,tol,0,f,fH,cost,cos2,sin2,h_NE,f_NE);
dhvE = h0+qE1; 
%===========================================================
function dhvV = HvV_o_tabl(f,fH,cost,cos2,sin2,h0,hoE,h_NE,f_NE)

tol = 1e-4;
dhvV = quadl(@Mvh_o_tabl,h0,hoE,tol,0,f,fH,cost,cos2,sin2,h_NE,f_NE);
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




