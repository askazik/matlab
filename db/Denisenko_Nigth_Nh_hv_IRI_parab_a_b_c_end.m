function [fig_h, sig_E] = Denisenko_Nigth_Nh_hv_IRI_parab_a_b_c_end(f, hv, f_N, h_N, ks_title, key_sig_E)
% key_sig_E = 0 - F область не считается

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

% ===== Поиск части таблицы ниже высоты отражения частоты fF_min = f(1);
fF_min = f(1);
for i=1:n-1
    if f_N(i)<=fF_min&&fF_min<=f_N(i+1)
        nF_min = i+1;
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

% ===== Поиск fval и fval по таблице.
for i=nE_min:nF_min-2
    if f_N(i)>=f_N(i+1)&&f_N(i+1)<=f_N(i+2)
        nE_val = i+1;
        hval = h_N(i+1);
        fval = f_N(i+1);
        break
    end
end

h_NE = h_N(1:nE_val+1); % Часть таблицы высот до глубины долины;
f_NE = f_N(1:nE_val+1); % Часть таблицы пл. частот до глубины долины;
h0 = h_N(1);

% ===== Расчет вклада в действ. высоты dhvE ниже hval. 
% ===== Определение вклада в действ. высоты dhvF выше hval.
for i=1:nf
    dhvE(i) = HvE_o_tabl(f(i),fH,cost,cos2,sin2,h0,hoE,hval,h_NE,f_NE);
    dhvF(i) = hv(i) - dhvE(i);
end

fE = f(1:nE);          % Частоты для расчета параболич. профиля.
dhv_E = dhvF(1:nE);    % Вклад в действ. высоты выше hval для расчета параболич. профиля. 
d_hvE = dhvE(1:nE);    % Вклад в действ. высоты ниже hval для параболич. профиля.
fF = f(nE+1:1:nf);       % Частоты для расчета модельного. профиля.
dhvF = dhvF(nE+1:1:nf);  % Вклад в действ. высоты выше параболич. профиля.
dhvEF = dhvE(nE+1:1:nf); % Вклад в действ. высоты ниже hval для модельного профиля.
nF = length(fF);

% ====== Расчет полутолщины Hval параболич. профиля.
tol = 1e-4;
for i=1:nE
    t0 = 0;
    tr = pi/2;
    A(1,i) = (fE(i)/fval)*quadl(@Mvh_o_parab,t0,tr,tol,0,fE(i),fH,cost,cos2,sin2,fval);
end
y = dhv_E';
B = A*A';
Hval = inv(B)*A*y;
dhv_parab = A*Hval;
Sval = (dhv_parab-dhv_E)*(dhv_parab-dhv_E)';
sig_E = sqrt(Sval/(nE));
% =========================================================================

if(key_sig_E==1)
    
    hvE_calc = d_hvE + dhv_parab; % Расчетные значения действ. высот для параболич. профиля.
    fpF = fE(nE);                             % Плазм. частота сшивки параболы и модели.
    hpF = hval + Hval*((fpF/fval)^2-1)^(1/2); % Высота сшивки параболы и модели.
    dif = (fval^2/fpF)*(hpF-hval)/Hval/Hval;  % Производная dh/dfN для сшивки параболы и модели.
    dh = (hpF-hval)/20;
    hp_E = hval:dh:hpF;
    fp_E = fval*sqrt((1+((hp_E-hval)./Hval).^2));
    
    % Расчет вклада параболы в действ. высоты.
    dhv_parab = [];
    for i=1:nF
        dhv_parab(i) = HvE_o_parab(fF(i),fH,cost,cos2,sin2,hval,hpF,fval,Hval);
        dhv_F(i) = dhvF(i) - dhv_parab(i);
    end
    
    % Расчет параметров модели.
    options = optimset('TolX',1e-4);
    m = 2; c = 0.75;
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
    hvF_calc = dhv_calc + dhvEF + dhv_parab; % Расчетные значения действ. высот для модельного профиля.
    hv_calc = [hvE_calc,hvF_calc];
    fd = [fE,fF];
    dh = (hmF-hpF)/20;
    hp_F = hpF+dh:dh:hmF;
    fp_F = fN_model(hp_F,c,m,hpF,hmF,fpF,foF,dif);
    hp = [hp_E,hp_F];
    fp = [fp_E,fp_F];
    
    fig_h = figure('NumberTitle','off','Name',['nE=',num2str(nE),'  Hv=',num2str(Hval),'  sE=',num2str(sig_E),...
        '  sF=',num2str(sig_F),'  foE=',num2str(foE),'  fv=',num2str(fval),'  fpF=',num2str(fpF),...
        '  foF=',num2str(foF),'  hmE=',num2str(hoE),'  hv=',num2str(hval),'  hpF=',num2str(hpF),...
        '  hmF=',num2str(hmF),'  a=',num2str(a),'  b=',num2str(b),'  c=',num2str(c)]);
    plot (f_N,h_N,'-b',fp,hp,'-r',f,hv,'-*k',fd,hv_calc,'-^r','LineWidth',2);
    ylabel('h & h_v, km')
    xlabel('f & f_N, MHz')
    legend('h(IRI)','h_c_o_r(h_v)','h_v(exper)','h_v(calc)','Location','southeast')
    title(ks_title)
    grid on
    
end % if key_sig_F

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
function dhval = HvE_o_parab(f,fH,cost,cos2,sin2,hval,hm,fval,Hval)

t0 = 0;
xb = (hm-hval)/Hval;
a = sqrt((f/fval)^2-1);
tb = asin(xb/a);
dhval = Hv_o_parab(t0,tb,f,fH,cost,cos2,sin2,fval,Hval);
%===========================================================
function [dhv_parab] = HV_theor_parab(f,fH,cost,cos2,sin2,fval,Hval)

t0 = 0;
tr = pi/2;
for i=1:length(f)
    dhv_parab(i) = Hv_o_parab(t0,tr,f(i),fH,cost,cos2,sin2,fval,Hval);
end
%====================================================
function dhv2 = SUM_dhv_parab(par,f,hv,fval,fH,cost,cos2,sin2)

Hval = par;
t0 = 0;
tr = pi/2;
n = length(f);
s = 0;
for i=1:n
    hvt(i) = Hv_o_parab(t0,tr,f(i),fH,cost,cos2,sin2,fval,Hval);
    s = s + ((hv(i)-hvt(i)))^2;
end
dhv2 = s;
%====================================================
function hv = Hv_o_parab(t0,tr,f,fH,cost,cos2,sin2,fval,Hval)

tol = 1e-4;
hv = (f/fval)*Hval*quadl(@Mvh_o_parab,t0,tr,tol,0,f,fH,cost,cos2,sin2,fval);
%===========================================================
function Mgh_o = Mvh_o_parab(t,f,fH,cost,cos2,sin2,fval)
% ядро гр. пок. прел. для интегрирования по высоте;

Y = fH/f;
YL = Y*cost;
a2 = (f/fval)^2-1;
X = (fval/f)^2*(1+a2*sin(t).^2);
Mgh_o = Mv_o(X,cost,cos2,sin2,Y,YL);
%============================================================




