function [Ne,h,outf,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                                        foF2, hmF2, foF1, hmF1, foE, hmE)
% Применение модели IRI2016
% Корректировка модели по имеющимся foF2, foF1, foE.
% Высоты не отрабатываются!!!!!!!!!!!!!!!

jmag = 0; % 0/1,geog/geom
%    JF switches to turn off/on (.true./.false.) several options
%
%    i       .true.                  .flase.          standard version
%    -----------------------------------------------------------------
%    1    Ne computed            Ne not computed                     t
%    2    Te, Ti computed        Te, Ti not computed                 t
%    3    Ne & Ni computed       Ni not computed                     t
%    4    B0 - Table option      B0 - Gulyaeva (1987)                t
%    5    foF2 - CCIR            foF2 - URSI                     false
%    6    Ni - DS-78 & DY-85     Ni - DS-95 & TTS-03             false
%    7    Ne - Tops: f10.7<188   f10.7 unlimited                     t            
%    8    foF2 from model        foF2 or NmF2 - user input           t
%    9    hmF2 from model        hmF2 or M3000F2 - user input        t
%   10    Te - Standard          Te - Using Te/Ne correlation        t
%   11    Ne - Standard Profile  Ne - Lay-function formalism         t
%   12    Messages to unit 6     no messages                         t
%   13    foF1 from model        foF1 or NmF1 - user input           t
%   14    hmF1 from model        hmF1 - user input (only Lay version)t
%   15    foE  from model        foE or NmE - user input             t
%   16    hmE  from model        hmE - user input                    t
%   17    Rz12 from file         Rz12 - user input                   t
%   18    IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
%   19    F1 probability model   critical solar zenith angle (old)   t
%   20    standard F1            standard F1 plus L condition        t
%   21    ion drift computed     ion drift not computed          false
%   22    ion densities in %     ion densities in m-3                t
%   23    Te_tops (Aeros,ISIS)   Te_topside (Intercosmos)        false
%   24    D-region: IRI-95       Special: 3 D-region models          t
%   25    F107D from AP.DAT      F107D user input (oarr(41))         t
%   26    foF2 storm model       no storm updating                   t
%   27    IG12 from file         IG12 - user input					 t
%   28    spread-F probability 	 not computed                    false
%   29    IRI01-topside          new options as def. by JF(30)   false
%   30    IRI01-topside corr.    NeQuick topside model   	     false 
%     (29,30) = (t,t) IRIold, (f,t) IRIcor, (f,f) NeQuick, (t,f) TTS   
%   ------------------------------------------------------------------
%
%  Depending on the jf() settings additional INPUT parameters may 
%  be required:
%
%       Setting              INPUT parameter
%    -----------------------------------------------------------------
%    jf(8)  =.false.     OARR(1)=user input for foF2/MHz or NmF2/m-3
%    jf(9)  =.false.     OARR(2)=user input for hmF2/km or M(3000)F2
%    jf(10 )=.false.     OARR(15),OARR(16)=user input for Ne(300km),
%       Ne(400km)/m-3. Use OARR()=-1 if one of these values is not 
%       available. If jf(23)=.false. then Ne(300km), Ne(550km)/m-3.
%    jf(13) =.false.     OARR(3)=user input for foF1/MHz or NmF1/m-3 
%    jf(14) =.false.     OARR(4)=user input for hmF1/km
%    jf(15) =.false.     OARR(5)=user input for foE/MHz or NmE/m-3 
%    jf(16) =.false.     OARR(6)=user input for hmE/km
%    jf(17) =.flase.     OARR(33)=user input for Rz12
%    jf(21) =.true.      OARR(41)=user input for daily F10.7 index
%    jf(23) =.false.     OARR(41)=user input for daily F10.7 index
%    jf(24) =.false.     OARR(41)=user input for daily F10.7 index
%          optional for jf(21:24); default is F10.7D=COV
%    jf(25) =.false.     OARR(41)=user input for daily F10.7 index
%          if oarr(41).le.0 then 12-month running mean is 
%          taken from internal file]
%    jf(27) =.flase.     OARR(39)=user input for IG12
%    jf(21) =.true.      OARR(41)=user input for daily F10.7 index
% --------------------------------------------------------------------
% Standard: t,t,t,t,f,f,t,t,t,t,t,t,t,t,t,t,t,t,t,t,f,t,f,t,t,t,t,f,f,f
% --------------------------------------------------------------------
%jf = logical([1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,0,0,0]);

% Без корректировки
jf = [1,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,0,0,0,...
     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% Тестируем на наличие данных для корректировки.
is_foF2 = ~isnan(foF2);
is_hmF2 = ~isnan(hmF2); 
is_foF1 = ~isnan(foF1); 
is_hmF1 = ~isnan(hmF1); 
is_foE = ~isnan(foE);
is_hmE = ~isnan(hmE);
% Корректируем вход по наличным данным
% c =======================================================================
% c Programs using subroutine IRI_WEB need to include (see IRITEST.FOR):
% c
% c        do i=1,100
% c          oarr(i,1)=-1.0
% c        end do
% c =======================================================================
oar = -ones(100,1);
% fN = sqrt(Ne/(1.24*10^10));
if is_foF2
    jf(8) = 0;
    oar(1,1)= foF2;
end
if is_hmF2
    jf(9) = 0;
    oar(2,1)= hmF2;    
end
if is_foF1
    jf(13) = 0;
    oar(3,1)= foF1;
end
if is_hmF1
    jf(14) = 0;
    oar(4,1)= hmF1;     
end
if is_foE
    jf(15) = 0;
    oar(5,1)= foE;      
end
if is_hmE
    jf(16) = 0;
    oar(6,1)= hmE;     
end
% alati = 30.; % lati/deg
% along = 30.; % long/deg
% iyyyy = 2010; % year(yyyy)
% mmdd = 101; % mmdd(or -ddd)
iut = 1; % iut(=0/1,LT/UT)
% dhour = 14.; % hour
% begin, end, and stepsize for the selected variable
% vbeg = 30.; % начальная высота
% vend = 1000.; % конечная высота
% vstp = 10.;  % шаг по высоте

% [outf,oarr] = iri_web(jmag,jf,alati,along,iyyyy,mmdd,iut,dhour,...
%                         vbeg,vend,vstp);
        
% height/km
% (enter  0 for list of peak heights and densities)
% (enter -1 for plasma frequencies, B0, M3000, valley, width and depth,)
% (         F1 probability, equatorial vertical ion drift, and)
% (         foF2 storm/quiet ratio, or 3 parameter of your choice)
height = 0;
h_tec_max = 0.; %upper height [km] for TEC integration (0 for no TEC)                    
ivar = 1; % variable? (1/2/../8 for height/lat/long/year/month/day/day of year/hour)

[outf,oarr] = mexIRI2016(jmag,jf,alati,along,iyyyy,mmdd,iut,dhour,...
                        height,h_tec_max,ivar,vbeg,vend,vstp,oar);
                    
% =========================================================================
% C-----------------------------------------------------------------
% C
% C INPUT:  JF(1:50)      true/false switches for several options
% C         JMAG          =0 geographic   = 1 geomagnetic coordinates
% C         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
% C         IYYYY         Year as YYYY, e.g. 1985
% C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
% C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL 
% C                          HOURS
% C         HEIBEG,       HEIGHT RANGE IN KM; maximal 100 heights, i.e.
% C          HEIEND,HEISTP        int((heiend-heibeg)/heistp)+1.le.100
% C
% C    JF switches to turn off/on (.true./.false.) several options
% C
% C    i       .true.                  .false.          standard version
% C    -----------------------------------------------------------------
% C    1    Ne computed            Ne not computed                     t
% C    2    Te, Ti computed        Te, Ti not computed                 t
% C    3    Ne & Ni computed       Ni not computed                     t
% C    4    B0,B1 - Bil-2000       B0,B1 - other models jf(31)     false
% C    5    foF2 - CCIR            foF2 - URSI                     false
% C    6    Ni - DS-1995 & DY-1985 Ni - RBV-2010 & TTS-2005        false
% C    7    Ne - Tops: f10.7<188   f10.7 unlimited                     t            
% C    8    foF2 from model        foF2 or NmF2 - user input           t
% C    9    hmF2 from model        hmF2 or M3000F2 - user input        t
% C   10    Te - Standard          Te - Using Te/Ne correlation        t
% C   11    Ne - Standard Profile  Ne - Lay-function formalism         t
% C   12    Messages to unit 6     to messages.txt on unit 11          t
% C   13    foF1 from model        foF1 or NmF1 - user input           t
% C   14    hmF1 from model        hmF1 - user input (only Lay version)t
% C   15    foE  from model        foE or NmE - user input             t
% C   16    hmE  from model        hmE - user input                    t
% C   17    Rz12 from file         Rz12 - user input                   t
% C   18    IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
% C   19    F1 probability model   critical solar zenith angle (old)   t
% C   20    standard F1            standard F1 plus L condition        t
% C   21    ion drift computed     ion drift not computed          false
% C   22    ion densities in %     ion densities in m-3                t
% C   23    Te_tops (Bil-1985)     Te_topside (TBT-2012)           false
% C   24    D-region: IRI-1990     FT-2001 and DRS-1995                t
% C   25    F107D from APF107.DAT  F107D user input (oarr(41))         t
% C   26    foF2 storm model       no storm updating                   t
% C   27    IG12 from file         IG12 - user                         t
% C   28    spread-F probability 	 not computed                    false
% C   29    IRI01-topside          new options as def. by JF(30)   false
% C   30    IRI01-topside corr.    NeQuick topside model   	     false 
% C (29,30) = (t,t) IRIold, (f,t) IRIcor, (f,f) NeQuick
% C   31    B0,B1 ABT-2009	     B0 Gulyaeva-1987 h0.5               t   
% C (4,31) = (t,t) Bil-00, (f,t) ABT-09, (f,f) Gul-87, (t,f) not used
% C   32    F10.7_81 from file     F10.7_81 - user input (oarr(46))    t
% C   33    Auroral boundary model on/off  true/false	             false
% C   34    Messages on            Messages off                        t
% C   35    foE storm model        no foE storm updating           false
% C   36    hmF2 w/out foF2_storm  with foF2-storm                     t
% C   37    topside w/out foF2-storm  with foF2-storm                  t
% C   38    turn WRITEs off in IRIFLIP   turn WRITEs on                t
% C   39    hmF2 (M3000F2)         new models                      false
% C   40    hmF2 AMTB-model        Shubin-COSMIC model                 t
% C   41    Use COV=F10.7_365      COV=f(IG12) (IRI before Oct 2015)   t
% C   42    Te with PF10.7 dep.	 w/o PF10.7 dependance               t
% C   43    B0 from model          B0 user input                       t
% C   44    B1 from model          B1 user input                       t
% C      ....
% C   50    
% C   ------------------------------------------------------------------
% C
% C  Depending on the jf() settings additional INPUT parameters may 
% c  be required:
% C
% C       Setting              INPUT parameter
% C    -----------------------------------------------------------------
% C    jf(8)  =.false.     OARR(1)=user input for foF2/MHz or NmF2/m-3
% C    jf(9)  =.false.     OARR(2)=user input for hmF2/km or M(3000)F2
% C    jf(10 )=.false.     OARR(15),OARR(16)=user input for Ne(300km),
% C       Ne(400km)/m-3. Use OARR()=-1 if one of these values is not 
% C       available. If jf(23)=.false. then Ne(300km), Ne(550km)/m-3.
% C    jf(13) =.false.     OARR(3)=user input for foF1/MHz or NmF1/m-3 
% C    jf(14) =.false.     OARR(4)=user input for hmF1/km
% C    jf(15) =.false.     OARR(5)=user input for foE/MHz or NmE/m-3 
% C    jf(16) =.false.     OARR(6)=user input for hmE/km
% C    jf(17) =.flase.     OARR(33)=user input for Rz12
% C    jf(25) =.false.     OARR(41)=user input for daily F10.7 index
% C    jf(27) =.false.     OARR(39)=user input for IG12
% C    jf(43) =.false.     OARR(10)=user input for B0
% C    jf(44) =.false.     OARR(87)=user input for B1
% C
% C
% C  OUTPUT:  OUTF(1:20,1:1000)
% C               OUTF(1,*)  ELECTRON DENSITY/M-3
% C               OUTF(2,*)  NEUTRAL TEMPERATURE/K
% C               OUTF(3,*)  ION TEMPERATURE/K
% C               OUTF(4,*)  ELECTRON TEMPERATURE/K
% C               OUTF(5,*)  O+ ION DENSITY/% or /M-3 if jf(22)=f 
% C               OUTF(6,*)  H+ ION DENSITY/% or /M-3 if jf(22)=f
% C               OUTF(7,*)  HE+ ION DENSITY/% or /M-3 if jf(22)=f
% C               OUTF(8,*)  O2+ ION DENSITY/% or /M-3 if jf(22)=f
% C               OUTF(9,*)  NO+ ION DENSITY/% or /M-3 if jf(22)=f
% C                 AND, IF JF(6)=.FALSE.:
% C               OUTF(10,*)  CLUSTER IONS DEN/% or /M-3 if jf(22)=f
% C               OUTF(11,*)  N+ ION DENSITY/% or /M-3 if jf(22)=f
% C               OUTF(12,*)  
% C               OUTF(13,*)  
% C  if(jf(24)    OUTF(14,1:11) standard IRI-Ne for 60,65,..,110km 
% C     =.false.)        12:22) Friedrich (FIRI) model at these heights 
% C                      23:33) standard Danilov (SW=0, WA=0) 
% C                      34:44) for minor Stratospheric Warming (SW=0.5) 
% C                      45:55) for major Stratospheric Warming (SW=1) 
% C                      56:66) weak Winter Anomaly (WA=0.5) conditions
% C                      67:77) strong Winter Anomaly (WA=1) conditions
% C               OUTF(15-20,*)  free
% c
% C            OARR(1:100)   ADDITIONAL OUTPUT PARAMETERS         
% C
% C      #OARR(1) = NMF2/M-3           #OARR(2) = HMF2/KM
% C      #OARR(3) = NMF1/M-3           #OARR(4) = HMF1/KM
% C      #OARR(5) = NME/M-3            #OARR(6) = HME/KM
% C       OARR(7) = NMD/M-3             OARR(8) = HMD/KM
% C       OARR(9) = HHALF/KM           #OARR(10) = B0/KM
% C       OARR(11) =VALLEY-BASE/M-3     OARR(12) = VALLEY-TOP/KM
% C       OARR(13) = TE-PEAK/K          OARR(14) = TE-PEAK HEIGHT/KM
% C      #OARR(15) = TE-MOD(300KM)     #OARR(16) = TE-MOD(400KM)/K
% C       OARR(17) = TE-MOD(600KM)      OARR(18) = TE-MOD(1400KM)/K
% C       OARR(19) = TE-MOD(3000KM)     OARR(20) = TE(120KM)=TN=TI/K
% C       OARR(21) = TI-MOD(430KM)      OARR(22) = X/KM, WHERE TE=TI
% C       OARR(23) = SOL ZENITH ANG/DEG OARR(24) = SUN DECLINATION/DEG
% C       OARR(25) = DIP/deg            OARR(26) = DIP LATITUDE/deg
% C       OARR(27) = MODIFIED DIP LAT.  OARR(28) = Geographic latitude
% C       OARR(29) = sunrise/dec. hours OARR(30) = sunset/dec. hours
% C       OARR(31) = ISEASON (1=spring) OARR(32) = Geographic longitude
% C      #OARR(33) = Rz12               OARR(34) = Covington Index
% C      #OARR(35) = B1                 OARR(36) = M(3000)F2
% C      $OARR(37) = TEC/m-2           $OARR(38) = TEC_top/TEC*100.
% C      #OARR(39) = gind (IG12)        OARR(40) = F1 probability 
% C      #OARR(41) = F10.7 daily        OARR(42) = c1 (F1 shape)
% C       OARR(43) = daynr              OARR(44) = equatorial vertical 
% C       OARR(45) = foF2_storm/foF2_quiet         ion drift in m/s
% C      #OARR(46) = F10.7_81           OARR(47) = foE_storm/foE_quiet 
% C       OARR(48) = spread-F probability          
% C       OARR(49) = Geomag. latitude   OARR(50) = Geomag. longitude  
% C       OARR(51) = ap at current time OARR(52) = daily ap
% C       OARR(53) = invdip/degree      OARR(54) = MLT-Te
% C       OARR(55) = CGM-latitude       OARR(56) = CGM-longitude
% C       OARR(57) = CGM-MLT            OARR(58) = CGM lat eq. aurl bodry
% C       OARR(59) = CGM-lati(MLT=0)    OARR(60) = CGM-lati for MLT=1
% C       OARR(61) = CGM-lati(MLT=2)    OARR(62) = CGM-lati for MLT=3
% C       OARR(63) = CGM-lati(MLT=4)    OARR(64) = CGM-lati for MLT=5
% C       OARR(65) = CGM-lati(MLT=6)    OARR(66) = CGM-lati for MLT=7
% C       OARR(67) = CGM-lati(MLT=8)    OARR(68) = CGM-lati for MLT=9
% C       OARR(69) = CGM-lati(MLT=10)   OARR(70) = CGM-lati for MLT=11
% C       OARR(71) = CGM-lati(MLT=12)   OARR(72) = CGM-lati for MLT=13
% C       OARR(73) = CGM-lati(MLT=14)   OARR(74) = CGM-lati for MLT=15
% C       OARR(75) = CGM-lati(MLT=16)   OARR(76) = CGM-lati for MLT=17
% C       OARR(77) = CGM-lati(MLT=18)   OARR(78) = CGM-lati for MLT=19
% C       OARR(79) = CGM-lati(MLT=20)   OARR(80) = CGM-lati for MLT=21
% C       OARR(81) = CGM-lati(MLT=22)   OARR(82) = CGM-lati for MLT=23
% C       OARR(83) = Kp at current time OARR(84) = magnetic declination 
% C       OARR(85) = L-value            OARR(86) = dipole moment 
% C                # INPUT as well as OUTPUT parameter
% C                $ special for IRIWeb (only place-holders)
% c-----------------------------------------------------------------------  

h = vbeg:vstp:vend;
Ne = double(outf(1,1:length(h)));
outf = double(outf);
oarr = double(oarr);
% plot(Ne,h);
%     grid on
%     title('N_e(h), 1/m^{-3}')
%     
% disp(['Geomagnetic inclination in degrees (для конечной высоты vend) = ' mat2str(double(oarr(25)))])