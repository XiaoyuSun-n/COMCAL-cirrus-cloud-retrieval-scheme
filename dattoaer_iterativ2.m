function ok = dattoaer_iterativ2( DatStr, InputDATDir, InputDNSDir, InputOZODir, OutputDir, ueko)

% Skript nach dattoaer_tropo. Dieses Skript iteriert LR um unter Wolken
% homogene Bedingungen zu finden
% gezielt Rechnungen zu wiederholen insb. sind LR355, LR532 usw hier
% Matrizen

% 
%
ok = 0;
MFile = [upper(mfilename) ': '];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameterdefinitionen
%%% kommen unten!
%%% erst Einlesen erforderlich
%%% find "Parameterdefinitionen"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LRiter = 1; %LR Iteration in Wolken ja / nein
LRiterS = 0;

FitRangeC=[15000 18000];   % 23000 - 26000
FitRangeCS=[11000 13000];   % fuer senkr pol insb. 355s
FitRangeA=[10000 11000]; % Randbedingungshoehen Klett 15000 16500
PolFitRange= FitRangeC; % [11000 11500]; % Randbedingungshoehe Depol
Angexpminus=1.2; %alpha propor Wvl^-Angexpminus
Hcalcrange = 32000;% Hoehenbereich - 32000
HcalcRaman = 20000; % fuer Raman N2
Pschwelle=1e-8;

BSR355schwelle = 2.0;
BSR355Sschwelle = 5.0;
BSR532schwelle = 2.8;
BSR532Sschwelle = 10.0;
BSR1064schwelle = 8;


val_LR355 = 35; val_LR532 = 36; val_LR1064=45; % Lidarverhaeltnisse Aer
val_LR355S = 35; val_LR532S = 35; 
val_LR355Wo = 20; val_LR532Wo = 28; val_LR1064Wo=32; % Lidarverhaeltnisse Wolke
val_LR355SWo = 20; val_LR532SWo = 28;
LRobergr = 100;    LRSobergr= 8000;
LRuntergr =  8;    LRSuntergr= 8;


eps = 0.05; % clear sky Hintergrund
itmax=50;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  wenn als selbstaendiges Programm: laden und vorbereiten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LvR2=[0.9 1 1.05 1.075 1.1 1.125 1.15 1.2 1.25 1.3 1.5 1.75 2 2.5 3 5 10];
LvR3=[0.7 1 1.3 1.4 1.5 1.6 1.75 2 2.5 3 4 5 7.5 10 20 25 99];
%LvR3=[0.7 1 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.5 1.6 1.7 2 2.5 3 4 5 7.5 10 25 50 99];

%Verzeichnisse
 OzoInDir=InputOZODir;                                 %['/lidar4/lidar/tropo/matlab/ozo/'];     
 PtuInDir=[InputOZODir(1:end-4) 'ptu/'];

%%% falls Ueberlapp
if ueko >0,
    numdate=datenum(str2num(DatStr(1:2)),str2num(DatStr(3:4)),str2num(DatStr(5:6)));
    
    %load /home/critter/matlab2/mein/matlab/ueberlappwerte2013_iris4_zpos3
    if numdate > 5186, %(13.3.2014)    
     str = ['load /home/critter/matlab2/mein/matlab/ueberlappwerte2014Maerz_iris3_zpos3neu'];
    end
    if numdate > 5512, %(02.2.2015)    
     str = ['load /home/critter/matlab2/mein/matlab/ueberlappwerte2015_Feb_iris2_zpos3'];
    end
    
    %oder nehme explizites File - dann length(ueko) >2
    if length(ueko) >2
        str =['load ' ueko];
    end
    
    eval(str)
end


%%% Einlesen der Daten
P355=[]; P532=[]; P532S=[]; P532A=[]; P532C=[]; P532SA=[]; P355SC=[];
P387=[]; P407=[]; P607=[]; P1064=[]; P607A=[]; P607C=[];
P355C=[]; P355A=[]; P355S=[]; P355SA=[]; P355SC=[]; P387A=[]; P387C=[];

%
%  load DAT file
%  
  DATInFile = [ InputDATDir DatStr '_komb.mat'];
  filestr = 1;
  
  if ~exist( DATInFile),
    DATInFile = [ InputDATDir DatStr '.mat'];
    filestr = 2;
    if ~exist( DATInFile),
    disp([MFile 'file ' DATInFile ' not found']);
    return
    end
  end

    load( DATInFile); disp(DATInFile);
    


  DNSInFile = [ InputDNSDir DatStr '.mat'];
  
  if ~exist(DNSInFile)
    error([MFile 'file ' DNSInFile ' not found']);
  else

    load(DNSInFile);
  end


RAWInFile = DATInFile;
s=strfind(RAWInFile,'/dat/');
RAWInFile(s:s+4)='/raw/';
s=strfind(RAWInFile,'_komb');
if ~isempty(s), RAWInFile=[RAWInFile(1:s-1) '.mat']; end
 
load(RAWInFile);

  
  
 
%%% Definitionen `

Sel532C=find(DATCounting==1 & abs(DATWavelength - 532e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'p');
Sel532A=find(DATCounting==0 & abs(DATWavelength - 532e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'p'); % sucht die entsprechenden Modes aus DAT-Files
Sel532SC=find(DATCounting==1 & abs(DATWavelength - 532e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 's');
Sel532SA=find(DATCounting==0 & abs(DATWavelength - 532e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 's'); % sucht die entsprechenden Modes aus DAT-Files
LenA = size(Sel532A,2); % Anzahl der DAT-Profile

Sel1064=find(DATCounting==0 & abs(DATWavelength - 1064e-9) < 1e-9 & DATHighLow =='H'); % sucht die entsprechenden Modes aus DAT-Files

Sel355C=find(DATCounting==1 & abs(DATWavelength - 355e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'p'); 
Sel355A=find(DATCounting==0 & abs(DATWavelength - 355e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'p'); 
% old KARL
if isempty(Sel355C)
 Sel355C=find(DATCounting==1 & abs(DATWavelength - 355e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'u');
end
if isempty(Sel355A)
 Sel355A=find(DATCounting==0 & abs(DATWavelength - 355e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'u');
end

Sel355SC=find(DATCounting==1 & abs(DATWavelength - 355e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 's'); 
Sel355SA=find(DATCounting==0 & abs(DATWavelength - 355e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 's'); 

Sel607C=find(DATCounting==1 & abs(DATWavelength - 607e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'u');
Sel607A=find(DATCounting==0 & abs(DATWavelength - 607e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'u');
Sel387C=find(DATCounting==1 & abs(DATWavelength - 387e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'u');
Sel387A=find(DATCounting==0 & abs(DATWavelength - 387e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'u');
Sel407C= find(DATCounting==1 & abs(DATWavelength - 407e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'u');
Sel407A= find(DATCounting==0 & abs(DATWavelength - 407e-9) < 1e-9 & DATHighLow =='H' & DATRecPol == 'u');
Sel660 = find(DATCounting==1 & abs(DATWavelength - 660e-9) < 1e-9 & DATHighLow =='H' );


% alte Version wenn noetig
Sel532L=find(DATCounting==0 & abs(DATWavelength - 532e-9) < 1e-9 & DATHighLow =='L' );
 if length(Sel532L) >0, 
    P532L=DATSignal(:,Sel532L); 
%     woh=80:121;
%     for j=1:length(Sel532L);
%      if  length(P532A) ==0; 
%        P532A=DATSignal(:,Sel532A);
%        for jj=1:length(Sel532A), P532A(:,j)=P532A(:,j)-mymean(P532A(800:end,j)); end
%      end
%       ff=mean(P532A(woh,j))./mean(P532L(woh,j));
%       P532L(:,j)=P532L(:,j).*ff;
%       [mini,wo]=min(abs(P532A(50:65,j)-P532L(50:65,j)));
%       P532A(1:50+wo-1,j)=P532L(1:50+wo-1,j);
%     end
 end
Sel607L=find(DATCounting==1 & abs(DATWavelength - 607e-9) < 1e-9 & DATHighLow =='L' );
if isempty(P607C),  P607C=DATSignal(:,Sel607C); end
for j=1:length(Sel607C), P607C(:,j)=P607C(:,j)-mymean(P607C(440:end,j)); end

if length(Sel607L) >0, P607L=DATSignal(:,Sel607L); end

% Defintion des Hoehenvektors
DATHeight = startsteptoheight(size(DATSignal,1),DATHeightStart,DATHeightStep);
Height = DATHeight(:,Sel532A);
H = Height(:,1);  Hsize=length(H);    dH=H(3)-H(2);  dHv=H; dHv(2:end)=diff(H);
HStat = 0; %%%%% H(1)-Hminsoll;
Selnocalc = find(H > Hcalcrange);

% Defintion der Zeitvektoren
DN532 = DATDateNumStart(Sel532A);     DM=DN532;  lDN532 = length(DN532);
DN355 = DATDateNumStart(Sel355A);                lDN355 = length(DN355);
DN1064 = DATDateNumStart(Sel1064);               lDN1064 = length(DN1064);
DN532S  = DATDateNumStart(Sel532SA);             lDN532S = length(DN532S);
DN355S  = DATDateNumStart(Sel355SA);             lDN355S = length(DN355S);
DN407  = DATDateNumStart(Sel407C);               lDN407 = length(DN407);
DN660  = DATDateNumStart(Sel660);                lDN660 = length(DN660);
DN607 = DATDateNumStart(Sel607A);                lDN607 = length(DN607);
DN387 = DATDateNumStart(Sel387A);                lDN387 = length(DN387);
%[datlen, wo] = max([lDN355 lDN355S lDN387 lDN407 lDN532 lDN532S lDN607 lDN660 lDN1064]);



%%%%%% zu den RAW-Signalen

if ~isempty(RAWSgnBck),
    R355A=RAWSgnBck(:,Sel355A);
    R355C=RAWSgnBck(:,Sel355C);
    R355SA=RAWSgnBck(:,Sel355SA);
    R355SC=RAWSgnBck(:,Sel355SC);
    R532A=RAWSgnBck(:,Sel532A);
    R532C=RAWSgnBck(:,Sel532C);
    R532SA=RAWSgnBck(:,Sel532SA);
    R532SC=RAWSgnBck(:,Sel532SC);
    R387A=RAWSgnBck(:,Sel387A);
    R387C=RAWSgnBck(:,Sel387C);
    R607A=RAWSgnBck(:,Sel607A);
    R607C=RAWSgnBck(:,Sel607C);
    R407A=RAWSgnBck(:,Sel407A);
    R407C=RAWSgnBck(:,Sel407C);
    R660C=RAWSgnBck(:,Sel660);
    R1064=RAWSgnBck(:,Sel1064);
end
    
    

% Definition der Blendenparamter
FieldStop=[];ZPosition=[];
if size(DATFieldStop(Sel532A)) > 0, FieldStop=DATFieldStop(Sel532A); end
if exist ('DATZPosition','var') ==1, 
    ZPosition=DATZPosition(Sel532A);  
else  ZPosition=0;
end


% 
% DNSHeight=zeros(size(DNSDensity));
% DNSsteps=size(DNSDensity,1);
% for jj=1:length(DNSHeightStep),
% DNSHeight(1,jj)=DNSHeightStart(jj);  
% for ji=2:DNSsteps,
%     DNSHeight(ji,jj)=DNSHeight(ji-1,jj)+DNSHeightStep(jj); 
% end
% end
% 
% DNSHeightgut=zeros(size(DATSignal,1),1); DNSDensitygut=DNSHeightgut; DNSTemperaturegut=DNSHeightgut;
% for jj=1: size(DNSHeightgut,2),
%     wo = find(~isnan(DNSDensity(:,jj)) & ~isnan(DNSTemperature(:,jj)));
%     DNSDensitygut=interp1(DNSHeight(wo,jj), DNSDensity(wo,jj),H);
%     DNSTemperaturegut=interp1(DNSHeight(wo,jj), DNSTemperature(wo,jj),H);
% end
% 
% DATDensity =zeros(size(DATSignal)); DATTemperature=DATDensity;
% dimen=size(DATSignal);
% for jj=1:dimen(2),
%     [q,w]=min(abs(DATDateNumStart(jj) - DNSDateNum)); 
%     if length(w) >1, w=w(1); end
%     DATDensity(:,jj)=DNSDensitygut(1:Hsize,1);   %DNSDensity(1:Hsize,w);
%     DATTemperature(:,jj)=DNSTemperaturegut(1:Hsize,1);
% end
% clear DNS*


% Defintion der Radiosondengroessen
dimen=size(P355);    
feu = zeros(size(P355)); mischvv=zeros(size(P355)); massmisch=zeros(size(P355)); ballonzeit=zeros(1, dimen(2));
Temp = zeros(size(P355)); Pressure = zeros(size(P355));  Density = zeros(size(P355));
for j=1 :dimen(2),
[feu(:,j),mischvv(:,j),massmisch(:,j),ballonzeit(j), Temp(:,j), Pressure(:,j), Density(:,j)] = ...
    hoppsonde(H,DN532(j),PtuInDir);

% Dichtekontrolle
   q=real(log(Density(:,j)));
   ee=find(~isnan(Density(:,j))); ee1=ee(end)+1;
   if H(ee1) < Hcalcrange,
    if j==1, disp('Radiosonde bricht zu fr?h ab - interpolieren!'); end
    qm=mymean(diff(q(ee1-11:ee1-1,1)));
    for j1=ee1:dimen(1), q(j1)=q(j1-1)+qm; end
    Temp(ee1:dimen(1),j)=Temp(ee1-1,j);
    Density(:,j)=exp(q);
   end
end


% Saettigungsdampfdruck ueber Eis
% Goff Gratch Gleichung 
TT=Temp; af=log(10);
q1=273.16./TT;  q2=TT./273.16;
  lp1=-9.09718.*(q1 -1);
  lp2=-3.56654.*(log(q1)./af);
  lp3=0.876793.*(1-TT./273.16);
  lp4=log(6.1071)./af;
  
  lp=lp1+lp2+lp3+lp4;
  Sattdruckeis=(10.^(lp)).*100;  %%% in Pa
  rfeis=mischvv.*Pressure./Sattdruckeis;

% Ozone
O3Density=zeros(Hsize,1); O3DensityAE=zeros(Hsize,1); 
tpos=floor(length(DATDateNumStart)/2);
[O3Density,O3DensityAE]=hoppozon(H,DATDateNumStart(tpos),OzoInDir);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ueberlapp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('HUeberlapp'),
    P355Aur = P355A; P355Cur = P355C; P355SAur = P355SA; P355SCur = P355SC;
    P532Aur = P532A; P532Cur = P532C; P532SAur = P532SA; P532SCur = P532SC;
    P607Aur = P607A; P607Cur = P607C; P387Aur = P387A;   P387Cur = P387C;
    P407Aur = P407A; P407Cur = P407C; P660ur = P660; P1064ur = P1064;
    U532 = interp1(HUeberlapp, Ueberlapp532,H);
    U532S = interp1(HUeberlapp, Ueberlapp532S,H);
    U355 = interp1(HUeberlapp, Ueberlapp355,H);
    U355S = interp1(HUeberlapp, Ueberlapp355S,H);
    U1064 = interp1(HUeberlapp, Ueberlapp1064,H);
    U387A = interp1(HUeberlapp, Ueberlapp387A,H);
    U407A = interp1(HUeberlapp, Ueberlapp407A,H);
    U607A = interp1(HUeberlapp, Ueberlapp607A,H);
    %U607A = interp1(HUeberlapp, Ueberlapp607A_3,H); % fur 2014!
%     U387A=(U387A-1).*0.65+1;  %F?r Jan 2019
%     U607A=(U607A-1).*1.2+1;
    if exist('Ueberlapp660'),
        U660 = interp1(HUeberlapp, Ueberlapp660,H);
    else 
        U660=ones(length(H),1);
    end
    
    for j=1:length(Sel532C),
     P532C(:,j) = P532C(:,j)./U532;
    end
    for j=1:length(Sel532A),
     P532A(:,j) = P532A(:,j)./U532;
    end
    for j=1:length(Sel532C),
     P532(:,j) = P532(:,j)./U532;
    end
    for j=1:length(Sel532SC),
     P532SC(:,j) = P532SC(:,j)./U532S;
    end
    for j=1:length(Sel532SA),
     P532SA(:,j) = P532SA(:,j)./U532S;
    end
    for j=1:length(Sel532SC),
     P532S(:,j) = P532S(:,j)./U532S;
    end
    
    for j=1:length(Sel355C),
     P355C(:,j) = P355C(:,j)./U355;
    end
    for j=1:length(Sel355C),
     P355A(:,j) = P355A(:,j)./U355;
    end
    for j=1:length(Sel355C),
     P355(:,j) = P355(:,j)./U355;
    end
    for j=1:length(Sel355SC),
     P355SC(:,j) = P355SC(:,j)./U355S;
    end
    for j=1:length(Sel355SA),
     P355SA(:,j) = P355SA(:,j)./U355S;
    end
    for j=1:length(Sel355SC),
     P355S(:,j) = P355S(:,j)./U355S;
    end
    
    for j=1:length(Sel1064),
     P1064(:,j) = P1064(:,j)./U1064;
    end
    
    for j=1:length(Sel387A),
     P387A(:,j) = P387A(:,j)./U387A;
    end
    
    for j=1:length(Sel387C),
     P387C(:,j) = P387C(:,j)./U387A;
     P387(:,j) = P387(:,j)./U387A;
    end
    
    for j=1:length(Sel407A),
     P407A(:,j) = P407A(:,j)./U407A;
     P407(:,j) = P407(:,j)./U407A;
    end
     
     for j=1:length(Sel607A),
      P607A(:,j) = P607A(:,j)./U607A;
      P607(:,j) = P607(:,j)./U607A;
     end
    
    for j=1:length(Sel660),
     P660(:,j) = P660(:,j)./U660;
    end
    
else      % kein Ueberlapp? Dann mu? wenigstens wasserfaktoren bekannt
    uvwasserfaktor = NaN;
    viswasserfaktor = NaN;
    
end   % exist HUeberlapp
%%%%%%%%%%%%%%%%%%%% Ende Ueberlapp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  weiter dimensionsvariable Parameterdefinitionen 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LR532arr=ones(size(P532)).*val_LR532;
LR532Sarr=ones(size(P532S)).*val_LR532S;
LR355arr=ones(size(P355)).*val_LR355;
LR355Sarr=ones(size(P355S)).*val_LR355S;
LR1064arr=ones(size(P1064)).*val_LR1064;
  

BSRAtFit532=1+eps;
BSRAtFit355 = 1+(2./3).^(4-Angexpminus).*eps;
BSRAtFit1064 = 1+(2).^(4-Angexpminus).*eps;
BSRAtFit532arr=zeros(dimen(2),1)+BSRAtFit532;
BSRAtFit532Sarr=zeros(dimen(2),1)+BSRAtFit532;
BSRAtFit355arr=zeros(dimen(2),1)+BSRAtFit355;
BSRAtFit355Sarr=zeros(dimen(2),1)+BSRAtFit355;
BSRAtFit1064arr=zeros(dimen(2),1)+BSRAtFit1064;
BSRAtFit355err= eps./2; % + zeros(dimen(2),1);
BSRAtFit532err= eps; % +zeros(dimen(2),1);
BSRAtFit1064err= eps.*2; %+ zeros(dimen(2),1);


Hmin=find(H>650); Hmin1=Hmin(1);
wohint=find(H>1200 & H< 1800);  % wird bei 532P endg?ltig fixiert

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Klett Rechnungen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rechnung fuer 532nm
Wvl532=5.3207e-7; WvlR532= ramshift(Wvl532,'N2');
LR532arrerr=ones(size(LR532arr))+15;    


Wvlfaktor= 1+(Wvl532./WvlR532).^Angexpminus;
Wvlfaktor2= 1-(Wvl532./WvlR532).^Angexpminus;
anfang=find(H>750); anfang=anfang(1);


%%% Berechnung von Alpha532 mit der Raman-Methode

AlRay532 = Density .* raytotwq( Wvl532, Temp, Density);
AlRay607 = Density .* raytotwq( WvlR532, Temp, Density);
WvlFct=-Angexpminus; HSel=find(H>0 & H<HcalcRaman); HSell=length(HSel);
Frqtreat = [0.2, 0.2, 60, 10];
Ptest = P607;
Mixx = [Wvl532, WvlR532]; 
ff=(Wvl532./WvlR532).^(-WvlFct)+1;
AlphaAer532Raman=zeros(size(P607A));  
AlphaAer532Ramanerr=zeros(size(P607A));
Alm532x=1e-9.*ones(dimen(1),1);
AOD532=zeros(size(Ptest));
AODVIS = zeros(size(Ptest)); %532+607
iter=8;

for jj=1:length(Sel607A),
TauRay532=qdrupvar(H,AlRay532(:,jj)+AlRay607(:,jj));    
ARayad=O3Density.*o3abswq(Wvl532,Temp(:,jj)); ARamad=O3Density.*o3abswq(WvlR532,Temp(:,jj));
AlRay532(:,jj)=AlRay532(:,jj)+ARayad; AlRay607(:,jj)=AlRay607(:,jj)+ARamad;     
Aterm2 = exp(qdrupvar(H,AlRay532(:,jj)+AlRay607(:,jj)));
XX= Density(:,jj)./H.^2 ./Aterm2;
%P607(:,jj) = P607(:,jj);   % -mymean(P607A(Hgross,jj));   
q1=Ptest(:,jj).*H.^2./Density(:,jj);
q2=poissonmittel(q1,5);
q3=q2.*0.6 + 0.4.*q1;
PR(:,jj)=q3./H.^2.*Density(:,jj); 
PRerr=real(sqrt(PR(:,jj)))./50+1e-4;
PRx(:,jj) = PR(:,jj); 

for j=1:iter,
    
[ Alm532, Alm532AE, dAdSgn, dAdRho] = alphaaer( H(HSel), HStat, ...        
   PR(HSel,jj), PRerr(HSel), Mixx, [], Density(HSel,jj), ...
   Temp(HSel, jj),Frqtreat, WvlFct, ...
   O3Density(HSel), [], 'ableit3'); %diffquotneu savitzkygolay
Alm532x(HSel) = Alm532;
wo=find(isnan(Alm532x)); Alm532x(wo)=1e-9;
PR(:,jj) = XX.*exp(-qdrupvar(H,ff.*Alm532x));
q2=find((~isnan(PR(:,jj)))); %% NaNs am Anfang herausschneiden
%anfang=q2(1); anfang=100;
PR(:,jj)=PR(:,jj)./PR(anfang,jj).*PRx(anfang,jj); 
%keyboard
end % iter
AlphaAer532Raman(HSel,jj)=Alm532;
AlphaAer532Ramanerr(HSel,jj)=Alm532AE;
AOD532(:,jj) = (real(log(Density(:,jj)./H.^2./PR(:,jj)))-TauRay532)./ff;
AODVIS(:,jj) = (real(log(Density(:,jj)./H.^2./PR(:,jj)))-TauRay532);
end
P607final = PR;


% PR=P607A;
% PRAE=real(sqrt(PR))./50*1e-4;
% s0 = 600; st=60;
% [ Alpha532schaetz, Alpha532schaetzMin, Alpha532schaetzMax, P607x, Hx607] = ...
%  alphaschaetz( H, HStat, PR, PRAE, Mixx, Density, Temp, O3Density, WvlFct, s0, st);
%  



%%% 532nm parallel polarisiert

EmitPol='p';RecPol='p'; 
BeRa532=Density.*raybckwq (Wvl532,'p','p', Temp, Density);  %H2=H-HStat;
phsoll =0.125./4094; 

% Rueckstreuverhaeltnis nach Quotientenmethode
BSR532quot=P532A; BSR532urquot=P532A;
for j=1: size(P532A,2),
    [m,w]=min(abs(DATDateNumStart(Sel607C) - DM(j)));
    BSR532quot(:,j)=P532(:,j)./PR(:,w);
end

% Signaldefinition
PC=P532C;
PA=P532A;  
P = PC; % oder P532
% ber=151:180;
nn=3; %length(ber);
f1=((1:nn)./nn)';
f2=((nn:-1:1)./nn)';
for j=1: length(Sel532A)
  %hwo=150:300; %find(PC(:,j) <10 & PC(:,j)>0.1);
  hwo=find(PC(:,j) <10 & PC(:,j)>0.1 & PA(:,j) >0.1 & H>1500);
  [c,k]=minanpass(PA(:,j),PC(:,j),hwo,H,5);
  if k > 0;
  PA2=(PA(:,j)-c)./k;
  [x,wo]=min(abs(PA2(hwo)-PC(hwo,j)));
  hwoend=hwo(wo); ber=hwoend+1:hwoend+3;
  P(1:hwoend,j)=PA2(1:hwoend);
  P(ber,j)=(f2.*PA2(ber)+f1.*PC(ber,j))./(f1+f2); 
  else
  P(:,j)=P532(:,j); %nehme Profil aus rawtodat
  end
end
P=P532;

[PCerr,SNR532C] = estinoiseprofile(PC,H,phsoll);
%PCerr=real(sqrt(PC))./50+1e-4; 
wo=find(PC<Pschwelle); PC(wo)=Pschwelle;
[PAerr,SNR532A] = estinoiseprofile(PA,H,phsoll);
%PAerr=real(sqrt(PA))./50+1e-4; 
wo=find(PA<Pschwelle); PA(wo)=Pschwelle;
[Perr,SNR532K] = estinoiseprofile(P,H,phsoll);
%Perr=real(sqrt(P))./50+1e-4; 
wo=find(P<Pschwelle); P(wo)=Pschwelle;
P532final = P;   P532finalerr=Perr;
wofit = find( H >FitRangeA(1) & H< FitRangeA(2));

los =1; bis =size(PC,2); %%%length(Sel532C)

Btemp=zeros(Hsize,1); Btemperr=Btemp; 
Betaaer=zeros(Hsize,1); Betaaererr=Btemp; 
BetaAer532C=zeros(size(PC)); BetaAer532Cerr=zeros(size(PC)); 
dBeta532CdLR=zeros(size(PC));
dBeta532CdR=zeros(size(PC)); 
dBeta532CdP=zeros(size(PC));
BSR532CKlett=zeros(size(PC)); BSR532CKletterr=zeros(size(PC));
attenu532P=zeros(size(PC));
attenu532PAE=zeros(size(PC));

BetaAer532A=zeros(size(PA)); BetaAer532Aerr=zeros(size(PA)); 
dBeta532AdLR=zeros(size(PA));
dBeta532AdR=zeros(size(PA)); 
dBeta532AdP=zeros(size(PA));
BSR532AKlett=zeros(size(PA)); BSR532AKletterr=zeros(size(PA)); 
BSRFit532Avektor=zeros(Hsize,1);


% Schleife ueber alle Zeitschritte: erst counting, dann analog
for i=los:bis,  
    
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
PC(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

aoint=exp(-1.*qdrupvar(H(Sel), AlRay532(Sel,i)));
arint=exp(-1.*qdrupvar(H(Sel), AlRay607(Sel,i)));
BSR532quot(Sel,i)=BSR532quot(Sel,i).*arint./aoint;
BSR532urquot(:,i)=BSR532quot(:,i);

% Rechnung counting

[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit532, FitRangeC, H(Sel), PC(Sel,i), PCerr(Sel,i), ...
    LR532arr(Sel,i), AlRay532(Sel,i), BeRa532(Sel,i));
% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2( BSRAtFit532, FitRangeC, DN532(i),H(Sel), HStat, PC(Sel,i), PCerr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR532arr(Sel,i), LR532arrerr(Sel,i), Wvl532, EmitPol, RecPol, ...
% Temp(Sel,i));
Betaaer(Sel)=Beta-BeRa532(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit532err)+abs(dBdLR.*LR532arrerr(Sel,i))+abs(dBdP.*PCerr(Sel,i));
Btemp(Sel)=Beta./BeRa532(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532(Sel,i); 

BSR532CKlett(Sel,i)=Btemp(Sel); BSR532CKletterr(Sel,i)=Btemperr(Sel); 
BetaAer532C(Sel,i)=Betaaer(Sel); BetaAer532Cerr(Sel,i)=Betaaererr(Sel,i); 

dBeta532CdP(Sel,i)=dBdP;
dBeta532CdR(Sel,i)=dBdR;
dBeta532CdLR(Sel,i)=dBdLR;

% Rechnung analog

Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
PA(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end


BSRAtFit532A=mymean(BSR532CKlett(wofit,i));   
if isnan(BSRAtFit532A), BSRAtFit532A=BSRAtFit532; end
BSRFit532Avektor(i)=BSRAtFit532A;

[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit532A, FitRangeA, H(Sel), PA(Sel,i), PAerr(Sel,i), ...
    LR532arr(Sel,i), AlRay532(Sel,i), BeRa532(Sel,i));
% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2( BSRAtFit532A, FitRangeA, DN532(i),H(Sel), HStat, PA(Sel,i), PAerr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR532arr(Sel,i), LR532arrerr(Sel,i), Wvl532, 'p', 'p', ...
% Temp(Sel,i));
if ~isempty(Sel)
Betaaer(Sel)=Beta-BeRa532(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit532err)+abs(dBdLR.*LR532arrerr(Sel,i))+abs(dBdP.*PAerr(Sel,i));
Btemp(Sel)=Beta./BeRa532(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532(Sel,i); 

BSR532AKlett(Sel,i)=Btemp(Sel); BSR532AKletterr(Sel,i)=Btemperr(Sel); 
BetaAer532A(Sel,i)=Betaaer(Sel); BetaAer532Aerr(Sel,i)=Betaaererr(Sel,i); 

dBeta532AdP(Sel,i)=dBdP;
dBeta532AdR(Sel,i)=dBdR;
dBeta532AdLR(Sel,i)=dBdLR;
end % if - anderenfalls bleiben BSR532A ... und alle =0
end

% Anpassung von analog und counting

BSR532KlebKlett = BSR532CKlett;  BSR532KlebKletterr = BSR532CKletterr;
hwo=find(H>8000 & H <12000);
nn=length(hwo);
f1=((1:nn)./nn)';
f2=((nn:-1:1)./nn)';
BSR532KlebKlett(1:hwo(1)-1,:)=BSR532AKlett(1:hwo(1)-1,:);
BSR532KlebKletterr(1:hwo(1)-1,:)=BSR532AKletterr(1:hwo(1)-1,:);
for j=los:bis,
  BSR532KlebKlett(hwo,j)=(f2.*BSR532AKlett(hwo,j)+f1.*BSR532CKlett(hwo,j))./(f1+f2);
  BSR532KlebKletterr(hwo,j)=(f2.*BSR532AKletterr(hwo,j)+f1.*BSR532CKletterr(hwo,j))./(f1+f2);
end


% Rechnung mit Kombinierten Profilen - P
BetaAer532Klett=zeros(size(P)); BetaAer532Kletterr=zeros(size(P)); 
dBeta532dLR=zeros(size(P));
dBeta532dR=zeros(size(P)); 
dBeta532dP=zeros(size(P));
BSR532Klett=zeros(size(P)); BSR532Kletterr=zeros(size(P)); 
C532Lidar=zeros(size(P));
%BSRFit532Avektor=zeros(Hsize,1);
Perr= real(sqrt(P))./50+1e-4;
for i= los:bis

Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end


% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2(BSRAtFit532arr(i), FitRangeC, DN532(i),H(Sel), HStat, P(Sel,i), Perr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR532arr(Sel,i), LR532arrerr(Sel,i), Wvl532, 'p', 'p', ...
% Temp(Sel,i));
[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit532arr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR532arr(Sel,i), AlRay532(Sel,i), BeRa532(Sel,i));
Betaaer(Sel)=Beta-BeRa532(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit532err)+abs(dBdLR.*LR532arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa532(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532(Sel,i); 

BSR532Klett(Sel,i)=Btemp(Sel); BSR532Kletterr(Sel,i)=Btemperr(Sel); 
BetaAer532Klett(Sel,i)=Betaaer(Sel); BetaAer532Kletterr(Sel,i)=Betaaererr(Sel,i); 
attenu532P(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);
attenu532PAE(Sel,i) = abs(Perr(Sel,i).*H(Sel).^2 ./ CLidar(Sel(1))) + abs(attenu532P(Sel,i) ./ CLidar(Sel(1)).*0.1.*CLidar(Sel(1)));

dBeta532dP(Sel,i)=dBdP;
dBeta532dR(Sel,i)=dBdR;
dBeta532dLR(Sel,i)=dBdLR;
C532Lidar(Sel,i) =CLidar;

% Korrektur der Aerosolextinktion muss aus analog erfolgen wegen Saettigung
% von Counting
AlphaAer532=BetaAer532Klett(:,i).*LR532arr(:,i);  aaerint=exp(-qdrupvar(H(Sel),Wvlfaktor2.*AlphaAer532(Sel)));
BSR532quot(Sel,i)=BSR532urquot(Sel,i)./aaerint;
end


%
% die entscheidende Iteration
%
if LRiter ~=0,

% iterativ
% wohint=40:60;
% wowolke=78:170;
z=zeros(1,dimen(2)); 
hwo=find(H>1500 & H< FitRangeC(1));
for j=1:dimen(2), z(j)=mymean(BSR532Klett(hwo,j)); end
listvalidmeasure=find(DATHV532p>650 & z > 1.05);


FitSel=find(H>FitRangeC(1) & H<FitRangeC(2));
gut = 1; % Prozent-Abweichung
q=zeros(bis,1).*NaN; q2=zeros(bis,1); liste=zeros(bis,1)+FitSel(1);
for j=1:bis, 
    q2(j)=max(BSR532Klett(1:FitSel(1),j));
    wo = find(BSR532Klett(1:FitSel(1),j) > BSR532schwelle);
    if ~isempty(wo), liste(j)=wo(1)-1; end 
end
[BSR00, wo] = min(q2(listvalidmeasure));
wohint = (Hmin1:min(liste)-1);
if isempty(wohint), 
   wo2=find(liste>Hmin1);
   liste = liste(wo2);
%    [m,w]=min(liste);
%    liste(w)=NaN;
end
wohint = (Hmin1:min(liste)-1);
for j=1:bis,
   q(j)=mymean(BSR532Klett(wohint,j));
end

BSR00 = q(wo);
HX=0:300:30000; 

Btemp2 = Btemp;

% LR532arr=ones(size(P)).*val_LR532;
% wowolke = find(BSR532Klett > BSR532schwelle);
% LR532arr(wowolke) = val_LR532Wo;
% while BSR00 unglaubw?rdig
iterp=0; condip=1; itpmax=5;
while condip
iterp=iterp+1;
LR532arr(:) = val_LR532;
BSR00problem=0;

for i= los:bis
    
wowolke = find(BSR532Klett(1:FitSel(1),i) > BSR532schwelle);
LR532arr(wowolke,i) = val_LR532Wo;
if ~isempty(wowolke)

Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

if length(Sel) >10,
    
iter=0; condi=1;  hz=0; tz=0;  % hz / tz Zaehler fuer hohe - tiefe Werte
refer=zeros(itmax,1);  refer2=zeros(itmax,1);
while condi
iter = iter+1;    
% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2( BSRAtFit532arr(i), FitRangeC, DN532(i),H(Sel), HStat, P(Sel,i), Perr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR532arr(Sel,i), LR532arrerr(Sel,i), Wvl532, 'p','p', ...
% Temp(Sel,i));
[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit532arr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR532arr(Sel,i), AlRay532(Sel,i), BeRa532(Sel,i));
Betaaer(Sel)=Beta-BeRa532(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit532err)+abs(dBdLR.*LR532arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa532(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532(Sel,i); 


% [Beta2,dB2dR, dB2dLR, dB2dP] = ...
% klettinv_ableit2( BSRAtFit532arr(i), FitRangeC, DN532(i),H(Sel), HStat, P(Sel,i), Perr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR532arr(Sel,i)+1, LR532arrerr(Sel,i), Wvl532, 'p','p', ...
% Temp(Sel,i));
[Beta2, dB2dR, dB2dLR, dB2dP, CLidar2] = klettinv_ableit4( BSRAtFit532arr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR532arr(Sel,i)+1, AlRay532(Sel,i), BeRa532(Sel,i));
Btemp2(Sel)=Beta2./BeRa532(Sel,i); % Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532(Sel,i); 

refer(iter) = mymean(Btemp(wohint));
refer2(iter) = mymean(Btemp2(wohint));
dxdlr = refer2(iter)-refer(iter);
diffi = BSR00 - refer(iter);
deltax=diffi./dxdlr;


if refer(iter) ./ BSR00 > 1+gut ./ 100,          % BSR zu gro?
    LR532arr(wowolke,i) = LR532arr(wowolke,i) +deltax;
    if mymean(LR532arr(wowolke,i)) > LRobergr, LR532arr(wowolke,i) = LRobergr; hz=hz+1; end
    if mymean(LR532arr(wowolke,i)) < LRuntergr, LR532arr(wowolke,i) = LRuntergr; tz=tz+1; end
elseif refer(iter) ./ BSR00 < 1-gut ./ 100,      % BSR zu klein
    LR532arr(wowolke,i) = LR532arr(wowolke,i) +deltax;
    if mymean(LR532arr(wowolke,i)) > LRobergr, LR532arr(wowolke,i) = LRobergr; hz=hz+1; end
    if mymean(LR532arr(wowolke,i)) < LRuntergr, LR532arr(wowolke,i) = LRuntergr; tz=tz+1; end
else
    condi=0;
end

if iter >=itmax, condi=0; end     % disp('keine Konvergenz'); i, end
if hz>3, condi = 0; end      % disp('hz > 3'); i, end
if tz>3, condi = 0; end      % disp('tz > 3'); i, end

end % while

% dieses Ergebnis kann falsch sein da BSR00 falsch
BX = interp1(H,Btemp,HX,'spline');
if min(BX(2:60)) < BSRAtFit532arr(i), BSR00problem=BSR00problem+1; end


% fertig und abspeichern
BSR532Klett(Sel,i)=Btemp(Sel); BSR532Kletterr(Sel,i)=Btemperr(Sel); 
BetaAer532Klett(Sel,i)=Betaaer(Sel); BetaAer532Kletterr(Sel,i)=Betaaererr(Sel,i); 
attenu532P(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);

dBeta532dP(Sel,i)=dBdP;
dBeta532dR(Sel,i)=dBdR;
dBeta532dLR(Sel,i)=dBdLR;
C532Lidar(Sel,i) = CLidar;

end % if Daten vorhanden

else
    % wowolke,i,
end % wowolke nicht leer

end % for Zeitschritte

if BSR00problem==1, BSR00=BSR00+0.05;
elseif BSR00problem >=2, BSR00=BSR00+0.1;
elseif BSR00problem ==0, condip=0;
end
if iterp > itpmax, condip=0; end

end %while BSR00

end % LRiter


% der Low Bereich aus dem Quotienten der analog Signale
BSR532P=BSR532Klett;   BSR532Perr=BSR532Kletterr;     %nicht KlebKlett
BSR532PKlett=BSR532KlebKlett;   BSR532PKletterr=BSR532KlebKletterr;
hrc=find(H >2600 & H < 4200); % Heightrangecompare
hnp=find(H>2000); hnp=hnp(1); %Height Null Position
hmp=find(H>3500); hnp=hnp(1); %Height End Position
for jj=1:bis, 
    tmp=mymean(BSR532P(hrc,jj));
    BSR532quot(:,jj)=BSR532quot(:,jj)./mymean(BSR532quot(hrc,jj)).*tmp;
    [m,w]=min(abs(BSR532P(hnp:hmp,jj)-BSR532quot(hnp:hmp,jj)));
    BSR532P(1:hnp-1+w,jj)=BSR532quot(1:hnp-1+w,jj); 
end
% wenn 607 schwach
% BSR532P=BSR532Klett; 
BSR532P=BSR532Klett; % empirisch
BSR532Perr=BSR532KlebKletterr+BSRAtFit532err;


% Gesamtgroessen parallel
[BetaAer532P, BetaAer532Perr]=  bsr2bsc(BSR532P, BSR532Perr, Density, Wvl532, Temp, 'p', 'p');  
[BetaAer532Plang, DN532lang] = maketimegaps(BetaAer532P, DN532, 1/24/2);  
[BSR532Plang, DN532lang] = maketimegaps(BSR532P, DN532, 1/24/2);  
[LR532arrl, DN532lang] = maketimegaps(LR532arr, DN532, 1/24/2);  
[BetaAer532Pquot, BetaAer532Pquoterr]=  bsr2bsc(BSR532quot, BSR532Perr, Density, Wvl532, Temp, 'p', 'p');  

% optische Dicke
q=size(BetaAer532Plang); optdivisp=zeros(q(2),1); bb1=BetaAer532Plang(Sel,:); bx=find(bb1<0); bb1(bx)=0;
bx=find(isnan(bb1)); bb1(bx)=0;
for jj=1:q(2), optdivisp(jj)=sum(bb1(:,jj).*LR532arrl(Sel,jj).*dH)+LR532arrl(1,jj).*(H(1)-HStat).*bb1(1,jj); end

AlphaAer532P = BetaAer532P;
Tau532P=zeros(size(BetaAer532P));
Hdiff=H-HStat; Hdiff(2:end)=diff(H);

for j=1:bis, AlphaAer532P(:,j)=BetaAer532P(:,j).*LR532arr(:,i); 
Tau532P(4:end,j)=cumsum(AlphaAer532P(4:end,j).*Hdiff(4:end));
end
%Alpha532Ptot = AlphaAer532P+AlRay532;


%
%%% 532nm senkrecht polarisiert
%

LR532Sarr = LR532arr; % extin aus "P"
EmitPol='p';RecPol='s'; 
BeRa532S=Density.*raybckwq (Wvl532,'p','s', Temp, Density);  

BSR532Squot=P532S; 
for j=1: size(P532S,2),
    [m,w]=min(abs(DATDateNumStart(Sel607C) - DM(j)));
    BSR532Squot(:,j)=P532S(:,j)./P607final(:,w);
end
BSR532Surquot=BSR532Squot;

% Signaldefintion
PC=P532SC;
PA=P532SA;
P = PC;
%ber=151:180;
nn=3;   %length(ber);
f1=((1:nn)./nn)';
f2=((nn:-1:1)./nn)';
for j=1: length(Sel532A)
  hwo=150:300; %find(PC(:,j) <10 & PC(:,j)>0.1);
  hwo=find(PC(:,j) <4 & PC(:,j)>0.08 & H > 1500);
  [c,k]=minanpass(PA(:,j),PC(:,j),hwo,H,5);
  if k > 0,
  PA(:,j)=(PA(:,j)-c)./k;
  if length(hwo) > 40, hwo=hwo(1:40); end 
  [x,wo]=min(abs(PA(hwo,j)-PC(hwo,j)));
  hwoend=hwo(wo); ber=hwoend+1:hwoend+3;
  P(1:hwoend,j)=PA(1:hwoend,j);
  P(ber,j)=(f2.*PA(ber,j)+f1.*PC(ber,j))./(f1+f2);
  else
  P(:,j)=P532S(:,j);
  end
end
P = P532S;

[PCerr,SNR532SC] = estinoiseprofile(PC,H,phsoll);
%PCerr=real(sqrt(PC))./50+1e-4; 
wo=find(PC<Pschwelle); PC(wo)=Pschwelle;
[PAerr,SNR532SA] = estinoiseprofile(PA,H,phsoll);
%PAerr=real(sqrt(PA))./50+1e-4; 
wo=find(PA<Pschwelle); PA(wo)=Pschwelle;
[Perr,SNR532SK] = estinoiseprofile(P,H,phsoll);
%Perr=real(sqrt(P))./50+1e-4; 
wo=find(P<Pschwelle); P(wo)=Pschwelle;

P532Sfinal = P;  P532Sfinalerr=Perr;
wofit = find( H >FitRangeA(1) & H< FitRangeA(2));

los =1; bis =size(PC,2); %%%length(Sel532C)

Btemp=zeros(Hsize,1); Btemperr=Btemp; 
Betaaer=zeros(Hsize,1); Betaaererr=Btemp; 
BetaAer532SC=zeros(size(PC)); BetaAer532SCerr=zeros(size(PC)); 
dBeta532SCdLR=zeros(size(PC));
dBeta532SCdR=zeros(size(PC)); 
dBeta532SCdP=zeros(size(PC));
BSR532SCKlett=zeros(size(PC)); BSR532SCKletterr=zeros(size(PC));
attenu532S = zeros(size(PC));
AlphaAer532S = zeros(size(PC));


BetaAer532SA=zeros(size(PA)); BetaAer532SAerr=zeros(size(PA)); 
dBeta532SAdLR=zeros(size(PA));
dBeta532SAdR=zeros(size(PA)); 
dBeta532SAdP=zeros(size(PA));
BSR532SAKlett=zeros(size(PA)); BSR532SAKletterr=zeros(size(PA)); 
BSRFit532SAvektor=zeros(Hsize,1);


% Schleife ueber alle Zeitschritte: erst counting, dann analog
for i=los:bis,  
    
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
PC(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

aoint=exp(-1.*qdrupvar(H(Sel), AlRay532(Sel,i)));
arint=exp(-1.*qdrupvar(H(Sel), AlRay607(Sel,i)));
BSR532Squot(Sel,i)=BSR532Squot(Sel,i).*arint./aoint; BSR532Surquot(Sel,i) = BSR532Squot(Sel,i);


[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit532, FitRangeC, H(Sel), PC(Sel,i), PCerr(Sel,i), ...
    LR532Sarr(Sel,i), AlphaAer532P(Sel,i), AlRay532(Sel,i), BeRa532S(Sel,i));
Betaaer(Sel)=Beta-BeRa532S(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit532err)+abs(dBdLR.*LR532arrerr(Sel,i))+abs(dBdP.*PCerr(Sel,i));
Btemp(Sel)=Beta./BeRa532S(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532S(Sel,i); 

BSR532SCKlett(Sel,i)=Btemp(Sel); BSR532SCKletterr(Sel,i)=Btemperr(Sel); 
BetaAer532SC(Sel,i)=Betaaer(Sel); BetaAer532SCerr(Sel,i)=Betaaererr(Sel,i); 

dBeta532SCdP(Sel,i)=dBdP;
dBeta532SCdR(Sel,i)=dBdR;
dBeta532SCdLR(Sel,i)=dBdLR;

% Rechnung analog
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
PA(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

BSRAtFit532SA=mymean(BSR532SCKlett(wofit,i));   
if isnan(BSRAtFit532SA), BSRAtFit532SA = BSRAtFit532; end
BSRFit532SAvektor(i)=BSRAtFit532SA;

[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit532SA, FitRangeA, H(Sel), PA(Sel,i), PAerr(Sel,i), ...
    LR532Sarr(Sel,i), AlphaAer532P(Sel,i), AlRay532(Sel,i), BeRa532S(Sel,i));
if ~isempty(Sel),
Betaaer(Sel)=Beta-BeRa532S(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit532err)+abs(dBdLR.*LR532arrerr(Sel,i))+abs(dBdP.*PAerr(Sel,i));
Btemp(Sel)=Beta./BeRa532S(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532S(Sel,i); 

BSR532SAKlett(Sel,i)=Btemp(Sel); BSR532SAKletterr(Sel,i)=Btemperr(Sel); 
BetaAer532SA(Sel,i)=Betaaer(Sel); BetaAer532SAerr(Sel,i)=Betaaererr(Sel,i); 

dBeta532SAdP(Sel,i)=dBdP;
dBeta532SAdR(Sel,i)=dBdR;
dBeta532SAdLR(Sel,i)=dBdLR;
end

end

% Anpassung von analog und counting

BSR532SKlebKlett = BSR532SCKlett;  BSR532SKlebKletterr = BSR532SCKletterr;
hwo=find(H>8000 & H <12000);
nn=length(hwo);
f1=((1:nn)./nn)';
f2=((nn:-1:1)./nn)';
BSR532SKlebKlett(1:hwo(1)-1,:)=BSR532SAKlett(1:hwo(1)-1,:);
BSR532SKlebKletterr(1:hwo(1)-1,:)=BSR532SAKletterr(1:hwo(1)-1,:);
for j=los:bis,
  BSR532SKlebKlett(hwo,j)=(f2.*BSR532SAKlett(hwo,j)+f1.*BSR532SCKlett(hwo,j))./(f1+f2);
  BSR532SKlebKletterr(hwo,j)=(f2.*BSR532SAKletterr(hwo,j)+f1.*BSR532SCKletterr(hwo,j))./(f1+f2);
end


% Rechnung mit Kombinierten Profilen - P
BetaAer532SKlett=zeros(size(P)); BetaAer532SKletterr=zeros(size(P)); 
dBeta532SdLR=zeros(size(P));
dBeta532SdR=zeros(size(P)); 
dBeta532SdP=zeros(size(P));
BSR532SKlett=zeros(size(P)); BSR532SKletterr=zeros(size(P)); 
C532SLidar=zeros(size(P));
aaerint = 1;
Perr= real(sqrt(P))./50+1e-4;

for i= los:bis

Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

if length(Sel) >10,

% [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit532Sarr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
%     LR532Sarr(Sel,i), AlRay532(Sel,i), BeRa532S(Sel,i));
[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit532Sarr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR532Sarr(Sel,i), AlphaAer532P(Sel,i), AlRay532(Sel,i), BeRa532S(Sel,i));
Betaaer(Sel)=Beta-BeRa532S(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit532err)+abs(dBdLR.*LR532arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa532S(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532S(Sel,i); 

BSR532SKlett(Sel,i)=Btemp(Sel); BSR532SKletterr(Sel,i)=Btemperr(Sel); 
BetaAer532SKlett(Sel,i)=Betaaer(Sel); BetaAer532SKletterr(Sel,i)=Betaaererr(Sel,i); 
attenu532S(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);

dBeta532SdP(Sel,i)=dBdP;
dBeta532SdR(Sel,i)=dBdR;
dBeta532SdLR(Sel,i)=dBdLR;
C532SLidar(Sel,i) = CLidar;

end % length Sel

end % for j von bis

%
% die entscheidende Iteration
%

if LRiterS ~=0,

% iterativ
%wohint=40:60;
%wowolke=78:170;
z=zeros(1,dimen(2)); 
hwo=find(H>1500 & H< FitRangeC(1));
for j=1:dimen(2), z(j)=mymean(BSR532SKlett(hwo,j)); end
listvalidmeasure=find(DATHV532s>650 & z > 1.05);

FitSel=find(H>FitRangeC(1) & H<FitRangeC(2));
gut = 1; % Prozent-Abweichung
q=zeros(bis,1).*NaN; q2=zeros(bis,1); liste=zeros(bis,1)+FitSel(1);
for j=1:bis, 
    q2(j)=max(BSR532SKlett(1:FitSel(1),j));
    wo = find(BSR532SKlett(1:FitSel(1),j) > BSR532Sschwelle);
    if ~isempty(wo), liste(j)=wo(1)-1; end 
end
[BSR00, wo] = min(q2(listvalidmeasure));
wohint = (Hmin1:min(liste)-1);
if isempty(wohint), 
   wo2=find(liste>Hmin1);
   liste = liste(wo2);
%    [m,w]=min(liste);
%    liste(w)=NaN;
end
wohint = (Hmin1:min(liste));
for j=1:bis,
   q(j)=mymean(BSR532SKlett(wohint,j));
end

BSR00 = max([q(wo),1.1.*BSRAtFit532Sarr(wo)]);


Btemp2 = Btemp;

% LR532Sarr=ones(size(P)).*val_LR532S;
% wowolke = find(BSR532S > BSR532Sschwelle);
% LR532Sarr(wowolke) = val_LR532SWo;
iterp=0; condip=1; itpmax=2;
while condip
iterp=iterp+1;
LR532Sarr(:) = val_LR532S;
BSR00problem = 0;

for i= los:bis
% if i==70,
%     q=3;
% end
wowolke = find(BSR532SKlett(1:FitSel(1),i) > BSR532Sschwelle);
LR532Sarr(wowolke,i) = val_LR532SWo;
if ~isempty(wowolke)
    
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

if length(Sel) >10,
    
iter=0; condi=1;  hz=0; tz=0;  % hz / tz Zaehler fuer hohe - tiefe Werte
refer=zeros(itmax,1);  refer2=zeros(itmax,1);
while condi
iter = iter+1;    

[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit532Sarr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR532Sarr(Sel,i), AlphaAer532P(Sel,i), AlRay532(Sel,i), BeRa532S(Sel,i));
Betaaer(Sel)=Beta-BeRa532S(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit532err)+abs(dBdLR.*LR532arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa532S(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532S(Sel,i); 


[Beta2, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit532Sarr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    (LR532Sarr(Sel,i)+1), AlphaAer532P(Sel,i), AlRay532(Sel,i), BeRa532S(Sel,i));

Btemp2(Sel)=Beta2./BeRa532S(Sel,i); % Btemperr(Sel)=Betaaererr(Sel,i)./BeRa532S(Sel,i); 

refer(iter) = mymean(Btemp(wohint));
refer2(iter) = mymean(Btemp2(wohint));
dxdlr = refer2(iter)-refer(iter);
diffi = BSR00 - refer(iter);
deltax=diffi./dxdlr;


if refer(iter) ./ BSR00 > 1+gut ./ 100,          % BSR zu gro?
    LR532Sarr(wowolke,i) = LR532Sarr(wowolke,i) +deltax;
    if mymean(LR532Sarr(wowolke,i)) > LRSobergr, LR532Sarr(wowolke,i) = LRSobergr; hz=hz+1; end
    if mymean(LR532Sarr(wowolke,i)) < LRSuntergr, LR532Sarr(wowolke,i) = LRSuntergr; tz=tz+1; end
elseif refer(iter) ./ BSR00 < 1-gut ./ 100,      % BSR zu klein
    LR532Sarr(wowolke,i) = LR532Sarr(wowolke,i) +deltax;
     if mymean(LR532Sarr(wowolke,i)) > LRSobergr, LR532Sarr(wowolke,i) = LRSobergr; hz=hz+1; end
    if mymean(LR532Sarr(wowolke,i)) < LRSuntergr, LR532Sarr(wowolke,i) = LRSuntergr; tz=tz+1; end
else
    condi=0;
end

if iter >=itmax, condi=0; end     % disp('keine Konvergenz'); i, end
if hz>3, condi = 0; end      % disp('hz > 3'); i, end
if tz>3, condi = 0; end      % disp('tz > 3'); i, end

end % while

% dieses Ergebnis kann falsch sein da BSR00 falsch
BX = interp1(H,Btemp,HX,'spline');
if min(BX(2:60)) < BSRAtFit532Sarr(i), BSR00problem=BSR00problem+1; end
    

BSR532SKlett(Sel,i)=Btemp(Sel); BSR532SKletterr(Sel,i)=Btemperr(Sel); 
BetaAer532SKlett(Sel,i)=Betaaer(Sel); BetaAer532SKletterr(Sel,i)=Betaaererr(Sel,i); 
attenu532S(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);


dBeta532SdP(Sel,i)=dBdP;
dBeta532SdR(Sel,i)=dBdR;
dBeta532SdLR(Sel,i)=dBdLR;
C532SLidar(Sel,i) = CLidar;

end % if Daten vorhanden

end % wowolke nicht leer


AlphaAer532S(:,i)=BetaAer532SKlett(:,i).*LR532Sarr(:,i);  
aaerint=exp(-qdrupvar(H(Sel),AlphaAer532S(Sel,i)- (1).*AlphaAer532P(Sel,i)));
BSR532Squot(Sel,i)=BSR532Surquot(Sel,i)./aaerint;

end % j Zeitschritte

if BSR00problem==1, BSR00=BSR00+0.05;
elseif BSR00problem >=2, BSR00=BSR00+0.1;
elseif BSR00problem ==0, condip=0;
end
if iterp > itpmax, condip=0; end

end %while BSR00

%Q= BSR532SKlett(1:FitSel(1),:); wo=find(Q < 1.002); Q(wo)= 1.002; BSR532SKlett(1:FitSel(1),:)=Q;

end  %%LRiter


% der Low Bereich aus dem Quotienten der analog Signale
BSR532S=BSR532SKlett;   BSR532Serr=BSR532SKletterr;
hrc=find(H >2600 & H < 4200); % Heightrangecompare
hnp=find(H>2000); hnp=hnp(1); %Height Null Position
hmp=find(H>3500); hnp=hnp(1); %Height End Position
for jj=1:bis, 
    tmp=mymean(BSR532S(hrc,jj));
    BSR532Squot(:,jj)=BSR532Squot(:,jj)./mymean(BSR532Squot(hrc,jj)).*tmp;
    [m,w]=min(abs(BSR532S(hnp:hmp,jj)-BSR532Squot(hnp:hmp,jj)));
    BSR532S(1:hnp-1+w,jj)=BSR532Squot(1:hnp-1+w,jj); 
end
% wenn  607 schwach:
% BSR532S=BSR532SKlett; 
BSR532S=BSR532SKlett; 
BSR532Serr=BSR532SKletterr+BSRAtFit532err; 


% Defintion Gesamtgroessen
[BetaAer532S, BetaAer532Serr]=  bsr2bsc(BSR532S, BSR532Serr, Density, Wvl532, Temp, 's', EmitPol );  
[BetaAer532Slang, DN532Slang] = maketimegaps(BetaAer532S, DN532S, 1/24/2);  
[BSR532Slang, DN532Slang] = maketimegaps(BSR532S, DN532S, 1/24/2);  
[BetaAer532Squot, BetaAer532Squoterr]=  bsr2bsc(BSR532Squot, BSR532Serr, Density, Wvl532, Temp, 's', EmitPol );

% optische Dicke
q=size(BetaAer532S); optdivisS=zeros(q(2),1); bb1=BetaAer532S(Sel,:); bx=find(bb1<0); bb1(bx)=0;
bx=find(isnan(bb1)); bb1(bx)=0;
for jj=1:q(2), optdivisS(jj)=sum(bb1(:,jj).*LR532Sarr(Sel,jj).*dH)+LR532Sarr(1,jj).*(H(1)-HStat).*bb1(1,jj); end

AlphaAer532S = BetaAer532S;
Tau532S=zeros(size(BetaAer532S));


for j=1:bis, AlphaAer532S(:,j)=AlphaAer532S(:,j).*LR532Sarr(:,i); 
Tau532S(:,j)=cumsum(AlphaAer532S(:,j).*Hdiff);
aaerint=exp(-qdrupvar(H(Sel),AlphaAer532S(Sel,i)- (1-Wvlfaktor).*AlphaAer532P(Sel,i)));
BSR532Squot(Sel,i)=BSR532Surquot(Sel,i)./aaerint;
end

%%% Zusammenfassung von 532 parallel und 532 senkrecht

% 532 total zusammenfassen p und s koennen unterschiedlich lang sein
BetaAer532tot = zeros(size(P532A));
BetaAer532toterr = zeros(size(P532A));
BetaRay532tot = zeros(size(P532A));
BSR532tot = zeros(size(P532A));
BSR532toterr = zeros(size(P532A));
AeroDep532 = zeros(size(P532A));
AeroDep532err = zeros(size(P532A));
los =1; bis =length(Sel532A);
for j=los:bis,
k= find(DATDateNumStart(Sel532SA) == DATDateNumStart(Sel532A(j)));
if length(k) ==1,  
   BetaAer532tot(:,j) = BetaAer532S(:,k) + BetaAer532P(:,j);
   BetaAer532toterr(:,j) = BetaAer532Serr(:,k) + BetaAer532Perr(:,j);
   BetaRay532tot(:,j)   = BeRa532S(:,k) + BeRa532(:,j);
   BSR532tot(:,j)     = BetaAer532tot(:,j) ./ BetaRay532tot(:,j) +1;
   BSR532toterr(:,j)     = BetaAer532toterr(:,j) ./ BetaRay532tot(:,j);
   AeroDep532(:,j) = BetaAer532S(:,k) ./ BetaAer532P(:,j);
   AeroDep532err(:,j) = BetaAer532Serr(:,k) ./ BetaAer532P(:,j) + ...
                        AeroDep532(:,j).*BetaAer532Perr(:,j) ./ BetaAer532P(:,j);
end
end



LR532Raman = zeros(size(P607C));   LR532Ramanerr=zeros(size(P607C)); 
los =1; bis =length(Sel607C);
for j=los:bis,
k= find(DATDateNumStart(Sel532C) == DATDateNumStart(Sel607C(j)));
if length(k) ==1,  
   LR532Raman(:,j)= AlphaAer532Raman(:,j) ./ BetaAer532tot(:,k);
   LR532Ramanerr(:,j)= abs(AlphaAer532Ramanerr(:,j) ./BetaAer532tot(:,k))+ ...
       abs(AlphaAer532Raman(:,j) .*BetaAer532toterr(:,k) ./BetaAer532tot(:,k).^2);
end
end


% LR532schaetz = zeros(size(P607x));   LR355schaetzerr=zeros(size(P607x)); 
% los =1; bis =length(Sel607C);
% for j=los:bis,
% k= find(DATDateNumStart(Sel532C) == DATDateNumStart(Sel607C(j)));
% if length(k) ==1,  
%    BetaAer532totx(:,k) = interp1(H,BetaAer532tot(:,k), Hx607);
%    BetaAer532totxerr(:,k) = interp1(H,BetaAer532toterr(:,k), Hx607);
%    h1 = Alpha532schaetz(:,j) - Alpha532schaetzMin(:,j);
%    h2 = Alpha532schaetzMax(:,j) - Alpha532schaetz(:,j);
%    Alpha532schaetzerr =  max(h1,h2);
%    LR532schaetz(:,j)= Alpha532schaetz(:,j) ./ BetaAer532totx(:,k);
%    LR532schaetzerr(:,j)= abs(Alpha532schaetzerr ./BetaAer532totx(:,k))+ ...
%        abs(Alpha532schaetz(:,j) .*BetaAer532totxerr(:,k) ./BetaAer532totx(:,k).^2);
% end
% end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Rechnung fuer 355nm
Wvl355=Wvl532.*2./3; WvlR355= ramshift(Wvl355,'N2');
WvlfaktorUV= 1+(Wvl355./WvlR355).^Angexpminus;
WvlfaktorUV2= 1-(Wvl355./WvlR355).^Angexpminus;

LR355arrerr=ones(size(LR355arr))+15;   
 

EmitPol='p';RecPol='p'; 
BeRa355=Density.*raybckwq (Wvl355,'p','p', Temp, Density);  H2=H-HStat;
AlRay355 = Density .* raytotwq( Wvl355, Temp, Density);
AlRay387 = Density .* raytotwq( WvlR355, Temp, Density);



% Signaldefinition
PC=P355C;
PA=P355A;
P = P355;
% %ber=151:180;
% nn= 3; %length(ber);
% f1=((1:nn)./nn)';
% f2=((nn:-1:1)./nn)';
% for j=1: length(Sel355A)
%   %hwo=150:300; %find(PC(:,j) <10 & PC(:,j)>0.1);
%   hwo =find(PC(:,j) <10 & PC(:,j)>0.1 & PA(:,j) >0.1);
%   [c,k]=minanpass(PA(:,j),PC(:,j),hwo,H,5);
%   if k > 0,
%   PA(:,j)=(PA(:,j)-c)./k;
%   [x,wo]=min(abs(PA(hwo,j)-PC(hwo,j)));
%   hwoend=hwo(wo); ber=hwoend+1:hwoend+3;
%   P(1:hwoend,j)=PA(1:hwoend,j); 
%   P(ber,j)=(f2.*PA(ber,j)+f1.*PC(ber,j))./(f1+f2);
%   hram=find(P387C(:,j) <10 & P387C(:,j)>0.03);
%   [c,k]=minanpass(P387A(:,j),P387C(:,j),hram,H,5);
%   P387A(:,j)=(P387A(:,j)-c)./k;
%   else
%   P(:,j)=P355(:,j);
%   end
% end

[PCerr,SNR355C] = estinoiseprofile(PC,H,phsoll);
%PCerr=real(sqrt(PC))./50+1e-4; 
wo=find(PC<Pschwelle); PC(wo)=Pschwelle;
[PAerr,SNR355A] = estinoiseprofile(PA,H,phsoll);
%PAerr=real(sqrt(PA))./50+1e-4; 
wo=find(PA<Pschwelle); PA(wo)=Pschwelle;
[Perr,SNR355K] = estinoiseprofile(P,H,phsoll);
%Perr=real(sqrt(P))./50+1e-4; 
wo=find(P<Pschwelle); P(wo)=Pschwelle;

wo=find(PC<Pschwelle); PC(wo)=Pschwelle;
%PCerr=real(sqrt(PC))./50+1e-4; 
[PCerr,SNR355C] = estinoiseprofile(PC,H);
wo=find(PA<Pschwelle); PA(wo)=Pschwelle;
%PAerr=real(sqrt(PA))./50+1e-4; 
[PAerr,SNR355A] = estinoiseprofile(PA,H);
wo=find(P<Pschwelle); P(wo)=Pschwelle;
%Perr=real(sqrt(P))./50+1e-4; 
[Perr,SNR355K] = estinoiseprofile(P,H);
P355final = P; P355finalerr= Perr;
wofit = find( H >FitRangeA(1) & H< FitRangeA(2));

los =1; bis =size(PC,2); %%%length(Sel355C)


%%% Berechnung von Alpha532 mit der Raman-Methode

for jj=1:bis;
  ARayad=O3Density.*o3abswq(Wvl355,Temp(:,jj)); ARamad=O3Density.*o3abswq(WvlR355,Temp(:,jj)); 
  AlRay355(:,jj)=AlRay355(:,jj)+ARayad; AlRay387(:,jj)=AlRay387(:,jj)+ARamad;   
end 

Ptest = P387;
Mixx = [Wvl355, WvlR355]; 
ff=(Wvl355./WvlR355).^(-WvlFct)+1;
AlphaAer355Raman=zeros(size(Ptest));
AlphaAer355Ramanerr=zeros(size(Ptest));
Alm355x =1e-9.*ones(dimen(1),1);
AOD355=zeros(size(Ptest));
AODUV=zeros(size(Ptest));
iter = 6; % Rechenzeit, aber bringt Verbesserung
anfang=find(H>700 & H< 900);

for jj=1:dimen(2),   
TauRay355=qdrupvar(H,AlRay355(:,jj)+AlRay387(:,jj));    
Aterm2 = exp(qdrupvar(H,AlRay355(:,jj)+AlRay387(:,jj)));
XX= Density(:,jj)./H.^2 ./Aterm2;
q1=Ptest(:,jj).*H.^2./Density(:,jj);
q2=poissonmittel(q1,5);
q3=q2.*0.6 + 0.4.*q1;
PR(:,jj)=q3./H.^2.*Density(:,jj);
PRerr=real(sqrt(PR(:,jj)))./50+1e-4;
PRx(:,jj) = PR(:,jj); 

for j=1:iter,

[ Alm355, Alm355AE, dAdSgn, dAdRho] = alphaaer( H(HSel), HStat, ...        
   PR(HSel,jj), PRerr(HSel), Mixx, [], Density(HSel,jj), ...
   Temp(HSel, jj),Frqtreat, WvlFct, ...
   O3Density(HSel), [], 'ableit3');    %diffquotneu  ableit2
Alm355x(HSel) = Alm355;
wo=find(isnan(Alm355x)); Alm355x(wo)=1e-9;
%Alm355x(1:400) = (Alm355x(1:400)+1e-5)./2;
PR(:,jj) = XX.*exp(-qdrupvar(H,ff.*Alm355x));
q2=find((~isnan(PR(:,jj)))); %% NaNs am Anfang herausschneiden
PR(:,jj)=PR(:,jj)./mymean(PR(anfang,jj)).*mymean(PRx(anfang,jj)); 
%keyboard
end % iter
AlphaAer355Raman(HSel,jj)=Alm355;
AlphaAer355Ramanerr(HSel,jj)=Alm355AE;
AOD355(:,jj) = (real(log(Density(:,jj)./H.^2./PR(:,jj)))-TauRay355) ./ff;
AODUV(:,jj) = (real(log(Density(:,jj)./H.^2./PR(:,jj)))-TauRay355);
end % length(Sel355A),
P387final = PR; 





% PR=P387A;
% PRAE=real(sqrt(PR))./50*1e-4;
% s0 = 1600; st=60;
% [ Alpha355schaetz, Alpha355schaetzMin, Alpha355schaetzMax, P387x, Hx387] = ...
%  alphaschaetz( H, HStat, PR, PRAE, Mixx, Density, Temp, O3Density, WvlFct, s0, st);
%  
    





%%%355nm parallel polarisiert


% Rueckstreuverhaeltnis nach Quotientenmethode
BSR355quot=P355;
for j=1: bis,
    [m,w]=min(abs(DATDateNumStart(Sel387A) - DN355(j)));
    BSR355quot(:,j)=P355(:,j)./PR(:,w);
end
BSR355urquot = BSR355quot;


Btemp=zeros(Hsize,1); Btemperr=Btemp; 
Betaaer=zeros(Hsize,1); Betaaererr=Btemp; 
BetaAer355C=zeros(size(PC)); BetaAer355Cerr=zeros(size(PC)); 
dBeta355CdLR=zeros(size(PC));
dBeta355CdR=zeros(size(PC)); 
dBeta355CdP=zeros(size(PC));
BSR355CKlett=zeros(size(PC)); BSR355CKletterr=zeros(size(PC));

BetaAer355A=zeros(size(PA)); BetaAer355Aerr=zeros(size(PA)); 
dBeta355AdLR=zeros(size(PA));
dBeta355AdR=zeros(size(PA)); 
dBeta355AdP=zeros(size(PA));
BSR355AKlett=zeros(size(PA)); BSR355AKletterr=zeros(size(PA)); 
BSRFit355Avektor=zeros(Hsize,1);
attenu355P = zeros(size(PC));


% Schleife ueber alle Zeitschritte: erst counting, dann analog
for i=los:bis,  
    
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
PC(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

aoint=exp(-1.*qdrupvar(H(Sel), AlRay355(Sel,i)));
arint=exp(-1.*qdrupvar(H(Sel), AlRay387(Sel,i)));
BSR355quot(Sel,i)=BSR355quot(Sel,i).*arint./aoint;
BSR355urquot(Sel,i)=BSR355quot(Sel,i);

% Rechnung counting

[Beta, dBdR, dBdLR, dBdP, ~] = klettinv_ableit4( BSRAtFit355, FitRangeC, H(Sel), PC(Sel,i), PCerr(Sel,i), ...
    LR355arr(Sel,i), AlRay355(Sel,i), BeRa355(Sel,i));
% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2(BSRAtFit355, FitRangeC, DN355(i),H(Sel), HStat, PC(Sel,i), PCerr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR355arr(Sel,i), LR355arrerr(Sel,i), Wvl355, EmitPol, RecPol, ...
% Temp(Sel,i));
Betaaer(Sel)=Beta-BeRa355(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit355err)+abs(dBdLR.*LR355arrerr(Sel,i))+abs(dBdP.*PCerr(Sel,i));
Btemp(Sel)=Beta./BeRa355(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa355(Sel,i); 

BSR355CKlett(Sel,i)=Btemp(Sel); BSR355CKletterr(Sel,i)=Btemperr(Sel); 
BetaAer355C(Sel,i)=Betaaer(Sel); BetaAer355Cerr(Sel,i)=Betaaererr(Sel,i); 

dBeta355CdP(Sel,i)=dBdP;
dBeta355CdR(Sel,i)=dBdR;
dBeta355CdLR(Sel,i)=dBdLR;

% Rechnung analog

Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
PA(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

BSRAtFit355A=mymean(BSR355CKlett(wofit,i));   
if isnan(BSRAtFit355A), BSRAtFit355A=1.015; end
BSRFit355Avektor(i)=BSRAtFit355A;

[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit355A, FitRangeA, H(Sel), PA(Sel,i), PAerr(Sel,i), ...
    LR355arr(Sel,i), AlRay355(Sel,i), BeRa355(Sel,i));
% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2( BSRAtFit355A, FitRangeA, DN355(i),H(Sel), HStat, PA(Sel,i), PAerr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR355arr(Sel,i), LR355arrerr(Sel,i), Wvl355, EmitPol, RecPol, ...
% Temp(Sel,i));
if ~isempty(Sel),
Betaaer(Sel)=Beta-BeRa355(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit355err)+abs(dBdLR.*LR355arrerr(Sel,i))+abs(dBdP.*PAerr(Sel,i));
Btemp(Sel)=Beta./BeRa355(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa355(Sel,i); 

BSR355AKlett(Sel,i)=Btemp(Sel); BSR355AKletterr(Sel,i)=Btemperr(Sel); 
BetaAer355A(Sel,i)=Betaaer(Sel); BetaAer355Aerr(Sel,i)=Betaaererr(Sel,i); 

dBeta355AdP(Sel,i)=dBdP;
dBeta355AdR(Sel,i)=dBdR;
dBeta355AdLR(Sel,i)=dBdLR;
end % if Daten da

end

% Anpassung von analog und counting
BSR355KlebKlett = BSR355CKlett;  BSR355KlebKletterr = BSR355CKletterr;
hwo=find(H>8000 & H <12000);
nn=length(hwo);
f1=((1:nn)./nn)';
f2=((nn:-1:1)./nn)';
BSR355KlebKlett(1:hwo(1)-1,:)=BSR355AKlett(1:hwo(1)-1,:);
BSR355KlebKletterr(1:hwo(1)-1,:)=BSR355AKletterr(1:hwo(1)-1,:);
for j=los:bis,
  BSR355KlebKlett(hwo,j)=(f2.*BSR355AKlett(hwo,j)+f1.*BSR355CKlett(hwo,j))./(f1+f2);
  BSR355KlebKletterr(hwo,j)=(f2.*BSR355AKletterr(hwo,j)+f1.*BSR355CKletterr(hwo,j))./(f1+f2);
end

% Rechnung mit Kombinierten Profilen - P
BetaAer355Klett=zeros(size(P)); BetaAer355Kletterr=zeros(size(P)); 
dBeta355dLR=zeros(size(P));
dBeta355dR=zeros(size(P)); 
dBeta355dP=zeros(size(P));
BSR355Klett=zeros(size(P)); BSR355Kletterr=zeros(size(P)); 
C355Lidar = zeros(size(P));
Perr= real(sqrt(P))./50+1e-4;
for i= los:bis

Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end


% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2( BSRAtFit355arr(i), FitRangeC, DN355(i),H(Sel), HStat, P(Sel,i), Perr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR355arr(Sel,i), LR355arrerr(Sel,i), Wvl355, 'p','p', ...
% Temp(Sel,i));
% 
[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit355arr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR355arr(Sel,i), AlRay355(Sel,i), BeRa355(Sel,i));

Betaaer(Sel)=Beta-BeRa355(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit355err)+abs(dBdLR.*LR355arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa355(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa355(Sel,i); 

BSR355Klett(Sel,i)=Btemp(Sel); BSR355Kletterr(Sel,i)=Btemperr(Sel); 
BetaAer355Klett(Sel,i)=Betaaer(Sel); BetaAer355Kletterr(Sel,i)=Betaaererr(Sel,i); 
attenu355P(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);


dBeta355dP(Sel,i)=dBdP;
dBeta355dR(Sel,i)=dBdR;
dBeta355dLR(Sel,i)=dBdLR;
C355Lidar(Sel,i) = CLidar;

AlphaAer355=BetaAer355Klett(:,i).*LR355arr(:,i);  aaerint=exp(-qdrupvar(H(Sel),WvlfaktorUV2.*AlphaAer355(Sel)));
BSR355quot(Sel,i)=BSR355urquot(Sel,i)./aaerint;
end

BSR355keineiter = BSR355Klett; 
BetaAer355keineiter = BetaAer355Klett;

%
% die entscheidende Iteration
%

if LRiter ~=0,

% iterativ
% wohint=40:60;
% wowolke=78:170;
z=zeros(1,dimen(2)); 
hwo=find(H>1500 & H< FitRangeC(1));
for j=1:dimen(2), z(j)=mymean(BSR355Klett(hwo,j)); end
listvalidmeasure=find(DATHV355p>650 & z > 1.015);


FitSel=find(H>FitRangeC(1) & H<FitRangeC(2));
gut = 1; % Prozent-Abweichung
q=zeros(bis,1).*NaN; q2=zeros(bis,1); liste=zeros(bis,1)+FitSel(1);
for j=1:bis, 
    q2(j)=max(BSR355Klett(1:FitSel(1),j));
    wo = find(BSR355Klett(1:FitSel(1),j) > BSR355schwelle);
    if ~isempty(wo), liste(j)=wo(1)-1; end 
end
[BSR00, wo] = min(q2(listvalidmeasure));
wohint = (Hmin1:min(liste)-1);
if isempty(wohint), 
   wo2=find(liste>Hmin1);
   liste = liste(wo2); 
%    [m,w]=min(liste);
%    liste(w)=NaN;
end
wohint = (Hmin1:min(liste));
for j=1:bis,
   q(j)=mymean(BSR355Klett(wohint,j));
end

BSR00 = q(wo);

Btemp2 = Btemp;

% wowolke = find(BSR355 > BSR355schwelle);
% LR355arr(wowolke) = val_LR355Wo;
iterp=0; condip=1; itpmax=5;
while condip
iterp=iterp+1;
LR355arr(:) = val_LR355;
BSR00problem=0;

for i= los:bis

wowolke = find(BSR355Klett(1:FitSel(1),i) > BSR355schwelle);
LR355arr(wowolke,i) = val_LR355Wo;
if ~isempty(wowolke)
    
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

if length(Sel) >10,
    
iter=0; condi=1;  hz=0; tz=0;  % hz / tz Zaehler fuer hohe - tiefe Werte
refer=zeros(itmax,1);  refer2=zeros(itmax,1);
while condi
iter = iter+1;    
[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit355arr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR355arr(Sel,i), AlRay355(Sel,i), BeRa355(Sel,i));
Betaaer(Sel)=Beta-BeRa355(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit355err)+abs(dBdLR.*LR355arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa355(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa355(Sel,i); 


[Beta2, dB2dR, dB2dLR, dB2dP, CLidar2] = klettinv_ableit4( BSRAtFit355arr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR355arr(Sel,i)+1, AlRay355(Sel,i), BeRa355(Sel,i));

Btemp2(Sel)=Beta2./BeRa355(Sel,i); % Btemperr(Sel)=Betaaererr(Sel,i)./BeRa355S(Sel,i); 

refer(iter) = mymean(Btemp(wohint));
refer2(iter) = mymean(Btemp2(wohint));
dxdlr = refer2(iter)-refer(iter);
diffi = BSR00 - refer(iter);
deltax=diffi./dxdlr;


if refer(iter) ./ BSR00 > 1+gut ./ 100,          % BSR zu gro?
    LR355arr(wowolke,i) = LR355arr(wowolke,i) +deltax;
    if mymean(LR355arr(wowolke,i)) > LRobergr, LR355arr(wowolke,i) = LRobergr; hz=hz+1; end
    if mymean(LR355arr(wowolke,i)) < LRuntergr, LR355arr(wowolke,i) = LRuntergr; tz=tz+1; end
elseif refer(iter) ./ BSR00 < 1-gut ./ 100,      % BSR zu klein
    LR355arr(wowolke,i) = LR355arr(wowolke,i) +deltax;
    if mymean(LR355arr(wowolke,i)) > LRobergr, LR355arr(wowolke,i) = LRobergr; hz=hz+1; end
    if mymean(LR355arr(wowolke,i)) < LRuntergr, LR355arr(wowolke,i) = LRuntergr; tz=tz+1; end
else
    condi=0;
end

if iter >=itmax, condi=0; end    % disp('keine Konvergenz'); i, end
if hz>3, condi = 0; end    % disp('hz > 3'); i, end
if tz>3, condi = 0; end    % disp('tz > 3'); i, end

end % while

% dieses Ergebnis kann falsch sein da BSR00 falsch
BX = interp1(H,Btemp,HX,'spline');
if min(BX(2:60)) < BSRAtFit355arr(i), BSR00problem=BSR00problem+1; end

   

BSR355Klett(Sel,i)=Btemp(Sel); BSR355Kletterr(Sel,i)=Btemperr(Sel); 
BetaAer355Klett(Sel,i)=Betaaer(Sel); BetaAer355Kletterr(Sel,i)=Betaaererr(Sel,i); 
attenu355P(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);


dBeta355dP(Sel,i)=dBdP;
dBeta355dR(Sel,i)=dBdR;
dBeta355dLR(Sel,i)=dBdLR;
C355Lidar(Sel,i) = CLidar;

end % if Daten vorhanden

end % if wowolke nicht leer

end  % for Zeitschritte

if BSR00problem==1, BSR00=BSR00+0.05;
elseif BSR00problem >=2, BSR00=BSR00+0.1;
elseif BSR00problem ==0, condip=0;
end
if iterp > itpmax, condip=0; end

end %while BSR00

end   % LRiter


% der Low Bereich aus dem Quotienten der analog Signale
BSR355P=BSR355Klett;   BSR355Perr=BSR355Kletterr;
BSR355PKlett = BSR355Klett; BSR355PKletterr= BSR355Kletterr;
hrc=find(H >2600 & H < 4200); % Heightrangecompare
hnp=find(H>2000); hnp=hnp(1); %Height Null Position
hmp=find(H>3500); hnp=hnp(1); %Height End Position
for jj=1:bis, 
    tmp=mymean(BSR355P(hrc,jj));
    BSR355quot(:,jj)=BSR355quot(:,jj)./mymean(BSR355quot(hrc,jj)).*tmp;
    [m,w]=min(abs(BSR355P(hnp:hmp,jj)-BSR355quot(hnp:hmp,jj)));
    BSR355P(1:hnp-1+w,jj)=BSR355quot(1:hnp-1+w,jj); 
end
% wenn 387 schwach
% BSR355P=BSR355Klett; 
BSR355P=BSR355Klett; 
BSR355Perr=BSR355Kletterr+BSRAtFit355err;

% Gesamtgroessen parallel
[BetaAer355P, BetaAer355Perr]=  bsr2bsc(BSR355P, BSR355Perr, Density, Wvl355, Temp, 'p', EmitPol );  
[BetaAer355Pquot, BetaAer355Pquoterr]=  bsr2bsc(BSR355quot, BSR355Perr, Density, Wvl355, Temp, 'p', EmitPol );  
[BetaAer355Plang, DN355lang] = maketimegaps(BetaAer355P, DN355, 1/24/2);  
[BSR355Plang, DN355lang] = maketimegaps(BSR355P, DN355, 1/24/2);  
[LR355arrl, DN355lang] = maketimegaps(LR355arr, DN355S, 1/24/2);  

% optische Dicke
q=size(BetaAer355P); optdivisp=zeros(q(2),1); bb1=BetaAer355P(Sel,:); bx=find(bb1<0); bb1(bx)=0;
bx=find(isnan(bb1)); bb1(bx)=0;
for jj=1:q(2), optdivisp(jj)=sum(bb1(:,jj).*LR355arr(Sel,jj).*dH)+LR355arr(1,jj).*(H(1)-HStat).*bb1(1,jj); end

AlphaAer355P = BetaAer355P;
Tau355P=zeros(size(BetaAer355P));


for j=1:bis, AlphaAer355P(:,j)=AlphaAer355P(:,j).*LR355arr(:,j); 
Tau355P(:,j)=cumsum(AlphaAer355P(:,j).*Hdiff);
end
%Alpha355Ptot = AlphaAer355P + AlRay355;


%
%
%%%  355nm senkrecht polarisiert
%

if ~isempty(Sel355SC)
    
LR355Sarr = LR355arr; % extin aus "P"
EmitPol='p';RecPol='s'; 
BeRa355S=Density.*raybckwq (Wvl355,'p','s', Temp, Density);  
LR355Sarrerr=ones(size(LR355Sarr))+15;   

BSR355Squot=P355S;
for j=1: size(P355S,2),
    [m,w]=min(abs(DATDateNumStart(Sel387A) - DN355(j)));
    BSR355Squot(:,j)=P355S(:,j)./P387final(:,w);
end
BSR355Surquot = BSR355Squot;

% Signaldefintion
PC=P355SC;
PA=P355SA;
P = PC;
%ber=151:180;
nn=3;  %length(ber);
f1=((1:nn)./nn)';
f2=((nn:-1:1)./nn)';
for j=1: length(Sel355A)
  hwo=150:300; %find(PC(:,j) <10 & PC(:,j)>0.1);
  hwo=find(PC(:,j) <10 & PC(:,j)>0.1 & PA(:,j) >0.1);
  [c,k]=minanpass(PA(:,j),PC(:,j),hwo,H,5);
  if k > 0,
  PA(:,j)=(PA(:,j)-c)./k;
  [x,wo]=min(abs(PA(hwo,j)-PC(hwo,j)));
  hwoend=hwo(wo); ber=hwoend+1:hwoend+3;
  P(1:hwoend,j)=PA(1:hwoend,j);
  P(ber,j)=(f2.*PA(ber,j)+f1.*PC(ber,j))./(f1+f2);
  else
  P(:,j)=P355S(:,j);
  end
end
P = P355S;

[PCerr,SNR355SC] = estinoiseprofile(PC,H,phsoll);
%PCerr=real(sqrt(PC))./50+1e-4; 
wo=find(PC<Pschwelle); PC(wo)=Pschwelle;
[PAerr,SNR355SA] = estinoiseprofile(PA,H,phsoll);
%PAerr=real(sqrt(PA))./50+1e-4; 
wo=find(PA<Pschwelle); PA(wo)=Pschwelle;
[Perr,SNR355SK] = estinoiseprofile(P,H,phsoll);
%Perr=real(sqrt(P))./50+1e-4; 
wo=find(P<Pschwelle); P(wo)=Pschwelle;
P355Sfinal = P; P355Sfinalerr=Perr;
wofit = find( H >FitRangeA(1) & H< FitRangeA(2));

los =1; bis =size(PC,2); %%%length(Sel355C)

Btemp=zeros(Hsize,1); Btemperr=Btemp; 
Betaaer=zeros(Hsize,1); Betaaererr=Btemp; 
BetaAer355SC=zeros(size(PC)); BetaAer355SCerr=zeros(size(PC)); 
dBeta355SCdLR=zeros(size(PC));
dBeta355SCdR=zeros(size(PC)); 
dBeta355SCdP=zeros(size(PC));
BSR355SCKlett=zeros(size(PC)); BSR355SCKletterr=zeros(size(PC));

BetaAer355SA=zeros(size(PA)); BetaAer355SAerr=zeros(size(PA)); 
dBeta355SAdLR=zeros(size(PA));
dBeta355SAdR=zeros(size(PA)); 
dBeta355SAdP=zeros(size(PA));
BSR355SAKlett=zeros(size(PA)); BSR355SAKletterr=zeros(size(PA)); 
BSRFit355SAvektor=zeros(Hsize,1);
attenu355S = zeros(size(PC));

% Schleife ueber alle Zeitschritte: erst counting, dann analog
for i=los:bis,  
    
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
PC(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

aoint=exp(-1.*qdrupvar(H(Sel), AlRay355(Sel,i)));
arint=exp(-1.*qdrupvar(H(Sel), AlRay387(Sel,i)));
BSR355Squot(Sel,i)=BSR355Squot(Sel,i).*arint./aoint;
BSR355Surquot(Sel,i)=BSR355Squot(Sel,i);

% Rechnung counting

[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit355, FitRangeC, H(Sel), PC(Sel,i), PCerr(Sel,i), ...
     LR355Sarr(Sel,i), AlphaAer355P(Sel,i), AlRay355(Sel,i), BeRa355S(Sel,i));
% [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit355S, FitRangeC, H(Sel), PC(Sel,i), PCerr(Sel,i), ...
%     LR355Sarr(Sel,i), AlRay355(Sel,i), BeRa355S(Sel,i));

Betaaer(Sel)=Beta-BeRa355S(Sel,i);
Betaaererr(Sel)=abs(dBdR.*BSRAtFit355err)+abs(dBdLR.*LR355arrerr(Sel,i))+abs(dBdP.*PCerr(Sel,i));
Btemp(Sel)=Beta./BeRa355S(Sel,i); Btemperr(Sel)=Betaaererr(Sel)./BeRa355S(Sel,i); 

BSR355SCKlett(Sel,i)=Btemp(Sel); BSR355SCKletterr(Sel,i)=Btemperr(Sel); 
BetaAer355SC(Sel,i)=Betaaer(Sel); BetaAer355SCerr(Sel,i)=Betaaererr(Sel); 

dBeta355SCdP(Sel,i)=dBdP;
dBeta355SCdR(Sel,i)=dBdR;
dBeta355SCdLR(Sel,i)=dBdLR;

% Rechnung analog
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
PA(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

if ~isempty(Sel)
BSRAtFit355SA=mymean(BSR355SCKlett(wofit,i));   
if isnan(BSRAtFit355SA), BSRAtFit355SA=BSRAtFit355; end
BSRFit355SAvektor(i)=BSRAtFit355SA;

[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit532SA, FitRangeA, H(Sel), PA(Sel,i), PAerr(Sel,i), ...
    LR355Sarr(Sel,i), AlphaAer355P(Sel,i), AlRay355(Sel,i), BeRa355S(Sel,i));
% [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit532SA, FitRangeA, H(Sel), PA(Sel,i), PAerr(Sel,i), ...
%     LR355Sarr(Sel,i), AlRay355(Sel,i), BeRa355S(Sel,i));
% [Beta,dBdR, dBdLR, dBdP] = ...
Betaaer(Sel)=Beta-BeRa355S(Sel,i);
Betaaererr(Sel)=abs(dBdR.*BSRAtFit355err)+abs(dBdLR.*LR355arrerr(Sel,i))+abs(dBdP.*PAerr(Sel,i));
Btemp(Sel)=Beta./BeRa355S(Sel,i); Btemperr(Sel)=Betaaererr(Sel)./BeRa355S(Sel,i); 

BSR355SAKlett(Sel,i)=Btemp(Sel); BSR355SAKletterr(Sel,i)=Btemperr(Sel); 
BetaAer355SA(Sel,i)=Betaaer(Sel); BetaAer355SAerr(Sel,i)=Betaaererr(Sel); 

dBeta355SAdP(Sel,i)=dBdP;
dBeta355SAdR(Sel,i)=dBdR;
dBeta355SAdLR(Sel,i)=dBdLR;
end  % if, Daten vorhanden

end

% Anpassung von analog und counting

BSR355SKlebKlett = BSR355SCKlett;  BSR355SKlebKletterr = BSR355SCKletterr;
hwo=find(H>8000 & H <12000);
nn=length(hwo);
f1=((1:nn)./nn)';
f2=((nn:-1:1)./nn)';
BSR355SKlebKlett(1:hwo(1)-1,:)=BSR355SAKlett(1:hwo(1)-1,:);
BSR355SKlebKletterr(1:hwo(1)-1,:)=BSR355SAKletterr(1:hwo(1)-1,:);
for j=los:bis,
  BSR355SKlebKlett(hwo,j)=(f2.*BSR355SAKlett(hwo,j)+f1.*BSR355SCKlett(hwo,j))./(f1+f2);
  BSR355SKlebKletterr(hwo,j)=(f2.*BSR355SAKletterr(hwo,j)+f1.*BSR355SCKletterr(hwo,j))./(f1+f2);
end

% Rechnung mit Kombinierten Profilen - P
BetaAer355SKlett=zeros(size(P)); BetaAer355SKletterr=zeros(size(P)); 
dBeta355SdLR=zeros(size(P));
dBeta355SdR=zeros(size(P)); 
dBeta355SdP=zeros(size(P));
BSR355SKlett=zeros(size(P)); BSR355SKletterr=zeros(size(P)); 
C355SLidar=zeros(size(P)); 
Perr= real(sqrt(P))./50+1e-4;
for i= los:bis

Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

if length(Sel) > 10
[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit355Sarr(i), FitRangeCS, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR355Sarr(Sel,i), AlphaAer355P(Sel,i), AlRay355(Sel,i), BeRa355S(Sel,i));
Betaaer(Sel)=Beta-BeRa355S(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit355err)+abs(dBdLR.*LR355Sarrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa355S(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa355S(Sel,i); 

BSR355SKlett(Sel,i)=Btemp(Sel); BSR355SKletterr(Sel,i)=Btemperr(Sel); 
BetaAer355SKlett(Sel,i)=Betaaer(Sel); BetaAer355SKletterr(Sel,i)=Betaaererr(Sel,i); 
attenu355S(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);

dBeta355SdP(Sel,i)=dBdP;
dBeta355SdR(Sel,i)=dBdR;
dBeta355SdLR(Sel,i)=dBdLR;
C355SLidar(Sel,i)=CLidar;


end % keine Daten

end




%die entscheidende Iteration


if LRiterS ~=0,

% iterativ
%wohint=40:60;
%wowolke=78:170;z=zeros(1,dimen(2)); 
hwo=find(H>1500 & H< FitRangeC(1));
for j=1:dimen(2), z(j)=mymean(BSR355SKlett(hwo,j)); end
listvalidmeasure=find(DATHV355s>650 & z > 1.015);


FitSel=find(H>FitRangeC(1) & H<FitRangeC(2));
gut = 1; % Prozent-Abweichung
q=zeros(bis,1).*NaN; q2=zeros(bis,1); liste=zeros(bis,1)+FitSel(1);
for j=1:bis, 
    q2(j)=max(BSR355SKlett(1:FitSel(1),j));
    wo = find(BSR355SKlett(1:FitSel(1),j) > BSR355Sschwelle);
    if ~isempty(wo), liste(j)=wo(1)-1; end 
end
[BSR00, wo] = min(q2(listvalidmeasure));
wohint = (Hmin1:min(liste)-1);
if isempty(wohint), 
   wo2=find(liste>Hmin1);
   liste = liste(wo2); 
%    [m,w]=min(liste);
%    liste(w)=NaN;
end
wohint = (Hmin1:min(liste));
for j=1:bis,
   q(j)=mymean(BSR355SKlett(wohint,j));
end

BSR00 = q(wo);

Btemp2 = Btemp;

% wowolke = find(BSR355SKlett > BSR355Sschwelle);
% LR355Sarr(wowolke) = val_LR355SWo;
iterp=0; condip=1; itpmax=2;
while condip
iterp=iterp+1;
LR355Sarr(:) = val_LR355S;
BSR00problem=0;

for i= 1:bis,

wowolke = find(BSR355SKlett(1:FitSel(1),i) > BSR355Sschwelle);
LR355Sarr(wowolke,i) = val_LR355SWo;
if ~isempty(wowolke)
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

if length(Sel) >10,
    
iter=0; condi=1;  hz=0; tz=0;  % hz / tz Zaehler fuer hohe - tiefe Werte
refer=zeros(itmax,1);  refer2=zeros(itmax,1);
while condi
iter = iter+1;    

% [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit355Sarr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
%     LR355Sarr(Sel,i), AlRay355(Sel,i), BeRa355S(Sel,i));

[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4s( BSRAtFit355Sarr(i), FitRangeCS, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR355Sarr(Sel,i), AlphaAer355P(Sel,i), AlRay355(Sel,i), BeRa355S(Sel,i));
Betaaer(Sel)=Beta-BeRa355S(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit355err)+abs(dBdLR.*LR355arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa355S(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa355S(Sel,i); 

% [Beta2, dB2dR, dB2dLR, dB2dP, CLidar2] = klettinv_ableit4( BSRAtFit355Sarr(i), FitRangeC, H(Sel), P(Sel,i), Perr(Sel,i), ...
%     LR355Sarr(Sel,i)+1, AlRay355(Sel,i), BeRa355S(Sel,i));
[Beta2, dB2dR, dB2dLR, dB2dP, CLidar2] = klettinv_ableit4s( BSRAtFit355Sarr(i), FitRangeCS, H(Sel), P(Sel,i), Perr(Sel,i), ...
    (LR355Sarr(Sel,i)+1), AlphaAer355P(Sel,i), AlRay355(Sel,i), BeRa355S(Sel,i));
Btemp2(Sel)=Beta2./BeRa355S(Sel,i); % Btemperr(Sel)=Betaaererr(Sel,i)./BeRa355S(Sel,i); 

refer(iter) = mymean(Btemp(wohint));
refer2(iter) = mymean(Btemp2(wohint));
dxdlr = refer2(iter)-refer(iter);
diffi = BSR00 - refer(iter);
deltax=diffi./dxdlr;


if (refer(iter) ./ BSR00 > 1+gut ./ 100),          % BSR zu gro?
    LR355Sarr(wowolke,i) = LR355Sarr(wowolke,i) +deltax;
    if mymean(LR355Sarr(wowolke,i)) > LRSobergr, LR355Sarr(wowolke,i) = LRSobergr; hz=hz+1; end
    if mymean(LR355Sarr(wowolke,i)) < LRSuntergr, LR355Sarr(wowolke,i) = LRSuntergr; tz=tz+1; end
elseif (refer(iter) ./ BSR00 < 1-gut ./ 100),      % BSR zu klein
    LR355Sarr(wowolke,i) = LR355Sarr(wowolke,i) +deltax;
    if mymean(LR355Sarr(wowolke,i)) > LRSobergr, LR355Sarr(wowolke,i) = LRSobergr; hz=hz+1; end
    if mymean(LR355Sarr(wowolke,i)) < LRSuntergr, LR355Sarr(wowolke,i) = LRSuntergr; tz=tz+1; end
else
    condi=0;
end

if iter >=itmax, condi=0; end   % disp('keine Konvergenz'); i, end
if hz>3, condi = 0; end   % disp('hz > 3'); i, end
if tz>3, condi = 0; end   % disp('tz > 3'); i, end

end % while
  
% dieses Ergebnis kann falsch sein da BSR00 falsch
BX = interp1(H,Btemp,HX,'spline');
if min(BX(2:60)) < BSRAtFit355Sarr(i), BSR00problem=BSR00problem+1; end

BSR355SKlett(Sel,i)=Btemp(Sel); BSR355SKletterr(Sel,i)=Btemperr(Sel); 
BetaAer355SKlett(Sel,i)=Betaaer(Sel); BetaAer355SKletterr(Sel,i)=Betaaererr(Sel,i); 
attenu355S(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);


dBeta355SdP(Sel,i)=dBdP;
dBeta355SdR(Sel,i)=dBdR;
dBeta355SdLR(Sel,i)=dBdLR;
C355SLidar(Sel,i)=CLidar;

end % if Daten vorhanden

end % wowolke nicht leer


end  % for Zeitschritte

if BSR00problem==1, BSR00=BSR00+0.05;
elseif BSR00problem >=2, BSR00=BSR00+0.1;
elseif BSR00problem ==0, condip=0;
end
if iterp > itpmax, condip=0; end

end %while BSR00


end   %LRiter


% der Low Bereich aus dem Quotienten der analog Signale
BSR355S=BSR355SKlett;   BSR355Serr=BSR355SKletterr;
hrc=find(H >2600 & H < 4200); % Heightrangecompare
hnp=find(H>2000); hnp=hnp(1); %Height Null Position
hmp=find(H>3500); hnp=hnp(1); %Height End Position
for jj=1:bis, 
    tmp=mymean(BSR355S(hrc,jj));
    BSR355Squot(:,jj)=BSR355Squot(:,jj)./mymean(BSR355Squot(hrc,jj)).*tmp;
    [m,w]=min(abs(BSR355S(hnp:hmp,jj)-BSR355Squot(hnp:hmp,jj)));
    BSR355S(1:hnp-1+w,jj)=BSR355Squot(1:hnp-1+w,jj); 
end
% wenn  387 schwach:
% BSR355S=BSR355SKlett; 
BSR355S=BSR355SKlett; 
BSR355Serr=BSR355SKletterr+BSRAtFit355err;

% Defintion Gesamtgroessen
[BetaAer355S, BetaAer355Serr]=  bsr2bsc(BSR355S, BSR355Serr, Density, Wvl355, Temp, 's', EmitPol );  
[BetaAer355Slang, DN355Slang] = maketimegaps(BetaAer355S, DN355S, 1/24/2);  
[BSR355Slang, DN355Slang] = maketimegaps(BSR355S, DN355S, 1/24/2);  
[BetaAer355Squot, BetaAer355Squoterr]=  bsr2bsc(BSR355Squot, BSR355Serr, Density, Wvl355, Temp, 's', EmitPol );  


% optische Dicke
q=size(BetaAer355S); optdiUVS=zeros(q(2),1); bb1=BetaAer355S(Sel,:); bx=find(bb1<0); bb1(bx)=0;
bx=find(isnan(bb1)); bb1(bx)=0;
for jj=1:q(2), optdiUVS(jj)=sum(bb1(:,jj).*LR355Sarr(Sel,jj).*dH)+LR355Sarr(1,jj).*(H(1)-HStat).*bb1(1,jj); end

AlphaAer355S = BetaAer355S;
Tau355S=zeros(size(BetaAer355S));


for j=1:bis, AlphaAer355S(:,j)=AlphaAer355S(:,j).*LR355Sarr(:,j); 
Tau355S(:,j)=cumsum(AlphaAer355S(:,j).*Hdiff);
aaerint=exp(-qdrupvar(H(Sel),AlphaAer355S(Sel,i)- (1).*AlphaAer355P(Sel,i)));
BSR355Squot(Sel,i)=BSR355Surquot(Sel,i)./aaerint;
end

%%%Zusammenfassung von 355 parallel und 355 senkrecht

% 355 total zusammenfassen p und s koennen unterschiedlich lang sein

BetaAer355tot = zeros(size(P355A));
BetaAer355toterr = zeros(size(P355A));
BetaRay355tot = zeros(size(P355A));
BSR355tot = zeros(size(P355A));
BSR355toterr = zeros(size(P355A));
AeroDep355 = zeros(size(P355A));
AeroDep355err = zeros(size(P355A));
los =1; bis =length(Sel355A);
for j=los:bis,
k= find(DATDateNumStart(Sel355SA) == DATDateNumStart(Sel355A(j)));
if length(k) ==1,  
   BetaAer355tot(:,j) = BetaAer355S(:,k) + BetaAer355P(:,j);
   BetaAer355toterr(:,j) = BetaAer355Serr(:,k) + BetaAer355Perr(:,j);
   BetaRay355tot(:,j)   = BeRa355S(:,k) + BeRa355(:,j);
   BSR355tot(:,j)     = BetaAer355tot(:,j) ./ BetaRay355tot(:,j)+1;
   BSR355toterr(:,j)  = BetaAer355toterr(:,j) ./ BetaRay355tot(:,j);
   AeroDep355(:,j) = BetaAer355S(:,k) ./ BetaAer355P(:,j);
   AeroDep355err(:,j) = BetaAer355Serr(:,k) ./ BetaAer355P(:,j) + ...
                        AeroDep355(:,j).*BetaAer355Perr(:,j) ./ BetaAer355P(:,j);
end
end
else % no depol data 355
    BetaAer355tot=[]; BetaAer355S=[]; BSR355tot=[]; BSR355SKlett=[];
    BetaAer355toterr=[]; BSR355toterr=[]; AeroDep355=[]; AeroDep355err=[];
end %if isempty Sel355SC




LR355Raman = zeros(size(P387C));   LR355Ramanerr=zeros(size(P387C)); 
los =1; bis =length(Sel387C);
for j=los:bis,
k= find(DATDateNumStart(Sel355C) == DATDateNumStart(Sel387C(j)));
if length(k) ==1,  
   LR355Raman(:,j)= AlphaAer355Raman(:,j) ./ BetaAer355tot(:,k);
   LR355Ramanerr(:,j)= abs(AlphaAer355Ramanerr(:,j) ./BetaAer355tot(:,k))+ ...
       abs(AlphaAer355Raman(:,j) .*BetaAer355toterr(:,k) ./BetaAer355tot(:,k).^2);
end
end


% LR355schaetz = zeros(size(P387x));   LR355schaetzerr=zeros(size(P387x)); 
% los =1; bis =length(Sel387C);
% for j=los:bis,
% k= find(DATDateNumStart(Sel355C) == DATDateNumStart(Sel387C(j)));
% if length(k) ==1,  
%    BetaAer355totx(:,k) = interp1(H,BetaAer355tot(:,k), Hx387);
%    BetaAer355totxerr(:,k) = interp1(H,BetaAer355toterr(:,k), Hx387);
%    h1 = Alpha355schaetz(:,j) - Alpha355schaetzMin(:,j);
%    h2 = Alpha355schaetzMax(:,j) - Alpha355schaetz(:,j);
%    Alpha355schaetzerr =  max(h1,h2);
%    LR355schaetz(:,j)= Alpha355schaetz(:,j) ./ BetaAer355totx(:,k);
%    LR355schaetzerr(:,j)= abs(Alpha355schaetzerr ./BetaAer355totx(:,k))+ ...
%        abs(Alpha355schaetz(:,j) .*BetaAer355totxerr(:,k) ./BetaAer355totx(:,k).^2);
% end
% end
% 

crk=zeros(size(BSR532tot)).*NaN;
wo=find(BSR532tot > 1.08 & BSR355tot > 1.04);
crk(wo)=BetaAer355tot(wo) ./ BetaAer532tot(wo);
wo = find(crk >5); crk(wo) =NaN;

hstartmeter=700;
hstart=find(H>=hstartmeter);   hstart=hstart(1);
RamFitRange=[14000 17000];
schwelle=1.0e3;
ph387=1./4094./8;
[TauVIStot, TauVISglatt,TgminVIS, TgmaxVIS, hschritteVIS] = tauschritte3(P607,H,hstart, Density(:,:), ...
    AlRay532(:,:), AlRay607(:,:), RamFitRange, schwelle, ph387, 0);

[taufit, alphafit, gradarr] = calctauundalpha2(hschritteVIS,TauVISglatt,33);
[tauminfit, alphaminfit, gradarr] = calctauundalpha2(hschritteVIS,TgminVIS,33);
[taumaxfit, alphamaxfit, gradarr] = calctauundalpha2(hschritteVIS,TgmaxVIS,33);
AngarrVIS=WvlFct.*ones(size(taufit));
atest=WvlFct.*ones(size(crk));
aw=0.3; % Angstroem Wolke
 for k=1:dimen(2)
   wo=find( crk(:,k) < 1.5 & BSR532tot(:,k) >8);
   atest(wo,k)=aw;
   AngarrVIS(:,k)=interp1(H,atest(:,k),hschritteVIS(:,k));
 end
f1a=1+(Wvl532./WvlR532).^AngarrVIS;
f2a=1+(Wvl532./WvlR532).^(-AngarrVIS);
Tau607final=taufit./f2a;
Tau532final=taufit./f1a;
Alpha607final=alphafit./f2a;
Alpha532final=alphafit./f1a;
Tau607minfinal=tauminfit./f2a;
Tau532minfinal=tauminfit./f1a;
Alpha607minfinal=alphaminfit./f2a;
Alpha532minfinal=alphaminfit./f1a;
Tau607maxfinal=taumaxfit./f2a;
Tau532maxfinal=taumaxfit./f1a;
Alpha607maxfinal=alphamaxfit./f2a;
Alpha532maxfinal=alphamaxfit./f1a;



[TauUVtot, TauUVglatt,TgminUV, TgmaxUV, hschritteUV] = tauschritte3(P387(:,:),H, hstart, Density(:,:), AlRay355(:,:), ...
    AlRay387(:,:), RamFitRange, schwelle, ph387, 0);

[taufit, alphafit, gradarr] = calctauundalpha2(hschritteUV,TauUVglatt,33);
[tauminfit, alphaminfit, gradarr] = calctauundalpha2(hschritteUV,TgminUV,33);
[taumaxfit, alphamaxfit, gradarr] = calctauundalpha2(hschritteUV,TgmaxUV,33);
AngarrUV=WvlFct.*ones(size(taufit));
atest=WvlFct.*ones(size(crk));
aw=0.3; % Angstroem Wolke
for k=1:dimen(2)
   wo=find( crk(:,k) < 1.5 & BSR532tot(:,k) >8);
   atest(wo,k)=aw;
   AngarrUV(:,k)=interp1(H,atest(:,k),hschritteUV(:,k));
 end
f1=1+(Wvl355./WvlR355).^AngarrUV;
f2=1+(Wvl355./WvlR355).^(-AngarrUV);
Tau387final=taufit./f2;
Tau355final=taufit./f1;
Alpha387final=alphafit./f2;
Alpha355final=alphafit./f1;
Tau387minfinal=tauminfit./f2;
Tau355minfinal=tauminfit./f1;
Alpha387minfinal=alphaminfit./f2;
Alpha355minfinal=alphaminfit./f1;
Tau387maxfinal=taumaxfit./f2;
Tau355maxfinal=taumaxfit./f1;
Alpha387maxfinal=alphamaxfit./f2;
Alpha355maxfinal=alphamaxfit./f1;










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Rechnung fuer 1064nm
Wvl1064=Wvl532.*2; 
LR1064arrerr=zeros(size(LR1064arr))+15;   



%%%1064nm 

EmitPol='p';RecPol='u'; 
BeRa1064=Density.*raybckwq (Wvl1064,'p','u', Temp, Density);  
AlRay1064ur = Density .* raytotwq( Wvl1064, Temp, Density);
AlRay1064 = zeros(size(AlRay1064ur));
phsollIR = phsoll ./ 0.08;

% Signaldefinition
P=P1064;
los =1; bis =size(P,2); %%%length(Sel1064C)
ber=find(H>39000 & H<60000); ber = ber';
for j=los:bis
 [c,k]=minanpass(P(:,j),P(:,j),ber, H, 2);
 P(:,j) = P(:,j) -k.*H-c;
end

[Perr,SNR1064] = estinoiseprofile(P,H,phsollIR);
%Perr=real(sqrt(P))./20+0.05; 
wo=find(P<Pschwelle); P(wo)=Pschwelle;

P1064final = P; P1064finalerr=Perr;


Btemp=zeros(Hsize,1); Btemperr=Btemp; 
Betaaer=zeros(Hsize,1); Betaaererr=Btemp; 
BetaAer1064=zeros(size(P)); BetaAer1064err=zeros(size(P)); 
dBeta1064dLR=zeros(size(P));
dBeta1064dR=zeros(size(P)); 
dBeta1064dP=zeros(size(P));
BSR1064Klett=zeros(size(P)); BSR1064Kletterr=zeros(size(P));
C1064Lidar=zeros(size(P));
attenu1064 = zeros(size(P));

% Schleife ueber alle Zeitschritte    1. Durchgang: Test

%%%%%% brauchen Test nicht mehr da Randbeding. 1064 aus 532 genommen wird

% for i=los:bis,  
%     
% Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
% P(:,i) > 0 & ...
% Density(:,i) > 0, 1);
% if length(Sel) > 1,
%         Sel = (Sel(1):Sel(2));
%     else
%         Sel = [];
% end
% 
% ARayad=O3Density.*o3abswq(Wvl1064,Temp(:,i)); 
% AlRay1064(:,i)=AlRay1064ur(:,i)+ARayad; 
% 
% % Rechnung 
% 
% % [Beta,dBdR, dBdLR, dBdP] = ...
% % klettinv_ableit2(  BSRAtFit1064arr(i), FitRangeA, [],H(Sel), HStat, P(Sel,i), Perr(Sel,i), ...
% % Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR1064arr(Sel,i), LR1064arrerr(Sel,i), Wvl1064, EmitPol, RecPol, ...
% % Temp(Sel,i));
% [Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit1064arr(i), FitRangeA, H(Sel), P(Sel,i), Perr(Sel,i), ...
%     LR1064arr(Sel,i), AlRay1064(Sel,i), BeRa1064(Sel,i));
% Betaaer(Sel)=Beta-BeRa1064(Sel,i);
% Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit1064err)+abs(dBdLR.*LR1064arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
% Btemp(Sel)=Beta./BeRa1064(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa1064(Sel,i); 
% 
% BSR1064Klett(Sel,i)=Btemp(Sel); BSR1064Kletterr(Sel,i)=Btemperr(Sel); 
% BetaAer1064(Sel,i)=Betaaer(Sel); BetaAer1064err(Sel,i)=Betaaererr(Sel,i); 
% 
% dBeta1064dP(Sel,i)=dBdP;
% dBeta1064dR(Sel,i)=dBdR;
% dBeta1064dLR(Sel,i)=dBdLR;
% C1064Lidar(Sel,i)=CLidar;
% 
% end

%%%% die Randbedingung k?nnte falsch sein, da da? IR Signal schon in H?hen,
%%%% die noch aerosolhaltig sind verrauscht wird, muessen also Rechnung
%%%% wiederholen mit Beding, da? R>RFit in Troposphaere
 wo=find(H>FitRangeA(1) & H<FitRangeA(2));
 BSRfitfaktor=2.^(4+WvlFct);
 for j=los:bis,
   q(j)=mymean(BSR532tot(wo,j));
   if q(j) > 1.08, BSRAtFit1064arr(j) = 1+(q(j)-1).*BSRfitfaktor; end
 end
 
%%% und nu 2. Durchgang, richtiges BSRAtFit

for i=los:bis,  
    
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

ARayad=O3Density.*o3abswq(Wvl1064,Temp(:,i)); 
AlRay1064(:,i)=AlRay1064ur(:,i)+ARayad; 
% Rechnung 

% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2( BSRAtFit1064arr(i), FitRangeA, [],H(Sel), HStat, P(Sel,i), Perr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR1064arr(Sel,i), LR1064arrerr(Sel,i), Wvl1064, EmitPol, RecPol, ...
% Temp(Sel,i));
[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit1064arr(i), FitRangeA, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR1064arr(Sel,i), AlRay1064(Sel,i), BeRa1064(Sel,i));
Betaaer(Sel)=Beta-BeRa1064(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit1064err)+abs(dBdLR.*LR1064arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa1064(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa1064(Sel,i); 

BSR1064Klett(Sel,i)=Btemp(Sel); BSR1064Kletterr(Sel,i)=Btemperr(Sel); 
BetaAer1064(Sel,i)=Betaaer(Sel); BetaAer1064err(Sel,i)=Betaaererr(Sel,i); 
attenu1064(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);


dBeta1064dP(Sel,i)=dBdP;
dBeta1064dR(Sel,i)=dBdR;
dBeta1064dLR(Sel,i)=dBdLR;
C1064Lidar(Sel,i)=CLidar;

end


%
% die entscheidende Iteration
%

if LRiter ~=0,

% iterativ
%wohint=40:60;
%wowolke=78:170;
z=zeros(1,dimen(2)); 
hwo=find(H>1500 & H< FitRangeC(1));
for j=1:dimen(2), z(j)=mymean(BSR1064Klett(hwo,j)); end
listvalidmeasure=find(DATHV1064>200 & z > 1.3);

FitSel=find(H>FitRangeC(1) & H<FitRangeC(2));
gut = 1; % Prozent-Abweichung
q=zeros(bis,1).*NaN; q2=zeros(bis,1); liste=zeros(bis,1)+FitSel(1);
for j=1:bis, 
    q2(j)=max(BSR1064Klett(1:FitSel(1),j));
    wo = find(BSR1064Klett(1:FitSel(1),j) > BSR1064schwelle);
    if ~isempty(wo), liste(j)=wo(1)-1; end 
end
[BSR00, wo] = min(q2(listvalidmeasure));
wohint = (Hmin1:min(liste)-1);
if isempty(wohint), 
   wo2=find(liste>Hmin1);
   liste = liste(wo2); 
%    [m,w]=min(liste);
%    liste(w)=NaN;
end
wohint = (Hmin1:min(liste));
for j=1:bis,
   q(j)=mymean(BSR1064Klett(wohint,j));
end

BSR00 = max([q(wo),1.4]);

Btemp2 = Btemp;

% wowolke = find(BSR1064Klett > BSR1064schwelle);
% LR1064arr(wowolke) = val_LR1064Wo;
iterp=0; condip=1; itpmax=5;
while condip
iterp=iterp+1;
LR1064arr(:) = val_LR1064;
BSR00problem=0;

for i= los:bis

wowolke = find(BSR1064Klett(1:FitSel(1),i) > BSR1064schwelle);
LR1064arr(wowolke,i) = val_LR1064Wo;
if ~isempty(wowolke)    
Sel = connrnge( H >= 0 & H <= Hcalcrange & ...
P(:,i) > 0 & ...
Density(:,i) > 0, 1);
if length(Sel) > 1,
        Sel = (Sel(1):Sel(2));
    else
        Sel = [];
end

if length(Sel) >10,
    
iter=0; condi=1;  hz=0; tz=0;  % hz / tz Zaehler fuer hohe - tiefe Werte
refer=zeros(itmax,1);  refer2=zeros(itmax,1);
while condi
iter = iter+1;    
% [Beta,dBdR, dBdLR, dBdP] = ...
% klettinv_ableit2( BSRAtFit1064arr(i), FitRangeC, DN1064(i),H(Sel), HStat, P(Sel,i), Perr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR1064arr(Sel,i), LR1064arrerr(Sel,i), Wvl1064, 'p','u', ...
% Temp(Sel,i));
[Beta, dBdR, dBdLR, dBdP, CLidar] = klettinv_ableit4( BSRAtFit1064arr(i), FitRangeA, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR1064arr(Sel,i), AlRay1064(Sel,i), BeRa1064(Sel,i));
Betaaer(Sel)=Beta-BeRa1064(Sel,i);
Betaaererr(Sel,i)=abs(dBdR.*BSRAtFit1064err)+abs(dBdLR.*LR1064arrerr(Sel,i))+abs(dBdP.*Perr(Sel,i));
Btemp(Sel)=Beta./BeRa1064(Sel,i); Btemperr(Sel)=Betaaererr(Sel,i)./BeRa1064(Sel,i); 

% [Beta2,dB2dR, dB2dLR, dB2dP] = ...
% klettinv_ableit2( BSRAtFit1064arr(i), FitRangeC, DN1064(i),H(Sel), HStat, P(Sel,i), Perr(Sel,i), ...
% Density(Sel,i), [], O3Density(Sel), O3DensityAE(Sel), LR1064arr(Sel,i)+1, LR1064arrerr(Sel,i), Wvl1064, 'p','u', ...
% Temp(Sel,i));
[Beta2, dB2dR, dB2dLR, dB2dP, CLidar2] = klettinv_ableit4( BSRAtFit1064arr(i), FitRangeA, H(Sel), P(Sel,i), Perr(Sel,i), ...
    LR1064arr(Sel,i)+1, AlRay1064(Sel,i), BeRa1064(Sel,i));

Btemp2(Sel)=Beta2./BeRa1064(Sel,i); % Btemperr(Sel)=Betaaererr(Sel,i)./BeRa1064(Sel,i); 

refer(iter) = mymean(Btemp(wohint));
refer2(iter) = mymean(Btemp2(wohint));
dxdlr = refer2(iter)-refer(iter);
diffi = BSR00 - refer(iter);
deltax=diffi./dxdlr;


if refer(iter) ./ BSR00 > 1+gut ./ 100,          % BSR zu gro?
    LR1064arr(wowolke,i) = LR1064arr(wowolke,i) +deltax;
    if mymean(LR1064arr(wowolke,i)) > LRobergr, LR1064arr(wowolke,i) = LRobergr; hz=hz+1; end
    if mymean(LR1064arr(wowolke,i)) < 5, LR1064arr(wowolke,i) = 5; tz=tz+1; end
elseif refer(iter) ./ BSR00 < 1-gut ./ 100,      % BSR zu klein
    LR1064arr(wowolke,i) = LR1064arr(wowolke,i) +deltax;
    if mymean(LR1064arr(wowolke,i)) > LRobergr, LR1064arr(wowolke,i) = LRobergr; hz=hz+1; end
    if mymean(LR1064arr(wowolke,i)) < LRuntergr, LR1064arr(wowolke,i) = LRuntergr; tz=tz+1; end
else
    condi=0;
end

if iter >=itmax, condi=0; end     % disp('keine Konvergenz'); i, end
if hz>3, condi = 0; end      % disp('hz > 3'); i, end
if tz>3, condi = 0; end      % disp('tz > 3'); i, end

end % while
    
% dieses Ergebnis kann falsch sein da BSR00 falsch
BX = interp1(H,Btemp,HX,'spline');
if min(BX(2:60)) < BSRAtFit1064arr(i), BSR00problem=BSR00problem+1; end


BSR1064Klett(Sel,i)=Btemp(Sel); BSR1064Kletterr(Sel,i)=Btemperr(Sel); 
BetaAer1064(Sel,i)=Betaaer(Sel); BetaAer1064err(Sel,i)=Betaaererr(Sel,i); 
attenu1064(Sel,i) = P(Sel,i).*H(Sel).^2 ./ CLidar(1);


dBeta1064dP(Sel,i)=dBdP;
dBeta1064dR(Sel,i)=dBdR;
dBeta1064dLR(Sel,i)=dBdLR;
C1064Lidar(Sel,i)=CLidar;

end % if Daten vorhanden

end  % wowolke nicht leer

end % for Zeitschritte

if BSR00problem==1, BSR00=BSR00+0.05;
elseif BSR00problem >=2, BSR00=BSR00+0.1;
elseif BSR00problem ==0, condip=0;
end
if iterp > itpmax, condip=0; end

end %while BSR00

end    %LRiter

BSR1064 = BSR1064Klett;  BSR1064err= BSR1064Kletterr;
 
[BetaAer1064lang, DN1064lang] = maketimegaps(BetaAer1064, DN1064, 1/24/2);  
[BSR1064lang, DN1064lang] = maketimegaps(BSR1064Klett, DN1064, 1/24/2);

% optische Dicke
q=size(BetaAer1064); optdiIR=zeros(q(2),1); bb1=BetaAer1064(Sel,:); bx=find(bb1<0); bb1(bx)=0;
bx=find(isnan(bb1)); bb1(bx)=0;
for jj=1:q(2), optdiIR(jj)=sum(bb1(:,jj).*LR1064arr(Sel,jj).*dH)+LR1064arr(1,jj).*(H(1)-HStat).*bb1(1,jj); end

AlphaAer1064 = BetaAer1064;
Tau1064=zeros(size(BetaAer1064));


for j=1:bis, AlphaAer1064(:,j)=AlphaAer1064(:,j).*LR1064arr(:,j); 
Tau1064(:,j)=cumsum(AlphaAer1064(:,j).*Hdiff);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Volumendepolarisation

% wo ist Aerosoldepol+0
FitSelPol=find(H>PolFitRange(1) & H<PolFitRange(2));

% 532nm
% analog
VD532A = zeros(size(P532A));
los =1; bis =length(Sel532A);
for j=los:bis,
k= find(DATDateNumStart(Sel532SA) == DATDateNumStart(Sel532A(j)));
if length(k) ==1,  
   VD532A(:,j)= P532SA(:,k) ./ P532A(:,j);
   VD532A(:,j) = VD532A(:,j)./mymean(VD532A(FitSelPol,j)).* 0.0143;
end
end

% counting
VD532C = zeros(size(P532C));
los =1; bis =length(Sel532C);
for j=los:bis,
k= find(DATDateNumStart(Sel532SC) == DATDateNumStart(Sel532C(j)));
if length(k) ==1,  
   VD532C(:,j)= P532SC(:,k) ./ P532C(:,j);
   VD532C(:,j) = VD532C(:,j)./mymean(VD532C(FitSelPol,j)).* 0.0143;
end
end

% zusammensetzen
VD532 = VD532A;
for j=los:bis,
 xx=find(H >7000 & H<12000);
 hilf=VD532C(xx,j);
 tmp=(abs(hilf-VD532A(xx,j))./hilf);
 [m,w] = min(tmp); 
 VD532(xx(1)+w:end,j) = VD532C(xx(1)+w:end,j);    
end

% aus den kombinierten Profilen:
VD532 = zeros(size(P532final));
VD532err = zeros(size(P532final));
los =1; bis =length(Sel532C);
for j=los:bis,
k= find(DATDateNumStart(Sel532SC) == DATDateNumStart(Sel532C(j)));
if length(k) ==1,  
   VD532(:,j)= P532Sfinal(:,k) ./ P532final(:,j);
   tmp = mymean(VD532(FitSelPol,j))./ 0.0143;
   VD532(:,j) = VD532(:,j)./tmp;
   VD532err(:,j) = 0.1.*VD532(:,j)./tmp.^2 + VD532(:,j)./tmp./P532final(:,j).*P532finalerr(:,j) ...
       + 1./tmp./P532Sfinal(:,j).*P532Sfinalerr(:,k);
end
end




% 355nm
% analog
VD355A = zeros(size(P355A));
los =1; bis =length(Sel355A);
for j=los:bis,
k= find(DATDateNumStart(Sel355SA) == DATDateNumStart(Sel355A(j)));
if length(k) ==1,  
   VD355A(:,j)= P355SA(:,k) ./ P355A(:,j);
   VD355A(:,j) = VD355A(:,j)./mymean(VD355A(FitSelPol,j)).* 0.0143;
end
end

% counting
VD355C = zeros(size(P355C));
los =1; bis =length(Sel355C);
for j=los:bis,
k= find(DATDateNumStart(Sel355SC) == DATDateNumStart(Sel355C(j)));
if length(k) ==1,  
   VD355C(:,j)= P355SC(:,k) ./ P355C(:,j);
   VD355C(:,j) = VD355C(:,j)./mymean(VD355C(FitSelPol,j)).* 0.0143;
end
end

% zusammensetzen
VD355 = VD355A;
for j=los:bis,
 xx=find(H >7000 & H<12000);
 hilf=VD355C(xx,j);
 tmp=(abs(hilf-VD355A(xx,j))./hilf);
 [m,w] = min(tmp); 
 VD355(xx(1)+w:end,j) = VD355C(xx(1)+w:end,j);    
end

% aus den kombinierten Profilen:
VD355 = zeros(size(P355final));
VD355err = zeros(size(P355final));
los =1; bis =length(Sel355C);
for j=los:bis,
k= find(DATDateNumStart(Sel355SC) == DATDateNumStart(Sel355C(j)));
if length(k) ==1,  
   VD355(:,j)= P355Sfinal(:,k) ./ P355final(:,j);
   tmp = mymean(VD355(FitSelPol,j))./ 0.0143;
   VD355(:,j) = VD355(:,j)./tmp;
   VD355err(:,j) = 0.1.*VD355(:,j)./tmp.^2 + VD355(:,j)./tmp./P355final(:,j).*P355finalerr(:,j) ...
       + 1./tmp./P355Sfinal(:,j).*P355Sfinalerr(:,k);
end
end




% wasser
%
% DNwUV   = DATDateNumStart(Sel407C);
% DNwVis  = DATDateNumStart(Sel660);
% DNRUV   = DATDateNumStart(Sel387C);
% DNRVis  = DATDateNumStart(Sel607C);


no=size(P407,2); wasserUV=zeros(size(P407));wasserUVprop=wasserUV;
T1=exp(-1.*qdrupvar(H(HSel),AlRay387(HSel,1)));
Al407 = Density .* raytotwq( 4.07e-7, Temp, Density);
T2=exp(-1.*qdrupvar(H(HSel),Al407(HSel,1)));
eouep = 5; %end of ueberlapppos
for j=1:no,
    q=0;  [m,q]=min(abs(DN407(j) - DN387));  if length(q) >1, q=q(1); disp(j);  end
     if q~=0,
       wasserUV(HSel,j) = P407(HSel,j) ./ P387(HSel,q).*T1./T2;  
       aaer_387=BetaAer355tot(HSel,q).*LR355arr(HSel,q).*(387./355).^(WvlFct);
       aaer_407=BetaAer355tot(HSel,q).*LR355arr(HSel,q).*(407./355).^(WvlFct);
       wo = find(isnan(aaer_387)); aaer_387(wo)=0;
       wo = find(isnan(aaer_407)); aaer_407(wo)=0;
       T1a=exp(-1.*qdrupvar(H(HSel),aaer_387(HSel)));
       T2a=exp(-1.*qdrupvar(H(HSel),aaer_407(HSel)));
       tvv= T1a./T2a; tvv=tvv./tvv(eouep);
       wasserUVprop(HSel,j) = wasserUV(HSel,j).*tvv;   % proportional - also unkalibriert
     end
   if q ~=j,  disp(j); disp(q);  end
end

no=size(P660,2); wasservis=zeros(size(P660));wasservisprop=wasservis;
T1=exp(-1.*qdrupvar(H(HSel),AlRay607(HSel,1)));
Al660 = Density .* raytotwq( 6.6e-7, Temp, Density);
T2=exp(-1.*qdrupvar(H(HSel),Al660(HSel,1)));
eouep = 5; %end of ueberlapppos
for j=1:no,
    q=0;  [m,q]=min(abs(DN660(j) - DN607));  if length(q) >1, q=q(1); disp(j);  end
     if q~=0,
       wasservis(HSel,j) = P660(HSel,j) ./ P607(HSel,q).*T1./T2;  
       aaer_607=BetaAer532tot(HSel,q).*LR532arr(HSel,q).*(607./532).^(WvlFct);
       aaer_660=BetaAer532tot(HSel,q).*LR532arr(HSel,q).*(660./532).^(WvlFct);
       wo = find(isnan(aaer_607)); aaer_607(wo)=0;
       wo = find(isnan(aaer_660)); aaer_660(wo)=0;
       T1a=exp(-1.*qdrupvar(H(HSel),aaer_607(HSel)));
       T2a=exp(-1.*qdrupvar(H(HSel),aaer_660(HSel)));
       tvv= T1a./T2a; tvv=tvv./tvv(eouep);
       wasservisprop(HSel,j) = wasservis(HSel,j).*tvv;  % proportional
     end
   if q ~=j,  disp(j); disp(q);  end
end


% ist Kalibration Wasser moeglich?
% 
wasserviskali = NaN; wasserUVkali = NaN; wo=[];
FitSelC=find(H>FitRangeC(1) & H<FitRangeC(2));
% UV
diffi=24.*(abs(ballonzeit - DN407)); % in Stunden
wo=find(diffi<1);
if length(wo) >=3 % ok mind. 3 Lidarprofile innerh. 1h zu Sonde
    norm = zeros(length(wo),1);
    for j = 1: length(wo);
     SNR407 = P407C(:,wo(j))./real(sqrt(P407C(:,wo(j))));
     wok = find(SNR407(FitSelC) >0);
     rausch = mymean(SNR407(FitSelC(wok)));
     gut=find(H>400 & SNR407>10.*rausch);
     if length(gut) >5  % ok mind. 5 Posit. zum Vgl
         norm(j) = sum(mischvv(gut)) ./sum(wasserUVprop(gut,wo(j)));
     end
    end   % for wo 
if all(norm) % alle Elemente ~=0
   var=max(norm)./min(norm);
   if var < 1.25  % alles das selbe, nun normieren
     wasserUVkali = mymean(norm);
     disp( 'normieren UV')
   end % homog Daten
end % Normierung zu allen Zeiten
end % length(wo)

% vis 
diffi=24.*(abs(ballonzeit - DN660)); % in Stunden
wo=find(diffi<1);
if length(wo) >=3 % ok mind. 3 Lidarprofile innerh. 1h zu Sonde
    norm = zeros(length(wo),1);
    for j = 1: length(wo);
     SNR660 = P660(:,wo(j))./real(sqrt(P660(:,wo(j))));
     wok = find(SNR660(FitSelC) >0);
     rausch = mymean(SNR660(FitSelC(wok)));
     gut=find(H>400 & SNR660>10.*rausch);
     if length(gut) >5  % ok mind. 5 Posit. zum Vgl
         norm(j) = sum(mischvv(gut)) ./sum(wasservisprop(gut,wo(j)));
     end
    end   % for wo 
if all(norm) % alle Elemente ~=0
   var=max(norm)./min(norm);
   if var < 1.25  % alles das selbe, nun normieren
     wasserviskali = mymean(norm);
     disp( 'normieren vis')
   end % homog Daten
end % Normierung zu allen Zeiten
end % length(wo)

% wenn das nicht fruchtet suche Kalibrierwerte aus Ueberlappfile
if isnan(wasserUVkali) & ~isnan(uvwasserfaktor)
    wasserUVkali = uvwasserfaktor; 
end
if isnan(wasserviskali) & ~isnan(viswasserfaktor)
    wasserviskali = viswasserfaktor; 
end

% das finale Ergebnis ist nun:
wasserUV = wasserUVprop .* wasserUVkali;
wasservis = wasservisprop .* wasserviskali;





%%%Datenqualitaet testen

LvRUV=[0.7 1 1.05 1.1 1.2 1.3 1.4 1.5 1.75 2 2.25 2.5 3 4 5 10 99];
LvRvis=[0.7 1 1.07 1.1 1.15 1.2 1.3 1.5 1.75 2 2.5 3 5 10  100];
LvRIR=[0.7 1 1.2 1.3 1.5 1.75 2 2.5 3 4 5 7 10 20 30 60 99];
LvVD=[6 10 14 15 16 17 20 25 30 50 75 100 150 200 300 600 1000].*1e-1;
LvR2=[0.9 1 1.05 1.075 1.1 1.125 1.15 1.2 1.25 1.3 1.5 1.75 2 2.5 3 5 10];
LvR3=[0.7 1 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.5 1.6 1.7 2 2.5 3 4 5 7.5 10 25 50 99];

% Abspeichern vorbereiten

% BSR532flag = 1;
% BSR355flag = 1;
% BSR1064flag = 1;

BetaAer532flag = 1;
BetaAer355flag = 1;
BetaAer1064flag = 1;

VD532flag = 1; % nach Umbau fehlerhaft
VD355flag = 1; % nach Umbau fehlerhaft


%Dnsp=Density(:,wasserzeitpos); Tp=Temp(:,wasserzeitpos); %Dichte- und Tempprof 1-dim


%%%Abspeichern in drei Formaten:
% a: Minimales Ergebnisfile
% b: erweiterter Datensatz
% c: Strato Datensatz

%%%a: Minimales Ergebnisfile
% 
% STRmindata=['save ' OutputDir DatStr '_mindata.mat' ...
%     ' H DN532 DN355 DN1064 BetaAer532tot BetaAer355tot BetaAer1064' ...
%     ' BetaAer532toterr BetaAer355toterr BSR532toterr BSR355toterr' ...
%     ' BetaAer532flag BetaAer355flag BetaAer1064flag LR532 LR355 LR1064' ...
%     ' VD532 VD532err VD355 VD355err VD532flag VD355flag BSR532tot BSR355tot ' ...
%     ' BSR1064 BSR532flag BSR355flag BSR1064flag Density Temp' ...
%     ' AlphaAer355Raman AlphaAer532Raman FieldStop ZPosition' ... %aus Gesamtprofil
%     ' BetaAer355 BetaAer355S BetaAer532 BetaAer532S BSR355Klett BSR355SKlett ' ...
%     ' BSR532Klett BSR532SKlett BSR355quot BSR532quot P532A P532C P1064 P607 P387' ...
%     ' P355A P355C P355SA P355SC P532SA P532SC P660 P407 P607 P532 P532final' ...
%     ' P532Sfinal P532finalerr P532Sfinalerr P387final P607final P1064final P1064finalerr' ...
%     ' feu rfeis O3Density BSR532gesKlett BSR532quot BSR532urquot' ...
%     ' wasservisprop wasserUVprop wasserviskali wasserUVkali ballonzeit'];
% 






STRmindata=['save ' OutputDir DatStr '_mindata_iterativ.mat' ...
 ' H DN355 DN355S DN532 DN532S DN1064 DN387 DN407 DN607 DN660 FieldStop ZPosition' ...
 ' DATHV355s DATHV355p DATHV532s DATHV532p DATHV387 DATHV407 DATHV660 DATHV1064' ...
 ' BSRAtFit355 BSRAtFit532 BSRAtFit1064 BSRFit355Avektor BSRFit355SAvektor' ...
 ' BSRFit532Avektor BSRFit532SAvektor LR355arr LR355Sarr LR532arr LR532Sarr LR1064arr ueko' ...
 ' R355A R355C R355SA R355SC R532A R532C R532SA R532SC R607A R607C' ...
 ' R387A R387C R407A R407C R660C R1064 P355 P355S' ...
 ' P355A P355C P355final P355finalerr P355SA P355SC P355Sfinal P355Sfinalerr P387C' ...
 ' P387 P387A P387final P407 P407C P407A P532A P532C P532 P532final P532finalerr P532SA P532SC' ...
 ' P532S P532Sfinal P532Sfinalerr P607 P607C P607A P607final P660 P1064 P1064final' ...
 ' P1064finalerr FitRangeC FitRangeA Density Temp feu rfeis mischvv O3Density ballonzeit' ...
 ' BSR355tot BSR355toterr BSR355quot BSR355Klett BSR355SKlett BSR355PKlett'...
 ' BSR532tot BSR532toterr BSR532quot BSR532Klett BSR532SKlett BSR532PKlett' ...
 ' BetaAer355P BetaAer355Perr BetaAer355S BetaAer355Serr BetaAer355tot' ...
 ' BetaAer355toterr BetaAer532P BetaAer532Perr BetaAer532S BetaAer532Serr' ...
 ' BetaAer532tot BetaAer532toterr BSR1064 BSR1064err BetaAer1064 BetaAer1064err' ...
 ' VD355 VD355err AeroDep355 AeroDep355err VD532 VD532err AeroDep532 AeroDep532err' ...
 ' AlphaAer355Raman AlphaAer355Ramanerr AlphaAer532Raman AlphaAer532Ramanerr' ...
 ' LR355Raman LR355Ramanerr LR532Raman LR532Ramanerr FitRangeC FitRangeA' ...
 ' wasservisprop wasserUVprop wasserviskali wasserUVkali' ...
 ' BetaAer355flag BetaAer532flag  BetaAer1064flag VD355flag VD532flag' ...
 ' BSRAtFit355arr BSRAtFit355Sarr BSRAtFit532arr BSRAtFit532Sarr BSRAtFit1064arr' ...
 ' attenu355P attenu355S attenu532P attenu532S attenu1064 SNR355A SNR355C SNR355K SNR532A SNR532C SNR532K' ...
 ' SNR355SK SNR532SK SNR1064' ...
 ' hschritteUV Alpha355final Alpha387final AngarrUV Tau355final Tau387final' ...
 ' hschritteVIS Alpha532final Alpha607final Tau532final Tau607final AngarrVIS' ...
 ' Alpha355minfinal Alpha355maxfinal Alpha387minfinal Alpha387maxfinal Tau355minfinal Tau355maxfinal' ...
 ' Tau387minfinal Tau387maxfinal Alpha532minfinal Alpha532maxfinal Alpha607minfinal Alpha607maxfinal' ...
 ' Tau532minfinal Tau532maxfinal Tau607minfinal Tau607maxfinal'];...
  

    


eval(STRmindata);
% ss=['save /atm_meas/awipev/lidar/critter/kampagnen/nacht2019/' DatStr 'alldata']
% eval(ss)


dHround=round(dH);
dT=[];
if length(DN532)>2, dT=86400.*(DN532(3)-DN532(2)); end
if dHround==60 
fs=14;
figure(1); clf
[R532l, DNRvis] = maketimegaps(BSR532tot, DN532, 1/24/2);
mycontourf( DNRvis, H(1:400)./1e3, R532l(1:400,:), LvR2), hold on
title(['\bf KARL: Ny-Alesund, ' datestr((DN532(1)+DN532(end))/2,1)], 'fontsize',fs+1);       
ylabel(['\bf altitude [km]'],'fontsize',fs)
xlabel('\bf time (UT)','fontsize',fs)
datetick('x',15), hold on,
xcolorbar2017( 'vert',LvR2  ,2, [], 0, 'N',LvR2(2))
text( 1.15, 0.95, ['\bf BS-Ratio'],'unit','normaliz', ...    %%1.15, 0.958
   'horizontalalign','center','fontsize',fs-3);
text( 1.15, 0.90, ['\bf'  '532 nm' ],'unit','normaliz', ...            %%1.15, 0.92 
    'horizontalalign','center','fontsize',fs-3);
PlotDir='/atm_meas/awipev/lidar/critter/kampagnen/moon2020/weiterebilder/';
Str1=['print -djpeg -r100 ' PlotDir 'karl' DatStr '.png']; 
eval(Str1)
end


% STRalldata=['save ' OutputDir DatStr '_alldata_ueko2.mat' ];  eval(STRalldata);


% savedir=['/lidar4/critter/projekte/langzeitdaten/2009/'];
% savestr=['save ' savedir DatStr 'alldata'];
% eval(savestr)
%  


% 
% PlotDir = '/lidar4/critter/speziellemessungen/Jul2015/'
% 
% str= ['save ' PlotDir 'KARL30m' DatStr '_alldata.mat '];
% eval(str)
% BSRAtFit355

return

















%f?r Datenbanken Speichern unter Normalbedingungen
% und nur wenn genug Daten
if length(DN355) >=3,
dT=86400.*(DN355(3)-DN355(2));
if ((dH >58 && dH< 62) && ( dT> 500 && dT < 720) ),
    

% nicht g?ltige Daten herausschmei?en! 
tropo=find(H>600 & H<12000);
validtimes = [];
bis=size(P1064,2);
for j=1:bis,
   if ((DN532(j) == DN1064(j)) &&  (max(BSR532P(tropo,j))< 50)  && (min(BSR532P(tropo,j)> 1)) ...
           && (max(AeroDep532(tropo,j))< 1)  && (min(AeroDep532(tropo,j))> 0) )
       validtimes=concat(1,validtimes,j);
   end
end


% welche H?he? abschneiden? mehr als 30km nie
hwo=find(H > 500 & H<=14000);    lhwo=length(hwo);
%Nur eine Zeit
DN = DN532(validtimes);
H=H(hwo); vv=validtimes;
R355A=R355A(hwo,vv);  R355C=R355C(hwo,vv);  R355SA=R355SA(hwo,vv);  R355SC=R355SC(hwo,vv); 
R532A=R532A(hwo,vv);  R532C=R532C(hwo,vv);  R532SA=R532SA(hwo,vv);  R532SC=R532SC(hwo,vv);  
R1064=R1064(hwo,vv);  R387A=R387A(hwo,vv);  R387C =R387C(hwo,vv);   R607A =R607A(hwo,vv); 
R607C=R607C(hwo,vv);  R407A=R407A(hwo,vv);  R407C =R407C(hwo,vv);   R660C = R660C(hwo,vv);
P355A=P355A(hwo,vv);  P355C=P355C(hwo,vv);  P355SA=P355SA(hwo,vv);  P355SC=P355SC(hwo,vv); 
P532A=P532A(hwo,vv);  P532C=P532C(hwo,vv);  P532SA=P532SA(hwo,vv);  P532SC=P532SC(hwo,vv);  
P1064=P1064(hwo,vv);  P387A=P387A(hwo,vv);  P387C =P387C(hwo,vv);   P607A =P607A(hwo,vv); 
P607C=P607C(hwo,vv);  P407A=P407A(hwo,vv);  P407C =P407C(hwo,vv);   P660 =P660(hwo,vv);
P355=P355(hwo,vv); P355S=P355S(hwo,vv);  P532=P532(hwo,vv);   P532S=P532S(hwo,vv);
attenu355P= attenu355P(hwo,vv); attenu355S=attenu355S(hwo,vv); attenu532P=attenu532P(hwo,vv);
attenu532S=attenu532S(hwo,vv); attenu1064=attenu1064(hwo,vv);
BetaAer532P = BetaAer532P(hwo,vv);  BetaAer532S = BetaAer532S(hwo,vv);  BetaAer355P = BetaAer355P(hwo,vv);  
BetaAer355S = BetaAer355S(hwo,vv);  BetaAer1064 = BetaAer1064(hwo,vv); BetaAer355tot= BetaAer355tot(hwo,vv);
BetaAer532tot= BetaAer532tot(hwo,vv);
BetaAer355toterr= BetaAer355toterr(hwo,vv);  BetaAer532toterr= BetaAer532toterr(hwo,vv);  BetaAer1064err = BetaAer1064err(hwo,vv);
LR355arr = LR355arr(hwo,vv);  LR355Sarr = LR355Sarr(hwo,vv);  LR532arr = LR532arr(hwo,vv);  
LR532Sarr = LR532Sarr(hwo,vv); LR1064arr = LR1064arr(hwo,vv);  BSRAtFit355arr=BSRAtFit355arr(vv);
BSRAtFit355Sarr=BSRAtFit355Sarr(vv); BSRAtFit532arr=BSRAtFit532arr(vv);
BSRAtFit532Sarr=BSRAtFit532Sarr(vv); BSRAtFit1064arr=BSRAtFit1064arr(vv);
SNR355C= SNR355C(hwo,vv); SNR355A= SNR355A(hwo,vv); SNR355K= SNR355K(hwo,vv); SNR355SK= SNR355SK(hwo,vv);
SNR532C= SNR532C(hwo,vv); SNR532A= SNR532A(hwo,vv); SNR532K= SNR532K(hwo,vv); SNR532SK= SNR532SK(hwo,vv); SNR1064=SNR1064(hwo,vv);
AeroDep355= AeroDep355(hwo,vv); AeroDep532= AeroDep532(hwo,vv); VD355= VD355(hwo,vv); VD532= VD532(hwo,vv);
AeroDep355err= AeroDep355err(hwo,vv); AeroDep532err= AeroDep532err(hwo,vv); VD355err= VD355err(hwo,vv); VD532err= VD532err(hwo,vv);


% Auff?llen mit NaN wenn Fehler zu gross
for j=1:length(vv)
    schlecht355=find(SNR355K(:,j) <2 | BetaAer355tot(:,j) <=0);
    BetaAer355tot(schlecht355,j) = NaN;
    schlecht355=find(SNR355K(:,j) <2 | SNR355SK(:,j) <2 | AeroDep355(:,j) >1 | AeroDep355(:,j) < 0);
    AeroDep355(schlecht355,j) = NaN;
    schlecht532=find(SNR532K(:,j) <2 | BetaAer532tot(:,j) <=0);
    BetaAer532tot(schlecht532,j) = NaN;
    schlecht532=find(SNR532K(:,j) <2 | SNR532SK(:,j) <2 | AeroDep532(:,j) >1 | AeroDep532(:,j) < 0);
    AeroDep532(schlecht532,j) = NaN;
    schlecht1064=find(SNR1064(:,j) <2 | BetaAer1064(:,j) <=0);
    BetaAer1064(schlecht1064,j) = NaN;
end

liste = ['DN H R355A R355C R355SA R355SC R532A R532C R532SA R532SC R1064 R387A R387C R607A R607C R407A R407C R660C' ...
' P355A P355C P355SA P355SC P532A P532C P532SA P532SC P1064 P387A P387C P607A P607C P407A P407C P660' ... 
' attenu355P attenu355S attenu532P attenu532S attenu1064 BetaAer532tot BetaAer532toterr BetaAer532P BetaAer532Perr' ...
' BetaAer532S BetaAer532Serr BetaAer355P BetaAer355Perr BetaAer355S BetaAer355Serr BetaAer1064 BetaAer1064err' ...
' LR355arr LR355Sarr LR532arr LR532Sarr LR1064arr' ...
' SNR355C SNR355A SNR355K SNR532C SNR532A SNR532K SNR1064' ...
];
str = ['save /atm_meas/awipev/lidar/karl/matlab/datenbanken/' DatStr '.mat ' liste];
eval(str)

jahr=['20' DatStr(1:2)];

for j=1:length(vv)
dataascii = zeros(lhwo,19);                    
dataascii(:,1)=H;
dataascii(:,2)=P355(:,j); dataascii(:,3)=P355S(:,j); 
dataascii(:,4)=BetaAer355tot(:,j); dataascii(:,5)=BetaAer355toterr(:,j); dataascii(:,6)=LR355arr(:,j);
dataascii(:,7)=AeroDep355(:,j); dataascii(:,8)=AeroDep355err(:,j); 
dataascii(:,9)=P532(:,j); dataascii(:,10)=P532S(:,j); 
dataascii(:,11)=BetaAer532tot(:,j); dataascii(:,12)=BetaAer532toterr(:,j); dataascii(:,13)=LR532arr(:,j);
dataascii(:,14)=AeroDep532(:,j); dataascii(:,15)=AeroDep532err(:,j); 
dataascii(:,16)=P1064(:,j);
dataascii(:,17)=BetaAer1064(:,j); dataascii(:,18)=BetaAer1064err(:,j);  dataascii(:,19)=LR1064arr(:,j);

q=datestr(DN(j));
q1=findstr(q,' ');
if ~isempty(q1), q(q1)='-'; end
fname = ['KARL-ASCII-' q]; 

str = ['save /atm_meas/awipev/lidar/karl/matlab/datenbanken/' jahr '/' fname  ' dataascii' ' -ascii'];
eval(str)

end

outdir=[OutputDir 'hdf/'];
hdffilename=['karl' DatStr '.h5'];
%hdf5write(hdffilename, outdir, liste)

end % einlesen

end % mind. 3 Zeiten an Daten

return

% 1Photon:
for j=1:49, Q10=sort(diff(P407C(500:end,j))); wo10=find(diff(Q10)>1e-10); sQ10=min(diff(Q10(wo10))), end


finaldimen=size(P532);
 
ncfilename=['karl' DatStr '.nc'];
ncfullfilename=['/atm_meas/awipev/lidar/karl/matlab/datenbanken/' jahr '/' ncfilename];
ncid = netcdf.create(ncfullfilename ,'NC_WRITE');
heightsteps = netcdf.defDim(ncid,'rows',finaldimen(1));
timesteps = netcdf.defDim(ncid,'length',finaldimen(2));
stringlength = netcdf.defDim(ncid,'stringlength',20);
varidTime = netcdf.defVar(ncid,'Time','NC_CHAR',[timesteps stringlength]);
varidMatlabTime = netcdf.defVar(ncid,'MatlabTime','NC_DOUBLE',[1 1]);
varidH = netcdf.defVar(ncid,'Height','NC_DOUBLE',[heightsteps 1]);
varidP355 = netcdf.defVar(ncid,'P355','NC_DOUBLE',[heightsteps timesteps]);
varidP355S = netcdf.defVar(ncid,'P355S','NC_DOUBLE',[heightsteps timesteps]);
varidB355 = netcdf.defVar(ncid,'BetaAer355tot','NC_DOUBLE',[heightsteps timesteps]);
varidB355err = netcdf.defVar(ncid,'BetaAer355toterr','NC_DOUBLE',[heightsteps timesteps]);
varidLR355 = netcdf.defVar(ncid,'LR355','NC_DOUBLE',[heightsteps timesteps]);
varidAeroDep355 = netcdf.defVar(ncid,'AeroDep355','NC_DOUBLE',[heightsteps timesteps]);
varidAeroDep355err = netcdf.defVar(ncid,'AeroDep355err','NC_DOUBLE',[heightsteps timesteps]);
varidP532 = netcdf.defVar(ncid,'P532','NC_DOUBLE',[heightsteps timesteps]);
varidP532S = netcdf.defVar(ncid,'P532S','NC_DOUBLE',[heightsteps timesteps]);
varidB532 = netcdf.defVar(ncid,'BetaAer532tot','NC_DOUBLE',[heightsteps timesteps]);
varidB532err = netcdf.defVar(ncid,'BetaAer532toterr','NC_DOUBLE',[heightsteps timesteps]);
varidLR532 = netcdf.defVar(ncid,'LR532','NC_DOUBLE',[heightsteps timesteps]);
varidAeroDep532 = netcdf.defVar(ncid,'AeroDep532','NC_DOUBLE',[heightsteps timesteps]);
varidAeroDep532err = netcdf.defVar(ncid,'AeroDep532err','NC_DOUBLE',[heightsteps timesteps]);
varidP1064 = netcdf.defVar(ncid,'P1064','NC_DOUBLE',[heightsteps timesteps]);
varidB1064 = netcdf.defVar(ncid,'BetaAer1064','NC_DOUBLE',[heightsteps timesteps]);
varidB1064err = netcdf.defVar(ncid,'BetaAer1064err','NC_DOUBLE',[heightsteps timesteps]);
varidLR1064 = netcdf.defVar(ncid,'LR1064','NC_DOUBLE',[heightsteps timesteps]);

netcdf.putAtt(ncid,varidMatlabTime,'description','readable time');
netcdf.putAtt(ncid,varidMatlabTime,'units','dd-mmm-yyyy-hh:mm:sec');
netcdf.putAtt(ncid,varidMatlabTime,'description','Time in matlab format');
netcdf.putAtt(ncid,varidMatlabTime,'units','[days after 0.Jan.0000]');
netcdf.putAtt(ncid,varidH,'description','altitude above sea level');
netcdf.putAtt(ncid,varidH,'units','[m]');
netcdf.putAtt(ncid,varidP355,'description','lidar profile @ 355nm "parallel" polarized');
netcdf.putAtt(ncid,varidP355,'units','arbitrary');
netcdf.putAtt(ncid,varidP355S,'description','lidar profile @ 355nm "perpendicular" polarized');
netcdf.putAtt(ncid,varidP355S,'units','arbitrary');
netcdf.putAtt(ncid,varidB355,'description','aerosol backscatter coefficient @ 355nm');
netcdf.putAtt(ncid,varidB355,'units','[1/(m*sr)]');
netcdf.putAtt(ncid,varidB355err,'description','maximal error of aerosol backscatter coefficient @ 355nm');
netcdf.putAtt(ncid,varidB355err,'units','[1/(m*sr)]');
netcdf.putAtt(ncid,varidLR355,'description','lidar ratio @ 355nm');
netcdf.putAtt(ncid,varidLR355,'units','[sr]');
netcdf.putAtt(ncid,varidP532,'description','lidar profile @ 532nm "parallel" polarized');
netcdf.putAtt(ncid,varidP532,'units','arbitrary');
netcdf.putAtt(ncid,varidP532S,'description','lidar profile @ 532nm "perpendicular" polarized');
netcdf.putAtt(ncid,varidP532S,'units','arbitrary');
netcdf.putAtt(ncid,varidB532,'description','aerosol backscatter coefficient @ 532nm');
netcdf.putAtt(ncid,varidB532,'units','[1/(m*sr)]');
netcdf.putAtt(ncid,varidB532err,'description','maximal error of aerosol backscatter coefficient @ 532nm');
netcdf.putAtt(ncid,varidB532err,'units','[1/(m*sr)]');
netcdf.putAtt(ncid,varidLR532,'description','lidar ratio @ 532nm');
netcdf.putAtt(ncid,varidLR532,'units','[sr]');
netcdf.putAtt(ncid,varidP1064,'description','lidar profile @ 1064nm  unpolarized');
netcdf.putAtt(ncid,varidP1064,'units','arbitrary');
netcdf.putAtt(ncid,varidB1064,'description','aerosol backscatter coefficient @ 1064nm');
netcdf.putAtt(ncid,varidB1064,'units','[1/(m*sr)]');
netcdf.putAtt(ncid,varidB1064err,'description','maximal error of aerosol backscatter coefficient @ 1064nm');
netcdf.putAtt(ncid,varidB1064err,'units','[1/(m*sr)]');
netcdf.putAtt(ncid,varidLR1064,'description','lidar ratio @ 1064nm');
netcdf.putAtt(ncid,varidLR1064,'units','[sr]');

netcdf.endDef(ncid);


netcdf.putVar(ncid,varidTime,q);
netcdf.putVar(ncid,varidMatlabTime,DN);
netcdf.putVar(ncid,varidH,H);
netcdf.putVar(ncid,varidP355,P355);
netcdf.putVar(ncid,varidP355,P355S);
netcdf.putVar(ncid,varidB355,BetaAer355tot);
netcdf.putVar(ncid,varidB355err,BetaAer355toterr);
netcdf.putVar(ncid,varidLR355,LR355arr);
netcdf.putVar(ncid,varidAeroDep355,AeroDep355);
netcdf.putVar(ncid,varidAeroDep355err,AeroDep355err);
netcdf.putVar(ncid,varidP532,P532);
netcdf.putVar(ncid,varidP532,P532S);
netcdf.putVar(ncid,varidB532,BetaAer532tot);
netcdf.putVar(ncid,varidB532err,BetaAer532toterr);
netcdf.putVar(ncid,varidLR532,LR532arr);
netcdf.putVar(ncid,varidAeroDep532,AeroDep532);
netcdf.putVar(ncid,varidAeroDep532err,AeroDep532err);
netcdf.putVar(ncid,varidP1064,P1064);
netcdf.putVar(ncid,varidB1064,BetaAer1064);
netcdf.putVar(ncid,varidB1064err,BetaAer1064err);
netcdf.putVar(ncid,varidLR1064,LR1064arr);


netcdf.close(ncid);





