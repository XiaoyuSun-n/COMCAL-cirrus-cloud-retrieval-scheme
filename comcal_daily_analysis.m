function comcal_daily_analysis( DTVon)

% Auswerteprogramm basierend auf nya_astar_tropo_all_n, ruft nacheinander verschiedene Funktionen aus \nya_astar auf
%
% DTVon: '220309'   Datum
 dh= 60;    % Hoehenaufloesung
 dt = 600;  % Zeitaufloesung

%glaetten = 5;

jahr =(DTVon(1:2));
numjahr = str2num(jahr);

status = 0;
if numjahr >=9, status=1; end % Fall 1 Teleskop neues KARL


if nargin < 3,
   dt=60; dh=60; %Standardeinstellungen
end


% Ueberlapp?
ueko=[];
%ueko = determineoverlapfile(DTVon);


%ComcalLidarDataDir = ['/lidar4/lidar/comcal/raw/2018/'];		% '/lidar4/lidar/comcal/tests/'
ComcalLidarDataDir = ['D:/Lidar/palau/ComCAL/20' DTVon(1:2) '_palau/' DTVon(1:4) '00/'];	

ComcalLidarDataDirZIP =  ['D:/Lidar/palau/ComCAL/20' DTVon(1:2) '_palau/' DTVon(1:4) '00/'];
% RohDaten gezippt


% wo liegen denn wohl die tollen Sondendaten?

OzoneInDir  = 'D:/Lidar/palau/Soundings/NASAAmes/';     %'/bsrn1/bsrnuser/2013/sonde/w_b/'; 
PtuInDir  =  'D:/Lidar/palau/Soundings/NASAAmes/';
  
AlternativSondenDir = '/bsrn1/bsrnuser/2017/sonde/w_b/';


DeZipDir = ['D:/Lidar/palau/ComCAL/20' DTVon(1:2) '_palau/entpackt/' DTVon];


  ZipMode = 1; % 1 = Entzippen 
  del = 0; % 1 = Rohdaten wieder aus DeZipDir loeschen   0 = nicht loeschen
  
  mtrdir = 'D:/Lidar/palau/ComCAL/matlab/';
  
  RAWDir =      [ mtrdir 'raw/']; % 
  DATDir =      [ mtrdir 'dat/'];
  DNSDir =      [ mtrdir 'dns/'];
  OZODir =      [ mtrdir 'ozo/'];
  AERDir =      [ mtrdir 'aer/'];
  PTUDir =      [ mtrdir 'ptu/'];
  
  DTBis = DTVon;
  DTptuVon=DTVon([1 2 3 4]);
  DTptuBis=DTBis([1 2 3 4]);
  
  AddTimeLicel2RAW=dt; %% passiert in licel2mat
  HeightRes = dh; %% passiert in licel2mat !!!


 %---------------------------------------------------------------------------------------------------------
 %---------------------------------------------------------------------------------------------------------
  
  
   comcal_ozo (DTptuVon, DTptuBis, OzoneInDir, OZODir);  % Ozone
  % comcal_ptu (DTptuVon, DTptuBis, PtuInDir, PTUDir);    % PTU      gibt
  % es nicht
  

  % AUSWERTEVERLAUF wird hier entschieden !!!
  % bitte entsprechend auskommentieren  
  
 
 HAVGRAWDir = [RAWDir 'havg' sprintf( '%04.0f', dh) '/'];
 HAVGDATDir = [DATDir 'havg' sprintf( '%04.0f', dh) '/'];
 HAVGAERDir = [AERDir 'havg' sprintf( '%04.0f', dh) '/'];
 HAVGDNSDir = [DNSDir 'havg' sprintf( '%04.0f', dh) '/'];
 
 AVGTime=dt; %passiert in AVG
 HTAVGRAWDir = [RAWDir 'havg' sprintf( '%04.0f',dh) '/tavg' sprintf( '%05.0f', AVGTime) '/'];
 HTAVGDATDir = [DATDir 'havg' sprintf( '%04.0f',dh) '/tavg' sprintf( '%05.0f', AVGTime) '/'];
 HTAVGAERDir = [AERDir 'havg' sprintf( '%04.0f',dh) '/tavg' sprintf( '%05.0f', AVGTime) '/'];
 HTAVGDNSDir = [DNSDir 'havg' sprintf( '%04.0f',dh) '/tavg' sprintf( '%05.0f', AVGTime) '/'];
    

%  if dt ~=600, %Abweichung vom Standard-Fall
%      HAVGRAWDir = HTAVGRAWDir;
%      HAVGDATDir = HTAVGDATDir;
%      HAVGAERDir = HTAVGAERDir;
%      HAVGDNSDir = HTAVGDNSDir;
%  end
%      




if ZipMode
  comcal_zip( DTVon, DTBis, ComcalLidarDataDirZIP, DeZipDir, HAVGRAWDir, AddTimeLicel2RAW, HeightRes, del, status); 
else
  comcal_raw( DTVon, DTBis, ComcalLidarDataDir, HAVGRAWDir, AddTimeLicel2RAW,HeightRes, status); % LICEl --> RAW
end
%  
%  %%nya_astar_tropo_avg( DTVon, DTBis, HAVGRAWDir, HTAVGRAWDir, AVGTime);
% 
 comcal_dat( DTVon, DTBis, HAVGRAWDir, HAVGDATDir, OZODir, PTUDir, HAVGDNSDir);  % RAW -> DAT/DNS


%%% comcal_aer( DTVon, DTBis, HAVGDATDir, HAVGDNSDir, OZODir, HAVGAERDir, ueko, 5);  % DAT --> AER
% 
str =['dattoaer_comcal_iterativ' DTVon];
list=['HAVGDATDir, HAVGDNSDir, OZODir, HAVGAERDir, ueko, 5'];
%str2=num2str(glaetten)
% 
str3 =[str '(' list ')']
% 

eval(str3)
 %dattoaer_comcal_iterativ220806( HAVGDATDir, HAVGDNSDir, OZODir, HAVGAERDir, ueko, glaetten);


 % 
 %comcal_plot(DTVon, dh)





