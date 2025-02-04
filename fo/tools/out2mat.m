%----------------------------------------------------------------
% Interface to run MEMLS forward opeartor using MSHM ascii output
% P. Shrestha
% Log: Oct 5 2021
%----------------------------------------------------------------
clc;clear;

%----------------------------------------------------------------
%% User Settings
%----------------------------------------------------------------
% Matrix Subset for each MSHM runs
np=1; %1;
nl=30;
nt=64*1150; %48; 17035; %48; %17035 ; %48;
%
it0 = 5750;  %336; %0
it1 = 6900; % ;12624; %12624; %0

% MEMLS parameters
m=0.1; q=0.2; sh=0; sv=0; ssh=0; ssv=0;
iangle = 35.0 ; % Incidence Angle
f=1;            % 1-X band, 2-Ku band   

diri ='__diri__/' ; %/home/ps98/scratch/mpdaf_CCI_run/out/ens__X__/' ;
addpath('/u/ps98/mpdaf/fo/MEMLS'); % Add MEMLS to path

%---------------------------------------------------------------
%% STEP 1: BIG TIME LOOP Convert MSHM out to mat files
%----------------------------------------------------------------

% MEMLS setup
if f==1
  band='X'; fre=9.6; UdiH=1; UdiV=1; diH=0.1; diV=0.1;
else
  band='Ku'; fre=17.2; UdiH=0; UdiV=0; diH=0.1; diV=0.1;
end

%CPSGM for b=it0:nt:it1      %BIG TIME LOOP
  for b=it0:1150:it1      %BIG TIME LOOP
  % CPSGM
  if b == 6900
   nt = 64*1149
  end

  mtime = sprintf('%05d',b); 
  %%
  zongSD = load(strcat(diri,'zongSD_',mtime,'.out'));
  save(strcat('zongSD_',mtime),'zongSD','-v7.3');
  clearvars zongSD;

  %%
  zongSWE = load(strcat(diri,'zongSWE_',mtime,'.out'));
  save(strcat('zongSWE_',mtime),'zongSWE','-v7.3');
  clearvars zongSWE;

  %%
  Tsoil = load(strcat(diri,'Tsoil_',mtime,'.out'));
  Tsoil=reshape(Tsoil,np,nt);
  save(strcat('Tsoil_',mtime),'Tsoil','-v7.3');
  clearvars Tsoil

  %%%§%
  Depth = load(strcat(diri,'Depth_',mtime,'.out'));
  Depth=reshape(Depth,np,nt,nl);
  save(strcat('Depth_',mtime),'Depth','-v7.3');
  clearvars Depth;

  %%
  Density = load(strcat(diri,'Density_',mtime,'.out'));
  Density=reshape(Density,np,nt,nl);
  save(strcat('Density_',mtime),'Density','-v7.3');
  clearvars Density;

  %%
  LW = load(strcat(diri,'LW_',mtime,'.out'));
  LW=reshape(LW,np,nt,nl);
  save(strcat('LW_',mtime),'LW','-v7.3');
  clearvars LW;

  %%
  Dsnow = load(strcat(diri,'Dsnow_',mtime,'.out'));
  Dsnow=reshape(Dsnow,np,nt,nl);
  save(strcat('Dsnow_',mtime),'Dsnow','-v7.3'); 
  clearvars Dsnow;

  %%
  Tsnow = load(strcat(diri,'Tsnow_',mtime,'.out'));
  Tsnow=reshape(Tsnow,np,nt,nl);
  save(strcat('Tsnow_',mtime),'Tsnow','-v7.3');
  clearvars Tsnow;

%----------------------------------------------------------------
%% STEP 2  Data subset loop, call MEMLS Interface
%----------------------------------------------------------------

% Load data
  foutname = strcat('memls_',band,'_',mtime,'.dat')

  load(strcat('Tsnow_',mtime,'.mat'));
  load(strcat('LW_',mtime,'.mat'));
  load(strcat('Depth_',mtime,'.mat'));
  load(strcat('Density_',mtime,'.mat'));
  load(strcat('Dsnow_',mtime,'.mat'));
  load(strcat('Tsoil_',mtime,'.mat'));

  %Initialze backscatter
  nt1 = 64*1150    % CPSGM
  sigma=NaN(nt1,3);

  for it =1:nt
  for p =1:np

  Tsnow1   = (reshape(Tsnow(p,it,:),1,nl))';
  LW1      = (reshape(LW(p,it,:),1,nl))';
  Depth1   = (reshape(Depth(p,it,:),1,nl))';
  Density1 = (reshape(Density(p,it,:),1,nl))';
  Dsnow1   = (reshape(Dsnow(p,it,:),1,nl))';
  Tsoil1   = Tsoil(p,it)';

  % Check number of snow layers
  num1 = find(Depth1>0); N1 = length(num1);
  if N1 == 0
    % if no snow layer, ground backscatter?
    sigma(it,:) = [-0.1 -0.1 -0.1];
  else
%  [~, sigma, ~] = amemls(fre, iangle, 0, 0, 0, 0, 273., 0.1, 0.2, 300., 0.1, 0., 280., 11, m, q);
    [~, sigma(it,:), ~] = amemls(fre, iangle, 0, 0, 0, 0, Tsnow1, LW1, Depth1,Density1, Dsnow1, 0, Tsoil1, 11, m, q);
  end

  end
  end
 
  % disp(sigma) 
  disp(foutname)
  save(foutname,"sigma","-ascii")
  delete *.mat

end  % BIG TIME LOOP
%----------------------------------------------------------------

