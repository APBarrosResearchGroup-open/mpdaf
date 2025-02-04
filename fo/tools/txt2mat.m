%----------------------------------------------------------------
% Interface to run MEMLS forward opeartor using Obseration ascii output
% P. Shrestha
% Log: Oct 8 2021
%----------------------------------------------------------------
clc;clear;

%----------------------------------------------------------------
%% User Settings
%----------------------------------------------------------------

% MEMLS parameters
m=0.1; q=0.2; sh=0.085; sv=0.06; ssh=0.064; ssv=0.045;
iangle = 35.0 ; % Incidence Angle
f=1;            % 1-X band, 2-Ku band   

fil  ='__fil__' ;
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

%%
%newstr   = split(fil,'_');
%if (size(newstr,1)>3)
%  newstr1 = split(newstr(4,1),'.')
%  foutname = strcat('memls_',band,'_',string(newstr(2,1)),'_',string(newstr(3,1)),'_',string(newstr1(1,1)),'.dat');
%  varname  = strcat(string(newstr(2,1)),'_',string(newstr(3,1)),'_',string(newstr1(1,1)));
%else
%  newstr1  = split(newstr(3,1),'.')
%  foutname = strcat('memls_',band,'_',string(newstr(2,1)),'_',string(newstr1(1,1)),'.dat');
%  varname  = strcat(string(newstr(2,1)),'_',string(newstr1(1,1)));
%end
newstr = split(fil,'/')
ndim = size(newstr)
newstr1 = split(newstr(ndim(1),1),'.');

foutname = strcat('memls_',band,string(newstr1(1,1)),'.dat');
varname  = strcat('Band1_',band,string(newstr1(1,1))); 
%
obsdata = load(fil);
save(varname,'obsdata','-v7.3');
clearvars obsdata;

%----------------------------------------------------------------
%% STEP 2  Data subset loop, call MEMLS Interface
%----------------------------------------------------------------

% Load data

load(strcat(varname,'.mat'));
nl = size(obsdata,1);

Tsnow1   = (reshape(obsdata(:,3),1,nl))';
LW1      = zeros(nl,1);
%Depth1   = zeros(nl,1) + 0.1; % 10 cm thick layers
Depth1   = (reshape(obsdata(:,1),1,nl))';
Density1 = (reshape(obsdata(:,2),1,nl))';
Dsnow1   = (reshape(obsdata(:,4),1,nl))';
Tsoil1   = (obsdata(1,5))';

% Check number of snow layers
 num1 = find(Depth1>0); N1 = length(num1);
 [~, sigma, ~] = amemls(fre, iangle, sh, sv, ssh, ssv, Tsnow1, LW1, Depth1,Density1, Dsnow1, 0, Tsoil1, 11, m, q);
 
% disp(sigma) 
disp(foutname)
save(foutname,"sigma","-ascii")
%delete *.mat

%----------------------------------------------------------------

