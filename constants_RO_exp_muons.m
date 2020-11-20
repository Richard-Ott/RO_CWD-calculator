function [fs]=constants_RO_exp_muons(lat,altitude)
% this function only works for one sample at a time! If you have more
% samples, then run this function in a loop!
% Richard Ott, 2020

% Rock density in kg/m3
fs.density = 2600;


% Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)to account for diminished cosmic ray attenuation with decreased atmospheric pressue 
p=1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*altitude)))); 

%Half lives ---------------
H10=1.387e6; % Chemeleff2010
fs.L10=log(2)/H10;
H26=0.705e6; % Niishizumi 2004
fs.L26=log(2)/H26;
H14=5700;
fs.L14=log(2)/H14;
fs.L21=0;

%Attenuation lengths
fs.att_l_spal=1600; %[kg/m2], Equal to 160g/cm2
fs.att_l_fm=43200; %[kg/m2], Equal to 4320g/cm2
fs.att_l_nm=15000; %[kg/m2], Equal to 1500g/cm2

% p_scaling_spal=1*exp((1013.25-p)/242); % WHERE DOES 242 COME FROM? LAMBDA?

% apply Stone scaling
addpath 'C:\Users\r_ott\Dropbox\Richard\PhD_ETH\matlab\CRONUS cosmo calculation\cronusearth-2.0\common\production'
scaling_stone = stone2000(lat,p,1);     % only for spallation
% fs.P10_top_spal=p_scaling_spal*4.01e3* scaling_stone;  %SLHL Phillips 2016
% fs.P26_top_spal=p_scaling_spal*27.93e3* scaling_stone; %SLHL Phillips 2016
% fs.P14_top_spal=p_scaling_spal*12.24e3* scaling_stone; %SLHL Phillips 2016
fs.P10_top_spal=4.01e3* scaling_stone;  %SLHL Phillips 2016
fs.P26_top_spal=27.93e3* scaling_stone; %SLHL Phillips 2016
fs.P14_top_spal=12.24e3* scaling_stone; %SLHL Phillips 2016

%%
%Muons
mindepth = 0;    % g/cm2
maxdepth = 9800; % g/cm2     was at 7800
%maxdepth = 2600;
%Muons 10Be
f_star=0.00191;     % Model 1A, alpha=1;
Natoms = 2.006e22;  % Oxygen atoms pr gram Quartz
sigma0 = 0.280e-30; % model 1A, alpha=1;

p_muons = p_rate_calc2(f_star,Natoms,sigma0,p,mindepth,maxdepth);

fs.P10_att_l_fm=p_muons.L(1)*10; %Converted to kg m-2
fs.P10_att_l_nm=p_muons.L(2)*10; %Converted to kg m-2
fs.P10_top_fm=p_muons.P(1)*1000; %Converted to atoms kg-1 yr-1
fs.P10_top_nm=p_muons.P(2)*1000; %Converted to atoms kg-1 yr-1

%Muons 14C
f_star=0.137; %Model 1A, alpha=1;
Natoms = 2.006e22; %Oxygen atoms pr gram Quartz
sigma0 = 2.37e-30; % model 1A, alpha=1;

p_muons=p_rate_calc2(f_star,Natoms,sigma0,p,mindepth,maxdepth);

fs.P14_att_l_fm=p_muons.L(1)*10; %Converted to kg m-2
fs.P14_att_l_nm=p_muons.L(2)*10; %Converted to kg m-2
fs.P14_top_fm=p_muons.P(1)*1000; %Converted to atoms kg-1 yr-1
fs.P14_top_nm=p_muons.P(2)*1000; %Converted to atoms kg-1 yr-1

%Muons 26Al
f_star=0.0133; %Model 1A, alpha=1;
Natoms = 1.003e22; %Si atoms pr gram Quartz
sigma0 = 3.89e-30; % model 1A, alpha=1;

p_muons=p_rate_calc2(f_star,Natoms,sigma0,p,mindepth,maxdepth);

fs.P26_att_l_fm=p_muons.L(1)*10; %Converted to kg m-2
fs.P26_att_l_nm=p_muons.L(2)*10; %Converted to kg m-2
fs.P26_top_fm=p_muons.P(1)*1000; %Converted to atoms kg-1 yr-1
fs.P26_top_nm=p_muons.P(2)*1000; %Converted to atoms kg-1 yr-1
return