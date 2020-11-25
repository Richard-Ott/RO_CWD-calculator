% Calculator for catchment average denudation rates. Denundation rates are
% calculated pixel-by-pixel and then averaged. 
% Uncertaintties are propagated following a Monte Carlo approach and
% include uncertainties on production rate and scaling model.
% Ignores decay (for now)
% Calculates the attenuation length for spallation following CRONUS and
% Sato (2008) with a dependency on time and air-pressure (for standard
% samples I use a rigidity cutoff of 4, CHECK BEFORE PUBLICATION)
%
% Input: - x,y of sampling location in lat, lon
%        - nuclide concentration and uncertainty
% Output: denudation rate with uncertainty in mm/a
%
% Richard Ott, 2020
clc
clear
close all

addpath('C:\Users\richa\Dropbox\Richard\PhD_ETH\GIS\DEM\filled_dem')
addpath('C:\Users\richa\Dropbox\Richard\PhD_ETH\data\geochronology\cosmo')

%% USER CHOICE ------------------------------------------------------------
nuclide = 10;    % 10 - 10Be, 14 - 14C, 26 - 26Al, 36 - 36Cl  (for now only with 10Be) 
mode    = 2;     % 1 - provide x,y coordinate of sampling location and nuclide concentrations manually
                 % 2 - provide a table of input data nx4 (Lat,Lon,N,dN)
utmzone = 35;    % utm zone of DEM
toposh = 0;      % do you want to include a topographic shielding correction?
density = 2.6;   % density in g/cm³

% production rate uncertainties as fractions 
Pn_uncert = 0.025;
Pms_uncert = 0.5;
Pmf_uncert = 0.5;
Sc_uncert = 0.09;% uncertainty scaling scheme

nMC = 5e3;         % Monte Carlo iterations
                 
%% LOAD DATA --------------------------------------------------------------
DEM = GRIDobj('crete_clipped_utm.tif');            % provide DEM

% load sample data 
if mode == 1
    x = input('Please enter the longitude of the sampling location. ');
    y = input('Please enter the latitude of the sampling location. ');
    N = input('Please enter the measured nuclide concentration ');
    dN = input('Please enter the measured nuclide concentration uncertainty');
    sample = input('Please enter the the sample name' ,'s');
elseif mode == 2
    rows = 1:5;
    [data,~,raw] = xlsread('cosmo_samples.xlsx');  % provide sample datas
    sample = raw(rows+2,1);
%     x = data(:,1);
%     y = data(:,2);
%     N = data(:,16)*1e4;
%     dN = data(:,17)*1e4;
    x = data(rows,1);
    y = data(rows,2);
    N = data(rows,16)*1e4;
    dN = data(rows,17)*1e4;
else
    error('You were supposed to set mode to either 1 or 2')
end

%% GET DRAINAGE BASIN -----------------------------------------------------
% compute flow of DEM
FD = FLOWobj(DEM,'preprocess','carve');
A = flowacc(FD);
S = STREAMobj(FD,'minarea',5e6,'unit','map');

% snap sampling location 2 stream
[x_utm,y_utm] = ll2utm(x,y,utmzone);
coords = [x_utm,y_utm]; 
[IX,dist] = snap2stream(A > 1e3,coords);
dist/1e3                                          % check distances in km to find where snapping might have gone wrong

% get drainage basins
DB = drainagebasins(FD,IX);

%% LOAD COSMO PARAMETERS --------------------------------------------------

% Scaling Equation Constants by Latitude (Stone,2000)
Lat = [0,10,20,30,40,50,60,90];
a = [31.8518,34.3699,40.3153,42.0983,56.7733,69.0720,71.8733,71.8733];
b = [250.3193,258.4759,308.9894,512.6857,649.1343,832.4566,863.1927,863.1927];
c = [-0.083393,-0.089807,-0.106248,-0.120551,-0.160859,-0.199252,-0.207069,-0.207069];
d = [7.4260e-5,7.9457e-5,9.4508e-5, 1.1752e-4,1.5463e-4,1.9391e-4,2.0127e-4,2.0127e-4];
e = [-2.2397e-8,-2.3697e-8,-2.8234e-8,-3.8809e-8,-5.0330e-8,-6.3653e-8,-6.6043e-8,-6.6043e-8];
m = [0.587,0.600,0.678,0.833,0.933,1.000,1.000,1.000];

[P,lambda,decayL,rigiditycutoff]=constants_RO_exp_muons(nuclide);    % get constants for production and attenuation, and radioactive decay

%% CALCULATE EROSION RATES ------------------------------------------------
% -------------------------------------------------------------------------
nBasins = length(IX);

% assign output variables
E          = zeros(nBasins,1);        % mean simulated rates
E_Std      = zeros(nBasins,1);        % uncertainty standard deviation
E_Med      = zeros(nBasins,1);        % median simluated rates
E_UpQuant  = zeros(nBasins,1);        % the quantiles are set to mimik a stdanard deviation (0.83 & 0.17)
E_LowQuant = zeros(nBasins,1);

for i = 1:nBasins
    % extract basins DEM
    Basin = DEM*nan;
    Basin.Z(DB.Z == i) = DEM.Z(DB.Z == i);
    Basin = crop(Basin);
    
    % topographic shielding
    if toposh
        Tshd = toposhielding(Basin,10,5);
        Tshd_mean = nanmean(Tshd.Z(Tshd.Z~=1));
    end
    %% CALCULATE PRODUCTION RATE FOR EVERY PIXEL AND AVERAGE
    
    % calculate Stone scaling parameters for every pixel
    [A,B,C,D,Es,M] = calc_Stone_ABDCE_catchment(Basin,utmzone,Lat,a,b,c,d,e,m);

    % Perform air pressure correction for elevation (Eq. 1 from Stone, 2000)to account for diminished cosmic ray attenuation with decreased atmospheric pressue 
    pres=1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*Basin.Z)))); 
    
    % calculate mean production
    [Pn,Pmf,Pms,Leff] = mean_production_catchment(Basin,0,P,lambda,A,B,C,D,Es,pres,rigiditycutoff);
    
    %% Calculate MC erosion rates
    
    % pre-sample random parameters for speed
    % generate random concentration with truncated normal distribution
    C_rand = truncnormrnd([nMC,1],N(i),dN(i),0,inf);

    % Uncertainty of neutrons: on the prod rate: 100*(0.1/3.9)=2.5% + 9% the scaling scheme
    PnStd = sqrt((Pn * Pn_uncert).^2 + (Pn* Sc_uncert).^2);
    Pn_rand = truncnormrnd([nMC,1],Pn,PnStd,0,inf);    % random n prod rate
    % Uncertainty of muon: on the prod rate: 50% + 9% the scaling scheme
    PmsStd = sqrt((Pms* Pms_uncert).^2 + (Pms* Sc_uncert).^2);
    Pms_rand = truncnormrnd([nMC,1],Pms,PmsStd,0,inf); % random ms prod rate
    PmfStd=sqrt((Pmf*Pmf_uncert).^2 + (Pmf* Sc_uncert).^2);
    Pmf_rand=truncnormrnd(nMC,Pmf,PmfStd,0,inf);       % random mf prod rate
    
    samples = nan(nMC,1);
    muN = density/Leff; muMs = density/lambda.Lms; muMf = density/lambda.Lmf;
    for j = 1:nMC
        if C_rand(j)~=0
            samples(j) = (1/C_rand(j))*((Pn_rand(j)/muN)+(Pms_rand(j)/muMs)+(Pmf_rand(j)/muMf));  % E in cm/a
        end    
    end
    samples = samples *10;   % convert to mm/a
    
    % calculate rate and uncertainties from MC
    E(i)          = mean(samples);        % mean simulated rates
    E_Std(i)      = std(samples);         % uncertainty standard deviation
    E_Med(i)      = median(samples);      % median simluated rates
    E_UpQuant(i)  = quantile(samples,0.83);        
    E_LowQuant(i) = quantile(samples,0.17);    
    
end

% rough estimate of integration time-scale only reliant on spallation in yrs
TimeScale=round(1./(muN*(E/10)));

vars = {'Sample','N','dN','E','dE','E_QuantUp','E_QuantLow','TimeScale'};
out_table = table(sample,N,dN,E,E_Std,E_UpQuant,E_LowQuant,TimeScale,'VariableNames',vars);

%% EXPORT 
export = input('Do you want to export the data? ');
if export
    tag = input('What filename? ','s');
    writetable(out_table,[tag '.xlsx'])
end




