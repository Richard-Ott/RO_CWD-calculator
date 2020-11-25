function [Pn,Pmf,Pms,Leff] = mean_production_catchment(DEM,Tshd,P,lambda,A,B,C,D,E,pressure,rigiditycutoff)
% Calculates production rates for catchments using expoential
% approximations for muons. If you do not want shielding then Tshd = 0
% Input:  all rasters as GRIDobj
%           - DEM
%           - Tshd (GRIDobj)
%           - lambda, structure with fields Ln, Lms, Lmf
%           - A,B,C,D,E (Gribdobjs)
%           - pres - air pressure 
% Richard Ott, 2020

if ~Tshd
    Tshd = DEM/DEM;
end

% calculate effective attenuation length for neutrons following CRONUS and
% Sato (2008). Using rigidity cutoff of 4 for speed and simplicity.
% Check before publication!

pressure_vector = pressure(:);
Leff = nan(length(pressure_vector),1);
for i = 1:length(pressure_vector)
    if isnan(pressure_vector(i))
        continue
    else
        Leff(i) = rawattenuationlength(pressure_vector(i),single(rigiditycutoff));
    end
end
Leff = nanmean(Leff);    

%% Calculation production from spallation and muons

% Production from Spallation, assuming neutron attenuation length in air of 150 g/cm2
Pn_grid = Tshd.Z.*P.Pn_SLHL.*(A.Z+B.Z.*exp(-pressure/Leff)+C.Z.*pressure+D.Z.*pressure.^2+E.Z.*pressure.^3); % Stone, 2000

% Muons caluclated using exponential approximation which should be fine
% for moderate to fast erosion rates (see Balco, 2017)

% Production from Slow Muons, assuming (1) sea level pressure of 1013.25 mbar 
Pms_grid = Tshd.Z.*P.Pms_SLHL.*exp((1013.25-pressure)/lambda.Lms); 

% Production from Fast Muons, assuming (1) sea level pressure of 1013.25 mbar
Pmf_grid = Tshd.Z.*P.Pmf_SLHL.*exp((1013.25-pressure)/lambda.Lmf); 


% take means
Pn = nanmean(Pn_grid(:));
Pms = nanmean(Pms_grid(:));
Pmf = nanmean(Pmf_grid(:));


end

