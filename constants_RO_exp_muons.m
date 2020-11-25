function [P,lambda,decayL,rigiditycutoff]=constants_RO_exp_muons(nuclide)

switch nuclide
    case 10
        % production rates SLHL
        P.Pn_SLHL = 4.01;       % neutron spallation SLHL production rate (10Be = 4.01 for Lal1991/Stone2000 scaling - according to Phillips et al., 2016) (at/g/yr)
        P.Pms_SLHL = 0.012;     % Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        P.Pmf_SLHL = 0.039;     % Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        % attenuation lengths
        rigiditycutoff =4;      % for calculation of effective attenuation length. I think 4 is fine for most applicatons but check before publicaiton
%         lambda.Ln  =  160;       % g/cm², in the mean_production_catchment file the Stone eqn for Pn uses 150. Not sure if that is attenuation length and which is best to use. Difference is quite dramatic
        lambda.Lmf =  1500;      % apparently Braucher 2011 but cant find it in there
        lambda.Lms =  4320;      % apparently Braucher 2011 but cant find it in there
        % decay constant
        decayL = log(2)/1.39e6;	 
    case 14
        % production rates SLHL
        P.Pn_SLHL = 12.24;		% neutron spallation SLHL production rate (26Al = 27.93 for Lal1991/Stone2000 scaling - according to Borchers et al., 2016) (at/g/yr)
        %P.Pms_SLHL = ;     
        %P.Pmf_SLHL = ;     
        % attenuation lengths
        rigiditycutoff =4;      % for calculation of effective attenuation length. I think 4 is fine for most applicatons but check before publicaiton
%         lambda.Ln  =  160;       % g/cm² 
        lambda.Lmf =  1500;      %  
        lambda.Lms =  4320;      % 
        % decay constant
        decayL = log(2)/5730;	 
    case 26
        % production rates SLHL
        P.Pn_SLHL = 27.93;      % neutron spallation SLHL production rate (26Al = 27.93 for Lal1991/Stone2000 scaling - according to Borchers et al., 2016) (at/g/yr)
        P.Pms_SLHL = 0.84;      % Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        P.Pmf_SLHL = 0.081;     % Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        % attenuation lengths
        rigiditycutoff =4;      % for calculation of effective attenuation length. I think 4 is fine for most applicatons but check before publicaiton
%         lambda.Ln  =  160;      % g/cm² 
        lambda.Lmf =  1500;      % check these
        lambda.Lms =  4320;      % check
        % decay constant
        decayL = log(2)/7.17e5;	 
    case 36
        
end

return