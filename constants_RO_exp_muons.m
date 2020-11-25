function [P,lambda,decayL]=constants_RO_exp_muons(nuclide)

switch nuclide
    case 10
        % production rates SLHL
        P.Pn_SLHL = 4.01;       % neutron spallation SLHL production rate (10Be = 4.01 for Lal1991/Stone2000 scaling - according to Phillips et al., 2016) (at/g/yr)
        P.Pms_SLHL = 0.012;     % Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        P.Pmf_SLHL = 0.039;     % Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        % attenuation lengths
        lambda.Ln  =  160;       % g/cm² 
        lambda.Lmf =  1500;      % 1500
        lambda.Lms =  4320;      % 4320 
        % decay constant
        decayL = log(2)/1.39e6;	 
    case 14
        % production rates SLHL
        P.Pn_SLHL = 12.24;		% neutron spallation SLHL production rate (26Al = 27.93 for Lal1991/Stone2000 scaling - according to Borchers et al., 2016) (at/g/yr)
        P.Pms_SLHL = 0.012;     % Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        P.Pmf_SLHL = 0.039;     % Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        % attenuation lengths
        lambda.Ln  =  160;       % g/cm² 
        lambda.Lmf =  1500;      %  
        lambda.Lms =  4320;      % 
        % decay constant
        decayL = log(2)/5730;	 
    case 26
        % production rates SLHL
        P.Pn_SLHL = 27.93;      % neutron spallation SLHL production rate (26Al = 27.93 for Lal1991/Stone2000 scaling - according to Borchers et al., 2016) (at/g/yr)
        P.Pms_SLHL = 0.012;     % Slow Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        P.Pmf_SLHL = 0.039;     % Fast Muon-induced Production Rate of 10Be (at/g/yr) scaled to sea level (Braucher et al. 2011) 
        % attenuation lengths
        lambda.Ln  =  160;      % g/cm² 
        lambda.Lmf =  1500;      % check these
        lambda.Lms =  4320;      % check
        % decay constant
        decayL = log(2)/7.17e5;	 
    case 36
        
end

return