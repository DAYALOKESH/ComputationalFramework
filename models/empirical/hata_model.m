function PL = hata_model(Scenario)
%HATA_MODEL Calculates Path Loss using Okumura-Hata or COST 231 Hata models.
%
%   PL = hata_model(Scenario)
%
%   Input:
%       Scenario - Struct with the following fields:
%           Frequency_MHz      : Carrier frequency in MHz
%           Tx_Height_m        : Transmitter antenna height in meters
%           Rx_Height_m        : Receiver antenna height in meters
%           Distance_Vector_km : Vector of distances in km
%           Environment        : String ('urban', 'suburban', 'rural', 'urban_large')
%
%   Output:
%       PL - Vector of Path Loss values in dB corresponding to Distance_Vector_km.
%
%   Logic:
%       - For f <= 1500 MHz, uses standard Okumura-Hata model.
%       - For f > 1500 MHz, uses COST 231 Hata model.
%
%   References:
%       [1] M. Hata, "Empirical Formula for Propagation Loss in Land Mobile Radio Services," 
%           IEEE Trans. Veh. Technol., vol. 29, no. 3, pp. 317-325, Aug. 1980.
%       [2] COST Action 231, "Digital mobile radio towards future generation systems," 
%           European Commission, Brussels, Tech. Rep., 1999.

    %% 1. Unpack and Validate Inputs
    f = Scenario.Frequency_MHz;
    h_tx = Scenario.Tx_Height_m;
    h_rx = Scenario.Rx_Height_m;
    d = Scenario.Distance_Vector_km;
    env = lower(string(Scenario.Environment));

    % Basic Validation
    if any(d <= 0)
        % warning('Distances <= 0 detected. Clamping to 0.001 km (1m) to avoid singularity.');
        d(d <= 0) = 0.001; 
    end
    
    if h_rx > h_tx
        warning('Receiver is higher than Transmitter. Hata model assumes h_tx > h_rx.');
    end

    %% 2. Model Selection
    % COST 231 extension is typically for 1500-2000 MHz. 
    % Standard Hata is 150-1500 MHz.
    
    use_cost231 = (f > 1500);
    
    %% 3. Pre-calculation of Terms
    log_f = log10(f);
    log_ht = log10(h_tx);
    log_d = log10(d);

    % Mobile Station Antenna Height Correction Factor a(h_rx)
    % This depends on environment (Urban Large vs others) and frequency.
    
    a_hr = zeros(size(d)); % usually scalar, but could be vector if h_rx was vector (not here)
    
    if env == "urban_large"
        if f <= 200 % Note: Hata formula split at 200/300/400 MHz varies by source. 
                    % Standard Hata: Large City correction varies by f
            % f <= 200 MHz
             a_hr = 8.29 * (log10(1.54 * h_rx))^2 - 1.1;
        else
            % f >= 400 MHz (Interpolate or use this for >200)
            % COST 231 typically uses this form for urban
             a_hr = 3.2 * (log10(11.75 * h_rx))^2 - 4.97;
        end
    else
        % For Small/Medium City, Suburban, Rural
        a_hr = (1.1 * log_f - 0.7) * h_rx - (1.56 * log_f - 0.8);
    end

    %% 4. Path Loss Calculation
    
    if use_cost231
        % --- COST 231 Hata Model (1.5 - 2.0 GHz) ---
        % L = 46.3 + 33.9 log(f) - 13.82 log(ht) - a(hr) + (44.9 - 6.55 log(ht)) log(d) + Cm
        
        A = 46.3;
        B = 33.9;
        C = 13.82;
        D = 44.9;
        E = 6.55;
        
        if env == "urban" || env == "urban_large"
            C_m = 3;
        else
            C_m = 0; % Suburban/Rural start from basic loss
        end
        
        L_base = A + B * log_f - C * log_ht - a_hr + (D - E * log_ht) * log_d + C_m;
        
        % Apply Suburban/Rural corrections relative to COST 231 Urban?
        % COST 231 usually defined for Urban. 
        % Commonly, the old Hata corrections are applied to the COST 231 base for other environments.
        
        if env == "suburban"
             % Standard Hata Suburban correction
             % L_sub = L_urban - 2(log(f/28))^2 - 5.4
             L_base = L_base - 2 * (log10(f / 28))^2 - 5.4;
             
        elseif env == "rural"
             % Standard Hata Rural correction
             % L_rural = L_urban - 4.78(log(f))^2 + 18.33 log(f) - 40.94
             L_base = L_base - 4.78 * (log_f)^2 + 18.33 * log_f - 40.94;
        end
        
        PL = L_base;
        
    else
        % --- Standard Okumura-Hata Model (150 - 1500 MHz) ---
        % L = 69.55 + 26.16 log(f) - 13.82 log(ht) - a(hr) + (44.9 - 6.55 log(ht)) log(d)
        
        A = 69.55;
        B = 26.16;
        C = 13.82;
        D = 44.9;
        E = 6.55;
        
        L_urban = A + B * log_f - C * log_ht - a_hr + (D - E * log_ht) * log_d;
        
        if env == "urban" || env == "urban_large"
            PL = L_urban;
            
        elseif env == "suburban"
            % L_sub = L_urban - 2(log(f/28))^2 - 5.4
            PL = L_urban - 2 * (log10(f / 28))^2 - 5.4;
            
        elseif env == "rural"
            % L_rural = L_urban - 4.78(log(f))^2 + 18.33 log(f) - 40.94
            PL = L_urban - 4.78 * (log_f)^2 + 18.33 * log_f - 40.94;
            
        else
            error('Unknown environment: %s. Use urban, urban_large, suburban, or rural.', env);
        end
    end
    
    % Ensure column vector output
    PL = PL(:);
    
end
