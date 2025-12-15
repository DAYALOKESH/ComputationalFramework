function PL = hata_wrapper(Scenario, TerrainProfile)
%HATA_WRAPPER Wrapper for the Empirical Hata/COST231 Model.
%
%   PL = hata_wrapper(Scenario, TerrainProfile)
%
%   Input:
%       Scenario - Struct with Frequency_MHz, Tx_Height_m, Rx_Height_m, Environment.
%       TerrainProfile - Struct with Distance_m (Elevation ignored by Hata).
%
%   Output:
%       PL - Vector of Total Path Loss (dB).

    % Prepare input for the existing hata_model function
    % hata_model expects Distance_Vector_km in the Scenario struct
    
    Scenario.Distance_Vector_km = TerrainProfile.Distance_m(:) / 1000.0;
    
    % Call the underlying model
    PL = hata_model(Scenario);
    
end
