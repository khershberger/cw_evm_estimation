%function [ evm ] = Untitled( pavs, vref, vmeas )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  
% Reassign value
vmeas = vout;

% Constants
pdf_start = -20;
pdf_stop  = 9;
% pdf_step  = 1;
%pdf_step = pavs(2)-pavs(1);
pdf_step = 1;

%Constans
% coeff_am = log(10)/20;
% coeff_pm = pi/180 ;

min_pavs = min(pavs);
max_pavs = max(pavs);

pavsi  = [min_pavs:pdf_step:max_pavs];
vrefi  = interp1(pavs, vref, pavsi);
vmeasi = interp1(pavs, vmeas, pavsi);
% pavsi  = pavs;
% vrefi  = vref;
% vmeasi = vmeas;

pout = 10*log10( vmeasi.^2 / (2*50) ) + 30;

vref_ideal = sqrt(10.^((pavsi-30)/10)*2*50);
gain_mag   = 20*log10( abs(vmeasi./vref_ideal));
gain_phase = 180/pi *  angle(vmeasi./vrefi);

rdb = [pdf_start:pdf_step:pdf_stop];
rlin = 10.^(rdb./10);
drlin = rlin .* 10.^(pdf_step/10) - rlin;
pdf = drlin .* rlin .* exp(-rlin);

steps_per_db = 1/pdf_step;
num_pavsi_valid = length(pavsi) - (pdf_stop/pdf_step) + floor(steps_per_db);  % not sure why +1

evm      = 0*pavsi;
evm_amam = 0*pavsi;
evm_ampm = 0*pavsi;

for k1 = 1:num_pavsi_valid
    tmp_evm_amam = 0;
    tmp_evm_ampm = 0;
    for k2 = 1:length(rdb)
        rpavs = pavsi(k1) + rdb(k2);
        if ((rpavs >= min_pavs ) && ( rpavs <= max_pavs ))
            step_offset = rdb(k2)*steps_per_db;
            delta_mag   = gain_mag(k1+step_offset)   - gain_mag(k1);
            delta_phase = gain_phase(k1+step_offset) - gain_phase(k1);
            
            tmp_evm_amam = tmp_evm_amam + pdf(k2)*delta_mag^2;
            tmp_evm_ampm = tmp_evm_ampm + pdf(k2)*delta_phase^2;
            % disp(sprintf('drlin: %g\trlin: %g\tpdf: %g\tamam: %g\tampm: %g\tevmam: %g\tevmpm: %g', drlin(k2), rlin(k2), pdf(k2), delta_mag, delta_phase, tmp_evm_amam, tmp_evm_ampm));
        end
    end 
    tmp_evm_amam = log(10)/20 * sqrt(tmp_evm_amam);
    tmp_evm_ampm = pi/180     * sqrt(tmp_evm_ampm);
    evm(k1)      = sqrt(tmp_evm_amam^2 + tmp_evm_ampm^2);
    % disp(sprintf('pavs: %g\tpout: %g\tevmAM: %g\tevmPM: %g\tevmRMS: %g',pavsi(k1),pout(k1), tmp_evm_amam, tmp_evm_ampm, evm(k1))); 
end

%end

