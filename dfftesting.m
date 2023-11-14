%% Little script to test the df/f on habituation days for 2 set ups 

clear all
%% pre-processing
res = 64; % resampling factor
lowpass_cutoff = [3]; %low-pass cut-off for Hz, e.g. 2
filt_steepness = .95; %how steep the filter transition is (0.5-1, default 0.85)
db_atten = 90; %decibel attenuation of filters (default = 60)
% base = 7; %baseline for z-scoring

tankfolder = 'C:\TDT\Synapse\Tanks\Pavlovian_conditioning-211110-110100\IA_EtOH_07_R3_R4-211129-161013';

[data, ts, conversion] = extractdf(tankfolder); % extract data from both set ups
%% filterting and traces
for l = 1:2 %go through the rats
    if l == 1
        side = 'L';
        filt490 = cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '7'), 2)); %%changed for nic
        filt405 =cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '0'), 2));       
    else 
      side = 'R';
      filt490 = cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '7'), 2));
       filt405 =cell2mat(data(contains(data(:,1), side) & contains(data(:,1), '0'), 2));
    end 

% Fit control to signal and calculate DFF
cfFinal = controlfit(filt490, filt405);
normDat= deltaFF(filt490,cfFinal);

%Highpass filter and moving average...
[hp_normDat, mov_normDat] = hpFilter(ts, normDat);

%lowpassfilter 
lp_normDat = lpFilter(hp_normDat, conversion, lowpass_cutoff,...
    filt_steepness, db_atten);
% figure
% plot(lp_normDat)

sesdat.filt490 = filt490;
sesdat.filt405 = filt405;

% visual inspection of the signal and control channel
figure
a = plot(filt490);
hold on
a = plot(filt405);

figure
plot(lp_normDat)
end





