%% Isis Alonso 18-11-20 photometry trace making


clear all
close all
fclose all

diary AnalysisLog
% dbstop in extractDataIndividual at 89 if (size(raw490, 1) ~= size(raw405, 1))


%% define where the stuff is
exp = 'IA_EtOH_07'; 
ses_num = 1;
type = 'pavlovianConflict';
tankfolder = 'C:\Users\Isis\surfdrive\1_Experiments\IA_EtOH07\Extracted photometry data\it3\Reward';
% rat = 'R2_R7';
time = [10, 40];%REMEMBER TO CHANGE THIS DEPENDING ON THEE XP!!!!!
%limits for trace making - the first one is the time before epoch, second one is after (ALWAYS KEEP BOTH POSITIVE)
% if you want plots of all sessions at the end of the extraction

%% pre-processing
res = 64; % resampling factor
lowpass_cutoff = 3; %low-pass cut-off for Hz, e.g. 2
filt_steepness = .95; %how steep the filter transition is (0.5-1, default 0.85)
db_atten = 90; %decibel attenuation of filters (default = 60)
exc = 100; %signal to take out for pre-processing (for me is 300s = 5mins)


%% individual session data extraction 
files = dir(fullfile(tankfolder));
files(ismember({files.name}, {'.', '..'})) = [];
for i = 1:length(files) %iterate through experiment folder
    if strfind(files(i).name, exp)
         fprintf('Updating traces for %10s \n', files(i).name)
         load(fullfile(tankfolder,  [files(i).name]))
         names = strsplit(files(i).name, {'_', '-', ' ', '.'}); %divide the file name into separte character vectors

         ev = sesdat.traces;
         conversion = sesdat.conversion;
         lp_normDat = sesdat.lp_dFF;

         if ~isempty(ev) %if events variable is not empty
n = 1;
for m = 1:size(ev, 1)
    tmp = cell2mat(ev(m, 2));
    for k = 1:size(tmp, 1)
        adjts = ceil(conversion*tmp(k, 1)); % get adjusted timestamps based on the sampling rate
        try
        signal = lp_normDat(adjts-ceil(time(1)*conversion):adjts+ceil(time(end)*conversion))'; 
        tmp2(k, :) = signal;
        catch
        fprintf('Trace around the timestamp goes over length of session signal! Ommiting and deleting ts...\n')
            continue %jump to the iteration
        end
    end
    if exist('tmp2', 'var') 
        if size(tmp, 1) > size(tmp2, 1)
        tmp(end-(size(tmp, 1) - size(tmp2, 1))+1:end, :) = [];
        end 
        tmp(:,2) = ones*m; tmp(:,3:end) = [];
    ev{m,2} = [tmp , tmp2];
    clear tmp tmp2
    end
end
sesdat.traces = ev; %save

save([tankfolder '\' files(i).name(1:end-4) '.mat'], 'sesdat')
         end
    end
end



