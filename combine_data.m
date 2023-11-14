%% Isis Alonso-Lozares 04-04-22 
%combine data saved from each session from each rat in one big structure so
%then it can be analysed 

clear all
close all

%% TODO


%% define where the stuff is
% exp = 'Nic_20 Photom'; %experiment
tankfolder = 'C:\Users\Isis\surfdrive\1_Experiments\IA_EtOH07\Extracted photometry data\it3\AllSessions';
% type = 'FE RE';
alldata = [];

%% load individual session data 
files = dir(fullfile(tankfolder));
files(ismember({files.name}, {'.', '..'})) = [];
for i = 1:length(files) %iterate through experiment folder
    if files(i).isdir == 0 
           load(fullfile(tankfolder,  [files(i).name]))
           alldata = [alldata, dat];
    end
end
%% 
