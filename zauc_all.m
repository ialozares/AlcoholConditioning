%%PLOT AUC and zscore  - ISIS ALONSO 16-02
fclose all
clear all
close all

%% TODO
% make so it is all saved in cells- easier and then I need less lines of
% code
% maybe put all cs data onto one cell structure?

%% define where the stuff is
exp = 'IA_EtOH07'; 
tankfolder = 'C:\Users\Isis\surfdrive\1_Experiments\IA_EtOH07\Extracted photometry data\it3\Conflict';
type = 'Condlict';
ses_num = 1;
r = 'R3';
gree = [0.47,0.67,0.19];
gris = [0.65,0.65,0.65];
yel = [0.93,0.69,0.13];
dg = [0.71,0.76,0.65];
dy = [0.80,0.72,0.54];

  p = 0.01;
    thres = 8;


% %% variables to save
dat.zauc = [];
dat.zauc2 = [];
% dat.magcsp = [];
% dat.magcsm = [];
% % dat.zpeak = [];
% almag = [];
% dat.csp = [];
% dat.csm = [];
% dat.cspresp = [];
% dat.cspom = [];
% dat.csmresp = [];
% dat.csmom = [];
% dat.magall=[];
% dat.magout =[];


%% stuff for plotting
forbar = categorical({'CS+', 'CS-'}) ;
%define colours for plotting with brewermap
% map = brewermap(9,'Pastel1');
map = [[0.47,0.67,0.19]; [0.93,0.69,0.13]];
set(0, 'DefaultAxesColorOrder', map)
gree = [0.47,0.67,0.19];
yel = [0.93,0.69,0.13];

%% load invidivual session data
files = dir(fullfile(tankfolder));
files(ismember({files.name}, {'.', '..'})) = [];
for i = 1:length(files) %iterate through experiment folder
  if strfind(files(i).name, r)
%   names = strsplit(files(i).name, {'-' , '_', ' ', '.'}); %divide the file name into separte character vectors
%   if names{5} == num2str(ses_num) || names{6} == num2str(ses_num) && names{3} == type
      load(fullfile(tankfolder,  [files(i).name]))
%    if isfile(fullfile(tankfolder, [ exp ' ' type ' session ' num2str(ses_num) ' data ' r '.mat']))
%            load(fullfile(tankfolder, [exp ' ' type ' session ' num2str(ses_num) ' data ' r '.mat']))

  names = strsplit(files(i).name, {'-' , '_', ' ', '.'}); %divide the file name into separte character vectors
  
%get variables
    data = sesdat.cond(:, 5:end);
    times = sesdat.cond(:, 1:4);
    mag = sesdat.mag;
   

%define time, baseline etc
time = linspace(-10, 40, size(data, 2)); %time vector the size of trace
base = (time >= -10) &(time <= -2);
auclims = (time >= -2) & (time <= 0); %limits of auc baseline
auclims3 = (time >= 0) & (time <= 2); %limits of auc comparison
auclims2 = (time >= 10) & (time <= 12); %limits of auc
 
%% zscore standardisation on df/f
zdata = zeros(size(data));
zbase = zeros(size(data));
tmp = 0;
    for m = 1:size(data, 1)
        zb = mean(data(m, base));
        zsd = std(data(m,base));
        for j = 1:size(data,2)
            tmp = tmp+1;
            zbase(m, tmp) = (data(m,j) -zb);
            zdata(m,tmp) = (data(m,j) - zb)/zsd;
        end
        tmp = 0;
    end
    zdata = [times, zdata];  
    
    %calculate auc during cs onset period
   
    cdat = zdata(zdata(:,2) == 1, 5:end); %cdat(cdat<0) = 0; %all cs+
    zauc.csp = trapz(cdat(:, auclims), 2);
    zauc.cspmean = mean(zauc.csp, 1);
    zauc.cspsem = std(zauc.csp, 0, 1)./sqrt(size(zauc.csp, 1));
    zauc2.csp = trapz(cdat(:, auclims2), 2);
    zauc2.cspmean = mean(zauc2.csp, 1);
    zauc2.cspsem = std(zauc2.csp, 0, 1)./sqrt(size(zauc2.csp, 1));
    dat.csp{1,ses_num} = cdat;
    
    cdat =zdata(zdata(:,2) == 2, 5:end); %cdat(cdat<0) = 0;
    zauc.csm = trapz(cdat(:, auclims), 2); %all cs-
    zauc.csmmean = mean(zauc.csm, 1);
    zauc.csmsem = std(zauc.csm, 0, 1)./sqrt(size(zauc.csm, 1));
    zauc2.csm = trapz(cdat(:, auclims2), 2); %all cs-
    zauc2.csmmean = mean(zauc2.csm, 1);
    zauc2.csmsem = std(zauc2.csm, 0, 1)./sqrt(size(zauc2.csm, 1));
     dat.csm{1,ses_num} = cdat;
    
    cdat = zdata(zdata(:,2) == 3, 5:end); %cdat(cdat<0) = 0;
    zauc.cspun = trapz(cdat(:, auclims), 2); %cs+ responded
    zauc.cspunmean = mean(zauc.cspun, 1);
    zauc.cspunpsem = std(zauc.cspun, 0, 1)./sqrt(size(zauc.cspun, 1));
    zauc2.cspun = trapz(cdat(:, auclims2), 2); %cs+ responded
    zauc2.cspunmean = mean(zauc2.cspun, 1);
    zauc2.cspunsem = std(zauc2.cspun, 0, 1)./sqrt(size(zauc2.cspun, 1));
     dat.cspun{1,ses_num} = cdat;
    
%      
%     cdat = zdata(zdata(:,2) == 1 & zdata(:,3) == 0, 5:end); %cdat(cdat<0) = 0;
%     zauc.cspom = trapz(cdat(:, auclims), 2); %cs+ not responded
%     zauc.cspommean = mean(zauc.cspom, 1);
%     zauc.cspomsem = std(zauc.cspom, 0, 1)./sqrt(size(zauc.cspom, 1));
%     zauc2.cspom = trapz(cdat(:, auclims2), 2); %cs+ not responded
%     zauc2.cspommean = mean(zauc2.cspom, 1);
%     zauc2.cspomsem = std(zauc2.cspom, 0, 1)./sqrt(size(zauc2.cspom, 1));
%      dat.cspom{1,ses_num} = cdat;

%          cdat = zdata(zdata(:,2) == 2 & zdata(:,3) ~= 0, 5:end); %cdat(cdat<0) = 0;
%     zauc.csmresp = trapz(cdat(:, auclims), 2); %cs+ responded
%     zauc.csmrespmean = mean(zauc.csmresp, 1);
%     zauc.csmrespsem = std(zauc.csmresp, 0, 1)./sqrt(size(zauc.csmresp, 1));
%     zauc2.csmresp = trapz(cdat(:, auclims2), 2); %cs+ responded
%     zauc2.csmrespmean = mean(zauc2.csmresp, 1);
%     zauc2.csmrespsem = std(zauc2.csmresp, 0, 1)./sqrt(size(zauc2.csmresp, 1));
%      dat.csmresp{1,ses_num} = cdat;
%     
     
%     cdat = zdata(zdata(:,2) == 2 & zdata(:,3) == 0, 5:end); %cdat(cdat<0) = 0;
%     zauc.csmom = trapz(cdat(:, auclims), 2); %cs+ not responded
%     zauc.csmommean = mean(zauc.csmom, 1);
%     zauc.csmomsem = std(zauc.csmom, 0, 1)./sqrt(size(zauc.csmom, 1));
%     zauc2.csmom = trapz(cdat(:, auclims2), 2); %cs+ not responded
%     zauc2.csmommean = mean(zauc2.csmom, 1);
%     zauc2.csmomsem = std(zauc2.csmom, 0, 1)./sqrt(size(zauc2.csmom, 1));
%      dat.csmom{1,ses_num} = cdat;

if ~isempty(mag)
    zmag = zscore(mag(:, 5:end));

    cdat = zmag(zmag(:,2) == 1, 5:end);
    dat.magcsp{1,ses_num} = cdat;

    cdat = zmag(zmag(:,2) == 2, 5:end);
    dat.magcsm{1,ses_num} = cdat;

    dat.magall{1, ses_num} = zmag;

    cdat = zmag(zmag(:,2) == 3 & zmag(:, 4) == 0, 5:end);
    dat.magout{1,ses_num} = cdat;
else
    fprintf('No magazine entries for sessio  %10s \n', files(i).name)
end




   

    %save variables...
    dat.sesdat{1, ses_num}= sesdat; 
    dat.zauc = [dat.zauc, zauc]; %allauc from cs
    dat.zauc2 = [dat.zauc2, zauc2];
    sesdat.phase = 'Reward';
  


    ses_num = ses_num + 1;
      end
  end
% end

dat.rat = r;
dat.phase = type;
dat.sessions = ses_num;
dat.exp = exp;
save(['C:\Users\Isis\surfdrive\1_Experiments\IA_EtOH07\Extracted photometry data\it3\'  exp ' All ' type ' sessions data ' r '.mat'], 'dat')
time = linspace(-10, 40, size(dat.csp{1}, 2));
%% -----------------------------------------------------------NOSEPOKES------------------------------------------------------------------------
n = figure;
subplot(1,3,1);
%% 
imagesc(time, 1:size(vertcat(dat.csp{:}), 1), (vertcat(dat.csp{:})));
    colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
     vline([0, 4], 'k')
    title('All CS+ individual trial acitivty')
subplot(1,3,2);
imagesc(time, 1:size(vertcat(dat.csm{:}), 1), vertcat(dat.csm{:}));
    colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
     vline([0, 4], 'k')
    title('All CS- individual trial acitivty')
    subplot(1,3,3);
imagesc(time, 1:size(vertcat(dat.cspun{:}), 1), vertcat(dat.cspun{:}));
    colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
     vline([0, 4], 'k')
    title('CS+ responded')
% subplot(2,3,2);
% imagesc(time, 1:size(vertcat(dat.csmresp{:}), 1), vertcat(dat.csmresp{:}));
%     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
%      vline([0, 4], 'k')
%     title('CS- responded')
%             subplot(2,3,3);
% imagesc(time, 1:size(vertcat(dat.cspom{:}), 1), vertcat(dat.cspom{:}));
%     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
%      vline([0, 4], 'k')
%     title('CS+ not responded')
% subplot(2,3,6);
% imagesc(time, 1:size(vertcat(dat.csmom{:}), 1), vertcat(dat.csmom{:}));
%     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
%      vline([0, 4], 'k')
%     title('CS- not responded')   
   set(findobj(gcf,'type','axes'), 'FontName', 'Helvetica',  'clim', [-3 3])
% %----------------------------------------------------individual nosepoke activity----------------------------------------------------------------

      saveas(n, [tankfolder '\photometryFigs\' exp ' All ' type ' cue onset plots ' r '.fig'])
      saveas(n, [tankfolder '\photometryFigs\' exp ' All ' type ' cue onset plots ' r '.png'])

%-----------------------------------------------------FIRST VS LAST-----------------------------------------------------------------------

%bootstrap 
%% 
tmp = bootstrap_data(dat.csp{1}, 5000, 0.001);
btsrp.cspF = CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.csp, 1),2);
clear tmp
tmp = bootstrap_data(dat.csp{end}, 5000, 0.001);
btsrp.cspL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.csp, 1),2);
clear tmp
tmp = bootstrap_data(dat.csm{1}, 5000, 0.001);
btsrp.csmF =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.csm, 1),2);
clear tmp
tmp = bootstrap_data(dat.csm{end}, 5000, 0.001);
btsrp.csmL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.csm, 1),2);
clear tmp

%permutation tests
[perm.cspFL, ~] = permTest_array(dat.csp{1}, dat.csp{end}, 1000);
[perm.csmFL, ~]= permTest_array(dat.csm{1}, dat.csm{end}, 1000);
[perm.csF, ~]= permTest_array(dat.csp{1}, dat.csm{1}, 1000);
[perm.csL, ~]= permTest_array(dat.csp{end}, dat.csm{end}, 1000);

% plot lines
a = figure
jbfill(time, (mean(dat.csp{1}, 1)-(std((dat.csp{1}), 0, 1)./sqrt(size(dat.csp{1}, 1)))),(std((dat.csp{1}), 0, 1)./sqrt(size(dat.csp{1}, 1))+mean(dat.csp{1}, 1)), dg, 'none', 0, 0.2)
hold on
jbfill(time, (mean(dat.csm{1}, 1)-(std((dat.csm{1}), 0, 1)./sqrt(size(dat.csm{1}, 1)))),(std((dat.csm{1}), 0, 1)./sqrt(size(dat.csm{1}, 1))+mean(dat.csm{1}, 1)), dy, 'none', 0, 0.2)
hold on
jbfill(time, (mean(dat.cspun{2}, 1)-(std((dat.cspun{2}), 0, 1)./sqrt(size(dat.cspun{2}, 1)))),(std((dat.cspun{2}), 0, 1)./sqrt(size(dat.cspun{2}, 1))+mean(dat.cspun{2}, 1)), 'k', 'none', 0, 0.2)
hold on
plot(time, mean(dat.csp{1}, 1), 'Color', dg)
hold on
plot(time,mean(dat.csm{1}, 1), 'Color', dy)
hold on
plot(time,mean(dat.cspun{2}, 1), 'Color', 'k')
a = figure
hold on
jbfill(time, (mean(dat.csp{end}, 1)-(std((dat.csp{end}), 0, 1)./sqrt(size(dat.csp{end}, 1)))),(std((dat.csp{end}), 0, 1)./sqrt(size(dat.csp{end}, 1))+mean(dat.csp{end}, 1)), gree, 'none')
hold on
jbfill(time, (mean(dat.csm{end}, 1)-(std((dat.csm{end}), 0, 1)./sqrt(size(dat.csm{end}, 1)))),(std((dat.csm{end}), 0, 1)./sqrt(size(dat.csm{end}, 1))+mean(dat.csm{end}, 1)), yel, 'none')
hold on
jbfill(time, (mean(dat.cspun{end}, 1)-(std((dat.cspun{end}), 0, 1)./sqrt(size(dat.cspun{end}, 1)))),(std((dat.cspun{end}), 0, 1)./sqrt(size(dat.cspun{end}, 1))+mean(dat.cspun{end}, 1)), 'r', 'none')
hold on
plot(time,mean(dat.csp{end}, 1),  'Color', gree)
plot(time,mean(dat.csm{end}, 1),  'Color', yel)
plot(time,mean(dat.cspun{end}, 1),  'Color', 'r')

%  f = get(gca, 'Children')
%   legend([f(1), f(2), f(5), f(6)], 'TEST CS+','TEST CS-', 'SAL CS+', 'SAL CS-')

ylinemin = max(mean(dat.csp{end}, 1)+2);
ylinemax = 2*ylinemin;

%significance bars for bootstrap
tmp = find(btsrp.cspF(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemin*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.cspF(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemax*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.cspL(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree) 
clear tmp id
tmp = find(btsrp.cspL(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(btsrp.csmF(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(btsrp.csmF(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(btsrp.csmL(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id
tmp = find(btsrp.csmL(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)

%plot significance bars for permutation
clear tmp id
tmp = find(perm.cspFL(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.cspFL(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.93)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csmFL(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csmFL(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id
tmp = find(perm.csL(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csL(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.33)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id
tmp = find(perm.csF(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-3)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.csF(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-3.33)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)

%some lines 
line([time(1) time(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
box off
set(gca, 'color', 'none')
  vline([0, 4], {'k'}, {'Cue', 'Alcohol'})
%   title(+ type +  " "  + r + "")

        saveas(a, [tankfolder '\photometryFigs\' exp ' FirstVsLast ' type ' cue onset plots' r '.fig'])
      saveas(a, [tankfolder '\photometryFigs\' exp ' FirstVsLast ' type ' cue onset plots' r '.png'])




        %% ----------------------------------------------------------------MAGAZINE ENTRIES-----------------------------------------------------------

  %------------------------------------------------------------------Individual plots------------------------------------------------------

 a = figure;
subplot(2,2,1);
imagesc(time, 1:size(vertcat(dat.magcsp{:}), 1), vertcat(dat.magcsp{:}));
    colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
     vline([0, 4], 'k')
    title('All CS+ magazine entry')
subplot(2,2,2);
imagesc(time, 1:size(vertcat(dat.magcsm{:}), 1), vertcat(dat.magcsm{:}));
    colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
     vline([0, 4], 'k')
    title('All CS- magazine entry')
    subplot(2,2,3);
imagesc(time, 1:size(vertcat(dat.magall{:}), 1), vertcat(dat.magall{:}));
    colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
     vline([0, 4], 'k')
    title('All magazine entries')
subplot(2,2,4);
imagesc(time, 1:size(vertcat(dat.magout{:}), 1), vertcat(dat.magout{:}));
    colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
     vline([0, 4], 'k')
    title('Magazine entries outside CS periods') 
   set(findobj(gcf,'type','axes'), 'FontName', 'Helvetica',  'clim', [-3 3])
      saveas(a, [tankfolder '\photometryFigs\' exp ' All ' type ' Magazine entries' r '.fig'])
      saveas(a, [tankfolder '\photometryFigs\' exp ' All ' type ' Magazine entries' r '.png'])


% %-----------------------------------------------------FIRST VS LAST-----------------------------------------------------------------------
% 
%bootstrap 
tmp = bootstrap_data(dat.magcsp{1}, 5000, 0.001);
btsrp.MagCspF = CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magcsp{1}, 1),2);
clear tmp
tmp = bootstrap_data(dat.magcsp{end}, 5000, 0.001);
btsrp.MagCspL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magcsp{end}, 1),2);
clear tmp
% % tmp = bootstrap_data(dat.magcsm{1}, 5000, 0.001);
% % btsrp.MagCsmF =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magcsm{1}, 1),2);
% % clear tmp
% % tmp = bootstrap_data(dat.magcsm{end}, 5000, 0.001);
% % btsrp.MagCsmL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magcsm{end}, 1),2);
% % clear tmp
tmp = bootstrap_data(dat.magout{1}, 5000, 0.001);
btsrp.MagOutF =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magout{1}, 1),2);
clear tmp
tmp = bootstrap_data(dat.magout{end}, 5000, 0.001);
btsrp.MagOutL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magout{end}, 1),2);
clear tmp

% %permutation tests
% [perm.MagCspFL, ~] = permTest_array(dat.magcsp{1}, dat.magcsp{end}, 1000);
% % [perm.MagCsmFL, ~]= permTest_array(dat.magcsm{1}, dat.magcsm{end}, 1000);
% % [perm.MagCsF, ~]= permTest_array(dat.magcsp{1}, dat.magcsm{1}, 1000);
% % [perm.MagCsL, ~]= permTest_array(dat.magcsp{end}, dat.magcsm{end}, 1000);
[perm.MagOFL,  ~] = permTest_array(dat.magout{1}, dat.magout{end}, 1000);
[perm.MagOCF,  ~] = permTest_array(dat.magout{1}, dat.magcsp{1}, 1000);
[perm.MagOCL,  ~] = permTest_array(dat.magout{end}, dat.magcsp{end}, 1000);

% % 
% plot lines
a = figure
jbfill(time, (mean(dat.magout{1}, 1)-(std((dat.magout{1}), 0, 1)./sqrt(size(dat.magout{1}, 1)))),(std((dat.magout{1}), 0, 1)./sqrt(size(dat.magout{1}, 1))+mean(dat.magout{1}, 1)), 'k', 'none', 0, 0.2)
hold on
jbfill(time, (mean(dat.magout{end}, 1)-(std((dat.magout{end}), 0, 1)./sqrt(size(dat.magout{end}, 1)))),(std((dat.magout{end}), 0, 1)./sqrt(size(dat.magout{end}, 1))+mean(dat.magout{end}, 1)), gris, 'none', 0, 0.2)

hold on
plot(time, mean(dat.magout{1}, 1), 'Color', 'k')
hold on
plot(time,mean(dat.magout{end}, 1), 'Color', gris)
hold on
jbfill(time, (mean(dat.magcsp{1}, 1)-(std((dat.magcsp{1}), 0, 1)./sqrt(size(dat.magcsp{1}, 1)))),(std((dat.magcsp{1}), 0, 1)./sqrt(size(dat.magcsp{1}, 1))+mean(dat.magcsp{1}, 1)), gree, 'none', 0, 0.2)
hold on
jbfill(time, (mean(dat.magcsp{end}, 1)-(std((dat.magcsp{end}), 0, 1)./sqrt(size(dat.magcsp{end}, 1)))),(std((dat.magcsp{end}), 0, 1)./sqrt(size(dat.magcsp{end}, 1))+mean(dat.magcsp{end}, 1)), dg, 'none', 0, 0.2)

hold on
plot(time, mean(dat.magcsp{1}, 1), 'Color', gree)
hold on
plot(time,mean(dat.magcsp{end}, 1), 'Color', dg)
% hold on
% % jbfill(time, (mean(dat.magcsm{1}, 1)-(std((dat.magcsm{1}), 0, 1)./sqrt(size(dat.magcsm{1}, 1)))),(std((dat.magcsm{1}), 0, 1)./sqrt(size(dat.magcsm{1}, 1))+mean(dat.magcsm{1}, 1)), yel, 'none')
% % hold on
% % jbfill(time, (mean(dat.magcsm{end}, 1)-(std((dat.magcsm{end}), 0, 1)./sqrt(size(dat.magcsm{end}, 1)))),(std((dat.magcsm{end}), 0, 1)./sqrt(size(dat.magcsm{end}, 1))+mean(dat.magcsm{end}, 1)), dy, 'none')
% % hold on
% % plot(time,mean(dat.magcsm{1}, 1),  'Color', yel)
% % plot(time,mean(dat.magcsm{end}, 1),  'Color', dy)
% 
% %  f = get(gca, 'Children')
% %   legend([f(1), f(2), f(5), f(6)], 'TEST CS+','TEST CS-', 'SAL CS+', 'SAL CS-')
% 
ylinemin = max(mean(dat.magcsp{1}, 1)+2);
ylinemax = 2*ylinemin;

%significance bars for bootstrap
tmp = find(btsrp.MagOutF(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemin*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','k','Color', 'k')
clear tmp id
tmp = find(btsrp.MagOutF(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemax*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','k','Color', 'k')
clear tmp id
tmp = find(btsrp.MagOutL(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris) 
clear tmp id
tmp = find(btsrp.MagOutL(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
clear tmp id
tmp = find(btsrp.MagCspF(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(btsrp.MagCspF(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(btsrp.MagCspL(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.MagCspL(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
% % tmp = find(btsrp.MagCsmF(2,:)<0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemin-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
% % clear tmp id
% % tmp = find(btsrp.MagCsmF(1,:)>0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
% % clear tmp id
% % tmp = find(btsrp.MagCsmL(2,:)<0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemin-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
% % clear tmp id
% % tmp = find(btsrp.MagCsmL(1,:)>0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
% 
% %plot significance bars for permutation
% clear tmp id
% tmp = find(perm.MagCspFL(1,:)<p);
% id = tmp(consec_idx(tmp, thres));
% plot(time(id), (ylinemax-1.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
% clear tmp id
% tmp = find(perm.MagCspFL(1,:)<p);
% id = tmp(consec_idx(tmp, thres));
% plot(time(id), (ylinemax-1.93)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% clear tmp id
% tmp = find(perm.MagOFL(1,:)<p);
% id = tmp(consec_idx(tmp, thres));
% plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','k','Color', 'k')
% clear tmp id
% tmp = find(perm.MagOFL(1,:)<p);
% id = tmp(consec_idx(tmp, thres));
% plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
% clear tmp id
% tmp = find(perm.MagOCF(1,:)<p);
% id = tmp(consec_idx(tmp, thres));
% plot(time(id), (ylinemax-2.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
% clear tmp id
% tmp = find(perm.MagOCF(1,:)<p);
% id = tmp(consec_idx(tmp, thres));
% plot(time(id), (ylinemax-2.33)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','k','Color', 'k')
% clear tmp id
% tmp = find(perm.MagOCL(1,:)<p);
% id = tmp(consec_idx(tmp, thres));
% plot(time(id), (ylinemax-3)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% clear tmp id
% tmp = find(perm.MagOCL(1,:)<p);
% id = tmp(consec_idx(tmp, thres));
% plot(time(id), (ylinemax-3.33)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
% 
%some lines 
line([time(1) time(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
box off
set(gca, 'color', 'none')
  vline(0, 'k', 'Magazine entry')
  title( "First vs second test day "  + r + "")

        saveas(a, [tankfolder '\photometryFigs\' exp ' FirstVsLast ' type ' magazine entry plots' r '.fig'])
      saveas(a, [tankfolder '\photometryFigs\' exp ' FirstVsLast ' type ' magazine entry plots' r '.png'])    
% 
% 
% 
% 
% %   %% ----------------------------------------------------------------MAGAZINE ENTRIES-----------------------------------------------------------
% % 
% %   %------------------------------------------------------------------Individual plots------------------------------------------------------
% % 
% %  a = figure;
% % subplot(2,2,1);
% % imagesc(time, 1:size(dat.magcsp{1}, 1), dat.magcsp{1});
% %     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
% %      vline([0, 4], 'k')
% %     title('All CS+ magazine entry')
% % subplot(2,2,2);
% % imagesc(time, 1:size(dat.magcsm{1}, 1), dat.magcsm{1});
% %     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
% %      vline([0, 4], 'k')
% %     title('All CS- magazine entry')
% %     subplot(2,2,3);
% % imagesc(time, 1:size(dat.magall{1}, 1), dat.magall{1});
% %     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
% %      vline([0, 4], 'k')
% %     title('All magazine entries')
% % subplot(2,2,4);
% % imagesc(time, 1:size(dat.magout{1}, 1), dat.magout{1});
% %     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
% %      vline([0, 4], 'k')
% %     title('Magazine entries outside CS periods') 
% %    set(findobj(gcf,'type','axes'), 'FontName', 'Helvetica',  'clim', [-3 3])
% %       saveas(a, [tankfolder '\photometryFigs\' exp ' All ' type ' Magazine entries' r '.fig'])
% %       saveas(a, [tankfolder '\photometryFigs\' exp ' All ' type ' Magazine entries' r '.png'])
% % 
% % 
% % %-----------------------------------------------------FIRST VS LAST-----------------------------------------------------------------------
% % 
% % %bootstrap 
% % tmp = bootstrap_data(dat.magall{1}, 5000, 0.001);
% % btsrp.MagCspF = CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magall{1}, 1),2);
% % clear tmp
% % tmp = bootstrap_data(dat.magall{2}, 5000, 0.001);
% % btsrp.MagCspL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magall{2}, 1),2);
% % clear tmp
% % tmp = bootstrap_data(dat.magcsm{1}, 5000, 0.001);
% % btsrp.MagCsmF =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magcasm{1}, 1),2);
% % clear tmp
% % tmp = bootstrap_data(dat.magcsm{2}, 5000, 0.001);
% % btsrp.MagCsmL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(dat.magcsm{2}, 1),2);
% % clear tmp
% % 
% % %permutation tests
% % [perm.MagCspFL, ~] = permTest_array(dat.magall{1}, dat.magall{end}, 1000);
% % [perm.MagCsmFL, ~]= permTest_array(dat.magcsm{1}, dat.magcsm{2}, 1000);
% % [perm.MagCsF, ~]= permTest_array(dat.magcsp{1}, dat.magcsm{1}, 1000);
% % [perm.MagCsL, ~]= permTest_array(dat.magcsp{2}, dat.magcsm{2}, 1000);
% % % 
% % % plot lines
% % a = figure
% % jbfill(time, (mean(dat.magall{1}, 1)-(std((dat.magall{1}), 0, 1)./sqrt(size(dat.magall{1}, 1)))),(std((dat.magall{1}), 0, 1)./sqrt(size(dat.magall{1}, 1))+mean(dat.magall{1}, 1)), 'k', 'none', 0, 0.2)
% % hold on
% % jbfill(time, (mean(dat.magall{end}, 1)-(std((dat.magall{end}), 0, 1)./sqrt(size(dat.magall{end}, 1)))),(std((dat.magall{end}), 0, 1)./sqrt(size(dat.magall{end}, 1))+mean(dat.magall{end}, 1)), gris, 'none', 0, 0.2)
% % % jbfill(time, (mean(dat.magcsm{1}, 1)-(std((dat.magcsm{1}), 0, 1)./sqrt(size(dat.magcsm{1}, 1))),(std((dat.magcsm{1}), 0, 1)./sqrt(size(dat.magcsm{1}, 1))+mean(dat.magcsm{1}, 1)), yel, 'none', 0, 0.2)
% % hold on
% % plot(time, mean(dat.magall{1}, 1), 'Color', 'k')
% % hold on
% % plot(time,mean(dat.magall{end}, 1), 'Color', dg)
% % hold on
% % % jbfill(time, (mean(dat.magcsp{2}, 1)-(std((dat.magcsp{2}), 0, 1)./sqrt(size(dat.magcsp{2}, 1)))),(std((dat.magcsp{2}), 0, 1)./sqrt(size(dat.magcsp{2}, 1))+mean(dat.magcsp{2}, 1)), dg, 'none')
% % % hold on
% % % jbfill(time, (mean(dat.magcsm{2}, 1)-(std((dat.magcsm{2}), 0, 1)./sqrt(size(dat.magcsm{2}, 1)))),(std((dat.magcsm{2}), 0, 1)./sqrt(size(dat.magcsm{2}, 1))+mean(dat.magcsm{2}, 1)), dy, 'none')
% % % hold on
% % % plot(time,mean(dat.magcsp{2}, 1),  'Color', dg)
% % % plot(time,mean(dat.magcspm{2}, 1),  'Color', dy)
% % 
% % %  f = get(gca, 'Children')
% % %   legend([f(1), f(2), f(5), f(6)], 'TEST CS+','TEST CS-', 'SAL CS+', 'SAL CS-')
% % 
% % ylinemin = max(mean(dat.magcsp{1}, 1));
% % ylinemax = 2*ylinemin;
% % 
% % %significance bars for bootstrap
% % tmp = find(btsrp.MagCspF(2,:)<0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), ylinemin*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','k','Color', 'k')
% % clear tmp id
% % tmp = find(btsrp.MagCspF(1,:)>0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), ylinemax*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','k','Color', 'k')
% % clear tmp id
% % tmp = find(btsrp.MagCspL(2,:)<0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemin-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris) 
% % clear tmp id
% % tmp = find(btsrp.MagCspL(1,:)>0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
% % clear tmp id
% % tmp = find(btsrp.MagCsmF(2,:)<0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemin-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
% % clear tmp id
% % tmp = find(btsrp.MagCsmF(1,:)>0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
% % clear tmp id
% % tmp = find(btsrp.MagCsmL(2,:)<0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemin-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
% % clear tmp id
% % tmp = find(btsrp.MagCsmL(1,:)>0);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
% % 
% % %plot significance bars for permutation
% % clear tmp id
% % tmp = find(perm.MagCspFL(1,:)<p);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-1.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor','k','Color', 'k')
% % clear tmp id
% % tmp = find(perm.MagCspFL(1,:)<p);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-1.93)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gris,'Color', gris)
% % clear tmp id
% % tmp = find(perm.MagCsmFL(1,:)<p);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
% % clear tmp id
% % tmp = find(perm.MagCsmFL(1,:)<p);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
% % clear tmp id
% % tmp = find(perm.MagCsF(1,:)<p);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-2.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
% % clear tmp id
% % tmp = find(perm.MagCsF(1,:)<p);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-2.33)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
% % clear tmp id
% % tmp = find(perm.MagCsL(1,:)<p);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-3)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
% % clear tmp id
% % tmp = find(perm.MagCsL(1,:)<p);
% % id = tmp(consec_idx(tmp, thres));
% % plot(time(id), (ylinemax-3.33)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
% % 
% % %some lines 
% % line([time(1) time(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
% % box off
% % set(gca, 'color', 'none')
% %   vline(0, 'k', 'Magazine entry')
% %   title( "First vs second test day "  + r + "")
% % 
% %         saveas(a, [tankfolder '\photometryFigs\' exp ' FirstVsLast ' type ' magazine entry plots' r '.fig'])
% %       saveas(a, [tankfolder '\photometryFigs\' exp ' FirstVsLast ' type ' magazine entry plots' r '.png'])    