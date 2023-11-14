fclose all
close all

% define stuff for plots
 xl = ([-5 20]);
  yl = ([-2 8]);
    xt = ([-5 -4 -3 -2 -1 -0 5 10 15 20]);
    gree = [0.47,0.67,0.19];
gris = [0.65,0.65,0.65];
yel = [0.93,0.69,0.13];
dg = [0.71,0.76,0.65];
dy = [0.80,0.72,0.54];

%% define stuff

csp_first = [];
csp_last = [];
csm_first = [];
csm_last = [];




exp = 'IA_DrPhoto02';
type = 'EXT vs RE';
for i = 1:size(dat, 2)
    if strfind(dat(i).Exp, exp)  
        if strfind(dat(i).phase, 'Extinction')    
     csp_first = [csp_first; cell2mat(dat(i).csp(end))];
     csm_first = [csm_first; cell2mat(dat(i).csm(end))];
        elseif strfind(dat(i).phase, 'FE RE')  
      csp_last = [csp_last; cell2mat(dat(i).csp(2))];
     csm_last = [csm_last; cell2mat(dat(i).csm(2))];
        end
    end
end 

% for i = 1:size(AllPhotom, 2)
%     if strfind(AllPhotom(i).Exp, 'IA_DrPhoto02')  
%             if strfind(AllPhotom(i).phase, 'FIRST TESTS')      
%     csp_first = [csp_first; cell2mat(AllPhotom(i).csp(1))];
%     csm_first = [csm_first; cell2mat(AllPhotom(i).csm(1))];
%             elseif strfind(AllPhotom(i).phase, 'Conditioning')
%      csp_last = [csp_last; cell2mat(AllPhotom(i).csp(end))];
%      csm_last = [csm_last; cell2mat(AllPhotom(i).csm(end))];
%             end
%          elseif strfind(AllPhotom(i).Exp, 'IA_EtOH06')
%              if strfind(AllPhotom(i).phase, 'Conditioning')  
%      csp_first = [csp_first; cell2mat(AllPhotom(i).csp(1))];
%     csm_first = [csm_first; cell2mat(AllPhotom(i).csm(1))];
% %              elseif strfind(AllPhotom(i).phase, 'Extinction')  
%      csp_last = [csp_last; cell2mat(AllPhotom(i).csp(end))];
%      csm_last = [csm_last; cell2mat(AllPhotom(i).csm(end))];
%              end
%         end
% end




time = linspace(-10, 20, size(csp_first, 2));

%% -----------------------------------------------------------CS ONSET------------------------------------------------------------------------

% %----------------------------------------------------individual trial activity----------------------------------------------------------------
% n = figure
% % for i = 1:size(dat.cs, 2)
% subplot(2,3,1);
% imagesc(time, 1:size(vertcat(csp_first, 1), (vertcat(csp_first)));
%     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
%      vline([0, 4], 'k')
%     title('All CS+ individual trial acitivty')
% subplot(2,3,4);
% imagesc(time, 1:size(vertcat(csm_first, 1), vertcat(csm_first));
%     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
%      vline([0, 4], 'k')
%     title('All CS- individual trial acitivty')
%     subplot(2,3,2);
% imagesc(time, 1:size(vertcat(csp_last, 1), vertcat(dat.cspresp{:}));
%     colormap(brewermap([], '*RdBu')); box off; colorbar; ylabel('Trial number')
%      vline([0, 4], 'k')
%     title('CS+ responded')
% subplot(2,3,5);
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
%    set(findobj(gcf,'type','axes'), 'FontName', 'Helvetica',  'clim', [-3 3])
% % %       saveas(n, [tankfolder '\photometryFigs\' exp ' All ' type ' cue onset plots' r '.fig'])
% % %       saveas(n, [tankfolder '\photometryFigs\' exp ' All ' type ' cue onset plots' r '.png'])

%-----------------------------------------------------FIRST VS LAST-----------------------------------------------------------------------



  p = 0.01;
    thres = 8;
%bootstrap 
%% 
tmp = bootstrap_data(csp_first, 5000, 0.001);
btsrp.cspF = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csp_first, 1),2);
clear tmp
tmp = bootstrap_data(csp_last, 5000, 0.001);
btsrp.cspL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(csp_last, 1),2);
clear tmp
tmp = bootstrap_data(csm_first, 5000, 0.001);
btsrp.csmF =CIadjust(tmp(1,:),tmp(2,:),tmp,size(csm_first, 1),2);
clear tmp
tmp = bootstrap_data(csm_last, 5000, 0.001);
btsrp.csmL =CIadjust(tmp(1,:),tmp(2,:),tmp,size(csm_last, 1),2);
clear tmp
save(['C:\Users\Isis\surfdrive\2_Papers\LH-GABApaper\photometryExperiments\Data\' ' ' exp ' ' type ' btsrp.mat'], 'btsrp')
%permutation tests
[perm.cspFL, ~] = permTest_array(csp_first, csp_last, 1000);
[perm.csmFL, ~]= permTest_array(csm_first, csm_last, 1000);
[perm.csF, ~]= permTest_array(csp_first, csm_first, 1000);
[perm.csL, ~]= permTest_array(csp_last, csm_last, 1000);
save(['C:\Users\Isis\surfdrive\2_Papers\LH-GABApaper\photometryExperiments\Data\' ' ' exp ' ' type ' perm.mat'], 'perm')
%% plots

% plot lines
a = figure
jbfill(time, (mean(csp_first, 1)-(std((csp_first), 0, 1)./sqrt(size(csp_first, 1)))),(std((csp_first), 0, 1)./sqrt(size(csp_first, 1))+mean(csp_first, 1)), dg, 'none', 0, 0.2)
hold on
jbfill(time, (mean(csm_first, 1)-(std((csm_first), 0, 1)./sqrt(size(csm_first, 1)))),(std((csm_first), 0, 1)./sqrt(size(csm_first, 1))+mean(csm_first, 1)), dy, 'none', 0, 0.2)
hold on
plot(time, mean(csp_first, 1), 'Color', dg)
hold on
plot(time,mean(csm_first, 1), 'Color', dy)
hold on
jbfill(time, (mean(csp_last, 1)-(std((csp_last), 0, 1)./sqrt(size(csp_last, 1)))),(std((csp_last), 0, 1)./sqrt(size(csp_last, 1))+mean(csp_last, 1)), gree, 'none')
hold on
jbfill(time, (mean(csm_last, 1)-(std((csm_last), 0, 1)./sqrt(size(csm_last, 1)))),(std((csm_last), 0, 1)./sqrt(size(csm_last, 1))+mean(csm_last, 1)), yel, 'none')
hold on
plot(time,mean(csp_last, 1),  'Color', gree)
plot(time,mean(csm_last, 1),  'Color', yel)

%  f = get(gca, 'Children')
%   legend([f(1), f(2), f(5), f(6)], 'TEST CS+','TEST CS-', 'SAL CS+', 'SAL CS-')

ylinemin = max(mean(csp_last, 1)+1);
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
plot(time(id), (ylinemax-1.5)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
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
     xlim(xl)
  ylim(yl)
  xticks(xt)
%   title('First vs last ' + type + 'all rats')

        saveas(a, ['C:\Users\Isis\surfdrive\2_Papers\LH-GABApaper\photometryExperiments\Figures\' exp ' ' type ' cue onset.fig'])
      saveas(a, ['C:\Users\Isis\surfdrive\2_Papers\LH-GABApaper\photometryExperiments\Figures\' exp ' ' type ' cue onset.png'])

