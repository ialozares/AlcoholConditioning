% %% load data
%     if strfind(files(i).name, type)
%            load(fullfile(tankfolder,  [files(i).name]))
% %            alldata = [alldata, dat];
% %     end
fclose all
close all

%% define variables
exp = 'IA_DrPhoto02';
type = 'DREADDs Ext Test all';
%% First test vs last conditioning
csp_vta_d = [];
csm_vta_d = [];
csp_vta_s = [];
csm_vta_s = [];

csp_CT_d = [];
csm_CT_d = [];
csp_CT_s = [];
csm_CT_s = [];

csp_NAc_d = [];
csm_NAc_d = [];
csp_NAc_s = [];
csm_NAc_s = [];



% for i = 1:size(dat, 2)
%     if strfind(dat(i).Exp, 'IA_DrPhoto02')  
%         if strfind(dat(i).phase, 'FIRST TESTS')  
%             if strfind(dat(i).group, 'VTA')    
%                 if strfind(dat(i).test, 'd-clz')
%     csp_vta_d = [csp_vta_d; cell2mat(dat(i).csp(1))];
%     csm_vta_d = [csm_vta_d; cell2mat(dat(i).csm(1))];
%      csp_vta_s = [csp_vta_s; cell2mat(dat(i).csp(2))];
%     csm_vta_s = [csm_vta_s; cell2mat(dat(i).csm(2))];
%                 elseif strfind(dat(i).test, 'sal')
%     csp_vta_s = [csp_vta_s; cell2mat(dat(i).csp(1))];
%     csm_vta_s = [csm_vta_s; cell2mat(dat(i).csm(1))];
%         csp_vta_d = [csp_vta_d; cell2mat(dat(i).csp(2))];
%     csm_vta_d = [csm_vta_d; cell2mat(dat(i).csm(2))];
%                 end
%             elseif strfind(AllPhotom(i).group, 'CT')
%                 if strfind(dat(i).test, 'd-clz')
%     csp_CT_d = [csp_CT_d; cell2mat(dat(i).csp(1))];
%     csm_CT_d = [csm_CT_d; cell2mat(dat(i).csm(1))];
%         csp_CT_s = [csp_CT_s; cell2mat(dat(i).csp(2))];
%     csm_CT_s = [csm_CT_s; cell2mat(dat(i).csm(2))];
%                 elseif strfind(dat(i).test, 'sal')
%     csp_CT_s = [csp_CT_s; cell2mat(dat(i).csp(1))];
%     csm_CT_s = [csm_CT_s; cell2mat(dat(i).csm(1))];
%      csp_CT_d = [csp_CT_d; cell2mat(dat(i).csp(2))];
%     csm_CT_d = [csm_CT_d; cell2mat(dat(i).csm(2))];
%                 end
%                 elseif strfind(AllPhotom(i).group, 'NAcS')
%                 if strfind(dat(i).test, 'd-clz')
%     csp_NAc_d = [csp_NAc_d; cell2mat(dat(i).csp(1))];
%     csm_NAc_d = [csm_NAc_d; cell2mat(dat(i).csm(1))];
%           csp_NAc_s = [csp_NAc_s; cell2mat(dat(i).csp(end))];
%     csm_NAc_s = [csm_NAc_s; cell2mat(dat(i).csm(end))];
%                 elseif strfind(dat(i).test, 'sal')
%     csp_NAc_s = [csp_NAc_s; cell2mat(dat(i).csp(1))];
%     csm_NAc_s = [csm_NAc_s; cell2mat(dat(i).csm(1))];
%         csp_NAc_d = [csp_NAc_d; cell2mat(dat(i).csp(2))];
%     csm_NAc_d = [csm_NAc_d; cell2mat(dat(i).csm(2))];
% 
%                 end
%             end
% %         elseif strfind(dat(i).phase, 'Extinction')
% %             if strfind(AllPhotom(i).group, 'NAcS')
% %       csp_NAc_s = [csp_NAc_s; cell2mat(dat(i).csp(1))];
% %     csm_NAc_s = [csm_NAc_s; cell2mat(dat(i).csm(1))];
% %             end
%         end
%     end
% end

for i = 1:size(dat, 2)
    if strfind(dat(i).Exp, 'IA_DrPhoto02')  
        if strfind(dat(i).phase, 'EXT TESTS')  
%             if strfind(dat(i).group, 'VTA')    
%                 if strfind(dat(i).test, 'd-clz')
    csp_vta_d = [csp_vta_d; cell2mat(dat(i).csp(1))];
    csm_vta_d = [csm_vta_d; cell2mat(dat(i).csm(1))];
     csp_vta_s = [csp_vta_s; cell2mat(dat(i).csp(2))];
    csm_vta_s = [csm_vta_s; cell2mat(dat(i).csm(2))];
%                 elseif strfind(dat(i).test, 'sal')
%     csp_vta_s = [csp_vta_s; cell2mat(dat(i).csp(1))];
%     csm_vta_s = [csm_vta_s; cell2mat(dat(i).csm(1))];
%         csp_vta_d = [csp_vta_d; cell2mat(dat(i).csp(2))];
%     csm_vta_d = [csm_vta_d; cell2mat(dat(i).csm(2))];
%                 end
%             elseif strfind(AllPhotom(i).group, 'CT')
% %                 if strfind(dat(i).test, 'd-clz')
%     csp_CT_d = [csp_CT_d; cell2mat(dat(i).csp(1))];
%     csm_CT_d = [csm_CT_d; cell2mat(dat(i).csm(1))];
%         csp_CT_s = [csp_CT_s; cell2mat(dat(i).csp(2))];
%     csm_CT_s = [csm_CT_s; cell2mat(dat(i).csm(2))];
% %                 elseif strfind(dat(i).test, 'sal')
% %     csp_CT_s = [csp_CT_s; cell2mat(dat(i).csp(1))];
% %     csm_CT_s = [csm_CT_s; cell2mat(dat(i).csm(1))];
% %      csp_CT_d = [csp_CT_d; cell2mat(dat(i).csp(2))];
% %     csm_CT_d = [csm_CT_d; cell2mat(dat(i).csm(2))];
% %                 end
%                 elseif strfind(AllPhotom(i).group, 'NAcS')
% %                 if strfind(dat(i).test, 'd-clz')
%     csp_NAc_d = [csp_NAc_d; cell2mat(dat(i).csp(1))];
%     csm_NAc_d = [csm_NAc_d; cell2mat(dat(i).csm(1))];
%           csp_NAc_s = [csp_NAc_s; cell2mat(dat(i).csp(2))];
%     csm_NAc_s = [csm_NAc_s; cell2mat(dat(i).csm(2))];
% %                 elseif strfind(dat(i).test, 'sal')
% %     csp_NAc_s = [csp_NAc_s; cell2mat(dat(i).csp(1))];
% %     csm_NAc_s = [csm_NAc_s; cell2mat(dat(i).csm(1))];
% %         csp_NAc_d = [csp_NAc_d; cell2mat(dat(i).csp(2))];
% %     csm_NAc_d = [csm_NAc_d; cell2mat(dat(i).csm(2))];
% 
% %                 end
%             end
%         elseif strfind(dat(i).phase, 'Extinction')
%             if strfind(AllPhotom(i).group, 'NAcS')
%       csp_NAc_s = [csp_NAc_s; cell2mat(dat(i).csp(1))];
%     csm_NAc_s = [csm_NAc_s; cell2mat(dat(i).csm(1))];
%             end
        end
%     end
        end
    end
% end
% end
% end




time = linspace(-10, 20, size(csp_vta_d, 2));

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

%-----------------------------------------------------d-clz VS -----------------------------------------------------------------------

gree = [0.47,0.67,0.19];
gris = [0.65,0.65,0.65];
yel = [0.93,0.69,0.13];
dg = [0.71,0.76,0.65];
dy = [0.80,0.72,0.54];


  p = 0.01;
    thres = 8;
%bootstrap 
%% 
tmp = bootstrap_data(csp_vta_d, 5000, 0.001);
btsrp.cspVTA_d = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csp_vta_d, 1),2);
clear tmp
tmp = bootstrap_data(csp_vta_s, 5000, 0.001);
btsrp.cspVTA_s = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csp_vta_s, 1),2);
clear tmp
tmp = bootstrap_data(csm_vta_d, 5000, 0.001);
btsrp.csmVTA_d = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csm_vta_d, 1),2);
clear tmp
tmp = bootstrap_data(csm_vta_s, 5000, 0.001);
btsrp.csmVTA_s = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csm_vta_s, 1),2);
clear tmp

% tmp = bootstrap_data(csp_NAc_d, 5000, 0.001);
% btsrp.cspNAc_d = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csp_NAc_d, 1),2);
% clear tmp
% tmp = bootstrap_data(csp_NAc_s, 5000, 0.001);
% btsrp.cspNAc_s = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csp_NAc_s, 1),2);
% clear tmp
% tmp = bootstrap_data(csm_NAc_d, 5000, 0.001);
% btsrp.csmNAc_d = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csm_NAc_d, 1),2);
% clear tmp
% tmp = bootstrap_data(csm_NAc_s, 5000, 0.001);
% btsrp.csmNAc_s = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csm_NAc_s, 1),2);
% clear tmp

% 
% tmp = bootstrap_data(csp_CT_d, 5000, 0.001);
% btsrp.cspCT_d = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csp_CT_d, 1),2);
% clear tmp
% tmp = bootstrap_data(csp_CT_s, 5000, 0.001);
% btsrp.cspCT_s = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csp_CT_s, 1),2);
% clear tmp
% tmp = bootstrap_data(csm_CT_d, 5000, 0.001);
% btsrp.csmCT_d = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csm_CT_d, 1),2);
% clear tmp
% tmp = bootstrap_data(csm_CT_s, 5000, 0.001);
% btsrp.csmCT_s = CIadjust(tmp(1,:),tmp(2,:),tmp,size(csm_CT_s, 1),2);
% clear tmp


clear tmp
save(['C:\Users\Isis\surfdrive\2_Papers\LH-GABApaper\photometryExperiments\Data\'  exp ' ' type ' btsrp.mat'], 'btsrp')
%permutation tests

[perm.cspVTA, ~] = permTest_array(csp_vta_d, csp_vta_s, 1000);
[perm.csmVTA, ~]= permTest_array(csm_vta_d, csm_vta_s, 1000);
[perm.csVTA_d, ~] = permTest_array(csp_vta_d, csm_vta_d, 1000);
[perm.csVTA_s, ~]= permTest_array(csp_vta_s, csm_vta_s, 1000);

% [perm.cspNAc, ~]= permTest_array(csp_NAc_d, csp_NAc_s, 1000);
% [perm.csmNAc, ~]= permTest_array(csp_NAc_d, csm_NAc_s, 1000);
% [perm.csNAc_d, ~]= permTest_array(csp_NAc_d, csm_NAc_d, 1000);
% [perm.csNAc_s, ~]= permTest_array(csp_NAc_s, csm_NAc_s, 1000);
% 
% [perm.cspCT, ~]= permTest_array(csp_CT_d, csp_CT_s, 1000);
% [perm.csmCT, ~]= permTest_array(csp_CT_d, csm_CT_s, 1000);
% [perm.csCT_d, ~]= permTest_array(csp_CT_d, csm_CT_d, 1000);
% [perm.csCT_s, ~]= permTest_array(csp_CT_s, csm_CT_s, 1000);
% 
% [perm.cspVTACT, ~]= permTest_array(csp_vta_d, csp_CT_d, 1000);
% [perm.cspNAcCT, ~]= permTest_array(csp_NAc_d, csp_CT_d, 1000);
% [perm.csmVTACT, ~]= permTest_array(csm_vta_d, csm_CT_d, 1000);
% [perm.csmNAcCT, ~]= permTest_array(csm_NAc_d, csm_CT_d, 1000);
save(['C:\Users\Isis\surfdrive\2_Papers\LH-GABApaper\photometryExperiments\Data\'  exp ' ' type ' perm.mat'], 'perm')

  xl = ([-5 20]);
  yl = ([-2 10]);
    xt = ([-5 -4 -3 -2 -1 -0 5 10 15 20]);
%% plots VTA

% plot lines
a = figure
jbfill(time, (mean(csp_vta_d, 1)-(std((csp_vta_d), 0, 1)./sqrt(size(csp_vta_d, 1)))),(std((csp_vta_d), 0, 1)./sqrt(size(csp_vta_d, 1))+mean(csp_vta_d, 1)), dg, 'none', 0, 0.2)
hold on
jbfill(time, (mean(csm_vta_d, 1)-(std((csm_vta_d), 0, 1)./sqrt(size(csm_vta_d, 1)))),(std((csm_vta_d), 0, 1)./sqrt(size(csm_vta_d, 1))+mean(csm_vta_d, 1)), dy, 'none', 0, 0.2)
hold on
plot(time, mean(csp_vta_d, 1), 'Color', dg)
hold on
plot(time,mean(csm_vta_d, 1), 'Color', dy)
hold on
jbfill(time, (mean(csp_vta_s, 1)-(std((csp_vta_s), 0, 1)./sqrt(size(csp_vta_s, 1)))),(std((csp_vta_s), 0, 1)./sqrt(size(csp_vta_s, 1))+mean(csp_vta_s, 1)), gree, 'none')
hold on
jbfill(time, (mean(csm_vta_s, 1)-(std((csm_vta_s), 0, 1)./sqrt(size(csm_vta_s, 1)))),(std((csm_vta_s), 0, 1)./sqrt(size(csm_vta_s, 1))+mean(csm_vta_s, 1)), yel, 'none')
hold on
plot(time,mean(csp_vta_s, 1),  'Color', gree)
plot(time,mean(csm_vta_s, 1),  'Color', yel)

%  f = get(gca, 'Children')
%   legend([f(1), f(2), f(5), f(6)], 'TEST CS+','TEST CS-', 'SAL CS+', 'SAL CS-')

ylinemin = max(mean(csp_vta_d, 1)+1);
ylinemax = 2*ylinemin;

%significance bars for bootstrap
tmp = find(btsrp.cspVTA_d(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemin*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.cspVTA_d(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemax*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.cspVTA_s(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree) 
clear tmp id
tmp = find(btsrp.cspVTA_s(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(btsrp.csmVTA_d(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(btsrp.csmVTA_d(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(btsrp.csmVTA_s(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id
tmp = find(btsrp.csmVTA_s(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)

%plot significance bars for permutation
clear tmp id
tmp = find(perm.cspVTA(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.5)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.cspVTA(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.63)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csmVTA(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csmVTA(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id

clear tmp id
tmp = find(perm.csVTA_d(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.50)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.csVTA_d(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.63)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csVTA_s(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-3)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csVTA_s(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-3.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id



%some lines 
line([time(1) time(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
box off
set(gca, 'color', 'none')
  vline([0, 4], {'k'}, {'Cue', 'Alcohol'})
  title('VTA all rats')
   xlim(xl)
  ylim(yl)
  xticks(xt)


%         saveas(a, ['C:\Users\Isis\surfdrive\2_Papers\LH-GABApaper\photometryExperiments\Figures\' exp ' ' type ' cue onset.fig'])
%       saveas(a, ['C:\Users\Isis\surfdrive\2_Papers\LH-GABApaper\photometryExperiments\Figures\' exp ' ' type ' cue onset.png'])


      %% plots NAC

% plot lines
a = figure
jbfill(time, (mean(csp_NAc_d, 1)-(std((csp_NAc_d), 0, 1)./sqrt(size(csp_NAc_d, 1)))),(std((csp_NAc_d), 0, 1)./sqrt(size(csp_NAc_d, 1))+mean(csp_NAc_d, 1)), dg, 'none', 0, 0.2)
hold on
jbfill(time, (mean(csm_NAc_d, 1)-(std((csm_NAc_d), 0, 1)./sqrt(size(csm_NAc_d, 1)))),(std((csm_NAc_d), 0, 1)./sqrt(size(csm_NAc_d, 1))+mean(csm_NAc_d, 1)), dy, 'none', 0, 0.2)
hold on
plot(time, mean(csp_NAc_d, 1), 'Color', dg)
hold on
plot(time,mean(csm_NAc_d, 1), 'Color', dy)
hold on
jbfill(time, (mean(csp_NAc_s, 1)-(std((csp_NAc_s), 0, 1)./sqrt(size(csp_NAc_s, 1)))),(std((csp_NAc_s), 0, 1)./sqrt(size(csp_NAc_s, 1))+mean(csp_NAc_s, 1)), gree, 'none')
hold on
jbfill(time, (mean(csm_NAc_s, 1)-(std((csm_NAc_s), 0, 1)./sqrt(size(csm_NAc_s, 1)))),(std((csm_NAc_s), 0, 1)./sqrt(size(csm_NAc_s, 1))+mean(csm_NAc_s, 1)), yel, 'none')
hold on
plot(time,mean(csp_NAc_s, 1),  'Color', gree)
plot(time,mean(csm_NAc_s, 1),  'Color', yel)

%  f = get(gca, 'Children')
%   legend([f(1), f(2), f(5), f(6)], 'TEST CS+','TEST CS-', 'SAL CS+', 'SAL CS-')

ylinemin = max(mean(csp_NAc_d, 1)+.1);
ylinemax = 2*ylinemin;

%significance bars for bootstrap
tmp = find(btsrp.cspNAc_d(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemin*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.cspNAc_d(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemax*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.cspNAc_s(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree) 
clear tmp id
tmp = find(btsrp.cspNAc_s(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(btsrp.csmNAc_d(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(btsrp.csmNAc_d(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(btsrp.csmNAc_s(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id
tmp = find(btsrp.csmNAc_s(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)

%plot significance bars for permutation
clear tmp id
tmp = find(perm.cspNAc(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.5)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.cspNAc(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.63)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csmNAc(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csmNAc(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id

clear tmp id
tmp = find(perm.csNAc_d(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.5)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.csNAc_d(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.93)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csNAc_s(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csNAc_s(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id



%some lines 
line([time(1) time(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
box off
set(gca, 'color', 'none')
  vline([0, 4], {'k'}, {'Cue', 'Alcohol'})
  title('NAcS all rats')
  xlim(xl)
  ylim(yl)
  xticks(xt)

        %% plots CT

% plot lines
a = figure
jbfill(time, (mean(csp_CT_d, 1)-(std((csp_CT_d), 0, 1)./sqrt(size(csp_CT_d, 1)))),(std((csp_CT_d), 0, 1)./sqrt(size(csp_CT_d, 1))+mean(csp_CT_d, 1)), dg, 'none', 0, 0.2)
hold on
jbfill(time, (mean(csm_CT_d, 1)-(std((csm_CT_d), 0, 1)./sqrt(size(csm_CT_d, 1)))),(std((csm_CT_d), 0, 1)./sqrt(size(csm_CT_d, 1))+mean(csm_CT_d, 1)), dy, 'none', 0, 0.2)
hold on
plot(time, mean(csp_CT_d, 1), 'Color', dg)
hold on
plot(time,mean(csm_CT_d, 1), 'Color', dy)
hold on
jbfill(time, (mean(csp_CT_s, 1)-(std((csp_CT_s), 0, 1)./sqrt(size(csp_CT_s, 1)))),(std((csp_CT_s), 0, 1)./sqrt(size(csp_CT_s, 1))+mean(csp_CT_s, 1)), gree, 'none')
hold on
jbfill(time, (mean(csm_CT_s, 1)-(std((csm_CT_s), 0, 1)./sqrt(size(csm_CT_s, 1)))),(std((csm_CT_s), 0, 1)./sqrt(size(csm_CT_s, 1))+mean(csm_CT_s, 1)), yel, 'none')
hold on
plot(time,mean(csp_CT_s, 1),  'Color', gree)
plot(time,mean(csm_CT_s, 1),  'Color', yel)

%  f = get(gca, 'Children')
%   legend([f(1), f(2), f(5), f(6)], 'TEST CS+','TEST CS-', 'SAL CS+', 'SAL CS-')

ylinemin = max(mean(csp_CT_d, 1)+1);
ylinemax = 2*ylinemin;

%significance bars for bootstrap
tmp = find(btsrp.cspCT_d(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemin*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.cspCT_d(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), ylinemax*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(btsrp.cspCT_s(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree) 
clear tmp id
tmp = find(btsrp.cspCT_s(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.4)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(btsrp.csmCT_d(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(btsrp.csmCT_d(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-0.8)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(btsrp.csmCT_s(2,:)<0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemin-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id
tmp = find(btsrp.csmCT_s(1,:)>0);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)

%plot significance bars for permutation
clear tmp id
tmp = find(perm.cspCT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.5)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.cspCT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.93)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csmCT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csmCT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id

clear tmp id
tmp = find(perm.csCT_d(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.5)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.csCT_d(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.93)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csCT_s(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csCT_s(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id



%some lines 
line([time(1) time(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
box off
set(gca, 'color', 'none')
  vline([0, 4], {'k'}, {'Cue', 'Alcohol'})
  title('Control all rats')
  xlim(xl)
  ylim(yl)
  xticks(xt)

  %% CT vs VTA 

  % plot lines
a = figure
jbfill(time, (mean(csp_vta_d, 1)-(std((csp_vta_d), 0, 1)./sqrt(size(csp_vta_d, 1)))),(std((csp_vta_d), 0, 1)./sqrt(size(csp_vta_d, 1))+mean(csp_vta_d, 1)), dg, 'none', 0, 0.2)
hold on
jbfill(time, (mean(csm_vta_d, 1)-(std((csm_vta_d), 0, 1)./sqrt(size(csm_vta_d, 1)))),(std((csm_vta_d), 0, 1)./sqrt(size(csm_vta_d, 1))+mean(csm_vta_d, 1)), dy, 'none', 0, 0.2)
hold on
plot(time, mean(csp_vta_d, 1), 'Color', dg)
hold on
plot(time,mean(csm_vta_d, 1), 'Color', dy)
hold on
jbfill(time, (mean(csp_CT_d, 1)-(std((csp_CT_d), 0, 1)./sqrt(size(csp_CT_d, 1)))),(std((csp_CT_d), 0, 1)./sqrt(size(csp_CT_d, 1))+mean(csp_CT_d, 1)), gree, 'none')
hold on
jbfill(time, (mean(csm_CT_d, 1)-(std((csm_CT_d), 0, 1)./sqrt(size(csm_CT_d, 1)))),(std((csm_CT_d), 0, 1)./sqrt(size(csm_CT_d, 1))+mean(csm_CT_d, 1)), yel, 'none')
hold on
plot(time,mean(csp_CT_d, 1),  'Color', gree)
plot(time,mean(csm_CT_d, 1),  'Color', yel)


ylinemin = max(mean(csp_CT_d, 1)+1);
ylinemax = 2*ylinemin;

%plot significance bars for permutation
clear tmp id
tmp = find(perm.cspVTACT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.5)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.cspVTACT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.93)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csmVTACT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csmVTACT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id



%some lines 
line([time(1) time(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
box off
set(gca, 'color', 'none')
  vline([0, 4], {'k'}, {'Cue', 'Alcohol'})
  title('Control s VTA all rats')
  xlim(xl)
  ylim(yl)
  xticks(xt)

   %% CT vs NAc

  % plot lines
a = figure
jbfill(time, (mean(csp_NAc_d, 1)-(std((csp_NAc_d), 0, 1)./sqrt(size(csp_NAc_d, 1)))),(std((csp_NAc_d), 0, 1)./sqrt(size(csp_NAc_d, 1))+mean(csp_NAc_d, 1)), dg, 'none', 0, 0.2)
hold on
jbfill(time, (mean(csm_NAc_d, 1)-(std((csm_NAc_d), 0, 1)./sqrt(size(csm_NAc_d, 1)))),(std((csm_NAc_d), 0, 1)./sqrt(size(csm_NAc_d, 1))+mean(csm_NAc_d, 1)), dy, 'none', 0, 0.2)
hold on
plot(time, mean(csp_NAc_d, 1), 'Color', dg)
hold on
plot(time,mean(csm_NAc_d, 1), 'Color', dy)
hold on
jbfill(time, (mean(csp_CT_d, 1)-(std((csp_CT_d), 0, 1)./sqrt(size(csp_CT_d, 1)))),(std((csp_CT_d), 0, 1)./sqrt(size(csp_CT_d, 1))+mean(csp_CT_d, 1)), gree, 'none')
hold on
jbfill(time, (mean(csm_CT_d, 1)-(std((csm_CT_d), 0, 1)./sqrt(size(csm_CT_d, 1)))),(std((csm_CT_d), 0, 1)./sqrt(size(csm_CT_d, 1))+mean(csm_CT_d, 1)), yel, 'none')
hold on
plot(time,mean(csp_CT_d, 1),  'Color', gree)
plot(time,mean(csm_CT_d, 1),  'Color', yel)


ylinemin = max(mean(csp_CT_d, 1)+1);
ylinemax = 2*ylinemin;

%plot significance bars for permutation
clear tmp id
tmp = find(perm.cspNAcCT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.5)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dg,'Color', dg)
clear tmp id
tmp = find(perm.cspNAcCT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-1.93)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',gree,'Color', gree)
clear tmp id
tmp = find(perm.csmNAcCT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',dy,'Color', dy)
clear tmp id
tmp = find(perm.csmNAcCT(1,:)<p);
id = tmp(consec_idx(tmp, thres));
plot(time(id), (ylinemax-2.13)*ones(size(time(id),2), 2), 's', 'MarkerSize', 7, 'MarkerFaceColor',yel,'Color', yel)
clear tmp id



%some lines 
line([time(1) time(end)], [0 0], 'LineStyle', '--', 'Color', 'k');
box off
set(gca, 'color', 'none')
  vline([0, 4], {'k'}, {'Cue', 'Alcohol'})
  title('Control s NAc all rats')
  xlim(xl)
  ylim(yl)
  xticks(xt)