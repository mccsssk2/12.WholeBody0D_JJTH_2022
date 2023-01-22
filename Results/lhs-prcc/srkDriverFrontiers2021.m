%
% Sanjay R. Kharche. Jan 2021.
%
% Call the Denise K. PRCC function the best you can. See the comments. for FIMH 2021.
%
%
clear all
clear all
close all
close all
%
%
wid = 6;
fs = 48;
%
% load data.
% first 82 columns are inputs. Next 43 columns are outputs.
datat035 = load('../Frontiers2020Dialysis_PRCC_Jan7_2020/rundir0/dir_0_35.5/dir035.dat');
datat037 = load('../Frontiers2020Dialysis_PRCC_Jan7_2020/rundir0/dir_0_37.5/dir037.dat');
datat135 = load('../Frontiers2020Dialysis_PRCC_Jan7_2020/rundir0/dir_1_35.5/dir135.dat');
datat137 = load('../Frontiers2020Dialysis_PRCC_Jan7_2020/rundir0/dir_1_37.5/dir137.dat');
	 
LHSmatrix035 = datat035(:,1:82);
LHSmatrix037 = datat037(:,1:82);
LHSmatrix135 = datat135(:,1:82);
LHSmatrix137 = datat137(:,1:82);
%
%
% now run thru' outputs.
% 
for i=1:1:8

if(i==3||i==6) continue; end;

output037 = datat037(:,82+i)';
output035 = datat035(:,82+i)';
output137 = datat137(:,82+i)';
output135 = datat135(:,82+i)';
%
% Parameter Labels. This is an array you may need if making figures inside the PRCC function. I will make my figures here in the driver.
PRCC_var=1:1:82;
[prcc037 signn sign_label]=PRCC(LHSmatrix037,output037,1,PRCC_var, 100);
[prcc035 signn sign_label]=PRCC(LHSmatrix035,output035,1,PRCC_var, 100);
[prcc137 signn sign_label]=PRCC(LHSmatrix137,output137,1,PRCC_var, 100);
[prcc135 signn sign_label]=PRCC(LHSmatrix135,output135,1,PRCC_var, 100);
%
%
% set all NaN to 0.
prcc037(isnan(prcc037)) = 0;
prcc035(isnan(prcc035)) = 0;
prcc137(isnan(prcc137)) = 0;
prcc135(isnan(prcc135)) = 0;
%
% rank all of them, pick first 10 from prcc037 (your control.).
[~,idx] = sort(abs(prcc037));
idx_10 = idx(end:-1:end-4);
% idx_10
prcc037ranked = prcc037(idx_10);
prcc035ranked = prcc035(idx_10);
prcc137ranked = prcc137(idx_10);
prcc135ranked = prcc135(idx_10);
%
%
figure('rend','painters','pos',[1 1 2000 800]);
bar(prcc037ranked, 'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'barWidth',0.5, 'linewidth',3);
hold on;
bar(prcc035ranked, 'EdgeColor','none','FaceColor','blue','barWidth',0.3);
hold on;
bar(prcc137ranked, 'EdgeColor','none','FaceColor','red','barWidth',0.175);
hold on;
bar(prcc135ranked, 'EdgeColor','none','FaceColor','green','barWidth',0.075);
% set(gca,'xticklabel',{idx_10(1), idx_10(2), idx_10(3), idx_10(4), idx_10(5)});
if(i==2)
set(gca,'xticklabel',{      'E_{dias,rv}'  ,    'IHR'  ,      'R_{up1}'  ,     'E_{sys,rv}'  ,        'G'         });
elseif(i==1)
set(gca,'xticklabel',{      'IHR'  ,    'G'  ,      'R_{sp1}'  ,     'R_{ll1}'  ,        'E_{dias,rv}'         });
elseif(i==4)
set(gca,'xticklabel',{      'G'  ,    'E_{dias,rv}'  ,      'R_{sp1}'  ,     'C_{a}'  ,        'R_{ll1}'         });
elseif(i==5)
set(gca,'xticklabel',{      'G'  ,    'IHR'  ,      'R_{sp1}'  ,     'C_{a}'  ,        'R_{ll1}'         });
elseif(i==7)
set(gca,'xticklabel',{      'IHR'  ,    'G'  ,      'E_{dias,rv}'  ,     'R_{up1}'  ,        'R_{sp1}'         });
elseif(i==8)
set(gca,'xticklabel',{      'R_{inlet}'  ,    'IHR'  ,      'R_{sp1}'  ,     'E_{dias,rv}'  ,        'R_{ll1}'         });
end;
yline(0,'linewidth',wid);
ax			= gca;
ax.LineWidth 	= wid;
ax.TickDir 	= 'out';
set(gca,'TickDir','out', 'FontSize', fs, 'linewidth', wid);
if(i==1)
legend('baseline','TH','HD','HD&TH', 'Location','north', 'Orientation', 'horizontal');
legend boxoff;
else
legend off;
end;

if(i==2||i==4||i==7) ylabel('PRCC.'); end;

clear idx;
clear idx_10;

text(11.25, -1.2,'.');
text(11.25,  1.25,'.');
if(i==1)
text(3,-0.8,'B. Heart rate.', 'FontSize', fs);
elseif(i==2)
text(3,-0.8,'A. Cardiac output.', 'FontSize', fs);
elseif(i==4)
text(3,-0.8,'C. Systolic pressure.', 'FontSize', fs);
elseif(i==5)
text(3,-0.8,'D. Diastolic pressure.', 'FontSize', fs);
elseif(i==7)
text(3,-0.8,'E. Aortic shear.', 'FontSize', fs);
elseif(i==8)
text(3,-0.8,'F. Renal shear.', 'FontSize', fs);
end;

ylim([-1.1 1.1]);
yticks([-0.75 0.75]);
box off;
fname = sprintf('Frontiers_PRCC%d.png', i);
saveas(gcf, fname);
clear fname;
hold off;
close all;

end;
