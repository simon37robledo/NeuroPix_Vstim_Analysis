 cluster_info = readtable('\\132.66.45.127\data\Large_scale_mapping_NP\Awake_exp\SA6\SA6_experiment_4_9_23\Insertion4\catgt_SA6_experiment_4_9_23_4_g0\cluster_info.tsv',...
     "FileType","text",'Delimiter', '\t');


  GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

  %%
  GNa = [138, 48, 35, 31, 67, 65, 72, 122, 28,40];
  TNa = [395, 281, 232, 228, 233, 207, 236, 337, 145,275];
  INa = [1, 1, 2, 3, 1 , 1, 2, 3, 4, 1];

  meanAnG = mean(GNa);

  meanAnT = mean(TNa);

  GNi = [45,61,91,110,128,116,128,213,236,258,267,201,233,200,252,315,99,229,137,174];
  TNi = [208,237,359,345,347,303,299,519,510,540,509,509,631,484,518,659,192,485,342,429];
  INi = [1, 2,3,4,5,6,7,1,2,1,2,1,2,3,4,5,6,7];
  Anim_i = ['PV103','PV103','PV103','PV103','PV103','PV103','PV103','PV102','PV102','PV67', 'PV67','PV67','PV67','PV27',...
      'PV27','PV27','PV27','PV27','PV27','PV27'];

  meanIG = mean(GNi);

  meanAnIT = mean(TNi);

  GNaw = [54,223,183,205,246];
  TNaw = [254,435,423,402,456];
  INaw = [1,1,2,3,1];
  Anim_aw = ['SA5','SA6','SA6','SA6','SA6'];


% % Calculate means and standard errors for each group
% mean_G = [mean(G_anesthetized), mean(G_immobilized), mean(G_awake)];
% std_G = [std(G_anesthetized), std(G_immobilized), std(G_awake)];
data = [GNa, GNi, GNaw];
datac = {GNa, GNi, GNaw};

grp = [zeros(1,length(GNa)),ones(1,length(GNi)), ones(1,length(GNaw))+1];

anL = sprintf('AN: %d I. & 5 E.',length(GNa));
imL = sprintf('IM: %d I. & 4 E.',length(GNi));
awL = sprintf('AW: %d I. & 3 E.',length(GNaw));

row1 = {'Anesthetized', 'Immobilized', 'Awake'};
row2 = {anL,imL,awL};

labelArray = [row1; row2];

tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));

% Create violin plots

h = boxplot(data, grp, 'Notch', 'on','Labels',{anL,imL,awL}, 'ExtremeMode','clip',...
    'Colors',[0.5 0.7 0.7]);
set(h, 'LineWidth', 2); 
grid on
hold on
f = plotSpread(datac);
set(findall(1,'type','line'),'markerSize',16);
ylabel('# Good units per insertion')
set(f{1},'color',[0 0.3 0.6]);
set(gcf, 'Color', 'white');


prettify_plot

%exportToPPTX(['num_GoddUnitsperInsertionBoxPlot.jpg'],)

sum(GNi)+sum(GNaw)
  