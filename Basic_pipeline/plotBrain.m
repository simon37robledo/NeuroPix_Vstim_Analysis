function plotBrain()
close all

path_BrainPV103_45 = '\\132.66.45.127\data\Large_scale_mapping_NP\Model_brain\PV_103_imaged_inclined_45TA.stl';

figure(1);
axes('Parent', gcf);
hold on;

brain = stlread(path_BrainPV103_45);


patch("Faces",brain.ConnectivityList,"Vertices",brain.Points,'FaceLighting',   'gouraud', ...
    'EdgeColor',       'none', ...
    'FaceColor',       [0.8 0.8 1.0], ...
    'AmbientStrength', 0.15, 'FaceAlpha', 0.3);

camlight('headlight');
material('dull');
% Fix the axes scaling, and set a nice view angle
axis('image');
sideview = [90 0];
topview = [-90 90];
frontview = [0 0];
view(sideview)

Ylims = ylim;
Xlims = xlim;
Zlims = zlim;

ymax = 80;
zmin =-65;

ylim([Ylims(1)+15 ymax])
zlim([zmin Zlims(2)])

% yticks([5:5:ymax-30])
% zticks([zmin:5:round(Zlims(2))])
xticks([round(Xlims(1)):5:round(Xlims(2))])

grid on
xlabel('Y (mm)')
ylabel('X (mm)')
zlabel('Z (mm)')
xticklabels(0:0.5:Xlims(2)/10);
yticklabels(0.5:0.5:ymax/10);
zticklabels(round((zmin - Zlims(2)))/10:0.5:0);
sizeDot = 50;
sizeText = 13;
offsetTextX = -1.3;
offsetTextY = 0;
offsetTextZ = 1;
end