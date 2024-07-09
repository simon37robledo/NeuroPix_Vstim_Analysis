%% plot insertions on brain


%%
close all
plotInser = 0;
stimPlot = 1;
path_brainB = '\\132.66.45.127\data\Large_scale_mapping_NP\Model_brain\PV_Brain_6mm_DVR_.stl';
path_brainO = '\\132.66.45.127\data\Large_scale_mapping_NP\Model_brain\PV_Brain_original_DVR_.stl';
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
%%
%%%%%%%%%PV103 No micromanipulator - imaging:

Ins1PV103 = [39.74 49 -21.44];

Ins2PV103 = [41.72 53.4 -20.94];

Ins3PV103 =[42.8 48.82 -22.05];

Ins4PV103 = [36.6 49.9 -21.4];

Ins5PV103 = [41.85 59.2 -21.8];

Ins6PV103 = [34.51 51.4 -21.38];

Ins7PV103 = [35.21 44.18 -23.85];

InscoorPV103 = [Ins1PV103;Ins2PV103;Ins3PV103;Ins4PV103;Ins5PV103;Ins6PV103;Ins7PV103];

%%%%%%%PV67 Microomanipulator coor

%1
Ins1PV67 = [44.24 50.51 -21.98]; %Ref

%2
xc = Ins1PV67(1)-853/100;
yc = Ins1PV67(2);

zcoor =  surfaceZ(xc,yc, brain);

Ins2PV67 = [xc yc zcoor];

%3
xc = Ins1PV67(1)-938/100;
yc = Ins1PV67(2)+994/100;

zcoor =  surfaceZ(xc,yc, brain);

Ins3PV67 = [xc yc zcoor];

%4
xc =Ins1PV67(1);
yc = Ins1PV67(2)+994/100;

zcoor =  surfaceZ(xc,yc, brain);

Ins4PV67 = [xc yc zcoor];


InscoorPV67 = [Ins1PV67; Ins2PV67; Ins3PV67; Ins4PV67];

%%%%%%%%%%%PV27 Microomanipulator coor: (ins1 = insertion 2 in lab notebook

InsRefPV27 = [39.18 57.17 -21.18];

%1
[disx,disy] = convMan2coor(394,8.36); %z and x from manipulator

xc = InsRefPV27(1)+disx;
yc = InsRefPV27(2)-disy;

zcoor =  surfaceZ(xc,yc, brain);

Ins1PV27 = [xc, yc, zcoor];

%2
Ins2PV27 = InsRefPV27;

%3
[disx,disy] = convMan2coor(-28.1,222.25); %z and x from manipulator

xc = InsRefPV27(1)+disx;
yc = InsRefPV27(2)-disy;

zcoor =  surfaceZ(xc,yc, brain);

Ins3PV27 = [xc, yc, zcoor];

%4
[disx,disy] = convMan2coor(-298.77,-1350.36); %z and x from manipulator

xc = InsRefPV27(1)+disx;
yc = InsRefPV27(2)-disy;

zcoor =  surfaceZ(xc,yc, brain);

Ins4PV27 = [xc, yc, zcoor];

%5
[disx,disy] = convMan2coor(-94.7,-264.18); %z and x from manipulator

xc = InsRefPV27(1)+disx;
yc = InsRefPV27(2)-disy;

zcoor =  surfaceZ(xc,yc, brain);

Ins5PV27 = [xc, yc, zcoor];

%6
[disx,disy] = convMan2coor(-94.7,319.15); %z and x from manipulator

xc = InsRefPV27(1)+disx;
yc = InsRefPV27(2)-disy;

zcoor =  surfaceZ(xc,yc, brain);

Ins6PV27 = [xc, yc, zcoor];

%7
[disx,disy] = convMan2coor(265.26,41); %z and x from manipulator

xc = InsRefPV27(1)+disx;
yc = InsRefPV27(2)-disy;

zcoor =  surfaceZ(xc,yc, brain);

Ins7PV27 = [xc, yc, zcoor];


InscoorPV27 = [Ins1PV27;Ins2PV27;Ins3PV27;Ins4PV27;Ins5PV27;Ins6PV27;Ins7PV27]; 

%%%%%%%%%%%%

animals = {InscoorPV67,InscoorPV103,InscoorPV27};

colorA = {[0 0.9 0.2],[0.9 0.2 0],[0 0.2 0.9]};
in=0;

if plotInser ==1

    for a = 1:length(animals)
        scatter3(animals{a}(:,1),animals{a}(:,2),animals{a}(:,3), 'bo', 'filled', 'MarkerEdgeColor','w', 'MarkerFaceColor',colorA{a}, 'LineWidth',2, 'SizeData',sizeDot);
        text(animals{a}(:,1)+offsetTextY,animals{a}(:,2)-offsetTextX,animals{a}(:,3)+offsetTextZ, string((1:length(animals{a}))+in),...
            'Color', 'k','FontSize',sizeText);

        in = in+length(animals{a});
    end

    prettify_plot

end
%%%%%%%%%%%%%% Depth and angle to plot probe trajectory


InsAnglesPV67 = [88 88 88 88];
 
InsDepthsPV67 = [3956 3907 4000 3934]; 

InsAnglesPV103= [70.5 70.5 70.5 70.5 70.5 70.5 70.5 ];
 
InsDepthsPV103= [4104.14 3964.18 4066.5 4123.87 4175.86 4225.55 4027.42]; 

InsAnglesPV27= [72.5 72.5 72.5  72.5 72.5 72.5 72.5];
 
InsDepthsPV27= [3913.9 3904.34 3525.7 3900.1 3914.34 3739.3 3906.8];

InsAnglesPV139 = [89 89];

InsDepthsPV139 = [3910 3907.23];

InsAnglesPV59 = [81 81 81];

InsDepthsPV59 = [2845.89 2858.91+500 3050.36];

InsAnglesPV32 = [69 72 72];

InsDepthsPV32 = [2400 2600 2200];


%%
animalD = {InsDepthsPV67,InsDepthsPV103,InsDepthsPV27};

animalA = {InsAnglesPV67,InsAnglesPV103,InsAnglesPV27};
plotInser =1;
if plotInser ==1

    for a = 1:length(animalA) %X doesn't change


        point1 = [animals{a}(:,1),animals{a}(:,2),animals{a}(:,3)];

        verticalDepth = sin(deg2rad(animalA{a})).*(animalD{a}/100); %depth of unit along vertical axis

        YDist = cos(deg2rad(animalA{a})).*(animalD{a}/100); %X distance of unit from insertion

        point2 = [animals{a}(:,1),animals{a}(:,2)+YDist',animals{a}(:,3)-verticalDepth'];

        for i = 1:length(animals{a})
           

            %plot3([point1(i,1),point2(i,1)],[point1(i,2),point2(i,2)],[point1(i,3),point2(i,3)],'LineWidth', 2) %'Color',colorA{a});
            x= linspace(point1(i,1),point2(i,1),size(rLFP,1)/4);
            y= linspace(point1(i,2),point2(i,2),size(rLFP,1)/4);
            z = linspace(point1(i,3),point2(i,3),size(rLFP,1)/4);
            c = mean(rLFP,2);
            %patch(x,y,z,c,'FaceColor','none','EdgeColor','interp');colorbar
             i =2;
            surface([x(:), x(:)], [y(:), y(:)], [z(:), z(:)], ...
                [c(:), c(:)], 'EdgeColor','flat', 'FaceColor','none','LineWidth',6);cb=colorbar;
            cb.Label.String = 'microVolts';
             set(gcf,'Color','white');
            hold on
        end

    end

end
%%%%% Plot units across depth

basic_pathPV102 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV102';
expPV102 = 'PV102_experiment_18_7_23';

basic_pathPV103 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV103';
expPV103 = 'PV103_Experiment_12_6_23';

basic_pathPV67 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV67';
expPV67= 'PV67_experiment_5_7_23';

basic_pathPV27 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV27';
expPVPV27 = 'PV27_Experiment_25_6_23';

basicPathA = {basic_pathPV67,basic_pathPV103,basic_pathPV27};

expA = {expPV67,expPV103,expPVPV27};

    num_colors = 64;
    valuesc = linspace(0, 1, num_colors); % Dark green to white
    darkNum = 0.2;
    invertOrder = linspace(64,1,64);
    colorComp = linspace(darkNum,1,num_colors);

    % Define the custom colormap from dark green to white
    color_mapPV67 = zeros(num_colors, 3); % Initialize the colormap

    % Create a gradient from dark green to white

    
    
    color_mapPV67(:, 2) = darkNum + (1-darkNum)*valuesc; % Set the green channel to the gradient
    color_mapPV67(:,1) = valuesc;
    color_mapPV67(:,3) = valuesc;


    color_mapPV67 = color_mapPV67(invertOrder,:);


    % Define the custom colormap from dark red to white
    color_mapPV103 = zeros(num_colors, 3); % Initialize the colormap

    % Create a gradient from dark red to white
    color_mapPV103(:, 1) = darkNum+(1-darkNum)*valuesc; % Set the red channel to the gradient
    color_mapPV103(:,2) = valuesc;
    color_mapPV103(:,3) = valuesc;

    color_mapPV103 = color_mapPV103(invertOrder,:);

     % Define the custom colormap from dark red to white
    color_mapPV27 = zeros(num_colors, 3); % Initialize the colormap

    % Create a gradient from dark blue to white
    
    color_mapPV27(:, 3) = darkNum + (1-darkNum)*valuesc; % Set the blue channel to the gradient
    color_mapPV27(:,1) = valuesc;
    color_mapPV27(:,3) = valuesc;

    color_mapPV27 = color_mapPV27(invertOrder,:);

    colormapsA = {color_mapPV67,color_mapPV103,color_mapPV27};

    figure(2)
    copyobj(findobj(1, 'Type', 'axes'), 2);

    colormapsGray = flipud(gray(64));


    % Define the custom colormap

    m = 64;
    % Colors
    color_palette = [1/2 0 0;   % Deep red
        1 0 0;     % Red
        1 1 0;     % Yellow
        1 1 1;     % White
        0 1 1;     % Cyan
        0 0 1;     % Blue
        0 0 1/2];  % Deep blue

    % Compute distributions along the samples
    color_dist = cumsum([0 1/10 1/5 1/5 1/5 1/5 1/10]);
    color_samples = round((m-1)*color_dist)+1;
    % Make the gradients
    J = zeros(m,3);
    J(color_samples,:) = color_palette(1:7,:);
    diff_samples = diff(color_samples)-1;
    for d = 1:1:length(diff_samples)
        if diff_samples(d)~=0
            color1 = color_palette(d,:);
            color2 = color_palette(d+1,:);
            G = zeros(diff_samples(d),3);
            for idx_rgb = 1:1:3
                g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
                g([1, length(g)]) = [];
                G(:,idx_rgb) = g';
            end
            J(color_samples(d)+1:color_samples(d+1)-1,:) = G;
        end
    end
    J = flipud(J);

    figPosition =  [-1356 329  985 518];

    if stimPlot == 1
        stimfigNumbers = 3:6;

        lims = Xlims;
        figure(3)
        copyobj(findobj(1, 'Type', 'axes'), 3);
        xlim([lims(2)/2 lims(2)])
        title('fullFieldFlash')

        prettify_plot


        figure(4)
        copyobj(findobj(1, 'Type', 'axes'), 4);
        title('linearlyMovingBall')
        xlim([lims(2)/2 lims(2)])

        prettify_plot

        figure(5)
        copyobj(findobj(1, 'Type', 'axes'), 5);
        title('rectGrid')
        xlim([lims(2)/2 lims(2)])
        
        prettify_plot

        figure(6)
        copyobj(findobj(1, 'Type', 'axes'), 6);
        title('StaticDriftingGrating')
        xlim([lims(2)/2 lims(2)])

        prettify_plot

        figure(7)
        copyobj(findobj(1, 'Type', 'axes'), 7);
        title('rectGrid Pos Y')
        xlim([lims(2)/2 lims(2)])

        prettify_plot

        figure(8)
        copyobj(findobj(1, 'Type', 'axes'), 8);
        title('rectGrid Pos X')
        xlim([lims(2)/2 lims(2)])

        prettify_plot

    end


itN =1;

RandNeurons =  cell(5,4);

sumEneurons =  cell(1,3);

randNselect = 0;

FFFrand = 3; %randi(4,1,1);

SGrand =  16; %randi(18,1,1);

MBrand =  17; %randi(18,1,1);

RGrand = 14; % randi(18,1,1);

CoorValues = cell(3,18);


exampleNeurons = {"PV67_experiment_5_7_23", "PV103_Experiment_12_6_23","PV27_Experiment_25_6_23";...
    1,3,1;
    120,11,108};


%%%%%%%%% Iterate through experiment, insertion, stimulus %%%%%%%%%%

 ResponsiveU =0;
 Totalunits = 0;


for a = 1:length(expA)
    RandNeurons =  cell(5,4);


    for in = 1:length(animalA{a})

        if a == 3
                
            path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
            (in+1)+string(filesep)+"catgt_"+string(expA{a})+"_"+(in+1)+"_g0");
        else

            path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
                in+string(filesep)+"catgt_"+string(expA{a})+"_"+in+"_g0");
        end

        NP = NPAPRecording(path);

        cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

        %Good units

        GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

        ShankDist = zeros(1, length(GoodU_or)); %distance of unit from surface along the shank
        verticalDepth = zeros(1, length(GoodU_or)); %depth of unit along vertical axis
        XDist = zeros(1, length(GoodU_or)); %X distance of unit from insertion

        InserDepth = animalD{a}(in);

        AngleInser = animalA{a}(in);


         for u = 1:length(GoodU_or)

             %Main unit channel position

             %change channel 0 to 1

             cluster_info.ch(cluster_info.ch ==0) = 1;

             ch = cluster_info.ch(cluster_info.cluster_id == GoodU_or(u));

             ShankDist(u) = InserDepth - (NP.chLayoutPositions(2,ch));

             if ShankDist(u) < 1
                 
                 ShankDist(u) = 1;

             end

             verticalDepth(u) = sin(deg2rad(AngleInser))*ShankDist(u); %depth of unit along vertical axis

             XDist(u) = cos(deg2rad(AngleInser))*ShankDist(u);


         end


         unitZ = animals{a}(in,3)-verticalDepth/100;
         randPlus = (rand(1,length(GoodU_or))-0.5)*1;
         unitY = (animals{a}(in,2)+XDist/100)+randPlus;
         unitX = animals{a}(in,1)+randPlus; %randomize so we can see several units

         unitXo =  repmat(animals{a}(in,1),1,length(unitY));
         unitYo = (animals{a}(in,2)+XDist/100);

         %Up down.
            

         if max(unitZ) > animals{a}(in,3)
             disp('outside neurons up')
             return
         end

         if min(unitZ) < -InserDepth/100+animals{a}(in,3)
             disp('outside neurons down')
             return
         end

         figure(2)
         
         scatter3(unitX,unitY,unitZ, ...
             'bo', 'filled','MarkerFaceColor',colorA{a}, 'LineWidth',2, 'SizeData',25);

         if stimPlot == 1

             thres = 2; %Z-score threshold 
             
             unitTrans = MinusRBt{2,a}{in}; %zscore position 2
             unitTune = MinusRBt{3,a}{in}; %tunning position 3

            

             CoorValues{1,itN} = unitTrans;
             CoorValues{2,itN} = [unitXo;unitYo;unitZ];
             CoorValues{3,itN} = sprintf('%s-%s',expA{a},in);

             %colorStim = {'fullFieldFlash','linearlyMovingBall','rectGrid','StaticDriftingGrating';...
             %[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],[0.4660 0.6740 0.1880]};
             %
             %                  if string(expA{a}) == string(exampleNeurons(1,3)) & in == cell2mat(exampleNeurons(2,3)) %PV27
             %
             %                      2+2
             %
             %                  end
             %
             %                   find(GoodU_or == 454)

             %Example neurons
             if a == 1 & in == 1

                 eNeuron = 116;

             elseif a == 2 & in == 1

                 eNeuron = 10;

             else a == 3 & in ==1

                 eNeuron = 108;

             end

             SizeN = 180;

             Totalunits = Totalunits+length(cell2mat(unitTrans(1,1)));

             [r, stimsIn] = size(unitTrans);

             vals = reshape(cell2mat(unitTrans(1,:)), length(cell2mat(unitTrans(1,1))), stimsIn);

             [rows, Cols] = find(vals>thres);

             ResponsiveU = ResponsiveU + length(unique(rows));

            % ResponsiveU = ResponsiveU + length(abs(vals(vals>thres)));

             %%%


             for s = 1:length(MinusRBt{2,a}{in})

          

                 stimName  = StimOrder{stimTTL{a}(2,s)};

                 indexStim = find(strcmp(StimOrder, stimName)==1);

                 %responseTOcolor = floor(normalize(cell2mat(unitTrans(1,s)),'range',[1 num_colors]));

                 %responseTOsize = round(normalize(cell2mat(unitTrans(1,s)),'range',[10 200]));

                 vals = cell2mat(unitTrans(1,s));

                 valsT = cell2mat(unitTune(1,s));

                 %vals(vals>0) = log10(vals(vals>0)+1);

                 %vals(vals<0) = -log10(abs(vals(vals<0)-1));

                 %valsN = (vals-mean(valStim))./std(valStim);
% 
                 maxval = 10; %0.5;
                 minval = -10; %-0.5;


                 vals(vals>maxval) =maxval;

                 vals(vals<minval) =minval;

                 indexR = abs(vals)>thres; %zscore



                 %%%% If RG or MV plot with special color code

                 if string(stimName) == "linearlyMovingBall"


                     % Calculate colormaps

                     [CM,CMP,h, cmapR_theta]=colormap2D('plotPolarExample',1); %create colormap

                     [nRC,~,nTC]=size(cmapR_theta);

                     rhoR = vals;

                     thetaR = valsT+deg2rad(90); %sum 90 degrees because 0 is north (90)


                     rhoR(rhoR == 0) = 0.001;

                     thetaR(thetaR == 0) = 0.001;

                     divFactor = max(unique(thetaR));


                    % thetaR(thetaR<0) = abs(thetaR(thetaR<0))+pi;

                     maxR =10;
                     rhoR(rhoR>maxR) = maxR;

%                      hf = figure('Position',[400*in 400 400 450]);

                     NumU =length(thetaR);


                     figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:)))) %Plot on specific figure


                     for i = 1:length(indexR)

                         if indexR(i) == 1

                             colorsN = cmapR_theta( ceil((thetaR(i)/(divFactor))*nRC), : , ceil(rhoR(i)/maxR*nTC));

                             scatter3(unitX(i),unitY(i),unitZ(i), ...
                                 'filled','CData',colorsN...
                                 ,'MarkerEdgeColor','k', 'LineWidth',0.5, 'SizeData',40)
                         else

                             scatter3(unitX(i),unitY(i),unitZ(i), ...
                                 'filled', 'MarkerFaceAlpha', 0 ,'MarkerEdgeColor','k', 'LineWidth',0.25, 'MarkerEdgeAlpha', 0.2,'SizeData',40)
                         end
                     end


                     %Plot example neuron bigger

                     if (a == 1 && in ==1) || (a == 2 && in ==1) || (a == 3 && in ==1)

                         xR = unitX(eNeuron); yR = unitY(eNeuron); zR = unitZ(eNeuron);

                         if abs(vals(eNeuron)) > 2

                             respC = cmapR_theta( ceil((thetaR(eNeuron)/(divFactor))*nRC), : , ceil(rhoR(eNeuron)/maxR*nTC));

                             figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                             scatter3(xR,yR,zR, ...
                                 'filled', 'MarkerFaceAlpha', 1 , 'MarkerFaceColor',respC, 'MarkerEdgeColor',colorE, ...
                                 'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)

                         else
                             figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                             scatter3(xR,yR,zR, ...
                                 'filled' ,'MarkerFaceColor','w','MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerEdgeAlpha', 0.8,'SizeData',SizeN)
                         end
                     end
                     

                 elseif string(stimName) == "rectGrid"

                     if string(expA{a}) == "PV67_experiment_5_7_23"
                         maxPos = 10; %grids for rect Grid
                        
                     else
                         maxPos = 4;     
                     end

                     C = colormap(jet(30));


                     %%%Plot Y
                     v = valsT(:,2);

                     ColorVars = linspace(0,1,length(C));

                     diff_S = abs(bsxfun(@minus, v, ColorVars));

                     % Find the index of the smallest difference for each element of Y
                     [~, idx] = min(diff_S, [], 2);

                     
               
                     figure(7) %Plot on specific figure


                     scatter3(unitX(indexR),unitY(indexR),unitZ(indexR), ...
                         'filled','CData',C(idx(indexR),:)...
                         ,'MarkerEdgeColor','k', 'LineWidth',0.5, 'SizeData',40)

                     figure(7)
                     scatter3(unitX(~indexR),unitY(~indexR),unitZ(~indexR), ...
                         'filled', 'MarkerFaceAlpha', 0 ,'MarkerEdgeColor','k', 'LineWidth',0.25, 'MarkerEdgeAlpha', 0.2,'SizeData',40)


                     colormap(C)
                     c = colorbar;
                     c.Title.String = "Y Pos";
                     caxis([0 1])
                     set(c, 'YDir', 'reverse');

                     %%%Plot X
                     v = valsT(:,1);

                     diff_S = abs(bsxfun(@minus, v, ColorVars));

                     % Find the index of the smallest difference for each element of Y
                     [~, idx] = min(diff_S, [], 2);
                     

                     figure(8) %Plot on specific figure
                     scatter3(unitX(indexR),unitY(indexR),unitZ(indexR), ...
                         'filled','CData',C(idx(indexR),:)...
                         ,'MarkerEdgeColor','k', 'LineWidth',0.5, 'SizeData',40)


                     figure(8)
                     scatter3(unitX(~indexR),unitY(~indexR),unitZ(~indexR), ...
                         'filled', 'MarkerFaceAlpha', 0 ,'MarkerEdgeColor','k', 'LineWidth',0.25, 'MarkerEdgeAlpha', 0.2,'SizeData',40)
                    
                     colormap(C)
                     c = colorbar;
                     c.Title.String = "X Pos";
                     caxis([0 1])


                     %Plot example neuron bigger
                     if (a == 1 && in ==1) || (a == 2 && in ==1) || (a == 3 && in ==1)

                         xR = unitX(eNeuron); yR = unitY(eNeuron); zR = unitZ(eNeuron);

                         if abs(vals(eNeuron)) > 2

                             respC = C(idx(eNeuron),:);

                             figure(7)
                             scatter3(xR,yR,zR, ...
                                 'filled', 'MarkerFaceAlpha', 1 , 'MarkerEdgeColor',respC, 'MarkerEdgeColor',colorE, ...
                                 'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)

                             figure(8)
                             scatter3(xR,yR,zR, ...
                                 'filled', 'MarkerFaceAlpha', 1 , 'MarkerEdgeColor',respC, 'MarkerEdgeColor',colorE, ...
                                 'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)

                         else
                             figure(7)
                             scatter3(xR,yR,zR, ...
                                 'filled' ,'MarkerFaceColor','w','MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerEdgeAlpha', 0.8,'SizeData',SizeN)
                             
                             figure(8)
                             scatter3(xR,yR,zR, ...
                                 'filled' ,'MarkerFaceColor','w','MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerEdgeAlpha', 0.8,'SizeData',SizeN)
                         end
                     end



                 else %%% FFF & SDG
                    
                     ColorVars = linspace(-maxval,maxval,64);

                     diff_S = abs(bsxfun(@minus, vals', ColorVars));

                     % Find the index of the smallest difference for each element of Y
                     [~, idx] = min(diff_S, [], 2);

                     figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:)))) %Plot on specific figure

                     scatter3(unitX(indexR),unitY(indexR),unitZ(indexR), ...
                         'filled','CData',J(idx(indexR),:)...
                         ,'MarkerEdgeColor','k', 'LineWidth',0.5, 'SizeData',40)

                     figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                     scatter3(unitX(~indexR),unitY(~indexR),unitZ(~indexR), ...
                         'filled', 'MarkerFaceAlpha', 0 ,'MarkerEdgeColor','k', 'LineWidth',0.25, 'MarkerEdgeAlpha', 0.2,'SizeData',40)

                     colormap(J)
                     %caxis([0])
                     c = colorbar;
                     c.TickLabels = linspace(-maxval,maxval,11);
                     c.Title.String = "Z-score";

                     %Plot example neuron bigger
                     if (a == 1 && in ==1) || (a == 2 && in ==1) || (a == 3 && in ==1)

                         xR = unitX(eNeuron); yR = unitY(eNeuron); zR = unitZ(eNeuron);

                         if abs(vals(eNeuron)) > 2

                             respC = J(idx(eNeuron));

                             figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                             scatter3(xR,yR,zR, ...
                                 'filled', 'MarkerFaceAlpha', 1 , 'MarkerEdgeColor',respC, 'MarkerEdgeColor',colorE, ...
                                 'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)

                         else
                             figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                             scatter3(xR,yR,zR, ...
                                 'filled' ,'MarkerFaceColor','w','MarkerEdgeColor','k', 'LineWidth',0.5, 'MarkerEdgeAlpha', 0.8,'SizeData',SizeN)
                         end
                     end

                 end



                 %                  if string(expA{a}) == string(exampleNeurons(1,2)) & in == cell2mat(exampleNeurons(2,2))
%   
%                      eNeuron = 120;
% 
%                      eN = GoodU_or(eNeuron);
% 
%                      xR = unitX(eNeuron); yR = unitY(eNeuron); zR = unitZ(eNeuron);
% 
% 
%                      if string(stimName) == "fullFieldFlash"
%                          c = 1;
%                          RandNeurons{1,c} = expA{a};
%                          RandNeurons{2,c} = in;
%                          RandNeurons{3,c} = eN;
%                          RandNeurons{4,c} = [vals(eNeuron)];
%                          RandNeurons{5,c} = stimName;
%                         
%                          figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
%                          scatter3(xR,yR,zR, ...
%                              'filled', 'MarkerFaceAlpha', 1 , 'CData',[respC; ], 'MarkerEdgeColor',colorE, ...
%                              'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)
% 
% 
%                      end
%                     
%                      if string(stimName) == "StaticDriftingGrating"
%                          c = 2;
%                          RandNeurons{1,c} = expA{a};
%                          RandNeurons{2,c} = in;
%                          RandNeurons{3,c} = eN;
%                          RandNeurons{4,c} = [vals(eNeuron)];
%                          RandNeurons{5,c} = stimName;
% 
%                          figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
%                          scatter3(xR,yR,zR, ...
%                              'filled', 'MarkerFaceAlpha', 1 , 'CData',[respC; 1 1 1], 'MarkerEdgeColor',colorE, ...
%                              'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)
% 
%                      end
% 
%                      if string(stimName) == "linearlyMovingBall"
% 
%                          figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
%                          scatter3(xR,yR,zR, ...
%                              'filled', 'MarkerFaceAlpha', 1 , 'CData',[respC; 1 1 1], 'MarkerEdgeColor',colorE, ...
%                              'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)
% 
%                      end
% 
%                      if string(stimName) == "rectGrid"
% 
%                          figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
%                          scatter3(xR,yR,zR, ...
%                              'filled', 'MarkerFaceAlpha', 1 , 'CData',[respC; 1 1 1], 'MarkerEdgeColor',colorE, ...
%                              'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)
% 
%                      end
% 
% 
% 
% 
% 
%                  end
% % 

                     if randNselect == 1

                         if length(stimName) == length('fullFieldFlash') &  ~isempty(FFFrand(FFFrand == itN)) %Check if stimulus if FFF and if the insertion is the random one for the stim.
                             c = 1;
                             RandNeurons{1,c} = expA{a};
                             RandNeurons{2,c} = in;

                             %rng(0, 'twister');

                             if ~isempty(find(indexR ==0))
                                 nr = find(indexR ==0); nrI = nr(randi([1, length(nr)])); randNR = vals(nrI); xNR = unitX(nrI); yNR = unitY(nrI); zNR = unitZ(nrI);
                                 Unr = GoodU_or(nrI);

                             else
                                 randNR = NaN;
                                 Unr = NaN;
                                 xNR = NaN;
                                 yNR = NaN;
                             end

                             %rng(1, 'twister');

                             if ~isempty(find(indexR ==1))
                                 r = find(indexR ==1); rI = r(randi([1, length(r)])); randR = vals(rI); xR = unitX(rI); yR = unitY(rI); zR = unitZ(rI);
                                 Ur = GoodU_or(rI);
                                 respC = J(indices(rI),:);
                             else
                                 randR = NaN;
                                 Ur = NaN;
                                 xR = NaN;
                                 yR = NaN;
                                 respC = [1 1 1];

                             end

                             RandNeurons{3,c} = [Ur Unr];
                             RandNeurons{4,c} = [randR randNR];

                             figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                             scatter3([xR;xNR],[yR;yNR],[zR;zNR], ...
                                 'filled', 'MarkerFaceAlpha', 1 , 'CData',[respC; 1 1 1], 'MarkerEdgeColor',colorE, 'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)
                             % text(Xlims(2),yR+1,zR, sprintf('%d-%d-r',itN,rI),"FontWeight","bold","Color",colorT,"FontSize",15)
                             % text(Xlims(2),yNR+1,zNR, sprintf('%d-%d-nr',itN,GoodU_or(nrI)),"FontWeight","bold","Color",colorT,"FontSize",15)
                         end


                         if  length(stimName) == length('rectGrid') & ~isempty(RGrand(RGrand == itN))
                             c = 2;
                             RandNeurons{1,c} = expA{a};
                             RandNeurons{2,c} = in;

                             rng(2, 'twister');

                             if ~isempty(find(indexR ==0))
                                 nr = find(indexR ==0); nrI = nr(randi([1, length(nr)])); randNR = vals(nrI); xNR = unitX(nrI); yNR = unitY(nrI); zNR = unitZ(nrI);
                                 Unr = GoodU_or(nrI);
                             else
                                 randNR = NaN;
                                 Unr = NaN;
                                 xNR = NaN;
                                 yNR = NaN;
                             end

                             rng(3, 'twister');

                             if ~isempty(find(indexR ==1))
                                 r = find(indexR ==1); rI = r(randi([1, length(r)])); randR = vals(rI); xR = unitX(rI); yR = unitY(rI); zR = unitZ(rI);
                                 Ur = GoodU_or(rI);
                                 respC = J(indices(rI),:);
                             else
                                 randR = NaN;
                                 Ur = NaN;
                                 xR = NaN;
                                 yR = NaN;
                                 respC = [1 1 1];
                             end

                             RandNeurons{3,c} = [Ur Unr];
                             RandNeurons{4,c} = [randR randNR];

                             figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                             scatter3([xR;xNR],[yR;yNR],[zR;zNR], ...
                                 'filled', 'MarkerFaceAlpha', 1 , 'CData',[respC; 1 1 1], 'MarkerEdgeColor',colorE, 'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)
                             % text(Xlims(2),yR+1,zR, sprintf('%d-%d-r',itN,rI),"FontWeight","bold","Color",colorT,"FontSize",15)
                             % text(Xlims(2),yNR+1,zNR, sprintf('%d-%d-nr',itN,GoodU_or(nrI)),"FontWeight","bold","Color",colorT,"FontSize",15)
                         end

                         if  length(stimName) == length('linearlyMovingBall')  & ~isempty(MBrand(MBrand == itN))
                             c =3;
                             RandNeurons{1,c} = expA{a};
                             RandNeurons{2,c} = in;

                             rng(4, 'twister');

                             if ~isempty(find(indexR ==0))
                                 nr = find(indexR ==0); nrI = nr(randi([1, length(nr)])); randNR = vals(nrI); xNR = unitX(nrI); yNR = unitY(nrI); zNR = unitZ(nrI);
                                 Unr = GoodU_or(nrI);
                             else
                                 randNR = NaN;
                                 Unr = NaN;
                                 xNR = NaN;
                                 yNR = NaN;
                             end

                             rng(5, 'twister');

                             if ~isempty(find(indexR ==1))
                                 r = find(indexR ==1); rI = r(randi([1, length(r)])); randR = vals(rI); xR = unitX(rI); yR = unitY(rI); zR = unitZ(rI);
                                 Ur = GoodU_or(rI);
                                 respC = J(indices(rI),:);
                             else
                                 randR = NaN;
                                 Ur = NaN;
                                 xR = NaN;
                                 yR = NaN;
                                 respC = [1 1 1];
                             end

                             RandNeurons{3,c} = [Ur Unr];
                             RandNeurons{4,c} = [randR randNR];

                             figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                             scatter3([xR;xNR],[yR;yNR],[zR;zNR], ...
                                 'filled', 'MarkerFaceAlpha', 1 , 'CData',[respC; 1 1 1], 'MarkerEdgeColor',colorE, 'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)
                             %  text(Xlims(2),yR+1,zR, sprintf('%d-%d-r',itN,rI),"FontWeight","bold","Color",colorT,"FontSize",15)
                             %  text(Xlims(2),yNR+1,zNR, sprintf('%d-%d-nr',itN,GoodU_or(nrI)),"FontWeight","bold","Color",colorT,"FontSize",15)
                         end

                         if length(stimName) == length('StaticDriftingGrating')  & ~isempty(SGrand(SGrand == itN))
                             c=4;
                             RandNeurons{1,c} = expA{a};
                             RandNeurons{2,c} = in;

                             rng(6, 'twister');

                             if ~isempty(find(indexR ==0))
                                 nr = find(indexR ==0); nrI = nr(randi([1, length(nr)])); randNR = vals(nrI); xNR = unitX(nrI); yNR = unitY(nrI); zNR = unitZ(nrI);
                                 Unr = GoodU_or(nrI);
                             else
                                 randNR = NaN;
                                 Unr = NaN;
                                 xNR = NaN;
                                 yNR = NaN;
                             end

                             rng(7, 'twister');

                             if ~isempty(find(indexR ==1))
                                 r = find(indexR ==1); rI = r(randi([1, length(r)])); randR = vals(rI); xR = unitX(rI); yR = unitY(rI); zR = unitZ(rI);
                                 Ur = GoodU_or(rI);
                                 respC = J(indices(rI),:);
                             else
                                 randR = NaN;
                                 Ur = NaN;
                                 xR = NaN;
                                 yR = NaN;
                                 respC = [1 1 1];
                             end

                             RandNeurons{3,c} = [Ur Unr];
                             RandNeurons{4,c} = [randR randNR];

                             figure(stimfigNumbers(strcmp(unitTrans(2,s),StimOrder(1,:))))
                             scatter3([xR;xNR],[yR;yNR],[zR;zNR], ...
                                 'filled', 'MarkerFaceAlpha', 1 , 'CData',[respC; 1 1 1], 'MarkerEdgeColor',colorE, 'LineWidth',2, 'MarkerEdgeAlpha', 1,'SizeData',SizeN)
                             % text(Xlims(2),yR+1,zR, sprintf('%d-%d-r',itN,rI),"FontWeight","bold","Color",colorT,"FontSize",15)
                             % text(Xlims(2),yNR+1,zNR, sprintf('%d-%d-nr',itN,GoodU_or(nrI)),"FontWeight","bold","Color",colorT,"FontSize",15)
                         end

                     end

                     %prettify_plot


             end

         end

         itN = itN+1;

    end

end
%%



%%
writerObj = VideoWriter('movBall.mp4','MPEG-4');
writerObj.FrameRate=30;
%writerObj.Quality=videoQuality;
open(writerObj);

%F=figure('position',[50 50 550 500],'color','w');h=axes;
%initial plot
figure(4)
F=gcf;
h=gca;
set(h,'nextplot','replacechildren');
set(F,'Renderer','zbuffer');
[V1,V2]=view;
view(h,[72,V2]);
axis off;
title('');
rotations = -5:0.5:360;
view(h,[rotations(1),V2]);
axis vis3d
for i=1:numel(rotations)
    view(h,[rotations(i),V2]);


    frame = getframe(F);
    writeVideo(writerObj,frame);

end
close(writerObj);


%%


%PV67
refIns = [44 34.6 Zlims(2)];

scatter3(refIns(1),refIns(2),refIns(3), 'bo', 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1], 'LineWidth',2, 'SizeData',sizeDot);
text(refIns(1),refIns(2)-1,refIns(3), '1', 'Color', 'blue','FontSize',sizeText);


%Y = X (fig) & -X = Y (fig)

Ins2 = [refIns(1)-853/100 refIns(2) refIns(3)];
scatter3(Ins2(1),Ins2(2),Ins2(3), 'bo', 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1], 'LineWidth',2, 'SizeData',sizeDot);
text(Ins2(1)+offsetTextY,Ins2(2)-offsetTextX,Ins2(3), '2', 'Color', 'blue','FontSize',sizeText);

Ins3 = [refIns(1)-938/100 refIns(2)+994/100 refIns(3)];
scatter3(Ins3(1),Ins3(2),Ins3(3), 'bo', 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1], 'LineWidth',2, 'SizeData',sizeDot);
text(Ins3(1)++offsetTextY,Ins3(2)-offsetTextX,Ins3(3), '3', 'Color', 'blue','FontSize',sizeText);

Ins4 = [refIns(1) refIns(2)+994/100  refIns(3)];
scatter3(Ins4(1),Ins4(2),Ins4(3), 'bo', 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 0 1], 'LineWidth',2, 'SizeData',sizeDot);
text(Ins4(1)+offsetTextY,Ins4(2)-offsetTextX,Ins4(3), '4', 'Color', 'blue','FontSize',sizeText);

prettify_plot

%PV102 %Not rand

refIns = [40.6 33.4 Zlims(2)];

Ins1 = [refIns(1) refIns(2) refIns(3)];
scatter3(refIns(1),refIns(2),refIns(3), 'bo', 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0], 'LineWidth',2, 'SizeData',sizeDot);
text(refIns(1)+offsetTextY,refIns(2)-offsetTextX,refIns(3), '1', 'Color', 'red','FontSize',sizeText);

Ins2 = [refIns(1)-240/100+415/100 refIns(2)+778/100 refIns(3)];
scatter3(Ins2(1),Ins2(2),Ins2(3), 'bo', 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0], 'LineWidth',2, 'SizeData',sizeDot);
text(Ins2(1)+offsetTextY,Ins2(2)-offsetTextX,Ins2(3), '2', 'Color', 'red','FontSize',sizeText);



prettify_plot



angle_degrees = 88; % Specify the desired angle in degrees

% Convert the angle to radians
angle_radians = deg2rad(angle_degrees);

% Define the starting and ending points for the lines
start_point = [0, 0, 0];
end_point = [cos(angle_radians), sin(angle_radians), 0];

% Plot the lines
plot3([start_point(1), end_point(1)], [start_point(2), end_point(2)], [start_point(3), end_point(3)], 'r', 'LineWidth', 2);

[vertices, faces, ~, ~] = stlread(path_brainO);




% convertorZ =
% convertorY =
% convertorX =

yticks([0:90])

%patch(brain,'FaceColor', 'none', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

%%%%%%%%%% FUNCTIONS %%%%%%%%%%%

%%%%%%% Micromanipulator coordinates conversion to matlab coordinates:
function [Xcoor,Ycoor]= convMan2coor(Zman,Xman)
Xcoor = Zman/100;
Ycoor = -Xman/100;
end


%%%%%%% Find max X (surface of brain) given x and y

function [closest_z]= surfaceZ(x_target,y_target,brain)

vertices = brain.Points;

% Calculate the squared distance between vertices and the target point
distances_squared = (vertices(:, 1) - x_target).^2 + (vertices(:, 2) - y_target).^2;

% Find vertices with distances within a small tolerance
tolerance = 1;  % Adjust the tolerance as needed
close_indices = find(distances_squared < tolerance);

if isempty(close_indices)
    % No vertices are very close
    [~, min_index] = min(distances_squared);
else
    % Among close vertices, find the one with the highest Z coordinate
    [~, max_z_index] = max(vertices(close_indices, 3));
    min_index = close_indices(max_z_index);
end

% Get the Z coordinate of the closest vertex
closest_z = vertices(min_index, 3);



end

