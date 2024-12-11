%%Test window

           A = [stimOn directions' offsets' sizes' speeds'];


        eNeuron = 75; %8
        %eNeuron = 1;

        orderS = [2 3 4 5;4 2 3 5;5 2 3 4;3 2 4 5];
        orderNames = {'dir_off_sizes_speeds';'sizes_dir_off_speeds';'speeds_dir_off_sizes';'off_dir_sizes_speeds'};


        cd(NP.recordingDir)
        NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;

        posX = squeeze(NeuronVals(:,:,3));
        posY = squeeze(NeuronVals(:,:,2));

        uDir = unique(directions);
        bin = 1;
        bin2 =50;
        trialsPerAngle = trialDivision*offsetN*speedN*sizeN*orientN;
        for k = 1

            [C sIndex2]= sortrows(A,orderS(k,:));

            %Sort directions:

            directimesSorted = C(:,1)';

            [Mr] = BuildBurstMatrix(goodU,round(p.t/bin2),round((directimesSorted-preBase)/bin2),round((stimDur+preBase*2)/bin2));

            Mr2 = [];

            for u = eNeuron

                j=1;

                mergeTrials = 1;

                if mergeTrials ~= 1 %Merge trials

                    for i = 1:mergeTrials:trials

                        meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end),u,:)),1);

                        Mr2(j,:) = meanb;

                        j = j+1;

                    end
                else
                    Mr2 = Mr(:,u,:);
                end

                [nT,nB] = size(Mr2);

                fig = figure;

                imagesc(squeeze(Mr2));colormap(flipud(gray(64)));
                %Plot stim start:
                xline(preBase/bin2,'k', LineWidth=1.5)
                %Plot stim end:
                xline(stimDur/bin2+preBase/bin2,'k',LineWidth=1.5)
                ylabel('Trials');xlabel('Time (ms)');
                title(sprintf('U.%d-Unit-phy-%d',u,GoodU_or(u)));

                %xticks([0.5 (preBase/bin):10:nB])
                %xticklabels([-preBase 0:10*bin:nB*bin])
                yticklabels([yticks]*mergeTrials)
                dirStart = C(1,2);
                offStart = C(1,3);
                for t = 1:nT
                    if dirStart ~= C(t,2)
                        yline(t,'r',string(C(t,2)),'LabelHorizontalAlignment', 'left','LabelVerticalAlignment', 'bottom',LineWidth=3);
                        dirStart = C(t,2);
                    end
                    if offStart ~= C(t,3)
                        yline(t,'b',string(C(t,3)),'LabelHorizontalAlignment', 'left','LabelVerticalAlignment', 'bottom',LineWidth=2);
                        offStart = C(t,3);
                    end
                               
                end
%                 %Directions
%                 v = nT/direcN:nT/direcN:nT-1;
%                 yline(v+0.5,'r', LineWidth=3);
%                 %Offsets
%                 v = nT/(direcN*offsetN):nT/(direcN*offsetN):nT-1;
%                 yline(v+0.5,'b', LineWidth=2);
%                 %sizes
%                 v = nT/(direcN*offsetN*sizeN):nT/(direcN*offsetN*sizeN):nT-1;
%                 yline(v+0.5, LineWidth=0.5);

                

                %                             hcb = colorbar();
                %                             title(hcb,'Spikes/sec');
                %caxis([0 max(0.2,max(max_mean_value(u)))])
                hold on
                %Plot rectangle:
                for d = 1:size(NeuronVals,2)
                rectangle('Position', [posX(u,d)/(bin2)+round(preBase/bin2), (posY(u,d)*trialDivision-trialDivision)/mergeTrials, window_size(2)/(bin2), trialDivision],...
                    'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
                end
                hold off
                %prettify_plot
                
                yyaxis right
                ylim([1,nT])
                yticks([trialsPerAngle:trialsPerAngle:nT])
                ax = gca;
                ax.YAxis(2).Color = 'k';
                yticklabels(sort(uDir,'descend'))
                ylabel('Angles')
                lims =xlim;
                xt = xticks;

                cd(NP.recordingDir)
                if ~exist(path+"\Figs",'dir')
                    mkdir Figs
                end
                cd(NP.recordingDir + "\Figs")
                set(fig,'Color','w')

            end

        end