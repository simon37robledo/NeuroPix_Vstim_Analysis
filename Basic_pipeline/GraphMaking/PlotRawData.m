

function [fig mx mn] = PlotRawData(fig,NP,chan,startTimes,window,freq,type,preBase,spikeTimes)
tic
raw_signal = squeeze(NP.getData(chan,startTimes,window));
toc
%raw_signal = squeeze(aaa(:,1,:));

%aaa = raw_signal;

fc = 300;

%%
% cd('D:\Mark_S13\Downloads')
% raw_signale = readNPY('Ch16_PV27_25_6_23_3.npy');
% raw_signale = squeeze(raw_signal);

if freq == "AP"
    [b, a] = butter(4, fc/(NP.samplingFrequency/2), 'high');
end

t = (0:length(raw_signal)-1) /NP.samplingFrequency;
%offset = 2.5*mean(std(raw_signal));
offset = mean(2*std(raw_signal));

figure(fig)

if type == "line"

    mx = -inf;
    mn = inf;

for tr = 1:min(size(raw_signal))

    if freq == "AP"

        if size(raw_signal,2) == 1
            raw_signal =  raw_signal'; 
        end
        signal = filtfilt(b, a, raw_signal(tr,:));

    else
        if size(raw_signal,2) == 1
            raw_signal =  raw_signal'; 
        end
        signal = raw_signal(tr,:);
    end

    j = size(raw_signal,1)-tr;
    plot(t, signal+offset*(j), 'LineWidth',0.5,'Color','k');

    if size(spikeTimes,2) > 0

        %Convert ms of raster to seconds
        spikeTimesIndex = find(spikeTimes(tr,:)>0)/1000;

        hold on
        %plot(spikeTimesIndex,repmat(min(signal-offset*(j)),1,length(spikeTimesIndex)),'.','Color','b','MarkerSize',7);

        %%Improve with plot function for several trials
        try
        xline(spikeTimesIndex,'LineWidth',1,'Color','b','Alpha',0.3) %Plot spikes.
        catch
        fprintf('Selected trial has no spikes')
        end
        
    end
    hold on

 
    if max(signal) > mx
        mx = max(signal);
    end

    if min(signal) < mn
        mn = min(signal);
    end
    hold on
end



%yticks([])

xlim([0 length(raw_signal)/NP.samplingFrequency]);
lims = xlim;

limsY = ylim;
ylim([limsY(1)-std(raw_signal)/2 limsY(2)+std(raw_signal)/2]);

stimDur = window-2*preBase;

xticks([0 preBase/1000:300/1000:(stimDur+preBase*2)/1000 lims(end)-0.00001])
xticklabels([-(preBase) 0:300:round((stimDur/100))*100 round((stimDur/100))*100 + preBase]./1000)

ylabel('10 Trials')
xlabel('Time (sec)')

ax = gca; % Get current axes
ax.YAxis.FontSize = 7; % Change font size of y-axis tick labels

%fig.Position = [  1027         246         415         241];
hold off

end
% cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
% print(gcf, sprintf('%s-10TRRawData-MovBall-%s-U%d-W%d-%dW-speed-500.pdf',NP.recordingName,orderNames{k},u,window_size(1),window_size(2)), '-dpdf', '-r300', '-vector');

end