

function fig = PlotRawData(fig,NP,chan,startTimes,window,freq,type,preBase)
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
offset = 2.5*mean(std(raw_signal));

figure(fig)

if type == "line"

for tr = 1:size(raw_signal,1)

    if freq == "AP"

        signal = filtfilt(b, a, raw_signal(tr,:));

    else

        signal = raw_signal(tr,:);

    end

    j = size(raw_signal,1)-tr;
    plot(t, signal+offset*(j), 'LineWidth',0.5,'Color','k');
    hold on
end

yticks([])

xlim([0 length(raw_signal)/NP.samplingFrequency]);
lims = xlim;

limsY = ylim;
ylim([limsY(1)-10 limsY(2)-10]);

stimDur = window-2*preBase;

xticks([0 preBase/1000:300/1000:(stimDur+preBase*2)/1000 lims(end)-0.00001])
xticklabels([-(preBase) 0:300:round((stimDur/100))*100 round((stimDur/100))*100 + preBase]./1000)
xline([preBase/1000 (stimDur+preBase)/1000],LineWidth=1.5)
ylabel('10 Trials')
xlabel('Time (sec)')
%fig.Position = [  1027         246         415         241];
hold off

end
% cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
% print(gcf, sprintf('%s-10TRRawData-MovBall-%s-U%d-W%d-%dW-speed-500.pdf',NP.recordingName,orderNames{k},u,window_size(1),window_size(2)), '-dpdf', '-r300', '-vector');

end