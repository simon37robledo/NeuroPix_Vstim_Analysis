%%%%PLot several trials one channel: 

chan = 16;

startTimes = directimesSorted-preBase;
window = stimDur+preBase*2;

freq = "AP"; %or "LFP"

type = "line"; %or heatmap

%function RawChan = PlotRawData(NP,chan,startTimes,window,freq,type)

raw_signal = squeeze(NP.getData(chan,startTimes,window));

aaa = raw_signal;

fc = 300;

raw_signal = squeeze(aaa(:,1,:));

cd('D:\Mark_S13\Downloads')
raw_signal = readNPY('Ch16_PV27_25_6_23_3.npy');
raw_signal = squeeze(raw_signal(1,1,:));

if freq == "AP"
    [b, a] = butter(4, fc/(NP.samplingFrequency/2), 'high');

    % Apply the filter using zero-phase filtering (to avoid phase shift)
    signal  = filtfilt(b, a, raw_signal);

else
    signal = raw_signal;
end

% Plot the raw and filtered signals
t = (0:length(raw_signal)-1) /NP.samplingFrequency;

subplot(2,1,1);
plot(t, raw_signal);
title('Raw Neural Signal');
xlabel('Time (s)');

subplot(2,1,2);
plot(t, signal);
title('High-Pass Filtered Signal (Spikes Only)');
xlabel('Time (s)');

%end