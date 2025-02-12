p = 'C:\Users\MarkS9\Documents\Spike_GLX_TempDATA\testBouncingBall_1000bals_g2\testBouncingBall_1000bals_g2_imec0';

NP = NPAPRecording(p);

a = NP.getAnalogData(1,0,NP.recordingDuration_ms-100);

figure;plot(squeeze(a)');