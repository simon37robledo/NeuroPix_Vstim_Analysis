NP = NPAPRecording('C:\Users\MarkS9\Documents\Spike_GLX_TempDATA\test_Diode_4_31_1_24_g0\test_Diode_4_31_1_24_g0_imec0');

A = NP.getAnalogData(1,0,NP.recordingDuration_ms);

plot(squeeze(A))