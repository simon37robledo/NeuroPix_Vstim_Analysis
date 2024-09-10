%Static and drifting grating to movbal conversion for size and speed

pixelsperCM = 32.5423728813559;%Acer oled screen with 1920*1080 res

BallSize = 120; %pix (diameter)
SpatialFreq1 = 1/(BallSize/pixelsperCM);%cylces per Cm
SpatialFreq2 = 1/(2*BallSize/pixelsperCM);

Speed = (([125 250 500 1000])/pixelsperCM)*% cm per sec

TempFreq1 = ([SpatialFreq2*Speed(3) SpatialFreq1*Speed(4)]); %cycle/sec

TempFreq2 = [SpatialFreq2*Speed(1) SpatialFreq1*Speed(2) ];

SpeedTest = SpatialFreq2/TempFreq1 SpatialFreq2/TempFreq2  SpatialFreq1/TempFreq1 SpatialFreq1/TempFreq2

%v = (SpatialFreq1/TempFreq1);