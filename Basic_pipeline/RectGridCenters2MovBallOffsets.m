%%% Rect center vs linearly moving ball offsets
rect = load('linearlyMovingBall_2024_7_10_12_34_40_554');

%%%Sreen size --> see screen resolution in display settings or prop name = rect):
offsets = 9;
coorRect = cell2mat(rect.VSMetaData.allPropVal(find(strcmp(rect.VSMetaData.allPropName,'rect'))));

%%%%Set up rect Grid = Divide screen side by number of offsets
rectSizeTR1 = coorRect(4)/offsets;

%%%Set up MovBall = Making sure that ball is within screen limits
divisionoff = floor(offsets/2);

halfScreen = coorRect(4)/2;

rBall = rectSizeTR1/2;

offsetsPositive = (halfScreen-rBall)/divisionoff:(halfScreen-rBall)/divisionoff:halfScreen-rBall;

offsetsNegative = sort(offsetsPositive*-1);

offsetsVector = [offsetsNegative 0 offsetsPositive];


cd('C:\\Users\MarkS9\Documents\GitHub\visual-stimulation-gui\stats\')

%%%%%%%%%%


figure
hold on; plot(rect.VSMetaData.allPropVal{21,1}.X1{1,1}(1,1),rect.VSMetaData.allPropVal{21,1}.Y1{1,1}(1,1),'.',MarkerSize=10)%,...
hold on; plot(rect.VSMetaData.allPropVal{21,1}.X2{1,1}(1,1),rect.VSMetaData.allPropVal{21,1}.Y2{1,1}(1,1),'.',MarkerSize=10)
hold on; plot(rect.VSMetaData.allPropVal{21,1}.X3{1,1}(1,1),rect.VSMetaData.allPropVal{21,1}.Y3{1,1}(1,1),'.',MarkerSize=10)%,...
hold on; plot(rect.VSMetaData.allPropVal{21,1}.X4{1,1}(1,1),rect.VSMetaData.allPropVal{21,1}.Y4{1,1}(1,1),'.',MarkerSize=10)%,...

Xc = (rect.VSMetaData.allPropVal{21,1}.X2{1,1}(1,1)-rect.VSMetaData.allPropVal{21,1}.X1{1,1}(1,1))/2+rect.VSMetaData.allPropVal{21,1}.X1{1,1}(1,1);
Yc = (rect.VSMetaData.allPropVal{21,1}.Y4{1,1}(1,1)-rect.VSMetaData.allPropVal{21,1}.Y1{1,1}(1,1))/2+rect.VSMetaData.allPropVal{21,1}.Y1{1,1}(1,1);

hold on; plot(Xc,Yc,'.',MarkerSize=10)%,...
xlim([1 700])
ylim([1 700])

rectangle('Position',rect.VSMetaData.allPropVal{21,1}.X1{1,1}(1,3),rect.VSMetaData.allPropVal{21,1}.Y1{1,1}(1,3),)




plot([rect.VSMetaData.allPropVal{21,1}.X1{1,1}(1,3),rect.VSMetaData.allPropVal{21,1}.X2{1,1}(1,3)],[rect.VSMetaData.allPropVal{21,1}.Y1{1,1}(1,3),rect.VSMetaData.allPropVal{21,1}.Y2{1,1}(1,3)],'.','Color','g')

%5x5 grid offsets
center = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,3) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,3))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,3);

rect1 = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,1) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,1))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,1);
Offest1 = rect1-center;

rect2 = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,2) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,2))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,2);
Offest2 = rect2-center;

rect4 = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,4) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,4))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,4);
Offest4 = rect4-center;

rect5 = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,5) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,5))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,5);
Offest5 = rect5-center;

%%Tilling ratios
rectSide = VSMetaData.allPropVal{22,1}(1);

ballSizes = [280 200 120 40];

TR1 = ballSizes(1)/rectSide;
TR2 = ballSizes(2)/rectSide;
TR3 = ballSizes(3)/rectSide;
TR4 = ballSizes(4)/rectSide;


