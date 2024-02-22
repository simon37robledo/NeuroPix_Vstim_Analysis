%%% Rect center vs linearly moving ball offsets

RG = VSMetaData;

figure
plot([VSMetaData.allPropVal{21,1}.X1{1,1}(1,3),VSMetaData.allPropVal{21,1}.X2{1,1}(1,3),...
    VSMetaData.allPropVal{21,1}.X3{1,1}(1,3),VSMetaData.allPropVal{21,1}.X4{1,1}(1,3)],...
    [VSMetaData.allPropVal{21,1}.Y1{1,1}(1,3),VSMetaData.allPropVal{21,1}.Y2{1,1}(1,3),...
    VSMetaData.allPropVal{21,1}.Y3{1,1}(1,3),VSMetaData.allPropVal{21,1}.Y4{1,1}(1,3)],'.')

xlim([400,600])
ylim([50,250])

plot([VSMetaData.allPropVal{21,1}.X1{1,1}(1,3),VSMetaData.allPropVal{21,1}.X2{1,1}(1,3)],[VSMetaData.allPropVal{21,1}.Y1{1,1}(1,3),VSMetaData.allPropVal{21,1}.Y2{1,1}(1,3)],'.','Color','g')

%5x5 grid
center = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,3) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,3))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,3);

rect1 = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,1) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,1))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,1);
Offest1 = rect1-center;

rect2 = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,2) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,2))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,2);
Offest2 = rect2-center;

rect4 = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,4) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,4))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,4);
Offest4 = rect4-center;

rect5 = (VSMetaData.allPropVal{21,1}.X2{1,1}(1,5) - VSMetaData.allPropVal{21,1}.X1{1,1}(1,5))/2+VSMetaData.allPropVal{21,1}.X1{1,1}(1,5);
Offest5 = rect5-center;