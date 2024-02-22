%% Plot scatter plots of responsivity and coordinates

M = [];

for in=1:itN-1

% Assuming C is your cell array
V = CoorValues{1,in};
unitN = length(CoorValues{1,in}{1,1});
[rows stimN]= size(V);

% Assuming your cell array is named S
S=V(2,:);

% Create a mapping from strings to numbers
mapping = containers.Map({'fullFieldFlash','linearlyMovingBall','rectGrid','StaticDriftingGrating'}, 1:4);

% Create a new cell array with two rows
new_V = cell(2, length(S));

% Assign the first row of new_S to be the same as S
new_V(1,:) = V(1,:);

% Assign the second row of new_S using the mapping
new_V(2,:) = values(mapping, S);


% Initialize M
Mval = zeros(unitN*stimN,2);

% Fill M
for i = 1:stimN
    Mval((i-1)*unitN+1:i*unitN,1) = new_V{1,i}';
    Mval((i-1)*unitN+1:i*unitN,2) =  new_V{2,i}';
end

C = CoorValues{2,in}';
Mcoor = repmat(C,length(S),1);

M = [M;[Mval, Mcoor]]; 

in
end

%%
Mf = M;%(M(:,1) > 0.1,:);
Mf(:,1) = log10(Mf(:,1)+1);
%Plot boxplots per stimulus
figure()
h = boxplot(Mf(:,1), Mf(:,2), 'Notch', 'on','Labels',{'fullFieldFlash','linearlyMovingBall','rectGrid','StaticDriftingGrating'},...
    'Colors',[0.5 0.7 0.7]); %,'OutlierSize', 0);
set(h, 'LineWidth', 2);
figure()
labels = {'fullFieldFlash','linearlyMovingBall','rectGrid','StaticDriftingGrating'};

f = violinplot(Mf(:,1), labels(Mf(:,2)));
grid on
