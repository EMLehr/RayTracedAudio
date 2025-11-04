fs = 48000;
% src = mono audio file
room = stlread("D:\Matlab_Projects\RayTracedAudio\room.stl");
src = audioread("Clave_Mono.wav");
% dimensions in meters
% 2.6416 = y = length = 104
% 2.9972 = z = height = 118
% 2.8666 = x = width = 112.85 in




%Effective absorption for entire wall instead of every object on wall



%%
% % Rescale room size to be correct
% % Step 1: Extract points and faces
% points = room.Points;            
% faces  = room.ConnectivityList;  
% 
% % Step 2: Current dimensions
% minCoords = min(points);
% maxCoords = max(points);
% dims = maxCoords - minCoords;   % [X, Y, Z]
% disp('Current dimensions [X, Y, Z]:')
% disp(dims)
% 
% % Step 3: Desired dimensions
% % Keep axis order the same: X, Y, Z
% targetDims = [2.1082, 2.6416, 2.9972];  % [X, Y, Z]
% 
% % Step 4: Compute scaling factors per axis
% scaleFactors = targetDims ./ dims;
% 
% % Step 5: Apply scaling
% points = points .* scaleFactors;
% 
% % Step 6: Rebuild triangulation
% roomFinal = triangulation(faces, points);
% 
% % Step 7: Visualize transparent
% figure
% trisurf(roomFinal, ...
%         'FaceColor', 'cyan', ...
%         'FaceAlpha', 0.3, ...
%         'EdgeColor', 'k');
% axis equal
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% title('Rescaled Room STL (No Reordering)')
% view(3)
% camlight
% lighting gouraud
% 
% % Save STL using triangulation object directly
% stlwrite(roomFinal, "room_rescaled_no_reorder.stl");
%disp("STL saved as 'room_rescaled_no_reorder.stl'");
%%

function visualizeGeneralRoom(tri,txinates,rxinates)
figure;
trisurf(tri, ...
    'FaceAlpha', 0.3, ...
    'FaceColor', [.5 .5 .5], ...
    'EdgeColor', 'none');
view(60, 30);
hold on; axis equal; grid off;
xlabel('x'); ylabel('y'); zlabel('z');
% Plot edges
fe = featureEdges(tri,pi/20);
numEdges = size(fe, 1);
pts = tri.Points;
a = pts(fe(:,1),:); 
b = pts(fe(:,2),:); 
fePts = cat(1, reshape(a, 1, numEdges, 3), ...
           reshape(b, 1, numEdges, 3), ...
           nan(1, numEdges, 3));
fePts = reshape(fePts, [], 3);
plot3(fePts(:, 1), fePts(:, 2), fePts(:, 3), 'k', 'LineWidth', .5); 

hold on
tx = txinates;
rx = rxinates;
scatter3(tx(1), tx(2), tx(3), 'sb', 'filled');
scatter3(rx(1,1), rx(1,2), rx(1,3), 'sr', 'filled');
end






x = 2.1082;
y = 2.6416;
z = 2.9972;
ft2 = 0.6096;

frontwall = x*z;
backwall = frontwall;
leftwall = y*z;
rightwall = leftwall;
ceiling = x*y;
floor = ceiling;

t = [0,0,0];
r = [x,y,z];
visualizeGeneralRoom(roomFinal,t,r);

% For 180 deg t = [x,y-ft2,z]

% Absorption Coefficients
freq = [125,250,500,1000,2000,4000];
% floor_A = [0.1,0.15,0.25,0.3,0.3,0.3];
% wall_panel_A = [0.11,0.4,0.7,0.74,0.88,0.89];
% wall_drywall_A = [0.31,0.33,0.14,0.1,0.1,0.12];
% wall_window_A = [0.12,0.06,0.04,0.03,0.02,0.02];
% wall_door_A = [0.35,0.39,0.44,0.49,0.54,0.57];
% ceiling_A = [0.45,0.7,0.8,0.8,0.65,0.45];

% Room dimensions [2.1082, 2.6416,2.9667]
% Left wall calcs: sound panel area 2.16294 * 2.6416 = 5.7136
% sound panel 2 area 0.807-0.045 * 2.69115-2.1715 = .7620 *.5196 = 0.3959
% Left wall drywall total area = 7.9174 - 5.7136 - .5196 = 1.6842
% Left wall sound panel total area = 6.2332
% Back wall calcs: Sound Panels =
% ((.809212-.144438)*(2.67391-.639982))+((2.05271-1.25454)*(2.16651-.208892))+((2.06766-1.39162)*(2.86026-2.16294))
% = 3.386
% Backwall drywall = 2.9327

% Effective absorption for back wall
backwallAreas = [3.386, 2.9327]; % Panels, Drywall

alpha = [
    0.11 0.4 0.7 0.74 0.88 0.89 % absorption coeff for sound panels
    0.31 0.33 0.14 0.1 0.1 0.12 % absorption coeff for drywall
    ];

numFreqs = length(freq);
effAbs_backWall = zeros(1,numFreqs);

for n = 1:numFreqs
    weightedSum = alpha(1,n) * backwallAreas(1) + alpha(2,n) * backwallAreas(2);
    effAbs_backWall(n) = weightedSum/backwall;
end


% Effective absorption for left wall
leftwallAreas = [6.2332, 1.6842]; % Panels, drywall
effAbs_leftWall = zeros(1,numFreqs);
for n = 1:numFreqs
    weightedSum = alpha(1,n) * leftwallAreas(1) + alpha(2,n) * leftwallAreas(2);
    effAbs_leftWall(n) = weightedSum/leftwall;
end
