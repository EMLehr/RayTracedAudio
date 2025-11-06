fs = 48000;
% src = mono audio file
room = stlread("D:\Matlab Projects\raytracingproject\RayTracedAudio\room.stl");
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
% Step 2: Current dimensions
% minCoords = min(points);
% maxCoords = max(points);
% dims = maxCoords - minCoords;   % [X, Y, Z]
% disp('Current dimensions [X, Y, Z]:')
% disp(dims)
% 
% Step 3: Desired dimensions
% Keep axis order the same: X, Y, Z
% targetDims = [2.1082, 2.6416, 2.9972];  % [X, Y, Z]
% 
% Step 4: Compute scaling factors per axis
% scaleFactors = targetDims ./ dims;
% 
% Step 5: Apply scaling
% points = points .* scaleFactors;
% 
% Step 6: Rebuild triangulation
% roomFinal = triangulation(faces, points);
% 
% Step 7: Visualize transparent
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

w = 0.0254; % distance from center head to ear

left_x = x - w; % left ear pos
right_x = x-w;  % right ear pos

r_x = 1.016;
r_y = 0.9398;
r_z = 1.4986;

frontwall = x*z;
backwall = frontwall;
leftwall = y*z;
rightwall = leftwall;
ceiling = x*y;
floor = ceiling;

% 180 deg
t = [r_x,r_y-ft2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);


% 0 deg
t = [r_x,r_y+ft2,r_z];
r = [r_x,r_y,r_z];
visualizeGeneralRoom(roomFinal,t,r);

% -22.5 deg
t = [r_x-ft2*sqrt(2-sqrt(2))/2,r_y+ft2*sqrt(2+sqrt(2))/2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);


% 45 deg
t = [r_x+ft2*sqrt(2)/2,r_y-ft2*sqrt(2)/2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);

% -67.5 deg
t = [r_x-ft2*sqrt(2+sqrt(2))/2,r_y+ft2*sqrt(2-sqrt(2))/2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);


% 112.5 deg
t = [r_x+ft2*(sqrt(2+sqrt(2))/2),r_y-ft2*sqrt(2-sqrt(2))/2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);

% 135 deg
t = [r_x+ft2*sqrt(2)/2,r_y-ft2*sqrt(2)/2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);

% -157.5 deg
t = [r_x-ft2*(sqrt(2-sqrt(2))/2),r_y-ft2*sqrt(2+sqrt(2))/2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);


% Absorption Coefficients
% thin carpet over thin felt on concrete floor_A = [0.1,0.15,0.25,0.3,0.3,0.3];
% acoustic banner wall_panel_A = [0.11,0.4,0.7,0.74,0.88,0.89];
% painted plaster surface on masonry wall wall_drywall_A = [0.02,0.02,0.02,0.02,0.02,0.02];
% 6mm glass wall_window_A = [0.12,0.06,0.04,0.03,0.02,0.02];
% acoustic door wall_door_A = [0.35,0.39,0.44,0.49,0.54,0.57];
% gypsum plaster tiles ceiling_A = [0.45,0.7,0.8,0.8,0.65,0.45];

% Scatter Coeff
% panel = [0.25, 0.35, 0.45, 0.55, 0.65, 0.75];
% door = [0.05 0.07 0.10 0.15 0.20 0.25];


% Room dimensions [2.1082, 2.6416,2.9667]
% Left wall calcs: sound panel area 2.16294 * 2.6416 = 5.7136
% sound panel 2 area 0.807-0.045 * 2.69115-2.1715 = .7620 *.5196 = 0.3959
% Left wall drywall total area = 7.9174 - 5.7136 - .5196 = 1.6842
% Left wall sound panel total area = 6.2332
% Back wall calcs: Sound Panels =
% ((.809212-.144438)*(2.67391-.639982))+((2.05271-1.25454)*(2.16651-.208892))+((2.06766-1.39162)*(2.86026-2.16294))
% = 3.386
% Backwall drywall = 2.9327
% Right wall calcs:
% Panel 1 = (1.4684-0.0118314)*(2.96562-2.2636) = 1.0225
% panel 2 = (1.4684-0.962617)*(2.2308-0.882232) = 0.6821
% Panel 3 =
% ((2.6416-0.962617)*(0.868334-.008933998))-((2.6416-2.10029)*(0.796715-.00893398))
% = 1.0165
% Right wall total sound panels area = 2.7211
% Right wall total area = 7.491
% Right wall door area = (0.962617-0.0118314)*2.2636 = 2.1522
%rightwall drywall area = 7.491 - 2.7211 - 2.1522 = 2.6177
% Frontwall area
% Box obstruction = 1.30145 * 0.796715 = 1.0369
% Front wall door area = (2.06766 - 1.58012) * (1.86924) = 0.9113
% Glass area = box area
% front wall drywall area = 6.3187 - (0.9113 + 1.0369*2) = 3.3336

freq = [125,250,500,1000,2000,4000];
numFreqs = length(freq);

% Calculates the effective absorption coefficients of a surface given the
% individual abs coeff, areas of the individual surfaces, total area of the
% surface, and num of freq that have a corresponding coeff
function effAbs = calcEffectiveAbsorption(alpha, areas, totalArea, numFreqs)
areas = areas(:);
[m, f] = size(alpha);

totalArea = sum(areas);

effAbs = sum(alpha .* areas, 1)./totalArea;
end


% Effective absorption for back wall
alpha = [
    0.11 0.4 0.7 0.74 0.88 0.89 % absorption coeff for sound panels
    0.31 0.33 0.14 0.1 0.1 0.12 % absorption coeff for drywall
    ];

backwallAreas = [3.386, 2.9327]; % Panels, Drywall
effAbs_backWall = calcEffectiveAbsorption(alpha,backwallAreas,backwall, numFreqs);

% Effective absorption for left wall
leftwallAreas = [6.2332, 1.6842]; % Panels, drywall
effAbs_leftWall = calcEffectiveAbsorption(alpha,leftwallAreas,leftwall,numFreqs);

% Eff Abs for right wall
rightwallAreas = [2.7211, 2.6177, 2.1522]; % Panels, drywall, door
alpha = [
    0.11 0.4 0.7 0.74 0.88 0.89 % absorption coeff for sound panels
    0.31 0.33 0.14 0.1 0.1 0.12 % absorption coeff for drywall
    0.35 0.39 0.44 0.49 0.54 0.57 % abs for acoustic door
    ];
rightwall = 7.491; % adjusted for box taking up space
effAbs_rightwall = calcEffectiveAbsorption(alpha,rightwallAreas,rightwall,numFreqs);

% Eff Abs for front wall
frontwallAreas = [1.0369, 0.9113, 3.3336]; % window, door, drywall
alpha = [
    [0.12 0.06 0.04 0.03 0.02 0.02]
    [0.35 0.39 0.44 0.49 0.54 0.57]
    [0.02 0.02 0.02 0.02 0.02 0.02]
    ];
frontwall = backwall - 1.0369; %Adjusted for box
effAbsfrontwall = calcEffectiveAbsorption(alpha,frontwallAreas,frontwall,numFreqs);

effAbsfloor = [0.1,0.15,0.25,0.3,0.3,0.3];

effAbsceiling = [0.45,0.7,0.8,0.8,0.65,0.45];