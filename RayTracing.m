% Initialize variables
freq = [125,250,500,1000,2000,4000]; % Freq vals for absoprtion and scatter coeffs
numFreqs = length(freq);
fs = 48000;
rng(0);
% Absorption Coefficients via ""
effAbsFloor = [0.1,0.15,0.25,0.3,0.3,0.3];
effAbsPan = [0.11,0.4,0.7,0.74,0.88,0.89];
effAbsWall = [0.02,0.02,0.02,0.02,0.02,0.02];
effAbsWindow = [0.12,0.06,0.04,0.03,0.02,0.02];
effAbsDoor = [0.35,0.39,0.44,0.49,0.54,0.57];
effAbsCeiling = [0.45,0.7,0.8,0.8,0.65,0.45];
effAbsBody = [0.15,0.25,0.35,0.45,0.55,0.65];

% Scatter Coeffs via LLM
effScFloor = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35];
effScPan = [0.25, 0.35, 0.45, 0.55, 0.65, 0.75];
effScWall = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05];
effScWindow = [0.05, 0.05, 0.05, 0.05, 0.05, 0.05];
effScDoor = [0.05, 0.07, 0.10, 0.15, 0.20, 0.25];
effScCeiling = [0.25, 0.35, 0.45, 0.55, 0.65, 0.75];
effScBody = [0.3, 0.35, 0.45, 0.55, 0.6, 0.65];

% STL files
room = stlread("STL_Files\room.stl");
ceiling = stlread("STL_Files\ceiling.stl");
floor = stlread("STL_Files\floor.stl");
northWall = stlread("STL_Files\northwall.stl");
southWall = stlread("STL_Files\southwall.stl");
eastWall = stlread("STL_Files\eastwall.stl");
westWall = stlread("STL_Files\westwall.stl");
box = stlread("STL_Files\boxthing.stl");
body = stlread("STL_Files\body.stl");
src = audioread("Audio_Files\Clave_Mono.wav");

% dimensions in meters
% 2.6416 = y = length = 104
% 2.9972 = z = height = 118
% 2.8666 = x = width = 112.85 in

% Only used as reference now that individual surfaces are being used
% Rescale room size to be correct
% % Step 1: Extract points and faces
% points = room.Points;            
% faces  = room.ConnectivityList;  
% 
% %Step 2: Current dimensions
% minCoords = min(points);
% maxCoords = max(points);
% dims = maxCoords - minCoords;   % [X, Y, Z]
% disp('Current dimensions [X, Y, Z]:')
% disp(dims)
% 
% %Step 3: Desired dimensions
% %Keep axis order the same: X, Y, Z
% targetDims = [2.1082, 2.6416, 2.9972];  % [X, Y, Z]
% 
% %Step 4: Compute scaling factors per axis
% scaleFactors = targetDims ./ dims;
% 
% %Step 5: Apply scaling
% points = points .* scaleFactors;
% 
% %Step 6: Rebuild triangulation
% roomFinal = triangulation(faces, points);
% 
% %Step 7: Visualize transparent
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

% Save STL using triangulation object directly
% stlwrite(roomFinal, "room_rescaled_no_reorder.stl");
%disp("STL saved as 'room_rescaled_no_reorder.stl'");

%% Dont touch this
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

%% scale, position, and rotate body appropriately
pts_body = body.Points;
faces_body = body.ConnectivityList;

%--- Step 1: Detect orientation and scale --------------------------------
dimsBody = max(pts_body) - min(pts_body);
disp("Original body dimensions [X Y Z] (in STL units):");
disp(dimsBody)

% Find the largest dimension (assumed to be height)
[~, tallestAxis] = max(dimsBody);

% Rotate if height isn't along Z
switch tallestAxis
    case 1 % height along X
        R = [0 0 1; 0 1 0; -1 0 0]; % rotate so X→Z
    case 2 % height along Y
        R = [1 0 0; 0 0 1; 0 -1 0]; % rotate so Y→Z
    otherwise
        R = eye(3); % already upright
end
pts_body = (R * pts_body')';
pts_body(:,3) = -pts_body(:,3);

%--- Step 2: Convert from mm→m if necessary ------------------------------
if max(dimsBody) > 3  % e.g., >3 meters → probably mm
    pts_body = pts_body / 1000;
end

%--- Step 3b & 4: Scale body, then align head while keeping feet at z=0 ---
minBody = min(pts_body);
maxBody = max(pts_body);
dimsBody = maxBody - minBody;

% Target height
targetHeight = 1.625; % meters
scaleFactor = targetHeight / dimsBody(3);
pts_body = pts_body * scaleFactor;

% Recompute min/max after scaling
minBody = min(pts_body);
maxBody = max(pts_body);

% Head center ~ear height below top
headZ = maxBody(3);
headCenter = [mean([minBody(1), maxBody(1)]), ...
              mean([minBody(2), maxBody(2)]), ...
              headZ - 0.09]; 

% Desired head position
r_x = 1.016; r_y = 0.9398; r_z = 1.4986;

% Compute translation for head
translation = [r_x, r_y, r_z] - headCenter;

% Apply translation
pts_body = pts_body + translation;

% Finally, shift feet to z=0 if slightly below
minZ = min(pts_body(:,3));
pts_body(:,3) = pts_body(:,3) - minZ; % ensure feet sit on floor


%--- Step 5: Visualize in room -------------------------------------------
% bodyFinal = triangulation(faces_body, pts_body);
% figure
% trisurf(roomFinal, 'FaceColor', 'cyan', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); hold on
% trisurf(bodyFinal, 'FaceColor', 'magenta', 'FaceAlpha', 0.7, 'EdgeColor', 'none');
% scatter3(r_x, r_y, r_z, 60, 'r', 'filled');
% axis equal
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('Human Body Oriented and Positioned in Room')
% view(3); camlight; lighting gouraud

%--- Step 6: Save transformed body as a new STL ---------------------------
 bodyFinal = triangulation(faces_body, pts_body);
% outputPath = "STL_Files/body_transformed.stl";
% stlwrite(bodyFinal, outputPath);
% disp("Saved transformed body STL as " + outputPath);

%%
function visualizeGeneralRoomMultiple(triList, txinates, rxinates)
% triList = cell array of triangulation objects
figure; hold on; axis equal; grid off;
xlabel('x'); ylabel('y'); zlabel('z');
view(60, 30);

faceColor = [0.5 0.5 0.5]; % same as visualizeGeneralRoom

for k = 1:length(triList)
    tri = triList{k};
    trisurf(tri, 'FaceAlpha', 0.3, 'FaceColor', faceColor, 'EdgeColor', 'none');
    
    % Plot edges
    fe = featureEdges(tri, pi/20);
    pts = tri.Points;
    a = pts(fe(:,1),:);
    b = pts(fe(:,2),:);
    numEdges = size(fe,1);
    fePts = cat(1, reshape(a,1,numEdges,3), reshape(b,1,numEdges,3), nan(1,numEdges,3));
    fePts = reshape(fePts,[],3);
    plot3(fePts(:,1), fePts(:,2), fePts(:,3), 'k', 'LineWidth', 0.5);
end

% Plot transmitter and receiver
scatter3(txinates(1), txinates(2), txinates(3), 'sb', 'filled');
scatter3(rxinates(1,1), rxinates(1,2), rxinates(1,3), 'sr', 'filled');

view(3); camlight; lighting gouraud;
end

%% Room visualization for testing
triList = {
    ceiling, floor, northWall, southWall, ...
    eastWall, westWall, box, bodyFinal
    };
x = 2.1082;
y = 2.6416;
z = 2.9972;
%ft2 = 0.6096;
ft2 = 0.7;
w = 0.0254 * 4; % distance from center head to ear


r_x = 1.016;
r_y = 0.955;
r_z = 1.4986;

headCenter = [r_x,r_y,r_z];
leftEar = headCenter + [-w,0,0];
rightEar = headCenter + [w,0,0];

frontwall = x*z;
backwall = frontwall;
leftwall = y*z;
rightwall = leftwall;
ceiling = x*y;
floor = ceiling;

% 180 deg
t = [r_x,r_y,r_z];
r = [leftEar,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);



% 0 deg
t = [r_x,r_y+ft2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);

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
angle = deg2rad(135);

t = headCenter+ft2*[cos(angle), sin(angle), 0];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);

% -157.5 deg
t = [r_x-ft2*(sqrt(2-sqrt(2))/2),r_y-ft2*sqrt(2+sqrt(2))/2,r_z];
r = [r_x,r_y,r_z];
%visualizeGeneralRoom(roomFinal,t,r);
%%


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


% Calculates the effective absorption coefficients of a surface given the
% individual abs coeff, areas of the individual surfaces, total area of the
% surface, and num of freq that have a corresponding coeff
function effAbs = calcEffective(alpha, areas, totalArea, numFreqs)
areas = areas(:);
[m, f] = size(alpha);

totalArea = sum(areas);

effAbs = sum(alpha .* areas, 1)./totalArea;
end


%% Effective absorption for back wall
alpha = [
    effAbsPan % absorption coeff for sound panels
    effAbsWall % absorption coeff for drywall
    ];

backwallAreas = [3.386, 2.9327]; % Panels, Drywall
effAbs_backWall = calcEffective(alpha,backwallAreas,backwall, numFreqs);


%% Effective absorption for left wall
leftwallAreas = [6.2332, 1.6842]; % Panels, drywall
effAbs_leftWall = calcEffective(alpha,leftwallAreas,leftwall,numFreqs);


%% Eff Abs for right wall
rightwallAreas = [2.7211, 2.6177, 2.1522]; % Panels, drywall, door
alpha = [
    effAbsPan % absorption coeff for sound panels
    effAbsWall % absorption coeff for drywall
    effAbsDoor % abs for acoustic door
    ];
rightwall = 7.491; % adjusted for box taking up space
effAbs_rightwall = calcEffective(alpha,rightwallAreas,rightwall,numFreqs);


%% Eff Abs for front wall
frontwallAreas = [1.0369, 0.9113, 3.3336]; % window, door, drywall
alpha = [
    effAbsWindow
    effAbsDoor
    effAbsWall
    ];
frontwall = backwall - 1.0369; %Adjusted for box
effAbs_frontwall = calcEffective(alpha,frontwallAreas,frontwall,numFreqs);

%%

% Via LLM
% Scattering coeffs for materials
% effScFloor = [0.1,0.15,0.2,0.25,0.3,0.35]
% effScCeiling = [0.25,0.35,0.45,0.55,0.65,0.75]
% Effective Scattering coeffs can be calculated same as absorption

%% Effective scattering for back wall
alpha = [
    effScPan % Scatter coeff for sound panels
    effScWall % Scatter coeff for drywall
    ];

backwallAreas = [3.386, 2.9327]; % Panels, Drywall
effSc_backWall = calcEffective(alpha,backwallAreas,backwall, numFreqs);

%% Effective scatter for left wall
effSc_leftWall = calcEffective(alpha,leftwallAreas,leftwall,numFreqs);

%% Eff scatter for right wall

alpha = [
    effScPan % absorption coeff for sound panels
    effScWall % absorption coeff for drywall
    effScDoor % abs for acoustic door
    ];

effSc_rightwall = calcEffective(alpha,rightwallAreas,rightwall,numFreqs);

%% Eff scatter for front wall
alpha = [
    effScWindow
    effScDoor
    effScWall
    ];

effSc_frontWall = calcEffective(alpha,frontwallAreas,frontwall,numFreqs);


%% Now for the simulation
absCoeff = {effAbs_frontwall, effAbs_rightwall, effAbs_leftWall, effAbs_backWall, effAbsCeiling, effAbsFloor, effAbsBody, effAbsWall};

%% 135 degrees
angle = deg2rad(-45);

t = headCenter+ft2*[cos(angle), sin(angle), 0];


r = [rightEar,r_y,r_z];
visualizeGeneralRoomMultiple(triList,t,r);


allPoints = [];
allConnectivity = [];
faceMaterialIndex = [];  % store which material each face uses
pointOffset = 0;

materials = [1,2,3,4,5,6,7,8]; % just IDs for each surface in triList

for k = 1:length(triList)
    tri = triList{k};
    pts = tri.Points;
    conn = tri.ConnectivityList + pointOffset;
    allPoints = [allPoints; pts];
    allConnectivity = [allConnectivity; conn];
    faceMaterialIndex = [faceMaterialIndex; repmat(materials(k), size(conn,1),1)];
    pointOffset = size(allPoints,1);
end

roomFinal = triangulation(allConnectivity, allPoints);

absCoeffMat = [effAbs_frontwall; effAbs_rightwall; effAbs_leftWall; ...
               effAbs_backWall; effAbsCeiling; effAbsFloor; ...
               effAbsBody; effAbsWall];  % 8x6 for 8 materials, 6 freqs

scCoeffMat = [effSc_frontWall; effSc_rightwall; effSc_leftWall; ...
               effSc_backWall; effScCeiling; effScFloor; ...
               effScBody; effScWall];

% Now map per-face using faceMaterialIndex
numFaces = size(roomFinal.ConnectivityList,1);
numFreqs = size(absCoeffMat,2);

MaterialAbsorption = zeros(numFaces, numFreqs);
MaterialScattering = zeros(numFaces, numFreqs);

for f = 1:numFaces
    matID = faceMaterialIndex(f);
    MaterialAbsorption(f,:) = absCoeffMat(matID,:);
    MaterialScattering(f,:) = scCoeffMat(matID,:);
end

[src, fs] = audioread("Audio_Files\Clave_Mono.wav");

% Ensure mono
if size(src,2) > 1
    src = mean(src,2); % average channels into mono
end


%% --- Convolve with right ear IR ---
irRight = acousticRoomResponse(roomFinal, t, rightEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
    ReceiverRadius = 0.02,...
    NumStochasticRays=1e5);

if size(irRight,2) > 1
    irRight = mean(irRight,2);
end

yRight = conv(src, irRight);

%% --- Convolve with left ear IR ---
irLeft = acousticRoomResponse(roomFinal, t, leftEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
    ReceiverRadius=0.02,...
    NumStochasticRays=1e5);

if size(irLeft,2) > 1
    irLeft = mean(irLeft,2);
end

yLeft = conv(src, irLeft);

%% --- Match lengths ---
len = max(length(yLeft), length(yRight));
yLeft(end+1:len) = 0;
yRight(end+1:len) = 0;

%% --- Combine into stereo ---
stereoSignal = [yLeft, yRight];

%% --- Normalize stereo output to avoid clipping ---
stereoSignal = stereoSignal / max(abs(stereoSignal(:)));

%% --- Play and save ---
sound(stereoSignal, fs);
audiowrite("Audio_Files\007meters\Clave135.wav", stereoSignal, fs);



%% 45 degrees
45
angle = deg2rad(45);

t = headCenter+ft2*[cos(angle), sin(angle), 0];


r = [rightEar,r_y,r_z];
visualizeGeneralRoomMultiple(triList,t,r);





%% --- Convolve with right ear IR ---
irRight = acousticRoomResponse(roomFinal, t, rightEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
    ReceiverRadius=0.02,...
    NumStochasticRays=1e5);

if size(irRight,2) > 1
    irRight = mean(irRight,2);
end

yRight = conv(src, irRight);

%% --- Convolve with left ear IR ---
irLeft = acousticRoomResponse(roomFinal, t, leftEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
   ReceiverRadius=0.02,...
   NumStochasticRays=1e5);

if size(irLeft,2) > 1
    irLeft = mean(irLeft,2);
end

yLeft = conv(src, irLeft);

%% --- Match lengths ---
len = max(length(yLeft), length(yRight));
yLeft(end+1:len) = 0;
yRight(end+1:len) = 0;

%% --- Combine into stereo ---
stereoSignal = [yLeft, yRight];

%% --- Normalize stereo output to avoid clipping ---
stereoSignal = stereoSignal / max(abs(stereoSignal(:)));

%% --- Play and save ---
sound(stereoSignal, fs);
audiowrite("Audio_Files\007meters\Clave45.wav", stereoSignal, fs);

%% -22.5 degrees
angle = deg2rad(112.5);

t = headCenter+ft2*[cos(angle), sin(angle), 0];


r = [rightEar,r_y,r_z];
visualizeGeneralRoomMultiple(triList,t,r);





%% --- Convolve with right ear IR ---
irRight = acousticRoomResponse(roomFinal, t, rightEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
    ReceiverRadius=0.02,...
    NumStochasticRays=1e5);

if size(irRight,2) > 1
    irRight = mean(irRight,2);
end

yRight = conv(src, irRight);

%% --- Convolve with left ear IR ---
irLeft = acousticRoomResponse(roomFinal, t, leftEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
   ReceiverRadius=0.02,...
   NumStochasticRays=1e5);

if size(irLeft,2) > 1
    irLeft = mean(irLeft,2);
end

yLeft = conv(src, irLeft);

%% --- Match lengths ---
len = max(length(yLeft), length(yRight));
yLeft(end+1:len) = 0;
yRight(end+1:len) = 0;

%% --- Combine into stereo ---
stereoSignal = [yLeft, yRight];

%% --- Normalize stereo output to avoid clipping ---
stereoSignal = stereoSignal / max(abs(stereoSignal(:)));

%% --- Play and save ---
sound(stereoSignal, fs);
audiowrite("Audio_Files\007meters\Clave22.wav", stereoSignal, fs);

%% 180 degrees
180
angle = deg2rad(-90);

t = headCenter+ft2*[cos(angle), sin(angle), 0];


r = [rightEar,r_y,r_z];
visualizeGeneralRoomMultiple(triList,t,r);





%% --- Convolve with right ear IR ---
irRight = acousticRoomResponse(roomFinal, t, rightEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
    ReceiverRadius=0.02,...
    NumStochasticRays=1e5);

if size(irRight,2) > 1
    irRight = mean(irRight,2);
end

yRight = conv(src, irRight);

%% --- Convolve with left ear IR ---
irLeft = acousticRoomResponse(roomFinal, t, leftEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
   ReceiverRadius=0.02,...
   NumStochasticRays=1e5);

if size(irLeft,2) > 1
    irLeft = mean(irLeft,2);
end

yLeft = conv(src, irLeft);

%% --- Match lengths ---
len = max(length(yLeft), length(yRight));
yLeft(end+1:len) = 0;
yRight(end+1:len) = 0;

%% --- Combine into stereo ---
stereoSignal = [yLeft, yRight];

%% --- Normalize stereo output to avoid clipping ---
stereoSignal = stereoSignal / max(abs(stereoSignal(:)));

%% --- Play and save ---
sound(stereoSignal, fs);
audiowrite("Audio_Files\007meters\Clave180.wav", stereoSignal, fs);

%% 112.5 degrees
112.5
angle = deg2rad(337.5);

t = headCenter+ft2*[cos(angle), sin(angle), 0];


r = [rightEar,r_y,r_z];
visualizeGeneralRoomMultiple(triList,t,r);





%% --- Convolve with right ear IR ---
irRight = acousticRoomResponse(roomFinal, t, rightEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
    ReceiverRadius=0.02,...
    NumStochasticRays=1e5);

if size(irRight,2) > 1
    irRight = mean(irRight,2);
end

yRight = conv(src, irRight);

%% --- Convolve with left ear IR ---
irLeft = acousticRoomResponse(roomFinal, t, leftEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
   ReceiverRadius=0.02,...
   NumStochasticRays=1e5);

if size(irLeft,2) > 1
    irLeft = mean(irLeft,2);
end

yLeft = conv(src, irLeft);

%% --- Match lengths ---
len = max(length(yLeft), length(yRight));
yLeft(end+1:len) = 0;
yRight(end+1:len) = 0;

%% --- Combine into stereo ---
stereoSignal = [yLeft, yRight];

%% --- Normalize stereo output to avoid clipping ---
stereoSignal = stereoSignal / max(abs(stereoSignal(:)));

%% --- Play and save ---
sound(stereoSignal, fs);
audiowrite("Audio_Files\007meters\Clave112.wav", stereoSignal, fs);


%% -67.5 degrees
-67.5
angle = deg2rad(157.5);

t = headCenter+ft2*[cos(angle), sin(angle), 0];


r = [rightEar,r_y,r_z];
visualizeGeneralRoomMultiple(triList,t,r);





%% --- Convolve with right ear IR ---
irRight = acousticRoomResponse(roomFinal, t, rightEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
    ReceiverRadius=0.02,...
    NumStochasticRays=1e5);

if size(irRight,2) > 1
    irRight = mean(irRight,2);
end

yRight = conv(src, irRight);

%% --- Convolve with left ear IR ---
irLeft = acousticRoomResponse(roomFinal, t, leftEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
   ReceiverRadius=0.02,...
   NumStochasticRays=1e5);

if size(irLeft,2) > 1
    irLeft = mean(irLeft,2);
end

yLeft = conv(src, irLeft);

%% --- Match lengths ---
len = max(length(yLeft), length(yRight));
yLeft(end+1:len) = 0;
yRight(end+1:len) = 0;

%% --- Combine into stereo ---
stereoSignal = [yLeft, yRight];

%% --- Normalize stereo output to avoid clipping ---
stereoSignal = stereoSignal / max(abs(stereoSignal(:)));

%% --- Play and save ---
sound(stereoSignal, fs);
audiowrite("Audio_Files\007meters\Clave67.wav", stereoSignal, fs);

%% -157.5 degrees
-157.5
angle = deg2rad(247.5);

t = headCenter+ft2*[cos(angle), sin(angle), 0];


r = [rightEar,r_y,r_z];
visualizeGeneralRoomMultiple(triList,t,r);





%% --- Convolve with right ear IR ---
irRight = acousticRoomResponse(roomFinal, t, rightEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
    ReceiverRadius=0.02,...
    NumStochasticRays=1e5);

if size(irRight,2) > 1
    irRight = mean(irRight,2);
end

yRight = conv(src, irRight);

%% --- Convolve with left ear IR ---
irLeft = acousticRoomResponse(roomFinal, t, leftEar, ...
    SampleRate=fs, ...
    Algorithm="stochastic ray tracing", ...
    BandCenterFrequencies=freq, ...
    MaterialAbsorption=MaterialAbsorption, ...
    MaterialScattering=MaterialScattering,...
   ReceiverRadius=0.02,...
   NumStochasticRays=1e5);

if size(irLeft,2) > 1
    irLeft = mean(irLeft,2);
end

yLeft = conv(src, irLeft);

%% --- Match lengths ---
len = max(length(yLeft), length(yRight));
yLeft(end+1:len) = 0;
yRight(end+1:len) = 0;

%% --- Combine into stereo ---
stereoSignal = [yLeft, yRight];

%% --- Normalize stereo output to avoid clipping ---
stereoSignal = stereoSignal / max(abs(stereoSignal(:)));

%% --- Play and save ---
sound(stereoSignal, fs);
audiowrite("Audio_Files\007meters\Clave157.wav", stereoSignal, fs);
