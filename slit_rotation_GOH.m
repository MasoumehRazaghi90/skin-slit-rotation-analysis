clear; close all; clc;

%% Plot settings
fontSize=20;
faceAlpha1=0.8; %transparency
markerSize=40; %For plotted points
markerSize2=10; %For nodes on patches
lineWidth1=1; %For meshes
cMap=spectral(250); %colormap


%% Control parameters

meshSearchTolerance=1e-3; %Tolerance for logic on mesh searching thresholds


%Material parameter set

c_ini     = 1.2099; % Initial slope
k1_ini    = 61.790; % "Modulus of fibers"
k2_ini    = 3.4157e-7; %  "non-linearity" or J-shape
kappa_ini = 0.2986; % [0,1/3], "distance between graphs"
k_ini     = 1.2099*1000; % Bulk modulus
gamma     = 41;
LangerAngle_perp = 0/180*pi; % Angle of Langer line with loading direction
LangerAngle_para = 90/180*pi; % Angle of Langer line with loading direction
angleSelection=1;%1 perpendicular to LL and 2 parallel
LangerAngles = [LangerAngle_perp, LangerAngle_para];


displacementMagnitude = 60;

% FEA control settings
% numTimeSteps=30; %Number of time steps desired
% max_refs=20; %Max reforms
% max_ups=0; %Set to zero to use full-Newton iterations
% opt_iter=10; %Optimum number of iterations
% max_retries=2; %Maximum number of retires
% dtmin=(1/numTimeSteps)/100; %Minimum time step size
% dtmax=(1/numTimeSteps); %Maximum time step size
% runMode='internal';

% % FEA control settings
numTimeSteps=60; %Number of time steps desired
max_refs=1000; %Max reforms
max_ups=20; %Set to zero to use full-Newton iterations
opt_iter=50; %Optimum number of iterations
max_retries=30; %Maximum number of retires
dtmin=(1/numTimeSteps)/10000; %Minimum time step size
dtmax=(1/numTimeSteps); %Maximum time step size
runMode='internal';


% Path names
% Path names
defaultFolder = fileparts(mfilename('fullpath'));
savePath = 'C:\Sahar-GIBBON-codes\paper3-final\final results-GOH\paper2-revisiong\GUPTA_2fibers\data\gupta_horizontal';


% Defining file names
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_stress_prin=[febioFebFileNamePart,'_stress_prin_out.txt']; %Log file name for exporting principal stress
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%% =========================================================
% Hybrid GIBBON-style code for slit geometry
% Practical goal:
%   0 deg  -> use global rotation of slit centers (works better)
%   45 deg -> use local rotation about each slit center (works better)
%   90 deg -> use global rotation of slit centers (works better)
%
% Valid intended angles for this script: 0, 45, 90
%% =========================================================

%% ----------------- Basic settings
pointSpacing = 2.5;   % target point spacing for boundary resampling 2 for 90, 45 works for 0 langer 90 work)
t            = 3;     % extrusion thickness (mm)

%% ----------------- Plate dimensions
plateWidth   = 60;    % mm
plateHeight  = 60;    % mm
plateCenter  = [plateWidth/2, plateHeight/2];

%% ----------------- Slit layout
numColumns   = 17;

% Top / bottom slit size
slitWidth    = 0.96;
slitHeight   = 4.8;

% Middle slit size
numVerticalPerColumn = 4;
vertSlitWidth        = 0.96;
vertSlitHeight       = 9.6;
vertGap              = 2.45;

% Margins
edgeMargin = 2.6;
topMargin  = 0;

%% ----------------- Choose angle: only 0,45,90
rotationAngle =90 ;   % set to 0 or 45 or 90
if rotationAngle == 0
    rotationLabel = 'Perpendicular';
elseif rotationAngle == 90
    rotationLabel = 'Parallel';
elseif rotationAngle == 45
    rotationLabel = '45 deg';
end


if ~ismember(rotationAngle,[0 45 90])
    error('This version is intentionally tuned only for 0, 45, and 90 degrees.');
end

%% ----------------- Spacing between columns
usableWidth = plateWidth - 2*edgeMargin;
gapBetweenColumns = (usableWidth - numColumns*slitWidth) / (numColumns-1);

fprintf('Horizontal gap between columns = %.3f mm\n', gapBetweenColumns);

%% ----------------- Create base plate as polyshape
outerX = [0 plateWidth plateWidth 0];
outerY = [0 0 plateHeight plateHeight];
pg = polyshape(outerX, outerY);

%% ----------------- Build slit list first
slitPolys = {};

% Top and bottom slits
for col = 1:numColumns

    xLeft = edgeMargin + (col-1)*(slitWidth + gapBetweenColumns);

    % top slit
    yBottomTop = plateHeight - slitHeight - topMargin;
    rectTop = [xLeft               yBottomTop;
        xLeft+slitWidth     yBottomTop;
        xLeft+slitWidth     yBottomTop+slitHeight;
        xLeft               yBottomTop+slitHeight];
    slitPolys{end+1} = rectTop; %#ok<SAGROW>

    % bottom slit
    yBottomBottom = topMargin;
    rectBottom = [xLeft               yBottomBottom;
        xLeft+slitWidth     yBottomBottom;
        xLeft+slitWidth     yBottomBottom+slitHeight;
        xLeft               yBottomBottom+slitHeight];
    slitPolys{end+1} = rectBottom; %#ok<SAGROW>
end

% Middle slits
yLowerLimit = topMargin + slitHeight;
yStartVert  = yLowerLimit + vertGap;

for col = 1:numColumns

    xLeft = edgeMargin + (col-1)*(slitWidth + gapBetweenColumns);

    for v = 0:(numVerticalPerColumn-1)

        yBottomVert = yStartVert + v*(vertSlitHeight + vertGap);

        rectMid = [xLeft                   yBottomVert;
            xLeft+vertSlitWidth     yBottomVert;
            xLeft+vertSlitWidth     yBottomVert+vertSlitHeight;
            xLeft                   yBottomVert+vertSlitHeight];

        slitPolys{end+1} = rectMid; %#ok<SAGROW>
    end
end

%% ----------------- Apply geometry rule depending on angle
switch rotationAngle
    case 45
        % For 45 degrees:
        % rotate each slit about its OWN center
        for i = 1:numel(slitPolys)
            coords = slitPolys{i};
            coordsRot = rotatePolygonAboutOwnCenter(coords, rotationAngle);
            slitPg = polyshape(coordsRot(:,1), coordsRot(:,2));
            pg = subtract(pg, slitPg);
        end

    case {0,90}
        % For 0 and 90 degrees:
        % rotate slit positions about PLATE center
        theta = deg2rad(rotationAngle);
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

        for i = 1:numel(slitPolys)
            coords = slitPolys{i};
            coordsC = coords - plateCenter;
            coordsR = (R * coordsC')';
            coordsF = coordsR + plateCenter;
            slitPg = polyshape(coordsF(:,1), coordsF(:,2));
            pg = subtract(pg, slitPg);
        end
end

%% ----------------- Show 2D polyshape before meshing
figure;
plot(pg,'FaceColor',[0.85 0.92 1.00],'EdgeColor','k');
axis equal;
xlabel('X (mm)');
ylabel('Y (mm)');
title(sprintf('Polyshape geometry, slit angle = %s', rotationLabel));

%% ----------------- Extract outer boundary + holes
N = numboundaries(pg);

outerCurve = [];
holeCurves = {};

for q = 1:N
    [xb,yb] = boundary(pg,q);

    if numel(xb)>1 && xb(1)==xb(end) && yb(1)==yb(end)
        xb(end) = [];
        yb(end) = [];
    end

    thisCurve = [xb(:) yb(:)];

    if ishole(pg,q)
        if signedArea2D(thisCurve(:,1),thisCurve(:,2)) > 0
            thisCurve = flipud(thisCurve);
        end
        holeCurves{end+1} = thisCurve;
    else
        if signedArea2D(thisCurve(:,1),thisCurve(:,2)) < 0
            thisCurve = flipud(thisCurve);
        end
        outerCurve = thisCurve;
    end
end

regionCell = [{outerCurve} holeCurves];

%% ----------------- Resample curves for GIBBON
interpMethod = 'linear';
closeLoopOpt = 1;
resampledCurves = cell(size(regionCell));

for k = 1:numel(regionCell)
    Vk = regionCell{k};
    indMust = 1:size(Vk,1);
    Vk = evenlySpaceCurve(Vk,pointSpacing,interpMethod,closeLoopOpt,indMust);
    resampledCurves{k} = Vk;
end

%% ----------------- 2D triangulation using GIBBON
resampleCurveOpt = 0;
[F,V] = regionTriMesh2D(resampledCurves,pointSpacing,resampleCurveOpt,0);
V(:,3) = 0;

%% ----------------- Show 2D mesh
cFigure; hold on;
gpatch(F,V,'w','k');
for k = 1:numel(resampledCurves)
    plotV(resampledCurves{k},'r-','LineWidth',1.5);
end
title(sprintf('2D triangular mesh, slit angle = %s', rotationLabel));
axisGeom;
camlight headlight;
drawnow;

%% ----------------- Extrude to 3D
numSteps = max(2,ceil(t/pointSpacing));
[E,VE,~,~] = patchThick(F,V,1,t,numSteps);
nElem = size(E,1);

if rotationAngle == 0
    rotName = 'vertical';
elseif rotationAngle == 45
    rotName = 'diagonal';
elseif rotationAngle == 90
    rotName = 'horizontal';
end

savePath = sprintf('%s-%d-ps%g', savePath, nElem, pointSpacing);

if ~exist(savePath,'dir')
    mkdir(savePath);
end

fprintf('Results folder = %s\n', savePath);

fprintf('\n');
fprintf('Number of nodes    = %d\n', size(VE,1));
fprintf('Number of elements = %d\n', size(E,1));
fprintf('Element type       = penta6\n');
fprintf('\n');

FE = element2patch(E,[],'penta6');
%% ----------------- Boundary faces
indb = tesBoundary(FE);

Fb = { ...
    FE{1}(indb{1},:), ...
    FE{2}(indb{2},:) ...
    };

%% ----------------- Show 3D mesh
cFigure; hold on;
gpatch(FE,VE,'bw','k',1);
title(sprintf('3D extruded mesh, slit angle = %s', rotationLabel));
axisGeom;
camlight headlight;
drawnow;



%% Defining the boundary conditions

%Prescribed displacement nodes
bcPrescribeList = find(VE(:,2) > max(VE(:,2)) - meshSearchTolerance);
bcSupportList   = find(VE(:,2) < min(VE(:,2)) + meshSearchTolerance);

%%
% Visualizing boundary conditions. Markers plotted on the semi-transparent
% model denote the nodes in the various boundary condition lists.

hf=cFigure;
title('Boundary conditions','FontSize',fontSize);
xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
hold on;

gpatch(Fb,VE,'kw','k',0.5);

hl(1)=plotV(VE(bcPrescribeList,:),'r.','MarkerSize',markerSize);
hl(2)=plotV(VE(bcSupportList,:),'g.','MarkerSize',markerSize);

legend(hl,{'Tensile Positive', 'Tensile Negative'});

axisGeom(gca,fontSize);
camlight headlight;
drawnow;


%% DEFINE FIBRE DIRECTIONS
nElem=length(E);
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%%
setenv('GIBBON_cFigure_close','0');

results = struct('lambda',[],'force',[],'energy_mean',[],'energy_sum',[],'damage',[]);
results(2) = results;

angles_deg = [0 90];
for ia = 1:numel(angles_deg)
    ang = angles_deg(ia);


    if ang == 0
        % mu_ini_local    = mu_perp_ini;
        % sigma_ini_local = sigma_perp_ini;
        % Dmax_ini_local  = Dmax_perp_ini;
        angLabel        = 'Perpendicular';
    elseif ang == 90
        % mu_ini_local    = mu_para_ini;
        % sigma_ini_local = sigma_para_ini;
        % Dmax_ini_local  = Dmax_para_ini;
        angLabel        = 'Parallel';
    else
        error('Only 0 and 90 degrees are supported in this setup.');
    end

    % Define fiber direction by rotating vector [1 0 0] around z-axis
    R = euler2DCM([0,0,-ang*pi/180]);
    e1 = (R*[1 0 0]')';  % Rotated fiber direction

    e1_dir = repmat(e1, nElem, 1);
    e2_dir = repmat([0 0 1], nElem, 1);  % Keep this fixed
    e3_dir = cross(e1_dir, e2_dir, 2);   % Orthogonal vector

    % Plot material directions
    V_cent = patchCentre(E, VE);  % Make sure VE is defined

    cFigure; hold on;
    title(sprintf('Fiber direction = %s', angLabel), 'FontSize', 20);
    gpatch(Fb, VE, 'kw', 'none', 0.25);
    quiverVec(V_cent, e1_dir, pointSpacing, 'r');
    quiverVec(V_cent, e2_dir, pointSpacing, 'g');
    quiverVec(V_cent, e3_dir, pointSpacing, 'b');
    axisGeom(gca, fontSize);
    camlight headlight;
    drawnow;



    %%
    % Visualizing material directions

    V_cent = patchCentre(E,VE);

    hf=cFigure; hold on;
    title(sprintf('Material directions  |  Fiber Direction = %s', angLabel), 'FontSize', fontSize);

    gpatch(Fb,VE,'kw','none',0.25);

    hf(1)=quiverVec(V_cent,e1_dir,pointSpacing,'r');
    hf(2)=quiverVec(V_cent,e2_dir,pointSpacing/2,'g');
    hf(3)=quiverVec(V_cent,e3_dir,pointSpacing/2,'b');

    legend(hf,{'e1-direction','e2-direction','e3-direction'});
    axisGeom(gca,fontSize);
    camlight headlight;


    legend(hf,{'e1-direction (fiber)','e2-direction','e3-direction'});
    axisGeom(gca,fontSize);
    camlight headlight;
    drawnow;


    %% Defining the FEBio input structure
    % See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
    % manual.

    %Get a template with default settings
    [febio_spec]=febioStructTemplate;

    %febio_spec version
    febio_spec.ATTR.version='4.0';

    %Module section
    febio_spec.Module.ATTR.type='solid';

    %Control section
    febio_spec.Control.analysis='STATIC';
    febio_spec.Control.time_steps=numTimeSteps;
    febio_spec.Control.step_size=1/numTimeSteps;
    febio_spec.Control.solver.max_refs=max_refs;
    febio_spec.Control.solver.qn_method.max_ups=max_ups;
    febio_spec.Control.time_stepper.dtmin=dtmin;
    febio_spec.Control.time_stepper.dtmax=dtmax;
    febio_spec.Control.time_stepper.max_retries=max_retries;
    febio_spec.Control.time_stepper.opt_iter=opt_iter;

    %Material section

    %% Material section: ONE GOH block
    %Material section

    febio_spec.Material.material{1}.ATTR.type = 'uncoupled solid mixture';
    materialName1='Material1_HGO';
    febio_spec.Material.material{1}.ATTR.name=materialName1;
    febio_spec.Material.material{1}.ATTR.type = 'solid mixture';
    febio_spec.Material.material{1}.ATTR.id=1;


    % GOH‑specific material parameters
    febio_spec.Material.material{1}.solid{1}.ATTR.type='HGO unconstrained';
    febio_spec.Material.material{1}.solid{1}.c     = c_ini/2;
    febio_spec.Material.material{1}.solid{1}.k1    = k1_ini;
    febio_spec.Material.material{1}.solid{1}.k2    = k2_ini;
    febio_spec.Material.material{1}.solid{1}.gamma = gamma;
    febio_spec.Material.material{1}.solid{1}.kappa = kappa_ini;
    febio_spec.Material.material{1}.solid{1}.k     = k_ini;


    % GOH‑specific material parameters
    febio_spec.Material.material{1}.solid{2}.ATTR.type='HGO unconstrained';
    febio_spec.Material.material{1}.solid{2}.c     = c_ini/2;
    febio_spec.Material.material{1}.solid{2}.k1    = k1_ini;
    febio_spec.Material.material{1}.solid{2}.k2    = k2_ini;
    febio_spec.Material.material{1}.solid{2}.gamma = -gamma;
    febio_spec.Material.material{1}.solid{2}.kappa = kappa_ini;
    febio_spec.Material.material{1}.solid{2}.k     = k_ini;
    % Mesh section
    % -> Nodes

    %%Area of interest
    febio_spec.Mesh.Nodes{1}.ATTR.name='Object1'; %The node set name
    febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(VE,1))'; %The node id's
    febio_spec.Mesh.Nodes{1}.node.VAL=VE; %The nodel coordinates

    % -> Elements
    partName1='KirigamiGripper';
    febio_spec.Mesh.Elements{1}.ATTR.name=partName1; %Name of this part
    febio_spec.Mesh.Elements{1}.ATTR.type='penta6'; %Element type
    febio_spec.Mesh.Elements{1}.elem.ATTR.id=(1:1:size(E,1))'; %Element id's
    febio_spec.Mesh.Elements{1}.elem.VAL=E; %The element matrix

    % -> NodeSets
    nodeSetName1='bcPrescribeList1';
    nodeSetName2='bcSupportList2';
    nodeSetName3='bcSupportList3';

    febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
    febio_spec.Mesh.NodeSet{1}.VAL=mrow(bcPrescribeList);

    febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
    febio_spec.Mesh.NodeSet{2}.VAL=mrow(bcSupportList);

    febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
    febio_spec.Mesh.NodeSet{3}.VAL=bcSupportList(1);

    %MeshDomains section
    febio_spec.MeshDomains.SolidDomain.ATTR.name=partName1;
    febio_spec.MeshDomains.SolidDomain.ATTR.mat=materialName1;

    %MeshData section
    % -> ElementData
    febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partName1;
    febio_spec.MeshData.ElementData{1}.ATTR.type='mat_axis';

    for q=1:1:size(E,1)
        febio_spec.MeshData.ElementData{1}.elem{q}.ATTR.lid=q;
        febio_spec.MeshData.ElementData{1}.elem{q}.a=e1_dir(q,:);
        febio_spec.MeshData.ElementData{1}.elem{q}.d=e2_dir(q,:);
    end

    %% Boundary condition section
    % -> Fix boundary conditions

    % Bottom surface: fully fixed (x,y,z)
    febio_spec.Boundary.bc{1}.ATTR.name = 'fix_bottom';
    febio_spec.Boundary.bc{1}.ATTR.type = 'zero displacement';
    febio_spec.Boundary.bc{1}.ATTR.node_set = nodeSetName2; % bottom nodes
    febio_spec.Boundary.bc{1}.x_dof = 1;
    febio_spec.Boundary.bc{1}.y_dof = 1;
    febio_spec.Boundary.bc{1}.z_dof = 1;

    % Top surface: fixed in x and z only
    febio_spec.Boundary.bc{2}.ATTR.name = 'fix_top_xz';
    febio_spec.Boundary.bc{2}.ATTR.type = 'zero displacement';
    febio_spec.Boundary.bc{2}.ATTR.node_set = nodeSetName1; % top nodes
    febio_spec.Boundary.bc{2}.x_dof = 1;
    febio_spec.Boundary.bc{2}.y_dof = 0;
    febio_spec.Boundary.bc{2}.z_dof = 1;

    % Top surface: prescribed displacement in y
    febio_spec.Boundary.bc{3}.ATTR.name = 'top_displacement_y';
    febio_spec.Boundary.bc{3}.ATTR.type = 'prescribed displacement';
    febio_spec.Boundary.bc{3}.ATTR.node_set = nodeSetName1;
    febio_spec.Boundary.bc{3}.dof = 'y';
    febio_spec.Boundary.bc{3}.value.ATTR.lc = 1;
    febio_spec.Boundary.bc{3}.value.VAL = displacementMagnitude;
    febio_spec.Boundary.bc{3}.relative = 0;

    %LoadData section
    % -> load_controller
    febio_spec.LoadData.load_controller{1}.ATTR.name='LC_1';
    febio_spec.LoadData.load_controller{1}.ATTR.id=1;
    febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
    febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
    %febio_spec.LoadData.load_controller{1}.extend='CONSTANT';
    febio_spec.LoadData.load_controller{1}.points.pt.VAL=[0 0; 1 1];

    %Output section
    % -> log file
    logNamePrefix = sprintf('ang_%03d', ang);
    febioLogFileName        = [logNamePrefix, '.txt'];
    febioLogFileName_disp   = [logNamePrefix, '_disp_out.txt'];
    febioLogFileName_stress_prin = [logNamePrefix, '_stress_prin_out.txt'];
    febioLogFileName_force  = [logNamePrefix, '_force_out.txt'];

    febio_spec.Output.logfile.ATTR.file = febioLogFileName;

    febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';

    febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
    febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
    febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';

    febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress_prin;
    febio_spec.Output.logfile.element_data{1}.ATTR.data='s1;s2;s3';
    febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';

    % strain energy density
    febioLogFileName_energy = [logNamePrefix, '_energy_out.txt'];
    febio_spec.Output.logfile.element_data{2}.ATTR.file = febioLogFileName_energy;
    febio_spec.Output.logfile.element_data{2}.ATTR.data = 'sed';
    febio_spec.Output.logfile.element_data{2}.ATTR.delim = ',';

    % damage
    febioLogFileName_damage = [logNamePrefix, '_damage_out.txt'];
    febio_spec.Output.logfile.element_data{3}.ATTR.file = febioLogFileName_damage;
    febio_spec.Output.logfile.element_data{3}.ATTR.data = 'D';
    febio_spec.Output.logfile.element_data{3}.ATTR.delim = ',';


    % Plotfile section
    febio_spec.Output.plotfile.compression=0;

    %% Quick viewing of the FEBio input file structure
    % The |febView| function can be used to view the xml structure in a MATLAB
    % figure window.

    %%
    %%|febView(febio_spec); %Viewing the febio file|

    %% Exporting the FEBio input file
    % Exporting the febio_spec structure to an FEBio input file is done using
    % the |febioStruct2xml| function.

    febFileName = fullfile(savePath, sprintf('%s_ang_%03d.feb', febioFebFileNamePart, ang));
    febioStruct2xml(febio_spec, febFileName);
    %Exporting to file and domNode
    %system(['gedit ',febioFebFileName,' &']);

    %% Running the FEBio analysis
    % To run the analysis defined by the created FEBio input file the
    % |runMonitorFEBio| function is used. The input for this function is a
    % structure defining job settings e.g. the FEBio input file name. The
    % optional output runFlag informs the user if the analysis was run
    % succesfully.

    febioAnalysis.run_filename=febFileName; %The input file name
    febioAnalysis.run_logname=fullfile(savePath,febioLogFileName); %The name for the log file
    febioAnalysis.disp_on=0; %Display information on the command window
    febioAnalysis.runMode=runMode;
    febioAnalysis.maxLogCheckTime=10; %Max log file checking time

    [runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

    %% Import FEBio results

    if runFlag==1 %i.e. a succesful run
        %%
        % Importing nodal forces from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),0,1);

        %Access data
        N_force_mat=dataStruct.data; %Displacement
        timeVec=dataStruct.time; %Time

        %%
        % Importing nodal displacements from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),0,1);

        %Access data
        N_disp_mat=dataStruct.data; %Displacement

        %Create deformed coordinate set
        V_DEF=N_disp_mat+repmat(VE,[1 1 size(N_disp_mat,3)]);

        %%
        %% True Poisson's Ratio vs. Axial Stretch

        % Reference dimensions
        w0 = max(VE(:,1)) - min(VE(:,1));  % Initial width (x-direction)
        h0 = max(VE(:,2)) - min(VE(:,2));  % Initial height (y-direction)
        A0 = w0 * h0;                    % Initial area

        nSteps = size(V_DEF,3);
        lambda_x_vec = zeros(nSteps,1);
        lambda_y_vec = zeros(nSteps,1);
        nu_true      = zeros(nSteps,1);
        J_vec        = zeros(nSteps,1);

        for i = 1:nSteps
            V_temp = V_DEF(:,:,i);

            w = max(V_temp(:,1)) - min(V_temp(:,1)); % Current width
            h = max(V_temp(:,2)) - min(V_temp(:,2)); % Current height
            A = w * h;

            lambda_x_vec(i) = w / w0;
            lambda_y_vec(i) = h / h0;
            J_vec(i)        = A / A0;%area ratio

            % True Poisson’s ratio (log-based definition)
            if lambda_y_vec(i) > 0 && lambda_x_vec(i) > 0
                nu_true(i) = -log(lambda_x_vec(i)) / log(lambda_y_vec(i));
            else
                nu_true(i) = NaN;
            end
        end

        %% Compute total reaction force in loading direction

        nSteps = size(N_force_mat,3);
        Fy_vec = zeros(nSteps,1);

        for i = 1:nSteps
            Fy_vec(i) = sum(N_force_mat(:,2,i)); % sum of reaction forces in y
        end

        %% Poisson’s ratio vs. Axial Stretch

        % Plot true Poisson’s ratio
        cFigure; hold on;
        plot(lambda_y_vec, nu_true, 'b-', 'LineWidth', 5);
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('\nu_{\rm true} (Poisson''s Ratio)', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');
        title(sprintf('Poisson''s Ratio vs. Stretch  |  Fiber Direction = %s', angLabel), 'FontSize', 20,'FontWeight', 'bold');
        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]); % restrict x-axis to start at 1
        axis square;
        box on;grid on;

        %% Area Ratio J vs. Axial Stretch

        cFigure; hold on;
        plot(lambda_y_vec, J_vec, 'k-', 'LineWidth', 5);
        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('J (Area Ratio)', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');
        title(sprintf('Area Ratio J vs. Axial Stretch  |  Fiber Direction = %s', angLabel), 'FontSize', 20,'FontWeight', 'bold');
        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]); % restrict x-axis to start at 1
        axis square;
        box on;grid on;
        %% Force vs Axial Stretch

        cFigure; hold on;
        plot(lambda_y_vec, Fy_vec, 'k-', 'LineWidth', 5);

        xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', 'FontSize', 20,'FontWeight', 'bold');
        ylabel('Reaction Force F_y', 'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');

        title(sprintf('Force vs Axial Stretch  |  Fiber Direction = %s', angLabel), ...
            'FontSize',20,'FontWeight','bold');

        set(gca, 'FontSize', 16,'FontWeight', 'bold');
        xlim([1 max(lambda_y_vec)]);
        axis square;
        box on; grid on;

        %%
        % Importing element stress from a log file
        dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress_prin),0,1);

        %Access data
        E_stress_mat=dataStruct.data;
        % Import strain energy density
        dataStruct = importFEBio_logfile(fullfile(savePath,febioLogFileName_energy),0,1);
        E_energy_mat = dataStruct.data;

        % Import damage
        dataStruct = importFEBio_logfile(fullfile(savePath,febioLogFileName_damage),0,1);
        E_damage_mat = dataStruct.data;
        % ===== STRAIN ENERGY vs STRETCH =====

        % nSteps = size(E_energy_mat,3);

        % energy_total_vec = zeros(nSteps,1);
        %
        % for i = 1:nSteps
        %     E_step = E_energy_mat(:,1,i);
        %     energy_total_vec(i) = mean(E_step);
        % end

        % cFigure; hold on;
        % plot(lambda_y_vec, energy_total_vec, 'k-', 'LineWidth', 5);
        %
        % xlabel('\lambda_y (Axial Stretch)', 'Interpreter', 'tex', ...
        %     'FontSize', 20,'FontWeight', 'bold');
        %
        % ylabel('Strain Energy Density', 'Interpreter', 'tex', ...
        %     'FontSize', 20, 'FontWeight', 'bold');
        %
        % title(sprintf('Strain Energy vs Stretch  |  Fiber Angle = %d°', ang), ...
        %     'FontSize',20,'FontWeight','bold');
        %
        % set(gca, 'FontSize', 16,'FontWeight', 'bold');
        % xlim([1 max(lambda_y_vec)]);
        % axis square;
        % box on; grid on;

        E_stress_mat_VM=sqrt(( (E_stress_mat(:,1,:)-E_stress_mat(:,2,:)).^2 + ...
            (E_stress_mat(:,2,:)-E_stress_mat(:,3,:)).^2 + ...
            (E_stress_mat(:,1,:)-E_stress_mat(:,3,:)).^2  )/2); %Von Mises stress

        %% Compute p90 stress using same rule as first code

        A = patchArea(E, V_DEF(:,:,1));
        A_tot = sum(A);

        p90_vec = zeros(size(N_disp_mat,3),1);

        for qt = 1:size(N_disp_mat,3)

            sigmaVM = E_stress_mat_VM(:,:,qt);

            % Same rule as first code
            Q = (sigmaVM(:) .* A) ./ A_tot;

            p90_vec(qt) = prctile(Q,90);

        end

        cFigure; hold on;
        plot(lambda_y_vec,p90_vec,'r-o','LineWidth',3,'MarkerFaceColor','r');
        xlabel('\lambda_y');
        ylabel('p90 stress');
        title(sprintf('p90 stress vs stretch | %s',angLabel));
        axis square; grid on; box on;

        % ===== DAMAGE CONTOUR =====

        CVd = faceToVertexMeasure(E, VE, E_damage_mat(:,:,end));

        hf = cFigure;
        title(sprintf('Damage contour | Fiber Direction = %s', angLabel), 'FontSize', 20)

        hp = gpatch(Fb, V_DEF(:,:,end), CVd, 'none', 1, lineWidth1);

        for qp = 1:numel(hp)
            hp(qp).FaceColor = 'interp';
        end

        axisGeom(gca,fontSize);
        colormap(parula);
        colorbar;

        Dmax = max(E_damage_mat(:),[],'omitnan');

        if isempty(Dmax) || ~isfinite(Dmax) || Dmax <= 0
            Dmax = 1;   % fallback when no damage developed
        end

        caxis([0 Dmax]);

        view(140,30);
        camlight headlight;

        % ===== DAMAGE vs STRETCH =====

        nSteps = size(E_energy_mat,3);

        energy_total_vec   = zeros(nSteps,1);
        energy_damage_vec  = zeros(nSteps,1);
        energy_intact_vec  = zeros(nSteps,1);
        damage_vec         = zeros(nSteps,1);

        for i = 1:nSteps

            E_step = E_energy_mat(:,1,i);
            D_step = E_damage_mat(:,1,i);

            % --- Split regions ---
            damaged = D_step > 0.01;   % threshold (important!)
            intact  = ~damaged;

            % --- Total energy ---
            energy_total_vec(i) = mean(E_step);

            % --- Damaged region energy ---
            if any(damaged)
                energy_damage_vec(i) = mean(E_step(damaged));
            else
                energy_damage_vec(i) = NaN;
            end

            % --- Intact region energy ---
            if any(intact)
                energy_intact_vec(i) = mean(E_step(intact));
            else
                energy_intact_vec(i) = NaN;
            end

            % --- Average damage ---
            damage_vec(i) = mean(D_step);

        end

        energy_total_mean   = energy_total_vec;
        energy_damage_mean  = energy_damage_vec;
        energy_intact_mean  = energy_intact_vec;

        energy_total_sum   = zeros(nSteps,1);
        energy_damage_sum  = zeros(nSteps,1);
        energy_intact_sum  = zeros(nSteps,1);

        for i = 1:nSteps

            E_step = E_energy_mat(:,1,i);
            D_step = E_damage_mat(:,1,i);

            damaged = D_step > 0.01;
            intact  = ~damaged;

            energy_damage_sum(i) = sum(E_step(damaged));
            energy_intact_sum(i) = sum(E_step(intact));

            energy_total_sum(i)  = energy_damage_sum(i) + energy_intact_sum(i);

        end


        %% Energy split vs Axial Stretch

        cFigure; hold on;

        plot(lambda_y_vec, energy_total_mean,  'k-',  'LineWidth', 5);
        plot(lambda_y_vec, energy_damage_mean, 'r-',  'LineWidth', 5);
        plot(lambda_y_vec, energy_intact_mean, 'b--', 'LineWidth', 5);

        legend({'Total (mean)','Damaged (mean)','Intact (mean)'}, 'FontSize',16);

        xlabel('\lambda_y (Axial Stretch)', 'FontSize',20,'FontWeight','bold');
        ylabel('Mean Strain Energy Density', 'FontSize',20,'FontWeight','bold');

        title(sprintf('Mean strain energy density (damaged vs intact regions) | %s', angLabel), ...
            'FontSize',20,'FontWeight','bold');

        set(gca,'FontSize',16,'FontWeight','bold');
        xlim([1 max(lambda_y_vec)]);
        axis square; box on; grid on;

        cFigure; hold on;

        plot(lambda_y_vec, energy_total_sum,  'k-',  'LineWidth', 5);
        plot(lambda_y_vec, energy_damage_sum, 'r-',  'LineWidth', 5);
        plot(lambda_y_vec, energy_intact_sum, 'b--', 'LineWidth', 5);

        legend({'Total','Damaged','Intact'}, 'FontSize',16);

        xlabel('\lambda_y (Axial Stretch)', 'FontSize',20,'FontWeight','bold');
        ylabel('Total Strain Energy', 'FontSize',20,'FontWeight','bold');

        title(sprintf('Energy decomposition (global) | %s', angLabel), ...
            'FontSize',20,'FontWeight','bold');

        set(gca,'FontSize',16,'FontWeight','bold');
        xlim([1 max(lambda_y_vec)]);
        axis square; box on; grid on;

        % ===== DAMAGE ANIMATION =====

        hf = cFigure;
        title(sprintf('Damage evolution | Fiber Direction = %s', angLabel), 'FontSize', 20)

        CVd = faceToVertexMeasure(E, VE, E_damage_mat(:,:,end));

        hp = gpatch(Fb, V_DEF(:,:,end), CVd, 'none', 1, lineWidth1);

        for qp = 1:numel(hp)
            hp(qp).FaceColor = 'interp';
        end

        axisGeom(gca,fontSize);
        colormap(parula);
        colorbar;

        Dmax = max(E_damage_mat(:),[],'omitnan');

        if isempty(Dmax) || ~isfinite(Dmax) || Dmax <= 0
            Dmax = 1;
        end

        caxis([0 Dmax]);

        view(140,30);
        camlight headlight;

        % Animation
        animStruct = [];
        animStruct.Time = timeVec;

        for qt = 1:size(N_disp_mat,3)

            CVd = faceToVertexMeasure(E, VE, E_damage_mat(:,:,qt));

            animStruct.Handles{qt} = hp;
            animStruct.Props{qt}   = {'Vertices','CData'};
            animStruct.Set{qt}     = {V_DEF(:,:,qt),CVd};

        end

        anim8(hf,animStruct);
        drawnow;

        % %%
        % p90_vec = zeros(size(N_disp_mat,3),1);
        %
        % for qt = 1:size(N_disp_mat,3)
        %     sigmaVM = E_stress_mat_VM(:,:,qt);
        %     p90_vec(qt) = prctile(sigmaVM(:),90);
        % end


        %% ================== STRESS ANIMATION ==================

        % Initial stress field
        CV = faceToVertexMeasure(E,VE,E_stress_mat_VM(:,:,end));


        hf = cFigure;
        gtitle([febioFebFileNamePart,': Press play to animate']);
        title(sprintf('$\\sigma_{vm}$ [MPa] | Fiber Direction = %s', angLabel), ...
            'Interpreter','latex','FontSize',fontSize);

        hp = gpatch(Fb,V_DEF(:,:,end),CV,'none',1,lineWidth1);

        for qp = 1:numel(hp)
            hp(qp).FaceColor = 'interp';
        end

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([0 max(E_stress_mat_VM(:))]);
        axis(axisLim(V_DEF));
        view(140,30);
        camlight headlight;

        % Animation
        animStruct = [];
        animStruct.Time = timeVec;

        for qt = 1:size(N_disp_mat,3)

            sigmaVM = E_stress_mat_VM(:,:,qt);
            CV = faceToVertexMeasure(E,VE,sigmaVM);


            animStruct.Handles{qt} = hp;
            animStruct.Props{qt}   = {'Vertices','CData'};
            animStruct.Set{qt}     = {V_DEF(:,:,qt),CV};
        end

        anim8(hf,animStruct);
        drawnow;


        %% ================== P90 CLIPPED STRESS ==================

        hf = cFigure;
        title(sprintf('$\\sigma_{vm}$ ($p_{90}$) [MPa] | Fiber Direction = %s', angLabel), ...
            'Interpreter','latex','FontSize',18);

        p90_plot_max = max(p90_vec);

        CV = faceToVertexMeasure(E,VE,E_stress_mat_VM(:,:,end));
        CV = min(CV, p90_vec(end));

        hp = gpatch(Fb, V_DEF(:,:,end), CV, 'none', 1, lineWidth1);

        for qp = 1:numel(hp)
            hp(qp).FaceColor = 'interp';
        end

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;
        caxis([0 p90_plot_max]);
        view(140,30);
        camlight headlight;

        % Animation
        animStruct = [];
        animStruct.Time = timeVec;

        for qt = 1:size(N_disp_mat,3)

            sigmaVM = E_stress_mat_VM(:,:,qt);
            CV = faceToVertexMeasure(E,VE,sigmaVM);

            % Apply percentile clipping
            CV = min(CV, p90_vec(qt));

            animStruct.Handles{qt} = hp;
            animStruct.Props{qt}   = {'Vertices','CData'};
            animStruct.Set{qt}     = {V_DEF(:,:,qt),CV};
        end

        anim8(hf,animStruct);
        drawnow;


        %% ================== DISPLACEMENT ANIMATION ==================

        % Initial displacement field
        CV = abs(N_disp_mat(:,1,end)); % x-displacement

        hf = cFigure;
        gtitle([febioFebFileNamePart,': Press play to animate']);
        title(sprintf('$u_x$ [mm] | Fiber Direction = %s', angLabel), ...
            'Interpreter','latex','FontSize',fontSize);

        hp = gpatch(Fb,V_DEF(:,:,end),CV,'none',1,lineWidth1);

        for qp = 1:numel(hp)
            hp(qp).FaceColor = 'interp';
        end

        axisGeom(gca,fontSize);
        colormap(cMap); colorbar;

        % ✅ Correct displacement color scaling
        caxis([0 max(abs(N_disp_mat(:,1,:)),[],'all')]);

        axis(axisLim(V_DEF));
        view(140,30);
        camlight headlight;

        % Animation
        animStruct = [];
        animStruct.Time = timeVec;

        for qt = 1:size(N_disp_mat,3)

            CV = abs(N_disp_mat(:,1,qt)); % displacement ONLY

            animStruct.Handles{qt} = hp;
            animStruct.Props{qt}   = {'Vertices','CData'};
            animStruct.Set{qt}     = {V_DEF(:,:,qt),CV};
        end

        anim8(hf,animStruct);
        drawnow;
        results(ia).lambda = lambda_y_vec;
        results(ia).force  = Fy_vec;
        results(ia).energy_mean = energy_total_mean;
        results(ia).energy_sum  = energy_total_sum;
        results(ia).damage = damage_vec;


    else
        warning('FEBio failed for angle %d', ang);
        disp(results(1).lambda)
        disp(results(2).lambda)

    end

end

plotComparison(results, 'force', ...
    'Reaction Force F_y', ...
    'Force vs Axial Stretch (Perpendicular vs Parallel)');

plotComparison(results, 'energy_sum', ...
    'Total Strain Energy', ...
    'Energy vs Axial Stretch (Perpendicular vs Parallel)');

plotComparison(results, 'damage', ...
    'Damage', ...
    'Damage vs Axial Stretch (Perpendicular vs Parallel)');

%% =========================================================
% Helper functions
%% =========================================================
function Vr = rotatePolygonAboutOwnCenter(V,angDeg)
c = mean(V,1);
R = [cosd(angDeg) -sind(angDeg); ...
    sind(angDeg)  cosd(angDeg)];
Vr = (R*(V-c)')' + c;
end

function A = signedArea2D(x,y)
x = x(:);
y = y(:);
A = 0.5 * sum(x.*y([2:end 1]) - y.*x([2:end 1]));
end

function plotComparison(results, fieldName, yLabelText, titleText)

cFigure; hold on;

plot(results(1).lambda, results(1).(fieldName), ...
    '-', 'Color', [0.85 0.15 0.15], 'LineWidth', 5);

plot(results(2).lambda, results(2).(fieldName), ...
    '-', 'Color', [0.15 0.35 0.85], 'LineWidth', 5);

xlabel('\lambda_y (Axial Stretch)', ...
    'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');

ylabel(yLabelText, ...
    'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');

title(titleText, ...
    'Interpreter', 'tex', 'FontSize', 20, 'FontWeight', 'bold');

legend({'Perpendicular','Parallel'}, ...
    'FontSize', 16, 'Location', 'best');

set(gca, 'FontSize', 16, 'FontWeight', 'bold');

lambdaAll = [results(1).lambda(:); results(2).lambda(:)];
lambdaAll = lambdaAll(isfinite(lambdaAll));
if isempty(lambdaAll)
    xlim([1 1.01]);
else
    xlim([1 max([lambdaAll; 1.01])]);
end

axis square;
box on;
grid on;

end
