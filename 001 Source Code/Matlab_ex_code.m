% 시간 설정 
sc = satelliteScenario;
sc.StartTime = datetime(2024,4,29,09,01,57); % In datetime
sc.StopTime = sc.StartTime + hours(10);       % In datetime
sc.SampleTime = 60;                          % In s

% 3가지 위성 궤도 설정
constellation = struct;
constellation(1).Name = "Walker Star";
constellation(1).ShellNames = ["Altitude 1", "Altitude 2"]; % 2개의 쉘 구성
% 2개의 쉘로 구성함으로서 다양한 커버리지를 동시에 확보 가능 
% 관심 지역을 대전이 아닌 극 지방 혹은 지구 전체로 하면 Mixed 방식의 경우 delta 방식보다 star방식과 더 유사하게
% 출력될 것 으로 예상됨 
constellation(1).ShellAltitude = [1015 1325]; % 각 쉘의 고도
constellation(1).ShellInclination = [98.98 50.88]; % 각 쉘의 경사 각 | 전 지구 커버리지 / 중위도 지역 위주
constellation(1).NumOrbits = [6 6]; % 2 2  33 / 각 쉘의 궤도 수 
constellation(1).NumSatellitesPerOrbit = [8 8]; % 22 44 / 각 궤도당 위성의 수
constellation(1).Type = 1;

constellation(2).Name = "Walker Delta";
constellation(2).ShellNames = ["Altitude 1", "Altitude 2"];
constellation(2).ShellAltitude = [1015 1325];
constellation(2).ShellInclination = [98.98 50.88];
constellation(2).NumOrbits = [6 6];
constellation(2).NumSatellitesPerOrbit = [8 8];
constellation(2).Type = 2;

constellation(3).Name = "Mixed Geometry";
constellation(3).ShellNames = ["Altitude 1", "Altitude 2"];
constellation(3).ShellAltitude = [1015 1325];
constellation(3).ShellInclination = [98.98 50.88];
constellation(3).NumOrbits = [6 6];
constellation(3).NumSatellitesPerOrbit = [8 8];
constellation(3).Type = 3;

% 송신기 설정
txConfig = struct;
txConfig.Frequency = 2e9;
txConfig.Power = 20;
txConfig.BitRate = 10;
txConfig.SystemLoss = 0;
txConfig.Bandwidth = 10e6;

% 관심지역 설정 (대전 중심)
lat = [36.23, 36.49];     % 대전 남북 방향
lon = [127.21, 127.49];   % 대전 동서 방향
numUEs = 200; % 사용자 장비 수 
minElevAngle = 30;

% 수신기 설정
rxConfig = struct;
rxConfig.MaxGByT = -5;
rxConfig.SystemLoss = 0;
rxConfig.PreReceiverLoss = 0;
rxConfig.RequiredEbNo = 11;

[latSpacing,lonSpacing] = findClosestFactors(numUEs);
latPts = linspace(lat(1),lat(2),latSpacing);
lonPts = linspace(lon(1),lon(2),lonSpacing);
[latMesh,lonMesh] = meshgrid(latPts,lonPts);
latCoord = latMesh(:);
lonCoord = lonMesh(:);
ue = groundStation(sc,latCoord,lonCoord);
[ue.MinElevationAngle] = minElevAngle;
rx = receiver(ue,Antenna=arrayConfig(Size=[1 1]), SystemLoss=rxConfig.SystemLoss, ...
    PreReceiverLoss=rxConfig.PreReceiverLoss, GainToNoiseTemperatureRatio=rxConfig.MaxGByT);

s = satelliteScenarioViewer(sc,ShowDetails=false); delete(s)
timeSteps = sc.StartTime:seconds(sc.SampleTime):sc.StopTime;
numTimeSteps = numel(timeSteps);

resultPerConstellation = struct("SatelliteVisibility",[], "LinkAvailability",[], "AvgLinkAvailability",[], ...
    "AvgSatelliteVisibility",[], "MaxEbNo",[], "MaxCNR",[], "Capacity",[], "NumTotalSatellites",[]);
results = repmat(resultPerConstellation,numel(constellation),1);

for constIdx = 1:numel(constellation)
    fprintf("Computing coverage statistics for constellation " + constIdx)
    [linkAvailability,satelliteVisibility,maxEbNo] = deal(zeros(numUEs,numTimeSteps));
    sat = addSatellites(sc,constellation(constIdx));
    tx = transmitter(sat, Frequency=txConfig.Frequency, Power=txConfig.Power, ...
        SystemLoss=txConfig.SystemLoss, BitRate=txConfig.BitRate, Antenna=arrayConfig(Size=[1 1]));

    for ueIdx = 1:numUEs
        if mod(ueIdx,ceil(0.2*numUEs)) == 0
            fprintf(".")
        end
        [~,el] = aer(ue(ueIdx),sat);
        elIdx = el >= ue(ueIdx).MinElevationAngle;
        satelliteVisibility(ueIdx,:) = any(elIdx,1);
        validTx = tx(any(elIdx,2));
        if ~isempty(validTx)
            links = link(validTx,rx(ueIdx));
            ebByNo = ebno(links);
            maxEbNo(ueIdx,:) = max(ebByNo);
            linkAvailability(ueIdx,:) = any(linkStatus(links));
            delete(links)
        end
    end

    cno = maxEbNo + pow2db(tx(1).BitRate) + 60;
    cnrdB = cno - pow2db(txConfig.Bandwidth);
    cnr = db2pow(cnrdB);
    capacity = txConfig.Bandwidth*log2(1+cnr);

    results(constIdx).SatelliteVisibility = satelliteVisibility;
    results(constIdx).AvgSatelliteVisibility = mean(mean(satelliteVisibility,2));
    results(constIdx).LinkAvailability = linkAvailability;
    results(constIdx).AvgLinkAvailability = mean(mean(linkAvailability,2));
    results(constIdx).MaxEbNo = maxEbNo;
    results(constIdx).MaxCNR = cnr;
    results(constIdx).Capacity = capacity;
    results(constIdx).NumTotalSatellites = numel(sat);
    delete(sat)
    fprintf(newline + "Computed coverage statistics for constellation " + constIdx + newline)
end

satelliteVisibility = [results.AvgSatelliteVisibility]'*100;
linkAvailability = [results.AvgLinkAvailability]'*100;
coveragePercent = zeros(numel(results),1);
threshold = 1;
for i = 1:numel(results)
    cTemp = results(i).Capacity;
    coveragePercent(i) = mean(cTemp > threshold*1e6,[1 2])*100;
end

numTotalSatellites = [results.NumTotalSatellites];
constellationNames = string({constellation.Name});
resultTable = table(constellationNames',numTotalSatellites',satelliteVisibility,linkAvailability,coveragePercent, ...
    'VariableNames',{'Constellation_Name','Num_Satellites','Satellite_Visibility','Link_Availability','Capacity_Over_1Mbps'});
disp(resultTable)
%% 
% 성능 비교 그래프 출력
figure;
data = [coveragePercent, linkAvailability, satelliteVisibility];
bar(categorical(constellationNames), data);
legend("Capacity > 1Mbps (%)", "Link Availability (%)", "Satellite Visibility (%)");
ylabel("Coverage Metric (%)");
title("LEO 궤도 방식별 커버리지 성능 비교");
grid on
colormap(turbo);

function [factor1,factor2] = findClosestFactors(n)
for i = round(sqrt(n)):-1:1
    if mod(n,i) == 0
        factor1 = i; factor2 = n/i; break
    end
end
end


function sat = addSatellites(sc, const)
% Add the satellites to the scenario based on constellation name

% 궤도 방식 선택 
switch const.Type
    case 1
        orbitFunc = @getOrbitalElements_Star;
    case 2
        orbitFunc = @getOrbitalElements;
    case 3
        orbitFunc = @getOrbitalElements_Mixed;
    otherwise
        error("Unknown constellation name: " + const.Name)
end

% 위성 추가
sat = [];
for i = 1:numel(const.ShellAltitude)
    [a,e,iang,raan,aop,ta,n] = orbitFunc( ...
        const.ShellAltitude(i), const.ShellInclination(i), ...
        const.NumOrbits(i), const.NumSatellitesPerOrbit(i));

    satShell = satellite(sc, a,e,iang,raan,aop,ta, ...
        Name=(const.Name + " " + const.ShellNames(i) + " " + string(1:n)'));
    
    sat = [sat(:); satShell(:)];
end

end

% walkerdelta
function [semimajoraxis,eccentricity, ...
    inclination,RAAN,argofperiapsis, ...
    trueanomaly,numTotSat] = getOrbitalElements(shellAltitude,shellInclination, ...
    numOrbits,numSatellitesPerPlane)
% Get the orbital elements based on the shell altitude, shell inclination,
% number of orbits and number of satellites per orbital plane

% Number of satellites in the constellation
numTotSat = numOrbits*numSatellitesPerPlane;

% Get orbital elements
orbitIdx = repelem(1:numOrbits,1,numSatellitesPerPlane);
planeIdx = repmat(1:numSatellitesPerPlane,1,numOrbits);
RAAN = 180*(orbitIdx-1)/numOrbits;
trueanomaly = 360*(planeIdx-1 + 0.5*(mod(orbitIdx,2)-1)) ...
    /numSatellitesPerPlane;
semimajoraxis = repmat((6371 + shellAltitude)*1e3,1,numTotSat); % meters
inclination = repmat(shellInclination,1,numTotSat);             % degrees
eccentricity = zeros(1,numTotSat);                              % degrees
argofperiapsis = zeros(1,numTotSat);                            % degrees

end

% walkerstar
function [semimajoraxis,eccentricity, ...
    inclination,RAAN,argofperiapsis, ...
    trueanomaly,numTotSat] = getOrbitalElements_Star(shellAltitude,shellInclination, ...
    numOrbits,numSatellitesPerPlane)

% 총 위성 수
numTotSat = numOrbits * numSatellitesPerPlane;

% 각 위성 인덱스 계산
orbitIdx = repelem(1:numOrbits, 1, numSatellitesPerPlane);
satIdx = repmat(1:numSatellitesPerPlane, 1, numOrbits);

% Walker Star 방식 특징:
% 1) 모든 위성의 true anomaly는 동일 간격
% 2) RAAN은 0~360도 전체에 균일하게 분포

RAAN = 360 * (orbitIdx - 1) / numOrbits;
trueanomaly = 360 * (satIdx - 1) / numSatellitesPerPlane;

% 기타 궤도 요소
semimajoraxis = repmat((6371 + shellAltitude) * 1e3, 1, numTotSat); % m
inclination = repmat(shellInclination, 1, numTotSat);               % deg
eccentricity = zeros(1, numTotSat);                                 % 원형 궤도
argofperiapsis = zeros(1, numTotSat);

end

%mixed
function [semimajoraxis, eccentricity, inclination, ...
          RAAN, argofperiapsis, trueanomaly, numTotSat] = ...
          getOrbitalElements_Mixed(shellAltitude, shellInclination, ...
                                   numOrbits, numSatellitesPerPlane)
% 각 shell마다 walker 방식 다르게 설정
% Shell 1: Delta 방식
% Shell 2: Star 방식
% 필요시 shell 개수 늘리면 확장 가능

nShells = numel(shellAltitude);
semimajoraxis = [];
eccentricity = [];
inclination = [];
RAAN = [];
argofperiapsis = [];
trueanomaly = [];
numTotSat = 0;

for i = 1:nShells
    if i == 1
        % Delta 방식
        [a, e, inc, raan, aop, ta, n] = getOrbitalElements( ...
            shellAltitude(i), shellInclination(i), ...
            numOrbits(i), numSatellitesPerPlane(i));
    else
        % Star 방식
        [a, e, inc, raan, aop, ta, n] = getOrbitalElements_Star( ...
            shellAltitude(i), shellInclination(i), ...
            numOrbits(i), numSatellitesPerPlane(i));
    end

    semimajoraxis = [semimajoraxis, a];
    eccentricity = [eccentricity, e];
    inclination = [inclination, inc];
    RAAN = [RAAN, raan];
    argofperiapsis = [argofperiapsis, aop];
    trueanomaly = [trueanomaly, ta];
    numTotSat = numTotSat + n;
end

end


function validateConstellation(constellation)
% Validates the dimensions of all fields in each constellation

for constIdx = 1:numel(constellation)
    validateattributes(constellation(constIdx).Name,["char","string"],"scalartext")
    numShellName = numel(constellation(constIdx).ShellNames);
    validateattributes(constellation(constIdx).ShellNames,"string","nonempty")
    validateattributes(constellation(constIdx).ShellAltitude,"double", ...
        {'positive','integer','numel',numShellName})
    validateattributes(constellation(constIdx).ShellInclination,"double", ...
        {'numel',numShellName})
    validateattributes(constellation(constIdx).NumOrbits,"double", ...
        {'positive','integer','numel',numShellName})
    validateattributes(constellation(constIdx).NumSatellitesPerOrbit,"double", ...
        {'positive','integer','numel',numShellName})
end

end