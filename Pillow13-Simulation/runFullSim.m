clear all;close all
cd 'C:\Users\nikpa_000\SkyDrive\Documents\Stanford\Grad\Research\SpikeSort\src\SpikeOverlapDetect\Pillow13-Simulation\';
seeds = [1,3,5,25,50];
for i = 1:length(seeds)
    cd 'C:\Users\nikpa_000\SkyDrive\Documents\Stanford\Grad\Research\SpikeSort\src\SpikeOverlapDetect\Pillow13-Simulation\';

    seed = seeds(i);
    num_el = 513;
    num_neur = 2;
    filen = ['data',num2str(seed),'_',num2str(num_el),'_',num2str(num_neur)];
    script0_simulateDataForTesting;
    cd('C:\Users\nikpa_000\SkyDrive\Documents\Stanford\Grad\Research\SpikeSort\src\mvision');

    %runPC;
    %generate vision bin file from raw Y data
    generateRDF;
    % dataPath = strcat('C:\Users\nikpa_000\SkyDrive\Documents\Stanford\Grad\Research\SpikeSort\src\SpikeOverlapDetect\Pillow13-Simulation\dat\mvision-dat\',filen,'\data111');
    % config = mVisionConfig();
    % noiseConfig = config.getNoiseConfig();
    % 
    % % Set up parameters and data source
    % time = noiseConfig.time; % Compute noise from 5 secs of data
    % timeToSkip = noiseConfig.timeToSkip; % Skip the first 5 sec of the file
    % 
    % dataSource = DataFileUpsampler(dataPath);
    % samplingRate = dataSource.samplingRate;
    % nSamples = time * samplingRate;
    % 
    % % get data to process noise on
    % dataSource.loadRandomBuffer(timeToSkip * samplingRate + dataSource.startSample, nSamples, false);
    % 
    % 
    % % noise calculation
    % [rmsNoise,~] = getNoise(dataSource.rawData(2:end,:)); % Remove TTL
    % rmsNoise = [0;rmsNoise]; 
    %test;
end


