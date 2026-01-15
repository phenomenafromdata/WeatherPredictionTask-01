% Run Experiment
rng('shuffle')
%mac desktop
tasks_path='/Users/neuroudp/Google Drive/LabCode/Magister_NeuroSc/Aspe';

tasks_path='/Users/daniel.rojas/Google Drive/LabCode/Magister_NeuroSc/Aspe';

participantInfo=GetSubjectInfo_Weather();

cd([tasks_path,filesep,'WeatherPredictionTask_UDP']);

WeatherPredictionTask_UDP(participantInfo);