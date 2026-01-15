% Initializes Eyelink, starts recording and executes calibration routine

% the script requires two variables:
% winPointer: the output of Psychtoolbox's Screen('OpenWindow')
% edfFile: a string indicating the filename of the edf data file (where eye
% data is going to be saved)

% DRL, Jan 4 2019


% disp('Waiting for user signal to start Eyelink...');
% KbStrokeWait(device2useID);
disp('Initializing Eyelink...')


% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).
el=EyelinkInitDefaults(winPointer);

% el.devicenumber=device2useID;


% Initialization of the connection with the Eyelink Gazetracker.
% exit program if this fails.
if ~EyelinkInit(0, 1)
    fprintf('Eyelink Init aborted.\n');
    Eyelink('Shutdown');
    sca;
    return;
else
    [v, vs]=Eyelink('GetTrackerVersion');
    fprintf('EYELINK CONNECTION SUCCESFUL. Running experiment on a ''%s'' tracker.\n', vs );
    disp('*************************')
end


% make sure that we get gaze data from the Eyelink
Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');

% open file to record data to
Eyelink('Openfile', edfFile);

disp('when done calibrating, press ESC to proceed to task')

% Calibrate the eye tracker
EyelinkDoTrackerSetup(el);

% do a final check of calibration using driftcorrection
% EyelinkDoDriftCorrection(el);

% start recording eye position
Eyelink('StartRecording');
% record a few samples before we actually start displaying
WaitSecs(0.1);
% mark zero-plot time in data file
Eyelink('Message', 'SYNCTIME');

% if ~dummyMode_EEG
%     NetStation('Event','EL01');
% end