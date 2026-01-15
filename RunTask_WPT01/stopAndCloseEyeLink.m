% stops Eyelink recording, saves data file, and closes eyelink connection

% DRL, Jan 4 2019

% if ~dummyMode_EEG
%     NetStation('Event','EL00');
% end


Eyelink('StopRecording');


Eyelink('CloseFile');

try  
    fprintf('Receiving data file ''%s''\n', edfFile);
    status=Eyelink('ReceiveFile');
    if status > 0
        fprintf('ReceiveFile status %d\n', status);
    end
    if 2==exist(edfFile, 'file')
        fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
    end
catch ME
    fprintf('Problem receiving data file ''%s''\n', edfFile );
    ME;
end


Eyelink('Shutdown');