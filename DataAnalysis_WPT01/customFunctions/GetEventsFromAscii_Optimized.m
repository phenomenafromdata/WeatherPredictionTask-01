function time_events=GetEventsFromAscii_Optimized(filename, eventnames)

% reads eyelink ascii file and pulls events timestamps


% May 2023: optimized from GetEventsFromAscii
% optimization involved 
% 1) Removing all uses of 'eval'
% 2) Replacing the function 'regexp' by 'strfind'

% by DRL, 2018




fid = fopen(filename);

delimiter = '\r';   % \r means return key
formatSpec = '%s';
out = textscan(fid, formatSpec, 'Delimiter', delimiter);
out = out{1};
fclose(fid);

eventData = cell(numel(eventnames), 1);
eventCount = zeros(numel(eventnames), 1);

for ii = 1:numel(out) % Loop through file rows
    for jj = 1:numel(eventnames) % Loop through events of interest
        line = out{ii};
        eventIdx = strfind(line, eventnames{jj});
        if ~isempty(eventIdx) % If event name is found in the line
            eventCount(jj) = eventCount(jj) + 1;
            eventData{jj}(eventCount(jj), 1) = str2double(line(4:eventIdx-1)); % First 3 letters are 'MSG' so remove them
            eventData{jj}(eventCount(jj), 2) = jj;
            
            if strcmpi(eventnames{jj}, 'START') 
                leftIdx = strfind(line, 'LEFT');
                if ~isempty(leftIdx)
                    t = line(eventIdx+numel(eventnames{jj}):leftIdx-1);
                    eventData{jj}(eventCount(jj), 1) = str2double(t);
                end
            elseif strcmpi(eventnames{jj}, 'END')
                samplesIdx = strfind(line, 'SAMPLES');
                if ~isempty(samplesIdx)
                    t = line(eventIdx+numel(eventnames{jj}):samplesIdx-1);
                    eventData{jj}(eventCount(jj), 1) = str2double(t);
                end
            end
        end
    end
end

time_events=cell2mat(eventData);
