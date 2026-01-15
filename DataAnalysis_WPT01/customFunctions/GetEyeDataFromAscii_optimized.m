function [EyeDataMatrix, Fs, info]=GetEyeDataFromAscii_optimized(filename)
% Reads eyelink-produced ascii file and imports it into Matlab.

% inputs:
% filename -  a string containg the name of the ascii file

%outputs:
% EyeDataMatrix - a data matrix of dimensions samples x variables
% columns are: 
% 1) time (in ms), 2) gaze position X - EYE1, 3) gaze position Y - EYE1, 
% 4) pupil diameter - EYE1, 5) gaze position X - EYE2, 6) gaze position Y - EYE2, 
% 7) pupil diameter - EYE2, 8) A constant value

% if only one eye was recorded, then there will be less columns


% New version based on GetEyeDataFromAscii
% DRL,May 2023

fid=fopen(filename);
delimiter='^\r';   %\r means return key
formatSpec = '%s';
out=textscan(fid,formatSpec,'Delimiter',delimiter); out=out{1,1};
fclose(fid);

% get sampling frequency, in samples/s
Fs=NaN;
targetString1='RATE';
for ii=1:numel(out)
    currString=out{ii,1};
    %search for sampling rate
    idx=strfind(currString,targetString1);
    if idx>0 %rate found, pull data
        Fs=str2double(currString(idx+numel(targetString1):idx+10));
        break
    end
end


% get the full matrix of eye tracking data

nCols2prealloc=8; %max n of data cols in the matrix

EyeDataMatrix=nan(size(out,1),nCols2prealloc);
current=0;
for ii=1:numel(out)
    currString=out{ii,1};
    
    if isstrprop(currString(1), 'digit')  % check it's a digit, meaning is a data row
        
        tline   = strrep(currString, ' . ', ' NaN '); % replace missing values
        tmp     = sscanf(tline, '%f');
        ncols   = numel(tmp);
        current = current + 1;
        
        EyeDataMatrix(current, 1:ncols) = tmp; 
    end
end
% remove extra rows
EyeDataMatrix=EyeDataMatrix(1:current,:);

%remove all-nan columns
cols2keep=sum(isnan(EyeDataMatrix))~=size(EyeDataMatrix,1);
EyeDataMatrix=EyeDataMatrix(:,cols2keep);

% get file header
info=out(1:10);
