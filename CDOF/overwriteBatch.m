function overwriteBatch()
%OVERWRITEBATCH Overwrites the existing batch file, so that the name of the
%GDF file is automatically provided without requiring user input. Designed
%for use when defmod.exe is used to set up some generalised modes.

warning('overwriteBatch currently reads the batch file as created by WamitRunCondition...')
warning('...when the defmod approach to generalised modes is used, and overwrites it.')

% Set read, write and compile locations
readLoc = [myMatPath '\snl_fbWecCntrlFork_CDOF\CDOF\wamrun\wam_run'];
writeLoc = [myMatPath '\snl_fbWecCntrlFork_CDOF\CDOF\wamrun\wam_run'];

% Read entire file as character array
txt = fileread([readLoc '.bat']);

% Split into lines
txt = regexp(txt, '\r\n', 'split');

% Rewrite these two lines of the file
txt{15} = '(echo n';
txt{16} = ['echo %fN%1) | "C:\wamitv7\defmod"' ];

% Write to new file
fid = fopen([writeLoc '.bat'], 'wt');
fprintf(fid, '%s\n', txt{:});
fclose(fid);

end

