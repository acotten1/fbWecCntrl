function writeDefmod(zLoc)
%WRITEDEFMOD Writes and compiles the DEFMOD.F file for use with a single heave CDOF
%generalised mode. The generalised mode is set over all panels at a certain
%z-location (e.g. the bottom of the object).

warning('writeDefmod currently relies upon a premade defmod file from which the structure is read before rewriting the file.')

% Set read, write and compile locations
readLoc = [myMatPath '\snl_fbWecCntrlFork_CDOF\CDOF\defmod_heaveCDOF'];
writeLoc = [myMatPath '\snl_fbWecCntrlFork_CDOF\CDOF\wamrun\newDefmod'];
compileLoc = 'C:\wamitv7\defmod';

% Read entire file as character array
txt = fileread([readLoc '.f']);

% Split into lines
txt = regexp(txt, '\r\n', 'split');

% Rewrite this line of the file using the new panel z-location
txt{495} = ['      IF(ABS(Z+' num2str(zLoc) ').LT.TOL) THEN']; % line 495 contains the z-location of the CDOF panels that will need replacing

% Write to new file
fid = fopen([writeLoc '.f'], 'wt');
fprintf(fid, '%s\n', txt{:});
fclose(fid);

% Compile into and executable and store in the WAMIT directory
command1 = ['gfortran "' writeLoc '.f" -o "' compileLoc '.exe"'];
system(command1);

end

