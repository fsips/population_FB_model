function setDirs

fullpath = mfilename('fullpath');
[dum, fullpath] = strtok( fliplr(fullpath), '/\' );
fullpath = fliplr(fullpath);

addpath(genpath(fullpath));