function DL = SMALL_learn(Problem,DL)
%% SMALL Dictionary Learning
%
%   Function gets as input Problem and Dictionary Learning (DL) structures
%   In Problem structure field b with the training set needs to be defined
%   In DL fields with name of the toolbox and solver, and parameters file
%   for particular dictionary learning technique needs to be present.
%
%   Outputs are Learned dictionary and time spent as a part of DL structure

%
%   Centre for Digital Music, Queen Mary, University of London.
%   This file copyright 2009 Ivan Damnjanovic.
%
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License as
%   published by the Free Software Foundation; either version 2 of the
%   License, or (at your option) any later version.  See the file
%   COPYING included with this distribution for more information.
%%

% SMALLboxInit

if (DL.profile)
    fprintf('\nStarting Dictionary Learning %s... \n', DL.name);
end

start=cputime;
tStart=tic;

%% toolbox configuration
% test if there is a locally modified version of the config
% otherwise reads the "default" config file
run('SMALL_learn_config_local.m');

%%
%   Dictionary Learning time
tElapsed=toc(tStart);
DL.time = cputime - start;
if (DL.profile)
    fprintf('\n%s finished task in %2f seconds (cpu time). \n', DL.name, DL.time);
    fprintf('\n%s finished task in %2f seconds (tic-toc time). \n', DL.name, tElapsed);
end
DL.time=tElapsed;
%   If dictionary is given as a sparse matrix change it to full

DL.D = full(D);

end

