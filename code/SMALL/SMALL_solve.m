function solver = SMALL_solve(Problem, solver)
%% SMALL sparse solver caller function
%
%   Function gets as input SMALL structure that contains SPARCO problem to
%   be solved, name of the toolbox and solver, and parameters file for
%   particular solver.
%
%   Outputs are solution, reconstructed signal and time spent

%   Centre for Digital Music, Queen Mary, University of London.
%   This file copyright 2009 Ivan Damnjanovic.
%
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License as
%   published by the Free Software Foundation; either version 2 of the
%   License, or (at your option) any later version.  See the file
%   COPYING included with this distribution for more information.
%
%%

% SMALLboxInit

if isa(Problem.A,'float')
    A = Problem.A;
    SparseLab_A=Problem.A;
    m = size(Problem.A,1);      % m is the no. of rows.
    n = size(Problem.A,2);      % n is the no. of columns.
else
    A  = @(x) Problem.A(x,1);  % The operator
    AT = @(y) Problem.A(y,2);  % and its transpose.
    SparseLab_A =@(mode, m, n, x, I, dim) SL_A(Problem.A, mode, m, n, x, I, dim);
    m = Problem.sizeA(1);      % m is the no. of rows.
    n = Problem.sizeA(2);      % n is the no. of columns.
    
end
% if signal that needs to be represented is different then training set for
% dictionary learning it should be stored in Problem.b1 matix
if isfield(Problem, 'b1')
    b = Problem.b1;
else
    b = Problem.b;            % The right-hand-side vector.
end
%%
if (solver.profile)
    fprintf('\nStarting solver %s... \n', solver.name);
end

start=cputime;
tStart=tic;

%% solvers configuration
% test if there is a locally modified version of the config
% otherwise reads the "default" config file
run('SMALL_solve_config.m');

%%
%   Sparse representation time
tElapsed=toc(tStart);
solver.time = cputime - start;
if (solver.profile)
    fprintf('Solver %s finished task in %2f seconds (cpu time). \n', solver.name, solver.time);
    fprintf('Solver %s finished task in %2f seconds (tic-toc time). \n', solver.name, tElapsed);
end
solver.time=tElapsed;
% geting around out of memory problem when converting big matrix from
% sparse to full...

if isfield(Problem, 'sparse')&&(Problem.sparse==1)
    solver.solution = y;
else
    solver.solution = full(y);
end
if isfield(Problem,'reconstruct')
    %   Reconstruct the signal from the solution
    solver.reconstructed  = Problem.reconstruct(solver.solution);
end
end
