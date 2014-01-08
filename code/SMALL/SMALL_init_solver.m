function solver = SMALL_init_solver(toolbox, name, param, profile)
%%   Function initialise SMALL structure for sparse representation.
%   Optional input variables:
%       toolbox - name of Dictionary Learning toolbox you want to use
%       name    - name of the algorithm from above toolbox
%       param   - parameters you want to set

%
%   Centre for Digital Music, Queen Mary, University of London.
%   This file copyright 2010 Ivan Damnjanovic.
%
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License as
%   published by the Free Software Foundation; either version 2 of the
%   License, or (at your option) any later version.  See the file
%   COPYING included with this distribution for more information.
%
%%

if ~ exist( 'toolbox', 'var' ) || isempty(toolbox) 
    solver.toolbox = []; 
else
    solver.toolbox = toolbox;
end
if ~ exist( 'name', 'var' ) || isempty(name) 
    solver.name = [];
else
    solver.name = name;
end
if ~ exist( 'param', 'var' ) || isempty(param) 
    solver.param = [];
else
    solver.param = param;
end
if ~ exist( 'profile', 'var' ) || isempty(profile) 
    solver.profile = 1;
else
    solver.profile = profile;
end
solver.add_constraints = 0;
solver.solution = [];
solver.reconstructed = [];
solver.time = [];

end