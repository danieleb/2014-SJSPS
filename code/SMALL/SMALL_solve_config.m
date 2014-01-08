%% Configuration file used in SMALL_solve
%
%   Please DO NOT use this file to change the solvers used in SMALLBox
%   If you want to change the solvers create a copy 
%     of this file named 'SMALL_learn_config_local.m'
%
%   Please refer to the documentation for further information

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


if strcmpi(solver.toolbox,'sparselab')
    y = eval([solver.name,'(SparseLab_A, b, n,',solver.param,');']);
elseif strcmpi(solver.toolbox,'sparsify')
    if isa(Problem.A,'float')
        y = eval([solver.name,'(b, A, n,',solver.param,');']);
    else
        y = eval([solver.name,'(b, A, n, ''P_trans'', AT,',solver.param,');']);
    end
elseif (strcmpi(solver.toolbox,'spgl1')||strcmpi(solver.toolbox,'gpsr'))
    y = eval([solver.name,'(b, A,',solver.param,');']);
elseif (strcmpi(solver.toolbox,'SPAMS'))
    y = eval([solver.name,'(b, A, solver.param);']);
elseif (strcmpi(solver.toolbox,'SMALL'))
    if isa(Problem.A,'float')
        y = eval([solver.name,'(A, b, n,',solver.param,');']);
    else
        y = eval([solver.name,'(A, b, n,',solver.param,',AT);']);
    end
elseif (strcmpi(solver.toolbox, 'ompbox'))
    G=A'*A;
    epsilon=solver.param.epsilon;
    maxatoms=solver.param.maxatoms;
    y = eval([solver.name,'(A, b, G,epsilon,''maxatoms'',maxatoms,''checkdict'',''off'');']);
    % danieleb: added call to omp functions with fast implementation.
elseif (strcmpi(solver.toolbox, 'ompbox_fast'))
    DtX=A'*b;
    XtX = sum(b.*b);
    G=A'*A;
    epsilon=solver.param.epsilon;
    maxatoms=solver.param.maxatoms;
    y = eval([solver.name,'(DtX, XtX, G,epsilon,''maxatoms'',maxatoms,''checkdict'',''off'');']);
elseif (strcmpi(solver.toolbox, 'ompsbox'))
    basedict = Problem.basedict;
    if issparse(Problem.A)
        A = Problem.A;
    else
        A = sparse(Problem.A);
    end
    G = dicttsep(basedict,A,dictsep(basedict,A,speye(size(A,2))));
    epsilon=solver.param.epsilon;
    maxatoms=solver.param.maxatoms;
    y = eval([solver.name,'(basedict, A, b, G,epsilon,''maxatoms'',maxatoms,''checkdict'',''off'');']);
    Problem.sparse=1;
elseif (strcmpi(solver.toolbox, 'ALPS'))
    if ~isa(Problem.A,'float')
        % ALPS does not accept implicit dictionary definition
        A = opToMatrix(Problem.A, 1);
    end
    [y, numiter, time, y_path] = wrapper_ALPS_toolbox(b, A, solver.param);
elseif (strcmpi(solver.toolbox, 'MMbox'))
    if ~isa(Problem.A,'float')
        % MMbox does not accept implicit dictionary definition
        A = opToMatrix(Problem.A, 1);
    end
    
    [y, cost] = wrapper_mm_solver(b, A, solver.param);

elseif (strcmpi(solver.toolbox, 'UNLocBox'))
    if ~isa(Problem.A,'float')
        % MMbox does not accept implicit dictionary definition
        A = opToMatrix(Problem.A, 1);
    end
    
    y = unloc_solver(b, A, solver.param,solver.name);
    
    %%
    %   Please do not make any changes to the 'SMALL_solve_config.m' file
    %   All the changes should be done to your local configuration file
    %    named 'SMALL_solve_config_local.m'
    %
    %   To introduce new sparse representation algorithm put the files in
    %   your Matlab path. Next, unique name <TolboxID> for your toolbox and
    %   prefferd API <Preffered_API> needs to be defined.
    %
    % elseif strcmpi(solver.toolbox,'<ToolboxID>')
    %
    %     % - Evaluate the function (solver.name - defined in the main) with
    %     %   parameters given above
    %
    %     y = eval([solver.name,'(<Preffered_API>);']);
    
else
    printf('\nToolbox has not been registered. Please change SMALL_solver file.\n');
    return
end
