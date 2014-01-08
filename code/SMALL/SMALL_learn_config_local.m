%% Configuration file used in SMALL_learn
%
%   Please DO NOT use this file to change the dictionary learning algorithms in SMALLBox
%   If you want to change the dictionary learning algorithms
%     create a copy of this file named 'SMALL_learn_config_local.m'
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

if strcmpi(DL.toolbox,'KSVD')
    param=DL.param;
    param.data=Problem.b;
    
    D = eval([DL.name,'(param)']);%, ''t'', 5);']);
elseif strcmpi(DL.toolbox,'KSVDS')
    param=DL.param;
    param.data=Problem.b;
    
    D = eval([DL.name,'(param, ''t'', 5);']);
elseif strcmpi(DL.toolbox,'SPAMS')
    
    X  = Problem.b;
    param=DL.param;
    
    D = eval([DL.name,'(X, param);']);
    %   As some versions of SPAMS does not produce unit norm column
    %   dictionaries, we need to make sure that columns are normalised to
    %   unit lenght.
    
    for i = 1: size(D,2)
        D(:,i)=D(:,i)/norm(D(:,i));
    end
elseif strcmpi(DL.toolbox,'SMALL')
    
    X  = Problem.b;
    param=DL.param;
    
    D = eval([DL.name,'(X, param);']);
    %   we need to make sure that columns are normalised to
    %   unit lenght.
    
    for i = 1: size(D,2)
        D(:,i)=D(:,i)/norm(D(:,i));
    end
    
elseif strcmpi(DL.toolbox,'TwoStepDL')
    
    DL=SMALL_two_step_DL(Problem, DL);
    
    %   we need to make sure that columns are normalised to
    %   unit lenght.
    
    for i = 1: size(DL.D,2)
        DL.D(:,i)=DL.D(:,i)/norm(DL.D(:,i));
    end
    D = DL.D;
    
elseif strcmpi(DL.toolbox,'SMALL_incoherentDL')
    
    DL=SMALL_incoherentDL(Problem, DL);
    
    %   we need to make sure that columns are normalised to
    %   unit lenght.
    
    for i = 1: size(DL.D,2)
        DL.D(:,i)=DL.D(:,i)/norm(DL.D(:,i));
    end
    D = DL.D;
    
elseif strcmpi(DL.toolbox,'MMbox')
    
    DL = wrapper_mm_DL(Problem, DL);
    
    %   we need to make sure that columns are normalised to
    %   unit lenght.
    
    for i = 1: size(DL.D,2)
        DL.D(:,i)=DL.D(:,i)/norm(DL.D(:,i));
    end
    D = DL.D;

%%    
%   Please do not make any changes to the 'SMALL_learn_config.m' file
%   All the changes should be done to your local configuration file
%    named 'SMALL_learn_config_local.m'
%
%   To introduce new dictionary learning technique put the files in
%   your Matlab path. Next, unique name <TolboxID> for your toolbox needs
%   to be defined and also prefferd API for toolbox functions <Preffered_API>
%
% elseif strcmpi(DL.toolbox,'<ToolboxID>')
%     % This is an example of API that can be used:
%     % - get training set from Problem part of structure
%     % - assign parameters defined in the main program
%
%     X  = Problem.b;
%     param=DL.param;
%
%     % - Evaluate the function (DL.name - defined in the main) with
%     %   parameters given above
%
%     D = eval([DL.name,'(<Preffered_API>);']);

else
    printf('\nToolbox has not been registered. Please change SMALL_learn file.\n');
    return
end
