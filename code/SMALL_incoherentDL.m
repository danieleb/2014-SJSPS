function DL=SMALL_incoherentDL(Problem, DL)

% determine which solver is used for sparse representation %

solver = DL.param.solver;

% determine which type of udate to use ('KSVD', 'MOD', 'ols' or 'mailhe') %

typeUpdate = DL.name;

sig = Problem.b;

% determine dictionary size %

if (isfield(DL.param,'initdict'))
    if (any(size(DL.param.initdict)==1) && all(iswhole(DL.param.initdict(:))))
        dictsize = length(DL.param.initdict);
    else
        dictsize = size(DL.param.initdict,2);
    end
end
if (isfield(DL.param,'dictsize'))    % this superceedes the size determined by initdict
    dictsize = DL.param.dictsize;
end

if (size(sig,2) < dictsize)
    error('Number of training signals is smaller than number of atoms to train');
end


% initialize the dictionary %

% todo: check second if statement
if (isfield(DL.param,'initdict')) && ~isempty(DL.param.initdict);
    if (any(size(DL.param.initdict)==1) && all(iswhole(DL.param.initdict(:))))
        dico = sig(:,DL.param.initdict(1:dictsize));
    else
        if (size(DL.param.initdict,1)~=size(sig,1) || size(DL.param.initdict,2)<dictsize)
            error('Invalid initial dictionary');
        end
        dico = DL.param.initdict(:,1:dictsize);
    end
else
    data_ids = find(colnorms_squared(sig) > 1e-6);   % ensure no zero data elements are chosen
    perm = randperm(length(data_ids));
    dico = sig(:,data_ids(perm(1:dictsize)));
end

% flow: 'sequential' or 'parallel'. If sequential, the residual is updated
% after each atom update. If parallel, the residual is only updated once
% the whole dictionary has been computed. Sequential works better, there
% may be no need to implement parallel. Not used with MOD.

if isfield(DL.param,'flow')
    flow =  DL.param.flow;
else
    flow = 'sequential';
end

% learningRate. If the type is 'ols', it is the descent step of
% the gradient (typical choice: 0.1). If the type is 'mailhe', the
% descent step is the optimal step*rho (typical choice: 1, although 2
% or 3 seems to work better). Not used for MOD and KSVD.

if isfield(DL.param,'learningRate')
    learningRate = DL.param.learningRate;
else
    learningRate = 0.1;
end

% number of iterations (default is 40) %

if isfield(DL.param,'nIte')
    nIte = DL.param.nIte;
else
    nIte = 40;
end

% show dictonary every specified number of iterations

if isfield(DL.param,'show_dict')
    show_dictionary=1;
    show_iter=DL.param.show_dict;
else
    show_dictionary=0;
    show_iter=0;
end

% This is a small patch that needs to be resolved in dictionary learning we
% want sparse representation of training set, and in Problem.b1 in this
% version of software we store the signal that needs to be represented
% (for example the whole image)


tmpTraining = Problem.b1;
Problem.b1 = sig;
if isfield(Problem,'reconstruct')
    Problem = rmfield(Problem, 'reconstruct');
end
solver.profile = 0;

% main loop %
Problem.A = normcols(dico);
solver = SMALL_solve(Problem, solver);
% userData = struct('initDic',normcols(dico),'initX',solver.solution,'dicDis',[],'snr',[]);
% set(0,'UserData',userData);
for i = 1:nIte
%     userData.oldX = solver.solution;            %store old X
    
    Problem.A = normcols(dico);
    solver = SMALL_solve(Problem, solver);
    
%     userData.newX = solver.solution;            %store new X
%     
%     userData.oldDic = dictionary(normcols(dico));           %store old dictionary
    %DICTIONARY UPDATE STEP
    if strcmpi(typeUpdate,'mocod')			%if update is MOCOD create parameters structure
        mocodParams = struct('zeta',DL.param.zeta,...	%coherence regularization factor
            'eta',DL.param.eta,...		%atoms norm regularization factor
            'Dprev',dico);				%previous dictionary
        % dico = dico_update(dico,sig,solver.solution,typeUpdate,flow,learningRate,mocodParams);
        if ~isfield(DL.param,'decFcn'), DL.param.decFcn = 'none'; end
        
        dico = dico_update_mocod(dico,sig,solver.solution,typeUpdate,flow,learningRate,mocodParams);
        
    else
        [dico, solver.solution] = dico_update(dico, sig, solver.solution, ...
            typeUpdate, flow, learningRate);
        dico = normcols(dico);
    end
    mu = DL.param.mu;
    DL.param.vis = false;
    switch lower(DL.param.decFcn)
        case 'ink-svd'
            dico = dico_decorr_symetric(dico,mu,solver.solution);
        case 'grassmannian'
            [n, m] = size(dico);
            dico = grassmannian(n,m,[],0.9,0.99,dico);
        case 'shrinkgram'
            dico = shrinkgram(dico,mu);
        case 'iterproj'
            dico = iterativeprojections(dico,mu,Problem.b1,solver.solution);
        case 'iprdis'
            dico = IPRDis(dico,mu,Problem.b,solver.solution,Problem.cat,DL.param);
        otherwise
    end
%     userData.newDic = dictionary(normcols(dico));             %store new dictionary
%     userData.dicDis = [userData.dicDis; distance(userData.oldDic,userData.newDic,solver)];
%     userData.snr = [userData.snr; snr(Problem.b,normcols(dico)*solver.solution)];
%     set(0,'UserData',userData);
    
    %     [dico, solver.solution] = dico_update(dico, sig, solver.solution, ...
    %         typeUpdate, flow, learningRate);
    %     if (decorrelate)
    %         dico = dico_decorr(dico, mu, solver.solution);
    %     end
    
    if ((show_dictionary)&&(mod(i,show_iter)==0))
        dictimg = SMALL_showdict(dico,[8 8],...
            round(sqrt(size(dico,2))),round(sqrt(size(dico,2))),'lines','highcontrast');
        figure(2); imagesc(dictimg);colormap(gray);axis off; axis image;
        pause(0.02);
    end
end

Problem.b1 = tmpTraining;
DL.D = dico;

end

function Y = colnorms_squared(X)

% compute in blocks to conserve memory
Y = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
    blockids = i : min(i+blocksize-1,size(X,2));
    Y(blockids) = sum(X(:,blockids).^2);
end

end
