function [dico, amp] = dico_update(dico, sig, amp, type, flow, rho)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [dico, amp] = dico_update(dico, sig, amp, type, flow, rho)
    %
    % perform one iteration of dictionary update for dictionary learning
    %
    % parameters:
    % - dico: the initial dictionary with atoms as columns
    % - sig: the training data
    % - amp: the amplitude coefficients as a sparse matrix
    % - type: the algorithm can be one of the following
    %   - ols: fixed step gradient descent, as described in Olshausen &
    %   Field95
    %   - opt: optimal step gradient descent, as described in Mailhe et
    %   al.08
    %   - MOD: pseudo-inverse of the coefficients, as described in Engan99
    %   - KSVD: PCA update as described in Aharon06. For fast applications,
    %   use KSVDbox rather than this code.
    %   - LGD: large step gradient descent. Equivalent to 'opt' with
    %          rho=2.
    % - flow: 'sequential' or 'parallel'. If sequential, the residual is
    % updated after each atom update. If parallel, the residual is only
    % updated once the whole dictionary has been computed.
    % Default: Sequential (sequential usually works better). Not used with
    % MOD.
    % - rho: learning rate. If the type is 'ols', it is the descent step of
    % the gradient (default: 0.1). If the type is 'opt', the
    % descent step is the optimal step*rho (default: 1, although 2 works
    % better. See LGD for more details). Not used for MOD and KSVD.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~ exist( 'flow', 'var' ) || isempty(flow)
        flow = 'sequential';
    end
    
    res = sig - dico*amp;
    nb_pattern = size(dico, 2);
    
    % if the type is random, then randomly pick another type
    switch type
        case 'rand'
            x = rand();
            if x < 1/3
                type = 'MOD';
            elseif type < 2/3
                type = 'opt';
            else
                type = 'KSVD';
            end
    end
    
    % set the learning rate to default if not provided
    if ~ exist( 'rho', 'var' ) || isempty(rho)
        switch type
            case 'ols'
                rho = 0.1;
            case 'opt'
                rho = 1;
        end
    end
    
    switch type
        case 'MOD'
            G = amp*amp';
            dico2 = sig*amp'*G^-1;
            for p = 1:nb_pattern
                n = norm(dico2(:,p));
                % renormalize
                if n > 0
                    dico(:,p) = dico2(:,p)/n;
                    amp(p,:) = amp(p,:)*n;
                end
            end
        case 'ols'
            for p = 1:nb_pattern
                grad = res*amp(p,:)';
                if norm(grad) > 0
                    pat = dico(:,p) + rho*grad;
                    pat = pat/norm(pat);
                    if nargin >5 && strcmp(flow, 'sequential')
                        res = res + (dico(:,p)-pat)*amp(p,:); %#ok<*NASGU>
                    end
                    dico(:,p) = pat;
                end
            end
        case 'opt'
            for p = 1:nb_pattern
                index = find(amp(p,:)~=0);
                vec = amp(p,index);
                grad = res(:,index)*vec';
                if norm(grad) > 0
                    pat = (vec*vec')*dico(:,p) + rho*grad;
                    pat = pat/norm(pat);
                    if nargin >5 && strcmp(flow, 'sequential')
                        res(:,index) = res(:,index) + (dico(:,p)-pat)*vec;
                    end
                    dico(:,p) = pat;
                end
            end
        case 'LGD'
            for p = 1:nb_pattern
                index = find(amp(p,:)~=0);
                vec = amp(p,index);
                grad = res(:,index)*vec';
                if norm(grad) > 0
                    pat = (vec*vec')*dico(:,p) + 2*grad;
                    pat = pat/norm(pat);
                    if nargin >5 && strcmp(flow, 'sequential')
                        res(:,index) = res(:,index) + (dico(:,p)-pat)*vec;
                    end
                    dico(:,p) = pat;
                end
            end
        case 'KSVD'
            for p = 1:nb_pattern
                index = find(amp(p,:)~=0);
                if ~isempty(index)
                    patch = res(:,index)+dico(:,p)*amp(p,index);
                    [U,~,V] = svd(patch);
                    if U(:,1)'*dico(:,p) > 0
                        dico(:,p) = U(:,1);
                    else
                        dico(:,p) = -U(:,1);
                    end
                    dico(:,p) = dico(:,p)/norm(dico(:,p));
                    amp(p,index) = dico(:,p)'*patch;
                    if nargin >5 && strcmp(flow, 'sequential')
                        res(:,index) = patch-dico(:,p)*amp(p,index);
                    end
                end
            end
    end
end

