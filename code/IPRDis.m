function [dic, mus, srs] = IPRDis(dic,mu,Y,X,cat,par)
%% Unit test
if nargin==1        %if input argument is a single number, run appropriate unit test
    feval(strcat('UnitTest',num2str(dic)));
    return;
end

%% Parameters and defaults
if ~exist('par','var') || isempty(par), par = struct; end   %parameters

def.nIte  = 20;     %number of iterations
def.eps   = 1e-6;   %tolerance level
def.alpha = .7;     %mix weight
def.vis   = false;  %visualize results (only for 2 or 3 dimensional datasets)
def.ver   = false;  %verose output

par = setdefaultoptions(def,par);   %set default parameters

nDim = size(dic,1);  %dimension of atoms of the dictionary

if isnumeric(cat), cat = num2str(cat); end
if ~isempty(cat)
    atoCat = ComputeAtomsCategories(X,cat);
    MU  = @(dic) MUCat(dic,atoCat); %mutual coherence function
else
    atoCat = [];
    MU  = @(dic) max(max(abs((dic'*dic)-diag(diag(dic'*dic))))); %mutual coherence function
end
SNR = @(dic) snr(Y,dic*omp(dic'*Y,dic'*dic,par.Tdata));

%% Main algorithm
dic = normc(dic);               %normalise columns

mus = zeros(par.nIte,1);		%coherence at each iteration
srs = zeros(par.nIte,1);		%signal to noise ratio at each iteration
iIter = 1;

if par.vis
        fprintf(1,'\n--------------------------------\n');
end
while iIter<=par.nIte && MU(dic)>mu
    % calculate snr and coherence
    mus(iIter) = MU(dic);
    srs(iIter) = SNR(normc(dic));
    
    if par.vis
        fprintf(1,'Iteration number %u - mu: %f, snr: %f\n', iIter,mus(iIter),srs(iIter));
    end
    
    % calculate gram matrix
    G = dic'*dic;
    
    %project onto the structural constraint set and cumpute convex combination
    H = ProjectStruct(G,mu,atoCat);
    alpha = par.alpha;
    G = alpha*H + (1-alpha)*G;  %convex combination (soft projection)
    
    % project into spectral constraint set and compute dictionary
    dic  = ProjectSpec(G,nDim);
    
    if par.vis
%         pause;
        scale = 3;
        if ~isempty(cat)
            for iCat=1:length(unique(cat))
                set(par.ploHan,'UData',scale*dic(1,atoCat==iCat));
                set(par.ploHan,'VData',scale*dic(2,atoCat==iCat));
            end
        else
            set(par.ploHan,'UData',scale*dic(1,:));
            set(par.ploHan,'VData',scale*dic(2,:));
        end
        drawnow
%         set(par.texHan,'String','Rotation');
    end
    
    
    % rotate dictionary
    C = Y*(dic*X)';
    [U, ~, V] = svd(C);
    W = V*U';
    dic = W*dic;
    
    if par.vis
%         pause(0.05)
        if ~isempty(cat)
            for iCat=1:length(unique(cat))
                set(par.ploHan,'UData',scale*dic(1,atoCat==iCat));
                set(par.ploHan,'VData',scale*dic(2,atoCat==iCat));
            end
        else
            set(par.ploHan,'UData',scale*dic(1,:));
            set(par.ploHan,'VData',scale*dic(2,:));
        end
        drawnow
    %set(par.texHan,'String','Decorrelation');
    end
    
    iIter = iIter+1;
end
if iIter<par.nIte
    mus(iIter:end) = mus(iIter);
    srs(iIter:end) = srs(iIter);
end
end

function dic  = ProjectSpec(G,nDim)
[Q, Lambda] = eig(G);           %evd decomposition of Gram matrix
[thrLam, idx] = ThreshLambda(Lambda,nDim);
Q = Q(:,idx);   %permute columns of Q to correspond to sorted eigenvalues
dic = sqrt(thrLam(1:nDim,:))*Q';
assert(rank(dic)<=nDim);
end

function [thrLam,idx] = ThreshLambda(Lambda,n)
[v, idx] = sort(diag(Lambda),'descend');    %sort eigenvalues from larger to smaller
v(v<0) = 0;                               %threshold negative eigenvalues
v(n+1:end) = 0;                    %threshold last nAto-nDim eigenvalues
thrLam = diag(v);                   %thresholded lambda
end

function H = ProjectStruct(G,mu,atoCat)
H = G;				%initialise matrix
nAto = size(G,1);
if ~isempty(atoCat)
    for iAtom=1:nAto
        for jAtom=1:nAto
            if iAtom==jAtom
                H(iAtom,jAtom) = 1;
            elseif abs(G(iAtom,jAtom)) > mu && atoCat(iAtom)~=atoCat(jAtom)
                H(iAtom,jAtom) = sign(G(iAtom,jAtom))*mu;
            end
        end
    end
else
    ind1 = find(abs(G(:))<=mu);		%find elements smaller than mu
    ind2 = find(abs(G(:))>mu);		%find elements bigger than mu
    H(ind1) = G(ind1);				%copy elements belonging to ind1
    H(ind2) = mu*sign(G(ind2));		%threshold elements belonging to ind2
    H(1:nAto+1:end) = 1;				%set diagonal to one
end
end


function UnitTest1 %#ok<DEFNU>
clear, clc, close all
stream = RandStream.getGlobalStream;    %set random number generator to reproduce results
stream.reset;
nData = 30000;						%number of data
theta = [pi/16, pi/8, pi/3, 4*pi/6];			%angles
nAngles = length(theta);			%number of angles
Q	  = [cos(theta); sin(theta)];	%rotation matrix
Y	  = Q*randmog(nAngles,nData);	%training data

figure, hold on
scale = 3;							%scale factor for plots
scatter(Y(1,:), Y(2,:),5,'k');		%scatter training data
theta1 = [pi/16,pi/8,4*pi/6];
phi = [cos(theta1); sin(theta1)];
O = zeros(size(phi));					%origin
hQ = quiver(O(1,:),O(2,:),scale*phi(1,:),scale*phi(2,:),'LineWidth',4,'Color','r');		%plot atoms
par.Tdata = 1;
X = omp(phi,Y,[],par.Tdata);
% hT = text(-2,5.5,'De-correlation','FontSize',20,'Color','k');
par.ploHan = hQ;
% par.texHan = hT;
xlim([-6,6]), ylim([-6,6]);
par.nIte = 5;
par.alpha = 1;
pause
IPRDis(phi,cos(30),Y,X,[],par);
end

function UnitTest2 %#ok<DEFNU>
clear, clc, close all
stream = RandStream.getGlobalStream;    %set random number generator to reproduce results
stream.reset;
nData = 10000;						%number of data
theta = [pi/16, pi/8, pi/3, 4*pi/6];			%angles
nAngles = length(theta);			%number of angles
Q	  = [cos(theta); sin(theta)];	%rotation matrix
Y	  = Q*randmog(nAngles,nData);	%training data
A = Y(1,:)./Y(2,:);
B = Q(1,:)./Q(2,:);
for i=1:nData
    for j=1:nAngles
        d(i,j) = abs(A(i)-B(j));
    end
    [~, traCat(i)] = min(d(i,:));
end
traCat(traCat==2) = 1;
traCat(traCat==3) = 2;
traCat(traCat==4) = 2;
figure
subplot(1,2,1), hold on
gscatter(Y(1,:),Y(2,:),traCat,'kr','o',4);
xlim([-4,4]), ylim([-4,4]);

theta1 = [pi/16,pi/8,4*pi/6];
phi = [cos(theta1); sin(theta1)];
O = zeros(size(phi));%origin
uniCat = unique(traCat);
nCat = length(uniCat);
X = omp(phi,Y,[],1);


C = zeros(size(phi,2),nCat);
for iCat=1:nCat
    C(:,iCat) = sum(abs(X(:,traCat==uniCat(iCat))'));
end
[~, atoCat] = max(C,[],2);  %cluster atoms
catCol = {'k','r'};
hQ = nan(length(unique(traCat)),1);
scale = 3;
for iCat=1:length(unique(traCat))
    idx = (atoCat==uniCat(iCat));
    hQ(iCat) = quiver(O(1,idx),O(2,idx),scale*phi(1,idx),scale*phi(2,idx),'LineWidth',2,'Color',catCol{iCat});		%plot atoms
end
hT = text(-2,3.5,'De-correlation','FontSize',12,'Color','k');

par.ploHan = hQ;
par.texHan = hT;
Ypro = phi*X;
subplot(1,2,2)
gscatter(Ypro(1,:),Ypro(2,:),traCat,'kr','o',4);



pause
phi = IPRDis(phi,cos(30),Y,X,traCat,par);
pause

X = omp(normcols(phi),Y,[],1);
Ypro = phi*X;
gscatter(Ypro(1,:),Ypro(2,:),traCat,'kr','o',4);
end

function X = randmog(m, n)
% RANDMOG - Generate mixture of Gaussians
s = [0.2 2];
% Choose which Gaussian
G1 = (rand(m, n) < 0.9);
% Make them
X = (G1.*s(1) + (1-G1).*s(2)) .* randn(m,n);
end

function [mu, mask] = MUCat(dic,cat)
nAto = size(dic,2);
mask = nan(nAto,nAto);
for iAto=1:nAto
    for jAto=1:nAto
        if cat(iAto)==cat(jAto)
            mask(iAto,jAto) = 0;
        else
            mask(iAto,jAto) = 1;
        end
    end
end
G = (dic'*dic).*mask;
mu = max(max(abs(G)));
end