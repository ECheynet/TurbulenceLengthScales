function L = Lx(u,meanU,dt,varargin)
% L = Lx(u,meanU,dt,varargin) calculates the streamwise integral length
% scale using the autocovariance function. It either compute the area under
% the curve from a zero lag to the frist zero crossing (method 1) or use an
% exponential fit (method 2)
%% Input
%   u: fluctuating component [1 x N] double
%   dt: time step in seconds [1x1] double
%   meanU : mean wind speed [1x1] double
% Optional parameter
%   -   method:  2 (exponential fit) or 1 (direct integration) 
%   -   tmax: [1x1] double
% 
% Output
% L = [1 x 1] double : Integral length scale
% 
% 
% 
% Author: E. Cheynet - UiB - last modified: 27-03-2022

%% Inputparseer
p = inputParser();
p.CaseSensitive = false;
p.addOptional('method',1);
p.addOptional('tmax',[]);
p.parse(varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%
method = p.Results.method ; % Turbulence factor, taken as 1.0 by default.
tmax = p.Results.tmax ;

%% AUTOCOVARIANCE

u = detrend(u); % remove the linear trend
N = numel(u);
if isempty(tmax), tmax = dt.*N;end



[dummy, lags] = xcov(u,u,N,'coef'); % is 2-sided and normalized

% Get the 1-sided autocovariance function
R = dummy(N+1:end);
tLag =lags(N+1:end)*dt;
        

% get the indice of the first zero-crossing
ind = zerocross(R);

if numel(ind)>=1
    % get the first zero-crossing
    ind = ind(1);
else 
    warning('No zero crossing found')
    L=NaN;
    return
end


%% non linear exponential fit if method==2
if method==2
    % Exponential fit
    expoFit =  @(coeff,x) exp(-coeff.*x);
    guess = 0.01; % first guess
    % fitting process
    try
        
        [~,indThres] = min(abs(tLag-tmax));
        indThres = min(indThres+1,ind);
        coefEsts = nlinfit(tLag(1:indThres),R(1:indThres), expoFit, guess);
        %             plot(tLag,autocov,tLag(1:ind),expoFit(coefEsts,tLag(1:ind)));xlim([0,100])
    catch exception % if statistic toolbox is not available, try with optimization toolbox'
        try
            option= optimset('Display','off');
            [~,indThres] = min(abs(tLag-tmax));
            indThres = min(indThres+1,ind);
            coefEsts = lsqcurvefit(@(L,t) expoFit(L,t),guess,tLag(1:indThres),R(1:indThres),0,500,option);
        catch exception
            disp exception
            warning('No toolbox available. Method is set to 1')
            method = 1;
        end
    end
end

%% Get the integral length scale
% Taylor hypothesis of frozen turbulence is applied
if method==1
    T = trapz(tLag(1:ind),R(1:ind));
    L = T*meanU;
elseif method==2
    % get integral time scale
    T = trapz(tLag(1:ind),expoFit(coefEsts,tLag(1:ind)));
    L = T*meanU;
else
    error(['method: "',method,...
        '" is unknown. Please, choose between:',...
        ' "expoDecay" or "DirectInt"'])
end

%% Nested functions
    function z=zerocross(v)
        z=find(diff(v>0)~=0)+1;
    end

%% Previous test code

% y = fft(detrend(u,'constant'));
% h0 = ifft(y.*conj(y));
% autocov = h0(1:round(numel(u)/2));
% autocov = autocov./nanmax(autocov);
% tLag = [0:numel(autocov)-1].*dt;
end

