function [L,newS,newF,Su,f] = fitvK(u,meanU,dt,component,varargin)
%  [L,newS,newF,Su,f] = fitvK(u,meanU,dt,component,varargin) computes the
%  integral length scale using a least-square fit of the von karman
%  spectrum to the estimated power spectral densities (PSDs).
%% Input
%   u: fluctuating component [1xN] float
%   meanU: mean wind speed [1x1] float
%   dt: time step in seconds
%   component : 'u','v' or 'w'
% Optional parameter: only used for the fitting procedure (see lsqcurvefit)
%% Output
% L = [1 x 1] double : Integral length scale
% newS: [1xNf] double: fitted von karman spectrum
% newF: [1xNf] double:frequency associated with newS
% Su: [1xM] double: Estimated PSD 
% f: [1xM] double:frequency associated with newS
% 
% Author: E. Cheynet - UiB - last modified: 27-03-2022

%% Inputparser
p = inputParser();
p.CaseSensitive = false;
p.addOptional('dispIter','off');
p.addOptional('tolX',1e-5);
p.addOptional('tolFun',1e-5);
p.addOptional('MaxFunEvals',600);
p.addOptional('Nfreq',60);
p.parse(varargin{:});
% shorthen the variables name
tolX = p.Results.tolX ;
tolFun = p.Results.tolFun ;
dispIter = p.Results.dispIter ;
MaxFunEvals = p.Results.MaxFunEvals ;
Nfreq = p.Results.Nfreq ;

assert(license('test','Curve_Fitting_Toolbox')==1,'The cohFit function requires Matlab''s Curve Fitting Toolbox.')
options=optimset('TolX',tolX,'TolFun',tolFun,'Display',dispIter,'MaxFunEvals',MaxFunEvals);

%% Compute the PSD with the welch method and 3 overlapping segments
N = numel(u);
fs = 1/dt;

nfft = round(600*fs);
% [Su0,f0] = pburg(detrend(u),30,nfft,fs);
[Su0,f0] = pwelch(detrend(u),nfft,round(nfft/2),N,fs);
Su0 = Su0./trapz(f0,Su0);
% Remove frequencies above 1 Hz (not useful for the spectral fit)
Su0 = Su0(f0>0 & f0<=1);
f0 = f0(f0>0 & f0<=1);

%% bin average the psd estimates
f3 = logspace(log10(f0(1)),log10(f0(end)),Nfreq);
[Su,f] = binAveraging(Su0,f0,'newX',f3);
%% Fit the von karman model
modelFun = @vK;
guess = 50;

L = lsqcurvefit(@(L,f) modelFun(L,f),guess,f,Su,1,500,options);
newF = logspace(log10(f(1)),log10(f(end)),Nfreq);
newS = vK(L,newF);


    function [S] = vK(L,f)
        % Normalized von karman spectrum
       
        
        fr = L.*f./meanU;
        
        
        if strcmpi(component,'u')
            S = 4*fr./(1+71.*fr.^2).^(5/6);
        elseif strcmpi(component,'v')
            S=  4*fr.*(1+755*fr.^2)./(1+283.*fr.^2).^(11/6);
        elseif strcmpi(component,'w')
            S=  4*fr.*(1+755*fr.^2)./(1+283.*fr.^2).^(11/6);
        else
            error('component unknown')
        end
        
        S = S./f;
        varU = trapz(f,S);
        S = S./varU;
        
        
    end








    function [y,x,varargout] = binAveraging(y0,x0,varargin)
        % function [y,x,varargout] = binAveraging(Y,x0,varargin) computes
        % non-overlapping bin-averaged quantities from vectors.
        %
        % Input:
        % x0: vector [1 x N] of absiccsa values, e.g. a time vector
        % y0: vector [1 x N] of ordinate values in each x0, e.g. a time series
        % varargin
        %
        %
        % options for Input:
        % newX: vector [1 x M]: new abcissa
        % Nbin: scalar [1 x 1]: number of bins (overwritten by newX, if specified)
        % binWidth: scalar [1 x 1]: width of bins (overwritten by newX, if specified)
        % BMAX: scalar: Max value of the bin (x0<=BMAX)
        % BMIN: scalar: Min value of the bin (x0>=BMIN)
        % averaging: string: 'mean" or 'median'
        % dispersion: string: 'std', 'decile_10_90' or 'decile_25_75'
        %
        % Output:
        % y: vector [1 x M]: binned values of y0
        % x: vector [1 x M]: new absiccsa
        % varargout: if specified, provides a [ 1x M] vector for the dispersion
        %
        % Author: E. Cheynet  - UiB - Norway
        % LAst modified 06/12/2019
        %
        %  see also histcounts accumarray
        
        %%  Inputparser
        
        p = inputParser();
        p.CaseSensitive = false;
        p.addOptional('newX',[]);
        p.addOptional('binMethod','auto');
        p.addOptional('BMIN',nanmin(x0));
        p.addOptional('BMAX',nanmax(x0));
        p.addOptional('binWidth',[]);
        p.addOptional('Nbin',[]);
        p.addOptional('averaging','mean');
        p.addOptional('dispersion','std');
        p.parse(varargin{:});
        % shorthen the variables name
        BMIN = p.Results.BMIN ;
        BMAX = p.Results.BMAX ;
        binMethod = p.Results.binMethod ;
        binWidth = p.Results.binWidth ;
        Nbin = p.Results.Nbin;
        newX = p.Results.newX ;
        averaging = p.Results.averaging ;
        dispersion = p.Results.dispersion ;
        y0 = y0(:);
        x0 = x0(:);
        
        %%
        if ~isempty(newX)
            [~,~,bin] = histcounts(x0,newX);
        elseif ~isempty(Nbin)
            [~,~,bin] = histcounts(x0,Nbin);
        elseif isempty(newX) && isempty(binWidth)
            [~,~,bin] = histcounts(x0,'BinMethod',binMethod,'BinLimits',[BMIN,BMAX]);
        elseif isempty(newX) && ~isempty(binWidth)
            [~,~,bin] = histcounts(x0,'BinMethod',binMethod,'BinLimits',[BMIN,BMAX],'binWidth',binWidth);
        else
            error('unspecfied options')
        end
        % bin(bin==0)=1;
        
        
        %% Mean values
        if strcmpi(averaging,'mean')
            y = accumarray(bin(bin>0), y0(bin>0), [], @nanmean, nan);
            x = accumarray(bin(bin>0), x0(bin>0), [], @nanmean, nan);
        elseif strcmpi(averaging,'median')
            y = accumarray(bin(bin>0), y0(bin>0), [], @nanmedian, nan);
            x = accumarray(bin(bin>0), x0(bin>0), [], @nanmedian, nan);
        else
            error(' ''averaging'' should be ''mean'' or ''median'' ')
        end
        
        
        %% Dispersion
        if strcmpi(dispersion,'std'),
            stdY = accumarray(bin(bin>0), y0(bin>0), [], @nanstd, nan);
        elseif strcmpi(dispersion,'decile_10_90'),
            stdY(1,:) = accumarray(bin(bin>0), y0(bin>0), [], @(x) quantile(x,0.1), nan);
            stdY(2,:) = accumarray(bin(bin>0), y0(bin>0), [], @(x) quantile(x,0.9), nan);
        elseif strcmpi(dispersion,'decile_25_75'),
            stdY(1,:) = accumarray(bin(bin>0), y0(bin>0), [], @(x) quantile(x,0.25), nan);
            stdY(2,:) = accumarray(bin(bin>0), y0(bin>0), [], @(x) quantile(x,0.75), nan);
        else
            error(' ''std'' should be ''decile_10_90'' or ''decile_25_75'' ')
        end
        
        stdY = stdY';
        
        
        y(isnan(x))=[];
        x(isnan(x))=[];
        
        
        if nargout ==3,    varargout = {stdY};end
        
        
        
    end



end

