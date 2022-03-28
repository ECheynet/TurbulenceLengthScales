function [L,data] = Lyz(u,pos)
% [Luv] = Lyz(u,v,dt,U) computes the crosswind turbulence length scale
% using an exponential fit to the correlation coefficients expressed as a
% function of the crosswind separation distance
%
% Input
%   u:  [Mxn] double: velocity turbulent component, can be u, v or w
%   pos: [Mx1] position of the sensors in a linear array (either pure
%   lateral or vertical separation)
% Ouput
%  L: [1x1] double:  crosswind turbulence length scale
%  data: structure variable useful to plot the correlation coefficients as
%  a function of the crosswind distance
% Author: E Cheynet - uiB - last modified 28/03/2022

%% Compute the correlation coefficient


[M,N]=size(u);
if N<M
    u=u';
    [M,N]=size(u);
elseif M ==1 && N >1
    error(' Correlation matrix needs more than 1 sensor');
end

% detrend the data
for ii=1:M,
    u(ii,:) = detrend(u(ii,:));
end


[R,p]=corrcoef(u');
R = R(:);

d = (pos(:)-pos(:)');
d = d(:);

R = R(d>=0);
d = d(d>=0);


%% Get the median of the distances when these are within +/- 10cm
[du,~,ind]=unique(round(d*10)/10);

clear newD newR err
for ii=1:numel(du),
    newR(ii) = nanmedian(R(ind==ii));
    newD(ii) = nanmedian(d(ind==ii));
    err(ii) = nanstd(R(ind==ii));
end

%% non linear fit
expoFit =  @(L,x) exp(-x./L(1)) + L(2); % L(1) is the  length scale and L(2) is the random error
guess = [50,0];


options=optimset('TolX',1e-4,'TolFun',1e-4,'Display','off');
coeff = lsqcurvefit(@(L,f) expoFit(L,f),guess,newD,newR,[1,0],[1000,0.9],options);
% coeff = nlinfit(newD,newR, expoFit, guess);
L = coeff(1);

data.R = newR;
data.d = newD;
data.err = err;
data.fun = expoFit;
data.coeff = coeff;
end

