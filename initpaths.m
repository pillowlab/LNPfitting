% initpaths.m
%
% Initialize paths needed for MIDlnp code

addpath tools_iSTAC   % directory for iSTAC (used to initialize filter estimates)
addpath tools_LNPfitting  % directory for MID code
addpath nlfuns;      % directory for nonlinearities

if ~exist('simdatadir','dir')
    mkdir simdatadir;
end
