% initpaths.m
%
% Initialize paths needed for MIDlnp code

addpath code_iSTAC   % directory for iSTAC (used to initialize filter estimates)
addpath code_MIDlnp  % directory for MID code
addpath nlfuns;      % directory for nonlinearities

if ~exist('simdata','dir')
    mkdir simdata;
end
