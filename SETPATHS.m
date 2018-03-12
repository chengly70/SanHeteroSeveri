% SETPATHS.m - for SA code repository
%
% This simple script sets the path to include relevant directories for the
% You must 'cd' into this directory in order
% to evaluate it.
%


basedir = pwd;  % The directory where this script lives

% Add a bunch sub-directories (with absoluate path names)
addpath([basedir]);
addpath([basedir '/thry/']);
addpath([basedir '/xpp_files/']);

