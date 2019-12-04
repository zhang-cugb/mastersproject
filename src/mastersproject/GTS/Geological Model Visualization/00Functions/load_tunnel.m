
function S = load_S1S3(BH_file, BH_dir)
% Export format: X, Y, Z
%
% Input:
%   [BH_file]  - optional: geology data filename
%   [BH_dir]   - optional: directory name where datafile is located
%
% Output:
%   geology data from OPTV
%       -   .x_m    : x-coordinate 
%       -   .y_m    : y-coordinate 
%       -   .z_m    : y-coordinate
%
% Nathan Dutler, 10.11.2017
% =========================================================================

%% get elevation datafilestring =================================================
if nargin == 1
    disp('Please provide data directory')
    return
elseif nargin < 1
	[BH_file,BH_dir] = uigetfile({ ...
        '*.txt', 'elevation(*.txt)'});
end

%% read elevation data ----------------------------------------------------

fidm     = fopen([BH_dir BH_file], 'r');
ft       = '%f %f  %f %f %f %s \r \n';
%            line X  Y  Z  LINE MARK  X   Y   Z
file_out = textscan(fidm, ft, 'Headerlines', 2);

S.x_m    = file_out{1};    % x position [m]
S.y_m    = file_out{2};    % y position [m]
S.z_m    = file_out{3};    % z position [m]
S.length = file_out{4};    % azimuth [m]
S.extra_l= file_out{5};    % dip [m]
S.text   = file_out{6};


fclose(fidm);