function tunnel = excav_points(BH_file, BH_dir, header)
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
ft       = '%f  %f   %f \r \n';
%            line X  Y  Z  LINE MARK  X   Y   Z
file_out = textscan(fidm, ft, 'Headerlines', header);

tunnel.x_m    = file_out{1};    % x position [m]
tunnel.y_m    = file_out{2};    % y position [m]
tunnel.z_m    = file_out{3};    % z position [m]



fclose(fidm);