%% User input for joccam: joint or single, static or time-lapse inversion

%% Grid parameters
d.dim=3;                     % dimension should be 2 or 3
d.grid.dx=2;              % pixel sizes for inverse solution
d.grid.dy=2;
d.grid.dz=2;
d.grid.mx=33;                % number of pixels in effective grid
d.grid.my=31;                % set to 1 for 2D; set to 5 or 7 for 2.5D
d.grid.mz=24;
%% Regularization
d.inv.regularization=1;      % 1 is covariance and 2 smoothing regularization operator
d.inv.IRLSd=0;               % If IRLSd=1, perform iteratively reweighted least-squares for data
d.inv.IRLSm=0;               % If IRLSm=1, perform iteratively reweighted least-squares for model
d.inv.IRLSp=1;               % exponent p in perturbed lp norm of Ekblom
d.inv.IRLSe=0.1;             % epsilon in the perturbed lp norm of Ekblom
% Model covariance description if regularization==1
d.inv.Ix=5;                % Integral scale in x-direction (assuming stationarity)
d.inv.Iy=5;              	 % Integral scale in y-direction
d.inv.Iz=5;              	 % Integral scale in z-direction
d.inv.cutoff=0.01;        	 % cut-off criteriion for the stochastic regularization
                    	     % where 0.01 is recommended (1% of maxvalue)
% Anisotropy in smoothing operator with regards to z-direction;
d.inv.anix=1;                 
d.inv.aniy=1;
%% Weights and step sizes used in the inversion
d.inv.slambda=3;             % Initial trade-off parameter
d.inv.dlambda=0.2;           % Delta lambda: step size of trade-off parameter
d.inv.mlambda=0.1;           % Minimum lambda
d.inv.jlambda=10000;         % Cross gradients weight (order 10000)
d.inv.cg_ind=[1 2; 1 3];     % 
d.inv.sRMS=1.2;              % RMS at which to start searching for a good lambda (instead of decreasing it monotoneously)
d.inv.fRMS=1;                % Final target RMS
d.inv.maxit=10;              % Number of iterations
%% General parameters
d.logfile='tomo3d.log';     % Logfile
d.inv.synthetic=0;           % If synthetic==1, data will be contaminate with noise
d.inv.psfi=10;                % Iteration where to save information to calculate the psf using jpsf.m

%% Data sets
i=1; % number of data set
d.ds(i).name        = 'SEISMIC'; % arbitrary name for the method / data set
d.ds(i).method      = 'tt'; % use 'tt' for traveltime inversion or 'er' for electrical resistivity 
d.ds(i).weight      = 1;    % weight of this data set
d.ds(i).folder      = 'seismic';  % folder for the forward calculation
d.ds(i).abs_error   = 0.02;     % absolute error (in units used)
d.ds(i).rel_error   = 0.005;     % relative error
d.ds(i).init        = 0.2;      % initial model value (slowness or resistivity)
d.ds(i).var         = 0.0025;     % model variance
d.ds(i).mtrans      = '';   % invert for logarithm  of the model parameters
d.ds(i).dtrans      = '';   % use logarithm  of the input data
% optional
%d.ds(i).init_1d     = '';    % file with 1d starting model (length must be d.grid.mz)
d.ds(i).dataf       = 'tt_all.dat'; % data file
%d.ds(i).tl_dataf    = 'seismic';    % beginning of the time-lapse file name. %i.dat will be added at the end (e.g. seismic1.dat)

