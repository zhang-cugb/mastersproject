%% User input for joccam: joint or single, static or time-lapse inversion

%% Grid parameters
d.dim=2;                     % dimension should be 2 or 3
d.grid.dx=1;              % pixel sizes for inverse solution
d.grid.dy=1;
d.grid.dz=1;
d.grid.mx=64;                % number of pixels in effective grid
d.grid.my=1;                % set to 1 for 2D; set to 5 or 7 for 2.5D
d.grid.mz=62;
%% Regularization
d.inv.regularization=1;      % 1 is covariance and 2 smoothing regularization operator
d.inv.IRLSd=0;               % If IRLSd=1, perform iteratively reweighted least-squares for data
d.inv.IRLSm=0;               % If IRLSm=1, perform iteratively reweighted least-squares for model
d.inv.IRLSp=1;               % exponent p in perturbed lp norm of Ekblom
d.inv.IRLSe=0.1;             % epsilon in the perturbed lp norm of Ekblom
% Model covariance description if regularization==1
d.inv.Ix=10;                % Integral scale in x-direction (assuming stationarity)
d.inv.Iy=10;              	 % Integral scale in y-direction
d.inv.Iz=10;              	 % Integral scale in z-direction
d.inv.cutoff=0.01;        	 % cut-off criteriion for the stochastic regularization
                    	     % where 0.01 is recommended (1% of maxvalue)
% Anisotropy in smoothing operator with regards to z-direction;
d.inv.anix=1;                 
d.inv.aniy=1;
%% Weights and step sizes used in the inversion
d.inv.slambda=2.6;             % Initial trade-off parameter
d.inv.dlambda=0.2;           % Delta lambda: step size of trade-off parameter
d.inv.mlambda=1.0;           % Minimum lambda
d.inv.jlambda=0;          % Cross gradients weight (order 10000)
d.inv.sRMS=0.5;              % RMS at which to start searching for a good lambda (instead of decreasing it monotoneously)
d.inv.fRMS=0.4;                % Final target RMS
d.inv.maxit=10;              % Number of iterations
%% General parameters
d.logfile='9_aniso.log';       % Logfile
d.inv.synthetic=0;           % If synthetic==1, data will be contaminate with noise
d.inv.psfi=10;               % Iteration where to save information to calculate the psf using jpsf.m

%% Data sets
% Include as many data sections as needed. Joint inversion currently works
% with max. 3 data sets.
i=1; % number of data set P-wave velocity
d.ds(i).name        = 'aniso_v'; % arbitrary name for the method / data set
d.ds(i).method      = 'tt_aniso_v'; % use 'tt' for traveltime inversion or 'er' for electrical resistivity 
d.ds(i).weight      = 1;    % weight of this data set
d.ds(i).folder      = 'seis_aniso';  % folder for the forward calculation
d.ds(i).abs_error   = 0.0001;     % absolute error (in units used)
d.ds(i).rel_error   = 0.0;     % relative error
d.ds(i).init        = 1/4689;      % initial model value (slowness or resistivity)
d.ds(i).var         = (1/4689)^2;     % model variance
d.ds(i).mtrans      = '';   % invert for logarithm  of the model parameters
d.ds(i).dtrans      = '';   % use logarithm  of the input data
%d.ds(i).eps         = 0.072; % only needed, if epsilon is not inverted for
%d.ds(i).delta       = 0.081; % only needed, if delta is not inverted for
%d.ds(i).eps         = 0.065; % only needed, if epsilon is not inverted for
%d.ds(i).delta       = 0.112; % only needed, if delta is not inverted for
d.ds(i).theta       = 133/180*pi; % only needed, if theta is not inverted for
d.ds(i).sym_axis    = [.5999 0.225 -.7678];
i=i+1; % number of data set anisotropic epsilon
d.ds(i).name        = 'aniso_e'; % arbitrary name for the method / data set
d.ds(i).method      = 'tt_aniso_e'; % use 'tt' for traveltime inversion or 'er' for electrical resistivity 
d.ds(i).weight      = 0.04;    % weight of this data set
d.ds(i).folder      = 'seis_aniso';  % folder for the forward calculation
d.ds(i).init        = 0.073;      % initial model value (slowness or resistivity)
d.ds(i).var         = 0.073*0.1;     % model variance
d.ds(i).mtrans      = '';   % invert for logarithm  of the model parameters
d.ds(i).dtrans      = '';   % use logarithm  of the input data
i=i+1; % number of data set anisotropic delta
d.ds(i).name        = 'aniso_d'; % arbitrary name for the method / data set
d.ds(i).method      = 'tt_aniso_d'; % use 'tt' for traveltime inversion or 'er' for electrical resistivity 
d.ds(i).weight      = 0.04;    % weight of this data set
d.ds(i).folder      = 'seis_aniso';  % folder for the forward calculation
d.ds(i).init        = 0.106;      % initial model value (slowness or resistivity)
d.ds(i).var         = 0.106*0.1;     % model variance
d.ds(i).mtrans      = '';   % invert for logarithm  of the model parameters
d.ds(i).dtrans      = '';   % use logarithm  of the input data
% i=i+1; % number of data set anisotropic theta
% d.ds(i).name        = 'aniso_t'; % arbitrary name for the method / data set
% d.ds(i).method      = 'tt_aniso_t'; % use 'tt' for traveltime inversion or 'er' for electrical resistivity 
% d.ds(i).weight      = 0.1;    % weight of this data set
% d.ds(i).folder      = 'seis_aniso';  % folder for the forward calculation
% d.ds(i).init        = 185/180*pi;      % initial model value (slowness or resistivity)
% d.ds(i).var         = pi*0.1;     % model variance
% d.ds(i).mtrans      = '';   % invert for logarithm  of the model parameters
% d.ds(i).dtrans      = '';   % use logarithm  of the input data
i=i+1; % number of data set station correction
d.ds(i).name        = 'aniso_st'; % arbitrary name for the method / data set
d.ds(i).method      = 'tt_aniso_st'; %station correction
d.ds(i).weight      = 4;    % for station correction: damping parameter (higher weight, more damping)
d.ds(i).folder      = 'seis_aniso';  % folder for the forward calculation
d.ds(i).init        = 0;      % initial model value (slowness or resistivity)
d.ds(i).mtrans      = '';   % invert for logarithm  of the model parameters
d.ds(i).dtrans      = '';   % use logarithm  of the input data

%% Time lapse option - if no time-lapse inversion desired, remove this part
% of set d.tl.steps=0;
d.tl.steps = 0;                 % time lapse steps (set to 0 if only static inversion)
d.tl.slambda = d.inv.slambda;
d.tl.maxit = 5;             % maximum iteration in each time lapse step
d.tl.mreffile = '';               % reference model (file with reference models (as written by occ_write_models with option ref)). If mref='', mref is inverted for.
%% jcinv parameters: parameters for inversion of the best parameter within
%  each zone. Just append this part to the input file you used for the
%  joccam inversion.
d.jcinv.calpha  = 0.9;          % update stepsize (1 for full update, 0 for none)
d.jcinv.dalpha  = 0;          % difference in alpha to test
d.jcinv.clambda = 0;          % damping (~100)
d.jcinv.cfile   = 'classes.txt';% text file containing the classes of each pixel
