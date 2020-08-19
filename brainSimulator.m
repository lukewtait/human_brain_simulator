function out = brain_simulator(pars)
% Brain simulator, written by Luke Tait
% 
% Outputs: A structure 'out', which contains the following fields: 
% out.time: time axis
% out.source: source space simulation
% out.sensor: sensor space simulation
% out.pars: Full parameters structure (see below). 
% Additionally, if pars.options.log_output is true (default), a text file
% will be created logging all output to the console, details of the
% simulations, and relevant citations. 
% 
% Inputs: a structure 'pars'
% If pars is empty, then all default settings are run. These are as
% outlined in the draft of my manuscript 'Network Substrates of Cognitive
% Impairment in Alzheimer's Disease', using the optimized connectome fit
% for the healthy older adults with natural frequencies drawn uniformly in
% the 6-13 Hz combined theta-alpha band (see manuscript). 
% 
% Alternatively, deviations from the default can be chosen. These are: 
% pars.atlas: 
% -> pars.atlas.atlas
%    Defines the atlas being used. Implemented are: 
%    -> 'brainnetome40' (default). The downsampled brainnetome atlas using 40
%        cortical regions as described in the manuscript and at http://atlas.brainnetome.org/
%    -> 'desikankilliany'. https://doi.org/10.1016/j.neuroimage.2006.01.021
%    -> 'destrieux'. https://doi.org/10.1016/j.neuroimage.2010.06.010
%    -> 'AAL'. https://doi.org/10.1006/nimg.2001.0978
%    -> 'custom'. Requires a custom leadfield for custom regions. 
% -> pars.atlas.regions
%    Defines the regions of the atlas to be used. This should be either
%    empty (pars.atlas.region = []) to use all regions in the atlas
%    (default) or a column vector of region indices (pars.atlas.region =
%    [1;3;5] to use the 1st, 3rd, and 5th regions). If pars.atlas.atlas is
%    custom, regions are simply labeled numerically. 
% 
% pars.connectome:
%    Defines the connectome. This can either be empty, which then uses a
%    default connectome (default; only implemented for desikankilliany or
%    brainnetome40), or an NxN matrix (where N is the number of regions
%    being used). Default connectomes are described in Tait et al (in press) for
%    brainnetome40, or https://doi.org/10.1016/j.neuroimage.2015.03.055 for
%    the desikankilliany. 
% 
% pars.dynamic_model: 
%    Defines the dynamic model on each node. Can be either 'hopf'
%    (default), as described in my manuscript, or 'nmm', which uses a
%    jansen-rit style neural mass model with addition of inhibitory self
%    connections (https://doi.org/10.1016/j.neuroimage.2007.05.032). 
%    Alternatively, pars.dynamic_model can be a custom ode function of the 
%    form
%    pars.dynamic_model = @(t,x,pars.model_parameters,W) odefunc(inputs), 
%    where pars.model_parameters is described below and W is a connectivity
%    matrix between regions (usually
%    pars.connectome*pars.model_parameters.G). 
%    At present, all simulations use the Euler-Maruyama stochastic
%    numerical integration scheme. 
% 
% pars.model_parameters: 
%    Defines the parameters of the model. If a custom ode is used, these
%    are largely dependent on the ode. However, model_parameters.dim must
%    be included, which is the number of dimensions for each region. This
%    will be set by default if a default ode is used. 
%    -> If pars.dynamic_model = 'hopf': 
%       -> pars.model_parameters.a: bifurcation parameter, scalar or vector
%          Nx1 (value for each node). Default -0.05 for all nodes.
%       -> pars.model_parameters.f: frequency of oscillations in Hz. 
%          Vector Nx1. Default is U(6,13). 
%       -> pars.model_parameters.std_noise: standard deviation of noise.
%          Scalar or vector Nx1. Default is 0.02 mV. 
%    -> If pars.dynamic_model = 'nmm':
%       -> pars.model_parameters.He: default 4mV
%       -> pars.model_parameters.te: tau_e in paper. default 4ms. 
%       -> pars.model_parameters.Hi: default 32mV
%       -> pars.model_parameters.ti: tau_i in paper. default 16ms. 
%       -> pars.model_parameters.gpe,gep,gpi,gip,gii: gamma 1,2,3,4,5 in 
%          paper. default 128,128,64,64,16 
%       -> pars.model_parameters.r: rho_1 in paper. default 2
%       -> pars.model_parameters.v0: rho_2 in paper. default 1
%       -> pars.model_parameters.u: input to system. default 0
% 
% pars.eeg_electrodes: 
%    Defines the EEG electrodes to be forward mapped to. This should be 
%    either empty (pars.eeg_electrodes = []) to do no forward mapping
%    or a cell array of electrode labels (e.g. {'Fp1','Fp2'}).
%    Default is the 19 clinical channels in the 10-20 system (e.g. our San
%    Marino data)
% 
% pars.simulation_options: 
% -> pars.simulation_options.t_init: amount of time to initialise (default
%    500ms)
% -> pars.simulation_options.t_step: time step for simulation (default 0.05
%    ms for nmm, or 0.5 ms for Hopf)
% -> pars.simulation_options.t_simulate: amount of time for simulation
%    (default 2000 ms)
% 
% pars.noise_function:
% -> Noise is automatically drawn as Gaussian. If set, pars.noise_function
%    can be used to transform this noise. It should have one input and one
%    output. eg. to draw uniform noise, pars.noise_function = @(in)
%    rand(size(in)).
% 
% pars.options: 
% -> pars.options.keep_full_simulation: If true, the full simulation is 
%    output. Otherwise, just the simulated EEG (observed variable) is 
%    output. Default false. 
% -> pars.options.plot_simulation: If true, plot simulation. Default is
%    true. 
% -> pars.options.simulation_video: If true, plot a video of simulation.
%    Extremely time consuming, and doesn't work for custom atlases. Default
%    is false. 
% -> pars.options.log_output: If true, console output will be logged to a
%    file. File name will be 'brainSimulator_YYYY_MM_DD_HH_mm_SS' where
%    YYYY is the year, MM is the month, DD is the day, HH is the hour, mm
%    is the minute, and SS is the seconds at which simulation began. 
% -> pars.options.simulation_video: If true, make video of simulation in 
%    source space. Default is false. NOTE: This is VERY time consuming, so
%    I would generally stick to false wherever possible.
% -> pars.options.rngseed: If you wish to repeat simulations with identical
%    noise inputs, you need to control the seed of the random number
%    generator. This can be a non-negative integer < 2^32. Default []
%    (don't use a specific seed). 
% -> pars.options.plot_eeg: If true, plot EEG. Default is true. 
% 
% pars.leadfield: a NchannelxNregion matrix mapping source space time
%    series to sensor space. Do not set this value to compute the leadfield
%    matrix. Usually only set if a custom atlas is included, as will be
%    automatically calculated for implemented atlases.

%% Versions
% v2: updated 27-03-2019 to address issues with custom ode functions
%     - addressed errors with validating input custom ode
%     - addressed errors with parameter and connectome input
%     - allow pars.model_parameters.observed_variables to be a function
%       handle so that custom observations can be used. e.g. Jansen-Rit
%       model needs x2-x3 to be observed EEG, so update allows this to be
%       the case.
% 
% v3: updated 03-05-2019
%     - fixed bug for pars.options.log_output=false
%     - added option pars.simulation_options.noise_function which takes in
%       Gaussian noise with variance 1, and outputs the noise to be used in
%       the simulation. This can also be used to overwrite Gaussian noise.
% 
% v4: updated 10-05-2019
%     - Updated to allow custom atlases and leadfield matrices, e.g.
%       defined externally through ICA analysis. 
%% Main code script

funcName = 'brainSimulator' ; % global variable just used for output - ignore
out = struct ; % initialize output structure 
if nargin < 1 % no input, make empty pars structure
    pars = [] ; 
end

% Check inputs are valid
pars = validate_pars(pars) ; 

% --- MAIN SCRIPT TO RUN ---
% Simulate source activity
pars = generate_atlas(pars) ; % build brain atlas
pars = generate_connectome(pars) ; % build connectome
pars = generate_odepars(pars) ; % build ode parameters
pars = generate_odefun(pars) ; % get ode function
[out.time,out.source,pars] = simulate_ode(pars) ; % simulate
% Visual output - plots and video
plot_simulation(out.time,out.source,pars) ; % plot the simulation (if pars.options.plot_simulation is true)
simulation_video(out.time,out.source,pars) ; % plot the simulation (if pars.options.simulation_video is true)
% Forward map to sensor activity
pars = generate_leadfield(pars) ; % make leadfield
out.sensor = map_to_eeg(out.source,pars) ; % output sensor space activity
% Visual output - plot
plot_eeg(out.time,out.sensor,pars) ; % plot the simulation (if pars.options.plot_eeg is true)
% output to log and close
pars = make_log_output(pars) ; 
% Save new pars to output
out.pars = pars ; 
% --- END OF MAIN SCRIPT ---

%% ----------------------- Main script functions  -------------------------
% The following sections contain the main functions for running the main
% script, largely in the order in which they are called during the script.
% The various plotting/video functions, validate_pars and log_output are 
% left until the end as these don't play a key role in calculations - the
% code could run without these. 

%% generate_atlas
% Here we load the atlas that was chosen, and extract the regions that we
% wish to simulate. Hence the resulting atlas will have the same labels in
% the same order as the regions chosen in the input. 
% 
% An atlas contains 3 entries: 
% labels: labels of the regions
% atlas2voxels: each region has an entry with a list of all voxels within
% that region
% voxels2atlas: a vector of length number_of_voxels, where each voxel has a
% number that corresponds to a region. 
function pars = generate_atlas(pars)
        
    if ~strcmp(pars.atlas.atlas,'custom')
        % output
        fprintf('Generating atlas... ')
        
        % load default atlases in
        load('default_data_mni','atlas')        
        % select the chosen atlas
        atlas = atlas.(pars.atlas.atlas) ; 

        if isempty(pars.atlas.regions) % default setting, all
            pars.atlas.regions = (1:length(atlas.labels))' ; 
        end
        if isa(pars.atlas.regions,'logical') % convert from logical to numeric
            validateattributes(pars.atlas.regions,{'logical'},{'numel',length(atlas.labels)},funcName,'pars.atlas.regions') % validate format
            pars.atlas.regions = find(pars.atlas.regions) ; % update pars
        end
        validateattributes(pars.atlas.regions,{'numeric'},{'<=',length(atlas.labels)},funcName,'pars.atlas.regions') % validate format
        atlas.labels = atlas.labels(pars.atlas.regions) ; % reorder labels
        atlas.atlas2voxels = atlas.atlas2voxels(pars.atlas.regions) ; % reorder atlas2voxels
        v2a = nan(size(atlas.voxels2atlas)) ; % temporary voxels2atlas
        for i = 1:length(atlas.atlas2voxels) % relabel voxels2atlas
            v2a(atlas.atlas2voxels{i}) = i ; 
        end
        atlas.voxels2atlas = v2a ; % update atlas

        pars.atlas.regions = atlas.labels ; % update pars
        pars.atlas.voxels2atlas = atlas.voxels2atlas ; 
        pars.atlas.atlas2voxels = atlas.atlas2voxels ; 
        
        fprintf('Done \n')
        
    else
        % output
        fprintf('Using user supplied atlas\n ')
    end
    
end

%% generate_connectome
% What we do here depends on whether the input connectome was a default
% (i.e. isa(pars.connectome,'char')=true) or a matrix. If a matrix, we just
% check that it is an NxN square matrix where N is the number of regions of
% the brain - other than that we just trust the user that the ordering of
% regions matches their input. 
% 
% If default, we load the connectome from the default data for the correct
% atlas, then ensure that the labels match those of the atlas regions. 
function pars = generate_connectome(pars)
    
    fprintf('Generating connectome... ')
    
    if isa(pars.connectome,'char') % default connectome
        load('default_data_mni','connectome') % load
        connectome = connectome.(pars.connectome) ; % choose default
        W = connectome.connectome ; % connectivity matrix
        lbl = connectome.labels ; % labels
        
        idx = zeros(length(pars.atlas.regions),1) ; % index to align to atlas regions
        for i = 1:length(pars.atlas.regions)
            try
                idx(i) = find(strcmp(pars.atlas.regions{i},lbl)) ; 
            catch
                error('Atlas and connectome do not align')
            end
        end
        W = W(idx,idx) ; % rearrange connectivity matrix
        pars.connectome = W ; % save as connectome
    else
        % check valid (square) input
        validateattributes(pars.connectome,{'numeric'},{'size',[length(pars.atlas.regions),length(pars.atlas.regions)],'finite','nonnan'},funcName,'pars.connectome')
    end
    
    fprintf('Done \n')
end

%% generate_odepars
% This function generates the parameters for a given dynamic model. For
% each model, it first deals with "general parameters", which are
% parameters that apply to all models. These include "dim", which is the
% dimensionality (number of variable) of a single oscillator, "N" which
% is the number of regions, and "G" which is a global scaling of the
% connectome. "Model specific" parameters are ones that are specific to the
% given dynamic model. 
% 
% Whilst this section looks long and scary, pretty much all parameters
% follow the same format: 
% if ~isfield(odepars,'X')                          <--- if you haven't included X in input, make default
%    odepars.X = defaultX*ones(odepars.N,1);        <--- set it to default, which is always heterogeneous, hence multiplying by ones(N,1) (same value for N regions)
%    else                                           <--- if you have included X in input, we need to check you have a value for each region and this is valid
%       if isscalar(odepars.X)                      <--- if you only put one value in, we are going to assume this value applies to all nodes heterogeneously
%          odepars.X=repmat(odepars.X,odepars.N,1); <--- repeat the value N times, i.e. once for each region
%       end
%       validateattributes(odepars.X,...)           <--- make sure the input is valid, e.g. if it was a string or had a length not equal to N, it is not valid
%    end
% end

function pars = generate_odepars(pars)
    
    fprintf('Generating ODE parameters... ')
    
    odepars = pars.model_parameters ;
    if ~isa(pars.dynamic_model,'function_handle')
        switch pars.dynamic_model
            case 'hopf'
                % general parameters
                odepars.dim = 2 ; % number of dimensions per region
                odepars.N = length(pars.atlas.regions) ; % number of regions
                if ~isfield(odepars,'G')
                    odepars.G = 1 ; % global coupling
                end
                validateattributes(odepars.G,{'numeric'},{'scalar'},funcName,'pars.model_parameters.G')
                % model specific
                % bifurcation parameter
                if ~isfield(odepars,'a') 
                    odepars.a = -0.05*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.a)
                        odepars.a = repmat(odepars.a,odepars.N,1) ; 
                    end
                    validateattributes(odepars.a,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.a')
                end
                % oscillator frequency
                if ~isfield(odepars,'f') 
                    if ~isempty(pars.options.rngseed) % random number generator seed if used
                        rng(pars.options.rngseed) ; 
                    end
                    odepars.f = 6+7*rand(odepars.N,1) ; % generate uniform in alpha (6-13 Hz) frequency by default
                else
                    validateattributes(odepars.f,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.f')
                end
                pars.model_parameters = odepars; % overwrite model parameters
                % noise variance
                if ~isfield(odepars,'std_noise') 
                    odepars.std_noise = [0.02*ones(odepars.N,1) ; zeros(odepars.N,1)] ; % generate uniform in alpha (8-12 Hz) frequency by default
                else
                    if isscalar(odepars.std_noise)
                        odepars.std_noise = repmat(odepars.std_noise,odepars.N*odepars.dim,1) ; 
                    end
                    validateattributes(odepars.std_noise,{'numeric'},{'size',[odepars.N*odepars.dim,1]},funcName,'pars.model_parameters.std_noise')
                end
                pars.model_parameters = odepars; % overwrite model parameters
                % observed variables
                pars.model_parameters.observed_variables = (1:odepars.N)' ; % overwrite if given

            case 'nmm'
                % general parameters
                odepars.dim = 12 ; 
                odepars.N = length(pars.atlas.regions) ; % number of regions
                if ~isfield(odepars,'G')
                    odepars.G = 1 ; % global coupling
                end
                validateattributes(odepars.G,{'numeric'},{'scalar'},funcName,'pars.model_parameters.G')

                % model specific
                % population gains
                if ~isfield(odepars,'He') 
                    odepars.He = 4*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.He)
                        odepars.He = repmat(odepars.He,odepars.N,1) ; 
                    end
                    validateattributes(odepars.He,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.He')
                end
                if ~isfield(odepars,'Hi') 
                    odepars.Hi = 32*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.Hi)
                        odepars.Hi = repmat(odepars.Hi,odepars.N,1) ; 
                    end
                    validateattributes(odepars.Hi,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.Hi')
                end

                % population time constants
                if ~isfield(odepars,'te') 
                    odepars.te = 4*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.te)
                        odepars.te = repmat(odepars.te,odepars.N,1) ; 
                    end
                    validateattributes(odepars.te,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.te')
                end
                if ~isfield(odepars,'ti') 
                    odepars.ti = 16*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.ti)
                        odepars.ti = repmat(odepars.ti,odepars.N,1) ; 
                    end
                    validateattributes(odepars.ti,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.ti')
                end

                % population connectivity
                if ~isfield(odepars,'gpe') 
                    odepars.gpe = 128*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.gpe)
                        odepars.gpe = repmat(odepars.gpe,odepars.N,1) ; 
                    end
                    validateattributes(odepars.gpe,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.gpe')
                end
                if ~isfield(odepars,'gep') 
                    odepars.gep = 128*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.gep)
                        odepars.gep = repmat(odepars.gep,odepars.N,1) ; 
                    end
                    validateattributes(odepars.gep,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.gep')
                end
                if ~isfield(odepars,'gpi') 
                    odepars.gpi = 64*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.gpi)
                        odepars.gpi = repmat(odepars.gpi,odepars.N,1) ; 
                    end
                    validateattributes(odepars.gpi,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.gpi')
                end
                if ~isfield(odepars,'gip') 
                    odepars.gip = 64*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.gip)
                        odepars.gip = repmat(odepars.gip,odepars.N,1) ; 
                    end
                    validateattributes(odepars.gip,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.gip')
                end
                if ~isfield(odepars,'gii') 
                    odepars.gii = 16*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.gii)
                        odepars.gii = repmat(odepars.gii,odepars.N,1) ; 
                    end
                    validateattributes(odepars.gii,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.gii')
                end

                % sigmoid parameters
                if ~isfield(odepars,'r') 
                    odepars.r = 2*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.r)
                        odepars.r = repmat(odepars.r,odepars.N,1) ; 
                    end
                    validateattributes(odepars.r,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.r')
                end
                if ~isfield(odepars,'v0') 
                    odepars.v0 = 1*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.v0)
                        odepars.v0 = repmat(odepars.v0,odepars.N,1) ; 
                    end
                    validateattributes(odepars.v0,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.v0')
                end

                % input
                if ~isfield(odepars,'u') 
                    odepars.u = 0*ones(odepars.N,1) ; 
                else
                    if isscalar(odepars.u)
                        odepars.u = repmat(odepars.u,odepars.N,1) ; 
                    end
                    validateattributes(odepars.u,{'numeric'},{'size',[odepars.N,1]},funcName,'pars.model_parameters.u')
                end

                % noise
                if ~isfield(odepars,'std_noise')
                    odepars.std_noise = [zeros(odepars.N*7,1) ; % v1-v7
                                         (1/16)*ones(odepars.N,1) ; % z1
                                         zeros(odepars.N*4,1)]; % z2-z5
                end

                pars.model_parameters = odepars ; 

                % observed variables
                pars.model_parameters.observed_variables = (5*odepars.N+1:6*odepars.N)' ; % overwrite if given
        end
    else % custom function
        pars.model_parameters.N = length(pars.atlas.regions) ; % number of regions
        if isa(pars.model_parameters.observed_variables,'function_handle')
            if nargin(pars.model_parameters.observed_variables) ~= 1 && nargout(pars.model_parameters.observed_variables) ~= 1
                error('If pars.dynamic_model.observed_variables is a function handle, it must have one input (the full simulation) and one output (the observation)')
            end
        else
            validateattributes(pars.model_parameters.observed_variables,{'numeric'},{'positive','integer','size',[pars.model_parameters.N,1],'<=',pars.model_parameters.N*pars.model_parameters.dim},funcName,'pars.model_parameters.observed_variables')
        end
    end
    
    fprintf('Done \n')
end

%% generate_odefun
% This is pretty simple, it just maps a string input of 'hopf' or 'nmm' to
% a function handle that gives the actual dynamic model, where the
% functions are in the next section of code. 
% 
% Below this are the functions containing ODEs for each of the default 
% models. They take input "t", which is the time, "X" which is the vector 
% of variables, "p" which is the structure of parameters, and "W" which is 
% the connectivity, and output "Xdot" which is the time derivative of "X". 
function pars = generate_odefun(pars)
    
    fprintf('Generating ODE function... ')
    if ~isa(pars.dynamic_model,'function_handle')
        switch pars.dynamic_model
            case 'hopf'
                pars.odefun = @(t,x) hopfODE(t,x,pars.model_parameters,pars.model_parameters.G*pars.connectome) ; % need to make ode function with parameters
            case 'nmm'
                pars.odefun = @(t,x) nmmODE(t,x,pars.model_parameters,pars.model_parameters.G*pars.connectome) ; 
        end
    else
        pars.odefun = @(t,x) pars.dynamic_model(t,x,pars.model_parameters,pars.connectome) ; 
    end
    
    fprintf('Done \n')
end

% HOPF ---------------
function Xdot = hopfODE(t,X,p,W)
    N = length(X) ; % dimensionality of full system
    x = X(1:N/2) ; % x - simulated activity
    y = X(N/2+1:N) ; % y - recovery currents

    xdot = (p.a-x.^2-y.^2).*x - 2*pi*p.f.*y + sum(W.*(x'-x),2) ; 
    ydot = (p.a-x.^2-y.^2).*y + 2*pi*p.f.*x ; 

    Xdot = [xdot ; ydot];
end

% NMM --------------
function Xdot = nmmODE(t,X,p,W)
    
    % Initialize Xdot (12 dimensions for each region as there are a total
    % of 12 v and x equations)
    Xdot = zeros(12*p.N,1) ; 
    
    % Get variables - v and x (which is vot) in the Moran paper, but we
    % call 'x' z here. 
    v1 = X(1:p.N) ; % PSP on stellate cells
    v2 = X(p.N+1:2*p.N) ; % PSP on pyramidal cells due to stellate cells
    v3 = X(2*p.N+1:3*p.N) ; % PSP on pyramidal cells due to inhibitory cells
    v4 = X(3*p.N+1:4*p.N) ; % PSP on inhibitory cells due to pyramidal cells
    v5 = X(4*p.N+1:5*p.N) ; % PSP on inhibitory cells due to inhibitory cells
    v6 = X(5*p.N+1:6*p.N) ; % Sum of PSPs on pyramidal (EEG)
    v7 = X(6*p.N+1:7*p.N) ; % Sum of PSPs on inhibitory (EEG)
    
    z1 = X(7*p.N+1:8*p.N) ; % derivative of v1
    z2 = X(8*p.N+1:9*p.N) ; % derivative of v2
    z3 = X(9*p.N+1:10*p.N) ; % derivative of v3
    z4 = X(10*p.N+1:11*p.N) ; % derivative of v4
    z5 = X(11*p.N+1:12*p.N) ; % derivative of v5
    
    % ODE for PSPs - this is the set of equations vdot=x(t) in Moran paper
    Xdot(1:p.N) = z1 ; % v1dot
    Xdot(p.N+1:2*p.N) = z2 ; % v2dot
    Xdot(2*p.N+1:3*p.N) = z3 ; % v3dot
    Xdot(3*p.N+1:4*p.N) = z4 ; % v4dot
    Xdot(4*p.N+1:5*p.N) = z5 ; % v5dot
    Xdot(5*p.N+1:6*p.N) = z2-z3 ; % v6dot
    Xdot(6*p.N+1:7*p.N) = z4-z5; % v7dot
    
    % Define sigma function
    Sig = @(v) 1./(1+exp(-p.r.*(v-p.v0))) - 1./(1+exp(p.r.*p.v0)) ; 
    
    % ODE for derivatives of PSPs - this is the set of equations
    % xdot=He/i*ke/i*u(t)-2*ke/i*x(t)-ke/i^2*v(t) in Moran paper
    Xdot(7*p.N+1:8*p.N) = (p.He./p.te).*(p.gpe.*Sig(v6)+ p.u + W*v6)-(2./p.te).*z1-((1./p.te).^2).*v1 ; % z1dot
    Xdot(8*p.N+1:9*p.N) = (p.He./p.te).*p.gep.*Sig(v1)-(2./p.te).*z2-((1./p.te).^2).*v2  ; % z2dot
    Xdot(9*p.N+1:10*p.N) = (p.Hi./p.ti).*p.gip.*Sig(v7)-(2./p.ti).*z3-((1./p.ti).^2).*v3 ; % z3dot
    Xdot(10*p.N+1:11*p.N) = (p.He./p.te).*p.gpi.*Sig(v6)-(2./p.te).*z4-((1./p.te).^2).*v4  ; % z4dot
    Xdot(11*p.N+1:12*p.N) = (p.Hi./p.ti).*p.gii.*Sig(v7)-(2./p.ti).*z5-((1./p.ti).^2).*v5 ; % z5dot
    
    % For a single region (p.N=1), Xdot(1:7) is vdot in Moran paper and
    % Xdot(8:12) is xdot in the Moran paper. 
end

%% simulate_ode
% This function is the core of running the simulation. It extracts the time
% axis, the initializes all variables at 0. Then a stochastic
% Euler-Maryuama solver is run to solve the system at each of the time
% points. Finally, if pars.options.keep_full_simulation is false, only the
% observation variables (corresponding to the LFP/EEG) are kept.  
function [t,x,pars] = simulate_ode(pars)
    
    %%% --- Get simulation options --- %%% 
    fprintf('Generating simulation parameters... ')
    dt = pars.simulation_options.t_step ; % time step
    t = -pars.simulation_options.t_init:dt:pars.simulation_options.t_simulate ; % time axis
    idxt = find(t>=0) ; % index of points for t>=0 (i.e. the simulation time, not the initialization time)
    Nt = length(t) ; % number of time points
    odefun = pars.odefun ; % function to simulate
    % in nmm model, noise input is scaled by He and te
    if strcmp(pars.dynamic_model,'nmm')
        noise_scaling = [zeros(pars.model_parameters.N*7,1) ; % v1-v7, no noise input
                         pars.model_parameters.He./pars.model_parameters.te ; % z1 
                         pars.model_parameters.He./pars.model_parameters.te ; % z2
                         pars.model_parameters.Hi./pars.model_parameters.ti ; % z3 
                         pars.model_parameters.He./pars.model_parameters.te ; % z4 
                         pars.model_parameters.Hi./pars.model_parameters.ti] ; % z5 
    else
        noise_scaling = ones(pars.model_parameters.N*pars.model_parameters.dim,1) ; 
    end
    % for hopf model, frequency is in Hz so we need time step in seconds
    if strcmp(pars.dynamic_model,'hopf')
        dt = dt/1000 ; 
    end
    fprintf('Done \n')
    
    %%% --- Initialize variables --- %%% 
    fprintf('Initializing variables... ')
    x = zeros(pars.model_parameters.N*pars.model_parameters.dim,Nt) ; 
    x(:,1) = 0.1*randn(pars.model_parameters.N*pars.model_parameters.dim,1) ; 
    fprintf('Done \n')
    
    %%% --- Simulate --- %%%
    fprintf('Simulating... '), tic
    pc_old = 0 ; % percentage complete
    if ~isempty(pars.options.rngseed) % set random number generator seed
        rng(pars.options.rngseed)
    end
    eta = randn(size(x)) ; 
    if ~isempty(pars.simulation_options.noise_function)
        if nargin(pars.simulation_options.noise_function) == 1
            eta = pars.simulation_options.noise_function(eta) ; 
        elseif nargin(pars.simulation_options.noise_function) == 2
            eta = pars.simulation_options.noise_function(eta,pars) ; 
        end
    end
    eta = eta.*noise_scaling.*pars.model_parameters.std_noise ; 
    msg = [] ;
    
    if pars.options.log_output
        diary off % turn log output off for percent complete
    end
    for i = 1:Nt
        % Output to console
        pc_new = round(100*(i/Nt)) ; 
        if pc_new>pc_old
            fprintf(repmat('\b',1,length(msg)))
            msg = sprintf('%d percent complete',pc_new) ; 
            fprintf(msg)
        end
        pc_old = pc_new ; 
        
        % Euler-Maruyama SDE solver
        f0 = odefun(t(i),x(:,i)) ; 
        x(:,i+1) = x(:,i) + dt*f0 + eta(:,i) ; % note dt/1000 as f in Hz 
    end
        
    % remove initialization time
    t = t(idxt) ; 
    x = x(:,idxt) ;
    
    fprintf(repmat('\b',1,length(msg)))
    if pars.options.log_output
        diary(['simulation_log_' pars.run_datetime '.txt']) % turn console output back on
    end
    msg = sprintf('Done\n',pc_new) ; 
    fprintf(msg)
    
    % Output simulation time
    fprintf('Simulation time %.2f seconds\n',toc)
    
    % Observation function
    if ~pars.options.keep_full_simulation
        if isa(pars.model_parameters.observed_variables,'function_handle')
            x = pars.model_parameters.observed_variables(x) ; 
        else
            x = x(pars.model_parameters.observed_variables,:) ; 
        end
    end
    
end

%% generate_leadfield
% The leadfield matrix is the matrix that maps source activity onto EEG
% electrodes as a linear mixture, i.e. leadfield(i,j) is the influence
% source j has on electrode i. 
% 
% This function is a little bit misnamed - it doesn't really calculate a
% new leadfield matrix. The leadfield matrix was precalculated in
% /unused_scripts/0_make_leadfield for each atlas and over 300 electrode
% locations in the 10-05 system. Since leadfield(i,j) is independent of all
% other entries in the leadfield (i.e. adding or taking other sources or
% electrodes will not influence that entry), this function essentially
% takes this precomputed leadfield for all regions and all electrodes and
% simply selects the rows and columns corresponding to the chosen regions
% and electrodes for simulation. 
function pars = generate_leadfield(pars)
        
    if isempty(pars.eeg_electrodes) % case where no forward modelling is to be done
        pars.leadfield = 'none' ; 
        return
    end
    
    if isfield(pars,'leadfield')
        validateattributes(pars.leadfield,{'numeric'},{'real','finite','nonnan','nonempty','nonsparse','size',[length(pars.eeg_electrodes),length(pars.atlas.regions)]})
        fprintf('Using user-supplied leadfield matrix \n ')
        return
    end
    
    fprintf('Making leadfield matrix... ')
    % load electrodes, atlases, and leadfield
    load('default_data_mni','electrode','leadfield','atlas')
    lf = leadfield.(pars.atlas.atlas) ; 
    atlas = atlas.(pars.atlas.atlas) ; 
    
    % only keep simulated regions for leadfield
    lbl = atlas.labels ; 
    ord = zeros(length(pars.atlas.regions),1) ;
    for i = 1:length(pars.atlas.regions)
        ord(i) = find(strcmp(pars.atlas.regions{i},lbl)) ; 
    end
    lf = lf(:,ord) ; 
    
    % only keep chosen eeg electrodes for leadfield
    lbl = electrode.elec.label ; 
    ord = zeros(length(pars.eeg_electrodes),1) ; 
    for i = 1:length(pars.eeg_electrodes)
        ord(i) = find(strcmp(pars.eeg_electrodes{i},lbl)) ; 
    end
    lf = lf(ord,:) ; 
    
    % save leadfield matrix to output
    pars.leadfield = lf ; 
    fprintf('Done\n')
end
    
%% map_to_eeg
% This section is very simple since the forward model is stated as: 
% X=L*S 
% where X is the EEG, L is the leadfield, and S is the source activity.
% This mapping is performed in this section. 
function x = map_to_eeg(s,pars)
    
    if isempty(pars.eeg_electrodes) % if not mapping to EEG
        x = 'No electrodes input, so no EEG generated' ; 
        return
    end
    
    if pars.options.keep_full_simulation
        if isa(pars.model_parameters.observed_variables,'function_handle')
            s = pars.model_parameters.observed_variables(s) ; 
        else
            s = s(pars.model_parameters.observed_variables,:) ; 
        end
    end
    fprintf('Forward mapping to EEG... ')
    x = pars.leadfield*s ; % forward model
    fprintf('Done\n')
    
end

%% Visual output
% This section is a collection of functions for plotting including: 
% plot_simulation: plots the source space data
% plot_eeg: plots the sensor space data
% simulation_video: makes a video of the source space data
% print_figure: used by the two plot_xyz functions to make a publication
% quality (500dpi, 10pt font, fit to most A4 pages) figures. 

% PLOT_SIMULATION -------------
function plot_simulation(t,x,pars)
    if ~pars.options.plot_simulation
        return
    end
    
    fprintf('Plotting simulation... ')
    figure
    if pars.options.keep_full_simulation
        x = x(pars.model_parameters.observed_variables,:) ; 
    end
    
    x = x-mean(x,2) ; % zero mean each trace for plotting purposes
    plot(t,x' + 5*std(x(:))*(1:pars.model_parameters.N),'k') ; % plot
    
    ax = gca ; 
    ax.YTick = 5*std(x(:))*(1:pars.model_parameters.N) ; 
    ax.YTickLabels = pars.atlas.regions ; 
    xlabel('Time [ms]')
    ylim([0 5*std(x(:))*(pars.model_parameters.N+1)])
    
    print_figure(['simulation_plot_' pars.run_datetime])
    close
    
    fprintf('Done\n')
end

% PLOT_EEG -------------
function plot_eeg(t,x,pars)
    if ~pars.options.plot_eeg
        return
    end
    
    fprintf('Plotting eeg... ')
    figure
    
    x = x-mean(x,2) ; % zero mean each trace for plotting purposes
    plot(t,x' + 5*std(x(:))*(1:length(pars.eeg_electrodes)),'k') ; % plot
    
    ax = gca ; 
    ax.YTick = 5*std(x(:))*(1:length(pars.eeg_electrodes)) ; 
    ax.YTickLabels = pars.eeg_electrodes ; 
    xlabel('Time [ms]')
    try
        ylim([0 5*std(x(:))*(length(pars.eeg_electrodes)+1)])
    end
    
    print_figure(['eeg_plot_' pars.run_datetime])
    close
    fprintf('Done\n')
end

% SIMULATION_VIDEO --------------
function simulation_video(t,x,pars)
    if ~pars.options.simulation_video
        return
    end
    
    fprintf('Making video (this may be time consuming)... ')
    
    fname = ['video_' pars.run_datetime] ; 
    % get video writer object
    vid = VideoWriter(fname,'MPEG-4') ; 
    vid.FrameRate = 50/(pars.simulation_options.t_step) ; %1/20th
    open(vid)
    
    % get figure
    h = figure ; 
    h.Units = 'centimeters' ; 
    h.Position = [0,0,15,12] ; 
    
    % make view and colormap
    vw = [-99.8 , -4] ; 
    cmap = [ones(32,1) , linspace(0,1,32)' , linspace(0,1,32)' ; flipud(linspace(0,1,32)') , flipud(linspace(0,1,32)') , ones(32,1)] ;
    
    % load head
    load('default_data_mni','headmodel') ; 
    
    % max and minimum
    range = 3*std(x(:)) ; 
    
    % make video
    for i = 1:length(t)
        clf % clear current figure
        
        % plot head
        patch('faces',headmodel.head.tri,'vertices',headmodel.head.pos,'linestyle','none','facealpha',0.1)
        
        % plot brain 
        X = zeros(length(pars.atlas.voxels2atlas),1) ;
        for j = 1:length(pars.atlas.atlas2voxels)
            X(pars.atlas.atlas2voxels{j}) = x(j,i) ; 
        end
        patch('faces',headmodel.cortex.tri,'vertices',headmodel.cortex.pos,'linestyle','none','facecolor','interp','facevertexCdata',X) ; 
        
        % visual
        ax = gca ; ax.View = vw ;
        camlight
        caxis([-range,range])
        colormap(cmap)
        title(sprintf('t = %.2f ms',t(i)))
        axis off
        drawnow
        
        % writevideo
        writeVideo(vid,getframe(h))
    end
    
    close(vid)
end

% PRINT_FIGURE
function print_figure(fname)

    fontsize = 10 ; 

    h = gcf ;
    % window appearance\\\
    h.MenuBar = 'figure' ; 
    h.ToolBar = 'figure' ; 
    % h.DockControls = 'off' ; 
    h.Color = 'white' ; 
    h.WindowStyle = 'normal' ; 

    % window position
    h.Units = 'centimeters' ;
    size = [15,12] ; 
    h.Position = [0,0,size] ; 
    h.Resize = 'on' ; 

    % check font sizes
    ax = h.Children ; 
    for i = 1:length(ax)
        axi = ax(i) ; 
        if strcmp(class(axi),'matlab.graphics.axis.Axes')
            axi.FontName = 'Helvetica' ; 
            axi.FontUnits = 'points' ; 
            axi.FontSize = fontsize ; 
        end
    end

    % set paper size
    h.PaperUnits = 'centimeters' ; 
    h.PaperSize = size ; 
    h.PaperPositionMode = 'manual' ; 
    h.PaperPosition = [0,0,size] ;
    h.InvertHardcopy = 'off' ; 

    % get format type - png with 500 dpi resolution
    format = '-dpng' ; 
    resolution = '-r500' ;

    print(fname,format,resolution)
    savefig(fname)

end

%% Function to validate parameters input
% Outline of this section: 
% - Make log so that console output is printed. This is mainly useful for
%   tracking errors etc.
% - Check that inputs are valid. This is basically to try and ensure errors 
%   are caught sooner rather than later after computation time. 
% - Make default parameters. All inputs are required for the code to run,
%   so those that aren't included in the input need to be set to default.
% ESSENTIALLY YOU CAN IGNORE THIS SECTION
function pars = validate_pars(pars)

    % Make log
    fname = char(datetime('now','Format','yyyy_MM_dd_HH_mm_ss')) ; 
    diary(['simulation_log_' char(fname) '.txt'])
    
    % OUTPUT: 
    fprintf('--- Running Brain Simulator ---\n\n')
    fprintf('Validating input parameters... ')

    % Check pars structure
    if isempty(pars) % make default
        pars = struct ; 
    end
    validateattributes(pars,{'struct'},{'scalar'},funcName,'pars') ; % validate format

    % Check inputs to pars

    % pars.atlas
    if ~isfield(pars,'atlas') % make default
        pars.atlas = struct ; 
    end
    validateattributes(pars.atlas,{'struct'},{'scalar'},funcName,'pars.atlas') % validate format
    
    % pars.atlas.atlas
    if ~isfield(pars.atlas,'atlas')
        pars.atlas.atlas = 'brainnetome40' ; 
    end
    validateattributes(pars.atlas.atlas,{'string','char'},{'scalartext'},funcName,'pars.atlas.atlas') % validate format
    if ~(any(strcmp(pars.atlas.atlas,{'destrieux','desikankilliany','AAL','brainnetome40','custom'})))
        error("pars.atlas.atlas must be 'custom', 'destrieux', 'desikankilliany', 'AAL', 'brainnetome40' (default)") % check valid inputs
    end
    switch pars.atlas.atlas
        case 'destrieux'
            pars.temp.atlascitation = 'Destrieux et al. (2010) https://dx.doi.org/10.1016%2Fj.neuroimage.2010.06.010' ; 
        case 'desikankilliany'
            pars.temp.atlascitation = 'Desikan et al. (2006) https://doi.org/10.1016/j.neuroimage.2006.01.021' ; 
        case 'AAL'
            pars.temp.atlascitation = 'Tzourio-Mazoyer et al. (2002) https://doi.org/10.1006/nimg.2001.0978' ; 
        case 'brainnetome40'
            pars.temp.atlascitation = 'Tait et al. (2019) Under review at Clinical Neurophysiology' ; 
        case 'custom'
            pars.temp.atlascitation = 'Custom atlas used' ; 
    end
    
    % pars.atlas.regions
    if ~isfield(pars.atlas,'regions')
        pars.atlas.regions = [] ; 
    end
    if ~isempty(pars.atlas.regions)
        validateattributes(pars.atlas.regions,{'numeric','logical'},{'vector','integer','positive'},funcName,'pars.atlas.regions') % validate format
    end
    if isempty(pars.atlas.regions) && strcmp(pars.atlas.atlas,'custom')
        error('When using a custom atlas user must define regions of interest, defaults not available')
    end
      
    % pars.connectome
    if ~isfield(pars,'connectome')
        pars.connectome = [] ; 
    end
    % make defaults
    if isempty(pars.connectome)
        switch pars.atlas.atlas
            case 'desikankilliany'
                pars.connectome = 'desikankilliany' ;
                pars.temp.connectomedescription = 'The connectome used was the default "tvb_connectivity_68" from The Virtual Brain. It was derived using the tractography data and the automated pipeline described in [Schirner et al. (2015) https://doi.org/10.1016/j.neuroimage.2015.03.055]' ; 
            case 'brainnetome40'
                pars.connectome = 'brainnetome40' ; 
                pars.temp.connectomedescription = 'The connectome used was the average healthy connectome described in [Tait et al. (2019) Under review with Clinical Neurophysiology]. This connectome was anatomically constrained by (binary) knowledge of known connections in the brain [Fan et al. (2016) https://dx.doi.org/10.1093%2Fcercor%2Fbhw157], with weights adjusted such that the "hopf" model gave realistic functional connectivity patterns.' ; 
            otherwise
                error('Default connectomes are only included with atlases desikankilliany, or brainnetome40')
        end
        pars.temp.isdefaultconnectome = true ; 
    % custom
    else
        pars.temp.isdefaultconnectome = false ; 
        validateattributes(pars.connectome,{'numeric'},{'square'},funcName,'pars.connectome') ; 
    end
      
    % pars.dynamic_model
    if ~isfield(pars,'dynamic_model')
        pars.dynamic_model = 'hopf' ; 
    end
    if isa(pars.dynamic_model,'function_handle') % custom function
        if nargin(pars.dynamic_model) ~= 4
            error('Custom dynamic models must have four inputs, being of the form xdot = func(t,x,model_parameters,connectivity_matrix)')
        end
    else
        validateattributes(pars.dynamic_model,{'string','char'},{'scalartext'},funcName,'pars.dynamic_model') % validate format
        if ~(any(strcmp(pars.dynamic_model,{'hopf','nmm'})))
            error("pars.dynamic_model must be 'hopf' (default), 'nmm' or a function handle") % check valid inputs
        end
        switch pars.dynamic_model
            case 'hopf'
                pars.temp.modeldescription = 'Each node was modelled as the normal form of a supercritical Hopf bifurcation. The model is implemented as described in [Tait et al. (2019) Under review with Clinical Neurophysiology], and was based on that of [Deco et al. (2017) https://dx.doi.org/10.1038%2Fs41598-017-03073-5]. For parameter "a" <= 0, the dynamics of the system are a steady state with resonant frequency "f" [Hz], meaning in the stochastic system driven by white noise the spectrum has a peak at "f" Hz. For "a" > 0, there are deterministic oscillations at "f" Hz.' ; 
            case 'nmm'
                pars.temp.modeldescription = 'Each node was modelled as a neural mass model consisting of three populations, implemented as described in [Moran et al. (2013) https://dx.doi.org/10.3389%2Ffncom.2013.00057; "LFP model"]. It is based on the Jansen-Rit formulation [Jansen and Rit (1995) https://doi.org/10.1007/BF00199471]. Pyramidal cells gained excitatory input from L4 stellate cells and supragranular inhibitory interneurons. The stellate populations and inhibitory populations both gain input from the pyramidal population, with no stellate <-> inhibitory connectivity. Following the modifications made by [Moran et al. (2007) https://doi.org/10.1016/j.neuroimage.2007.05.032], recurrent inhibition was also included in the model.' ; 
        end
    end
    
    
    % pars.model_parameters
    if ~isfield(pars,'model_parameters')
        pars.model_parameters = struct ; 
    end
    validateattributes(pars.model_parameters,{'struct'},{'scalar'},funcName,'pars.model_parameters') % validate format
    
    % pars.model_parameters - custom ode
    if isa(pars.dynamic_model,'function_handle')
        if ~isfield(pars.model_parameters,'dim') || ~isfield(pars.model_parameters,'observed_variables')
            error('When the dynamic model is a function handle, pars.model_parameters must have a field dim which is the number of dimensions per region in the model.')
        end
        validateattributes(pars.model_parameters.dim,{'numeric'},{'scalar','integer','positive'},funcName,'pars.model_parameters.dim')
    end
    
    % pars.simulation_options
    if ~isfield(pars,'simulation_options')
        pars.simulation_options = struct ; 
    end
    validateattributes(pars.simulation_options,{'struct'},{'scalar'},funcName,'pars.simulation_options') % validate format
    
    % pars.simulation_options.t_init
    if ~isfield(pars.simulation_options,'t_init')
        pars.simulation_options.t_init = 500 ; 
    end
    validateattributes(pars.simulation_options.t_init,{'numeric'},{'scalar','nonnegative'},funcName,'pars.simulation_option.t_init')
    
    % pars.simulation_options.t_step
    if ~isfield(pars.simulation_options,'t_step')
        if strcmp(pars.dynamic_model,'nmm')
            pars.simulation_options.t_step = 0.05 ; 
        else
            pars.simulation_options.t_step = 0.5 ; 
        end
    end
    validateattributes(pars.simulation_options.t_step,{'numeric'},{'scalar','nonnegative','<',pars.simulation_options.t_init},funcName,'pars.simulation_option.t_step')
    
    % pars.simulation_options.t_simulate
    if ~isfield(pars.simulation_options,'t_simulate')
        pars.simulation_options.t_simulate = 2000 ; 
    end
    validateattributes(pars.simulation_options.t_simulate,{'numeric'},{'scalar','nonnegative','>',pars.simulation_options.t_step},funcName,'pars.simulation_option.t_simulate')
    
    % pars.simulation_options.noise_function
    if ~isfield(pars.simulation_options,'noise_function')
        pars.simulation_options.noise_function = [] ; 
    end
    if ~isempty(pars.simulation_options.noise_function)
        validateattributes(pars.simulation_options.noise_function,{'function_handle'},{},funcName,'pars.simulation_option.noise_function')
    end
    
    % pars.options
    if ~isfield(pars,'options')
        pars.options = struct ; 
    end
    validateattributes(pars.options,{'struct'},{'scalar'},funcName,'pars.options') % validate format
    
    % pars.options.keep_full_simulation
    if ~isfield(pars.options,'keep_full_simulation')
        pars.options.keep_full_simulation = false ; 
    end
    validateattributes(pars.options.keep_full_simulation,{'logical'},{'scalar'},funcName,'pars.options.keep_full_simulation')
    
    % pars.options.plot_simulation
    if ~isfield(pars.options,'plot_simulation')
        pars.options.plot_simulation = true ; 
    end
    validateattributes(pars.options.plot_simulation,{'logical'},{'scalar'},funcName,'pars.options.plot_simulation')
    
    % pars.options.simulation_video
    if ~isfield(pars.options,'simulation_video')
        pars.options.simulation_video = false ; 
    end
    validateattributes(pars.options.simulation_video,{'logical'},{'scalar'},funcName,'pars.options.simulation_video')
    
    % pars.options.log_output
    if ~isfield(pars.options,'log_output')
        pars.options.log_output = true ; 
    end
    validateattributes(pars.options.log_output,{'logical'},{'scalar'},funcName,'pars.options.log_output')
    if ~pars.options.log_output
        diary off
        delete(['simulation_log_' char(fname) '.txt'])
    end
    pars.run_datetime = fname ; % file diary/video output
    
    % pars.options.rngseed
    if ~isfield(pars.options,'rngseed')
        pars.options.rngseed = [] ; 
    end
    validateattributes(pars.options.rngseed,{'numeric'},{'nonnegative','integer','<',2^32},funcName,'pars.options.rngseed')
    
    % pars.eeg_electrodes
    if ~isfield(pars,'eeg_electrodes')
        pars.eeg_electrodes = {'Fp1';'Fp2';'F7';'F3';'Fz';'F4';'F8';'T3';'C3';'Cz';'C4';'T4';'T5';'P3';'Pz';'P4';'T6';'O1';'O2'} ; 
    end
    if ~isempty(pars.eeg_electrodes) % check all valid labels
        validateattributes(pars.eeg_electrodes,{'cell'},{},funcName,'pars.eeg_electrodes')
        pars.eeg_electrodes = pars.eeg_electrodes(:) ; % ensure column vector
        
        if ~isfield(pars,'leadfield') % electrodes not in the 10-20 system allowed, but a custom leadfield must be supplied
            load('default_data_mni','electrode')
            lbl = electrode.elec.label ; clear electrode
            for i = 1:length(pars.eeg_electrodes)
                isgood = any(strcmp(pars.eeg_electrodes{i},lbl),2) ; 
                if ~isgood
                    error('Electrode label %s not valid',pars.eeg_electrodes{i})
                end
            end
        end
    end
    
    % note - if atlas is custom + we wish to forward map, user must supply
    % a leadfield matrix
    if strcmp(pars.atlas.atlas,'custom') && ~isempty(pars.eeg_electrodes) && ~isfield(pars,'leadfield')
        error('If custom atlas is used, a leadfield matrix must be supplied to forward map to EEG')
    end
        
    % pars.options.plot_eeg
    if ~isfield(pars.options,'plot_eeg')
        pars.options.plot_eeg = true ; 
    end
    validateattributes(pars.options.plot_eeg,{'logical'},{'scalar'},funcName,'pars.options.plot_eeg')
    if isempty(pars.eeg_electrodes)
        pars.options.plot_eeg = false ; % no forward model, so cannot plot
    end
   
    % output
    fprintf('Done\n')
end

%% Log output

function pars = make_log_output(pars)
    if ~pars.options.log_output
        pars = rmfield(pars,'temp') ; 
        return
    end
    pause(1)
    
    fprintf('\n\n--- Overview of simulation options ---\n\n')
    
    % Atlas
    fprintf('- ATLAS - \n')
    fprintf('Atlas: %s\n',pars.atlas.atlas)
    fprintf('More info: %s\n',pars.temp.atlascitation)
    if isa(pars.dynamic_model,'function_handle')
        fprintf('%d regions used in simulations:\n',pars.model_parameters.N)
        disp(pars.atlas.regions)
    else
        fprintf('%d regions used in simulations (outlined in table below).\n',pars.model_parameters.N)
    end
    
    % Connectome
    fprintf('\n- CONNECTOME - \n')
    if pars.temp.isdefaultconnectome
        fprintf('Connectome: %s\n',pars.atlas.atlas)
        fprintf('More info: %s\n',pars.temp.connectomedescription)
    else
        fprintf('Connectome: Custom\n')
    end
    
    % Dynamic model
    fprintf('\n - DYNAMIC MODEL - \n')
    if isa(pars.dynamic_model,'function_handle')
        fprintf('Dynamic model: Custom\n')
    else
        fprintf('Dynamic model: %s\n',pars.dynamic_model)
        fprintf('More info: %s\n',pars.temp.modeldescription)
    end
    
    % Model parameters
    if ~isa(pars.dynamic_model,'function_handle')
        switch pars.dynamic_model
            case 'hopf'
                c = cell(0) ; 
                if iscell(pars.atlas.regions)
                    c(:,1) = pars.atlas.regions ; 
                else
                    c(:,1) = num2cell(pars.atlas.regions) ;
                end
                for i = 1:pars.model_parameters.N
                    c{i,2} = pars.model_parameters.a(i) ; 
                    c{i,3} = pars.model_parameters.f(i) ; 
                end
                c = [{'Region Label','a','f [Hz]'} ; c] ; 
            case 'nmm' 
                c = cell(0) ; 
                if iscell(pars.atlas.regions)
                    c(:,1) = pars.atlas.regions ; 
                else
                    c(:,1) = num2cell(pars.atlas.regions) ;
                end
                for i = 1:pars.model_parameters.N
                    c{i,2} = pars.model_parameters.He(i) ; 
                    c{i,3} = pars.model_parameters.Hi(i) ;
                    c{i,4} = pars.model_parameters.te(i) ;
                    c{i,5} = pars.model_parameters.ti(i) ;
                    c{i,6} = pars.model_parameters.gpe(i) ;
                    c{i,7} = pars.model_parameters.gep(i) ; 
                    c{i,8} = pars.model_parameters.gpi(i) ;
                    c{i,9} = pars.model_parameters.gip(i) ;
                    c{i,10} = pars.model_parameters.gii(i) ;
                    c{i,11} = pars.model_parameters.r(i) ;
                    c{i,12} = pars.model_parameters.v0(i) ;
                    c{i,13} = pars.model_parameters.u(i) ;
                end
                c = [{'Region Label','He [mV]','Hi [mV]','tau_e [ms]','tau_i [ms]','g_pe','g_ep','g_pi','g_ip','g_ii','r','v0','u'} ; c] ;
        end
        fprintf('Model parameters:   \n')
        disp(c)
    end
    
    % Simulation 
    fprintf('\n- NUMERICAL INTEGRATION - \n')
    fprintf('Dynamics were simulated for %d ms. Numerical integration used a stochastic Euler-Maryuama method with step size %.3f ms. Simulations were initialized with all variables equal to zero, and %d ms of simulation was ran and disregarded to avoid effects of transients.\n',pars.simulation_options.t_simulate,pars.simulation_options.t_step,pars.simulation_options.t_init) 
    
    % Forward model
    fprintf('\n- FORWARD MODELLING TO EEG - \n')
    fprintf('Dynamics were mapped to %d EEG electrodes (described in the table below). The cortical surface and 3 boundary element method (BEM) shells (scalp, skull, brain) were extracted from the Brainstorm software package [https://neuroimage.usc.edu/brainstorm/], using the ICBM152_2016 template model, which is a non-linear average MRI over 152 subjects [Fonov et al. (2009) http://dx.doi.org/10.1016/S1053-8119(09)70884-5]. 8004 sources were distributed along the cortex, and oriented normally to the surface. For each region of the brain, the centroid and average orientation of all source points within the region were used to construct a leadfield matrix from the region to each electrode using the OpenMEEG BEM algorithm [Gramfort et al. (2011) http://dx.doi.org/10.1155/2011/923703] implemented in the Fieldtrip software [Oostenveld et al.(2011)http://dx.doi.org/10.1155/2011/156869].\n',length(pars.eeg_electrodes))
    fprintf('EEG electrodes:\n')
    disp(pars.eeg_electrodes)
    
    diary off
    pars = rmfield(pars,'temp') ; 
end

%% End function
end