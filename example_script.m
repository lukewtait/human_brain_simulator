%% Outline
% This script contains 4 example scripts, described below. 
% 1) How to run with just the default settings. 
% 2) How to run using prebuilt models and atlases, but using custom
% settings.
% 3) How to simulate sensor space dynamics (instead of source dynamics)
% 4) A more complicated example. Using entirely custom settings, including
% a custom system of ODEs (the Jansen-Rit equations), and a custom atlas,
% we simulated the evoked response of the visual cortex. 
addpath(pwd) ; % ensure directory with brainSimulator.m is on path

%% Example script 1: All default settings
clear, clc, close all

cd('example_run/example1_default')
pars = struct ; % here, we are making an empty pars structure. Whilst pars has lots of fields, any empty field is set to its default value. Therefore by not adding fields to pars, we are entirely default. 
out = brainSimulator(pars) ; % runs brain simulator with all default settings. 
save('out','-struct','out') % save output
cd('../..') % return to main directory

% note the structure out has 4 fields: 
% time: the time axis of the simulation
% source: source space activity
% sensor: sensor space (EEG) activity
% pars: the full 'default' pars field. This is updated with defaults etc as
% the script runs, so this is the final version and has all information on
% the simulations you could want. 

% There were a number of other outputs: 
% simulation_log_[datetime].txt gives a log of time simulation and a number
% of key details. This is a really useful summary of all the setting etc in
% the modelling. IMPORTANT: This also has relevant citations for if you use
% this in a paper. 
% simulation_plot_[datetime].png/fig are figures plotting the source space
% simulation time series
% eeg_plot_[datetime].png/fig are figures plotting the sensor (EEG) space
% time series

%% Example script 2: Custom settings using prebuilt models
clear, clc, close all
% To fully simulate, we need a handful of things: 
% 1) Define where in the brain we want to simulate. To do so, we must 
% define an atlas of brain regions, and choose the regions we want to use. 
% 2) Define model dynamics to simulate. This tells us how each region acts 
% independently.  To do so, we must choose an ODE formulation (dynamic 
% model) and model parameters.
% 3) Define a connectome. This tells us how each region influences the
% others. 
% 4) Choose simulation options like time steps and length of simulation
% 5) Choose EEG electrodes to record the activity on the scalp

% (1) Define atlas
% atlas structure - describes the atlas, has 2 fields
atlas = struct ; 
atlas.atlas = 'desikankilliany' ; % use the 'desikankilliany' atlas, this could also be 'destrieux','AAL', or 'brainnetome40' (default). 
atlas.regions = randperm(68,10) ; % here, we are choosing a random 10 of the 68 regions in the desikankilliany atlas to simulate. If you want to do all, you can either just not include this field (default), or use atlas.regions = [] ; 

% (2) Model dynamics
% dynamic model
dynamic_model = 'nmm' ; % here, we are using the neural mass model dynamic model instead of the default hopf model. Note dynamics are much faster than before - default parameters are set to generate gamma rhythms

% model parameters structure - described the model parameters. Number of
% fields depends on the model. Here we will change the time constants and
% the global coupling. 
model_parameters = struct ; 
model_parameters.te = 16 ; % slow down the excitatory time constant (default is 4)
model_parameters.ti = 16+16*rand(10,1) ; % slow down the inhibitory time constant. Note we have set this to uniformly distrubted between 16 and 32. Most parameters can be heterogeneous for nodes. If you choose a scalar, it assumes homogeneous parameter values.  (default is 16)
model_parameters.G = 20 ; % global coupling constant, set very high - we expect hypersynchrony (default is 1)

% (3) Connectome 
% connectome
connectome = rand(10).*~eye(10) ; % here we are generating a random 10x10 matrix to use as our connectome. You can input any custom connectome, or if using the 'desikankilliany' or 'brainnetome40' atlases you can also leave this blank or set to string 'default'

% (4) Simulation options
% simulation options - describes the simulation times. Has 3 fields.
simulation_options = struct ; 
simulation_options.t_init = 1000 ; % simulate for 1000 ms to initialize variables, this is thrown away (default is 500ms)
simulation_options.t_step = 0.01 ; % time step of 0.01 ms (default is 0.05 ms for nmm, 0.5ms for hopf)
simulation_options.t_simulate = 5000 ; % simulate for 5 seconds (default is 2 s)

% (5) EEG electrodes
eeg_electrodes = {'Fp1';'Fp2';'F7';'F3';'Fz';'F4';'F8';'T3';'C3';'Cz';'C4';'T4';'T5';'P3';'Pz';'P4';'T6';'O1';'O2'} ; % this is the default setting, which is the 19 channels commonly used in clinic. You can use a combination of over 300 channels in the 10-05 system. 

% additional options
options = struct ; 
options.keep_full_simulation = false ; % this is default, and it means we only keep the EEG-like variable. If set to true, we keep all variables.
options.plot_simulation = false ; % by default this is true, but we are setting to false. This means it will not plot the source space data
options.simulation_video = false ; % you can set this to true, but it's very time consuming (default is false)
options.log_output = true ; % this gives a lot of informaiton about the simulations in the text file, including citations etc (default is true)
options.rngseed = 123 ; % we can choose to set a seed for the random number generator for reproducibility. Leave this blank or set to empty to not choose a seed. 
options.plot_eeg = true ; % plot the EEG. (default true)

% compile into the pars structure
pars = struct ; 
pars.atlas = atlas ; 
pars.dynamic_model = dynamic_model ; 
pars.model_parameters = model_parameters ; 
pars.connectome = connectome ; 
pars.simulation_options = simulation_options ; 
pars.eeg_electrodes = eeg_electrodes ; 
pars.options = options ; 

% run the code
cd('example_run/example2_custom')
out = brainSimulator(pars) ; 
save('out','-struct','out') % save output
cd('../..') % return to main directory

%% Example script 3: Sensor space simulation
clear, clc, close all
% Quick example for if you want to simulate directly in sensor space, as
% opposed to simulating in source space and mapping to EEG. 

% make pars structure for brain simulator
pars = struct ; 
 
% make custom atlas
pars.atlas = struct ; 
pars.atlas.atlas = 'custom' ; % use a custom atlas, not a prebuilt one
pars.atlas.regions = [1:19]' ; % use 19 electrodes

% make connectome - example based on timecourses
pars.connectome = rand(19).*(rand(19)<0.4).*(~eye(19)) ; % make a random connectome for example (with 40% probability of connection between regions)

% leadfield should be 1-to-1
pars.leadfield = eye(19) ; 
% Essentially, we are simulating 19 'regions of the brain', then mapping
% them directly to an electrode (i.e. region 1 is equivalent to electrode
% 1, etc). 

% only plot EEG
pars.options.plot_simulation = false ; 
pars.options.plot_EEG = true ; 

cd('example_run/example3_sensor_only')
out = brainSimulator(pars) ; 
save('out','-struct','out') % save output
cd('../..') % return to main directory

%% Example script 4: Simulate a visual evoked potential using no prebuilt options
clear, clc, close all
% This section simulates a visual evoked potential using ICA derived ROIs
% (not atlas based) and the Jansen-Rit model. 

% make pars structure for brain simulator
pars = struct ; 
 
% make custom atlas
pars.atlas = struct ; 
pars.atlas.atlas = 'custom' ; % use a custom atlas, not a prebuilt one
pars.atlas.regions = [1:3]' ; % use 4 ROIs

% make connectome - example based Fig 11 of Jansen and Rit (1995) Biol
% Cybern 73:357-366
pars.connectome = [0    100 0 ; 
                   4000 0   0 ; 
                   0    0   0 ];
% this one is just a toy model where ROI 1 (primary visual) strongly
% influences ROI 2 (e.g. extrastriate), whilst ROI 2 weakly influences ROI
% 1. ROI 3 (e.g. non-visual areas) is unconnected.
                
% dynamic model - Jansen Rit
addpath('example_run/example4_vep') % path to where custom ODE files are stored
pars.dynamic_model = @(t,x,model_pars,W) JRode(t,x,model_pars,W) ; 
% Jansen Rit ODE model. Custom ODE models need to be in the format
% fun(t,x,model_parameters,W), where t is time, x is a vector of state
% variables, model_parameters is pars.model_parameters, and W is your
% connectome. 
 
% model parameters
pars.model_parameters = struct ; 
pars.model_parameters.dim = 6 ; % JR model has 6 dimensions (equations) per ROI
pars.model_parameters.A = 3.25 ; 
pars.model_parameters.a = 100 ; 
pars.model_parameters.B = 22 ; 
pars.model_parameters.b = 50 ; 
pars.model_parameters.C = 125 ; % chosen so system is in non-oscillating regime
pars.model_parameters.v0 = 6 ; 
pars.model_parameters.e0 = 2.5 ; 
pars.model_parameters.r = 0.56 ;
pars.model_parameters.tdelay = 55/1000 ; % delay until stimulus reaches ROI 3
pars.model_parameters.p = @(t) JRstim_VEP(t,pars.model_parameters.tdelay) ;  
pars.model_parameters.std_noise = 0 ; % deterministic
N = 3 ; 
pars.model_parameters.observed_variables = @(x) x(N+1:2*N,:) - x(2*N+1:3*N,:); % in JR model, it is (x_ep - x_ip) that is the simulated EEG. 
 
% simulation options
pars.simulation_options = struct ; 
pars.simulation_options.t_init = 4 ; % note that unlike previous examples (which were in milliseconds), JR model is in seconds
pars.simulation_options.t_step = 0.1e-3 ; 
pars.simulation_options.t_simulate = 0.5 ;
 
% options
pars.options = struct ; 
pars.options.log_output = true ; 
pars.options.plot_eeg = true ; 
pars.options.plot_simulation = true ; 
 
% eeg electrodes
pars.eeg_electrodes = {'vep'} ; % this is the average over all occipital electrodes, but this averaging is accounted for in the leadfield 
 
% leadfield
pars.leadfield = sort(rand(1,N),'descend')  ; 
% this leadfield is chosen to simulate an electrode over primary visual
% cortex. Since we have set parameters such that ROI 1 is primary visual,
% ROI 2 is extrastriate, and ROI3 is non-visual areas, it makes sense that
% ROI 1 has the strongest leadfield, ROI 2 is next strongest, and ROI 3 is
% weakest, hence the descending sort. 

% simulate 
cd('example_run/example4_vep')
out = brainSimulator(pars) ; 
save('out','-struct','out') % save output
cd('../..') % return to main directory
 

