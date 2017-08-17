% this is the new one - edited 8 May 2013
% cs_adj.m
% This script computes the adjoint solutions to the TMI steady-advection
% problem for the 8 core site locations (9 core sites total). Saves every
% 100 years.

%%%%%% driver for TMI transient tracer simulations.
%%%
%%% G. Jake Gebbie, ggebbie@whoi.edu, WHOI, 27 July 2011.
%%% Based off: G. Gebbie & P. Huybers, "The mean age of ocean waters
%%% inferred from radiocarbon observations: sensitivity to surface
%%% sources and accounting for mixing histories", submitted, JPO.
%%% TMI = Total Matrix Intercomparison, Gebbie & Huybers, JPO, 2010.

clear

%% Load the TMI tendency matrix, At,
%  which satisfies dc/dt = At*c, where c is a global tracer distribution.
% load the feb file to get the variable Nfield
load ../../tmi/At_4deg_27july2011.mat
Nfield = length(it);

% Load inmixlyr, which specifies ML locations; these are the only that must
% be saved.
load ../../tmi/inmixlyr

%% List the years for which the global tracer is output.
years = [0:1:100 110:10:1000 1025:25:2000 2100:100:5000];
yearsmid = (years(1:end-1)+years(2:end))./2;
NY = length(years);

%% Loop through all core locations
path = '../../data/benthic/';
%files = dir([path '*.mat']);

files = {
    'GeoB1711.mat'
    'GeoB9526_4.mat'
    'M35003_4.mat'
    'MD07_3076.mat'
    'MD98_2165.mat'
    'MD99_2334K.mat'
    'NIOP_905.mat'
    'TR163_31b.mat'
        }

lf = length(files);

for cind = 1:lf
%    eval(['load ' path files(cind).name ';']);
%    name = strrep(files(cind).name,'.mat','');
eval(['load ' path char(files(cind)) ';']);
name = strrep(char(files(cind)),'.mat','');
    eval(['latind=' name '.tmi_latind' ';'])
    eval(['lonind=' name '.tmi_lonind' ';'])
    eval(['depthind=' name '.tmi_depthind' ';'])
    
    expname = name;
    
    outputfile = ['cadjs/',expname,'adj'];
    disp(outputfile)
    % Define the initial conditions.
    c0 = (jt==latind & it==lonind & kt==depthind);
    
    % Set the time tendency at the data location to zero
    Atc = At';
    Atc(c0,:)=0;
    %Atc(c0,c0) = 1;
    
    % Pre-allocate arrays.
    C = nan(NY,Nfield);
    T = nan(NY,1);
    
    % set the initial conditions.
    C(1,:) = c0;
    
    % options for the ODE solver.
    options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:Nfield,'Jacobian',Atc);
    
    %for nn = 1:NY
    for nn = 2:2:(NY-1)
        disp(nn)
        ind = nn-1:nn+1;
        yrs = years(ind);
        Cin  = squeeze(C(ind(1),:));
        tic
        [Ttmp,Cout] = ode15s(@(t,x) Atc*x,yrs,Cin,options);
        toc
        C(ind,:) = squeeze(Cout);
        T(ind) = squeeze(Ttmp);
    end
    
    C = C(:,inmixlyr);
    eval(['save ',outputfile,' C T']) % save intermediate values if necessary
    disp(['saved ' name])
end

