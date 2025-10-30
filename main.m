%% Parameter Settings for earthquake location inversion using a 2D uniform model
close all;clear;
tic;
% study area
x_interval=[0,38.3];z_interval=[0,12.1];
dx=0.1; dz=0.1; % to plot

dimparameter=3;

% receivers and hypocenters location and origin time
receiver=[
15.3716	0.01
20.8944	0.01
26.4712	0.01
15.3716	5.5328
26.4712	5.5328
15.3716	11.0556
20.8944	11.0556
26.4712	11.0556];

source=[20.8943943943944	5.53284966525435];

origintime_true=20;

source_indx=round((source-[x_interval(1),z_interval(1)])./[dx,dz]+1);

% velocity
v0=3;

% observational data
timeobs=origintime_true+sqrt( ( source(1)-receiver(:,1) ).^2+( source(2)-receiver(:,2)).^2 )/v0;  %台站到时

% grid
x=x_interval(1):dx:x_interval(2);z=z_interval(1):dz:z_interval(2);
[xx,zz]=ndgrid(x,z);
v=v0*ones(size(xx));

% initial sampling settings
close all;nsamples=60000;
theta_true=[source,origintime_true*ones(size(source,1),1)]';
% pool = parpool('local', min(size(source, 1),10));
%% Random Walk Metropolis(MH)
mh_time=0;initial_point=[20;0;0];
mh_burn_in=1000;
mh_burn_in=mh_burn_in*ones(1,size(source,1));
mh_results=cell(size(source,1),1);mh_potential=mh_results;
mh_accept=zeros(size(source,1),1);

tic;
% parfor i=1:size(source,1)
for i=1:size(source,1)  
    [mh_results{i},mh_potential{i},mh_accept(i)]=MH_Sampler_2D([0;0;0],[38.3;12.1;30], [01;01;01], ...
    initial_point, nsamples,timeobs(:,i),receiver,v0);
    mh_burn_in(1,i)=burned_estimate(mh_results{i},0.01);  
    mh_ess{i}=multiESS(mh_results{i}(:,mh_burn_in(1,i):end)');
end
toc;
mh_time=mh_time+toc;

for i=1:size(source,1)
[~, mh_maxIndex(i)] = max(mh_potential{i}(:, mh_burn_in(:,i):end), [], 2);
mh_map_res(:, i) = mh_results{i}(:, mh_burn_in(:,i)-1+mh_maxIndex(i));
mh_map_error(:, i) = mh_map_res(:,i)-theta_true(:,i);
mh_mean_res(:,i)=mean(mh_results{i}(:, mh_burn_in(:,i):end),2);
mh_mean_error(:,i)=mh_mean_res(:,i)-theta_true(:,i);
end

%ESS per second
mh_time_cell = repmat({mh_time}, size(mh_ess)); 
per_mh_ess= cellfun(@(e,t) e ./ t, mh_ess, mh_time_cell, 'UniformOutput', false);

% Gelman-Rubin test    
initial_point=[10;0;0];
% parfor i=1:size(source,1)
for i=1:size(source,1)
[mh_results2{i},~,~]=MH_Sampler_2D([0;0;0],[38.3;12.1;30], [1;1;1], initial_point, nsamples,...
    timeobs(:,i),receiver,v0);
     mh_burn_in2(1,i)=burned_estimate(mh_results2{i},0.01);
     mh_GR(:,i)=gelman_rubin(mh_results{i}, mh_results2{i},max(mh_burn_in(1,i),mh_burn_in2(1,i)));
end

%% HMC
hmc_time=0;initial_point=[20;0.0;0];
hmc_results=cell(size(source,1),1);hmc_potential=mh_results;
hmc_accept=zeros(size(source,1),1);
M=[1e0,0,0;0,1e0,0;0,0,1e0 ]/0.01^2;  %mass
hmc_burn_in=1000;
hmc_burn_in=hmc_burn_in*ones(1,size(source,1));
tic;
for i=1:size(source,1)
    [hmc_results{i},hmc_accept(i),hmc_potential{i},~,~,~]=Hamilmc(initial_point,0.1,20,3, ...
    nsamples,timeobs(:,i),receiver,M,v0); 
    hmc_burn_in(1,i)=burned_estimate(hmc_results{i},0.01);
    hmc_ess{i}=multiESS(hmc_results{i}(:,hmc_burn_in(1,i):end)');
end
toc;
hmc_time=hmc_time+toc;

for i=1:size(source,1)
    [~, hmc_minIndex(i)] = min(hmc_potential{i}(:, hmc_burn_in(:,i):end), [], 2);
    hmc_map_res(:, i) = hmc_results{i}(:, hmc_burn_in(:,i)-1+hmc_minIndex(i));
    hmc_map_error(:, i) = hmc_map_res(:,i)-theta_true(:,i);
    hmc_mean_res(:,i)=mean(hmc_results{i}(:, hmc_burn_in(:,i):end),2);
    hmc_mean_error(:,i)=hmc_mean_res(:,i)-theta_true(:,i);
end

% ESS per second
hmc_time_cell = repmat({hmc_time}, size(hmc_ess)); 
per_hmc_ess= cellfun(@(e,t) e ./ t, hmc_ess, hmc_time_cell, 'UniformOutput', false);

% Gelman-Rubin test    
initial_point=[10;00.0;0];
for i=1:size(source,1)
[hmc_results2{i},~,~,~,~,~]=Hamilmc(initial_point,0.1,20,3,nsamples,timeobs(:,i), ...
    receiver,M,v0);
     hmc_burn_in2(1,i)=burned_estimate(hmc_results2{i},0.01);
     hmc_GR(:,i)=gelman_rubin(hmc_results{i}, hmc_results2{i},max(hmc_burn_in(1,i),hmc_burn_in2(1,i)));
end
%% NUTS
clc;
% negative potential energy
logp = @(x) -Potential_energy(x,timeobs,receiver,v0);  
% the gradient of the negative potential energy
grad = @(x) -gradient_P(x,timeobs,receiver,v0);  

f = @(x) deal(logp(x), grad(x));


nuts_time=0;initial_point=[20;00.0;0];
theta0 = initial_point;
nuts_results=cell(size(source,1),1);nuts_potential=nuts_results;
nuts_accept=zeros(size(source,1),1);

nuts_burn_in=1000;
nuts_burn_in=nuts_burn_in*ones(1,size(source,1));
tic;
for i=1:size(source,1)
    [nuts_results{i}, nuts_potential{i}, epsilon] = runNUTS(f, theta0, ...
    'delta', 0.8, 'n_warmup', nuts_burn_in, 'n_mcmc', nsamples, 'max_tree_depth', 10,'n_updates', 1000);

    nuts_results{i}=[theta0,nuts_results{i}];

    nuts_potential{i}=-[0;nuts_potential{i}]';
    nuts_burn_in(1,i)=burned_estimate(nuts_results{i},0.01);
    nuts_ess{i}=multiESS(nuts_results{i}(:,nuts_burn_in(1,i):end)');
end
toc;
nuts_time=nuts_time+toc;

for i=1:size(source,1)
    [~, nuts_minIndex(i)] = min(nuts_potential{i}(:, nuts_burn_in(:,i):end), [], 2);
    nuts_map_res(:, i) = nuts_results{i}(:, nuts_burn_in(:,i)-1+nuts_minIndex(i));
    nuts_map_error(:, i) = nuts_map_res(:,i)-theta_true(:,i);
    nuts_mean_res(:,i)=mean(nuts_results{i}(:, nuts_burn_in(:,i):end),2);
    nuts_mean_error(:,i)=nuts_mean_res(:,i)-theta_true(:,i);
end

% ESS per second
nuts_time_cell = repmat({nuts_time}, size(nuts_ess)); 
per_nuts_ess= cellfun(@(e,t) e ./ t, nuts_ess, nuts_time_cell, 'UniformOutput', false);

% Gelman-Rubin test
initial_point=[10; 0.0;0];
nuts_burn_in2=nuts_burn_in;
for i=1:size(source,1)
    [nuts_results2{i},~ , ~] = runNUTS(f, initial_point, ...
    'delta', 0.8, 'n_warmup', nuts_burn_in, 'n_mcmc', nsamples, 'max_tree_depth', 10,'n_updates', 1000);

    nuts_results2{i}=[theta0,nuts_results2{i}];
    nuts_burn_in2(1,i)=burned_estimate(nuts_results2{i},0.01);
    nuts_GR(:,i)=gelman_rubin(nuts_results{i}, nuts_results2{i},max(nuts_burn_in(1,i),nuts_burn_in2(1,i)));
end

%% Zig-Zag

zz_time=0;initial_point=[20;0;0];
zz_results=cell(size(source,1),1);zz_potential=zz_results;zz_Time=zz_potential;
zz_flip=zeros(size(source,1),1);
zz_burn_in=1000;
zz_burn_in=zz_burn_in*ones(1,size(source,1));
tic;
for i=1:size(source,1)
[zz_Time{i},zz_results{i},zz_v,zz_potential{i},zz_flip(i),m,M,i0,taui0]=zig_zag_method(initial_point,[0;0;0],[38.3;12.1;30],dimparameter,nsamples, ...
    timeobs(:,i),0.01,[0.0;0.0;0.0],receiver,v0);
zz_burn_in(1,i)=burned_estimate(zz_results{i},0.01);
zz_ess{i}=multiESS(zz_results{i}(:,zz_burn_in(1,i):end)');
end
toc;
zz_time=zz_time+toc;

for i=1:size(source,1)
[~, zz_minIndex(i)] = min(zz_potential{i}(:, zz_burn_in(:,i):end), [], 2);
zz_map_res(:, i) = zz_results{i}(:, zz_burn_in(:,i)-1+zz_minIndex(i));
zz_map_error(:, i) = zz_map_res(:,i)-theta_true(:,i);
zz_mean_res(:,i)=mean(zz_results{i}(:, zz_burn_in(:,i):end),2);
zz_mean_error(:,i)=zz_mean_res(:,i)-theta_true(:,i);
end

% ESS per second
zz_time_cell = repmat({zz_time}, size(zz_ess)); 
per_zz_ess= cellfun(@(e,t) e ./ t, zz_ess, zz_time_cell, 'UniformOutput', false);

% Gelman-Rubin test
initial_point=[10;0;0];
for i=1:size(source,1)
[~,zz_results2{i},~,~,~]=zig_zag_method(initial_point,[0;0;0],[200;200;30],dimparameter,nsamples, ...
    timeobs(:,i),0.01,[0.0;0.0;0.0],receiver,v0);
     zz_burn_in2(1,i)=burned_estimate(zz_results2{i},0.01);
     zz_GR(:,i)=gelman_rubin(zz_results{i}, zz_results2{i},max(zz_burn_in(1,i),zz_burn_in2(1,i)));
end