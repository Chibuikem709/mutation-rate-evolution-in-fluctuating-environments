%General Parameters
N = 1000; %population size
num_loci = 1000;%number of loci in accessory genome
m = [100]; %m = [1, 10,50, 100, 1000, 10000];
mu_1 =  10^(-2);%mutation rate
rev_mut = .001;
num_essential = round(.25*num_loci);

%scalar values of some things
s_baggage = .01;%costs
s_gains = .1;%weights
frac_needed_scalar = .1; %fraction_needed
alpha= 1.96;
%vector values of the same things things
% weights = .1 * ones(1, num_loci);%vector of length "num_loci" denoting fitness contribution
%             %of each functional locus when needed
% costs = .001 * ones(1, num_loci);%vector denoting fitness cost of functional loci when not needed
% fraction_needed = .1 * ones(1, num_loci);%vector denoting fraction of environments each locus is needed


%simulation parameters
num_reps = 1; %number of replicate runs
time = 1000;
%initialize data structures
avg_hamming_dist = nan(num_reps,time);
m_array_index = 0;





%for the figure
close all;
figure(1);
fontSize = 15;






change_rate = m;%number of gens btwn environmental changes


for j = 1:num_reps
    pop = [N,ones(1,num_loci),mu_1]; %starts all individuals with all functional loci
    rand_vect = rand(1,num_loci); %initializes with a random environment
    env = rand_vect<=frac_needed_scalar;
    k=1;
    counter = 0;
    avg_hamming_dist(j,1) = sum(pop(:,2:end-1) ~= repmat(env, [length(pop(:,1)) 1]),2);
    point_in_time = 0;
    
    while k == 1
        
        counter = counter + 1;
        point_in_time = point_in_time + 1;
        if mod(counter,change_rate)==0 %determines if environment should change this generation
            rand_vect = rand(1,num_loci); %initializes with a random environment
            env = rand_vect<=frac_needed_scalar; %resets each environment
        end
        
        pop = mutate3(pop, rev_mut, num_essential); %mutate population
        %%%%%%%%%%% deals with extinction %%%%%%%%%%%%%%%%%%%
        if (size(pop,1) == 1) && (isnan(pop(:,2)))
            break
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fitnesses = compute_fitness3(pop,s_gains,env,s_baggage); %computes fitnesses for each unique genotype
        
        pop = wright_fisher2(fitnesses); %performs wright_fisher sampling to update population to next generation
        avg_hamming_dist(j,point_in_time) = mean(sum(pop(:,2:end-1) ~= repmat(env, [length(pop(:,1)) 1]),2));
        
        
        
        if point_in_time == time
            k = 0;
        end
        
        
        
    end
    
    
end
x2=1:length(avg_hamming_dist);
y2=avg_hamming_dist;
plot(x2,y2)
title([' Replicate =' num2str(num_reps) ], 'FontSize', fontSize);
ylabel('Avg. Hamming Distance', 'FontSize', fontSize);
xlabel('Time (Generation)', 'FontSize', fontSize);
