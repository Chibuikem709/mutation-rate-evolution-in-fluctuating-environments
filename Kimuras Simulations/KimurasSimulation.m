%model parameters
    %General Parameters
        N = 1000; %population size 
        num_loci = 2;%number of loci in accessory genome
        mu_1 =  0;%mutation rate
        mu_2 = 100*mu_1;
        rev_mut = .001;
        num_essential = round(.25*num_loci);

    %scalar values of some things
        s_baggage = .1;%costs
        s_gains = .1;%weights
        frac_needed_scalar = .1; %fraction_needed 
    %vector values of the same things things
        % weights = .1 * ones(1, num_loci);%vector of length "num_loci" denoting fitness contribution 
        %             %of each functional locus when needed
        % costs = .001 * ones(1, num_loci);%vector denoting fitness cost of functional loci when not needed
        % fraction_needed = .1 * ones(1, num_loci);%vector denoting fraction of environments each locus is needed
    


%simulation parameters
    num_reps = 10000; %number of replicate runs
    
%initialize data structures
        is_fixed = zeros(1, num_reps);
        p_fix = zeros(1, num_reps);
        s = zeros(1, num_reps);
        rel_fit = zeros(1, num_reps);
        g1 = [0,1];
        g2 = [0,0];
        env = [0,0];
        stats = cell(2,2);
%for the figure
    close all;
    figure(1);
    fontSize = 15;


for j = 1:num_reps
    pop = [N-1,g1,mu_1;1,g2,mu_2]; %starts all individuals with all functional loci
    rand_vect = rand(1,num_loci); %initializes with a random environment
    
    k=1;
   
    while k == 1
      
    
        pop = mutate3(pop, rev_mut, num_essential); %mutate population
        %%%%%%%%%%% deals with extinction %%%%%%%%%%%%%%%%%%%
        if (size(pop,1) == 1) && (isnan(pop(:,2)))
            break
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fitnesses = compute_fitness3(pop,s_gains,env,s_baggage); %computes fitnesses for each unique genotype
        pop = wright_fisher2(fitnesses); %performs wright_fisher sampling to update population to next generation
        
        focual_class = fitnesses(ismember(fitnesses(:,3:end-1), g2, 'rows'),:);
        other_class = fitnesses(ismember(fitnesses(:,3:end-1), g1, 'rows'),:);
        rel_fit(j)= focual_class(1)/other_class(1);
        s(j)= rel_fit(j)-1;
        p_fix(j)= (1-exp(-s(j)))/(1-exp(-s(j)*N));
        

        if size(pop,1) == 1
            k = 0;
            if sum(pop(2:end-1) == g2) == length(g2)
                is_fixed(j) = 1;
            else
                 is_fixed(j) = 0;
            end
        end
    end
            

end
fig = scatter(1:num_reps,is_fixed);
stats{1,1} = 'Predict P_fix';
stats{2,1} = 'Observed P_fix';
stats{1,2} = num2str(mean(p_fix));
stats{2,2} = num2str(sum(is_fixed,2)/num_reps);

