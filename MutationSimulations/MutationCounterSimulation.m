%This uses the Graphical driver template to investigate the effect of 
%mutation rate on the average number of loci in the genome.

%model parameters
    %General Parameters
        N = 1000; %population size 
        num_loci = 1000;%number of loci in accessory genome
        mu_1 =  10^(-7);%mutation rate
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
    num_reps = 10; %number of replicate runs
    num_gens = 100;%number of generations to run simulation 

%initialize data structures
        avg_loci = nan(num_reps, num_gens); %avg # of loci of each indv. in each generation
        entropies = avg_loci;
        freq_mutator = avg_loci;
        
        per_gen_KO = zeros(2,num_gens);
        per_gen_GOF = zeros(2,num_gens);
        per_gen_Lethal = zeros(2,num_gens);
        
        per_rep_KO =zeros(2,num_reps);
        per_rep_GOF =zeros(2,num_reps);
        per_rep_Lethal =zeros(2,num_reps);

%for the figure
    close all;
    
    fontSize = 15;
   
    


    
    for j = 1:num_reps
        pop = [N,ones(1,num_loci),mu_1]; %starts all individuals with all functional loci
        rand_vect = rand(1,num_loci); %initializes with a random environment
        env = rand_vect<=frac_needed_scalar;
        for k = 1:num_gens

            output = mutate4(pop, rev_mut, num_essential); %mutate population
            %%%%%%%%%%% deals with extinction %%%%%%%%%%%%%%%%%%%
            if (size(pop,1) == 1) && (isnan(pop(:,2)))
                break
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Data Collection and Statistics
            %Top row is Observed Number and bottom row is Expected Number
            
            per_gen_KO(:,k) = output(:,1);
            per_gen_GOF(:,k) = output(:,2);
            per_gen_Lethal(:,k) = output(:,3);
                

        end
        per_rep_KO(:,j) = mean(per_gen_KO,2);
        per_rep_GOF(:,j) = mean(per_gen_GOF,2);
        per_rep_Lethal(:,j) = mean(per_gen_Lethal,2);
    end
    %[num_KO_mut, num_GOF_mut,num_lethal_mut]
    stat = [mean(per_rep_KO,2), mean(per_rep_GOF,2), mean(per_rep_Lethal,2)];
 %Graphics
                
