function new_pop = compute_fitness3(pop,s_gains,env,s_baggage)
%"classes" is a cell array derived from pop using the unique_genotypes vector
%"weights" is the fitness weight of each function (fitness gain when function is
%   needed), should be a vector of length=num_functions
%"env" specifies the functions needed in the current environment, the
    %number of columns in "env" is the number of functions, the number of rows
    %is 1 and the number of layers is the number of replicates
%"costs" is a vector specifying the cost of carrying unneeded functions 
    %(fitness loss when function is present but not needed
    

     %extract information from the matrix 
     
    num_in_class = pop(:,1); %first column of the matrix is the number of individuals with that genotype 
                                    %(due to output of unique_genotypes function)
    genotypes = pop(:,2:end-1); %remaining columns of the matrix are the genotypes
    %mutation_rate = pop(:,end);
    num_genotypes = length(num_in_class); %total number of unique genotypes
    
    %weight_mat = repmat(weights,num_genotypes,1); % repeats vector "weights" across each row for easy matrix operations
    %cost_mat = repmat(costs,num_genotypes,1); %same as prev line but for "costs" vector
    
    
    env_mat = repmat(env,num_genotypes,1);
    
    gets_contribution = (genotypes.*env_mat==1); %determines which loci result in a fitness contribution
    gets_cost = (genotypes.*(~env_mat)==1); %determines which loci result in a fitness cost
    
    %fitness_effects = (gets_contribution.*s_gains)-(gets_cost.*s_baggage); %determines the combined fitness effects at each site in each genotype
    
    total_fitnesses = (1+s_gains).^(sum(gets_contribution,2)).*(1-s_baggage).^(sum(gets_cost,2));
    
    
    new_pop = [total_fitnesses,pop]; %outputs the number of each unique genotype and its fitness
    
    new_pop(isnan(pop(:,2)),1) = 0;