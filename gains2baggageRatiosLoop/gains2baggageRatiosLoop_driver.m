   %General Parameters
        N = 1000; %population size 
        num_loci = 1000;%number of loci in accessory genome
        m = [100]; %m = [1, 10,50, 100, 1000, 10000];
        mu_1 =  10^(-7);%mutation rate
        mu_2 = 100*mu_1;
        rev_mut = .001;
        num_essential = round(.25*num_loci);

    %scalar values of some things
        baggage2gains_ratio = [.005,.05;.01, .1;.05,.5; .1, 1]; 
        frac_needed_scalar = .1; %fraction_needed 
        alpha= 1.96;
    %vector values of the same things things
        % weights = .1 * ones(1, num_loci);%vector of length "num_loci" denoting fitness contribution 
        %             %of each functional locus when needed
        % costs = .001 * ones(1, num_loci);%vector denoting fitness cost of functional loci when not needed
        % fraction_needed = .1 * ones(1, num_loci);%vector denoting fraction of environments each locus is needed
    
        
    %simulation parameters
        num_reps = 10; %number of replicate runs
    %initialize data structures
        is_fixed = zeros(length(baggage2gains_ratio(:,1)), num_reps);
        stop_time = zeros(length(baggage2gains_ratio(:,1)), num_reps);
        m_array_index = 0;
       



    
%for the figure
    close all;
   
    
    

for i = 1:length(baggage2gains_ratio(:,1)) 
    
    change_rate = m;%number of gens btwn environmental changes
    s_baggage = baggage2gains_ratio(i,1) ;
    s_gains = baggage2gains_ratio(i,2);
    
    for j = 1:num_reps
        pop = [N-1,ones(1,num_loci),mu_1;1,ones(1,num_loci),mu_2]; %starts all individuals with all functional loci
        rand_vect = rand(1,num_loci); %initializes with a random environment
        env = rand_vect<=frac_needed_scalar;
        k=1;
        counter = round((rand)*(change_rate - 1));
        while k == 1

            counter = counter + 1;
            
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
            
            
            
            mu_poly = unique(pop(:, end));
            if length(mu_poly) == 1
                k = 0;
                if mu_poly == mu_1
                    is_fixed(i,j) = 0;
                    stop_time(i,j)= counter;
                elseif mu_poly == mu_2
                    is_fixed(i,j) = 1;
                    stop_time(i,j)= counter;

                else
                     error('Value other than 0 or 1 are not permited')
                end
            end

        
               
        end


    end
end
    


pfix_r1 = sum(is_fixed(1,:),2)/num_reps;
pfix_r2 = sum(is_fixed(2,:),2)/num_reps;
pfix_r3 = sum(is_fixed(3,:),2)/num_reps;
pfix_r4 = sum(is_fixed(4,:),2)/num_reps;
phat =[pfix_r1;pfix_r2; pfix_r3; pfix_r4]% pfix_cr10; pfix_cr1000; pfix_cr10000
se = sqrt(((1-phat).*(phat))./N);
ci_plus = phat + alpha*se;
ci_minus = phat - alpha*se;
ci = [ci_minus, ci_plus]
baggage2gains_ratio
figure
scatter([1:length(baggage2gains_ratio(:,1))],phat)
errorbar([1:length(baggage2gains_ratio(:,1))],phat,ci_minus, ci_plus,'.')
a = [ 'P_{fix} for m =',num2str(m),''];
title(a, 'FontSize', 10);
xlabel('Ratio', 'FontSize', fontSize);
ylabel('P_{fix}', 'FontSize', fontSize);