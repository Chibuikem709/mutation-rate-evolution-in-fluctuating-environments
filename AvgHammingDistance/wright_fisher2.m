function pop = wright_fisher2(fitnesses)

N=sum(fitnesses(:,2)); %compute population size

genotypes = fitnesses(:,3:end); %pulls out matrix of genotypes
weighted_fitnesses = fitnesses(:,1).*fitnesses(:,2); %weight fitness values by multiplying by number of individuals with that fitness
rel_fitnesses = weighted_fitnesses/sum(weighted_fitnesses); %compute relative fitness values (p parameters for multinomial)

samples = fast_multinomial(N,rel_fitnesses); %perform multinomial sampling

extinct = samples==0; %determines which genotypes went extinct (were not sampled)
genotypes = genotypes(~extinct,:); %removes extinct genotypes from genotypes matrix
samples = samples(~extinct); %removes zero elements from samples

pop = [samples,genotypes];


    