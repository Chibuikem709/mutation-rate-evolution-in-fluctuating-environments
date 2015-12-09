function [outputs] = mutate4(pop, rev_mu, num_essential)


new_pop = [];
dims = size(pop);
num_KO_mut = zeros(2,dims(1));
num_GOF_mut = zeros(2,dims(1));
num_lethal_mut = zeros(2,dims(1));


for k = 1:dims(1)
    current_class = pop(k,:);
    genotype = current_class(2:dims(2)); %genotype of current class
    
    num_loci = dims(2)-2;
    num_ind = current_class(1);
    mu = current_class(dims(2)); %mutation rate of current_genotype
    num_func = num_ind*(sum(genotype(1:num_loci),2)); %total number of sites that can mutate
    num_KO = num_ind*(num_loci-(num_func/num_ind)); %total number of knocked out sites that can back mutate
    
    num_muts = fast_binornd(num_func,mu);  %determine total number of mutations
    num_rev_muts = fast_binornd(num_KO, mu*rev_mu);%determine total number of reverse mutations
   
    
    expanded_class = repmat(genotype,num_ind,1); %expanded form where each row is a genotype
    func_sites = find_functionals(genotype, num_ind, 1); %sites where mutations can occur
    KO_sites= find_functionals(genotype, num_ind, 0);
  
    mut_sites = randperm(num_func,num_muts); %sites where new mutations are assigned
    rev_sites = randperm(num_KO, num_rev_muts);
    
    expanded_class(func_sites(mut_sites))=0; %change mutated values to 0
    expanded_class(KO_sites(rev_sites))=1;%change back mutated values to 1
    
    is_lethal = (poissrnd(mu*num_essential, 1, num_ind)>0);
    expanded_class(is_lethal,:) = -1;
    
    %Observed Number
    num_KO_mut(1,k) = num_muts;
    num_GOF_mut(1,k) = num_rev_muts;
    num_lethal_mut(1,k) = sum(is_lethal);
    
    %Expected Number
    num_KO_mut(2,k) = num_func*mu;
    num_GOF_mut(2,k)= num_KO*mu*rev_mu;
    num_lethal_mut(2,k)= mu*num_essential;
    
    
    
    if k ==1
        new_pop =expanded_class;
    else
    old_dims = size(new_pop);
    temp_pop = zeros(old_dims(1)+num_ind,dims(2)-1); %initializes matrix to accomodate new rows
    temp_pop(1:old_dims(1),:)=new_pop; 
    temp_pop((old_dims(1)+1):end,:) = expanded_class; 
    new_pop = temp_pop;
    end
    
end

%Averages
num_KO_mut = mean(num_KO_mut,2);
num_GOF_mut = mean(num_GOF_mut,2);
num_lethal_mut = mean(num_lethal_mut,2);

[genotypes,~,indices] = unique(new_pop,'rows'); %finds row indices of unique rows pop
num_classes = max(indices); %determines the total number of unique rows
num_in_classes = hist(indices,num_classes)'; %determines how many rows are in each fitness class
genotypes(genotypes(:,1) == -1, :) = nan;

outputs = [num_KO_mut, num_GOF_mut,num_lethal_mut]; %assigns the classes matrix for the current 


end 
function indices = find_functionals(genotype, num_ind, func_or_KO)

    start_indices = ((0:length(genotype)-1)*num_ind)+1;

    func_starts = start_indices(genotype == func_or_KO);

    temp = repmat(func_starts,num_ind,1);
    temp2 = repmat((0:(num_ind-1))',1, length(func_starts));

    indices = temp + temp2;

    indices = reshape(indices, [1, num_ind*length(func_starts)]);

end
    
    
































