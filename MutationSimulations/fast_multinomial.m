function rv = fast_multinomial(N,p)

%uses iterative calls to "fast_binornd" function to compute multinomial
%random variable. "N" is a scalar parameter representing the number of
%outcomes. "p" is a vector corresponding to the probability of each outcome
%and must sum to 1. The output "rv" is a column vector equal in length to "p" and
%with a sum equal to "N".

length_p = length(p); 
rv = zeros(length_p,1); %initializes output vector rv

%the spirit of the loop below is to split the multinomial sample into a
%series of binomial draws. Starting with p(1) it will perform a binomial
%sampling with probability of success of p(1) and all other outcomes
%being failure 1-p(1). Then it will move to p(2), reweight the remaining values of
%p such that they sum to 1, and perform another binomial. The process
%iterates until N successful multinomial draws are completed

if abs(sum(p)-1)>= 10^(-10)
    sprintf(num2str(abs(sum(p)-1)))
    error('input vector p must sum to 1')
    
end

counter =1; %keeps count of loop interation


while ((sum(rv)<N) && (counter<=length_p-1)) %keeps sampling until all N samples are generated or all binomials are computed
    remaining_N = N-sum(rv); %determines how many more successful draws are needed to reach N
    focal_p = p(counter); %moves to next p in loop interation
    remaining_p = sum(p((counter+1):end)); %determines total probability remaining
    effective_p = focal_p/(focal_p+remaining_p); %reweights probabilities to sum to 1
    rv(counter) = fast_binornd(remaining_N,effective_p); %stores number of successful outcomes for focal p
    counter = counter+1; %updates counter
end


rv(end) = N-sum(rv); %any remaining samples must correspond to last outcome
