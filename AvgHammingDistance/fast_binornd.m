function rv = fast_binornd(N,p)

%Computes a binomial random variable using much faster draws from Poisson
%and Normal distributions where these approximations are appropriate. N must be a constant
%scalar integer value (convenient for Wright-Fisher calculations) p can be any 
%multidimensional numerical array of values between 0 and 1



large_p = p>0.5; %find values where p is greater than 0.5
p(large_p) = 1-p(large_p); %flop large values so that they are less than 0.5

rv = zeros(size(p)); %draws a binomial RV for each value of p in the matrix p

if N < 25 %if N is small than the full calculation can be done without much computation

    rv = binornd(N,p);
    rv(large_p) = N-rv(large_p);
    
else %if N is large than use approximations of the binomial
    

    symmetric = ((N*p)>5); %if N*p is greater than 5 than binomial dist is approximately symmetric, can use gaussian

    rv(~symmetric) = poissrnd(N*p(~symmetric));

    rv(symmetric) = round(normrnd(N.*p(symmetric),sqrt(N.*p(symmetric).*(1-p(symmetric)))));
                        %rounded to produce integer output from normal dist
    neg_vals = find(rv(symmetric)<0);

    %deals with possibility of negative values from norm by redrawing until
    %all values are positive (should be a very rare occurrence for N*p>5
    %but could still happen
    
    counter = 1;

        while  ~isempty(neg_vals)
            rv(neg_vals) = round(normrnd(N.*p(neg_vals),sqrt(N.*p(neg_vals).*(1-p(neg_vals)))));
            neg_vals = find(rv(symmetric)<0);
            counter = counter+1;
            if counter>100
                error('lots of negatives')
            end
        end
        
    
    
    %undo flipping to make p <0.5
    rv(large_p) = N-rv(large_p);
end


    
    
    
    
    


