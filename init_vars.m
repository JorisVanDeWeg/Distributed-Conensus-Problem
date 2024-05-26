function Extras = init_vars(name)
    % Extras: Probailities of 'Special' Occurences
    % 1) Transmission Error Current Node; [, Probability of Error]
    % 2) Short term Dropout (per Iteration); [, Prob. of Dropout, Stopping Iteration #] 
    % 3) Long term Dropout; [, Iteration of Dropout, # of nodes in dropout ]
    % 4) Random Addition; [, Prob. of Addition, Stopping Iteration #]
    % 5) One time Addition; [, Iteration of Addition, # of Additions]
    % Ideal if it all scenarios are set to false
    
    % Hardcoded cases
    if name == "ADMM"
        Extras(1, 1:2) = [false, 0.4];
        Extras(2, :) = [false, 1, 1.5e5]; 
        Extras(3, :) = [false, 1.5e5, 20]; 
        Extras(4, :) = [false, 1.3333e-4 , 1.5e5];
        Extras(5, :) = [false, 1.5e5, 20]; 
    end
    if name == "PDMM"
        Extras(1, 1:2) = [false, 0.4];
        Extras(2, :) = [false, 1, 3e4]; 
        Extras(3, :) = [false, 1e4, 20];
        Extras(4, :) = [false, 2e-3, 1e4]; 
        Extras(5, :) = [false, 1e4, 20]; 
    end
    if name == "RG"
        Extras(1, 1:2) = [false, 0.4];
        Extras(2, :) = [false, 1, 1*10^5]; 
        Extras(3, :) = [false, 1*10^5, 20];
        Extras(4, :) = [false, 2*10^-4, 1*10^5]; 
        Extras(5, :) = [false, 1*10^5, 20]; 
    end
    if name == "RGRW"
        Extras(1, 1:2) = [false, 0.4];
        Extras(2, :) = [false, 1, 1*10^5]; 
        Extras(3, :) = [false, 5*10^4, 20];
        Extras(4, :) = [false, 4*10^-4, 5*10^4]; 
        Extras(5, :) = [false, 5*10^4, 20]; 
    end
end
