function fixed_X = fixed_states(x0)

fixed_X = ones(size(x0))*NaN;

% % SCFA - from 100 in co1 to 40 in co5
% power = (40/100)^(1/4);
% fixed_X(16:20) = [100 100*power 100*power^2 100*power^3 100*power^4];
% 
% % Carbs
% fixed_X(21:25) = [0.025/(60*24) 0 0 0 0];
% 
% % RS
% fixed_X(26:30) = [0.075/(60*24) 0 0 0 0];