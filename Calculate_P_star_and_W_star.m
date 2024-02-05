function [PA_optimal , PB_optimal , w_optimal ]=Calculate_P_star_and_W_star(L , Parameters ,HA , HB , HCA , HCB ,W)


sigma_a = 10^-10;
sigma_b = 10^-10;
sigma_c = 10^-10;

% initialize

[PA_star , PB_star]= Optimal_PA_PB(Parameters , HA , HB , HCB , HCA , W);

[ w_star_one_it , W_iteration] = Optimal_phase_IRS(L ,Parameters , HA , HB , HCB , HCA , W , PA_star ,PB_star);


% Main Loop

for j=1:7
    for i=1:2

    New_Optimal_W = W_iteration;

    [ w_star , W_iteration] = Optimal_phase_IRS(L ,Parameters , HA , HB , HCB , HCA , New_Optimal_W , PA_star ,PB_star);


    end

    [PA_star , PB_star]= Optimal_PA_PB(Parameters , HA , HB , HCB , HCA , New_Optimal_W);

    F_W_PA_PB = sigma_c + PA_star*real(trace(HCA*W_iteration)) + (PB_star)*real(trace(HCB*W_iteration));




end



% Optimal Values

for i=1:L+1

    w_optimal_L(i) = exp(1i*w_star(i)/w_star(L+1));

end

PA_optimal = PA_star;
PB_optimal = PB_star;
w_optimal = w_optimal_L;

end

