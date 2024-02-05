function [  G_W_PA_PB]  = AO_convergence(L , Parameters ,HA , HB , HCA , HCB ,W, HA_dagger , HB_dagger ,HCA_dagger , HCB_dagger )


sigma_a = 10^-10;
sigma_b = 10^-10;
sigma_c = 10^-10;

% initialize

[PA_star , PB_star]= Optimal_PA_PB(Parameters , HA , HB , HCB , HCA , W);

[ w_star_one_it , W_iteration] = Optimal_phase_IRS(L ,Parameters , HA , HB , HCB , HCA , W , PA_star ,PB_star);


% Main Loop

for j=1:100
    for i=1:3

    New_Optimal_W = W_iteration;

    [ w_star , W_iteration] = Optimal_phase_IRS(L ,Parameters , HA , HB , HCB , HCA , New_Optimal_W , PA_star ,PB_star);


    end

    [PA_star , PB_star]= Optimal_PA_PB(Parameters , HA , HB , HCB , HCA , New_Optimal_W);

    F_W_PA_PB = sigma_c + PA_star*real(trace(HCA*W_iteration)) + (PB_star)*real(trace(HCB*W_iteration));

    G_W_PA_PB(j) = log(sigma_b + real(PA_star*trace(HB*W_iteration))) + log(sigma_a + real((PB_star)*trace(HA*W_iteration)))...
            -log(F_W_PA_PB);




end



end