function [PA_optimal , PB_star]= Optimal_PA_PB(Parameters , HA , HB , HCB , HCA , W)

syms PA;

sigam_a = Parameters(1);
sigam_b = Parameters(2);
sigam_c = Parameters(3);
P_max   = Parameters(4);
PA_min   = Parameters(5);
PB_min   = Parameters(6);

F_W_PA_PB = sigam_c + PA*real(trace(HCA*W)) + (P_max-PA)*real(trace(HCB*W));

G_W_PA_PB = log(sigam_b + real(PA*trace(HB*W))) + log(sigam_a + real((P_max-PA)*trace(HA*W)))...
            -log(F_W_PA_PB);

diffrential  = diff(G_W_PA_PB);

PA_diff = double(solve(diffrential == 0));

for i=1:length(PA_diff)

   if PA_diff(i)<P_max && PA_diff(i) >PA_min

      PA_diff_feasible = PA_diff(i);

   end
end

PB_diff = P_max - PA_diff_feasible;


% caluculate objective function for and find Maximum

F_W_PA_PB_min = sigam_c + PA_min*real(trace(HCA*W)) + (PB_min)*real(trace(HCB*W));

G_W_PA_PB_min = log(sigam_b + real(PA_min*trace(HB*W))) + log(sigam_a + real((PB_min)*trace(HA*W)))...
               -log(F_W_PA_PB_min);


F_W_PA_PB_min_max_A = sigam_c + PA_min*real(trace(HCA*W)) + (P_max-PA_min)*real(trace(HCB*W));


G_W_PA_PB_min_max_A  = log(sigam_b + real(PA_min*trace(HB*W))) + log(sigam_a + real((P_max-PA_min)*trace(HA*W)))...
               -log(F_W_PA_PB_min_max_A);


F_W_PA_PB_min_max_B = sigam_c + PA_min*real(trace(HCA*W)) + (P_max-PA_min)*real(trace(HCB*W));


G_W_PA_PB_min_max_B  = log(sigam_b + real(PB_min*trace(HB*W))) + log(sigam_a + real((P_max-PB_min)*trace(HA*W)))...
               -log(F_W_PA_PB_min_max_B);


F_W_PA_PB_diff  = sigam_c + PA_diff_feasible*real(trace(HCA*W)) + (PB_diff)*real(trace(HCB*W));

G_W_PA_PB_diff  = log(sigam_b + real(PA_diff_feasible*trace(HB*W))) + log(sigam_a + real((P_max-PB_diff)*trace(HA*W)))...
               -log(F_W_PA_PB_diff);


maximum = max([G_W_PA_PB_min , G_W_PA_PB_min_max_B , G_W_PA_PB_min_max_A , G_W_PA_PB_diff]);

if maximum == G_W_PA_PB_diff

    PA_optimal = PA_diff_feasible;
    PB_star    = PB_diff;

elseif maximum == G_W_PA_PB_min_max_A

    PA_optimal = PA_min;
    PB_star    =  P_max-PA_min;

elseif  maximum == G_W_PA_PB_min_max_B

    PA_optimal = PB_min;
    PB_star    = P_max - PB_min;

else
      PA_optimal = PA_min;
     PB_star    = PB_min;

end