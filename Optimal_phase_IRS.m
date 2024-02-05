function [ W_star , optimal_W]= Optimal_phase_IRS(L,Parameters , HA , HB , HCB , HCA , W_hat , PA_star ,PB_star)

Size_W = [L+1,L+1];

sigam_a = Parameters(1);
sigam_b = Parameters(2);
sigam_c = Parameters(3);

cvx_begin

   variable W( Size_W );

   minimize (-( (1/log(2) )*(log(sigam_b + real(PA_star*(trace(HB*W)))))  + (log(sigam_a + real(PB_star*trace(HA*W))))...
             -log( sigam_c + abs(real(PA_star*trace(HCA*W_hat))) + abs(real(PB_star*trace(HCB*W_hat)))) - ...
              (sigam_c + abs(real(PA_star*trace(HCA*(W-W_hat)))) + abs(real(PB_star*trace(HCB*(W-W_hat)))) - sigam_c)/(sigam_c + abs(real(PA_star*trace(HCA*(W_hat)))) +abs(real(PB_star*trace(HCB*(W_hat))))) ))  ;

    subject to 

        W >= 1;
        diag(W) == 1;
        norm(W-W_hat) <=30;


cvx_end


optimal_W = W;

W_star = generate_rank_one_vector(W);



end

