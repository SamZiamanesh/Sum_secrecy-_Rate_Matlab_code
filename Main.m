clc;
clear;
close all;

L= linspace(10,40,4);

sigma_a = 10^-10;
sigma_b = 10^-10;
sigma_c = 10^-10;

P_max  = 10^(1.5);

PA_min = 1;
PB_min = 1;


PA_max = 0.2*P_max;
PB_max = 0.8*P_max;

Gij_Los = randn(50 ,1)+1j*randn(50,1); % line of sight component of channel


%% IRS



for i=1:length(L)


% create channels

Bob_IRS_channel = Channel_user_IRS(L(i),[40,0,0],[0,10,0] , Gij_Los);

Alice_IRS_channel = Channel_user_IRS(L(i),[-30,0,0],[0,10,0] , Gij_Los) ; 

charli_IRS_channel = Channel_user_IRS(L(i),[2,0,0],[0,10,0] , Gij_Los);


Alice_Bob_channel =  Channel_each_user([-30,0,0] , [40,0,0] , Gij_Los);

Alice_charli_channel = Channel_each_user([-30,0,0],[0,0,0] , Gij_Los);

Bob_charli_channel   =  Channel_each_user([40,0,0],[0,0,0] , Gij_Los);

% 

HA_dagger = [Bob_IRS_channel'*diag(Alice_IRS_channel) ,Alice_Bob_channel ] ;
HA = HA_dagger'*HA_dagger;


HB_dagger = [Alice_IRS_channel'*diag(Bob_IRS_channel) , Alice_Bob_channel ] ;
HB = HB_dagger'*HB_dagger;

HCA_dagger = [Alice_IRS_channel'*diag(charli_IRS_channel) ,Alice_charli_channel] ;
HCA = HCA_dagger'*HCA_dagger;

HCB_dagger = [Bob_IRS_channel'*diag(charli_IRS_channel) , Bob_charli_channel ] ;
HCB = HCB_dagger'*HCB_dagger;


W_init_dagger = 1*exp(1j*unifrnd(-1 ,1 , [1,L(i)+1]));
W = W_init_dagger'*W_init_dagger;


Parameters = [sigma_a , sigma_b , sigma_c , P_max , PA_min ,PB_min ];


[PA_optimal , PB_optimal , w_optimal ] = Calculate_P_star_and_W_star(L(i) , Parameters ,HA , HB , HCA , HCB ,W);

I_y_a_P_W_star = log(1+real(PA_optimal*abs(w_optimal*HA_dagger')^2)/sigma_a);
I_y_b_P_W_star = log(1+real(PB_optimal*abs(w_optimal*HB_dagger')^2)/sigma_b);
I_y_c_P_W_star = log(1+(real(PA_optimal*abs(w_optimal*HCA_dagger')^2)+real(PB_optimal*abs(w_optimal*HCB_dagger')^2))/sigma_c);


R_sum(i)=I_y_a_P_W_star + I_y_b_P_W_star-I_y_c_P_W_star ; 



I_y_a_P_star = log(1+real(PA_optimal*abs(W_init_dagger*HA_dagger')^2)/sigma_a);
I_y_b_P_star = log(1+real(PB_optimal*abs(W_init_dagger*HB_dagger')^2)/sigma_b);
I_y_c_P_star = log(1+(real(PA_optimal*abs(W_init_dagger*HCA_dagger')^2)+real(PB_optimal*abs(W_init_dagger*HCB_dagger')^2))/sigma_c);

R_sum_W_out(i)=I_y_a_P_star + I_y_b_P_star-I_y_c_P_star ; 



end

%% NO IRS


for i=1:length(L)


Alice_Bob_channel =  Channel_each_user([-30,0,0] , [40,0,0] , Gij_Los);

Alice_charli_channel = Channel_each_user([-30,0,0],[0,0,0] , Gij_Los);

Bob_charli_channel   =  Channel_each_user([40,0,0],[0,0,0] , Gij_Los);



HA_dagger = Alice_Bob_channel ;
HA = HA_dagger'*HA_dagger;


HB_dagger =  Alice_Bob_channel  ;
HB = HB_dagger'*HB_dagger;

HCA_dagger = Alice_charli_channel;
HCA = HCA_dagger'*HCA_dagger;

HCB_dagger = Bob_charli_channel;
HCB = HCB_dagger'*HCB_dagger;


W_init_dagger = 1*exp(1j*unifrnd(0 ,1 , [1,1]));
W = W_init_dagger'*W_init_dagger;


Parameters = [sigma_a , sigma_b , sigma_c , P_max , PA_min ,PB_min ];


[PA_star , PB_star]= Optimal_PA_PB(Parameters , HA , HB , HCB , HCA , W);

I_y_a_P_star = log(1+real(PA_star*abs(W_init_dagger*HA_dagger')^2)/sigma_a);
I_y_b_P_star = log(1+real(PB_star*abs(W_init_dagger*HB_dagger')^2)/sigma_b);
I_y_c_P_star = log(1+(real(PA_star*abs(W_init_dagger*HCA_dagger')^2)+real(PB_star*abs(W_init_dagger*HCB_dagger')^2))/sigma_c);


R_sum_no_IRS(i)=I_y_a_P_star + I_y_b_P_star-I_y_c_P_star ; 

end


%% AO algrothim convergence

for i=1:1

% create channels

Bob_IRS_channel = Channel_user_IRS(5,[40,0,0],[0,10,0] , Gij_Los);

Alice_IRS_channel = Channel_user_IRS(5,[-30,0,0],[0,10,0] , Gij_Los) ; 

charli_IRS_channel = Channel_user_IRS(5,[2,0,0],[0,10,0] , Gij_Los);


Alice_Bob_channel =  Channel_each_user([-30,0,0] , [40,0,0] , Gij_Los);

Alice_charli_channel = Channel_each_user([-30,0,0],[0,0,0] , Gij_Los);

Bob_charli_channel   =  Channel_each_user([40,0,0],[0,0,0] , Gij_Los);

% 

HA_dagger = [Bob_IRS_channel'*diag(Alice_IRS_channel) ,Alice_Bob_channel ] ;
HA = HA_dagger'*HA_dagger;


HB_dagger = [Alice_IRS_channel'*diag(Bob_IRS_channel) , Alice_Bob_channel ] ;
HB = HB_dagger'*HB_dagger;

HCA_dagger = [Alice_IRS_channel'*diag(charli_IRS_channel) ,Alice_charli_channel] ;
HCA = HCA_dagger'*HCA_dagger;

HCB_dagger = [Bob_IRS_channel'*diag(charli_IRS_channel) , Bob_charli_channel ] ;
HCB = HCB_dagger'*HCB_dagger;


W_init_dagger = 1*exp(1j*unifrnd(0 ,1 , [1,5+1]));
W = W_init_dagger'*W_init_dagger;


Parameters = [sigma_a , sigma_b , sigma_c , P_max , PA_min ,PB_min ];


 G_W_PA_PB= AO_convergence(5 , Parameters ,HA , HB , HCA , HCB ,W , HA_dagger , HB_dagger ,HCA_dagger , HCB_dagger );


end






%% W* and P* for different Pmax

diff_P_max = [10^0.5 , 10^1 , 10^2];
L_1= linspace(10,30,3);

for j=1:length(diff_P_max)

  for i=1:length(L_1)


%create channels

     Bob_IRS_channel = Channel_user_IRS(L_1(i),[40,0,0],[0,10,0] , Gij_Los);

     Alice_IRS_channel = Channel_user_IRS(L_1(i),[-30,0,0],[0,10,0] , Gij_Los) ; 

     charli_IRS_channel = Channel_user_IRS(L_1(i),[2,0,0],[0,10,0] , Gij_Los);

     Alice_Bob_channel =  Channel_each_user([-30,0,0] , [40,0,0] , Gij_Los);

     Alice_charli_channel = Channel_each_user([-30,0,0],[0,0,0] , Gij_Los);

     Bob_charli_channel   =  Channel_each_user([40,0,0],[0,0,0] , Gij_Los);



     HA_dagger = [Bob_IRS_channel'*diag(Alice_IRS_channel) ,Alice_Bob_channel ] ;
     HA = HA_dagger'*HA_dagger;


     HB_dagger = [Alice_IRS_channel'*diag(Bob_IRS_channel) , Alice_Bob_channel ] ;
     HB = HB_dagger'*HB_dagger;

     HCA_dagger = [Alice_IRS_channel'*diag(charli_IRS_channel) ,Alice_charli_channel] ;
     HCA = HCA_dagger'*HCA_dagger;

     HCB_dagger = [Bob_IRS_channel'*diag(charli_IRS_channel) , Bob_charli_channel ] ;
     HCB = HCB_dagger'*HCB_dagger;


     W_init_dagger = 1*exp(1j*unifrnd(-1 ,1 , [1,L_1(i)+1]));
     W = W_init_dagger'*W_init_dagger;


    diff_Parameters = [sigma_a , sigma_b , sigma_c , diff_P_max(j) , PA_min ,PB_min ];


    [PA_optimal , PB_optimal , w_optimal ] = Calculate_P_star_and_W_star(L_1(i) , diff_Parameters ,HA , HB , HCA , HCB ,W);

    I_y_a_P_W_star = log(1+real(PA_optimal*abs(w_optimal*HA_dagger')^2)/sigma_a);
    I_y_b_P_W_star = log(1+real(PB_optimal*abs(w_optimal*HB_dagger')^2)/sigma_b);
    I_y_c_P_W_star = log(1+(real(PA_optimal*abs(w_optimal*HCA_dagger')^2)+real(PB_optimal*abs(w_optimal*HCB_dagger')^2))/sigma_c);


    diff_R_sum(j,i)=I_y_a_P_W_star + I_y_b_P_W_star-I_y_c_P_W_star ; 



   I_y_a_P_star = log(1+real(PA_optimal*abs(W_init_dagger*HA_dagger')^2)/sigma_a);
   I_y_b_P_star = log(1+real(PB_optimal*abs(W_init_dagger*HB_dagger')^2)/sigma_b);
   I_y_c_P_star = log(1+(real(PA_optimal*abs(W_init_dagger*HCA_dagger')^2)+real(PB_optimal*abs(W_init_dagger*HCB_dagger')^2))/sigma_c);

   diff_R_sum_W_out(j,i)=I_y_a_P_star + I_y_b_P_star-I_y_c_P_star ; 



  end


  figure(3)

  subplot(2,2,j);
  plot(L_1, diff_R_sum(j,:) ,'-squre','Linewidth',2,"color","blue");
  hold on;
  plot(L_1, diff_R_sum_W_out(j,:)  ,'-o','Linewidth',2,"color","red");

  legend(["S:(W* ,P*)" , "S:(W,P*)"]);
  title("Sum Secrecy Rate for different P max(dBm)= ",10*log(diff_P_max(j))/log(10));
  xlabel("L");
  ylabel("Sum Secrecy Rate(bps/Hz) ");
  grid on;
  set(gca , "Fontsize",15);



end




%% plot Results

figure(1);

plot(diff(G_W_PA_PB) ,'Linewidth',2);
title("Error");
xlabel("Iteration");
ylabel("G(W,P_A,P_B)");
grid on;
set(gca , "Fontsize",15);

% figure(2);
% 
plot(G_W_PA_PB,'-squre','Linewidth',2);
title("AO algrothim convergence");
xlabel("Iteration");
ylabel("G(W,P_A,P_B)");
grid on;
set(gca , "Fontsize",15);
% 

figure(2);

load("R_sum_W_star_P_star_1.mat");

plot(L, R_sum ,'-squre','Linewidth',2,"Color","blue");
hold on;
plot(L, R_sum_W_out ,'-o','Linewidth',2,"Color","red");
hold on;
plot(L, R_sum_no_IRS ,'-*','Linewidth',2,"Color","green");


legend(["S:(W* ,P*)" , "S:(W,P*)" , "S:(NO-IRS,p*)"]);
xlabel("L");
ylabel("Sum Secrecy Rate(bps/Hz) ");
grid on;
set(gca , "Fontsize",15);

