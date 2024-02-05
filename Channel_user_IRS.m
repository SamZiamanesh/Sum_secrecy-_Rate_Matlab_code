function hij_user_IRS = Channel_user_IRS(L,d1,d2 , Gij_Los)


L0 = 10^(-3);
dij = sqrt( (d1(1) -d2(1))^2 + (d1(2) -d2(2))^2 + (d1(3) -d2(3))^2 );
cij = 2;
bij = 0;

Gij_NLos = Gij_Los;

hij_user_IRS = sqrt(L0 * dij^(-cij)) * (sqrt(bij/(1+bij))*Gij_Los(1:L) + sqrt(1/(1+bij))*Gij_NLos(1:L) );


end
