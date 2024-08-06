function [dydt] = Evolution_Probability_FFSP_2(t,y,C,Hidden_Species,Y,S_bis,mu,MAK,W_star,New,index_state,network)


A=Evolution_Matrix_FFSP_2(C,Hidden_Species,Y,S_bis,mu,MAK,W_star,New,index_state,network);
dydt=A*y;

end