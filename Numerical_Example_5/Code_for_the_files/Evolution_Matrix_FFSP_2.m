function [A] = Evolution_Matrix_FFSP_2(C,Hidden_Species,Y,S_bis,mu,MAK,W_star,New,index_state,network)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
Shape_State_Space=size(Hidden_Species);
Rows_A=Shape_State_Space(1);
A=zeros(Rows_A,Rows_A);
Columns=1:Shape_State_Space(2);
for i=1:Rows_A
    if mu > 0
    for v=1:mu
       
        if index_state(i,v) ~= 0
           if MAK == 1 && network == 1
           %[W] = propensity(New(i,v),Y,C,S_bis);
           W = propensity_1_bimolecular(New(i,Columns),Y,C,S_bis,v);
           A(i,index_state(i,v))=W;
           elseif MAK == 1 && network == 0
               W = propensity_1(New(i,Columns),Y,C,S_bis,v);
               A(i,index_state(i,v))=W;    
           else
               %W=zeros(1,length(W_star));
               %for k=1:length(W_star)
               %W(k)=W_star{k}(New(i,v),Y);
               %A(i,index_state(i,v))=W(v);
               %end
               A(i,index_state(i,v))=W_star{v}(New(i,Columns),Y);
           end
           %A(i,index_state(i,v))=W(v);
           %A(i,index_state(i,v))=W;
        
        end
        Columns=Columns+repmat(Shape_State_Space(2),1,Shape_State_Space(2));   
    end
    end
    Columns=1:Shape_State_Space(2);
    if MAK == 1 && network == 1 
        [W] = propensity_bimolecular(Hidden_Species(i,:),Y,C,S_bis);
    elseif MAK == 1 && network == 0
        [W] = propensity_G(Hidden_Species(i,:),Y,C,S_bis);
    else
        W=zeros(1,length(W_star));
        for k=1:length(W_star)
        W(k)=W_star{k}(Hidden_Species(i,:),Y);
        end
    end
    A(i,i)=-sum(W);
end

        
        
end
