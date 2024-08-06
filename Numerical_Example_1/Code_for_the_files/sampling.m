function [V0] = sampling(p0,State_Space,N_tot)
%p0 initial probability distribution of the chemical system
%State_Space of the hidden species
Start=randsmpl(p0, 1, N_tot);
Shape=size(State_Space);
if length(unique(Start)) == 1
Start=unique(Start);
Start=State_Space(Start,:);
V0=repmat(Start', [1, N_tot]);
else
    V0=zeros(Shape(2),N_tot);
    for i=1:N_tot
        V0(:,i)=State_Space(Start(i),:);
    end
end

end

