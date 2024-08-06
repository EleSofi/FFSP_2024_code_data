function [New,index_state] = indeces_FFSP(Hidden_Species,mu,S_u)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here
Shape_State=size(Hidden_Species);
index_state=zeros(Shape_State(1),mu);
New=zeros(Shape_State(1),Shape_State(2)*mu);
Columns=1:Shape_State(2);
for i=1:Shape_State(1)
for v=1:mu
        New(i,Columns)=Hidden_Species(i,:)-S_u(:,v)';
        [index]=ismember(Hidden_Species,New(i,Columns),'rows');
        if isempty(index(index == 1)) ~= 1
            index_state(i,v)=find(index == 1);
        end
        Columns=Columns+repmat(Shape_State(2),1,Shape_State(2));
end
Columns=1:Shape_State(2);
end
end

