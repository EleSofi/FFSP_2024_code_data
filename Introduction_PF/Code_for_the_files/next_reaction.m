function [I] = next_reaction(W)
%Function which achieves the realizations of the index of the next reaction
%random variable with probability mass function P(I_next=i)=(w(i)/w0), with
%i=1,2,3..M (assuming we have a system with M reactions).
%W is a vector whose components w(i) are the propensity functions for each
%reaction of the system.
%w0 is the sum of all the propensity function given the current state of
%the system.
u=rand;
w0=sum(W);
s=W(1)/w0;
count=2;
I=1;
while (u > s) && (count <= length(W))
    I=I+1;
    s=s+(W(count)/w0);
    count=count+1;
end
end

