function [T] = Exponential(lambda)
%This is a function for gaining an exponential random variable realization. 
%The input is
%the lambda parameter of the Exponential Distribution. The output T is the
%realization of the Exponential random variable obtained applying the
%inverse exponential distribution function to a standard uniform random
%variable realization. 
u=rand;
T=((-log(1-u))/lambda);
end

