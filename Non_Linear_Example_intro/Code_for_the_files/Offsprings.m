function [V_tilda,w] = Offsprings(V,w)
w_tilda=w/sum(w);
o=zeros(length(w),1);

u=rand(length(w)-1,1);

g=length(w);
h=length(w);
for i=1:length(w)-1
    if fraction(length(w)*w_tilda(i))+fraction(g-(length(w)*w_tilda(i))) < 1
        if u(i) < 1 - (fraction(length(w)*w_tilda(i))/fraction(g))
            o(i)=fix(length(w)*w_tilda(i));
        else
            o(i)=fix(length(w)*w_tilda(i))+(h-fix(g));
        end
    else 
        if u(i) < 1 - ((1-fraction(length(w)*w_tilda(i)))/(1-fraction(g)))
            o(i)=fix(length(w)*w_tilda(i))+1;
        else
            o(i)=fix(length(w)*w_tilda(i))+(h-fix(g));
        end
    end
    g=g-(length(w)*w_tilda(i));
    h=h-o(i);
end
o(end)=h;

V_tilda=zeros(size(V));
j=1;
for i=1:length(w)
    w(i)=1;
    for l=1:o(i)
        V_tilda(:,j)=V(:,i);
        j=j+1;
    end
end
        
end
            



