function hrf = glover_hrf(t,p)

%	calculate hrf of the form
%	hrf = a*t.^n1.*exp(-t/t1) - b*t.^n2.*exp(-t/t2)
%	see Glover, NI 9:416 (1999).

a = p(1); b = p(2); t1 = p(3); t2 = p(4); n1 = p(5); n2 = p(6);

c1 = max(t.^n1.*exp(-t/t1));
c2 = max(t.^n2.*exp(-t/t2));
hrf = a*t.^n1.*exp(-t/t1)/c1 - b*t.^n2.*exp(-t/t2)/c2;
end
