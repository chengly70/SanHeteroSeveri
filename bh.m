function b =bh(v) % Beta for inactivation gating variable h in Na current
b=2000./(320*exp((-0.1)*(v+75))+1);
end