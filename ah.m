function a=ah(v)     % Alpha for inactivation gating variable h in Na current
a=20*exp(-0.125*(v+75));
end