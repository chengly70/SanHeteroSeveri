function a=adL(v) %Alpha value for variable h
a=(-0.02839*(v+41.8))./(exp(-(v+41.8)./2.5)-1)-(0.0849*(v+6.8))./(exp(-(v+6.8)./4.8)-1);
end