function a=an(v)%Alpha for variable n
a=0.01*(v+50)./(1-exp(-(v+50)./10));
end