function ss=yfss(v)    % steady state for gating variable y in funny current
ss=1./(1+exp((v+52.5)./9));
end