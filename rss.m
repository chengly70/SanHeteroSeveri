function ss=rss(v)  %  steady state for inactivation gating variable r in transcient outward K current
ss=1./(1+exp(-(v-19.3)./15));
end