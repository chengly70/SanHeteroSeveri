function ss=ass(v,ach)      % steady state for activation gating variable a in Ach-activated K current
ss=aa(ach)./(aa(ach)+ba(v));
end