function tau=taur(v)    % time constant for activation gating variable q in transient outward K current
part1=exp(0.09*(v+30.61));
part2=exp(-0.12*(v+23.84));
tau=0.001*0.66*1.4*(15.59./(1.037*part1+0.369*part2)+2.98);
end