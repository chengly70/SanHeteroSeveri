function tau=tauq(v) % time constant for activation gating variable q in transient outward K current
part1=exp(-0.08*(v+44));
part2=exp(0.1*(v+45.93));
tau=0.001*0.6*(65.17./(0.57*part1+0.065*part2)+10.1);
end