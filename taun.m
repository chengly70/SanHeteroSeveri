function tau=taun(v)        % time constant for activation gating variable n in Ks current
part1=28./(1+exp(-(v-40)./3));
part2=exp(-(v-5)./25);
scale=1;
tau=scale*1./(part1+part2);
end