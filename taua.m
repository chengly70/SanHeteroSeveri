function tau=taua(v,ach)    % time constant for activation gating variable a in Ach-activated K current
tau=1./(aa(ach)+ba(v));
end