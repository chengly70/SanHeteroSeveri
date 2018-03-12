function tau=taufT(v)       % time constant for inactivation gating variable fT in T type Ca current
tau=1./(16.67*exp(-(v+75)./83.3)+16.67*exp((v+75)./15.38));
end