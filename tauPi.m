function tau=tauPi(v)       % time constant for inactivation gating variable Pi in Kr current
tau=1./(100*exp(-v./54.645)+656*exp(v./106.157));
end