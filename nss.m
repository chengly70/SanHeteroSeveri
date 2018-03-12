function ss = nss(v)
ss = (14./(1+exp(-(v-40)./9)))./(14./(1+exp(-(v-40)./9)) + exp(-v./45));