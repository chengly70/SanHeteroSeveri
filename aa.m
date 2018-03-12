function a=aa(ach)      % for activation gating variable a in Ach-activated K current
a=(3.5988-0.0256)./(1+(0.0000012155./(ach.^1.6951)))+0.0256;
end