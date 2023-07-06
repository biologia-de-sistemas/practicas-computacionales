function  [Pss_off, Bss_off]=StableFocalPoints_off(P_env ,gamma_B, kappa_P, delta_P)

Pss_off=(P_env*kappa_P)/(delta_P + delta_P*gamma_B);
Bss_off=1;

end

