function [Ds_est] = DTP2D_lookup_Dsa(theta0a,T)
data = [0.954433145	2.71235E-15
        0.507974605	2.14359E-14
        0.123300993	2.55E-13];
Ds0 = interp1(data(:,1),data(:,2),theta0a,'linear','extrap');
Ea = 30e3;
R = 8.314;
Ds_est = Ds0*exp(-Ea/R*(1/T - 1/298.15));
end
