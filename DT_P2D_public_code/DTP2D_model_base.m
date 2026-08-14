function [output,paras,sep_elem] = DTP2D_model_base(f_vector,factors,multi_soc_range,T,type_anode)

    omegagen= f_vector*(2*pi());
    N = length(omegagen);

    A_coat = 12.6*2/10000;
    x_1 = 0.8781;
    x_0 = 0.0216;
    y_0 = 0.9319;
    y_1 = 0.3532;
    c_imp = zeros(N,length(multi_soc_range));
    a_imp = zeros(N,length(multi_soc_range));
    s_imp = zeros(N,length(multi_soc_range));
    fc_imp = zeros(N,length(multi_soc_range));

    for i = 1:length(multi_soc_range)
        soc = multi_soc_range(i)*0.01;
        x = x_0 + (x_1 - x_0)*soc;
        y = y_0 + (y_1 - y_0)*soc;
        L1 = factors(6,i)*1e-11;
        L2 = factors(6,i+length(multi_soc_range))*1e-11;
        RL1 = factors(7,i)*1e-4;
        RL2 = factors(7,i+length(multi_soc_range))*1e-4;
        R_itsca(i) = factors(1,i)*0.0255;
        R_itscc(i) = factors(1,i+length(multi_soc_range))*0.0255;
        R_filma = 0;                            R_filmc = 0;
        C_filma = 0;                        C_filmc = 0;
        if type_anode == 0
            Rpa =  (17.8)*1e-6;              Rpc = (10)*1e-6;
        elseif type_anode == 1
            Rpa =  (17.8)*1e-6;              Rpc = (10)*1e-6;
        elseif type_anode == 2
            Rpa =  (16.5)*1e-6;              Rpc = (10)*1e-6;
        end

        Dsa = factors(4,i)*DTP2D_lookup_Dsa(x,T);      Dsc = factors(4,i+length(multi_soc_range))*DTP2D_lookup_Dsc(y,T);
        Cdla =  factors(3,i)*0.2;             Cdlc = factors(3,i+length(multi_soc_range))*0.2;
        cta = 29626;                        ctc = 48786;
        ka = factors(2,i)*DTP2D_lookup_ka(x,T,cta);
        kc = factors(2,i+length(multi_soc_range))*DTP2D_lookup_kc(y,T,ctc);
        c_e_ref = 1000;
        if type_anode == 0
            La = 79e-6;                         Lc = 60.0e-6;
        elseif type_anode == 1
            La = 79e-6;                         Lc = 60.0e-6;
        elseif type_anode == 2
            La = 77e-6;                         Lc = 60.0e-6;
        end

        bruga = 1.44;                       brugc = 1.44;
        epsla =   0.237;                    epslc = 0.234;
        n_am1_vf =    0.977;                p_am1_vf = 0.9792;
        epssa =   (1-epsla)*n_am1_vf;       epssc = (1-epslc)*p_am1_vf;
        taua = epsla^(-bruga)*factors(8,i);              tauc = epslc^(-brugc)*factors(8,i+length(multi_soc_range));
        if type_anode == 0
            sigmaa = 5.6022e+03;                     sigmac = 24.2718;
        elseif type_anode == 1
            sigmaa = 5.6022e+03;                     sigmac = 24.2718;
        elseif type_anode == 2
            sigmaa = 4.2088e+03;                     sigmac = 24.2718;
        end

        cs0a = x*cta;                       cs0c = y*ctc;
        alphaa = 0.5;
        alphac = 0.5;

        c_e = 1120;
        Di0 = factors(1,2*length(multi_soc_range)+1)*DTP2D_lookup_De(c_e/1000,T);
        epsls = 0.5;
        Lsep = 12.5e-6;
        F = 96487;
        R = 8.314;
        iapp = 1;
        tplus = 0.363;
        nu=1;
        taus = 1.8;

        dlnfdlnc = (0.601-0.24*(c_e/(1000))^0.5+0.982*(1-0.0052*(T-298.15))*(c_e/(1000))^1.5)/(1-0.363)-1;
        kappa = factors(2,2*length(multi_soc_range)+1);
        aa =factors(5,i)*3*epssa/Rpa;
        ac =factors(5,i+length(multi_soc_range))*3*epssc/Rpc;
        A_tota = aa*A_coat*La;
        A_totc = ac*A_coat*Lc;

        i0a = F*ka*((c_e/c_e_ref)^alphaa)*((cta-cs0a)^alphaa)*cs0a^alphac;
        i0c = F*kc*((c_e/c_e_ref)^alphaa)*((ctc-cs0c)^alphaa)*cs0c^alphac;

        sigmaeffa =(epssa/taua)*sigmaa;
        sigmaeffc =(epssc/tauc)*sigmac;
        kappaeffa = (epsla/taua)*kappa;
        kappaeffc = (epslc/tauc)*kappa;
        kappaeffs = (epsls/taus)*kappa;
        Dieffa = (epsla/taua)*Di0;
        Dieffc = (epslc/tauc)*Di0;
        Dieffs = (epsls/taus)*Di0;

        dx = 0.0001;
        chg = 0.5;
        dUdcc = -((1/ctc)*(DTP2D_lookup_Uc(soc+dx,type_anode) - DTP2D_lookup_Uc(soc-dx,type_anode))/(2*dx));
        dUdca = ((1/cta)*(DTP2D_lookup_Ua(soc+dx,type_anode) - DTP2D_lookup_Ua(soc-dx,type_anode))/(2*dx));

        for k = 1:N
            sbar_a=1i*omegagen(k)*epsla*La^2/Dieffa;
            sbar_c=1i*omegagen(k)*epslc*Lc^2/Dieffc;
            sbar_s=1i*omegagen(k)*epsls*Lsep^2/Dieffs;
            s_a=1i*omegagen(k);
            s_c=1i*omegagen(k);

            Rcta=R*T/i0a/F/(alphaa+alphac);
            Rctc=R*T/i0c/F/(alphaa+alphac);
            Rdifa=-dUdca*Rpa/Dsa/F;
            Rdifc=-dUdcc*Rpc/Dsc/F;
            Ysa=(sqrt(s_a*Rpa^2/Dsa)-tanh(sqrt(s_a*Rpa^2/Dsa)))/tanh(sqrt(s_a*Rpa^2/Dsa));
            Ysc=(sqrt(s_c*Rpc^2/Dsc)-tanh(sqrt(s_c*Rpc^2/Dsc)))/tanh(sqrt(s_c*Rpc^2/Dsc));
            zetaa = s_a*Cdla + 1/((Rcta)+Rdifa/Ysa);
            zetac = s_c*Cdlc + 1/((Rctc)+Rdifc/Ysc);
            betaa = 1/F*(s_a * C_filma + 1/(1/zetaa+R_filma));
            betac = 1/F*(s_c * C_filmc + 1/(1/zetac+R_filmc));

            Theta1a = (1-tplus)*La/F/Dieffa*(2*R*T*(1-tplus)/F*(1/c_e)*(1+dlnfdlnc))/(1/La/aa/F/betaa);
            Theta1c = (1-tplus)*Lc/F/Dieffc*(2*R*T*(1-tplus)/F*(1/c_e)*(1+dlnfdlnc))/(1/Lc/ac/F/betac);
            Theta2a = La*(1/sigmaeffa+1/kappaeffa)/(1/La/aa/F/betaa);
            Theta2c = Lc*(1/sigmaeffc+1/kappaeffc)/(1/Lc/ac/F/betac);
            Theta3 = 2*kappaeffs*R*T*(1-tplus)/Dieffs/F*(1/c_e)*(1+dlnfdlnc);

            lambda1a = 1/2*(sbar_a + Theta1a + Theta2a + sqrt(sbar_a^2 + 2*Theta1a*sbar_a - 2*sbar_a*Theta2a + Theta1a^2 + 2*Theta1a*Theta2a + Theta2a^2));
            lambda1c = 1/2*(sbar_c + Theta1c + Theta2c + sqrt(sbar_c^2 + 2*Theta1c*sbar_c - 2*sbar_c*Theta2c + Theta1c^2 + 2*Theta1c*Theta2c + Theta2c^2));
            lambda2a = 1/2*(sbar_a + Theta1a + Theta2a - sqrt(sbar_a^2 + 2*Theta1a*sbar_a - 2*sbar_a*Theta2a + Theta1a^2 + 2*Theta1a*Theta2a + Theta2a^2));
            lambda2c = 1/2*(sbar_c + Theta1c + Theta2c - sqrt(sbar_c^2 + 2*Theta1c*sbar_c - 2*sbar_c*Theta2c + Theta1c^2 + 2*Theta1c*Theta2c + Theta2c^2));

            Lambda1a = -La^3*aa*iapp*(1-tplus)*betaa/sigmaeffa/nu/Dieffa/(lambda1a-lambda2a)*((1/sqrt(lambda2a)/sinh(sqrt(lambda2a))-1/sqrt(lambda1a)/sinh(sqrt(lambda1a)))...
                +sigmaeffa/kappaeffa*(1/sqrt(lambda2a)/tanh(sqrt(lambda2a))-1/sqrt(lambda1a)/tanh(sqrt(lambda1a))));
            Lambda1c = -Lc^3*ac*iapp*(1-tplus)*betac/sigmaeffc/nu/Dieffc/(lambda1c-lambda2c)*((1/sqrt(lambda2c)/sinh(sqrt(lambda2c))-1/sqrt(lambda1c)/sinh(sqrt(lambda1c)))...
                +sigmaeffc/kappaeffc*(1/sqrt(lambda2c)/tanh(sqrt(lambda2c))-1/sqrt(lambda1c)/tanh(sqrt(lambda1c))));
            Lambda3 = Lsep/(sqrt(sbar_s)*Dieffs*sinh(sqrt(sbar_s)));
            Lambda2a = La/Dieffa/(lambda1a-lambda2a)*((sbar_a-lambda2a+Theta1a)/(sqrt(lambda1a))/tanh(sqrt(lambda1a))-(sbar_a-lambda1a+Theta1a)/(sqrt(lambda2a))/tanh(sqrt(lambda2a))) + Lambda3*cosh(sqrt(sbar_s));
            Lambda2c = Lc/Dieffc/(lambda1c-lambda2c)*((sbar_c-lambda2c+Theta1c)/(sqrt(lambda1c))/tanh(sqrt(lambda1c))-(sbar_c-lambda1c+Theta1c)/(sqrt(lambda2c))/tanh(sqrt(lambda2c))) + Lambda3*cosh(sqrt(sbar_s));

            xi_star_a = (Lambda1a + (Lambda1c/Lambda2c)*Lambda3)/(Lambda2a - Lambda3^2/Lambda2c);
            xi_star_c = (Lambda1c + (Lambda1a/Lambda2a)*Lambda3)/(Lambda2c - Lambda3^2/Lambda2a);

            C_1_s = Theta3/(sqrt(sbar_s)*sinh(sqrt(sbar_s)))*(xi_star_c - xi_star_a*cosh(sqrt(sbar_s)));
            C_2_s = Theta3*xi_star_a/sqrt(sbar_s);
            cbar_s_x0 = C_1_s;
            cbar_s_x1 = Theta3/sqrt(sbar_s)*(xi_star_c/tanh(sqrt(sbar_s))-xi_star_a/sinh(sqrt(sbar_s)));
            phi2_s_x0 = Lsep*iapp/kappaeffs*((cbar_s_x0-cbar_s_x1)/iapp+1);
            phi2_s_x1 = 0;

            C_2_c = Theta1c*iapp/sqrt(lambda1c)/(lambda1c-lambda2c)*(sigmaeffc/kappaeffc-(sbar_c+Theta1c-lambda2c)*sigmaeffc*xi_star_c/iapp/Lc^2/ac/(1-tplus)/betac);
            C_1_c = -Theta1c*iapp/sqrt(lambda1c)/(lambda1c-lambda2c)/sinh(sqrt(lambda1c))-C_2_c/tanh(sqrt(lambda1c));
            C_4_c = -Theta1c*iapp/sqrt(lambda2c)/(lambda1c-lambda2c)*(sigmaeffc/kappaeffc-(sbar_c+Theta1c-lambda1c)*sigmaeffc*xi_star_c/iapp/Lc^2/ac/(1-tplus)/betac);
            C_3_c = Theta1c*iapp/sqrt(lambda2c)/(lambda1c-lambda2c)/sinh(sqrt(lambda2c))-C_4_c/tanh(sqrt(lambda2c));
            C_7_c = -Lc*iapp/sigmaeffc*(1-sbar_c*Lc^2*ac*F*betac/sigmaeffc/lambda1c/lambda2c);
            C_8_c = -Lc/sigmaeffc*((sbar_c-lambda1c)/Theta1c*C_1_c*(1-Lc^2*ac*F*betac/sigmaeffc/lambda1c)+(sbar_c-lambda2c)/Theta1c*C_3_c*(1-Lc^2*ac*F*betac/sigmaeffc/lambda2c));
            chi_lambda1_c = -Theta1c*iapp/sqrt(lambda1c)/(lambda1c-lambda2c)/tanh(sqrt(lambda1c))-C_2_c/sinh(sqrt(lambda1c));
            chi_lambda2_c =  Theta1c*iapp/sqrt(lambda2c)/(lambda1c-lambda2c)/tanh(sqrt(lambda2c))-C_4_c/sinh(sqrt(lambda2c));
            phi1_c_x1 = -Lc^3*ac*F*betac/sigmaeffc^2*((sbar_c-lambda1c)/Theta1c/lambda1c*chi_lambda1_c+(sbar_c-lambda2c)/Theta1c/lambda2c*chi_lambda2_c)+C_7_c+C_8_c;

            C_1_a = Theta1a*iapp/sqrt(lambda1a)/(lambda1a-lambda2a)*(1/tanh(sqrt(lambda1a))+sigmaeffa/kappaeffa/sinh(sqrt(lambda1a))-(sbar_a+Theta1a-lambda2a)*nu*sigmaeffa*xi_star_a/La^2/aa/iapp/betaa/(1-tplus)/sinh(sqrt(lambda1a)));
            C_2_a = -Theta1a*iapp/sqrt(lambda1a)/(lambda1a-lambda2a);
            C_3_a = -Theta1a*iapp/sqrt(lambda2a)/(lambda1a-lambda2a)*(1/tanh(sqrt(lambda2a))+sigmaeffa/kappaeffa/sinh(sqrt(lambda2a))-(sbar_a+Theta1a-lambda1a)*nu*sigmaeffa*xi_star_a/La^2/aa/iapp/betaa/(1-tplus)/sinh(sqrt(lambda2a)));
            C_4_a = Theta1a*iapp/sqrt(lambda2a)/(lambda1a-lambda2a);
            chi_lambda1_a =  Theta1a*iapp/sqrt(lambda1a)/(lambda1a-lambda2a)/sinh(sqrt(lambda1a))+Theta1a*iapp/sqrt(lambda1a)/(lambda1a-lambda2a)/tanh(sqrt(lambda1a))*(sigmaeffa/kappaeffa-(sbar_a+Theta1a-lambda2a)*nu*sigmaeffa*xi_star_a/La^2/aa/iapp/(1-tplus)/betaa);
            chi_lambda2_a = -Theta1a*iapp/sqrt(lambda2a)/(lambda1a-lambda2a)/sinh(sqrt(lambda2a))-Theta1a*iapp/sqrt(lambda2a)/(lambda1a-lambda2a)/tanh(sqrt(lambda2a))*(sigmaeffa/kappaeffa-(sbar_a+Theta1a-lambda1a)*nu*sigmaeffa*xi_star_a/La^2/aa/iapp/(1-tplus)/betaa);
            C_7_a = -La*iapp/sigmaeffa*(1-La^2*aa*F*betaa/iapp/sigmaeffa*((sbar_a-lambda1a)/Theta1a/sqrt(lambda1a)*C_2_a+(sbar_a-lambda2a)/Theta1a/sqrt(lambda2a)*C_4_a));
            C_8_a = -La/sigmaeffa*((sbar_a-lambda1a)/Theta1a*chi_lambda1_a*(1-La^2*aa*F*betaa/sigmaeffa/lambda1a)+(sbar_a-lambda2a)/Theta1a*chi_lambda2_a*(1-La^2*aa*F*betaa/sigmaeffa/lambda2a))+Lsep*iapp/kappaeffs*(1+(cbar_s_x0-cbar_s_x1)/iapp)-C_7_a;
            phi1_a_x0 = -La^3*aa*F*betaa/sigmaeffa^2*((sbar_a-lambda1a)*C_1_a/Theta1a/lambda1a+(sbar_a-lambda2a)*C_3_a/Theta1a/lambda2a)+C_8_a;

            c_imp(k,i) = (RL2*L2*s_c)/(RL2 + L2*s_c) -(phi1_c_x1-phi2_s_x1)/iapp;
            a_imp(k,i) = (RL1*L1*s_a)/(RL1 + L1*s_a) -(phi2_s_x0-phi1_a_x0)/iapp;
            s_imp(k,i) = -(phi2_s_x1-phi2_s_x0)/iapp;
            fc_imp(k,i) = -(phi1_c_x1-phi1_a_x0)/iapp;

            Z_cta_tot(k,i) = Rcta/A_tota/(1 + s_a*Cdla*Rcta);
            Z_ctc_tot(k,i) = Rctc/A_totc/(1 + s_c*Cdlc*Rctc);
            Z_cta(k,i) = Rcta/A_tota/(1 + s_a*Cdla*A_tota*Rcta/A_tota);
            Z_ctc(k,i) = Rctc/A_totc/(1 + s_c*Cdlc*A_totc*Rctc/A_totc);
            Z_difa(k,i) = Rdifa/A_tota/Ysa;
            Z_difc(k,i) = Rdifc/A_totc/Ysc;
        end

        parasa(:,i) = [R_itsca(i), i0a, Cdla, Dsa, aa, L1, RL1, taua, Di0, kappa]';
        parasc(:,i) = [R_itscc(i), i0c, Cdlc, Dsc, ac, L2, RL2, tauc, Di0, kappa]';
    end

    for i = 1:length(multi_soc_range)*2
        if i <= length(multi_soc_range)
            output(:,2*i-1) = R_itsca(i) + (1/A_coat)*real(a_imp(:,i));
            output(:,2*i) = (1/A_coat)*imag(a_imp(:,i));
            paras(:,i) = [parasa(:,i)];
        elseif i > length(multi_soc_range)
            output(:,2*i-1) = R_itscc(i - length(multi_soc_range)) + (1/A_coat)*real(c_imp(:,(i - length(multi_soc_range))));
            output(:,2*i) = (1/A_coat)*imag(c_imp(:,(i - length(multi_soc_range))));
            paras(:,i) = [parasc(:,(i-length(multi_soc_range)))];
        end
    end

    sep_elem.Z_cta = Z_cta;

    sep_elem.Z_ctc = Z_ctc;

    sep_elem.Z_difa = Z_difa;

    sep_elem.Z_difc = Z_difc;

end
