function [R_ohma_dist, Z_Kela_dist, Z_cta_dist, Z_Dela_dist, Z_Dsa_dist, ...
          R_ohmc_dist, Z_Kelc_dist, Z_ctc_dist, Z_Delc_dist, Z_Dsc_dist, ...
          R_ohm_dist, Z_Kel_dist, Z_ct_dist, Z_Del_dist, Z_Ds_dist] = ...
    DTP2D_element_separation( ...
          f_data, z_model_dist, ...
          Z_cta_dist_ini, Z_ctc_dist_ini, ...
          Z_difa_dist_ini, Z_difc_dist_ini, ...
          Z_cta_tot_dist_ini, Z_ctc_tot_dist_ini, ...
          Z_Kel_ct_dist_LCS, multi_soc_range)

    N = length(multi_soc_range);
    nFreq = length(f_data);

    Z_Kela_dist = zeros(nFreq,N);
    Z_cta_dist  = zeros(nFreq,N);
    Z_Dela_dist_ini = zeros(nFreq,N);
    Z_Dela_dist = zeros(nFreq,N);
    Z_Dsa_dist  = zeros(nFreq,N);

    Z_Kelc_dist = zeros(nFreq,N);
    Z_ctc_dist  = zeros(nFreq,N);
    Z_Delc_dist_ini = zeros(nFreq,N);
    Z_Delc_dist = zeros(nFreq,N);
    Z_Dsc_dist  = zeros(nFreq,N);

    Z_Kel_dist = zeros(nFreq,N);
    Z_ct_dist  = zeros(nFreq,N);
    Z_Del_dist_ini = zeros(nFreq,N);
    Z_Del_dist = zeros(nFreq,N);
    Z_Ds_dist  = zeros(nFreq,N);

    R_ohma_dist = zeros(nFreq,N);
    R_ohmc_dist = zeros(nFreq,N);
    R_ohm_dist  = zeros(nFreq,N);

    for i = 1:N
        Z_modela = complex(z_model_dist(:,2*i-1), z_model_dist(:,2*i));
        Z_modelc = complex(z_model_dist(:,2*N + 2*i-1), z_model_dist(:,2*N + 2*i));
        Z_model = Z_modela + Z_modelc;

        Za_Kel_ct_now = complex(Z_Kel_ct_dist_LCS(:,2*i-1), Z_Kel_ct_dist_LCS(:,2*i));
        Zc_Kel_ct_now = complex(Z_Kel_ct_dist_LCS(:,2*N + 2*i-1), Z_Kel_ct_dist_LCS(:,2*N + 2*i));

        Z_Kela_dist(:,i) = Za_Kel_ct_now - Z_cta_tot_dist_ini(:,i);
        Z_Kelc_dist(:,i) = Zc_Kel_ct_now - Z_ctc_tot_dist_ini(:,i);
        Z_Kel_dist(:,i)  = Z_Kela_dist(:,i) + Z_Kelc_dist(:,i);

        Z_Dela_dist_ini(:,i) = Z_modela ...
            - Z_Kela_dist(:,i) ...
            - Z_cta_tot_dist_ini(:,i) ...
            - Z_difa_dist_ini(:,i);

        Z_Delc_dist_ini(:,i) = Z_modelc ...
            - Z_Kelc_dist(:,i) ...
            - Z_ctc_tot_dist_ini(:,i) ...
            - Z_difc_dist_ini(:,i);

        Z_Del_dist_ini(:,i) = Z_Dela_dist_ini(:,i) + Z_Delc_dist_ini(:,i);

        idxA = find(imag(Z_modela) <= 0, 1);
        idxC = find(imag(Z_modelc) <= 0, 1);
        idxT = find(imag(Z_model) <= 0, 1);
        if isempty(idxA), idxA = 1; end
        if isempty(idxC), idxC = 1; end
        if isempty(idxT), idxT = 1; end

        R_ohma_dist(:,i) = real(Z_modela(idxA));
        R_ohmc_dist(:,i) = real(Z_modelc(idxC));
        R_ohm_dist(:,i)  = real(Z_model(idxT));

        Z_cta_dist(:,i) = Z_cta_dist_ini(:,i) ...
            + real(Z_Kela_dist(end,i) - Z_cta_dist_ini(1,i));

        Z_ctc_dist(:,i) = Z_ctc_dist_ini(:,i) ...
            + real(Z_Kelc_dist(end,i) - Z_ctc_dist_ini(1,i));

        Z_ct_dist(:,i) = Z_cta_dist(:,i) + Z_ctc_dist(:,i);

        Z_Dela_dist(:,i) = Z_Dela_dist_ini(:,i) ...
            + real(Z_cta_dist(end,i) - Z_Dela_dist_ini(1,i));

        Z_Delc_dist(:,i) = Z_Delc_dist_ini(:,i) ...
            + real(Z_ctc_dist(end,i) - Z_Delc_dist_ini(1,i));

        Z_Del_dist(:,i) = Z_Del_dist_ini(:,i) ...
            + real(Z_ct_dist(end,i) - Z_Del_dist_ini(1,i));

        Z_Dsa_dist(:,i) = Z_difa_dist_ini(:,i) ...
            + Z_Dela_dist(end,i) - Z_difa_dist_ini(1,i);

        Z_Dsc_dist(:,i) = Z_difc_dist_ini(:,i) ...
            + Z_Delc_dist(end,i) - Z_difc_dist_ini(1,i);

        Z_Ds_dist(:,i) = Z_Dsa_dist(:,i) + Z_Dsc_dist(:,i);
    end
end
