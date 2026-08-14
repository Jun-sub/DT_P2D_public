%% run_DT_P2D_MSCEF_example.m
% Supplementary MSCEF fitting example for the DT-P2D model.
% Each EIS CSV file contains Frequency [Hz], Real(Z) [Ohm], and Imag(Z) [Ohm].
% Set root_folder and save_folder to "AUTO" to use ./Example_Data and ./Save_example.
% The script performs preliminary P2D fitting, DT-P2D fitting, impedance-element separation,
% and exports the summary table, figures, and MATLAB workspace.

clear; clc; close all;

root_folder = "G:\공유 드라이브\Battery Software Group (2025)\Papers\2026_EIS (DT-P2D)\EES\EES_batteries\1st_revision\추가 피규어 및 테이블, 데이터, 공개용 코드\DT_P2D_public_code\변경 버전\Example_Data";
save_folder = "G:\공유 드라이브\Battery Software Group (2025)\Papers\2026_EIS (DT-P2D)\EES\EES_batteries\1st_revision\추가 피규어 및 테이블, 데이터, 공개용 코드\DT_P2D_public_code\변경 버전\Save_example";

DOE = 'NG_PoorE';
set_name = 'Set1';
multi_soc_anode = [20 60];
multi_soc_cathode = [20 60];

num_iter = 100;
num_iter_dist = 30;
start_points = 5;
use_parallel = 1;
use_relative_weight = true;
save_results = true;

T = 298.15;

script_folder = string(fileparts(mfilename('fullpath')));

if strlength(string(root_folder)) == 0 || strcmpi(string(root_folder), "AUTO")

    root_folder = fullfile(script_folder, "Example_Data");
else

    root_folder = string(root_folder);
end

if strlength(string(save_folder)) == 0 || strcmpi(string(save_folder), "AUTO")
    save_folder = fullfile(script_folder, "Save_example");
else
    save_folder = string(save_folder);
end

if ~isfolder(root_folder)
    error(['EIS data folder was not found:\n%s\n\n' ...
        'For AUTO mode, create the following folder and place the CSV files inside it:\n' ...
        '<code folder>\\Example_Data'], root_folder);
end

if save_results
    ensure_folder(save_folder);
end

fprintf('Data folder   : %s\n', root_folder);
fprintf('Output folder : %s\n', save_folder);

multi_soc_range = unique([multi_soc_anode multi_soc_cathode]);
type_anode = get_anode_type(DOE);
N = length(multi_soc_range);

case_tag = string(DOE);
if strcmpi(DOE, 'NG_PoorE') && strlength(string(set_name)) > 0
    case_tag = case_tag + "_" + string(set_name);
end

[f_data, z_integ_data, data_loaded] = load_DT_P2D_3E_data(root_folder, DOE, set_name, multi_soc_range, multi_soc_anode, multi_soc_cathode);
fprintf('Loaded %d selected EIS datasets. Common frequency points: %d\n', height(data_loaded), numel(f_data));
disp(data_loaded);

weight = zeros(length(f_data), 4*N);
for i = 1:2*N
    if i <= N
        col = 2*i-1;
    else
        col = 2*i-1;
    end
    znow_abs = sqrt(z_integ_data(:,col).^2 + z_integ_data(:,col+1).^2);
    if use_relative_weight
        weight(:,col) = 1 ./ max(znow_abs, eps);
        weight(:,col+1) = weight(:,col);
    else
        weight(:,col:col+1) = 1;
    end
end

weight_idx = ones(1,4*N);
weight_idx(1:2*N) = repelem(ismember(multi_soc_range,multi_soc_anode),2);
weight_idx(2*N+1:4*N) = repelem(ismember(multi_soc_range,multi_soc_cathode),2);
weight(:,weight_idx == 0) = 0;
weighted_data = z_integ_data .* weight;

fprintf('\nStart preliminary P2D fitting\n');

factors_ini = ones(1,10);
factors_integ_ini = ones(10,2*N+1);
factors_integ_ini(1,2*N+1) = 1;
factors_integ_ini(2,2*N+1) = 1;

lb = factors_integ_ini*0.01;
ub = factors_integ_ini*50;
lb(8,1:2*N) = 0.7;  ub(8,1:2*N) = 1.3;
lb(9,1:2*N) = 0.15; ub(9,1:2*N) = 3;
lb(3:10,2*N+1) = 1; ub(3:10,2*N+1) = 1;

bound_idx = ones(1,2*N);
bound_idx(1:N) = ismember(multi_soc_range,multi_soc_anode);
bound_idx(N+1:2*N) = ismember(multi_soc_range,multi_soc_cathode);
lb(:,bound_idx == 0) = factors_ini(1);
ub(:,bound_idx == 0) = factors_ini(1);

model_base_weighted = @(factors,freq) DTP2D_model_base(freq,factors,multi_soc_range,T,type_anode).*weight;
options_base = optimoptions('lsqcurvefit', 'Display','iter', 'MaxIterations',num_iter, ...
    'MaxFunctionEvaluations',1e8, 'FunctionTolerance',1e-6, 'StepTolerance',1e-6, ...
    'FiniteDifferenceType','central', 'Algorithm','trust-region-reflective');

[factors_integ_hat, ~] = DTP2D_multistart_lsqcurvefit(model_base_weighted, factors_integ_ini, f_data, weighted_data, lb, ub, options_base, start_points, use_parallel);
[z_model, para_base] = DTP2D_model_base(f_data,factors_integ_hat,multi_soc_range,T,type_anode);

fprintf('\nStart DT-P2D distributed fitting\n');

factors_integ_ini_dist = ones(12,2*N+1);
factors_integ_ini_dist(1:10,1:2*N+1) = factors_integ_hat;
factors_integ_ini_dist(11,1:2*N) = 0.5;
factors_integ_ini_dist(12,1:N) = 5;
factors_integ_ini_dist(12,N+1:2*N) = 0.5;

lb_dist = factors_integ_ini_dist*0.1;
ub_dist = factors_integ_ini_dist*5.0;
lb_dist(1,2*N+1) = 0.01; ub_dist(1,2*N+1) = 50;
lb_dist(8,1:2*N) = 0.7;   ub_dist(8,1:2*N) = 1.3;
lb_dist(:,bound_idx == 0) = factors_ini(1);
ub_dist(:,bound_idx == 0) = factors_ini(1);
lb_dist(11,1:2*N) = 0.01; ub_dist(11,1:2*N) = 5;
lb_dist(12,1:2*N) = 0.01; ub_dist(12,1:2*N) = 15;
lb_dist(3:12,2*N+1) = 1; ub_dist(3:12,2*N+1) = 1;
lb_dist(11:12,bound_idx == 0) = 0.5;
ub_dist(11:12,bound_idx == 0) = 0.5;

model_dist_weighted = @(factors,freq) DTP2D_model_distributed(freq,factors,multi_soc_range,T,type_anode).*weight;
options_dist = optimoptions('lsqcurvefit', 'Display','iter', 'MaxIterations',num_iter_dist, ...
    'MaxFunctionEvaluations',1e8, 'FunctionTolerance',1e-8, 'StepTolerance',1e-8, ...
    'FiniteDifferenceType','central', 'Algorithm','trust-region-reflective');

if num_iter_dist == 0
    factors_integ_hat_dist = factors_integ_ini_dist;
else
    factors_integ_hat_dist = lsqcurvefit(model_dist_weighted,factors_integ_ini_dist,f_data,weighted_data,lb_dist,ub_dist,options_dist);
end

factors_kel = factors_integ_hat_dist;
factors_kel(4,1:2*N) = 1e8;
factors_kel(1,2*N+1) = 1e8;
factors_kel(11,2*N+1) = 6;
[Z_Kel_ct_dist_LCS, ~, ~] = DTP2D_model_distributed(f_data,factors_kel,multi_soc_range,T,type_anode);
[z_model_dist, para_dist, sep_elem_dist] = DTP2D_model_distributed(f_data,factors_integ_hat_dist,multi_soc_range,T,type_anode);

[R_ohma_dist, Z_Kela_dist, Z_cta_dist, Z_Dela_dist, Z_Dsa_dist, ...
    R_ohmc_dist, Z_Kelc_dist, Z_ctc_dist, Z_Delc_dist, Z_Dsc_dist, ...
    ~, ~, ~, ~, ~] = ...
    DTP2D_element_separation(f_data, z_model_dist, ...
    sep_elem_dist.Z_cta_dist, sep_elem_dist.Z_ctc_dist, ...
    sep_elem_dist.Z_difa_dist, sep_elem_dist.Z_difc_dist, ...
    sep_elem_dist.Z_cta_tot_dist, sep_elem_dist.Z_ctc_tot_dist, ...
    Z_Kel_ct_dist_LCS, multi_soc_range);

if save_results

    summary_table = make_DT_P2D_summary_table( ...
        case_tag, multi_soc_range, multi_soc_anode, multi_soc_cathode, ...
        para_dist, z_integ_data, z_model_dist, ...
        R_ohma_dist, Z_Kela_dist, Z_cta_dist, Z_Dela_dist, Z_Dsa_dist, ...
        R_ohmc_dist, Z_Kelc_dist, Z_ctc_dist, Z_Delc_dist, Z_Dsc_dist);

    writecell(summary_table, fullfile(save_folder, case_tag + "_summary_table.csv"));

    save(fullfile(save_folder, case_tag + "_DT_P2D_workspace.mat"));
end

if save_results

    for socNow = multi_soc_anode(:).'
        figA = plot_electrode_soc_grid(f_data, multi_soc_range, socNow, 'Anode', ...
            z_integ_data, z_model, z_model_dist, 0, ...
            R_ohma_dist, Z_Kela_dist, Z_cta_dist, Z_Dela_dist, Z_Dsa_dist, ...
            save_folder, false);

        if ~isempty(figA)
            savefig(figA, fullfile(save_folder, ...
                case_tag + "_Anode_SOC" + string(socNow) + "_DT_P2D_element_separation.fig"));
        end
    end

    for socNow = multi_soc_cathode(:).'
        figC = plot_electrode_soc_grid(f_data, multi_soc_range, socNow, 'Cathode', ...
            z_integ_data, z_model, z_model_dist, 2*N, ...
            R_ohmc_dist, Z_Kelc_dist, Z_ctc_dist, Z_Delc_dist, Z_Dsc_dist, ...
            save_folder, false);

        if ~isempty(figC)
            savefig(figC, fullfile(save_folder, ...
                case_tag + "_Cathode_SOC" + string(socNow) + "_DT_P2D_element_separation.fig"));
        end
    end
else
    plot_electrode_soc_grid(f_data, multi_soc_range, multi_soc_anode, 'Anode', ...
        z_integ_data, z_model, z_model_dist, 0, ...
        R_ohma_dist, Z_Kela_dist, Z_cta_dist, Z_Dela_dist, Z_Dsa_dist, ...
        save_folder, false);

    plot_electrode_soc_grid(f_data, multi_soc_range, multi_soc_cathode, 'Cathode', ...
        z_integ_data, z_model, z_model_dist, 2*N, ...
        R_ohmc_dist, Z_Kelc_dist, Z_ctc_dist, Z_Delc_dist, Z_Dsc_dist, ...
        save_folder, false);
end

fprintf('\nDone. Output folder: %s\n', save_folder);

function C = make_DT_P2D_summary_table(case_tag, multi_soc_range, multi_soc_anode, multi_soc_cathode, ...
        para_dist, z_data, z_model_dist, ...
        R_ohma_dist, Z_Kela_dist, Z_cta_dist, Z_Dela_dist, Z_Dsa_dist, ...
        R_ohmc_dist, Z_Kelc_dist, Z_ctc_dist, Z_Delc_dist, Z_Dsc_dist)

    N = length(multi_soc_range);

    headers = {'Paras','Unit'};
    colMeta = struct('electrode',{},'soc',{},'colPara',{},'colZ',{});

    for socNow = multi_soc_anode(:).'
        idx = find(multi_soc_range == socNow, 1);
        if isempty(idx)
            continue;
        end
        headers{end+1} = sprintf('Anode_SOC%d', socNow);
        colMeta(end+1).electrode = 'Anode';
        colMeta(end).soc = socNow;
        colMeta(end).colPara = idx;
        colMeta(end).colZ = 2*idx - 1;
    end

    for socNow = multi_soc_cathode(:).'
        idx = find(multi_soc_range == socNow, 1);
        if isempty(idx)
            continue;
        end
        headers{end+1} = sprintf('Cathode_SOC%d', socNow);
        colMeta(end+1).electrode = 'Cathode';
        colMeta(end).soc = socNow;
        colMeta(end).colPara = N + idx;
        colMeta(end).colZ = 2*N + 2*idx - 1;
    end

    C = cell(0, numel(headers));
    C(end+1,:) = headers;

    C(end+1,:) = make_summary_row('i0', '[A/m^2]', 2, para_dist, colMeta);
    C(end+1,:) = make_summary_row('Cdl', '[F/m^2]', 3, para_dist, colMeta);
    C(end+1,:) = make_summary_row('Ds', '[m^2/s]', 4, para_dist, colMeta);
    C(end+1,:) = make_summary_row('av', '[m^2/m^3]', 5, para_dist, colMeta);
    C(end+1,:) = make_summary_row('tau', '[-]', 8, para_dist, colMeta);
    C(end+1,:) = make_summary_row('Del', '[m^2/s]', 9, para_dist, colMeta);
    C(end+1,:) = make_summary_row('Kel', '[S/m]', 10, para_dist, colMeta);
    C(end+1,:) = make_summary_row('sigma_DRT', '[-]', 11, para_dist, colMeta);
    C(end+1,:) = make_summary_row('sigma_DDT', '[-]', 12, para_dist, colMeta);

    C(end+1,:) = make_resistance_row('R_ohm', '[Ohm]', colMeta, multi_soc_range, ...
        R_ohma_dist, R_ohmc_dist, 'ohm');
    C(end+1,:) = make_resistance_row('R_Kel', '[Ohm]', colMeta, multi_soc_range, ...
        Z_Kela_dist, Z_Kelc_dist, 'width');
    C(end+1,:) = make_resistance_row('R_ct', '[Ohm]', colMeta, multi_soc_range, ...
        Z_cta_dist, Z_ctc_dist, 'width');
    C(end+1,:) = make_resistance_row('R_Del', '[Ohm]', colMeta, multi_soc_range, ...
        Z_Dela_dist, Z_Delc_dist, 'width');
    C(end+1,:) = make_resistance_row('R_Ds', '[Ohm]', colMeta, multi_soc_range, ...
        Z_Dsa_dist, Z_Dsc_dist, 'width');

    C(end+1,:) = make_rmse_row('RMSE', '[Ohm]', colMeta, z_data, z_model_dist, 'ohm');
    C(end+1,:) = make_rmse_row('RMSE', '[%]', colMeta, z_data, z_model_dist, 'percent');

    titleRow = cell(1, numel(headers));
    titleRow(:) = {''};
    titleRow{1} = char(case_tag);
    C = [titleRow; C];
end

function row = make_summary_row(parameter_name, unit_text, rowIdx, para_matrix, colMeta)

    row = cell(1, numel(colMeta) + 2);
    row{1} = parameter_name;
    row{2} = unit_text;

    for k = 1:numel(colMeta)
        row{2+k} = para_matrix(rowIdx, colMeta(k).colPara);
    end
end

function row = make_resistance_row(parameter_name, unit_text, colMeta, multi_soc_range, ...
        anode_component, cathode_component, mode)

    row = cell(1, numel(colMeta) + 2);
    row{1} = parameter_name;
    row{2} = unit_text;

    for k = 1:numel(colMeta)
        idx = find(multi_soc_range == colMeta(k).soc, 1);

        if strcmpi(colMeta(k).electrode, 'Anode')
            component = anode_component(:,idx);
        else
            component = cathode_component(:,idx);
        end

        if strcmpi(mode, 'ohm')
            row{2+k} = extract_ohmic_value(component);
        else
            row{2+k} = extract_component_width(component);
        end
    end
end

function row = make_rmse_row(parameter_name, unit_text, colMeta, z_data, z_model, mode)

    row = cell(1, numel(colMeta) + 2);
    row{1} = parameter_name;
    row{2} = unit_text;

    for k = 1:numel(colMeta)
        col = colMeta(k).colZ;

        Zexp = complex(z_data(:,col), z_data(:,col+1));
        Zfit = complex(z_model(:,col), z_model(:,col+1));

        rmse_ohm = sqrt(mean(abs(Zfit - Zexp).^2, 'omitnan'));

        if strcmpi(mode, 'percent')
            denom = sqrt(mean(abs(Zexp).^2, 'omitnan'));
            if denom <= 0 || ~isfinite(denom)
                row{2+k} = NaN;
            else
                row{2+k} = 100 * rmse_ohm / denom;
            end
        else
            row{2+k} = rmse_ohm;
        end
    end
end

function val = extract_ohmic_value(component)

    component = component(:);
    idx = find(isfinite(real(component)), 1, 'first');

    if isempty(idx)
        val = NaN;
    else
        val = real(component(idx));
    end
end

function val = extract_component_width(component)

    r = real(component(:));
    r = r(isfinite(r));

    if isempty(r)
        val = NaN;
    else
        val = max(r) - min(r);
    end
end
