function fig = plot_electrode_soc_grid(f_data, multi_soc_range, soc_list, electrode, ...
    z_data, z_model, z_model_dist, colBase, ...
    R_ohm_e, Z_Kel_e, Z_ct_e, Z_Del_e, Z_Ds_e, ...
    fig_folder, save_results)

    if isempty(soc_list)
        fig = [];
        return;
    end

    soc_list = soc_list(:).';
    nSOC = numel(soc_list);

    fig = figure('Name',[electrode ' DT-P2D element separation']);
    set(gcf,'position',[100 100 1000 1000])

    tiledlayout(nSOC, 2, 'TileSpacing','compact', 'Padding','compact');

    for row = 1:nSOC
        socNow = soc_list(row);
        plot_idx = find(multi_soc_range == socNow, 1);
        if isempty(plot_idx)
            continue;
        end

        cols = colBase + (2*plot_idx-1 : 2*plot_idx);

        nexttile;
        plot_one_panel(f_data, false, electrode, socNow, ...
            z_data(:,cols), z_model(:,cols), z_model_dist(:,cols), ...
            R_ohm_e(:,plot_idx), Z_Kel_e(:,plot_idx), Z_ct_e(:,plot_idx), ...
            Z_Del_e(:,plot_idx), Z_Ds_e(:,plot_idx));

        nexttile;
        plot_one_panel(f_data, true, electrode, socNow, ...
            z_data(:,cols), z_model(:,cols), z_model_dist(:,cols), ...
            R_ohm_e(:,plot_idx), Z_Kel_e(:,plot_idx), Z_ct_e(:,plot_idx), ...
            Z_Del_e(:,plot_idx), Z_Ds_e(:,plot_idx));
    end

    if save_results
        savefig(fig, fullfile(fig_folder, string(electrode) + "_DT_P2D_element_separation.fig"));
    end
end

function plot_one_panel(f_data, isZoom, electrode, socNow, zData, zBase, zDist, ...
    R_ohm, Z_Kel, Z_ct, Z_Del, Z_Ds)

    c_map = [0 0.450980392156863 0.760784313725490;
             0.937254901960784 0.752941176470588 0;
             0.803921568627451 0.325490196078431 0.298039215686275;
             0.125490196078431 0.521568627450980 0.305882352941177;
             0.572549019607843 0.368627450980392 0.623529411764706;
             0.882352941176471 0.529411764705882 0.152941176470588;
             0.301960784313725 0.733333333333333 0.835294117647059;
             0.933333333333333 0.298039215686275 0.592156862745098;
             0.494117647058824 0.380392156862745 0.282352941176471;
             0.454901960784314 0.462745098039216 0.470588235294118];

    c_ohm = c_map(10,:);
    c_kel = c_map(6,:);
    c_ct  = c_map(8,:);
    c_del = c_map(5,:);
    c_ds  = c_map(4,:);

    if strcmpi(electrode, 'Cathode')
        c_dt = [0.9330 0.7450 0.3010];
    else
        c_dt = [0.3010 0.7450 0.9330];
    end

    if isZoom
        f_zoom_lb = 10;
        idxAxis = f_data > f_zoom_lb;
        if ~any(idxAxis)
            idxAxis = true(size(f_data));
        end
    else
        idxAxis = true(size(f_data));
    end

    % Use the zoom frequency range only for the axis limit. The curves and
    % markers are plotted using all frequency points so that the fitted model
    % traces are not artificially truncated in the zoom-in panels.
    idxPlot = true(size(f_data));

    if strcmpi(electrode, 'Cathode')
        c_p2d = [1 0 0];
    else
        c_p2d = [0 0 1];
    end

    hold on; box on;

    plot(zData(idxPlot,1), -zData(idxPlot,2), 'ok', ...
        'LineWidth', 1.1, 'MarkerSize', 6.0, 'DisplayName','Exp Data');
    plot(zBase(idxPlot,1), -zBase(idxPlot,2), 'o', ...
        'LineWidth', 1.1, 'MarkerSize', 5.8, 'Color', c_p2d, 'DisplayName','P2D');
    plot(zDist(idxPlot,1), -zDist(idxPlot,2), 'o', ...
        'LineWidth', 1.1, 'MarkerSize', 6.0, 'Color', c_dt, 'DisplayName','DT-P2D');

    rOhm = R_ohm(find(isfinite(R_ohm),1,'first'));
    if isempty(rOhm)
        rOhm = 0;
    end

    plot([0 real(rOhm)], [0 0], '-', ...
        'LineWidth', 2.5, 'Color', c_ohm, 'DisplayName','Z_{ohm}');
    plot(real(Z_Kel), -imag(Z_Kel), '-', ...
        'LineWidth', 2.5, 'Color', c_kel, 'DisplayName','Z_{Kel}');
    plot(real(Z_ct), -imag(Z_ct), '-', ...
        'LineWidth', 2.5, 'Color', c_ct, 'DisplayName','Z_{ct}');
    plot(real(Z_Del), -imag(Z_Del), '-', ...
        'LineWidth', 2.5, 'Color', c_del, 'DisplayName','Z_{Del}');
    plot(real(Z_Ds), -imag(Z_Ds), '-', ...
        'LineWidth', 2.5, 'Color', c_ds, 'DisplayName','Z_{Ds}');

    axis_limit = axis_limit_from_eis_data(zData(idxAxis,:));

    set(gca,'Box','on', ...
        'PlotBoxAspectRatio',[1 1 1], ...
        'FontUnits','points', 'FontSize',10, 'FontName','Times New Roman', ...
        'XLim',[0 axis_limit], 'YLim',[0 axis_limit])

    grid off;
    xlabel('Z_{re} [Ohm]');
    ylabel('-Z_{im} [Ohm]');

    if isZoom
        title(sprintf('%s SOC%d P2D + dist Zoom-in', electrode, socNow));
    else
        title(sprintf('%s SOC%d P2D + dist', electrode, socNow));
    end

    legend('Location','best');
    hold off;
end

function axis_limit = axis_limit_from_eis_data(zData)
    vals = abs(zData(:));
    vals = vals(isfinite(vals));

    if isempty(vals)
        axis_limit = 1;
        return;
    end

    axis_limit = 1.1 * max(vals);

    if axis_limit <= 0
        axis_limit = 1;
    end
end
