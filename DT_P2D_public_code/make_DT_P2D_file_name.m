function fileName = make_DT_P2D_file_name(DOE,set_name,electrode,soc)
    if strcmpi(DOE,'NG_PoorE')
        if strlength(string(set_name)) == 0
            error('NG_PoorE requires set_name, e.g., ''Set1''.');
        end
        fileName = sprintf('%s_%s_%s_SOC%d.csv', DOE, set_name, electrode, soc);
    else
        fileName = sprintf('%s_%s_SOC%d.csv', DOE, electrode, soc);
    end
end
