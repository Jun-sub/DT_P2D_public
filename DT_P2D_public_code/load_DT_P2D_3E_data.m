function [f_data,z_integ_data,loaded_table] = load_DT_P2D_3E_data(root_folder, DOE, set_name, multi_soc_range, multi_soc_anode, multi_soc_cathode)
% z_integ_data has 4*N columns:
%   anode   SOC i: columns 2*i-1 and 2*i
%   cathode SOC i: columns 2*N+2*i-1 and 2*N+2*i
% Each pair is [real(Z), imag(Z)].

    N = length(multi_soc_range);
    f_data = [];
    z_integ_data = [];
    rows = {};

    for pass = 1:2
        if pass == 1
            electrode = 'Anode';
            selected = multi_soc_anode;
            colOffset = 0;
        else
            electrode = 'Cathode';
            selected = multi_soc_cathode;
            colOffset = 2*N;
        end

        for i = 1:N
            soc = multi_soc_range(i);
            cols = colOffset + (2*i-1:2*i);

            if ismember(soc, selected)
                fileName = make_DT_P2D_file_name(DOE,set_name,electrode,soc);
                filePath = fullfile(root_folder,fileName);

                if ~isfile(filePath)
                    error('Raw EIS file not found: %s', filePath);
                end

                raw = readmatrix(filePath);
                if size(raw,2) < 3
                    error('Expected [frequency, real, imag] columns in: %s', filePath);
                end

                freq = raw(:,1);
                Z = raw(:,2:3);
                [freq,ord] = sort(freq,'descend');
                Z = Z(ord,:);

                if isempty(f_data)
                    f_data = freq(:);
                    z_integ_data = zeros(numel(f_data),4*N);
                else
                    tol = max(1e-9,1e-8*max(abs(f_data)));
                    if numel(freq) ~= numel(f_data) || any(abs(freq(:)-f_data(:)) > tol)
                        error('Frequency grid mismatch in file: %s', filePath);
                    end
                end

                z_integ_data(:,cols) = Z;
                rows(end+1,:) = {string(DOE), string(set_name), string(electrode), soc, string(fileName), numel(freq)};
            end
        end
    end

    if isempty(f_data)
        error('No selected raw EIS data were loaded.');
    end

    loaded_table = cell2table(rows,'VariableNames',{'DOE','Set','Electrode','SOC','FileName','Npoints'});
end
