function type_anode = get_anode_type(DOE)
    if strcmpi(DOE,'NG_PoorE')
        type_anode = 0;
    elseif strcmpi(DOE,'NG_OptE')
        type_anode = 1;
    elseif strcmpi(DOE,'BG_OptE')
        type_anode = 2;
    else
        error('Unknown DOE: %s', DOE);
    end
end
