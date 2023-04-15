obj = plasma_ST_profiles

%Requested parameter names
query_param_names_column = cell(0);
query_param_names_column{end+1,1} = 'T_e_estimate';
query_param_names_column{end+1,1} = 'n_e_estimate';
query_param_names_column{end+1,1} = 'V_plasma_estimate';
query_param_names_column = query_param_names_column';
qpn_length = length(query_param_names_column);
MD_arrays_column = cell(1,qpn_length);

MD_arrays_ca = [query_param_names_column;MD_arrays_column]';

switch obj.st_profile_type
    case 'stationary'
        array_size_vector = [obj.Nx,obj.Ny,2];
    case 'timeseries'
        array_size_vector = [obj.Nx,obj.Ny,obj.Nt,2];
end

for ind = 1:qpn_length

    MD_arrays_ca{ind,2} = cell(1,obj.Nr);
    for nr = 1:obj.Nr
        MD_arrays_ca{ind,2}{1,nr} = nan(array_size_vector);
    end

        switch obj.st_profile_type
            case 'stationary'
                for nx = 1:obj.Nx
                    for ny = 1:obj.Ny
                        for nr = 1:obj.Nr
                            for nrf = 1:2
                                if obj.Nr > 1    
                                    var_str = ['obj.pp_ca{',num2str(nx),',',num2str(ny),',',num2str(nr),'}.fit_iv_curve_out_ca{1,',num2str(nrf),'}.',MD_arrays_ca{ind,1}];
                                else
                                    var_str = ['obj.pp_ca{',num2str(nx),',',num2str(ny),'}.fit_iv_curve_out_ca{1,',num2str(nrf),'}.',MD_arrays_ca{ind,1}];
                                end
                                try
                                    disp(eval(var_str));
                                    MD_arrays_ca{ind,2}{1,nr}(nx,ny,nrf) = eval(var_str);
                                catch err
                                    disp(err)
                                end
                            end
                        end
                    end
                end
            case 'timeseries'
                for nx = 1:obj.Nx
                    for ny = 1:obj.Ny
                        for nt = 1:obj.Nt
                            for nr = 1:obj.Nr
                                for nrf = 1:2
                                    if obj.Nr > 1
                                        var_str = ['obj.pp_ca{',num2str(nx),',',num2str(ny),',',num2str(nt),',',num2str(nr),'}.fit_iv_curve_out_ca{1,',num2str(nrf),'}.',MD_arrays_ca{ind,1}];
                                    else
                                        var_str = ['obj.pp_ca{',num2str(nx),',',num2str(ny),'}.fit_iv_curve_out_ca{1,',num2str(nrf),'}.',MD_arrays_ca{ind,1}];
                                    end
                                    try
                                        disp(eval(var_str));
                                        MD_arrays_ca{ind,2}{1,nr}(nx,ny,nt,nrf) = eval(var_str);
                                    catch err
                                        disp(err)
                                    end
                                end
                            end
                        end
                    end
                end
        end

end


    
obj.MD_arrays_ca = MD_arrays_ca;