classdef PlasmaSTProfiles < handle
    %PLASMASTPROFILES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        %ST profile type
        st_profile_type
        %Values admitted:
        %1) 'stationary'
        %2) 'timeseries'
        
        source_folder
        source_folder_V_float
        external_V_float_is_present
        
        x_axis
        y_axis
        x_begin
        y_begin
        delta_x
        delta_y
        x_end
        y_end
        
        Nx
        Ny
        Nr
        
        meas_number_STR
        
        %Plasma Parameter cell array
        pp_ca
        %If st_profile_type = 

        time_axis_rise
        time_axis_fall
        Nt
        
        MD_arrays_ca
        
        %Diagnostics
        h_rect_ca
        save_flag_array
        time_info_loaded = false
    end
    
    methods
        function obj = PlasmaSTProfiles(varargin)
            %Class constructor with variable input arguments
            %One of the following formats is expected:
            %PlasmaSTProfiles(st_profile_type,source_folder)
            %PlasmaSTProfiles(st_profile_type,source_folder,source_folder_V_float)
            
            if isempty(varargin)
                %No input arguments -> error
                fprintf('Error: No input data!\n');
                return
            end
            
            if nargin >= 1
                obj.st_profile_type = varargin{1};
            else
                fprintf('ST profile type not specified: ''stationary'' or ''timeseries''?\n');
                return
            end
            
            if nargin >= 2
                obj.source_folder = varargin{2};
            else
                fprintf('No IV source folder input!\n');
                return
            end
            
            if nargin >= 3
                obj.external_V_float_is_present = true;
                obj.source_folder_V_float = varargin{3};
            else
                fprintf('No V_float source folder: V_float determination from IV signals\n');
                obj.external_V_float_is_present = false;
            end
            
        end
        
        function obj = readMapFile_parallel(obj,N_repetitions,varargin)
            
            %Optional argument: map file directory, if none is provided,
            %the function will use the source folder provided at the
            %creation of the instance of the class
            map_file = 'map.map';
            if isempty(varargin)
                full_file_path = strcat(obj.source_folder,map_file);
            end
            
            if nargin >= 3
                full_file_path = varargin{2};
            end

            %Read the map file
            B = tdfread(full_file_path);
            ind = B.x0x25_idx;
            x_pos = B.x;
            y_pos = B.y;
            ind = ind(:,1);

            obj.x_axis = squeeze(unique(x_pos(:,1)));
            obj.y_axis = squeeze(unique(y_pos(:,1)));
            
            obj.x_begin = min(obj.x_axis);
            obj.y_begin = min(obj.y_axis);
            
            obj.delta_x = min(abs(diff(obj.x_axis)));
            obj.delta_y = min(abs(diff(obj.y_axis)));
            
            obj.x_end = max(obj.x_axis);
            obj.y_end = max(obj.y_axis);
            
            
            obj.Nx = length(obj.x_axis);
            obj.Ny = length(obj.y_axis);
            obj.Nr = N_repetitions;
            
            %Organize data in multidimensional array
            %Measurement number 2Dspace*1Drepetition array
            obj.meas_number_STR = NaN(obj.Nx,obj.Ny,obj.Nr);
            
            %Prepare plasma parameter timeseries array
            obj.pp_ca = cell(obj.Nx,obj.Ny,obj.Nr);
            
            %Prepare the save flag array
            obj.save_flag_array = zeros(obj.Nx,obj.Ny,obj.Nr);

            p = gcp();
            for nx = 1:obj.Nx
                for ny = 1:obj.Ny
                    for nr = 1:obj.Nr
                        %file_name = strcat('meas_',num2str(k+(repetition-1),'%04.f'),'.h5');
                        %file_names_STR(k_x,k_y,repetition) = file_name;
                        k = nx + (ny-1)*obj.Nx;
                        obj.meas_number_STR(nx,ny,nr) = obj.Nr*(k-1)+(nr-1);
                        %Build the measurement file name
                        file_number = obj.meas_number_STR(nx,ny,nr);
                        iv_signals_source_full_file_path = strcat(obj.source_folder,'meas_',num2str(file_number,'%04.f'),'.h5');
                        %Check the ST profile type (stationary or timeseries)
                        switch obj.st_profile_type
                            case 'stationary'
                                %Create a PlsmaParamTimeseries object and include it in the cell array
                                if obj.external_V_float_is_present
                                    v_float_signal_source_full_file_path = strcat(obj.source_folder_V_float,'meas_',num2str(file_number,'%04.f'),'.h5');
                                    f(nx,ny,nr) = parfeval(p,@PlasmaParamStationary,1,iv_signals_source_full_file_path,v_float_signal_source_full_file_path);
                                else
                                    f(nx,ny,nr) = parfeval(p,@PlasmaParamStationary,1,iv_signals_source_full_file_path);
                                end
                            case 'timeseries'
                                %Create a PlsmaParamTimeseries object and include it in the cell array
                                if obj.external_V_float_is_present
                                    v_float_signal_source_full_file_path = strcat(obj.source_folder_V_float,'meas_',num2str(file_number,'%04.f'),'.h5');
                                    f(nx,ny,nr) = parfeval(p,@PlasmaParamTimeseries,1,iv_signals_source_full_file_path,v_float_signal_source_full_file_path);
                                else
                                    f(nx,ny,nr) = parfeval(p,@PlasmaParamTimeseries,1,iv_signals_source_full_file_path);
                                end
                            otherwise
                                fprintf(strcat('Unknown ST profile type specification: ',obj.st_profile_type));
                        end

                        
                    end
                end
            end
            

            
            for nx = 1:obj.Nx
                for ny = 1:obj.Ny
                    for nr = 1:obj.Nr
                      [~,value] = fetchNext(f);
                      %fprintf('Got result with index: %d.\n', completedIdx);
                      fprintf('Coordinates:\n');
                      fprintf(strcat('nx = ',num2str(nx),'\n'));
                      fprintf(strcat('ny = ',num2str(ny),'\n'));
                      fprintf(strcat('nr = ',num2str(nr),'\n'));
                      obj.pp_ca{nx,ny,nr} = value;
                    end
                end
            end
            
        end
        
        function obj = extractPlasmaSTProfiles(obj,S_probe,A_ion,langmuir_iv_fit_version,options)
            %Prepare the folder for the plasma parameters timeseries
            mkdir([obj.source_folder,'plasma_parameters_timeseries_list']);

            
            fprintf('Start of plasma parameter extraction for the whole profile.\n');
            %Visual diagnostics
            f_diagn = figure;
            set(f_diagn,'Renderer','painters');
            ax_diagn = axes;
            obj.h_rect_ca = cell(obj.Ny,obj.Nx);
            for nx = 1:obj.Nx
                for ny = 1:obj.Ny
                    h_rect = rectangle('Position',[obj.x_axis(nx)-obj.delta_x/2,obj.y_axis(ny)-obj.delta_y/2,obj.delta_x,obj.delta_y],'FaceColor','red');
                    obj.h_rect_ca{ny,nx} = h_rect;
                end
            end
            axis(ax_diagn,'equal');
            
            %Prepare to copy the timing information
            time_info_loaded = false;
            
            for nx = 1:obj.Nx
                for ny = 1:obj.Ny
                    for nr = 1:obj.Nr
                        
                        %Diagnostics
                        obj.h_rect_ca{ny,nx}.FaceColor = 'yellow';
                        fprintf(['X_index = ',num2str(nx),'\n']);
                        fprintf(['Y_index = ',num2str(ny),'\n']);
                        fprintf(['Repetition_index = ',num2str(nr),'\n\n']);
                         try
                            %Build the measurement file name
                            file_number = obj.meas_number_STR(nx,ny,nr);
                            source_full_file_path = strcat(obj.source_folder,'meas_',num2str(file_number,'%04.f'),'.h5');
                            %Create a PlsmaParamTimeseries object and perform plasma parameter extraction on it
                            pp_timeseries = PlasmaParamTimeseries(source_full_file_path);
                            pp_timeseries.extractPlasmaParameters(S_probe,A_ion,langmuir_iv_fit_version,options);
                            %Save the object in the post-processing folder
                            filename = strcat('plasma_param_timeseries_',num2str(file_number,'%04.f'),'.mat');
                            save([obj.source_folder,'plasma_parameters_timeseries_list\',filename],'pp_timeseries');
                            fprintf('Plasma Parameter Timeseries saved.\n');
                            obj.h_rect_ca{ny,nx}.FaceColor = 'green';
                            %If the time information has not yet been loaded
                            if ~time_info_loaded
                                %Take the time information from the first successfully extracted plasma parameter timeseries
                                obj.time_axis_rise = pp_timeseries.ivt_source.total_time_axis_rise;
                                obj.time_axis_fall = pp_timeseries.ivt_source.total_time_axis_fall;
                                obj.Nt = pp_timeseries.ivt_source.NC_tot;
                                time_info_loaded = true;
                            end
                            obj.pp_timeseries_ca{nx,ny,nr} = pp_timeseries;
                         catch err_extract
                             fprintf('Plasma Parameter extraction failed: jumping to next interation...\n');
                             obj.h_rect_ca{ny,nx}.FaceColor = 'red';
                             continue
                         end 
                        
                    end
                end
            end
        end
        
        function obj = extractPlasmaSTProfiles_parallel(obj,S_probe,A_ion,langmuir_iv_fit_version,options)

            fprintf('Start of plasma parameter extraction for the whole profile.\n');
            %Visual diagnostics
            f_diagn = figure;
            set(f_diagn,'Renderer','opengl');
            ax_diagn = axes;
            title(ax_diagn,'Diagnostics');
            obj.h_rect_ca = cell(obj.Ny,obj.Nx);
            for nx = 1:obj.Nx
                for ny = 1:obj.Ny
                    h_rect = rectangle('Position',[obj.x_axis(nx)-obj.delta_x/2,obj.y_axis(ny)-obj.delta_y/2,obj.delta_x,obj.delta_y],'FaceColor','red');
                    obj.h_rect_ca{ny,nx} = h_rect;
                    obj.h_rect_ca{ny,nx}.FaceColor = 'red';
                end
            end
            axis(ax_diagn,'equal');
            
            %Prepare parfor index
            Nx = obj.Nx;
            Ny = obj.Ny;
            Nr = obj.Nr;
            N = int16(Nx * Ny * Nr);
            
            
            %Also the unused ones must be preapared in order for the PARFOR
            %to work (I don't know why)
            
            
            
            switch obj.st_profile_type
                case 'stationary'
                    %Prepare data queue
                    pp_stationary_dataqueue = parallel.pool.DataQueue;
                    %Prepare the folder for the plasma parameters stationary
                    mkdir([obj.source_folder,'plasma_parameters_stationary_list']);
                    %Call the nested function "assignPpstatToPosition" every time that a pp_stationary has been created
                    afterEach(pp_stationary_dataqueue, @assignPpstatToPosition);
                    parfor n = 0:N-1
                        nr = int16(mod(n,Nr)+1); %#ok<PROPLC>
                        aux = int16(idivide(n,Nr,'floor')); %#ok<PROPLC>
                        ny = int16(mod(aux,Ny)+1); %#ok<PROPLC>
                        nx = int16(idivide(aux,Ny,'floor')+1); %#ok<PROPLC>
                        try
                            %Build the measurement file name
                            file_number = obj.meas_number_STR(nx,ny,nr);
                            source_full_file_path = strcat(obj.source_folder,'meas_',num2str(file_number,'%04.f'),'.h5');
                            if obj.external_V_float_is_present
                                source_V_float_full_file_path = strcat(obj.source_folder_V_float,'meas_',num2str(file_number,'%04.f'),'.h5');
                                pp_stationary = PlasmaParamStationary(source_full_file_path,source_V_float_full_file_path);
                            else
                                pp_stationary = PlasmaParamStationary(source_full_file_path);
                            end
                            args = {pp_stationary,nx,ny,nr};
                            send(pp_stationary_dataqueue,args);
                        catch err_create
                            fprintf('Plasma Parameter Stationary object creation failed.\n');
                        end
                    end
                case 'timeseries'
                    %Prepare data queue
                    pp_timeseries_dataqueue = parallel.pool.DataQueue;
                    %Prepare the folder for the plasma parameters timeseries
                    mkdir([obj.source_folder,'plasma_parameters_timeseries_list']);
                    %Call the nested function "assignPptsToPosition" every time that a pp_timeseries has been created
                    afterEach(pp_timeseries_dataqueue, @assignPptsToPosition);
                    parfor n = 0:N-1
                        nr = int16(mod(n,Nr)+1); %#ok<PROPLC>
                        aux = int16(idivide(n,Nr,'floor')); %#ok<PROPLC>
                        ny = int16(mod(aux,Ny)+1); %#ok<PROPLC>
                        nx = int16(idivide(aux,Ny,'floor')+1); %#ok<PROPLC>
                        try
                            %Build the measurement file name
                            file_number = obj.meas_number_STR(nx,ny,nr);
                            source_full_file_path = strcat(obj.source_folder,'meas_',num2str(file_number,'%04.f'),'.h5');
                            if obj.external_V_float_is_present
                                source_V_float_full_file_path = strcat(obj.source_folder_V_float,'meas_',num2str(file_number,'%04.f'),'.h5');
                                pp_timeseries = PlasmaParamTimeseries(source_full_file_path,source_V_float_full_file_path);
                            else
                                pp_timeseries = PlasmaParamTimeseries(source_full_file_path);
                            end
                            args = {pp_timeseries,nx,ny,nr};
                            send(pp_timeseries_dataqueue,args);
                        catch err_create
                            fprintf('Plasma Parameter Timeseries object creation failed.\n');
                        end
                    end
                otherwise
                    fprintf(strcat('Unknown st_profile_type sepcification: ',obj.st_profile_type,'\n'));
                    return
            end
            
            function assignPptsToPosition(args)
                pp_timeseries_to_save = args{1};
                nx_to_save = args{2};
                ny_to_save = args{3};
                nr_to_save = args{4};
                %Diagnostics
                fprintf(['Assigning pp_timeseries number = ',num2str(nx_to_save*ny_to_save*nr_to_save-1),'\n']);
                fprintf(['X_index = ',num2str(nx_to_save),'\n']);
                fprintf(['Y_index = ',num2str(ny_to_save),'\n']);
                fprintf(['Repetition_index = ',num2str(nr_to_save),'\n\n']);
                 try
                      obj.pp_timeseries_ca{nx_to_save,ny_to_save,nr_to_save} = pp_timeseries_to_save;
                 catch err_assign_position
                      fprintf('Plasma Parameter Timeseries position assigment failed.\n');

                 end 
            end
            
            function assignPpstatToPosition(args)
                pp_stationary_to_save = args{1};
                nx_to_save = args{2};
                ny_to_save = args{3};
                nr_to_save = args{4};
                %Diagnostics
                fprintf(['Assigning pp_stationary number = ',num2str(nx_to_save*ny_to_save*nr_to_save-1),'\n']);
                fprintf(['X_index = ',num2str(nx_to_save),'\n']);
                fprintf(['Y_index = ',num2str(ny_to_save),'\n']);
                fprintf(['Repetition_index = ',num2str(nr_to_save),'\n\n']);
                 try
                      obj.pp_ca{nx_to_save,ny_to_save,nr_to_save} = pp_stationary_to_save;
                 catch err_assign_position
                      fprintf('Plasma Parameter Stationary position assigment failed.\n');

                 end 
            end
            
            
            %Prepare parfor index
            Nx = obj.Nx;
            Ny = obj.Ny;
            Nr = obj.Nr;
            N = int16(Nx * Ny * Nr);
            
            switch obj.st_profile_type
                case 'stationary'
                    %Prepare data queue
                    pp_stationary_dataqueue = parallel.pool.DataQueue;
                    %Call the nested function "savePpstatToFile" every time that a pp_stationary has finished with plasma parameter extraction
                    afterEach(pp_stationary_dataqueue, @savePpstatToFile);
                    %parfor
                    parfor n = 0:N-1
                        nr = int16(mod(n,Nr)+1); %#ok<PROPLC>
                        aux = int16(idivide(n,Nr,'floor')); %#ok<PROPLC>
                        ny = int16(mod(aux,Ny)+1); %#ok<PROPLC>
                        nx = int16(idivide(aux,Ny,'floor')+1); %#ok<PROPLC>
                        try
                            pp_stationary = obj.pp_ca{nx,ny,nr}.extractPlasmaParameters(S_probe,A_ion,langmuir_iv_fit_version,options);
                        catch err_extract
                            fprintf('Plasma Parameter Stationary data analysis failed: jumping to next task...\n');
                            pp_stationary = obj.pp_ca{nx,ny,nr};
                        end
                        args = {pp_stationary,nx,ny,nr};
                        send(pp_stationary_dataqueue,args);
                    end
                case 'timeseries'
                    %Prepare data queue
                    pp_timeseries_dataqueue = parallel.pool.DataQueue;
                    %Call the nested function "savePptsToFile" every time that a pp_timeseries has finished with plasma parameter extraction
                    afterEach(pp_timeseries_dataqueue, @savePptsToFile);
                    parfor n = 0:N-1
                        nr = int16(mod(n,Nr)+1); %#ok<PROPLC>
                        aux = int16(idivide(n,Nr,'floor')); %#ok<PROPLC>
                        ny = int16(mod(aux,Ny)+1); %#ok<PROPLC>
                        nx = int16(idivide(aux,Ny,'floor')+1); %#ok<PROPLC>

                        try
                            pp_timeseries = obj.pp_timeseries_ca{nx,ny,nr}.extractPlasmaParameters(S_probe,A_ion,langmuir_iv_fit_version,options);
                        catch err_extract
                            fprintf('Plasma Parameter Timeseries data analysis failed: jumping to next task...\n');
                            pp_timeseries = obj.pp_timeseries_ca{nx,ny,nr};
                        end
                        args = {pp_timeseries,nx,ny,nr};
                        fprintf('Here');
                        send(pp_timeseries_dataqueue,args);
                    end
            end
            
            %Nested function which saves the pp_stationary with extracted
            %plasma parameters to a file
            function savePpstatToFile(args)
                pp_stationary_to_save = args{1};
                nx_to_save = args{2};
                ny_to_save = args{3};
                nr_to_save = args{4};
                obj.pp_ca{nx_to_save,ny_to_save,nr_to_save} = pp_stationary_to_save; 
                %Diagnostics
                fprintf(['Saving file number = ',num2str(nx_to_save*ny_to_save*nr_to_save-1),'\n']);
                fprintf(['X_index = ',num2str(nx_to_save),'\n']);
                fprintf(['Y_index = ',num2str(ny_to_save),'\n']);
                fprintf(['Repetition_index = ',num2str(nr_to_save),'\n\n']);
                 try
                    %Save the object in the post-processing folder
                    file_number_to_save = obj.meas_number_STR(nx_to_save,ny_to_save,nr_to_save);
                    filename = strcat('plasma_param_stationary_',num2str(file_number_to_save,'%04.f'),'.mat');
                    save([obj.source_folder,'plasma_parameters_stationary_list\',filename],'pp_stationary_to_save');
                    fprintf('Plasma Parameter Stationary saved.\n');
                    obj.save_flag_array(nx_to_save,ny_to_save,nr_to_save) = true;
                    %If the time information has not yet been loaded
                    %NO TIME INFORMATION NEEDED IN THE STATIONARY MEASUREMENT
%                     if ~obj.time_info_loaded
%                         %Take the time information from the first successfully extracted adn saved plasma parameter s
%                         obj.time_axis_rise = pp_timeseries_to_save.ivt_source.total_time_axis_rise;
%                         obj.time_axis_fall = pp_timeseries_to_save.ivt_source.total_time_axis_fall;
%                         obj.Nt = pp_timeseries_to_save.ivt_source.NC_tot;
%                         obj.time_info_loaded = true;
%                     end
                    if all(obj.save_flag_array(nx_to_save,ny_to_save,:))
                        obj.h_rect_ca{ny_to_save,nx_to_save}.FaceColor = 'green';
                        drawnow
                    else
                        obj.h_rect_ca{ny_to_save,nx_to_save}.FaceColor = 'yellow';
                        drawnow
                    end
                 catch err_save
                     fprintf('Plasma Parameter Stationary data analysis failed: jumping to next task...\n');
                     obj.h_rect_ca{ny_to_save,nx_to_save}.FaceColor = 'red';
                     drawnow
                 end 

            end
            
            %Nested function which saves the pp_timeseries with extracted
            %plasma parameters to a file
            function savePptsToFile(args)
                pp_timeseries_to_save = args{1};
                nx_to_save = args{2};
                ny_to_save = args{3};
                nr_to_save = args{4};
                obj.pp_ca{nx_to_save,ny_to_save,nr_to_save} = pp_timeseries_to_save;
                %Diagnostics
                fprintf(['Saving file number = ',num2str(nx_to_save*ny_to_save*nr_to_save-1),'\n']);
                fprintf(['X_index = ',num2str(nx_to_save),'\n']);
                fprintf(['Y_index = ',num2str(ny_to_save),'\n']);
                fprintf(['Repetition_index = ',num2str(nr_to_save),'\n\n']);
                 try
                    %Save the object in the post-processing folder
                    file_number_to_save = obj.meas_number_STR(nx_to_save,ny_to_save,nr_to_save);
                    filename = strcat('plasma_param_timeseries_',num2str(file_number_to_save,'%04.f'),'.mat');
                    save([obj.source_folder,'plasma_parameters_timeseries_list\',filename],'pp_timeseries_to_save');
                    fprintf('Plasma Parameter Timeseries saved.\n');
                    obj.save_flag_array(nx_to_save,ny_to_save,nr_to_save) = true;
                    %If the time information has not yet been loaded
                    if ~obj.time_info_loaded
                        %Take the time information from the first successfully extracted adn saved plasma parameter timeseries
                        obj.time_axis_rise = pp_timeseries_to_save.ivt_source.total_time_axis_rise;
                        obj.time_axis_fall = pp_timeseries_to_save.ivt_source.total_time_axis_fall;
                        obj.Nt = pp_timeseries_to_save.ivt_source.NC_tot;
                        obj.time_info_loaded = true;
                    end
                    if all(obj.save_flag_array(nx_to_save,ny_to_save,:))
                        obj.h_rect_ca{ny_to_save,nx_to_save}.FaceColor = 'green';
                        drawnow
                    else
                        obj.h_rect_ca{ny_to_save,nx_to_save}.FaceColor = 'yellow';
                        drawnow
                    end
                 catch err_save
                     fprintf('Plasma Parameter Timeseries data analysis failed: jumping to next task...\n');
                     obj.h_rect_ca{ny_to_save,nx_to_save}.FaceColor = 'red';
                     drawnow
                 end 

            end
            
            
        end
        
        function obj = buildMDArrays(obj)            
            %Prepare the folder for the plasma parameters timeseries
            mkdir([obj.source_folder,'plasma_parameters_multidim_arrays']);
            
            fprintf('Start of plasma parameters multidimensional array reconstruction.\n');
            obj.T_e_MDa = nan(obj.Nx,obj.Ny,obj.Nr,obj.Nt,2);
            obj.I_sat_i_MDa = nan(obj.Nx,obj.Ny,obj.Nr,obj.Nt,2);
            obj.V_plasma_MDa = nan(obj.Nx,obj.Ny,obj.Nr,obj.Nt,2);
            for nx = 1:obj.Nx
                for ny = 1:obj.Ny
                    for nr = 1:obj.Nr
                        fprintf(['X_index = ',num2str(nx),'\n']);
                        fprintf(['Y_index = ',num2str(ny),'\n']);
                        fprintf(['Repetition_index = ',num2str(nr),'\n\n']);
                        %Access the corresponding pp_timeseries
                        file_number = obj.meas_number_STR(nx,ny,nr);
                        pp_ts_filename = strcat('plasma_param_timeseries_',num2str(file_number,'%04.f'),'.mat');
                        pp_ts_file_path = strcat(obj.source_folder,'plasma_parameters_timeseries_list\',pp_ts_filename);
                        if exist(pp_ts_file_path, 'file')
                            S = load(pp_ts_file_path);
                            try 
                                pp_timeseries = S.pp_timeseries_to_save;
                            catch
                                pp_timeseries = S.pp_timeseries;
                            end
                            %PROVISIONAL FIX-----
                            if isempty(obj.Nt)
                                obj.Nt = pp_timeseries.ivt_source.NC_tot;
                            end
                            %--------------------
                            for nc = 1:obj.Nt
                                for rf = 1:2
                                    obj.T_e_MDa(nx,ny,nr,nc,rf) = pp_timeseries.T_e_a(nc,rf);
                                    obj.I_sat_i_MDa(nx,ny,nr,nc,rf) = pp_timeseries.I_sat_i_a(nc,rf);
                                    obj.V_plasma_MDa(nx,ny,nr,nc,rf) = pp_timeseries.V_plasma_a(nc,rf);
                                end
                            end
                            clear pp_timeseries
                        else
                            fprintf('File not found.\n');
                        end
                    end
                end
            end
            T_e_MDa_filepath = strcat(obj.source_folder,'plasma_parameters_multidim_arrays\','T_e_MDa');
            T_e_MDa = obj.T_e_MDa; %#ok<NASGU,PROP>
            save(T_e_MDa_filepath,'T_e_MDa');
            
            I_sat_i_MDa_filepath = strcat(obj.source_folder,'plasma_parameters_multidim_arrays\','I_sat_i_MDa');
            I_sat_i_MDa = obj.I_sat_i_MDa; %#ok<NASGU,PROP>
            save(I_sat_i_MDa_filepath,'I_sat_i_MDa');
            
            V_plasma_MDa_filepath = strcat(obj.source_folder,'plasma_parameters_multidim_arrays\','V_plasma_MDa');
            V_plasma_MDa = obj.V_plasma_MDa; %#ok<NASGU,PROP>
            save(V_plasma_MDa_filepath,'V_plasma_MDa');
            
            fprintf('Plasma parameters multidimensional array reconstruction has finished.\n');
        end
        
        function obj = buildMDArrays_query(obj)
            
            %Prepare the folder for the plasma parameters timeseries
            mkdir([obj.source_folder,'plasma_parameters_multidim_arrays']);
            
                %Requested parameter names
                query_param_names_column = cell(0);
                query_param_names_column{end+1,1} = 'T_e_estimate';
                query_param_names_column{end+1,1} = 'n_e_estimate';
                query_param_names_column{end+1,1} = 'V_plasma_estimate';
                query_param_names_column{end+1,1} = 'V_float_estimate';
                query_param_names_column = query_param_names_column';
                qpn_length = length(query_param_names_column);
                MD_arrays_column = cell(1,qpn_length);

                obj.MD_arrays_ca = [query_param_names_column;MD_arrays_column]';

                switch obj.st_profile_type
                    case 'stationary'
                        array_size_vector = [obj.Nx,obj.Ny,2];
                    case 'timeseries'
                        array_size_vector = [obj.Nx,obj.Ny,obj.Nt,2];
                end

                for ind = 1:qpn_length

                    obj.MD_arrays_ca{ind,2} = cell(1,obj.Nr);
                    for nr = 1:obj.Nr
                        obj.MD_arrays_ca{ind,2}{1,nr} = nan(array_size_vector);
                    end

                        switch obj.st_profile_type
                            case 'stationary'
                                for nx = 1:obj.Nx
                                    for ny = 1:obj.Ny
                                        for nr = 1:obj.Nr
                                            for nrf = 1:2
                                                if obj.Nr > 1    
                                                    var_str = ['obj.pp_ca{',num2str(nx),',',num2str(ny),',',num2str(nr),'}.fit_iv_curve_out_ca{1,',num2str(nrf),'}.',obj.MD_arrays_ca{ind,1}];
                                                else
                                                    var_str = ['obj.pp_ca{',num2str(nx),',',num2str(ny),'}.fit_iv_curve_out_ca{1,',num2str(nrf),'}.',obj.MD_arrays_ca{ind,1}];
                                                end
                                                try
                                                    obj.MD_arrays_ca{ind,2}{1,nr}(nx,ny,nrf) = eval(var_str);
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
                                                        var_str = ['obj.pp_ca{',num2str(nx),',',num2str(ny),',',num2str(nt),',',num2str(nr),'}.fit_iv_curve_out_ca{1,',num2str(nrf),'}.',obj.MD_arrays_ca{ind,1}];
                                                    else
                                                        var_str = ['obj.pp_ca{',num2str(nx),',',num2str(ny),'}.fit_iv_curve_out_ca{1,',num2str(nrf),'}.',obj.MD_arrays_ca{ind,1}];
                                                    end
                                                    try
                                                        obj.MD_arrays_ca{ind,2}{1,nr}(nx,ny,nt,nrf) = eval(var_str);
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
        end
        
        function obj = loadPpTimeseriesCaFromFileList(obj)

            fprintf('Start of loadPpTimeserieCaFromFileList.\n');
            for nx = 1:obj.Nx
                for ny = 1:obj.Ny
                    for nr = 1:obj.Nr
                        fprintf(['X_index = ',num2str(nx),'\n']);
                        fprintf(['Y_index = ',num2str(ny),'\n']);
                        fprintf(['Repetition_index = ',num2str(nr),'\n\n']);
                        %Access the corresponding pp_timeseries
                        file_number = obj.meas_number_STR(nx,ny,nr);
                        pp_ts_filename = strcat('plasma_param_timeseries_',num2str(file_number,'%04.f'),'.mat');
                        pp_ts_file_path = strcat(obj.source_folder,'plasma_parameters_timeseries_list\',pp_ts_filename);
                        if exist(pp_ts_file_path, 'file')
                            load(pp_ts_file_path,'pp_timeseries');
                            obj.pp_timeseries_ca{nx,ny,nr} = pp_timeseries;
                        else
                            fprintf('File not found.\n');
                        end
                    end
                end
            end
        end
        
        function obj = loadMDArraysFromFiles(obj)
            load(strcat(obj.source_folder,'plasma_parameters_multidim_arrays\','T_e_MDa.mat'),'T_e_MDa');
            obj.T_e_MDa = T_e_MDa;
            load(strcat(obj.source_folder,'plasma_parameters_multidim_arrays\','I_sat_i_MDa.mat'),'I_sat_i_MDa');
            obj.I_sat_i_MDa = I_sat_i_MDa;
            load(strcat(obj.source_folder,'plasma_parameters_multidim_arrays\','V_plasma_MDa.mat'),'V_plasma_MDa');
            obj.V_plasma_MDa = V_plasma_MDa;
        end
        
        function pp_timeseries = getPpTimeseries(obj,x_pos,y_pos,nr)
            
            cond_x = (x_pos >= obj.x_begin && x_pos <= obj.x_end);
            cond_y = (y_pos >= obj.y_begin && y_pos <= obj.y_end);
            cond_rep = (nr >= 1 && nr <= obj.Nr);
            
            if(cond_x && cond_y && cond_rep)
                [~,nx] = min(abs(obj.x_axis - x_pos));
                [~,ny] = min(abs(obj.y_axis - y_pos));
                pp_timeseries = obj.pp_timeseries_ca{nx,ny,nr};
                return
            else
                fprintf('Error: requested position or repetition index is out of range.\n')
            end
        end
        
    end
    
end

