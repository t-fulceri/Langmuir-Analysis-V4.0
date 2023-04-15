classdef VSignal < handle
    %Contains the V_float signal timetrace
    
    properties
        init = false;
        daq_settings = DaqSettings;
        date
        time
        time_step
        f_sampling
        signal_length
        time_axis
        V_signal
        V_mean
        V_std
        V_ADC_step
        
    end
    
    methods
        
        function obj = VSignal(varargin)
            %The class constructor accepts the following optional arguments:
            %varargin{1} = full_file_path
            
            %If the file path information is present
            if nargin >= 1
                %Retrieve the file path information
                full_file_path = varargin{1};
                %Initialize the object with the data read from the file
                assignPropertiesFromHD5File( obj, full_file_path);
                %If arrived here, set the object as initialized
                obj.init = true;
            else
                %Initialize all the properties to default values
                obj.date = 0;
                obj.time = 0;
                obj.time_step = 0;
                obj.f_sampling = 0;
                obj.signal_length = 0;
                obj.time_axis = 0;
                obj.V_signal = 0;
                obj.V_ADC_step = 0;
            end

        end
        
        function obj = assignPropertiesFromHD5File( obj, full_file_path)
            
            assignDAQ_settingsFromFile( obj, full_file_path);
            
            [~,name,ext] = fileparts(full_file_path);
            fprintf(['assignPropertiesFromHD5File\n',name,ext,'\n']);
            fprintf('Settings:\n');
            fprintf(['voltage_sign = ',num2str(obj.daq_settings.v_sign),'\n']);
            fprintf(['current_sign = ',num2str(obj.daq_settings.i_sign),'\n']);
            fprintf(['voltage_deamplification = ',num2str(obj.daq_settings.v_deampl),'\n\n']);
            fprintf(['current_divider = ',num2str(obj.daq_settings.i_divider),'\n\n']);

            %Collect all the relevant information from the HD5 file
            %Check metadata
            metadata = h5info(full_file_path);
            %Extract timestep from metadata [s]
            obj.time_step = metadata.Groups.Attributes(2,1).Value(1);
            %Extract date and time from metadata
            obj.date = metadata.Datasets.Attributes(1,1).Value;
            obj.time = metadata.Datasets.Attributes(2,1).Value;

            %Read data
            data = h5read(full_file_path,'/dataset');
            %Signal length
            obj.signal_length = length(data);
            %Time axis
            obj.time_axis = squeeze(obj.time_step.*(0:obj.signal_length-1))';
            %Extract voltage
            %NOTE: PROBE VOLTAGE IS IN CHANNEL 1 OF THE HD5 FILE PRODUCED
            %BY THE OSCILLOSCOPE
            voltage = squeeze(data(1,:))';
            %Optional inversion of voltage signal (depends on the daq settings)
            obj.V_signal = double(obj.daq_settings.v_deampl).*double(obj.daq_settings.v_sign).*voltage;
            %Extract uncertainty coming from ADC discretization;
            obj.V_ADC_step = signalToADCStep(voltage);
            %Mean of the signal
            obj.V_mean = mean(obj.V_signal);
            %Standard deviation of the signal
            obj.V_std = std(obj.V_signal);
            %Calculate sampling frequency [Hz]
            obj.f_sampling = 1/obj.time_step;
            
        end
        
        function obj = assignDAQ_settingsFromFile( obj, full_file_path)
            
            [source_folder,~,~] = fileparts(full_file_path);
            daq_settings_path = strcat(source_folder,'\daq_settings.txt');
            obj.daq_settings.readDaqSettingsFromFile(daq_settings_path);
        end
        
    end
    
end

