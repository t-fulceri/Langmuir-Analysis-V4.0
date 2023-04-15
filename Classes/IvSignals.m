classdef IvSignals < handle
    %Contains the voltage/current signals with minimal metadata
    
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
        I_signal
        I_smooth
        I_smooth_exp_fit
        I_CM_corrected
        V_ADC_step
        I_ADC_step
        
        
    end
    
    methods
        
        function obj = IvSignals(varargin)
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
                obj.I_signal = 0;
                obj.I_smooth = 0;
                obj.I_smooth_exp_fit = 0;
                obj.I_CM_corrected = 0;
                obj.V_ADC_step = 0;
                obj.I_ADC_step = 0;
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
            voltage = squeeze(data(1,:))';
            %Optional inversion of voltage signal (depends on the daq settings)
            obj.V_signal = double(obj.daq_settings.v_deampl).*double(obj.daq_settings.v_sign).*voltage;
            %Extract current
            current = squeeze(data(2,:))';
            %Optional inversion of current signal (depends on the daq settings)
            obj.I_signal = (1/double(obj.daq_settings.i_divider)).*double(obj.daq_settings.i_sign).*current;

            %Extract uncertainty coming from ADC discretization;
            obj.V_ADC_step = signalToADCStep(voltage);
            obj.I_ADC_step = signalToADCStep(current);

            %Calculate sampling frequency [Hz]
            obj.f_sampling = 1/obj.time_step;
            
        end
        
        function obj = assignDAQ_settingsFromFile( obj, full_file_path)
            
            [source_folder,~,~] = fileparts(full_file_path);
            daq_settings_path = strcat(source_folder,'\daq_settings.txt');
            obj.daq_settings.readDaqSettingsFromFile(daq_settings_path);
        end
        
        function obj = calculateCMEffect( obj )
            
            %CM time constant (extracted from specific measurements)
            tau = 1.4e-3;
            smooth_param = 1000;
            
            obj.I_smooth = smooth(obj.I_signal,smooth_param,'sgolay');
            I_max = max(obj.I_smooth);

            [xData, yData] = prepareCurveData( obj.time_axis, obj.I_smooth );
            % Set up fittype and options.
            ft = fittype( 'A*exp(-t/1.4e-3)', 'independent', 't', 'dependent', 'I' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = I_max;

            % Fit model to data.
            [fitresult, ~] = fit( xData, yData, ft, opts );
            A = fitresult.A;
            obj.I_smooth_exp_fit = A.*exp(-obj.time_axis./tau);

            obj.I_CM_corrected = obj.I_signal + A - obj.I_smooth_exp_fit;
        end
        
        function I_vacuum_effect = extractVacuumEffect( obj )
            %The vacuum effect is a constant offset of about 1 mA plus an
            %a current oscillating at the same frequency as the drive due
            %to a hysteresis effect
            %ONLY USE THIS METHOD ID THE SIGNAL IS COMING FROM A VACUUM
            %MEASUREMENT
            I_mean = mean(obj.I_signal);
            [ main_freq, main_phase ] = extractMainFreqAndPhase( obj.I_signal - I_mean, obj.f_sampling ,false);
            
            [xData, yData] = prepareCurveData( obj.time_axis, obj.I_signal - I_mean );
            ft = fittype( 'A*cos(2*pi*f*t + phi)', 'independent', 't', 'dependent', 'I_osc' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [0.001 main_freq main_phase];
            % Fit model to data.
            [fitresult, ~] = fit( xData, yData, ft, opts );
            
            I_vacuum_effect = I_mean + (fitresult.A).*cos((2*pi*fitresult.f).*obj.time_axis + fitresult.phi);
        end
        
    end
    
end

