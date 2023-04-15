classdef DaqSettings < handle
    %DAQSETTINGS Contains the information about the data acquisition sytem
    %used to record the V and I signals from the oscilloscope
    
    properties
        v_sign
        i_sign
        v_deampl
        i_divider
    end
    
    methods
        function obj = DaqSettings(varargin)
            %The class constructor accepts the following optional arguments:
            %varargin{1} = v_sign
            %varargin{2} = i_sign
            %varargin{3} = v_deampl
%             obj.v_sign = 0;
%             obj.i_sign = 0;
%             obj.v_deampl = 0;
%             if nargin >= 1
%                 v_sign = varargin{1};
%                 obj.v_sign = v_sign;
%             else
%                 obj.v_sign = 0;
%             end
%             if nargin >= 2
%                 i_sign = varargin{2};
%                 obj.i_sign = i_sign;
%             else
%                 obj.i_sign = 0;
%             end
%             if nargin >= 3
%                 v_deampl = varargin{3};
%                 obj.v_deampl = v_deampl;
%             else
%                 obj.v_deampl = 0;
%             end
        end
        
        function obj = readDaqSettingsFromFile(obj,full_file_path)
            
            fileID = fopen(full_file_path,'r');
            C = textscan(fileID,'%s %5.3f');
            fclose(fileID);
            values = C{:,2};
            
            obj.v_sign = values(1);
            obj.i_sign = values(2);
            obj.v_deampl = values(3);
            obj.i_divider = values(4);

        end
    end
    
end

