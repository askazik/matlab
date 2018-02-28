classdef classAmplitude < memmapfile
% classAmplitude - ������ � �������� ����� ������ ��������.
% ����������� �������� �������� � ��������� ����� ������, � ������� 
% ����������� ������� ������������ �������� ���������� �� ��������� 
% ������� ��� ���������� ������ �����������.
% ��� ���������� ������ ������ ������������ ���������� ����������� ������ 
% � ������ (Memory Mapped Files).
% ��� �������� ������������ �������� ���������� ������������ ��� ������
% inputParser class.
%
% See also: <a href="matlab:doc memmapfile">memmapfile class</a>, 
%<a href="matlab:doc inputParser">inputParser class</a>.
%
% Created by Alexy Skazik (c)                                   Feb-2018
    
    properties (GetAccess = private, SetAccess = private)
        fileProperties
    end
        
    properties (GetAccess = public, SetAccess = private)
        version
        frqs
        heights
        date
        dt
    end
    properties
        position
    end
    
    methods
        %% Class constructor
        function this = classAmplitude(fileName, varargin)
            % ����������� ������ ������� � ����� ���������.
            
            % Verify input arguments
            p = inputParser;
            isfrq = @(str) exist('str','var') & ischar(str) & strcmp(str(end-3:end),'.frq'); 
            p.addRequired('fileName', isfrq);
            p.parse(fileName, varargin{:});
            
            % Object Initialization
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            this = this@memmapfile(p.Results.fileName);
         
            % Post Initialization
            % Any code, including access to object
            this.fileProperties = this.getHeader();
            this.version = this.fileProperties.ver;
            this.frqs = this.fileProperties.frqs;
            this.heights = this.getRealHeights();
            this.date = datenum(this.fileProperties.year + 1900, ...
                                this.fileProperties.mon, ...
                                this.fileProperties.mday, ...
                                this.fileProperties.hour, ...
                                this.fileProperties.min, ...
                                this.fileProperties.sec);
            this.dt = length(this.frqs) * 1. / this.fileProperties.pulse_frq;
            
            % ��������� ��������� � ������.
            this.offset = this.fileProperties.data_begin_position;
            % ������ ����� ������.
            myFormat = cell(length(this.frqs),3);
            j = 1;
            for i = 1:2:2*length(this.frqs)
                % ���� ������ ���������
                myFormat{i,1} = 'int16';
                myFormat{i,2} = [this.fileProperties.count_height,2];
                myFormat{i,3} = ['frq_', int2str(j)];
                % ���������� �� �������� ��������� (�������� � �����)
                myFormat{i+1,1} = 'int16';
                myFormat{i+1,2} = [1,1];
                myFormat{i+1,3} = ['gain_', int2str(j)];
                
                j = j + 1;
            end
            this.format = myFormat;
        end
        
        %% Disp overload
        function disp(obj)
            disp@memmapfile(obj);

            str{1}=sprintf('Version: (%d)', obj.version);
            str{2}=sprintf('Frequensies, kHz: (%s)', sprintf('%d, ',obj.frqs));
            str{3}=sprintf('Heights, km: (%.3f : %.3f : %.3f)', ...
                obj.heights(1), ...
                obj.heights(2)-obj.heights(1), ...
                obj.heights(end));
            str{4}=sprintf('Date: (%s) UT', datestr(obj.date));
            str{5}=sprintf('dt: (%.3f), s', obj.dt);
            str = sprintf('%s\n',str{:});
            fprintf(str)
        end

        %% Get version of the data file.
        function value = get.version(this)
            % �������� ������ ����� ������.
            value = this.version;
        end
        
        %% Get frequensies of the data file.
        function value = get.frqs(this)
            % �������� ������� ������������.
            value = this.frqs;
        end
        
        %% Get heights of the data file.
        function value = get.heights(this)
            % �������� ������ ������������.
            value = this.heights;
        end
        
        %% Get date of the data file.
        function value = get.date(this)
            % �������� ���� ������������ � �������� �������.
            value = this.date;
        end
        
        %% Get time step of the data file.
        function value = get.dt(this)
            % �������� ����� ����� ���������� ������������.
            value = this.dt;
        end
                      
        %% Get data from data unit.
        function out = getUnit(this, i)
            % �������� ���� ������ �� ����� ���������� ���������.
            %d = this.Data;
            d = get(this, 'Data');
            cur = d(i);
            fnames = fieldnames(cur);
            
            j = 1;
            for i = 1:2:length(fnames)
                tmp = double(cur.(fnames{i}));
                % ��������� ����������.
                out.amplitudes(:,j) = complex(tmp(:,1),tmp(:,2));
                out.gain(:,j) = cur.(fnames{i+1});
                
                j = j + 1;
            end
        end
    end
    
    methods (Access = private, Sealed)
        %% Get real heights file version.
        function heights = getRealHeights(this)
            % ����������� ����������� ����� ��������� �������� ������ �����
            % ������.
            
            heights = this.fileProperties.height_min : ...
                        this.fileProperties.height_step / 1000. : ...
                        this.fileProperties.height_max;
            switch this.version
                case 0
                    % �������� �����
                case 1
                    % ��� ��������� ������
                    heights = (0 : ...
                        this.fileProperties.height_step : ...
                        (this.fileProperties.count_height - 1) * this.fileProperties.height_step) / 1000.;
                case 2
                    % �������� �����
                otherwise
                    % �������� �����
            end
            
        end
        
        %% Read the file header.
        function properties = getHeader(this)
            % ������ ��������� ����� ������.
            %
            % struct dataHeader { 	    // === ��������� ����� ������ ===
            %   unsigned ver; // ����� ������
            %   struct tm time_sound; // GMT ����� ��������� ������������
            %   unsigned height_min; // ��������� ������, �� (��, ��� ���� ��� ��������� �������������)
            % 	unsigned height_max; // �������� ������, �� (��, ��� ���� ��� ��������� �������������)
            %   unsigned height_step; // ��� �� ������, � (�������� ���, ����������� �� ������� ���)
            %   unsigned count_height; // ����� ����� (������ ��������� ������ ��� ��� ������������, fifo ������ ��� 4��. �.�. �� ������ 1024 �������� ��� ���� ������������ �������)
            %   unsigned count_modules; // ���������� �������/������ ������������
            % 	unsigned pulse_frq; // ������� ����������� ���������, ��
            % };
            %
            % struct tm
            % Member	Type	Meaning	Range
            % tm_sec	int	seconds after the minute	0-60*
            % tm_min	int	minutes after the hour	0-59
            % tm_hour	int	hours since midnight	0-23
            % tm_mday	int	day of the month	1-31
            % tm_mon	int	months since January	0-11
            % tm_year	int	years since 1900
            % tm_wday	int	days since Sunday	0-6
            % tm_yday	int	days since January 1	0-365
            % tm_isdst	int	Daylight Saving Time flag
            % =========================================================================
            
            % 1. ��������� ���� �� ������.
            [fid, message] = fopen(this.filename, 'r');
            if ( fid == -1)
                error(message);
            end
            
            % 2. ��������� ��������� ������������.
            properties.ver = fread(fid,1,'uint32');
            properties.sec = fread(fid,1,'int32');
            properties.min = fread(fid,1,'int32');
            properties.hour = fread(fid,1,'int32');
            properties.mday = fread(fid,1,'int32');
            properties.mon = fread(fid,1,'int32');
            properties.year = fread(fid,1,'int32');
            properties.wday = fread(fid,1,'int32');
            properties.yday = fread(fid,1,'int32');
            properties.isdst = fread(fid,1,'int32');
            properties.height_min = fread(fid,1,'uint32');
            properties.height_max = fread(fid,1,'uint32');
            properties.height_step = fread(fid,1,'uint32');
            properties.count_height = fread(fid,1,'uint32');
            properties.count_modules = fread(fid,1,'uint32');
            properties.pulse_frq = fread(fid,1,'uint32');
            
            % ���������������� ������ ������ ������������
            for i=1:properties.count_modules
                properties.frqs(i) = fread(fid,1,'uint32');
            end
            
            properties.data_begin_position = ftell(fid);
            
            % 3. ��������� ���� ������.
            if ~isempty(fid)
                fclose(fid);
            end
        end
        
    end
end

