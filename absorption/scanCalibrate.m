function scanCalibrate

step_accumulate = 1; % ��� ���������� ������
fname = 'd:\!data\2018-02-13\20180216083227.frq';
h = figure('NumberTitle','off');
out = getData(fname, 0);

if out.properties.ver == 1
    i = 0;
    while 1
        
        a = zeros(out.properties.count_height,5);
        b = a;
        c = a;
        for i = i:(i+step_accumulate)
            out = getData(fname, i); % ����� �����������, ��� ���������� �������� ��������!!!            
            a = a + abs(real(out.data));
            b = b + abs(imag(out.data));
            c = c + abs(out.data);
        end
        frqs = out.properties.frqs';
        
        % f = ', num2str(out.properties.frqs(frq_num)),' kHz'
        tmp = ['time offset = ',num2str(out.time),' s, Ku = ',strcat(int2str(out.gain))];
        set(h,'Name',tmp)

        subplot(1,3,1);
        plot(a/step_accumulate,out.h);
        legend(int2str(frqs))
        grid on
        xlim([-8500 8500])

        subplot(1,3,2);
        plot(b/step_accumulate,out.h);
        legend(int2str(frqs))
        grid on
        xlim([-8500 8500])

        subplot(1,3,3);
        plot(c/step_accumulate,out.h);
        legend(int2str(frqs))
        grid on
        xlim([0 12000])
        
        pause(0.01);
    end
else
    close(h)
    error('scanCalibrate:UnsupportedFileVersion', ...
    'Specified version of the ionogramm file %d is not supported.', ...
    out.properties.ver)
end

function out = getData(fname, idx_position)
%==========================================================================
% ������� ��������� ���� ������ �� ����� ������������ ������������ ��
% ���������� ��������. (�������� "�����".)
% �������������� ������������ ��� ����������� ����� ��������� �- � �-
% ��������� �� ����� �� ��������� � ���������� ��������������� ��
% ����������� �������� �������.
% === ����: ===============================================================
% fname - ������ ���� � ����� ������
% idx_position - ����� ��������� ������ ����� ��� ��������� �������
% ������������ (�� 0)
%==========================================================================

% �������� ���������
%fname = '20161103070242.frq';
%idx_position = 100;

% 1. �������� ���� �� �������������.
if ~ischar(fname)
    error('��� ����� - ������ ������!'); 
end
if ~exist(fname, 'file') == 2
   error('��������� ���� �� ���������� ��� ��������� � ������ �����.'); 
end

% 2. ��������� ���� �� ������ � ��������� ��������� ������������.
[fid, message] = fopen(fname, 'r');
if ( fid == -1)
    error(message); 
end
properties = getHeader(fid);

out.properties = properties;
out.time = idx_position * ( properties.count_modules / double(properties.pulse_frq)); % � ��������, � ������ ������������
height_max = (properties.height_step/1000.) * (properties.count_height-1);
out.h = (0 : properties.height_step/1000. : height_max)'; % ������ ������� ���������

% 3. ��������� ���� ������, ��������� �������� ���������.
tmp = getDataUnit(properties, idx_position);
out.data = tmp.data;
out.gain = tmp.gain;

% 4. ��������� ���� ������.
if ~isempty(fid)
    fclose(fid);
end

function properties = getHeader(fid)
% =========================================================================
% ������ ���������.
% struct dataHeader { 	    // === ��������� ����� ������ ===
%   unsigned ver; // ����� ������
%   struct tm time_sound; // GMT ����� ��������� ������������
%   unsigned height_min; // ��������� ������, � (��, ��� ���� ��� ��������� �������������)
% 	unsigned height_max; // �������� ������, � (��, ��� ���� ��� ��������� �������������)
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
properties.ver = fread(fid,1,'uint32');

properties.sec = fread(fid,1,'int32');
properties.min = fread(fid,1,'int32');
properties.hour = fread(fid,1,'int32');
properties.mday = fread(fid,1,'int32');
properties.mon = fread(fid,1,'int32');
properties.year = fread(fid,1,'int32')+1900;
properties.wday = fread(fid,1,'int32');
properties.yday = fread(fid,1,'int32');
properties.isdst = fread(fid,1,'int32');

properties.height_min = fread(fid,1,'uint32');
properties.height_max = fread(fid,1,'uint32');
properties.height_step = fread(fid,1,'uint32');
properties.count_height = fread(fid,1,'uint32');
properties.count_modules = fread(fid,1,'uint32');
properties.pulse_frq = fread(fid,1,'uint32');

for i=1:properties.count_modules % ���������������� ������ ������ ������������
    properties.frqs(i) = fread(fid,1,'uint32');
end

properties.data_begin_position = ftell(fid);
properties.fid = fid;
    
function out = getDataUnit(properties, idx_position)
% =========================================================================
% ������ ������ ��� ���� ��������� ������ � ���� ������ �������.
% properties - ��������� ����� ���������
% heights_count - ����� ����� � ����� ���������
% idx_position - ����� ����� (������� � 0)
% =========================================================================
frq_unit_size = properties.count_height * 4 + 2;
full_unit_size = properties.count_modules * frq_unit_size;
offset = properties.data_begin_position + idx_position * full_unit_size; % � ������
fseek(properties.fid, offset, 'bof');

gain = zeros(1,properties.count_modules);
data = zeros(properties.count_height,properties.count_modules);
for i = 1:properties.count_modules
    dirty_data = fread(properties.fid, frq_unit_size/2-1,'int16'); % ������ ����������� � �������
    gain(i) = fread(properties.fid, 1,'int16'); % �������� ��������� ��� ���������, � ��

    % ��������� ����������.
    a = dirty_data(1:2:end) * 1.0;
    b = dirty_data(2:2:end) * 1.0;
    
    data(:,i) = complex(a,b);
end

out.data = data;
out.gain = gain;