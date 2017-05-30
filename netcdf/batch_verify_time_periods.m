function batch_verify_time_periods
% ������ ������ 26.05.2017.
% 1. �������� ������������ ������ ������� �������� � ��������� ����� (�� �������).
% 2. ����������� ���������� ����� ������� � ��������� ������ (� ��������. ����� ��������� ������.)

clc

dt = 0.0002505;     % ����� ����� ��������� � ��������
workdir = uigetdir('d:/!data/kuleshoff', '����������� �����...');

dirlist = dir([workdir, '\T*.*']);
s_in_day = 24 * 3600;
N = 2621440;

prevEnd = 0;
prevFile = '';
for i = 1:length(dirlist)
    unit = dirlist(i);
    fullname = [workdir, '\', unit.name];
       
    if unit.bytes/2 < N % ������� ������ ���� ������
        disp(['���� <',unit.name,'> (',num2str(unit.bytes),' bytes) - �������� �������!'])
    end
    
    begTime = getTimeFromFName(unit.name);
    dT = dt * (unit.bytes/2 - 1) / s_in_day; % ������������ ����� � ������
    endTime = begTime + dT;
    
    if prevEnd > 0
        cur_delay = (begTime - prevEnd) * s_in_day;
        if abs(cur_delay) > 10
            disp([num2str(i),': ',prevFile,' -> ',unit.name,', ', num2str(cur_delay/60), ' ���, ������������������ ���������!'])
        end
    end
    prevEnd = endTime;
    prevFile = unit.name;
end

function begTime = getTimeFromFName(fname)
% ���������� ���� ������� �� ����� �����

% ���
year = 2000 + str2num(fname(2:3));
% �����
caseSwitch = fname(4);
month = 0;
switch (caseSwitch)
    case 'A'
        month = 10;
    case 'B'
        month = 11;
    case 'C'
        month = 12;
    otherwise
        month = str2num(caseSwitch);
end

day = str2num(fname(5:6));
hour = str2num(fname(7:8));
min = str2num(fname(10:11));
% ������� ������� ������
sec = 6 * str2num(fname(12));
                    
begTime = datenum(year, month, day, hour, min, sec);
         