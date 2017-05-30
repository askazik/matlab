function batch_verify_time_periods
% Начало работы 26.05.2017.
% 1. Проверка правильности файлов формата Кулешова в указанной папке (по размеру).
% 2. Определение промежутка между данными в отдельных файлах (В отсчётах. Нужно заполнять нулями.)

clc

dt = 0.0002505;     % Время между отсчётами в секундах
workdir = uigetdir('d:/!data/kuleshoff', 'Проверяемая папка...');

dirlist = dir([workdir, '\T*.*']);
s_in_day = 24 * 3600;
N = 2621440;

prevEnd = 0;
prevFile = '';
for i = 1:length(dirlist)
    unit = dirlist(i);
    fullname = [workdir, '\', unit.name];
       
    if unit.bytes/2 < N % Сохранён полный файл данных
        disp(['Файл <',unit.name,'> (',num2str(unit.bytes),' bytes) - меньшего размера!'])
    end
    
    begTime = getTimeFromFName(unit.name);
    dT = dt * (unit.bytes/2 - 1) / s_in_day; % длительность файла в сутках
    endTime = begTime + dT;
    
    if prevEnd > 0
        cur_delay = (begTime - prevEnd) * s_in_day;
        if abs(cur_delay) > 10
            disp([num2str(i),': ',prevFile,' -> ',unit.name,', ', num2str(cur_delay/60), ' мин, последовательность разорвана!'])
        end
    end
    prevEnd = endTime;
    prevFile = unit.name;
end

function begTime = getTimeFromFName(fname)
% Извлечение даты выборки из имени файла

% Год
year = 2000 + str2num(fname(2:3));
% Месяц
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
% Выделим десятки секунд
sec = 6 * str2num(fname(12));
                    
begTime = datenum(year, month, day, hour, min, sec);
         