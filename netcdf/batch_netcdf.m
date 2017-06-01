function batch_netcdf
% Начало работы 26.05.2017.
% 1. Проверка правильности файлов формата Кулешова в указанной папке (по размеру).
% 2. Определение промежутка между данными в отдельных файлах.
% 3. Если промежуток больше 10 с. (время неопределенности формата),
% записывается новый временной массив.

dt = 0.0002505;     % Время между отсчётами в секундах

clc
workdir = uigetdir('d:\!data\kuleshov\', 'Проверяемая папка...');

% Создание netCDF файла.
try
    cmode = bitor(netcdf.getConstant('CLOBBER'), netcdf.getConstant('NETCDF4'));
    ncid = netcdf.create([workdir, '\directory.nc'], cmode);

    % Create an global attribute.
    netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'dt',dt);
    
    dirlist = dir([workdir, '\T*.*']);
    s_in_day = 24 * 3600;
    
    % Перебор по файлам в рабочей папке и внесение данных в хранилище.
    prevEnd = 0;
    prevFile = '';
    key_new_group = true;
    for i = 1:length(dirlist)
        unit = dirlist(i);
        if unit.bytes == 0 % пропускаем файлы нулевого размера
            continue
        end
        
        begTime = getTimeFromFName(unit.name);
        dT = dt * (unit.bytes/2 - 1) / s_in_day; % длительность файла в сутках
        endTime = begTime + dT;
        
        if prevEnd > 0 % определение разрыва оси времени
            cur_delay = (begTime - prevEnd) * s_in_day;
            if abs(cur_delay) > 10
                disp([num2str(i),': ',prevFile,' -> ',unit.name,', ', num2str(cur_delay/60), ' мин, последовательность разорвана!'])
                key_new_group = true;
            end
        end
        
        % создание новой группы при наличие разрыва времени
        if key_new_group
            group_name = datestr(begTime,'yyyy-mm-dd_HH-MM-SS');
            childGroupId = netcdf.defGrp(ncid,group_name);
            % define dimensions
            time_dimid = netcdf.defDim(childGroupId,'time',netcdf.getConstant('NC_UNLIMITED'));
            % define variables
            varid = netcdf.defVar(childGroupId,'amplitude','short',time_dimid);
            netcdf.defVarDeflate(childGroupId,varid,true,true,9);
            % Leave define mode and enter data mode to write data.
            netcdf.endDef(childGroupId);
            
            key_new_group = false;
            start = 0;
        end
        
        % Извлечём данные из текущего файла
        fullname = [workdir, '\', unit.name];
        fid = fopen(fullname,'rb');
        if fid
            ampl = fread(fid,inf,'int16=>int16');
            fclose(fid);
        end
        count = length(ampl);
        netcdf.putVar(childGroupId,varid,start,count,ampl)
        start = start + count;
        
        prevEnd = endTime;
        prevFile = unit.name;
        
    end
    
catch ME
    warning(ME.message)
end
netcdf.close(ncid)
ncdisp([workdir, '\directory.nc'])

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