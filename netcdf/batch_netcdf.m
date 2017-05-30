function batch_netcdf
% Начало работы 26.05.2017.
% 1. Проверка правильности файлов формата Кулешова в указанной папке (по размеру).
% 2. Определение промежутка между данными в отдельных файлах.
% 3. Если промежуток больше 10 с. (время неопределенности формата),
% записывается новый временной массив.

clc

dt = 0.0002505;     % Время между отсчётами в секундах
workdir = uigetdir('d:\!data\kuleshov\', 'Проверяемая папка...');

dirlist = dir([workdir, '\T*.*']);
s_in_day = 24 * 3600;
N = 2621440;

% Создание netCDF файла.
nc_file = [workdir, '\data.nc'];
% nccreate(nc_file,'v_amplitude',...
%           'Dimensions',{'d_time',inf},...
%           'Datatype','int16',...
%           'DeflateLevel',9,...
%           'Format','netcdf4_classic');
% ncdisp(nc_file);

% cmode = bitor(netcdf.getConstant('CLOBBER'), ...
%     netcdf.getConstant('64BIT_OFFSET'));
cmode = bitor(netcdf.getConstant('CLOBBER'), netcdf.getConstant('NETCDF4'));
%cmode = bitor(cmode, netcdf.getConstant('CLASSIC_MODEL'));
ncid = netcdf.create(nc_file, cmode);

% define dimensions
time_dimid = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));

% define variables
varid = netcdf.defVar(ncid,'amplitude','NC_SHORT',[time_dimid]);
netcdf.defVarDeflate(ncid,varid,true,true,9);

% Create an attribute associated with the variable.
netcdf.putAtt(ncid,0,'dt',dt);

netcdf.endDef(ncid);
ncdisp(nc_file)

% Перебор по файлам в рабочей папке и внесение данных в хранилище.
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
netcdf.close(ncid)

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
         