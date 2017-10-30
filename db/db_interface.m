function varargout = db_interface(varargin)
% 
% DB_INTERFACE M-file for db_interface.fig
%      DB_INTERFACE, by itself, creates a new DB_INTERFACE or raises the
%      existing
%      singleton*.
%
%      H = DB_INTERFACE returns the handle to a new DB_INTERFACE or the handle to
%      the existing singleton*.
%
%      DB_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in DB_INTERFACE.M with the given input
%      arguments.
%
%      DB_INTERFACE('Property','Value',...) creates a new DB_INTERFACE or
%      raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before db_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to db_interface_OpeningFcn via
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help db_interface

% Last Modified by GUIDE v2.5 23-Oct-2017 13:26:18

%% Инициализации интерфейса

% Begin initialization code - DO NOT EDIT
global dbName user password
global connection data data_proc
global last_saved_mat % последний сохраненный файл MAT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @db_interface_OpeningFcn, ...
                   'gui_OutputFcn',  @db_interface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before db_interface is made visible.
function db_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to db_interface (see VARARGIN)

% Choose default command line output for db_interface
handles.output = hObject;
handles.keyButtonDown = false;
handles.faceColor = 'r';
handles.faceMarker = 'o';

% Update handles structure
guidata(hObject, handles);

OpenMenuItem_Callback(hObject, eventdata, handles)

% UIWAIT makes db_interface wait for user response (see UIRESUME)
% uiwait(handles.db_interface);


% --- Outputs from this function are returned to the command line.
function varargout = db_interface_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
%
%% Начало блока рукописных функций
% --------------------------------------------------------------------

function cur = getRowByKey(key)
% Извлекает строку данных из таблицы parus.
% key - значение ключа iono_key
global connection

if ~isempty(connection)
% SQL query to get all fields from Table.
    query = ['SELECT * FROM parus WHERE iono_key = ',num2str(key)];
    data = myblob_command (connection, query);    
end %if

% Сохраняем данные в структуру
cur = struct(...
    'time_sound',data{2},...
    'station_id',data{3},...
    'latitude',data{4},...
    'longitude',data{5},...
    'height_min',data{11},...
    'height_max',data{12},...
    'height_step',data{13},...
    'freq_min',data{14},...
    'freq_max',data{15},...
    'count_freq',data{16},...
    'signal_data',data{22}...
);

function plotRow(key)
% Отрисовывает ионограмму по заданной строке данных из таблицы parus.
% key - значение ключа iono_key

cur = getRowByKey(key);
data = cur.signal_data;
n = length(data);

if ~strcmp(data,'null')
i = 1; % устанавливаем индекс на начало
frq_i = 1;
h_n = 1 + (cur.height_max - cur.height_min)/cur.height_step;
mat_o = nan(h_n, cur.count_freq);
mat_x = nan(h_n, cur.count_freq);
while i < n
    FrequencyData = readFrequencyData(data,i);
    i = i + 13;
    
    % Чтение откликов обыкновенных волн
    for i_o = 1:FrequencyData.count_o  
        SignalResponse = readSignalResponse(data,i);
        i = i + 6;
        
        if i + SignalResponse.count_samples-1 > n % если блок битый?
            break;
        end
        arr_o = FrequencyData.thereshold_o + typecast(data(i:i+SignalResponse.count_samples-1), 'uint8');
        % Определим номер начальной высоты
        j = 1 + (SignalResponse.height_begin - uint32(cur.height_min))/uint32(cur.height_step);
        mat_o(j:j+uint32(SignalResponse.count_samples)-1,frq_i) = arr_o;
        
        i = i + SignalResponse.count_samples;
    end
    
    % Чтение откликов необыкновенных волн
    for i_x = 1:FrequencyData.count_x        
        SignalResponse = readSignalResponse(data,i);
        i = i + 6;
        
        if i + SignalResponse.count_samples-1 > n % если блок битый?
            break;
        end
        arr_x = FrequencyData.thereshold_x + typecast(data(i:i+SignalResponse.count_samples-1), 'uint8');
        % Определим номер начальной высоты
        j = 1 + (SignalResponse.height_begin - uint32(cur.height_min))/uint32(cur.height_step);
        mat_x(j:j+uint32(SignalResponse.count_samples)-1,frq_i) = arr_x;
        
        i = i + SignalResponse.count_samples;
    end
    frq_i = frq_i + 1;
end

% Определение режима рисования
mnu_key_o = 0;
mnu_key_x = 0;
if strcmp(get(findobj('Tag','itmnu_o'),'Checked'),'on')
   mnu_key_o = 1; 
end
if strcmp(get(findobj('Tag','itmnu_x'),'Checked'),'on')
   mnu_key_x = 2; 
end
mnu_key = mnu_key_o + mnu_key_x;

clm = jet;
% Assign white (all 1's) to black (the first row in myColorMap).
clm(1, :) = [1 1 1];
switch mnu_key
   case {1,2}
        if mnu_key_o == 1
            mat = mat_o;
        else
            mat = mat_x;
        end
        max_a = max(max(mat));
        min_a = min(min(mat));
        indNaN = find(isnan(mat));
        mat(indNaN) = 0;
        min_a = min(min(mat));
        img = mat2gray((flipud(mat)-min_a)/max_a);
        
    case {3} % совместный рисунок о- и х- компонент
        % наложение о- на х- компоненту
        ind_o = find(~isnan(mat_o));
        ind_x = find(~isnan(mat_x));
        mat = mat_x;
        %mat(ind_o) = mat_o(ind_o);
        mat(ind_x) = 120;
        mat(ind_o) = 230;
        ind_x = find(isnan(mat));
        mat(ind_x) = 0;
        img = mat2gray(flipud(mat));
        
    otherwise % компоненты не выбраны
        mat = ones(h_n, cur.count_freq);
        img = mat2gray(flipud(mat));
        clm = white;
end

    iptsetpref('ImshowAxesVisible','off')
    h = imshow(img,'Colormap',clm,'XData',[cur.freq_min cur.freq_max]/1000,'YData',[cur.height_min, cur.height_max]/1000);    
    
    set(h,'Tag','imgMain');
    axis normal
    
    pos = get(gca,'Position');
    ax = axes('XLim',[cur.freq_min cur.freq_max]/1000,'YLim',[cur.height_min, cur.height_max]/1000,...
                'Tag','axMain','Position',pos,'Color','none','Box','on',...
                'FontSize',16,'FontWeight','bold');
    set(gca,'FontName','FixedWidth','FontSize',12);
    grid on
      
    xlabel('Частота, МГц');
    
    switch cur.station_id
        case '34502'
            station_name = 'Москва';
        case '34506'
            station_name = 'Ростов-на-Дону';
        otherwise
            station_name = 'неизвестно';
    end
    title(['UTC = ',cur.time_sound, ', ',station_name],'Color','k')     
end

function FrequencyData = readFrequencyData(indata,idx)
% Чтение заголовка строки
% indata - данные
% idx - индех начала строки

data = indata(idx:idx+12);
FrequencyData = struct(...
    'frequency',typecast(data(1:2), 'uint16'),...
    'gain_control',typecast(data(3:4), 'uint16'),...
    'pulse_time',typecast(data(5:6), 'uint16'),...
    'pulse_length',typecast(data(7), 'uint8'),...
    'band',typecast(data(8), 'uint8'),...
    'type',typecast(data(9), 'uint8'),...
    'thereshold_o',typecast(data(10), 'uint8'),...
    'thereshold_x',typecast(data(11), 'uint8'),...
    'count_o',typecast(data(12), 'uint8'),...
    'count_x',typecast(data(13), 'uint8')...
);

function SignalResponse = readSignalResponse(indata,idx)
% Чтение заголовка строки для отдельной волны
% indata - данные
% idx - индех начала строки

data = indata(idx:idx+6);
SignalResponse = struct(...
    'height_begin',typecast(data(1:4), 'uint32'),...
    'count_samples',typecast(data(5:6), 'uint16')...
);

function loadDBContour()
% Загрузка обводки из БД
global connection data data_proc

index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
key = data.Data{index_selected,1};

if ~isempty(connection)
% SQL query to get all fields from Table.
    query = ['SELECT * FROM ionosphere_bottom WHERE iono_key = ', num2str(key)];
    data_proc = myblob_command (connection, query);    
end %if

% Сохраняем данные в структуру
data_proc = struct(...  
    'iono_key',data_proc{1},...      % bigint(11) unsigned NOT NULL
    'fmin',data_proc{2},...          % float
    'fmin_q',data_proc{3},...        % char(2)
    'foF2',data_proc{4},...          % float
    'foF2_q',data_proc{5},...        % char(2)
    'fxF2',data_proc{6},...          % float
    'fxF2_q',data_proc{7},...        % char(2)
    'hmF2',data_proc{8},...          % float
    'fxI',data_proc{9},...           % float
    'fxI_q',data_proc{10},...        % char(2)
    'hpF2',data_proc{11},...         % float
    'hpF2_q',data_proc{12},...       % char(2)
    'hvF2',data_proc{13},...         % float
    'hvF2_q',data_proc{14},...       % char(2)
    'hvF',data_proc{15},...          % float
    'hvF_q',data_proc{16},...        % char(2)
    'foF1',data_proc{17},...         % float
    'foF1_q',data_proc{18},...       % char(2)
    'hmF1',data_proc{19},...         % float
    'foE',data_proc{20},...          % float
    'foE_q',data_proc{21},...        % char(2)
    'hvE',data_proc{22},...          % float
    'hvE_q',data_proc{23},...        % char(2)
    'hmE',data_proc{24},...          % float
    'foEs',data_proc{25},...         % float
    'foEs_q',data_proc{26},...       % char(2)
    'fbEs',data_proc{27},...         % float
    'fbEs_q',data_proc{28},...       % char(2)
    'ftEs',data_proc{29},...         % float
    'ftEs_q',data_proc{30},...       % char(2)
    'hvEs',data_proc{31},...         % float
    'hvEs_q',data_proc{32},...       % char(2)
    'MUF3000F1',data_proc{33},...    % float
    'MUF3000F1_q',data_proc{34},...  % char(2)
    'MUF3000F2',data_proc{35},...    % float
    'MUF3000F2_q',data_proc{36},...  % char(2)
    'type_Es',data_proc{37},...      % varchar(5)
    'type_SpreadF',data_proc{38},... % varchar(5)
    'trace_Eo',data_proc{39},...     % varbinary(240)
    'trace_F1o',data_proc{40},...    % varbinary(240)
    'trace_F2o',data_proc{41},...    % varbinary(240)
    'trace_Ex',data_proc{42},...     % varbinary(240)
    'trace_F1x',data_proc{43},...    % varbinary(240)
    'trace_F2x',data_proc{44},...    % varbinary(240)
    'profile',data_proc{45}...       % blob
);

% Разберём данные обрисовок по массивам
trace_Eo = getArrayFromBLOB(data_proc.trace_Eo);
    plotPoints(trace_Eo, 'o', 5, 'g','E_o');
trace_F1o = getArrayFromBLOB(data_proc.trace_F1o);
    plotPoints(trace_F1o, 'o', 5, 'b','F1_o');
trace_F2o = getArrayFromBLOB(data_proc.trace_F2o);
    plotPoints(trace_F2o, 'o', 5, 'r','F2_o');
trace_Ex = getArrayFromBLOB(data_proc.trace_Ex);
    plotPoints(trace_Ex, 'x', 10, 'g','E_x');
trace_F1x = getArrayFromBLOB(data_proc.trace_F1x);
    plotPoints(trace_F1x, 'x', 10, 'b','F1_x');
trace_F2x = getArrayFromBLOB(data_proc.trace_F2x);
    plotPoints(trace_F2x, 'x', 10, 'r','F2_x');
    
% Отрисовка предельных частот
plotLine(data_proc.foF2, 1, 'r','foF2');
plotLine(data_proc.fxF2, 1, 'g','fxF2');
plotLine(data_proc.foF1, 1, 'b','foF1');
plotLine(data_proc.foE, 1, 'm','foE');
plotLine(data_proc.foEs, 1, 'c','foEs');

% Отрисовка минимальных действующих высот
plotLineX(data_proc.hmF2, 1, 'r','hmF2');
plotLineX(data_proc.hmF1, 1, 'b','hmF1');
plotLineX(data_proc.hmE, 1, 'm','hmE');


function pt = getArrayFromBLOB(blob)
% Извлечение данных из БЛОБ массива БД

pt = 0; % в случае отсутствия данных
if ~ischar(blob) && ~isempty(blob)
    n = length(blob);
    npt = n / 8; 
    if ~rem(n,8)
        for i=1:npt
            indata = blob((i-1)*8+1:i*8);
            pt(i,1) = typecast(indata(1:4), 'single');
            pt(i,2) = typecast(indata(5:8), 'single');
        end
    end
end

function plotPoints(pts, marker, msize, mcolor,tag)
% Отрисовка точек
n = length(pts);
if n > 1
    for i=1:n
        line(pts(i,1),pts(i,2),'Marker',marker,...
            'LineWidth',3,...
            'MarkerEdgeColor',mcolor,...
            'MarkerFaceColor',mcolor,...
            'MarkerSize',msize,...
            'Tag','DBPoints',...
            'DisplayName',[tag, ' из БД'],...
            'LineStyle','none');
    end
end


function plotLine(f, msize, mcolor, tag)
% Отрисовка линий критических частот
if isnan(f) % запись в БД отсутствует
    f = get(gca,'XLim');
    f = f(1);
end
line([f,f],get(gca,'YLim'),...
         'Marker','o',...
         'LineWidth',msize,...
         'Color',mcolor,...
         'Tag','DBLine',...
         'DisplayName',tag);
     
function plotLineX(hm, msize, mcolor, tag)
% Отрисовка линий минимальных действующих высот
if isnan(hm) % запись в БД отсутствует
    hm = get(gca,'YLim');
    hm = hm(1);
end
line(get(gca,'XLim'),[hm,hm],...
         'Marker','o',...
         'LineWidth',msize,...
         'Color',mcolor,...
         'Tag','DBLine',...
         'DisplayName',tag);     

function deleteDBContour()
% Очистка ионограммы от точек обводки из БД

delete(findobj('Tag','DBPoints'));
delete(findobj('Tag','DBLine'));


function plotIRI(varargin)
% Отрисовка IRI
% Если задать функции один аргумент, IRI рассчитывается по значениям,
% скорректированного рабочего графика.
global data data_proc 

index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
dat = data.Data(index_selected,:);

alati = dat{3};
along = dat{4};
vbeg = 80; % double(dat{5})/1000.;
vend = double(dat{6})/1000.;
vstp = 5; % double(dat{7})/1000.;
timesound = dat{2};
iyyyy = str2num(timesound(1:4));
mmdd = str2num([timesound(6:7),timesound(9:10)]);

k = strfind(timesound, ':');
if k(1) == 13
    dhour = str2num(timesound(12)) + str2num(timesound(14:15))/60. + str2num(timesound(17:18))/3600.;
else
    dhour = str2num(timesound(12:13)) + str2num(timesound(15:16))/60. + str2num(timesound(18:19))/3600.;
end
 
%[Ne,h,outf,oarr] = iri2011(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp);
foF2 = data_proc.foF2;
hmF2 = data_proc.hmF2;
foF1 = data_proc.foF1;
hmF1 = data_proc.hmF1;
foE = data_proc.foE;
hmE = data_proc.hmE;
if (nargin == 1) % данные для расчёта берутся из графика
    yl = get(gca,'YLim');
    xl = get(gca,'XLim');
    
    xlim = get(findobj('DisplayName','foF2'),'XData');
        if xlim(1) == xl(1)
            foF2 = 0;
        else
            foF2 = xlim(1);
        end
    ylim = get(findobj('DisplayName','hmF2'),'YData');
        if ylim(1) == yl(1)
            hmF2 = 0;
        else
            hmF2 = ylim(1);
        end    
    xlim = get(findobj('DisplayName','foF1'),'XData');
        if xlim(1) == xl(1)
            foF1 = 0;
        else
            foF1 = xlim(1);
        end        
    ylim = get(findobj('DisplayName','hmF1'),'YData');
        if ylim(1) == yl(1)
            hmF1 = 0;
        else
            hmF1 = ylim(1);
        end    
    xlim = get(findobj('DisplayName','foE'),'XData');    
        if xlim(1) == xl(1)
            foE = 0;
        else
            foE = xlim(1);
        end    
    ylim = get(findobj('DisplayName','hmE'),'YData');
        if ylim(1) == yl(1)
            hmE = 0;
        else
            hmE = ylim(1);
        end
end

% Не используются в расчётах.
hmF2 = NaN;
foF1 = NaN;
hmF1 = NaN;
foE = NaN;
hmE = NaN;
[Ne,h,outf,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, foE, hmE);
% C      #OARR(1) = NMF2/M-3           #OARR(2) = HMF2/KM
fmF2 = sqrt(oarr(1)/(1.24*10^10));
hmF2 = oarr(2);

Ne(find(Ne<0)) = 0;

fN = round(sqrt(Ne/(1.24*10^10))' * 100)/100;

keyobj = findobj('Tag','IRI2016');
if ~isempty(keyobj)
    delete(keyobj)
end

keyobj = findobj('Tag','axNe');
if ~isempty(keyobj)
    delete(keyobj)
end

% Отрисовка плазменной частоты
line(double(fN),double(h),'Tag','IRI2016','Color','k','LineWidth',2,'DisplayName','IRI-2016');
% Отрисовка точки максимума
line(double(fmF2),double(hmF2),'Tag','IRI2016-F2max','Color','k','LineWidth',2,'Marker','o','DisplayName','IRI-2016 F2m');


function plotEDP()
% Отрисовка профиля ИПГ из базы
global data_proc

% ...старый мир разрушим до основанья, а затем...
keyobj = findobj('Tag','EDP');
if ~isempty(keyobj)
    delete(keyobj)
end

% выберем данные
profile = getArrayFromBLOB(data_proc.profile);

if length(profile)>1 % проверка на наличие данных
    % Отрисовка плазменной частоты
    line(profile(:,1),profile(:,2),'Tag','EDP','Color','g','LineWidth',2,'DisplayName','ИПГ из БД');
end
    

function key = isSelectedGroup2
% Возвращает 1 если выбрана кнопка из группы 2
state(1) = get(findobj('Tag','toggle_Es'),'State');
state(2) = get(findobj('Tag','toggle_E'),'State');
state(3) = get(findobj('Tag','toggle_F1'),'State');
state(4) = get(findobj('Tag','toggle_F2'),'State');
key = any(strcmp(state,'on'));

% Определяет, когда указатель находится рядом с линией.
function yesno = isIntoLine(hline)
 currentPoint = get(gca,'CurrentPoint');
    
    x = currentPoint(1,1);
    xline = get(hline,'XData');                   
    Dx = abs(xline(1) - x(1));
    
    yesno = 0;

    x_lim = get(gca,'XLim');
    dx = (x_lim(2) - x_lim(1))*0.02; % 2% от разброса
    
if Dx <= dx
    yesno = 1;
end    

% Определяет, когда указатель находится рядом с линией высот.
function yesno = isIntoLineX(hline)
 currentPoint = get(gca,'CurrentPoint');
    
    y = currentPoint(1,2);
    yline = get(hline,'YData');                   
    Dy = abs(yline(1) - y(1));
    
    yesno = 0;

    y_lim = get(gca,'YLim');
    dy = (y_lim(2) - y_lim(1))*0.02; % 2% от разброса
    
if Dy <= dy
    yesno = 1;
end  

% Определяет выбранную линию критических частот.
function hline = getSelectedLine(handles)

hline = [];
if strcmp(get(handles.toggle_fc,'State'),'on')
    % компонента
    if strcmp(get(handles.toggle_o,'State'),'on')
        str1 = 'fo';
    end
    if strcmp(get(handles.toggle_x,'State'),'on')
        str1 = 'fx';
    end

    % слой ионосферы
    if strcmp(get(handles.toggle_Es,'State'),'on')
        str2 = 'Es';
    end
    if strcmp(get(handles.toggle_E,'State'),'on')
        str2 = 'E';
    end
    if strcmp(get(handles.toggle_F1,'State'),'on')
        str2 = 'F1';
    end
    if strcmp(get(handles.toggle_F2,'State'),'on')
        str2 = 'F2';
    end
    str = [str1, str2];
    hline = findobj('DisplayName',str);
end

% Определяет выбранную линию минимальных действующих высот.
function hline = getSelectedLineX(handles)

hline = [];
if strcmp(get(handles.toggle_hm,'State'),'on')

    % слой ионосферы
    if strcmp(get(handles.toggle_hmE,'State'),'on')
        str = 'hmE';
    end
    if strcmp(get(handles.toggle_hmF1,'State'),'on')
        str = 'hmF1';
    end
    if strcmp(get(handles.toggle_hmF2,'State'),'on')
        str = 'hmF2';
    end
    hline = findobj('DisplayName',str);
end

% Delete points into any click neighbourhood.
function deletePoint(pt)
x0 = pt(1,1);
y0 = pt(1,2);
    dx = x0/20;
    dy = y0/20;
h = findobj('Type','line');
for i = 1:size(h)
    str = get(h(i),'Tag');
    key1 = strfind(str, 'trace_');
    if isempty(key1)
        key1 = 0;
    end
    if key1 == 1 || strcmp(str,'DBPoints')
        x = get(h(i),'XData');
        y = get(h(i),'YData');
        if ((x0-dx) < x(1)) & ((x0+dx) > x(1)) & ...
           ((y0-dy) < y(1)) & ((y0+dy) > y(1))
                
            delete(h(i));
            break
        end
    end
end

% --------------------------------------------------------------------
% Конец блока рукописных функций
% --------------------------------------------------------------------

%% Начало блока генерированных функций событий
% --------------------------------------------------------------------

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dbName user password connection data

% Open DB connection here.
try
    %connection = database(dbName,user,password,drv,url);
    connection = database(dbName,user,password,...
        'Vendor','MySQL',...
        'Server','localhost');
catch ME
    disp('Удаленная база не найдена.');
    ParamMenuItem_Callback(hObject, eventdata, handles);
    connection = database(dbName,user,password,...
        'Vendor','MySQL',...
        'Server','localhost');
end

% Найдем все года в БД
query = 'select distinct(year(time_sound)) from ionosphere.parus ORDER BY 1';
curs = exec(connection, query);
data = get(fetch(curs));
% Заполним выпадающий список
lst_years = data.Data(:,1);
n_lst_years = length(lst_years);
set(handles.popupmenuYear, 'String', lst_years);
set(handles.popupmenuYear, 'Value', n_lst_years);
close(curs);

% Найдём все месяцы в последнем годе
query = ['select distinct MONTHNAME(time_sound), month(time_sound) from ionosphere.parus where year(time_sound) = ',...
    num2str(lst_years{n_lst_years}),...
    ' ORDER BY 2'];
curs = exec(connection, char(query));
data = get(fetch(curs));
lst_months = data.Data(:,1);
numbers_months = data.Data(:,2);
n_lst_months = length(lst_months);
set(handles.popupmenuMonth, 'String', lst_months);
set(handles.popupmenuMonth, 'Value', n_lst_months);
close(curs);

query = [...
    'SELECT iono_key, time_sound, latitude, longitude, height_min, height_max, height_step, station_id FROM parus ',...
    'WHERE YEAR(time_sound) = ', num2str(lst_years{n_lst_years}), ' AND ',...
    'MONTH(time_sound) = ', num2str(numbers_months{n_lst_months}), ' ORDER BY time_sound'];
curs = exec(connection, char(query));
data = get(fetch(curs));
close(curs);

set(handles.listboxDate, 'String', data.Data(:,2)');
key = data.Data{get(handles.listboxDate, 'Value'),1};
plotRow(key);
loadDBContour();
plotIRI();
plotEDP();


% --------------------------------------------------------------------
function ParamMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to ParamMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dbName user password

prompt={'Database name:',...
        'Username:',...
        'Password:'};

name='Connection properties';
numlines=1;
defaultanswer={dbName,user,password};
answer=inputdlg(prompt,name,numlines,defaultanswer);

dbName = answer{1};
user = answer{2};
password = answer{3};

save('ini', 'dbName', 'user', 'password', '-append');


% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selection = questdlg(['Close ' get(handles.db_interface,'Name') '?'],...
                     ['Close ' get(handles.db_interface,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

db_interface_CloseRequestFcn(hObject, eventdata, handles)

delete(handles.db_interface)


% --- Executes on selection change in listboxDate.
function listboxDate_Callback(hObject, eventdata, handles)
% hObject    handle to listboxDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listboxDate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxDate
global data

ccursor = get(gcf,'Pointer');
set(gcf,'Pointer','watch');

index_selected = int32(get(hObject,'Value'));
key = data.Data{index_selected,1};
plotRow(key);

if strcmp(get(findobj('Tag','itmnu_FromDB'),'Checked'),'on')
    loadDBContour();
end

if strcmp(get(findobj('Tag','itmnu_IRI'),'Checked'),'on')
    plotIRI();
    plotEDP();
end

% Отработка включенной легенды
if strcmp(get(findobj('Tag','uitoggletoolLegend'),'State'),'on')
    uitoggletoolLegend_OnCallback(hObject, eventdata, handles);
end

set(gcf,'Pointer',ccursor);


% --- Executes during object creation, after setting all properties.
function listboxDate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function db_interface_CreateFcn(hObject, eventdata, handles)
% hObject    handle to db_interface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global dbName user password url drv

% Загрузим инициализационные данные
load 'ini.mat'


% --- Executes when user attempts to close db_interface.
function db_interface_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to db_interface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global connection

close(connection);
delete(hObject);

% --------------------------------------------------------------------
function itmnu_FromDB_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_FromDB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

isChecked = get(hObject,'Checked');
if strcmp(isChecked,'on')
    set(hObject,'Checked','off');
    deleteDBContour();
else
    set(hObject,'Checked','on');
    loadDBContour();
end


% --------------------------------------------------------------------
function itmnu_ox_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isChecked = get(hObject,'Checked');
if strcmp(isChecked,'on')
    set(hObject,'Checked','off');
else
    set(hObject,'Checked','on');
end
OpenMenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function itmnu_IRI_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_IRI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

isChecked = get(hObject,'Checked');
if strcmp(isChecked,'on')
    set(hObject,'Checked','off');
    
    keyobj = findobj('Tag','IRI2016');
    if ~isempty(keyobj)
        delete(keyobj)
    end
    keyobj = findobj('Tag','IRI2016-F2max');
    if ~isempty(keyobj)
        delete(keyobj)
    end
else
    set(hObject,'Checked','on');
    plotIRI(1);
end


% --------------------------------------------------------------------
function itmnu_EDP_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_EDP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isChecked = get(hObject,'Checked');
if strcmp(isChecked,'on')
    set(hObject,'Checked','off');
    
    keyobj = findobj('Tag','EDP');
    if ~isempty(keyobj)
        delete(keyobj)
    end
else
    set(hObject,'Checked','on');
    plotEDP();
end

% --------------------------------------------------------------------
function itmnu_NIIF_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_NIIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isChecked = get(hObject,'Checked');
if strcmp(isChecked,'on')
    set(hObject,'Checked','off');
    
    keyobj = findobj('Tag','NIIF');
    if ~isempty(keyobj)
        delete(keyobj)
    end

%     keyobj = findobj('Tag','axNe');
%     if keyobj
%         delete(keyobj)
%     end
else
    set(hObject,'Checked','on');
    %plotNIIF();
end

% --------------------------------------------------------------------
function itmnu_Erase_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_Erase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function itmnu_Save_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data data_proc

% Сохраним данные из БД
%itmnu_Save_DB_Callback(hObject, eventdata, handles);

% Сохраним данные (быть может измененные) с графика
% Данные зондирования
index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
alati = data.Data{index_selected,3};
along = data.Data{index_selected,4};
timesound = data.Data{index_selected,2};

% Компоненты имени файла
% Название файла
switch data.Data{index_selected,8}
    case '34502'
        station_name = 'm';
    case '34506'
        station_name = 'r';
    otherwise
        station_name = '';
end
fname = strcat(timesound(1:10),'_',timesound(12:13),'-',timesound(15:16),station_name);

% Поиск линий
hlines = findobj('Tag','DBLine');
    tags = get(hlines,'DisplayName');
n = length(hlines);
for i = 1:n
    h = hlines(i);
    switch tags{i}
        case 'foF2'   
            foF2 = getF(h);  
        case 'fxF2'   
            fxF2 = getF(h);  
        case 'foF1'   
            foF1 = getF(h);  
        case 'foE'   
            foE = getF(h);            
        case 'foEs'   
            foEs = getF(h);               
        case 'hmF2'
            hmF2 = getH(h); 
        case 'hmF1'
            hmF1 = getH(h); 
        case 'hmE'
            hmE = getH(h);           
        otherwise
            disp('Неизвестный таг.')
    end
end

% Получаем оцифровки следов. Если отсутствуют на графике, используем
% сохранённые последовательности из БД.
trace_Eo = getTrace('trace_Eo');
    if isempty(trace_Eo)
        trace_Eo = getArrayFromBLOB(data_proc.trace_Eo);
    end
trace_Ex = getTrace('trace_Ex');
    if isempty(trace_Ex)
        trace_Ex = getArrayFromBLOB(data_proc.trace_Ex);
    end
trace_F1o = getTrace('trace_F1o');
    if isempty(trace_F1o)
        trace_F1o = getArrayFromBLOB(data_proc.trace_F1o);
    end
trace_F1x = getTrace('trace_F1x');
    if isempty(trace_F1x)
        trace_F1x = getArrayFromBLOB(data_proc.trace_F1x);
    end
trace_F2o = getTrace('trace_F2o');
    if isempty(trace_F2o)
        trace_F2o = getArrayFromBLOB(data_proc.trace_F2o);
    end
trace_F2x = getTrace('trace_F2x');
    if isempty(trace_F2x)
        trace_F2x = getArrayFromBLOB(data_proc.trace_F2x);
    end
profile = getArrayFromBLOB(data_proc.profile);

% Формируем структуру для записи в БД
foF2 = ifempty(foF2);
fxF2 = ifempty(fxF2);
hmF2 = ifempty(hmF2);
foF1 = ifempty(foF1);
hmF1 = ifempty(hmF1);
foE = ifempty(foE);
hmE = ifempty(hmE);
foEs = ifempty(foEs);
colnames = {'foF2', 'fxF2', 'hmF2', 'foF1', 'hmF1', 'foE', 'hmE', 'foEs', ...
            'trace_Eo', 'trace_F1o', 'trace_F2o', 'trace_Ex', 'trace_F1x', 'trace_F2x'};
exdata = {foF2, fxF2, hmF2, foF1, hmF1, foE, hmE, foEs, ...
          getBLOBFromArray(trace_Eo), getBLOBFromArray(trace_F1o), getBLOBFromArray(trace_F2o),...
          getBLOBFromArray(trace_Ex), getBLOBFromArray(trace_F1x), getBLOBFromArray(trace_F2x)};
updateDB(colnames,exdata);

save(fname, 'timesound', 'index_selected', 'alati', 'along', ...
    'trace_Eo', 'trace_F1o', 'trace_F2o', 'trace_Ex', 'trace_F1x', 'trace_F2x', 'profile', ...
    'foF2', 'fxF2', 'hmF2', 'foF1', 'hmF1', 'foE', 'hmE', 'foEs');

function out = ifempty(x)

out = x;
if isempty(x)
    out = 0;
else
    out = x;
end

% --------------------------------------------------------------------
function updateDB(colnames,exdata)
% Сохранение изменений графика в БД
global data connection

table = 'ionosphere_bottom';

index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
key = data.Data{index_selected,1};

if ~isempty(connection)
    %set(conn,'AutoCommit','off');
    whereclause = strcat('WHERE iono_key = ''',num2str(key),'''');
    % MySQL syntax:
    % UPDATE table_name SET field1=new-value1, field2=new-value2 [WHERE Clause]
    for i = 1:length(colnames)
        tmpdata = exdata{i};
        if tmpdata(1) ~= 0,
            if length(tmpdata) ~= 1,
                tmpdata = mysql_real_escape_string(tmpdata);
                % Send a Matlab object to the database;
                % Therefore:
                % Generate a Matlab matrix.
                % Ask Matlab for a temporary MAT-file to buffer the matrix.
                % Save the matrix to the MAT-file using Matlab's SAVE command.
                % Define where (table, column, row, ...) the MAT-file/matrix/blob will be stored.
                % Send the blob to the database.  
                tmpfile = [tempname, '.mat'];
                save (tmpfile, 'tmpdata');
                myblob_to_db (connection, table, colnames{i}, whereclause, tmpfile);
            else              
                % Insert a new row with a description in the string column
                % and leave the blob column empty
                command = ['UPDATE ionosphere_bottom SET ',colnames{i},'=''',num2str(tmpdata),'''  ',whereclause];
                myblob_command (connection, command);
            end
            %query = strcat('UPDATE ionosphere_bottom SET ',colnames(i),'=''',tmpdata,'''  ',whereclause);
            %update(conn, 'ionosphere_bottom', colnames(i), tmpdata, whereclause);
            %curs = exec(conn, query);
        end
    end
    %commit(conn);   
    % Close connection to mySQL database
    % Don't forget to do this... ;-)
    myblob_close (connection);
end %if

% --------------------------------------------------------------------
function rstr = mysql_real_escape_string(data)
% библиотечная функция MySQL mysql_real_escape_string, 
% добавляет обратную косую черту к следующим символам: 
% \x00, \n, \r, \, ', " и \x1a.
% ['0x00', '0x0A', '0x0D', '0x5C', '0x27', '0x22', '0x1A']
n = length(data);
rstr = [];
escape_arr = [0, 10, 13, 92, 39, 34, 26];
add_sumb = 92; % '0x5C';

if n > 0
    for i=1:n
        must_add = 0;
        for j = 1:7
            if data(i) == escape_arr(j)
                must_add = 1;
            end
        end
        if must_add
            rstr = horzcat(rstr,add_sumb,data(i));
        else
            rstr = horzcat(rstr,data(i));
        end
    end
end

% --------------------------------------------------------------------
function blob = getBLOBFromArray(arr)
% Преобразование двумерных массивов данных в БЛОБ для сохранения в БД

blob = 0; % в случае отсутствия данных
if ~isempty(arr)
    a = single(arr);
    blob = reshape(a',[],1)'; 
    blob = typecast(blob, 'uint8');
end

% --------------------------------------------------------------------
function h = getH(handle)
% Получает значение высоты
hLim = get(gca,'YLim');

h = get(handle,'YData');
h = h(1);

if h <= hLim(1) || h >= hLim(2)
    h = [];
end

% --------------------------------------------------------------------
function f = getF(handle)
% Получает значение частоты
fLim = get(gca,'XLim');

f = get(handle,'XData');
f = f(1);

if f <= fLim(1) || f >= fLim(2)
    f = [];
end

% --------------------------------------------------------------------
function trace = getTrace(tag)
% Получает двумерный массив оцифровки
hlines = findobj('Tag',tag);
n = length(hlines);
trace = [];
for i = 1:n
    f = get(hlines(i),'XData');
    h = get(hlines(i),'YData');
    trace(i,1) = f;
    trace(i,2) = h;
end

    
% --------------------------------------------------------------------
function itmnu_Profile_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_Profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function toggle_o_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_x,'State','off');
enableToggleGroup2(handles);

handles.faceMarker = 'o';
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function toggle_x_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_o,'State','off');
enableToggleGroup2(handles);

handles.faceMarker = 'x';
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toggle_pt_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_fc,'State','off');
set(handles.toggle_hm,'State','off');
set(handles.toggle_delete,'State','off');
enableToggleGroup1(handles);
enableToggleGroup2(handles);


% --------------------------------------------------------------------
function toggle_fc_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_pt,'State','off');
set(handles.toggle_hm,'State','off');
set(handles.toggle_delete,'State','off');
enableToggleGroup1(handles);
enableToggleGroup2(handles);

function enableToggleGroup1(handles)
% Разрешение выбора компонент.
if strcmp(get(handles.toggle_o,'Enable'),'off')
    set(handles.toggle_o,'Enable','on');
    set(handles.toggle_x,'Enable','on');
end

function enableToggleGroup2(handles)
% Разрешение выбора компонент.
if strcmp(get(handles.toggle_Es,'Enable'),'off')
    set(handles.toggle_Es,'Enable','on');
    set(handles.toggle_E,'Enable','on');
    set(handles.toggle_F1,'Enable','on');
    set(handles.toggle_F2,'Enable','on');
end

function enableToggleGroup3(handles)
% Разрешение выбора компонент.
if strcmp(get(handles.toggle_hmF2,'Enable'),'off')
    set(handles.toggle_hmF2,'Enable','on');
    set(handles.toggle_hmF1,'Enable','on');
    set(handles.toggle_hmE,'Enable','on');
end

% --------------------------------------------------------------------
function toggle_Group1_OffCallback(hObject, eventdata, handles)
% hObject    handle to toggle_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_Es,'Enable','off');
set(handles.toggle_E,'Enable','off');
set(handles.toggle_F1,'Enable','off');
set(handles.toggle_F2,'Enable','off');

set(handles.toggle_o,'Enable','off');
set(handles.toggle_x,'Enable','off');

% --------------------------------------------------------------------
function toggle_Group2_OffCallback(hObject, eventdata, handles)
% hObject    handle to toggle_o (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_Es,'Enable','off');
set(handles.toggle_E,'Enable','off');
set(handles.toggle_F1,'Enable','off');
set(handles.toggle_F2,'Enable','off');

% --------------------------------------------------------------------
function toggle_Es_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_Es (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_E,'State','off');
set(handles.toggle_F1,'State','off');
set(handles.toggle_F2,'State','off');

handles.faceColor = 'c';
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toggle_E_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_Es,'State','off');
set(handles.toggle_F1,'State','off');
set(handles.toggle_F2,'State','off');

handles.faceColor = 'm';
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toggle_F1_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_F1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_Es,'State','off');
set(handles.toggle_E,'State','off');
set(handles.toggle_F2,'State','off');

handles.faceColor = 'b';
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function toggle_F2_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_F2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_Es,'State','off');
set(handles.toggle_E,'State','off');
set(handles.toggle_F1,'State','off');

handles.faceColor = 'r';
% Update handles structure
guidata(hObject, handles);


% --- Executes on mouse press over axes background.
function db_interface_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to db_interface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.keyButtonDown = true;
% Тип следа
switch handles.faceColor
    case 'r'   
        type = 'F2';  
    case 'b'   
        type = 'F1';  
    case 'm'   
        type = 'E';  
    case 'c'   
        type = 'Es';
    otherwise
        disp('Неизвестный тип следа.')
end

currentPoint = get(gca,'CurrentPoint');

% Обработка обрисовки следа точками
if strcmp(get(handles.toggle_pt,'State'),'on')
    x = currentPoint(1,1);
    y = currentPoint(1,2);
    yLim = get(gca,'YLim');
    xLim = get(gca,'XLim');
    if (x<xLim(2)) & (x>xLim(1)) & (y<yLim(2)) & (y>yLim(1)),
        if strcmp(get(handles.toggle_o,'State'),'on')
            % 'trace_Eo', 'trace_F1o', 'trace_F2o'
            line(x,y,'Marker',handles.faceMarker,...
                    'LineWidth',2,...
      				'MarkerEdgeColor','k',...
                    'MarkerFaceColor',handles.faceColor,...
                    'MarkerSize',7,...
                    'Tag',strcat('trace_',type,'o'),...
                    'DisplayName',strcat(type,'o',' вручную'));
        end
        if strcmp(get(handles.toggle_x,'State'),'on')
            % 'trace_Ex', 'trace_F1x', 'trace_F2x'
            line(x,y,'Marker',handles.faceMarker,...
                    'LineWidth',2,...
      				'MarkerEdgeColor',handles.faceColor,...
                    'MarkerSize',10,...
                    'Tag',strcat('trace_',type,'x'),...
                    'DisplayName',strcat(type,'x',' вручную'));
        end        
    end
end   

% Отработка удаления точек
if strcmp(get(gcf,'Pointer'),'circle')
    deletePoint(currentPoint);
end
   
% Update handles structure
guidata(hObject, handles);



% --- Executes on mouse motion over figure - except title and menu.
function db_interface_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to db_interface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 currentPoint = get(gca,'CurrentPoint');
    x = currentPoint(1,1);
    y = currentPoint(1,2);
   
    hline = getSelectedLine(handles);
    xline = get(hline,'XData');
   
    if ~isempty(xline)
        if isIntoLine(hline)
            set(gcf,'Pointer','fleur');
            if handles.keyButtonDown                   
                set(hline,'XData',[x x]);
            end
        else
            set(gcf,'Pointer','arrow');            
        end
    else
        hline = getSelectedLineX(handles);
        yline = get(hline,'YData');
        if ~isempty(yline)
            if isIntoLineX(hline)
                set(gcf,'Pointer','fleur');
                if handles.keyButtonDown                   
                    set(hline,'YData',[y y]);
                end
            else
                set(gcf,'Pointer','arrow');            
            end
        end
    end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function db_interface_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to db_interface (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.keyButtonDown = false;

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function itmnu_Save_DB_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_Save_DB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data data_proc

% Сохраняем данные из БД в mat-файл для дальнейшей обработки.

% Компоненты имени файла
% Название файла
switch data.Data{index_selected,8}
    case '34502'
        station_name = 'm';
    case '34506'
        station_name = 'r';
    otherwise
        station_name = '';
end
fname = strcat(timesound(1:10),'_',timesound(12:13),'-',timesound(15:16),station_name);

% Разберём данные обрисовок по массивам
trace_Eo = getArrayFromBLOB(data_proc.trace_Eo);
trace_F1o = getArrayFromBLOB(data_proc.trace_F1o);
trace_F2o = getArrayFromBLOB(data_proc.trace_F2o);
trace_Ex = getArrayFromBLOB(data_proc.trace_Ex);
trace_F1x = getArrayFromBLOB(data_proc.trace_F1x);
trace_F2x = getArrayFromBLOB(data_proc.trace_F2x);
profile = getArrayFromBLOB(data_proc.profile);

foF2 = data_proc.foF2;
fxF2 = data_proc.fxF2;
hmF2 = data_proc.hmF2;
foF1 = data_proc.foF1;
hmF1 = data_proc.hmF1;
foE = data_proc.foE;
hmE = data_proc.hmE;
foEs = data_proc.foEs;

% Данные зондирования
index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
alati = data.Data{index_selected,3};
along = data.Data{index_selected,4};
timesound = data.Data{index_selected,2};

% Сохранение
save(fname, 'timesound', 'index_selected', 'alati', 'along', ...
    'trace_Eo', 'trace_F1o', 'trace_F2o', 'trace_Ex', 'trace_F1x', 'trace_F2x', 'profile', ...
    'foF2', 'fxF2', 'hmF2', 'foF1', 'hmF1', 'foE', 'hmE', 'foEs');


% --------------------------------------------------------------------
function toggle_hm_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_hm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_fc,'State','off');
set(handles.toggle_pt,'State','off');
set(handles.toggle_delete,'State','off');

enableToggleGroup3(handles);


% --------------------------------------------------------------------
function toggle_hmE_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_hmE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_hmF1,'State','off');
set(handles.toggle_hmF2,'State','off');

handles.faceColor = 'm';
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function toggle_hmF1_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_hmF1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_hmE,'State','off');
set(handles.toggle_hmF2,'State','off');

handles.faceColor = 'b';
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function toggle_hmF2_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_hmF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_hmE,'State','off');
set(handles.toggle_hmF1,'State','off');

handles.faceColor = 'r';
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function toggle_hm_OffCallback(hObject, eventdata, handles)
% hObject    handle to toggle_hm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.toggle_hmE,'Enable','off');
set(handles.toggle_hmF1,'Enable','off');
set(handles.toggle_hmF2,'Enable','off');


% --------------------------------------------------------------------
function uipushtool_calculate_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool_calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotIRI(1);


% --------------------------------------------------------------------
function SaveFigMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFigMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data

% Данные зондирования
index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
timesound = deblank(data.Data{index_selected,2});
i1 = findstr(timesound, ':');
i2 = findstr(timesound, ' ');
timesound(i1) = '-';
timesound(i2) = '_';
fname = timesound(1:19);

saveas(gcf,fname,'fig');



% --------------------------------------------------------------------
function SaveIonMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveIonMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data

% Данные зондирования
index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
timesound = deblank(data.Data{index_selected,2});
i1 = findstr(timesound, ':');
i2 = findstr(timesound, ' ');
timesound(i1) = '-';
timesound(i2) = '_';
fname = timesound(1:19);

ppm = get(gcf,'PaperPositionMode');
set(gcf,'PaperPositionMode','auto');
    print('-noui','-dtiff','-r300',fname);
set(gcf,'PaperPositionMode',ppm);


% --------------------------------------------------------------------
function toggle_delete_OffCallback(hObject, eventdata, handles)
% hObject    handle to toggle_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'Pointer','arrow');


% --------------------------------------------------------------------
function toggle_delete_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggle_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Запрещение выбора компонент.
set(handles.toggle_o,'Enable','off');
set(handles.toggle_x,'Enable','off');
set(handles.toggle_Es,'Enable','off');
set(handles.toggle_E,'Enable','off');
set(handles.toggle_F1,'Enable','off');
set(handles.toggle_F2,'Enable','off');
set(handles.toggle_hmE,'Enable','off');
set(handles.toggle_hmF1,'Enable','off');
set(handles.toggle_hmF2,'Enable','off');

set(gcf,'Pointer','circle');
set(handles.toggle_fc,'State','off');
set(handles.toggle_hm,'State','off');
set(handles.toggle_pt,'State','off');

%%
function data = myblob_command (connection, query)

curs = exec(connection, query);
ddata = fetch(curs);
data = ddata.Data;


function changeYearMonthDay()
global connection data

years = get(findobj('Tag','popupmenuYear'), 'String');
cur_year = years{get(findobj('Tag','popupmenuYear'),'Value')};
months = get(findobj('Tag','popupmenuMonth'), 'String');
cur_month = months{get(findobj('Tag','popupmenuMonth'),'Value')};

query = [...
    'SELECT iono_key, time_sound, latitude, longitude, height_min, height_max, height_step, station_id FROM parus ',...
    'WHERE YEAR(time_sound) = ', cur_year, ' AND ',...
    'MONTHNAME(time_sound) = ', strcat('''',cur_month,''''), ' ORDER BY time_sound'];

curs = exec(connection, char(query));
data = get(fetch(curs));
close(curs);

set(findobj('Tag','listboxDate'), 'String', data.Data(:,2)');
key = data.Data{get(findobj('Tag','listboxDate'), 'Value'),1};
plotRow(key);
loadDBContour();


% --- Executes during object creation, after setting all properties.
function popupmenuYear_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuYear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenuMonth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuYear.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuYear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuYear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuYear


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuYear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuMonth.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuMonth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuMonth


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxDate.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listboxDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxDate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxDate


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxDate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuYear.
function popupmenuYear_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuYear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global connection
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuYear contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuYear

lst_years = get(hObject, 'String');
n_lst_years = get(hObject, 'Value');

% Найдём все месяцы в последнем годе
query = join(string({'select distinct MONTHNAME(time_sound), month(time_sound) from ionosphere.parus where year(time_sound) = ',...
    num2str(lst_years{n_lst_years}),...
    ' ORDER BY 2'}));
curs = exec(connection, char(query));
data = get(fetch(curs));
close(curs);
lst_months = data.Data(:,1);
n_lst_months = length(lst_months);
set(handles.popupmenuMonth, 'String', lst_months);
set(handles.popupmenuMonth, 'Value', n_lst_months);

changeYearMonthDay();


% --- Executes on selection change in popupmenuMonth.
function popupmenuMonth_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuMonth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuMonth contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuMonth
changeYearMonthDay();


% --------------------------------------------------------------------
function itmnu_Save_Real_to_MAT_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_Save_Real_to_MAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global data data_proc
global last_saved_mat

% Сохраним данные из БД
%itmnu_Save_DB_Callback(hObject, eventdata, handles);

% Сохраним данные (быть может измененные) с графика
% Данные зондирования
index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
alati = data.Data{index_selected,3};
along = data.Data{index_selected,4};
timesound = data.Data{index_selected,2};

% Компоненты имени файла
% Название файла
switch data.Data{index_selected,8}
    case '34502'
        station_name = 'm';
    case '34506'
        station_name = 'r';
    otherwise
        station_name = '';
end
fname = strcat(timesound(1:10),'_',timesound(12:13),'-',timesound(15:16),station_name);

% Поиск линий
hlines = findobj('Tag','DBLine');
    tags = get(hlines,'DisplayName');
n = length(hlines);
for i = 1:n
    h = hlines(i);
    switch tags{i}
        case 'foF2'   
            foF2 = getF(h);  
        case 'fxF2'   
            fxF2 = getF(h);  
        case 'foF1'   
            foF1 = getF(h);  
        case 'foE'   
            foE = getF(h);            
        case 'foEs'   
            foEs = getF(h);               
        case 'hmF2'
            hmF2 = getH(h); 
        case 'hmF1'
            hmF1 = getH(h); 
        case 'hmE'
            hmE = getH(h);           
        otherwise
            disp('Неизвестный таг.')
    end
end

% Получаем оцифровки следов. Если отсутствуют на графике, используем
% сохранённые последовательности из БД.
trace_Eo = getTrace('trace_Eo');
    if isempty(trace_Eo)
        trace_Eo = getArrayFromBLOB(data_proc.trace_Eo);
    end
trace_Ex = getTrace('trace_Ex');
    if isempty(trace_Ex)
        trace_Ex = getArrayFromBLOB(data_proc.trace_Ex);
    end
trace_F1o = getTrace('trace_F1o');
    if isempty(trace_F1o)
        trace_F1o = getArrayFromBLOB(data_proc.trace_F1o);
    end
trace_F1x = getTrace('trace_F1x');
    if isempty(trace_F1x)
        trace_F1x = getArrayFromBLOB(data_proc.trace_F1x);
    end
trace_F2o = getTrace('trace_F2o');
    if isempty(trace_F2o)
        trace_F2o = getArrayFromBLOB(data_proc.trace_F2o);
    end
trace_F2x = getTrace('trace_F2x');
    if isempty(trace_F2x)
        trace_F2x = getArrayFromBLOB(data_proc.trace_F2x);
    end
profile = getArrayFromBLOB(data_proc.profile);

% Получение профиля по IRI
tmp = findobj('Tag','IRI2016');
iri_profile = [get(tmp,'XData')',get(tmp,'YData')'];

save(fname, 'timesound', 'index_selected', 'alati', 'along', ...
    'trace_Eo', 'trace_F1o', 'trace_F2o', 'trace_Ex', 'trace_F1x', 'trace_F2x', 'profile', 'iri_profile', ...
    'foF2', 'fxF2', 'hmF2', 'foF1', 'hmF1', 'foE', 'hmE', 'foEs');

% Сохраним имя файла для истории
last_saved_mat = fname;

% --------------------------------------------------------------------
function uitoggletoolLegend_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
legend('off');


% --------------------------------------------------------------------
function uitoggletoolLegend_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletoolLegend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[lgd,icons,plots,txt] = legend('show');

% Поиск несовпадающих (уже отсортированных!) элементов
n = length(txt);
numbers = zeros(1,n);
numbers(1) = 1;
j = 1;
for i=2:n
    if ~strcmp(txt{i}, txt{numbers(1,j)})
       j = j + 1;
       numbers(1,j) = i;
    end
end

% Заполнение легенды
hs = zeros(1,j);
for i=1:j
    tmp = findobj('DisplayName',txt{numbers(i)});
    hs(i) = tmp(1);
end
legend(hs, txt{numbers(1:j)});


% --- Executes on selection change in popupmenuDay.
function popupmenuDay_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuDay contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuDay
changeYearMonthDay()

% --- Executes during object creation, after setting all properties.
function popupmenuDay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuDay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function itmnu_Save_MAT_and_Work_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_Save_MAT_and_Work (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global last_saved_mat

% Вначале сохраним MAT-файл для последующей независимой обработки
itmnu_Save_Real_to_MAT_Callback();
S = load(last_saved_mat);

% Внешняя процедура обработки
tmp = flipud(S.trace_F2o);
f = tmp(:,1)';
hv = tmp(:,2)';

% Подготовим IRI профиль
tmp = S.iri_profile;
f_N_tmp = tmp(:,1)';
h_N_tmp = tmp(:,2)';
% Определение точки максимума (должна быть отображена на графике)
obj = findobj('Tag','IRI2016-F2max');
fmF2 = get(obj,'XData');
hmF2  = get(obj,'YData');
% Обрезание по точке максимума
ind = find(h_N_tmp < hmF2);
f_N = f_N_tmp(ind);
h_N = h_N_tmp(ind);
% Дополним точкой максимума
f_N = [f_N, fmF2];
h_N = [h_N, hmF2];

H = Denisenko_Nigth_Nh_from_hv_IRI_real(f, hv, f_N, h_N, S.timesound);
savefig(H,last_saved_mat);


% --------------------------------------------------------------------
function itmnu_Save_MAT_and_Work_Cor_fmE_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_Save_MAT_and_Work_Cor_fmE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global last_saved_mat
global data data_proc 

index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
dat = data.Data(index_selected,:);

alati = dat{3};
along = dat{4};
vbeg = 80; % double(dat{5})/1000.;
vend = double(dat{6})/1000.;
vstp = 5; % double(dat{7})/1000.;
timesound = dat{2};
iyyyy = str2num(timesound(1:4));
mmdd = str2num([timesound(6:7),timesound(9:10)]);

k = strfind(timesound, ':');
if k(1) == 13
    dhour = str2num(timesound(12)) + str2num(timesound(14:15))/60. + str2num(timesound(17:18))/3600.;
else
    dhour = str2num(timesound(12:13)) + str2num(timesound(15:16))/60. + str2num(timesound(18:19))/3600.;
end

% Вначале сохраним MAT-файл для последующей независимой обработки
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% В процессе корректировки изменяется профиль IRI, но это не отображается
% на графике!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
itmnu_Save_Real_to_MAT_Callback();
S = load(last_saved_mat);

% Получим оцифровку
tmp = flipud(S.trace_F2o);
f = tmp(:,1)';
hv = tmp(:,2)';

% Получим подлежащий корректировке IRI профиль
tmp = S.iri_profile;
f_ = tmp(:,1);
for i=1:length(f_)-2 % Поиск локального минимума, соответствующего E-слою
    if f_(i)<=f_(i+1)&&f_(i+1)>f_(i+2)
        fmE = f_(i+1);
        break
    end
end

% Определяем частоты fmE построения корректировок
foE = fmE:0.05:f(1)-0.01;
foF2 = get(findobj('DisplayName','foF2'),'XData');
foF2 = foF2(1);
% Не используются в расчётах.
hmF2 = NaN;
foF1 = NaN;
hmF1 = NaN;
hmE = NaN;

SE = zeros(length(foE),1);
SF = zeros(length(foE),1);
calc_F = 1;
for i = 1:length(foE)
    [Ne,h,outf,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, foE(i), hmE);

    % Подготавливаем скорректированный профиль
    Ne(find(Ne<0)) = 0;
    fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;    
    % Обрезание по точке максимума, Дополним точкой максимума
    hoF2 = oarr(2);
    ind = find( h < hoF2);
    f_N = [fN(ind), foF2];
    h_N = [h(ind), hoF2];

    % Вызываем внешнюю процедуру, учитывающую E-слой
    [H, SE(i), SF(i)] = Denisenko_Nigth_Nh_from_hv_IRI_real_fmE(f, hv, f_N, h_N, S.timesound, calc_F);
    if calc_F
        set(H, 'Name', [get(H,'Name'), ', fmE = ', num2str(foE(i)), ' MHz']);
    end
end

% Выбираем минимальное SE.
ind_E = find(SE == min(SE));

    [Ne,h,outf,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, foE(ind_E), hmE);

    % Подготавливаем скорректированный профиль
    Ne(find(Ne<0)) = 0;
    fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;    
    % Обрезание по точке максимума, Дополним точкой максимума
    hoF2 = oarr(2);
    ind = find( h < hoF2);
    f_N = [fN(ind), foF2];
    h_N = [h(ind), hoF2];
    
calc_F = 1;
[H, SEx, SFx] = Denisenko_Nigth_Nh_from_hv_IRI_real_fmE(f, hv, f_N, h_N, S.timesound, calc_F);
set(H, 'Name', [get(H,'Name'), ', fmE = ', num2str(foE(ind_E)), ' MHz']);

% Выбираем минимальное SE+SF.
SESF = SE + SF;
ind_E = find(SESF == min(SESF));

    [Ne,h,outf,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, foE(ind_E), hmE);

    % Подготавливаем скорректированный профиль
    Ne(find(Ne<0)) = 0;
    fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;    
    % Обрезание по точке максимума, Дополним точкой максимума
    hoF2 = oarr(2);
    ind = find( h < hoF2);
    f_N = [fN(ind), foF2];
    h_N = [h(ind), hoF2];
    
calc_F = 1;
[H, SEx, SFx] = Denisenko_Nigth_Nh_from_hv_IRI_real_fmE(f, hv, f_N, h_N, S.timesound, calc_F);
set(H, 'Name', [get(H,'Name'), ', fmE = ', num2str(foE(ind_E)), ' MHz. Минимум SE + SF.']);

% Перезапись измененных параметров в файл.
foE_ = foE; % копируем т.к. затирается из МАТ файла
load(last_saved_mat);
foE = foE_(ind_E);
% Получение профиля по IRI
iri_profile = [f_N',h_N'];
save(last_saved_mat, 'timesound', 'index_selected', 'alati', 'along', ...
    'trace_Eo', 'trace_F1o', 'trace_F2o', 'trace_Ex', 'trace_F1x', 'trace_F2x', 'profile', 'iri_profile', ...
    'foF2', 'fxF2', 'hmF2', 'foF1', 'hmF1', 'foE', 'hmE', 'foEs');


% --------------------------------------------------------------------
function itmnu_Save_MAT_and_Work_Cor_fmE_fval_Callback(hObject, eventdata, handles)
% hObject    handle to itmnu_Save_MAT_and_Work_Cor_fmE_fval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global last_saved_mat
global data data_proc 

index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
dat = data.Data(index_selected,:);

alati = dat{3};
along = dat{4};
vbeg = 80; % double(dat{5})/1000.;
vend = double(dat{6})/1000.;
vstp = 5; % double(dat{7})/1000.;
timesound = dat{2};
iyyyy = str2num(timesound(1:4));
mmdd = str2num([timesound(6:7),timesound(9:10)]);

k = strfind(timesound, ':');
if k(1) == 13
    dhour = str2num(timesound(12)) + str2num(timesound(14:15))/60. + str2num(timesound(17:18))/3600.;
else
    dhour = str2num(timesound(12:13)) + str2num(timesound(15:16))/60. + str2num(timesound(18:19))/3600.;
end

% Вначале сохраним MAT-файл для последующей независимой обработки
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% В процессе корректировки изменяется профиль IRI, но это не отображается
% на графике!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
itmnu_Save_Real_to_MAT_Callback();
S = load(last_saved_mat);

% Получим оцифровку
tmp = flipud(S.trace_F2o);
f = tmp(:,1)';
hv = tmp(:,2)';

% Получим подлежащий корректировке IRI профиль
tmp = S.iri_profile;
f_ = tmp(:,1);
for i=1:length(f_)-2 % Поиск локального минимума, соответствующего E-слою
    if f_(i)<=f_(i+1)&&f_(i+1)>f_(i+2)
        fmE = f_(i+1);
        break
    end
end

% Определяем частоты fmE построения корректировок
foE = fmE:0.05:min(f)-0.01;
foF2 = get(findobj('DisplayName','foF2'),'XData');
foF2 = foF2(1);
% Не используются в расчётах.
hmF2 = NaN;
foF1 = NaN;
hmF1 = NaN;
hmE = NaN;

SF1 = zeros(length(foE),1);
SF2 = zeros(length(foE),1);
for i = 1:length(foE)
    [Ne,h,outf,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, foE(i), hmE);

    % Подготавливаем скорректированный профиль
    Ne(find(Ne<0)) = 0;
    fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;    
    % Обрезание по точке максимума, Дополним точкой максимума
    hoF2 = oarr(2);
    ind = find( h < hoF2);
    f_N = [fN(ind), foF2];
    h_N = [h(ind), hoF2];

    % Вызываем внешнюю процедуру, учитывающую E-слой
    [H, SF1(i)] = Denisenko_Nigth_Nh_from_hv_IRI_fval(f, hv, f_N, h_N, S.timesound, 2, 1);
    close(H);
    [H, SF2(i)] = Denisenko_Nigth_Nh_from_hv_IRI_fval(f, hv, f_N, h_N, S.timesound, 3, 2);
    %set(H, 'Name', [get(H,'Name'), ', fmE = ', num2str(foE(i)), ' MHz']);
    close(H);
end

% Выбираем минимальное SF1.
ind_E = find(SF1 == min(SF1));

    [Ne,h,outf,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, foE(ind_E), hmE);

    % Подготавливаем скорректированный профиль
    Ne(find(Ne<0)) = 0;
    fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;    
    % Обрезание по точке максимума, Дополним точкой максимума
    hoF2 = oarr(2);
    ind = find( h < hoF2);
    f_N = [fN(ind), foF2];
    h_N = [h(ind), hoF2];
    
[H, SFx] = Denisenko_Nigth_Nh_from_hv_IRI_fval(f, hv, f_N, h_N, S.timesound, 2, 1);
set(H, 'Name', [get(H,'Name'), ', fmE = ', num2str(foE(ind_E)), ' MHz']);
savefig(H,strcat(last_saved_mat,'m2p1'));

% Выбираем минимальное SF2.
ind_E = find(SF2 == min(SF2));

    [Ne,h,outf,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, foE(ind_E), hmE);

    % Подготавливаем скорректированный профиль
    Ne(find(Ne<0)) = 0;
    fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;    
    % Обрезание по точке максимума, Дополним точкой максимума
    hoF2 = oarr(2);
    ind = find( h < hoF2);
    f_N = [fN(ind), foF2];
    h_N = [h(ind), hoF2];
    
[H, SFx] = Denisenko_Nigth_Nh_from_hv_IRI_fval(f, hv, f_N, h_N, S.timesound, 3, 2);
set(H, 'Name', [get(H,'Name'), ', fmE = ', num2str(foE(ind_E)), ' MHz']);
savefig(H,strcat(last_saved_mat,'m3p2'));

% % Перезапись измененных параметров в файл.
% foE_ = foE; % копируем т.к. затирается из МАТ файла
% load(last_saved_mat);
% foE = foE_(ind_E);
% % Получение профиля по IRI
% iri_profile = [f_N',h_N'];
% save(last_saved_mat, 'timesound', 'index_selected', 'alati', 'along', ...
%     'trace_Eo', 'trace_F1o', 'trace_F2o', 'trace_Ex', 'trace_F1x', 'trace_F2x', 'profile', 'iri_profile', ...
%     'foF2', 'fxF2', 'hmF2', 'foF1', 'hmF1', 'foE', 'hmE', 'foEs');


% --------------------------------------------------------------------
function Denisenko_Nigth_Nh_hv_IRI_parab_a_b_c_end_Callback(hObject, eventdata, handles)
% hObject    handle to Denisenko_Nigth_Nh_hv_IRI_parab_a_b_c_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global last_saved_mat
global data data_proc 

index_selected = int32(get(findobj('Tag','listboxDate'),'Value'));
dat = data.Data(index_selected,:);

alati = dat{3};
along = dat{4};
vbeg = 80; % double(dat{5})/1000.;
vend = double(dat{6})/1000.;
vstp = 5; % double(dat{7})/1000.;
timesound = dat{2};
iyyyy = str2num(timesound(1:4));
mmdd = str2num([timesound(6:7),timesound(9:10)]);

k = strfind(timesound, ':');
if k(1) == 13
    dhour = str2num(timesound(12)) + str2num(timesound(14:15))/60. + str2num(timesound(17:18))/3600.;
else
    dhour = str2num(timesound(12:13)) + str2num(timesound(15:16))/60. + str2num(timesound(18:19))/3600.;
end

% Вначале сохраним MAT-файл для последующей независимой обработки
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% В процессе корректировки изменяется профиль IRI, но это не отображается
% на графике!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
itmnu_Save_Real_to_MAT_Callback();
S = load(last_saved_mat);

% Получим оцифровку
tmp = flipud(S.trace_F2o);
f = tmp(:,1)';
hv = tmp(:,2)';

% Получим подлежащий корректировке IRI профиль
tmp = S.iri_profile;
f_ = tmp(:,1);
for i=1:length(f_)-2 % Поиск локального минимума, соответствующего E-слою
    if f_(i)<=f_(i+1)&&f_(i+1)>f_(i+2)
        fmE = f_(i+1);
        break
    end
end

% Определяем частоты fmE построения корректировок
df = (min(f)-0.01 - fmE)/10;
foE = fmE:df:(min(f)-0.01);
foF2 = get(findobj('DisplayName','foF2'),'XData');
foF2 = foF2(1);
% Не используются в расчётах.
hmF2 = NaN;
foF1 = NaN;
hmF1 = NaN;
hmE = NaN;

if hv(3)<hv(1)
   % nE = 3;
    
    SE = zeros(length(foE),1);
    for i = 1:length(foE)
        [Ne,h,~,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
            foF2, hmF2, foF1, hmF1, foE(i), hmE);
        
        % Подготавливаем скорректированный профиль
        Ne(find(Ne<0)) = 0;
        fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;
        % Обрезание по точке максимума, Дополним точкой максимума
        hoF2 = oarr(2);
        ind = find( h < hoF2);
        f_N = [fN(ind), foF2];
        h_N = [h(ind), hoF2];
        
        % Вызываем внешнюю процедуру, учитывающую E-слой
        [~, SE(i)] = Denisenko_Nigth_Nh_hv_IRI_parab_a_b_c_end(f, hv, f_N, h_N, S.timesound, 0);
        %set(H, 'Name', [get(H,'Name'), ', fmE = ', num2str(foE(i)), ' MHz']);
        %close(H);
    end
    
    % Выбираем минимальное SF.
    ind_E = find(SE == min(SE));
    new_foEm = foE(ind_E(1));
else
    % nE = 1;
    new_foEm = fmE;
end

    [Ne,h,~,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, new_foEm, hmE);

    % Подготавливаем скорректированный профиль
    Ne(find(Ne<0)) = 0;
    fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;    
    % Обрезание по точке максимума, Дополним точкой максимума
    hoF2 = oarr(2);
    ind = find( h < hoF2);
    f_N = [fN(ind), foF2];
    h_N = [h(ind), hoF2];
    
[H, ~, hmF2] = Denisenko_Nigth_Nh_hv_IRI_parab_a_b_c_end(f, hv, f_N, h_N, S.timesound, 1);
% set(H, 'Name', [get(H,'Name'), ', fmE = ', num2str(new_foEm), ' MHz']);
savefig(H,strcat(last_saved_mat,'_a_b_c'));

% Добавление ещё одного скоректированного профиля..
% Корректируем IRI еще и на hmF2
[Ne,h,~,oarr] = iri2016cor(alati,along,iyyyy,mmdd,dhour,vbeg,vend,vstp,...
                    foF2, hmF2, foF1, hmF1, new_foEm, hmE);
                
    % Подготавливаем скорректированный профиль
    Ne(find(Ne<0)) = 0;
    fN = round(sqrt(Ne/(1.24*10^10)) * 100)/100;    
    % Обрезание по точке максимума, Дополним точкой максимума
    hoF2 = oarr(2);
    ind = find( h < hoF2);
    f_N = [fN(ind), foF2];
    h_N = [h(ind), hoF2];
% Получение профиля по IRI
% iri_profile_new = [f_N',h_N'];
plot(f_N',h_N','-g','LineWidth',2,'DisplayName','IRI with hmF2');
legend('Location','southeast')

