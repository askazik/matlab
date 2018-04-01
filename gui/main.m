classdef main < matlab.mixin.SetGet
% main - класс управляющего окна для вызова специализированных окон.
%
% Created by Alexy Skazik (c)                                   Feb-2018
    
    properties (SetAccess = private)
        handles; % хэндлы графических объектов класса
        FileIni = 'ini.mat'; % инициализационный файл
        StructIni; % инициализационная структура
        FileTag; % тэг файлов данных
    end
    
    methods (Access = protected)
        
        function gui(hObject)
        % Отрисовываем окно с меню.
        
            hObject.handles.main = figure('MenuBar', 'None', ...
                'Name', 'Обработка экспериментальных данных', ...
                'NumberTitle', 'off',...
                'Resize','off',...
                'DockControls','off');
            set(hObject.handles.main, ...
                'CloseRequestFcn', ...
                {@(eventdata, handles)Close(hObject, eventdata, handles)});
        end
        
        function CreateMenuBar(obj)
        % Отрисовываем окно с меню.
            fig_handle = obj.handles.main;
            
            % Open ini file
            obj.StructIni = lab705_Ini(obj.FileIni);
            
            % Меню источников данных
            % TODO: Place menu items (label/callback) into repository.
            mSource = uimenu(fig_handle,'Label','Источники');
                uimenu(mSource,'Label','Настройка доступа к БД mySQL...');
                uimenu(mSource,'Label','Папка с IRI...');
                uimenu(mSource,'Label','Ионограммы ИПГ (mySQL)',...
                    'Separator','on','Tag','0');
                uimenu(mSource,'Label','Ионограммы (картинки)...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','1');
                uimenu(mSource,'Label','Ионограммы НИИФ (.ion)...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','2');
                uimenu(mSource,'Label','Амплитуды НИИФ (.frq)...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','3');
                uimenu(mSource,'Label','Записи RWM...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','4');
                uimenu(mSource,'Label','Записи внешнего зондирования...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','5');
                uimenu(mSource,'Label','Выход',...
                    'Separator','on','Accelerator','Q',...
                    'Callback', {@(src, event)Close(obj, src, event)},...
                    'Tag','6');
            
            % Меню помощи
            mHelp = uimenu(fig_handle,'Label','Помощь');
                uimenu(mHelp,'Label','Оглавление');
                uimenu(mHelp,'Label','О программе',...
                    'Separator','on',...
                    'Callback', {@(src, event)About(obj, src, event)});
        end 
        
        % Функции реакции на события интерфейса.
        function Close(hObject, ~, ~)
            tmp = hObject.StructIni;
            PathWork = tmp.PathWork;
            PathIRI = tmp.PathIRI;
            PathHelp = tmp.PathHelp;
            dbName = tmp.dbName;
            password = tmp.password;
            user = tmp.user; 

            % Сохранение ini-файла.
            save hObject.FileIni PathWork PathIRI PathHelp ...
                dbName password user '-append';
            
            % Уничтожение окна - выход из программы.
            delete(hObject.handles.main)
        end
        
        % Функции реакции на события интерфейса.
        function OpenDir(hObject, src, ~)
            Dir = hObject.StructIni.PathWork;
            Tag = get(src,'Tag');
            Title = get(src,'Label');
            
            Dir = uigetdir(Dir,Title);
            if Dir ~= 0
                if ~isfield(hObject.handles,'children') % дочерние окна
                    % первое дочернее окно
                    hObject.handles.children = list(Dir,Tag,Title);
                else
                    hObject.handles.children(end+1) = list(Dir,Tag,Title);
                end
            end
        end
        
        function ini(hObject)
            % Чтение инициализационного файла. Настройка соединения с БД.
            S = load(hObject.FileIni);
            
            if ~isfield(S,'PathWork') % путь к папке с файлами измерений
                S.PathWork = '';
            end
            if ~isfield(S,'PathIRI') % путь к папке с файлами IRI
                S.PathIRI = '.\iri';
            end
            if ~isfield(S,'PathHelp') % путь к папке с файлами помощи
                S.PathHelp = '.\help';
            end
            if ~isfield(S,'dbName')
                S.dbName = 'ionosphere';
            end
            if ~isfield(S,'password')
                S.password = '27670977';
            end
            if ~isfield(S,'user')
                S.user = 'root';
            end
            
            hObject.StructIni = S;
        end
        
        function About(~, ~, ~)
            % Краткая информация о программе.
            
            text = {'Программа предназначена для визуализации и обработки',...
                'разнородных данных ионосферных измерений,',...
                'сохранённых в базах данных и различных файловых',...
                'форматах.', ...
                '', ...
                'Сказик А.И. (2018)'};
            hAbout = msgbox(text, ...
                'О программе', ...
                'help', ...
                'modal');
            uiwait(hAbout);
        end
    end
       
    methods (Access = public)
        
        function obj = main()
        % Конструктор
            %% Object Initialization %%
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            obj = obj@matlab.mixin.SetGet();
            
            %% Post Initialization %%
            obj.ini();
            obj.gui();
            obj.CreateMenuBar();
            
            % Установка размеров окна.
            OutPos = get(obj.handles.main,'OuterPosition');
            InPos = get(obj.handles.main,'InnerPosition');
            pos = [30,...
                OutPos(2)+OutPos(4),...
                OutPos(3),OutPos(4)-InPos(4)];
            set(obj.handles.main,'OuterPosition', pos);
        end
    end
end