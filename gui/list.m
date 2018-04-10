classdef list < matlab.mixin.SetGet
% list - класс окна для представления списка результатов эксперимента.
%
% Created by Alexei Skazik (c)                                   Feb-2018
    
    properties (SetAccess = private)
        handles; % хэндлы графических объектов класса
        Dir; % рабочая папка
        Tag = '1'; % ключ отображаемых файлов
        Title; % заголовок окна
        FilesList = {}; % список файлов для работы
        
        % Files extensions for current Tag
        ExtForTag = struct(...
            'f1', {'*.jpg','*.tif','*.png','*.bmp','*.gif'},...
            'f2', '*.ion',...
            'f3', '*.frq',...
            'f4', '*.cat',...
            'f5', '*.dat');
    end
    
    properties (SetAccess = protected)
        Ext; % расширения отображаемых файлов
    end
    
    methods (Access = public)
        
        function obj = list(Dir, Tag, Title)
        % Конструктор
            %% Object Initialization %%
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            obj = obj@matlab.mixin.SetGet();
        
            obj.Dir = Dir; % рабочая папка
            obj.Tag = Tag; % ключ отображаемых файлов
            obj.Title = Title; % заголовок окна
        
            obj.retriveFilesList();
            obj.gui();
        end
        
        function delete(obj)
            %remove the closerequestfcn from the figure, this prevents an
            %infitie loop with the following delete command
            set(obj.handles.main,  'closerequestfcn', '');
            % Уничтожение окна - выход из программы.
            delete(obj.handles.main);
            %clear out the pointer to the figure - prevents memory leaks
            obj.handles = [];
        end
    end
        
    methods (Access = protected)    
        
        function gui(obj)
            % Отрисовываем окно.     
            obj.handles.main = figure('MenuBar', 'None', ...
                'Name', obj.Title, ...
                'NumberTitle', 'off',...
                'Resize','on',...
                'DockControls','off');
            set(obj.handles.main, ...
                'CloseRequestFcn', ...
                {@(eventdata, handles)Close(obj, eventdata, handles)});
                      
            % Установка размеров окна.
            r = groot;
            sc = r.ScreenSize();
            dx = sc(3)/4;
            pos = [sc(3) - dx - 10,...
                50, dx, sc(4) - 50];
            set(obj.handles.main,'OuterPosition', pos);
            
            % Стандартные
%            OutPos = get(obj.handles.main,'OuterPosition');
            InPos = get(obj.handles.main,'InnerPosition');
%             pos = [30,...
%                 OutPos(2)+OutPos(4),...
%                 OutPos(3),OutPos(4)-InPos(4)];
%             set(obj.handles.main,'OuterPosition', pos);
            
            % uicontrols
            % Add a text uicontrol to label the slider.
            dy = 10;
            dx = InPos(3);
            y = dy;
            obj.handles.btn = uicontrol('Style','pushbutton',...
                'Position',[10 y dx-20 30],...
                'String',obj.Dir,...
                'Callback', ...
                {@(eventdata, handles)BtnCallback(obj, eventdata, handles)});
            % Create pop-up menu
            y = y + 30 + dy;
            obj.handles.listbox = uicontrol('Style', 'listbox',...
                'String', obj.FilesList,...
                'Position', ...
                [10 y dx-20 InPos(4) - y - 10], ...
                'Callback', ...
                {@(eventdata, handles)ListboxCallback(obj, eventdata, handles)});
    
            % Make figure visble after adding all components
            obj.handles.main.Visible = 'on';
        end
              
        % Функции реакции на события интерфейса.
        function Close(hObject, ~, ~)
            % hObject.save_gui_ini(hObject.IniGUI)            
            % Уничтожение окна - выход из программы.
            delete(hObject)
        end
        
        function ListboxCallback(obj, ~, ~)
            
            listbox_items = obj.handles.listbox.String;
            idx = obj.handles.listbox.Value;
            fname = listbox_items(idx,:);
            
            fullname = fullfile(obj.Dir, fname);
            disp(fullname)
        end
        
        function BtnCallback(obj, ~, ~)

            title = 'Choose new working directory...';
            selDir = uigetdir(obj.Dir,title);
            if ~isequal(selDir,0)
                obj.Dir = selDir;
            end
            obj.retriveFilesList();
            obj.handles.listbox.Value = 1;
            obj.handles.listbox.String = obj.FilesList;
            obj.handles.btn.String = obj.Dir;
        end           
        
        function retriveFilesList(obj)
            
            obj.Ext = obj.getExtensionsForTag(obj.Tag);
            fileslist = {};
            for i = 1:length(obj.Ext)
                d = dir(fullfile(obj.Dir,'\',obj.Ext{i})); % Список файлов
                if ~isempty(d)
                    for j = 1:length(d)
                        if isempty(fileslist)
                            fileslist{1} = d(j).name;
                        else
                            fileslist{end+1} = d(j).name;
                        end
                    end
                end
            end
            
            if ~isempty(fileslist)
                set(findobj('Tag','listbox'), 'String', fileslist);
            else
                msg = ['В папке <',obj.Dir,...
                    '> отсутствуют соответствующие задаче файлы данных.'];
                msgbox(msg, 'Error','error');
            end
            obj.FilesList = fileslist;
        end
        
        function ext = getExtensionsForTag(hObject, tag)
            field = ['f',tag];
            if isfield(hObject.ExtForTag, field)
                n = length(hObject.ExtForTag);
                for i=1:n
                    ext{i} = hObject.ExtForTag(i).(field);
                end
            else
                ext = '*.*';
                msg = strcat('Set extension <',ext,'> for Tag <',tag,'>.');
                warning(msg)
            end

        end
    end
end

