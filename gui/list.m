classdef list < matlab.mixin.SetGet
% list - класс окна для представления списка результатов эксперимента.
%
% Created by Alexy Skazik (c)                                   Feb-2018
    
    properties (SetAccess = private)
        Dir; % рабочая папка
        Tag = '1'; % ключ отображаемых файлов
        Title; % заголовок окна
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
        
            obj.Ext = obj.getExtensionsForTag(Tag);
            fileslist = {};
            for i = 1:length(obj.Ext)
                d = dir(fullfile(Dir,'\',obj.Ext{i})); % Список файлов
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
                ME = MException('list:noFiles', ...
                'В папке <%s> отсутствуют соответствующие задаче файлы данных.',Dir);
                throw(ME);
            end
        end
        
        function ext = getExtensionsForTag(hObject, tag)
            switch tag
                case '1'
                    ext = {'*.jpg','*.tif','*.png','*.bmp','*.gif'};
                case '2'
                    ext = {'*.ion'};
                case '3'
                    ext = {'*.frq'};
                case '4'
                    ext = {'*.cat'};
                case '5'
                    ext = {'*.dat'};
                otherwise
                    ME = MException('list:noTag', ...
                        'Tag <%s> не может быть обработан.',tag);
                    throw(ME);
            end

        end
    end
end

