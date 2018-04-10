classdef db_ini < matlab.mixin.SetGet
% db_ini - класс инициализационного файла для подключения к базе данных.
% Если файл отсутствует попытаемся корректно создать его.
%
% Created by Alexei Skazik (c)                                   Feb-2018
    
    properties (SetAccess=private, AbortSet=true)
        IniName; % имя инициализационного файла
        IniStructure % параметры инициализации
        Tag % database tag ('mysql', 'oracle',...)
    end
    
    %% Object methods
    methods (Access = public)    
        
        function obj = db_ini(fname, tag)    
        % Construct CheckDatabase object with <ini_sources.mat> initialisation
        % file by default.
        % fname - name of initialization file,
        % tag - database tag ('mysql', 'oracle',...).
        
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            obj = obj@matlab.mixin.SetGet();
            
            % Check database tag
            if ~exist('tag','var')|| strcmp(tag,'') || isempty(tag)
                tag = 'my';
            end
            
            % Check of ini file existence
            is_exist = (exist(fname,'file') == 2);
            if ~is_exist
                % Select/create/exit              
                % Construct a questdlg with three options
                choice = questdlg(...
                    strcat('File <',fname,'> not found!'), ...
                    'Error', ...
                    'Select file','Create file by default [mySQL]','Exit', ...
                    'Select file');
                % Handle response
                [pathname,~,~] = fileparts(fname);
                switch choice
                    case 'Select file'
                        [filename, pathname] = uigetfile('*.mat', ...
                            'Select a MATLAB MAT file',...
                            pathname);
                        if isequal(filename,0)
                            ME = MException('db_ini:noInitDBFile', ...
                            'The configuration file for DB was not selected.');
                            throw(ME)
                        else
                            fname = fullfile(pathname, filename);
                        end
                    case 'Create file by default [mySQL]'
                        mysql.jar = '.\jar\mysql-connector-java-5.1.42-bin.jar';
                        mysql.url = 'jdbc:mysql://localhost:3306';
                        mysql.user ='root';
                        mysql.pwd = '2760977';
                        mysql.dbname = 'ionosphere';
                        fname = fullfile(pathname, fname);
                        save(fname,tag);
                    case 'Exit'
                        ME = MException('db_ini:noInitDBFile', ...
                        'The configuration file for DB was not selected.');
                        throw(ME)
                end
            end
            obj.IniName = fname;
            obj.Tag = tag;
            obj.IniStructure = load(obj.IniName,obj.Tag);
        end
        
        function disp(obj)    
        % Display of object parameters
            first = {'File';'Tag'};
            last = {obj.IniName;obj.Tag};
            names = fieldnames(obj.IniStructure.(obj.Tag));
            for i = 1:length(names)
                value = obj.IniStructure.(obj.Tag).(names{i});
                first(i+2) = names(i);
                last{i+2} = value;
            end
            s = strcat(first, {': '}, {sprintf('\t')}, last);
            disp(s)
        end
        
        function save(obj)    
        % Save of object parameters
            eval([obj.Tag, ' = obj.IniStructure.', obj.Tag]);
            save(obj.IniName, obj.Tag);
        end
        
        function setIniStructureForTag(obj,tag,val)
            obj.IniStructure.(tag) = val;
            obj.save();
        end
        
        function value = getIniStructureForTag(obj,tag)
            value = obj.IniStructure.(tag);
        end
        
    end
    
    %% Properties set/get
    methods       
        function value = get.IniName(obj)
            value = obj.IniName;
        end
        
        function obj = set.IniName(obj,val)
            obj.IniName=val;
        end
    end
    
end

