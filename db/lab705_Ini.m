classdef lab705_Ini < matlab.mixin.SetGet
% lab705_Ini - класс инициализационного файла.
% Если файл отсутствует попытаемся корректно создать его.
%
% Created by Alexy Skazik (c)                                   Feb-2018
    
    properties (SetAccess=private, AbortSet=true)
        IniName; % имя инициализационного файла
        IniStructure % параметры инициализации
    end
    
    %% Object methods
    methods (Access = public)    
        
        function obj = lab705_Ini(fname)    
        % Construct CheckDatabase object with <ini.mat> initialisation
        % file by default.
        % fname - name of initialisation file.
        
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            obj = obj@matlab.mixin.SetGet();
                       
            % Check of ini file existence
            try
                S = load(fname);
            catch ME
                % Select/create/exit              
                % Construct a questdlg with three options
                choice = questdlg(...
                    ME.message, ...
                    'Error', ...
                    'Select file','Create file','Exit', ...
                    'Select file');
                % Handle response
                [pathstr,name,~] = fileparts(fname);
                switch choice
                    case 'Select file'
                        [filename, pathname] = uigetfile('*.mat', ...
                            'Select a MATLAB MAT file');
                        if isequal(filename,0)
                            return
                        else
                            fname = fullfile(pathname, filename);
                        end
                    case 'Create file'
                        jar = 'C:\Program Files (x86)\MySQL\Connector.J 5.1\mysql-connector-java-5.1.42-bin.jar';
                        url = 'jdbc:mysql://localhost:3306';
                        user ='root';
                        password = '2760977';
                        dbName = 'ionosphere';
                        
                        save(fullfile(pathstr, name),...
                            'jar','url','user','password','dbName');
                    case 'Exit'
                        return
                end
            end
            obj.IniName = fname;
            obj.IniStructure = load(obj.IniName);
        end
        
        function disp(obj)    
        % Display of object parameters
            disp(['The initialization file: ',obj.IniName]); % имя инициализационного файла
            disp(obj.IniStructure); % параметры инициализации
        end
        
        function save(obj, in_structure, fname)    
        % Save of object parameters
            if ~exist('fname','var')
                fname = obj.IniName; 
            end
            
            if exist('in_structure','var')
                obj.IniStructure = in_structure; 
                names = fieldnames(obj.IniStructure);
                for i = 1:length(names)
                    tmp_field = names{i};
                    eval([tmp_field, ' = ''', obj.IniStructure.(tmp_field),'''']);

                    % Create file if not exist
                    if i == 1
                        save(fname, tmp_field);
                    end
                    save(fname, tmp_field, '-append');
                end
            end
        end
        
    end
    
    %% Properties set/get
    methods
        function obj = set.IniStructure(obj,val)
            obj.IniStructure = val;
        end
        
        function value = get.IniStructure(obj)
            value = obj.IniStructure;
        end
        
        function value = get.IniName(obj)
            value = obj.IniName;
        end
    end
    
end

