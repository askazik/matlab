classdef lab705_MySQL < matlab.mixin.SetGet
% lab705_MySQL - класс соединения с БД mySQL.
% Если соединение отсутствует попытаемся настроить его.
%
% Created by Alexy Skazik (c)                                   Feb-2018
    
    properties (SetAccess = private)
        IniObject; % инициализационный объект
        java_connection % объект соединения с БД
    end
    
    methods (Access = public)
        
        %% Object Initialization %%
        function obj = lab705_MySQL(fname)    
        % Construct CheckDatabase object with <ini.mat> initialisation
        % file by default.
        % fname - name of initialisation file.
        
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            obj = obj@matlab.mixin.SetGet();
        
            % Set of ini file name
            if ~exist('fname','var')
                fname = 'ini.mat'; 
            end
            
            % Check of ini file existence
            obj.IniObject = lab705_Ini(fname);
            S = obj.IniObject.IniStructure;
            
            % set mySQL jar
            if ~isfield(S, 'jar')
                [filename, pathname] = uigetfile('*.jar', ...
                    'Select mySQL-jar file');
                if isequal(filename,0)
                    return
                else
                    fname = fullfile(pathname, filename);
                end
                S.jar = fname;
            end
            
            % set mySQL url
            if ~isfield(S, 'url')
                answer  = inputdlg(...
                    'Standart is (jdbc:mysql://localhost:3306):',...
                    'Enter mySQL url', ...
                    [1 50],...
                    {'jdbc:mysql://localhost:3306'});
                if isempty(answer)
                    return
                else
                    S.url = answer{:};
                end
            end
            
            % set user
            if ~isfield(S, 'user')
                answer  = inputdlg(...
                    'Standart is (root):',...
                    'Enter mySQL username', ...
                    [1 25],...
                    {'root'});
                if isempty(answer)
                    return
                else
                    S.user = answer{:};
                end
            end

            % set password
            if ~isfield(S, 'password')
                answer  = inputdlg(...
                    '',...
                    'Enter mySQL password for username', ...
                    [1 25]);
                if isempty(answer)
                    return
                else
                    S.password = answer{:};
                end
            end
            
            % select schema(dbName)
            if ~isfield(S, 'dbName')
                schema = chooseSchemaDialog(obj);
                if isempty(schema)
                    return
                else
                    S.dbName = schema;
                end
            end
            s = orderfields(S);
            obj.IniObject.save(s);
        end
        
        %% Setup mySQL jar
        function setupJar(obj)
            
            dynamic_jars = javaclasspath;
            n_jar = length(dynamic_jars);

            % Check is jar in dynamic path
            i = 1;
            while ~(i <= n_jar && ...
                    strcmp(dynamic_jars{i},obj.IniObject.IniStructure.jar))
                i = i + 1;
            end
            % Add to dynamic path if jar is'nt in dynamic path.
            if i > n_jar
                javaaddpath(obj.IniObject.IniStructure.jar)
            end            
        end
        
        %% Open JDBC Driver with current USER/PASSWORD paar.
        function conn = openConnection(obj)
        % Username and password you chose when installing mySQL
            conn = [];
            S = obj.IniObject.IniStructure;
            try
                props=java.util.Properties;
                props.setProperty('user', S.user);
                props.setProperty('password', S.password);

                % Create the database connection
                driver = com.mysql.jdbc.Driver;
                url = S.url;
                conn = driver.connect(url, props);

            catch ME
                disp(ME.ExceptionObject)
            end
        end
        
        %% List existing schemas into mySQL.
        function schemas = getSchemas(obj)
            schemas = {};
            conn = obj.java_connection;
            if ~isempty(conn)
                sql = 'show databases'; % Gets all records

                % Statement 
                stmt = conn.createStatement();
                % ResultSet 
                rs = stmt.executeQuery(sql);

                % Get rows count into ResultSet
                count = -1;
                if rs.last()
                   count = rs.getRow();
                end

                % Read the results into an array of result structs
                rs.beforeFirst();
                schemas = cell(count,1);

                i=1;
                while rs.next()
                    schemas{i} = (char(rs.getString(1)));
                    i=i+1;
                end
            end
        end
        
        %% Open choice schema dialog.
        function choice = chooseSchemaDialog(obj)
            
            % Get all schemas in current database
            setupJar(obj);
            obj.java_connection = openConnection(obj);
            schemas = getSchemas(obj);

            d = dialog('Position',[300 300 220 150],'Name','Schema...');
            txt = uicontrol('Parent',d,...
                   'Style','text',...
                   'Position',[10 80 210 40],...
                   'String','Select a database(schema)');

            popup = uicontrol('Parent',d,...
                   'Style','popup',...
                   'Position',[10 70 200 25],...
                   'String',schemas,...
                   'Callback',@popup_callback);

            btn = uicontrol('Parent',d,...
                   'Position',[89 20 70 25],...
                   'String','Close',...
                   'Callback','delete(gcf)');

            choice = [];

            % Wait for d to close before running to completion
            uiwait(d);

               function popup_callback(popup,event)
                  idx = popup.Value;
                  popup_items = popup.String;
                  % This code uses dot notation to get properties.
                  % Dot notation runs in R2014b and later.
                  % For R2014a and earlier:
                  % idx = get(popup,'Value');
                  % popup_items = get(popup,'String');
                  choice = char(popup_items(idx,:));
               end
        end
    end
end