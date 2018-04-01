classdef db_mysql < matlab.mixin.SetGet
% db_mysql - класс соединения с БД mySQL.
% Если соединение отсутствует попытаемся настроить его.
%
% Created by Alexy Skazik (c)                                   Feb-2018
    
    properties (SetAccess = private)
        IniStructure; % инициализационная структура
        java_connection % объект соединения с БД
    end
    
    %% Object methods
    methods (Access = private)    
        %% Setup mySQL jar
        function setupJar(obj)
            
            dynamic_jars = javaclasspath;
            n_jar = length(dynamic_jars);

            % Check is jar in dynamic path
            i = 1;
            while i <= n_jar && ...
                    ~strcmp(dynamic_jars{i},obj.IniStructure.jar)
                i = i + 1;
            end
            % Add to dynamic path if jar is'nt in dynamic path.
            if i > n_jar
                javaaddpath(obj.IniStructure.jar)
            end            
        end
    end
        
    methods (Access = public)
        
        function obj = db_mysql(IniName)    
        % Construct Database object.
        % IniName - name of initialisation file.
        
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            obj = obj@matlab.mixin.SetGet();
        
            iniobj = db_ini(IniName);
            obj.IniStructure = iniobj.getIniStructureForTag('mysql');
            
            % set mySQL jar
            if ~isfield(obj.IniStructure, 'jar')
                [filename, pathname] = uigetfile('*.jar', ...
                    'Select mySQL-jar file');
                if isequal(filename,0)
                    ME = MException('db_mysql:noJarSelected', ...
                    'You can specify jar file for JDBC.');
                    throw(ME)
                else
                    fname = fullfile(pathname, filename);
                end
                obj.IniStructure.jar = fname;
            end
            
            % set mySQL url
            if ~isfield(obj.IniStructure, 'url')
                answer  = inputdlg(...
                    'Standart is (jdbc:mysql://localhost:3306):',...
                    'Enter mySQL url', ...
                    [1 50],...
                    {'jdbc:mysql://localhost:3306'});
                if isempty(answer)
                    ME = MException('db_mysql:noUrlSelected', ...
                    'You can specify url for JDBC.');
                    throw(ME)
                else
                    obj.IniStructure.url = answer{:};
                end
            end
            
            % set user
            if ~isfield(obj.IniStructure, 'user')
                answer  = inputdlg(...
                    'Standart is (root):',...
                    'Enter mySQL username', ...
                    [1 25],...
                    {'root'});
                if isempty(answer)
                    ME = MException('db_mysql:noUserSelected', ...
                    'You can specify user for JDBC.');
                    throw(ME)
                else
                    obj.IniStructure.user = answer{:};
                end
            end

            % set password
            if ~isfield(obj.IniStructure, 'pwd')
                answer  = inputdlg(...
                    '',...
                    'Enter mySQL password for username', ...
                    [1 25]);
                if isempty(answer)
                    ME = MException('db_mysql:noPwdSelected', ...
                    'You can specify pwd for JDBC.');
                    throw(ME)
                else
                    obj.IniStructure.pwd = answer{:};
                end
            end
            
            % select schema(dbName)
            if ~isfield(obj.IniStructure, 'dbname')
                schema = chooseSchemaDialog(obj);
                if isempty(schema)
                    return
                else
                    obj.IniStructure.dbname = schema;
                end
            end
            obj.IniStructure = orderfields(obj.IniStructure);
            iniobj.setIniStructureForTag('mysql',obj.IniStructure);
            
            obj.setupJar();
        end
        
        %% Open JDBC Driver with current USER/PASSWORD paar.
        function conn = openConnection(obj)
        % Username and password you chose when installing mySQL
            obj.java_connection = [];
            S = obj.IniStructure;
            try
                props=java.util.Properties;
                props.setProperty('user', S.user);
                props.setProperty('password', S.pwd);

                % Create the database connection
                driver = com.mysql.jdbc.Driver;
                url = S.url;
                obj.java_connection = driver.connect(url, props);
                conn = obj.java_connection;

            catch ME
                disp(ME)
                rethrow(ME);
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
        
        function delete(obj)
        % Destructor
            if ~isempty(obj.java_connection)
                obj.java_connection.close();
                disp('mySQL was disconnected from JDBC.')
            else
                disp('mySQL was not connected with JDBC.')
            end
        end
    end
end