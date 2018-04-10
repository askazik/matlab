classdef main < matlab.mixin.SetGet
% main - ����� ������������ ���� ��� ������ ������������������ ����.
%
% Created by Alexei Skazik (c)                                   Feb-2018
    
    properties (SetAccess = private)
        handles; % ������ ����������� �������� ������
        IniGUI = 'ini_gui.mat'; % ����������������� ����
        IniSources = 'ini_sources.mat'; % ����������������� ���� ����������
        StructIniGUI; % ����������������� ���������
        FileTag; % ��� ������ ������
        
        % Default GUI
        DefaultGUI = struct(...
            'PathWork', '',...
            'PathIRI', '.\iri',...
            'PathHelp', '.\help');
    end
    
    methods (Access = protected)
        
        function gui(hObject)
        % ������������ ���� � ����.
        
            hObject.handles.main = figure('MenuBar', 'None', ...
                'Name', '��������� ����������������� ������', ...
                'NumberTitle', 'off',...
                'Resize','off',...
                'DockControls','off');
            set(hObject.handles.main, ...
                'CloseRequestFcn', ...
                {@(eventdata, handles)Close(hObject, eventdata, handles)});
        end
        
        function CreateMenuBar(obj)
        % ������������ ���� � ����.
            fig_handle = obj.handles.main;
            
            % ���� ���������� ������
            % TODO: Place menu items (label/callback) into repository.
            mSource = uimenu(fig_handle,'Label','���������');
                uimenu(mSource,'Label','���������� ������������...');
                uimenu(mSource,'Label','����� � IRI...',...
                    'Callback', {@(src, event)OpenIRIDir(obj, src, event)});
                uimenu(mSource,'Label','���������� ��� (mySQL)',...
                    'Separator','on','Tag','0');
                uimenu(mSource,'Label','���������� (��������)...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','1');
                uimenu(mSource,'Label','���������� ���� (.ion)...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','2');
                uimenu(mSource,'Label','��������� ���� (.frq)...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','3');
                uimenu(mSource,'Label','������ RWM...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','4');
                uimenu(mSource,'Label','������ �������� ������������...',...
                    'Callback', {@(src, event)OpenDir(obj, src, event)},...
                    'Tag','5');
                uimenu(mSource,'Label','�����',...
                    'Separator','on','Accelerator','Q',...
                    'Callback', {@(src, event)Close(obj, src, event)});
            
            % ���� ������
            mHelp = uimenu(fig_handle,'Label','������');
                uimenu(mHelp,'Label','����������');
                uimenu(mHelp,'Label','� ���������',...
                    'Separator','on',...
                    'Callback', {@(src, event)About(obj, src, event)});
        end 
        
        % ������� ������� �� ������� ����������.
        function Close(hObject, ~, ~)
            hObject.save_gui_ini(hObject.IniGUI);
            children = hObject.handles.children;
            for i=1:length(children)
                delete(children(i))
            end
            delete(hObject);
        end
               
        % ������� ������� �� ������� ����������.
        function OpenDir(hObject, src, ~)
            Dir = hObject.StructIniGUI.PathWork;
            Tag = get(src,'Tag');
            Title = get(src,'Label');
            
            Dir = uigetdir(Dir,Title);
            if Dir ~= 0
                if ~isfield(hObject.handles,'children') % �������� ����
                    % ������ �������� ����
                    hObject.handles.children = list(Dir,Tag,Title);
                else
                    hObject.handles.children(end+1) = list(Dir,Tag,Title);
                end
            end
        end
        
        function OpenIRIDir(hObject, ~, ~)
            Dir = hObject.StructIniGUI.PathIRI;
            Title = 'Choose directory of the IRI files...';
            Dir = uigetdir(Dir,Title);
            if ~isequal(Dir,0)
                hObject.StructIniGUI.PathIRI = Dir;
            end
        end        
        
        % Create default GUI init file.
        function set_default_gui(hObject)
            S = hObject.DefaultGUI;
            fields = fieldnames(S);
            for i = 1:numel(fields)
                hObject.StructIniGUI.(fields{i}) = S.(fields{i});
            end
            msg = strcat('Set default GUI init parameters.');
            warning(msg)
        end
        
        % Verify init file structure.
        function verify_gui_ini(hObject)
            S = hObject.DefaultGUI;
            S_cur = hObject.StructIniGUI;
            fields = fieldnames(S);
            for i = 1:numel(fields)
                if ~isfield(S_cur, fields{i})
                    hObject.StructIniGUI.(fields{i}) = S.(fields{i});
                    msg = strcat('Create and fill not existed init GUI field. <',...
                        fields{i}, '> = ', S.(fields{i}));
                    warning(msg)
                end
            end
        end
        
        % Save GUI structure.
        function save_gui_ini(hObject, fname)
            S = orderfields(hObject.StructIniGUI);
            fields = fieldnames(S);
            for i = 1:numel(fields)
                evalc([fields{i}, ' = ''', S.(fields{i}),'''']);
                if i == 1
                    save(fname, fields{i});
                else
                    save(fname, fields{i},'-append');
                end
            end
        end
        
        function ini_gui(hObject)
            % Open GUI init file if exist.
            try
                hObject.StructIniGUI = load(hObject.IniGUI);
            catch ME
                [filename, pathname] = uigetfile('*.mat', ME.message);
                if isequal(filename,0)
                    hObject.set_default_gui()
                else
                    hObject.IniGUI = fullfile(pathname, filename);
                end
            end
            hObject.verify_gui_ini();
            hObject.save_gui_ini(hObject.IniGUI);
        end
        
        function About(~, ~, ~)
            % ������� ���������� � ���������.
            
            text = {'��������� ������������� ��� ������������ � ���������',...
                '����������� ������ ����������� ���������,',...
                '���������� � ����� ������ � ��������� ��������',...
                '��������.', ...
                '', ...
                '������ �.�. (2018)'};
            hAbout = msgbox(text, ...
                '� ���������', ...
                'help', ...
                'modal');
            uiwait(hAbout);
        end
    end
       
    methods (Access = public)
        
        function obj = main()
        % �����������
            %% Object Initialization %%
            % Call superclass constructor before accessing object
            % You cannot conditionalize this statement
            obj = obj@matlab.mixin.SetGet();
            
            %% Post Initialization %%
            obj.ini_gui();
            obj.gui();
            obj.CreateMenuBar();
            
            % ��������� �������� ����.
            OutPos = get(obj.handles.main,'OuterPosition');
            InPos = get(obj.handles.main,'InnerPosition');
            pos = [30,...
                OutPos(2)+OutPos(4),...
                OutPos(3),OutPos(4)-InPos(4)];
            set(obj.handles.main,'OuterPosition', pos);
        end
        
        function delete(obj)
            %remove the closerequestfcn from the figure, this prevents an
            %infitie loop with the following delete command
            set(obj.handles.main,  'closerequestfcn', '');
            % ����������� ���� - ����� �� ���������.
            delete(obj.handles.main);
            %clear out the pointer to the figure - prevents memory leaks
            obj.handles = [];
                        
            disp('The successful completion of the program.')
        end
    end
end