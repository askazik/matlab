classdef main < matlab.mixin.SetGet
% main - ����� ������������ ���� ��� ������ ������������������ ����.
%
% Created by Alexy Skazik (c)                                   Feb-2018
    
    properties (SetAccess = private)
        handles; % ������ ����������� �������� ������
        FileIni = 'ini.mat'; % ����������������� ����
        StructIni; % ����������������� ���������
        FileTag; % ��� ������ ������
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
            
            % Open ini file
            obj.StructIni = lab705_Ini(obj.FileIni);
            
            % ���� ���������� ������
            % TODO: Place menu items (label/callback) into repository.
            mSource = uimenu(fig_handle,'Label','���������');
                uimenu(mSource,'Label','��������� ������� � �� mySQL...');
                uimenu(mSource,'Label','����� � IRI...');
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
                    'Callback', {@(src, event)Close(obj, src, event)},...
                    'Tag','6');
            
            % ���� ������
            mHelp = uimenu(fig_handle,'Label','������');
                uimenu(mHelp,'Label','����������');
                uimenu(mHelp,'Label','� ���������',...
                    'Separator','on',...
                    'Callback', {@(src, event)About(obj, src, event)});
        end 
        
        % ������� ������� �� ������� ����������.
        function Close(hObject, ~, ~)
            tmp = hObject.StructIni;
            PathWork = tmp.PathWork;
            PathIRI = tmp.PathIRI;
            PathHelp = tmp.PathHelp;
            dbName = tmp.dbName;
            password = tmp.password;
            user = tmp.user; 

            % ���������� ini-�����.
            save hObject.FileIni PathWork PathIRI PathHelp ...
                dbName password user '-append';
            
            % ����������� ���� - ����� �� ���������.
            delete(hObject.handles.main)
        end
        
        % ������� ������� �� ������� ����������.
        function OpenDir(hObject, src, ~)
            Dir = hObject.StructIni.PathWork;
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
        
        function ini(hObject)
            % ������ ������������������ �����. ��������� ���������� � ��.
            S = load(hObject.FileIni);
            
            if ~isfield(S,'PathWork') % ���� � ����� � ������� ���������
                S.PathWork = '';
            end
            if ~isfield(S,'PathIRI') % ���� � ����� � ������� IRI
                S.PathIRI = '.\iri';
            end
            if ~isfield(S,'PathHelp') % ���� � ����� � ������� ������
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
            obj.ini();
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
    end
end