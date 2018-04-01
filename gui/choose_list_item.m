classdef choose_list_item < handle
%CHOOSE_LIST_ITEM Draw dialog with possibility of an item choose.
%   Detailed explanation goes here
%
% Created by Alexy Skazik (c)                                   Apr-2018

    properties (SetAccess=private)
        d % dialog handle
        txt % label handle
        popup % popup list handle
        btn % button handle
        items_list % items list of string cells
        choosed_item
        choosed_item_number
        Label = 'Choose an item from popup list:';
    end
    
    methods (Access = public)
        
        function choice = getChoice(obj)
            choice = obj.choosed_item;
        end
        
        function obj = choose_list_item(list_of_cells, preferred_item, label)
            
            obj.items_list = list_of_cells;
            
            obj.d = dialog('Position',[300 300 220 150],...
                'Name','Choose dialog...');
            obj.txt = uicontrol('Parent',obj.d,...
                'Style','text',...
                'Position',[10 80 210 40],...
                'String',obj.Label);
            obj.popup = uicontrol('Parent',obj.d,...
                'Style','popup',...
                'Position',[10 70 200 25],...
                'String',obj.items_list,...
                'Callback',@(src, event)popup_callback(obj, src, event));
            obj.btn = uicontrol('Parent',obj.d,...
                'Position',[89 20 70 25],...
                'String','Close',...
                'Callback',@(src, event)close(obj, src, event));
            obj.choosed_item_number = 1; % by default
            
            % Check preferred_item
            if exist('preferred_item','var') && ...
                    ~strcmp(preferred_item,'') && ...
                    ~isempty(preferred_item)
                obj.setPreferredItem(preferred_item);
            end
            
            % Check label
            if exist('label','var') && ...
                    ~strcmp(label,'') && ...
                    ~isempty(label)
                set(obj.txt,'String',label);
            end                   
            
            % Wait for d to close before running to completion
            uiwait(obj.d);
        end
    end
            
    methods (Access = private)       
        %class deconstructor - handles the cleaning up of the class &
        %figure. Either the class or the figure can initiate the closing
        %condition, this function makes sure both are cleaned up
        function delete(obj)
            %delete the figure
            delete(obj.d);
            %clear out the pointer to the figure - prevents memory leaks
            obj.d = [];
        end
        
        function obj = close(obj, src, event)        
            obj = popup_callback(obj, src, event);
            delete(gcf);
        end
        
        function obj = popup_callback(obj, ~, ~)
            idx = obj.popup.Value;
            popup_items = obj.popup.String;
            % This code uses dot notation to get properties.
            % Dot notation runs in R2014b and later.
            % For R2014a and earlier:
            % idx = get(popup,'Value');
            % popup_items = get(popup,'String');
            obj.choosed_item = char(popup_items(idx,:));
            obj.choosed_item_number = idx;
        end

        function setPreferredItem(obj, preferred_item)
            popup_items = obj.popup.String;
            for idx = 1:length(popup_items)
                choice = char(popup_items(idx,:));
                if strcmp(choice,preferred_item)
                    obj.popup.Value = idx;
                    break;
                end
            end
        end
    end
    
    %% Properties set/get
    methods              

    end
    
end

