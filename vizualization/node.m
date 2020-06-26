classdef node < handle
    properties
        data;
        children;
    end
    
    methods
        function obj = node(data, varargin)
            if nargin > 0
                obj.data = data;
                
                numvarargs = length(varargin);
                if numvarargs == 1
                    obj.children = varargin{1};
                else
                    obj.children = [];
                end 
            end
        end
        
        function insert_child(obj, n, index)
            nl = numel(obj.children);
            new_children(nl+1) = node();
            
            for ii = 1 : index-1
                new_children(ii) = obj.children(ii);
            end
            new_children(index) = n;
            
            for ii = index : nl
                new_children(ii+1) = obj.children(ii);
            end
            obj.children = new_children;
        end
        
        function add_child(obj, n)
            obj.children = [obj.children, n];
        end
            
    end
end