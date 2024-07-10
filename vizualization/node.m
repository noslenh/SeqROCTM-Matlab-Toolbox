%Class to model a node of a tree
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2024

classdef node < handle
    properties
        data;
        children;
    end
    
    methods
        %constructor
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
        
        %insert a son at a given position
        function insert_son(obj, n, index)
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
        
        %add a son
        function add_son(obj, n)
            obj.children = [obj.children, n];
        end
        
        %print to console
        function textstr = print(obj)
            nl = numel(obj.children);
            whiteSpaces = repmat(' ', 1, 2*length(obj.data));
            if nl==0
                textstr = [whiteSpaces, '--(', obj.data, ')', '-*', newline]; 
            else
                textstr = [whiteSpaces, '--(', obj.data, ')', newline];
                for i = 1 : nl
                    textstr = [textstr, obj.children(i).print];
                end
            end     
        end
            
    end
end