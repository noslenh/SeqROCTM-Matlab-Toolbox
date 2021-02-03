%Class to model a node of a tree

%Author : Noslen Hernandez (noslenh@gmail.com), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2021

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
            
    end
end