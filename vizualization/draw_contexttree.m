function draw_contexttree(contexts, A, varargin)
%DRAW_CONTEXTTREE Draws a context tree
%
% Inputs
%
%   contexts    : set of contexts
%   A           : Alphabet
%   varargin(1) : color for drawing the tree
%   varargin(2) : height that the figure will reserve for the design
%
%Author : Noslen Hernandez (noslen.hernandez-gonzalez@inrae.fr), Aline Duarte (alineduarte@usp.br)
%Date   : 02/2019

% convert the set of contexts to a tree class
tree = contexts_to_tree(contexts, A);

% create an object of the class to draw
dt = DrawTree(tree);
buchheim(dt);

switch length(varargin)
    case 1
        dt.do_graphic(varargin{1});
    case 2
        dt.do_graphic(varargin{1}, varargin{2});
    case 0
        dt.do_graphic();
end

end