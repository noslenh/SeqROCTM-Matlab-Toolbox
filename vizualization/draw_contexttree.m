function draw_contexttree(contexts, alphabet, varargin)
%DRAW_CONTEXTTREE draw the tree defined by the contexts
%
% Inputs
%
% contexts    : set of contexts
% Alphabet    : Alphabet
% varargin(1) : color for drawing the tree
% varargin(2) : height that the figure will reserve for the design

%Author : Aline Duarte, Noslen Hernandez
%Date   : 02/2019


% convert the set of contexts to a tree class
tree = contexts_to_tree(contexts, alphabet);

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