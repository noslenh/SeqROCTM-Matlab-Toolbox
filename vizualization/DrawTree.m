classdef DrawTree < handle
    properties
        x;
        y;
        tree;
        children;
        parent;
        thread;
        offset;
        pancestor;
        change; shift;
        lmost_sibling;
        number;
        mod;
        
        max_x;
        height;
    end
    
    methods
        function obj = DrawTree(tree, varargin) % parent, depth, number
            
            % only want 3 optional inputs at most
            numvarargs = length(varargin);
            if numvarargs > 3
                error('DrawTree:TooManyInputs', ...
                    'requires at most 3 optional inputs');
            end
            
            % set defaults for optional inputs (parent, depth, number)
            optargs = {[], 0, 1};
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [obj.parent, depth, obj.number] = optargs{:};
            
            obj.x = -1;
            obj.y = depth;
            obj.tree = tree;
            for i = 1 : numel(tree.children)
                obj.children = [obj.children, ...
                    DrawTree(tree.children(i), obj, depth + 0.5, i)]; % aqui estoy poniendo 0.5 entre las Y
            end
            
            obj.thread = [];
            obj.offset = 0;
            obj.pancestor = obj;
            obj.lmost_sibling = [];
            obj.change = 0;
            obj.shift = 0;
            obj.mod = 0;
        end
        
        function l = left(obj)
            if numel(obj.children) > 0
                l = obj.children(1);
            else
               l = obj.thread;
            end
        end
        
        function l = right(obj)
            nc = numel(obj.children);
            if nc > 0
                l = obj.children(nc);
            else
               l = obj.thread;
            end
        end
        
        function n = left_brother(obj)
            n = [];
            if ~isempty(obj.parent)
                for ch = obj.parent.children
                    if ch == obj
                        return;
                    else
                        n = ch;
                    end
                end
            end
        end
        
        function n = get_lmostsibling(obj)
            if isempty(obj.lmost_sibling)&&(~isempty(obj.parent))&&(~isequal(obj,obj.parent.children(1)))
                obj.lmost_sibling = obj.parent.children(1);
            end
            n = obj.lmost_sibling;
        end
        
        function dt = buchheim(dt)
            firstwalk(dt);
            [min, max, h] = second_walk(dt);
            if min < 0 
                third_walk(dt, -min);
            end
            dt.max_x = max;
            dt.height = h;
        end
        
        function third_walk(dt, m)
            dt.x = dt.x + m;
            for c = dt.children
                third_walk(c, m);
            end
        end
        
        function firstwalk(v, varargin) % distance
            % only want 1 optional inputs at most
            numvarargs = length(varargin);
            if numvarargs > 1
                error('DrawTree:firstwalk:TooManyInputs', ...
                    'requires at most 1 optional inputs');
            end
            % set defaults for optional inputs (distance)
            optargs = {1.2};
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            distance = optargs{1};
            
            if isempty(v.children)
                if ~isempty(v.get_lmostsibling)
                    v.x = v.left_brother().x + distance;
                else
                    v.x = 0;
                end
            else
                default_ancestor = v.children(1);
                for w = v.children
                    firstwalk(w);
                    default_ancestor = apportion(w, default_ancestor, distance);
                end
                execute_shifts(v);
                
                midpoint = (v.children(1).x + v.children(numel(v.children)).x) / 2;

                ell = v.children(1);
                arr = v.children(numel(v.children));
                w = v.left_brother();
                if ~isempty(w)
                    v.x = w.x + distance;
                    v.mod = v.x - midpoint;
                else
                    v.x = midpoint;
                end
            end
        end
        
        function default_ancestor = apportion(v, default_ancestor, distance)
            w = v.left_brother();
            if ~isempty(w) 
                % in buchheim notation
                % i == inner; o == outer; r == right; l == left;
                vir = v; vor = v;
                vil = w;
                vol = v.get_lmostsibling();
                sir = v.mod; sor = v.mod;
                sil = vil.mod;
                sol = vol.mod;
                while ~isempty(vil.right())&&(~isempty(vir.left()))
                    vil = vil.right();
                    vir = vir.left();
                    vol = vol.left();
                    vor = vor.right();
                    vor.pancestor = v;
                    shift1 = (vil.x + sil) - (vir.x + sir) + distance;
                    if shift1 > 0
                        a = ancestor(vil, v, default_ancestor);
                        move_subtree(a, v, shift1);
                        sir = sir + shift1;
                        sor = sor + shift1;
                    end
                    sil = sil + vil.mod;
                    sir = sir + vir.mod;
                    sol = sol + vol.mod;
                    sor = sor + vor.mod;
                end
                if ~isempty(vil.right()) && isempty(vor.right())
                    vor.thread = vil.right();
                    vor.mod = vor.mod + sil - sor;
                else
                    if ~isempty(vir.left()) && isempty(vol.left())
                        vol.thread = vir.left();
                        vol.mod = vol.mod + sir - sol;
                    end
                    default_ancestor = v;
                end
            end
        end
        
        function move_subtree(wl, wr, shift)
            subtrees = wr.number - wl.number;
            wr.change = wr.change - shift / subtrees;
            wr.shift = wr.shift + shift;
            wl.change = wl.change + shift / subtrees;
            wr.x = wr.x + shift;
            wr.mod = wr.mod + shift;
        end
        
        function execute_shifts(v)
            shift1 = 0; change1 = 0;
            nc = numel(v.children);
            for ii = nc:-1:1
                w = v.children(ii);
                w.x = w.x + shift1;
                w.mod = w.mod + shift1;
                change1 = change1 + w.change;
                shift1 = w.shift + change1;
            end
        end
        
        function an = ancestor(vil, v, default_ancestor)
            found = false;
            ii = 1;
            while ~found && (ii <= numel(v.parent.children))  %vil.ancestor in v.parent.children
                found = (v.parent.children(ii) == vil.pancestor);
                ii = ii+1;
            end
            if found
                an = vil.pancestor;
            else
                an = default_ancestor;
            end
        end

        function [min, max, h] = second_walk(v, varargin) % m = 0, depth = 0, min = [], max = [], h = []
            
            numvarargs = length(varargin);
            if numvarargs > 5
                error('DrawTree:secondwalk:TooManyInputs', ...
                    'requires at most 3 optional inputs');
            end
            % set defaults for optional inputs (m, depth)
            optargs = {0, 0, [], [], []};
            optargs(1:numvarargs) = varargin;
            % Place optional args in memorable variable names
            [m, depth, min, max, h] = optargs{:};

            v.x = v.x + m;
            v.y = depth;
            
            if isempty(min) || v.x < min
                min = v.x;
            end
            
            % to known the grater x-coord and the height
            if isempty(max) || v.x > max
                max = v.x;
            end
            
            if isempty(h) || v.y > h
                h = v.y;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for w = v.children
                [min, max, h] = second_walk(w, m + v.mod, depth+0.5, min, max, h);
            end
        end
        
        function do_plot(obj)
            for c = obj.children
                plot([obj.x, c.x], [-obj.y, -c.y]);
                do_plot(c);
            end
            if isempty(obj.children)
                txt = mat2str(obj.tree.data);
                th = text(obj.x, -obj.y, txt(2:end-1)); 
                pp = get(th,'Extent');
                set(th,'position',[pp(1)-pp(3)/2, pp(2), 0]);
            end
        end
        
        function do_graphic(obj, varargin)
            hold on;

            switch length(varargin)
                case 0  %default
                    ylim([-obj.height - 0.3 0.3]);
                    axis off;
                    set(gca,'colororder',[1 0 0]);
                case 1  % set the specified color
                    ylim([-obj.height - 0.3 0.3]);
                    axis off;
                    set(gca,'colororder',varargin{1});
                case 2 % set the specified color and height
                    ylim([-varargin{2}*0.5-0.3 0.3]);
                    axis off;
                    set(gca,'colororder',varargin{1});
            end
            do_plot(obj);
        end
        
    end  
end