classdef Domain
    
    properties
        nodes
        edges
    end
    
    methods
        
        % Builder
        function obj=Domain(input1,input2)
            w = whos('input1');
            switch w.class               
                case 'char'
                    switch input1
                        case 'square'
                            if nargin == 2
                                r = input2;
                            else
                                r = 1;
                            end
                            input1 = [-r -r;r -r;r r;-r r];
                            input2 = {{1,2};{2,3};{3,4};{4,1}};
                            obj = Domain(input1,input2);
                            
                            
                        case 'disc'
                            if nargin == 2
                                r = input2;
                            else
                                r = 1;
                            end
                            input1 = [r 0; -r 0];
                            input2 = {{1,2,[num2str(r) '*cos(t)'],[num2str(r) '*sin(t)'],0,pi},{2,1,[num2str(r) '*cos(t)'],[num2str(r) '*sin(t)'],pi,2*pi}};
                            obj = Domain(input1,input2);
                    end
                      
                case 'double'
                    switch nargin
                        case 1
                            obj.nodes = input1;
                            N =  size(input1,1);
                            edge = cell(N,1);
                            for i = 1 : N-1
                                edge{i} = {i,i+1};
                            end
                            edge{N} = {N,1};
                            obj.edges = edge;
                            
                        case 2
                            obj.nodes = input1;
                            obj.edges = input2;                       
                    end
                    
                otherwise
                    error('input1 class must be double or char')
                     obj.nodes = [];
                     obj.edges = [];   
            end
                    
    end
        
        % Methods
        function disp(obj)
            % DOMAIN.DISP:
            %       obj.disp display informations contained in the Domain
            %       object.
            disp('Geometry object')
            disp('Nodes :')
            disp(obj.nodes)
            for i = 1 : length(obj.edges)
                edge = obj.edges{i};
                a = num2str(edge{1});
                b = num2str(edge{2});
                switch length(edge)
                    case {2,3}
                        disp(['Edge ' num2str(i) ' : from ' a ' to ' b])
                    case {6,7}
                        expx = edge{3};
                        expy = edge{4};
                        tmin = num2str(edge{5});
                        tmax = num2str(edge{6});
                        disp(['Edge ' num2str(i) ' : from ' a ' to ' b ', on the curve (' expx ',' expy '),  t' char(1028) '[' tmin ',' tmax ']' ])
                end
                
            end
            disp(' ')
        end
        function plot(obj,I)
            
            switch nargin
                case 1
                    obj.plot(1:length(obj.edges));
                case 2
                    Ne = length(I);
                    for i = 1 : Ne
                        obj.plotedge(obj.edges{i});
                        hold on
                    end
                    axis image
                    hold off
            end
        end
        
        % Tools
        function b = eq(obj,domain)
            if size(obj.nodes,1) == size(domain.nodes,1)
            s = sum(sum(abs(obj.nodes-domain.nodes)));
            else 
                s = 1;
            end
            b = ~(s>0);
            e = obj.edges;
            if b
                for i = 1 : length(e)
                    a = (e{i}{1} ~= domain.edges{i}{1}) || (e{i}{2} ~= domain.edges{i}{2});
                    if a
                        b=0;
                        break
                    end
                end
            end
        end            
        function plotedge(obj,edge)
            le = length(edge);
            
            switch le
                case {3,7}
                    switch edge{length(edge)}
                        case 'l'
                            c = 'r';
                        case 'r'
                            c = 'b';
                            
                        case {'rl','lr'}
                            c = 'g';
                    end
                otherwise
                    c = 'r';
            end
            switch le
                
                case {2,3}
                    x1 = obj.nodes(edge{1},:);
                    x2 = obj.nodes(edge{2},:);
                    plot([x1(1) x2(1)],[x1(2) x2(2)],c,'linewidth',2);
                case {6,7}
                    N = 1000;
                    expx = edge{3};
                    expy = edge{4};
                    tmin = edge{5};
                    tmax = edge{6};
                    t = linspace(tmin,tmax,N);
                    x = eval(expx);
                    y = eval(expy);
                    h = (1e-10)*(tmax-tmin);
                    t = t+h;
                    xh = eval(expx);
                    yh = eval(expy);
                    
                    xp = (xh - x)/h;
                    yp = (yh - y)/h;
                    
                    v = sqrt(xp.^2+yp.^2);
                    v = (v(2:end) + v(1:end-1)/2);
                    
                    cv =[0 cumsum(1./v)];
                    cv = cv/max(cv);
                    t = cv*(tmax-tmin)+tmin;
                    
                    x = eval(expx);
                    y = eval(expy);
                    
                    
                    plot(x,y,c,'linewidth',2)
            end
            
            
            
        end
        function [M,Ibound]=matrix(obj,dx)
            
            % DOMAIN.MATRIX:
            %       [M,Ibound]=obj.matrix compile the data contained in the Domain
            %       object in a matrix M to be used by initmesh. Ibound is
            %       the number of the initial edge.
            E = obj.edges;
            M = [];
            Ibound = [];
            for i = 1 : length(E)
                edge = E{i};
                [x,y] = obj.discretize_edge(edge,dx);
                N = length(x);
                switch length(edge)
                    case {2,6}
                        m = [2*ones(1,N-1) ; x(1:end-1) ; x(2:end) ;  y(1:end-1) ; y(2:end) ; ones(1,N-1) ; zeros(6,N-1)];
                    case {3,7}
                        switch edge{length(edge)}
                            case 'l'
                                side = [ones(1,N-1) ; zeros(1,N-1)];
                            case 'r'
                                side = [zeros(1,N-1) ; ones(1,N-1)];
                                
                            case {'rl','lr'}
                                side = ones(2,N-1);
                        end
                        m = [2*ones(1,N-1) ; x(1:end-1) ; x(2:end) ;  y(1:end-1) ; y(2:end) ; side ; zeros(5,N-1)];
                        
                        
                end
                M = [M m];
                Ibound = [Ibound ; i*ones(N-1,1)];
            end
        end  
        function [x,y] = discretize_edge(obj,edge,dx)
            le = length(edge);
            switch le
                case {2,3}
                    x1 = obj.nodes(edge{1},:);
                    x2 = obj.nodes(edge{2},:);
                    l = sqrt((x2(1)-x1(1)).^2 + (x2(2)-x1(2)).^2);
                    N = max(floor(l/dx)+1,2);
                    x = linspace(x1(1),x2(1),N);
                    y = linspace(x1(2),x2(2),N);
                case {6,7}
                    N = 1000;
                    expx = edge{3};
                    expy = edge{4};
                    tmin = edge{5};
                    tmax = edge{6};
                    t = linspace(tmin,tmax,N);
                    x = eval(expx);
                    y = eval(expy);
                    h = (1e-10)*(tmax-tmin);
                    t = t+h;
                    xh = eval(expx);
                    yh = eval(expy);
                    
                    xp = (xh - x)/h;
                    yp = (yh - y)/h;
                    
                    v = sqrt(xp.^2+yp.^2);
                    v = (v(2:end) + v(1:end-1))/2;
                    l = sum(v)/(N-1)*(tmax-tmin);
                    N = max(floor(l/dx)+2,3);
                    t = linspace(tmin,tmax,N);
                    x = eval(expx);
                    y = eval(expy);
                    h = (1e-10)*(tmax-tmin);
                    t = t+h;
                    xh = eval(expx);
                    yh = eval(expy);
                    
                    xp = (xh - x)/h;
                    yp = (yh - y)/h;
                    
                    v = sqrt(xp.^2+yp.^2);
                    v = (v(2:end) + v(1:end-1)/2);
                    cv =[0 cumsum(1./v)];
                    cv = cv/max(cv);
                    t = cv*(tmax-tmin)+tmin;
                    
                    x = eval(expx);
                    y = eval(expy);
            end
            
            
        end
    end
end