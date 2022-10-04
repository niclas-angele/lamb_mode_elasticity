classdef Mesh < handle
    properties (SetAccess = private)
        domain
        nodes
        triangles
        edges
        tricenters

        Iboundnodes
        Iintnodes
        
        Nnodes
        Ntriangles
        Nedges
           
        int1
        intx
        intxx
        basis
        stats
        
        hexa
        Ihexacenters
        Nhexa
    end
    methods (Access = public)
        
        
        % Constructors
        function obj = Mesh(varargin)
            % Class MESH (this help is obsolete)
            %
            % Constructor :
            %
            % obj = Mesh
            % Create a defaut test mesh object of the unit square with 2 triangles.
            %
            % obj = Mesh(dx)
            % Create a mesh of the unit square of resolution dx.
            % 
            % obj = Mesh(p,dx)
            % Create a mesh from a polygonal simply connected domain given
            % - p is (2xN) double.
            %
            % obj = Mesh(domain,dx)
            % Create a mesh from a domain object (see Domain)
            %
            % obj = Mesh(p,dx,'save') or obj = Mesh(domain,dx,'save') asks
            % Matlab to save and load the object to build it only one time
            % in case of repetitive call.
            %
            % METHOD LIST :
            %
            %
            % obj.disp
            % obj.plot
            %
            % u = obj.P0(input,order)
            % u = obj.P1(input,order)
            % u = obj.P0bound(input,order)
            % z = obj.P1togrid(u,x,y)
            % u = obj.gridtoP1(u_grid,x,y)
            % uP0 = obj.P1toP0(uP1)
            % p = obj.P0prod(u,v)
            % n = obj.normL(u)
            % v = obj.norm(u,p)
            %
            % surf(u)
            % contour(u)
            % quiver(u)
            %
            % gu = obj.grad(u)
            % gsu = obj.grads(u)
            % divu = obj.div(u)
            %
            % u = obj.solvesystem(A,F)
            %
            % u = obj.solve(a,c,f,bclist)
            % u = obj.solvepde(a,c,f,bclist)
            % K = obj.stiffness(a)
            % M = obj.mass(c)
            % F = obj.source(f)
            % Fb = obj.boundsource(E,f)
            % [Mb,Fb,Id,Vd] = obj.bcmatrices(bclist)
            % [A,F,Id,Vd] = obj.pde(a,c,f,bclist)
            % [A,F] = obj.system(a,c,f,bclist)
            %
            % u = obj.solve_vect(a,c,f,bclist)
            % u = obj.solvepde_vect(a,c,f,bclist)
            % K = obj.stiffness_vect(a)
            % A = obj.mass_vect(a)
            % F = obj.source_vect(f)
            % F = obj.boundsource_vect(E,f)
            % [Mb,Fb,Id,Vd] = obj.bcmatrices_vect(bclist)
            % [A,F,Id,Vd] = obj.pde_vect(a,c,f,bclist)
            % [A,F] = obj.system_vect(a,c,f,bclist)
            %
            % E = obj.buildbasis
            % [int1,intx,intxx] = obj.integrals
            % [A,F] = obj.assemble(A,F,Id,Vd)
            % u = obj.eval(input,x,y,order)
            % m = obj.displace(u)
            if ~license('test','pde_toolbox')
                error('Mesh class needs the pde_toolbox')
            end
           
            
            % Mesh;
            if nargin == 0
                D = Domain([-1 -1;1 -1;1 1;-1 1]);
                    obj = Mesh(D,0.9);
                    return
            end
            
            arg1 = varargin{1};
            arg1origin = arg1;
            
            % Mesh('square',...) shortcuts
            if ischar(arg1)
                switch arg1
                    case {'square','Square'}
                        arg1 = Domain([-1 -1;1 -1;1 1;-1 1]);
                    case {'disc','Disc'} % Mesh('disc',dx)
                        arg1 = Domain('disc');
                end
            elseif isnumeric(arg1) && isscalar(arg1) % Mesh(dx,...)
                D = Domain([-1 -1;1 -1;1 1;-1 1]);
                cellarg = [{D} varargin];
                obj = Mesh(cellarg{:});
                return
            elseif isnumeric(arg1) && (size(arg1,1) == 2 || size(arg1,2) == 2) % Mesh(P,...)
                arg1 = Domain(arg1);
            end
            
            
            switch nargin
                case 1 % Mesh(Domain),
                    obj = Mesh(arg1,0.9);
                    return
                case 2
                    arg2 = varargin{2};
                    if isscalar(arg2) && isnumeric(arg2) % Effective builder Mesh(D,dx)
                        disp('buiding mesh v2.4_beta')
                        time = cputime;
                        [gm,Ibound]=arg1.matrix(arg2);
                        [p,e,t]=initmesh(gm,'Hmax',arg2,'Hgrad',1.5,'Jiggle','mean','JiggleIter',10);
                        p=p';
                        e=e';
                        e = [e(:,[1 2]) Ibound(e(:,5))];
                        t=t';
                        dx = arg2;
                        
                        obj.domain = arg1;
                        obj.nodes=p;
                        obj.edges=e;
                        obj.triangles=t(:,1:3);
                        obj.meshpropbuilder(dx);
                        obj.stats(8)= cputime-time;
                        disp 'mesh done'
                         
                    else
                        disp(arg2);
                        error('Wrong argument')                       
                    end
                        
                        
                case 3
                    arg2 = varargin{2};
                    arg3 = varargin{3};
                    if ischar(arg3)
                        switch arg3
                            case 'save'
                                a = dir('mesh*.mat');
                                N = length(a);
                                for i = 1 : N
                                    load(a(i).name,'saved_dx');
                                    load(a(i).name,'saved_domain');
                                    if saved_dx == arg2 && saved_domain == arg1
                                        disp('loading existing mesh')
                                        load(a(i).name,'saved_mesh');
                                        obj = saved_mesh;
                                        return
                                    end
                                end
                                obj = Mesh(arg1,arg2);
                                saved_mesh = obj;
                                saved_dx = arg2;
                                saved_domain = arg1;
                                name = ['mesh' num2str(N+1) '.mat' ];
                                save(name,'saved_mesh','saved_dx','saved_domain')
                                return
                            case 'hexa' % Mesh(P,dx,'hexa') Hexagonal Mesh
                                disp('buildind hexagonal mesh')
                                time = cputime;
                                [p,e,t,hex] = hexagonalmesh(obj,arg1origin,arg2);
                                
                                obj = Mesh(p,e,t);
                                obj.hexa = hex;
                                obj.Nhexa = size(hex,1);
                                obj.Ihexacenters = obj.triangles( hex(:,1),3);
                                obj.stats(8) =  cputime-time;
                        end
                    else % Mesh(p,e,t)
                            time = cputime;
                            p = arg1origin;
                            e = arg2;
                            t = arg3;
                            dx = -1;
                            obj.domain = [];
                            obj.nodes = p;
                            obj.edges = e;
                            obj.triangles = t;
                            obj.meshpropbuilder(dx);
                            obj.stats(8)= cputime-time;
                    end
            end
        end
        function m = displace(obj,u)
            if ~obj.isP1(u,1)
                error('Needs P1 vectorial function')
            end
            p = obj.nodes + u;
            e = obj.edges;
            t = obj.triangles;
            disp 'builing displaced Mesh'
            m = Mesh(p,e,t);
            disp 'mesh done'
            m.domain = 'displaced';
   
        end
        function u = solve(obj,a,c,f,bclist)
            if nargin == 3
                A = a;
                F = c;
            else
                [A,F] = obj.system(a,c,f,bclist);
            end
            u = A\F;
        end
        
        
        % Display
        function disp(obj)
            % MESH.DISP
            %       obj.disp displays data of the Mesh object.
            disp('Mesh object')
            disp(['Nodes :            ' num2str(obj.Nnodes)])
            disp(['Triangles :        ' num2str(obj.Ntriangles)])
            disp(['dx :               ' num2str(obj.stats(1))])
            disp(['Min resolution :   ' num2str(obj.stats(2))])
            disp(['Max resolution :   ' num2str(obj.stats(3))])
            disp(['Mean resolution :  ' num2str(obj.stats(4))])
            disp(['Uniformity :       ' num2str(obj.stats(5))])
            disp(['Mean quality :     ' num2str(obj.stats(6))])
            disp(['Memory used(kB) :  ' num2str(obj.stats(7))])
            disp(['Time(s) :          ' num2str(obj.stats(8))])
            disp(' ')
        end
        function plot(obj)
            % MESH.PLOT
            %       obj.plot plots the Mesh object.
            t = obj.triangles;
            t = [t ones(size(t,1),1)];
            pdemesh(obj.nodes',obj.edges',t');
            axis image
        end
        
        
        % Function spaces
        function u = P0(obj,input,order)
            % P0(mesh)
            % P0(mesh,input)
            % P0(mesh,input,order)
            
            switch nargin
                case 1
                    u = obj.P0(0,0);
                    return
                case 2
                    u = obj.P0(input,-1);
                    return
            end
            
            
            switch class(input)
                case  'double'
                    Nn = obj.Ntriangles;
                    [Nl,Nc] = size(input);
                    switch Nl
                        case 1
                            val = ones(Nn,1)*input;
                        case 2
                            val  = [input(1,:) input(2,:)];
                        case Nn
                            val = input;
                        otherwise
                            if Nc == 1
                                val = input.';
                                u = obj.P0(mesh,val,order);
                                return
                            else
                                error('Wrong row number in input')
                            end
                    end
                    u = obj.tensor(val,order);
                    
                case {'char','cell','function_handle'}
                    val = obj.eval(input,'P0');
                    u = obj.tensor(val,order);
                    return
                    
                otherwise
                    class(input)
                    error('Wrong input class')
            end
            
        end
        function u = P1(obj,input,order)
            % P1(mesh)
            % P1(mesh,input)
            % P1(mesh,input,order)
            
            switch nargin
                case 1
                    u = obj.P1(0,0);
                    return
                case 2
                    u = obj.P1(input,-1);
                    return
            end
            
            
            switch class(input)
                case  'double'
                    Nn = obj.Nnodes;
                    [Nl,Nc] = size(input);
                    switch Nl
                        case 1
                            val = ones(Nn,1)*input;
                        case 2
                            val  = [input(1,:) input(2,:)];
                        case Nn
                            val = input;
                        otherwise
                            if Nc == 1
                                val = input.';
                                u = obj.P1(mesh,val,order);
                                return
                            else
                                error('Wrong row number in input')
                            end
                    end
                    u = obj.tensor(val,order);
                    
                case {'char','cell','function_handle'}
                    val = obj.eval(input,'P1');
                    u = obj.tensor(val,order);
                    return
                    
                otherwise
                    class(input)
                    error('Wrong input class')
            end
            
        end
        function u = P0bound(obj,input,order)
            % P0bound(mesh)
            % P0bound(mesh,input)
            % P0bound(mesh,input,order)
            
            switch nargin
                case 1
                    u = obj.P0bound(0,0);
                    return
                case 2
                    u = obj.P0bound(input,-1);
                    return
            end
            
            
            switch class(input)
                case  'double'
                    e1 = obj.edges(:,1);
                    e2 = obj.edges(:,2);
                    x = (obj.nodes(e1,1) + obj.nodes(e2,1))/2;
                    y = (obj.nodes(e1,2) + obj.nodes(e2,2))/2;
                    Nn = length(x);
                    [Nl,Nc] = size(input);
                    switch Nl
                        case 1
                            val = ones(Nn,1)*input;
                        case 2
                            val  = [input(1,:) input(2,:)];
                        case Nn
                            val = input;
                        otherwise
                            if Nc == 1
                                val = input.';
                                u = obj.P0bound(mesh,val,order);
                                return
                            else
                                error('Wrong row number in input')
                            end
                    end
                    u = obj.tensor(val,order);
                    
                case {'char','cell','function_handle'}
                    val = obj.eval(input,'P0bound');
                    u = obj.tensor(val,order);
                    return
                    
                otherwise
                    class(input)
                    error('Wrong input class')
            end
            
        end
        function b = isP0(obj,u,order)
            b = true;
            [Nl,Nc] = size(u);
            Nt = obj.Ntriangles;
            if Nl ~= Nt
                b = false;
                return
            end
            if nargin == 3 && Nc ~= 2^order
                b = false;
            end
        end
        function b = isP1(obj,u,order)
            b = true;
            [Nl,Nc] = size(u);
            Nt = obj.Nnodes;
            if Nl ~= Nt
                b = false;
                return
            end
            if nargin == 3 && Nc ~= 2^order
                b = false;
            end
        end
          
        
        % Converters
        function varargout = P1togrid(obj,u,x,y)
             if ~obj.isP1(u)
                 error('Needs P1 input Function')
             end
             N = size(u,2);
             varargout = cell(1,N);
             for n = 1 : N
                varargout{n} = tri2grid(obj.nodes',obj.triangles',u(:,n),x,y);
             end
        end
           
        
        % Operators
        function gu = grad(obj,u)
            if ~obj.isP1(u)
                error('grad needs P1 Function')
            end
            Nc = size(u,2);
            
            Nt = obj.Ntriangles;
            T = obj.triangles;
            a1 = obj.basis(:,1:2);
            a2 = obj.basis(:,4:5);
            a3 = obj.basis(:,7:8);
            g = zeros(Nt,2);
            gu = zeros(Nt,2*Nc);
            for k = 1 : Nc
                for n = 1 : Nt
                    t = T(n,:);
                    g(n,:) = u(t,k).'*[a1(n,:) ; a2(n,:) ; a3(n,:)];
                end
                gu(:,[k Nc+k]) = g;
            end 
        end
        function gsu = strain(obj,u)
             if ~obj.isP1(u,1)
                 error('gradsym needs P1 vector Function')
             end
             gu = obj.grad(u);
             s = (gu(:,2)+gu(:,3))/2;
             gsu = [gu(:,1) s s gu(:,4)];
        end       
        function nu = norm(obj,u,p)
             if nargin == 2
                 p = 2;
             end
             au = abs(u);
             if strcmp(p,'inf')
                 nu = max(au,[],2);
             else
                 nu = (sum(au.^p,2)).^(1/p);
             end
        end
        function detu = det(obj,u)
            Nc = size(u,2);
            if Nc ~= 4
                error('Needs an order 2 tensor Function')
            end
            detu = u(:,1).*u(:,4) - u(:,2).*u(:,3);
            
        end
        
        
        % Function plots
        function surf(obj,u)
            Nc = size(u,2);
            if isreal(u)
                if Nc == 1
                    obj.surfrs(u);
                else
                    N1 = floor(sqrt(Nc));
                    N2 = ceil(Nc/N1);
                    for n = 1 : Nc
                        subplot(N1,N2,n)
                        obj.surfrs(u(:,n));
                        title(n)
                    end
                end
            else
                if Nc == 1
                    subplot 121
                    obj.surfrs(real(u));
                    title('Re(u)')
                    subplot 122
                    obj.surfrs(imag(u));
                    title('Im(u)')
                else
                    for n = 1 : Nc
                        subplot(Nc,2,2*n-1)
                        obj.surfrs(real(u(:,n)));
                        title(['Re(u_{' num2str(n) '})'])
                        subplot(Nc,2,2*n)
                        obj.surfrs(imag(u(:,n)));
                        title(['Im(u_{' num2str(n) '})'])
                    end
                end
            end
        end    
        function contour(obj,u)
            Nc = size(u,2);
            if isreal(u)
                if Nc == 1
                    obj.contourrs(u);
                else
                    N1 = floor(sqrt(Nc));
                    N2 = ceil(Nc/N1);
                    for n = 1 : Nc
                        subplot(N1,N2,n)
                        obj.contourrs(u(:,n));
                        title(n)
                    end
                end
            else
                if Nc == 1
                    subplot 121
                    obj.surfrs(real(u));
                    title('Re(u)')
                    subplot 122
                    obj.surfrs(imag(u));
                    title('Im(u)')
                else
                    for n = 1 : Nc
                        subplot(Nc,2,2*n-1)
                        obj.surfrs(real(u(:,n)));
                        title(['Re(u_{' num2str(n) '})'])
                        subplot(Nc,2,2*n)
                        obj.surfrs(imag(u(:,n)));
                        title(['Im(u_{' num2str(n) '})'])
                    end
                end
            end
        end
        function quiver(obj,u,varargin)
            if obj.isP0(u,1)
                x = obj.tricenters(:,1);
                y = obj.tricenters(:,2);
            elseif obj.isP1(u,1)
                x = obj.nodes(:,1);
                y = obj.nodes(:,2);
            else
                error('Needs a P0 or P1 vector function')
            end
            
            varargin = [{x,y,u(:,1),u(:,2)} varargin];
            quiver(varargin{:});
        end
        function streamline(obj,u,xstart,ystart,varargin)
            if obj.isP0(u,1)
                x = obj.tricenters(:,1);
                y = obj.tricenters(:,2);
            elseif obj.isP1(u,1)
                x = obj.nodes(:,1);
                y = obj.nodes(:,2);
            else
                error('Needs a P0 or P1 vector function')
            end
            dx = obj.mesh.stats(1)/2;
            xmin = min(x);
            xmax = max(x);
            ymin = min(y);
            ymax = max(y);
            Nx = floor((xmax-xmin)/dx)+1;
            Ny = floor((ymax-ymin)/dx)+1;
            
            X = linspace(xmin,xmax,Nx);
            Y = linspace(ymin,ymax,Ny);
            [Xgrid,Ygrid] = meshgrid(X,Y);
            
            [u,v] = obj.grid(X,Y);
            
            varargin = [{Xgrid,Ygrid,u,v,xstart,ystart} varargin];
            streamline(varargin{:});
        end
        
 
        % Tools
        function val = eval(obj,input,space) 
            switch space
                case 'P1'
                    N = obj.Nnodes;
                    x = obj.nodes(:,1);
                    y = obj.nodes(:,2);
                case 'P0'
                    N = obj.Ntriangles;
                    x = obj.tricenters(:,1);
                    y = obj.tricenters(:,2);
                case 'P0bound'
                    e1 = obj.edges(:,1);
                    e2 = obj.edges(:,2);
                    x = (obj.nodes(e1,1) + obj.nodes(e2,1))/2;
                    y = (obj.nodes(e1,2) + obj.nodes(e2,2))/2;
                    N = length(x);
                otherwise
                    error('Space not treated')
            end

            switch class(input)
                case 'double'
                    Nl = size(input,1);
                    switch Nl
                        case 1   
                            val = input*ones(N,1);
                        case N
                            val = input;
                        otherwise
                            error('Wrong input size')
                    end 
                case 'function_handle'
                    val = feval(input,x,y);
                case 'char'
                    val = eval(input);
                case 'cell'
                    Ncell = length(input);
                    val = zeros(N,Ncell);
                    for n = 1 : Ncell
                        val(:,n) = obj.eval(input{n},space);
                    end 
                otherwise
                    error('input type not recognized')
                    
            end
            
        end
        
        % Submesh
        function m = submesh(obj)
            time = cputime;
            m = Mesh;
            m.domain = obj.domain;
            
            Nn = obj.Nnodes;
            Nt = obj.Ntriangles;
            Ne = obj.Nedges;
            p = obj.nodes;
            t = obj.triangles;
            e = obj.edges;
            dx = obj.stats(1)/2;
            
            % Find Edges numbering
            I = zeros(3*Nt,1);
            J = I;
            V = ones(3*Nt,1);
            for n = 1 : Nt
                ip1 = t(n,1);
                ip2 = t(n,2);
                ip3 = t(n,3);
                I(3*n-2:3*n) = [min(ip1,ip3);min(ip2,ip1);min(ip3,ip2)];
                J(3*n-2:3*n) = [max(ip1,ip3);max(ip2,ip1);max(ip3,ip2)];
            end
            
            mat = sparse(I,J,V,Nn,Nn);
            mat = (mat>0);
            nie = sum(mat(:));
            
            mat2 = sparse(Nn,Nn);
            mat2(mat) = 1:nie;
            mat = mat2+mat2.';
            
            % Create new nodes and new triangles
            newp = [p ; zeros(nie,2)];
            newt = zeros(4*Nt,3);

            for n = 1 : Nt
                ip1 = t(n,1);
                ip2 = t(n,2);
                ip3 = t(n,3);
                
                p12 = (p(ip1,:) + p(ip2,:))/2;
                p23 = (p(ip2,:) + p(ip3,:))/2;
                p31 = (p(ip3,:) + p(ip1,:))/2;
                
                ip12 = Nn + mat(ip1,ip2);
                newp(ip12,:) = p12;
                ip23 = Nn + mat(ip2,ip3);
                newp(ip23,:) = p23;
                ip31 = Nn + mat(ip3,ip1);
                newp(ip31,:) = p31;
                
                newt(4*n-3:4*n,:)  = [ip1 ip12 ip31
                    ip12 ip2 ip23
                    ip23 ip3 ip31
                    ip12 ip23 ip31];
                
            end
            
            m.nodes = newp;
            m.triangles = newt;
            
            % Create new edges
            newe = zeros(2*Ne,3);
            for n = 1 : Ne
                ip1 = e(n,1);
                ip2 = e(n,2);
                lab = e(n,3);
                ip12 = Nn+mat(ip1,ip2);
                newe(2*n-1:2*n,:) = [ip1 ip12 lab ; ip12 ip2 lab];
            end
            
            m.edges = newe;
            m.meshpropbuilder(dx);
            m.stats(8) = cputime-time;
        end
        function u = P0subtoP0(obj,usub)
            N = size(usub,2);
            Ntsub = obj.Ntriangles;
            Nt = Ntsub/4;
            ur = reshape(usub,4,N*Nt);
            s = obj.int1*ones(1,N);
            s = s(:);
            Ssub = reshape(s,4,N*Nt);
            u = sum(ur.*Ssub)./sum(Ssub);
            u = reshape(u,Nt,N);
        end
        function usub = P0toP0sub(obj,u)
            Nt = obj.Ntriangles;
            N = size(u,2);
            u = reshape(u,1,Nt*N);
            usub = ones(4,1)*u;
            usub = reshape(usub,4*Nt,N);
        end
        
        
        
        %%%%%%%% After this line it s only old functions from v2.2
        
        % Refiners
        function m = refine(obj,domain)
            disp('refining mesh')
            m = obj;
            
            dx = obj.stats(1);
            time=cputime;
            [gm,Ibound]=domain.matrix(dx);
            
            p = obj.nodes;
            e = obj.edges;
            t = obj.triangles;
            
            ne = size(e,1);
            nt = size(t,1);
            e2 = [e(:,1:2) zeros(ne,1) ones(ne,1) e(:,1) ones(ne,1) zeros(ne,1)];
            t2 = [t ones(nt,1)];
            
            [p,e,t]=refinemesh(gm,p',e2',t2');
            p=jigglemesh(p,e,t,'Iter',2);
            p=p';
            e=e';
            t=t';
            m.nodes=p;
            m.edges=[e(:,[1 2]) Ibound(e(:,5))];
            m.triangles=t(:,1:3);
            
            m.tricenters=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
            
            m.Iboundnodes = union(m.edges(:,1),m.edges(:,2));
            m.Iintnodes = setdiff(1:size(p,1), m.Iboundnodes)';
            
            
            q=mean(pdetriq(p',t'));
            a=2*sqrt(pdetrg(p',t'))/(5^(1/4));
            r=mean(a);
            rmin=min(a);
            rmax=max(a);
            [i1,ix,ixx] = m.integrals;
            
            m.Nnodes = size(p,1);
            m.Ntriangles = size(t,1);
            m.Nedges = size(e,1);
            
            m.int1 = i1;
            m.intx = ix;
            m.intxx = ixx;
            m.basis = m.buildbasis;
            w=whos('obj');
            mem=w.bytes/1000;
            
            time=cputime-time;
            m.stats=[dx rmin rmax r rmax/rmin q mem time];
            
            
        end
        
    
        % Interpolation
        function u = gridtoP0(obj,X,Y,u_grid)
            u = interp2(X,Y,u_grid,obj.tricenters(:,1),obj.tricenters(:,2));
        end
        function uP0 = P1toP0(obj,uP1)
            N=size(uP1,2);
            nT=obj.Ntriangles;
            uP0=zeros(nT,N);
            for i = 1 : N
                uP0(:,i)=pdeintrp(obj.nodes',obj.triangles',uP1(:,i));
            end
        end
        function uP1 = P0toP1(obj,uP0)
            N=size(uP0,2);
            nN=obj.Nnodes;
            uP1=zeros(nN,N);
            for i = 1 : N
                uP1(:,i)=pdeprtni(obj.nodes',obj.triangles',uP0(:,i).');
            end
        end
        function uP0h = P0toP0hex(obj,uP0)
           Nh = obj.Nhexa;
           dim = size(uP0,2);
           uP0h = reshape(uP0,6,Nh*dim);
           uP0h = sum(uP0h)/6;
           uP0h = reshape(uP0h,Nh,dim);
        end
        function uP0 = P0hextoP0(obj,uP0h)
            Nh = obj.Nhexa;
            dim = size(uP0h,2);
            uP0 =  reshape(uP0h,1,Nh*dim);
            uP0 = ones(6,1)*uP0;
            uP0 = reshape(uP0,6*Nh,dim);
        end
        
        function p = P0prod(obj,u,v)
            a = obj.int1;
            p = sum(u.*v.*a);
        end 
        function n = normL2(obj,u)
            % n = obj.normL2(u) returns le L^2-norm of the P1 function u.
           Nt = obj.Ntriangles;
                M = obj.mass(ones(Nt,1));
                n = sqrt(u'*(M*u));   
        end
         
      
        % Differential operators
        function du = div(obj,u)
            gu1 = obj.grad(u(:,1));
            gu2 = obj.grad(u(:,2));
            du = gu1(:,1)+gu2(:,2);
        end
        
        % Solvers
        function u = solvesystem(obj,A,F)
            Nn = obj.Nnodes;
           U = A\F;
           Nu = length(U);
           switch Nu
               case Nn
                   u = U;
               case 2*Nn
                   u = reshape(U,Nn,2);
               case Nn+1
                   u = U(1:Nn);
               case 2*Nn+3
                   u = reshape(U(1:2*Nn,:),Nn,2);
           end
               
        end

        % Scalar PDE matrices
        function K = stiffness(obj,a)
            Nt = obj.Ntriangles;
            if size(a,2) ~= 4 || size(a,1)~= Nt
                disp('MESH.intgugv error : input must be an P0 order 2 tensor !');
                K = [];
                return
            end
            
            tri = obj.triangles;
            
            tp=tri';
            
            rep3 = floor((3:9*Nt+2)/3);
            rep9 = floor((9:9*Nt+8)/9);
            rep33 = 3*(rep9-1)+1 + mod(0:9*Nt-1,3);
            
            Inodes = tp(:);
            I = Inodes(rep3);
            J = Inodes(rep33);
            
            alphamat = obj.basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            alphai = alpha(rep3,:);
            alphaj = alpha(rep33,:);
            
            a = a(rep9,:);
            
            i1 = obj.int1(rep9);
            
            agej = [a(:,1).*alphaj(:,1) + a(:,2).*alphaj(:,2)      a(:,3).*alphaj(:,1) + a(:,4).*alphaj(:,2)];
            geiagej = i1.*(alphai(:,1).*agej(:,1) + alphai(:,2).*agej(:,2));
            K = sparse(I,J,geiagej);
        end
        function A = intugv(obj,b)
            Nt = obj.Ntriangles;
            if size(b,2) ~= 2 || size(b,1)~= Nt
                disp('MESH.intgugv error : input must be an P0 vector function !');
                A = [];
                return
            end
            
            tri = obj.triangles;
            
            tp=tri';
            
            rep3 = floor((3:9*Nt+2)/3);
            rep9 = floor((9:9*Nt+8)/9);
            rep33 = 3*(rep9-1)+1 + mod(0:9*Nt-1,3);

            Inodes = tp(:);
            I = Inodes(rep3);
            J = Inodes(rep33);
            
            alphamat = obj.basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            alphai = alpha(rep3,:);
            alphaj = alpha(rep33,:);
            
            betamat = obj.basis(:,[3 6 9]);
            beta = reshape(betamat',1,3*Nt)';
            betaj = beta(rep33,:);
            
            b = b(rep9,:);
            
            i1 = obj.int1(rep9);
            ix = obj.intx(rep9,:);
            
            
            V1 = (alphai(:,1).*b(:,1) + alphai(:,2).*b(:,2)).*(alphaj(:,1).*ix(:,1) + alphaj(:,2).*ix(:,2));
            V2 = (alphai(:,1).*b(:,1) + alphai(:,2).*b(:,2)).*betaj.*i1;
           
            
        
            aeiej = V1+V2;
           
            A = sparse(I,J,aeiej);
        end
        function M = mass(obj,c)
            if nargin == 1
                M  = obj.mass(obj.P0(1));
                return
                
            end
            Nt = obj.Ntriangles;
            if size(c,2) ~= 1 || size(c,1)~= Nt
                disp('MESH.mass error : input must be an P0 scalar function !');
                M = [];
                return
            end
            
            tri = obj.triangles;
            
            tp=tri';
            
            rep3 = floor((3:9*Nt+2)/3);
            rep9 = floor((9:9*Nt+8)/9);
            rep33 = 3*(rep9-1)+1 + mod(0:9*Nt-1,3);

            Inodes = tp(:);
            I = Inodes(rep3);
            J = Inodes(rep33);
            
            alphamat = obj.basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            alphai = alpha(rep3,:);
            alphaj = alpha(rep33,:);

            betamat = obj.basis(:,[3 6 9]);
            beta = reshape(betamat',1,3*Nt)';
            betai = beta(rep3,:);
            betaj = beta(rep33,:);
            
            c = c(rep9,:);
            
            i1 = obj.int1(rep9);
            ix = obj.intx(rep9,:);
            ixx = obj.intxx(rep9,:);
            
            V1 = alphai(:,1).*(ixx(:,1).*alphaj(:,1) + ixx(:,2).*alphaj(:,2)) +... 
                 alphai(:,2).*(ixx(:,3).*alphaj(:,1) + ixx(:,4).*alphaj(:,2));
            V2 = (alphai(:,1).*ix(:,1) + alphai(:,2).*ix(:,2)).*betaj +...
                 (alphaj(:,1).*ix(:,1) + alphaj(:,2).*ix(:,2)).*betai;
            V3 = betai.*betaj.*i1;
            
        
            aeiej = c.*(V1 + V2 + V3);
           
            M = sparse(I,J,aeiej);
        end 
        function Mb = boundmass(obj,E,g)
             Ne = obj.Nedges;
             Nn = obj.Nnodes;
            if size(g,2) ~= 1 || size(g,1)~= Ne
                error('input 2 must be a P0bound scalar function !');
            end
            
            ebound = obj.edges;

            im = ismember(ebound(:,3),E);
            g = g(im,:);
            e = ebound(im,1:2);
            
            x = obj.nodes(:,1);
            y = obj.nodes(:,2);
            L = sqrt((x(e(:,2)) - x(e(:,1))).^2 + (y(e(:,2)) - y(e(:,1))).^2);
            e = e';
            I = e([1 1 2 2],:);
            J = e([1 2 1 2],:);
            
            V = 1/6*[2;1;1;2]*(L.*g)';
            Mb = sparse(I,J,V,Nn,Nn);
        end
        
        function F = source(obj,f)
            Nt = obj.Ntriangles;
            if size(f,2) ~= 1 || size(f,1)~= Nt
                disp('MESH.source error : input must be a P0 scalar function !');
                F = [];
                return
            end
            tri = obj.triangles;
            tp=tri';
            rep3 = floor((3:3*Nt+2)/3);
            
            I = tp(:);
            
            alphamat = obj.basis(:,[1 2 4 5 7 8]);
            alpha = reshape(alphamat',2,3*Nt)';
            
            betamat = obj.basis(:,[3 6 9]);
            beta = reshape(betamat',1,3*Nt)';

            f = f(rep3,:);
            
            i1 = obj.int1(rep3);
            ix = obj.intx(rep3,:);
            
            fei = f.*(alpha(:,1).*ix(:,1) + alpha(:,2).*ix(:,2) + beta.*i1);
            F = full(sparse(I,ones(3*Nt,1),fei));

        end  
        function Fb = boundsource(obj,E,f)
            Ne = obj.Nedges;
            Nn = obj.Nnodes;
            if size(f,2) ~= 1 || size(f,1)~= Ne
                disp('MESH.boundsource error : input must be a P0bound scalar function !');
                Fb = [];
                return
            end
            ebound = obj.edges;

            im = ismember(ebound(:,3),E);
            f = f(im,:);
            e = ebound(im,1:2);
            Ne = length(e);
            x = obj.nodes(:,1);
            y = obj.nodes(:,2);
            L = sqrt((x(e(:,2)) - x(e(:,1))).^2 + (y(e(:,2)) - y(e(:,1))).^2);
            ep = e';
            
            I = ep(:);
            J = ones(2*Ne,1);

            rep2 = floor((2:2*Ne+1)/2);

            f = f(rep2,:);
            
            fei = f.*L(rep2)/2;
            Fb = full(sparse(I,J,fei,Nn,1));
        end
        
        % Vectorial PDE matrices
        function K = stiffness_vect(obj,a)
            Nt = obj.Ntriangles;
            if size(a,2) ~= 16 || size(a,1)~= Nt
                disp('MESH.intgugv error : input must be an P0 order 4 tensor');
                K = [];
                return
            end
            
            A1 = obj.stiffness(a(:,1:4));
            A2 = obj.stiffness(a(:,5:8));
            A3 = obj.stiffness(a(:,9:12));
            A4 = obj.stiffness(a(:,13:16));
            K = [A1 A2 ; A3 A4];
            
        end
        function A = mass_vect(obj,a)
            Nt = obj.Ntriangles;
            if size(a,2) ~= 4 || size(a,1)~= Nt
                disp('MESH.intgugv error : input must be an P0 order 2 tensor !');
                A = [];
                return
            end
            
            A1 = obj.mass(a(:,1));
            A2 = obj.mass(a(:,2));
            A3 = obj.mass(a(:,3));
            A4 = obj.mass(a(:,4));
            A = [A1 A2 ; A3 A4];
            
        end       
        function F = source_vect(obj,f)
            Nt = obj.Ntriangles;
            if size(f,2) ~= 2 || size(f,1)~= Nt
                disp('MESH.intv_vect error : input must be a P0 vector function !');
                F = [];
                return
            end
            F1 = obj.source(f(:,1));
            F2 = obj.source(f(:,2));
            F = [F1;F2];

        end 
        function F = boundsource_vect(obj,E,f)
            Ne = obj.Nedges;
            if size(f,2) ~= 2 || size(f,1)~= Ne
                disp('MESH.intboundv_vect error : input must be a P0 vector function !');
                F = [];
                return
            end
            F1 = obj.boundsource(E,f(:,1));
            F2 = obj.boundsource(E,f(:,2));
            F = [F1;F2];
         end
        
        % Boundary Conditions applyers
        function [Mb,Fb,Id,Vd] = bcmatrices(obj,bclist)
            Nbclist = length(bclist);
            Nn = obj.Nnodes;
            e = obj.edges(:,3);
            if Nbclist == 0
                Id = -1;
                Vd = [];
                Mb = sparse(Nn,Nn);
                Fb = zeros(Nn,1);
                return
            end
            onlyneumann = 1;
            Id = [];
            Vd = zeros(Nn,1);
            Mb = sparse(Nn,Nn);
            Fb = zeros(Nn,1);
            for k = 1 : Nbclist
                bc = bclist{k};
                if bc{2} == 0
                    onlyneumann = 0;
                    ed=bc{1};
                    id=[];
                    for i = 1 : length(ed)
                        id=[id ; obj.edges(e==ed(i),1) ; obj.edges(e==ed(i),2)];
                    end
                    id=unique(id);
                    Id=union(Id,id);
                    Vdk=obj.P1(bc{4},0);
                    Vd(id)=Vdk(id);
                end
                if bc{2} == 1
                    g = obj.P0bound(bc{3});
                    if sum(abs(g))~=0
                        onlyneumann = 0;
                        Mb = Mb + obj.boundmass(bc{1},g);
                    end
                    h = obj.P0bound(bc{4},0);
                    Fb = Fb + obj.boundsource(bc{1},h);
                    
                        
                end    
                
            end
            if onlyneumann
                    Id = -1;
            end
        end   
         function [Mb,Fb,Id,Vd] = bcmatrices_vect(obj,bclist)
            Nbclist = length(bclist);
            Nn = obj.Nnodes;
            e = obj.edges(:,3);
            if Nbclist == 0
                Id = -1;
                Vd = [];
                Mb = sparse(2*Nn,2*Nn);
                Fb = zeros(2*Nn,1);
                return
            end
            onlyneumann = 1;
            Id = [];
            Vd = zeros(2*Nn,1);
            Mb = sparse(2*Nn,2*Nn);
            Fb = zeros(2*Nn,1);
            for k = 1 : Nbclist
                bc = bclist{k};
                if bc{2} == 0
                    onlyneumann = 0;
                    ed=bc{1};
                    id=[];
                    for i = 1 : length(ed)
                        id=[id ; obj.edges(e==ed(i),1) ; obj.edges(e==ed(i),2)];
                    end
                    id=unique(id);
                    Id=union(Id,id);
                    h = obj.P1(bc{4},1);
                    Vd(id)=h(id,1);
                    Vd(Nn+id)=h(id,2);
                end
                
                if bc{2} == 1
                    g = obj.P0bound(bc{3},2);
                    if sum(abs(g))~=0
                        onlyneumann = 0;
                        Mb = Mb + obj.boundmass_vect(bc{1},g);
                    end
                    h = obj.P0bound(bc{4},1);
                    Fb = Fb + obj.boundsource_vect(bc{1},h);
                    
                        
                end    
                
            end
            Id = [Id ; Nn+Id];
            
            if onlyneumann
                    Id = -1;
            end
        end   
        
        % Scalar pde builders and solvers
        function [A,F,Id,Vd] = pde(obj,a,c,f,bclist)
            Nn = obj.Nnodes;
            K = sparse(Nn,Nn);
            M = sparse(Nn,Nn);
            S = zeros(Nn,1);
            
            
            if nnz(a)
                K = obj.stiffness(obj.P0(a,2));
            end
            if nnz(c)
                M = obj.mass(obj.P0(c));
            end
             if nnz(f)
                S = obj.source(obj.P0(f));
             end
             
             [Mb,Sb,Id,Vd] = obj.bcmatrices(bclist);
             A = K+M+Mb;
             F = S + Sb;
        end     
        function [A,F] = system(obj,a,c,f,bclist)
            [A,F,Id,Vd] = obj.pde(a,c,f,bclist);
            [A,F] = obj.assemble(A,F,Id,Vd);
        end

        
        % Vectorial pde builders and solvers
        function [A,F,Id,Vd] = pde_vect(obj,a,c,f,bclist)
            a = obj.P0(a,4);
            c = obj.P0(c,2);
            f = obj.P0(f,1);
            Nn = obj.Nnodes;
            K = sparse(2*Nn,2*Nn);
            M = sparse(2*Nn,2*Nn);
            S = zeros(2*Nn,1);
            
            if norm(a)
                K = obj.stiffness_vect(a);
            end
            if norm(c)
                M = obj.mass_vect(c);
            end
             if norm(f)
                S = obj.source_vect(f);
             end
             
             [Mb,Sb,Id,Vd] = obj.bcmatrices_vect(bclist);
             A = K+M+Mb;
             F = S + Sb;
        end
        function [A,F] = system_vect(obj,a,c,f,bclist)
            [A,F,Id,Vd] = obj.pde_vect(a,c,f,bclist);
            [A,F] = obj.assemble(A,F,Id,Vd);
        end
        function u = solve_vect(obj,a,c,f,bclist)
            [A,F] = obj.system_vect(a,c,f,bclist);
            Nn = obj.Nnodes;
            U = A\F;
            u = [U(1:Nn) U(Nn+1:end)];
        end

        % Tools
        function E = buildbasis(obj)
            Node1=obj.nodes(obj.triangles(:,1),:);
            Node2=obj.nodes(obj.triangles(:,2),:);
            Node3=obj.nodes(obj.triangles(:,3),:);
            
            E = zeros(obj.Ntriangles,9);
            for j = 1 : obj.Ntriangles
                M=[[Node1(j,:) ; Node2(j,:) ; Node3(j,:)] [1;1;1]];
                N=inv(M);
                E(j,:)=N(:)';
            end
        end
        function [int1,intx,intxx] = integrals(obj)
            p=obj.nodes;
            t=obj.triangles;
            
            p1=p(t(:,1),:)';
            p2=p(t(:,2),:)';
            p3=p(t(:,3),:)';
            M=[p2-p1;p3-p1];
            det=abs(M(1,:).*M(4,:)-M(2,:).*M(3,:));
            
            int1 = det'/2;
            
            intx=[det.*(1/2*p1(1,:)+1/6*(M(1,:)+M(3,:))) ; det.*(1/2*p1(2,:)+1/6*(M(2,:)+M(4,:)))]';
            
            I1=1/2*[p1(1,:).^2;p1(1,:).*p1(2,:);p1(1,:).*p1(2,:);p1(2,:).^2];
            I2=1/6*[p1(1,:).*(M(1,:)+M(3,:));p1(1,:).*(M(2,:)+M(4,:));p1(2,:).*(M(1,:)+M(3,:));p1(2,:).*(M(2,:)+M(4,:))];
            I3=[I2(1,:);I2(3,:);I2(2,:);I2(4,:)];
            
            L1=1/24*(2*M(1,:).^2+2*M(1,:).*M(3,:)+2*M(3,:).^2);
            L2=1/24*(2*M(2,:).*M(1,:)+M(4,:).*M(1,:)+M(3,:).*M(2,:)+2*M(4,:).*M(3,:));
            L4=1/24*(2*M(2,:).^2+2*M(2,:).*M(4,:)+2*M(4,:).^2);
            
            I4=[L1;L2;L2;L4];
            I=I1+I2+I3+I4;
            intxx=[det.*I(1,:);det.*I(2,:);det.*I(3,:);det.*I(4,:)]';
            
        end
        function [A,F] = assemble(obj,A,F,Id,Vd)
            Nn = obj.Nnodes;
            if Id == -1
                if length(F) == Nn
                    if abs(sum(F))> 1e-12
                        disp('Mesh.assemble warning : sum(F) must be zero for full Neumann boundary conditions !')
                        F = F - sum(F)/Nn;
                    end
                    A(1,:)=0;
                    A(1,1)=1;
                    F(1)=0;
                elseif length(F) == 2*Nn
                    F1 = F(1:Nn);
                    F2 = F(Nn+1:2*Nn);
                    if abs(sum(F))> 2e-12
                        F1 = F1 - sum(F1)/Nn;
                        F2 = F2 - sum(F2)/Nn;
                    end
                    A([1 2 Nn+1],:)=0;
                    A(1,1) = 1;
                    A(2,2) = 1;
                    A(Nn+1,Nn+1) = 1;
                    F1([1 2 Nn+1]) = 0;
                    
                end
            elseif ~isempty(Id)
                NId = length(Id);
                A(Id,:) = 0;
                A(Id,Id) = speye(NId,NId);
                F(Id) = Vd(Id);
            end
        end 
        function surfrs(obj,u)
            N = size(u,1);
            t = obj.triangles;
            t = [t ones(length(t),1)];
            
            if obj.Ntriangles == N
                pdesurf(obj.nodes',t',u');
                view([0 90])
                colormap jet
                colorbar
                axis image
            elseif  obj.Nnodes == N
                pdesurf(obj.nodes',t',u);
            else
                error('Wrong data format')
            end
            view([0 90])
            colormap jet
            colorbar
            axis image
        end
        function contourrs(obj,u)
            N = size(u,1);
            t = obj.triangles;
            t = [t ones(length(t),1)];
            
            if obj.Ntriangles == N
                pdesurf(obj.nodes',t',u');
                view([0 90])
                colormap jet
                colorbar
                axis image
            elseif  obj.Nnodes == N
                pdecont(obj.nodes',t',u);
            else
                error('Wrong data format')
            end
            view([0 90])
            colormap jet
            colorbar
            axis image
        end
        function valout = tensor(obj,valin,order)
            N = size(valin,2);
            n = size(valin,1);
            switch order
                case -1
                    valout = valin;
                case 0
                    switch N
                        case 1
                            valout = valin;
                        otherwise
                            error('input and order 0 are not compatible')
                    end
                case 1
                    switch N
                        case 1
                            valout = valin*[1 1];
                        case 2
                            valout = valin;
                        otherwise
                            error('input and order 1 are not compatible')
                    end
                    
                case 2
                    switch N
                        case 1
                            valout = valin*[1 0 0 1];
                        case 2
                            valout = [valin(:,1) zeros(n,2) valin(:,2)];
                        case 3
                            valout = [valin(:,1) valin(:,3) valin(:,3) valin(:,2)];
                        case 4
                            valout = valin;
                        otherwise
                            error('input and order 2 are not compatible')
                    end
                    
                case 4
                    switch N
                        case 1
                            I = [1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1];
                            valout = valin(:,1)*I;
                        case 2
                             I = [1 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1];
                            T = [1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1];
                            D = [1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1];
                            valout = (valin(:,1))*(I+T) + valin(:,2)*D;
                        case 4
                            s11 = valin(:,1);
                            s12 = valin(:,2);
                            s21 = valin(:,3);
                            s22 = valin(:,4);
                            valout = [s11.*s11 s11.*s21  s21.*s11 s21.*s21 ...
                                s11.*s12 s11.*s22  s21.*s12 s21.*s22 ...
                                s12.*s11 s12.*s21  s22.*s11 s22.*s21 ...
                                s12.*s12 s12.*s22  s22.*s12 s22.*s22 ...
                                ];
                        case 8
                            s11 = valin(:,1);
                            s12 = valin(:,2);
                            s21 = valin(:,3);
                            s22 = valin(:,4);
                            t11 = valin(:,5);
                            t12 = valin(:,6);
                            t21 = valin(:,7);
                            t22 = valin(:,8);
                            
                            valout = [s11.*t11 s11.*t21  s21.*t11 s21.*t21 ...
                                s11.*t12 s11.*t22  s21.*t12 s21.*t22 ...
                                s12.*t11 s12.*t21  s22.*t11 s22.*t21 ...
                                s12.*t12 s12.*t22  s22.*t12 s22.*t22 ...
                                ];
                            
                        case 16
                            valout = valin;
                        otherwise
                            error('input and order 2 are not compatible')
                    end
                    
                otherwise
                    disp(order)
                    error('Invalid order')
                    
            end
        end
        function I = boundary(obj,L)
            if nargin == 1
                e1 = obj.edges(:,1);
                e2 = obj.edges(:,2);
                I = union(e1,e2);
            else
                I = [];
                for n = 1 : length(L)
                    ind = obj.edges(:,3) == L(n);
                    e1 = obj.edges(ind,1);
                    e2 = obj.edges(ind,2);
                    I = union(I,union(e1,e2));
                end               
            end
        end
            
        
        % Hexagonal meshes
        function [p,e,t,hex] = hexagonalmesh(obj,box,dx)
            xmin = box(1);
            xmax = box(2);
            ymin = box(3);
            ymax = box(4);
            Nx = floor(((xmax-xmin)/(3*dx)));
            Ny = floor(((ymax-ymin)/(sqrt(3)/2*dx)))-1;
            hmesh = [];
            
            for j = 1 : Ny
                yj = ymin+j*dx*sqrt(3)/2;
                if mod(j,2)
                    for i = 1 : Nx
                        xi = xmin +dx+3*(i-1)*dx;
                        hmesh = addhexa(obj,hmesh,[xi yj],dx);
                    end
                else
                    for i = 1 : Nx
                        xi = xmin +2.5*dx+3*(i-1)*dx;
                        hmesh = addhexa(obj,hmesh,[xi yj],dx);
                    end
                end
                
            end
            
            nodes = hmesh.nodes;
            [a,I,J] = uniquetol(nodes,dx/10,'ByRows',true);
            newnodes = nodes(I,:);
            
            triangles = hmesh.triangles;
            Ntriangles = hmesh.Ntriangles;
            t = triangles(:);
            t = J(t);
            newtriangles = reshape(t,Ntriangles,3);
            
            hmesh.triangles = newtriangles;
            hmesh.nodes = newnodes;
            
            cedges = [];
            for i = 1 : hmesh.Nhex
                h = hmesh.hex(i,:);
                for k = 1 : 6
                    e = hmesh.triangles(h(k),1:2);
                    e = sort(e);
                    ce = e(1)+1i*e(2);
                    ind = find(abs(cedges - ce)<dx/10,1,'first');
                    if isempty(ind)
                        cedges = [cedges ; ce];
                    else
                        cedges = [cedges(1:ind-1) ; cedges(ind+1:end)];
                    end
                end
            end
            Nedges = length(cedges);
            edges = [floor(real(cedges)) floor(imag(cedges))];
            hmesh.edges = [edges ones(Nedges,1)];
            
            p = hmesh.nodes;
            e = hmesh.edges;
            t = hmesh.triangles;
            hex = hmesh.hex;
            
        end
        function hmesh = addhexa(obj,hmesh,c,dx)
            if isempty(hmesh)
                p = zeros(7,2);
                p(7,1:2) = c;
                for i = 1:6
                    p(i,1:2) = c + dx*[cos((i-1)*pi/3) sin((i-1)*pi/3)];
                end
                t = [1 2 7
                    2 3 7
                    3 4 7
                    4 5 7
                    5 6 7
                    6 1 7];
                
                hmesh.hex =  [1 2 3 4 5 6];
                hmesh.triangles = t;
                hmesh.nodes = p;
                
                hmesh.Nhex =  1;
                hmesh.Ntriangles = 6;
                hmesh.Nnodes = 7;
                return
            end
            
            Nn = size(hmesh.nodes,1);
            Nt = size(hmesh.triangles,1);
            hm = addhexa(obj,[],c,dx);
            hmesh.hex = [ hmesh.hex ; Nt + hm.hex];
            hmesh.triangles = [ hmesh.triangles ; Nn + hm.triangles];
            hmesh.nodes = [ hmesh.nodes ; hm.nodes];
            
            hmesh.Nhex = hmesh.Nhex + 1;
            hmesh.Ntriangles = hmesh.Ntriangles + 6;
            hmesh.Nnodes = hmesh.Nnodes + 7;
        end
        
    end
    
    methods (Access = private)
        
        % Private tools
        function meshpropbuilder(obj,dx)
            p = obj.nodes;
            e = obj.edges;
            t = obj.triangles;
            obj.tricenters=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
            
            obj.Iboundnodes = union(obj.edges(:,1),obj.edges(:,2));
            obj.Iintnodes = setdiff(1:size(p,1), obj.Iboundnodes)';
            
            
            q=mean(pdetriq(p',t'));
            a=2*sqrt(pdetrg(p',t'))/(5^(1/4));
            r=mean(a);
            rmin=min(a);
            if dx == -1
                dx = rmin;
            end
            rmax=max(a);
            [i1,ix,ixx] = obj.integrals;
            
            obj.Nnodes = size(p,1);
            obj.Ntriangles = size(t,1);
            obj.Nedges = size(e,1);
            
            obj.int1 = i1;
            obj.intx = ix;
            obj.intxx = ixx;
            obj.basis = obj.buildbasis;
            
            props = properties(obj);
            
            total_mem = 0;
            for ii=1:length(props)
                curr_prop = obj.(props{ii});
                s = whos('curr_prop');
                total_mem = total_mem + s.bytes;
            end
            mem = total_mem/1024;
            
            obj.stats=[dx rmin rmax r rmax/rmin q mem 0];
        end
    end
end

% methods to be repaired
function m=refine(obj,g,f,level)
disp('Refine is obsolete')
% MESH.REFINE
%       m=obj.refine(g) refines the Mesh object. g is the
%       initial Geometry object.
%
%       m=obj.refine(g,f) refines only the area specified by
%       the boolean function f. Matrix f is 1xnT, 1xnP or expr.
%
%       m=obj.refine(g,f,level) refines around the highest
%       value of the non negative map f. Matrix f is 1xnT, 1xnP
%       or expr. level is a double in [0,1].
time=cputime;
nt=obj.stats(3);
m=obj;
GM=g.matrix;
p=m.nodes;
e=m.edges;
t=m.triangles;

switch nargin
    case 2
        F=ones(1,nt);
        s=1/2;
    case 3
        F=obj.tensor(f,0,1);
        s=1/2;
    case 4
        F=obj.tensor(f,0,1);
        s=level;
end
maxf=max(F);
it=find(F>((1-s)*maxf));
[p,e,t]=refinemesh(GM,p',e',t');
p=jigglemesh(p,e,t,'Iter',2);

m.nodes=p;
m.edges=e;
m.triangles=t;
w=whos('obj');
mem=w.bytes/1000;
q=mean(pdetriq(p,t));
a=pdetrg(p,t);
r=mean(a);
rmin=min(a);
rmax=max(a);
time=cputime-time;
obj.Nnodes = size(p,2);
obj.Ntriangles = size(t,2);
m.stats=[dx rmin rmax r rmax/rmin q mem time+obj.stats(10)];
end






