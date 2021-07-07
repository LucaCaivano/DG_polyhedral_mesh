%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Modification of PolyMesher      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [region] = generate_mesh_plate(Dati,dom,Domain,NElem,MaxIter,P)
if ~exist('P','var'), P=PolyMshr_RndPtSet(NElem,Domain,dom); end

BdBox = Domain(dom,'BdBox');

NElem = size(P,1);
Tol=5e-3; It=0; Err=1; c=1.5;%c=1.5;
BdBox = Domain(dom,'BdBox');
Area = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
Pc = P;
while(It<=MaxIter && Err>Tol)
  Alpha = c*sqrt(Area/NElem);
  P = Pc;                                       %Lloyd's update
  R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha,dom);   %Generate the reflections 
  [Node,Element] = voronoin([P;R_P],{'Qbb','Qz'});           %Construct Voronoi diagram
  [Pc,A] = PolyMshr_CntrdPly(Element,Node,NElem);
  Area = sum(abs(A));
  Err = sqrt(sum((A.^2).*sum((Pc-P).*(Pc-P),2)))*NElem/Area^1.5;
%   fprintf('It: %3d   Error: %1.3e\n',It,Err); It=It+1;
end
[Node,Element] = PolyMshr_ExtrNds(NElem,Node,Element);  %Extract node list
[Node,Element] = PolyMshr_CllpsEdgs(Node,Element,0.1);  %Remove small edges
[Node,Element] = PolyMshr_RsqsNds(Node,Element);        %Reoder Nodes
BC=Domain(dom,'BC',Node); Supp=BC{1}; Load=BC{2};           %Recover BC arrays
% poly_mesh_plot(Node,Element,NElem,'k');


ne = length(Element);
elem_area = zeros(ne,1);
nedge = zeros(ne,1);
BBox = zeros(ne,4);

coords_element = cell(1,ne);
max_kb = cell(1,ne);

coord = Node;



% distortion_factor = 0.15;
% x0 = BdBox(1); x1 = BdBox(2);
% y0 = BdBox(3); y1 = BdBox(4);
% for ip = 1 : length(coord)
%    x = coord(ip,1);
%    y = coord(ip,2);
% 
%    %dx = (-1).^ip.*(0.25.*distortion_factor).*rand(1);
%    %dy = (-1).^ip.*(0.25.*distortion_factor).*rand(1);
%  
%    %if (abs(x-x0) > 1.e-1 && abs(x-x1)>1.e-1 &&  abs(y-y0)>1.e-1 && abs(y-y1)>1.e-1); x = x + dx; y = y + dy; end
% 
%    % boundary points
%   % if ((y  == y0 || y  == y1) && (x  ~= x0 && x ~= x1)); x = x + dx; end
% 
%    % boundary point
%   % if ((x  == x0 || x  == x1) && (y  ~= y0 && y ~= y1)); y = y + dy; end
% 
%   %if (abs(x-x0) > 1.e-1 && abs(x-x1)>1.e-1 &&  abs(y-y0)>1.e-1 && abs(y-y1)>1.e-1)
%   %    x = x + sin(10*pi*x)/50;
%   %    y = y + sin(5*pi*y)/90;
%      
%   %end
%    coord(ip,1) = x;
%    coord(ip,2) = y;
% 
% end

%%


for i = 1:length(Element)
%     disp(i)
    nedge(i) = length(Element{i});
    coords_element{i} = coord(Element{i},:);
    x_min = min( coords_element{i}(:,1));x_max = max( coords_element{i}(:,1));
    y_min = min( coords_element{i}(:,2));y_max = max( coords_element{i}(:,2));
    BBox(i,:)=[x_min x_max y_min y_max];

    elem_area(i) = polyarea(coords_element{i}(:,1),coords_element{i}(:,2));
    max_kb{i} = zeros(nedge(i),1);
    for j = 1:nedge(i)
        if j<nedge(i)
            v1 = coords_element{i}(j,:); v2 = coords_element{i}(j+1,:);
            ch_j = j;
            ch_j_1 = j +1;
        else
            v1 = coords_element{i}(j,:); v2 = coords_element{i}(1,:);
            ch_j = j;
            ch_j_1 = 1;
        end
        for k = 1:nedge(i)
            
            
%             if k~=j && k~=(j+1)
            if k~= ch_j && k~= ch_j_1

%                 disp([j,k,nedge(i)])

                
                v3 = coords_element{i}(k,:);
                [x_tria,y_tria]=poly2cw([v1(1) v2(1) v3(1)],[v1(2) v2(2) v3(2)]);
                [x1,y1] = polybool('intersection',coords_element{i}(end:-1:1,1),coords_element{i}(end:-1:1,2),x_tria,y_tria);
                
                area = polyarea(x_tria,y_tria);
                
                if 1-any(isnan(x1)) && abs(polyarea(x1,y1)- area) < 1e-6
                    
%                     disp([area,max_kb{i}(j)]);
                    if area>max_kb{i}(j)
                        max_kb{i}(j) = area;
                    end
                end
                
            end
        end

    end
    
end

region=struct('nedges', nedge',...
    'BBox',BBox,...
    'ne',length(Element),...
    'coord',coord);
region.coords_element = coords_element;
region.connectivity = Element;
region.area = elem_area;
region.max_kb = max_kb;
region.degree = Dati.fem;


%%% POLYMESHER FUNCTIONS %%%%

%------------------------------------------------- GENERATE RANDOM POINTSET
function P = PolyMshr_RndPtSet(NElem,Domain,dom)
P=zeros(NElem,2); BdBox=Domain(dom,'BdBox'); Ctr=0;
while Ctr<NElem  
  Y(:,1) = (BdBox(2)-BdBox(1))*rand(NElem,1)+BdBox(1);
  Y(:,2) = (BdBox(4)-BdBox(3))*rand(NElem,1)+BdBox(3);
  d = Domain(dom,'Dist',Y);
  I = find(d(:,end)<0);                 %Index of seeds inside the domain
  NumAdded = min(NElem-Ctr,length(I));  %Number of seeds that can be added
  P(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
  Ctr = Ctr+NumAdded;
end
%--------------------------------------------------------- REFLECT POINTSET
function R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha,dom)
eps=1e-8; eta=0.9;
d = Domain(dom,'Dist',P);  
NBdrySegs = size(d,2)-1;          %Number of constituent bdry segments
n1 = (Domain(dom,'Dist',P+repmat([eps,0],NElem,1))-d)/eps;
n2 = (Domain(dom,'Dist',P+repmat([0,eps],NElem,1))-d)/eps;
I = abs(d(:,1:NBdrySegs))<Alpha;  %Logical index of seeds near the bdry
P1 = repmat(P(:,1),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,1)
P2 = repmat(P(:,2),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,2)
R_P(:,1) = P1(I)-2*n1(I).*d(I);  
R_P(:,2) = P2(I)-2*n2(I).*d(I);
d_R_P = Domain(dom,'Dist',R_P);
J = abs(d_R_P(:,end))>=eta*abs(d(I)) & d_R_P(:,end)>0;
R_P=R_P(J,:); R_P=unique(R_P,'rows');
%---------------------------------------------- COMPUTE CENTROID OF POLYGON
function [Pc,A] = PolyMshr_CntrdPly(Element,Node,NElem)
Pc=zeros(NElem,2); A=zeros(NElem,1);
for el = 1:NElem
  vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(Element{el}); 
  vxS=vx([2:nv 1]); vyS=vy([2:nv 1]); %Shifted vertices
  temp = vx.*vyS - vy.*vxS;
  A(el) = 0.5*sum(temp);
  Pc(el,:) = 1/(6*A(el,1))*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
end
%------------------------------------------------------- EXTRACT MESH NODES
function [Node,Element] = PolyMshr_ExtrNds(NElem,Node0,Element0)
map = unique([Element0{1:NElem}]);
cNode = 1:size(Node0,1);
cNode(setdiff(cNode,map)) = max(map);
[Node,Element] = PolyMshr_RbldLists(Node0,Element0(1:NElem),cNode);
%----------------------------------------------------- COLLAPSE SMALL EDGES
function [Node0,Element0] = PolyMshr_CllpsEdgs(Node0,Element0,Tol)
while(true)
  cEdge = [];
  for el=1:size(Element0,1)
    if size(Element0{el},2)<4, continue; end;  %Cannot collapse triangles
    vx=Node0(Element0{el},1); vy=Node0(Element0{el},2); nv=length(vx);
    beta = atan2(vy-sum(vy)/nv, vx-sum(vx)/nv);
    beta = mod(beta([2:end 1]) -beta,2*pi);
    betaIdeal = 2*pi/size(Element0{el},2);
    Edge = [Element0{el}',Element0{el}([2:end 1])'];
    cEdge = [cEdge; Edge(beta<Tol*betaIdeal,:)];
  end
  if (size(cEdge,1)==0), break; end
  cEdge = unique(sort(cEdge,2),'rows');
  cNode = 1:size(Node0,1);
  for i=1:size(cEdge,1)
    cNode(cEdge(i,2)) = cNode(cEdge(i,1));
  end
  [Node0,Element0] = PolyMshr_RbldLists(Node0,Element0,cNode);
end
%--------------------------------------------------------- RESEQUENSE NODES
function [Node,Element] = PolyMshr_RsqsNds(Node0,Element0)
NNode0=size(Node0,1); NElem0=size(Element0,1);
ElemLnght=cellfun(@length,Element0); nn=sum(ElemLnght.^2); 
i=zeros(nn,1); j=zeros(nn,1); s=zeros(nn,1); index=0;
for el = 1:NElem0
  eNode=Element0{el}; ElemSet=index+1:index+ElemLnght(el)^2;
  i(ElemSet) = kron(eNode,ones(ElemLnght(el),1))';
  j(ElemSet) = kron(eNode,ones(1,ElemLnght(el)))';
  s(ElemSet) = 1;
  index = index+ElemLnght(el)^2;
end
K = sparse(i,j,s,NNode0, NNode0);
p = symrcm(K);
cNode(p(1:NNode0))=1:NNode0;
[Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode);
%------------------------------------------------------------ REBUILD LISTS
function [Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode)
Element = cell(size(Element0,1),1);
[foo,ix,jx] = unique(cNode);
if size(Node0,1)>length(ix), ix(end)=max(cNode); end;
Node = Node0(ix,:); 
for el=1:size(Element0,1)
  Element{el} = unique(jx(Element0{el}));
  vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(vx);
  [foo,iix] = sort(atan2(vy-sum(vy)/nv,vx-sum(vx)/nv));
  Element{el} = Element{el}(iix);
end