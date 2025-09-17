function [DEFL,REACT,ELE_FOR,AFLAG] = ud_3d1el(...
		nnodes,coord,concen,fixity,nele,ends,A,Izz,Iyy,J,Cw,IsSym,Ysc,Zsc,Betay,Betaz,Betaw,Zzz,Zyy,Ayy,Azz,...
		E,v,Fy,YldSurf,Wt,webdir,beta_ang,w,thermal,truss,anatype);

	DEFL=[]; REACT=[]; ELE_FOR=[];

clc;

% creating matrix of nodal DOFs
nodeDOFs = zeros(length(fixity(:,1)),6);
for i = 1:length(fixity(:,1))
    for j = 1:6
        nodeDOFs(i,j) = 6*i-(6-j);
    end
end

% Initialising global structure stiffness matrix
K_global = zeros(size(nodeDOFs,1)*6,size(nodeDOFs,1)*6);

% creating matrix of element internal DOFs
memb_id = zeros(length(ends(:,1)),12);

for i = 1:size(memb_id,1)
    for j = 1:2
        memb_id(i,[6*j-5:6*j]) = nodeDOFs(ends(i,j),:);
    end
end

% creating nodal DOF load vector 
dof_loads = zeros(size(concen,1)*6,1);
FeF_g = zeros(size(concen,1)*6,1);

for i = 1:size(concen,1)
    dof_loads(6*i-5:6*i,1) = concen(i,:);
end

% identifying the free, known and supported DOFs 
ctr1=1;ctr2=1;ctr5=1;
delta_n = [];id_delta_n = [];
for i = 1:size(fixity,1)
    for j = 1:6
        if isnan(fixity(i,j)) == 1
            id_delta_f(ctr1,1) = nodeDOFs(i,j);
            delta_f(ctr1,1) = NaN;
            ctr1 = ctr1+1;
        elseif fixity(i,j) == 0 
            id_delta_s(ctr2,1) = nodeDOFs(i,j);
            delta_s(ctr2,1) = 0;
            ctr2 = ctr2+1;
        else
            id_delta_n(ctr5,1) = nodeDOFs(i,j);
            delta_n(ctr5,1) = fixity(i,j);
            ctr5 = ctr5+1;
        end
    end
end

delta = [delta_f;delta_n;delta_s]; % creating the initial nodal displacement vector

% Assembling of member stiffness matrices into global structure stiffness
% matrix
for e = 1:size(ends,1)

    n1 = ends(e,1); % start node 
    n2 = ends(e,2); % end node
    
    % start node coordinates
    x1 = coord(n1,1);
    y1 = coord(n1,2);
    z1 = coord(n1,3);
    
    % end node coordinates
    x2 = coord(n2,1);
    y2 = coord(n2,2);
    z2 = coord(n2,3);

    L = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2); % Length

    % if ends(e,3) == 1 % condition for flexural release at the ith node 
    % call the stiffness function with flag
    elk = BDLT_estiff(A(e),Izz(e),Iyy(e),J(e),Ayy(e),Azz(e),E(e),v(e),L,ends(e,[3:4]));
  
    % call the transformation matrix with flag
    gamma = BDLT_etran([x1 y1 z1],[x2 y2 z2],webdir(e,:));

    % call the memberFeFs function to get the FeF vector for this member
    FeF = memberFeFs(w(e,:),L,ends(e,[3:4]));

    index = memb_id(e,:); % extracting the member internal DOFs 
    
    % Assembling the element stiffnesses into the global stiffness matrix
    K_global(index,index) = K_global(index,index) + gamma'*elk*gamma;

    % Assembling element FeFs into the global FeF vector
    FeF_g(index,1) = FeF_g(index,1) + gamma'*FeF;

end

% Solving for the free DOF displacements
if isempty(delta_n) == 1
    delta_f = inv(K_global(id_delta_f,id_delta_f))*(dof_loads(id_delta_f,1)...
              -FeF_g(id_delta_f,1)-K_global(id_delta_f,id_delta_s)*delta_s);
else
    delta_f = inv(K_global(id_delta_f,id_delta_f))*(dof_loads(id_delta_f,1)...
              -FeF_g(id_delta_f,1)-K_global(id_delta_f,id_delta_n)*delta_n-K_global(id_delta_f,id_delta_s)*delta_s);
end
delta(1:length(id_delta_f),1) = delta_f;

% Rearranging nodal displacement vector for MASTAN2 readability
DEFL = fixity;ctr3 = 1;
for i=1:size(DEFL,1)
    for j = 1:6
        if any(id_delta_f(:)==nodeDOFs(i,j)) == 1
            DEFL(i,j) = delta(ctr3);
            ctr3 = ctr3+1;
        end
    end
end

% Calculating nodal reactions
REACT = zeros(length(fixity(:,1)),6);ctr4=1;ctr6=1;
if isempty(delta_n) == 1
    R_support = K_global(id_delta_s,id_delta_f)*delta_f+K_global(id_delta_s,id_delta_s)*delta_s+FeF_g(id_delta_s);
else
    R_support = K_global(id_delta_s,id_delta_f)*delta_f+K_global(id_delta_s,id_delta_n)*delta_n+FeF_g(id_delta_s);
    R_known = K_global(id_delta_n,id_delta_f)*delta_f+K_global(id_delta_n,id_delta_n)*delta_n+FeF_g(id_delta_n);
end
for i=1:size(REACT,1) % Rearranging for MASTAN2 readability
    for j = 1:6
        if any(id_delta_s(:)==nodeDOFs(i,j)) == 1
            REACT(i,j) = R_support(ctr4);
            ctr4 = ctr4+1;
        elseif any(id_delta_n(:) == nodeDOFs(i,j)) == 1
            REACT(i,j) = R_known(ctr6);
            ctr6 = ctr6+1;
        end
    end
end

% Computing member end forces in local coordinate system
for i = 1:size(DEFL,1)
    DEFL_vector(6*i-5:6*i,1) = DEFL(i,:);
end
for e = 1:size(ends,1)

    n1 = ends(e,1); 
    n2 = ends(e,2);

    x1 = coord(n1,1);
    y1 = coord(n1,2);
    z1 = coord(n1,3);

    x2 = coord(n2,1);
    y2 = coord(n2,2);
    z2 = coord(n2,3);

    L = sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);

    % if ends(e,3) == 1 % condition for flexural release at the ith node 
    % call the stiffness function with flag
    elk = BDLT_estiff(A(e),Izz(e),Iyy(e),J(e),Ayy(e),Azz(e),E(e),v(e),L,ends(e,[3:4]));
  
    % call the transformation matrix with flag
    gamma = BDLT_etran([x1 y1 z1],[x2 y2 z2],webdir(e,:));

    % call the memberFeFs function to get the FeF vector for this member
    FeF = memberFeFs(w(e,:),L,ends(e,[3:4]));

    index = memb_id(e,:); 

    del_elem = gamma*DEFL_vector(index,1); % converting the member global displacements into local coordinates
   
    ELE_FOR(e,:) = elk*del_elem + FeF; % Calculating member end forces in local coordinates
end

    % Condition for unstable system
    if isnan(DEFL(:))==1
	    AFLAG = inf;
    else
        AFLAG = 1;
    end

end
