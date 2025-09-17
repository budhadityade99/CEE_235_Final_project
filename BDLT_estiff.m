function [elk] = BDLT_estiff(A,Izz,Iyy,J,Ayy,Azz,E,v,L,release)
        % condition when Iyy or J not provided to prevent unstable system
        % error
        if Iyy==0 || J==0
            Iyy=1e-7*Izz;
            J=1e-7*Izz;
        end
        G = E/(2*(1+v)); % Shear Modulus
        elk = zeros(12,12); % initialize local element stiffness matrix with zeros
        %Populating local element stifness matrix with axial stiffness terms
        elk(:,1) = [E*A/L;0;0;0;0;0;-E*A/L;0;0;0;0;0];
        elk(:,7) = -[E*A/L;0;0;0;0;0;-E*A/L;0;0;0;0;0];
        %Populating local element stifness matrix with torsional stiffness terms
        elk(:,4) = [0;0;0;G*J/L;0;0;0;0;0;-G*J/L;0;0;];
        elk(:,10) = -[0;0;0;G*J/L;0;0;0;0;0;-G*J/L;0;0;];
        %Initialize eta terms to incorporate shear deformations along local y and z axes
        nz = E*Izz/(Ayy*G); ny = E*Iyy/(Azz*G);
        %Create bending stiffness terms with modifications for shear deformations
        bnd_stiffzz = (E*Izz/(L*(L^2/12+nz)))*[1 L/2 -1 L/2; L/2 (L^2/3+nz) -L/2 (L^2/6-nz); -1 -L/2 1 -L/2; L/2 (L^2/6-nz) -L/2 (L^2/3+nz);];
        bnd_stiffyy = (E*Iyy/(L*(L^2/12+ny)))*[1 -L/2 -1 -L/2; -L/2 (L^2/3+ny) L/2 (L^2/6-ny); -1 L/2 1 L/2; -L/2 (L^2/6-ny) L/2 (L^2/3+ny);];
        %Populating local element stifness matrix with bending stiffness terms
        elk([2 6 8 12],[2 6 8 12]) = bnd_stiffzz;
        elk([3 5 9 11],[3 5 9 11]) = bnd_stiffyy;
    if release == [0 0]
        % Rigid member
        elk = elk;
    elseif release == [1 0]
        % Flexural release at starting node (ith node) using static condensation
        % row and column operations to send the 5th and 6th rows and columns
        % (My and Mz) to the end of the matrix for static condensation
        elk_reshaped1(1:4,:) = elk(1:4,:);
        elk_reshaped1(5:10,:) = elk(7:12,:);
        elk_reshaped1(11:12,:) = elk(5:6,:);
        elk_reshaped2(:,1:4) = elk_reshaped1(:,1:4);
        elk_reshaped2(:,5:10) = elk_reshaped1(:,7:12);
        elk_reshaped2(:,11:12) = elk_reshaped1(:,5:6);
        elk1 = elk_reshaped2(1:10,1:10)-elk_reshaped2(1:10,11:12)*inv(elk_reshaped2(11:12,11:12))*elk_reshaped2(11:12,1:10);
        % Populating the condensed matrix in the original 12x12 form
        elk1_r = [elk1(1:4,:);zeros(2,size(elk1,2));elk1(5:10,:)];
        elk = [elk1_r(:,1:4) zeros(size(elk1_r,1),2) elk1_r(:,5:10)];
    elseif release == [0 1]
        % Flexural release at ending node (jth node) using static condensation        
        elk2 = elk(1:10,1:10)-elk(1:10,11:12)*inv(elk(11:12,11:12))*elk(11:12,1:10);
        % Populating the condensed matrix in the original 12x12 form
        elk2_r = [elk2;zeros(2,size(elk2,2))];
        elk = [elk2_r zeros(size(elk2_r,1),2)];
    else
        % Flexural releases at both nodes (truss member)
        elk = zeros(12,12);
        elk([1,7],[1,7]) = (A*E/L)*[1 -1;-1 1];
        elk([4,10],[4,10]) = (G*J/L)*[1 -1;-1 1];
    end
end