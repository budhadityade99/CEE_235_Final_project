function [gamma] = BDLT_etran(coordi,coordj,webdir)

    L = sqrt((coordi(1)-coordj(1))^2+(coordi(2)-coordj(2))^2+(coordi(3)-coordj(3))^2); % element length
    cosiney = webdir; %direction cosine of local y axis with respect to global coordinate system
    cosinex = [(coordj(1)-coordi(1))/L (coordj(2)-coordi(2))/L (coordj(3)-coordi(3))/L]; %direction cosine of local x axis with respect to global coordinate system
    cosinez = cross(cosinex,cosiney); %direction cosine of local z axis with respect to global coordinate system
    gamma1 = [cosinex; cosiney; cosinez]; % gamma matrix
    gamma =zeros(12,12); % initialize transformation matrix with zeros
    %Populate transformation matrix with appropriate direction cosines
    gamma([1 2 3],[1 2 3]) = gamma1;
    gamma([4 5 6],[4 5 6]) = gamma1;
    gamma([7 8 9],[7 8 9]) = gamma1;
    gamma([10 11 12],[10 11 12]) = gamma1;
end