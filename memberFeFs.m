function [FeF] = memberFeFs(w,L,release)

    FeF = zeros(12,1);
    if release == [0 0] % No flexural release
        % FeFs for axially distributed load
        FeF(1,1) = -w(1)*L/2; FeF(7,1) = -w(1)*L/2;
        % FeFs for distributed load along local y-axis
        FeF(2,1) = -w(2)*L/2; FeF(6,1) = -w(2)*L^2/12; FeF(8,1) = -w(2)*L/2; FeF(12,1) = w(2)*L^2/12;
        % FeFs for distributed load along local z-axis
        FeF(3,1) = -w(3)*L/2; FeF(5,1) = w(3)*L^2/12; FeF(9,1) = -w(3)*L/2; FeF(11,1) = -w(3)*L^2/12;
    elseif release == [1 0] % Flexural release at ith (starting node)
        % FeFs for axially distributed load
        FeF(1,1) = -w(1)*L/2; FeF(7,1) = -w(1)*L/2;
        % FeFs for distributed load along local y-axis
        FeF(2,1) = -3*w(2)*L/8; FeF(8,1) = -5*w(2)*L/8; FeF(12,1) = w(2)*L^2/8;
        % FeFs for distributed load along local z-axis
        FeF(3,1) = -3*w(3)*L/8; FeF(9,1) = -5*w(3)*L/8; FeF(11,1) = -w(3)*L^2/8;
    elseif release == [0 1] % Flexural release at jth node (end node)
        % FeFs for axially distributed load
        FeF(1,1) = -w(1)*L/2; FeF(7,1) = -w(1)*L/2;
        % FeFs for distributed load along local y-axis
        FeF(2,1) = -5*w(2)*L/8; FeF(6,1) = -w(2)*L^2/8; FeF(8,1) = -3*w(2)*L/8;
        % FeFs for distributed load along local z-axis
        FeF(3,1) = -5*w(3)*L/8; FeF(5,1) = w(3)*L^2/8; FeF(9,1) = -3*w(3)*L/8;
    else % Flexural release at both nodes (Truss member)
        FeF(1,1) = -w(1)*L/2; FeF(2,1) = -w(1)*L/2;
    end
end