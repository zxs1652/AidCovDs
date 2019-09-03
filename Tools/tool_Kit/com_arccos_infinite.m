function [k_Later_IJ,k_Later_II,k_Later_JJ] = com_arccos_infinite(k_Prior_IJ,k_Prior_II,k_Prior_JJ,temp_R)
    theta_Prior = acos(k_Prior_IJ/sqrt(k_Prior_II*k_Prior_JJ));
    sin_Theta = sin(theta_Prior);
    cos_Theta = cos(theta_Prior);
    switch temp_R
        case 0
            j_Res = pi - theta_Prior;
        case 1
            j_Res = (pi - theta_Prior)*cos_Theta - sin_Theta;
        case 2
            j_Res = (pi - theta_Prior)*(1 + 2*cos_Theta^2) + 3*sin_Theta*cos_Theta;
        case 3
            j_Res = (pi - theta_Prior)*(9*sin_Theta^2*cos_Theta+15*cos_Theta^3) + 4*sin_Theta^3 + 15*sin_Theta*cos_Theta^2;
    end            
    k_Later_IJ = 1/pi*(k_Prior_II*k_Prior_JJ)^(temp_R/2)*j_Res;
    k_Later_II = 1/pi*(k_Prior_II*k_Prior_II)^(temp_R/2)*j_Res;
    k_Later_JJ = 1/pi*(k_Prior_JJ*k_Prior_JJ)^(temp_R/2)*j_Res;

end