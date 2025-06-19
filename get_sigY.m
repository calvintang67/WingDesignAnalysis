function z = get_sigY(xDist,zDist,Mx,Mz,Ix,Iz,Ixz)

denom = (Ix*Iz-Ixz^2);
z = -xDist*(Ix*Mz-Ixz*Mx)/denom - zDist*(Iz*Mx-Ixz*Mz)/denom;
