%wing shear flow
clear all;
close all;


Vx = 1; Vz = 1; My = 1;  %test loads will be applied individually

%bending moment test loads
Mx = 1; Mz = 1;

%define a few 
numTopStringers = 8;
numBottomStringers = 5;
numNoseTopStringers = 2;
numNoseBottomStringers = 3;

t_upper = 0.02/12;
t_lower = 0.02/12;
t_upper_front = 0.02/12; 
t_lower_front = 0.02/12;
t_frontSpar = 0.04/12;
t_rearSpar = 0.04/12;

frontSpar = 0.3;
backSpar = 0.75;
chord = 5;


sparCaps(1).posX = frontSpar*chord;
sparCaps(2).posX = frontSpar*chord;
sparCaps(3).posX = backSpar*chord;
sparCaps(4).posX = backSpar*chord;

sparCaps(1).posZ = get_z(frontSpar,1)*chord;
sparCaps(2).posZ = get_z(frontSpar,0)*chord;
sparCaps(3).posZ = get_z(backSpar,1)*chord;
sparCaps(4).posZ = get_z(backSpar,0)*chord;

sparCaps(1).area = .1;
sparCaps(2).area = .1;
sparCaps(3).area = .1;
sparCaps(4).area = .1;

upperStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numTopStringers + 1);
lowerStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numBottomStringers + 1);
upperNoseStringerGap = (sparCaps(1).posX - 0)/(numNoseTopStringers + 1);
lowerNoseStringerGap = (sparCaps(1).posX - 0)/(numNoseBottomStringers + 1);


%set stringers spaced evenly along X axis betwen Spars
%top Stringers
for i=1:numTopStringers
    topStringers(i).posX = sparCaps(1).posX + upperStringerGap*i;
    topStringers(i).posZ = get_z(topStringers(i).posX/chord,1)*chord;
    topStringers(i).area = .1;
end

%bottom Stringers
for i=1:numBottomStringers
    bottomStringers(i).posX = sparCaps(4).posX - lowerStringerGap*i;
    bottomStringers(i).posZ = get_z(bottomStringers(i).posX/chord,0)*chord;
    bottomStringers(i).area = .1;

end

%nose bottom Stringers
for i=1:numNoseBottomStringers
    noseBottomStringers(i).posX = sparCaps(2).posX - lowerNoseStringerGap*i;
    noseBottomStringers(i).posZ = get_z(noseBottomStringers(i).posX/chord,0)*chord;
    noseBottomStringers(i).area = .1;
end

%nose top Stringers
for i=1:numNoseTopStringers
    noseTopStringers(i).posX = upperNoseStringerGap*i;
    noseTopStringers(i).posZ = get_z(noseTopStringers(i).posX/chord,1)*chord;
    noseTopStringers(i).area = .1;
end


centroid.posX = sum([sparCaps.posX].*[sparCaps.area]) + ...
    sum([topStringers.posX].*[topStringers.area]) + ...
    sum([bottomStringers.posX].*[bottomStringers.area]) + ...
    sum([noseTopStringers.posX].*[noseTopStringers.area]) + ...
    sum([noseBottomStringers.posX].*[noseBottomStringers.area]);

centroid.posX = centroid.posX / ( sum([sparCaps.area]) + sum([topStringers.area]) + ...
    sum([bottomStringers.area]) + sum([noseTopStringers.area]) + sum([noseBottomStringers.area]));

centroid.posZ = sum([sparCaps.posZ].*[sparCaps.area]) + ...
    sum([topStringers.posZ].*[topStringers.area]) + ...
    sum([bottomStringers.posZ].*[bottomStringers.area]) + ...
    sum([noseTopStringers.posZ].*[noseTopStringers.area]) + ...
    sum([noseBottomStringers.posZ].*[noseBottomStringers.area]);
centroid.posZ = centroid.posZ / ( sum([sparCaps.area]) + sum([topStringers.area]) + ...
    sum([bottomStringers.area]) + sum([noseTopStringers.area]) + sum([noseBottomStringers.area]));




%summing contributions for inertia terms
Ix = 0; Iz = 0; Ixz = 0;

for i=1:4 %spar caps
    Ix = Ix + sparCaps(i).area*(sparCaps(i).posZ-centroid.posZ)^2;
    Iz = Iz + sparCaps(i).area*(sparCaps(i).posX-centroid.posX)^2;
    Ixz = Ixz + sparCaps(i).area*(sparCaps(i).posX-centroid.posX)*(sparCaps(i).posZ-centroid.posZ);
end



for i=1:numTopStringers %top stringers
    Ix = Ix + topStringers(i).area*(topStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + topStringers(i).area*(topStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + topStringers(i).area*(topStringers(i).posX-centroid.posX)*(topStringers(i).posZ-centroid.posZ);
end
for i=1:numBottomStringers %bottom stringers
    Ix = Ix + bottomStringers(i).area*(bottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + bottomStringers(i).area*(bottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + bottomStringers(i).area*(bottomStringers(i).posX-centroid.posX)*(bottomStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseTopStringers %nose top stringers
    Ix = Ix + noseTopStringers(i).area*(noseTopStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + noseTopStringers(i).area*(noseTopStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + noseTopStringers(i).area*(noseTopStringers(i).posX-centroid.posX)*(noseTopStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseBottomStringers %nose bottom stringers
    Ix = Ix + noseBottomStringers(i).area*(noseBottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + noseBottomStringers(i).area*(noseBottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + noseBottomStringers(i).area*(noseBottomStringers(i).posX-centroid.posX)*(noseBottomStringers(i).posZ-centroid.posZ);
end

%Ixz = -Ixz;

%define webs

%% web cell 1
web = [];
%upper webs
numStringers = numTopStringers;
stringerGap = upperStringerGap;
webThickness = t_upper;
tempStringers = topStringers;

for i=1:(numStringers+1)
    web(i).xStart = sparCaps(1).posX + stringerGap*(i-1);
    web(i).xEnd = sparCaps(1).posX + stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart/chord,1)*chord;
    web(i).zEnd = get_z(web(i).xEnd/chord,1)*chord;
    if i==1
        web(i).dp_area = sparCaps(1).area;
        web(i).dP_X = 0;
        web(i).dP_Z = 0;
        web(i).qPrime_X = 0;
        web(i).qPrime_Z = 0;
    else
        web(i).dp_area = tempStringers(i-1).area;
        dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);  %just Vx
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);  %just Vz
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z;
    end
    tempInt = get_int(web(i).xStart/chord,web(i).xEnd/chord,1)*chord^2;  %integral of airfoil function
    triangle1 = abs( (web(i).xStart - sparCaps(1).posX)*web(i).zStart/2);
    triangle2 = abs((web(i).xEnd - sparCaps(1).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart/chord,web(i).xEnd/chord,1)*chord;
    web(i).dS_over_t = web(i).ds / web(i).thickness;
    
    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
    web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);
end
webTop = web;
web = [];

%rear spar
i=1;
web(i).xStart = sparCaps(3).posX;
web(i).xEnd = sparCaps(4).posX;
web(i).thickness = t_rearSpar;
web(i).zStart = sparCaps(3).posZ;
web(i).zEnd = sparCaps(4).posZ;
web(i).dp_area = sparCaps(3).area;
dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime_X = webTop(numTopStringers+1).qPrime_X - web(i).dP_X;
web(i).qPrime_Z = webTop(numTopStringers+1).qPrime_Z - web(i).dP_Z;

web(i).Area = (sparCaps(3).posX-sparCaps(1).posX)*sparCaps(3).posZ/2 + ...
    abs((sparCaps(3).posX-sparCaps(1).posX)*sparCaps(4).posZ/2);
web(i).ds = abs(sparCaps(3).posZ - sparCaps(4).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;

web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

webRearSpar = web;
web = [];


%lower webs
numStringers = numBottomStringers;
stringerGap = lowerStringerGap;
webThickness = t_lower;
tempStringers = bottomStringers;

for i=1:(numStringers+1)
    web(i).xStart = sparCaps(4).posX - stringerGap*(i-1);
    web(i).xEnd = sparCaps(4).posX - stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart/chord,0)*chord;
    web(i).zEnd = get_z(web(i).xEnd/chord,0)*chord;
    dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
    if i==1
        web(i).dp_area = sparCaps(4).area;
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime_X = webRearSpar.qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = webRearSpar.qPrime_Z - web(i).dP_Z;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP_X = get_dp(dx,dz, Vx,0,Ix,Iz,Ixz,web(i).dp_area);
        web(i).dP_Z = get_dp(dx,dz, 0,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z;
    end
    
    tempInt = get_int(web(i).xEnd/chord,web(i).xStart/chord,0)*chord^2;  %integral of airfoil function
    triangle2 = abs((web(i).xStart - sparCaps(1).posX)*web(i).zStart/2);
    triangle1 = abs((web(i).xEnd - sparCaps(1).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart/chord,web(i).xEnd/chord,0)*chord;
    web(i).dS_over_t = web(i).ds / web(i).thickness;

    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X*(web(i).xEnd-web(i).xStart);
    web(i).qp_dx_Z = web(i).qPrime_Z*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X*(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z*(web(i).zEnd-web(i).zStart);

    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webBottom = web;
web = [];

%front Spar
i=1;
web(i).xStart = sparCaps(2).posX;
web(i).xEnd = sparCaps(1).posX;
web(i).thickness = t_frontSpar;
web(i).zStart = sparCaps(2).posZ;
web(i).zEnd = sparCaps(1).posZ;
web(i).dp_area = sparCaps(2).area;
dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime_X = webBottom(numBottomStringers+1).qPrime_X - web(i).dP_X;
web(i).qPrime_Z = webBottom(numBottomStringers+1).qPrime_Z - web(i).dP_Z;
web(i).Area = 0;
web(i).ds = abs(sparCaps(2).posZ - sparCaps(1).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;

web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

webFrontSpar = web;
web = [];




%% web cell 2

%lower nose webs
numStringers = numNoseBottomStringers;
stringerGap = lowerNoseStringerGap;
webThickness = t_lower_front;
tempStringers = noseBottomStringers;

for i=1:(numStringers+1)
    web(i).xStart = sparCaps(2).posX - stringerGap*(i-1);
    web(i).xEnd = sparCaps(2).posX - stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart/chord,0)*chord;
    web(i).zEnd = get_z(web(i).xEnd/chord,0)*chord;
    dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;

    if i==1
        web(i).dp_area = sparCaps(2).area;
        web(i).dP_X = 0;
        web(i).dP_Z = 0;
        web(i).qPrime_X = 0;
        web(i).qPrime_Z = 0;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z;
    end
    tempInt = get_int(web(i).xEnd/chord,web(i).xStart/chord,0)*chord^2;  %integral of airfoil function
    triangle1 = abs((web(i).xStart - sparCaps(2).posX)*web(i).zStart/2);
    triangle2 = abs((web(i).xEnd - sparCaps(2).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart/chord,web(i).xEnd/chord,0)*chord;
    web(i).dS_over_t = web(i).ds / web(i).thickness;

    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
    web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webLowerNose = web;
web = [];

%upper nose webs
numStringers = numNoseTopStringers;
stringerGap = upperNoseStringerGap;
webThickness = t_upper_front;
tempStringers = noseTopStringers;

for i=1:(numStringers+1)
    web(i).xStart = stringerGap*(i-1);
    web(i).xEnd = stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart/chord,1)*chord;
    web(i).zEnd = get_z(web(i).xEnd/chord,1)*chord;
    dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;
    if i==1
        web(i).dp_area = 0;
        web(i).dP_X = 0;
        web(i).dP_Z = 0;
        web(i).qPrime_X = webLowerNose(numNoseBottomStringers+1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = webLowerNose(numNoseBottomStringers+1).qPrime_Z - web(i).dP_Z;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
        web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime_X = web(i-1).qPrime_X - web(i).dP_X;
        web(i).qPrime_Z = web(i-1).qPrime_Z - web(i).dP_Z;
    end
    tempInt = get_int(web(i).xStart/chord,web(i).xEnd/chord,1)*chord^2;  %integral of airfoil function
    triangle2 = abs((web(i).xStart - sparCaps(2).posX)*web(i).zStart/2);
    triangle1 = abs((web(i).xEnd - sparCaps(2).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart/chord,web(i).xEnd/chord,1)*chord;
    web(i).dS_over_t = web(i).ds / web(i).thickness;

    web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
    web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
    web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
    web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
    web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
    web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
    web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
    web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

end
webUpperNose = web;
web = [];


%front Spar
i=1;
web(i).xStart = sparCaps(1).posX;
web(i).xEnd = sparCaps(2).posX;
web(i).thickness = t_frontSpar;
web(i).zStart = sparCaps(1).posZ;
web(i).zEnd = sparCaps(2).posZ;
web(i).dp_area = sparCaps(1).area;
dx = web(i).xStart-centroid.posX; dz = web(i).zStart-centroid.posZ;

web(i).dP_X = get_dp(dx,dz,Vx,0,Ix,Iz,Ixz,web(i).dp_area);
web(i).dP_Z = get_dp(dx,dz,0,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime_X = webUpperNose(numNoseTopStringers+1).qPrime_X - web(i).dP_X;
web(i).qPrime_Z = webUpperNose(numNoseTopStringers+1).qPrime_Z - web(i).dP_Z;
web(i).Area = 0;
web(i).ds = abs(sparCaps(1).posZ - sparCaps(2).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;
web(i).q_dS_over_t_X = web(i).qPrime_X * web(i).dS_over_t;
web(i).q_dS_over_t_Z = web(i).qPrime_Z * web(i).dS_over_t;
web(i).two_A_qprime_X = 2*web(i).Area*web(i).qPrime_X;
web(i).two_A_qprime_Z = 2*web(i).Area*web(i).qPrime_Z;
web(i).qp_dx_X = web(i).qPrime_X *(web(i).xEnd-web(i).xStart);
web(i).qp_dx_Z = web(i).qPrime_Z *(web(i).xEnd-web(i).xStart);
web(i).qp_dz_X = web(i).qPrime_X *(web(i).zEnd-web(i).zStart);
web(i).qp_dz_Z = web(i).qPrime_Z *(web(i).zEnd-web(i).zStart);

webFrontSparCell2 = web;
web = [];


%check that q'*dx sums up to Vx
disp(' Checking force due to Vx = 1:')
Fx = sum([webTop.qp_dx_X])+webRearSpar.qp_dx_X+ sum([webBottom.qp_dx_X])+webFrontSpar.qp_dx_X;  %cell 1
Fx = Fx + sum([webLowerNose.qp_dx_X])+ sum([webUpperNose.qp_dx_X]);  %cell 2

Fz = sum([webTop.qp_dz_X])+webRearSpar.qp_dz_X+ sum([webBottom.qp_dz_X])+webFrontSpar.qp_dz_X;  %cell 1
Fz = Fz + sum([webLowerNose.qp_dz_X])+ sum([webUpperNose.qp_dz_X]);  %cell 2
disp([' Fx = ' num2str(Fx) '  Fz = ' num2str(Fz)]) 

%check that q'*dz sums up to Vz


disp(' Checking force due to Vz = 1:')
Fx = sum([webTop.qp_dx_Z])+webRearSpar.qp_dx_Z+ sum([webBottom.qp_dx_Z])+webFrontSpar.qp_dx_Z;  %cell 1
Fx = Fx + sum([webLowerNose.qp_dx_Z])+ sum([webUpperNose.qp_dx_Z]);  %cell 2
Fz = sum([webTop.qp_dz_Z])+webRearSpar.qp_dz_Z+ sum([webBottom.qp_dz_Z])+webFrontSpar.qp_dz_Z;  %cell 1
Fz = Fz + sum([webLowerNose.qp_dz_Z])+ sum([webUpperNose.qp_dz_Z]);  %cell 2
disp([' Fx = ' num2str(Fx) '  Fz = ' num2str(Fz)]) 

%%

% sum up the ds/t and  q*ds/t to solve 2 equations, 2 unknowns

% [A]*[q1s q2s] = B

A11 = sum([webTop.dS_over_t])+webRearSpar.dS_over_t+ sum([webBottom.dS_over_t])+webFrontSpar.dS_over_t;
A22 = sum([webLowerNose.dS_over_t])+ sum([webUpperNose.dS_over_t])+webFrontSparCell2.dS_over_t;
A12 = -webFrontSpar.dS_over_t;
A21 = -webFrontSparCell2.dS_over_t;

B1_X = sum([webTop.q_dS_over_t_X])+webRearSpar.q_dS_over_t_X+ sum([webBottom.q_dS_over_t_X])+webFrontSpar.q_dS_over_t_X;
B2_X = sum([webLowerNose.q_dS_over_t_X])+ sum([webUpperNose.q_dS_over_t_X])+webFrontSparCell2.q_dS_over_t_X;
B1_Z = sum([webTop.q_dS_over_t_Z])+webRearSpar.q_dS_over_t_Z+ sum([webBottom.q_dS_over_t_Z])+webFrontSpar.q_dS_over_t_Z;
B2_Z = sum([webLowerNose.q_dS_over_t_Z])+ sum([webUpperNose.q_dS_over_t_Z])+webFrontSparCell2.q_dS_over_t_Z;

Amat = [A11 A12; A21 A22];
Bmat_X = -[B1_X;B2_X];
Bmat_Z = -[B1_Z;B2_Z];

qs_X = inv(Amat)*Bmat_X;
qs_Z = inv(Amat)*Bmat_Z;



sum_2_a_q_X = sum([webTop.two_A_qprime_X])+webRearSpar.two_A_qprime_X+ sum([webBottom.two_A_qprime_X]);  %cell 1 qprimes
sum_2_a_q_X = sum_2_a_q_X + sum([webLowerNose.two_A_qprime_X])+ sum([webUpperNose.two_A_qprime_X]);   %cell 2 qprimes
sum_2_a_q_X = sum_2_a_q_X +  2*qs_X(1)*(sum([webTop.Area])+webRearSpar.Area+ sum([webBottom.Area]));
sum_2_a_q_X = sum_2_a_q_X +  2*qs_X(2)*(sum([webLowerNose.Area])+ sum([webUpperNose.Area]));

sum_2_a_q_Z = sum([webTop.two_A_qprime_Z])+webRearSpar.two_A_qprime_Z+ sum([webBottom.two_A_qprime_Z]);  %cell 1 qprimes
sum_2_a_q_Z = sum_2_a_q_Z + sum([webLowerNose.two_A_qprime_Z])+ sum([webUpperNose.two_A_qprime_Z]);   %cell 2 qprimes
sum_2_a_q_Z = sum_2_a_q_Z +  2*qs_Z(1)*(sum([webTop.Area])+webRearSpar.Area+ sum([webBottom.Area]));
sum_2_a_q_Z = sum_2_a_q_Z +  2*qs_Z(2)*(sum([webLowerNose.Area])+ sum([webUpperNose.Area]));

%shear center
sc.posX =  sum_2_a_q_Z / Vz + frontSpar*chord;
sc.posZ = - sum_2_a_q_X / Vx;

disp(['shear center: x = ' num2str(sc.posX) ' z = ' num2str(sc.posZ)]);

% now consider the torque representing shifting the load from the quarter
% chord to the SC  (need to check signs on these moments)
 
torque_Z = Vz*(sc.posX - 0.25*chord);
torque_X = -Vx*sc.posZ;
torque_Y = My;


Area1 = sum([webTop.Area]) + webRearSpar.Area + sum([webBottom.Area]);
%check area
Area1_check = get_int(frontSpar,backSpar,1)*chord^2 + get_int(frontSpar,backSpar,0)*chord^2;

Area2 = sum([webLowerNose.Area]) + sum([webUpperNose.Area]);
Area2_check = get_int(0,frontSpar,1)*chord^2 + get_int(0,frontSpar,0)*chord^2;


%for twist equation  (see excel spreadsheet example)

q1t_over_q2t = (A22/Area2 + webFrontSpar.dS_over_t/Area1)/(A11/Area1 + webFrontSpar.dS_over_t/Area2);

q2t = torque_X/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;
qt_X = [q1t;q2t];

q2t = torque_Z/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;
qt_Z = [q1t;q2t];

q2t = torque_Y/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;
qt_Y = [q1t;q2t];



% --- - add up all shear flows:  qtot = (qPrime + qs) + qt

webNames = {'webTop' 'webBottom' 'webRearSpar' 'webUpperNose' 'webLowerNose' 'webFrontSpar'};

for j=1:length(webNames)
    webTemp = [];
    webTemp = eval(webNames{j});
    for i=1:length(webTemp)
        if j<4  %cell 1 webs
            webTemp(i).qtot_Z = webTemp(i).qPrime_Z + qs_Z(1) + qt_Z(1);
            webTemp(i).qtot_X = webTemp(i).qPrime_X + qs_X(1) + qt_X(1);
            webTemp(i).qtot_Y = qt_Y(1);
        elseif j<6  %cell 2 webs
            webTemp(i).qtot_Z = webTemp(i).qPrime_Z + qs_Z(2) + qt_Z(2);
            webTemp(i).qtot_X = webTemp(i).qPrime_X + qs_X(2) + qt_X(2);
            webTemp(i).qtot_Y = qt_Y(2);
            
        else  %shared web
            webTemp(i).qtot_Z = webTemp(i).qPrime_Z + qs_Z(1) + qt_Z(1) - qs_Z(2) - qt_Z(2);
            webTemp(i).qtot_X = webTemp(i).qPrime_X + qs_X(1) + qt_X(1) - qs_X(2) - qt_X(2);
            webTemp(i).qtot_Y = qt_Y(1)-qt_Y(2);
        end
        webTemp(i).qtot_dz_Z = webTemp(i).qtot_Z*(webTemp(i).zEnd-webTemp(i).zStart);
        webTemp(i).qtot_dx_Z = webTemp(i).qtot_Z*(webTemp(i).xEnd-webTemp(i).xStart);
        webTemp(i).qtot_dz_X = webTemp(i).qtot_X*(webTemp(i).zEnd-webTemp(i).zStart);
        webTemp(i).qtot_dx_X = webTemp(i).qtot_X*(webTemp(i).xEnd-webTemp(i).xStart);
        webTemp(i).qtot_dz_Y = webTemp(i).qtot_Y*(webTemp(i).zEnd-webTemp(i).zStart);
        webTemp(i).qtot_dx_Y = webTemp(i).qtot_Y*(webTemp(i).xEnd-webTemp(i).xStart);
    end
    eval([webNames{j} ' = webTemp;']);
end



check_dz_Z = sum([webTop.qtot_dz_Z]) + sum([webBottom.qtot_dz_Z]) + sum([webUpperNose.qtot_dz_Z]) + ...
    sum([webLowerNose.qtot_dz_Z]) + webRearSpar.qtot_dz_Z + webFrontSpar.qtot_dz_Z
check_dx_Z = sum([webTop.qtot_dx_Z]) + sum([webBottom.qtot_dx_Z]) + sum([webUpperNose.qtot_dx_Z]) + ...
    sum([webLowerNose.qtot_dx_Z]) + webRearSpar.qtot_dx_Z + webFrontSpar.qtot_dx_Z

check_dz_X = sum([webTop.qtot_dz_X]) + sum([webBottom.qtot_dz_X]) + sum([webUpperNose.qtot_dz_X]) + ...
    sum([webLowerNose.qtot_dz_X]) + webRearSpar.qtot_dz_X + webFrontSpar.qtot_dz_X
check_dx_X = sum([webTop.qtot_dx_X]) + sum([webBottom.qtot_dx_X]) + sum([webUpperNose.qtot_dx_X]) + ...
    sum([webLowerNose.qtot_dx_X]) + webRearSpar.qtot_dx_X + webFrontSpar.qtot_dx_X

check_dz_Y = sum([webTop.qtot_dz_Y]) + sum([webBottom.qtot_dz_Y]) + sum([webUpperNose.qtot_dz_Y]) + ...
    sum([webLowerNose.qtot_dz_Y]) + webRearSpar.qtot_dz_Y + webFrontSpar.qtot_dz_Y
check_dx_Y = sum([webTop.qtot_dx_Y]) + sum([webBottom.qtot_dx_Y]) + sum([webUpperNose.qtot_dx_Y]) + ...
    sum([webLowerNose.qtot_dx_Y]) + webRearSpar.qtot_dx_Y + webFrontSpar.qtot_dx_Y

%adding bending stress to stringer structures

for i=1:length(topStringers)
    dx = topStringers(i).posX-centroid.posX; dz = topStringers(i).posZ-centroid.posZ;
    topStringers(i).sigY_Mx = get_sigY(dx,dz,Mx,0,Ix,Iz,Ixz);
    topStringers(i).sigY_Mz = get_sigY(dx,dz,0,Mz,Ix,Iz,Ixz);
end
for i=1:length(bottomStringers)
    dx = bottomStringers(i).posX-centroid.posX; dz = bottomStringers(i).posZ-centroid.posZ;
    bottomStringers(i).sigY_Mx = get_sigY(dx,dz,Mx,0,Ix,Iz,Ixz);
    bottomStringers(i).sigY_Mz = get_sigY(dx,dz,0,Mz,Ix,Iz,Ixz);
end
for i=1:length(noseTopStringers)
    dx = noseTopStringers(i).posX-centroid.posX; dz = noseTopStringers(i).posZ-centroid.posZ;
    noseTopStringers(i).sigY_Mx = get_sigY(dx,dz,Mx,0,Ix,Iz,Ixz);
    noseTopStringers(i).sigY_Mz = get_sigY(dx,dz,0,Mz,Ix,Iz,Ixz);
end
for i=1:length(noseBottomStringers)
    dx = noseBottomStringers(i).posX-centroid.posX; dz = noseBottomStringers(i).posZ-centroid.posZ;
    noseBottomStringers(i).sigY_Mx = get_sigY(dx,dz,Mx,0,Ix,Iz,Ixz);
    noseBottomStringers(i).sigY_Mz = get_sigY(dx,dz,0,Mz,Ix,Iz,Ixz);
end




%plotting airfoil cross-section

xChord = 0:.01:1;
xChord = xChord*chord;
upperSurface = zeros(1,length(xChord));
lowerSurface = zeros(1,length(xChord));

for i=1:length(xChord)
    upperSurface(i) = get_z(xChord(i)/chord,1)*chord;
    lowerSurface(i) = get_z(xChord(i)/chord,0)*chord;
end

figure; hold on; axis equal; grid on;
%plot(xChord,z_camber,'-')
plot(xChord,upperSurface,'-k','linewidth',2)
plot(xChord,lowerSurface,'-k','linewidth',2)
plot([0 1],[0 0],'--k','linewidth',1)


for i = 1:length(webTop)
   vecX = [frontSpar*chord webTop(i).xStart webTop(i).xEnd];
   vecZ = [0 webTop(i).zStart webTop(i).zEnd];
   fill(vecX,vecZ,[0.9 0.9 0.9])
end

for i = 1:length(webBottom)
   vecX = [frontSpar*chord webBottom(i).xStart webBottom(i).xEnd];
   vecZ = [0 webBottom(i).zStart webBottom(i).zEnd];
   fill(vecX,vecZ,[0.9 0.9 0.9])
end

for i = 1:length(webUpperNose)
   vecX = [frontSpar*chord webUpperNose(i).xStart webUpperNose(i).xEnd];
   vecZ = [0 webUpperNose(i).zStart webUpperNose(i).zEnd];
   fill(vecX,vecZ,[0.7 0.9 1.0])
end

for i = 1:length(webLowerNose)
   vecX = [frontSpar*chord webLowerNose(i).xStart webLowerNose(i).xEnd];
   vecZ = [0 webLowerNose(i).zStart webLowerNose(i).zEnd];
   fill(vecX,vecZ,[0.7 0.9 1.0])
end

   vecX = [frontSpar*chord sparCaps(3).posX sparCaps(4).posX];
   vecZ = [0 sparCaps(3).posZ sparCaps(4).posZ];
   fill(vecX,vecZ,[0.9 0.9 0.9])


sparCapSize = 18;
stringerSize = 18;
plot([sparCaps(1).posX sparCaps(2).posX],[sparCaps(1).posZ sparCaps(2).posZ],'-k','linewidth',2)
plot([sparCaps(3).posX sparCaps(4).posX],[sparCaps(3).posZ sparCaps(4).posZ],'-k','linewidth',2)
plot([sparCaps.posX],[sparCaps.posZ],'.b','markersize',sparCapSize)
plot([topStringers.posX],[topStringers.posZ],'.r','markersize',stringerSize)
plot([bottomStringers.posX],[bottomStringers.posZ],'.r','markersize',stringerSize)
plot([noseTopStringers.posX],[noseTopStringers.posZ],'.r','markersize',stringerSize)
plot([noseBottomStringers.posX],[noseBottomStringers.posZ],'.r','markersize',stringerSize)
plot(centroid.posX,centroid.posZ,'.k','markerSize',18)
plot(sc.posX,sc.posZ,'.g','markersize',18)





