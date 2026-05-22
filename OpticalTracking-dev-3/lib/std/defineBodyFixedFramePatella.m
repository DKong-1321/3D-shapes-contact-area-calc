function gTp0 = defineBodyFixedFramePatella(patella,right)
med = patella.medial.map(@(x) x.translations_mean);
lat = patella.lateral.map(@(x) x.translations_mean);
dist = patella.distal.map(@(x) x.translations_mean);
inf = patella.inferior.map(@(x) x.translations_mean);
sup = patella.superior.map(@(x) x.translations_mean);
if med.is_none || lat.is_none
    gTp0 = [];
    return
end
med = med.unwrap();
lat = lat.unwrap();
try
    dist = dist.unwrap();
catch
    inf = inf.unwrap();
    sup = sup.unwrap();
end

%defines a body fixed frame for the patella.  
%my notation, _ implies it is a vector e.g. xa_ is the direction vector of
%the x-axis

angle = @(u_,v_) acosd(dot(u_,v_)/(norm(u_,2)*norm(v_,2))); %define a function to calculate the angle between two vectors
ucross=@(u_,v_) cross(u_,v_)/norm(cross(u_,v_),2); %define function to find unit cross product
uvector=@(a,b) (b-a)/norm(b-a,2); %define a function to find a unit vector from a to b

%x, y and z notation of grood and suntay 1983, however the body frame is
%one decided upon by RvA for a particular experiment in Nov 2018 rather
%than the body frame defined by Grood and Sunray (as the points necessary
%to contructs G&S body frame such as the ankle joint centre were not
%available)

%x = Medial lateral = epicondylar axis, positive to the right, thus is
%orientated laterally in right knee, and medially in the left knee
%y = Anterior posterior, anterior is positive
%z = Proximal Distal, proximal is positive
%i,j,k correspond to unit base vectors in the x, y and z direactions

origin = (med+lat)/2;
if right
    i_ = uvector(med,lat)'; %RIGHT KNEE, x Axis
else
    i_ = uvector(lat,med)'; %LEFT KNEE, x Axis
end

if dist
    tempk_= uvector(dist,origin)'; %the distal point is approximate and thus this axis is not necessarily perpendicular to epicondylar axis
else
    tempk_= uvector(inf,sup)'; %the inferior (distal) point is approximate and thus this axis is not necessarily perpendicular to epicondylar axis
end
j_ = ucross(tempk_,i_); % y-axis
k_ = ucross(i_,j_); %recalculate k so perpendicular to give orthogonal coordinate system.

rot=[i_,j_,k_];

rot=[rot,[0 0 0]';0 0 0 1];


trans=[1 0 0 origin(1);
       0 1 0 origin(2);
       0 0 1 origin(3)
       0 0 0 1];


gTp0=trans*rot; %same as gTp0=[rot,origin';0 0 0 1];

%check for orthogonality
% angle(i_,j_)
% angle(j_,k_)
% angle(i_,k_)


end
