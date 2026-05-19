function [tf, pf] = grood_and_suntay(femur, tibia, patella, fTt, fTp, right)
    tf = tibiofemoral(fTt, femur, tibia, right);
    pf = patellofemoral(fTp, femur, patella, right);
end

function result = tibiofemoral(relative_position, femur, tibia, right)
    if isempty(relative_position)
        result = array2table([], 'VariableNames', {'flexion', 'varus', 'external', 'lateral', 'anterior', 'superior'});
        return;
    end

    result = table();

    fTt = relative_position;
    e1_=femur.i;%Femoral X axis in global reference frame, Grood and Suntay definition
    e3_=tibia.k;%Tibial z axis in global reference frame, Grood and Suntay definition
    e2_=ucross(e3_,e1_);%Floating axis in global reference frame, Grood and Suntay definition

    flexion = asind(dot(-e2_,femur.k));
    beta = acosd(dot(femur.i,tibia.k));
    if right
        external = asind(dot(-e2_,tibia.i));
        varus = 90-beta;
    else
        external = asind(dot(e2_,tibia.i));
        varus = -(90-beta);
    end

    %     H_=tibiaOrigin-femurOrigin;
    H_ = tibia.origin(1:3, :) - femur.origin(1:3, :);%translation vector of from femoral origin to tibial origin (in global coordinate frame)
    if right
        lateral = dot(H_,e1_);%projected onto the medial lateral axis e1
    else
        lateral = dot(H_,-e1_);
    end
    anterior = dot(H_,e2_);%projected onto the anterior posterior axis e2
    distal = -dot(H_,e3_);%projected onto the compression distraction axis e3, minus sign to make distraction +ve

    result.flexion = flexion';
    result.varus = varus';
    result.external = external';
    result.lateral = lateral';
    result.anterior = anterior';
    result.superior = distal';

    %Code to check between different notations
    [angles,trfi]=rotationsAndTranslations(fTt,right);

    trfi=-trfi;%which also equals gTti{i,1}\(grfi{i,1}-grti{i,1})

    B=fTt;%B notation used in equations 18-20 in grood and suntay, note they use a matrix with translation rows 2-4 in column 1, here it is column 4 rows 1-3, and rotations are similarly mapped
    q=B(1,4,:);
    if ~right
        q=-q;
    end
    q(2,:)=squeeze(B(2,4,:)).*cosd(angles(:,1)) - squeeze(B(3,4,:)) .* sind(angles(:,1));
    q(3,:)=-squeeze(pagemtimes(B(1,3,:), B(1,4,:)) + pagemtimes(B(2,3,:), B(2,4,:)) + pagemtimes(B(3,3,:), B(3,4,:)));

    q = squeeze(q)';

    % fRtCheck{i,1}=[dot(femur.i,tibia.i) dot(femur.i,tibia.j) dot(femur.i,tibia.k)
    %     dot(femur.j,tibia.i) dot(femur.j,tibia.j) dot(femur.j,tibia.k)
    %     dot(femur.k,tibia.i) dot(femur.k,tibia.j) dot(femur.k,tibia.k)];
    %
    % fRt{i,1}-fRtCheck{i,1};

end

function result = patellofemoral(relative_position, femur, patella, right)
    result = table();
    if isempty(relative_position)
        var_type = {'int8', 'int8', 'int8', 'int8', 'int8', 'int8'};
        result = table('Size', [0 6], 'VariableTypes', var_type, 'VariableNames', {'flexion', 'medial_rotation', 'lateral_tilt', 'medial_shift', 'anterior', 'superior'});
        return;
    end
    e1_=femur.i;%Femoral X axis in global reference frame, Grood and Suntay definition
    pe3_=patella.k;
    pe2_=ucross(pe3_,e1_);
    fTp = relative_position;

    [angles, ~] = rotationsAndTranslations(fTp, right);
    %                         Flexion(i,1)=asind(-dot(e2_,K_));
    %                     Flexion(i,1)=acosd(dot(J_,e2_));
    beta = acosd(dot(femur.i,patella.k));
    if right
        external = asind(dot(-pe2_,patella.i));
        %                         ExtRotation(i,1)=acosd(dot(j_,e2_));
        varus=90-beta;
    else
        external=asind(dot(pe2_,patella.i));
        %                         ExtRotation=acosd(dot(j_,e2_));
        varus=-(90-beta);
    end

    pH_ = [patella.origin(1:3,:) - femur.origin(1:3,:)];%translation vector of from femoral origin to tibial origin (in global coordinate frame)
    if right
        lateral=dot(pH_,e1_);%projected onto the medial lateral axis e1
    else
        lateral=dot(pH_,-e1_);
    end
    anterior=dot(pH_,pe2_);%projected onto the anterior posterior axis e2
    superior=-dot(pH_,pe3_);%projected onto the compression distraction axis e3, minus sign to make distraction +ve

    result.flexion = angles(:,1);
    result.medial_rotation = varus'; % Equivalent to varus
    result.lateral_tilt = external'; % Equivalent to external
    result.medial_shift = lateral'; % Equivalent to lateral
    result.anterior = anterior';
    result.superior = superior';
end
