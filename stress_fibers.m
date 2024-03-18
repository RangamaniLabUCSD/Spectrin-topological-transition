function [r_s,edges_s,edge_type,P] = stress_fibers(r_s,edges_s,edge_type,P)
    % create nodes for the stress fibers at the bottom
    aux = find(r_s(:,2) == min(r_s(:,2)));
    aux = aux(1:2:end);
    r_s_aux = r_s(aux,:) - [0 1*P.y0/2];

    r_s_aux2 = r_s(r_s(:,2) ==  P.y0/2,:) - [0 2*P.y0/2];
    r_s_aux3 = sort([r_s_aux;r_s_aux2]);
    aux_stress = size(r_s,1)+(1:size(r_s_aux3,1))';
    P.stress =aux_stress;
%     create edges connecting the stress fibers
    edges_s = [edges_s;aux_stress(1:end-1) aux_stress(2:end)];
    edge_type = [edge_type;2*ones(size(r_s_aux3,1)-1,1)];
    r_s = [r_s;r_s_aux3];
    
    % create nodes for the stress fibers at the top
    r_s_aux = r_s(r_s(:,2) == max(r_s(:,2)),:) + [0 P.y0/2];

    r_s_aux2 = r_s(r_s(:,2) == max(r_s(:,2))-P.y0/2,:) + [0 2*P.y0/2];
    r_s_aux3 = sort([r_s_aux;r_s_aux2]);
    aux_stress = size(r_s,1)+(1:size(r_s_aux3,1))';
    P.stress =[P.stress;aux_stress];
    %     create edges connecting the stress fibers
    edges_s = [edges_s;aux_stress(1:end-1) aux_stress(2:end)];
    edge_type = [edge_type;2*ones(size(r_s_aux3,1)-1,1)];
    r_s = [r_s;r_s_aux3];

    % create nodes corresponding to focal adhesion
    auxx = round(P.max_x0/P.x0);
    auxy =  ceil(max(r_s(:,2))/P.y0);
    r_s_aux = [min(r_s(:,1))-10*P.x0/2 min(r_s(:,2));min(r_s(:,1))-10*P.x0/2 max(r_s(:,2));
        max(r_s(:,1))+10*P.x0/2 min(r_s(:,2));max(r_s(:,1))+10*P.x0/2 max(r_s(:,2))];
     % create edges connecting focal adhesion with stress fibers
    edges_s = [edges_s;auxx*auxy+1 size(r_s,1)+1;auxx*(auxy+1) size(r_s,1)+3;
        auxx*(auxy+1)+1 size(r_s,1)+2;auxx*(auxy+2) size(r_s,1)+4];
    edge_type = [edge_type;2;2;2;2];
    P.adhesion = size(r_s,1)+(1:4)';
    r_s = [r_s;r_s_aux];
    P.actin = setdiff(unique(edges_s),[P.stress;P.adhesion]);

end
