function [r_s,edges_s,edge_type,myosin,P] = initial_myosin(r_s,T_s,edges_s,edge_type,P)
%     find the triangle formed by spectrin bundles in the bottom left of
%     the network
    aux = find(r_s(P.actin,1) == min(r_s(P.actin,1)) & r_s(P.actin,2) == min(r_s(P.actin,2)));
    aux_T = [find(T_s(:,1) == P.actin(aux));find(T_s(:,2) == P.actin(aux));find(T_s(:,3) == P.actin(aux))];
    aux_T = aux_T(1);
    P.myosin_T = aux_T;
%     add a node in the middle of that tringle 
    r_s = [r_s;sum(r_s(T_s(aux_T,:),:))./3];
%    connect that node to the closest stress fiber node and add the
%    corresponding edge
    aux2 = find(r_s(P.stress,1) == min(r_s(P.stress,1)) & r_s(P.stress,2) == min(r_s(P.stress,2)));
    edges_s = [edges_s;P.stress(aux2) size(r_s,1)];
    % % % % % % do the same for the other corners of the mesh
    aux = find(r_s(P.actin,1) == min(r_s(P.actin,1)) & r_s(P.actin,2) == max(r_s(P.actin,2)));
    aux_T = [find(T_s(:,1) == P.actin(aux));find(T_s(:,2) == P.actin(aux));find(T_s(:,3) == P.actin(aux))];
    aux_T = aux_T(1);
    P.myosin_T = [P.myosin_T;aux_T];
    r_s = [r_s;sum(r_s(T_s(aux_T,:),:))./3];
    aux2 = find(r_s(P.stress,1) == min(r_s(P.stress,1)) & r_s(P.stress,2) == max(r_s(P.stress,2)));
    edges_s = [edges_s;P.stress(aux2) size(r_s,1)];
    % % % % % % 
    aux = find(r_s(P.actin,1) == max(r_s(P.actin,1)) & r_s(P.actin,2) == min(r_s(P.actin,2)));
    aux_T = [find(T_s(:,1) == P.actin(aux));find(T_s(:,2) == P.actin(aux));find(T_s(:,3) == P.actin(aux))];
    aux_T = aux_T(1);
    P.myosin_T = [P.myosin_T;aux_T];
    r_s = [r_s;sum(r_s(T_s(aux_T,:),:))./3];
    aux2 = find(r_s(P.stress,1) == max(r_s(P.stress,1)) & r_s(P.stress,2) == min(r_s(P.stress,2)));
    edges_s = [edges_s;P.stress(aux2) size(r_s,1)];
    % % % % % 
    aux = find(r_s(P.actin,1) == max(r_s(P.actin,1)) & r_s(P.actin,2) == max(r_s(P.actin,2)));
    aux_T = [find(T_s(:,1) == P.actin(aux));find(T_s(:,2) == P.actin(aux));find(T_s(:,3) == P.actin(aux))];
    aux_T = aux_T(2);
    P.myosin_T = [P.myosin_T;aux_T];
    r_s = [r_s;sum(r_s(T_s(aux_T,:),:))./3];
    aux2 = find(r_s(P.stress,1) == max(r_s(P.stress,1)) & r_s(P.stress,2) == max(r_s(P.stress,2)));
    edges_s = [edges_s;P.stress(aux2) size(r_s,1)];

%     set type 3 for the myosin linkers
    edge_type = [edge_type;3*ones(size(P.myosin_T))];
    myosin = ((size(r_s,1)-length(P.myosin_T)+1):size(r_s,1))';
end