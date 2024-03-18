function [r_s,T_s,edges_s,edge_type] = create_mesh(P)
%     put nodes in position
    aux_x = P.x0:P.x0:P.max_x0;
    aux_y = [0:P.y0:P.max_y0 P.max_y0+P.y0/2];
    aux_y2 = [0 P.y0/2:P.y0:P.max_y0+P.y0/2];
    x_s = [];
    y_s = [];
    for l=1:length(aux_x)
        if mod(l,2) == 1
            y_s = [y_s;aux_y'];
            x_s = [x_s;aux_x(l)*ones(size(aux_y,2),1)];
        else
            y_s = [y_s;aux_y2'];
            x_s = [x_s;aux_x(l)*ones(size(aux_y2,2),1)];
        end
    end
    x_s = x_s(:);
    y_s = y_s(:);
    r_s = [x_s y_s];
%     create mesh
    T_s_delaunay = delaunayTriangulation(x_s,y_s);
    T_s = T_s_delaunay.ConnectivityList;%tringulation
%     get edges
    edges_s = edges(T_s_delaunay);%edges list
    edge_type = zeros(size(edges_s,1),1);%type of edges
% take away connecting edges that have different length than P.y0
    aux_edges = [];
    for l = 1:size(edges_s,1)
        if (y_s(edges_s(l,1)) == min(y_s) && y_s(edges_s(l,2)) == P.y0/2) &&...
            (x_s(edges_s(l,1)) == x_s(edges_s(l,2)))
            aux_edges = [aux_edges;l];
        elseif (y_s(edges_s(l,2)) == min(y_s) && y_s(edges_s(l,1)) == P.y0/2) &&...
            (x_s(edges_s(l,1)) == x_s(edges_s(l,2)))
            aux_edges = [aux_edges;l];
        elseif (y_s(edges_s(l,1)) == max(y_s) && y_s(edges_s(l,2)) == max(y_s)-P.y0/2) &&...
            (x_s(edges_s(l,1)) == x_s(edges_s(l,2)))
            aux_edges = [aux_edges;l];
        elseif (y_s(edges_s(l,2)) == max(y_s) && y_s(edges_s(l,1)) == max(y_s)-P.y0/2) &&...
            (x_s(edges_s(l,1)) == x_s(edges_s(l,2)))
            aux_edges = [aux_edges;l];
        elseif y_s(edges_s(l,1)) == min(y_s) && y_s(edges_s(l,2)) == min(y_s)
            aux_edges = [aux_edges;l];
        elseif y_s(edges_s(l,1)) == max(y_s) && y_s(edges_s(l,2)) == max(y_s)
            aux_edges = [aux_edges;l];
        end
    end

    edges_s(aux_edges,:) = [];
    edge_type(aux_edges) = [];
end
    