clear all
close all

% Sup Fig 7, show how spectrin recovers - shrink


% setting up the initial configuration 
P.y0 = 180*1.2;
P.x0 = sqrt(P.y0^2 - (P.y0/2)^2);
P.max_y0 = 17*P.y0/2;
P.max_x0 = P.x0*15;

[r_s,T_s,edges_s,edge_type] = create_mesh(P);
P.actin = 1:size(r_s,1);

P.k0 = 1;%spectrin spring constant
P.d00 = 180;%spectrin resting length

P.delta_t = 0.0001*20;
P.zeta = 0.025*50;%drag

P.t_ini=0;
P.t_end = 60;
aux_t = P.t_ini:P.delta_t:P.t_end;

save_aux = 10;
edges_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);


ll=1;
r_save{ll,1} = r_s;

for k=2:length(aux_t)
    f = zeros(size(r_s));

    for l = 1:size(edges_s,1) 
        r_ij = r_s(edges_s(l,1),:) - r_s(edges_s(l,2),:);
        d = sqrt(dot(r_ij,r_ij,2));
        if edge_type(l) == 0
            f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.k0*r_ij*(d-P.d00)/d;
            f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.k0*r_ij*(d-P.d00)/d;
        end
    end
    r_s(P.actin,:) = P.delta_t*f(P.actin,:)./P.zeta + r_s(P.actin,:);
    
    if mod(k,save_aux) == 1
        ll=ll+1;
        r_save{ll,1} = r_s;

    end
end


