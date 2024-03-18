clear
close all
rng(040522)


% Supplementary Figure 8L, faster focal adhesions

P.y0 = 180;%nm, corresponds with spectrin bundle length
P.x0 = sqrt(P.y0^2 - (P.y0/2)^2);%nm; height of the triangles
P.max_y0 = 17*P.y0/2;
P.max_x0 = P.x0*15;

[r_s,T_s,edges_s,edge_type] = create_mesh(P);
% r_s:position of the nodes, edges_s: edge connectivity matrix, edge_type :
% type of edge
[r_s,edges_s,edge_type,P] = stress_fibers(r_s,edges_s,edge_type,P);
P = check_free_T(T_s,edges_s,edge_type,P);

P.max_r = 5*P.y0/2;%nm, maximum lenght of myosins linkers
P.min_r = 3*P.y0/4;%nm, minimum lenght of myosins linkers
% create myosin linkers
[r_s,edges_s,edge_type,myosin,P] = initial_myosin(r_s,T_s,edges_s,edge_type,P);
% myosin: entries of the positin vector r_s corresponding to
% myosin linkers nodes

P.d_inter = P.y0/3;%auxiliar variable to check myosin intersections
P.max_r2 = 5*P.y0/2;%nm, maximum lenght of myosins rods
P.min_r2 = 3*P.y0/4;%nm, minimum lenght of myosins rods
P.myosin_T2 = []; %P.myosin_T2: entries of the triangulation T_s corresponding to the triangles that have myosin rods
P.myosin2_ini = 10;%number of initial myosin rodsX2

% create myosin rods
myosin2 = [];
while size(myosin2,1) < P.myosin2_ini
    [r_s,edges_s,edge_type,myosin2,P] = create_myosin2(r_s,T_s,edges_s,edge_type,myosin2,P);
end
% myosin2: entries of the positin vector r_s corresponding to
% myosin rods nodes

% create a plot to check if the mesh is correct
figure;
trimesh(T_s,r_s(:,1),r_s(:,2),[],'edgecolor',[0.75 0.75 0.75])
hold on 
% plot spectrin bundles
aux_e = find(edge_type == 0);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'b')
end
% plot stress fibers
aux_e = find(edge_type == 2);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'g')
end
% plot myosin linkers
aux_e = find(edge_type == 3);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'m')
end
% plot myosin rods
aux_e = find(edge_type == 4);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'r')
end
axis('equal')
view(2)
plot(r_s(P.adhesion,1),r_s(P.adhesion,2),'ko')%plot focal adhesions
plot(r_s(P.actin,1),r_s(P.actin,2),'ms')%plot actin short filaments 
plot(r_s(P.stress,1),r_s(P.stress,2),'c*')
plot(r_s(myosin,1),r_s(myosin,2),'m^')
plot(r_s(myosin2,1),r_s(myosin2,2),'b^')


P.k0 = 1;%pN/nm,spectrin spring constant
P.d00 = P.y0;%nm, spectrin resting length

P.k2 = 4;%pN/nm,stress fiber spring constant
P.gamma2 = 40/P.y0;%pN/nm,stress fiber cable constant
P.d02 = P.x0;%nm,stress fiber resting length

P.gamma3 = 0.25;%pN/nm,myosin cable constant
P.gamma4 = P.gamma3*.3/.7;%pN/nm,myosin2 cable constant

P.psi2_a = 1/100;%1/s,myosin2 addition rate
P.psi2_r = 1/160;%1/s,myosin2 removal rate

P.delta_t = 0.0001*20;%s, time step size
P.zeta = 0.025*50;%pN s/nm, drag
P.zeta_a = 5*P.zeta/16;%3*P.zeta/4;%pN s/nm,adhesion drag

P.t_ini=0;%s, initial time
P.t_end = 600;%s,end time
P.t_stress = 300;%s, time to stop moving stress fibers
aux_t = P.t_ini:P.delta_t:P.t_end;

P.th = 5e-2;%pN, force threshold

save_aux = 10;
edges_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
type_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
mr_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
mr2_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
r_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);


ll=1;
edges_save{ll,1} = edges_s;
type_save{ll,1} = edge_type;
mr_save{ll,1} = myosin;
mr2_save{ll,1} = myosin2;
r_save{ll,1} = r_s;

for k=2:length(aux_t)
    aux_m = 1;
    aux_m2 = 1;
    f = zeros(size(r_s));
    edges_s1 = [];
    edges_s2 = [];
    edges_s3 = [];
    
    for l = 1:size(edges_s,1) 
%         calculate the distance of edges
        r_ij = r_s(edges_s(l,1),:) - r_s(edges_s(l,2),:);
        d = sqrt(dot(r_ij,r_ij,2));
%         depending of the type of edge compute the force and store the
%         value
        if edge_type(l) == 0
            if P.k0*(d-P.d00)/d > -P.th 
                f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.k0*r_ij*(d-P.d00)/d;
                f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.k0*r_ij*(d-P.d00)/d;
            else
%                 save the spectrin bundles that get detach
                edges_s1 = [edges_s1;l];
            end  
        elseif edge_type(l) == 2 
            f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.k2*r_ij*(d-P.d02)/d - P.gamma2*r_ij;
            f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.k2*r_ij*(d-P.d02)/d + P.gamma2*r_ij;
        elseif edge_type(l) == 3 
            f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.gamma3.*r_ij;
%             distribute the force among the three nodes of the spectrin
%             triangle that the myosin linker attaches to 
            aux_a = T_s(P.myosin_T(aux_m),:);
            f(aux_a,:) =  f(aux_a,:) + P.gamma3*r_ij/3;
            aux_m = aux_m + 1;
        elseif edge_type(l) == 4
            aux_a = T_s(P.myosin_T2(2*aux_m2-1),:);
%             distribute the force among the three nodes of the pectrin
%             triangles that the myosin rods attaches to 
            f(aux_a,:) =  f(aux_a,:) - P.gamma4*r_ij/3;
            aux_b = T_s(P.myosin_T2(2*aux_m2),:);
            f(aux_b,:) =  f(aux_b,:) + P.gamma4.*r_ij/3;
            aux_m2 = aux_m2 + 1;
        end
    end
%     update the nodes position 
    r_s(P.actin,:) = P.delta_t*f(P.actin,:)./P.zeta + r_s(P.actin,:);
    r_s(P.stress,:) = P.delta_t*f(P.stress,:)./P.zeta + r_s(P.stress,:);
    if aux_t(k) < P.t_stress
        r_s(P.adhesion,1) = r_s(P.adhesion,1);
        r_s(P.adhesion([1 3]),2) = r_s(P.adhesion([1 3]),2) + P.delta_t/P.zeta_a;
        r_s(P.adhesion([2 4]),2) = r_s(P.adhesion([2 4]),2) - P.delta_t/P.zeta_a;
    else
        r_s(P.adhesion,:) = r_s(P.adhesion,:);    
    end
%  for the spectrin edges that were detached, get the corresponding triangles  
    if ~isempty(edges_s1)
        for l = 1:length(edges_s1)
            aux_T = intersect([find(T_s(:,1) == edges_s(edges_s1(l),1));
                    find(T_s(:,2) == edges_s(edges_s1(l),1));
                    find(T_s(:,3) == edges_s(edges_s1(l),1))],...
                    [find(T_s(:,1) == edges_s(edges_s1(l),2));
                    find(T_s(:,2) == edges_s(edges_s1(l),2));
                    find(T_s(:,3) == edges_s(edges_s1(l),2))]);
            P.myosin_Tfree = setdiff(P.myosin_Tfree,aux_T);
        end
        
        for l = 1:length(edges_s1)
%             check if the myosin linker in the spectrin triangle that lost an
%             edge can attach to another spectrin triangle 
            [r_s,myosin,edges_s2,P] = check_myosin(T_s,r_s,edges_s,edge_type,myosin,edges_s1(l),edges_s2,P);
%             check if the myosin rod in the spectrin triangle that lost an
%             edge can attach to another spectrin triangle 
            [r_s,myosin2,edges_s3,P] = check_myosin2(T_s,r_s,edges_s,edge_type,myosin2,edges_s1(l),edges_s3,P);
        end
    end
%     randomly add a myosin rod
    if rand <= P.psi2_a*P.delta_t
       [r_s,edges_s,edge_type,myosin2,P] = create_myosin2(r_s,T_s,edges_s,edge_type,myosin2,P); 
    end
    
%     update position of the myosins
    for l = 1:length(myosin)
        r_s(myosin(l),:) = sum(r_s(T_s(P.myosin_T(l),:),:))./3;
    end
    
    for l = 1:length(myosin2)
        r_s(myosin2(l),:) = sum(r_s(T_s(P.myosin_T2(l),:),:))./3;
    end
    edges_s([edges_s1;edges_s2;edges_s3],:) = [];
    edge_type([edges_s1;edges_s2;edges_s3],:) = [];
    
%     remove small myosin
    [edges_s,edge_type,myosin,P] = remove_small_myosin(r_s,edges_s,edge_type,myosin,P); 
    [edges_s,edge_type,myosin2,P] = remove_small_myosin2(r_s,edges_s,edge_type,myosin2,P);
%     randomly remove a myosin rod
    if rand <= P.psi2_r*P.delta_t
        [edges_s,edge_type,myosin2,P] = remove_myosin2(edges_s,edge_type,myosin2,P);
    end
%     save data 
    if mod(k,save_aux) == 1
        ll=ll+1;
        edges_save{ll,1} = edges_s;
        type_save{ll,1} = edge_type;
        mr_save{ll,1} = myosin;
        mr2_save{ll,1} = myosin2;
        r_save{ll,1} = r_s;
    end
end

k = fix(find(aux_t == P.t_end)/save_aux)+1;
ff_vec = linspace(-.17,0.17);
col_map = turbo(100);
figure
set(gcf, 'Position',  [50, 50, 600, 300])
axes('Position', [0.1 0.15 0.85 0.8]);
hold on 
edges_s = edges_save{k,1};
edge_type = type_save{k,1};
r_s = r_save{k,1};
myosin2 = mr2_save{k,1};
myosin = mr_save{k,1};
aux2 = r_s(edges_s(edge_type == 0,1),:) - r_s(edges_s(edge_type == 0,2),:);
s_length_end = sqrt(dot(aux2,aux2,2));

aux_e = find(edge_type == 0);
r_ij = r_s(edges_s(aux_e,1),:) - r_s(edges_s(aux_e,2),:);
d = sqrt(dot(r_ij,r_ij,2)); 
ff = - P.k0.*(d-P.d00)./d;%- P.k0.*(d-P.d00)./d;
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color',col_map(find(ff_vec>=ff(l),1),:),'linewidth',1.5)
end

aux_e = find(edge_type == 2);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color','k','linewidth',1.5)
end

aux_e = find(edge_type == 3);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color','m','linewidth',1.5)
end

aux_e = find(edge_type == 4);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color','m','linewidth',1.5)
end
plot(r_s(P.adhesion,1),r_s(P.adhesion,2),'ko','markerfacecolor','k')
plot(r_s(P.stress,1),r_s(P.stress,2),'o','MarkerEdgeColor','k')
plot(r_s(myosin,1),r_s(myosin,2),'o','MarkerEdgeColor','m')
plot(r_s(myosin2,1),r_s(myosin2,2),'o','MarkerEdgeColor','m')
set(gca,'FontSize',14)
xlabel('[nm]')
ylabel('[nm]')
axis('equal')
colormap('turbo')
c = colorbar;
set(c.XLabel,{'String','Rotation','Position'},{'Force [pN]',0,[0.5 max(ff_vec)*1.15]})
caxis([min(ff_vec) max(ff_vec)])
title(['time = ' num2str(aux_t((k-1)*save_aux+1),'%.4f') ' s  '])
