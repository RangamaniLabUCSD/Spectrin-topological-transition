clear all
close all
tic
% Fig 8F, spectrin detachment on


% setting up the initial configuration 
P.y0 = 180;
P.x0 = sqrt(P.y0^2 - (P.y0/2)^2);
P.max_y0 = 17*P.y0/2;
P.max_x0 = P.x0*15;

[r_s,T_s,edges_s,edge_type] = create_mesh(P);
P.actin = 1:size(r_s,1);

aux = find(r_s(:,1) == min(r_s(:,1)) & r_s(:,2) <= max(r_s(:,2)));
r_s_aux = r_s(aux,:) - [P.y0*2 0];
edge_type = [edge_type;5*ones(size(aux))];
edges_s = [edges_s;aux size(r_s,1)+(1:length(aux))'];
r_s = [r_s;r_s_aux];


aux = find(r_s(:,1) == max(r_s(:,1)) & r_s(:,2) <= max(r_s(:,2)));
r_s_aux = r_s(aux,:) + [P.y0*2 0];
edge_type = [edge_type;5*ones(size(aux))];
edges_s = [edges_s;aux size(r_s,1)+(1:length(aux))'];
r_s = [r_s;r_s_aux];

P.adhesion = [find(r_s(:,1) == min(r_s(:,1)));
    find(r_s(:,1) == max(r_s(:,1)))];

P.k0 = 1;%spectrin spring constant
P.d00 = P.y0;%spectrin resting length

P.gamma5 = 150;

P.th = 5e-2;%force threshold

P.delta_t = 0.0001*20;
P.zeta = 0.025*50;%drag

P.t_ini=0;
P.t_end = 120;
aux_t = P.t_ini:P.delta_t:P.t_end;

figure;
trimesh(T_s,r_s(:,1),r_s(:,2),[],'edgecolor',[0.75 0.75 0.75])
hold on 
aux_e = find(edge_type == 0);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'b')
end
aux_e = find(edge_type == 5);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'r')
end
axis('equal')
view(2)


save_aux = 10;
r_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
edges_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
type_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);

ll=1;
edges_save{ll,1} = edges_s;
r_save{ll,1} = r_s;
type_save{ll,1} = edge_type;

for k=2:length(aux_t)
    f = zeros(size(r_s));
    edges_s1 = [];
    
    for l = 1:size(edges_s,1) 
        r_ij = r_s(edges_s(l,1),:) - r_s(edges_s(l,2),:);
        d = sqrt(dot(r_ij,r_ij,2));
        if edge_type(l) == 0
            if P.k0*(d-P.d00)/d > -P.th 
                f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.k0*r_ij*(d-P.d00)/d;
                f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.k0*r_ij*(d-P.d00)/d;
            else
                edges_s1 = [edges_s1;l];
            end
        elseif edge_type(l) == 5
           f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.gamma5.*r_ij;
           f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.gamma5.*r_ij;
        end
    end
    f(P.adhesion,:) = 0;
    r_s = P.delta_t*f./P.zeta + r_s;
    
    edges_s(edges_s1,:) = [];
    edge_type(edges_s1,:) = [];
 
    if mod(k,save_aux) == 1
        ll=ll+1;
        edges_save{ll,1} = edges_s;
        r_save{ll,1} = r_s;
        type_save{ll,1} = edge_type;
    end
end

cols = lines(4);

r_s = r_save{1,1};
edges_s = edges_save{1,1};
edge_type = type_save{1,1};
aux2 = r_s(edges_s(edge_type == 0,1),:) - r_s(edges_s(edge_type == 0,2),:);
s_length_ini = sqrt(dot(aux2,aux2,2));
h = figure;
set(gcf, 'Position',  [50, 50, 600, 300])
axes('Position', [0.1 0.15 0.85 0.8]);
hold on 
col_map = turbo(100);
aux_e = find(edge_type == 0);
r_ij = r_s(edges_s(aux_e,1),:) - r_s(edges_s(aux_e,2),:);
d = sqrt(dot(r_ij,r_ij,2));  
ff = - P.k0.*(d-P.d00)./d;
ff_vec = linspace(-.17,0.17);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color',col_map(find(ff_vec>=ff(l),1),:),'linewidth',1.5)
end
plot(r_s(edges_s(aux_e,:),1),r_s(edges_s(aux_e,:),2),'o','markerfacecolor',cols(4,:),'MarkerEdgeColor',cols(4,:))
aux_e = find(edge_type == 5);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color','k','linewidth',1.5)
end
plot(r_s(P.adhesion,1),r_s(P.adhesion,2),'ko','markerfacecolor','k')
set(gca,'FontSize',14)
xlabel('[nm]')
ylabel('[nm]')
axis('equal')
colormap('turbo')
c = colorbar;
set(c.XLabel,{'String','Rotation','Position'},{'Force [pN]',0,[0.5 max(ff_vec)*1.15]})
caxis([min(ff_vec) max(ff_vec)])
ax_vec = [min(r_s(:,1))-100 max(r_s(:,1))+100 min(r_s(:,2))-100 max(r_s(:,2))+100];
axis(ax_vec)


%%%%%
edges_s = edges_save{end,1};
edge_type = type_save{end,1};
r_s = r_save{end,1};
aux2 = r_s(edges_s(edge_type == 0,1),:) - r_s(edges_s(edge_type == 0,2),:);
s_length_end = sqrt(dot(aux2,aux2,2));
h = figure;
set(gcf, 'Position',  [50, 50, 600, 300])
axes('Position', [0.1 0.15 0.85 0.8]);
hold on 
col_map = turbo(100);
aux_e = find(edge_type == 0);
r_ij = r_s(edges_s(aux_e,1),:) - r_s(edges_s(aux_e,2),:);
d = sqrt(dot(r_ij,r_ij,2)); 
ff = - P.k0.*(d-P.d00)./d;
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color',col_map(find(ff_vec>=ff(l),1),:),'linewidth',1.5)
end
plot(r_s(edges_s(aux_e,:),1),r_s(edges_s(aux_e,:),2),'o','markerfacecolor',cols(4,:),'MarkerEdgeColor',cols(4,:))
aux_e = find(edge_type == 5);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color','k','linewidth',1.5)
end

plot(r_s(P.adhesion,1),r_s(P.adhesion,2),'ko','markerfacecolor','k')
set(gca,'FontSize',14)
xlabel('[nm]')
ylabel('[nm]')
axis('equal')
colormap('turbo')
c = colorbar;
set(c.XLabel,{'String','Rotation','Position'},{'Force [pN]',0,[0.5 max(ff_vec)*1.15]})
caxis([min(ff_vec) max(ff_vec)])
axis(ax_vec)

figure;
set(gcf, 'Position',  [50, 50, 300, 300])
axes('Position', [0.17 0.2 0.8 0.75]);
hold on 
histogram(s_length_ini)
histogram(s_length_end)
set(gca,'FontSize',14)
xlabel('[nm]')
ylabel('Count')
legend('initial','final')