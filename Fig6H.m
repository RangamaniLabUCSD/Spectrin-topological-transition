clear all
close all

% Fig 8H, tooked strees fibers
P.y0 = 180;
P.x0 = sqrt(P.y0^2 - (P.y0/2)^2);
P.max_y0 = 17*P.y0/2;
P.max_x0 = P.x0*15;

[r_s,T_s,edges_s,edge_type] = create_mesh(P);

aux = find(r_s(:,2) == min(r_s(:,2)));
P.stress =aux;
edges_s = [edges_s;aux(1:end-1) aux(2:end)];

aux = find(r_s(:,2) == max(r_s(:,2)));
r_s_aux2 = r_s(r_s(:,2) == max(r_s(:,2))-P.y0/2,:) + [0 P.y0/2];
aux_stress = size(r_s,1)+(1:size(r_s_aux2,1))';
P.stress =[P.stress;aux;aux_stress];
edges_s = [edges_s;aux(1:end-1) aux_stress;aux_stress aux(2:end)];
r_s = [r_s;r_s_aux2];

for l = 1:size(edges_s,1)
    if r_s(edges_s(l,1),2) == min(r_s(:,2)) && r_s(edges_s(l,2),2) == min(r_s(:,2))
        edge_type(l) = 2;
    elseif r_s(edges_s(l,1),2) == max(r_s(:,2)) && r_s(edges_s(l,2),2) == max(r_s(:,2))
        edge_type(l) = 2;
    end
end


    % focal adhesion
auxx = round(P.max_x0/P.x0);
auxy =  ceil(max(r_s(:,2))/P.y0);
aux_adh = 5*P.x0;
r_s_aux = [min(r_s(:,1))-aux_adh min(r_s(:,2));min(r_s(:,1))-aux_adh max(r_s(:,2));
    max(r_s(:,1))+aux_adh min(r_s(:,2));max(r_s(:,1))+aux_adh max(r_s(:,2))];
edges_s = [edges_s;1 size(r_s,1)+1;auxx*(auxy+1)-auxy size(r_s,1)+3;
    (auxy+1) size(r_s,1)+2;auxx*(auxy+1) size(r_s,1)+4];
edge_type = [edge_type;2;2;2;2];
P.adhesion = size(r_s,1)+(1:4)';
r_s = [r_s;r_s_aux];
P.actin = setdiff(unique(edges_s),[P.stress;P.adhesion]);

% setting up the initial configuration 

P.d_inter = P.y0/3;

figure;
trimesh(T_s,r_s(:,1),r_s(:,2),[],'edgecolor',[0.75 0.75 0.75])
hold on 
aux_e = find(edge_type == 0);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'b')
end
aux_e = find(edge_type == 2);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'g')
end

axis('equal')
view(2)
plot(r_s(P.adhesion,1),r_s(P.adhesion,2),'ko')
plot(r_s(P.actin,1),r_s(P.actin,2),'ms')
plot(r_s(P.stress,1),r_s(P.stress,2),'c*')

P.k0 = 1;%spectrin spring constant
P.d00 = P.y0;%spectrin resting length

P.k2 = 4;%stress fiber spring constant
P.gamma2 = 40/180;%stress fiber cable constant
P.d02 = P.x0;%P.y0*P.k2/(P.k2+P.gamma2);%P.x0 + P.gamma2/P.k2;%stress fiber resting length

P.delta_t = 0.0001*20;
P.zeta = 0.025*50;%drag
P.zeta_a = 3*P.zeta/4;%3*P.zeta/4;%adhesion drag

P.t_ini=0;
P.t_end = 600;
aux_t = P.t_ini:P.delta_t:P.t_end;

P.th = 5e-2;%1.01e-1;force threshold

save_aux = 10;
edges_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
type_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);
r_save = cell(P.t_end/(save_aux*P.delta_t)+1,1);


ll=1;
edges_save{ll,1} = edges_s;
type_save{ll,1} = edge_type;
r_save{ll,1} = r_s;

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
        elseif edge_type(l) == 2 
            f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.k2*r_ij*(d-P.d02)/d - P.gamma2*r_ij;
            f(edges_s(l,2),:) =  f(edges_s(l,2),:) + P.k2*r_ij*(d-P.d02)/d + P.gamma2*r_ij;
        end
    end
    r_s(P.actin,:) = P.delta_t*f(P.actin,:)./P.zeta + r_s(P.actin,:);
    r_s(P.stress,:) = P.delta_t*f(P.stress,:)./P.zeta + r_s(P.stress,:);
    r_s(P.adhesion,:) = r_s(P.adhesion,:);
   
    edges_s(edges_s1,:) = [];
    edge_type(edges_s1,:) = [];

    if mod(k,save_aux) == 1
        ll=ll+1;
        edges_save{ll,1} = edges_s;
        type_save{ll,1} = edge_type;
        r_save{ll,1} = r_s;

    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % /
cols = lines(4);

edges_s = edges_save{1,1};
edge_type = type_save{1,1};
r_s = r_save{1,1};
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
aux_e = find(edge_type == 2);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color','k','linewidth',1.5)
end
plot(r_s(P.adhesion,1),r_s(P.adhesion,2),'ko','markerfacecolor','k')
plot(r_s(P.stress,1),r_s(P.stress,2),'o','MarkerEdgeColor','k')
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
ff = -P.k0.*(d-P.d00)./d;%- P.k0.*(d-P.d00)./d;
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color',col_map(find(ff_vec>=ff(l),1),:),'linewidth',1.5)
end
plot(r_s(edges_s(aux_e,:),1),r_s(edges_s(aux_e,:),2),'o','markerfacecolor',cols(4,:),'MarkerEdgeColor',cols(4,:))
aux_e = find(edge_type == 2);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color','k','linewidth',1.5)
end
plot(r_s(P.adhesion,1),r_s(P.adhesion,2),'ko','markerfacecolor','k')
plot(r_s(P.stress,1),r_s(P.stress,2),'o','MarkerEdgeColor','k')
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

neighbors = zeros(1,length(P.actin));
l = size(edges_save,1);
for ll = 1:length(P.actin)
    aux1 = [find(edges_save{l,1}(:,1) == P.actin(ll));
        find(edges_save{l,1}(:,2) == P.actin(ll))];
    if ~isempty(aux1)
        aux1 = unique(aux1);
        neighbors(1,ll) = length(aux1);
    end  
end

h = figure;
set(gcf, 'Position',  [50, 50, 300, 300])
axes('Position', [0.17 0.2 0.8 0.75]);
hold on
histogram(neighbors(1,:),'Normalization','probability')
