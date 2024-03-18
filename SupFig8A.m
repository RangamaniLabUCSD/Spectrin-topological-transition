clear all
close all
tic
rng(040522)

% Supplementary Figure 8A - Single myosin


% setting up the initial configuration 
P.y0 = 180;
P.x0 = sqrt(P.y0^2 - (P.y0/2)^2);
P.max_y0 = 17*P.y0/2;
P.max_x0 = P.x0*15;

[r_s,T_s,edges_s,edge_type] = create_mesh(P);
[r_s,edges_s,edge_type,P] = stress_fibers(r_s,edges_s,edge_type,P);
P = check_free_T(T_s,edges_s,edge_type,P);

P.d_inter = P.y0/3;

P.max_r = 5*P.y0/2;
P.min_r = 3*P.y0/4;
[r_s,edges_s,edge_type,myosin,P] = initial_myosin(r_s,T_s,edges_s,edge_type,P);



P.max_r2 = 5*P.y0/2;
P.min_r2 = 3*P.y0/4;
P.myosin_T2 = [];
P.myosin2_ini = 10;

myosin2 = [];
while size(myosin2,1) < P.myosin2_ini
    [r_s,edges_s,edge_type,myosin2,P] = create_myosin2(r_s,T_s,edges_s,edge_type,myosin2,P);
end



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
aux_e = find(edge_type == 3);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'m')
end
aux_e = find(edge_type == 4);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],'r')
end
axis('equal')
view(2)
plot(r_s(P.adhesion,1),r_s(P.adhesion,2),'ko')
plot(r_s(P.actin,1),r_s(P.actin,2),'ms')
plot(r_s(P.stress,1),r_s(P.stress,2),'c*')
plot(r_s(myosin,1),r_s(myosin,2),'m^')
plot(r_s(myosin2,1),r_s(myosin2,2),'b^')


P.k0 = 1;%spectrin spring constant
P.d00 = P.y0;%spectrin resting length

P.k2 = 4;%stress fiber spring constant
P.gamma2 = 40/P.y0;%stress fiber cable constant
P.d02 = P.x0;%P.y0*P.k2/(P.k2+P.gamma2);%P.x0 + P.gamma2/P.k2;%stress fiber resting length

P.gamma3 = 0.25;%0.25;%75/320;%myosin cable constant
P.gamma4 = P.gamma3*.3/.7;%myosin2 cable constant

P.psi2_a = 1/100;%myosin2 addition rate
P.psi2_r = 1/160;%myosin2 removal rate

P.delta_t = 0.0001*20;
P.zeta = 0.025*50;%drag
P.zeta_a = 3*P.zeta/4;%3*P.zeta/4;%adhesion drag

P.t_ini=0;
P.t_end = 600;
P.t_stress = 300;
aux_t = P.t_ini:P.delta_t:P.t_end;

P.th = 5e-2;%1.01e-1;force threshold

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
        elseif edge_type(l) == 3 
            f(edges_s(l,1),:) =  f(edges_s(l,1),:) - P.gamma3.*r_ij;
            aux_a = T_s(P.myosin_T(aux_m),:);
            f(aux_a,:) =  f(aux_a,:) + P.gamma3*r_ij/3;
            aux_m = aux_m + 1;
        elseif edge_type(l) == 4
            aux_a = T_s(P.myosin_T2(2*aux_m2-1),:);
            f(aux_a,:) =  f(aux_a,:) - P.gamma4*r_ij/3;
            aux_b = T_s(P.myosin_T2(2*aux_m2),:);
            f(aux_b,:) =  f(aux_b,:) + P.gamma4.*r_ij/3;
            aux_m2 = aux_m2 + 1;
        end
    end
    r_s(P.actin,:) = P.delta_t*f(P.actin,:)./P.zeta + r_s(P.actin,:);
    r_s(P.stress,:) = P.delta_t*f(P.stress,:)./P.zeta + r_s(P.stress,:);
    if aux_t(k) < P.t_stress
        r_s(P.adhesion,1) = r_s(P.adhesion,1);
        r_s(P.adhesion([1 3]),2) = r_s(P.adhesion([1 3]),2) + P.delta_t/P.zeta_a;
        r_s(P.adhesion([2 4]),2) = r_s(P.adhesion([2 4]),2) - P.delta_t/P.zeta_a;
    else
        r_s(P.adhesion,:) = r_s(P.adhesion,:);    
    end
    
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
            [r_s,myosin,edges_s2,P] = check_myosin(T_s,r_s,edges_s,edge_type,myosin,edges_s1(l),edges_s2,P);
            [r_s,myosin2,edges_s3,P] = check_myosin2(T_s,r_s,edges_s,edge_type,myosin2,edges_s1(l),edges_s3,P);
        end
    end
    
    
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
    
    
    if mod(k,save_aux) == 1
        ll=ll+1;
        edges_save{ll,1} = edges_s;
        type_save{ll,1} = edge_type;
        mr_save{ll,1} = myosin;
        mr2_save{ll,1} = myosin2;
        r_save{ll,1} = r_s;
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
cols = lines(4);

edges_s = edges_save{1,1};
edge_type = type_save{1,1};
r_s = r_save{1,1};
myosin = mr_save{1,1};
myosin2 = mr2_save{1,1};
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
ax_vec = [min(r_s(:,1))-100 max(r_s(:,1))+100 min(r_s(:,2))-100 max(r_s(:,2))+100];
axis(ax_vec)


%%%%%
edges_s = edges_save{end,1};
edge_type = type_save{end,1};
r_s = r_save{end,1};
myosin = mr_save{end,1};
myosin2 = mr2_save{end,1};
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
% ff = - P.k0.*(d-P.d00)./d;
% ff_vec = linspace(-.17,0.17);
ff = (- P.k0.*(d-P.d00)./d)./s_length_end;
ff_vec = linspace(-.001,0.001);
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
% set(c.XLabel,{'String','Rotation','Position'},{'Force [pN]',0,[0.5 max(ff_vec)*1.15]})
% caxis([min(ff_vec) max(ff_vec)])
% ax_vec = [min(r_s(:,1))-100 max(r_s(:,1))+100 min(r_s(:,2))-100 max(r_s(:,2))+100];
% axis(ax_vec)
set(c.XLabel,{'String','Rotation','Position'},{'Tension [pN/nm]',0,[0.5 max(ff_vec)*1.15]})
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


total_edge = zeros(size(edges_save));
myosin2_number = zeros(size(edges_save));
myosin2_time = [];
myosin2_lifetime = [];
myosin_number = zeros(size(edges_save));
myosin_time = [];
myosin_lifetime = [];


for l = 1:size(edges_save,1)
    myosin2_number(l,1) = length(mr2_save{l,1});
    myosin_number(l,1) = length(mr_save{l,1});
    aux = type_save{l,1};
    total_edge(l) = sum(aux==0);
    
    
    if l == 1
        myosin2_time = ones(size(mr2_save{l,1}));
        myosin_time = ones(size(mr_save{l,1}));
    else
        if ~isempty(setdiff(mr2_save{l,1},mr2_save{l-1,1}))
            myosin2_time = [myosin2_time;
                l*ones(size(setdiff(mr2_save{l,1},mr2_save{l-1,1})))];
%             disp(l)
        end
        if ~isempty(setdiff(mr2_save{l-1,1},mr2_save{l,1}))
            aux = setdiff(mr2_save{l-1,1},mr2_save{l,1});
            aux_ind = zeros(size(aux));
            for k = 1:length(aux)
                aux_ind(k) = find(mr2_save{l-1,1} == aux(k));
            end
            myosin2_lifetime = [myosin2_lifetime;myosin2_time(aux_ind(1:2:end)) l*ones(size(aux_ind(1:2:end)))];
            myosin2_time(aux_ind) = [];
        end
        if ~isempty(setdiff(mr_save{l,1},mr_save{l-1,1}))
            myosin_time = [myosin_time;
                l*ones(size(setdiff(mr_save{l,1},mr_save{l-1,1})))];
          
        end
        if ~isempty(setdiff(mr_save{l-1,1},mr_save{l,1}))
            aux = setdiff(mr_save{l-1,1},mr_save{l,1});
            aux_ind = zeros(size(aux));
            for k = 1:length(aux)
                aux_ind(k) = find(mr_save{l-1,1} == aux(k));
            end
            myosin_lifetime = [myosin_lifetime;myosin_time(aux_ind) l*ones(size(aux_ind))];
            myosin_time(aux_ind) = [];
        end
        
    end
        
end

h = figure;
set(gcf, 'Position',  [50, 50, 500, 200])
axes('Position', [0.1 0.2 0.85 0.75]);
set(gca,'FontSize',14)
hold on
plot(aux_t(1:save_aux:end),myosin2_number./2,'linewidth',1.5)
xlabel('time [s]')
ylabel('# Myosin rods')



figure;
set(gcf, 'Position',  [50, 50, 250, 200])
axes('Position', [0.2 0.2 0.75 0.75]);
x1 = (aux_t((myosin2_lifetime(:,2)-1)*save_aux+1)-aux_t((myosin2_lifetime(:,1)-1)*save_aux+1))';
g1 = repmat({'Myosin'},size(x1,1),1);
hold on
h =boxplot(x1,g1,'Widths',0.5);
set(h,{'linew'},{2})
set(gca,'FontSize',14)
ylabel('Myosin rod lifetime [s]')


neighbors = zeros(1,length(P.actin));
l = size(edges_save,1);
aux = type_save{l,1};
for ll = 1:length(P.actin)
    aux1 = [find(edges_save{l,1}(:,1) == P.actin(ll));
        find(edges_save{l,1}(:,2) == P.actin(ll))];
    if ~isempty(aux1)
        aux1 = unique(aux1);
        neighbors(ll) = length(aux1);
    end  
end


figure
col_neighbor = turbo(10);
col_neighbor = [col_neighbor(7:10,:);col_neighbor(1:3,:)];
r_s = r_save{end,1};
r_s_a = r_s(P.actin,:);
figure;
set(gcf, 'Position',  [50, 50, 600, 300])
axes('Position', [0.1 0.15 0.85 0.8]);
hold on 
for l = 0:6
    scatter(r_s_a(neighbors==l,1),r_s_a(neighbors==l,2),50,col_neighbor(l+1,:),'filled')
end
alpha 0.5
legend('0-neighbor','1-neighbor','2-neighbors','3-neighbors','4-neighbors','5-neighbors','6-neighbors')
set(gca,'FontSize',14)
xlabel('[nm]')
ylabel('[nm]')
axis('equal')
colorbar
