close all 
clear all 


cols = lines(4);

r_s = r_save{1,1};
ax_vec =[87.0614872174388,2905.92230826158,-100,2044];%[min(r_s(:,1))-100 max(r_s(:,1))+100 min(r_s(:,2))-100 max(r_s(:,2))+100];% -300+[87.0614872174388,2905.92230826158,-100,2044];-300+[87.0614872174388,2905.92230826158,-100,2044];%
aux2 = r_s(edges_s(edge_type == 0,1),:) - r_s(edges_s(edge_type == 0,2),:);
s_length_ini = sqrt(dot(aux2,aux2,2));
figure;
set(gcf, 'Position',  [50, 50, 600, 300])
axes('Position', [0.1 0.15 0.85 0.8]);
hold on 
col_map = turbo(100);
aux_e = find(edge_type == 0);
r_ij = r_s(edges_s(aux_e,1),:) - r_s(edges_s(aux_e,2),:);
d = sqrt(dot(r_ij,r_ij,2));  
ff = - P.k0.*(d-P.d00)./d;
ff_vec = linspace(-.25,0.25);
for l=1:size(aux_e,1)
    plot([r_s(edges_s(aux_e(l),1),1),r_s(edges_s(aux_e(l),2),1)],...
       [r_s(edges_s(aux_e(l),1),2),r_s(edges_s(aux_e(l),2),2)],...
       'color',col_map(find(ff_vec<=ff(l),1,'last'),:),'linewidth',1.5)
end
plot(r_s(edges_s(aux_e,:),1),r_s(edges_s(aux_e,:),2),'o','markerfacecolor',cols(4,:),'MarkerEdgeColor',cols(4,:))
set(gca,'FontSize',14)
xlabel('[nm]')
ylabel('[nm]')
axis('equal')
colormap('turbo')
c = colorbar;
set(c.XLabel,{'String','Rotation','Position'},{'Force [pN]',0,[0.5 max(ff_vec)*1.15]})
caxis([min(ff_vec) max(ff_vec)])
axis(ax_vec)


%%%%%

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
