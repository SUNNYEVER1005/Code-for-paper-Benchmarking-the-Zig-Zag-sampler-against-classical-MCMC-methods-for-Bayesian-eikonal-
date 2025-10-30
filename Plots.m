%% MH
results      = mh_results;   % MH results
mean_error   = mh_mean_error;
map_error    = mh_map_error;
GR           = mh_GR;

%% HMC
results      = hmc_results;
mean_error   = hmc_mean_error;
map_error    = hmc_map_error;
GR           = hmc_GR;

%% NUTS
results      = nuts_results;
mean_error   = nuts_mean_error;
map_error    = nuts_map_error;
GR           = nuts_GR;

%% Zig-Zag
results      = zz_results;
mean_error   = zz_mean_error;
map_error    = zz_map_error;
GR           = zz_GR;

%%  (a)
% Velocity-model plot with station and hypocentre locations
% Low-velocity in red, high-velocity in blue
fontsizenum = 14;

% Subplot 1: geological velocity map with source & receivers
figure;
imagesc(x_interval, z_interval, v.');
set(gca, 'YDir', 'reverse');          % flip Y axis
hold on;

jet_map = colormap('jet');
colormap(flipud(jet_map));            % reverse jet so red = slow

% hide tick marks
set(gca, 'TickLength', [0 0]);

size0 = 5;
% receivers – black downward triangles
plot(receiver(:,1), receiver(:,2), 'v', ...
     'MarkerSize', size0, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

% true hypocentre – black pentagram
plot(source(:,1), source(:,2), 'p', ...
     'MarkerSize', size0, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

hcb = colorbar('southoutside');
hcb.Orientation = 'horizontal';
hcb.FontSize    = fontsizenum;

axis equal; axis tight;

legend({'Receiver', 'Hypocentre'}, ...
       'Location', 'best', 'FontSize', fontsizenum, 'FontName', 'Times New Roman');

ylabel('Depth (km)', 'FontSize', fontsizenum, 'FontName', 'Times New Roman');
xlabel('Distance (km)', 'FontSize', fontsizenum, 'FontName', 'Times New Roman');
ax = gca;
ax.XAxisLocation = 'top';

text(mean(xlim), min(ylim)-0.5, '(a)', ...
     'HorizontalAlignment','center', 'VerticalAlignment','top', ...
     'FontSize',fontsizenum, 'FontName','Times New Roman', 'Units','data');
text(mean(xlim), min(ylim)-1, 'km/s', ...
     'HorizontalAlignment','center', 'VerticalAlignment','top', ...
     'FontSize',fontsizenum, 'FontName','Times New Roman', 'Units','data');

set(gca, 'FontSize', fontsizenum);
hold off;

%% (b) Gelman–Rubin  plots  (threshold = 2)
GR_TOTAL=[mh_GR, hmc_GR, nuts_GR, zz_GR];
[gr_row, gr_col] = size(GR_TOTAL);

figure;
colors = [0.9647 0.3255 0.0784;   % red
          1.0000 0.7333 0.0000;   % yellow
          0.0000 0.6314 0.9451];  % blue

b = bar(GR', 'grouped');
for kk = 1:min(gr_row,3)
    b(kk).FaceColor = colors(kk,:);
end

ax = gca;
ax.XAxisLocation = 'bottom';
ax.FontSize = fontsizenum;

if gr_col == 4
    ax.XTick      = 1:4;
    ax.XTickLabel = {'MH','HMC','NUTS','Zig-Zag'};
end

h = legend({'$x_s$','$z_s$','$\tau$'}, ...
           'Interpreter','latex', 'Location','best', ...
           'FontSize',fontsizenum, 'FontName','Times New Roman');

ylabel('Value of $\hat{R}$', 'Interpreter','latex');
grid off; box on;

text(mean(xlim), min(ylim), '(b)', ...
     'HorizontalAlignment','center', 'VerticalAlignment','top', ...
     'FontSize',fontsizenum, 'FontName','Times New Roman', 'Units','data');

ylim([0 2]);

% font check
fonts = listfonts;
if any(strcmpi(fonts,'Helvetica'))
    ax.YAxis.FontName = 'Helvetica';
else
    ax.YAxis.FontName = get(groot,'DefaultAxesFontName');
end
ax.YAxis.FontSize = fontsizenum;

%% (c) ESS per second plot
mh_vec   = cell2mat(per_mh_ess(:));
hmc_vec  = cell2mat(per_hmc_ess(:));
nuts_vec = cell2mat(per_nuts_ess(:));
zz_vec   = cell2mat(per_zz_ess(:));

x_ess = 1:numel(mh_vec);

figure; hold on; grid on;
plot(x_ess, log(mh_vec),  '-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','MH');
plot(x_ess, log(hmc_vec), '-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','HMC');
plot(x_ess, log(nuts_vec),'-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','NUTS');
plot(x_ess, log(zz_vec), '-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Zig-Zag');

xlabel('Hypocentre index');
ylabel('log(ESS)');
h = legend('Location','best');
set(h, 'FontSize',fontsizenum, 'FontName','Times New Roman');

grid off;
ax = gca;
ax.XAxisLocation = 'top';
box on;
ax.FontSize = fontsizenum;

%% (d) Error norms
% mean error
L2_loc_eps   = sum(sqrt(mean_error(1,:).^2 + mean_error(2,:).^2)) / size(mean_error,2);
L2_time_eps  = sum(sqrt(mean_error(3,:).^2)) / size(mean_error,2);
Linf_loc_eps = max(max(abs(mean_error(1,:)), abs(mean_error(2,:)))); % L∞
Linf_time_eps= max(mean_error(3,:));

% MAP error
L2_loc_map   = sum(sqrt(map_error(1,:).^2 + map_error(2,:).^2)) / size(map_error,2);
L2_time_map  = sum(sqrt(map_error(3,:).^2)) / size(map_error,2);
Linf_loc_map = max(max(abs(map_error(1,:)), abs(map_error(2,:))));
Linf_time_map= max(map_error(3,:));