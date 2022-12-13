%% 绘图
%% 测试例子

X1 = 1:0.1:10;
Y1 = sin(X1);
Y2 = cos(X1);
figure(1)
% 设置图窗属性
set(gcf, 'unit','centimeters', 'position', [5,5,9,6]) % 屏幕上的起始坐标5，5； 大小9*6

linewidth_line = 1.2; %提前设置好线形大小
markersize = 2.5;
linewidth_gca = 7;
linewidth_label = 9;
fontsize_ledgend = 7;
fontsize_label = 7;
% 绘制第一个图， 设置线形 颜色 线宽 节点大小
% - 实线 -- 虚线 : -. 点线 
plot(X1, Y1, '-', 'color', 'r', 'linewidth', linewidth_line, 'markersize', markersize);
plot(X1, Y2, '-.', 'color', 'b', 'linewidth', linewidth_line, 'markersize', markersize);

hold on 
grid on; % 打开网格
% grid off;
grid minor; % 设置网格大小

legend('Y1', 'Y2'); % 设置图例
legend('orientation', 'horizontal') % 将图例排布设为横向

xlim([1,10]);
ylim([-1,1]); % 设置绘图范围

% 设置图像属性
set(gca, 'linewidth', linewidth_line, 'Fontsize', fontsize_label); % 设置线宽， 字体大小
set(gca, 'GridlineStyle', '--');

xlabel('Xlabel');
ylabel('Ylabel');
title('Title');
%% Figure 1. 数据结构示意图

% 使用函数titlelayout创建分区布局图，并调整参数
figure(1)
set(gcf, 'unit', 'centimeters', 'position', [10, 10, 19, 10]); % 中等画幅宽度
% set(gcf, 'fontname', 'Timesnewroman')
tiledlayout(2,2, 'TileSpacing', 'compact');  % 图间距设置为最低，none为无间距   , 'TileSpacing', 'compact'

nexttile
% X1 = turn_to_pointcloud(MAT);
% X1(:,1:2) = X1(:,1:2)/6.4;
% X1(:,3) = X1(:,3)*1e6;
X = linspace(0, 80, 512);
[X, Y] = meshgrid(X);
mesh(X, Y, MAT*1e6 +4);
set(gca, 'xtick', [40, 80], 'ytick', [0,40,80]);
xlabel('x \mum'); ylabel('y \mum'); zlabel('height \mum');
xlim([0, 80]); ylim([0, 80]);
% grid off
% title('(a) Full 3D Image of Depth Figure')

nexttile
X = linspace(0, 8, 51);
[X, Y] = meshgrid(X);
mesh(X, Y, MAT(30:80, 10:60)'*1e6+4);
set(gca, 'xtick', [4, 8], 'ytick', [0,4,8]);
xlabel('x \mum'); ylabel('y \mum'); zlabel('height \mum');
xlim([0,8]); ylim([0,8]);
% grid off

nexttile
X = linspace(0, 8, 51);
[X, Y] = meshgrid(X);
mesh(X, Y, MAT(90:140, 90:140)*1e6+4);
set(gca, 'xtick', [4, 8], 'ytick', [0,4,8], 'ztick', [2, 3]);
xlabel('x \mum'); ylabel('y \mum'); zlabel('height \mum');
xlim([0,8]); ylim([0,8]); zlim([2,3.3])
title('(c)');
% grid off

nexttile
X = linspace(0, 8, 51);
[X, Y] = meshgrid(X);
mesh(X, Y, MAT(226:276,226:276)*1e6+4);
set(gca, 'xtick', [4, 8], 'ytick', [0,4,8], 'ztick', [4, 5]);
xlabel('x \mum'); ylabel('y \mum'); zlabel('height \mum');
xlim([0,8]); ylim([0,8]);
title('(d)');

set(gcf, 'paperunits', 'centimeters');
set(gcf,'PaperSize',[19, 12]);
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperPosition',[0,0,19 12]);
set(gcf,'Renderer','painters');
% 
exportgraphics(gcf,'peaks.tiff','Resolution',600);%输出分辨率为300的PNG图片
% exportgraphics(gcf,'peaks.pdf','ContentType','image');%输出矢量pdf图片2.1
% exportgraphics(gcf,'peaks.eps');%输出矢量eps图片

% print 1.tif -dtiff -r600

% 调整xlabel和ylabel角度，分别为-20 和 10
% 添加平面颜色， 调整字体


% grid off

% (a) mesh figure of the depth image
% (b) two adjacent planes
% (c) two adjacent planes
% (d) Planes and adjacent surface

%% Fig 3  对三种类型数据的标注结果
% 图片竖向排版，全都从设计模型里截取
% 
% 排版顺序：
% 
% 目标平面 倾斜平面 曲面

figure(2)
set(gcf, 'unit', 'centimeters', 'position', [10, 10, 19, 8]); % 中等画幅宽度
% set(gcf, 'fontname', 'Timesnewroman')
tiledlayout(2,3, 'TileSpacing', 'compact');  % 图间距设置为最低，none为无间距   , 'TileSpacing', 'compact'

target = 100*(s.*D);
tiledlayout(2,3);

nexttile
X = linspace(1,5,32);
[X, Y] = meshgrid(X);
% Z1 = MAT(246:277,236:267);
Z1 = MAT(246:277,233:264);
mesh(X,Y,Z1);
set(gca, 'xtick', [0, 5], 'ytick', [0,5], 'ztick', [12.8]);
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel('x/\mum'); ylabel('y/\mum'); zlabel('z/\mum');
set(get(gca, 'xlabel'), 'Rotation', 20); % 'Position',[1,1,1]
set(get(gca, 'ylabel'), 'Rotation', -35);
xlim([0,5]); ylim([0,5]);

nexttile
X = linspace(1,5,32);
[X, Y] = meshgrid(X);
Z3 = MAT(1:32,1:32)/6.4;
mesh(X,Y,Z3);
set(gca, 'xtick', [0, 5], 'ytick', [0,5], 'ztick', [-5.0,-4.5, -4.0]);
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel('x/\mum'); ylabel('y/\mum'); zlabel('z/\mum');
set(get(gca, 'xlabel'), 'Rotation', 20);
set(get(gca, 'ylabel'), 'Rotation', -35);
xlim([0,5]); ylim([0,5]);

nexttile
X = linspace(1,5,32);
[X, Y] = meshgrid(X);
Z7 = MAT(191:222,191:222)/6.4;
mesh(X,Y,Z7);
set(gca, 'xtick', [0, 5], 'ytick', [0,5], 'ztick', [-1,-0.5]);
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel('x/\mum'); ylabel('y/\mum'); zlabel('z/\mum');
set(get(gca, 'xlabel'), 'Rotation', 20);
set(get(gca, 'ylabel'), 'Rotation', -35);
xlim([0,5]); ylim([0,5]);

nexttile
X = linspace(1,5,32);
[X, Y] = meshgrid(X);
T1 = target(241:272,241:272);
mesh(X,Y,T1);
set(gca, 'xtick', [0, 5], 'ytick', [0,5], 'ztick', [0, 0.1]);
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel('x/\mum'); ylabel('y/\mum'); zlabel('Eigenvalues');
set(get(gca, 'xlabel'), 'Rotation', 20);
set(get(gca, 'ylabel'), 'Rotation', -35);
xlim([0,5]); ylim([0,5]);zlim([-0.25, 0.25])

nexttile
X = linspace(1,5,32);
[X, Y] = meshgrid(X);
% T3 = target(2:33,2:33)/6.4e4;
T3 = target(2:33,22:53)*1e3;
mesh(X,Y,T3);
set(gca, 'xtick', [0, 5], 'ytick', [0,5]); % , 'ztick', [1.5, 1.7]
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel('x/\mum'); ylabel('y/\mum'); zlabel('Eigenvalues');
set(get(gca, 'xlabel'), 'Rotation', 20);
set(get(gca, 'ylabel'), 'Rotation', -35);
xlim([0,5]); ylim([0,5]);zlim([4,5.5]);

nexttile
X = linspace(1,5,32);
[X, Y] = meshgrid(X);
T5 = target(191:222,191:222)*0.25e3;
mesh(X,Y,T5);
set(gca, 'xtick', [0, 5], 'ytick', [0,5], 'ztick', [0, 1, 2.5]);
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel('x/\mum'); ylabel('y/\mum'); zlabel('Eigenvalues');
set(get(gca, 'xlabel'), 'Rotation', 20);
set(get(gca, 'ylabel'), 'Rotation', -35);
xlim([0,5]); ylim([0,5]);

set(gcf, 'paperunits', 'centimeters');
set(gcf,'PaperSize',[19, 12]);
set(gcf,'PaperPositionMode','manual','PaperPosition',[0,0,19 12])
set(gcf,'Renderer','painters');
% exportgraphics(gcf,'peaks.tiff','Resolution',600);%输出分辨率为300的PNG图片
exportgraphics(gcf, 'opeaks.pdf', 'ContentType', 'vector'); % 输出pdf图片
% 保存后使用pdf编辑器添加子标题，之后再使用AI导出为tiff格式
%% Figure 4 提取到的特征的分布
% 这个图，每次运行主程序之后，将P数据保存出来，然后放在一起统计画出来，全部画直方图
% 
% 格式：标准模型 数据1 数据2

figure(3)

% tiledlayout(3);
set(gcf, 'unit', 'centimeters', 'position', [10, 10, 9, 14]); % 中等画幅宽度
tiledlayout(2,3, 'TileSpacing', 'compact'); 

subplot(3,1,1);
histogram(P_0(:,4)*5);
% set(gca, 'xtick', [0, 5], 'ytick', [0,5], 'ztick', [0, 1, 2.5]);
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel({'Eigenvalues'; '(a)'}); ylabel('Points Quantity');
xlim([-0.5,4]);


subplot(3,1,2);
histogram(P_2(:,4)*5)
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel({'Eigenvalues'; '(b)'}); ylabel('Points Quantity');
xlim([-0.5,4]);

subplot(3,1,3);
histogram(P_3(:,4)*7);
set(gca, 'Fontname', 'Times New Roman', 'Fontsize', 9);
xlabel({'Eigenvalues'; '(c)'}); ylabel('Points Quantity');
xlim([-0.5,4]);

% exportgraphics(gcf, 'Figure4.tiff', 'Resolution', 600); 
exportgraphics(gcf, 'Figure4.pdf', 'ContentType', 'vector'); 
%% Fig 5 对数据的整体标注和分割提取效果
% 数据格式（3x3）：
% 
% 标准模型 标注图 提取图
% 
% 实测模型 标注图 提取图
% 
% 实测模型 标注图 提取图
% 
% Figure 3 阈值取33画出来的图好看

f5 = figure(5);
set(gcf, 'unit', 'centimeters', 'position', [10,10,14,16]);
tiledlayout(3,3, 'TileSpacing', 'compact'); 

% group 1
nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,MAT_0/6.4+5);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
xlabel({'Length/\mum'}); ylabel('Width/\mum');
view(2);
xlim([0,80]); ylim([0,80]);

nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,sD_0);
set(gca,'ytick', []);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
xlabel({'Length/\mum';'(a)'}); ylabel('');
view(2);
xlim([0,80]); ylim([0,80]);

nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,XY_0);
set(gca,'ytick', []);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
grid off;
xlabel({'Length/\mum'}); ylabel('');
view(2);
xlim([0,80]); ylim([0,80]);

% group 2
nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,MAT_2/6.4+5);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
xlabel({'Length/\mum'}); ylabel({'Width/\mum'});
view(2);
xlim([0,80]); ylim([0,80]);

nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,sD_2);
set(gca,'ytick', []);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
xlabel({'Length/\mum';'(b)'}); ylabel({''});
view(2);
xlim([0,80]); ylim([0,80]);

nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,XY_2);
set(gca,'ytick', []);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
grid off;
xlabel({'Length/\mum'}); ylabel('');
view(2);
xlim([0,80]); ylim([0,80]);


% group 3
nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,MAT_3/6.4+5);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
xlabel({'Length/\mum'}); ylabel('Width/\mum');
view(2);
xlim([0,80]); ylim([0,80]);
nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,sD_3);
set(gca,'ytick', []);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
xlabel({'Length/\mum'; '(c)'}); ylabel('');
view(2);
xlim([0,80]); ylim([0,80]);

nexttile
X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,XY_3);
set(gca,'ytick', []);
set(gca, 'fontname', 'Times New Roman', 'FontSize', 9);
grid off;
xlabel({'Length/\mum'}); ylabel('');
view(2);
xlim([0,80]); ylim([0,80]);

exportgraphics(f5, 'seg_results.tiff', 'Resolution', 600);
% exportgraphics(f5, 'seg_results.pdf', 'ContentType', 'vector');
%% Fig 6 测量对象标注

f6 = figure();
set(gcf, 'unit', 'centimeters', 'position', [10,10,9,6]);

X = linspace(0,80,512);
[X,Y] = meshgrid(X);
mesh(X,Y,MAT/6.4+5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel({'x/\mum'}); ylabel('y/\mum'); zlabel('z/\mum');
ylim([0, 40]);
set(gca, 'xtick', [0, 80], 'ytick', [40,80], 'ztick', [0,2,4]);
set(get(gca, 'xlabel'), 'Rotation', 3);
set(get(gca, 'ylabel'), 'Rotation', -80);

view([175.94 52.67])

% h1 = annotation('line',[0.2675 0.5086],[0.4769 0.4901], 'color', 'r', 'LineStyle', '--', 'LineWidth', 1);
% set(h1,'X',[0.5409 0.538],'Y',[0.5606 0.2259]);
% annotation('line',[0.2675 0.5086],[0.4769 0.4901], 'color', 'r', 'LineStyle', '--', 'LineWidth', 1);
% annotation('doublearrow',[0.4322 0.4292],[0.4857 0.3712]);
% h2 = annotation('line',[0.2705 0.538],[0.2127 0.2347]);
% set(h2,'X',[0.3616 0.3587],'Y',[0.4901 0.2215], 'LineStyle', '--', 'LineWidth', 1);
% h3 = annotation('line',[0.2675 0.2646],[0.4681 0.2127]);
% set(h3,'X',[0.3645 0.538],'Y',[0.2259 0.2347], 'LineStyle', '--', 'LineWidth', 1);
% annotation('doublearrow',[0.5027 0.5762],[0.4901 0.4945], 'color', 'r', 'LineStyle', '--', 'LineWidth', 0.25)

annotation('line',[0.5409 0.5409],[0.565 0.2743], 'color', 'r', 'LineStyle', '-.', 'LineWidth', 1);

text(45, 30, 6.5, 'S_{1}', 'FontName', 'Times New Roman', 'FontSize', 11);
% text(45.5, 37, 4.2, 'D_{1}', 'color', 'r', 'FontName', 'Times New Roman', 'FontSize', 11);
text(20, 30, 5.5, 'S_{2}', 'FontName', 'Times New Roman', 'FontSize', 11);
text(52.5, 35, 4.5, 'C', 'FontName', 'Times New Roman', 'FontSize', 11);
text(43, 35, 3, 'Axis l', 'FontName', 'Times New Roman', 'FontSize', 11);

exportgraphics(f6, 'Figure6.tiff', 'Resolution', 600);
% exportgraphics(f6, 'Figure6.pdf', 'ContentType', 'vector');
%% Fig 7 展示最小二乘平面的计算和回转轴的计算

f7 = figure("Color", [1,1,1]);
tiledlayout(2,3, 'TileSpacing', 'compact'); 
set(gcf, 'Units', 'centimeters', 'Position', [10,10,14,7]);

nexttile
pcshow(XY_1/6.4);
hold on;
[sf1, ~] = size(Flat_1);
pcshow(Flat_1/6.4, repmat([240,100,73]/255, sf1,1));
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 9 );
xlabel('x(\mum)', 'color', [0,0,0], 'position', [17,-17, 1]); 
ylabel('y(\mum)', 'color', [0,0,0], 'position', [-18,5, 1]);
% xlabel('x(\mum)', 'position', [17,-17, 1]); ylabel('y(\mum)', 'position', [-18,5, 1]); 
set(get(gca, 'xlabel'), 'Rotation', 22);
set(get(gca, 'ylabel'), 'Rotation', -37);


nexttile
pcshow(XY_2/6.4);
hold on;
[sf2, ~] = size(Flat_2);
pcshow(Flat_2/6.4, repmat([240,100,73]/255, sf2,1));
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 9 );
xlabel('x(\mum)', 'color', [0,0,0], 'position', [17,-17, 1]); 
ylabel('y(\mum)', 'color', [0,0,0], 'position', [-18,5, 1]);
% xlabel('x(\mum)', 'position', [17,-17, 1]); ylabel('y(\mum)', 'position', [-18,5, 1]); 
set(get(gca, 'xlabel'), 'Rotation', 22);
set(get(gca, 'ylabel'), 'Rotation', -37);

nexttile
pcshow(XY_3/6.4);
hold on;
[sf3, ~] = size(Flat_3);
pcshow(Flat_3/6.4, repmat([240,100,73]/255, sf3,1));
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 9 );
xlabel('x(\mum)', 'color', [0,0,0], 'position', [17,-22, 1]); 
ylabel('y(\mum)', 'color', [0,0,0], 'position', [-19,5, 1]);
% xlabel('x(\mum)', 'position', [17,-17, 1]); ylabel('y(\mum)', 'position', [-18,5, 1]); 
set(get(gca, 'xlabel'), 'Rotation', 22);
set(get(gca, 'ylabel'), 'Rotation', -37);

nexttile
pcshow(Col_1/6.4);
hold on;
line_z = linspace(min(Col_1(:,3)), max(Col_1(:,3)), 100);
line = [la_1(1,1)*line_z' + repmat(la_1(1,2),100,1), ...
        la_1(2,1)*line_z' + repmat(la_1(2,2),100,1), ...
        line_z'];
pcshow(line/6.4, repmat([240,100,73]/255, 100,1));
xlim([20,60]); ylim([20,60]);
set(gca, 'XTick', [20,40,60], 'YTick', [20,40,60]);
annotation(f7,'arrow','Position', [0.229 0.23 0 0.13],...
    'Color', [240,100,73]/255,'LineStyle','none',...
    'HeadStyle','vback1', ...
    'HeadWidth', 5);
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 9 );
xlabel('x(\mum)', 'color', [0,0,0], 'position', [17,-5, 1]); 
ylabel('y(\mum)', 'color', [0,0,0], 'position', [-5,5, 1]);
% xlabel('x(\mum)', 'position', [17,-17, 1]); ylabel('y(\mum)', 'position', [-18,5, 1]); 
set(get(gca, 'xlabel'), 'Rotation', 22);
set(get(gca, 'ylabel'), 'Rotation', -37);

nexttile
pcshow(Col_2/6.4);
hold on;
line_z = linspace(min(Col_2(:,3)), max(Col_2(:,3)), 100);
line = [la_2(1,1)*line_z' + repmat(la_2(1,2),100,1), ...
        la_2(2,1)*line_z' + repmat(la_2(2,2),100,1), ...
        line_z'];
pcshow(line/6.4, repmat([240,100,73]/255, 100,1));
xlim([20,60]); ylim([20,60]);
set(gca, 'XTick', [20,40,60], 'YTick', [20,40,60]);
annotation(f7,'arrow','Position', [0.512 0.23 0 0.13],...
    'Color', [240,100,73]/255,'LineStyle','none',...
    'HeadStyle','vback1', ...
    'HeadWidth', 5);
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 9 );
xlabel('x(\mum)', 'color', [0,0,0], 'position', [20,-3, 1]); 
ylabel('y(\mum)', 'color', [0,0,0], 'position', [-3,10, 1]);
% xlabel('x(\mum)', 'position', [17,-17, 1]); ylabel('y(\mum)', 'position', [-18,5, 1]); 
set(get(gca, 'xlabel'), 'Rotation', 22);
set(get(gca, 'ylabel'), 'Rotation', -37);

nexttile
pcshow(Col_3/6.4);
hold on;
line_z = linspace(min(Col_3(:,3)), max(Col_3(:,3)), 100);
line = [la_3(1,1)*line_z' + repmat(la_3(1,2),100,1), ...
        la_3(2,1)*line_z' + repmat(la_3(2,2),100,1), ...
        line_z'];
pcshow(line/6.4, repmat([240,100,73]/255, 100,1));
xlim([20,60]); ylim([20,60]);
set(gca, 'XTick', [20,40,60], 'YTick', [20,40,60]);
annotation(f7,'arrow','Position', [0.811 0.23 0 0.13],...
    'Color', [240,100,73]/255,'LineStyle','none',...
    'HeadStyle','vback1', ...
    'HeadWidth', 5);
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 9 );
xlabel('x(\mum)', 'color', [0,0,0], 'position', [17,-5, 1]); 
ylabel('y(\mum)', 'color', [0,0,0], 'position', [-5,5, 1]);
% xlabel('x(\mum)', 'position', [17,-17, 1]); ylabel('y(\mum)', 'position', [-18,5, 1]); 
set(get(gca, 'xlabel'), 'Rotation', 22);
set(get(gca, 'ylabel'), 'Rotation', -37);

% text(0, 0, 0, '(a)',  'FontName', 'Times New Roman', 'FontSize', 11);
% text(0, 0, '(b)',  'FontName', 'Times New Roman', 'FontSize', 11);

exportgraphics(f7, 'Figure7.tiff', 'Resolution', 600);
% exportgraphics(f7, 'Figure7.pdf', 'ContentType', 'vector');
%% Fig 8 曲面划分示意图

f8 = figure('color', 'w');
set(gcf, 'unit', 'centimeters', 'position', [10,10,9,9]);

% P(:,1) = P(:,1)*6.4;
% P(:,2) = P(:,2)*6.4;
P_ = P;
P_(:,3) = P_(:,3)+6.8;
pcshow(P_);
hold on;

% 画切分示意图，半径R_t/6.4
h_t = 3.4;
theta2 = linspace(0,2*pi,512);
circle = [R_t*cos(theta2')/6.4, R_t*sin(theta2')/6.4, h_t*ones(512,1)];
circle = [circle; [0,0,h_t]];
pcshow(circle, repmat([240,100,73]/255, 513, 1));
hold on;

R_n1 = linspace(0,R_t/6.4,512);
line_n1 = [zeros(512,1), -R_n1', h_t*ones(512,1)];
pcshow(line_n1, repmat([255,170,50]/255, 512,1));
hold on;
line_n2 = [R_n1'*sin(25*pi/180), -R_n1'*cos(25*pi/180), h_t*ones(512,1)];
pcshow(line_n2, repmat([255,170,50]/255, 512,1));
hold on;
theta1 = linspace(0,25*pi/180,30);
circle_t = [R_t*sin(theta1')/12, -R_t*cos(theta1')/12, h_t*ones(30,1)];
pcshow(circle_t, repmat([48,15,164]/255, 30, 1));
hold on;
theta3 = linspace(25*pi/180, 2*pi, 512);
circle_n3 = [R_t*sin(theta3')/24, -R_t*cos(theta3')/24, h_t*ones(512,1)];
pcshow(circle_n3, repmat([240,100,73]/255, 512, 1));

text(1,-0.5,3.3,'\theta_{1}', 'FontName', 'Times New Roman', 'FontSize', 11);
annotation('arrow',[0.170588235294118 0.491087962962963],...
    [0.691176470588235 0.482268518518518],...
    'Color',[240,100,73]/255,'LineStyle','none',...
    'HeadStyle','vback1');

set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11 );
xlabel('x/\mum', 'color', [0,0,0]); ylabel('y/\mum', 'color', [0,0,0]); zlabel('z/\mum', 'color', [0,0,0]);
set(get(gca, 'xlabel'), 'Rotation', 24);
set(get(gca, 'ylabel'), 'Rotation', -40);
set(gca, 'ztick', [0,3.4]);

exportgraphics(f8, 'Figure8.tiff', 'Resolution', 600);
exportgraphics(f8, 'Figure8.pdf', 'ContentType', 'vector');

%% Fig 9 旋转半径的计算

f9 = figure();
set(gcf, 'Units', 'centimeters', 'Position', [10,10,19,14]);

subplot(2,2,1);
plot(rol_2(:,1)/6.4, rol_2(:,2)/6.4+0.9, '.');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel({'r/\mum' ;'(a)'}); ylabel('h/\mum');
ylim([0.3, 0.45]);

subplot(2,2,2);
yyaxis left; % 双坐标画图
plot(r_2/6.4, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 1);
xlabel({'Groups'; '(b)'}); ylabel('R/\mum');
ylim([0,120/6.4]);
yyaxis right; 
plot(sr_2, '.', 'Color', 'r', "LineWidth", 1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
ylabel('Points Quantity');
% legend('R', 'Quantity')

subplot(2,2,3);
plot(rol_3(:,1)/6.4, rol_3(:,2)/6.4, ".");
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel({'r/\mum' ;'(c)'}); ylabel('h/\mum');
ylim([0.3, 0.4]);

subplot(2,2,4);
yyaxis left;
plot(r_3/6.4, 'Color', 'b', 'LineStyle', '-', 'LineWidth', 1);
xlabel({'Groups'; '(d)'}); ylabel('R/\mum');
ylim([0,80/6.4]);
yyaxis right;
plot(sr_3, '.', 'Color', 'r', "LineWidth", 1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
ylabel('Points Quantity');
% legend('R', 'Quantity');

exportgraphics(f9, 'Figure9.tiff', 'Resolution', 600);
exportgraphics(f9, 'Figure9.pdf', 'ContentType', 'vector');
%% 设计肩面倾斜的模型
% 注释详情参考model_full.m

phai = 0.00873;
% phai = 0;
Col=Column(phai); 
Col = Col+zeros(512,512);
% mesh(Col)
% 注：由于锥体是对称结构，只能修改曲面的单位来生成合适的模型
% 此时模型的整体单位应与旋转体的单位一致，默认高度为 3um。
% Col = ColumnWtR();
Pyr=11.6676*Pyramid(16)-2;
% if phai ~= 0
%     incf = sin(0.2*phai)*repmat(-255:256, 512,1);
%     Pyr=Pyr - incf';
% end
% mesh(Pyr);

i=Col<(2); j=Col>Pyr;
Z=Pyr.*i.*j+Col.*(ones(512,512)-i.*j);
% Z = Pyr.*j + Col*512/80;
% Z = Z*512/80;
% mesh(Z)

% 以上步骤结束后顶面、曲面、肩面都是倾斜的，倾斜0.00873rad，0.5度
% 将顶面变为水平，操作方法是用一个平面切削一下顶面，稍微改变两面距离
c = 12.5*ones(512,512)/6.4;
k = c < Z;
Z = c.*k + Z.*(ones(512,512)-k);

% 啊哈哈哈噪声来喽
Z = noise(6.4*Z);

[x,y] = size(Col);

X=[]; Y=[]; Z_=[];
% Z=ceil(Z*64)+64;
mesh(Z)
%% 设计顶面尺寸改变的模型
% 顶面直径： 9
% 
% 相对高度：1.5

phai = 0;
r = 4.5;

% 这里不用改
x_in=[4 7.7 8.32 9.72 12 10.23 16.59 17.72 20];%(5.63,-1)
y_in=[-0.5 -2.3 -2.6 -2.8 -3 -2.9 -2.6 -2.4 -2];

x_out=r:(80/512):20; 
y_out=spline(x_in,y_in,x_out)+2; % 这样做出来的图就可以和论文里一摸一样了

% 将剩余部分用 0 和 1.5 填满
x_f = 0:(80/512):(r-(80/512));
y_f = 1.5*ones(size(x_f));
[~, syf] = size(y_f); [~, syo] = size(y_out);
x = 1:256;
y = zeros(1,256);
y(1:syf) = y_f;
y(1+syf:syf+syo) = y_out;
% plot(y)
% rotation
theta=linspace(-pi,pi,100);
[Rr, Rt] = meshgrid(x, theta);
x=Rr.*cos(Rt);
yy=Rr.*sin(Rt);
if phai == 0
    z=repmat(y,100,1); %根据mesh函数的输入格式，扩充矩阵
% figure(2); mesh(x,yy,z);

% plot inclination col
% phai = 0.01; % rads
else
    [x, yy, z] = incli_col([1:256],y,phai);
end
% grid
% [X,Y]=meshgrid(-35:70/255:35,-35:70/255:35);
[X,Y]=meshgrid(-256:513/512:256,-256:513/512:256); % 这里可以将旋转曲面采样为任意大小的网格
Z=griddata(x,yy,z,X,Y,'cubic');

% 将补充的平面也进行倾斜
finc = repmat(sin(phai)*(-256:255), 512,1);
isn = isnan(Z);
Z(isn) = 0;
Z = Z + isn.*finc';

Col = Z;

Col = Col+zeros(512,512);
% mesh(Col)
% 注：由于锥体是对称结构，只能修改曲面的单位来生成合适的模型
% 此时模型的整体单位应与旋转体的单位一致，默认高度为 3um。
% Col = ColumnWtR();
Pyr=11.6676*Pyramid(16)-2;

i=Col<(2); j=Col>Pyr;
Z=Pyr.*i.*j+Col.*(ones(512,512)-i.*j);

% 啊哈哈哈噪声来喽
% Z = noise(6.4*Z);

[x,y] = size(Col);

X=[]; Y=[]; Z_=[];
% Z=ceil(Z*64)+64;
mesh(Z)
%% Fig 补 不同方法对的标注结果

f10 = figure('color', 'w');
tiledlayout(2,2, 'TileSpacing', 'none'); 
set(gcf, 'Units', 'centimeters', 'Position', [10,10,9,10]);

nexttile
pcshow(ptc_2);
hold on;
[s_2, ~] = size(XY_2);
pcshow(XY_2, repmat([240,100,73]/255, s_2,1));
view(2);
set(gca, 'XTick', [0,80], 'YTick', [0,80]);
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11 );
xlabel({'x(\mum)'; '(a)'}); ylabel({'y(\mum)'});
% xlabel({'x(Pixels)'; ' '}, 'color', [0,0,0]); ylabel({'y(Pixels)'}, 'color', [0,0,0]);
% set(get(gca, 'xlabel'), 'Rotation', 24);
% set(get(gca, 'ylabel'), 'Rotation', -40);

nexttile
pcshow('figure_2_rg.ply')
view(2);
set(gca, 'XTick', [0,80], 'YTick', []);
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11 );
xlabel({'x(\mum)'; '(d)'});
% xlabel({'x(Pixels)'; ' '}, 'color', [0,0,0]); ylabel({'y(Pixels)'}, 'color', [0,0,0]);
% set(get(gca, 'xlabel'), 'Rotation', 24);
% set(get(gca, 'ylabel'), 'Rotation', -40);

nexttile
pcshow(ptc_3);
hold on;
[s_3, ~] = size(XY_3);
pcshow(XY_3, repmat([240,100,73]/255, s_3,1));
view(2);
set(gca, 'XTick', [0,80], 'YTick', [0,80]);
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11 );
xlabel({'x(\mum)'; '(c)'}); ylabel({'y(\mum)'});
% xlabel({'x(Pixels)'; ' '}, 'color', [0,0,0]); ylabel({'y(Pixels)'}, 'color', [0,0,0]);
% set(get(gca, 'xlabel'), 'Rotation', 24);
% set(get(gca, 'ylabel'), 'Rotation', -40);

nexttile
pcshow('figure_3_rg.ply')
view(2);
set(gca, 'XTick', [0,80], 'YTick', []);
set(gca, 'color', 'w', 'xcolor', [0,0,0], 'ycolor', [0,0,0], 'zcolor', [0,0,0]);
set(gca, 'FontName', 'Times New Roman', 'Fontsize', 11 );
xlabel({'x(\mum)'; '(d)'});
% xlabel({'x(Pixels)'; ' '}, 'color', [0,0,0]); ylabel({'y(Pixels)'}, 'color', [0,0,0]);
% set(get(gca, 'xlabel'), 'Rotation', 24);
% set(get(gca, 'ylabel'), 'Rotation', -40);

% exportgraphics(f10, 'Figure10.tiff', 'Resolution', 600);
exportgraphics(f10, 'Figure10.emf', 'ContentType', 'vector');