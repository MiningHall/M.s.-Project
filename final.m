%% 计算回转面的回转轴方向
% 注：这里只有方向是准确的
% 
% 对于方程为：(x-b)/a = (y-d)/c = z/1 的回转轴
% 
% 方程的最小二乘解计算式为：
% 
% [a b; c d] = [sum(x_i*z_i) sum(x_i); sum(y_i*z_i) sum(y_i)] * [sum(z_i*z_i) 
% sum(z_i); sum(z_i) n] ^ -1

% [size_col,~] = size(col);
% 这种好像不对
% col_x = col(:,1)'; col_y = col(:,2)'; col_z = [col(:,3)';ones(1,size_col)];
% abcd = [col_x;col_y]*col_z'/(col_z*col_z');

%% 
% 按如下方法，优化点到线的距离
% 
% 方程：
% 
% [x_i^2 x_i*y_i -x_i; x_i*y_i y_i^2 -y_i; x_i y_i -n] [m; n; D] = [x_i*z_i; 
% y_i*z_i; z_i]

% 最小二乘
c_xi2 = dot(col(:,1),col(:,1)); c_yi2 = dot(col(:,2),col(:,2));
c_xiyi = dot(col(:,1),col(:,2)); 
c_xizi = dot(col(:,1),col(:,3)); c_yizi = dot(col(:,2),col(:,3));
c_xi = sum(col(:,1)); c_yi = sum(col(:,2)); c_zi = sum(col(:,3));
mnD =  [c_xi2 c_xiyi -c_xi;...
        c_xiyi c_yi2 -c_yi;...
        c_xi c_yi -size_col]\...
       [-c_xizi; -c_yizi; -c_zi]
% 绘制回转轴
col_z = 30:60;
% col_x = abcd(1,2)+abcd(1,1)*col_z;
% col_y = abcd(2,2)+abcd(2,1)*col_z;
col_x = mnD(1)*col_z;
col_y = 3470 + mnD(2)*col_z;
col_xyz = [col_x;col_y;col_z]'; % 转换成n*3矩阵
col = [col;col_xyz];
pcshow(col);
%% 在轴向面上做投影，再求圆心                       

% 在平面上找一点定为x轴，再做叉乘确定y轴
z = [mnD(1),mnD(2),1];
% 假设我们的平面过原点，取点 [1,1,z] 用作计算
x = [1 0 -mnD(1)];
y = cross(x,z);
% 给三个坐标全都归一化
x = x/dot(x,x); y = y/dot(y,y); z = z/dot(z,z);
[x;y;z]*[x;y;z]'
% 计算旋转矩阵，用_来表示
col_ = ([x;y;z] * col')';

pcshow(col_)
%% 
% 这个地方算出来怎么看都不像是个对称轴

Carl = col;
%% 
% 重要：这里Carl还是用了Z来算

% 计算新坐标系的坐标轴向量
% x_为平面上一随机向量
% x_ = [x(1), y(1), C_z(1)] - [x(101), y(101), C_z(s)]
% x_ = [0,0,abc(3)] - [x(1),y(1),C_z(1)];
x_ = x_ / sqrt(dot(x_, x_));
z_ = [abc(1),abc(2),1]; z_ = z_/sqrt(dot(z_,z_));
y_ = cross(x_, z_); y_ = y_/sqrt(dot(y_, y_));
% 计算旋转矩阵，用_来表示
an = [x_;y_;z_] * [x_;y_;z_]'
Carl_ = ([x_;y_;z_] * Carl')';
% Cm = [min(Carl_(:,1)), min(Carl_(:,2)), min(Carl_(:,3))];
Cm = [0, 0, abc(3)];
Carl_ = Carl_ - repmat(Cm, s, 1);
% Carl_ = ([x_;y_;z_] * Carl')'
% 计算最小二乘圆心，仅仅考虑Carl_在x、y上的投影
xi2 = Carl_(:,1).*Carl_(:,1); yi2 = Carl_(:,2).*Carl_(:,2);
J = [2*sum(xi2), 2*dot(Carl_(:,1),Carl_(:,2)), -sum(Carl_(:,1)); ...
     2*dot(Carl_(:,1), Carl_(:,2)), 2*sum(yi2), -sum(Carl_(:,2)); ...
     2*sum(Carl_(:,1)), 2*sum(Carl_(:,2)), -s]
xy2 = xi2+yi2;
S = [dot(Carl_(:,1),xy2); dot(Carl_(:,2),xy2); sum(xy2)];
xyM = J\S;
R = sqrt(xyM(1)^2 + xyM(2)^2 - xyM(3));
xyR = [xyM(1), xyM(2), R]
% 原坐标系下画圆
circle = R*[cos([0:2*pi/500:2*pi]);sin([0:2*pi/500:2*pi]); zeros(1,501)] + repmat([xyR(1); xyR(2); 0], 1,501);
Circle = [x_;y_;z_]'*(circle + repmat(Cm', 1,501));
% 原坐标系下的圆心
xyzC = [x_;y_;z_]'*([xyR(1); xyR(2); 0] + Cm');
figure
pcshow(pointCloud(Carl))
hold on
pcshow(Carl_,[0.5,0,0])
hold on
pcshow(xyzC', [1,1,1])
hold on
pcshow(Circle', [1,1,0])
%%

% for i=2:628
%     col_test = [col_test;[10*sin((1:628)/100)',10*cos((1:628)/100)', 0.1*i*ones(628,1)];
% end

%% 
% 这里还需再把点云转换成平面，使用平面轮廓提取计算圆心和半径

out_line = zeros(512,512);
out_line(Y(:,1),Y(:,2)) = 1;
edge(out_line);
%% 
% 记录一下指标要求：
% 
% 顶面的圆度、柱面的圆柱度；（柱面本身也很难界定，这里再议）

% 角度条件

% MAT = MAT';
S = 4*pi*ones(512,512);
for i = 1:511
    for j = 1:511
        a = [1 0 MAT(i+1,j) - MAT(i,j)];
        b = [0 1 MAT(i+1,j+1) - MAT(i+1,j)];
        c = [-1 0 MAT(i,j+1) - MAT(i+1,j+1)];
        d = [0 -1 MAT(i,j) - MAT(i,j+1)];
%         sigma =  acos(dot(a,b)/(norm(a)*norm(b))) + acos(dot(b,c)/(norm(b)*norm(c))) +...
%                 acos(dot(c,d)/(norm(c)*norm(d))) + acos(dot(d,a)/(norm(d)*norm(a)));
        sigma = dot(a,b)/(norm(a)*norm(b)) + dot(b,c)/(norm(b)*norm(c)) +...
                dot(c,d)/(norm(c)*norm(d)) + dot(d,a)/(norm(d)*norm(a));
        S(i:i+1,j:j+1) = min(S(i:i+1,j:j+1),sigma*ones(2,2));
    end

end
%%
% 对角线距离条件
D = 1000*ones(512,512);
for i = 1:511
    for j = 1:511
        ac = sqrt(2+(MAT(i,j)-MAT(i+1,j+1))^2);
        bd = sqrt(2+(MAT(i+1,j)-MAT(i,j+1))^2);
        ed = abs(ac-bd);
        D(i:i+1,j:j+1) = min(D(i:i+1,j:j+1),ed*ones(2,2));
    end
end

D_ = ones(512,512)-(D-min(min(D)))/(max(max(D))-min(min(D)));

figure()
mesh(D_)
%%
% 四边距离条件
Dl = 100*ones(512,512);
for i = 1:511
    for j = 1:511
        d = max(max(MAT(i:i+1,j:j+1))) - min(min(MAT(i:i+1,j:j+1)));
        Dl(i:i+1,j:j+1) = min(Dl(i:i+1,j:j+1),d*ones(2,2));
    end
end
%%
% 主要用这个
% 方向条件
% MAT = MAT*6.4e6;
% MAT = MAT*1e6;
min_MAT = min(min(MAT));

MAT2 = zeros(514,514);
MAT2(2:513,2:513) = MAT;
MAT2(:,1) = MAT2(:,2);
MAT2(:,514) = MAT2(:,513);
MAT2(1,:) = MAT2(2,:);
MAT2(514,:) = MAT2(513,:);

% conv
s = zeros(512,512);
D = zeros(512,512);
for i = 2:513
    for j = 2:513
        D(i-1,j-1) = max(max(MAT2(i-1:i+1,j-1:j+1))) - min(min(MAT2(i-1:i+1,j-1:j+1)));
        a = MAT2(i-1,j); b = MAT2(i,j+1);c = MAT2(i+1,j); d = MAT2(i,j-1);
        C = MAT2(i,j);
        s1 = cross([1,0,C-a], [0,-1,C-b]);
        s2 = cross([0,-1,C-b], [-1,0,C-c]);
        s3 = cross([-1,0,C-c], [0,1,C-d]);
        s4 = cross([0,1,C-d], [1,0,C-a]);
        ss = [dot(s1,s1),dot(s2,s2),dot(s3,s3),dot(s4,s4)];
        s(i-1,j-1) = min(ss);
    end
end

mesh(s)
mesh(D)
mesh(s.*D)
%%
S = 2*pi*ones(512,512) -S;
S_= ones(512,512) - (S-min(min(S)))/(max(max(S))-min(min(S)));

figure()
mesh(S_)
%%
SD = Dl.*D_;
mesh(SD)
%%
% 按照比例，这里可以将z坐标乘6.4

[x,y]=size(MAT);
% MAT = MAT*6.4e6;

% 根据欧氏距离将 *顶面* 和 *肩面* 初步提取出来
t = s.*D;
% 给所有点按照欧式距离大小打上标签
P = zeros(x*y,4);
ptr = 1;
for i = 1:x
    for j = 1:y
        P(ptr,:) = [i, j, MAT(i,j), t(i,j)]; 
        ptr = ptr+1;
    end
end

% 将标签按位数取整，保留小数点后十位
P(:,4) = roundn(P(:,4),-20);
% pcshow(pointCloud(P(:,1:3)))

% 使用 tabulate 函数对标签进行统计，得到如：
%       Value Count Percent
%         0.9    20 100.00%
% 格式的数据；
% 且标签 tab 按 Value 的升序排列。
tab = tabulate(P(:,4));

% 将初始阈值设定为 30%
% 获得基础平面
% threshold = 25.76;
threshold = 28;
per= 0; ptr=0; % per 为百分比统计值；ptr 为计数器指针
while(per+tab(ptr+1,3)<threshold)
    ptr = ptr+1;
    per = per+tab(ptr,3);
end
P1 = P(:,4) < tab(ptr,1); % 提取到所有符合阈值的点
P1 = repmat(P1,1,4).*P;
Pl = P - P1;
P1(P1(:,1)==0,:)=[]; % 删除空行

% 从上面的三团点云中分离出顶面和肩面
P1cloud = pointCloud(P1(:,1:3));
figure()
pcshow(P1cloud)
dist_threshold = 10
[labels,numClusters] = pcsegdist(P1cloud, dist_threshold); % 对三组点云进行分割
tab_cluster = tabulate(labels);% 给阈值按数量大小排序
[~,I]=sort(tab_cluster(:,2),'descend');% 给点云标签数量按降序排列
tab_cluster=tab_cluster(I,:);

% 找出数量前两名标签所对应的坐标，第一名 X 代表肩面，第三名 Y 代表顶面
[x,~]=find(labels==tab_cluster(1,1)); [y,~]=find(labels==tab_cluster(3,1));
[z,~]=find(labels==tab_cluster(2,1));
X = P1(x,1:3);
Y = P1(y,1:3);
% 这个 Z 作为底面的圆环留备后用
Z = P1(z,1:3);
[s_z, ~] = size(Z);

% % 找出纵数量大小前两名的标签对应坐标，纵坐标第一名Y代表顶面，纵坐标第二名X代表肩面
% [x,~]=find(labels==tab_cluster(1,1)); [y,~]=find(labels==tab_cluster(3,1));
% [z,~]=find(labels==tab_cluster(2,1));

% 计算现两团点云在包括的点的个数，并更新 per
% 点云中共有262144个点
[s_x,~] = size(X);
[s_y, ~] = size(Y);
s = s_x+s_y;
per = 100*s/262144;
per_x = 100*s_x/262144
per_y = 100*s_y/262144

% 使用最小二乘方法，分别计算两面的倾斜角
abc_x = least_square_3d(X);
arg_x = atan(dot(abc_x(1:2),abc_x(1:2)));
abc_y = least_square_3d(Y);
arg_y = atan(dot(abc_y(1:2),abc_y(1:2)));

% 计算这两面所占的比例，修正threshold
% 修正使用的计算方法为：按照面积比来计算倾斜区域内的应有点数
% 区域生长迭代的方向为占比 25.76% 的点云数据，其中：
%   top 包含 1593 点，占比 0.61%
%   middle 包含 65926 点，占比 25.15%
%   两者共占比 25.76%
% 倾斜后的面积与原面积比为 cos(arg_):1
threshold = cos(arg_x)*25.15 + cos(arg_y)*0.61;
% threshold_x = cos(arg_x)*25.15;
% threshold_y = cos(arg_y)*0.61;
threshold_x = cos(arg_x)*22.1;
threshold_y = cos(arg_y)*0.6;
[tab_size,~] = size(tab)

% 以这两个面的点云为中心，开始做区域生长。
% 开始进行迭代，每轮迭代后根据方向向量更新阈值
while(per < threshold)
    c_x = sum(X)/s_x;
    c_y = sum(Y)/s_y;
    if ptr == tab_size
        break
    end
    % 找出需要生长的下一个尺度内的点，以总数的 1% 为尺度
    n = 0; % n 为百分比计数器
    while(ptr < tab_size && (n+tab(ptr+1,3))<0.5) % 找出下个 1% 的点
        ptr = ptr+1;
        n = n+tab(ptr,3);
    end
    p = Pl(:,4) < tab(ptr); % 找出下个1%的点
    p = repmat(p,1,4).*Pl;
    Pl = Pl-p; % 保存剩下的点
    
    p = p(:,1:3); % p 只取坐标，并删除 p 和 Pl 的空行
    p(p(:,1)==0,:)=[];
    Pl(Pl(:,1)==0,:)=[];
    
    
    % 确认这些点是否与 X 和 Y 邻接
    % 这里将 x 与 y 分开进行判断，获得更准确的顶面结果
    [s_n,~] = size(p);
    if per_x < threshold_x
        d_x = p - repmat(c_x, s_n, 1);
        d_x = sqrt(sum(d_x.*d_x,2)); % 按行求和计算距离
        % 如果距离小于阈值，则把该点加入点云集合中
        % 这里将到圆心的距离变为唯一距离指标
        dx_p = repmat((d_x > 12.5*dist_threshold).*(d_x < 20*dist_threshold), 1, 3) .* p;
        p = p - dx_p;
        
        dx_p(dx_p(:,1)==0,:) = [];
        p(p(:,1)==0,:) = [];
        X = [X;dx_p];
        [s_n,~] = size(p);
%         for i = 1:s_n
%         % 分别计算距离的平方
%         d_x = repmat(p(i,:),s_x,1) - X;
%         d_X = sqrt(dot(d_x',d_x')); % 这里必须要用转置才可以
% 
%         % 如果距离小于阈值，则把该点加入点云集合中
% %         if min(d_X) < dist_threshold
% %             s_x = s_x + 1;
% %             X(s_x,:) = p(i,:);
% %         elseif min(d_Y) < dist_threshold
% %             s_y = s_y + 1;
% %             Y(s_y,:) = p(i,:);
% %         end
%             if min(d_X) < dist_threshold && max(d_X) < 30*dist_threshold
%                 s_x = s_x + 1;
%                 X(s_x,:) = p(i,:);
%             end
%         end
    end
    if per_y < threshold_y
        d_y = p - repmat(c_x, s_n, 1);
        d_y = sqrt(sum(d_y.*d_y,2)); % 按行求和计算距离
        % 如果距离小于阈值，则把该点加入点云集合中
        dy_p = repmat(d_y < 3*dist_threshold, 1, 3) .* p;
        p = p - dy_p;
        
        dy_p(dy_p(:,1)==0,:) = [];
        p(p(:,1)==0,:) = [];
        Y =[Y; dy_p];
        
%         for i = 1:s_n
%             d_y = repmat(p(i,1:3),s_y,1) - Y;
%             d_Y = sqrt(dot(d_y',d_y'));
%             % 如果距离小于阈值，则把该点加入点云集合中
%     %         if min(d_X) < dist_threshold
%     %             s_x = s_x + 1;
%     %             X(s_x,:) = p(i,:);
%     %         elseif min(d_Y) < dist_threshold
%     %             s_y = s_y + 1;
%     %             Y(s_y,:) = p(i,:);
%     %         end
%             if min(d_Y) < dist_threshold && max(d_Y) < 5*dist_threshold
%                 s_y = s_y + 1;
%                 Y(s_y,:) = p(i,:);
%             end
%         end
    end
    
    % 使用最小二乘方法，分别计算更新后 X 和 Y 的倾斜角
    abc_x = least_square_3d(X);
    arg_x = atan(sqrt(dot(abc_x(1:2),abc_x(1:2))));
    abc_y = least_square_3d(Y);
    arg_y = atan(sqrt(dot(abc_y(1:2),abc_y(1:2))));
    
    % 根据倾斜角修正阈值
    threshold_x = cos(arg_x)*threshold_x;
    threshold_y = cos(arg_y)*threshold_y;
    threshold = threshold_x + threshold_y
    
    % 更新 X 和 Y 的大小和占比，更新X，Y点集占比
    [s_x,~] = size(X);
    [s_y, ~] = size(Y);
    
    per_x = 100*s_x/262144;
    per_y = 100*s_y/262144;
    
    s = s_x+s_y;    
    per = per_x + per_y;
%     n_per = 100*s/262144 - per
%     if (100*s/262144 - per) < 0.0005
%         break
%     end
end
per
% 整理点云，填充点云中出现的空隙

% 画图
figure()
pcshow(X);
figure()
pcshow(Y);
% 算得最小二乘平面方程和平面的方向向量
% 算x
Carl = X;
[s,~] = size(Carl);
x = Carl(:,1); y = Carl(:,2); z = Carl(:,3);
x2 = dot(x,x); y2 = dot(y,y);
xy = dot(x,y); xz = dot(x,z); yz = dot(y,z);
abc_s = [x2, xy, -sum(x); ...
       xy, y2, -sum(y); ...
       -sum(x), -sum(y), s] \[-xz; -yz; sum(z)]
S_z = repmat(abc_s(3), s, 1) - abc_s(1)*x - abc_s(2)*y;
figure
pcshow(pointCloud(Carl))
hold on
flat_X = [x,y,S_z];
pcshow(flat_X,[0.5,0,0])

% 算y
Carl = Y;
[s,~] = size(Carl);
x = Carl(:,1); y = Carl(:,2); z = Carl(:,3);
x2 = dot(x,x); y2 = dot(y,y);
xy = dot(x,y); xz = dot(x,z); yz = dot(y,z);
abc_t = [x2, xy, -sum(x); ...
       xy, y2, -sum(y); ...
       -sum(x), -sum(y), s] \[-xz; -yz; sum(z)]
T_z = repmat(abc_t(3), s, 1) - abc_t(1)*x - abc_t(2)*y;
figure
pcshow(pointCloud(Carl))
hold on
flat_Y = [x,y,T_z];
pcshow(flat_Y,[0.5,0,0])
%% 划分得到旋转曲面
% matalb 距离分割函数：pcsegdist(pointCloud(col_pyr),5);

col_mat = zeros(512,512);
for i=1:s_x
   col_mat(X(i,1),X(i,2)) = X(i,3); 
end
for i=1:s_y
    col_mat(Y(i,1),Y(i,2)) = Y(i,3); 
end
% for i=1:s_z
%     col_mat(Z(i,1),Z(i,2)) = Z(i,3); 
% end
col_mat = MAT - col_mat;

% 提取到剩余的点
col_pyr = turn_to_pointcloud(col_mat);

% 使用距离进行分割
cp_labels = pcsegdist(pointCloud(col_pyr),5);

% 选出第一个和第二个最大子集
col_tab = tabulate(cp_labels);
[~,I]=sort(col_tab(:,3),'descend');
pyr = col_pyr.*repmat(cp_labels==I(1),1,3);
col = col_pyr.*repmat(cp_labels==I(2),1,3);
col(col(:,1)==0,:)=[];
pyr(pyr(:,1)==0,:)=[];
pcshow(col);
%% 使用切片迭代的方法计算回转轴
% 计算流程：
% 
% 切片->平面投影->计算最小二乘圆心.->拟合直线->旋转->重新切片投影并计算圆心拟合直线
% 
% 终止条件(未设置)：
% 
% 循环10次

% input col
max_height = max(col(:,3));
min_height = min(col(:,3));
mn = (max_height-min_height)/20;
start_height = min_height + mn;

% set a vector to save center points
[col_size,~] = size(col);
col_center = [];

while(start_height < max_height)
    % got one chip and the points left
    chip = col(:,3) < start_height;
    chip = chip.*col(:,3) > (start_height-mn);
    chip = repmat(chip,1,3) .* col;
    % remove zeros, calculate size
    chip(chip(:,3)==0,:) = []; [chip_size,~] = size(chip);
    % calculate center point of chip and get its 3d location
    abc = least_square_3d(chip);
    chip_flatz = abc(3)*ones(chip_size,1) - abc(1)*chip(:,1) - abc(2)*chip(:,2); % points on flat
    chip_flat = [chip(:,1), chip(:,2), chip_flatz];
    [chip_center, chip_R] = least_square_2dcircle(chip_flat(:,1:2));
    % store the location in vector
    col_center = [col_center; chip_center(1), chip_center(2), start_height];
    % set threshold of chips
    start_height = start_height + (max_height-min_height)/20;
end
% col_center(:,3) = col_center(:,3)*1e7;

% 去掉偏心太多的点
center_density = [sum(col(:,1))/col_size sum(col(:,2))/col_size];
col_center(abs(col_center(:,1)-center_density(1))>10,:) = [];
col_center(abs(col_center(:,2)-center_density(2))>10,:) = [];
[ccs, ~] = size(col_center);
col_center(col_center(:,3)<sum(col_center(:,3))/ccs, :) = [];

% calculate the line and the direction
line_args = least_square_3dline(col_center);
line_z = -25:0.125:25;
line = [line_args(1,1)*line_z' + repmat(line_args(1,2),401,1), ...
        line_args(2,1)*line_z' + repmat(line_args(2,2),401,1), ...
        line_z'];
% rotate pointcloud into the drection of line 
% [axis, col_] = rotate_3d([line_args(1,1), line_args(2,1), 0], col);

pcshow(line,[0,0,1]);
hold on
% pcshow(col);
% hold on
pcshow(col_center,[1,0,0])
%% 计算轴与平面的夹角

line_direction = [line_args(1,1), line_args(2,1), 1]
% 顶面夹角
top_direction = [abc_t(1),abc_t(2), 1];
sigma_lt = acos(dot(line_direction,top_direction)/(norm(line_direction)*norm(top_direction)))
sigma_lt_d = sigma_lt*180/pi
% 肩面夹角
sho_direction = [abc_s(1),abc_s(2), 1];
sigma_ls = acos(dot(line_direction,sho_direction)/(norm(line_direction)*norm(sho_direction)))
sigma_ls_d = sigma_ls*180/pi
%% 计算两平面之间的夹角

% 平面夹角
sigma_ts = acos(dot(top_direction,sho_direction)/(norm(top_direction)*norm(sho_direction)))
sigma_ts_d = sigma_ts*180/pi
%%
% 循环两次的实际轴向：
r*[line_args(1,1),line_args(2,1),1]'
%% 计算如下参数，并将数据保存用以画图
% 两个平面的轮廓半径 || 平面质心之间的距离
% 
% 两平面夹角 || 平面与曲面回转轴的夹角
%% 利用点云的 x、y 坐标提取出轮廓，并划出顶面外圆

[y_size,~]=size(Y);
out_line = zeros(512,512);
for i=1:y_size
    out_line(Y(i,1),Y(i,2)) = 1;
end
I = edge(out_line);
Y_outline = turn_to_pointcloud(I.*MAT);
[outline_Size,~]=size(Y_outline);
Y_outline = MAT.*I;

Y_outline = turn_to_pointcloud(Y_outline);
pcshow(Y_outline)
% a*x + b*y + z = c
abc = least_square_3d(Y_outline);
Y_outline_z = abc(3)-abc(1)*Y_outline(:,1)-abc(2)*Y_outline(:,2);
outline_flat = [Y_outline(:,1), Y_outline(:,2), Y_outline_z];
% 旋转，函数输入为平面的三个参数
[axis_t, outline_flat_rotate] = rotate_3d(abc, outline_flat);
h = outline_flat_rotate(1,3);
outline_flat_rotate(:,3) = zeros(outline_Size,1);
% 计算圆及圆心

Carl_ = outline_flat_rotate;
% 计算最小二乘圆心，仅仅考虑Carl_在x、y上的投影
[center, R_t] = least_square_2dcircle(outline_flat_rotate(:,1:2))
% 修整算法，在第一个圆内的点统统删掉，再重新计算
% 用这种方法筛出内部的闭环轮廓
o_outline = sum((outline_flat_rotate(:,1:2)-repmat(center',outline_Size,1))...
          .*(outline_flat_rotate(:,1:2)-repmat(center',outline_Size,1)) , 2) > R_t^2;
o_outline = repmat(o_outline, 1, 3).*outline_flat_rotate;
o_outline(o_outline(:,2) == 0, :) = [];

[center, R_t] = least_square_2dcircle(o_outline)
R_t/3.2
% % 原坐标系下画圆
circle = R_t*[cos([0:2*pi/500:2*pi]);sin([0:2*pi/500:2*pi]); zeros(1,501)]' + repmat([center(1), center(2), h], 501,1);
Circle = axis_t*circle';
% 绘制原坐标系下的圆心
Center_t = axis_t'*([center(1); center(2); h])
figure
pcshow(pointCloud(Y))
hold on
pcshow(outline_flat,[0.5,0,0])
hold on
pcshow(Circle', [1,1,0])
%% 使用几乎相同的方法计算肩面
% 注意微调数据并将用作画图的轮廓保存起来

[x_size,~]=size(X);
out_line = zeros(512,512);
for i=1:x_size
    out_line(X(i,1),X(i,2)) = 1;
end
I = edge(out_line,'Canny');
X_outline = MAT.*I;
X_outline = turn_to_pointcloud(X_outline);
X_outline_tab = pcsegdist(pointCloud(X_outline), 4);
XO_tab = tabulate(X_outline_tab);
[~,XO_I]=sort(XO_tab(:,3),'descend');
% XO_I(2)是可调参数
X_outline = X_outline.*repmat(X_outline_tab==XO_I(1),1,3);

X_outline(X_outline(:,1)==0,:) = [];
% X_outline(X_outline(:,1)>400,:) = [];
[outline_Size,~]=size(X_outline);
pcshow(X_outline)
% a*x + b*y + z = c
abc = least_square_3d(X_outline);
X_outline_z = abc(3)-abc(1)*X_outline(:,1)-abc(2)*X_outline(:,2);
outline_flat = [X_outline(:,1), X_outline(:,2), X_outline_z];
% 旋转，函数输入为平面的三个参数
[axis, outline_flat_rotate] = rotate_3d(abc, outline_flat);
h = outline_flat_rotate(1,3);
outline_flat_rotate(:,3) = zeros(outline_Size,1);
% 计算圆及圆心
Carl_ = outline_flat_rotate(:, 1:2);
% 计算最小二乘圆心，仅仅考虑Carl_在x、y上的投影
[center, R_s] = least_square_2dcircle(outline_flat_rotate(:,1:2))
% % 圆周滤波，舍弃掉圆周上突变的点（均值滤波）
[c_size,~] = size(Carl_);
r_c = sqrt(sum((Carl_ - repmat(center', c_size, 1)).*(Carl_ - repmat(center', c_size, 1)), 2));
r_a = sum(r_c)/c_size; % 平均半径
outline_flat_rotate(r_c > r_a, :) = [];
[center, R_s] = least_square_2dcircle(outline_flat_rotate(:,1:2))
R_s/3.2
% 原坐标系下画圆
circle = R_s*[cos([0:2*pi/500:2*pi]);sin([0:2*pi/500:2*pi]); zeros(1,501)]' + repmat([center(1), center(2), h], 501,1);
Circle = axis*circle';

% 绘制原坐标系下的圆心
Center_s = axis*([center(1), center(2), h]')
figure
pcshow(pointCloud(X))
hold on
pcshow(outline_flat,[0.5,0,0])
hold on
pcshow(Circle', [1,1,0])
%% 计算两平面质心之间的距离

dist_1 = Center_s - Center_t;
d1 = sqrt(dot(dist_1,dist_1))/6.4
%% 计算顶面到曲面底部的距离
% 按照顶面和肩面的平均方向旋转曲面，得到其自顶向下的距离

[~, col_d] = rotate_3d(abc, col);
% col_tab = tabulate(col_d(:,3));
% [col_x, col_y] = find(col_d ==min(col_d(:,3)));
% col_l = col_d(col_x,:);
% dist_2 = col_l' - Center_t;
% d2 = sqrt(dot(dist_2,dist_2))
% OC = line_direction/sqrt(dot(line_direction,line_direction));
% dc = cross(dist_2,OC);
% d3 = sqrt(dot(dc,dc));
% 
% d2 = sqrt(d2^2 - d3^2)/6.4
d2 = (max(col_d(:,3)) - min(col_d(:,3)))/6.4
%% 计算平面和曲面的过渡半径
% 这里只计算顶面和曲面的夹角
% 
% 计算原理（参考老师给的文献）：
% 
% 根据圆的圆心和范围，划定小段旋转曲面作为计算范围

% 划分过渡部分的点云
mat_y = turnto_mat(Y);
mat_f = turnto_mat(flat_Y);
mat_b = (mat_y < mat_f).*mat_y;
Blunt = turn_to_pointcloud(mat_b);
[rol, r,sr, P] = Blunt_R(Blunt, Center_t, axis_t, 0.125*pi);
% r是过渡半径，sr是得到半径处的采样点数
r(r==0)=[]; sr(sr==0)=[];
figure()
plot(r)
hold on
plot(sr,'.')
%% 设计不同模型并计算画图

x_in=[4 7.2 8.32 9.72 12 10.23 16.59 17.72 20];%(5.63,-1)
y_in=[0 -2 -2.4 -2.8 -3 -2.9 -2.6 -2.4 -2];

x_out=4:(80/512):20; % 这里，生成按比例放大的效果
y_out=spline(x_in,y_in,x_out)+2; % 这样做出来的图就可以和论文里一摸一样了


% 将剩余部分用 0 和 2 填满
x_f = 0:(80/512):(4-(80/512));
y_f = 2*ones(size(x_f));
[~, syf] = size(y_f); [~, syo] = size(y_out);
x = 1:256;
y = zeros(1,256);
y(1:syf) = y_f;
y(1+syf:syf+syo) = y_out;
figure
x = [-fliplr(x), x]/6.4;
y = [fliplr(y), y]+1;
y(449:512) = 1:-1/64:1/64;
y(1:64) = 1/64:1/64:1;
plot(x,y)
% annotation('doublearrow',[0.3125, 1/3],[25/80  1])