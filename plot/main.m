clear,clc
%% 原点

oriDK200GeV4fm = importfile('ori_DK_200GeV_4fm.dat', 7, 507);
oriAi200GeV4fm = importfile('ori_Ai_200GeV_4fm.dat', 7, 507);
oriDK200GeV6fm = importfile('ori_DK_200GeV_6fm.dat', 7, 507);
oriAi200GeV6fm = importfile('ori_Ai_200GeV_6fm.dat', 7, 507);
oriDK200GeV8fm = importfile('ori_DK_200GeV_8fm.dat', 7, 507);
oriAi200GeV8fm = importfile('ori_Ai_200GeV_8fm.dat', 7, 507);


figure
semilogy(oriDK200GeV4fm(:,1),oriDK200GeV4fm(:,2),oriAi200GeV4fm(:,1),oriAi200GeV4fm(:,2))
title('Au-Au $\sqrt{s} = 200\,\mathrm{GeV}$ $b = 4\,\mathrm{fm}$', 'Interpreter', 'latex')
xlabel('$t(\mathrm{fm})$', 'Interpreter', 'latex')
ylabel('$eB(\mathrm{MeV}^2$)', 'Interpreter', 'latex')
legend('Kharzeev''s method', 'new method')
ylim([10^1,10^6])
print(gcf,'ori_200GeV_4fm.jpg', '-djpeg', '-r300');

figure
semilogy(oriDK200GeV6fm(:,1),oriDK200GeV6fm(:,2),oriAi200GeV6fm(:,1),oriAi200GeV6fm(:,2))
title('Au-Au $\sqrt{s} = 200\,\mathrm{GeV}$ $b = 6\,\mathrm{fm}$', 'Interpreter', 'latex')
xlabel('$t(\mathrm{fm})$', 'Interpreter', 'latex')
ylabel('$eB(\mathrm{MeV}^2$)', 'Interpreter', 'latex')
legend('Kharzeev''s method', 'new method')
ylim([10^1,10^6])
print(gcf,'ori_200GeV_6fm.jpg', '-djpeg', '-r300');

figure
semilogy(oriDK200GeV8fm(:,1),oriDK200GeV8fm(:,2),oriAi200GeV8fm(:,1),oriAi200GeV8fm(:,2))
title('Au-Au $\sqrt{s} = 200\,\mathrm{GeV}$ $b = 8\,\mathrm{fm}$', 'Interpreter', 'latex')
xlabel('$t(\mathrm{fm})$', 'Interpreter', 'latex')
ylabel('$eB(\mathrm{MeV}^2$)', 'Interpreter', 'latex')
legend('Kharzeev''s method', 'new method')
ylim([10^1,10^6])
print(gcf,'ori_200GeV_8fm.jpg', '-djpeg', '-r300');

%% x-y 平面图
% 导入数据
xyAiAu200GeV6fm001 = importfile1('x-y_Ai_Au_200GeV_6fm_0.01fm.dat', 1, 40401);
xyAiAu200GeV6fm005 = importfile1('x-y_Ai_Au_200GeV_6fm_0.05fm.dat', 1, 40401);
xyAiAu200GeV6fm01 = importfile1('x-y_Ai_Au_200GeV_6fm_0.1fm.dat', 1, 40401);
xyAiAu200GeV6fm02 = importfile1('x-y_Ai_Au_200GeV_6fm_0.2fm.dat', 1, 40401);
zmax = 5.3e4;
zmin = -4e4;
% t = 0.01 fm
Y = reshape(xyAiAu200GeV6fm001(:,1), [201, 201]);
X = reshape(xyAiAu200GeV6fm001(:,2), [201, 201]);
Z = reshape(xyAiAu200GeV6fm001(:,3), [201, 201]);
figure
meshc(X,Y,Z)
title('Au-Au $\sqrt{s} = 200\,\mathrm{GeV}$ $b = 6\,\mathrm{fm}$ $t = 0.01\,\mathrm{fm}$', 'Interpreter', 'latex')
xlabel('$z(\mathrm{fm})$', 'Interpreter', 'latex');
ylabel('$x(\mathrm{fm})$', 'Interpreter', 'latex');
zlabel('$eB_y(\mathrm{MeV}^2)$', 'Interpreter', 'latex');
zlim([zmin, zmax])
print(gcf,'x-z_200GeV_001fm.jpg', '-djpeg', '-r300');
% t = 0.05 fm
Y = reshape(xyAiAu200GeV6fm005(:,1), [201, 201]);
X = reshape(xyAiAu200GeV6fm005(:,2), [201, 201]);
Z = reshape(xyAiAu200GeV6fm005(:,3), [201, 201]);
figure
meshc(X,Y,Z)
title('Au-Au $\sqrt{s} = 200\,\mathrm{GeV}$ $b = 6\,\mathrm{fm}$ $t = 0.05\,\mathrm{fm}$', 'Interpreter', 'latex')
xlabel('$z(\mathrm{fm})$', 'Interpreter', 'latex');
ylabel('$x(\mathrm{fm})$', 'Interpreter', 'latex');
zlabel('$eB_y(\mathrm{MeV}^2)$', 'Interpreter', 'latex');
zlim([zmin, zmax])
print(gcf,'x-z_200GeV_005fm.jpg', '-djpeg', '-r300');
% t = 0.1 fm 
Y = reshape(xyAiAu200GeV6fm01(:,1), [201, 201]);
X = reshape(xyAiAu200GeV6fm01(:,2), [201, 201]);
Z = reshape(xyAiAu200GeV6fm01(:,3), [201, 201]);
figure
meshc(X,Y,Z)
title('Au-Au $\sqrt{s} = 200\,\mathrm{GeV}$ $b = 6\,\mathrm{fm}$ $t = 0.1\,\mathrm{fm}$', 'Interpreter', 'latex')
xlabel('$z(\mathrm{fm})$', 'Interpreter', 'latex');
ylabel('$x(\mathrm{fm})$', 'Interpreter', 'latex');
zlabel('$eB_y(\mathrm{MeV}^2)$', 'Interpreter', 'latex');
zlim([zmin, zmax])
print(gcf,'x-z_200GeV_010fm.jpg', '-djpeg', '-r300');
% t = 0.2 fm
Y = reshape(xyAiAu200GeV6fm02(:,1), [201, 201]);
X = reshape(xyAiAu200GeV6fm02(:,2), [201, 201]);
Z = reshape(xyAiAu200GeV6fm02(:,3), [201, 201]);
figure
meshc(X,Y,Z)
title('Au-Au $\sqrt{s} = 200\,\mathrm{GeV}$ $b = 6\,\mathrm{fm}$ $t = 0.2\,\mathrm{fm}$', 'Interpreter', 'latex')
xlabel('$z(\mathrm{fm})$', 'Interpreter', 'latex');
ylabel('$x(\mathrm{fm})$', 'Interpreter', 'latex');
zlabel('$eB_y(\mathrm{MeV}^2)$', 'Interpreter', 'latex');
zlim([zmin, zmax])
print(gcf,'x-z_200GeV_020fm.jpg', '-djpeg', '-r300');

