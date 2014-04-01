close all;
clear all;
%%Setup simulations

f_prefix='2D_weak_1'

NE=[8,16,32,64];
NT=[8,16,32,64];
DELTA=[1];
DELTA_BC=[1, 100, 10000];

%%Set up the correct directories

%loztop
%exe_directory = '/home/loztop/Dropbox/Dphil/libmesh_git/weak_poro_paper/2D_convergence_sinh/'
%output_plot_folder='/home/loztop/Dropbox/Dphil/libmesh_git/weak_poro_paper/2D_convergence_sinh/matlab_files/converge/pdf_plots/'

%clpc59
exe_directory = '/users/lorenzb/Dphil/libmesh_git/weak_poro_paper/2D_convergence_sinh/'
output_plot_folder='/home/scratch/Dropbox/Dphil/libmesh_git/weak_poro_paper/2D_convergence_sinh/matlab_files/converge/pdf_plots/'

exe_filename = 'ex11-opt'

res_directory =strcat([exe_directory 'data/matfiles/'])
res_directory_plot =strcat([exe_directory 'data/plots/'])

str_nt='NT_';
str_ne='NE_';

format long;


for (i=1:length(DELTA_BC))
 for(j=1:length(NE))

error_file_name = strcat([res_directory f_prefix,...
     '_' num2str(DELTA_BC(i)) '_' str_nt int2str(NT(j)) '_' str_ne int2str(NE(j)) '_.mat']);

A = load('-ascii', error_file_name)


L2_u_h(i,j,1:2)=A(1:2);  
L2_u_h(i,j,3)=A(3);

H1_u_h(i,j,1:2)=A(1:2); 
H1_u_h(i,j,3)=A(4);

L2_w_h(i,j,1:2)=A(1:2);  
L2_w_h(i,j,3)=A(5);
Hd_w_h(i,j,1:2)=A(1:2);  
Hd_w_h(i,j,3)=A(6);

L2_p_h(i,j,1:2)=A(1:2);   %IS INF
L2_p_h(i,j,3)=A(7);

 end
end
 
 %%
%%Do actual plotting
close all;

%DISP_Linfty_H1
h=figure;
set([gca]             , ...
    'FontSize'   , 12           );

x_axis_NE=log2(NE(1:end));
H1_u_h_hand(1)=plot(x_axis_NE,log2(H1_u_h(1,1:end,3)),'--kd','LineWidth',2,'MarkerSize',8)
hold on;
H1_u_h_hand(2)=plot(x_axis_NE,log2(H1_u_h(2,1:end,3)),'-kx','LineWidth',2,'MarkerSize',10)
hold on;
H1_u_h_hand(3)=plot(x_axis_NE,log2(H1_u_h(3,1:end,3)),'--ko','LineWidth',2,'MarkerSize',8)
hold on;
H1_u_h_hand(4)=plot(x_axis_NE,1*(-x_axis_NE+x_axis_NE(1)) -5,'--b','LineWidth',2,'MarkerSize',8)
hold on;

title('2D convergence of $||\mathbf{u}-\mathbf{u}_{h}||_{l^{\infty}(H^{1})}$','interpreter','latex','FontSize',19);
xlabel('$\log(1/\Delta t)$, $\log(1/\Delta h)$ ','interpreter','latex','FontSize',19)
ylabel('$\log(\mbox{error})$ ','interpreter','latex','FontSize',19)
 grid on;
set(gca,'YTick',[-11:1:-4])
set(gca,'XTick',[2:1:8])
axis([2 8 -11 -4])
set(h, 'Position', [100 100 700 700])

hLegend = legend( ...
   [H1_u_h_hand(1),H1_u_h_hand(2),H1_u_h_hand(3),H1_u_h_hand(4)], ...
 strcat(['\gamma=',num2str(DELTA_BC(1))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(2))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(3))]) , ...
   '1st order ' , ...
   'location', 'SouthWest');
print(h,strcat([output_plot_folder,'disp_linf_h1']),'-dpdf','-r0')




%DISP_L2_L2
h=figure;
set([gca]             , ...
    'FontSize'   , 12           );

x_axis_NE=log2(NE(1:end));
L2_u_h_hand(1)=plot(x_axis_NE,log2(L2_u_h(1,1:end,3)),'--kd','LineWidth',2,'MarkerSize',8)
hold on;
L2_u_h_hand(2)=plot(x_axis_NE,log2(L2_u_h(2,1:end,3)),'-kx','LineWidth',2,'MarkerSize',10)
hold on;
L2_u_h_hand(3)=plot(x_axis_NE,log2(L2_u_h(3,1:end,3)),'--ko','LineWidth',2,'MarkerSize',8)
hold on;
L2_u_h_hand(4)=plot(x_axis_NE,2*(-x_axis_NE+x_axis_NE(1)) -8,'--b','LineWidth',2,'MarkerSize',8)
hold on;

title('2D convergence of $||\mathbf{u}-\mathbf{u}_{h}||_{l^{2}(L^{2})}$','interpreter','latex','FontSize',19);
xlabel('$\log(1/\Delta t)$, $\log(1/\Delta h)$ ','interpreter','latex','FontSize',19)
ylabel('$\log(\mbox{error})$ ','interpreter','latex','FontSize',19)
 grid on;
set(gca,'YTick',[-19:1:-6])
set(gca,'XTick',[2:1:8])
axis([2 8 -19 -6])
set(h, 'Position', [100 100 700 700])

hLegend = legend( ...
   [L2_u_h_hand(1),L2_u_h_hand(2),L2_u_h_hand(3),L2_u_h_hand(4)], ...
    strcat(['\gamma=',num2str(DELTA_BC(1))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(2))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(3))]) , ...
   '2nd order ' , ...
   'location', 'SouthWest');
print(h,strcat([output_plot_folder,'disp_l2_l2']),'-dpdf','-r0')






%P_Linf_L2
h=figure;
set([gca]             , ...
    'FontSize'   , 12           );

x_axis_NE=log2(NE(1:end));
L2_p_h_hand(1)=plot(x_axis_NE,log2(L2_p_h(1,1:end,3)),'--kd','LineWidth',2,'MarkerSize',8)
hold on;
L2_p_h_hand(2)=plot(x_axis_NE,log2(L2_p_h(2,1:end,3)),'-kx','LineWidth',2,'MarkerSize',10)
hold on;
L2_p_h_hand(3)=plot(x_axis_NE,log2(L2_p_h(3,1:end,3)),'--ko','LineWidth',2,'MarkerSize',8)
hold on;
L2_p_h_hand(4)=plot(x_axis_NE,1*(-x_axis_NE+x_axis_NE(1)) -1,'--b','LineWidth',2,'MarkerSize',8)
hold on;

title('2D convergence of $||{p}-{p}_{h}||_{l^{\infty}(L^{2})}$','interpreter','latex','FontSize',19);
xlabel('$\log(1/\Delta t)$, $\log(1/\Delta h)$ ','interpreter','latex','FontSize',19)
ylabel('$\log(\mbox{error})$ ','interpreter','latex','FontSize',19)
 grid on;
set(gca,'YTick',[-11:1:0])
set(gca,'XTick',[2:1:8])
axis([2 8 -11 0])
set(h, 'Position', [100 100 700 700])

hLegend = legend( ...
   [L2_p_h_hand(1),L2_p_h_hand(2),L2_p_h_hand(3),L2_p_h_hand(4)], ...
   strcat(['\gamma=',num2str(DELTA_BC(1))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(2))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(3))]) , ...
   '1st order ' , ...
   'location', 'SouthWest');
print(h,strcat([output_plot_folder,'p_linf_l2']),'-dpdf','-r0')









%FLUX_L2_div
h=figure;
set([gca]             , ...
    'FontSize'   , 12           );

x_axis_NE=log2(NE(1:end));
Hd_w_h_hand(1)=plot(x_axis_NE,log2(Hd_w_h(1,1:end,3)),'--kd','LineWidth',2,'MarkerSize',8)
hold on;
Hd_w_h_hand(2)=plot(x_axis_NE,log2(Hd_w_h(2,1:end,3)),'-kx','LineWidth',2,'MarkerSize',10)
hold on;
Hd_w_h_hand(3)=plot(x_axis_NE,log2(Hd_w_h(3,1:end,3)),'--ko','LineWidth',2,'MarkerSize',8)
hold on;
Hd_w_h_hand(4)=plot(x_axis_NE,1*(-x_axis_NE+x_axis_NE(1)) -1,'--b','LineWidth',2,'MarkerSize',8)
hold on;

title('2D convergence of $||\nabla \cdot( \mathbf{z}-\mathbf{z}_{h})||_{l^{2}(L^{2})}$','interpreter','latex','FontSize',19);
xlabel('$\log(1/\Delta t)$, $\log(1/\Delta h)$ ','interpreter','latex','FontSize',19)
ylabel('$\log(\mbox{error})$ ','interpreter','latex','FontSize',19)
 grid on;
set(gca,'YTick',[-9:1:2])
set(gca,'XTick',[2:1:8])
axis([2 8 -9 2])
set(h, 'Position', [100 100 700 700])

hLegend = legend( ...
   [Hd_w_h_hand(1),Hd_w_h_hand(2),Hd_w_h_hand(3),Hd_w_h_hand(4)], ...
  strcat(['\gamma=',num2str(DELTA_BC(1))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(2))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(3))]) , ...
   '1st order ' , ...
   'location', 'SouthWest');
print(h,strcat([output_plot_folder,'flux_l2_div']),'-dpdf','-r0')





%FLUX_L2_L2
h=figure;
set([gca]             , ...
    'FontSize'   , 12           );

x_axis_NE=log2(NE(1:end));
L2_w_h_hand(1)=plot(x_axis_NE,log2(L2_w_h(1,1:end,3)),'--kd','LineWidth',2,'MarkerSize',8)
hold on;
L2_w_h_hand(2)=plot(x_axis_NE,log2(L2_w_h(2,1:end,3)),'-kx','LineWidth',2,'MarkerSize',10)
hold on;
L2_w_h_hand(3)=plot(x_axis_NE,log2(L2_w_h(3,1:end,3)),'--ko','LineWidth',2,'MarkerSize',8)
hold on;
L2_w_h_hand(4)=plot(x_axis_NE,2*(-x_axis_NE+x_axis_NE(1)) -4,'--b','LineWidth',2,'MarkerSize',8)
hold on;

title('2D convergence of $|| \mathbf{z}-\mathbf{z}_{h}||_{l^{2}(L^{2})}$','interpreter','latex','FontSize',19);
xlabel('$\log(1/\Delta t)$, $\log(1/\Delta h)$ ','interpreter','latex','FontSize',19)
ylabel('$\log(\mbox{error})$ ','interpreter','latex','FontSize',19)
 grid on;
set(gca,'YTick',[-14:1:-1])
set(gca,'XTick',[2:1:8])
axis([2 8 -14 -1])
set(h, 'Position', [100 100 700 700])

hLegend = legend( ...
   [L2_w_h_hand(1),L2_w_h_hand(2),L2_w_h_hand(3),L2_w_h_hand(4)], ...
  strcat(['\gamma=',num2str(DELTA_BC(1))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(2))]) , ...
  strcat(['\gamma=',num2str(DELTA_BC(3))]) , ...
   '2nd order ' , ...
   'location', 'SouthWest');
print(h,strcat([output_plot_folder,'flux_l2_l2']),'-dpdf','-r0')


 