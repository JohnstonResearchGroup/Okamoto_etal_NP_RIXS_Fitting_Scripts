A = load('./Data/data_B.txt');
B = load('./Data/RIXS_exp_data_o.txt');


figure(1); hold on;
plot(A(:,1),'.')
plot(B(:,1)*0.001,B(:,2),'-r')
