clear all;

% load('graphArray_20Nodes_10000ConnectedDiGraph');
% load('graphArray_100Nodes_100ConnectedDiGraph');
load('graphArray_10Nodes_100ConnectedDiGraph');
Num_graphs = size(arr,1);
numberNodes = sqrt(size(arr,2));
numGraphs = 100;
numIniCond = 1000;

s = rng;
U = 2*rand(numberNodes,numIniCond) - 1;
save('s');


for iniCond = 1%:numIniCond
        X0 = [1,2,3,4,5,6,7,8,9,10]';
        Y0 = ones(numberNodes,1);
        mu = mean(X0);
        sum(X0);
        Diam = 5;
        numIterations= 20000*Diam;
        tol = 1e-3;
        
        X_nrpush = zeros(numberNodes,numIterations+1);
        X_we = zeros(numberNodes,numIterations+1);
        X_hm = zeros(numberNodes,numIterations+1);
        X_tn = zeros(numberNodes,numIterations+1);
        X_hzcsc = zeros(numberNodes,numIterations+1);
        X_km = zeros(numberNodes,numIterations+1);
                  
        X_nrpush(:,1) = X0;
        X_we(:,1) = zeros(numberNodes,1);
        X_hm(:,1) = X0;
        X_tn(:,1) = X0;
        X_hzcsc(:,1) = X0;
        X_km(:,1) = X0;
        
        Y_nrpush = zeros(numberNodes,numIterations+1);
        Y_we = zeros(numberNodes,numIterations+1);
        
        Y_nrpush(:,1)= Y0;
        Y_we(:,1) = zeros(numberNodes,1);
        

        Z = zeros(numberNodes,numIterations+1);
        Z(:,1)= X0./Y0;
        
        SS = [];
        W = eye(numberNodes);
   
        for graphNo = 70%:700+numGraphs
            
             currentG = arr(graphNo,:);  
             currentG = reshape(currentG,numberNodes,numberNodes)'+eye(numberNodes);
             Laplacian = PLaplacian(currentG,numberNodes);
             NumOutNeighbors = sum(currentG);
             inverseNumNeighbors = 1./NumOutNeighbors;
             Weight_Matrix_cons = currentG*diag(inverseNumNeighbors,0);
             W_hm = P_Huang_Manton(currentG, numberNodes);
             W_tn = PDoubleStochastic(currentG, numberNodes);
             W_hzcsc = PDoubleStochastic(currentG, numberNodes);
       

            for i = 1:numIterations
                %% sequence initialization for the algorithms
                if i > 200 % 500
%                    beta = 1/(i)^(1.2); 
                    beta = 10/(i)^(1.1);
                    beta_hm = 1/(i)^(1.2);
                    beta_tn = 1/(i)^(0.99);
                    beta_km = 1/(i)^(0.9);
                    b_we = 0.04;
%                      beta = 0.01;
%                      beta = 0;
                else
                    beta = 0.2;
                    beta_hm = 1/(i)^(1.2);
                    beta_tn = 1/(i)^(0.99);
                    beta_km = 1/(i)^(0.9);
                    b_we = 0.04;
                end

                 Beta1 = diag([beta*ones(3,1); beta*ones(5,1); beta*ones(2,1)]);
                 Beta2 = [beta*ones(3,1); beta*ones(5,1); beta*ones(2,1)]*ones(1,numberNodes);
                 Beta = (Beta2 - Beta1);   

                 Alpha = diag(Weight_Matrix_cons);
                 Alpha = (1 - beta*(1 - Alpha))./Alpha;
                 Alpha = diag(Alpha);

                 if i > 1 % 20 
    %                      theta1 = 0.7^i;
                     theta1 = 1/(i)^(1.1);
                 else
                     theta1 = 100*(40+10*zeta(1.1));
                 end
                 
                Theta = [theta1*ones(numberNodes/2,1); theta1*ones(numberNodes/2,1)];
                
                %% Noise realizations
                
                noise_x = rand([numberNodes,1]);
                noise_y = rand([numberNodes,1]);
%                 noise_x = -1 + 2*rand([numberNodes,1]);
%                 noise_y = -1 + 2*rand([numberNodes,1]);
                
%                 if (mod(i,50) == 0)
%                      noise_x = -200 + 400*rand([numberNodes,1]);
%                      noise_y = -200 + 400*rand([numberNodes,1]);
%                 end
                
%                 sigma = 0.02;
%                 noise_mean = 0;
%                 noise_x = normrnd(noise_mean,sigma,[numberNodes,1]);
%                 noise_y = normrnd(noise_mean,sigma,[numberNodes,1]);
                
                %% NR-PushSum updates   
                X_nrpush(:,i+1) = Alpha.*Weight_Matrix_cons*X_nrpush(:,i) + Beta.*Weight_Matrix_cons*X_nrpush(:,i) + (Theta).*X0 +  (Beta)*noise_x;

                Y_nrpush(:,i+1) = Alpha.*Weight_Matrix_cons*Y_nrpush(:,i) + Beta.*Weight_Matrix_cons*Y_nrpush(:,i) + (Theta).*Y0  + (Beta)*noise_y;
             
                Z(:,i+1) = X_nrpush(:,i+1)./Y_nrpush(:,i+1);

                %% WE updates
                Y_we(:,i) = X0 + Laplacian*X_we(:,i) + noise_y;
                
                X_we(:,i+1) = X_we(:,i) - b_we*Laplacian*Y_we(:,i) + noise_x;
                
                %% HM updates
                X_hm(:,i+1) = (1 - beta_hm)*X_hm(:,i) + beta_hm*W_hm*( X_hm(:,i) + noise_x );
                
                %% TN updates
                X_tn(:,i+1) = ( (1 - beta_tn)*eye(numberNodes) + beta_tn*W_tn )*X_hm(:,i) + beta_tn*noise_x;
         
                %% HZCSC updates
                X_hzcsc(:,i+1) = W_hzcsc*( X_hm(:,i) + noise_x );
                
                %% KM updates
                X_km(:,i+1) = X_km(:,i) - beta_km*( Laplacian*X_km(:,i) + noise_x );
                
            end
        end
end

iter_val = 1:510:i;
% figure(1);
% plot(iter_val,Z(1,iter_val)', 'LineWidth', 2, 'MarkerEdgeColor','k'); 
% hold on
% plot(iter_val,Y_we(1,iter_val)', 'LineWidth', 2, 'MarkerEdgeColor','r');
% hold on
% plot(iter_val,X_hm(1,iter_val)', 'LineWidth', 2, 'MarkerEdgeColor','g');
% hold on
% plot(iter_val,X_tn(1,iter_val)', 'LineWidth', 2, 'MarkerEdgeColor','b');
% hold on
% plot(iter_val,X_hzcsc(1,iter_val)', 'LineWidth', 2, 'MarkerEdgeColor','m');
% hold on
% plot(iter_val,X_km(1,iter_val)', 'LineWidth', 2, 'MarkerEdgeColor','c');
% hold on
% plot(iter_val,mu*ones(length(X_nrpush(1,iter_val)'),1), '--', 'LineWidth', 2, 'MarkerEdgeColor','y');
% xlabel('Algorithm Iterations');
% ylabel('Agent 1 State Trajectory');
% legend('NR-PushSum','WE Algorithm','HM Algorithm','TN Algorithm','HZCSC Algorithm','KM Algorithm','Initial Average');
% 
% figure(2); 
% plot(iter_val,abs(0.1*sum(Z(:,iter_val)) - mu), 'LineWidth', 2, 'MarkerEdgeColor','k'); 
% hold on
% plot(iter_val,abs(0.1*sum(Y_we(:,iter_val)) - mu), 'LineWidth', 2, 'MarkerEdgeColor','r');
% hold on
% plot(iter_val,abs(0.1*sum(X_hm(:,iter_val)) - mu), 'LineWidth', 2, 'MarkerEdgeColor','g');
% hold on
% plot(iter_val,abs(0.1*sum(X_tn(:,iter_val)) - mu), 'LineWidth', 2, 'MarkerEdgeColor','b');
% hold on
% plot(iter_val,abs(0.1*sum(X_hzcsc(:,iter_val)) - mu), 'LineWidth', 2, 'MarkerEdgeColor','m');
% hold on
% plot(iter_val,abs(0.1*sum(X_km(:,iter_val)) - mu), 'LineWidth', 2, 'MarkerEdgeColor','c');
% xlabel('Algorithm Iterations');
% ylabel('Consensus residual');
% legend('NR-PushSum','WE Algorithm','HM Algorithm','TN Algorithm','HZCSC Algorithm','KM Algorithm');

figure(3); 
plot(iter_val,ConsErr(Z,numberNodes,mu,iter_val), 'LineWidth', 2, 'MarkerEdgeColor','k'); 
hold on
plot(iter_val,ConsErr(Y_we,numberNodes,mu,iter_val), 'LineWidth', 2, 'MarkerEdgeColor','r');
hold on
plot(iter_val,ConsErr(X_hm,numberNodes,mu,iter_val), 'LineWidth', 2, 'MarkerEdgeColor','g');
hold on
plot(iter_val,ConsErr(X_tn, numberNodes, mu,iter_val), 'LineWidth', 2, 'MarkerEdgeColor','b');
hold on
plot(iter_val,ConsErr(X_hzcsc, numberNodes, mu, iter_val), 'LineWidth', 2, 'MarkerEdgeColor','m');
hold on
plot(iter_val,ConsErr(X_km, numberNodes,mu,iter_val), 'LineWidth', 2, 'MarkerEdgeColor','c');
xlabel('Algorithm Iterations');
ylabel('Consensus residual');
legend('NR-PushSum','WE Algorithm','HM Algorithm','TN Algorithm','HZCSC Algorithm','KM Algorithm');

figure(4);
plot(Z(:,1:i)','LineWidth', 2); 
hold on
plot(mu*ones(length(Z(:,1:i)'),1), '--', 'LineWidth', 2);
xlabel('Algorithm Iterations');
ylabel('Agent States');
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Agent 7','Agent 8','Agent 9','Agent 10', 'Initial Average');

  
