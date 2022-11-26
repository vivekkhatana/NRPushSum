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
%         X0 = U(:,iniCond);
        X0 = [1,2,3,4,5,6,7,8,9,10]';
%         X0 = [1,2,3,4,5,-1,-2,-3,-4,-5]';
        Y0 = ones(numberNodes,1);
        mu = mean(X0);
        sum(X0);
        Diam = 5;
        numIterations= 5000*Diam;
        tol = 1e-4;

       
%         a = 0.4; b = 0.6; c = 0.3; d = 0.7;
        
        X = zeros(numberNodes,numIterations+1);
        X_a = zeros(numberNodes,numIterations+1);
  
        Y = zeros(numberNodes,numIterations+1);
        Y_a = zeros(numberNodes,numIterations+1);
%                 
        X(:,1) = X0;
        Y(:,1)= Y0;
%         X_a(:,1) = X0;
%         Y_a(:,1)= Y0;
        Z = zeros(numberNodes,numIterations+1);
        Z(:,1)= X0./Y0;
        Z_tilde = zeros(numberNodes,numIterations+1);
              
        numUpdates = ceil(numIterations/Diam)+1;
        MXP = zeros(numUpdates,1);
        MNP = zeros(numUpdates,1);
        mxp_l = Z(:,1);
        mnp_l = mxp_l;        
   
        for graphNo = 70%:700+numGraphs
            
             currentG = arr(graphNo,:);   
             currentG = reshape(currentG,numberNodes,numberNodes)'+eye(numberNodes);
             NumOutNeighbors = sum(currentG);
             inverseNumNeighbors = 1./NumOutNeighbors;
             Weight_Matrix_cons = currentG*diag(inverseNumNeighbors,0);
%              Weight_Matrix_cons = PDoubleStochastic(currentG, numberNodes);
             countMXP = 1;

            for i = 1:numIterations
                 
                if i > 200 % 500
%                    beta = 1/(i)^(1.5); 
                    beta = 100/(i)^(1.1);
%                      beta = 0.01;
%                      beta = 0;
                else
                    beta = 0.2;
%                     beta = 0.1;
                end

                 Beta1 = diag([0.9*beta*ones(3,1); 0.5*beta*ones(5,1); beta*ones(2,1)]);
                 Beta2 = [0.9*beta*ones(3,1); 0.5*beta*ones(5,1); beta*ones(2,1)]*ones(1,numberNodes);
                 Beta = (Beta2 - Beta1);   

                 Alpha = diag(Weight_Matrix_cons);
                 Alpha = (1 - Beta1*(1 - Alpha))./Alpha;
                 Alpha = diag(Alpha);
                 
                 sigma = 0;
                 noise_mean = 0.1;
                 
                 if i > 40 % 20 
%                      theta1 = 0.7^i;
                     theta1 = 1/(i)^(1.1);
                     Theta = [theta1*ones(numberNodes/2,1); theta1*ones(numberNodes/2,1)];
                 else
                     theta1 = 100;
                     Theta = [theta1*ones(numberNodes/2,1); theta1*ones(numberNodes/2,1)];
                 end
                 
%                   Theta = [theta1*ones(numberNodes/2,1); theta1*ones(numberNodes/2,1)];

%              X(:,i+1) = Alpha.*Weight_Matrix_cons*X(:,i) + Beta.*Weight_Matrix_cons*X(:,i) + (Theta).*X0 +  (Beta)*normrnd(noise_mean,sigma,[numberNodes,1]);
% 
%              Y(:,i+1) = Alpha.*Weight_Matrix_cons*Y(:,i) + Beta.*Weight_Matrix_cons*Y(:,i) + (Theta).*Y0  + (Beta)*normrnd(noise_mean,sigma,[numberNodes,1]);

             X(:,i+1) = Alpha.*Weight_Matrix_cons*X(:,i) + Beta'.*Weight_Matrix_cons*X(:,i) + (Theta).*X0 + (Beta)'*(-1 + 2*rand([numberNodes,1]));

             Y(:,i+1) = Alpha.*Weight_Matrix_cons*Y(:,i) + Beta'.*Weight_Matrix_cons*Y(:,i) + (Theta).*Y0 + (Beta)'*(-1 + 2*rand([numberNodes,1]));
             
             Z(:,i+1) = X(:,i+1)./Y(:,i+1);
%               
             M = Alpha.*Weight_Matrix_cons + Beta.*Weight_Matrix_cons;
             M1 = (Beta.*Weight_Matrix_cons)/beta;
%              abs(eig(M))
%              M = [a*Weight_Matrix_cons a*Weight_Matrix_cons; c*Weight_Matrix_cons d*Weight_Matrix_cons];
             

                 mxp_l=max(currentG*diag(mxp_l,0),[],2);
                 tempd=currentG*diag(mnp_l,0);
                 tempd(tempd==0)=NaN;
                 mnp_l=min(tempd,[],2);
                

                 if(mod(i,Diam)==0)
                    if(max(mxp_l)==min(mxp_l))
                       MXP(countMXP)= max(mxp_l);
                       else
                       fprintf("Error with mxp calculation\n");
                    end
                    if(max(mnp_l)==min(mnp_l))
                       MNP(countMXP)= min(mnp_l);
                    else
                       fprintf("Error  with mnp calculation\n");
                    end
                    if(MXP(countMXP)-MNP(countMXP) < tol)
                       fprintf("Converged with consensus value =%f \n",(MXP(countMXP)+MNP(countMXP))/2);
%                        break;
                    else
                       fprintf("Resetting MXP and MNP \n");
                       countMXP = countMXP+1;
                       mxp_l = Z(:,i+1);
                       mnp_l = mxp_l;
                    end
                end
         
            end
        end
end

% figure(1);
% plot(Z(:,1:i)')   
% MXP_x=(1:countMXP)*Diam-Diam+1;
% hold on;stairs(MXP_x,MXP(1:countMXP),'-bx')
% hold on;stairs(MXP_x,MNP(1:countMXP),'-ko')
% hold off;
% legend('Node 1','Node 2','Node 3','Node 4','Node 5','Node 6','Node 7','Node 8','Node 9','Node 10','Node 11','Node 12','Node 13','Node 14','Node 15','Node 16','Node 17','Node 18','Node 19','Node 20','MXP','MNP')

figure(1);
plot(Z(:,1:i)'); 
hold on
plot(mu*ones(length(X(:,1:i)'),1), '--', 'LineWidth', 1);
xlabel('Algorithm Iterations');
ylabel('Agent States');
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Agent 7','Agent 8','Agent 9','Agent 10', 'Initial Average');


figure(2); 
plot(X(:,1:i)');
xlabel('Algorithm Iterations');
ylabel('Numerator States');
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Agent 7','Agent 8','Agent 9','Agent 10', 'Initial Average');

 
figure(3);
plot(Y(:,1:i)');
xlabel('Algorithm Iterations');
ylabel('Denominator States');
legend('Agent 1','Agent 2','Agent 3','Agent 4','Agent 5','Agent 6','Agent 7','Agent 8','Agent 9','Agent 10', 'Initial Average');

  
