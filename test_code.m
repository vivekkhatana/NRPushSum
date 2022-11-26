clear all

% load('graphArray_20Nodes_10000ConnectedDiGraph');
load('graphArray_100Nodes_100ConnectedDiGraph');
% load('graphArray_10Nodes_100ConnectedDiGraph');
Num_graphs = size(arr,1);
numberNodes = sqrt(size(arr,2));
numGraphs = 100;
numIniCond = 1000;

s = rng;
U = rand(numberNodes,numIniCond);
save('s');


for iniCond = 1%:numIniCond
        X0 = 10*U(:,iniCond);
        Y0 = 1*ones(numberNodes,1);
        mu = mean(X0);
        
        Diam = 6;
        numIterations = 200*Diam;
        tol = 1e-3;
       
        beta = 1e-4;
        
%         a = 0.4; b = 0.6; c = 0.3; d = 0.7;
        
        X = zeros(numberNodes,numIterations+1);
        X_a = zeros(numberNodes,numIterations+1);
  
        Y = zeros(numberNodes,numIterations+1);
        Y_a = zeros(numberNodes,numIterations+1);
%                 
        X(:,1) = X0;
        Y(:,1)= Y0;
        X_a(:,1) = X0;
        Y_a(:,1)= Y0;
        Z = zeros(numberNodes,numIterations+1);
        Z(:,1)= X0./Y0;
        
              
        numUpdates = ceil(numIterations/Diam)+1;
        MXP = zeros(numUpdates,1);
        MNP = zeros(numUpdates,1);
        mxp_l = Z(:,1);
        mnp_l = mxp_l;        
   
        for graphNo = 57%:700+numGraphs
            
             currentG = arr(graphNo,:);   
             currentG = reshape(currentG,numberNodes,numberNodes)'+eye(numberNodes);
             NumOutNeighbors = sum(currentG);
             inverseNumNeighbors = 1./NumOutNeighbors;
%              Weight_Matrix_cons = currentG*diag(inverseNumNeighbors,0);
             Weight_Matrix_cons = PDoubleStochastic(currentG, numberNodes);
             countMXP = 1;
            
            for i = 1:numIterations

             X_a(:,i+1) = (1-beta)*Weight_Matrix_cons*X_a(:,i) + (beta)*Weight_Matrix_cons*(X(:,i)) + X0  + Weight_Matrix_cons*normrnd(0,1,[numberNodes,1]);

             X(:,i+1) = beta*Weight_Matrix_cons*X_a(:,i+1)  +  beta*Weight_Matrix_cons*normrnd(0,1,[numberNodes,1]) ;

             Y_a(:,i+1) = (1-beta)*Weight_Matrix_cons*Y_a(:,i) + (beta)*Weight_Matrix_cons*(Y(:,i)) + Y0 + Weight_Matrix_cons*normrnd(0,1,[numberNodes,1]);

             Y(:,i+1) = beta*Weight_Matrix_cons*Y_a(:,i+1) + beta*Weight_Matrix_cons*normrnd(0,1,[numberNodes,1]) ;

             Z(:,i+1)= X(:,i+1)./Y(:,i+1);

%               
             M = [(1-beta)*Weight_Matrix_cons (1-beta)*beta*Weight_Matrix_cons^2; (beta)*Weight_Matrix_cons (beta)^2*Weight_Matrix_cons^2];
             
             eig(M);             

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
                       break;
                    else
%                        fprintf("Resetting MXP and MNP \n");
                       countMXP = countMXP+1;
                       mxp_l = Z(:,i+1);
                       mnp_l = mxp_l;
                    end
                end
         
            end
        end
end

figure(1);
plot(Z(:,1:i)')  
%  plot(Y(:,1:i)')  
% MXP_x=(1:countMXP)*Diam-Diam+1;
% hold on;stairs(MXP_x,MXP(1:countMXP),'-bx')
% MNP_x=(1:countMXP)*Diam-Diam+1;
% hold on;stairs(MNP_x,MNP(1:countMXP),'-ko')
hold off;
legend('Node 1','Node 2','Node 3','Node 4','Node 5','Node 6','Node 7','Node 8','Node 9','Node 10','Node 11','Node 12','Node 13','Node 14','Node 15','Node 16','Node 17','Node 18','Node 19','Node 20','MXP','MNP')


figure(2);
plot(X(:,1:i)') 
 
figure(3);
plot(Y(:,1:i)') 
  
