function PMatrix =  P_Huang_Manton(A, nodes)

% This is the matrix used in COORDINATION AND CONSENSUS OF NETWORKED AGENTS WITH NOISY MEASUREMENTS:
% STOCHASTIC ALGORITHMS AND ASYMPTOTIC BEHAVIOR


PMatrix = zeros(nodes,nodes);
   
for i = 1:nodes
    Ni = get_neighbors(A,i);
    for j = 1:length(Ni)
        
        PMatrix(i,Ni(j)) = 1/(length(Ni)-1);
        
    end
end

for i = 1:nodes
    for j = 1:nodes
        if(i==j)
            PMatrix(i,j) = 0;
        end        
    end
end
    
end