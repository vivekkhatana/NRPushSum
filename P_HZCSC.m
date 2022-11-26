function PMatrix =  P_HZCSC(A, nodes)

% This is the matrix used in Consensus Under Bounded Noise in Discrete Network Systems: An Algorithm With Fast
% Convergence and High Accuracy


PMatrix = zeros(nodes,nodes);
   
for i = 1:nodes
    Ni = get_neighbors(A,i);
    for j = 1:length(Ni)
        
        PMatrix(i,Ni(j)) = 1/length(Ni);
        
    end
end

for i = 1:nodes
    for j = 1:nodes
        if(i==j)
            PMatrix(i,j) = 1 - (sum(PMatrix(i,1:i-1)) + sum(PMatrix(i,i+1:end)));
        end        
    end
end
    
end