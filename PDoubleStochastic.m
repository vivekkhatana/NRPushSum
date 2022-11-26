function PMatrix =  PDoubleStochastic(A, nodes)
    PMatrix = zeros(nodes,nodes);
    P = zeros(1,nodes);
   
for i = 1:nodes
    Ni = get_neighbors(A,i);
    for j = 1:length(Ni)
    P(i,Ni(j)) = 1/length(Ni);
    end
end
    
    
   for i = 1:nodes
     for j = 1:nodes
        if (i~=j)
            PMatrix(i,j) = min(P(i,j),P(j,i));
            PMatrix(j,i) = PMatrix(i,j);
        end
     end
     PMatrix(i,i) = 1-sum(PMatrix(i,:));
    end
    
end