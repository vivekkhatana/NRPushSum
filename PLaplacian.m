function PMatrix =  PLaplacian(A, nodes)
    PMatrix = zeros(nodes,nodes);
   
for i = 1:nodes
    for j = 1:nodes
        if (i~= j)
            PMatrix(i,j) = -A(i,j);
        end
    end
    PMatrix(i,i) = -sum(PMatrix(i,:));        
end
    
end