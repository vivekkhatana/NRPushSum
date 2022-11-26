function PMatrix =  PRowStochastic(A, nodes)
    PMatrix = zeros(nodes,nodes);
   
   
for i = 1:nodes
    Ni = get_neighbors(A,i);
    for j = 1:length(Ni)
    PMatrix(i,Ni(j)) = 1/length(Ni);
    end
end
    
end