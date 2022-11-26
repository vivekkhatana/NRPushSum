function PMatrix =  PColStochastic_outdeg(A, nodes)
    PMatrix = zeros(nodes,nodes);
    No = zeros(1,nodes);
   
    
    for i = 1:nodes
       for j = 1:nodes
           if A(j,i) ~=0
               No(i) = No(i) + 1;
           end
       end
    end
    
     
    for i = 1:nodes
        count = 0; 
        for j = 1:nodes
            if count < No(i)
              if A(j,i) ~=0
                PMatrix(j,i) = (1/No(i));
                count = count + 1;
              end
            end
            if count == No(i)
                PMatrix(j,i) = 1-sum(PMatrix(1:j-1,i));
            end
        end
    end
   
 end