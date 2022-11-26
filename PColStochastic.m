function PMatrix =  PColStochastic(A, nodes)
    PMatrix = zeros(nodes,nodes);
    outdeg = zeros(1,nodes);
   
    
    for i = 1:nodes
       for j = 1:nodes
           if A(j,i) ~=0
               outdeg(i) = outdeg(i) + 1;
           end
       end
    end
    
     
    for i = 1:nodes
        count = 0; 
        for j = 1:nodes
            if count < outdeg(i)
              if A(j,i) ~=0
                PMatrix(j,i) = (1/outdeg(i)).*rand(1,1);
                count = count + 1;
              end
            end
            if count == outdeg(i)
                PMatrix(j,i) = 1-sum(PMatrix(1:j-1,i));
            end
        end
    end
   
 end