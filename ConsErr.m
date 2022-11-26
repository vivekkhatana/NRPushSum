function Consensus_Error_Traj_ret =  ConsErr(A, nodes, avg, iter_val)

[~,b] = size(A);
Consensus_Error_Traj = zeros(nodes,b);
   
for i = 1:nodes
  Consensus_Error_Traj(i,:) = (A(i,:) - avg).^2;
end
    
  Consensus_Error_Traj_ret =  Consensus_Error_Traj(iter_val);
end