function [PHI] = B_bar_matrix(A, B, Np, Nc)
  F = [];
  PHI = [];
  for j = 1:Nc
    for i = (1-j):(Np-j)

      if i < 0
        F = [F; 0*A^i*B];
      else
        F = [F; A^i*B];
      end

    end
    % Add to PHI
    PHI = [PHI F];
    % Clear F
    F = [];
  end
  PHI = cat(1,zeros(2,Np),PHI);
end