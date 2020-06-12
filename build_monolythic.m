function [K, F] = build_monolythic(Au, Av, T1, T2, G1, G2, B, bu, bv, bp, d)
    % Empty matrices of appropiate sizes
    zuv = sparse(size(Au,1), size(Av,1));
    zpp = sparse(size(G1,1), size(G1,1));
    zud = sparse(size(Au,1), size(B,1));
    zvd = sparse(size(Av,1), size(B,1));
    zpd = sparse(size(G1,1), size(B,1));
     
     K = [ Au    zuv     T1   zud
          zuv'    Av     T2   zvd
           G1     G2    zpp   zpd
          zud'   zvd'   zpd'    B  ];

    F = [ bu;     bv;    bp;   d];
end