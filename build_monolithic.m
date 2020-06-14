function [K, F] = build_monolithic(Au, Av, T1, T2, G1, G2, B, bu, bv, bp, d)
    % Empty matrices of appropiate sizes
    uv0 = sparse(size(Au,1), size(Av,1));
    pp0 = sparse(size(G1,1), size(G1,1));
    ud0 = sparse(size(Au,1), size(B,1));
    vd0 = sparse(size(Av,1), size(B,1));
    pd0 = sparse(size(G1,1), size(B,1));
     
     K = [ Au    uv0     T1   ud0
          uv0'    Av     T2   vd0
           G1     G2    pp0   pd0
          ud0'   vd0'   pd0'    B  ];

    F = [ bu;     bv;    bp;   d];
end