function [ x1,f1 ] = get_pdfs( x , parS, parW, parO )
     
%     tic;
    
    parS.gam    = x(1);
    parS.delta  = x(2);
    parS.zeta   = x(3);
    parS.lambda = x(4);


    [ f1, x1 ] = a2_WPI_function_opt( parS, parW, parO );

    
end

