function val = objPS_Pf( x00 , parOO, parS, parW, parO, chn )
    lb = parOO.lb;
    ub = parOO.ub;
    Pf_acceptable = parOO.Pf_acceptable;
    x_failure     = parOO.x_failure;
    
    x = lb + (ub - lb).*x00;

    global computed
    idx = find( (computed(:,1)==x(1) & computed(:,2)==x(2) &  computed(:,3)==x(3)  &  computed(:,4)==x(4)), 1);
    
    if isempty(idx) 
        fprintf('a = %f, d = %f, z = %f, l = %f, ', x(1),x(2),x(3),x(4));
        if x(4) > 2*sqrt(x(2))
            val = 0;
            Pf = 10;
            shft = 0;
        else
            if x_failure > 0
                parW.targetPDF = 1;
                [ x1,f1 ] = get_pdfs( x , parS, parW, parO );
                
                [Pf,shft] = get_Pf(x1,f1,parOO);

                if (Pf - Pf_acceptable) >= 10e-7
                    val = 0;
                else
                    parW.targetPDF = 3;
                    [ x2,f2 ] = get_pdfs( x , parS, parW, parO );
                    val = get_Ph(x2,f2,parOO,x(1));
                end
            else
                parW.targetPDF = 3;
                [ x2,f2 ] = get_pdfs( x , parS, parW, parO );
                val = get_Ph(x2,f2,parOO,x(1));
                Pf = 0;
                shft = 0;
            end
        end
        fprintf('P = %f, Pf = %f, dx = %f \n', -val, Pf, shft);
        computed = [computed ;[x(1) x(2) x(3) x(4) val Pf shft]];
    else
        val = computed(idx,5);
        Pf  = computed(idx,6);
        shft  = computed(idx,7);
    end
    
    save(['computed_CHN_' num2str(chn) '_xf' num2str(parOO.x_failure) '.mat'],'computed');
    
end

