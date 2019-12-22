function INITPs = choose_initP(fname,parOO,n_smpl,n_choo)

    load(fname)
    computed = zeros(size(Init_grid,1),7);
    for i = 1:size(Init_grid,1)
        if Init_grid(i,4) > 2*sqrt(Init_grid(i,2))
            val = 0;
            Pf = 10;
            shft = 0;
        else
            x1 = Init_x(i,:);
            x2 = Init_y(i,:);
            f1 = Init_pdfx(i,:);
            f2 = Init_pdfy(i,:);
            [Pf,shft] = get_Pf(x1,f1,parOO);
            if (Pf - parOO.Pf_acceptable) >= 10e-7 && parOO.x_failure > 0
                val = 0;
            else
                val = get_Ph(x2,f2,parOO,Init_grid(i,1));
            end            
        end
        computed(i,:) = [Init_grid(i,:),val,Pf,shft];
    end
    sorted_comp = sortrows(computed,5);
    s0 = sorted_comp(1:n_smpl,1:4);
    s = s0./max(s0,[],1);
    OUT = 2:n_smpl;
    IN  = [1];
    while length(IN) < n_choo
        no = length(OUT);
        ni = length(IN);
        dps = zeros(no,ni);
        for j = 1:no
            for k = 1:ni
                dps(j,k) = [0.1 0.4 0 0.5]*abs(s(OUT(j),:)-s(IN(k),:)).';
            end
        end
        score = min(dps,[],2);
        [ms,is] = max(score);
        IN = [IN,OUT(is)];
        OUT(is) = [];
    end
    disp(sorted_comp(IN,1:6))
    INITPs = sorted_comp(IN,:);
end

