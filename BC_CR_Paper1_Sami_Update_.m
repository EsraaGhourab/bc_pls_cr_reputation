clc; clear all ;close all;
%=========================================================================
%============================ Constants ==================================
%=========================================================================
% load('/Users/esraa/Documents/Career/downloads files/All Results/high_PU_load.mat')
%=========================================================================
%============================ Constants ==================================
%=========================================================================
osd = .8;
osr = 1;
ose = .2;
ord = 1;
ore = .1;
pout = 10^-6;
P_out = db2pow(pout)/2;
R = 7.5;
w = 1;Rc=1;
% No = db2pow(ys);
%=========================================================================

%translation to watt
ys = 20;  %SNR
No = P_out/ys;
P_No = P_out/No ;
M = 4;  % number of relays
L = 10;  % number of Eavesdropper
n = 1e1; % samples in channels
% Gsr =
% G= Csr / Grd ;
G1 = db2pow(0) ;
ofake = .2;
%=========================================================================
%=============================== Capacity and pe ===========================
%=========================================================================
SNR = 0:.2:10;
alfa = 0;
%=========================================================================
%=============================== Capacity and pe ===========================
%=========================================================================
MER = 0:1:30;
osd = .5;
osr = .1;
ord = 1;
ore = .5;

osd_su1 = .45;
osd_su2 =osd_su1;
opd = .2;
ofake = .5;
ose_su11 = .3;
ose_su22 = .3;
opd_su2 = opd;
%=========================================================================
%=============================== Capacity and pe ===========================
%=========================================================================
d = random('exp', 1.4,[1 11]) ;
load('/Users/esraaghourab/iCloud Drive (Archive) - 1/Documents/Career/master files/master/CR/sim_data/high_PU_load.mat')
% %=============================== Direct transmission =====================================
% %==========================================================================================
wCount = zeros(1,length(MER));wCount2 = zeros(1,length(MER));wCount3 = zeros(1,length(MER));
wCount_M = zeros(1,length(MER));wCount2_M = zeros(1,length(MER));wCount3_M = zeros(1,length(MER));
wcount = 0;Rcount = 0;wcount2 = 0;Rcount2 = 0;wcount3 = 0;Rcount3 = 0;
SU1 =  zeros(M,length(MER));SU2 =  SU1;SU3 =  SU1;
MU =  zeros(M,length(MER));MU2 =  MU;MU3 =  MU;
AU =  zeros(M,length(MER));AU2 =  AU;AU3 =  AU;
error = 0;Gc = 0;gc = 0;Gc2 = 0;gc2 = 0;Gc3 = 0;gc3 = 0;
%=========================================================================
% %=============================== Direct transmission =====================================
% %==========================================================================================
wcount_M = 0;Rcount_M = 0;wcount2_M = 0;Rcount2_M = 0;wcount3_M = 0;Rcount3_M = 0;
error_M = 0;Gc_M = 0;gc_M = 0;Gc2_M = 0;gc2_M = 0;Gc3_M = 0;gc3_M = 0;
RCount = zeros(1,length(MER));RCount2 = zeros(1,length(MER));RCount3 = zeros(1,length(MER));
RCount_M = zeros(1,length(MER));RCount2_M = zeros(1,length(MER));RCount3_M = zeros(1,length(MER));
SU1_M =  zeros(M,length(MER));SU2_M =  SU1_M;SU3_M =  SU1_M;
MU_M =  zeros(M,length(MER));MU2_M =  MU_M;MU3_M =  MU_M;
AU_M =  zeros(M,length(MER));AU2_M =  AU_M;AU3_M =  AU_M;
%%
for i = 2 : length(MER)
    mer = MER(i);
    mer = db2pow(mer);
    hs = random('exp', osr,[M n]) ;  %'rayl' , 'exp'   % no sqrt here
    hsd = random('exp', osd,[1 n]) ;
    hh_i = random('exp', .04,[1 n]) ;
    %---------------------------------------- SU class 1 data ---------------------------------------------------------------------------
    %=======================================================================
    Snesing_Matrix1 =  randi([0 1], [3 1e3]);  % spectrum sensing matrix
    Snesing_SU1_Matrix =  (SU1_matrix(:,1:3))';  % spectrum sensing matrix
    Snesing_SU1_Matrix = [Snesing_SU1_Matrix Snesing_Matrix1];
    Channel_su1 = zeros(size(Snesing_SU1_Matrix)); %el available channel el 3ndi
    nChannel_su1 = size(Snesing_SU1_Matrix, 1); % number of available channel
    Tmax_su1 = size(Snesing_SU1_Matrix, 2); %maximum time instant available
    f_su1 =nChannel_su1;
    t_su1 = Tmax_su1;
    n = t_su1;
    item = 1;
    for time = 1:Tmax_su1 % counter 3la el time lkol time instant

        %         lfsr = [xor(xor(xor(lfsr(16), lfsr(14)), lfsr(13)),lfsr(11)) lfsr(1:15)];
        timeSample_su1 = Snesing_SU1_Matrix(:, time);

        for ch = 1:nChannel_su1

            if(timeSample_su1(ch) == 1)
                transmitted_su1(item) =  random('exp', osd_su1 ,[1 1])  ; % real data
                Channel_su1(ch, time) = transmitted_su1(item) ;
                item = item+1;
            else
                Channel_su1(ch, time) =random('exp', ofake ,[1 1]) ; % fake data
            end
        end
    end
    %------------------------------------------------------------------------------------  Variance of SU1 to Eav and PU to Eav  --------------------------------------------------------------------------------------
    ose = osd/mer;
    ope = opd/mer;
    %------------------------------------------------------------------------------  Fading channel cooeffitionts for  SU1 with the AWGN -----------------------------------------------------------------------
    Channel_su1 = Channel_su1;
    Channel_su22 = random('exp', .2 ,[f_su1 n]);
    %------------------------------------------------------------------------------- rayleigh Fading channel cooeffitionts for PU and SU1 -----------------------------------------------------------------------
    %a3ml for loop 3la f 34an a2dr a3ml generation lkol el channels general
    hsd1 =  (Channel_su1(1, :));
    hsd2 =  (Channel_su1(2, :));
    hsd3 =  (Channel_su1(3, :));

    hd = random('exp', ord,[M length(hsd1)]) ;

    arg = min(hsd1,hd) ; % no square here
    arg2 = min(hsd2,hd) ; % no square here
    arg3 = min(hsd3,hd) ; % no square here
    %     arg = (hs.*hd)./(hs+hd) ; % no square here
    hsrd = max(arg);
    hsrd2 = 1.2 .* max(arg2);
    hsrd3 = 2.2 .* max(arg3);

    ose = osd/mer;
    hse = min(random('exp', ose,[L n]));

    alpha = ((2 ^ R) - 1) .\ mer;
    csd = log2( 2 + P_No .* hsrd );
    csd2 = log2( 2 + P_No .* hsrd2 );
    csd3 = log2( 2 + P_No .* hsrd3 );

    %     hse = random('exp', ose, [1 length(hsd1)]);
    hpe = random('exp', ope, [1 length(hsd1)]);
    cse = ( log2(1+ ( (hse * P_No  ) ./  ( ( hpe * P_No) +1 ))) ) ;

    cs = csd - cse;
    cs2 = csd2 - cse;
    cs3 = csd3 - cse;

    Cs_alfa(i) = 3 .* mean( max(cs,0)); % max and zero
    Cs2_alfa(i) = 3.2 .* mean( max(cs2,0)); % max and zero
    Cs3_alfa(i) = 3.42.* mean( max(cs3,0)); % max and zero

    for q = 1 : M
        k = i;
        if cs(i) > 0
            W_1_SU(i) = Cs_alfa(i);
        else
            W_1_SU(i) = 0;
        end

        if cs2(k) > 0
            W_2_SU(k) = Cs2_alfa(i);
        else
            W_2_SU(k) = 0;
        end

        if cs3(k) > 0
            W_3_SU(k) = Cs3_alfa(i);
        else
            W_3_SU(k) = 0;
        end

        x = round(max(Cs_alfa));
        Cost = randi([1,x],1,length(MER));

        if W_1_SU(k) >= Cost(k)
            SU1(k) = SU1(k-1) + 1;
            W_1_SU_(k+1) = W_1_SU(k) - Cost(k);
        elseif W_2_SU(k) >= Cost(k)
            SU2(k) = SU2(k-1) + 1;
            W_2_SU_(k+1) = W_2_SU(k) - Cost(k);
        else
            SU3(k) = SU3(k-1) + 1;
            W_3_SU(k+1) = W_3_SU(k) - Cost(k);
        end

        SU1_prob = sum(SU1)./length(SU1);
        SU2_prob = sum(SU2)./length(SU1);
        SU3_prob = sum(SU3)./length(SU1);


        Prob_int_direct(i) = sum(cs < 0)./length(cs);
        Prob2_int_direct(i) = sum(cs2 < 0)./length(cs);
        Prob3_int_direct(i) = sum(cs3 < 0)./length(cs);
        BER(i,q) = Prob_int_direct(i);

        z = 0.0009 * rand(1,length(Prob_int_direct));
        Prob_int_direct = z + Prob_int_direct;Prob2_int_direct = z + Prob2_int_direct;Prob3_int_direct = z + Prob3_int_direct;

        thr = (max(Prob_int_direct) + min(Prob_int_direct))./2;
        thr2 = (max(Prob2_int_direct) + min(Prob2_int_direct))./2;
        thr3 = (max(Prob3_int_direct) + min(Prob3_int_direct))./2;

        Rate(i) = 1 - Prob_int_direct(i)
        Rate2(i) = 1 - Prob2_int_direct(i)
        Rate3(i) = 1 - Prob3_int_direct(i)

        if SU1(k) > 0 && Rate (i) > thr
            Rate_(i) = Rate(i-1) +1
            W_1_SU(k+1) = W_1_SU_(i) + W_1_SU(k) + 1;
        else
            Rate_(i) = Rate(i-1);
        end

        if SU2(k) > 0 && Rate2 (i) > thr2
            Rate2_(i) = Rate2(i-1) +1
            W_2_SU(k+1) = W_2_SU_(k)+ W_2_SU(i) + 1;
        else
            Rate2_(i) = Rate2(i-1);
        end

        if SU3(k) > 0 && Rate3 (i) > thr3
            Rate3_(i) = Rate3(i-1) +1
            W_3_SU(k+1) = W_3_SU(k) + W_2_SU_(i) + 1;
        else
            Rate3_(i) = Rate3(i-1);
        end

        R_best(i) = (max(Rate_(i)))
        R_best2(i) = (max(Rate2_(i)))
        R_best3(i) = (max(Rate3_(i)))

        R_best_(i) = (max(Rate(i)))
        R_best2_(i) = (max(Rate2(i)))
        R_best3_(i) = (max(Rate3(i)))

        thrr = 2.6;
        x1(i) = (sum((Rate_(i)>thrr)));
        x2(i) = (sum((Rate_(i)<thrr)));

        x1_2(i) = (sum((Rate2_(i)>thrr)));
        x2_2(i) = (sum((Rate2_(i)<thrr)));

        x1_3(i) = (sum((Rate3_(i)>thrr)));
        x2_3(i) = (sum((Rate3_(i)<thrr)));

        G1(i) = mean(max (Rate_(i) , thrr));
        G2(i) = mean(min (Rate_(i) , thrr));

        G1_2(i) = mean(max (Rate2_(i) , thrr));
        G2_2(i) = mean(min (Rate2_(i) , thrr));

        G1_3(i) = mean(max (Rate3_(i) , thrr));
        G2_3(i) = mean(min (Rate3_(i) , thrr));

                if Rate(i) > thrr
                    Gc = Gc+1;
                    Gc_(q,i) = Gc;
                else
                    gc = gc + 1;
                    gc_(q,i) = gc;
                end
                Gt(q,i) = G1(i);
                gt(q,i) = G2(i);

        Gg(i) = G1(i) + G2(i);
        Ggt(q,i) = Gg(i);

        xg = ((x1(i)) .* G1) + ((x2(i)) .* G2);
        Xgt(q,i) = xg(i);

        Trust_value(i) = sum (xg(i) ./ Gg(i));
        Trust_value(i) = Trust_value(i)./length(Trust_value);

        if W_1_SU(k) > W_1_SU(k-1) && Rate(k) > Rate(k-1)
            Trust_value(i) = Trust_value(i-1) +2;
        elseif W_1_SU(k) > W_1_SU(k-1) || Rate(k) > Rate(k-1)
            Trust_value(i) = Trust_value(i-1) + 1;
        else
            Trust_value(i) = Trust_value(i-1) -1;
        end

        if Trust_value(i) < 0
            Trust_value(i) = Trust_value(i-1)
        end
        Trust_value_(i) = Trust_value(i) ./ length(Trust_value) ;

        Thr = 0.3;
        if Trust_value_(i) > Thr
            AU(i) = AU(i-1) + 1;
        else
            MU(i) = MU(i-1) + 1;
        end
        PD = sum(AU) ./ length(MER);
        PMD = sum(MU) ./ length(MER);

                if Rate2(i) > thrr
                    Gc2 = Gc2+1;
                    Gc2(q,i) = Gc2;
                else
                    gc2 = gc2 + 1;
                    gc_2(q,i) = gc2;
                end
                Gt2(q,i) = G1_2(i);
                gt2(q,i) = G2_2(i);

        Gg2(i) = G1_2(i) + G2_2(i);
        Ggt2(q,i) = Gg2(i);

        xg2 = ((x1_2(i)) .* G1_2) + ((x2_2(i)) .* G2_2);
        Xgt2(q,i) = xg2(i);

        Trust_value2(i) = sum (xg2(i) ./ Gg2(i)) ;
        Trust_value2(i) = Trust_value2(i) ./ length(Trust_value2) ;

        if W_2_SU(k) > W_2_SU(k-1) && Rate2(k) > Rate2(k-1)
            Trust_value2(i) = Trust_value2(i-1) +2;
        elseif W_2_SU(k) > W_2_SU(k-1) || Rate2(k) > Rate2(k-1)
            Trust_value2(i) = Trust_value2(i-1) + 1;
        else
            Trust_value2(i) = Trust_value2(i-1) -1;
        end

        if Trust_value2(i) < 0
            Trust_value2(i) = Trust_value2(i-1);
        end

        Trust_value2_(i) = Trust_value2(i) ./ length(Trust_value2) ;

        if Trust_value2_(i) > Thr
            AU2(i) = AU2(i-1) + 1;
        else
            MU2(i) = MU2(i-1) + 1;
        end
        PD2 = sum(AU2) ./ length(MER);
        PMD2 = sum(MU2) ./ length(MER);
        
                if Rate3(i) > thrr
                    Gc3 = Gc3+1;
                    Gc_3(q,i) = Gc3;
                else
                    gc3 = gc3 + 1;
                    gc_3(q,i) = gc3;
                end
                Gt3(q,i) = G1_3(i);
                gt3(q,i) = G2_2(i);

        Gg3(i) = G1_3(i) + G2_3(i);
        Ggt3(q,i) = Gg3(i);

        xg3 = ((x1_3(i)) .* G1_3) + ((x2_3(i)) .* G2_3);
        Xgt3(q,i) = xg3(i);

        Trust_value3(i) = sum (xg3(i) ./ Gg3(i)) ;
        Trust_value3(i) = Trust_value3(i) ./ length(Trust_value3) ;

        if W_3_SU(k) > W_3_SU(k-1) && Rate3(k) > Rate3(k-1)
            Trust_value3(i) = Trust_value3(i-1) +2;
        elseif W_3_SU(k) > W_3_SU(k-1) || Rate3(k) > Rate3(k-1)
            Trust_value3(i) = Trust_value3(i-1) + 1;
        else
            Trust_value3(i) = Trust_value3(i-1) - 1;
        end

        if Trust_value3(i) < 0
            Trust_value3(i) = Trust_value3(i-1)
        end

        Trust_value3_(i) = Trust_value3(i) ./ length(Trust_value3) ;

        if Trust_value3_(i) > Thr
            AU3(i) = AU3(i-1) + 1;
        else
            MU3(i) = MU3(i-1) + 1;
        end
        %             end
        PD3 = sum(AU3) ./ length(MER);
        PMD3 = sum(MU3) ./ length(MER);

    end

    TRV(q,i) = (Trust_value(i));
    TRV2(q,i) = (Trust_value2(i));
    TRV3(q,i) = (Trust_value3(i));
    TRVt = Xgt ./ Ggt;

end
Pd = AU / MER
Pd2 = AU2 / MER
Pd3 = AU3 / MER
Pmd = MU / MER
Pmd2 = MU2 / MER
Pmd3 = MU3 / MER

Cs_tot = [Cs_alfa; Cs2_alfa; Cs3_alfa];
Pint_tot = [Prob_int_direct;Prob2_int_direct;Prob3_int_direct];
Rate_tot = [Rate; Rate2;Rate3];
R_best = (max((Rate_tot)));
Cs_Best = max((Cs_tot));
Pint_best = min((Pint_tot));
%=========================================================================
%=========================================================================
%%
for i = 2 : length(MER)
    mer = MER(i);
    mer = db2pow(mer);
    %---------------------------------------- SU class 1 data ---------------------------------------------------------------------------
    %=======================================================================
    Snesing_Matrix1_M =  randi([0 1], [3 10]);  % spectrum sensing matrix
    Snesing_SU1_Matrix_M =  (SU2_matrix(1:length(MER),2:4))';  % spectrum sensing matrix
    Snesing_SU1_Matrix_M = [Snesing_SU1_Matrix_M Snesing_Matrix1_M];
    Channel_su1_M = zeros(size(Snesing_SU1_Matrix_M)); %el available channel el 3ndi
    nChannel_su1_M = size(Snesing_SU1_Matrix_M, 1); % number of available channel
    Tmax_su1_M = size(Snesing_SU1_Matrix_M, 2); %maximum time instant available
    f_su1_M =nChannel_su1_M;
    t_su1_M = Tmax_su1_M;
    n = t_su1_M;
    item = 1;
    for time_M = 1:Tmax_su1_M % counter 3la el time lkol time instant
        timeSample_su1_M = Snesing_SU1_Matrix_M(:, time_M);
        for ch_M = 1:nChannel_su1_M
            if(timeSample_su1_M(ch_M) == 1)
                transmitted_su1_M(item) =  random('exp', osd ,[1 1])  ; % real data
                Channel_su1_M(ch_M, time_M) = transmitted_su1_M(item) ;
                item = item+1;
            else
                Channel_su1_M(ch_M, time_M) =random('exp', ofake ,[1 1]) ; % fake data
            end
        end
    end
    %------------------------------------------------------------------------------------  Variance of SU1 to Eav and PU to Eav  --------------------------------------------------------------------------------------
    ose = osd./mer;
    ope = opd./mer;
    %------------------------------------------------------------------------------  Fading channel cooeffitionts for  SU1 with the AWGN -----------------------------------------------------------------------
    Channel_su1_M = Channel_su1_M;
    hpd_M = random('exp', .3 ,[1 n]);
    %------------------------------------------------------------------------------- rayleigh Fading channel cooeffitionts for PU and SU1 -----------------------------------------------------------------------
    %a3ml for loop 3la f 34an a2dr a3ml generation lkol el channels general
    hsd1_M =  (Channel_su1_M(1, :));
    hsd2_M =  (Channel_su1_M(2, :));
    hsd3_M =  (Channel_su1_M(3, :));

    hd_M = random('exp', ord,[M length(hsd1_M)]) ;
    hse_M = max(random('exp', ose,[M n])) ;

    arg_M = min(hsd1_M,hd_M) ; % no square here
    arg2_M = min(hsd2_M,hd_M) ; % no square here
    arg3_M = min(hsd3_M,hd_M) ; % no square here

    hsrd_M = max(arg_M);
    hsrd2_M = max(arg2_M);
    hsrd3_M = max(arg3_M);

    alpha_M = ((2 ^ R) - 1) .\ mer;
    csd_M = ( log2(1+( (  hsrd_M .* P_No  ) ./  ( ( hsrd_M .* P_No + hpd_M .* P_No) +1 ))) );
    csd2_M = ( log2(1+( ( hsrd2_M .* P_No  ) ./  ( ( hsrd2_M .* P_No + hpd_M .* P_No) +1 ))) );
    csd3_M = ( log2(1+( ( hsrd3_M .* P_No  ) ./  ( ( hsrd3_M .* P_No + hpd_M .* P_No) +1 ))) );

    hse_ = min(random('exp', ose, [M length(hsd1_M)]));

    hpe_M = random('exp', ope, [1 length(hsd1_M)]);
    cse_M = ( log2(1+( ( hse_ * (2 .* P_No)) ./  ( ( hpe_M .* P_No) +1 ))) );
    cse2_M = ( log2(1+( ( hse_ * (2 .* P_No)) ./  ( ( hpe_M .* P_No) +1 ))) );
    cse3_M = ( log2(1+( ( hse_ * (2 .* P_No)) ./  ( ( hpe_M .* P_No) +1 ))) );

    cs_M = csd_M - cse_M;
    cs2_M = csd2_M - cse2_M;
    cs3_M = csd3_M - cse3_M;

    Cs_alfa_M(i) = 3 .* mean( max(cs_M,0)); % max and zero
    Cs2_alfa_M(i) = 3.2 .* mean( max(cs2_M,0)); % max and zero
    Cs3_alfa_M(i) = 3.42.* mean( max(cs3_M,0)); % max and zero

    for q = 1 : M

        k = i;
        if cs_M(i) > 0
            W_1_SU_M(i) = Cs_alfa_M(i);
        else
            W_1_SU_M(i) = 0;
        end

        if cs2_M(k) > 0
            W_2_SU_M(k) = Cs2_alfa_M(i);
        else
            W_2_SU_M(k) = 0;
        end

        if cs3_M(k) > 0
            W_3_SU_M(k) = Cs3_alfa_M(i);
        else
            W_3_SU_M(k) = 0;
        end

        x_M = round(max(Cs_alfa_M));
        Cost_M = randi([0,x_M],1,length(MER));

        if W_1_SU_M(k) >= Cost_M(k)
            SU1_M(k) = SU1_M(k-1) + 1;
            W_1_SU_M_(k+1) = W_1_SU_M(k) - Cost_M(k);
        elseif W_2_SU_M(k) >= Cost_M(k)
            SU2_M(k) = SU2_M(k-1) + 1;
            W_2_SU_M_(k+1) = W_2_SU_M(k) - Cost_M(k);
        else
            SU3_M(k) = SU3_M(k-1) + 1;
            W_3_SU_M_(k+1) = W_3_SU_M(k) - Cost_M(k);
        end

        SU1_prob_M = sum(SU1_M)./length(SU1_M);
        SU2_prob_M = sum(SU2_M)./length(SU1_M);
        SU3_prob_M = sum(SU3_M)./length(SU1_M);

        Prob_int_direct_M(i) = sum(cs_M < 0)./length(cs_M);
        Prob2_int_direct_M(i) = sum(cs2_M < 0)./length(cs2_M);
        Prob3_int_direct_M(i) = sum(cs3_M < 0)./length(cs3_M);
        BER_(i,q) = Prob_int_direct_M(i);

        Prob_int_direct_M(i) = z(i) + Prob_int_direct_M(i);Prob2_int_direct_M(i) = z(i) + Prob2_int_direct_M(i);Prob3_int_direct_M(i) = z(i) + Prob3_int_direct_M(i);

        thr_M = (max(Prob_int_direct_M) + min(Prob_int_direct_M))./2;
        thr2_M = (max(Prob2_int_direct_M) + min(Prob2_int_direct_M))./2;
        thr3_M = (max(Prob3_int_direct_M) + min(Prob3_int_direct_M))./2;

        Rate_M(i) = 1 - (Prob_int_direct_M(i))
        Rate2_M(i) = 1 - (Prob2_int_direct_M(i))
        Rate3_M(i) = 1 - (Prob3_int_direct_M(i))

        if SU1_M(k) > 0 && Rate_M (i) > thr_M
            Rate_M_(i) = Rate_M(i-1) +1
            W_1_SU_M(k+1) = W_1_SU_M(k) + W_1_SU_M_(k) + 1;
        else
            Rate_M_(i) = Rate_M(i-1);
        end

        if SU2_M(k) > 0 && Rate2_M (i) > thr2_M
            Rate2_M_(i) = Rate2_M(i-1) +1
            W_2_SU_M(k+1) = W_2_SU_M(k) + W_2_SU_M_(k) + 1;
        else
            Rate2_M_(i) = Rate2_M(i-1);
        end

        if SU3_M(k) > 0 && Rate3_M (i) > thr3_M
            Rate3_M_(i) = Rate3_M(i-1) +1
            W_3_SU_M(k+1) = W_3_SU_M(k) + W_3_SU_M_(k) + 1;
        else
            Rate3_M_(i) = Rate3_M(i-1);
        end

        R_best_M(i) = (max(Rate_M_(i)))
        R_best2_M(i) = (max(Rate2_M_(i)))
        R_best3_M(i) = (max(Rate3_M_(i)))

        R_best_M_(i) = (max(Rate_M(i)))
        R_best2_M_(i) = (max(Rate2_M(i)))
        R_best3_M_(i) = (max(Rate3_M(i)))

        thrr_M = 2.6;
        x1_M(i) = (sum((Rate_M_(i)>thrr_M)));
        x2_M(i) = (sum((Rate_M_(i)<thrr_M)));

        x1_2_M(i) = (sum((Rate2_M_(i)>thrr_M)));
        x2_2_M(i) = (sum((Rate2_M_(i)<thrr_M)));

        x1_3_M(i) = (sum((Rate3_M_(i)>thrr_M)));
        x2_3_M(i) = (sum((Rate3_M_(i)<thrr_M)));

        G1_M(i) = mean(max (Rate_M_(i) , thrr_M));
        G2_M(i) = mean(min (Rate_M_(i) , thrr_M));

        G1_2_M(i) = mean(max (Rate2_M_(i) , thrr_M));
        G2_2_M(i) = mean(min (Rate2_M_(i) , thrr_M));

        G1_3_M(i) = mean(max (Rate3_M_(i) , thrr_M));
        G2_3_M(i) = mean(min (Rate3_M_(i) , thrr_M));

                if Rate_M(i) > thrr_M
                    Gc_M = Gc_M+1;
                    Gc__M(q,i) = Gc_M;
                else
                    gc_M = gc_M + 1;
                    gc__M(q,i) = gc_M;
                end
                Gt_M(q,i) = G1_M(i);
                gt_M(q,i) = G2_M(i);

        Gg_M(i) = G1_M(i) + G2_M(i);
        Ggt_M(q,i) = Gg_M(i);

        xg_M = ((x1_M(i)) .* G1_M) + ((x2_M(i)) .* G2_M);
        Xgt_M(q,i) = xg_M(i);

        Trust_value_M(i) = sum (xg_M(i) ./ Gg_M(i)) ;
        Trust_value_M(i) = Trust_value_M(i) ./length(Trust_value_M) ;

        if W_1_SU_M(k) > W_1_SU_M(k-1) && Rate_M(k) > Rate_M(k-1)
            Trust_value_M(i) = Trust_value_M(i-1) +2;
        elseif W_1_SU_M(k) > W_1_SU_M(k-1) || Rate_M(k) > Rate_M(k-1)
            Trust_value_M(i) = Trust_value_M(i-1) + 1;
        else
            Trust_value_M(i) = Trust_value_M(i-1) - 1;
        end

        if Trust_value_M(i) < 0
            Trust_value_M(i) = Trust_value_M(i-1)
        end
        Trust_value_M_(i) = Trust_value_M(i)./length(Trust_value_M);

        Thr_M = 0.6;
        if Trust_value_M_(i) > Thr_M
            AU_M(i) = AU_M(i-1) + 1;
        else
            MU_M(i) = MU_M(i-1) + 1;
        end

        PD_M = sum(AU_M) ./ length(MER);
        PMD_M = sum(MU_M) ./ length(MER);

                if Rate2_M(i) > thrr_M
                    Gc2_M = Gc2_M+1;
                    Gc2_M(q,i) = Gc2_M;
                else
                    gc2_M = gc2_M + 1;
                    gc_2M(q,i) = gc2_M;
                end
                Gt2_M(q,i) = G1_2_M(i);
                gt2_M(q,i) = G2_2_M(i);

        Gg2_M(i) = G1_2_M(i) + G2_2_M(i);
        Ggt2_M(q,i) = Gg2_M(i);

        xg2_M = ((x1_2_M(i)) .* G1_2_M) + ((x2_2_M(i)) .* G2_2_M);
        Xgt2_M(q,i) = xg2_M(i);

        Trust_value2_M(i) = sum (xg2_M(i) ./ Gg2_M(i)) ;
        Trust_value2_M(i) = Trust_value2_M(i)./length(Trust_value2_M);

        if W_2_SU_M(k) > W_2_SU_M(k-1) && Rate2_M(k) > Rate2_M(k-1)
            Trust_value2_M(i) = Trust_value2_M(i-1) +2;
        elseif W_2_SU_M(k) > W_2_SU_M(k-1) || Rate2_M(k) > Rate2_M(k-1)
            Trust_value2_M(i) = Trust_value2_M(i-1) + 1;
        else
            Trust_value2_M(i) = Trust_value2_M(i-1) - 1;
        end

        if Trust_value2_M(i) < 0
            Trust_value2_M(i) = Trust_value2_M(i-1)
        end

        Trust_value2_M_(i) = Trust_value2_M(i)./length(Trust_value2_M);

        if Trust_value2_M_(i) > Thr_M
            AU2_M(i) = AU2_M(i-1) + 1;
        else
            MU2_M(i) = MU2_M(i-1) + 1;
        end
        PD2_M = sum(AU2_M) ./ length(MER);
        PMD2_M = sum(MU2_M) ./ length(MER);

                if Rate3_M(i) > thrr_M
                    Gc3_M = Gc3_M+1;
                    Gc_3M(q,i) = Gc3_M;
                else
                    gc3_M = gc3_M + 1;
                    gc_3M(q,i) = gc3_M;
                end
                Gt3_M(q,i) = G1_3_M(i);
                gt3_M(q,i) = G2_2_M(i);
        
        Gg3_M(i) = G1_3_M(i) + G2_3_M(i);
        Ggt3M(q,i) = Gg3_M(i);

        xg3_M = ((x1_3_M(i)) .* G1_3_M) + ((x2_3_M(i)) .* G2_3_M);
        Xgt3_M(q,i) = xg3_M(i);

        Trust_value3_M(i) = sum (xg3_M(i) ./ Gg3_M(i)) ;
        Trust_value3_M(i) = Trust_value3_M(i)./length(Trust_value3_M);

        if W_3_SU_M(k) > W_3_SU_M(k-1) && Rate3_M(k) > Rate3_M(k-1)
            Trust_value3_M(i) = Trust_value3_M(i-1) +2;
        elseif W_3_SU_M(k) > W_3_SU_M(k-1) || Rate3_M(k) > Rate3_M(k-1)
            Trust_value3_M(i) = Trust_value3_M(i-1) + 1;
        else
            Trust_value3_M(i) = Trust_value3_M(i-1) - 1;
        end

        if Trust_value3_M(i) < 0
            Trust_value3_M(i) = Trust_value3_M(i-1)
        end

        Trust_value3_M_(i) = Trust_value3_M(i)./length(Trust_value3_M);

        if Trust_value3_M_(i) > Thr_M
            AU3_M(i) = AU3_M(i-1) + 1;
        else
            MU3_M(i) = MU3_M(i-1) + 1;
        end
        %             end
        PD3_M = sum(AU3_M) ./ length(MER);
        PMD3_M = sum(MU3_M) ./ length(MER);

    end

    TRV_M(q,i) = (Trust_value_M(i));
    TRV2_M(q,i) = (Trust_value2_M(i));
    TRV3_M(q,i) = (Trust_value3_M(i));
    TRVt_M = Xgt_M ./ Ggt_M;

end
Pd_M = AU_M / MER
Pd2_M = AU2_M / MER
Pd3_M = AU3_M / MER
Pmd_M = MU_M / MER
Pmd2_M = MU2_M / MER
Pmd3_M = MU3_M / MER

Cs_tot_M = [Cs_alfa_M; Cs2_alfa_M; Cs3_alfa_M];
Pint_tot_M = [Prob_int_direct_M;Prob2_int_direct_M;Prob3_int_direct_M];
Rate_tot_M = [Rate_M; Rate2_M;Rate3_M];
R_best_M = (max((Rate_tot_M)));
Cs_Best_M = max((Cs_tot_M));
Pint_best_M = min((Pint_tot_M));

%=========================================================================
%=========================================================================
%%
% for q = 1 : M
for j = 1: length(Trust_value)
    Thr = .6;
    if Trust_value_(j) < Thr
        wcount = wcount + 1;
        wCount(j) = wcount;

        if wCount(j) == 3
            Rcount = Rcount + 1;
            RCount(j) = Rcount;
            wcount = 0;
            wCount(j) = 0;
            if RCount(j) == 2
                Rcount = 0;
                RCount(j) = 0;
            end
        end
    end
end
%%
for j = 1: length(Trust_value)
    Thr = 0.6;
    if Trust_value2_(j) < Thr
        wcount2 = wcount2 + 1;
        wCount2(j) = wcount2;
        if wCount2(j) == 3
            Rcount2 = Rcount2 + 1;
            RCount2(j) = Rcount2;
            wcount2 = 0;
            wCount2(j) = 0;
            if RCount2(j) == 2
                Rcount2 = 0;
                RCount2(j) = 0;
            end
        end
    end
end
%%
for j = 1: length(Trust_value)
    Thr = 0.6;
    if Trust_value3_(j) < Thr
        wcount3 = wcount3 + 1;
        wCount3(j) = wcount3;
        if wCount3(j) == 3
            Rcount3 = Rcount3 + 1;
            RCount3(j) = Rcount3;
            wcount3 = 0;
            wCount3(j) = 0;
            if RCount3(j) == 2
                Rcount3 = 0;
                RCount3(j) = 0;
            end
        end
    end
end
%%
%         else
%             Trust_value(i) = 5;
% for q = 1 : M
for j = 1: length(Trust_value_M)
    Thr_M = 0.6;
    if Trust_value_M_(j) < Thr_M
        wcount_M = wcount_M + 1;
        wCount_M(j) = wcount_M;

        if wCount_M(j) == 3
            Rcount_M = Rcount_M + 1;
            RCount_M(j) = Rcount_M;
            wcount_M = 0;
            wCount_M(j) = 0;
            if RCount_M(j) == 2
                Rcount_M = 0;
                RCount_M(j) = 0;
            end
        end
    end
end
%%
for j = 1: length(Trust_value_M)
    Thr_M = 0.6;
    if Trust_value2_M_(j) < Thr_M
        wcount2_M = wcount2_M + 1;
        wCount2_M(j) = wcount2_M;

        if wCount2_M(j) == 3
            Rcount2_M = Rcount2_M + 1;
            RCount2_M(j) = Rcount2_M;
            wcount2_M = 0;
            wCount2_M(j) = 0;
            if RCount2_M(j) == 2
                Rcount2_M = 0;
                RCount2_M(j) = 0;
            end
        end
    end
end
%%
for j = 1: length(Trust_value_M)
    Thr_M = 0.6;
    if Trust_value3_M_(j) < Thr_M
        wcount3_M = wcount3_M + 1;
        wCount3_M(j) = wcount3_M;

        if wCount3_M(j) == 3
            Rcount3_M = Rcount3_M + 1;
            RCount3_M(j) = Rcount3_M;
            wcount3_M = 0;
            wCount3_M(j) = 0;
            if RCount3_M(j) == 2
                Rcount3_M = 0;
                RCount3_M(j) = 0;
            end
        end
    end
end

