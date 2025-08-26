#=
import Pkg; Pkg.add("Combinatorics")
import Pkg; Pkg.add("LinearAlgebra")
import Pkg; Pkg.add("Printf")
import Pkg; Pkg.add("DelimitedFiles")
import Pkg; Pkg.add("Plots")
import Pkg; Pkg.add("StatsBase")

import Pkg; Pkg.add("LaTeXStrings")
import Pkg; Pkg.add("Graphs")
import Pkg; Pkg.add("GraphRecipes")
import Pkg; Pkg.add("Colors")
import Pkg; Pkg.add("JLD2")
=#

using Plots
using LaTeXStrings
using Colors
using Printf
using StatsBase
using DelimitedFiles
using Combinatorics
using LinearAlgebra
using Graphs
using GraphRecipes
using JLD2

using Random: shuffle



if Sys.iswindows()
   SL="\\";
else
   SL="/";
end

function index(DIM, ORDER)
    B = ones(DIM^ORDER, ORDER)
    
    if DIM > 1
        for i = 2:(DIM^ORDER)
            if B[i-1, ORDER] < DIM
                B[i, ORDER] = B[i-1, ORDER] + 1
            end
            
            for i_DIM = 1:ORDER-1
                if round((i/DIM^i_DIM - floor(i/DIM^i_DIM))*DIM^i_DIM) == 1
                    if B[i-DIM^i_DIM, ORDER-i_DIM] < DIM
                        for j = 0:DIM^i_DIM-1
                            B[i+j, ORDER-i_DIM] = B[i+j-DIM^i_DIM, ORDER-i_DIM] + 1
                        end
                    end
                end
            end
        end
        
        i_BB = 1
        BB = Vector{Int}[]
        for i = 1:size(B,1)
            jn = 1
            for j = 2:ORDER
                if B[i, j] >= B[i, j-1]
                    jn += 1
                end
            end
            if jn == ORDER
                push!(BB, B[i, :])
                i_BB += 1
            end
        end
    else
        println("DIM=1!!!")
    end
    
    return hcat(BB...)
end

function monomial_list(nr_delays, order)
    # monomials
    P_ODE = index(nr_delays+1, order)'
    P_ODE = P_ODE .- ones(Int64,size(P_ODE))
    
    P_ODE = P_ODE[2:size(P_ODE,1),:];
    
    return P_ODE
end

function make_MOD_nr(SYST,NrSyst)

   DIM=length(unique(SYST[:,1]));
   order=size(SYST,2)-1;

   P=[[0 0]; monomial_list(DIM*NrSyst,order)];

   MOD_nr=fill(0,size(SYST,1)*NrSyst,2);
   for n=1:NrSyst
       for i=1:size(SYST,1)
           II=SYST[i,2:end]';
           II[II .> 0] .+= DIM*(n-1);

           Nr=i+size(SYST,1)*(n-1);
           MOD_nr[Nr,2]=findall( sum(abs.(repeat(II,size(P,1),1)-P),dims=2)' .== 0 )[1][2] - 1;
           MOD_nr[Nr,1]=SYST[i,1]+DIM*(n-1);
       end
       #P[MOD_nr[1:size(SYST,1),2].+1,1:2]
   end
   MOD_nr=reshape(MOD_nr',size(SYST,1)*NrSyst*2)';

   return MOD_nr,DIM,order,P

end


function integrate_ODE_general_BIG(MOD_nr,MOD_par,dt,L,DIM,ODEorder,X0,FNout,CH_list,DELTA,TRANS=nothing)
  if TRANS===nothing
     TRANS=0;
  end

  if Sys.iswindows()
     if !isfile("i_ODE_general_BIG.exe")
        cp("i_ODE_general_BIG","i_ODE_general_BIG.exe");
     end

     CMD=".\\i_ODE_general_BIG.exe";
  else
     CMD="./i_ODE_general_BIG";
  end

  MOD_NR = join(MOD_nr, " ");
  CMD = "$CMD -MODEL $MOD_NR";
  MOD_PAR = join(MOD_par, " ");
  CMD = "$CMD -PAR $MOD_PAR";
  ANF=join(X0," ");
  CMD = "$CMD -ANF $ANF";
  CMD = "$CMD -dt $dt";
  CMD = "$CMD -L $L";
  CMD = "$CMD -DIM $DIM";
  CMD = "$CMD -order $ODEorder";
  if TRANS>0
     CMD = "$CMD -TRANS $TRANS";
  end
  if length(FNout)>0
     CMD = "$CMD -FILE $FNout";
  end
  CMD = "$CMD -DELTA $DELTA";
  CMD = "$CMD -CH_list $(join(CH_list," "))";

  if length(FNout)>0
     if Sys.iswindows()
        run(Cmd(string.(split(CMD, " "))));
     else
        run(`sh -c $CMD`);
     end
  else
     if Sys.iswindows()
       X = read(Cmd(string.(split(CMD, " "))),String);
     else
       X = read(`sh -c $CMD`,String);
     end
     X = split(strip(X), '\n');
     X = hcat([parse.(Float64, split(row)) for row in X]...)';

     return X
  end
end

function dir_exist(DIR)
    if !isdir(DIR)
        mkdir(DIR)
    end
end

function number_to_string(n::Number)
   return @sprintf("%.15f", n);
end

function make_MODEL_new(MOD,SSYM,mm)
   MODEL=findall(x -> x == 1, MOD[mm,:]);
   L_AF=length(MODEL)+1;
   SYM=SSYM[mm,:];
   model = join([lpad(x, 2, '0') for x in MODEL], "_");

   return MODEL,SYM,model,L_AF
end

function make_MODEL(SYST)
   order=size(SYST,2);
   nr_delays=2;

   P_ODE=monomial_list(nr_delays,order);

   MODEL=fill(0,size(SYST,1))';
   for i=1:size(SYST,1)
       II=SYST[i,:]';

       MODEL[i] = findall( sum(abs.(repeat(II,size(P_ODE,1),1)-P_ODE),dims=2)' .== 0 )[1][2];
   end
   #P_ODE[MODEL',:]
   L_AF=length(MODEL)+1;

   return MODEL, L_AF, order
end

function make_MOD_new_new(N_MOD,nr_delays,order)

    # MOD,P_ODE,SSYM=make_MOD_new_new(N_MOD,nr_delays,order);

    # Claudia 06/16/2015 Matlab
    # Claudia 08/17/2023 Julia
    # Claudia 01/27/2024

    #nr_delays=2; order=3; N_MOD=2:3;
    #nr_delays=2; order=2; N_MOD=3;

    if nr_delays!=2
       println("only nr_delays=2 supported");
       global nr_delays
       nr_delays=2;
    end

    P_DDA=monomial_list(nr_delays,order); L=size(P_DDA,1);

    PP = -P_DDA;
    PP[PP.==-1].=2;
    PP[PP.==-2].=1;
    PP=sort(PP,dims=2);

    as_ints(a::AbstractArray{CartesianIndex{L}}) where L = reshape(reinterpret(Int, a), (L, size(a)...));

    f=fill(0,size(P_DDA,1),2);
    for k1=1:size(P_DDA,1)
        f[k1,1]=k1;
        ff = findall(x -> x ==0, sum(abs.(P_DDA-repeat(PP[k1,:],1,size(PP,1))'),dims=2));
        if length(ff)>0
           f[k1,2]=as_ints(ff)[1];
        end
    end

    MOD=fill(0,1,size(P_DDA,1))
    for n_N = 1:length(N_MOD)
        N = N_MOD[n_N];
        C = collect(combinations(1:L, N));
        C = reshape(collect(Iterators.flatten(C)),(size(C[1],1),size(C,1)))';
        M = zeros(Int, size(C, 1), L);

        for c = 1:size(C, 1)
            M[c, C[c, :]] .= 1
        end

        M1=sort(M.*reshape(repeat(1:L,size(M,1)),(size(M,2),size(M,1)))',dims=2)[:,end-N+1:end];
        M2=-M1;

        for k1=1:size(f,1)
            M2[M2.==-f[k1,1]].=f[k1,2];
        end
        M2=sort(M2,dims=2);

        f2=fill(0,size(M1,1),2);
        for k1=1:size(M1,1)
            f2[k1,1]=k1;
            ff = findall(x -> x ==0, sum(abs.(M1-repeat(M2[k1,:],1,size(M2,1))'),dims=2));
            if length(ff)>0
               f2[k1,2]=as_ints(ff)[1];
            end
        end
        f2=sort(f2,dims=2);
        f2=unique(f2,dims=1);
        f2=f2[f2[:,1].!=f2[:,2],2];
        f2=setdiff(1:size(M1,1),f2);
        #M1=M1[f2,:];

        MOD=[MOD;M[f2,:]];
    end
    MOD=MOD[2:end,:];

    SSYM = fill(-1, size(MOD, 1), 2);
    for n_M=1:size(MOD,1)
        p=P_DDA[findall(x->x==1,MOD[n_M,:]),:];

        SSYM[n_M,1]=length(unique([value for value in p if value > 0]));

        p=convert(Matrix{Float64},p); p[p .== 0] .= NaN;
        p1=mod.(p.+2,2).+1;
        p1[isnan.(p1)].=0;p1=Int.(p1);
        p[isnan.(p)].=0;p=Int.(p);

        p1=sort!(p1,dims=2);
        p1=sortslices(p1,dims=1);

        if sum(abs.(p-p1))==0
           SSYM[n_M,2]=1;
        else
           SSYM[n_M,2]=0;
        end
    end

    return MOD,P_DDA,SSYM
end


function add_noise(s,SNR)

   # check the length of the noise free signal
   N = length(s);

   # n is the noise realization, make it zero mean and unit variance
   n = randn(N);
   n .= (n.-mean(n))./std(n);
   # c is given  from SNR = 10*log10( var(s)/var(c*n) )
   c= sqrt( var(s)*10^-(SNR/10) );
   
   s_out = (s+c.*n);

   return s_out

end

function make_TAU_ALL(SSYM,DELAYS)

   uSYM=unique(SSYM,dims=1);
   for k=1:size(uSYM,1)
      s=uSYM[k,:]'; nr=s[1]; sym=s[2];
    
      FN=@sprintf("TAU_ALL__%d_%d",s[1],s[2]);
      if !isfile(FN)
        fid=open(FN,"w");        
        if nr==1
           for tau1 = 1:length(DELAYS)
               @printf(fid,"%d\n",DELAYS[tau1]);
           end
        elseif nr==2
           if sym==0
              for tau1 = 1:length(DELAYS)
                  for tau2 = 1:length(DELAYS)
                      if tau1!=tau2
                         @printf(fid,"%d %d\n",DELAYS[tau1],DELAYS[tau2]);
                      end
                  end
              end
           elseif sym==1
              for tau1 = 1:length(DELAYS)
                  for tau2 = 1:length(DELAYS)
                      if tau1<tau2
                         @printf(fid,"%d %d\n",DELAYS[tau1],DELAYS[tau2]);
                      end
                  end
              end
           end
        end
        close(fid);
      end
   end
end
    

function generate_TR(M,N)
   TR=shuffle(1:M);TR=repeat(sort(TR[1:N].-1),1,M)';
   for k=2:M
       TR[k,:]=TR[k-1,:].+1;
   end
   TR=mod.(TR,M).+1;
   TR=sortslices(sort(TR,dims=2),dims=1);

   TS=zeros(Int64,M,M-N);
   for k=1:M
       TS[k,:]=setdiff(1:M,TR[k,:]);
   end

   return TR, TS
end

function make__TR_TS__subj(NN,N_cv,TR_TS_FILE="TRTS.ascii",TR_NUM=nothing)

 if TR_NUM === nothing
    TR_NUM=Int.(ceil.(NN.*3/4));

    @printf("N_train is [%s]\n",(join([string(x) for x in TR_NUM]," ")))
 end    

 TRTS=fill(NaN,N_cv,maximum(NN));
 for n_NN=1:length(NN)
     M=NN[n_NN]; N=TR_NUM[n_NN];

     TR=zeros(Int64,N_cv,N);
     TS=zeros(Int64,N_cv,M-N);
     R=[];
     yn=1;
     while yn==1
        for k=1:floor(N_cv/M)
            if k==1
               (TR, TS)=generate_TR(M,N);
            else
               (TR1, TS1)=generate_TR(M,N);
               TR=[TR;TR1]; TS=[TS;TS1];
            end
        end
        if ceil(N_cv/M)>floor(N_cv/M)
           (TR1,TS1)=generate_TR(M,N);
        
           MinMax=1000;
           for k=1:1000
               r1=shuffle(1:M); r1=r1[1:(N_cv-size(TR,1))];
               r=TR1[r1,:];  
               #histogram(r[:], bins=1:M, legend=false)
               h1 = fit(Histogram, vec(r), collect(1:M + 1)); h2 = h1.weights;
               m=maximum(h2)-minimum(h2);
               if m<MinMax
                  MinMax=m; 
                  R=r1;
               end
           end
           TR=[TR;TR1[R,:]];TS=[TS;TS1[R,:]];
        end
        #display(histogram(TR[:], bins=1:M, legend=false))
        #if size(remove_duplicate_rows(TR),1)==N_cv
        #   yn=0;
        #end
        
        tr=sortslices(sort(TR[1:M:end,:],dims=2),dims=1);
        if length(findall(x->x==0,sum(abs.(diff(tr,dims=1)),dims=2)))==0
           yn=0;
        end
     end

     if n_NN==1
        TRTS=hcat(TR,TS);
     else
        TRTS=vcat(TRTS, hcat(TR,TS));
     end
 end

 fid=open(TR_TS_FILE,"w");
 for i1=1:size(TRTS,1)
     for i2=1:size(TRTS,2)
         @printf(fid,"%5d ",TRTS[i1,i2]);
     end
     @printf(fid,"\n");
 end
 close(fid);

 #writedlm(TR_TS_FILE, TRTS,' ');

 return TR_NUM

end



function run_CV_ASCII(AF_file,CV_file,TR_TS_FILE,N_TAU,N_train,N_test,N_trial)

  if Sys.iswindows()
     if !isfile("run_CV_ASCII.exe")
        cp("run_CV_ASCII","run_CV_ASCII.exe");
     end

     CMD=".\\run_CV_ASCII.exe";
  else
     CMD="./run_CV_ASCII";
  end

  CMD="./run_CV_ASCII";
  CMD = "$CMD -N_train $(join([string(x) for x in N_train]," "))";
  CMD = "$CMD -N_test $(join([string(x) for x in N_test]," "))";
  CMD = "$CMD -N_TRIAL $N_trial";
  CMD = "$CMD -N_TAU $N_TAU";
  CMD = "$CMD -TR_TS_file $TR_TS_FILE";
  CMD = "$CMD -IN_FN $CV_file";
  CMD = "$CMD -OUT_FN $CV_file";

  if Sys.iswindows()
     run(Cmd(string.(split(CMD, " "))));
  else
     run(`sh -c $CMD`);    
  end

  return CMD
end

function ROC_area(Y,L)

    n=size(Y,1);
    s = sortperm(Y,dims=1);
    f = findall(L[s] .== 1);
    s = length(f);
    A = (mean([index[1] for index in f]) - (s + 1) / 2) / (n - s);

    return A

end


function deriv_all(data, dm, order=nothing, dt=nothing)
   if order===nothing 
      order=2 
   end
   if dt===nothing 
      dt=1
   end

   t=collect(1+dm:length(data)-dm)
   L=length(t)

   #### second order:
   
   if order==2
      ddata = zeros(L)
      for n1=1:dm
          ddata += (data[t.+n1].-data[t.-n1])/n1
      end
      ddata /= (dm/dt);  
   end
   
   #### third order:
   
   if order==3
      ddata = zeros(L)
      d=0;
      for n1=1:dm
          for n2=n1+1:dm
              d+=1
              ddata -= (((data[t.-n2].-data[t.+n2])*n1^3- 
                      (data[t.-n1].-data[t.+n1])*n2^3)/(n1^3*n2-n1*n2^3))
          end
      end
      ddata /= (d/dt);  
   end
   
   return ddata
end


