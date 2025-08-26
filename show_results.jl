####  show_results.jl  ####

include("setup_DDA.jl");                                               # set of Julia functions

@load "b_values.jld2";                                                 # b-values
NN=[1 1].*size(bList,2);                                               # number of 'subjects'

N_Trial=50;                                                            # number of trials 

AP=fill(NaN,size(MOD,1),(N_MOD+1)+1+2);                                # A' matrix for CV
                                                                       # length of weight vector is L_AF+1 
                                                                       # L_AF=N_MOD+1

for mm=1:size(MOD,1)                                                   # loop over all models
    (MODEL,SYM,model,L_AF) = make_MODEL_new(MOD,SSYM,mm);              # model parameters
    TAU_name = @sprintf("TAU_ALL__%d_%d",SYM[1],SYM[2]);               # file with delay pairs
    TAU=readdlm(TAU_name,' ',Int64); N_TAU=size(TAU,1);                # number of delays

    OUT_DIR = @sprintf("%s/%s",DDA_DIR,model);                         # output folder
    CV_FN = @sprintf("%s/%s",CV_DIR,model);                            # CV file
    FN_JLD = @sprintf("%s/%s.jld2",CV_DIR,model);                      # jld2 file

    Aprime=load(FN_JLD,"Aprime");                                      # A' for one model and all delays
    W=load(FN_JLD,"W");                                                # weights

    tau=findall(Aprime[:,1].==maximum(Aprime))[1];                     # find delays with highest A'

    AP[mm,1]=Aprime[tau][1];                                           # highest A' and
    AP[mm,2]=tau;                                                      # corresponding delay index and
    AP[mm,3:end]=W[tau,:];                                             # weights

    Aprime=nothing; W=nothing;
end

mm=findall(AP[:,1].==maximum(AP[:,1]))[1];                             # find model with highest A'
Aprime=AP[mm,1]; tau=AP[mm,2]; W=AP[mm,3:end];                         # corresponding delay index and weights
(MODEL,SYM,model,L_AF) = make_MODEL_new(MOD,SSYM,mm);                  # model parameters
TAU_name = @sprintf("TAU_ALL__%d_%d",SYM[1],SYM[2]);                   # file with delay pairs
TAU=readdlm(TAU_name,Int64); TAU=TAU[Int(tau),:];                      # best delays

FN_results="best_SNR5.jld2";                                           # save parameters for best A'
@save FN_results Aprime W TAU mm tau                                   # in file

D=fill(NaN,2,N_Trial,maximum(NN));                                     # features for best DDA model
for n_COND = 1:size(COND,2)                                            # loop over conditions and
    for n_SUBJ = 1:size(bList,2)                                       # subjects
        @printf("%d %2d\n", n_COND, n_SUBJ);
        FN_DDA = @sprintf("%s/%s/%d_%02d.jld2",DDA_DIR,model,n_COND,n_SUBJ);
                                                                       # DDA file name
        Q=load(FN_DDA,"Q");Q=Q[:,:,Int(tau)];                          # load features
        Q=[ones(1,N_Trial);Q]'*W;                                      # distance from hyperplane

        D[n_COND,:,n_SUBJ]=Q;                                          # combined distances from hyperplane

        Q = nothing; GC.gc();
    end
end

Q=permutedims(D,[2,3,1]);                                              # plot results
Q=reshape(Q,N_Trial,length(bList));
b=reshape(bList',length(bList),1);
plot(b,Q',
     seriestype=:scatter,
     color=:red,markerstrokecolor=:red,
     legend=false,
     markersize=3)

plot!(b,reshape(mean(D,dims=2)[:,1,:]',length(bList),1),
     seriestype=:scatter,
     color=:blue,markerstrokecolor=:blue,
     legend=false,
     markersize=5)

plot!([.37 .49]',[0.5 0.5]',
      legend=false,
      color=:black,
      linewidth=3)

display(current());

Q=permutedims(D,[2,3,1]);
Q=reshape(Q,N_Trial,length(bList));
Q=reshape(Q,N_Trial*length(bList),1);
L=Int.(reshape([zeros(size(bList,2)*N_Trial,1);
                ones(size(bList,2)*N_Trial,1)],length(bList)*N_Trial,1));
ROC_area(Q,L)

