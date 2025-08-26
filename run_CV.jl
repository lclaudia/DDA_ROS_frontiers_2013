####  run_CV.jl  ####

include("setup_DDA.jl");                                               # set of Julia functions

@load "b_values.jld2"                                                  # 'subjects' from Roessler system

TR_TS_FILE = "TRTS.ascii";                                             # TRain TeSt file for CV
N_Trial=50;                                                            # 50 trials per subject
NN=[9 9];                                                              # 9 subjects in each group
N_cv=NN[1]*5;                                                          # repeat CV 45 times

N_train=[6 6];                                                         # 6 subjets for training
if !isfile(TR_TS_FILE)                                                 # make file if it does not exist
   make__TR_TS__subj(NN,N_cv,TR_TS_FILE,N_train); 
end
N_test=NN.-N_train;                                                    # 3 subjets for testing

for mm = shuffle(1:size(MOD,1))                                        # loop over models
    println(mm);

    (MODEL,SYM,model,L_AF) = make_MODEL_new(MOD,SSYM,mm);              # DDA settings
    TAU_name = @sprintf("TAU_ALL__%d_%d",SYM[1],SYM[2]);
    N_TAU=readdlm(TAU_name,Int64); N_TAU=size(N_TAU,1);                # delays
    OUT_DIR = @sprintf("%s/%s",DDA_DIR,model);                         # DDA ouput folder
    CV_FN = @sprintf("%s/%s",CV_DIR,model);                            # CV ouput folder
    FN_JLD = @sprintf("%s/%s.jld2",CV_DIR,model);                      # jld2 file

    if !isfile(FN_JLD)                                                 # check if it exists
       run(`touch $FN_JLD`);                                           # placeholder if it does not exist

       AF=fill(NaN,N_TAU,length(COND),L_AF*N_Trial,maximum(NN));       # matrix with all DDA outputs
                                                                       # code adds the 1 automatically!
       for n_COND = 1:length(COND)                                     # loop over conditions
           for n_SUBJ = 1:size(bList,2)                                # loop over subjects
               @printf("%d %2d\n", n_COND, n_SUBJ);
               FN_DDA = @sprintf("%s/%d_%02d.jld2",OUT_DIR,n_COND,n_SUBJ);
     
               Q=load(FN_DDA,"Q");                                     # load features
               Q=reshape(Q,L_AF*N_Trial,N_TAU)';                       # reshape

               AF[:,n_COND,:,n_SUBJ]=Q;                                # all DDA outputs in one matrix 
               Q = nothing; GC.gc();
           end
       end
       AF=reshape(AF,N_TAU*length(COND),L_AF*N_Trial*maximum(NN));     # reshape DDA matrix
       writedlm(CV_FN, map(number_to_string, AF),' ');                 # write to file
   
       CMD=run_CV_ASCII(CV_FN,CV_FN,TR_TS_FILE,N_TAU,
                        N_train,N_test,N_Trial);                       # run CV

       Aprime=readdlm(join([CV_FN,"_Aprime"]));                        # read file with A' for each delay pair
       W=readdlm(join([CV_FN,"_W"]));                                  # read file with corresponding weights
       #plot(Aprime,seriestype=:scatter,color=:blue,markerstrokecolor=:blue,legend=false,markersize=2)
         
       @save FN_JLD Aprime W                                           # save results in jld2 file
         
       run(`sh -c "rm $(CV_FN)_*"`);                                   # delete all
       run(`rm $CV_FN`);                                               #    ascii files
    end
end   

