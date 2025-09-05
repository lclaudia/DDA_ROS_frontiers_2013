####  run_DDA.jl  ####

include("setup_DDA.jl");                                               # setup parameters and 
                                                                       # define set of Julia functions

if isfile("b_values.jld2")
   @load "b_values.jld2"                                               # load the b-values if file exists
end

for mm = shuffle(1:size(MOD,1))                                        # loop over randomized model numbers
    println(mm);

    (MODEL,SYM,model,L_AF) = make_MODEL_new(MOD,SSYM,mm);              # parameters of specific model
    TAU_name = @sprintf("TAU_ALL__%d_%d",SYM[1],SYM[2]);               # delay file
    TAU = readdlm(TAU_name); N_TAU=size(TAU,1);
    OUT_DIR = @sprintf("%s/%s",DDA_DIR,model);                         # output folder
    dir_exist(OUT_DIR);                                                # create if it does not exist

    for n_COND = shuffle(1:size(COND,2))                               # loop over conditions
       for n_SUBJ = shuffle(1:size(bList,2))                           # loop over subjects
          FN_data = @sprintf("%s/%d_%02d.ascii",                       # file names
                              DATA_DIR, n_COND, n_SUBJ);
          FN_DDA = @sprintf("%s/%d_%02d.DDA",  
                             OUT_DIR, n_COND, n_SUBJ);
          FN_JLD = @sprintf("%s/%d_%02d.jld2",  
                             OUT_DIR, n_COND, n_SUBJ);
          if !isfile(FN_JLD)                                           # check if exists
             run(`touch $FN_JLD`);                                     # if not, create placeholder
          
             @printf("%d %2d: %s\n",n_COND,n_SUBJ,FN_JLD);

             if Sys.iswindows()
                if !isfile("run_DDA_AsciiEdf.exe")
                   cp("run_DDA_AsciiEdf","run_DDA_AsciiEdf.exe");
                end
           
                CMD=".\\run_DDA_AsciiEdf.exe";
             else
                CMD="./run_DDA_AsciiEdf";
             end
 
             CMD = "$CMD -ASCII";                                      # DDA executable
             CMD = "$CMD -MODEL $(join(MODEL," "))";                   # model
             CMD = "$CMD -TAU_file $TAU_name";                         # delay file         
             CMD = "$CMD -dm $dm -order $order -nr_tau $nr_delays";    # DDA parameters
             CMD = "$CMD -DATA_FN $FN_data -OUT_FN $FN_DDA";           # input and output files
             CMD = "$CMD -SELECT 1 0 0 0";

             if Sys.iswindows()
                run(Cmd(string.(split(CMD, " "))));
             else
                run(`sh -c $CMD`);
             end

             Q=readdlm(join([FN_DDA,"_ST"]));                          # read output file
             Q=reshape(Q[3:end],L_AF,50,N_TAU);                        # reshape

             @save FN_JLD Q                                            # save DDA features

             Q = nothing;GC.gc();                                      # clear variable
          end
       end
    end
end

