####  make_data.jl  ####

#include("DDAfunctions.jl");                                           # set of Julia functions

if !isfile("b_values.jld2")
   b=collect(LinRange(.37,.44,1000));b=b[randperm(1000)];b1=sort(b[1:9]);
   b=collect(LinRange(.46,.49,1000));b=b[randperm(1000)];b2=sort(b[1:9]);
                                                                       # randomly generate "subjects"
                                                                       # i.e.: 9 b-values for 2 conditions 
   
   COND=[1 2];bList=[b1 b2]';
   
   @save "b_values.jld2" COND bList                                    # save values if file does not exist 
else
   @load "b_values.jld2"                                               # load the b-values if file exists
end

a = 0.2; c = 5.7;                                                      # parameters in the Roessler system
SNR=5;                                                                 # signal to noise ratio in dB
DATA_DIR = @sprintf("DATA_%d", SNR);                                   # data location
dir_exist(DATA_DIR);                                                   # make folder if it does not exist

dt = .05; DIM=3; order=2;                                              # choice of parameters 

MOD_nr=[0 2 0 3 1 1 1 2 2 0 2 3 2 6];                                  # system encoding for Roessler

N_Trial=50;                                                            # 50 trials
L = 5 * (500 + 512) * N_Trial;                                         # integration length 
                                                                       # (500+512) data points for 50 trials
                                                                       # multiply with 5 for 
                                                                       # later downsampling of data

for n_COND = 1:size(COND,2)                                            # loop over conditions
    for n_SUBJ = 1:size(bList,2)                                       # loop over subjects
        @printf("%d %2d\n", n_COND, n_SUBJ);
        FN = @sprintf("%s/%d_%02d.ascii", DATA_DIR, n_COND, n_SUBJ);   # file name
        if !isfile(FN)                                                 # check if file exists
           X0 = rand(3,1) * 10;                                        # random initial conditions	
           MOD_par = [-1, -1, 1, a, bList[n_COND, n_SUBJ], -c, 1];     # model parameters
  
           X=integrate_ODE_general_BIG(MOD_nr,MOD_par,dt,L,DIM,order,X0,"",1,5,5000);
                                                                       # integrate with a transient of 5000 
           
           X = X[1:5:end,1];                                           # take only x and every 5th point
           X = reshape(X,(500 + 512),50);                              # reshape the data
  
           x = zeros(512,50);                                          # 50 time series, 512 data points each
           for kk = 1:50
               R = shuffle(1:500);                                     # randomized start of trial
               x[:,kk] = add_noise(X[R[1]:R[1] + 511,kk], SNR);        # add noise and radomize start of trial
           end;
  
           writedlm(FN, map(number_to_string, x),' ');                 # write to ascii file
  
           X = nothing; GC.gc();                                       # clear variable X
        end
   end; 
end; 

