####  setup_DDA.jl  ####

include("DDAfunctions.jl");                                            # set of Julia functions

nr_delays=2; order=2; N_MOD=3;                                         # DDA model parameters:
                                                                       # 2 delays, quadratic nonlinearity, 
                                                                       # 3 term models 
(MOD,P,SSYM) = make_MOD_new_new(N_MOD,nr_delays,order);                # make list of models MOD,
                                                                       # multicombination P,
                                                                       # and symmetry for each model

dm=4;                                                                  # data points for numerical derivative

TM=50;                                                                 # maximal delay
DELAYS=collect(1:TM);                                                  # all possible delays

make_TAU_ALL(SSYM,DELAYS);                                             # make delay files for each symmetry

N_Trial=50;

SNR=5;                                                                 # signal to noise ratio in dB

DDA_DIR = @sprintf("DDA_%d", SNR);                                     # output folder
dir_exist(DDA_DIR);
 
DATA_DIR = @sprintf("DATA_%d", SNR);                                   # data folder
dir_exist(DATA_DIR);

CV_DIR = @sprintf("CV_%d", SNR);                                       # folder for cross-validation
dir_exist(CV_DIR);


