#############################################################################################
Macros to do fast simulation of processes important for tuning of reconstruction.
Currently fast simulation of ionization digitization and cluster finder - AliTPCclusterFast
#############################################################################################


How to use it?
a) which macro to use (I know it was somewhere in AliRoot but with the GIT page I dont find it anymore),
b) what is the basic functionality of the functions
c) what do I need to run it (aliroot version?) 
d) how would I run it and extract space point resolution and dEdx resolution
    (best step by step for dummies)?




a) Which macro to use (I know it was somewhere in AliRoot

   Example case - submit 40 jobs with 100 tracks.
   
   source $ALICE_ROOT/TPCdev/TPC/fastSimul/simul.sh
   makeEnvLocal             #this is just example please setup your environmnet script to set env variables 
   makeSubmitRUN 40 100

c) what do I need to run it (aliroot version?) 
   Recent AliRoot


b) What is the basic functionality of the functions?
   Provides cluster and track functionitility - possible to modify parameters of the reconstruction/ resp. harware setup   See example usage simul.C:DrawdEdxResolExample()
   
   
   

