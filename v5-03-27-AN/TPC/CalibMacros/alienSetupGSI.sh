#example shell script to setup alien environmnet at GSI
# To be used in the proof and batch jobs
# The certificate will be automatically prolonged

# 1.) init alien environment
. .alienv217login 
#export LD_LIBRARY_PATH=/usr/local/pub/debian4.0/x86_64/gcc411-21/alice/alien/v2-17/lib64/:$LD_LIBRARY_PATH
#    or your update alien
#  source /d/alice06/wiechula/bin/alienv  head640812

# 2.) check your pasword file
#     if not htere create it
#     it should be read only for owner
ls -al $HOME/.globus/.pwd 

# 3.) enable alien
alien proxy-init -valid 24:00 -pwstdin < $HOME/.globus/.pwd 
alien-token-init
source /tmp/gclient_env_$UID

 
#
#
#
export alien_CLOSE_SE=ALICE::GSI::SE
#
# 4.) Do what you need
# e.g find the list of files stored at gsi
# 4.a) alien 
#  .b) setSElimit alice::gsi::se 
#  .c) find /alice/data/2009/LHC09c/ root_archive.zip > file://esdalien0.txt
#
# Or use the alien in the aliroot session:
# YOU HAVE TO SPECIFY THE STORAGE
# gSystem->Setenv("alien_CLOSE_SE","ALICE::GSI::SE")
  
