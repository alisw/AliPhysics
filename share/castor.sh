if [ "$CASTOR_USER_POOL" = "" ]
  then
    CASTOR_USER_POOL=$STAGE_POOL 
fi

if [ "$CASTOR_BIN" = "" ]
then
    nsls=`which nsls` 
    CASTOR_BIN=`dirname $nsls`
fi

if [ ! -d $CASTOR_BIN ] 
then
  printf "CASTOR is not supported on this platform.\n"
  exit 255
fi

###########################################################################
ALIFS_Usage()
{
   printf "CASTOR Implementation:\n\n"
   printf "Usage: alifs [-help][-p <pool>] <command [options]>   \n"
   printf "              ls [-cdilRTu] [--class] [--comment] path\n" 
   printf "              mv oldname newname...                   \n"
   printf "              rm [-f] [-i] [-r] dirname...            \n"
   printf "              mkdir [-m absolute_mode] [-p] dirname...\n"
   printf "              cp [-s maxsize] f1 f2                   \n" 
   printf "              cp f1 <dir2>                            \n"
   exit
}
###########################################################################
ALIFS_ls()
{
   $CASTOR_BIN/nsls $*
}
###########################################################################
ALIFS_mkdir()
{
   $CASTOR_BIN/nsmkdir $*
}
###########################################################################
ALIFS_mv()
{
   $CASTOR_BIN/nsrename $*
}
###########################################################################
ALIFS_rm()
{
   $CASTOR_BIN/nsrm $*
}
###########################################################################
ALIFS_cp()
{
   $CASTOR_BIN/rfcp $*
}
