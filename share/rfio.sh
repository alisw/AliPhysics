if [ "$RFIO_BIN" = "" ]
then
    rfdir=`which rfdir` 
    RFIO_BIN=`dirname $rfdir`
fi

if [ ! -d $RFIO_BIN ] 
then
  printf "HPSS/RFIO is not supported on this platform.\n"
  exit 255
fi

###########################################################################
ALIFS_Usage()
{
   printf "HPSS/RFIO Implementation:\n\n"
   printf "Usage: alifs [-help][-p <pool>] <command [options]>   \n"
   printf "              ls [-R] path                            \n" 
   printf "              mv oldname newname...                   \n"
   printf "              rm [-r] pathname...                     \n"
   printf "              mkdir [-m absolute_mode] [-p] dirname...\n"
   printf "              cp [-s maxsize] f1 f2                   \n" 
   printf "              cp f1 <dir2>                            \n"
   exit
}
###########################################################################
ALIFS_ls()
{
   $RFIO_BIN/rfdir $*
}
###########################################################################
ALIFS_mkdir()
{
   $RFIO_BIN/rfmkdir $*
}
###########################################################################
ALIFS_mv()
{
   $RFIO_BIN/rfrename $*
}
###########################################################################
ALIFS_rm()
{
   $RFIO_BIN/rfrm $*
}
###########################################################################
ALIFS_cp()
{
   $RFIO_BIN/rfcp $*
}
