#
# Hallo world -  Print AliRoot/Root/Alien system info
#

#
# HOST info
#
echo --------------------------------------
echo 
echo HOSTINFO
echo 
echo HOSTINFO HOSTNAME"      "$HOSTNAME
echo HOSTINFO DATE"          "`date`
echo HOSTINFO gccpath"       "`which gcc` 
echo HOSTINFO gcc version"   "`gcc --version | grep gcc`
echo --------------------------------------    

#
# ROOT info
#
echo --------------------------------------
echo
echo ROOTINFO
echo 
echo ROOTINFO ROOT"           "`which root`
echo ROOTINFO VERSION"        "`root-config --version`
echo 
echo --------------------------------------


#
# ALIROOT info
#
echo --------------------------------------
echo
echo ALIROOTINFO
echo 
echo ALIROOTINFO ALIROOT"        "`which aliroot`
echo ALIROOTINFO VERSION"        "`echo $ALICE_LEVEL`
echo ALIROOTINFO TARGET"         "`echo $ALICE_TARGET`
echo 
echo --------------------------------------

#
# Alien info
#
echo --------------------------------------
echo
echo ALIENINFO
for a in `alien --printenv`; do echo ALIENINFO $a; done 
echo
echo --------------------------------------


#
# Local Info
#
echo PWD `pwd`
echo Dir 
ls -al
