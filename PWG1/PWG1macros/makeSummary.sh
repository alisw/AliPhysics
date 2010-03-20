nSeg=`cat out*.log  | grep violation | grep -c egmentation`


cat  out*.C.log | grep SysInfoMem
cat  out*.C.log | grep SysInfoTime
echo LogSize:"     "`wc -l  < out*.log`
echo HasSegFault:  $nSeg


if [ $nSeg -gt 0 ] ; then
    echo SegFault `pwd`
    cat out*.log  | grep 0x00
fi

