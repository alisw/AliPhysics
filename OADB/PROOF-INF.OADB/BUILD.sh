#! /bin/sh

if [ "" = "clean" ]; then
   make clean
   exit 0
fi

make
rc=$?
echo "rc=$?"
if [ $? != "0" ] ; then
   exit 1
fi
exit 0

