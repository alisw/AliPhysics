#! /bin/sh

# dirty hack!
# touch libPWGUDselectors.pkg
# touch libPWGUDdep.pkg

if [ "" = "clean" ]; then
   make clean
   exit 0
fi

make libPWGUDbase.so
rc=$?
echo "rc=$?"
if [ $? != "0" ] ; then
   exit 1
fi
exit 0
