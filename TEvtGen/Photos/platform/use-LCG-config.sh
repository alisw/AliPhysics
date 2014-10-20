#!/bin/sh

if test ! -d ../src; then
  echo "WARNING! Script must be running from PHOTOS/platform directory!"
  return
fi

echo "This option will overwrite all configuration scripts/makefiles"
echo "and modify the configuration procedure to match LCG setup."
echo ""
echo "You will need autotools version 2.59 or higher."
echo ""
echo "Proceed? (Yes/No)"
read ANSWER

ANSWER=`echo $ANSWER | tr "[:upper:]" "[:lower:]"`

if test "$ANSWER" = "yes" || test "$ANSWER" = "y"; then
  echo "Removing previous installation scripts"
  rm -rf ../config* ../make* ../Make*
  rm -rf ../src/make.inc ../src/*/Makefile ../src/photos-fortran/make*
  rm -rf ../examples/config* ../examples/make* ../examples/Make*

  echo "Copying and configuring new scripts"
  cp -rf LCGCONFIG/* ../.
  cd ..
  autoreconf --install --force
	echo "Done."
else
	echo "Aborted."
fi
