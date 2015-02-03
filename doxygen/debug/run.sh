#!/bin/bash

cd "${ALICE_ROOT}"/../src/

[[ $EDITOR == '' ]] && EDITOR=vim

# TPC
#Files=$( find TPC/ -maxdepth 1 -name '*.h' -or -name '*.cxx' -or -name '*.C' )
#Files=$( find TPC/Attic/ -maxdepth 1 -name '*.h' -or -name '*.cxx' -or -name '*.C' )

Files=( 'TPC/Attic/AliTPCPid.cxx' )

while [[ $# -gt 0 ]] ; do
  case "$1" in
    -r) RestoreOnly=1 ;;
    -x) DoxygenOnly=1 ;;
    -d) Debug='--debug=debug' ;;
    -o) Stdout='-o' ;;
  esac
  shift 1
done

cd "${ALICE_ROOT}"/../src/doxygen

for F in ${Files[@]} ; do
  F="${ALICE_ROOT}/../src/${F}"
  if [[ $DoxygenOnly != 1 ]] ; then
    git checkout "$F"
    r=$?
  fi
  if [[ $RestoreOnly != 1 ]] ; then
    ./thtml2doxy.py $Stdout $Debug $( cat debug/includes.txt ) "$F"
    r=$?
  fi
  if [[ $r != 0 ]] ; then
    echo -e "\033[31mFatal error at ${F}: stopping here\033[m"
    break
  else
    while [[ 1 ]] ; do
      git diff "$F"
      echo ''
      FNorm=$( cd "`dirname "$F"`";pwd )/$( basename "$F" )
      AliRootLen=${#ALICE_ROOT}
      FNorm=${FNorm:$AliRootLen}
      echo -e "\033[35mFile \033[34m${FNorm}\033[m"
      echo -e -n "\033[35mWhat to do? \033[34medit, stage, continue\033[m \033[35m> \033[m"
      read ans
      case "$ans" in
        edit)
          "$EDITOR" "$F"
        ;;
        stage)
          git add "$F" || exit 1
          break
        ;;
        continue)
          break
        ;;
      esac
    done
  fi
done

echo -e "\033[32mAll OK\033[m"
