#!/bin/bash

cd "${ALICE_ROOT}"/../src/

[[ $EDITOR == '' ]] && EDITOR=vim

# List of files from external source
Files=$(
  cat "$(dirname "$0")"/list-of-files-to-process.txt | \
  sed -e '/^$/d' | \
  grep -v ^# | \
  sort )

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
    #./thtml2doxy.py $Stdout $Debug $( cat debug/includes.txt ) "$F"
    ./thtml2doxy.py $Stdout $Debug -I$ALICE_ROOT/include "$F"
    r=$?

    # remove trailing whitespace
    cat "$F" | sed -e 's/\s*$//' > "$F".0
    mv "$F".0 "$F"
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
      echo -e -n "\033[35mWhat to do? \033[34medit, stage, parts, continue, restore\033[m \033[35m> \033[m"
      read ans
      case "$ans" in
        edit)
          "$EDITOR" "$F"
        ;;
        stage)
          git add "$F" || exit 1
          break
        ;;
        parts)
          git add -p "$F" || exit 1
          git checkout "$F" || exit 1
          echo -e "\033[35mPart of \033[34m${FNorm}\033[35m was staged, the rest was restored. Press \033[34mEnter\033[35m to continue.\033[m"
          read
          break
        ;;
        continue)
          break
        ;;
        restore)
          git checkout "$F" || exit 1
          break
        ;;
      esac
    done
  fi
done

echo -e "\033[32mAll OK\033[m"
