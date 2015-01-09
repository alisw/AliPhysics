#!/bin/bash

cd "$( dirname "$0" )"

bin="${PWD}/../thtml2doxy.py"
includeslist="${PWD}/includes.txt"

function process_file() {

  src="$1"
  quiet="$2"  # set it to 1 to activate
  stop="$3"  # set it to 1 to activate

  quiet_flag='-d'
  [[ $quiet == 1 ]] && quiet_flag=''

  t=$( mktemp -d /tmp/thtml2doxy-XXXXX )
  mkdir -p $t/pass1
  mkdir -p $t/pass2

  out1="${t}/pass1/$( basename "${src}" )"
  out2="${t}/pass2/$( basename "${src}" )"

  # pass1
  cd "${t}/pass1"
  "${bin}" -o $quiet_flag $(cat "${includeslist}") "${src}" > "${out1}"
  r1=$?
  if [[ $stop == 1 || $r1 != 0 ]] ; then
    read -p "run1 terminated with ${r1}: type 'sh' to enter a shell..." ans
    [[ "$ans" == 'sh' ]] && bash
  fi

  if [[ $r1 != 0 ]] ; then
    # Abort here
    echo "run1 broken (exitcode ${r1}), makes no sense to continue"
    cd /
    rm -rf "${t}"
    return 1
  fi

  # pass2
  cd "${t}/pass2"
  "${bin}" -o $quiet_flag $(cat "${includeslist}") "${out1}" > "${out2}"
  r2=$?
  if [[ $stop == 1 || $r2 != 0 ]] ; then
    read -p "run2 terminated with ${r2}: type 'sh' to enter a shell..." ans
    [[ "$ans" == 'sh' ]] && bash
  fi

  if [[ $quiet != 1 ]] ; then

    # diff1
    echo '=== BEGIN DIFF1 ==='
    diff -rupN "${src}" "${out1}" | pygmentize -ldiff
    echo '=== END DIFF1 ==='
    echo ''

    # diff2
    if [[ $r2 == 0 ]] ; then
      echo '=== BEGIN DIFF2 ==='
      diff -rupN "${out1}" "${out2}" | pygmentize -ldiff
      echo '=== END DIFF2 ==='
      echo ''
    fi

    # report
    echo '=== BEGIN REPORT ==='
    echo "run1: ${r1}"
    echo "run2: ${r2}"
    echo '=== END REPORT ==='
    echo ''

  else

    # compact output
    cya='\033[36m'
    red='\033[31m'
    gre='\033[32m'
    non='\033[m'

    r1s="r1:${red}fail${non}($r1)"
    [[ $r1 == 0 ]] && r1s="r1:${gre}ok${non}"

    r2s="r2:${red}fail${non}($r2)"
    difs="diff:${red}fail${non}"
    ndif=0
    if [[ $r2 == 0 ]] ; then

      r2s="r2:${gre}ok${non}"

      diff -rupN "${out1}" "${out2}" > "${t}/diff2.txt"
      ndif=$( cat "${t}/diff2.txt" | wc -l )
      [[ $ndif == 0 ]] && difs="diff:${gre}ok${non}"

    fi

    echo -e "${cya}${src}${non} ${r1s} ${r2s} ${difs}"
    if [[ $ndif != 0 ]] ; then
      cat "${t}/diff2.txt" | pygmentize -ldiff

      read -p "unexpected diffs found between run1 and run2: type 'sh' to enter a shell..." ans
      [[ "$ans" == 'sh' ]] && bash

      # Abort
      cd /
      rm -rf "${t}"
      return 1

    fi

  fi

  # cleanup
  cd /
  rm -rf "${t}"
  return 0

}

if [[ ! -e "${includeslist}" ]] ; then
  echo 'no list of includes found'
  exit 1
fi

# e.g. all TPC files
exec 3< <( find "${PWD}/../../TPC/" -name '*.h' )
while read -r -u3 file ; do
  #ext="${file##*.}"
  #base="${file%.*}"
  process_file "${file}" 1 0 || break
done
exec 3<&-

# e.g. test single classes
# class="${PWD}/../../TPC/TPCsim/AliTPC.cxx"
# ext="${class##*.}"
# base="${class%.*}"
# process_file "${base}.cxx" 0 1
