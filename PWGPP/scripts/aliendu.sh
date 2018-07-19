# Code to obtain du information recursivelly. 
# *** Beta version (not optimized) to get something working as critically needed
#
# Future recursive du: 
#  du in alien -  similar (functionality and switehces) to GNU coreutils should be implemented directly in AliEn
#  export information into Elastic DB 
# 
# Example usage:
# 1.) source function   
#      source $AliPhysics_SRC/PWGPP/scripts/aliendu.sh; 
# 2.) Execute aliendu recursivelly:
#      aliendu <basedir> <recursiveSize> <verbosity> <maxNuberOfDiretories>
#      aliendu /alice/cern.ch/user/p/pwg_pp/JIRA/ 1000000000 1 100 
#      aliendu /alice/cern.ch/user/m/miranov/ 1000000000  1 100
#  author: marian.ivanov@cern.ch

 
aliendu(){
    #
    #prefix=/alice/cern.ch/user/m/miranov/
    #
    if [ -z "$ALILOG_HOST" ]; then
	source $ALICE_ROOT/libexec/alilog4bash.sh
    fi;
    if [ $# -ne 4 ] ; then
	printf "\n"
	alilog_error "aliendu: Illegal number of parameters $# - 4 expected"
	alilog_info "aliendu: Example usage"
	printf "aliendu <basedir> <recursiveSize> <verbosity> <maxNuberOfDiretories>\n"
        printf "aliendu /alice/cern.ch/user/m/miranov/ 10000000000  1   20\n\n"
        return;
    fi
    
    local prefix=$1;
    local recursiveSize=$2 
    local verbose=$3
    local maxDirs=$4
    if [ -z "$prefix" ]; then
       echo Empty path;
       return;
    fi;
    #aliensh -c "gbbox du $prefix" 
    local -i dSize
    dSize=`aliensh -c "gbbox du $prefix" | grep " uses " |  gawk '{ print $7}'`
    if [ -z "$dSize" ]; then
       echo Empty $prefix	
       return;	
    fi;	
    printf  "$dSize\t$prefix\n"  
    local myPrefix=$prefix
    if [ $dSize -gt $recursiveSize ] ; then
	#echo $dSize $recursiveSize
	local -i subDirs=`alien_ls $prefix | grep -c ""`
	if [ $subDirs -gt $maxDirs ] ; then  # do not folow long directories
	    return;
	fi;
        for a in `alien_ls  $prefix`; do 
	   aliendu $myPrefix/$a $recursiveSize $verbose $maxDirs
        done        
    fi;
}

makeDUTree(){
    input=$1
    if [ $# -ne 1 ] ; then
	    printf "\n"
	    alilog_error "makeDUtree: Illegal number of parameters $# - 1 expected"
	    alilog_infi 
	    return
    fi;
    echo "dsize/F:hName/C:depth/I:basename/C"
    cat $input| sort -k1 -r -n | while read CMD; do
	    line=`echo $CMD | tr -s /`
	    dName=`echo $line | gawk '{print $2}'`
	    bName=`basename $dName`
            nDirs=`echo $line|  awk 'BEGIN{FS="/"} {print NF?NF-1:0}'`
	    printf "$line\t$nDirs\t$bName\n"
    done
}

# example usage:
#   aliendu /alice/cern.ch/user/p/pwg_pp/JIRA/ 2000000000000 1  100 >   pwg_pp.JIRA.dutree
#   aliendu /alice/cern.ch/user/p/pwg_pp/triggeredRaw 2000000000000 1  100 >   pwg_pp.triggeredRaw.dutree
#   aliendu /alice/cern.ch/user/p/pwg_pp/ 2000000000000 1  100 >   pwg_pp.dutree
#   makeDUTree pwg_pp.dutree  > pwg_pp.tree
# TTree tree
# tree.ReadFile("pwg_pp.tree")
