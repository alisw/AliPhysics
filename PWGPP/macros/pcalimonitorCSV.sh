# source $AliPhysics_SRC/PWGPP/macros/pcalimonitorCSV.sh

alias helpCat='pygmentize -O style=borland,linenos=1 -l bash'


init(){
    [[ -z "$ALILOG_HOST" ]] && source ${ALICE_ROOT}/libexec/alilog4bash.sh
    alilog_info "List of functions:"
    helpCat<<HELP_USAGE
"List of functions:"
dumpCSV
printSummary
examleUsage
exampleCase
HELP_USAGE
}

dumpCSV(){
  [[ -z $2 ]] && helpCat <<HELP_USAGE
  "makeDirs"
  Input:
        $1 - production path
        $2 - output file
  Output:
       csv file
  Example usage:
       dumpCSV   https://alimonitor.cern.ch/prod/jobs.jsp?t=20831 table.csv
HELP_USAGE
    [[ $# -ne 2 ]] &&return
    inPath="$1"
    output=$2
    echo dumpCSV ${inPath} -- ${output}
    curl -Lk  --tlsv1 --cert $HOME/.globus/usercert.pem --key $HOME/.globus/userkey.pem -o table.mif   "${inPath}"'&res_path=mif'
    cat table.mif |sed  s_/I:_,_g| sed s_/D:_,_g  |  sed s_/C:_,_g | sed s_\"_,_g | sed s_:_,_g | sed s_/I__ > $output
}

printSummary(){
 [[ -z $2 ]] && helpCat <<HELP_USAGE
  "makeDirs"
  Input:
        $1 - production csv
        $2 - selection
        $3 - where sort statements
  Output:
       JIRA print
  Example usage:
     printSummary JobID,RunNo,input_events,wall_time,outputsize,outputsize/input_events table "where input_events>10000"
HELP_USAGE
    [[ $# -lt 2 ]] &&return
    echo csvsql --query "select $1 from '$2' '${3}' " $2.csv
    csvsql --query "select $1 from '$2' ${3} " $2.csv | csvlook -l

}

exampleUsage(){
    helpCat<<HELP_USAGE
    exampleUsage
    This is example usage of the pcalimonitorCSV.sh
HELP_USAGE
    source $AliPhysics_SRC/PWGPP/macros/pcalimonitorCSV.sh
    dumpCSV "https://alimonitor.cern.ch/prod/jobs.jsp?t=20831&outputdir=pass1" table.csv
    printSummary JobID,RunNo,input_events,wall_time,outputsize,outputsize/input_events table
    #
    dumpCSV      "https://alimonitor.cern.ch/prod/jobs.jsp?t=20831" table.csv
    printSummary JobID,RunNo,input_events,wall_time,outputsize,outputsize/input_events table "where input_events>40000 ORDER BY input_events DESC"

}

exampleCase(){
    curl -Lk  --tlsv1 --cert $HOME/.globus/usercert.pem --key $HOME/.globus/userkey.pem -o table.mif  "https://alimonitor.cern.ch/prod/jobs.jsp?t=20831&res_path=mif"
    cat table.mif |sed  s_/I:_,_g| sed s_/D:_,_g  |  sed s_/C:_,_g | sed s_\"_,_g | sed s_:_,_g | sed s_/I__ > filters.csv
    csvcut -n filters.csv
    csvcut -c JobID,RunNo,input_events,wall_time,outputsize  filters.csv |csvlook -l
    csvsql --query "select JobID,RunNo,input_events,wall_time,outputsize,outputsize/input_events  from 'filters'" filters.csv | csvlook -l
}


init

