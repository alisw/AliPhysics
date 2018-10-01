#!/usr/bin/env bash
# Create browsable web page like in QA TPC layout
#
#
# Example usage:
#  (  aliroot -b -q  $AliPhysics_SRC/PWGPP/QA/scripts/qaTrending.C+ '$AliPhysics_SRC/PWGPP/QA/detectorQAscripts/qaTrendingTRD.C("LHC15o","cpass1_pass1")' )
#  ( source $AliPhysics_SRC/PWGPP/scripts/makeHtml.sh table.html index.html 0);

## create run tabletable
## fill html web page
##
if [ "$#" -lt 3 ]; then
    echo "makeHtml.sh : illegal number of input parameters"
    echo "1. name of output html file"
    echo "2. input table"
    echo "3. height Draw"
    echo "4. comma separated list of static tables"
    exit;
fi
## Get parameters
source $ALICE_ROOT/libexec/alilog4bash.sh
{

    outputHtml=$1
    inputTable=$2
    htmlList=$4
    heightDraw=$3
    alilog_info "outputTable $outputHtml"
    alilog_info "inputTable  $inputTable"
    alilog_info "height draw $heightDraw"
    alilog_info "htmlList    $htmlList"
}
## include template
cat   $AliPhysics_SRC/PWGPP/scripts/html/tableBrowser.html  >   $outputHtml
# make period tabs
## Write Period Table
echo "<body class=\"content\">" >> $outputHtml
echo '<div class="w3-container">' >>$outputHtml
echo '</div>' >>$outputHtml
echo '<ul class="w3-navbar w3-black">' >>$outputHtml
for a in `ls periodTable*.html`; do
 tname=`echo $a | sed s_periodTable__|sed s_"\.html"__`;
 echo "<li><a href=\"#\" onclick=\"openQATab('$tname')\"\>$tname</a></li>" >>$outputHtml
done;
echo "</ul>" >>$outputHtml
##
echo "<div class=\"topWrapper\">" >> $outputHtml
echo "<div class=\"leftDiv\">" >> $outputHtml
#cat treePeriodTable.html  >> $outputHtml  ###### Period Table for the upper left corner
for a in `ls periodTableS*.html`; do
 tname=`echo $a | sed s_periodTable__|sed s_"\.html"__`;
 echo  '<div id="$tname" class="w3-container QATab"'>>$outputHtml
 cat $a >>$outputHtml
 echo "</div>" >>$outputHtml
done;
echo "</div>" >> $outputHtml

if [ $heightDraw -gt 0 ] ; then
## echo "<div class=\"rightDiv\" id=\"canvasDiv\">"  >> $outputHtml
    echo "  <canvas width=\"800\" height=\"$heightDraw\"  id=\"canvasDraw\"></canvas>"  >> $outputHtml
    echo "</div>"  >> $outputHtml
    echo ""  >> $outputHtml
fi;
echo "</div>"  >> $outputHtml
##
## Write Custom Query Box and Run Table
##
echo "<div class=\"lowerDiv\">"  >> $outputHtml
echo "<table border=\"0\" cellpadding=\"1\" cellspacing=\"2\">"  >> $outputHtml
echo "    <tbody>" >> $outputHtml
echo "        <tr>"  >> $outputHtml
echo "            <td>Custom query:</td>"  >> $outputHtml
echo "            <td><input id=\"globalSelectionMI\" class=\"globalSelectionMI\" name=\"globalSelectionMI\" type=\"text\" size=\"50\"></td>"  >> $outputHtml
echo "        </tr>"  >> $outputHtml
echo "    </tbody>"  >> $outputHtml
echo "</table>"  >> $outputHtml
echo '<table id="runTable" class="display" cellspacing="0" width="100%">' >>$outputHtml
cat $inputTable | grep -v "<table"   >> $outputHtml  ##### Run Table for the lower corner
echo "</div>" >> $outputHtml
#echo "</document>"       >> $outputHtml
