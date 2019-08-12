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

########################################################################################
#  0.) include web page template (CSS and java scripts)
########################################################################################
{
    cat   $AliPhysics_SRC/PWGPP/scripts/html/tableBrowser_v1.html  >   $outputHtml
    echo "<body class=\"content\">" >> $outputHtml
    echo "  <div class=\"topWrapper\">" >> $outputHtml
    echo "    <div class=\"leftDiv\">" >> $outputHtml
}

#######################################################################################
#  1.) Make period tabs -
#      code snippet https://www.w3schools.com/howto/tryit.asp?filename=tryhow_js_tabs_fade
#           a.) Make buttons
#           b.) Make tab content
#######################################################################################
{
    # 2.a) Make buttons ( intendation  6)
    printf '\n      <div class="tab">\n' >>$outputHtml
    for a in `ls tab*.html`; do
        tname=`echo $a |sed s_"\.html"__`;
        # echo "<li><a href=\"#\" onclick=\"openQATab('$tname')\"\>$tname</a></li>" >>$outputHtml
        # example for w3c
        #  <button class="tablinks" onclick="openCity(event, 'London')">London</button>
        echo '       <button class="tabLinks" onclick="openQATab(event,'"'"$tname"'"')">'$tname'</button>' >>$outputHtml
    done;
    echo "      </div>" >>$outputHtml
    # 2.b Make tab content

    #cat treePeriodTable.html  >> $outputHtml  ###### Period Table for the upper left corner
    for a in `ls tab*.html`; do
        tname=`echo $a |sed s_"\.html"__`;
        printf  '\n      <div id="'$tname'" class="QATab">\n'>>$outputHtml
        #printf  '<h2 style="margin-bottom:3px;margin-top:3px">TPC QA trending</h2>' >>$outputHtml
        cat $a  | sed s_"^"_"          "_  >>$outputHtml
        printf "        </div>\n\n" >>$outputHtml
    done;
    printf "      </div>\n" >> $outputHtml
}

###########################################################################################################
# 2.) Draw canvas
###########################################################################################################
{
    echo "    <div class=\"rightDiv\" id=\"canvasDiv\">"  >> $outputHtml
    echo "        <canvas width=\"800\" height=\"$heightDraw\"  id=\"canvasDraw\"></canvas>"  >> $outputHtml
    echo "    </div>"  >> $outputHtml
    echo ""  >> $outputHtml
    echo "</div>"  >> $outputHtml
}
###########################################################################################################
# 3.) Write Custom Query Box and table
###########################################################################################################
{
echo "<div class=\"lowerDiv\">"  >> $outputHtml
echo "<table border=\"0\" cellpadding=\"1\" cellspacing=\"2\">"  >> $outputHtml
echo "    <tbody>" >> $outputHtml
echo "        <tr>"  >> $outputHtml
echo "            <td>Custom query:</td>"  >> $outputHtml
echo "            <td><input id=\"globalSelectionMI\" class=\"globalSelectionMI\" name=\"globalSelectionMI\" type=\"text\" size=\"50\"></td>"  >> $outputHtml
echo "            <td>Custom draw:</td>"  >> $outputHtml
echo "            <td><input id=\"globalDrawMI\" class=\"globalDrawMI\" name=\"globalDrawMI\" type=\"text\", size=\"50\"></td>"  >> $outputHtml
echo "        </tr>"  >> $outputHtml
echo "    </tbody>"  >> $outputHtml
echo "</table>"  >> $outputHtml
echo '<table id="runTable" class="display" cellspacing="0" width="100%">' >>$outputHtml
cat $inputTable | grep -v "<table"   >> $outputHtml  ##### Run Table for the lower corner
echo "</div>" >> $outputHtml
#echo "</document>"       >> $outputHtml
}
