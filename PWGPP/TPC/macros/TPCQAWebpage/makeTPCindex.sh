#
# Example usage:
#  ( source  $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/makeTPCindex.sh trending.root );

# create run tabletable
fileName=$1

aliroot -l -b -q $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/dumpTable.C+\(\"${fileName}\",\"runTable\"\)

#
# fill html web page
#
cat $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/tpcTemplate.html     >   index.html

## Write Period Table
echo "<body class=\"content\">" >> index.html
echo "<div class=\"topWrapper\">" >> index.html
echo "<div class=\"leftDiv\">" >> index.html
cat treePeriodTable.inc  >> index.html ###### Period Table for the upper left corner
echo "</div>" >> index.html

## Write Canvas for Plots
echo "<div class=\"rightDiv\">"  >> index.html
echo "  <canvas width=\"800\" height=\"580\"  id=\"canvasDraw\"></canvas>"  >> index.html
echo "</div>"  >> index.html
echo "</div>"  >> index.html
echo ""  >> index.html

## Write Custom Query Box and Run Table
echo "<div class=\"lowerDiv\">"  >> index.html
echo "<table border=\"0\" cellpadding=\"1\" cellspacing=\"2\">"  >> index.html
echo "    <tbody>" >> index.html
echo "        <tr>"  >> index.html
echo "            <td>Custom query:</td>"  >> index.html
echo "            <td><input id=\"globalSelectionMI\" class=\"globalSelectionMI\" name=\"globalSelectionMI\" type=\"text\" size=\"50\"></td>"  >> index.html
echo "        </tr>"  >> index.html
echo "    </tbody>"  >> index.html
echo "</table>"  >> index.html
cat treeRunTable.inc     >> index.html ##### Run Table for the lower corner
echo "</div>" >> index.html
#echo "</document>"       >> index.html
