#!/bin/bash

if [ "0$1" == "0" ]; then
    echo "No bream type specified, assuming \"Pb-Pb\""
    echo "(Valid bram types: pp, Pb-Pb)"
    echo ""
    BEAMTYPE=Pb-Pb
elif [ $1 != "pp" ] && [ $1 != "Pb-Pb" ]; then
    echo "Invalid beam type specified: \"$1\""
    echo "(Valid bram types: pp, Pb-Pb)"
    echo ""
    BEAMTYPE=Pb-Pb
else
    BEAMTYPE=$1
fi

echo "#axis definitions (PbPb)"
echo "#name,nbins,low,high"
echo "resetIncludingDownstream=1"

#axis='nClustersSPD,100,0.,30e3'
#fAxes["nClustersSPD"].set( 100, 0., 800.,  &fnClustersSPD );

#histogram=',fHistSPDclusters_SPDrawSize,SPD clusters vs SPD raw size,rawSizeSPD,nClustersSPD'
#NewHistogram(",fHistSPDclusters_SPDrawSize,SPD clusters vs SPD raw size,rawSizeSPD,nClustersSPD");

cat AliHLTGlobalPromptRecoQAComponent.cxx | sed -n -e "/Start Axes for $BEAMTYPE/,/End Axes for $BEAMTYPE/ p" | grep "^ *fAxes\[\"" | sed "s/ //g" | sed "s/fAxes\[\"/axis='/" | sed "s/\"\]\.set(/,/" | sed "s/,*&[a-zA-Z0-9_]*);/'/"

cat AliHLTGlobalPromptRecoQAComponent.cxx | sed -n -e '/Start Common Axes/,/End Common Axes/ p' | grep "^ *fAxes\[\"" | sed "s/ //g" | sed "s/fAxes\[\"/axis='/" | sed "s/\"\]\.set(/,/" | sed "s/,*&[a-zA-Z0-9_]*);/'/"

echo "#histogram definitions"
echo "#triggerRegEx,histName,histTitle,xAxisName,yAxisName"

cat AliHLTGlobalPromptRecoQAComponent.cxx | sed -n -e '/Start Histograms/,/End Histograms/ p' | grep "^ *NewHistogram(" | sed "s/ //g" | sed "s/NewHistogram(\"/histogram='/" | sed "s/\");/'/"
