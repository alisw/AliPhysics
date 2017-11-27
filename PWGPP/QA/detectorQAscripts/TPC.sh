#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumer     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

runLevelQA()
{
  qaFile=$1

  cp $ALICE_PHYSICS/PWGPP/TPC/macros/MakeTrend.C .
  cp $ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/ConfigCalibTrain.C .

  echo $dataType/$year/$period/$pass/$runNumber
  echo $ocdbStorage
  export edataType=$dataType
  export eyear=$year
  export eperiod=$period
  export epass=$pass
  export erunNumber=$runNumber
  export eocdbStorage=$ocdbStorage
  echo aliroot -b -q -l "ConfigCalibTrain.C\($runNumber,\"$ocdbStorage\"\)  MakeTrend.C(\"$qaFile\",$runNumber)" 
  aliroot -b -q -l ConfigCalibTrain.C\($runNumber,\"$ocdbStorage\"\)  "MakeTrend.C(\"$qaFile\",$runNumber)" 

  cp $ALICE_PHYSICS/PWGPP/TPC/macros/drawPerformanceTPCQAMatch.C .
  echo aliroot -b -q -l " drawPerformanceTPCQAMatch.C(\"$qaFile\")"
  aliroot -b -q -l " drawPerformanceTPCQAMatch.C(\"$qaFile\")"
  makeHTMLindexPerRun
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGPP/TPC/macros/drawPerformanceTPCQAMatchTrends.C .
  cp $ALICE_PHYSICS/PWGPP/TPC/macros/qaConfig.C .
  cp $ALICE_PHYSICS/PWGPP/TPC/macros/makePeriodTrendingTree.C .
  aliroot -b -q -l "drawPerformanceTPCQAMatchTrends.C(\"trending.root\",\"PbPb\")"
  # aliroot -b -q -l "makePeriodTrendingTree.C(\"trending.root\",\"PbPb\")"
  if [[ ${dataType} =~ "sim" ]]; then 
    mcrddir="mcrd_com"
    mkdir -p $mcrddir
  fi

  makeHTMLindexPerPeriod

  if [[ ${dataType} =~ "sim" ]]; then  
    echo "running tpcMCValidation.C in " $PWD
    echo "MC period: $period;" 
    echo "make direcotry: $mcrddir"

        aliroot -q  "$ALICE_PHYSICS/PWGPP/TPC/macros/tpcMCValidation.C+(\"$period\",\"$mcrddir\")"
    cd - 
  fi 
}


makeHTMLindexPerRun()
{
  cat > index.html <<EOF
<div align="left"><br>
<h2>Run Data Quality</h2>
<a href="TPC_event_info.png">Event Information</a><br>
<a href="cluster_occupancy.png">Cluster Occupancy</a><br>
<a href="eta_phi_pt.png">#eta, #phi and pt</a><br>
<a href="cluster_in_detail.png">Number of clusters in #eta and #phi</a><br>
<a href="dca_in_detail.png">DCAs vs #eta</a><br>
<a href="TPC_dEdx_track_info.png">TPC dEdx</a><br>
<a href="dca_and_phi.png">DCAs vs #phi</a><br>
<a href="TPC-ITS.png">TPC-ITS matching</a><br>
<a href="dcar_pT.png">dcar vs pT</a><br>
<a href="pullPhiConstrain.png">Tracking parameter phi</a><br>
<a href="rawQAInformation.png">Raw QA Information</a><br>
<a href="canvasROCStatusOCDB.png">Canvas ROC Status OCDB</a><br>
<!--  <a href="res_pT_1overpT.png">Resolution vs pT and 1/pT</a><br>  //-->
<!--  <a href="eff_all+all_findable.png">Efficiency all charged + findable</a><br>  //-->
<!--  <a href="eff_Pi_K_P.png">Efficiency #pi, K, p</a><br>  //-->
<!--  <a href="eff_Pi_K_P_findable.png">Efficiency findable #pi, K, p</a><br>  //-->
</div>

EOF

}

makeHTMLindexPerPeriod()
{
  #
# Example usage:
#  ( source  $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/makeTPCindex.sh trending.root );

# create run tabletable
fileName="trending.root"

#make a local copy of the external dependences
copyFileFromRemote http://tablefilter.free.fr/TableFilter/tablefilter.js .
copyFileFromRemote http://methvin.com/splitter/splitter.js .

if [[ ${dataType} =~ "sim" ]]; then 
    aliroot -l -b -q $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/rootlogon.C $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/dumpTable.C+\(\"${fileName}\",\"runTable\"\,\"$period\",\"anchorper\",\"anchorpass\"\)
else
    aliroot -l -b -q $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/rootlogon.C $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/dumpTable.C+\(\"${fileName}\",\"runTable\"\)
fi

#
# fill html web page
#

index="index.html"

if [[ ${dataType} =~ "sim" ]]; then 
mcrdindex="${mcrddir}/index_mcrd.html"
else
mcrdindex=""
fi


cat $ALICE_PHYSICS/PWGPP/TPC/macros/TPCQAWebpage/tpcTemplate.html  | tee -a   >   $index $mcrdindex

## Write Period Table
echo "<body class=\"content\">" | tee -a $index $mcrdindex
echo "<div class=\"topWrapper\">" | tee -a $index $mcrdindex
echo "<div class=\"leftDiv\">" | tee -a $index $mcrdindex
cat treePeriodTable.inc  | tee -a >> $index
cat treePeriodTableMCRD.inc  | tee -a >> $mcrdindex ###### Period Table for the upper left corner
echo "</div>" | tee -a >> $index $mcrdindex

## Write Canvas for Plots
echo "<div class=\"rightDiv\">"  | tee -a >> $index $mcrdindex
echo "  <canvas width=\"800\" height=\"580\"  id=\"canvasDraw\"/>"  | tee -a >> $index $mcrdindex
echo "</div>"  | tee -a >> $index $mcrdindex
echo "</div>"  | tee -a >> $index $mcrdindex
echo ""  | tee -a >> $index $mcrdindex

## Write Custom Query Box and Run Table
echo "<div class=\"lowerDiv\">"  | tee -a >> $index $mcrdindex
echo "<table border=\"0\" cellpadding=\"1\" cellspacing=\"2\">"  | tee -a >> $index $mcrdindex
echo "    <tbody>" | tee -a >> $index $mcrdindex
echo "        <tr>"  | tee -a >> $index $mcrdindex
echo "            <td>Custom query:</td>"  | tee -a >> $index $mcrdindex
echo "            <td><input id=\"globalSelectionMI\" class=\"globalSelectionMI\" name=\"globalSelectionMI\" type=\"text\", size=\"50\"></td>"  | tee -a >> $index $mcrdindex
echo "        </tr>"  | tee -a >> $index $mcrdindex
echo "    </tbody>"  | tee -a >> $index $mcrdindex
echo "</table>"  | tee -a >> $index $mcrdindex
cat treeRunTable.inc     | tee -a >> $index $mcrdindex ##### Run Table for the lower corner
echo "</div>" | tee -a >> $index $mcrdindex
echo "</document>"       | tee -a >> $index $mcrdindex




#   cat > $index $mcrdindex <<EOF
# <div align="left"><br>
# <h2>Periodical Data Quality</h2>
# <a href="meanTPCncl_vs_run.png">Mean Number of TPC Clusters</a><br>
# <a href="meanTPCnclF_vs_run.png"># of Found Clusters / # of Findable Clusters</a><br>
# <br>
# <a href="meanMIP_vs_run.png">Mean of MIPs</a><br>
# <a href="meandEdxele_vs_run.png">Mean electron energy loss p(0.32,0.38)GeV/c, dEdx(70,100)</a><br>
# <br>
# <a href="meanVertX_vs_run.png">Mean of Vertex_X</a>, <a href="meanVertY_vs_run.png">Vertex_Y</a>, <a href="meanVertZ_vs_run.png">Vertex_Z</a>
# <br><br>
# <a href="meanMult_vs_run.png">Multiplicities of Primary Tracks</a><br>
# <a href="TPC-ITS-matching-efficiency_vs_run.png">TPC-ITS matching efficiency</a><br>
# <a href="ITS-TPC-matching-quality_vs_run.png">ITS-TPC matching quality</a><br>

# <a href="1overPt_vs_run.png">Delta 1/pt</a><br>
# <a href="DCAOffset_vs_run.png">DCAs</a><br>
# <br>
# <p><font size="4">Runs:</font></p>
# EOF

# local dir
# for dir in 000*; do
#   echo "<a href="${dir}">${dir}</a>" | tee -a >> $index $mcrdindex
# done

#   cat | tee -a >> $index $mcrdindex <<EOF
# <br>
# <p><font size="4">Additional plots</font></p>
# <a href="occ_AC_Side_IROC_OROC_vs_run.png">Nr of Chambers with lower gain (occupancy)</a><br>
# <br>
# <a href="dcar_fitting_run.png">DCAr fitting parameters</a><br>
# <a href="dcar_0_vs_run.png">DCAr fitting parameters (0)</a><br>
# <a href="dcar_1_vs_run.png">DCAr fitting parameters (1)</a><br>
# <a href="dcar_2_vs_run.png">DCAr fitting parameters (2)</a><br>
# <a href="dcaz_0_vs_run.png">DCAz fitting parameters (0)</a><br>
# <a href="dcaz_1_vs_run.png">DCAz fitting parameters (1)</a><br>
# <a href="dcaz_2_vs_run.png">DCAz fitting parameters (2)</a><br>
# <br>
# <a href="resolutionMIP_vs_run.png">Resolution of MIPs</a><br>
# <a href="resolutionMeandEdxEle_vs_run.png">Resolution of mean electron energy loss</a><br>
# <a href="ElectroMIPSeparation_vs_run.png">Separation between electon and MIPs energy loss</a><br>
# <br>
# <a href="MIPattachSlopes_vs_run.png">Attachment parameter p1, A side and C side</a><br>
# <br>
# <a href="pullPhiConstrain_vs_run.png">Tracking parameter phi</a><br>
# <br>
# <br>
# <a href="prodinfo">Production information</a><br>

# </div>

# EOF
}
