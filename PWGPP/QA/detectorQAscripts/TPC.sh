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
  aliroot -b -q -l "drawPerformanceTPCQAMatchTrends.C(\"trending.root\",\"PbPb\")"
  makeHTMLindexPerPeriod
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
<!--  <a href="res_pT_1overpT.png">Resolution vs pT and 1/pT</a><br>  //-->
<!--  <a href="eff_all+all_findable.png">Efficiency all charged + findable</a><br>  //-->
<!--  <a href="eff_Pi_K_P.png">Efficiency #pi, K, p</a><br>  //-->
<!--  <a href="eff_Pi_K_P_findable.png">Efficiency findable #pi, K, p</a><br>  //-->
</div>

EOF
}

makeHTMLindexPerPeriod()
{
  cat > index.html <<EOF
<div align="left"><br>
<h2>Periodical Data Quality</h2>
<a href="meanTPCncl_vs_run.png">Mean Number of TPC Clusters</a><br>
<a href="meanTPCnclF_vs_run.png"># of Found Clusters / # of Findable Clusters</a><br>
<br>
<a href="meanMIP_vs_run.png">Mean of MIPs</a><br>
<a href="meandEdxele_vs_run.png">Mean electron energy loss p(0.32,0.38)GeV/c, dEdx(70,100)</a><br>
<br>
<a href="meanVertX_vs_run.png">Mean of Vertex_X</a>, <a href="meanVertY_vs_run.png">Vertex_Y</a>, <a href="meanVertZ_vs_run.png">Vertex_Z</a>
<br><br>
<a href="meanMult_vs_run.png">Multiplicities of Primary Tracks</a><br>
<a href="TPC-ITS-matching-efficiency_vs_run.png">TPC-ITS matching efficiency</a><br>
<a href="ITS-TPC-matching-quality_vs_run.png">ITS-TPC matching quality</a><br>

<a href="1overPt_vs_run.png">Delta 1/pt</a><br>
<a href="DCAOffset_vs_run.png">DCAs</a><br>
<br>
<p><font size="4">Runs:</font></p>
EOF

local dir
for dir in 000*; do
  echo "<a href="${dir}">${dir}</a>" >> index.html
done

  cat >> index.html <<EOF
<br>
<p><font size="4">Additional plots</font></p>
<a href="occ_AC_Side_IROC_OROC_vs_run.png">Nr of Chambers with lower gain (occupancy)</a><br>
<br>
<a href="dcar_fitting_run.png">DCAr fitting parameters</a><br>
<a href="dcar_0_vs_run.png">DCAr fitting parameters (0)</a><br>
<a href="dcar_1_vs_run.png">DCAr fitting parameters (1)</a><br>
<a href="dcar_2_vs_run.png">DCAr fitting parameters (2)</a><br>
<a href="dcaz_0_vs_run.png">DCAz fitting parameters (0)</a><br>
<a href="dcaz_1_vs_run.png">DCAz fitting parameters (1)</a><br>
<a href="dcaz_2_vs_run.png">DCAz fitting parameters (2)</a><br>
<br>
<a href="resolutionMIP_vs_run.png">Resolution of MIPs</a><br>
<a href="resolutionMeandEdxEle_vs_run.png">Resolution of mean electron energy loss</a><br>
<a href="ElectroMIPSeparation_vs_run.png">Separation between electon and MIPs energy loss</a><br>
<br>
<a href="MIPattachSlopes_vs_run.png">Attachment parameter p1, A side and C side</a><br>
<br>
<a href="pullPhiConstrain_vs_run.png">Tracking parameter phi</a><br>
<br>
<br>
<a href="prodinfo">Production information</a><br>

</div>

EOF
}
