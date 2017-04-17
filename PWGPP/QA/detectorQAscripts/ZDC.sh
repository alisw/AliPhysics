runLevelQA()
{
  qaFile=$1

  cp $ALICE_PHYSICS/PWGPP/ZDC/macros/MakeTrendZDC.C .
  aliroot -b -q -l "MakeTrendZDC.C(\"$qaFile\",$runNumber)"

  cp $ALICE_PHYSICS/PWGPP/ZDC/macros/DrawPerformanceZDCQAMatch.C .
  aliroot -b -q -l "DrawPerformanceZDCQAMatch.C(\"trending.root\")"
  makeHTMLindexPerRun
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGPP/ZDC/macros/DrawPerformanceZDCQAMatchTrends.C .
  aliroot -b -q -l "DrawPerformanceZDCQAMatchTrends.C(\"prodQAhistos.root\")"
  makeHTMLindexPerPeriod
}

makeHTMLindexPerRun()
{
  cat > index.html<<EOF
  <div align="left"><br>
  <h1>ZDC Data Quality (single run checks)</h1>
  <br>
  <h2>Signals: spectra</h2>
  </div>
  <a href="cZNA_Spectra.png">ZNA spectrum</a><br>
  <a href="cZNC_Spectra.png">ZNC spectrum</a><br>
  <a href="cZPA_Spectra.png">ZPA spectrum</a><br>
  <a href="cZPC_Spectra.png">ZPC spectrum</a><br>
  <a href="cZEM1_Spectra.png">ZEM1 spectrum</a><br>
  <a href="cZEM2_Spectra.png">ZEM2 spectrum</a><br>
  <h2>Signals: mean values</h2>
  </div>
  <a href="cZNA_Mean_Values.png">ZNA mean signal</a><br>
  <a href="cZNC_Mean_Values.png">ZNC mean signal</a><br>
  <a href="cZPA_Mean_Values.png">ZPA mean signal</a><br>
  <a href="cZPC_Mean_Values.png">ZPC mean signal</a><br>
  <a href="cZEM1_Mean_Values.png">ZEM1 mean signal</a><br>
  <a href="cZEM2_Mean_Values.png">ZEM2 mean signal</a><br>
  <h2>Uncalibrated signals: spectra</h2>
  </div>
  <a href="cZNA_Spectra_Uncal.png">ZNA spectrum</a><br>
  <a href="cZNC_Spectra_Uncal.png">ZNC spectrum</a><br>
  <a href="cZPA_Spectra_Uncal.png">ZPA spectrum</a><br>
  <a href="cZPC_Spectra_Uncal.png">ZPC spectrum</a><br>
  <h2>Uncalibrated signals: mean values</h2>
  </div>
  <a href="cZNA_Mean_Uncalib.png">ZNA mean signal</a><br>
  <a href="cZNC_Mean_Uncalib.png">ZNC mean signal</a><br>
  <a href="cZPA_Mean_Uncalib.png">ZPA mean signal</a><br>
  <a href="cZPC_Mean_Uncalib.png">ZPC mean signal</a><br>
  <a href="cZEM1_Mean_Uncalib.png">ZEM1 mean signal</a><br>
  <a href="cZEM2_Mean_Uncalib.png">ZEM2 mean signal</a><br>
  <h2>Centroids</h2>
  </div>
  <a href="cZNA_X_centroid.png">ZNA centroid X coordinate</a><br>
  <a href="cZNA_Y_centroid.png">ZNA centroid Y coordinate</a><br>
  <a href="cZNC_X_centroid.png">ZNC centroid X coordinate</a><br>
  <a href="cZNC_Y_centroid.png">ZNC centroid Y coordinate</a><br>
  <h2>Timing</h2>
  </div>
  <a href="cTimingSum.png">Timing (ZNC TDC + ZNA TDC)</a><br>
  <a href="cTimingDiff.png">Timing (ZNC TDC - ZNA TDC)</a><br>
  </div>
EOF
}

makeHTMLindexPerPeriod()
{
  cat > index.html<<EOF
  <div align="left"><br>
  <h1>ZDC Data Quality (trending plots)</h1>
  <br>
  <h2>Signals: mean values</h2>
  </div>
  <a href="ZN_signals_trending.png">ZNA,ZNC signals trending</a><br>
  <a href="ZP_signals_trending.png">ZPA,ZPC signals trending</a><br>
  <a href="ZN_signals_trending.png">ZNA,ZNC uncalibrated signals trending</a><br>
  <a href="ZP_signals_trending.png">ZPA,ZPC uncalibrated signals trending</a><br>
  <a href="ZEM_signals_trending.png">ZEM1,ZEM2 signals trending</a><br>
  <h2>Centroids</h2>
  </div>
  <a href="ZNA_centroids_trending.png">ZNA centroids trending</a><br>
  <a href="ZNC_centroids_trending.png">ZNC centroids trending</a><br>
  <h2>Timing</h2>
  </div>
  <a href="ZN_timing_trending.png">ZN timing trending</a><br>
  <br><br>
  <h4>Analyzed runs:</h4>
  </div>
EOF

local dir
for dir in 000*; do
  echo "<a href="${dir}">${dir}</a>" >> index.html
done

}
