#available variables:
#  $dataType     e.g. data or sim
#  $year         e.g. 2011
#  $period       e.g. LHC13g
#  $runNumber     e.g. 169123
#  $pass         e.g. cpass1,pass1,passMC
#  #ocdbStorage  e.g. "raw://", "local://./OCDB"

runLevelQA()
{
  qaFile=$1
  export TOFqaFile=$qaFile
  export TOFrunNumber=$runNumber
  export TOFocdbStorage=$ocdbStorage
  cp $ALICE_PHYSICS/PWGPP/TOF/trending/MakeTrendingTOFQAv2.C .
aliroot -l -b << EOF
  gSystem->AddIncludePath("-I${ALICE_ROOT}/TOF ");
  gSystem->AddIncludePath("-I${ALICE_ROOT}/include ");
  gSystem->Load("libTOFbase");
  .L MakeTrendingTOFQAv2.C+g
  const TString qaFile = gSystem->Getenv("TOFqaFile");
  const TString runNumber = gSystem->Getenv("TOFrunNumber");
  const TString ocdbStorage = gSystem->Getenv("TOFocdbStorage");
  Printf("qaFile %s, runNumber %s, ocdbStorage %s", qaFile.Data(), runNumber.Data(), ocdbStorage.Data());
  MakeTrendingTOFQAv2(qaFile, runNumber.Atoi(), "", 0, 1, 1, -2., 2., 10., 100., ocdbStorage, 1, 1, 1, 1)
EOF
  rm MakeTrendingTOFQAv2_C.so
  rm MakeTrendingTOFQAv2_C.d
  unset TOFqaFile
  unset TOFrunNumber
  unset TOFocdbStorage
}

periodLevelQA()
{
  trendingFile=$1

  cp $ALICE_PHYSICS/PWGPP/TOF/trending/DrawTrendingTOFQA.C .
  aliroot -b -q -l "DrawTrendingTOFQA.C+g(\"trending.root\")"
  rm DrawTrendingTOFQA_C.so
  rm DrawTrendingTOFQA_C.d
}
