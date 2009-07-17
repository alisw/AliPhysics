#ifndef __CINT__
#endif
Bool_t RunOnProof(TString macro,Long64_t numberOfEvents, Long64_t skipEvents)
{

  fgMode = "proof";

  // adds standard ANALYSIS parfiles
  fgPARS = "STEERBase:ESD:AOD:ANALYSIS:ANALYSISalice";

  // adds standard RESONANCE parfiles
  fgPARS += ":PWG2resonances";

  fgPARSClean = "";
//   fgPARSClean = "all";

  fgUser = gSystem->ExpandPathName("$USER");

//     fReset="RESET";
//   fgProofToConnect = Form("%s@localhost:2093",fgUser.Data());
//     fgProofToConnect = Form("%s@localhost:3093",fgUser.Data());
  fgProofToConnect = Form("%s@alicecaf.cern.ch",fgUser.Data());
//    fgProofToConnect = Form("%s@alicepc100.jinr.ru",fgUser.Data());
//     fgPARSClean = "PWG2resonances:PWG2resonancesTEST";
//     fgPARSClean="ALL";

  fgDataType = (Int_t) AliRsnUtils::kDataSet;
  fgInputFileName = "/COMMON/COMMON/LHC08c12_0.9TeV_0.5T";
//     fgInputFileName = "/COMMON/COMMON/LHC08c11_10TeV_0.5T";
//   fgInputFileName = "/COMMON/COMMON/LHC09a4_10TeV";

//    fgDataType = (Int_t) AliRsnUtils::kTxt;
//    fgInputFileName = "ESD_alicepc100_PDC08_test_ssd.txt";
//    fgInputFileName = "ESD_alicepc100_PDC08_test_hdd.txt";


  // sets up macro
  fgMacro = macro;

  // number of events
  fgNumOfEvents = numberOfEvents;

  // number of events to skip
  fgNumOfEventsSkip = skipEvents;


  return runProof();
}

Bool_t RunLocaly(TString macro,Long64_t numberOfEvents, Long64_t skipEvents)
{

  fgMode = "local";
//   fgUseLocalLibs = kTRUE;
  fgAlirootLibPath = "$ALICE_ROOT/lib/tgt_$ALICE_TARGET";
//   fgAlirootLibPath = "$ALIMV/$ALICE_TARGET/lib";

  // adds standard ANALYSIS parfiles
  fgPARS = "STEERBase:ESD:AOD:ANALYSIS:ANALYSISalice";
  if (fgUseLocalLibs) {
//     Info("","Loading AliRoot Libs ...");
//     gROOT->Macro(Form("%s/macros/loadlibs.C",gSystem->ExpandPathName("$ALICE_ROOT/")));
//     Info("","Loading AliRoot Libs is done...");
    fgPARS = "";
    fgLIBS = "STEERBase:ESD:AOD";
//     fgLIBS += ":ANALYSIS:ANALYSISalice";
    fgPARS += "ANALYSIS:ANALYSISalice";
  }
  // adds standard RESONANCE parfiles
  fgPARS += ":PWG2resonances";

  // clean packages
//   fgPARSClean = "all";
//   fgPARSClean = "ANALYSISTest:PWG2resonances:PWG2resonancesMV";

  // input
  fgDataType = (Int_t) AliRsnUtils::kTxt;
  fgInputFileName = "ESD.txt";
  fgTreeName = "esdTree";
//     fgInputFileName = "AOD_vala_PDC08.txt";
//     fgTreeName = "aodTree";
//     fgDoMixing = kTRUE;

  // sets up macro
  fgMacro = macro;
  // number of events
  fgNumOfEvents = numberOfEvents;
  // number of events to skip
  fgNumOfEventsSkip = skipEvents;

  return runLocal();

}

Bool_t RunOnAliEn()
{

  TGrid::Connect("alien://");

  fgMode = "alien";

//   fgAlienShouldRun = kTRUE;
//   fgAlienShoudlCopy = kTRUE;

  // adds standard ANALYSIS parfiles
  fgPARS = "STEERBase:ESD:AOD:ANALYSIS:ANALYSISalice";

  // adds standard RESONANCE parfiles
  fgPARS += ":PWG2resonances";

  fgPARSClean = "";

  TString productionName = "PDC07f";
  Int_t tagCollNum = 160009;
  if (num>0)
    tagCollNum = num;
  fgAlienSplit=50;
  fgProjectDir = Form("/alice/cern.ch/user/m/mvala/RSNTASK/ANALYSIS/%s_50/ALL",productionName.Data());
  //     fProjectDir = Form("/alice/cern.ch/user/m/mvala/RSNTASK/ANALYSIS/%s/%d",productionName.Data(),tagCollNum);
  fgOutputDir = Form("/alice/cern.ch/user/m/mvala/DATA/RSN/%s_50/%d",productionName.Data(),tagCollNum);


//   fgDataType = AliRsnUtils::kXmlCollection;
//   fgCollName = Form("/alice/cern.ch/user/m/mvala/DATA/ESD/collections/PDC07f/%d.xml",tagCollNum);
  fgDataType = AliRsnUtils::kXmlCollectionTag;
  fgCollName = Form("/alice/cern.ch/user/m/mvala/DATA/ESD/collections/PDC07f/%d_tag.xml",tagCollNum);
  fgInputFileName = Form("alien://%s",fCollName.Data());


  // extra files
  fgExtraInputFiles = "PWG2resonancesUtils.C:runRSNFilterAdvanced.C:runPhiAnalysisAll.C:runPIDComparison.C";
  fgExtraInputFiles += ":runRSNProcessInfo.C:RsnConfig_PHIKK.C";

  fgJDL = (TAlienJDL*) gGrid->GetJDLGenerator();

  // jdl part (you can use mini.jdl or set up by yourself
  // utils->GetJDL()->... [take look in TAlienJDL in $ROOTSYS/net/alien/])
//   utils->GetJDL()->Parse ( "mini.jdl" );
  fgJDL->SetExecutable("root.sh");
//   fgJDL->SetExecutable("aliroot.sh");
  fgJDL->AddToPackages("VO_ALICE@APISCONFIG::V2.4");
  fgJDL->AddToPackages("VO_ALICE@ROOT::v5-21-01-alice");
//   fgJDL->AddToPackages("VO_ALICE@AliRoot::v4-15-Rev-04");
  fgJDL->AddToOutputArchive("log_archive.zip:stdout,stderr@ALICE::CERN::SE");
  fgJDL->AddToOutputArchive("root_archive.zip:*.root@ALICE::CERN::SE");
//   fgJDL->SetOutputDirectory(Form("%s/",pname.Data()));
  fgJDL->SetOutputDirectory(Form("%s/#alien_counter_03i#",fgOutputDir.Data()));
  fgJDL->SetSplitMode("se", fgAlienSplit);
  fgJDL->SetInputDataList("wn.xml");
  fgJDL->SetInputDataListFormat("xml-single");
//     fgJDL->SetInputDataListFormat(Form("merge:%s",fCollName.Data()));

  fgJDL->AddToInputDataCollection(Form("LF:%s,nodownload",fgCollName.Data()));
//   fgJDL->AddToMerge("histOut.root:/alice/jdl/mergerootfile.jdl:histOut.Merged.root");
//   fgJDL->SetMergedOutputDirectory(Form("%s",pname.Data()));

  return runAlien();
}
