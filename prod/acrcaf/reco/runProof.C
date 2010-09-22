void runProof(Int_t runNumber, Int_t nev = 10000, Int_t firstev = 0,TString proofCluster="aliprod@alice-caf.cern.ch",TString workers="")
{
  gEnv->SetValue("XSec.GSI.DelegProxy","2");

  // Select ROOT version
  TProof::Mgr(proofCluster.Data())->SetROOTVersion(Form("VO_ALICE@ROOT::%s",gSystem->GetFromPipe("echo $ROOTSYS | awk -F/ '{print $NF}'").Data()));

  // set ALIROOT mode to REC
  TList * list = new TList();
  list->Add(new TNamed("ALIROOT_MODE", "REC"));    
  list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
  // Login to CAF

  TString workersStr;
  if (!workers.IsNull()) workersStr=Form("workers=%s",workers.Data());
  
  // connecting to proof cluster
  TProof::Open(proofCluster.Data(),workersStr.Data());

  // seting up aliroot version
  TString alirootVer = Form("VO_ALICE@AliRoot::%s",gSystem->GetFromPipe("echo $ALICE_ROOT | awk -F/ '{print $NF}'").Data());
  
  // Enable AliRoot
  Printf("Setting aliroot version on caf to %s ...",alirootVer.Data());
  if (gProof->EnablePackage(alirootVer.Data(), list)) {
    Error("run.C",Form("Error enabling proof package %s !!!!",alirootVer.Data()));
    return;
  }


  gSystem->Load("libMonaLisa.so");
  TMonaLisaWriter monalisa("pcalishuttle02.cern.ch","SHIFTER_RECO_CAF");
  SendMonaLisaData(&monalisa,runNumber,"Started",0);

  // Temporary fix in order to avoid timeouts on the master
  Int_t numWorkers = gProof->GetParallel();
  gProof->SetParameter("PROOF_PacketizerStrategy", (Int_t)0);
  gProof->SetParameter("PROOF_PacketAsAFraction", (Int_t) nev/numWorkers);
//  gProof->SetParameter("PROOF_MinPacketTime", 8);
//  gProof->SetParameter("PROOF_MaxPacketTime", 5);

  // Set some verbosity
  // gProof->SetLogLevel(3); 


  // Run reconstruction
  gROOT->LoadMacro("rec.C");
  gROOT->ProcessLine(Form("rec(%d,%d,%d);",runNumber,nev,firstev));

  TProof::Mgr(proofCluster.Data())->GetSessionLogs()->Save("*",Form("log/run%d.log",runNumber));

  // Check the produced dataset
  TFileCollection *coll = gProof->GetDataSet(Form("run%d",runNumber));
  if (coll) {
    Int_t nEvents = coll->GetTotalEntries("/esdTree");
    if (nEvents > 0) {
      cout << "===========================================================================" << endl;
      cout << nEvents << " events reconstructed and stored in the dataset run" << runNumber << endl;
      cout << "===========================================================================" << endl;
      cout << "The dataset is:" << endl;
      coll->Print();
      cout << "===========================================================================" << endl;
      SendMonaLisaData(&monalisa,runNumber,"Done",nEvents);
    }
    else {
      SendMonaLisaData(&monalisa,runNumber,"No_Events",nEvents);
    }
  }
  else {
      SendMonaLisaData(&monalisa,runNumber,"Failure",0);
  }
}

void SendMonaLisaData(TMonaLisaWriter *monalisa, Int_t runNumber, const char* status, Int_t nEvents)
{
  TMonaLisaText mlStatus("Status",status);
  TMonaLisaValue mlEventCount("Event_count",nEvents);
  TList mlList;
  mlList.Add(&mlStatus);
  mlList.Add(&mlEventCount);
  TString mlID;
  mlID.Form("%d",runNumber);
  monalisa->SendParameters(&mlList, mlID.Data());
}

void LoadRecMacroOnClient() {
  gSystem->AddIncludePath(Form("-I\"%s/include\"", gSystem->Getenv("ALICE_ROOT")));
  gROOT->Macro("$ALICE_ROOT/macros/loadlibsrec.C");
}
