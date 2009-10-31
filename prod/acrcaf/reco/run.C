void run(Int_t runNumber,Int_t nev=10000, Int_t firstev=0)
{
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  // Select ROOT version
  TProof::Mgr("aliprod@alicecaf")->SetROOTVersion("v5-24-00b-caf");
  // Login to CAF
  TProof::Open("aliprod@alicecaf");

  // Enable AliRoot
  gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-17-Release.rec/AF-v4-17-rec.par");
  gProof->EnablePackage("AF-v4-17-rec.par");

  gSystem->Load("libMonaLisa.so");
  TMonaLisaWriter monalisa("pcalishuttle02.cern.ch",
			   "SHIFTER_RECO_CAF");
  SendMonaLisaData(&monalisa,runNumber,"Started",0);

  // Temporary fix in order to avoid timeouts on the master
  gProof->SetParameter("PROOF_PacketAsAFraction",20);

  // Run reconstruction
  gROOT->LoadMacro("rec.C");
  gROOT->ProcessLine(Form("rec(%d,%d,%d);",runNumber,nev,firstev));

  TProof::Mgr("aliprod@alicecaf")->GetSessionLogs()->Save("*",Form("log/run%d.log",runNumber));

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
