void initHists(Int_t fRunNo, Int_t fPass, TList* list); 
void fill()
{
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libOADB");

  // LHC10h parameters
  Float_t fV0MScaleFactor   = 1.0;
  Float_t fSPDScaleFactor   = 1.0;
  Float_t fTPCScaleFactor   = 1.0;
  Float_t fV0MScaleFactorMC = 0.8;
  Float_t fV0MSPDOutlierPar0 = -0.579579;
  Float_t fV0MSPDOutlierPar1 = 0.273949;
  Float_t fV0MTPCOutlierPar0 = -1.03873;
  Float_t fV0MTPCOutlierPar1 = 0.125691;
  Float_t fV0MSPDSigmaOutlierPar0  = 190.3;
  Float_t fV0MSPDSigmaOutlierPar1  = -3.3;
  Float_t fV0MSPDSigmaOutlierPar2  = 0.015;
  Float_t fV0MTPCSigmaOutlierPar0  = 96.8;
  Float_t fV0MTPCSigmaOutlierPar1  = 1.7;
  Float_t fV0MTPCSigmaOutlierPar2  = 0.008;
  Float_t fV0MZDCOutlierPar0       = 6200.0;
  Float_t fV0MZDCOutlierPar1       = - 0.25 ;
  Float_t fV0MZDCEcalOutlierPar0   = 235000.;
  Float_t fV0MZDCEcalOutlierPar1   = - 9.5;
  Float_t fZVCut = 10.0;
  Float_t fOutliersCut = 6.0;
  Bool_t fUseScaling  = kFALSE;
  Bool_t fUseCleaning = kTRUE;
  
  
  AliOADBContainer* con     = new AliOADBContainer("Centrality");

  
  AliOADBCentrality*  cent1  = new AliOADBCentrality();
  cent1->SetScaleFactors(fV0MScaleFactor,fSPDScaleFactor,fTPCScaleFactor,fV0MScaleFactorMC);
  cent1->SetOutlierV0MSPDFactors(fV0MSPDOutlierPar0,fV0MSPDOutlierPar1,fV0MSPDSigmaOutlierPar0,fV0MSPDSigmaOutlierPar1,fV0MSPDSigmaOutlierPar2);
  cent1->SetOutlierV0MTPCFactors(fV0MTPCOutlierPar0,fV0MTPCOutlierPar1,fV0MTPCSigmaOutlierPar0,fV0MTPCSigmaOutlierPar1,fV0MTPCSigmaOutlierPar2);
  cent1->SetOutlierV0MZDCFactors(fV0MZDCOutlierPar0,fV0MZDCOutlierPar1);
  cent1->SetOutlierV0MZDCEcalFactors(fV0MZDCEcalOutlierPar0,fV0MZDCEcalOutlierPar1);
  cent1->SetZVCut(fZVCut);
  cent1->SetOutliersCut(fOutliersCut);
  cent1->SetUseScaling(fUseScaling);
  cent1->SetUseCleaning(fUseCleaning);
  TList* list1  = new TList();
  TList* list12  = new TList();
  list1->SetName("AliCentralityBy1D_137161");
  list12->SetName("AliCentralityByFunction_137161");
  initHists(137161, 2, list1, list12); 
  cent1->SetHistReferences(list1, list12);
  con->AppendObject(cent1,  136851,  137165);
  
  AliOADBCentrality*  cent2  = new AliOADBCentrality();
  cent2->SetScaleFactors(fV0MScaleFactor,fSPDScaleFactor,fTPCScaleFactor,fV0MScaleFactorMC);
  cent2->SetOutlierV0MSPDFactors(fV0MSPDOutlierPar0,fV0MSPDOutlierPar1,fV0MSPDSigmaOutlierPar0,fV0MSPDSigmaOutlierPar1,fV0MSPDSigmaOutlierPar2);
  cent2->SetOutlierV0MTPCFactors(fV0MTPCOutlierPar0,fV0MTPCOutlierPar1,fV0MTPCSigmaOutlierPar0,fV0MTPCSigmaOutlierPar1,fV0MTPCSigmaOutlierPar2);
  cent2->SetOutlierV0MZDCFactors(fV0MZDCOutlierPar0,fV0MZDCOutlierPar1);
  cent2->SetOutlierV0MZDCEcalFactors(fV0MZDCEcalOutlierPar0,fV0MZDCEcalOutlierPar1);
  cent2->SetZVCut(fZVCut);
  cent2->SetOutliersCut(fOutliersCut);
  cent2->SetUseScaling(fUseScaling);
  cent2->SetUseCleaning(fUseCleaning);
  TList* list2  = new TList();
  TList* list22  = new TList();
  list2->SetName("AliCentralityBy1D_137366");
  list22->SetName("AliCentralityByFunction_137366");
  initHists(137366, 2, list2, list22); 
  cent2->SetHistReferences(list2, list22);
  con->AppendObject(cent2,  137230,  137531);
  
  
  AliOADBCentrality*  cent3  = new AliOADBCentrality();
  cent3->SetScaleFactors(fV0MScaleFactor,fSPDScaleFactor,fTPCScaleFactor,fV0MScaleFactorMC);
  cent3->SetOutlierV0MSPDFactors(fV0MSPDOutlierPar0,fV0MSPDOutlierPar1,fV0MSPDSigmaOutlierPar0,fV0MSPDSigmaOutlierPar1,fV0MSPDSigmaOutlierPar2);
  cent3->SetOutlierV0MTPCFactors(fV0MTPCOutlierPar0,fV0MTPCOutlierPar1,fV0MTPCSigmaOutlierPar0,fV0MTPCSigmaOutlierPar1,fV0MTPCSigmaOutlierPar2);
  cent3->SetOutlierV0MZDCFactors(fV0MZDCOutlierPar0,fV0MZDCOutlierPar1);
  cent3->SetOutlierV0MZDCEcalFactors(fV0MZDCEcalOutlierPar0,fV0MZDCEcalOutlierPar1);
  cent3->SetZVCut(fZVCut);
  cent3->SetOutliersCut(fOutliersCut);
  cent3->SetUseScaling(fUseScaling);
  cent3->SetUseCleaning(fUseCleaning);
  TList* list3  = new TList();
  TList* list32  = new TList();
  list3->SetName("AliCentralityBy1D_137722");
  list32->SetName("AliCentralityByFunction_137722");
  initHists(137722, 2, list3, list32); 
  cent3->SetHistReferences(list3, list32);
  con->AppendObject(cent3, 137539, 137848);
  
  
  AliOADBCentrality*  cent4  = new AliOADBCentrality();
  cent4->SetScaleFactors(fV0MScaleFactor,fSPDScaleFactor,fTPCScaleFactor,fV0MScaleFactorMC);
  cent4->SetOutlierV0MSPDFactors(fV0MSPDOutlierPar0,fV0MSPDOutlierPar1,fV0MSPDSigmaOutlierPar0,fV0MSPDSigmaOutlierPar1,fV0MSPDSigmaOutlierPar2);
  cent4->SetOutlierV0MTPCFactors(fV0MTPCOutlierPar0,fV0MTPCOutlierPar1,fV0MTPCSigmaOutlierPar0,fV0MTPCSigmaOutlierPar1,fV0MTPCSigmaOutlierPar2);
  cent4->SetOutlierV0MZDCFactors(fV0MZDCOutlierPar0,fV0MZDCOutlierPar1);
  cent4->SetOutlierV0MZDCEcalFactors(fV0MZDCEcalOutlierPar0,fV0MZDCEcalOutlierPar1);
  cent4->SetZVCut(fZVCut);
  cent4->SetOutliersCut(fOutliersCut);
  cent4->SetUseScaling(fUseScaling);
  cent4->SetUseCleaning(fUseCleaning);
  TList* list4  = new TList();
  TList* list42  = new TList();
  list4->SetName("AliCentralityBy1D_138150");
  list42->SetName("AliCentralityByFunction_138150");
  initHists(138150, 2, list4, list42); 
  cent4->SetHistReferences(list4, list42);
  con->AppendObject(cent4,  138125, 138154);


  AliOADBCentrality*  cent5  = new AliOADBCentrality();
  cent5->SetScaleFactors(fV0MScaleFactor,fSPDScaleFactor,fTPCScaleFactor,fV0MScaleFactorMC);
  cent5->SetOutlierV0MSPDFactors(fV0MSPDOutlierPar0,fV0MSPDOutlierPar1,fV0MSPDSigmaOutlierPar0,fV0MSPDSigmaOutlierPar1,fV0MSPDSigmaOutlierPar2);
  cent5->SetOutlierV0MTPCFactors(fV0MTPCOutlierPar0,fV0MTPCOutlierPar1,fV0MTPCSigmaOutlierPar0,fV0MTPCSigmaOutlierPar1,fV0MTPCSigmaOutlierPar2);
  cent5->SetOutlierV0MZDCFactors(fV0MZDCOutlierPar0,fV0MZDCOutlierPar1);
  cent5->SetOutlierV0MZDCEcalFactors(fV0MZDCEcalOutlierPar0,fV0MZDCEcalOutlierPar1);
  cent5->SetZVCut(fZVCut);
			   cent5->SetOutliersCut(fOutliersCut);
  cent5->SetUseScaling(fUseScaling);
  cent5->SetUseCleaning(fUseCleaning);
  TList* list5  = new TList();
  TList* list52  = new TList();
  list5->SetName("AliCentralityBy1D_138200");
  list52->SetName("AliCentralityByFunction_138200");
  initHists(138200, 2, list5, list52); 
  cent5->SetHistReferences(list5, list52);
  con->AppendObject(cent5,  138190, 138275);


  AliOADBCentrality*  cent6  = new AliOADBCentrality();
  cent6->SetScaleFactors(fV0MScaleFactor,fSPDScaleFactor,fTPCScaleFactor,fV0MScaleFactorMC);
  cent6->SetOutlierV0MSPDFactors(fV0MSPDOutlierPar0,fV0MSPDOutlierPar1,fV0MSPDSigmaOutlierPar0,fV0MSPDSigmaOutlierPar1,fV0MSPDSigmaOutlierPar2);
  cent6->SetOutlierV0MTPCFactors(fV0MTPCOutlierPar0,fV0MTPCOutlierPar1,fV0MTPCSigmaOutlierPar0,fV0MTPCSigmaOutlierPar1,fV0MTPCSigmaOutlierPar2);
  cent6->SetOutlierV0MZDCFactors(fV0MZDCOutlierPar0,fV0MZDCOutlierPar1);
  cent6->SetOutlierV0MZDCEcalFactors(fV0MZDCEcalOutlierPar0,fV0MZDCEcalOutlierPar1);
  cent6->SetZVCut(fZVCut);
			   cent6->SetOutliersCut(fOutliersCut);
  cent6->SetUseScaling(fUseScaling);
  cent6->SetUseCleaning(fUseCleaning);
  TList* list6  = new TList();
  TList* list62  = new TList();
  list6->SetName("AliCentralityBy1D_139172");
  list62->SetName("AliCentralityByFunction_139172");
  initHists(139172, 2, list6, list62); 
  cent6->SetHistReferences(list6, list62);
  con->AppendObject(cent6,  138358,  139517);  
  

  AliOADBCentrality*  defaultCent  = new AliOADBCentrality("oadbDefault");
  defaultCent->SetScaleFactors(fV0MScaleFactor,fSPDScaleFactor,fTPCScaleFactor,fV0MScaleFactorMC);
  defaultCent->SetOutlierV0MSPDFactors(fV0MSPDOutlierPar0,fV0MSPDOutlierPar1,fV0MSPDSigmaOutlierPar0,fV0MSPDSigmaOutlierPar1,fV0MSPDSigmaOutlierPar2);
  defaultCent->SetOutlierV0MTPCFactors(fV0MTPCOutlierPar0,fV0MTPCOutlierPar1,fV0MTPCSigmaOutlierPar0,fV0MTPCSigmaOutlierPar1,fV0MTPCSigmaOutlierPar2);
  defaultCent->SetOutlierV0MZDCFactors(fV0MZDCOutlierPar0,fV0MZDCOutlierPar1);
  defaultCent->SetOutlierV0MZDCEcalFactors(fV0MZDCEcalOutlierPar0,fV0MZDCEcalOutlierPar1);
  defaultCent->SetZVCut(fZVCut);
  defaultCent->SetOutliersCut(fOutliersCut);
  defaultCent->SetUseScaling(fUseScaling);
  defaultCent->SetUseCleaning(kFALSE);
  defaultCent->SetHistReferences(list6, list62);
  con->AddDefaultObject(defaultCent);

  
  con->WriteToFile("centrality.root");
}

void initHists(Int_t fRunNo, Int_t fPass, TList* list1, TList* list2) 
{
  TString fileName =(Form("%s/COMMON/CENTRALITY/data/pass%d/AliCentralityBy1D_%d.root", AliOADBContainer::GetOADBPath(), fPass, fRunNo));
  TString fileName2=(Form("%s/COMMON/CENTRALITY/data/pass%d/AliCentralityByFunction_%d.root", AliOADBContainer::GetOADBPath(), fPass, fRunNo));
  TDirectory *owd = gDirectory;
  // Check if the file is present
  TString path = gSystem->ExpandPathName(fileName.Data());
  if (gSystem->AccessPathName(path)) {
     AliError(Form("File %s does not exist", path.Data()));
     return;
  }
  fFile  = TFile::Open(fileName);
  fFile->ls();
  owd->cd();

  fHtempV0M  = (TH1F*) (fFile->Get("hmultV0_percentile"));
  fHtempTRK  = (TH1F*) (fFile->Get("hNtracks_percentile"));
  //fHtempTKL  = (TH1F*) (fFile->Get("hNtracklets_percentile"));
  //fHtempCL0  = (TH1F*) (fFile->Get("hNclusters0_percentile"));
  fHtempCL1  = (TH1F*) (fFile->Get("hNclusters1_percentile"));

  fHtempV0M  ->SetName("fHOutMultV0M_percentile");
  fHtempTRK  ->SetName("fHOutMultTRK_percentile");
  fHtempCL1  ->SetName("fHOutMultCL1_percentile");
  fHtempV0M  ->SetTitle("fHOutMultV0M_percentile");
  fHtempTRK  ->SetTitle("fHOutMultTRK_percentile");
  fHtempCL1  ->SetTitle("fHOutMultCL1_percentile");

  list1->Add(fHtempV0M);
  list1->Add(fHtempTRK);
  //list->Add(fHtempTKL);
  //list->Add(fHtempCL0);
  list1->Add(fHtempCL1);
  list1->ls();


  TDirectory *owd = gDirectory;
  TString path = gSystem->ExpandPathName(fileName2.Data());
  if (gSystem->AccessPathName(path)) {
     AliError(Form("File %s does not exist", path.Data()));
     return;
  }   
  fFile2  = TFile::Open(fileName2);
  owd->cd();
  fFile2->ls();
  //  fHtempV0MvsFMD =  (TH1F*) (fFile2->Get("hmultV0vsmultFMD_all_percentile"));
  //  fHtempTKLvsV0M  = (TH1F*) (fFile2->Get("hNtrackletsvsmultV0_all_percentile"));
  fHtempZEMvsZDC  = (TH2F*) (fFile2->Get("hEzemvsEzdc_all_percentile"));
  
  fHtempZEMvsZDC->SetName("fHOutMultZEMvsZDC");
  fHtempZEMvsZDC->SetTitle("fHOutMultZEMvsZDC");

  // list2->Add(fHtempV0MvsFMD);
  // list2->Add(fHtempTKLvsV0M);
  list2->Add(fHtempZEMvsZDC);
  owd->cd();

}
