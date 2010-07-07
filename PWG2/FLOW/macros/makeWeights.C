//=======================================================================//
// Macro makeWeights.C is used to make phi, pt and eta weights. 
// Before using the macro makeWeights.C you should already have 
// available the output .root files from various methods stored 
// in the common output file "AnalysisResults.root". 
// When calling this macro you must specify the analysis type 
// and the method whose output file you would like to use 
// to make the weights for the subsequent runs (for cumulants, GFC and QC,
// you must also specify the cumulant order).
//=======================================================================//

//========================= phi-weights ========================
// Phi-weights are obtained by inverting and normalizing the
// azimuthal acceptance profile. This procedure isn't applicable 
// if the detector has a gap in azimuthal acceptance (i.e. if 
// there exists a phi bin with no entries in the histogram for
// detector's azimuthal acceptance profile).
//========================= pt-weights =========================
// You can make pt-weights in three different ways:
// 1.) pt-weights are growing linearly as a function of pt for 
//     pt <= ptCutoff. For pt > ptCutoff the pt-weights are 
//     constant and equal to ptMax. To enable this option set
//     useLinearPtWeights = kTRUE;
// 2.) pt-weights are growing quadratically as a function of pt
//     for pt <= ptCutoff. For pt > ptCutoff the pt-weights are 
//     constant and equal to ptMax. To enable this option set
//     useQadraticPtWeights = kTRUE;
// 3.) pt-weights are simply v(pt) from the specified method's 
//     result for differential flow vs pt. To enable this option 
//     set useLinearPtWeights = kFALSE and useQadraticPtWeights = kFALSE.
Double_t ptCutoff = 2; // in GeV
Double_t ptWeightsMax = 0.2; 
Bool_t useLinearPtWeights = kTRUE;
Bool_t useQuadraticPtWeights = kFALSE;
//========================= eta-weights ========================
// Eta-weights are simply v(eta) from the specified method's
// result for differential flow vs eta.
//==============================================================

enum libModes {mLocal,mLocalSource};

//void makeWeights(TString type="", TString method="GFC", TString cumulantOrder="4th", Int_t mode=mLocalSource)
void makeWeights(TString type="ESD", TString method="", TString cumulantOrder="", Int_t mode=mLocal)
{ 
 // 1. type: "ESD", "AOD", "ESDMCkineESD", "ESDMCkineMC", for Monte Carlo and 'on the fly' use simply "", ;
 // 2. method: MCEP, LYZ1, LYZ2, LYZEP, SP, FQD, GFC or QC; 
 // 3. cumulantOrder: 2nd, 4th, 6th or 8th;             
 // 4. mode: if mode = mLocal -> analyze data on your computer using aliroot
 //          if mode = mLocalSource -> analyze data on your computer using root + source files 
 
 // Cross-check if the user's settings make sense:
 CrossCheckSettings();
 
 // Load needed libraries:
 LoadLibrariesMW(mode);

 // Name of the common output file:
 TString outputFileName = "AnalysisResults.root";
 TFile *outputFile = NULL;
 // Access the common output file:
 if(!(gSystem->AccessPathName(Form("%s%s%s",gSystem->pwd(),"/",outputFileName.Data()),kFileExists)))
 {
  outputFile = TFile::Open(outputFileName.Data(),"READ");
 } else
   {
    cout<<endl;
    cout<<"WARNING: Couldn't find the file "<<outputFileName.Data()<<" in "<<endl;
    cout<<"         directory "<<gSystem->pwd()<<" !!!!"<<endl;
    cout<<endl;
    exit(0);
   }
 
 // Access the output file of the specified method in the common output file:
 TString methodFileName = "output";
 ((methodFileName.Append(method.Data())).Append("analysis")).Append(type.Data()))); 
 TDirectoryFile *methodFile = (TDirectoryFile*)outputFile->FindObjectAny(methodFileName.Data());
 TList *methodList = NULL;
 if(methodFile)
 {
  methodFile->GetObject(Form("cobj%s",method.Data()),methodList);
  if(!methodList)
  {
   cout<<endl;
   cout<<"WARNING: Couldn't access the list "<<methodList->GetName()<<" !!!!"<<endl;
   cout<<endl;
   exit(0);  
  }   
 } else 
   {
    cout<<endl;
    cout<<"WARNING: Couldn't find the file "<<methodFileName.Data()<<".root in "<<endl;
    cout<<"         common output file "<<gSystem->pwd()<<"/"<<outputFileName.Data()<<" !!!!"<<endl;
    cout<<endl;
    exit(0);
   }
 
 // Accessing common control and results histograms from which phi, pt and eta weights will be made:
 AliFlowCommonHist *commonHist = NULL;
 AliFlowCommonHistResults *commonHistRes = NULL; 
 if(!(method=="GFC"||method=="QC"))
 {
  commonHist = dynamic_cast<AliFlowCommonHist*> methodList->FindObject(Form("AliFlowCommonHist%s",method.Data()));
  commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject(Form("AliFlowCommonHistResults%s",method.Data()));
 } else if(method=="GFC") // GFC has distinct common hist results for different cumulant orders (but control histos are the same for different orders)
   {
    commonHist = dynamic_cast<AliFlowCommonHist*> methodList->FindObject(Form("AliFlowCommonHist%s",method.Data()));
    commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject(Form("AliFlowCommonHistResults%sOrder%s",cumulantOrder.Data(),method.Data()));    
   } else // this is for sure QC - treated separately because it has distinct both common control and result histograms for different QC orders
     {
      commonHist = dynamic_cast<AliFlowCommonHist*> methodList->FindObject(Form("AliFlowCommonHist%sOrder%s",cumulantOrder.Data(),method.Data()));
      commonHistRes = dynamic_cast<AliFlowCommonHistResults*> methodList->FindObject(Form("AliFlowCommonHistResults%sOrder%s",cumulantOrder.Data(),method.Data()));
     }     
 if(!commonHist)
 {
  cout<<endl;
  cout<<"WARNING: commonHist is NULL !!!!"<<endl;
  cout<<endl;
  exit(0);
 }
 if(!commonHistRes)
 {
  cout<<endl;
  cout<<"WARNING: commonHistRes is NULL !!!!"<<endl;
  cout<<endl;
  exit(0);
 }

 // Making the output file "weights.root" in which list "weights" will be saved.
 // List "weights" will hold histograms phiWeights, ptWeights and etaWeights which
 // will hold the phi, pt and eta weights, respectively:
 TFile *weightsFile = new TFile("weights.root","RECREATE"); 
 TList *weightsList = new TList();
 gStyle->SetOptStat(0); // remove statistic box from all histograms

 // phi-weights:
 TH1F *phiWeights = (TH1F*)commonHist->GetHistPhiRP()->Clone("phi_weights"); // to be improved (transferred into TH1D eventually)
 phiWeights->SetTitle("#phi-weights: correcting for non-uniform acceptance");
 phiWeights->SetYTitle("w_{#phi}");
 phiWeights->SetXTitle("#phi"); 
 Int_t nBinsPhi = 0; // number of phi bins
 Int_t counterOfEmptyBinsPhi = 0; // number of empty phi bins
 Double_t nParticlesInBin = 0.; // number of particles in particular phi bin
 Double_t nParticlesPerBin = 0.; // average number of particles per phi bin
 Double_t nParticles = 0.; // number of particles in all phi bins 
 // calculate phi-weights:
 nBinsPhi = phiWeights->GetNbinsX();
 nParticles = phiWeights->Integral();
 if(nBinsPhi) nParticlesPerBin = nParticles/nBinsPhi; 
 for(Int_t b=1;b<=nBinsPhi;b++)
 {
  Double_t wPhi = 0.; // phi-weight for particular phi bin 
  nParticlesInBin = phiWeights->GetBinContent(b);
  if(nParticlesInBin) 
  {
   wPhi = nParticlesPerBin/nParticlesInBin;
  } else
    {
     counterOfEmptyBinsPhi++;
    }
  phiWeights->SetBinContent(b,wPhi);
 }
 if(!counterOfEmptyBinsPhi)
 {
  weightsList->Add(phiWeights);
  cout<<"Phi weights created."<<endl;
 } else
   {
    cout<<"WARNING: Couldn't create phi weights because "<<counterOfEmptyBinsPhi<<" phi bins were empty !!!!"<<endl;
   }
 
// phi-weights for eta subevents:
if (method=="SP"){
  //subevent 0
  TH1F *phiWeightsSub0 = (TH1F*)commonHist->
    GetHistPhiSub0()->Clone("phi_weights_sub0"); 
  phiWeightsSub0->SetTitle("#phi-weights for subevent 0");
  phiWeightsSub0->SetYTitle("w_{#phi}");
  phiWeightsSub0->SetXTitle("#phi"); 
  Int_t nBinsPhiSub0 = 0; // number of phi bins
  Int_t counterOfEmptyBinsPhiSub0 = 0; // number of empty phi bins
  Double_t nParticlesInBinSub0 = 0.; // number of particles in this phi bin
  Double_t nParticlesPerBinSub0 = 0.; // average number of particles/bin
  Double_t nParticlesSub0 = 0.; // number of particles in all phi bins 
  //subevent 1
  TH1F *phiWeightsSub1 = (TH1F*)commonHist->
    GetHistPhiSub1()->Clone("phi_weights_sub1"); 
  phiWeightsSub1->SetTitle("#phi-weights for subevent 0");
  phiWeightsSub1->SetYTitle("w_{#phi}");
  phiWeightsSub1->SetXTitle("#phi"); 
  Int_t nBinsPhiSub1 = 0; // number of phi bins
  Int_t counterOfEmptyBinsPhiSub1 = 0; // number of empty phi bins
  Double_t nParticlesInBinSub1 = 0.; // number of particles in this phi bin
  Double_t nParticlesPerBinSub1 = 0.; // average number of particles/bin
  Double_t nParticlesSub1 = 0.; // number of particles in all phi bins 

  // calculate phi-weights for subevent 0:
  nBinsPhiSub0 = phiWeightsSub0->GetNbinsX();
  nParticlesSub0 = phiWeightsSub0->Integral();
  if(nBinsPhiSub0) nParticlesPerBinSub0 = nParticlesSub0/nBinsPhiSub0; 
  for(Int_t b=1;b<=nBinsPhiSub0;b++) {
    Double_t wPhiSub0 = 0.; // phi-weight for particular phi bin 
    nParticlesInBinSub0 = phiWeightsSub0->GetBinContent(b);
    if(nParticlesInBinSub0) {
      wPhiSub0 = nParticlesPerBinSub0/nParticlesInBinSub0;
    } else {
      counterOfEmptyBinsPhiSub0++;
    }
    phiWeightsSub0->SetBinContent(b,wPhiSub0);
  }
  if(!counterOfEmptyBinsPhiSub0) {
    weightsList->Add(phiWeightsSub0);
    cout<<"Phi weights created."<<endl;
  } else {
    cout<<"WARNING: Couldn't create phi weights for subevent 0 because "<<counterOfEmptyBinsPhiSub0<<" phi bins were empty !!!!"<<endl;
  }

  // calculate phi-weights for subevent 1:
  nBinsPhiSub1 = phiWeightsSub1->GetNbinsX();
  nParticlesSub1 = phiWeightsSub1->Integral();
  if(nBinsPhiSub1) nParticlesPerBinSub1 = nParticlesSub1/nBinsPhiSub1; 
  for(Int_t b=1;b<=nBinsPhiSub1;b++) {
    Double_t wPhiSub1 = 0.; // phi-weight for particular phi bin 
    nParticlesInBinSub1 = phiWeightsSub1->GetBinContent(b);
    if(nParticlesInBinSub1) {
      wPhiSub1 = nParticlesPerBinSub1/nParticlesInBinSub1;
    } else {
      counterOfEmptyBinsPhiSub1++;
    }
    phiWeightsSub1->SetBinContent(b,wPhiSub1);
  }
  if(!counterOfEmptyBinsPhiSub1) {
    weightsList->Add(phiWeightsSub1);
    cout<<"Phi weights created."<<endl;
  } else {
    cout<<"WARNING: Couldn't create phi weights for subevent 1 because "<<counterOfEmptyBinsPhiSub1<<" phi bins were empty !!!!"<<endl;
  }

 }

 // pt-weights:  
 Double_t ptMin = AliFlowCommonConstants::GetMaster()->GetPtMin();
 Double_t ptMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
 Int_t nBinsPt  = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
 Double_t ptBinWidth = 0.;
 if(nBinsPt) ptBinWidth = (ptMax-ptMin)/nBinsPt;
 TH1D *ptWeights = new TH1D("pt_weights","",nBinsPt,ptMin,ptMax);
 ptWeights->SetXTitle("p_{t} [GeV]");
 ptWeights->SetYTitle("w_{p_{T}}");
 if(useLinearPtWeights) 
 {
  ptWeights->SetTitle("Linear p_{T}-weights: optimizing the flow signal");
 } else if(useQuadraticPtWeights)
   { 
    ptWeights->SetTitle("Quadratic p_{T}-weights: optimizing the flow signal");
   } else
     {
      ptWeights->SetTitle("Differential flow as p_{T}-weights: optimizing the flow signal");   
     }    
 // calculate pt weights:    
 for(Int_t b=1;b<=nBinsPt;b++)
 {
  if(useLinearPtWeights)
  {
   if(ptMin+b*ptBinWidth < ptCutoff)
   {
    ptWeights->SetBinContent(b,(ptMin+b*ptBinWidth)*(ptWeightsMax/ptCutoff)); 
   } else
     {
      ptWeights->SetBinContent(b,ptWeightsMax); 
     }  
     if(b==nBinsPt)
     {
      weightsList->Add(ptWeights);
      cout<<"Pt weights (linear) created."<<endl;
     }
  } else if(useQuadraticPtWeights)
    {
     if(ptMin+b*ptBinWidth < ptCutoff)
     {
      ptWeights->SetBinContent(b,pow(ptMin+b*ptBinWidth,2.)*(ptWeightsMax/pow(ptCutoff,2.))); 
     } else
       {
        ptWeights->SetBinContent(b,ptWeightsMax); 
       } 
       if(b==nBinsPt)
       {
        weightsList->Add(ptWeights);
        cout<<"Pt weights (quadratic) created."<<endl;
       }
    } else // differential flow result is used as a pt-weight: 
      {
       ptWeights->SetBinContent(b,commonHistRes->GetHistDiffFlowPtPOI()->GetBinContent(b));
       if(b==nBinsPt)
       {
        weightsList->Add(ptWeights);
        cout<<"Pt weights (from differential flow) created."<<endl;
       }
      } 
 } // end of for(Int_t b=1;b<=nBinsPt;b++)

 // eta-weights:
 TH1D *etaWeights = commonHistRes->GetHistDiffFlowEtaPOI()->Clone("eta_weights");
 etaWeights->SetXTitle("#eta [GeV]");
 etaWeights->SetYTitle("w_{#eta}");
 etaWeights->SetTitle("Differential flow as #eta-weights: optimizing the flow signal"); 
 if(etaWeights) 
 {
  weightsList->Add(etaWeights);
  cout<<"Eta weights (from differential flow) created."<<endl;
 } 
 
 // Save list holding histogram with weights:
 weightsFile->WriteObject(weightsList,"weights","SingleKey");
 
 cout<<"New file \"weights.root\" created to hold those phi, pt and eta weights."<<endl;

 delete weightsList;
 delete weightsFile; 
 
} // end of void makeWeights(TString type="", TString method="QC", TString cumulantOrder="4th", Int_t mode=mLocal)

void CrossCheckSettings() 
{
 // Check in this method if the settings make sense
 
 if(useLinearPtWeights && useQuadraticPtWeights)
 {
  cout<<endl;
  cout<<"WARNING: You cannot set useLinearPtWeights and useQuadraticPtWeights to kTRUE"<<endl;
  cout<<"          at the same time. Please make up your mind."<<endl;
  cout<<endl;
  exit(0);
 }
 if(ptCutoff<0.)
 { 
  cout<<endl;
  cout<<"WARNING: It doesn't make much sense to have ptCutoff < 0."<<endl;
  cout<<endl;
  exit(0);
 } 
 if(ptWeightsMax<0.)
 { 
  cout<<endl;
  cout<<"WARNING: It doesn't make much sense to have ptWeightsMax < 0."<<endl;
  cout<<endl;
  exit(0);
 } 
 
} // end of void CrossCheckSettings()

void LoadLibrariesMW(const libModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  //gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");
  gSystem->Load("libPhysics");
  
  //----------------------------------------------------------
  // >>>>>>>>>>> Local mode <<<<<<<<<<<<<< 
  //----------------------------------------------------------
  if (mode==mLocal) {
    //--------------------------------------------------------
    // If you want to use already compiled libraries 
    // in the aliroot distribution
    //--------------------------------------------------------
    
    //==================================================================================  
    //load needed libraries:
    gSystem->AddIncludePath("-I$ROOTSYS/include");
    //gSystem->Load("libTree");
    
    // for AliRoot
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libPWG2flowCommon");
    //cerr<<"libPWG2flowCommon loaded ..."<<endl;
    
  }
  
  else if (mode==mLocalSource) {
    
    // In root inline compile
   
   // Constants  
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZConstants.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCumuConstants.cxx+");
    
    // Flow event
    gROOT->LoadMacro("AliFlowCommon/AliFlowVector.cxx+"); 
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimple.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowEvent.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowEventSimple.cxx+");
    
    // Cuts
    gROOT->LoadMacro("AliFlowCommon/AliFlowTrackSimpleCuts.cxx+");    
    
    // Output histosgrams
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHist.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowCommonHistResults.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist1.cxx+");
    gROOT->LoadMacro("AliFlowCommon/AliFlowLYZHist2.cxx+");
    
    cout << "finished loading macros!" << endl;  
    
  }  
  
} // end of void LoadLibrariesMW(const libModes mode) 





