//==========================================================================================
// Before using the macro makeWeights.C you should already have available the output .root 
// files from various methods from the previous run over data (without any weights). When 
// calling this macro you must specify the analysis type and the method from which output 
// file you would like to make the weights for the next run (for the cumulants, GFC and QC,
// you must also specify the order): 
// 
// 1. type of analysis can be: ESD, AOD, MC, ESDMC0 or ESDMC1;
//
// 2. method can be: MCEP, LYZ1, LYZ2, LYZEP, FQD, GFC or QC; 
//
// 3. cumulant order can be: 2nd, 4th, 6th or 8th.                                                   
//==========================================================================================


enum libModes {mLocal,mLocalSource};
//mLocal: Analyze data on your computer using aliroot
//mLocalSource: Analyze data on your computer using root + source files


//void makeWeights(TString type="", TString method="GFC", TString cumulantOrder="4th", Int_t mode=mLocalSource)
void makeWeights(TString type="ESD", TString method="GFC", TString cumulantOrder="4th", Int_t mode=mLocal)
{
 // load needed libraries:
  LoadWeightLibraries(mode);


 // open the output file from the first run of the specified method:
 TString inputFileName = "output";
 TFile* file = NULL;
 file = TFile::Open(((((inputFileName.Append(method.Data())).Append("analysis")).Append(type.Data())).Append(".root")).Data(), "READ"); 
 
 // using pt weights linear or quadratic in pt:
 Bool_t useLinearPt    = kTRUE;
 Bool_t useQuadraticPt = kFALSE;
 if(useLinearPt && useQuadraticPt)
 {
  cout<<" WARNING: you can use pt weights linear or quadratic in pt, but not at the same time. "<<endl;
  exit(0);
 }
 Bool_t useFunctionOfPt = useLinearPt||useQuadraticPt; 
 
 // accessing the results:
 TString cobj = "cobj";
 TString afc  = "AliFlowCommonHist";
 TString afcr = "AliFlowCommonHistResults";
 
 TList *pList = NULL;
 AliFlowCommonHist *commonHist = NULL;
 AliFlowCommonHistResults *commonHistRes = NULL; 
 
 if(file) 
 {
  file->GetObject((cobj.Append(method.Data()).Data()),pList); 
  if(pList) 
  {
   if(!(method=="GFC"||method=="QC"))
   {
    commonHist    = dynamic_cast<AliFlowCommonHist*> (pList->FindObject((afc.Append(method.Data())).Data()));
    commonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pList->FindObject((afcr.Append(method.Data())).Data()));
   }else if(method=="GFC")
    {
     commonHist    = dynamic_cast<AliFlowCommonHist*> (pList->FindObject((afc.Append(method.Data())).Data()));
     commonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pList->FindObject((((afcr.Append(cumulantOrder.Data()).Append("Order")).Append(method.Data())).Data())));     
    }else
     {
      commonHist    = dynamic_cast<AliFlowCommonHist*> (pList->FindObject(((afc.Append(cumulantOrder.Data())).Append("Order")).Append(method.Data())));
      commonHistRes = dynamic_cast<AliFlowCommonHistResults*> (pList->FindObject((((afcr.Append(cumulantOrder.Data()).Append("Order")).Append(method.Data())).Data())));
     }     
  }//end of if(pList)  
 }//end of if(file) 

 //making the output file and creating the TList to hold the histograms with weights:
 TFile* outputFile = new TFile("weights.root","RECREATE"); 
 TList* listWeights = new TList();
 Int_t nBinsPhi = 0;
 Double_t nParticlesInBin = 0;  // number of particles in particular phi bin
 Double_t nParticlesPerBin = 0; // average number of particles per phi bin
 Double_t nParticles = 0;       // number of particles in all phi bins 
 Double_t wPhi = 0.;  
 
 // common control histos:
 if(commonHist)
 {
  // azimuthal acceptance:
  (commonHist->GetHistPhiRP())->SetName("phi_weights");
  (commonHist->GetHistPhiRP())->SetTitle("phi_weights: correction for non-uniform acceptance");
  (commonHist->GetHistPhiRP())->SetYTitle("weights");
  (commonHist->GetHistPhiRP())->SetXTitle("#phi"); 
  nBinsPhi = (commonHist->GetHistPhiRP())->GetNbinsX();
  nParticles = (commonHist->GetHistPhiRP())->Integral();
  if(nBinsPhi) nParticlesPerBin = nParticles/nBinsPhi; 
  // loop over phi bins:
  for(Int_t b=1;b<nBinsPhi+1;b++)
  {
   nParticlesInBin = (commonHist->GetHistPhiRP())->GetBinContent(b);
   // calculate the phi weight wPhi for each bin:
   if(nParticlesInBin) wPhi = nParticlesPerBin/nParticlesInBin;
   (commonHist->GetHistPhiRP())->SetBinContent(b,wPhi);
  }
  listWeights->Add(commonHist->GetHistPhiRP());
 }else{cout<<" WARNING: the common control histos from the 1st run were not accessed."<<endl;} 
 
 // common results histos:
 if(commonHistRes)
 {
  // diff. flow (pt):
  (commonHistRes->GetHistDiffFlowPtPOI())->SetName("pt_weights");
  if(!useFunctionOfPt) listWeights->Add(commonHistRes->GetHistDiffFlowPtPOI());
  // diff. flow (eta):
  (commonHistRes->GetHistDiffFlowEtaPOI())->SetName("eta_weights");
  listWeights->Add(commonHistRes->GetHistDiffFlowEtaPOI());
 }else{cout<<" WARNING: the common results histos from the 1st run were not accessed."<<endl;}  
 
 // pt weights linear and quadratic in pt:
 if(useFunctionOfPt)
 {
  Double_t ptMin = AliFlowCommonConstants::GetPtMin();
  Double_t ptMax = AliFlowCommonConstants::GetPtMax();
  Int_t nBinsPt  = AliFlowCommonConstants::GetNbinsPt();
  Double_t ptCutOff = 2.0; // for pt > ptCutOff use constant weights 
  if(nBinsPt==0) 
  { 
   cout<<" WARNING: number of pt bins is 0. "<<endl;
   exit(0);
  } 
  Double_t ptBinWidth = 1.*(ptMax-ptMin)/nBinsPt;
 
  TH1D *ptFunction = new TH1D("ptFunction", "ptFunction", nBinsPt, ptMin, ptMax);
  ptFunction->SetName("pt_weights");
  ptFunction->SetXTitle("p_{t} [GeV]");
  ptFunction->SetYTitle("weights");
  Int_t power = 0;
  if(useLinearPt) 
  {
   power = 1;
   ptFunction->SetTitle("p_{t} weights: #propto p_{t}");
  } 
  if(useQuadraticPt)
  { 
   power = 2;
   ptFunction->SetTitle("p_{t} weights: #propto p_{t}^{2}");
  }   
  for(Int_t b=1;b<nBinsPt+1;b++)
  {
   if(ptMin+b*ptBinWidth < ptCutOff)
   {
    ptFunction->SetBinContent(b, pow(ptMin+b*ptBinWidth, power)); 
   }else
    {
     ptFunction->SetBinContent(b, pow(ptCutOff, power)); 
    }  
  } 
  listWeights->Add(ptFunction); 
 }   
        
 outputFile->WriteObject(listWeights,"weights","SingleKey");

 delete listWeights;
 delete outputFile; 
}  


void LoadWeightLibraries(const libModes mode) {
  
  //--------------------------------------
  // Load the needed libraries most of them already loaded by aliroot
  //--------------------------------------
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libXMLIO.so");
  gSystem->Load("libPhysics.so");
  
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
    gSystem->Load("libTree.so");
    
    // for AliRoot
    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libPWG2flowCommon.so");
    cerr<<"libPWG2flowCommon.so loaded ..."<<endl;
    
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
  
}



