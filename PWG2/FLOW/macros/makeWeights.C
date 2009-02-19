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

void makeWeights(TString type="ESD", TString method="GFC", TString cumulantOrder="4th")
{
 //load needed libraries:
 gSystem->AddIncludePath("-I$ROOTSYS/include");
 gSystem->Load("libTree.so");

 //for AliRoot
 gSystem->AddIncludePath("-I$ALICE_ROOT/include");
 gSystem->Load("libANALYSIS.so");
 gSystem->Load("libPWG2flowCommon.so");
 cerr<<"libPWG2flowCommon.so loaded ..."<<endl;

 //open the output file from the first run of the specified method:
 TString inputFileName = "output";
 TFile* file = NULL;
 file = TFile::Open(((((inputFileName.Append(method.Data())).Append("analysis")).Append(type.Data())).Append(".root")).Data(), "READ"); 
   
 //accessing the results:
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
 //common control histos:
 if(commonHist)
 {
  //azimuthal acceptance:
  (commonHist->GetHistPhiInt())->SetName("phi_weights");
  (commonHist->GetHistPhiInt())->SetTitle("Phi_weights: correction for non-uniform acceptance");
  (commonHist->GetHistPhiInt())->SetYTitle("weights");
  (commonHist->GetHistPhiInt())->SetXTitle("#phi"); 
  nBinsPhi = (commonHist->GetHistPhiInt())->GetNbinsX();
  nParticles = (commonHist->GetHistPhiInt())->Integral();
  if(nBinsPhi) nParticlesPerBin = nParticles/nBinsPhi; 
  // loop over phi bins:
  for(Int_t b=1;b<nBinsPhi+1;b++)
  {
   nParticlesInBin = (commonHist->GetHistPhiInt())->GetBinContent(b);
   // calculate the phi weight wPhi for each bin:
   if(nParticlesInBin) wPhi = nParticlesPerBin/nParticlesInBin;
   (commonHist->GetHistPhiInt())->SetBinContent(b,wPhi);
  }
  listWeights->Add(commonHist->GetHistPhiInt());
 }else{cout<<" WARNING: the common control histos from the 1st run were not accessed."<<endl;} 
 //common results histos:
 if(commonHistRes)
 {
  //diff. flow (pt):
  (commonHistRes->GetHistDiffFlowPtPOI())->SetName("pt_weights");
  listWeights->Add(commonHistRes->GetHistDiffFlowPtPOI());
  //diff. flow (eta):
  (commonHistRes->GetHistDiffFlowEtaPOI())->SetName("eta_weights");
  listWeights->Add(commonHistRes->GetHistDiffFlowEtaPOI());
 }else{cout<<" WARNING: the common results histos from the 1st run were not accessed."<<endl;}  
 
 outputFile->WriteObject(listWeights,"weights","SingleKey");

 delete listWeights;
 delete outputFile; 
}  



