//=====================================================================================
// Before using the macro makeWeightsForTheSecondRun.C you should already have         
// available the output .root files from various methods from the first run over data 
// without any weights. When calling this macro you must specify the analysis type 
// and the method from which output file you would like to make the weights for the 
// second run (for the cumulants, GFC and QC, you must also specify the order): 
// 
// 1. type of analysis can be: ESD, AOD, MC, ESDMC0 or ESDMC1;
//
// 2. method can be: MCEP, LYZ1, LYZ2, LYZEP, FQD, GFC or QC; 
//
// 3. cumulant order can be: 2nd, 4th, 6th or 8th.                                                   
//=====================================================================================

void makeWeightsForTheSecondRun(TString type="ESD", TString method="GFC", TString cumulantOrder="4th")
{
 //load needed libraries:
 gSystem->AddIncludePath("-I$ROOTSYS/include");
 gSystem->Load("libTree.so");

 //for AliRoot
 gSystem->AddIncludePath("-I$ALICE_ROOT/include");
 gSystem->Load("libANALYSIS.so");
 gSystem->Load("libPWG2flow.so");
 cerr<<"libPWG2flow.so loaded ..."<<endl;

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

 //making the output file and storing the histograms needed for the weights:
 TFile* outputFile = new TFile("weightsForTheSecondRun.root","RECREATE"); 
 //common control histos:
 if(commonHist)
 {
  //azimuthal acceptance:
  (commonHist->GetHistPhiInt())->SetName("phi_weights");
  //normalizing:
  Double_t norm=(commonHist->GetHistPhiInt())->Integral();
  if(norm)
  {
   (commonHist->GetHistPhiInt())->Scale(1./norm);
  } 
  //writing the normalized histogram in output file: 
  (commonHist->GetHistPhiInt())->Write();
 }else
  {
   cout<<" WARNING: the common control histos from the 1st run were not accessed."<<endl;
  } 
 //common results histos:
 if(commonHistRes)
 {
  //diff. flow in pt from the first run:
  (commonHistRes->GetHistDiffFlowPtPOI())->SetName("pt_weights");
  (commonHistRes->GetHistDiffFlowPtPOI())->Write();
  //diff. flow in eta from the first run:
  (commonHistRes->GetHistDiffFlowEtaPOI())->SetName("eta_weights");
  (commonHistRes->GetHistDiffFlowEtaPOI())->Write();
 }else
  {
   cout<<" WARNING: the common results histos from the 1st run were not accessed."<<endl;
  }  
 
 delete outputFile;
}  
