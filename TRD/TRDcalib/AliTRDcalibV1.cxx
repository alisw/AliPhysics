#define AliTRDcalibV1_cxx
// The class definition in Calib.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("AliTRDcalibV1.cxx")
// Root > T->Process("AliTRDcalibV1.cxx","some options")
// Root > T->Process("AliTRDcalibV1.cxx+")
//

#include "AliTRDcalibV1.h"
#include <TTree.h>
#include <TObject.h>
#include <TH2.h>
#include <TStyle.h>
#include <TH2I.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <AliESD.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h> 
#include <AliESDtrack.h>
#include <AliESDfriendTrack.h>

#include <AliTRDgeometry.h>
#include <AliTRDtrack.h>
#include <AliCDBManager.h>
#include <AliTRDCalibraFillHisto.h>
#include <AliTRDCalibraVdriftLinearFit.h>


AliTRDcalibV1::AliTRDcalibV1(TTree *) : 
   TSelector(),
   fESD(0),
   fev(0),
   fevf(0),
   fo(0),
   ft(0),
   fesdTrack(0),
   ffriendTrack(0),
   fFileNo(0)     
 {
   //G__SetCatchException(0);     
 }  
//_____________________________________________________________________
void AliTRDcalibV1::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   //TString option = GetOption();
  
}
//______________________________________________________________________
void AliTRDcalibV1::SlaveBegin(TTree * tree)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

  //printf("Slave Begin\n");

   //TString option = GetOption();
  if(tree) Init(tree);

  fo = 0x0;
  ft = 0x0;
  fesdTrack = 0x0;
  ffriendTrack = 0x0;

  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //cdbManager->SetSpecificStorage("TRD/Calib/FEE","local:///u/bailhach/aliroot/database30head/");
  cdbManager->SetRun(0);

  // instance calibration
  fcalib = AliTRDCalibraFillHisto::Instance();
  fcalib->SetNz(0,0);
  fcalib->SetNrphi(0,0);
  fcalib->SetNz(1,0);
  fcalib->SetNrphi(1,0);
  fcalib->SetHisto2d();
  fcalib->SetVector2d();
  fcalib->SetLinearFitterOn();
  fcalib->SetLinearFitterDebugOn();
  fcalib->SetCH2dOn();
  fcalib->SetPH2dOn();
  fcalib->SetPRF2dOn();
  fcalib->Init2Dhistos();
  //fcalib->SetDebugLevel(1);
  fcalib->SetNumberClusters(14);

}
//__________________________________________________________________________
void   AliTRDcalibV1::CleanESD(){
  //
  //printf("CleanESD\n");
  if (fev!=0){
    delete fev;
    fev = 0;
  }
  if (fevf!=0){
    delete fevf;
    fevf =0;
  }
}
//_________________________________________________________________________________
Bool_t AliTRDcalibV1::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Calib::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
  //printf("process\n");
 
  if (!fChain) return kFALSE;  
  //printf("process1\n");
  Int_t nBytes;
  Int_t nTRDcls = 0;
  nBytes = fChain->GetTree()->GetEntry(entry);
  //printf("There are %d bytes for these event\n",nBytes);
  if (!fev || (nBytes == 0)) { 
    return kFALSE;
  }
  if(fev->GetAliESDOld()) fev->CopyFromOldESD();
  //printf("process2\n");
  Int_t ntr = fev->GetNumberOfTracks();
  //printf("Tracks new %d\n",ntr);
  
  if (!fevf || (fevf->GetNumberOfTracks()!=ntr)) {
    return kFALSE;
  }
  
  if(ntr>0){
    
    fev->SetESDfriend(fevf);
    
    //printf("Number of friends tracks %d\n",fevf->GetNumberOfTracks());
    
    
    for(int itrk=0; itrk<fev->GetNumberOfTracks(); itrk++){
     
      fesdTrack = fev->GetTrack(itrk);
      if(!(nTRDcls = fesdTrack->GetTRDncls())) continue;
      if(!(fesdTrack->GetFriendTrack())) continue;
      //ffriendTrack = new AliESDfriendTrack(*(fesdTrack->GetFriendTrack()));
      ffriendTrack = fevf->GetTrack(itrk);        
         
      Int_t icalib=0;
      while((fo = (TObject *)(ffriendTrack->GetCalibObject(icalib++)))){
	//printf("Name of calibObject %s\n",fo->IsA()->GetName());
	if(strcmp(fo->IsA()->GetName(), "AliTRDtrackV1") != 0) continue;
	//printf("\tfound %s @ 0x%x; calib object %d\n", fo->IsA()->GetName(), fo, icalib-1);
	ft = (AliTRDtrackV1 *)fo;
	
	fcalib->UpdateHistogramsV1(ft);
      }
    }
  }
  
  //CleanESD();
  return kTRUE;
}
//______________________________________________________________________________________________
void AliTRDcalibV1::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

  if(!fOutput)
    {
      printf("ERROR: Output list not initialized\n");
      return;
    }
  
  fCH2d = new TH2I(*(fcalib->GetCH2d()));
  fPH2d = new TProfile2D(*(fcalib->GetPH2d()));
  fPRF2d = new TProfile2D(*(fcalib->GetPRF2d()));

  AliTRDCalibraVdriftLinearFit *ju = fcalib->GetVdriftLinearFit();
  for(Int_t det = 0; det < 540; det++){
    fVdriftLinear[det] = new TH2F(*(ju->GetLinearFitterHisto(det,kTRUE)));
  }
  


  fOutput->Add(fCH2d);
  fOutput->Add(fPH2d);
  fOutput->Add(fPRF2d);
  for(Int_t det = 0; det < 540; det++){
    fOutput->Add(fVdriftLinear[det]);
  }

  fcalib->Destroy();
 
}
//____________________________________________________________________________________
void AliTRDcalibV1::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

 
  printf("InTerminate()\n");
  if (!fOutput) return;
 
  //fOutput->Print();
  
  fCH2d = dynamic_cast<TH2I*>(fOutput->FindObject("CH2d"));
  
  if (!fCH2d)
    {
      printf("Error: %s not returned\n","fCH2d");
      return;
    }

  fPH2d = dynamic_cast<TProfile2D*>(fOutput->FindObject("PH2d"));
  
  if (!fPH2d)
    {
      printf("Error: %s not returned\n","fPH2d");
      return;
    }

  fPRF2d = dynamic_cast<TProfile2D*>(fOutput->FindObject("PRF2d"));
  
  if (!fPRF2d)
    {
      printf("Error: %s not returned\n","fPRF2d");
      return;
    }


  const char * Name = 0x0;
  Name = "LFDV%dversion0";
  TObjArray array(540);

  for(Int_t det = 0; det < 540; det++){

    TString Namehisto (Form(Name, det));
    fVdriftLinear[det] = dynamic_cast<TH2F*>(fOutput->FindObject((const char *)Namehisto));
    if (!fVdriftLinear[det])
      {
	printf("Error: %s not returned\n",(const char *)Namehisto);
	return;
      }
    array.AddAt(fVdriftLinear[det],det);
  }

  AliTRDCalibraVdriftLinearFit vdriftlinearfit = AliTRDCalibraVdriftLinearFit(array);
  vdriftlinearfit.FillPEArray();

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.10);

 
  printf("There are %d files analysed\n",fFileNo);
  
  
  TCanvas *u = new TCanvas("u","",50,50,600,800);
  u->Divide(3,1);
  u->cd(1);
  fCH2d->DrawCopy("lego");
  u->cd(2);
  fPH2d->DrawCopy("lego");
  u->cd(3);
  fPRF2d->DrawCopy("lego");

  TObjArray *arrayE = vdriftlinearfit.GetEArray();

  Double_t totalsum = 0.0;
  for(Int_t k = 0; k < 540; k++){
    TVectorD *h = (TVectorD *) arrayE->UncheckedAt(k);
    if(h){
      totalsum += (*h)[2];
    }
  }

  TFile file("Output.root","recreate");
  fOutput->Write();  

}
//___________________________________________________________________________________________
void AliTRDcalibV1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  //if (counter>1) return;
  tree->SetBranchStatus("*",1);
  //
  // New AliESDevent format
  //
  if (!fChain->GetBranch("ESD")){
    //
    //
    //
    if (fev) delete fev;
    fev = new AliESDEvent();
    fev->ReadFromTree(tree); // Attach the branch with ESD friends
    fevf = (AliESDfriend*)fev->FindListObject("AliESDfriend");
    tree->SetBranchAddress("ESDfriend.",&fevf); 
    return;
  }
  
  fChain->SetBranchAddress("ESD",&fESD);
  Info("Init","Enter");
  Bool_t isOK=kFALSE;
  if (fChain->GetBranch("ESDfriend")) {
    fChain->SetBranchAddress("ESDfriend",&fevf);
    Info("Init","V0-ESDfriend.");
    isOK=kTRUE;
  }
  if (fChain->GetBranch("ESDfriend.")){
    Info("Init","V1-ESDfriend.");
    fChain->SetBranchAddress("ESDfriend.",&fevf);
    isOK=kTRUE;
  }
  return;

}
//___________________________________________________________________________________________
Bool_t AliTRDcalibV1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

  ++fFileNo;
  printf ("Processing file no %d\n",fFileNo);
  
  return kTRUE;
}
//__________________________________________________________________________________
Int_t AliTRDcalibV1::ReadEvent(Long64_t entry){
  //
  //
  //


  if (!fChain) return -1;  
  if (!fChain->GetTree()) return -1; 
  try {
    fChain->GetTree()->GetEntry(entry);
  } catch (std::bad_alloc) {
    printf("Pica vyjebana pojebany skurveny kokot piciak\n");
    return -1;
  }
  if (!fESD && !fev) { 
    return -2;
  }
  Int_t ntracks = (fESD) ? fESD->GetNumberOfTracks() : fev->GetNumberOfTracks();   

  if (fev){
    fev->SetESDfriend(fevf);   
  }

  
  if (!fevf || fevf->GetNumberOfTracks() != ntracks) {
    try {
      delete fESD;
    }
    catch (std::bad_alloc) {
      printf("Pica vyjebana pojebany skurveny kokot piciak\n");
      fESD =0;
      return -1;
    }
    return -3;
  }
  if (fESD) fESD->SetESDfriend(fevf);
  return 0;
}
