
//Root includes
#include "TNamed.h"
//AliRoot includes
#include "AliLoader.h" 
#include "AliMUONConstants.h"
#include "AliMUONData.h"
#include "AliMUONDigit.h"
#include "AliMUONHit.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONGlobalTrigger.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"

ClassImp(AliMUONData)
 
//_____________________________________________________________________________
AliMUONData::AliMUONData():TNamed()
{
  fLoader        = 0x0;
  fHits          = 0x0;    // One event in treeH per primary track
  fDigits        = 0x0;  // One event in treeH per detection plane
  fNdigits       = 0x0;
  fRawClusters   = 0x0; //One event in TreeR/RawclusterBranch per tracking detection plane
  fGlobalTrigger = 0x0; //! List of Global Trigger 1st event in TreeR/GlobalTriggerBranch
  fLocalTrigger  = 0x0;  //! List of Local Trigger, 1st event in TreeR/LocalTriggerBranch
  fRecTracks     = 0x0;       
//default constructor
}
//_____________________________________________________________________________
AliMUONData::AliMUONData(AliLoader * loader, const char* name, const char* title):
  TNamed(name,title)
{
  fLoader        = loader;
  fHits          = new TClonesArray("AliMUONHit",1000);
  fNhits         = 0;
  fDigits        = new TObjArray(AliMUONConstants::NCh());
  fNdigits       = new Int_t[AliMUONConstants::NCh()];
  for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
    fDigits->AddAt(new TClonesArray("AliMUONDigit",10000),iDetectionPlane); 
    fNdigits[iDetectionPlane]=0;
  }
  fRawClusters   = new TObjArray(AliMUONConstants::NTrackingCh());
  fNrawclusters  = new Int_t[AliMUONConstants::NTrackingCh()];
  for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NTrackingCh();iDetectionPlane++) {
    fRawClusters->AddAt(new TClonesArray("AliMUONRawCluster",10000),iDetectionPlane); 
    fNrawclusters[iDetectionPlane]=0;
  }
  fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1);    
  fNglobaltrigger =0;
  fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);   
  fNlocaltrigger = 0;
  fRecTracks     = new TClonesArray("AliMUONTrack", 10);
  fNrectracks    = 0; // really needed or GetEntriesFast sufficient ????
  //default constructor
}
//_____________________________________________________________________________
AliMUONData::AliMUONData(const AliMUONData& rMUONData):TNamed(rMUONData)
{
  // Dummy copy constructor
  ;
}
//_____________________________________________________________________________
AliMUONData::~AliMUONData()
{
  if (fHits) {
    fHits->Delete();
    delete fHits;
  }
  if (fDigits) {
    fDigits->Delete();
    delete fDigits;
  }
  if (fRawClusters) {
    fRawClusters->Delete();
    delete fRawClusters;
  }
  if (fGlobalTrigger){
    fGlobalTrigger->Delete();
    delete fGlobalTrigger;
  }  
  if (fLocalTrigger){
    fLocalTrigger->Delete();
    delete fLocalTrigger;
  }
  if (fRecTracks){
    fRecTracks->Delete();
    delete fRecTracks;
  }
  //detructor 
}
//_____________________________________________________________________________
void AliMUONData::AddDigit(Int_t id, Int_t *tracks, Int_t *charges, Int_t *digits)
{
  //
  // Add a MUON digit to the list of Digits of the detection plane id
  //
  TClonesArray &ldigits = * Digits(id) ; 
  new(ldigits[fNdigits[id]++]) AliMUONDigit(tracks,charges,digits);
}
//_____________________________________________________________________________
void AliMUONData::AddGlobalTrigger(Int_t *singlePlus, Int_t *singleMinus,
				   Int_t *singleUndef,
				   Int_t *pairUnlike, Int_t *pairLike)
{
  // add a MUON Global Trigger to the list (only one GlobalTrigger per event !)
  TClonesArray &globalTrigger = *fGlobalTrigger;
  new(globalTrigger[fNglobaltrigger++]) 
    AliMUONGlobalTrigger(singlePlus, singleMinus,  singleUndef, pairUnlike, pairLike);
}
//_____________________________________________________________________________
void AliMUONData::AddHit(Int_t fIshunt, Int_t track, Int_t iChamber, 
			 Int_t idpart, Float_t X, Float_t Y, Float_t Z, 
			 Float_t tof, Float_t momentum, Float_t theta, 
			 Float_t phi, Float_t length, Float_t destep)
{
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliMUONHit(fIshunt, track, iChamber, 
				  idpart, X, Y, Z, 
				  tof, momentum, theta, 
				  phi, length, destep);
}
//____________________________________________________________________________
void AliMUONData::AddLocalTrigger(Int_t *localtr)
{
  // add a MUON Local Trigger to the list
  TClonesArray &localTrigger = *fLocalTrigger;
  new(localTrigger[fNlocaltrigger++]) AliMUONLocalTrigger(localtr);
}
//_____________________________________________________________________________
void AliMUONData::AddRawCluster(Int_t id, const AliMUONRawCluster& c)
{
  //
  // Add a MUON rawcluster to the list in the detection plane id
  //
  TClonesArray &lrawcl = *((TClonesArray*) fRawClusters->At(id));
  new(lrawcl[fNrawclusters[id]++]) AliMUONRawCluster(c);
}
//_____________________________________________________________________________
void AliMUONData::AddRecTrack(const AliMUONTrack& track)
{
  //
  // Add a MUON rectrack
  //
  TClonesArray &lrectracks = *fRecTracks;
  new(lrectracks[fNrectracks++]) AliMUONTrack(track);
}
//____________________________________________________________________________
Bool_t   AliMUONData::IsRawClusterBranchesInTree()
{
  if (TreeR()==0x0) {
    Error("TreeR","No treeR in memory");
    return kFALSE;
  }
  else {
     char branchname[30];
     sprintf(branchname,"%sRawClusters1",GetName());
     TBranch * branch = 0x0;
     branch = TreeR()->GetBranch(branchname);
     if (branch)  return kTRUE;
     else return kFALSE;    
  }
}
//____________________________________________________________________________
Bool_t   AliMUONData::IsTriggerBranchesInTree()
{
 if (TreeR()==0x0) {
    Error("TreeR","No treeR in memory");
    return kFALSE;
  }
  else {
     char branchname[30];
     sprintf(branchname,"%sLocalTrigger",GetName());
     TBranch * branch = 0x0;
     branch = TreeR()->GetBranch(branchname);
     if (branch)  return kTRUE;
     else return kFALSE;    
  }
}
//____________________________________________________________________________
void AliMUONData::Fill(Option_t* option)
{
  // Method to fill the trees
  const char *cH   = strstr(option,"H");
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  //const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP
  
  char branchname[30];
  TBranch * branch = 0x0;

  //
  // Filling TreeH
  if ( TreeH() && cH ) {
    TreeH()->Fill();
  }  
  //
  // Filling TreeD
  if ( TreeD() && cD) {
    TreeD()->Fill();
  }

  //
  // filling rawclusters
  if ( TreeR()  && cRC ) {
    if ( IsTriggerBranchesInTree() ) {
      // Branch per branch filling
      for (int i=0; i<AliMUONConstants::NTrackingCh(); i++) {
	sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
	branch = TreeR()->GetBranch(branchname);
	branch->Fill();
      }
    }
    else  TreeR()->Fill();
  }
  
 //
  // filling trigger 
  if ( TreeR()  && cGLT) {
    if (IsRawClusterBranchesInTree()) {
      // Branch per branch filling
      sprintf(branchname,"%sLocalTrigger",GetName());
      branch = TreeR()->GetBranch(branchname); 
      branch->Fill();
      sprintf(branchname,"%sGlobalTrigger",GetName());
      branch = TreeR()->GetBranch(branchname);
      branch->Fill();
    }
    else  TreeR()->Fill();
  }
  //
  // filling tracks
  if ( TreeT() && cRT ) {
    sprintf(branchname,"%sTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    branch->Fill();
  }
}
//_____________________________________________________________________________
void AliMUONData::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the MUON.
  //
  const Int_t kBufferSize = 4000;
  char branchname[30];
  
  const char *cH   = strstr(option,"H");
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP
  
  TBranch * branch = 0x0;
  
  // Creating Branches for Hits
  if (TreeH() && cH) {
    if (fHits == 0x0)  fHits = new TClonesArray("AliMUONHit",1000);
    fNhits = 0;
    sprintf(branchname,"%sHits",GetName());  
    branch = TreeH()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return ;
    }
    branch = TreeH()->Branch(branchname,&fHits,kBufferSize);
    Info("MakeBranch","Making Branch %s for hits \n",branchname);
  }  
  
  //Creating Branches for Digits
  if (TreeD() && cD ) {
    // one branch for digits per chamber
    if (fDigits  == 0x0) {
      fDigits  = new TObjArray(AliMUONConstants::NCh());
      for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
	fDigits->AddAt(new TClonesArray("AliMUONDigit",10000),iDetectionPlane); 
      }
    }
    if (fNdigits == 0x0) {
      fNdigits = new Int_t[AliMUONConstants::NCh()];
      for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
	fNdigits[iDetectionPlane]=0;
      }
    }
    for (Int_t iDetectionPlane=0; iDetectionPlane<AliMUONConstants::NCh() ;iDetectionPlane++) {
      sprintf(branchname,"%sDigits%d",GetName(),iDetectionPlane+1);
      branch = 0x0;
      branch = TreeD()->GetBranch(branchname);
      if (branch) {  
	Info("MakeBranch","Branch %s is already in tree.",GetName());
	return;
      }
      TClonesArray * digits = Digits(iDetectionPlane); 
      branch = TreeD()->Branch(branchname, &digits, kBufferSize);
      Info("MakeBranch","Making Branch %s for digits in detection plane %d\n",branchname,iDetectionPlane+1);
      }
  }
  
  if (TreeR() && cRC ) {
    //  one branch for raw clusters per tracking detection plane
    //        
    Int_t i;
    if (fRawClusters == 0x0) {
      fRawClusters = new TObjArray(AliMUONConstants::NTrackingCh());
      for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
	fRawClusters->AddAt(new TClonesArray("AliMUONRawCluster",10000),i); 
      }
    }

    if (fNrawclusters == 0x0) {
      fNrawclusters= new Int_t[AliMUONConstants::NTrackingCh()];
      for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
	fNrawclusters[i]=0;
      }
    }
    
    for (i=0; i<AliMUONConstants::NTrackingCh() ;i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);	
      branch = 0x0;
      branch = TreeR()->GetBranch(branchname);
      if (branch) {  
	Info("MakeBranch","Branch %s is already in tree.",GetName());
	return;
      }
      branch = TreeR()->Branch(branchname, &((*fRawClusters)[i]),kBufferSize);
      Info("MakeBranch","Making Branch %s for rawcluster in detection plane %d\n",branchname,i+1);
    }
  }

  if (TreeR() && cGLT ) {
    //
    // one branch for global trigger
    //
    sprintf(branchname,"%sGlobalTrigger",GetName());
    branch = 0x0;
    
    if (fGlobalTrigger == 0x0) {
      fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1); 
      fNglobaltrigger = 0;
    }
    branch = TreeR()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return ;
    }
    branch = TreeR()->Branch(branchname, &fGlobalTrigger, kBufferSize);
    Info("MakeBranch", "Making Branch %s for Global Trigger\n",branchname);
    
    //
    // one branch for local trigger
    //  
    sprintf(branchname,"%sLocalTrigger",GetName());
    branch = 0x0;
    
    if (fLocalTrigger == 0x0) {
      fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);
      fNlocaltrigger = 0;
    }
    branch = TreeR()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return;
    }
    branch = TreeR()->Branch(branchname, &fLocalTrigger, kBufferSize);
    Info("MakeBranch", "Making Branch %s for Global Trigger\n",branchname);  
  }
  
  if (TreeT() && cRT ) {
    if (fRecTracks == 0x0)  fRecTracks = new TClonesArray("AliMUONTrack",10);
    fNrectracks = 0;
    sprintf(branchname,"%sTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) {  
      Info("MakeBranch","Branch %s is already in tree.",GetName());
      return ;
    }
    branch = TreeT()->Branch(branchname,&fRecTracks,kBufferSize);
    Info("MakeBranch","Making Branch %s for tracks \n",branchname);
  }  

  if (TreeP() && cRP ) {
    Info("MakeBranch","Making Branch for TreeP is not yet ready. \n");
  }
}

//____________________________________________________________________________
void AliMUONData::ResetDigits()
{
    //
    // Reset number of digits and the digits array for this detector
    //
    if (fDigits == 0x0) return;
    for ( int i=0;i<AliMUONConstants::NCh();i++ ) {
      if ((*fDigits)[i])    ((TClonesArray*)fDigits->At(i))->Clear();
      if (fNdigits)  fNdigits[i]=0;
    }
}
//______________________________________________________________________________
void AliMUONData::ResetHits()
{
  // Reset number of clusters and the cluster array for this detector
  fNhits   = 0;
  if (fHits) fHits->Clear();
}
//_______________________________________________________________________________
void AliMUONData::ResetRawClusters()
{
    // Reset number of raw clusters and the raw clust array for this detector
    //
  for ( int i=0;i<AliMUONConstants::NTrackingCh();i++ ) {
    if ((*fRawClusters)[i])    ((TClonesArray*)fRawClusters->At(i))->Clear();
    if (fNrawclusters)  fNrawclusters[i]=0;
  }
}
//_______________________________________________________________________________
void AliMUONData::ResetTrigger()
{
  //  Reset Local and Global Trigger 
  fNglobaltrigger = 0;
  if (fGlobalTrigger) fGlobalTrigger->Clear();
  fNlocaltrigger = 0;
  if (fLocalTrigger) fLocalTrigger->Clear();
}
//____________________________________________________________________________
void AliMUONData::ResetRecTracks()
{
  // Reset tracks information
  fNrectracks = 0;
  if (fRecTracks) fRecTracks->Clear();
}
//_____________________________________________________________________________
void AliMUONData::SetTreeAddress(Option_t* option)
{
  const char *cH   = strstr(option,"H");
  const char *cD   = strstr(option,"D");   // Digits branches in TreeD
  const char *cRC  = strstr(option,"RC");  // RawCluster branches in TreeR
  const char *cGLT = strstr(option,"GLT"); // Global and Local Trigger branches in TreeR
  const char *cRT  = strstr(option,"RT");  // Reconstructed Track in TreeT
  //const char *cRP  = strstr(option,"RP");  // Reconstructed Particle in TreeP
  
  // Set branch address for the Hits, Digits, RawClusters, GlobalTrigger and LocalTrigger Tree.
  char branchname[30];
  TBranch * branch = 0x0;

  //
  // Branch address for hit tree
  if ( TreeH() && cH ) {
    if (fHits == 0x0) fHits     = new TClonesArray("AliMUONHit",1000);
    fNhits =0;
  } 
  if (TreeH() && fHits && cH) {
    sprintf(branchname,"%sHits",GetName());  
    branch = TreeH()->GetBranch(branchname);
    if (branch) {
      Info("SetTreeAddress","(%s) Setting for Hits",GetName());
      branch->SetAddress(&fHits);
    }
    else { //can be invoked before branch creation
      Warning("SetTreeAddress","(%s) Failed for Hits. Can not find branch in tree.",GetName());
    }
  }
  
  //
  // Branch address for digit tree
  if ( TreeD() && cD) {
    if (fDigits == 0x0) { 
      fDigits = new TObjArray(AliMUONConstants::NCh());
      fNdigits= new Int_t[AliMUONConstants::NCh()];
      for (Int_t i=0; i<AliMUONConstants::NCh() ;i++) {
	fDigits->AddAt(new TClonesArray("AliMUONDigit",10000),i); 
	fNdigits[i]=0;
      }
    }
  }

  if (TreeD() && fDigits && cD) {
    for (int i=0; i<AliMUONConstants::NCh(); i++) {
      sprintf(branchname,"%sDigits%d",GetName(),i+1);
      branch = TreeD()->GetBranch(branchname);
      TClonesArray * digits = Digits(i);
      if (branch) branch->SetAddress( &digits );
      else Warning("SetTreeAddress","(%s) Failed for Digits Detection plane %d. Can not find branch in tree.",GetName(),i);
    }
  }
  
  //
  // Branch address for rawclusters, globaltrigger and local trigger tree
  if (TreeR() ) {
    if (fRawClusters == 0x0 && cRC) {
      fRawClusters = new TObjArray(AliMUONConstants::NTrackingCh());
      fNrawclusters= new Int_t[AliMUONConstants::NTrackingCh()];
      for (Int_t i=0; i<AliMUONConstants::NTrackingCh();i++) {
	fRawClusters->AddAt(new TClonesArray("AliMUONRawCluster",10000),i); 
	fNrawclusters[i]=0;
      }
    }
    if (fLocalTrigger == 0x0 && cGLT) {
      fLocalTrigger  = new TClonesArray("AliMUONLocalTrigger",234);
    }
    if (fGlobalTrigger== 0x0 && cGLT) {
        fGlobalTrigger = new TClonesArray("AliMUONGlobalTrigger",1); 
    }

  }
  if ( TreeR()  && fRawClusters && cRC) {
    for (int i=0; i<AliMUONConstants::NTrackingCh(); i++) {
      sprintf(branchname,"%sRawClusters%d",GetName(),i+1);
      if (fRawClusters) {
	branch = TreeR()->GetBranch(branchname);
	if (branch) branch->SetAddress(&((*fRawClusters)[i]));
	else Warning("SetTreeAddress","(%s) Failed for RawClusters Detection plane %d. Can not find branch in tree.",GetName(),i);
      }
    }
  }
  if ( TreeR()  && fLocalTrigger && cGLT) {
    sprintf(branchname,"%sLocalTrigger",GetName());
    branch = TreeR()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fLocalTrigger);
    else Warning("SetTreeAddress","(%s) Failed for LocalTrigger. Can not find branch in tree.",GetName());
  }
  if ( TreeR() && fGlobalTrigger && cGLT) {
    sprintf(branchname,"%sGlobalTrigger",GetName());
    branch = TreeR()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fGlobalTrigger);
    else Warning("SetTreeAddress","(%s) Failed for LocalTrigger. Can not find branch in tree.",GetName());
  }

  if ( TreeT() && fRecTracks && cRT ) {
    sprintf(branchname,"%sTrack",GetName());  
    branch = TreeT()->GetBranch(branchname);
    if (branch) branch->SetAddress(&fRecTracks);
    else Warning("SetTreeAddress","(%s) Failed for Tracks. Can not find branch in tree.",GetName());
  }
}
//_____________________________________________________________________________
