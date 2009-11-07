#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector3.h"


#include "AliESDVertex.h"
#include "AliESDEvent.h"

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"

#include "AliPerfAnalyzeInvPt.h"
#include "AliRecInfoCuts.h"
#include "AliMCInfoCuts.h"

#include "AliPerformancePtCalib.h"

using namespace std;

ClassImp(AliPerformancePtCalib)

//________________________________________________________________________
  AliPerformancePtCalib::AliPerformancePtCalib():
  AliPerformanceObject("AliPerformancePtCalib"),
  
 // option parameter for Analyse()
  fNThetaBins(0), 
  fNPhiBins(0),
  fRange(0),
  fExclRange(0),
  fFitGaus(0) ,
  // option for user defined 1/pt shift
  fShift(0),
  fDeltaInvP(0),
  
  //options for cuts
  fOptTPC(0),
  fESDcuts(0),
  fRefitTPC(0),
  fRefitITS(0),
  fDCAcut(0),

  fMinPt(0),
  fMaxPt(0),
  fMinNClustersTPC(0),
  fMaxChi2PerClusterTPC(0),
  fMaxDCAtoVertexXY(0),
  fMaxDCAtoVertexZ(0),

  fCutsRC(0),
  fCutsMC(0),

  fList(0),
  // histograms
  fHistInvPtTheta(0),
  fHistInvPtPhi(0),
  fHistPtTheta(0),
  fHistPtPhi(0),

  fHistPtShift0(0),
  fHistPrimaryVertexPosX(0),
  fHistPrimaryVertexPosY(0),
  fHistPrimaryVertexPosZ(0),
  fHistTrackMultiplicity(0),
  fHistTrackMultiplicityCuts(0),

  fHistTPCMomentaPosP(0),
  fHistTPCMomentaNegP(0),
  fHistTPCMomentaPosPt(0),
  fHistTPCMomentaNegPt(0),
  fHistUserPtShift(0),

  //esd track cuts
  fESDTrackCuts(0),
  
  // analysis folder 
  fAnalysisFolder(0)
{
   
  // Dummy constructor

  fShift = kFALSE;
   fDeltaInvP = 0.00;
  //options for cuts
  fOptTPC =  kTRUE;                      // read TPC tracks yes/no
  fESDcuts = kTRUE;                      // read ESD track cuts
  fRefitTPC = kFALSE;                    // require TPC refit
  fRefitITS = kFALSE;                    // require ITS refit
  fDCAcut = kTRUE;                       // apply DCA cuts
 

  fCutsRC = NULL;
  fCutsMC = NULL;

  fMinPt=0.15;  // GeV/c
  fMaxPt=1.e10; // GeV/c 
  fMinNClustersTPC = 50;
  fMaxChi2PerClusterTPC = 4.0;
  fMaxDCAtoVertexXY = 2.4; // cm
  fMaxDCAtoVertexZ  = 3.2; // cm
  
 // options for function Analyse()
  fFitGaus = kFALSE;
   fNThetaBins = 0;
   fNPhiBins =0  ;
   fRange =0;
   fExclRange =0;
   fFitGaus =0;
  
  Init();
}

//________________________________________________________________________
AliPerformancePtCalib::AliPerformancePtCalib(Char_t * name="AliPerformancePtCalib",Char_t* title ="AliPerformancePtCalib") :
   AliPerformanceObject(name,title),
   
   // option parameter for Analyse()
   fNThetaBins(0), 
   fNPhiBins(0),
   fRange(0),
   fExclRange(0),
   fFitGaus(0) ,
  fShift(0),
  fDeltaInvP(0),
   
   //options for cuts
   fOptTPC(0),
  fESDcuts(0),
   fRefitTPC(0),
  fRefitITS(0),
  fDCAcut(0),

  fMinPt(0),
  fMaxPt(0),
  fMinNClustersTPC(0),
  fMaxChi2PerClusterTPC(0),
  fMaxDCAtoVertexXY(0),
  fMaxDCAtoVertexZ(0),

  fCutsRC(0),
  fCutsMC(0),
   fList(0),
   
   // histograms
  fHistInvPtTheta(0),
  fHistInvPtPhi(0),
  fHistPtTheta(0),
  fHistPtPhi(0),

  fHistPtShift0(0),
  fHistPrimaryVertexPosX(0),
  fHistPrimaryVertexPosY(0),
  fHistPrimaryVertexPosZ(0),
  fHistTrackMultiplicity(0),
  fHistTrackMultiplicityCuts(0),

  fHistTPCMomentaPosP(0),
  fHistTPCMomentaNegP(0),
  fHistTPCMomentaPosPt(0),
  fHistTPCMomentaNegPt(0),
   fHistUserPtShift(0),	
   //esd track cuts														     
  fESDTrackCuts(0),
  
  // analysis folder 
  fAnalysisFolder(0)

  
{
  // Constructor
  fShift = kFALSE;
  fDeltaInvP = 0.00;
  //options for cuts
  fOptTPC =  kTRUE;                      // read TPC tracks yes/no
  fESDcuts = kTRUE;                      // read ESD track cuts
  fRefitTPC = kFALSE;                    // require TPC refit
  fRefitITS = kFALSE;                    // require ITS refit
  fDCAcut = kTRUE;                       // apply DCA cuts
 
  fCutsRC = NULL;
  fCutsMC = NULL;

  fMinPt=0.15;  // GeV/c
  fMaxPt=1.e10; // GeV/c 
  fMinNClustersTPC = 50;
  fMaxChi2PerClusterTPC = 4.0;
  fMaxDCAtoVertexXY = 2.4; // cm
  fMaxDCAtoVertexZ  = 3.2; // cm

  // options for Analyse()
 fFitGaus = kFALSE;
   fNThetaBins = 0;
   fNPhiBins =0  ;
   fRange =0;
   fExclRange =0;
   fFitGaus =0;
  
  Init();
}

AliPerformancePtCalib::~AliPerformancePtCalib() { 
//
// destructor
//
if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0; 
}

//________________________________________________________________________
void AliPerformancePtCalib::Init() 
{
  // Create histograms
  // Called once

   fList = new TList();
  // init folder
  fAnalysisFolder = CreateFolder("folderPt_TPC","Analysis Pt Resolution Folder");
  fList->Add(fAnalysisFolder);
  // Primary Vertex:
  fHistPrimaryVertexPosX       = new TH1F("fHistPrimaryVertexPosX", "Primary Vertex Position X;Primary Vertex Position X (cm);Events",100,-0.5,0.5);
  fList->Add(fHistPrimaryVertexPosX);
  fHistPrimaryVertexPosY       = new TH1F("fHistPrimaryVertexPosY", "Primary Vertex Position Y;Primary Vertex Position Y (cm);Events",100,-0.5,0.5);
  fList->Add(fHistPrimaryVertexPosY);
  fHistPrimaryVertexPosZ       = new TH1F("fHistPrimaryVertexPosZ", "Primary Vertex Position Z;Primary Vertex Position Z (cm);Events",200,-2.0,2.0);
  fList->Add(fHistPrimaryVertexPosZ);
  // Multiplicity:
  fHistTrackMultiplicity     = new TH1F("fHistTrackMultiplicity", "Multiplicity distribution;Number of tracks;Events", 250, 0, 250);
  fList->Add(fHistTrackMultiplicity);
  fHistTrackMultiplicityCuts = new TH1F("fHistTrackMultiplicityCuts", "Multiplicity distribution;Number of tracks after cuts;Events", 250, 0, 250);
  fList->Add(fHistTrackMultiplicityCuts);
 
  // momentum histos
  //pt shift 0 only needed if shift in 1/pt is applied
  fHistPtShift0 = new TH1F("fHistPtShift0","1/pt dN/pt vs. pt of ESD track  ",600,0.0,6.0);
  fList->Add(fHistPtShift0);
  fHistInvPtTheta = new TH2F("fHistInvPtTheta","#theta vs 1/pt ",900, -4.5, 4.5,300,0.0,3.0);
  fList->Add(fHistInvPtTheta);
  fHistInvPtPhi   = new TH2F("fHistInvPtPhi","#phi vs 1/pt",900, -4.5, 4.5,325,0.0,6.5);
  fList->Add(fHistInvPtPhi);
  fHistPtTheta    = new TH2F("fHistPtTheta"," #theta vs pt ",300, 0.0, 15.0,300,0.0,3.0);
  fList->Add(fHistPtTheta);
  fHistPtPhi      = new TH2F("fHistPtPhi"," #phi vs  pt ",300, 0.0,15.0,325,0.0,6.5);
  fList->Add(fHistPtPhi);
	     
   // mom test histos
  fHistTPCMomentaPosP  =  new TH2F("fHistTPCMomentaPosP","TPC p vs global esd track p pos",300,0.0,15.0,300,0.0,15.0);
  fList->Add(fHistTPCMomentaPosP);
  fHistTPCMomentaNegP  =  new TH2F("fHistTPCMomentaNegP","TPC p vs global esd track p neg",300,0.0,15.0,300,0.0,15.0);
  fList->Add(fHistTPCMomentaNegP);
  fHistTPCMomentaPosPt =  new TH2F("fHistTPCMomentaPosPt","TPC pt vs global esd track pt pos",300,0.0,15.0,300,0.0,15.0);
  fList->Add(fHistTPCMomentaPosPt);
  fHistTPCMomentaNegPt =  new TH2F("fHistTPCMomentaNegPt","TPC pt vs global esd track pt neg",300,0.0,15.0,300,0.0,15.0);
  fList->Add(fHistTPCMomentaNegPt);

  //user pt shift check
   fHistUserPtShift = new TH1F("fHistUserPtShift","user defined shift in 1/pt",100,-0.5,1.5);
   fList->Add(fHistUserPtShift);

  
  // esd track cuts  
  fESDTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  
  //fESDTrackCuts->DefineHistoqgrams(1);
  fESDTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDTrackCuts->SetRequireTPCRefit(kTRUE);
  fESDTrackCuts->SetAcceptKinkDaughters(kTRUE);
  fESDTrackCuts->SetMinNClustersTPC((Int_t)fMinNClustersTPC);
  fESDTrackCuts->SetMaxChi2PerClusterTPC(fMaxChi2PerClusterTPC);
  fESDTrackCuts->SetMaxDCAToVertexXY(fMaxDCAtoVertexXY);
  fESDTrackCuts->SetMaxDCAToVertexZ(fMaxDCAtoVertexZ);
  fESDTrackCuts->SetDCAToVertex2D(kTRUE);
  fESDTrackCuts->SetPtRange(fMinPt,fMaxPt); 
}

//________________________________________________________________________
void AliPerformancePtCalib::SetPtShift(Double_t shiftVal ) { 
   if(!(shiftVal==0)) { fShift=kTRUE; fDeltaInvP = shiftVal; } 
}

//________________________________________________________________________
void AliPerformancePtCalib::Exec(AliMCEvent*, AliESDEvent* const esdEvent, AliESDfriend*, Bool_t, Bool_t) 
{

 
  if (!esdEvent) {
    Printf("ERROR: Event not available");
    return;
  }

  if (!(esdEvent->GetNumberOfTracks())) {
    Printf(" PtCalib task: There is no track in this event");
    return;
  }

  if(fShift) fHistUserPtShift->Fill(fDeltaInvP);
  
  fHistTrackMultiplicity->Fill(esdEvent->GetNumberOfTracks());


  // read primary vertex info
  Double_t tPrimaryVtxPosition[3];
  //  Double_t tPrimaryVtxCov[3];
  const AliESDVertex *primaryVtx = esdEvent->GetPrimaryVertexTPC();
 
  tPrimaryVtxPosition[0] = primaryVtx->GetXv();
  tPrimaryVtxPosition[1] = primaryVtx->GetYv();
  tPrimaryVtxPosition[2] = primaryVtx->GetZv();
  
  fHistPrimaryVertexPosX->Fill(tPrimaryVtxPosition[0]);
  fHistPrimaryVertexPosY->Fill(tPrimaryVtxPosition[1]);
  fHistPrimaryVertexPosZ->Fill(tPrimaryVtxPosition[2]);


  //_fill histos for pt spectra and shift of transverse momentum
  Int_t count=0;
 
  for(Int_t j = 0;j<esdEvent->GetNumberOfTracks();j++){// track loop
    AliESDtrack *ESDTrack = esdEvent->GetTrack(j);
    if(!ESDTrack) continue;
    
    
    if(fESDcuts == kTRUE){
        if(!fESDTrackCuts->AcceptTrack(ESDTrack))continue;
    }
       
    //track cuts
    if(fRefitTPC) if(AddTPCcuts(ESDTrack)) continue;
    if(fRefitITS) if(AddITScuts(ESDTrack)) continue;
    if(fDCAcut)   if(AddDCAcuts(ESDTrack)) continue ;
      
    // fill histos
    if(fOptTPC == kTRUE){ //TPC tracks
      const AliExternalTrackParam *TPCTrack = ESDTrack->GetTPCInnerParam(); 
      if(!TPCTrack) continue;
      if(fabs(TPCTrack->Eta())> 0.8) continue;
      
      Double_t signedPt = TPCTrack->GetSignedPt();
      Double_t invPt = 0.0;
      if(!signedPt==0) {
	invPt = 1.0/signedPt;

	fHistPtShift0->Fill(fabs(signedPt));
	
	if(fShift == kTRUE ){
	  invPt += fDeltaInvP; //shift momentum for tests
	  if(!invPt==0) signedPt = 1.0/invPt;
	  else continue;
	}

	fHistInvPtTheta->Fill(invPt,TPCTrack->Theta());
	fHistInvPtPhi->Fill(invPt,TPCTrack->Phi());
	fHistPtTheta->Fill(fabs(signedPt),TPCTrack->Theta());
	fHistPtPhi->Fill(fabs(signedPt),TPCTrack->Phi());


	Double_t pTPC  = TPCTrack->GetP();
	Double_t pESD =	ESDTrack->GetP();
	Double_t ptESD  = ESDTrack->GetSignedPt();
	
	if(ESDTrack->GetSign()>0){//compare momenta ESD track and TPC track
	  fHistTPCMomentaPosP->Fill(fabs(pESD),fabs(pTPC));
	  fHistTPCMomentaPosPt->Fill(fabs(ptESD),fabs(signedPt));
 	}
 	else{
 	  fHistTPCMomentaNegP->Fill(fabs(pESD),fabs(pTPC));
 	  fHistTPCMomentaNegPt->Fill(fabs(ptESD),fabs(signedPt));
	}
	count++;
      }
      else continue;
    }
   
    else{// ESD tracks
      Double_t invPt = 0.0;
      Double_t signedPt = ESDTrack->GetSignedPt();
      if(!signedPt==0){
	invPt = 1.0/signedPt; 

	fHistPtShift0->Fill(fabs(signedPt));
	  
	if(fShift == kTRUE ){
	  invPt += fDeltaInvP;//shift momentum for tests
	  if(!invPt==0) signedPt = 1.0/invPt;
	  else continue;
	}
	fHistInvPtTheta->Fill(invPt,ESDTrack->Theta());
	fHistInvPtPhi->Fill(invPt,ESDTrack->Phi());
	fHistPtTheta->Fill(signedPt,ESDTrack->Theta());
	fHistPtPhi->Fill(signedPt,ESDTrack->Phi());
	
	count++;
      }
    }
  }
    
  fHistTrackMultiplicityCuts->Fill(count);
    
}    



//______________________________________________________________________________________________________________________
Bool_t AliPerformancePtCalib::AddTPCcuts(AliESDtrack *ESDTrack){
  
  Bool_t cut = kFALSE;
  
  if ((ESDTrack->GetStatus()&AliESDtrack::kTPCrefit)==0) cut=kTRUE; // TPC refit
  if (ESDTrack->GetTPCNcls()<50) cut=kTRUE; // min. nb. TPC clusters
  if(cut) return kTRUE;
  return kFALSE;
}
//______________________________________________________________________________________________________________________
Bool_t AliPerformancePtCalib::AddDCAcuts(AliESDtrack *ESDTrack){
  
  Bool_t cut = kFALSE;
  
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z and impact parameters:
  ESDTrack->GetImpactParameters(dca,cov);
  if(TMath::Abs(dca[0])>3. || TMath::Abs(dca[1])>3.) cut=kTRUE;
  if(ESDTrack->GetKinkIndex(0)>0) cut=kTRUE;
  if(cut) return kTRUE;
  return kFALSE;
}

//______________________________________________________________________________________________________________________
Bool_t AliPerformancePtCalib::AddITScuts(AliESDtrack *ESDTrack){

  Bool_t cut = kFALSE;
  
  if ((ESDTrack->GetStatus()&AliESDtrack::kITSrefit)==0) cut=kTRUE; // ITS refit
  Int_t clusterITS[200]; 
  if(ESDTrack->GetITSclusters(clusterITS)<2) cut=kTRUE;  // min. nb. ITS clusters //3
  
  if(cut) return kTRUE;
  return kFALSE;
}

//______________________________________________________________________________________________________________________

void AliPerformancePtCalib::Analyse()
{
  
  
  AliPerfAnalyzeInvPt *ana = new  AliPerfAnalyzeInvPt("AliPerfAnalyzeInvPt","AliPerfAnalyzeInvPt");
  
     
  TH1::AddDirectory(kFALSE);
 
  ana->SetProjBinsTheta(fThetaBins,fNThetaBins);
  ana->SetProjBinsPhi(fPhiBins,fNPhiBins);
  ana->SetMakeFitOption(fFitGaus,fExclRange,fRange);
  
  TObjArray *aFolderObj = new TObjArray;
  ana->StartAnalysis(fHistInvPtTheta,fHistInvPtPhi, aFolderObj);
  
  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
  if(ana) delete ana;
  
}

//______________________________________________________________________________________________________________________
TFolder* AliPerformancePtCalib::ExportToFolder(TObjArray * array) 
{
  // recreate folder every time and export objects to new one
  //
  AliPerformancePtCalib * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
    // get name and title from old folder
    name = folder->GetName();  
    title = folder->GetTitle();  

    // delete old one
    delete folder;

    // create new one
    newFolder = CreateFolder(name.Data(),title.Data());
    newFolder->SetOwner();

    // add objects to folder
    while(i < size) {
      newFolder->Add(array->At(i));
      i++;
    }
  }

  return newFolder;
}

//______________________________________________________________________________________________________________________
Long64_t AliPerformancePtCalib::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
    {
      AliPerformancePtCalib* entry = dynamic_cast<AliPerformancePtCalib*>(obj);
      if (entry == 0) continue; 
  
      fHistInvPtTheta->Add(entry->fHistInvPtTheta);
      fHistInvPtPhi->Add(entry-> fHistInvPtPhi);
      fHistPtTheta->Add(entry->fHistPtTheta);
      fHistPtPhi->Add(entry->fHistPtPhi);
  
      fHistPtShift0->Add(entry->fHistPtShift0);
      fHistPrimaryVertexPosX->Add(entry->fHistPrimaryVertexPosX);
      fHistPrimaryVertexPosY->Add(entry->fHistPrimaryVertexPosY);
      fHistPrimaryVertexPosZ->Add(entry->fHistPrimaryVertexPosZ);
      fHistTrackMultiplicity->Add(entry->fHistTrackMultiplicity);
      fHistTrackMultiplicityCuts->Add(entry->fHistTrackMultiplicityCuts);
  
      fHistTPCMomentaPosP->Add(entry->fHistTPCMomentaPosP);
      fHistTPCMomentaNegP->Add(entry->fHistTPCMomentaNegP);
      fHistTPCMomentaPosPt->Add(entry->fHistTPCMomentaPosPt);
      fHistTPCMomentaNegPt->Add(entry->fHistTPCMomentaNegPt);
  
      count++;
    }
  
  return count;
}

//______________________________________________________________________________________________________________________
TFolder* AliPerformancePtCalib::CreateFolder(TString name,TString title) { 
  // create folder for analysed histograms
  //
  TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
//______________________________________________________________________________________________________________________
void AliPerformancePtCalib::SetProjBinsPhi(const Double_t *phiBinArray,Int_t nphBins){
   
   fNPhiBins = nphBins;
  
   for(Int_t k = 0;k<fNPhiBins;k++){
      fPhiBins[k] = phiBinArray[k];
   }
   Printf("AliPerformancePtCalib: number of bins in phi set to %i",fNPhiBins);

}
//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalib::SetProjBinsTheta(const Double_t *thetaBinArray, Int_t nthBins){
  
   fNThetaBins = nthBins;
   for(Int_t k = 0;k<fNThetaBins;k++){
      fThetaBins[k] = thetaBinArray[k];
   }
   Printf("AliPerformancePtCalib: number of bins in theta set to %i",fNThetaBins);

}
//____________________________________________________________________________________________________________________________________________
void AliPerformancePtCalib::SetMakeFitOption(const Bool_t setGausFit, const Double_t exclusionR,const Double_t fitR ){

  
   fRange = fitR;
   fFitGaus = setGausFit;
   fExclRange  = exclusionR;
  
   if(fFitGaus) Printf("AliPerformancePtCalib: set MakeGausFit with fit range %2.3f and exclusion range in 1/pt: %2.3f",fRange,fExclRange);
   else  Printf("AliPerformancePtCalib: set standard polynomial fit with fit range %2.3f and exclusion range in 1/pt: %2.3f",fRange,fExclRange);
 
}
