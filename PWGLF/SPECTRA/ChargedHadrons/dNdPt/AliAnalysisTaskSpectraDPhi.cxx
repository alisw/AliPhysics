/*   This macro produces: Pt spectra for different Delta phi 
     Aditya Nath Mishra Wigner RCP, Budapest
     Please report bugs to: amishra@cern.ch  
     First version: 05/05/2020

*/

#include "AliAnalysisTaskSpectraDPhi.h"

// ROOT includes
#include <TList.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2D.h>
#include <TProfile.h>
#include <THnSparse.h>
#include <TFile.h>


// AliRoot includes
#include <AliAnalysisManager.h>
#include <AliAnalysisFilter.h>
#include <AliESDInputHandler.h>
#include <AliESDEvent.h>
#include <AliEventCuts.h>
#include <AliESDVertex.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrackCuts.h>
#include <AliESDVZERO.h>
#include <AliAODVZERO.h>
#include <TTreeStream.h>
#include <AliHeader.h>
#include <AliAnalysisUtils.h>
#include <AliMultiplicity.h>
#include <AliMultSelection.h>
#include <AliMultSelectionTask.h>
#include <AliAODInputHandler.h> 
#include <AliAODHandler.h> 
#include <AliAODVertex.h>
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliDataFile.h>


#include <iostream>
using namespace std;

const Double_t pi = 3.1415926535897932384626433832795028841971693993751058209749445;

const Int_t nNchBins = 200;
Double_t NchBins[nNchBins+1]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
			      21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
			      39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,
			      57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,
			      75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,
			      93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,
			      109,110,111,112,113,114,115,116,117,118,119,120,121,122,
			      123,124,125,126,127,128,129,130,131,132,133,134,135,136,
			      137,138,139,140,141,142,143,144,145,146,147,148,149,150,
			      151,152,153,154,155,156,157,158,159,160,161,162,163,164,
			      165,166,167,168,169,170,171,172,173,174,175,176,177,178,
			      179,180,181,182,183,184,185,186,187,188,189,190,191,192,
			      193,194,195,196,197,198,199,200};


const Int_t nPtBining = 76;
Double_t PtBining[nPtBining+1] = {0.,0.05,0.1,0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
			      0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,
			      1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
			      4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20,22.0,
			      24.0,26.0,28.0,32.0,36.0,42.0,50.0,60.0,80.0,100.0,130.0,
			      160.0,200.0,300,400,500};

const Int_t ndPhiBins =18;
Double_t dPhiBins[ndPhiBins] ={0.174533,0.349066,0.523599,0.698132,0.872665,1.0472,1.22173,1.39626,1.5708,1.74533,1.91986,2.0944,2.26893,2.44346,2.61799,2.79253,2.96706,3.14};

Double_t dPhiBins2[ndPhiBins+1] ={0,0.174533,0.349066,0.523599,0.698132,0.872665,1.0472,1.22173,1.39626,1.5708,1.74533,1.91986,2.0944,2.26893,2.44346,2.61799,2.79253,2.96706,3.14};


const Int_t nV0Mbins = 10;
Double_t V0Mbins[nV0Mbins+1]={0.00, 1.00, 5.00, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100};

ClassImp(AliAnalysisTaskSpectraDPhi)

//_____________________________________________________________________________
AliAnalysisTaskSpectraDPhi::AliAnalysisTaskSpectraDPhi():
AliAnalysisTaskSE(),
  fEventCuts(0x0),
  fESD(0x0),
  fAOD(0x0),
  fTrackFilter(0x0),
  fTrackFilterDCA(0x0),
// fTrackFilterMatchEff(0x0),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),
  ftrigBit(0x0),
  fRandom(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),  
  fEtaCut(0.8),  
  fLeadMin(0.0),
  fTriggeredEventMB(-999),
  fRun(-999),
  fEventId(-999),
  fdcaxy(-999),
  fdcaz(-999),
  fisPS(kFALSE),
  fListOfObjects(0),
  fEventCounter(0x0),
  fEvents(0x0),
  fVtxBeforeCuts(0x0), 
  fVtxAfterCuts(0x0),
  isINEL0Rec(kFALSE),
  hpT(0x0),
  hPhi(0x0),
  hEta(0x0),		
  hPtL(0x0),
  hEtaL(0x0),
  hPhiL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  hpTvsDphiOA(0x0),
  hpTvsDphiSA(0x0),
  fMultSelection(0x0),
  ftrackmult08(-999),   
  fv0mpercentile(-999)
{
  // Default constructor (should not be used)
   for(Int_t i=0;i<ndPhiBins;i++){
     hpTDphiBinsOA[i]  = 0;
     hpTDphiBinsSA[i]  = 0;
     hMultDphiBinsOA[i] = 0;
     hMultDphiBinsSA[i] = 0;
   }
}

//______________________________________________________________________________
AliAnalysisTaskSpectraDPhi::AliAnalysisTaskSpectraDPhi(const char *name):
  AliAnalysisTaskSE(name),
  fEventCuts(0x0),
  fESD(0x0),
  fAOD(0x0),
  fTrackFilter(0x0),
  fTrackFilterDCA(0x0),
  // fTrackFilterMatchEff(0x0),
  fAnalysisType("ESD"),
  fAnalysisMC(kFALSE),
  ftrigBit(0x0),
  fRandom(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),  
  fEtaCut(0.8),  
  fLeadMin(5.0),
  fTriggeredEventMB(-999),
  fRun(-999),
  fEventId(-999),
  fdcaxy(-999),
  fdcaz(-999),
  fisPS(kFALSE),
  fListOfObjects(0),
  fEventCounter(0x0),
  fEvents(0x0),
  fVtxBeforeCuts(0x0), 
  fVtxAfterCuts(0x0),
  isINEL0Rec(kFALSE),
  hpT(0x0),
  hPhi(0x0),
  hEta(0x0),		
  hPtL(0x0),
  hEtaL(0x0),
  hPhiL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  hpTvsDphiOA(0x0),
  hpTvsDphiSA(0x0),
  fMultSelection(0x0),
  ftrackmult08(-999),   
  fv0mpercentile(-999)
{ 
   for(Int_t i=0;i<ndPhiBins;i++){
     hpTDphiBinsOA[i]  = 0;
     hpTDphiBinsSA[i]  = 0;
     hMultDphiBinsOA[i] = 0;
     hMultDphiBinsSA[i] = 0;
   }
   
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskSpectraDPhi::Exit(const char *msg) {

  Printf("%s", msg);
  return;
}


//_____________________________________________________________________________
AliAnalysisTaskSpectraDPhi::~AliAnalysisTaskSpectraDPhi()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fListOfObjects;
    fListOfObjects = 0x0;
  }

}

//______________________________________________________________________________
void AliAnalysisTaskSpectraDPhi::UserCreateOutputObjects()
{

  fRandom = new TRandom(0); // 0 means random seed

  //OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  // Definition of trackcuts
  if(!fTrackFilter){	
    fTrackFilter = new AliAnalysisFilter("trackFilter2015");
    SetTrackCuts(fTrackFilter);
  }

  //
  // Histograms
  //  

  //1D hists
  fEvents        = 0;
  fEventCounter  = 0;
  fVtxBeforeCuts = 0;   
  fVtxAfterCuts  = 0;
  hpT            = 0;
  hPhi           = 0;
  hEta           = 0;
  hPtL           = 0;
  hPhiL          = 0;
  hEtaL          = 0;
  hDphi          = 0;
  hRefMult08     = 0;
  hV0Mmult       = 0;
  hpTvsDphiOA    = 0;
  hpTvsDphiSA    = 0;

  fEvents        = new TH1I("fEvents","Number of analyzed events; Events; Counts", 6, 0, 6);
  fEventCounter  = new TH1D("fEventCounter","fEventCounter",10,0,10);
  fVtxBeforeCuts = new TH1D("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 12, -30, 30);
  fVtxAfterCuts  = new TH1D("fVtxAfterCuts", "Vtx distribution (after cuts); Vtx z [cm]; Counts", 12, -30, 30);
  hpT            = new TH1D("hpT","p_{T} Spectra",nPtBining,PtBining);
  hPhi           = new TH1D("hPhi", ";#phi (rad); count", 64,0,2.0*TMath::Pi());
  hEta           = new TH1D("hEta", "#eta Distribution ;#eta; count", 30, -1.5, 1.5);
  hPtL           = new TH1D("hPtL",";#it{p}_{T}^{leading} (GeV/#it{c});counts",nPtBining,PtBining);
  hPhiL          = new TH1D("hPhiL","; #phi^{leading} (rad);counts",64,0,2.0*TMath::Pi());
  hEtaL          = new TH1D("hEtaL","; #eta^{leading};counts",30,-1.5,1.5);
  hDphi          = new TH1D("hDphi","; d#phi Distribution;counts",64,-(TMath::Pi())/2.0,3.0*(TMath::Pi())/2.0);
  hRefMult08     = new TH1D("hRefMult08","Multiplicity (-0.8 < #eta < 0.8);Reference Multiplicity;count",nNchBins,NchBins);
  hV0Mmult       = new TH1D("hV0Mmult","V0M ;V0M percentile;count",nNchBins,NchBins);
  hpTvsDphiOA    = new TH2D("hpTvsDphiOA","p_{T} vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",314,-3.14,3.14,nPtBining,PtBining);
  hpTvsDphiSA    = new TH2D("hpTvsDphiSA","p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",314,0,3.14,nPtBining,PtBining);

  fListOfObjects->Add(fEvents);
  fListOfObjects->Add(fVtxBeforeCuts);
  fListOfObjects->Add(fVtxAfterCuts);
  fListOfObjects->Add(fEventCounter);
  fListOfObjects->Add(hpT);
  fListOfObjects->Add(hPhi);
  fListOfObjects->Add(hEta);
  fListOfObjects->Add(hPtL);	
  fListOfObjects->Add(hPhiL);
  fListOfObjects->Add(hEtaL);
  fListOfObjects->Add(hDphi);
  fListOfObjects->Add(hRefMult08);
  fListOfObjects->Add(hV0Mmult);

  for(Int_t i=0;i<ndPhiBins;i++){
  hpTDphiBinsOA[i]  = 0;
  hpTDphiBinsOA[i] = new TH1D(Form("hpTdPhiBinsOA%d",i),"Charged particles (Opening #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBining,PtBining);
  fListOfObjects->Add(hpTDphiBinsOA[i]);
  
  hpTDphiBinsSA[i]  = 0;
  hpTDphiBinsSA[i] = new TH1D(Form("hpTdPhiBinsSA%d",i),"Charged particles (Sliding #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBining,PtBining);
  fListOfObjects->Add(hpTDphiBinsSA[i]);
  
  hMultDphiBinsOA[i] = 0;
  hMultDphiBinsOA[i] = new TH1D(Form("hMultdPhiBinsOA%d",i),"Charged particles  (Opening #Delta#phi);Multiplicity;count",1000,0,1000);
  fListOfObjects->Add(hMultDphiBinsOA[i]);
  
  hMultDphiBinsSA[i] = 0;
  hMultDphiBinsSA[i] = new TH1D(Form("hMultdPhiBinsSA%d",i),"Charged particles  (Sliding #Delta#phi);Multiplicity;count",1000,0,1000);
  fListOfObjects->Add(hMultDphiBinsSA[i]);
  }
  
  fEventCuts.AddQAplotsToList(fListOfObjects);

  PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskSpectraDPhi::UserExec(Option_t *)
{
  
  // -----------------------------------------------------
  //			 Main loop
  //    First we make sure that we have valid input(s)!
  // -----------------------------------------------------
  
  // Main loop
  //
  // First we make sure that we have valid input(s)!
  //


  AliVEvent *event = InputEvent();
  if (!event) 
    {
      Error("UserExec", "Could not retrieve event");
      return;
    }

  if (fAnalysisMC)
    {
      fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
      if (!fMCEvent)
	{ 
	  cout<<"Could not retrieve MC event"<<endl;
	  return;
	}

      AliHeader *header = fMCEvent->Header();
      if(!header) {AliDebug( AliLog::kError , "Header not avaible" ); return; }
			
    }


  fESD = dynamic_cast<AliESDEvent*>(event);
  if(!fESD)
    {
      Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }    

  /************ BEGINNING OF EVENT SELECTION *******************/

  // Get trigger decision
  fTriggeredEventMB = 0; //init    
  if (!fAnalysisMC){  // for data check event selection as well
    if((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & ftrigBit )
      {
	fTriggeredEventMB = 1;  //event triggered as minimum bias
      }
  }
  else {
    if (ftrigBit) {
      fTriggeredEventMB = 1;
    }
  }
  Bool_t SPDvsClustersBG = kFALSE;

  AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
  if (!AnalysisUtils)
    {
      cout<<"------- No AnalysisUtils Object Found --------"<<AnalysisUtils<<endl;
      return;
    }
  else
    SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD); // We want NO BG
	
  Bool_t isNotPileUp = !fESD->IsPileupFromSPD(5,0.8); 
  Bool_t IncompleteDAQ = fESD->IsIncompleteDAQ(); // we want is not incomplete DAQ

  // vertex 
  const AliESDVertex * vertex    =    fESD->GetPrimaryVertex(); // tracks vertex, if not -> spd vertex, if not TPC vertex
  Bool_t isVtxGood = vertex->GetStatus() && selectVertex2015pp( fESD ,kTRUE,kFALSE,kTRUE); // requires Tracks and spd vertex, and Zconsistency of 5mm
  double vertex_z = vertex->GetZ();
  Bool_t isVtxInZCut = (TMath::Abs(vertex_z)   <= fVtxCut); // Zvtx in +- 10

 
 
  // Implement INEL>0
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  Bool_t isINEL0 = kFALSE;
  for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i) {
    if (TMath::Abs(mult->GetEta(i)) < 1.) isINEL0 = kTRUE;
  }
  
  /********** IS PHYSICS SELECTION FLAG ****************************/
	
  fisPS = fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp;
  //fisPS = !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp; 
		

  // recontructed INEL > 0 is PS + vtx + Zvtx inside +-10 ------
  isINEL0Rec = kFALSE;
  if ( isINEL0 && fisPS && isVtxGood && isVtxInZCut) isINEL0Rec = kTRUE;

  if (fisPS) // VZ hack
    fVtxBeforeCuts->Fill(vertex_z);
  if (isINEL0Rec)
    fVtxAfterCuts->Fill(vertex_z);

  fEventCounter->Fill(0); // all events

  if(fTriggeredEventMB)  
    fEventCounter->Fill(1); // triggered events

  if ( fTriggeredEventMB && !IncompleteDAQ )  
    fEventCounter->Fill(2); // trigger + IsIncompleteDAQ

  if (fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG)   
    fEventCounter->Fill(3); // trigger + IsIncompleteDAQ + BG rejection

  if(fPileUpRej)
    {
      if(fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp) //     
	fEventCounter->Fill(4); // trigger + IsIncompleteDAQ + BG rejection + PileUp
    }
  if (fisPS)  
    fEventCounter->Fill(5); //PS: trigger + IsIncompleteDAQ + BG rejection + PileUp + 1 Tracklet in eta +-1
		
  if (fisPS && isVtxGood)
    fEventCounter->Fill(6); //PS + GetPrimaryVertex
	
  if (isINEL0Rec)
    fEventCounter->Fill(7); //PS + GetPrimaryVertex + isVtxInZCut

  // -------------------------------------- multiplcity estimators section ------------------------------------------ //

  ftrackmult08 = -999;
  fv0mpercentile = -999;

  ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);     //reference
  // ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);     //tracklets
  //ftrackmult08 = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD); //Combined estimator

  hRefMult08->Fill(ftrackmult08);

  fMultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection"); // Esto es para 13 TeV
  if (!fMultSelection)
    cout<<"------- No AliMultSelection Object Found --------"<<fMultSelection<<endl;
  else
    fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  hV0Mmult->Fill(fv0mpercentile);
  
  cout<<"------- V0M mult ==  "<<fv0mpercentile<<"--------"<<endl;
  
   // ------------------------------------------ end of mult estimators -------------------------------------------------//

   if (!fEventCuts.AcceptEvent(event)) {
    PostData(1, fListOfObjects);
    return;
  }

  if(isINEL0Rec && !fAnalysisMC) {
    AnalysisSpectra(fESD);
  }
   
  //  cout<<"hello!!!"<<endl;

  // Post output data.
  PostData(1, fListOfObjects);

}
//________________________________________________________________________
void AliAnalysisTaskSpectraDPhi::AnalysisSpectra(AliESDEvent* esdEvent){

  fEvents->Fill(0);
  
  fRun  = esdEvent->GetRunNumber();
  fEventId = 0;
  if(esdEvent->GetHeader()) fEventId = GetEventIdAsLong(esdEvent->GetHeader());

  const Int_t nESDTracks = fESD->GetNumberOfTracks();


  // selection on leading particle
  Double_t pt_leading     = 0;
  Double_t eta_leading    = 0;
  Double_t phi_leading    = 0;
  Int_t    i_leading      = 0;

  Double_t mult    = 0;
  Double_t multDphiBinsOA[100];
  Double_t multDphiBinsSA[100];

  for(Int_t iT = 0; iT < nESDTracks; iT++) {
    
    AliVParticle* part = fESD->GetTrack(iT);
    Double_t eta  = part->Eta();
    Double_t pt   = part->Pt();
    Double_t phi  = part->Phi();

    if( TMath::Abs(eta) > fEtaCut ) continue;
    // if( !TMath::Abs(pt) > 0.15 ) continue;
    if(pt<0.15) continue;
    if( !fTrackFilter->IsSelected(part) ) continue;
    mult++;
    if(pt>pt_leading){
      pt_leading      = pt;
      eta_leading     = eta;
      phi_leading     = phi;
      i_leading = iT;
    }
    
    hPhi->Fill(phi);
    hpT->Fill(pt);  
    hEta->Fill(eta);
  }

  if(pt_leading<fLeadMin) return;
  hPtL->Fill(pt_leading);
  hEtaL->Fill(eta_leading);
  hPhiL->Fill(phi_leading);

   for(Int_t iT = 0; iT < nESDTracks; iT++) {

     // exclude the auto-correlation
     //  if(iT==i_leading) continue;

    AliVParticle* part = fESD->GetTrack(iT);
    Double_t eta  = part->Eta();
    Double_t pt   = part->Pt();
    Double_t phi  = part->Phi();

    if( TMath::Abs(eta) > fEtaCut ) continue;
    // if( !TMath::Abs(pt) > 0.15 ) continue;
    if(pt<0.15) continue;
    if( !fTrackFilter->IsSelected(part) ) continue;

    Double_t DPhi = DeltaPhi( phi, phi_leading );
    Double_t DeltaPhi = TMath::Abs(DPhi);
    
    hDphi->Fill(DPhi);
    hpTvsDphiOA->Fill(DPhi,pt);
    hpTvsDphiSA->Fill(DeltaPhi,pt);
    
    for(Int_t j=0;j<ndPhiBins;j++){
      if(TMath::Abs(DPhi)<= dPhiBins[j]){
	multDphiBinsOA[j]++;
	hpTDphiBinsOA[j]->Fill(pt);
      }

      if(DPhi > dPhiBins2[j] && DPhi <= dPhiBins2[j+1]){
	multDphiBinsSA[j]++;
	hpTDphiBinsSA[j]->Fill(pt);
      }
    }
    
    for(Int_t k=0;k<ndPhiBins;k++){
      hMultDphiBinsOA[k]->Fill(multDphiBinsOA[k]);
      hMultDphiBinsSA[k]->Fill(multDphiBinsSA[k]);
    }
   }
   
   
}// end loop over tracks

   //____________________________________________________________
   void AliAnalysisTaskSpectraDPhi::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts;
  // TPC
  esdTrackCuts->SetMinNCrossedRowsTPC(70);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);

  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
  // 7*(0.0015+0.0050/pt^1.1)
  esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");

  esdTrackCuts->SetMaxDCAToVertexZ(2);
  esdTrackCuts->SetDCAToVertex2D(kFALSE);
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);

  esdTrackCuts->SetMaxChi2PerClusterITS(36);

  /*
    AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
    esdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);//
    esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
    esdTrackCuts->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);//
    esdTrackCuts->SetMaxChi2PerClusterTPC(4);//
    esdTrackCuts->SetAcceptKinkDaughters(kFALSE);//
    esdTrackCuts->SetRequireTPCRefit(kTRUE);//
    esdTrackCuts->SetRequireITSRefit(kTRUE);//
    esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
    AliESDtrackCuts::kAny);//
    esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");//
    esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);//
    esdTrackCuts->SetMaxDCAToVertexZ(2);//
    esdTrackCuts->SetDCAToVertex2D(kFALSE);//
    esdTrackCuts->SetRequireSigmaToVertex(kFALSE);//
    esdTrackCuts->SetMaxChi2PerClusterITS(36);//
  */
  fTrackFilter->AddCuts(esdTrackCuts);

}

Double_t AliAnalysisTaskSpectraDPhi::DeltaPhi(Double_t phia, Double_t phib)
{
  Double_t dphi = -999;
  Double_t pi = TMath::Pi();

  if (phia < 0)         phia += 2*pi;
  else if (phia > 2*pi) phia -= 2*pi;
  if (phib < 0)         phib += 2*pi;
  else if (phib > 2*pi) phib -= 2*pi;
  dphi = phib - phia;
  if (dphi > pi)        dphi -= 2*pi;
  else if (dphi < -pi)  dphi += 2*pi;
  
  return dphi;
}

//_____________________________________________________________________________
ULong64_t AliAnalysisTaskSpectraDPhi::GetEventIdAsLong(AliVHeader* header)
{
  // To have a unique id for each event in a run!
  // Modified from AliRawReader.h
  return ((ULong64_t)header->GetBunchCrossNumber()+
	  (ULong64_t)header->GetOrbitNumber()*3564+
	  (ULong64_t)header->GetPeriodNumber()*16777215*3564);
}

//______________________________________________________________________
Bool_t AliAnalysisTaskSpectraDPhi::selectVertex2015pp(AliESDEvent *esd,
						   Bool_t checkSPDres, //enable check on vtx resolution
						   Bool_t requireSPDandTrk, //ask for both trk and SPD vertex 
						   Bool_t checkProximity) //apply cut on relative position of spd and trk verteces 
{

  if (!esd) return kFALSE;
  
  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
  const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
  Bool_t hasSPD = spdVertex->GetStatus();
  Bool_t hasTrk = trkVertex->GetStatus();
 
  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;
  
  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD) return kFALSE;
    //on demand check the spd vertex resolution and reject if not satisfied
    if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
      if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE; 
    }
  }
  return kTRUE;
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskSpectraDPhi::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}
