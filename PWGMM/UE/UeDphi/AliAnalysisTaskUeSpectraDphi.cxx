/*   This macro produces: Delta phi correclation in different multiplicity classes and leading particle dE/dx vs leading particle momentum vs Multiplicity (near side, away side and transverse side) 
     Aditya Nath Mishra ICN-UNAM
     Please report bugs to: amishra@cern.ch / aditya.mishra@correo.nucleares.unam.mx 
     First version: 07/03/2019

 */

#include "AliAnalysisTaskUeSpectraDphi.h"

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
const Int_t ndPhiBins =18;
Double_t dPhiBinsOA[ndPhiBins] ={0.174533,0.349066,0.523599,0.698132,0.872665,1.0472,1.22173,1.39626,1.5708,1.74533,1.91986,2.0944,2.26893,2.44346,2.61799,2.79253,2.96706,3.14};
Double_t dPhiBinsSA[ndPhiBins+1] ={0,0.174533,0.349066,0.523599,0.698132,0.872665,1.0472,1.22173,1.39626,1.5708,1.74533,1.91986,2.0944,2.26893,2.44346,2.61799,2.79253,2.96706,3.14};

ClassImp(AliAnalysisTaskUeSpectraDphi)

//_____________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::AliAnalysisTaskUeSpectraDphi():
AliAnalysisTaskSE(),
  fESD(0x0),
  fEventCuts(0x0),
  fTrackFilter(0x0),
  fAnalysisType("ESD"),
  ftrigBit(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),  
  fEtaCut(0.8),
  fNcl(70),
  fTriggeredEventMB(-999),
  fListOfObjects(0x0),
  fisPS(kFALSE),
  isINEL0Rec(kFALSE),
  fZvtx_SPD(0x0),
  fisTracklet(kFALSE),
  fVtxBeforeCuts(0x0), 
  fVtxAfterCuts(0x0),
  fEvents(0x0),
  hpT(0x0),
  hEta(0x0),
  hPhi(0x0),
  hPtL(0x0),
  hPhiL(0x0),
  hEtaL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  hpTvsDphiOA(0x0),
  hpTvsDphiSA(0x0),
  fMultSelection(0x0),
  ftrackmult08(-999),   
  fv0mpercentile(-999)
{
  for(Int_t i=0; i<18; i++){
    hpTDphiBinsOA[i] = 0;
    hpTDphiBinsSA[i] = 0;
    hMultDphiBinsOA[i] = 0;
    hMultDphiBinsSA[i] = 0;
  }
  // Default constructor (should not be used)
}

//______________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::AliAnalysisTaskUeSpectraDphi(const char *name):
  AliAnalysisTaskSE(name),
  fESD(0x0),
  fEventCuts(0x0),
  fTrackFilter(0x0),
  fAnalysisType("ESD"),
  ftrigBit(0x0),
  fPileUpRej(kFALSE),
  fVtxCut(10.0),  
  fEtaCut(0.8),
  fNcl(70),
  fTriggeredEventMB(-999),
  fListOfObjects(0x0),
  fisPS(kFALSE),
  isINEL0Rec(kFALSE),
  fZvtx_SPD(0x0),
  fisTracklet(kFALSE),
  fVtxBeforeCuts(0x0), 
  fVtxAfterCuts(0x0),
  fEvents(0x0),
  hpT(0x0),
  hEta(0x0),
  hPhi(0x0),
  hPtL(0x0),
  hPhiL(0x0),
  hEtaL(0x0),
  hDphi(0x0),
  hRefMult08(0x0),
  hV0Mmult(0x0),
  hpTvsDphiOA(0x0),
  hpTvsDphiSA(0x0),
  fMultSelection(0x0),
  ftrackmult08(-999),   
  fv0mpercentile(-999)
{
  for(Int_t i=0; i<18; i++){
    hpTDphiBinsOA[i] = 0;
    hpTDphiBinsSA[i] = 0;
    hMultDphiBinsOA[i] = 0;
    hMultDphiBinsSA[i] = 0;
  }
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________

void AliAnalysisTaskUeSpectraDphi::Exit(const char *msg) {

  Printf("%s", msg);
  return;
}


//_____________________________________________________________________________
AliAnalysisTaskUeSpectraDphi::~AliAnalysisTaskUeSpectraDphi()
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
void AliAnalysisTaskUeSpectraDphi::UserCreateOutputObjects()
{

  const Int_t nPtBins      = 70; 
  Double_t PtBins[nPtBins+1] = {
    0.15,0.16,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,
    0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,
    1.6,1.7,1.8,1.9,2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4,
    4.5,5,5.5,6,6.5,7,8,9,10,11,12,13,14,15,16,18,20,22.0,
    24.0,26.0,28.0,32.0,36.0,42.0,50.0,60.0,80.0,100.0,130.0,
    160.0,200.0};

  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested 

  // Definition of trackcuts
  if(!fTrackFilter){	
    fTrackFilter = new AliAnalysisFilter("trackFilter2015");
    SetTrackCuts(fTrackFilter);
  }

  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();

  //
  // Histograms
  //  

  fEvents = 0;
  fEvents = new TH1D("fEvents","Number of events; Events; Counts", 10, 0, 10);
  fListOfObjects->Add(fEvents);
  
  fVtxBeforeCuts = 0;
  fVtxBeforeCuts = new TH1D("fVtxBeforeCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 300, -30, 30);
  fListOfObjects->Add(fVtxBeforeCuts);

  fVtxAfterCuts = 0;
  fVtxAfterCuts = new TH1D("fVtxAfterCuts", "Vtx distribution (before cuts); Vtx z [cm]; Counts", 300, -30, 30);
  fListOfObjects->Add(fVtxAfterCuts);

  hpT = new TH1D("hpT","",nPtBins,PtBins);
  fListOfObjects->Add(hpT);

  hEta = new TH1D("hEta","; #eta^{leading};counts",20,-1,1);
  fListOfObjects->Add(hEta);
  
  hPhi = new TH1D("hPhi", ";#phi (rad); count", 64,0,2.0*TMath::Pi());
  fListOfObjects->Add(hPhi);

  hPtL = 0;
  hPtL = new TH1D("hPtL",";#it{p}_{T}^{leading} (GeV/#it{c});counts",nPtBins,PtBins);
  fListOfObjects->Add(hPtL);
	
  hEtaL = 0;
  hEtaL = new TH1D("hEtaL","; #eta^{leading};counts",20,-1,1);
  fListOfObjects->Add(hEtaL);
	
  hPhiL = 0;
  hPhiL = new TH1D("hPhiL","; #phi^{leading} (rad);counts",64,0,2.0*TMath::Pi());
  fListOfObjects->Add(hPhiL);

  hDphi = 0;
  hDphi = new TH1D("hDphi","",64,-2*TMath::Pi(),2*TMath::Pi());
  fListOfObjects->Add(hDphi);

  hpTvsDphiOA = 0;
  hpTvsDphiOA = new TH2D("hpTvsDphiOA","p_{T} vs #Delta#phi (-#pi to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",314,-3.14,3.14,nPtBins,PtBins);
  fListOfObjects->Add(hpTvsDphiOA);
  
  hpTvsDphiSA = 0;
  hpTvsDphiSA = new TH2D("hpTvsDphiSA","p_{T} vs #Delta#phi  (0 to #pi);#Delta#phi (rad.);p_{T} (GeV/c)",314,0,3.14,nPtBins,PtBins);
  fListOfObjects->Add(hpTvsDphiSA);

  hRefMult08 = 0;
  hRefMult08 = new TH1D("hRefMult08","Multiplicity (-0.8 < #eta < 0.8);N_{ch};count",500,0,500);   
  fListOfObjects->Add(hRefMult08);
	
  hV0Mmult = 0;
  hV0Mmult = new TH1D("hV0Mmult","V0M ;V0M percentile;count",110,0,110);
  fListOfObjects->Add(hV0Mmult);

    for(Int_t i=0;i<18;i++){
     hpTDphiBinsOA[i] = 0;
     hpTDphiBinsOA[i] = new TH1D(Form("hpTdPhiBinsOA%d",i),"Charged particles (opening #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBins,PtBins);
     fListOfObjects->Add(hpTDphiBinsOA[i]);
     
     hpTDphiBinsSA[i] = 0;
     hpTDphiBinsSA[i] = new TH1D(Form("hpTdPhiBinsSA%d",i),"Charged particles (sliding #Delta#phi);#it{p}_{T} (GeV/c);count",nPtBins,PtBins);
     fListOfObjects->Add(hpTDphiBinsSA[i]);

      hMultDphiBinsOA[i] = 0;
      hMultDphiBinsOA[i] = new TH1D(Form("hMultdPhiBinsOA%d",i),"Charged particles (opening #Delta#phi);Multiplicity;count",500,0,500);
     fListOfObjects->Add(hMultDphiBinsOA[i]);
     
     hMultDphiBinsSA[i] = 0;
     hMultDphiBinsSA[i] = new TH1D(Form("hMultdPhiBinsSA%d",i),"Charged particles (sliding #Delta#phi);Multiplicity;count",500,0,500);
     fListOfObjects->Add(hMultDphiBinsSA[i]);
    }

  fEventCuts.AddQAplotsToList(fListOfObjects);

  PostData(1, fListOfObjects);

}

//______________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::UserExec(Option_t *)
{

  // -----------------------------------------------------
  //			 InputEvent
  // -----------------------------------------------------

  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  // -----------------------------------------------------
  //			 E S D
  // -----------------------------------------------------

  if (fAnalysisType == "ESD"){
    fESD = dynamic_cast<AliESDEvent*>(event);

    if(!fESD){
      Printf("%s:%d ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
      this->Dump();
      return;
    }
  }

  /************ BEGINNING OF EVENT SELECTION *******************/
  // Get trigger decision
  fTriggeredEventMB = 0; //init
  if((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & ftrigBit ){
    fTriggeredEventMB = 1;
  }
  
  Bool_t SPDvsClustersBG = kFALSE;
  
  AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
  if (!AnalysisUtils)
    {
      cout<<"------- No AnalysisUtils Object Found --------"<<AnalysisUtils<<endl;
      return;
    }
  else SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD); // We want NO BG

  Bool_t isNotPileUp = !fESD->IsPileupFromSPD(5,0.8); 
  Bool_t IncompleteDAQ = fESD->IsIncompleteDAQ(); // we want is not incomplete DAQ

  // vertex 
  const AliESDVertex * vertex    =    fESD->GetPrimaryVertex(); // tracks vertex, if not -> spd vertex, if not TPC vertex
  Bool_t isVtxGood = vertex->GetStatus() && selectVertex2015pp( fESD ,kTRUE,kFALSE,kTRUE); // requires Tracks and spd vertex, and Zconsistency of 5mm
  Double_t vertex_z = vertex->GetZ();
  Bool_t isVtxInZCut = (TMath::Abs(vertex_z) <= fVtxCut); // Zvtx in +- 10

  // Implement INEL>0
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  Bool_t isINEL0 = kFALSE;
  for (Int_t i = 0; i < mult->GetNumberOfTracklets(); ++i){ if (TMath::Abs(mult->GetEta(i)) < 1.) isINEL0 = kTRUE;}

  /********** IS PHYSICS SELECTION FLAG ****************************/	
  fisPS = fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp; 
  
  // recontructed INEL > 0 is PS + vtx + Zvtx inside +-10 ------
  isINEL0Rec = kFALSE;
  if ( isINEL0 && fisPS && isVtxGood && isVtxInZCut) isINEL0Rec = kTRUE;
  
  if (fisPS) fVtxBeforeCuts->Fill(vertex_z);              // VZ hack
  if (isINEL0Rec) fVtxAfterCuts->Fill(vertex_z);
  fEvents->Fill(0); // all events
  if(fTriggeredEventMB) fEvents->Fill(1); // triggered events
  if ( fTriggeredEventMB && !IncompleteDAQ ) fEvents->Fill(2); // trigger + IsIncompleteDAQ
  if (fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG) fEvents->Fill(3); // trigger + IsIncompleteDAQ + BG rejection
  if(fPileUpRej)
    {
      if(fTriggeredEventMB && !IncompleteDAQ && !SPDvsClustersBG && isNotPileUp) //     
	fEvents->Fill(4); // trigger + IsIncompleteDAQ + BG rejection + PileUp
    }
  if (fisPS)fEvents->Fill(5); //PS: trigger + IsIncompleteDAQ + BG rejection + PileUp + 1 Tracklet in eta +-1
  if (fisPS && isVtxGood) fEvents->Fill(6); //PS + GetPrimaryVertex
  if (isINEL0Rec) fEvents->Fill(7); //PS + GetPrimaryVertex + isVtxInZCut
	
  // -------------------------------------- multiplcity estimators section ------------------------------------------ //

  ftrackmult08 = -999;
  fv0mpercentile = -999;

  //ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8);     //reference
  ftrackmult08=AliESDtrackCuts::GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);     //tracklets
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

  if (isINEL0Rec) {MakeAnalysis(fESD);}
  cout<<"hello!!!"<<endl;

  // Post output data.
  PostData(1, fListOfObjects);

}
//________________________________________________________________________
void AliAnalysisTaskUeSpectraDphi::MakeAnalysis(AliESDEvent* fESD){


  // selection on leading particle
  Double_t pt_leading    = 0;
  Double_t p_leading    = 0;
  Double_t eta_leading    = 0;
  Double_t phi_leading    = 0;

  Int_t    i_leading = 0;
 
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    AliESDtrack* esdTrack = fESD->GetTrack(i);

    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t momentum = esdTrack->P();
    Double_t pt       = esdTrack->Pt();
	
    if(TMath::Abs(eta) > fEtaCut) continue;
    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;
    if(pt<0.15) continue;
    Short_t ncl = esdTrack->GetTPCsignalN();
    if(ncl<fNcl) continue;
	
    if(pt>pt_leading){
      pt_leading      = pt;
      p_leading       = momentum;
      eta_leading     = eta;
      phi_leading     = phi;
      i_leading = i;
    }

    hpT->Fill(pt);
    hEta->Fill(eta);
    hPhi->Fill(phi);

  }// end loop over tracks

  if(pt_leading<0.15)
    return;

  hPtL->Fill(pt_leading);
  hEtaL->Fill(eta_leading);
  hPhiL->Fill(phi_leading);

  Int_t mult = 0;
  Int_t multOA[18];
  Int_t multSA[18];

  for(Int_t j=0;j<ndPhiBins;j++){
    multOA[j] =0;
    multSA[j]=0;
  }
    
  for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    AliESDtrack* esdTrack = fESD->GetTrack(i);
    Double_t eta      = esdTrack->Eta();
    Double_t phi      = esdTrack->Phi();
    Double_t pt       = esdTrack->Pt();

    if(TMath::Abs(eta) > fEtaCut) continue;
    //quality cuts, standard 2015 track cuts
    if(!fTrackFilter->IsSelected(esdTrack)) continue;
    if(pt<0.15) continue;
    Short_t ncl = esdTrack->GetTPCsignalN();
    if(ncl<fNcl) continue;
    mult++;
		
    Double_t DPhiOA = DeltaPhi( phi, phi_leading );
    Double_t DPhiSA = TMath::Abs(DPhiOA);
  
    hDphi->Fill(DPhiOA);
    hpTvsDphiOA->Fill(DPhiOA,pt);
    hpTvsDphiSA->Fill(DPhiSA,pt);

    for(Int_t j=0;j<ndPhiBins;j++){
      if(TMath::Abs(DPhiOA)<= dPhiBinsOA[j]){
	hpTDphiBinsOA[j]->Fill(pt);
	 multOA[j]++;
      }
      if(DPhiSA > dPhiBinsSA[j] && DPhiSA <= dPhiBinsSA[j+1]){
	hpTDphiBinsSA[j]->Fill(pt);
	multSA[j]++;
      }
    }
  }// end loop over tracks
   for(Int_t k=0;k<ndPhiBins;k++){
      hMultDphiBinsOA[k]->Fill(multOA[k]);
      hMultDphiBinsSA[k]->Fill(multSA[k]);
    }
}

//______________________________________________________________________
Bool_t AliAnalysisTaskUeSpectraDphi::selectVertex2015pp(AliESDEvent *esd,
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
Bool_t AliAnalysisTaskUeSpectraDphi::IsGoodSPDvertexRes(const AliESDVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04 && spdVertex->GetZRes()<0.25)) return kFALSE;
  return kTRUE;
}

//____________________________________________________________
void AliAnalysisTaskUeSpectraDphi::SetTrackCuts(AliAnalysisFilter* fTrackFilter){

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
  
  fTrackFilter->AddCuts(esdTrackCuts);

}

Double_t AliAnalysisTaskUeSpectraDphi::DeltaPhi(Double_t phia, Double_t phib,
					    Double_t rangeMin, Double_t rangeMax)
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
