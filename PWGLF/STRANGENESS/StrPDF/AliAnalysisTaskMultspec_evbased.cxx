class TTree;
class TParticle;
class TVector3;

class AliESDVertex;
class AliESDv0;
//class AliESDcascade;
class AliAODVertex;
class AliAODv0;
class AliAODcascade;

#include <Riostream.h>
#include "TFile.h"
#include "TCanvas.h"
#include "THistManager.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliPID.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliESDcascade.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskMultspec.h"


#include "AliAnalysisTaskMultspec_evbased.h"


ClassImp(AliAnalysisTaskMultspec_evbased)

AliAnalysisTaskMultspec_evbased::AliAnalysisTaskMultspec_evbased() : AliAnalysisTaskSE(),
fHistos(nullptr),
fpart(666),
fPIDResponse(0), 
fEventCuts(0),
fRejectPileupEvts(kTRUE),
fRemovePythiaGen(kTRUE),
fisMC(kFALSE),
ffillV0(nullptr),
ffillCasc(nullptr)
{
  //default constructor
}

AliAnalysisTaskMultspec_evbased::AliAnalysisTaskMultspec_evbased(const char *name, int particle, TString lExtraOptions, bool ismc, bool removePythiaGen) : AliAnalysisTaskSE(name),
fHistos(nullptr),
fpart(particle),
fPIDResponse(0),
fEventCuts(0),
fRejectPileupEvts(kTRUE),
fRemovePythiaGen(kTRUE),
fisMC(kFALSE),
ffillV0(nullptr),
ffillCasc(nullptr)
{

  //SetAddCasc(AddCasc);
  SetIsMC(ismc);
  SetParticle(particle);
    
  //Standard output
  DefineOutput(1, TList::Class()); // Miscellaneous Histograms
  //DefineOutput(2, TTree::Class()); // Output Tree
}

AliAnalysisTaskMultspec_evbased::~AliAnalysisTaskMultspec_evbased()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    if (fHistos) {
        delete fHistos;
        fHistos = 0x0;
    }
 
}

//________________________________________________________________________
void AliAnalysisTaskMultspec_evbased::UserCreateOutputObjects()
{

  ffillV0 = new V0filler_Mult;
  ffillCasc = new Cascfiller_Mult;

  //miscellaneous histograms
  fHistos = new THistManager("fHistos");
  fHistos->CreateTH1("hNevts", "", 1, 0, 1);  // #events histo

  fHistos->CreateTH1("ImassK0sCheck","",1000,0.44,0.56, "s");
  fHistos->CreateTH1("ImassLamCheck","",1000,1.08,1.15, "s");
  fHistos->CreateTH1("ImassXiCheck","",1000,1.26,1.38, "s");
  fHistos->CreateTH1("ImassOmCheck","",1000,1.62,1.73, "s");
  fHistos->CreateTH1("RecoSpectrum","",20,0,20, "s");
  fHistos->CreateTH3("CorrMultNstrange", "CorrMultNstrange", 100, 0, 100, 100, 0, 100, 200, 0, 200, "s");  //event selection code, for debug

  // PID Setup
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  inputHandler->SetNeedField();
  
  //fEventCuts Setup
  if (fRejectPileupEvts) fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);

  //Output posting
  DataPosting();

}// end UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskMultspec_evbased::UserExec(Option_t *)
{
 
  //number of reconstructed V0
  int n_rec = 0;

  //get event from the input handler and cast it into the desired type of event
  AliAODEvent *lAODevent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!lAODevent) { 
    AliWarning("ERROR: event not available \n"); 
    DataPosting(); 
    return; 
  }

  // take track of number of analyzed events
  fHistos->FillTH1("hNevts", 0.5);
  
  // accepting events that pass internal selections in "AliEventCuts" 
  bool isEvtAccepted = fEventCuts.AcceptEvent(lAODevent);

  if (!isEvtAccepted) {
    DataPosting(); 
    return;
  }
  
  //get trigger information
  //fTriggerMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  // multiplicity Information 
  int lEvSelCode = -666;
  double lMultiplicity_v0 = -666.;
  double lMultiplicity_ch = -666.;
  AliMultiplicity* fMultiplicity =  (AliMultiplicity*)lAODevent->GetMultiplicity();
  AliMultSelection *MultSelection = (AliMultSelection*) lAODevent->FindListObject("MultSelection");
  if (!MultSelection) { 
    AliWarning("AliMultSelection object not found!"); 
    DataPosting(); 
    return; 
  } else {
    lMultiplicity_v0 = MultSelection->GetMultiplicityPercentile("V0M");
    // lMultiplicity_ch = MultSelection->GetNumberOfTracklets()->GetMultiplicityPercentile("RefMult08");
    // lMultiplicity_v0 = fMultiplicity->GetMultiplicityPercentile("V0M");
    lMultiplicity_ch = fMultiplicity->GetNumberOfTracklets();
    if(fpart<kxi_Mult) ffillV0->Mult = lMultiplicity_v0;
    else ffillCasc->Mult = lMultiplicity_v0;
    lEvSelCode = MultSelection->GetEvSelCode(); //==0 means event is good. Set by AliMultSelectionTask
  }
  if(lEvSelCode != 0) {   //skip everything if event selection code !=0
    DataPosting(); 
    return; 
  }

  //acquire best PV
  const AliVVertex *lBestPrimVtx = lAODevent->GetPrimaryVertex();
  double lBestPV[3] = {-666., -666., -666.};
  lBestPrimVtx->GetXYZ(lBestPV);
  float lfBestPV[3];
  for (int i=0; i<3; i++) lfBestPV[i] = (float)lBestPV[i];

  //acquire magnetic field
  double lMagField = -666;
  lMagField = lAODevent->GetMagneticField();
  
  //generator index
  int genindex = 0;
  enum indexGen {kPythia8_MultStrInj_pp5,kk0s_injector,kxi_injector};
  int indexRecP = 0;
  int indexRecN = 0;
  int indexRecB = 0;
  
  // start v0 part
  if(fpart<kxi_Mult){
    
    int nv0s = 0;
    nv0s = (int) lAODevent->GetNumberOfV0s();

    //start loop over V0 candidates
    for (int iV0 = 0; iV0 < nv0s; iV0++) {
        
      //get the iV0_th candidate in the event
      AliAODv0 *v0 = lAODevent->GetV0(iV0);
      if (!v0 || v0->GetOnFlyStatus()) continue;
      
      //get the 2D radius of the decay vertex
      ffillV0->V0Rad = v0->RadiusV0();
  
      //retrieve daughter AODTracks
      AliAODTrack *pTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(0);
      AliAODTrack *nTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(1);
      if (!pTrack || !nTrack) { AliWarning("ERROR: Could not retrieve one of the daughter tracks"); continue; }
     
      //pt and rapidity (cut the unwanted rapidities)
      ffillV0->Pt = v0->Pt();
      if( fpart==kk0s_Mult && TMath::Abs(v0->RapK0Short())>0.5 ) continue;
      if( (fpart==klam_Mult || fpart==kalam_Mult) && TMath::Abs(v0->RapLambda())>0.5 ) continue;

      //Daughters' Eta must be inside 0.8 to be true
      bool lV0_etaPos = TMath::Abs(pTrack->Eta())<0.8 ? kTRUE : kFALSE ;
      bool lV0_etaNeg = TMath::Abs(nTrack->Eta())<0.8 ? kTRUE : kFALSE ;
      if(!lV0_etaPos || !lV0_etaNeg) continue;
  
      //like-sign V0s (if any)
      if(((int) pTrack->GetSign()) == ((int) nTrack->GetSign()))  continue;
  
      //crossed raws
      double_t lCrosRawsPos = pTrack->GetTPCClusterInfo(2, 1);
      double_t lCrosRawsNeg = nTrack->GetTPCClusterInfo(2, 1);
      ffillV0->LeastCRaws = (int) TMath::Min(lCrosRawsPos, lCrosRawsNeg);
      
      //crossed raws / Findable clusters
      if(pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0) ffillV0->LeastCRawsOvF = 666;
      else{
        double_t lCrosRawsOvFPos = lCrosRawsPos/((double) (pTrack->GetTPCNclsF()));
        double_t lCrosRawsOvFNeg = lCrosRawsNeg/((double) (nTrack->GetTPCNclsF()));
        ffillV0->LeastCRawsOvF = TMath::Min(lCrosRawsOvFPos, lCrosRawsOvFNeg);
      }
      
      //dca info
      ffillV0->DcaPosToPV  = v0->DcaPosToPrimVertex();
      ffillV0->DcaNegToPV  = v0->DcaNegToPrimVertex();
      ffillV0->DcaV0Daught = v0->DcaV0Daughters();
  
      //cosPA
      ffillV0->V0CosPA = v0->CosPointingAngle(lBestPV);
  
      // invariant mass with different hypotheses
      if( fpart==kk0s_Mult ) {
          if ( v0->MassK0Short()<0.4||v0->MassK0Short()>0.6 ) continue;
          else {
            ffillV0->InvMass = v0->MassK0Short();
            ffillV0->CompMass = v0->MassLambda();
          }
      }
      else if( fpart==klam_Mult ){
          if ( v0->MassLambda()<1.08||v0->MassLambda()>1.16 ) continue;
          else {
            ffillV0->InvMass = v0->MassLambda();
            ffillV0->CompMass = v0->MassK0Short();
          }
      }
      else if( fpart==kalam_Mult ){
          if ( v0->MassAntiLambda()<1.08||v0->MassAntiLambda()>1.16 ) continue;
          else {
            ffillV0->InvMass = v0->MassAntiLambda();
            ffillV0->CompMass = v0->MassK0Short();
          }
      }
  
      //PID
      if( fpart==kk0s_Mult ) {
          ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
          ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
      }
      else if( fpart==klam_Mult ) {
          ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton);
          ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
      }
      else if( fpart==kalam_Mult ) {
          ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
          ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton);
      }
  
      //distance over total momentum
      ffillV0->DistOverTotP = v0->DecayLengthV0(lBestPV)/(v0->P()+1e-10);//avoid division by zero
  
      //TOF matching
      bool lNegTOFmatching = ( nTrack->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE );
      bool lPosTOFmatching = ( pTrack->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE );
      ffillV0->TOFmatch = (!lNegTOFmatching && !lPosTOFmatching) ? kFALSE : kTRUE;
      
      //ITS matching
      bool lNegITSmatching = ( !(nTrack->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE );
      bool lPosITSmatching = ( !(pTrack->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE );
      ffillV0->ITSmatch = (!lNegITSmatching && !lPosITSmatching) ? kFALSE : kTRUE;

      //ITSTOFtwo
      ffillV0->ITSTOFtwo = ( (lNegTOFmatching || lNegITSmatching) && (lPosTOFmatching || lPosITSmatching) ) ? kTRUE : kFALSE;
      
      //cut application
      if(ApplyCut(fpart)) {
        fHistos->FillTH1( fpart==kk0s_Mult ? "ImassK0sCheck" : "ImassLamCheck", ffillV0->InvMass );
        if( (ffillV0->InvMass>0.48 && ffillV0->InvMass<0.515) || (ffillV0->InvMass>1.11 && ffillV0->InvMass<1.12) ) n_rec++;
      }

    } // end of V0 loop

  } //end V0 part
  else {  // start cascade part

    int ncasc = 0;
    ncasc = lAODevent->GetNumberOfCascades();
    
    //start casc loop
    for (int i_casc = 0; i_casc < ncasc; i_casc++) {

      //get the iCasc_th candidate in the event
      AliAODcascade *casc = lAODevent->GetCascade(i_casc);
      if (!casc || casc->GetOnFlyStatus()) continue;
    
      //charge
      ffillCasc->charge = casc->ChargeXi()<0 ? kTRUE : kFALSE ;
    
      //cascade and V0 2D radii
      double lVtxCasc[3];
      lVtxCasc[0] = casc->DecayVertexXiX();
      lVtxCasc[1] = casc->DecayVertexXiY();
      lVtxCasc[2] = casc->DecayVertexXiZ();
      ffillCasc->CascRad = TMath::Sqrt(lVtxCasc[0]*lVtxCasc[0]+lVtxCasc[1]*lVtxCasc[1]);
      ffillCasc->V0Rad = casc->RadiusSecVtx();
    
      //get daughter tracks (positive, negative and bachelor)
      AliAODTrack *pTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDaughter(0));
      AliAODTrack *nTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDaughter(1));
      AliAODTrack *bTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDecayVertexXi()->GetDaughter(0));
      if (!pTrackCasc || !nTrackCasc || !bTrackCasc) { 
        AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...\n"); 
        continue; 
      }
      
      //transverse momentum
      ffillCasc->Pt = TMath::Sqrt(casc->Pt2Xi());
    
      //candidate's rapidity (mass hypothesis dependent)
      if(fpart==kxi_Mult && (TMath::Abs(casc->RapXi())>0.5)) continue;
      if(fpart==kom_Mult && (TMath::Abs(casc->RapOmega())>0.5)) continue;
    
    
      //daughters' etas must be inside 0.8 to be true
      bool letaPos = TMath::Abs(pTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      bool letaNeg = TMath::Abs(nTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      bool letaBac = TMath::Abs(bTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      if(!letaPos || !letaNeg || !letaBac) continue;
      
      //crossed raws
      double lCrosRawsPos = pTrackCasc->GetTPCClusterInfo(2, 1);
      double lCrosRawsNeg = nTrackCasc->GetTPCClusterInfo(2, 1);
      double lCrosRawsBac = bTrackCasc->GetTPCClusterInfo(2, 1);
      ffillCasc->LeastCRaws = (int) (lCrosRawsPos<lCrosRawsNeg ? std::min(lCrosRawsPos, lCrosRawsBac) : std::min(lCrosRawsNeg, lCrosRawsBac));
      
      //crossed raws / Findable clusters
      if(pTrackCasc->GetTPCNclsF()<=0 || nTrackCasc->GetTPCNclsF()<=0 || bTrackCasc->GetTPCNclsF()<=0) ffillCasc->LeastCRawsOvF = 0; //avoid division by zero
      else{    
        double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrackCasc->GetTPCNclsF()));
        double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrackCasc->GetTPCNclsF()));
        double lCrosRawsOvFBac = lCrosRawsBac / ((double)(bTrackCasc->GetTPCNclsF()));
        ffillCasc->LeastCRawsOvF = lCrosRawsOvFPos<lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos, lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg, lCrosRawsOvFBac);
      }
      
      //DCA info
      ffillCasc->DcaCascDaught = casc->DcaXiDaughters();
      ffillCasc->DcaBachToPV = casc->DcaBachToPrimVertex();
      ffillCasc->DcaPosToPV = casc->DcaPosToPrimVertex();
      ffillCasc->DcaNegToPV = casc->DcaNegToPrimVertex();
      ffillCasc->DcaV0Daught = casc->DcaV0Daughters();
      ffillCasc->DcaV0ToPV = casc->DcaV0ToPrimVertex();
    
      //cascade and V0 cosine of pointing angle
      ffillCasc->CascCosPA = casc->CosPointingAngleXi((const Double_t&) lBestPV[0], (const Double_t&) lBestPV[1], (const Double_t&) lBestPV[2]);
      ffillCasc->V0CosPA = casc->CosPointingAngle(lBestPV);
    
      //candidate's invariant mass
      if( fpart==kxi_Mult ) {
        if ( casc->MassXi()<1.26||casc->MassXi()>1.38 ) continue;
        else {
          ffillCasc->InvMass = casc->MassXi();
          ffillCasc->CompMass = casc->MassOmega();
        }
      }
      else if( fpart==kom_Mult ){
        if ( casc->MassOmega()<1.62||casc->MassOmega()>1.73 ) continue;
        else {
          ffillCasc->InvMass = casc->MassOmega();
          ffillCasc->CompMass = casc->MassXi();
        }
      }
    
      //PID
      if(casc->ChargeXi()<0){
        ffillCasc->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kProton);
        ffillCasc->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kPion);
      }
      else {
        ffillCasc->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kPion);
        ffillCasc->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kProton);
      }
      //check if bachelor's sign is always checked in compliance with Lambda/antilambda nature
      if(fpart==kxi_Mult) ffillCasc->NSigBac = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
      else ffillCasc->NSigBac = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);
      if( casc->ChargeXi()<0 && bTrackCasc->Charge()>0 ) printf("Wrong charge association\n");
    
      //distance over total momentum
      ffillCasc->DistOverTotP = (TMath::Sqrt(TMath::Power(lVtxCasc[0]-lBestPV[0], 2)+TMath::Power(lVtxCasc[1]-lBestPV[1], 2)+TMath::Power(lVtxCasc[2]-lBestPV[2], 2)))/(TMath::Sqrt(casc->Ptot2Xi())+1e-10);
    
      //TOF matching
      bool lNegTOFmatching = nTrackCasc->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE ;
      bool lPosTOFmatching = pTrackCasc->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE ;
      bool lBacTOFmatching = bTrackCasc->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE ;
      ffillCasc->TOFmatch = (!lNegTOFmatching && !lPosTOFmatching && !lBacTOFmatching) ? kFALSE : kTRUE;
      
      //ITS matching
      bool lNegITSmatching = !(nTrackCasc->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE ;
      bool lPosITSmatching = !(pTrackCasc->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE ;
      bool lBacITSmatching = !(bTrackCasc->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE ;
      ffillCasc->ITSmatch = (!lNegITSmatching && !lPosITSmatching && !lBacITSmatching) ? kFALSE : kTRUE;
      
      //ITSTOFtwo
      bool negcond = (lNegTOFmatching || lNegITSmatching) ? kTRUE : kFALSE;
      bool poscond = (lPosTOFmatching || lPosITSmatching) ? kTRUE : kFALSE;
      bool baccond = (lBacTOFmatching || lBacITSmatching) ? kTRUE : kFALSE;
      ffillCasc->ITSTOFtwo = ( (negcond && poscond) || (negcond && baccond) || (poscond && baccond) ) ? kTRUE : kFALSE;

      //V0 daughter mass
      double_t imassla = 0.;
      (casc->ChargeXi()<0 ? imassla=casc->MassLambda() : imassla=casc->MassAntiLambda());
      if (TMath::Abs(imassla-1.115683)</*0.008*/0.0035) ffillCasc->GoodInvMassLam = kTRUE;
      else ffillCasc->GoodInvMassLam = kFALSE;
    
      //calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
      ffillCasc->BacBarCosPA = casc->BachBaryonCosPA();
    
      if(ApplyCut(fpart)) {
        fHistos->FillTH1( fpart==kxi_Mult ? "ImassXiCheck" : "ImassOmCheck", ffillCasc->InvMass );
        if( (ffillCasc->InvMass>1.315 && ffillCasc->InvMass<1.33) || (ffillCasc->InvMass>1.665 && ffillCasc->InvMass<1.68) ) n_rec++;
      }

    } //end casc loop
    
  }
  
  fHistos->FillTH3("CorrMultNstrange",n_rec,lMultiplicity_v0,lMultiplicity_ch);
  fHistos->FillTH1("RecoSpectrum",n_rec);

  DataPosting();

}

bool AliAnalysisTaskMultspec_evbased::ApplyCut(int parttype){

  bool passcut = kFALSE;

  if(parttype==kk0s_Mult){
     if(((ffillV0->DcaV0Daught)<0.3) && ((ffillV0->DcaPosToPV)>0.13) && ((ffillV0->DcaNegToPV)>0.13) && ((ffillV0->V0CosPA)>0.999) && ((ffillV0->DistOverTotP)<(15./0.498115)) && ((ffillV0->V0Rad)>0.7) && ((ffillV0->LeastCRawsOvF)>=0.9) && ((ffillV0->NSigPos)>-1) && ((ffillV0->NSigPos)<1) && ((ffillV0->NSigNeg)>-1) && ((ffillV0->NSigNeg)<1) && ((ffillV0->ITSmatch)==kTRUE || (ffillV0->TOFmatch)==kTRUE)) passcut = kTRUE;
   }
  else if(parttype==klam_Mult || parttype==kalam_Mult){
     if(((ffillV0->DcaV0Daught)<0.5) && ((ffillV0->DcaPosToPV)>0.13) && ((ffillV0->DcaNegToPV)>0.13) && ((ffillV0->V0CosPA)>0.999) && ((ffillV0->DistOverTotP)<(20./1.11603)) && ((ffillV0->V0Rad)>1.05) && ((ffillV0->LeastCRawsOvF)>=0.9) && ((ffillV0->NSigPos)>-1) && ((ffillV0->NSigPos)<1) && ((ffillV0->NSigNeg)>-1) && ((ffillV0->NSigNeg)<1) && ((ffillV0->ITSmatch)==kTRUE || (ffillV0->TOFmatch)==kTRUE)) passcut = kTRUE;
  }
  else if(parttype==kxi_Mult){
     if(((ffillCasc->CascRad)>1) && ((ffillCasc->V0Rad)>2.5) && ((ffillCasc->DcaNegToPV)>0.45) && ((ffillCasc->DcaPosToPV)>0.45) && ((ffillCasc->DcaV0Daught)<0.5) && ((ffillCasc->DcaBachToPV)>0.18) && ((ffillCasc->DcaV0ToPV)>0.135) && ((ffillCasc->DcaCascDaught)<0.065) && ((ffillCasc->V0CosPA)>0.99) && ((ffillCasc->CascCosPA)>0.999) && ((ffillCasc->NSigPos)>-1) && ((ffillCasc->NSigPos)<1) && ((ffillCasc->NSigNeg)>-1) && ((ffillCasc->NSigNeg)<1) && ((ffillCasc->NSigBac)>-1) && ((ffillCasc->NSigBac)<1) && ((ffillCasc->DistOverTotP)<(10.5/1.32171)) && ((ffillCasc->LeastCRaws)>77)  && ((ffillCasc->ITSmatch)==kTRUE || (ffillCasc->TOFmatch)==kTRUE)) passcut = kTRUE;

  }
  else if(parttype==kom_Mult){
     if(((ffillCasc->CascRad)>0.9) && ((ffillCasc->V0Rad)>2.9) && ((ffillCasc->DcaNegToPV)>0.38) && ((ffillCasc->DcaPosToPV)>0.38) && ((ffillCasc->DcaV0Daught)<0.5) && ((ffillCasc->DcaBachToPV)>0.09) && ((ffillCasc->DcaV0ToPV)>0.122) && ((ffillCasc->DcaCascDaught)<0.4) && ((ffillCasc->V0CosPA)>0.99) && ((ffillCasc->CascCosPA)>0.999) && ((ffillCasc->NSigPos)>-1) && ((ffillCasc->NSigPos)<1) && ((ffillCasc->NSigNeg)>-1) && ((ffillCasc->NSigNeg)<1) && ((ffillCasc->NSigBac)>-1) && ((ffillCasc->NSigBac)<1) && ((ffillCasc->DistOverTotP)<(5.5/1.67245)) && ((ffillCasc->LeastCRaws)>81)  && ((ffillCasc->ITSmatch)==kTRUE || (ffillCasc->TOFmatch)==kTRUE)) passcut = kTRUE;
  }

  return passcut;

}

//________________________________________________________________________
void AliAnalysisTaskMultspec_evbased::Terminate(Option_t *)
{

}

//________________________________________________________________________
void AliAnalysisTaskMultspec_evbased::DataPosting() {

  PostData(1, fHistos->GetListOfHistograms());
  // PostData(2, fTree);
  
}
