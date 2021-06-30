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

#include "AliAnalysisTaskV0multspec.h"


ClassImp(AliAnalysisTaskV0multspec)

AliAnalysisTaskV0multspec::AliAnalysisTaskV0multspec() : AliAnalysisTaskSE(),
fHistos_misc(nullptr), 
fTreeV0(nullptr), 
fTreeCascade(nullptr),
fAddCasc(kFALSE),
fpart(666),
fisMC(kTRUE),
fPIDResponse(0), 
ffillV0(nullptr),
ffillCasc(nullptr)
{
  //default constructor
}

AliAnalysisTaskV0multspec::AliAnalysisTaskV0multspec(const char *name, int particle, bool AddCasc, TString lExtraOptions, bool ismc) : AliAnalysisTaskSE(name),
fHistos_misc(nullptr), 
fTreeV0(nullptr), 
fTreeCascade(nullptr),
fAddCasc(AddCasc),
fpart(particle),
fisMC(ismc),
fPIDResponse(0), 
ffillV0(nullptr),
ffillCasc(nullptr)
{

  //SetAddCasc(AddCasc);
  SetIsMC(ismc);
  SetParticle(particle);
    
  //Standard output
  DefineOutput(1, TList::Class()); // Miscellaneous Histograms
  DefineOutput(2, TTree::Class()); // Output Tree, V0s
  if(fAddCasc) DefineOutput(3, TTree::Class()); // Output Tree, Cascades
  
}

AliAnalysisTaskV0multspec::~AliAnalysisTaskV0multspec()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    if (fHistos_misc) {
        delete fHistos_misc;
        fHistos_misc = 0x0;
    }
    if (fTreeV0) {
        delete fTreeV0;
        fTreeV0 = 0x0;
    }
    if (fTreeCascade) {
        delete fTreeCascade;
        fTreeCascade = 0x0;
    }
 
 
}

//________________________________________________________________________
void AliAnalysisTaskV0multspec::UserCreateOutputObjects()
{

  //miscellaneous histograms
  fHistos_misc = new THistManager("fHistos_misc");
  fHistos_misc->CreateTH1("hNevts", "", 1, 0, 1);  //dummy #events histo
  fHistos_misc->CreateTH1("hV0multspec_nocut", "", 1000, 0, 1000, "s");  //dummy V0 multiplicity spectrum for checking
  fHistos_misc->CreateTH2("hV0multVSEvmult_nocut", "", 1000, 0, 1000, 100,0,100, "s");  //dummy V0 multiplicity VS Event multiplicity spectrum for checking
  fHistos_misc->CreateTH1("hCascmultspec_nocut", "", 100, 0, 100, "s");  //dummy Cascades multiplicity spectrum for checking
  fHistos_misc->CreateTH2("hCascmultVSEvmult_nocut", "", 100, 0, 100, 100,0,100, "s");  //dummy Cascades multiplicity VS Event multiplicity spectrum for checking
  fHistos_misc->CreateTH1("hEvSelCode", "", 1000, 0, 1000, "s");  //event selection code, for monitoring
  //MC-related TH
  if(fisMC) {
      fHistos_misc->CreateTH2("hR0gen", "", 100, 0, 100, 100, 0, 100, "s");  //histogram Gen_distrib for zero-Rec events for given part VS mult
      fHistos_misc->CreateTH2("hV0multspec_Gen", "", 100, 0, 100, 100, 0, 100, "s");  //histogram Gen_distrib given part VS mult
  }

  //V0 tree
  ffillV0 = new V0filler;
  fTreeV0 = new TTree("fTreeV0","V0s");
  fTreeV0->Branch("V0filler", ffillV0);
  
  //Cascade tree
  if(fAddCasc){
      ffillCasc = new Cascfiller;
      fTreeCascade = new TTree("fTreeCascade","Cascades");
      fTreeCascade->Branch("Cascfiller",ffillCasc);
  }
  
  // PID Setup
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  inputHandler->SetNeedField();

  //Output posting
  DataPosting();

}// end UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskV0multspec::UserExec(Option_t *)
{

  //get event from the input handler and cast it into the desired type of event
  AliAODEvent *lAODevent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!lAODevent) { 
    AliWarning("ERROR: event not available \n"); 
    DataPosting(); 
    return; 
  }

  //get MC event and fill histo
  AliMCEvent  *lMCev  = 0x0;
  if(fisMC){
    lMCev = MCEvent();
    if (!lMCev) {
      Printf("ERROR: Could not retrieve MC event in file %s\n",fInputHandler->GetTree()->GetCurrentFile()->GetName());
      DataPosting(); 
      return;
    }
  }
  //get MC header and MC array
  AliAODMCHeader* header = 0x0;
  TClonesArray* MCTrackArray = 0x0;
  if(fisMC) {
      header = static_cast<AliAODMCHeader*>(lAODevent->FindListObject(AliAODMCHeader::StdBranchName()));
      if (!header) {
        AliWarning("No header found.");
        DataPosting(); 
        return;
      } 
      MCTrackArray = dynamic_cast<TClonesArray*>(lAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (MCTrackArray == NULL){
        AliWarning("No MC track array found.");
        DataPosting(); 
        return;
      }  
  }
  
  //dumb histo for checking
  fHistos_misc->FillTH1("hNevts", 0.5);

  //get trigger information
  //fTriggerMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  // Multiplicity Information 
  int lEvSelCode = -666;
  double lMultiplicity = -666.;
  AliMultSelection *MultSelection = (AliMultSelection*) lAODevent->FindListObject("MultSelection");
  if (!MultSelection) { 
    AliWarning("AliMultSelection object not found!"); 
    DataPosting(); 
    return; 
  } else {
    lMultiplicity = MultSelection->GetMultiplicityPercentile("V0M");
    ffillV0->Mult = lMultiplicity;
    if(fAddCasc) ffillCasc->Mult = lMultiplicity;
    lEvSelCode = MultSelection->GetEvSelCode(); //==0 means event is good. Set by AliMultSelectionTask
    fHistos_misc->FillTH1("hEvSelCode", lEvSelCode);
  }

  //skip everything if event selection code !=0
  if(lEvSelCode != 0) { 
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

  // start v0 part
  int nv0s = 0;
  nv0s = (int) lAODevent->GetNumberOfV0s();
  fHistos_misc->FillTH1("hV0multspec_nocut",nv0s);
  fHistos_misc->FillTH2("hV0multVSEvmult_nocut",nv0s,lMultiplicity);
  bool newV0 = kTRUE;
  
  //loop over generated MC 
  int initgens = 253;
  ffillV0->nGen  = initgens;
  if(fisMC){
    int ngen = 0;
    for (int i_MCtrk = 0;  i_MCtrk < lMCev->GetNumberOfTracks(); i_MCtrk++){
      AliVParticle *lPart = (AliAODMCParticle*) lMCev->GetTrack(i_MCtrk);
      if(!lPart || lPart->Y()<-0.5 || lPart->Y()>0.5 || !lPart->IsPhysicalPrimary()) continue;
      if( fpart==kk0s && lPart->PdgCode()==310   ) ngen++;
      if( fpart==klam && lPart->PdgCode()==3122  ) ngen++;
      if( fpart==kalam && lPart->PdgCode()==-3122 ) ngen++;
    }
    ffillV0->nGen  = ngen;
    //fill special gen histo if nv0s is zero
    if(nv0s==0) {
        fHistos_misc->FillTH2("hR0gen", ngen, lMultiplicity);
    }
    //fill generated histos 
    fHistos_misc->FillTH2("hV0multspec_Gen", ngen, lMultiplicity);
  }

  
  for (int iV0 = 0; iV0 < nv0s; iV0++) {
      
    //set a fV0_IsNewEvt flag if this is the first V0 of the loop
    ffillV0->IsNewEvt = newV0;

    //get the iV0_th candidate in the event
    AliAODv0 *v0 = lAODevent->GetV0(iV0);
    if (!v0 || v0->GetOnFlyStatus()) continue;

    //get the 2D radius of the decay vertex
    ffillV0->V0Rad = v0->RadiusV0();

    //pt and rapidity (cut the unwanted rapidities)
    ffillV0->Pt = v0->Pt();
    if( fpart==kk0s && TMath::Abs(v0->RapK0Short())>0.5 ) continue;
    if( (fpart==klam || fpart==kalam) && TMath::Abs(v0->RapLambda())>0.5 ) continue;

    //retrieve daughter AODTracks
    AliAODTrack *pTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(0);
    AliAODTrack *nTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(1);
    if (!pTrack || !nTrack) { AliWarning("ERROR: Could not retrieve one of the daughter tracks"); continue; }

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
    if( fpart==kk0s ) {
        if ( v0->MassK0Short()<0.4||v0->MassK0Short()>0.6 ) continue;
        else ffillV0->InvMass = v0->MassK0Short();
    }
    else if( fpart==klam ){
        if ( v0->MassLambda()<1.08||v0->MassLambda()>1.16 ) continue;
        else ffillV0->InvMass = v0->MassLambda();
    }
    else if( fpart==kalam ){
        if ( v0->MassAntiLambda()<1.08||v0->MassAntiLambda()>1.16 ) continue;
        else ffillV0->InvMass = v0->MassAntiLambda();
    }

    //PID
    if( fpart==kk0s ) {
        ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
        ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
    }
    else if( fpart==klam ) {
        ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton);
        ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
    }
    else if( fpart==kalam ) {
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
    
    //tree filling 
    fTreeV0->Fill();
    newV0 = kFALSE;


  } // end of V0 loop

  // start cascade part
  int ncasc = 0;
  ncasc = lAODevent->GetNumberOfCascades();
  fHistos_misc->FillTH1("hCascmultspec_nocut",ncasc);
  fHistos_misc->FillTH2("hCascmultVSEvmult_nocut",ncasc,lMultiplicity);
  bool newCasc = kTRUE;

  if(fAddCasc) {
  
    for (int i_casc = 0; i_casc < ncasc; i_casc++) {
      
      //set a fCasc_IsNewEvt flag if this is the first V0 of the loop
      ffillCasc->IsNewEvt = newCasc;
  
      //get the iCasc_th candidate in the event
      AliAODcascade *casc = lAODevent->GetCascade(i_casc);
      if (!casc) continue;
  
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
  
      //daughters' etas
      bool letaPos = TMath::Abs(pTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      bool letaNeg = TMath::Abs(nTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      bool letaBac = TMath::Abs(bTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      if(!letaPos || !letaNeg || !letaBac) continue;
  
      //PID
      ffillCasc->NSigPosProton = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kProton);
      ffillCasc->NSigPosPion = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kPion);
      ffillCasc->NSigNegProton = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kProton);
      ffillCasc->NSigNegPion = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kPion);
      ffillCasc->NSigBacPion = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
      ffillCasc->NSigBacKaon = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);
  
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
      
      //candidate's rapidity (mass hypothesis dependent)
      ffillCasc->yXi = TMath::Abs(casc->RapXi())<0.5 ? kTRUE : kFALSE ;
      ffillCasc->yOm = TMath::Abs(casc->RapOmega())<0.5 ? kTRUE : kFALSE ;
  
      //charge
      ffillCasc->charge = casc->ChargeXi()>0 ? kTRUE : kFALSE ;
  
      //V0 daughter mass (later to be checked against nominal)
      if (casc->ChargeXi()<0) ffillCasc->InvMassLam = casc->MassLambda();
      else ffillCasc->InvMassLam = casc->MassAntiLambda();
  
      //transverse momentum
      ffillCasc->Pt = TMath::Sqrt(casc->Pt2Xi());
  
      //distance over total momentum
      ffillCasc->DistOverTotP = (TMath::Sqrt(TMath::Power(lVtxCasc[0]-lBestPV[0], 2)+TMath::Power(lVtxCasc[1]-lBestPV[1], 2)+TMath::Power(lVtxCasc[2]-lBestPV[2], 2)))/(TMath::Sqrt(casc->Ptot2Xi())+1e-10);
  
      //candidate's invariant mass
      ffillCasc->InvMassXi = casc->MassXi();
      ffillCasc->InvMassOm = casc->MassOmega();
  
      //calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
      ffillCasc->BacBarCosPA = casc->BachBaryonCosPA();
  
      //tree filling 
      fTreeCascade->Fill();
      newCasc = kFALSE;
      
    } //end casc loop
  }
  
  DataPosting();

}

//________________________________________________________________________
void AliAnalysisTaskV0multspec::Terminate(Option_t *)
{

}

//________________________________________________________________________
void AliAnalysisTaskV0multspec::DataPosting() {

  PostData(1, fHistos_misc->GetListOfHistograms());
  PostData(2, fTreeV0);
  if(fAddCasc) PostData(3, fTreeCascade);
  
}
