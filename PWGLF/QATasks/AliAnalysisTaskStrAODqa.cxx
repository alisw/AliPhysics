class TTree;
class TParticle;
class TVector3;
class AliAODVertex;
class AliAODv0;
class AliAODcascade;

#include <Riostream.h>
#include "TH3.h"
#include "TCanvas.h"
#include "THistManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliPID.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskStrAODqa.h"

ClassImp(AliAnalysisTaskStrAODqa)

AliAnalysisTaskStrAODqa::AliAnalysisTaskStrAODqa()
: AliAnalysisTaskSE(),
//outputs
fHistos_eve(0),
fHistos_V0(0),
fHistos_Casc(0),
//objects from the manager
fPIDResponse(0),
//variables for MC 
  fMCEvent(0),
  fReadMCTruth(0),
//variables for V0 cuts
fV0_DcaV0Daught(0),
fV0_DcaPosToPV(0),
fV0_DcaNegToPV(0),
fV0_V0CosPA(0),
fV0_V0Rad(0),
fV0_Pt(0),
fV0_yK0S(0),
fV0_yLam(0),
fV0_etaPos(0),
fV0_etaNeg(0),
fV0_InvMassK0s(0),
fV0_InvMassLam(0),
fV0_InvMassALam(0),
fV0_LeastCRaws(0),
fV0_LeastCRawsOvF(0),
fV0_NSigPosProton(0),
fV0_NSigPosPion(0),
fV0_NSigNegProton(0),
fV0_NSigNegPion(0),
fV0_DistOverTotP(0),
//variables for Cascade analysis
fCasc_isNotTPCRefit(0),
fCasc_DcaCascDaught(0),
fCasc_CascCosPA(0),
fCasc_CascRad(0),
fCasc_etaPos(0),
fCasc_etaNeg(0),
fCasc_etaBac(0),
fCasc_NSigPosProton(0),
fCasc_NSigPosPion(0),
fCasc_NSigNegProton(0),
fCasc_NSigNegPion(0),
fCasc_NSigBacPion(0),
fCasc_NSigBacKaon(0),
fCasc_LeastCRaws(0),
fCasc_LeastCRawsOvF(0),
fCasc_DcaV0Daught(0),
fCasc_V0CosPA(0),
fCasc_V0CosPAToXi(0),
fCasc_DcaV0ToPV(0),
fCasc_DcaBachToPV(0),
fCasc_yXi(0),
fCasc_yOm(0),
fCasc_charge(0),
fCasc_Pt(0),
fCasc_Ptot(0),
fCasc_DistOverTotP(0),
fCasc_V0DistOverTotP(0),
fCasc_CascCtauXi(0),
fCasc_CascCtauOmega(0),
fCasc_V0Ctau(0),
fCasc_InvMassXi(0),
fCasc_InvMassOm(0),
fCasc_V0Rad(0),
fCasc_DcaPosToPV(0),
  fCasc_DcaNegToPV(0),
  fCasc_InvMassLambda(0)

{

}

AliAnalysisTaskStrAODqa::AliAnalysisTaskStrAODqa(const char *name, TString lExtraOptions)
: AliAnalysisTaskSE(name),
//outputs
fHistos_eve(0),
fHistos_V0(0),
fHistos_Casc(0),
//objects from the manager
fPIDResponse(0),
//variables for MC 
fMCEvent(0),
fReadMCTruth(0),
//variables for V0 cuts
fV0_DcaV0Daught(0),
fV0_DcaPosToPV(0),
fV0_DcaNegToPV(0),
fV0_V0CosPA(0),
fV0_V0Rad(0),
fV0_Pt(0),
fV0_yK0S(0),
fV0_yLam(0),
fV0_etaPos(0),
fV0_etaNeg(0),
fV0_InvMassK0s(0),
fV0_InvMassLam(0),
fV0_InvMassALam(0),
fV0_LeastCRaws(0),
fV0_LeastCRawsOvF(0),
fV0_NSigPosProton(0),
fV0_NSigPosPion(0),
fV0_NSigNegProton(0),
fV0_NSigNegPion(0),
fV0_DistOverTotP(0),
//variables for Cascade analysis
fCasc_isNotTPCRefit(0),
fCasc_DcaCascDaught(0),
fCasc_CascCosPA(0),
fCasc_CascRad(0),
fCasc_etaPos(0),
fCasc_etaNeg(0),
fCasc_etaBac(0),
fCasc_NSigPosProton(0),
fCasc_NSigPosPion(0),
fCasc_NSigNegProton(0),
fCasc_NSigNegPion(0),
fCasc_NSigBacPion(0),
fCasc_NSigBacKaon(0),
fCasc_LeastCRaws(0),
fCasc_LeastCRawsOvF(0),
fCasc_DcaV0Daught(0),
fCasc_V0CosPA(0),
fCasc_V0CosPAToXi(0),
fCasc_DcaV0ToPV(0),
fCasc_DcaBachToPV(0),
fCasc_yXi(0),
fCasc_yOm(0),
fCasc_charge(0),
fCasc_Pt(0),
fCasc_Ptot(0),
fCasc_DistOverTotP(0),
fCasc_V0DistOverTotP(0),
fCasc_CascCtauXi(0),
fCasc_CascCtauOmega(0),
fCasc_V0Ctau(0),
fCasc_InvMassXi(0),
fCasc_InvMassOm(0),
fCasc_V0Rad(0),
fCasc_DcaPosToPV(0),
fCasc_DcaNegToPV(0),
fCasc_InvMassLambda(0)
{

    //Standard output
    DefineOutput(1, TList::Class()); // Event Histograms
    DefineOutput(2, TList::Class()); // V0 Histograms
    DefineOutput(3, TList::Class()); // Cascades Histograms

}


AliAnalysisTaskStrAODqa::~AliAnalysisTaskStrAODqa()
{

}

//________________________________________________________________________
void AliAnalysisTaskStrAODqa::UserCreateOutputObjects()
{

    //histograms for event variables
    fHistos_eve = new THistManager("histos_eve");
    //fHistos_eve->CreateTH1("hcent", "", 100, 0, 100, "s");  //storing #events in bins of centrality
    fHistos_eve->CreateTH1("henum", "", 1, 0, 1);  //storing total #events
    fHistos_eve->CreateTH3("GeneratedParticles", "", 7, 0, 7, 100, 0, 10, 200, -10, 10);  //storing generated particles
    //    for (int iP=1; iP<=kNParticles; iP++) ((TH2*)fHistos_eve->FindObject("GeneratedParticles"))->GetXaxis()->SetBinLabel(iP, kParticleNames[iP-1]);    

    fHistos_V0 = new THistManager("histos_V0");
    fHistos_V0->CreateTH1("CosPA", "", 100, 0.9, 1.);
    fHistos_V0->CreateTH1("Radius", "", 100, 0., 10.);
    fHistos_V0->CreateTH1("V0DCANegToPV",  "", 50, 0., 1.);
    fHistos_V0->CreateTH1("V0DCAPosToPV", "", 50, 0., 1.);
    fHistos_V0->CreateTH1("V0DCAV0Daughters",  "", 50, 0., 2.);

    fHistos_V0->CreateTH2("ResponsePionFromLambda", "", 80, 0., 10., 80, -4., 4.);
    fHistos_V0->CreateTH2("ResponseProtonFromLambda", "", 80, 0., 10., 80, -4., 4.);

    fHistos_V0->CreateTH2("ImassK0S", "", 100, 0., 10., 200, 0.4, 0.6);
    fHistos_V0->CreateTH2("ImassLam", "", 100, 0., 10., 200, 1.07, 1.17);
    fHistos_V0->CreateTH2("ImassALam", "", 100, 0., 10., 200, 1.07, 1.17);
    fHistos_V0->CreateTH2("ImassK0STrue", "", 100, 0., 10., 200, 0.4, 0.6);
    fHistos_V0->CreateTH2("ImassLamTrue", "", 100, 0., 10., 200, 1.07, 1.17);
    fHistos_V0->CreateTH2("ImassALamTrue", "", 100, 0., 10., 200, 1.07, 1.17);

    fHistos_Casc = new THistManager("histos_Casc");
    fHistos_Casc->CreateTH2("XiProgSelections","", 30, 0.5, 30.5, 2 , -2, 2);
    fHistos_Casc->CreateTH2("OmegaProgSelections","", 30, 0.5, 30.5, 2, -2, 2);
    fHistos_Casc->CreateTH2("CascCosPA","", 200,0.90,1.0, 2, -2, 2);
    fHistos_Casc->CreateTH2("V0CosPA","", 100,0.9,1.0, 2, -2, 2);
    fHistos_Casc->CreateTH2("V0CosPAToXi","", 100,0.9,1.0, 2, -2, 2);
    fHistos_Casc->CreateTH2("CascRadius","", 100,0.0,10.0, 2, -2, 2);
    fHistos_Casc->CreateTH2("V0Radius","", 100,0.0,10.0, 2, -2, 2);
    fHistos_Casc->CreateTH2("CascyXi","", 200,-2.0,2.0, 2, -2, 2);
    fHistos_Casc->CreateTH2("CascyOmega","", 200,-2.0,2.0, 2, -2, 2);
    fHistos_Casc->CreateTH2("CascCtauXi","", 100,0,100, 2, -2, 2);
    fHistos_Casc->CreateTH2("CascCtauOmega","", 100,0,100, 2, -2, 2);
    fHistos_Casc->CreateTH2("V0Ctau","", 100,0,100, 2, -2, 2);
    fHistos_Casc->CreateTH2("CascPt","", 100, 0, 25, 2, -2, 2);
    fHistos_Casc->CreateTH2("DcaV0Daughters","", 100,0,2, 2, -2, 2);
    fHistos_Casc->CreateTH2("DcaCascDaughters","", 100,0,2, 2, -2, 2);
    fHistos_Casc->CreateTH2("DcaV0ToPV","", 40,0,0.2, 2, -2, 2);
    fHistos_Casc->CreateTH2("DcaBachToPV","", 40,0,0.2, 2, -2, 2);
    fHistos_Casc->CreateTH2("DcaPosToPV","", 40,0,0.2, 2, -2, 2);
    fHistos_Casc->CreateTH2("DcaNegToPV","", 40,0,0.2, 2, -2, 2);
    fHistos_Casc->CreateTH2("InvMassLambdaDaughter", "",  100,1.1, 1.13, 2, -2, 2); 
    fHistos_Casc->CreateTH2("CascCfrPt","", 100, 0, 25, 100, 0, 25);
    fHistos_Casc->CreateTH2("CascCfrP","", 100, 0, 25, 100, 0, 25);

    fHistos_Casc->CreateTH2("ImassXiPlu","",100, 0, 10,80,1.28,1.36);
    fHistos_Casc->CreateTH2("ImassXiMin","",100, 0, 10,80,1.28,1.36);
    fHistos_Casc->CreateTH2("ImassOmPlu","",100, 0, 10,80,1.63,1.71);
    fHistos_Casc->CreateTH2("ImassOmMin","",100, 0, 10,80,1.63,1.71);
    fHistos_Casc->CreateTH2("ImassXiPluTrue","",100, 0, 10,80,1.28,1.36);
    fHistos_Casc->CreateTH2("ImassXiMinTrue","",100, 0, 10,80,1.28,1.36);
    fHistos_Casc->CreateTH2("ImassOmPluTrue","",100, 0, 10,80,1.63,1.71);
    fHistos_Casc->CreateTH2("ImassOmMinTrue","",100, 0, 10,80,1.63,1.71);

    // PID Setup
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();

    //Output
    PostData(1, fHistos_eve->GetListOfHistograms()    );
    PostData(2, fHistos_V0->GetListOfHistograms()    );
    PostData(3, fHistos_Casc->GetListOfHistograms()    );

}// end UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskStrAODqa::UserExec(Option_t *)
{

    //get event from the inpot handler and cast it into the desired type of event
    AliAODEvent *lAODevent = 0x0;
    lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
    if (!lAODevent) {
      AliWarning("ERROR: AODevent not available \n");
      PostData(1, fHistos_eve->GetListOfHistograms()    );
      PostData(2, fHistos_V0->GetListOfHistograms()    );
      PostData(3, fHistos_Casc->GetListOfHistograms()    );
      return;
    }

    // dumb histo for checking
    fHistos_eve->FillTH1("henum", 0.5);

    // Multiplicity Information
/*    AliMultSelection *MultSelection = (AliMultSelection*) lAODevent -> FindListObject("MultSelection");
    if( !MultSelection ) {
      AliWarning("AliMultSelection object not found!");
      PostData(1, fHistos_eve->GetListOfHistograms()    );
      PostData(2, fHistos_V0->GetListOfHistograms()    );
      PostData(3, fHistos_Casc->GetListOfHistograms()    );
     return;
   }
    else if(MultSelection->GetEvSelCode() != 0){
        PostData(1, fHistos_eve->GetListOfHistograms()    );
        PostData(2, fHistos_V0->GetListOfHistograms()    );
        PostData(3, fHistos_Casc->GetListOfHistograms()    );
        return;
    }
*/
    //acquire best PV
    AliAODVertex *lBestAODPrimVtx = lAODevent->GetPrimaryVertex();
    double lBestPV[3]  = {-666., -666., -666.};
    lBestAODPrimVtx->GetXYZ( lBestPV );


    //MC generated part 

    TClonesArray* AODMCTrackArraybis =0x0;
    if(fReadMCTruth){
      fMCEvent= MCEvent();
      if (fMCEvent){
	AODMCTrackArraybis = dynamic_cast<TClonesArray*>(lAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (AODMCTrackArraybis == NULL){
	  return;
	  Printf("ERROR: stack not available");
	}  
	for(Long_t i = 0; i < AODMCTrackArraybis->GetEntriesFast(); i++) {
	  AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArraybis->At(i));
	  if (!particle) continue;
	  if (!(particle->IsPhysicalPrimary()))continue; //we are mainly interested in the primaries because we want to see the effect of the injection of strange particles on the pT spectrum

	  if (particle->GetPdgCode()==  310)     fHistos_eve->FillTH3("GeneratedParticles", 0.5, particle->Pt(),particle->Y()); //K0s
	  if (particle->GetPdgCode()== 3122)     fHistos_eve->FillTH3("GeneratedParticles", 1.5, particle->Pt(),particle->Y()); //Lambda 
	  if (particle->GetPdgCode()==-3122)     fHistos_eve->FillTH3("GeneratedParticles", 2.5, particle->Pt(),particle->Y()); //AntiLambda
	  if (particle->GetPdgCode()== 3312)     fHistos_eve->FillTH3("GeneratedParticles", 3.5, particle->Pt(),particle->Y()); //Xi- 
	  if (particle->GetPdgCode()==-3312)     fHistos_eve->FillTH3("GeneratedParticles", 4.5, particle->Pt(),particle->Y()); //Xi+ 
	  if (particle->GetPdgCode()== 3334)     fHistos_eve->FillTH3("GeneratedParticles", 5.5, particle->Pt(),particle->Y()); //Omega- 
	  if (particle->GetPdgCode()==-3334)     fHistos_eve->FillTH3("GeneratedParticles", 6.5, particle->Pt(),particle->Y()); //Omega+ 
	}
      }
    }

    //MC variables useful for both V0s and Cascades
    AliAODMCParticle* particlePos;
    AliAODMCParticle* particleNeg;
    AliAODMCParticle* particleBach;
    Int_t PdgPos=0;
    Int_t PdgNeg=0;
    Int_t PdgBach=0;

    Int_t labelMotherPos=0;
    Int_t labelMotherNeg=0;
    Int_t labelGMotherPos=0;
    Int_t labelGMotherNeg=0;
    Int_t labelMotherBach=0;
    AliAODMCParticle* MotherPos;
    AliAODMCParticle* MotherNeg;
    AliAODMCParticle* GMotherPos;
    AliAODMCParticle* GMotherNeg;
    AliAODMCParticle* MotherBach;
    Int_t PdgMotherPos=0;
    Int_t PdgMotherNeg=0;
    Int_t PdgMotherLambda=0;
    Int_t PdgGMotherPos=0;
    Int_t PdgGMotherNeg=0;
    Int_t PdgMotherBach=0;
    Int_t V0PDGCode=0;

    // start of v0 part
    // ----------------
    int nv0s = 0;
    nv0s = lAODevent->GetNumberOfV0s();

    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {

        AliAODv0 *v0 = lAODevent->GetV0(iV0);
        if (!v0 || v0->GetOnFlyStatus()) continue;

        fV0_Pt = v0->Pt();
        fV0_yK0S = v0->RapK0Short();
        fV0_yLam  = v0->RapLambda();
        fV0_V0Rad = v0->RadiusV0();
        fV0_InvMassK0s = v0->MassK0Short();
        fV0_InvMassLam = v0->MassLambda();
        fV0_InvMassALam = v0->MassAntiLambda();

        //retrieve daughter AODTracks
        AliAODTrack *pTrack = (AliAODTrack*) lAODevent->GetTrack(0);
        AliAODTrack *nTrack = (AliAODTrack*) lAODevent->GetTrack(1);
        if (!pTrack || !nTrack) { AliWarning("ERROR: Could not retrieve one of the daughter tracks\n"); continue; }



	//-------------------------------------------------------                                                                                           
  	//---------MC information--------------------------------                                                                                            
 	//-------------------------------------------------------                                                                                       

	Bool_t isK0s=kFALSE;
	Bool_t isLambda=kFALSE;
	Bool_t isAntiLambda=kFALSE;

	Int_t labelPos    = pTrack->GetLabel();
	Int_t labelNeg    = nTrack->GetLabel();
	
	TClonesArray* AODMCTrackArray =0x0;
	if(fReadMCTruth){
	  fMCEvent= MCEvent();
	  if (fMCEvent){
	    AODMCTrackArray = dynamic_cast<TClonesArray*>(lAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
	    if (AODMCTrackArray == NULL){
	      return;
	      Printf("ERROR: stack not available");
	    }
	  
	    particlePos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelPos)));
	    particleNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelNeg)));
	   
	    PdgPos = particlePos->GetPdgCode(); 
	    PdgNeg = particleNeg->GetPdgCode();        
            
	    labelMotherPos=particlePos->GetMother();
	    labelMotherNeg=particleNeg->GetMother();

	    MotherPos = static_cast<AliAODMCParticle*>(AODMCTrackArray-> At(TMath::Abs(labelMotherPos)));
	    MotherNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray-> At(TMath::Abs(labelMotherNeg)));

	    PdgMotherPos = MotherPos->GetPdgCode();                                
	    PdgMotherNeg = MotherNeg->GetPdgCode();                   
	    
	    isK0s = (PdgPos==211 && PdgNeg==-211 && PdgMotherPos == 310 && PdgMotherNeg == 310 && labelMotherPos==labelMotherNeg  && MotherPos->IsPhysicalPrimary());
	    isLambda = (PdgPos==2212 && PdgNeg==-211 && PdgMotherPos == 3122 && PdgMotherNeg == 3122 && labelMotherPos==labelMotherNeg  && MotherPos->IsPhysicalPrimary());
	    isAntiLambda = (PdgPos==-2212 && PdgNeg==211 && PdgMotherPos == -3122 && PdgMotherNeg == -3122 && labelMotherPos==labelMotherNeg  && MotherPos->IsPhysicalPrimary());

	  }
	}
	//	cout << "\nMC info: "<< endl;
	//	cout << "isK0s " << isK0s << " isLambda " << isLambda << endl;

        //Daughters' Eta
        fV0_etaPos = pTrack->Eta();
        fV0_etaNeg = nTrack->Eta();
        if( pTrack->GetSign() == nTrack->GetSign()) { continue; } // remove like-sign V0s (if any)

        //crossed raws
        double_t lCrosRawsPos = pTrack->GetTPCClusterInfo(2,1);
        double_t lCrosRawsNeg = nTrack->GetTPCClusterInfo(2,1);
        fV0_LeastCRaws = (int) TMath::Min(lCrosRawsPos,lCrosRawsNeg);
        //crossed raws / Findable clusters
        double_t lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrack->GetTPCNclsF())+1e-10);
        double_t lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrack->GetTPCNclsF())+1e-10);
        fV0_LeastCRawsOvF = TMath::Min(lCrosRawsOvFPos,lCrosRawsOvFNeg);

        //dca info
        fV0_DcaPosToPV  = v0->DcaPosToPrimVertex();
        fV0_DcaNegToPV  = v0->DcaNegToPrimVertex();
        fV0_DcaV0Daught = v0->DcaV0Daughters();

        //cosPA
        fV0_V0CosPA = v0->CosPointingAngle(lBestAODPrimVtx);

        //PID
        fV0_NSigPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
        fV0_NSigPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
        fV0_NSigNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
        fV0_NSigNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );

        //distance over total momentum
        fV0_DistOverTotP = v0->DecayLengthV0(lBestPV)/(v0->P()+1e-10);//avoid division by zero

        //filling histos
        fHistos_V0->FillTH1("CosPA",fV0_V0CosPA);
        fHistos_V0->FillTH1("Radius",fV0_V0Rad);
        fHistos_V0->FillTH1("V0DCANegToPV",fV0_DcaNegToPV);
        fHistos_V0->FillTH1("V0DCAPosToPV",fV0_DcaPosToPV);
        fHistos_V0->FillTH1("V0DCAV0Daughters",fV0_DcaV0Daught);
        if(ApplyCuts(0,0,0)){
          fHistos_V0->FillTH2("ImassK0S", fV0_Pt, fV0_InvMassK0s);
	        if (isK0s)      fHistos_V0->FillTH2("ImassK0STrue", fV0_Pt, fV0_InvMassK0s);
        }
        if(ApplyCuts(1,0,0)){
          fHistos_V0->FillTH2("ImassLam", fV0_Pt, fV0_InvMassLam);
	        if (isLambda)      fHistos_V0->FillTH2("ImassLamTrue", fV0_Pt, fV0_InvMassLam);
          if(fV0_DcaV0Daught < 1.0 && fV0_V0CosPA > 0.999 && TMath::Abs(fV0_InvMassK0s-0.497614) > 0.012 && TMath::Abs(fV0_InvMassALam-1.115683) > 0.08 && TMath::Abs(fV0_InvMassLam-1.115683) < 0.002){ 
            fHistos_V0->FillTH2("ResponsePionFromLambda", fV0_Pt, fV0_NSigNegPion);
            fHistos_V0->FillTH2("ResponseProtonFromLambda", fV0_Pt, fV0_NSigPosProton);
          }
        }
        if(ApplyCuts(2, 0,0)){
	  fHistos_V0->FillTH2("ImassALam", fV0_Pt, fV0_InvMassALam);      
	  if (isAntiLambda)      fHistos_V0->FillTH2("ImassALamTrue", fV0_Pt, fV0_InvMassALam);
    }  
  } // end of V0 loop

    // start of cascades part
    //-----------------------


    int ncasc = 0;
    ncasc = lAODevent->GetNumberOfCascades();
    //    cout << "\nnumber of cascade candidates in this event : " << ncasc << endl;

    for (int i_casc = 0; i_casc < ncasc; i_casc++) {
        AliAODcascade *casc = lAODevent->GetCascade(i_casc);
        if (!casc) continue;
	//	cout << " I have a cascade n. " << i_casc <<  endl;
        //get daughter tracks (positive, negative and bachelor)
        AliAODTrack *pTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDaughter(0));
        AliAODTrack *nTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDaughter(1));
        AliAODTrack *bTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDecayVertexXi()->GetDaughter(0));
        if (!pTrackCasc || !nTrackCasc || !bTrackCasc ) { AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...\n"); continue; }
	

	//-------------------------------------------------------                                                                                           
  	//---------MC information--------------------------------                                                                                            
 	//-------------------------------------------------------                                                                                       

	Int_t labelPos    = pTrackCasc->GetLabel();
	Int_t labelNeg    = nTrackCasc->GetLabel();
	Int_t labelBach   = bTrackCasc->GetLabel();

	Bool_t isXiNeg=kFALSE;
	Bool_t isXiPos=kFALSE;
	Bool_t isOmegaNeg=kFALSE;
	Bool_t isOmegaPos=kFALSE;
	
	TClonesArray* AODMCTrackArray =0x0;
	if(fReadMCTruth){
	  fMCEvent= MCEvent();
	  if (fMCEvent){
	    AODMCTrackArray = dynamic_cast<TClonesArray*>(lAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
	    if (AODMCTrackArray == NULL){
	      return;
	      Printf("ERROR: stack not available");
	    }
	  
	    particlePos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelPos)));
	    particleNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelNeg)));
	    particleBach = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelBach)));
	   
	    PdgPos = particlePos->GetPdgCode(); 
	    PdgNeg = particleNeg->GetPdgCode();        
	    PdgBach = particleBach->GetPdgCode();     
            
	    labelMotherPos=particlePos->GetMother();
	    labelMotherNeg=particleNeg->GetMother();
	    labelMotherBach=particleBach->GetMother();
	    MotherPos = static_cast<AliAODMCParticle*>(AODMCTrackArray-> At(TMath::Abs(labelMotherPos)));
	    MotherNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray-> At(TMath::Abs(labelMotherNeg)));
	    MotherBach = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelMotherBach)));
	    PdgMotherPos = MotherPos->GetPdgCode();                                
	    PdgMotherNeg = MotherNeg->GetPdgCode();                   
	    PdgMotherBach = MotherBach->GetPdgCode();                 
	    
	    labelGMotherPos=MotherPos->GetMother();
	    labelGMotherNeg=MotherNeg->GetMother();
	    GMotherPos = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelGMotherPos)));
	    GMotherNeg = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TMath::Abs(labelGMotherNeg)));
	    PdgGMotherPos = GMotherPos->GetPdgCode(); 
	    PdgGMotherNeg = GMotherNeg->GetPdgCode();                                              

	    isXiNeg = (PdgPos==2212 && PdgNeg==-211 && PdgMotherPos == 3122 && PdgMotherNeg == 3122 && labelMotherPos==labelMotherNeg  && !MotherPos->IsPhysicalPrimary() && PdgBach==-211 && PdgMotherBach==3312 && PdgGMotherPos==3312 &&PdgGMotherNeg==3312 && GMotherPos->IsPhysicalPrimary() && labelGMotherPos==labelGMotherNeg && labelGMotherNeg==labelMotherBach);
	    isXiPos = (PdgPos==211 && PdgNeg==-2212 && PdgMotherPos == -3122 &&  PdgMotherNeg == -3122 && labelMotherPos==labelMotherNeg  && !MotherPos->IsPhysicalPrimary() && PdgBach==211 && PdgMotherBach==-3312 && PdgGMotherPos==-3312 && PdgGMotherNeg==-3312 && GMotherPos->IsPhysicalPrimary() && labelGMotherPos==labelGMotherNeg && labelGMotherNeg==labelMotherBach);
	    isOmegaNeg = (PdgPos==2212 && PdgNeg==-211 && PdgMotherPos == 3122 && PdgMotherNeg == 3122 && labelMotherPos==labelMotherNeg  && !MotherPos->IsPhysicalPrimary() && PdgBach==-321 && PdgMotherBach==3334 && PdgGMotherPos==3334 &&PdgGMotherNeg==3334 && GMotherPos->IsPhysicalPrimary() && labelGMotherPos==labelGMotherNeg  && labelGMotherNeg==labelMotherBach);
	    isOmegaPos = (PdgPos==211 && PdgNeg==-2212 && PdgMotherPos == -3122 &&  PdgMotherNeg == -3122 && labelMotherPos==labelMotherNeg  && !MotherPos->IsPhysicalPrimary() && PdgBach==321 && PdgMotherBach==-3334 && PdgGMotherPos==-3334 && PdgGMotherNeg==-3334 && GMotherPos->IsPhysicalPrimary() && labelGMotherPos==labelGMotherNeg && labelGMotherNeg==labelMotherBach);
	    
	  }
	}
	Bool_t isXi = (isXiPos || isXiNeg);
	Bool_t isOmega = (isOmegaPos || isOmegaNeg);

//	cout << "\nMC info: "<< endl;
//	cout << " isXiPos" << isXiPos << " isXiNeg " << isXiNeg << endl;
//	cout << " isOmegaPos" << isOmegaPos << " isOmegaNeg " << isOmegaNeg << endl;

	//----------------------------------------------
       	//---------------- end of MC information---------
       
	//daughters' etas
        fCasc_etaPos = pTrackCasc->Eta();
        fCasc_etaNeg = nTrackCasc->Eta();
        fCasc_etaBac = bTrackCasc->Eta();

        //PID
        fCasc_NSigPosProton = fPIDResponse->NumberOfSigmasTPC( pTrackCasc, AliPID::kProton );
        fCasc_NSigPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrackCasc, AliPID::kPion );
        fCasc_NSigNegProton = fPIDResponse->NumberOfSigmasTPC( nTrackCasc, AliPID::kProton );
        fCasc_NSigNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrackCasc, AliPID::kPion   );
        fCasc_NSigBacPion   = fPIDResponse->NumberOfSigmasTPC( bTrackCasc, AliPID::kPion );
        fCasc_NSigBacKaon   = fPIDResponse->NumberOfSigmasTPC( bTrackCasc, AliPID::kKaon );

	//info for TPC information
	ULong_t pStatus    = pTrackCasc->GetStatus();
	ULong_t nStatus    = nTrackCasc->GetStatus();
	ULong_t bachStatus = bTrackCasc->GetStatus();
	fCasc_isNotTPCRefit = ((pStatus&AliAODTrack::kTPCrefit) == 0 || (nStatus&AliAODTrack::kTPCrefit)  == 0 || (bachStatus&AliAODTrack::kTPCrefit)  == 0);
	//	cout <<"\n\n Refit " <<  fCasc_isNotTPCRefit << endl;
        //crossed raws
        double lCrosRawsPos = pTrackCasc->GetTPCClusterInfo(2,1);
        double lCrosRawsNeg = nTrackCasc->GetTPCClusterInfo(2,1);
        double lCrosRawsBac = bTrackCasc->GetTPCClusterInfo(2,1);
        fCasc_LeastCRaws = (int) ( lCrosRawsPos<lCrosRawsNeg ? std::min(lCrosRawsPos,lCrosRawsBac) : std::min(lCrosRawsNeg,lCrosRawsBac) );
        //crossed raws / Findable clusters
        double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrackCasc->GetTPCNclsF()+1e-10));
        double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrackCasc->GetTPCNclsF()+1e-10));
        double lCrosRawsOvFBac = lCrosRawsBac / ((double)(bTrackCasc->GetTPCNclsF()+1e-10));
        fCasc_LeastCRawsOvF = (float) ( lCrosRawsOvFPos<lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos,lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg,lCrosRawsOvFBac) );

        //DCA info
        fCasc_DcaCascDaught = casc->DcaXiDaughters();
        fCasc_DcaBachToPV = casc->DcaBachToPrimVertex();
        fCasc_DcaPosToPV = casc->DcaPosToPrimVertex();
        fCasc_DcaNegToPV = casc->DcaNegToPrimVertex();
        fCasc_DcaV0Daught = casc->DcaV0Daughters();
        fCasc_DcaV0ToPV = casc->DcaV0ToPrimVertex();

        //cascade and V0 cosine of pointing angle
        //fCasc_CascCosPA = casc->CosPointingAngleXi( lBestAODPrimVtx );
        fCasc_CascCosPA = casc->CosPointingAngleXi( (const Double_t&) lBestPV[0] , (const Double_t&) lBestPV[1] , (const Double_t&) lBestPV[2]);
        fCasc_V0CosPA   = casc->CosPointingAngle( lBestAODPrimVtx ); //check: should not take the secondary vertex insted? Depends on what I cut...
	fCasc_V0CosPAToXi = casc->CosPointingAngle(casc->GetDecayVertexXi()); //Cosine of pointing angle with respect to secondary vertex

        fCasc_yXi = casc->RapXi();
        fCasc_yOm = casc->RapOmega();
        fCasc_charge = (int)casc->ChargeXi();
        fCasc_Pt = sqrt(casc->Pt2Xi());
	fCasc_Ptot= sqrt(casc->Ptot2Xi());

	//Total V0 momentum (there might be an easier way which I'm not aware of)      
	Double_t lBMom[3]={0};
	Double_t  lNMom[3]={0};
	Double_t lPMom[3]={0};
	pTrackCasc->GetPxPyPz( lBMom);
	nTrackCasc->GetPxPyPz( lPMom);
	bTrackCasc->GetPxPyPz( lNMom);
	Float_t lV0TotMomentum = TMath::Sqrt(  TMath::Power( lNMom[0]+lPMom[0] , 2)
					       + TMath::Power( lNMom[1]+lPMom[1] , 2)
					       + TMath::Power( lNMom[2]+lPMom[2] , 2) );


        //distance over total momentum of cascade
	Double_t lPosXi[3]={-1000, -1000, -1000}; //decay vertex of Xi

	lPosXi[0] = casc->DecayVertexXiX();
	lPosXi[1] = casc->DecayVertexXiY();
	lPosXi[2] = casc->DecayVertexXiZ();

        //fCasc_DistOverTotP = casc->DecayLengthXi(lBestPV)/(casc->P()+1e-10);
	//        fCasc_DistOverTotP = casc->DecayLength( lBestAODPrimVtx )/(fCasc_Ptot+1e-10); //this method gives the worng values!!!!!!

	Double_t lXiDecayLength = TMath::Sqrt(
				     TMath::Power( lPosXi[0] - lBestPV[0] , 2) +
				     TMath::Power( lPosXi[1] - lBestPV[1] , 2) +
				     TMath::Power( lPosXi[2] - lBestPV[2] , 2)
				     );

	fCasc_DistOverTotP = lXiDecayLength/(fCasc_Ptot+1e-10);

	//	cout << lXiDecayLength << " def  " << casc->DecayLength( lBestAODPrimVtx) << " defxi " << endl;

	//distance over total momentum of V0 from cascade
	Double_t lPosV0Xi[3]={-1000, -1000, -1000}; //decay vertex of V0 from cascade
	casc->GetXYZ( lPosV0Xi ); 	
	Double_t Casc_V0Dist =  TMath::Sqrt(  TMath::Power( lPosV0Xi[0]-lPosXi[0] , 2)
					    + TMath::Power( lPosV0Xi[1]-lPosXi[1] , 2)
					    + TMath::Power( lPosV0Xi[2]-lPosXi[2] , 2) );
	fCasc_V0DistOverTotP =  Casc_V0Dist/(lV0TotMomentum +1e-10);

        //cascade and V0 radii
	//        fCasc_CascRad =  casc->RadiusSecVtx(); //warning: this method is identical to casc->RadiusV0()
        fCasc_V0Rad = casc->RadiusV0();
	fCasc_CascRad   = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );

        //candidate's invariant mass
        fCasc_InvMassXi = casc->MassXi();
        fCasc_InvMassOm = casc->MassOmega();

	//candidate ctau for cascades and for their V0 daughters
	fCasc_CascCtauXi=1.32171*fCasc_DistOverTotP; //1.32171 GeV/c is Xi+- mass
	fCasc_CascCtauOmega=1.67245*fCasc_DistOverTotP; //1.67245 GeV/c is Omega+- mass
	fCasc_V0Ctau=1.115683*fCasc_V0DistOverTotP; //1.115683 GeV/c is Lambda mass

	//Invmass of Lambda as cascade daughter
	if ( fCasc_charge < 0)
	  fCasc_InvMassLambda   = casc->MassLambda();
	else
	  fCasc_InvMassLambda    = casc->MassAntiLambda();


        //filling histos
	fHistos_Casc->FillTH2("CascCfrPt", fCasc_Pt, casc->Pt());
	fHistos_Casc->FillTH2("CascCfrP", fCasc_Ptot, casc->P());
	fHistos_Casc->FillTH2("CascCosPA", fCasc_CascCosPA,fCasc_charge);
	fHistos_Casc->FillTH2("V0CosPA", fCasc_V0CosPA,fCasc_charge);
	fHistos_Casc->FillTH2("V0CosPAToXi", fCasc_V0CosPAToXi,fCasc_charge);
	fHistos_Casc->FillTH2("CascRadius", fCasc_CascRad,fCasc_charge);
	fHistos_Casc->FillTH2("V0Radius", fCasc_V0Rad,fCasc_charge);
	fHistos_Casc->FillTH2("V0Ctau", fCasc_V0Ctau,fCasc_charge);
	fHistos_Casc->FillTH2("CascPt", fCasc_Pt,fCasc_charge);
	fHistos_Casc->FillTH2("DcaV0Daughters", fCasc_DcaV0Daught,fCasc_charge);
	fHistos_Casc->FillTH2("DcaCascDaughters", fCasc_DcaCascDaught,fCasc_charge);
	fHistos_Casc->FillTH2("DcaV0ToPV", fCasc_DcaV0ToPV,fCasc_charge);
	fHistos_Casc->FillTH2("DcaBachToPV", fCasc_DcaBachToPV,fCasc_charge);
	fHistos_Casc->FillTH2("DcaPosToPV", fCasc_DcaPosToPV,fCasc_charge);
	fHistos_Casc->FillTH2("DcaNegToPV", fCasc_DcaNegToPV,fCasc_charge);
	fHistos_Casc->FillTH2("InvMassLambdaDaughter", fCasc_InvMassLambda,fCasc_charge);

//	cout << "\nMC info after a while: "<< endl;
//	cout << " isXiPos" << isXiPos << " isXiNeg " << isXiNeg << endl;
//	cout << " isOmegaPos" << isOmegaPos << " isOmegaNeg " << isOmegaNeg << endl;
      if (isXi)          fHistos_Casc->FillTH1("XiProgSelections"   ,19, fCasc_charge);
      if (isOmega)       fHistos_Casc->FillTH1("OmegaProgSelections",19, fCasc_charge);


	if(ApplyCuts(3, isXi,0)) {
	    fHistos_Casc->FillTH2("CascyXi", fCasc_yXi, 1.);
	    fHistos_Casc->FillTH2("CascCtauXi", fCasc_CascCtauXi, 1.);
	    fHistos_Casc->FillTH2("ImassXiPlu", fCasc_Pt, fCasc_InvMassXi);
	    if (isXiPos)  fHistos_Casc->FillTH2("ImassXiPluTrue", fCasc_Pt, fCasc_InvMassXi);
    }

	if(ApplyCuts(4, isXi, 0)) {
	    fHistos_Casc->FillTH2("CascyXi", fCasc_yXi, -1.);
	    fHistos_Casc->FillTH2("CascCtauXi", fCasc_CascCtauXi, -1.);
	    fHistos_Casc->FillTH2("ImassXiMin", fCasc_Pt, fCasc_InvMassXi);
	    if (isXiNeg)  fHistos_Casc->FillTH2("ImassXiMinTrue", fCasc_Pt, fCasc_InvMassXi);
    }

	if(ApplyCuts(5, 0, isOmega)) {
	    fHistos_Casc->FillTH2("CascyOmega", fCasc_yOm, 1.);
	    fHistos_Casc->FillTH2("CascCtauOmega", fCasc_CascCtauOmega, 1.);
	    fHistos_Casc->FillTH2("ImassOmPlu", fCasc_Pt, fCasc_InvMassOm);
	    if (isOmegaPos)  fHistos_Casc->FillTH2("ImassOmPluTrue", fCasc_Pt, fCasc_InvMassOm);
    }

	if(ApplyCuts(6, 0, isOmega)) {
	    fHistos_Casc->FillTH2("CascyOmega", fCasc_yOm, -1.);
	    fHistos_Casc->FillTH2("CascCtauOmega", fCasc_CascCtauOmega, -1.);
	    fHistos_Casc->FillTH2("ImassOmMin", fCasc_Pt, fCasc_InvMassOm);
	    if (isOmegaNeg)  fHistos_Casc->FillTH2("ImassOmMinTrue", fCasc_Pt, fCasc_InvMassOm);
    }

  }

  PostData(1, fHistos_eve->GetListOfHistograms()    );
  PostData(2, fHistos_V0->GetListOfHistograms()    );
  PostData(3, fHistos_Casc->GetListOfHistograms()    );

}

//________________________________________________________________________
void AliAnalysisTaskStrAODqa::Terminate(Option_t *)
{

/*    TCanvas *can = new TCanvas("can","can");
    can->cd();
    (fHistos_eve->FindObject("henum"))->Draw();

    TCanvas *can2 = new TCanvas("can2","can2");
    can2->cd();
    (fHistos_K0S->FindObject("himass"))->Draw("E");

    TCanvas *can3 = new TCanvas("can3","can3");
    can3->cd();
    ((TH1D*)fHistos_Lam->FindObject("himass"))->Draw("E");
    ((TH1D*)fHistos_ALam->FindObject("himass"))->SetLineColor(kRed);
    ((TH1D*)fHistos_ALam->FindObject("himass"))->Draw("Esame");
*/

}

//________________________________________________________________________
bool AliAnalysisTaskStrAODqa::ApplyCuts(int part, Bool_t isXi, Bool_t isOmega)
{
    if(part<3){ //we are checking V0s

        // check candidate's rapidity (particle hypothesis' dependent)
        if( (part==0) && TMath::Abs(fV0_yK0S)>0.5 ) return kFALSE;
        if( (part > 0) && TMath::Abs(fV0_yLam)>0.5 ) return kFALSE;
        // check candidate daughters' pseudo-rapidity
        if( TMath::Abs(fV0_etaPos)>0.8 || TMath::Abs(fV0_etaNeg)>0.8 ) return kFALSE;
        // check candidate daughters' crossed TPC raws (note that the checked value is the lowest between the two daughter)
        if( fV0_LeastCRaws<70 ) return kFALSE;
        // check candidate daughters' crossed TPC raws over findable
        if( fV0_LeastCRawsOvF<0.8 ) return kFALSE;
        // check candidate daughters' DCA to Primary Vertex (needs to be large because V0 decay is far from the Primary Vertex)
        if( fV0_DcaPosToPV<0.1 || fV0_DcaNegToPV<0.1) return kFALSE;
        // check candidate daughters' DCA between them (needs to be small because they have to come from the same secondary vertex)
        if( fV0_DcaV0Daught>0.5 ) return kFALSE;
        // check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
        if( fV0_V0Rad<3.0 ) return kFALSE;
        // check the cosine of the Pointing Angle (angle between candidate's momentum and vector connecting Primary and secondary vertices)
        if( fV0_V0CosPA<0.998 ) return kFALSE;
	//reject Lambda candidates when considering K0s
        if( part == 0 && TMath::Abs(fV0_InvMassLam)<0.005 ) return kFALSE;
        // check PID for all daughters (particle hypothesis' dependent)
        if( (part==0) && (TMath::Abs(fV0_NSigPosPion)>3 || TMath::Abs(fV0_NSigNegPion)>3) ) return kFALSE;
        if( (part==1) && (TMath::Abs(fV0_NSigPosProton)>3 || TMath::Abs(fV0_NSigNegPion)>3) ) return kFALSE;
        if( (part==2) && (TMath::Abs(fV0_NSigNegProton)>3 || TMath::Abs(fV0_NSigPosPion)>3) ) return kFALSE;
        // check candidate's proper lifetime (particle hypothesis' dependent). Remember: c*tau = L*m/p
        if( (part==0) && (0.497*fV0_DistOverTotP >20) ) return kFALSE;
	      if( (part>0) && (1.115683*fV0_DistOverTotP >30) ) return kFALSE;

    }

    else { //we are checking cascades
      //check sign of the Cascade
      if (part == 3 && fCasc_charge<0) return kFALSE;
      if (part == 5 && fCasc_charge<0) return kFALSE;
      if (part == 4 && fCasc_charge>0) return kFALSE;
      if (part == 6 && fCasc_charge>0) return kFALSE;

      if (isXi && part<5)           fHistos_Casc->FillTH1("XiProgSelections"   ,1, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",1, fCasc_charge);

      //TPC refit for daughter tracks
            if (fCasc_isNotTPCRefit) return kFALSE;
            if (isXi && part<5)           fHistos_Casc->FillTH1("XiProgSelections"   ,20, fCasc_charge);
            if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",20, fCasc_charge);

      // check candidate's rapidity (particle hypothesis' dependent)
      if( (part<5) && TMath::Abs(fCasc_yXi)>0.5 ) return kFALSE;
      if( (part>=5) && TMath::Abs(fCasc_yOm)>0.5 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,2, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",2, fCasc_charge);

      // check candidate daughters' pseudo-rapidity
      if( TMath::Abs(fCasc_etaPos)>0.8 || TMath::Abs(fCasc_etaNeg)>0.8 || TMath::Abs(fCasc_etaBac)>0.8 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,3, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",3, fCasc_charge);

      // check candidate daughters' crossed TPC raws (note that the checked value is the lowest among the daughters)
      if( fCasc_LeastCRaws<70 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,4, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",4, fCasc_charge);

      // check candidate daughters' crossed TPC raws over findable
      if( fCasc_LeastCRawsOvF<0.8 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,5, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",5, fCasc_charge);

      // check candidate's 2D decay distance from PV (if it is too small, then it's not a weak decay)
      if( fCasc_CascRad<1.0 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,6, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",6, fCasc_charge);

      // check candidate V0 daughter's 2D decay distance from PV (if it is too small, then it's not a weak decay)
      if( fCasc_V0Rad<5.0 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,7, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",7, fCasc_charge);

      // check the cosine of the Pointing Angle for both cascade and V0 (angle between candidate's momentum and vector connecting Primary and secondary vertices)
      if( part<5 && fCasc_CascCosPA<0.9992 ) return kFALSE;
      if( part>=5 && fCasc_CascCosPA<0.9992 ) return kFALSE; //tighter selection for Omega
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,8, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",8, fCasc_charge);

      // if( fCasc_V0CosPA<0.95 ) return kFALSE;
      if( fCasc_V0CosPAToXi<0.99 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,9, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",9, fCasc_charge);

      // check candidate daughters' DCA to Primary Vertex (needs to be large because decay is far from the Primary Vertex)
      if( fCasc_DcaBachToPV<0.17 ) return kFALSE;
      //      if( fCasc_DcaBachToPV>3 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,10, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",10, fCasc_charge);

      if( fCasc_DcaV0ToPV<0.15 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,11, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",11, fCasc_charge);

      // check V0 daughters' DCA to Primary Vertex. Different cut for meson and baryon daughters, so different conditions for + and - candidates
      if( fCasc_charge>0 &&(fCasc_DcaPosToPV < 0.3 || fCasc_DcaNegToPV < 0.11)) return kFALSE;
      if( fCasc_charge<0 &&(fCasc_DcaPosToPV < 0.11 || fCasc_DcaNegToPV < 0.3)) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,12, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",12, fCasc_charge);

      // check V0 daughter's daughters DCA between them (needs to be small because they have to come from the same secondary vertex)
      if( fCasc_DcaV0Daught>1. ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,13, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",13, fCasc_charge);

      // check candidate daughter's DCA between them (needs to be small because they have to come from the same secondary vertex)
      if( fCasc_DcaCascDaught>0.3 ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,14, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",14, fCasc_charge);

      // check candidate V0 daughter's mass difference from nominal Lambda mass
      if( TMath::Abs(fCasc_InvMassLambda-1.115683)>0.005) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,15, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",15, fCasc_charge);

      //XI rejection (only for Omegas)
      if( (part > 4) && TMath::Abs(fCasc_InvMassXi-1.32171)<0.003) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,16, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",16, fCasc_charge);

      // check candidate's proper lifetime (particle hypothesis' dependent). Remember: c*tau = L*m/p
      if( (part<5) && (fCasc_CascCtauXi> (4.91*3)) ) return kFALSE;   //4.91 is the ctau of xi in cm
      if( (part>=5) && (fCasc_CascCtauOmega > (2.461*3)) ) return kFALSE;   //2.461 is the ctau of om in cm
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,17, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",17, fCasc_charge);

      // check PID for all daughters (particle hypothesis' dependent)
      if( (part==3) && (TMath::Abs(fCasc_NSigPosPion)>3 || TMath::Abs(fCasc_NSigNegProton)>3 || TMath::Abs(fCasc_NSigBacPion)>3) ) return kFALSE;
      if( (part==4) && (TMath::Abs(fCasc_NSigNegPion)>3 || TMath::Abs(fCasc_NSigPosProton)>3 || TMath::Abs(fCasc_NSigBacPion)>3) ) return kFALSE;
      if( (part==5) && (TMath::Abs(fCasc_NSigPosPion)>3 || TMath::Abs(fCasc_NSigNegProton)>3 || TMath::Abs(fCasc_NSigBacKaon)>3) ) return kFALSE;
      if( (part==6) && (TMath::Abs(fCasc_NSigNegPion)>3 || TMath::Abs(fCasc_NSigPosProton)>3 || TMath::Abs(fCasc_NSigBacKaon)>3) ) return kFALSE;
      if (isXi && part<5)          fHistos_Casc->FillTH1("XiProgSelections"   ,18, fCasc_charge);
      if (isOmega && part>=5)       fHistos_Casc->FillTH1("OmegaProgSelections",18, fCasc_charge);
    }

    return kTRUE; //survived!

}

//________________________________________________________________________
bool AliAnalysisTaskStrAODqa::ApplyCutsNSigmaTPC(int part)
{
  if( (part==0) && (TMath::Abs(fV0_NSigPosPion)>3 || TMath::Abs(fV0_NSigNegPion)>3) ) return kFALSE;
  if( (part==1) && (TMath::Abs(fV0_NSigPosProton)>3 || TMath::Abs(fV0_NSigNegPion)>3) ) return kFALSE;
  if( (part==2) && (TMath::Abs(fV0_NSigNegProton)>3 || TMath::Abs(fV0_NSigPosPion)>3) ) return kFALSE;
  if( (part==3) && (TMath::Abs(fCasc_NSigPosPion)>3 || TMath::Abs(fCasc_NSigNegProton)>3 || TMath::Abs(fCasc_NSigBacPion)>3) ) return kFALSE;
  if( (part==4) && (TMath::Abs(fCasc_NSigNegPion)>3 || TMath::Abs(fCasc_NSigPosProton)>3 || TMath::Abs(fCasc_NSigBacPion)>3) ) return kFALSE;
  if( (part==5) && (TMath::Abs(fCasc_NSigPosPion)>3 || TMath::Abs(fCasc_NSigNegProton)>3 || TMath::Abs(fCasc_NSigBacKaon)>3) ) return kFALSE;
  if( (part==6) && (TMath::Abs(fCasc_NSigNegPion)>3 || TMath::Abs(fCasc_NSigPosProton)>3 || TMath::Abs(fCasc_NSigBacKaon)>3) ) return kFALSE;

  return kTRUE;
}
