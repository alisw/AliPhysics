// Class AliRsnCutCascade
//
// General Cascade cut for Xi/Omega related Analysis.
// Based on AliRsnCutV0, it contains every cuts for the CutV0 as prefix "V0"
//
// authors: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
//          advised by Beomkyu Kim(kimb@cern.ch)
//

#include <Riostream.h>
#include <TFormula.h>
#include <TBits.h>

#include "AliLog.h"
#include "AliESDtrackCuts.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnCutCascade.h"

ClassImp(AliRsnCutCascade)

AliRsnCutCascade::AliRsnCutCascade(const char *name, Int_t hypothesis, AliPID::EParticleType pid, AliPID::EParticleType pid2, AliPID::EParticleType pid3) :
AliRsnCut(name, AliRsnTarget::kDaughter),
fHypothesis(0),
fMass(0.0),
fMassTolerance(0.01),
fMassToleranceVeto(0.01),
fSwitch(0),
fV0LowRadius(0),
fV0HighRadius(0),
fCascadeLowRadius(0),
fCascadeHighRadius(0),
fV0MinDCAVertex(0.3),
fCascadeMinDCAVertex(0.3),
fV0MaxDCAVertex(0.3),
fCascadeMaxDCAVertex(0.3),
fV0MinCosPointAngle(0.95),
fCascadeMinCosPointAngle(0.95),
fV0MaxDaughtersDCA(0.5),
fCascadeMaxDaughtersDCA(0.5),
fMinTPCcluster(70),
fMaxRapidity(0.8),
fMaxPseudorapidity(1e20),
fPID(pid),
fPID2(pid2),
fPID3(pid3),
fPIDCutV0Proton(0),
fPIDCutV0Pion(0),
fPIDCutBachelor(0),
fESDtrackCuts(0x0),
fCutQuality(Form("%sDaughtersQuality", name)),
fAODTestFilterBit(5)
{
    //
    // Default constructor.
    // Initializes all cuts in such a way that all of them are disabled.
    //
    
    SetHypothesis(hypothesis);
}

//_________________________________________________________________________________________________
AliRsnCutCascade::AliRsnCutCascade(const AliRsnCutCascade &copy) :
AliRsnCut(copy),
fHypothesis(copy.fHypothesis),
fMass(copy.fMass),
fMassTolerance(copy.fMassTolerance),
fMassToleranceVeto(copy.fMassToleranceVeto),
fSwitch(copy.fSwitch),
fV0LowRadius(copy.fV0LowRadius),
fV0HighRadius(copy.fV0HighRadius),
fCascadeLowRadius(copy.fCascadeLowRadius),
fCascadeHighRadius(copy.fCascadeHighRadius),
fV0MinDCAVertex(copy.fV0MinDCAVertex),
fCascadeMinDCAVertex(copy.fCascadeMinDCAVertex),
fV0MaxDCAVertex(copy.fV0MaxDCAVertex),
fCascadeMaxDCAVertex(copy.fCascadeMaxDCAVertex),
fV0MinCosPointAngle(copy.fV0MinCosPointAngle),
fCascadeMinCosPointAngle(copy.fCascadeMinCosPointAngle),
fV0MaxDaughtersDCA(copy.fV0MaxDaughtersDCA),
fCascadeMaxDaughtersDCA(copy.fCascadeMaxDaughtersDCA),
fMinTPCcluster(copy.fMinTPCcluster),
fMaxRapidity(copy.fMaxRapidity),
fMaxPseudorapidity(copy.fMaxPseudorapidity),
fPID(copy.fPID),
fPID2(copy.fPID2),
fPID3(copy.fPID3),
fPIDCutV0Proton(copy.fPIDCutV0Proton),
fPIDCutV0Pion(copy.fPIDCutV0Pion),
fPIDCutBachelor(copy.fPIDCutBachelor),
fESDtrackCuts(copy.fESDtrackCuts),
fCutQuality(copy.fCutQuality),
fAODTestFilterBit(copy.fAODTestFilterBit)
{
    //
    // Copy constructor.
    // Just copy all data member values.:IsSelected: Object is not a V0 (RESONANCES/AliRsnCutCascade.cxx:149)
    //
    
    fCutQuality.SetPtRange(0.15, 1E+20);
    fCutQuality.SetEtaRange(-0.8, 0.8);
    fCutQuality.SetSPDminNClusters(1);
    fCutQuality.SetITSminNClusters(0);
    fCutQuality.SetITSmaxChi2(1E+20);
    fCutQuality.SetTPCminNClusters(fMinTPCcluster);
    fCutQuality.SetTPCmaxChi2(4.0);
    fCutQuality.SetRejectKinkDaughters();
    fCutQuality.SetAODTestFilterBit(5);
    
}

//_________________________________________________________________________________________________
AliRsnCutCascade &AliRsnCutCascade::operator=(const AliRsnCutCascade &copy)
{
    //
    // Assignment operator.
    // Just copy all data member values.
    //
    
    if (this == &copy)
        return *this;
    fHypothesis = copy.fHypothesis;
    fMass = copy.fMass;
    fMassTolerance = copy.fMassTolerance;
    fMassToleranceVeto = copy.fMassToleranceVeto;
    fSwitch = copy.fSwitch;
    fV0LowRadius = copy.fV0LowRadius;
    fV0HighRadius  = copy.fV0HighRadius;
    fCascadeLowRadius = copy.fCascadeLowRadius;
    fCascadeHighRadius  = copy.fCascadeHighRadius;
    fV0MinDCAVertex = copy.fV0MinDCAVertex;
    fCascadeMinDCAVertex = copy.fCascadeMinDCAVertex;
    fV0MaxDCAVertex = copy.fV0MaxDCAVertex;
    fCascadeMaxDCAVertex = copy.fCascadeMaxDCAVertex;
    fV0MinCosPointAngle = copy.fV0MinCosPointAngle;
    fCascadeMinCosPointAngle = copy.fCascadeMinCosPointAngle;
    fV0MaxDaughtersDCA = copy.fV0MaxDaughtersDCA;
    fCascadeMaxDaughtersDCA = copy.fCascadeMaxDaughtersDCA;
    fMinTPCcluster = copy.fMinTPCcluster;
    fMaxRapidity = copy.fMaxRapidity;
    fMaxPseudorapidity = copy.fMaxPseudorapidity;
    fCutQuality = copy.fCutQuality;
    fPID = copy.fPID;
    fPID2 = copy.fPID2;
    fPID3 = copy.fPID3;
    fPIDCutV0Proton = copy.fPIDCutV0Proton;
    fPIDCutV0Pion = copy.fPIDCutV0Pion;
    fPIDCutBachelor = copy.fPIDCutBachelor;
    fESDtrackCuts = copy.fESDtrackCuts;
    fCutQuality = copy.fCutQuality;
    fAODTestFilterBit = copy.fAODTestFilterBit;
    
    return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutCascade::IsSelected(TObject *object)
{
    //:IsSelected: Object is not a V0 (RESONANCES/AliRsnCutCascade.cxx:149)
    
    // Cut checker.
    // Checks the type of object being evaluated
    // and then calls the appropriate sub-function (for ESD or AOD)
    //
    
    // coherence check
    if (!TargetOK(object)) return kFALSE;
    
    // check cast
    AliESDcascade *Xiesd = (AliESDcascade*) fDaughter->Ref2ESDcascade();
    AliAODcascade *Xiaod = (AliAODcascade*) fDaughter->Ref2AODcascade();
    AliDebugClass(1, Form("ClassName: %s, %s, %s",fDaughter->GetRef()->ClassName(),Xiesd,Xiaod));
    
    // operate depending on cast:IsSelected: Object is not a V0 (RESONANCES/AliRsnCutCascade.cxx:149)
    
    if (Xiesd) {
        return CheckESD(Xiesd);
    } else if (Xiaod) {
        return CheckAOD(Xiaod);
    } else {
        AliDebugClass(1, "Object is not a Cascade");
        return kFALSE;
    }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutCascade::CheckESD(AliESDcascade *Xi)
{
    //
    // Check an ESD V0.
    // This is done using the default track checker for ESD.
    // It is declared static, not to recreate it every time.
    //
    
    AliDebugClass(1, "Check ESD");
    if (Xi->GetOnFlyStatus()) {
        AliDebugClass(1, "Rejecting Xi/V0 in 'on fly' status");
        return kFALSE; // if kTRUE, then this V0 is recontructed
    }
    
    // retrieve pointer to owner event
    AliESDEvent *lESDEvent = fEvent->GetRefESD();
    Double_t xPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetX();
    Double_t yPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetY();
    Double_t zPrimaryVertex = lESDEvent->GetPrimaryVertex()->GetZ();
    AliDebugClass(2, Form("Primary vertex: %f %f %f", xPrimaryVertex, yPrimaryVertex, zPrimaryVertex));
    
    // retrieve the V0 daughters and Bachelor Pion/Kaon
    UInt_t lIdxPos      = (UInt_t) TMath::Abs(Xi->GetPindex());
    UInt_t lIdxNeg      = (UInt_t) TMath::Abs(Xi->GetNindex());
    UInt_t lIdxBPi      = (UInt_t) TMath::Abs(Xi->GetBindex());
    AliESDtrack *pTrack = lESDEvent->GetTrack(lIdxPos);
    AliESDtrack *nTrack = lESDEvent->GetTrack(lIdxNeg);
    AliESDtrack *bTrack = lESDEvent->GetTrack(lIdxBPi);
    
    Double_t v0q = 0;
    
    // filter like-sign V0
    if ( TMath::Abs( ((pTrack->GetSign()) - (nTrack->GetSign())) ) < 0.1) {
        AliDebugClass(2, "Failed like-sign V0 check");
        return kFALSE;
    }
    
    // check quality cuts
    if (fESDtrackCuts) {
        AliDebugClass(2, "Checking quality cuts");
        if (!fESDtrackCuts->IsSelected(pTrack)) {
            AliDebugClass(2, "Positive daughter failed quality cuts");
            return kFALSE;
        }
        if (!fESDtrackCuts->IsSelected(nTrack)) {
            AliDebugClass(2, "Negative daughter failed quality cuts");
            return kFALSE;
        }
        if (!fESDtrackCuts->IsSelected(bTrack)) {
            AliDebugClass(2, "Bachelor daughter failed quality cuts");
            return kFALSE;
        }
    }
    
    // topological checks
    if (TMath::Abs(Xi->GetDcaV0Daughters()) > fV0MaxDaughtersDCA) {
        AliDebugClass(2, "Failed check on DCA between V0 daughters");
        return kFALSE;
    }
    if (TMath::Abs(Xi->GetDcaXiDaughters()) > fCascadeMaxDaughtersDCA) {
        AliDebugClass(2, "Failed check on DCA between Cascade daughters");
        return kFALSE;
    }
    if ((TMath::Abs(Xi->GetD(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex)) > fV0MaxDCAVertex) || (TMath::Abs(Xi->GetD(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex)) < fV0MinDCAVertex)) {
        AliDebugClass(2, "Failed check on V0 DCA to primary vertex");
        return kFALSE;
    }
    if (TMath::Abs(Xi->GetDcascade(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex)) > fCascadeMaxDCAVertex) {
        AliDebugClass(2, "Failed check on Cascade DCA to primary vertex");
        return kFALSE;
    }
    
    Double_t XiPosition[3];
    Xi->GetXYZcascade(XiPosition[0], XiPosition[1], XiPosition[2]);
    Double_t V0CPA = TMath::Abs(Xi->GetV0CosineOfPointingAngle(XiPosition[0], XiPosition[1], XiPosition[2]));
    if ( (V0CPA < fV0MinCosPointAngle) || (V0CPA >= 1 ) ) {
        AliDebugClass(2, "Failed check on V0 cosine of pointing angle");
        return kFALSE;
    }
    
    Double_t XiCPA = TMath::Abs(Xi->GetCascadeCosineOfPointingAngle(xPrimaryVertex, yPrimaryVertex, zPrimaryVertex));
    if ( (XiCPA < fCascadeMinCosPointAngle) || (XiCPA >= 1 ) ) {
        AliDebugClass(2, "Failed check on Cascade cosine of pointing angle");
        return kFALSE;
    }
    
    if (TMath::Abs(Xi->Y(fHypothesis)) > fMaxRapidity) {
        AliDebugClass(2, "Failed check on Cascade rapidity");
        return kFALSE;
    }

    if (TMath::Abs(Xi->Eta()) > fMaxPseudorapidity){
        AliDebugClass(2, "Failed check on Cascade pseuorapidity");
        return kFALSE;
    }
    //
    Double_t v0Position[3]; // from $ALICE_ROOT/ANALYSIS/AliESDV0Cuts.cxx
    Xi->GetXYZ(v0Position[0],v0Position[1],v0Position[2]);
    
    Double_t V0radius = TMath::Sqrt(TMath::Power(v0Position[0],2) + TMath::Power(v0Position[1],2));
    Double_t Xiradius = TMath::Sqrt(TMath::Power(XiPosition[0],2) + TMath::Power(XiPosition[1],2));
    if ( ( V0radius < fV0LowRadius ) || ( V0radius > fV0HighRadius ) ) {
        AliDebugClass(2, "Failed V0 fiducial volume");
        return kFALSE;
    }
    if ( ( Xiradius < fCascadeLowRadius ) || ( Xiradius > fCascadeHighRadius ) ) {
        AliDebugClass(2, "Failed Cascade fiducial volume");
        return kFALSE;
    }
    
    // Lifetime cut for negative and positive track
    
    // Total Momentum
    // No way to import V0 momentum, change it to Xi momentum.
    
    Double_t tXimom[3];
    Xi->GetPxPyPz( tXimom[0],tXimom[1],tXimom[2] );
    
    // REMOVE LIFE TIME CUT
    /*
     Double_t fLength = TMath::Sqrt(TMath::Power(v0Position[0]- xPrimaryVertex,2) + TMath::Power(v0Position[1] - yPrimaryVertex,2)+ TMath::Power(v0Position[2]- zPrimaryVertex,2));
     
     if( TMath::Abs(fMass*fLength/lV0TotalMomentum) > fLife)
     {
     AliDebugClass(2, "Failed Lifetime Cut on positive track V0");
     return kFALSE;
     }
     */
    
    Xi->ChangeMassHypothesis(v0q, fHypothesis);
    if ((TMath::Abs(Xi->GetEffMassXi() - fMass)) > fMassTolerance) {
        AliDebugClass(2, Form("Cascade is not in the expected inv mass range, input mass: %f, cut mass: %f",TMath::Abs(Xi->GetEffMassXi() - fMass),fMass));
        return kFALSE;
    }
    
    //Set Switch to kTRUE to use Competing Cascade Rejection
    if(fSwitch){
        if(fHypothesis == kXiMinus) {
            Xi->ChangeMassHypothesis(v0q, kOmegaMinus);
            if ((TMath::Abs(Xi->GetEffMassXi() - 1.6725)) < fMassToleranceVeto) {
                Xi->ChangeMassHypothesis(v0q, kXiMinus);
                return kFALSE;
            }
            Xi->ChangeMassHypothesis(v0q, kXiMinus);
        }
        if(fHypothesis == kXiPlusBar) {
            Xi->ChangeMassHypothesis(v0q, kOmegaPlusBar);
            if ((TMath::Abs(Xi->GetEffMassXi() - 1.6725)) < fMassToleranceVeto) {
                Xi->ChangeMassHypothesis(v0q, kXiPlusBar);
                return kFALSE;
            }
            Xi->ChangeMassHypothesis(v0q, kXiPlusBar);
        }
        if(fHypothesis == kOmegaMinus) {
            Xi->ChangeMassHypothesis(v0q, kXiMinus);
            if ((TMath::Abs(Xi->GetEffMassXi() - 1.3217)) < fMassToleranceVeto) {
                Xi->ChangeMassHypothesis(v0q, kOmegaMinus);
                return kFALSE;
            }
            Xi->ChangeMassHypothesis(v0q, kOmegaMinus);
        }
        if(fHypothesis == kOmegaPlusBar) {
            Xi->ChangeMassHypothesis(v0q, kXiPlusBar);
            if ((TMath::Abs(Xi->GetEffMassXi() - 1.3217)) < fMassToleranceVeto) {
                Xi->ChangeMassHypothesis(v0q, kOmegaPlusBar);
                return kFALSE;
            }
            Xi->ChangeMassHypothesis(v0q, kOmegaPlusBar);
        }
    }
    
    // check PID on proton or antiproton from V0
    
    // check initialization of PID object
    AliPIDResponse *pid = fEvent->GetPIDResponse();
    if (!pid) {
        AliFatal("NULL PID response");
        return kFALSE;
    }
    
    Double_t posnsTPC  = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, fPID));
    Double_t posnsTPC2 = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, fPID2));
    
    Double_t negnsTPC  = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, fPID));
    Double_t negnsTPC2 = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, fPID2));
    
    Double_t bpinsTPC  = TMath::Abs(pid->NumberOfSigmasTPC(bTrack, fPID3));
    
    Double_t maxTPC  = fPIDCutV0Proton;
    Double_t maxTPC2 = fPIDCutV0Pion;
    Double_t maxTPC3 = fPIDCutBachelor;
    
    // applies the cut differently depending on the PID and the momentum
    
    if(fHypothesis==kXiMinus) {
        if (! ((posnsTPC <= maxTPC) && (negnsTPC2 <= maxTPC2) && (bpinsTPC <= maxTPC3)) ) {
            AliDebugClass(2, "Failed check on Cascade PID");
            return kFALSE;
        }
    }
    else if(fHypothesis==kXiPlusBar) {
        if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2) && (bpinsTPC <= maxTPC3)) ) {
            AliDebugClass(2, "Failed check on Cascade PID");
            return kFALSE;
        }
    }
    else if(fHypothesis==kOmegaMinus) {
        if (! ((posnsTPC <= maxTPC) && (negnsTPC2 <= maxTPC2) && (bpinsTPC <= maxTPC3)) ) {
            AliDebugClass(2, "Failed check on Cascade PID");
            return kFALSE;
        }
    }
    else if(fHypothesis==kOmegaPlusBar) {
        if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2) && (bpinsTPC <= maxTPC3)) ) {
            AliDebugClass(2, "Failed check on Cascade PID");
            return kFALSE;
        }
    }
    
    // if we reach this point, all checks were successful
    AliDebugClass(2, "Good Cascade");
    return kTRUE;
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutCascade::CheckAOD(AliAODcascade *Xi)
{
    //
    // Check an AOD Cascade.
    // This is done doing directly all checks, since there is not
    // an equivalent checker for AOD tracks
    //
    
    AliDebugClass(2, "Check AOD");
    if (Xi->GetOnFlyStatus()) {
        AliDebugClass(2, "Rejecting Cascasde in 'on fly' status");
        return kFALSE; // if kTRUE, then this V0 is recontructed
    }
    
    // retrieve pointer to owner event
    AliAODEvent *lAODEvent = fEvent->GetRefAOD();
    Double_t xPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetX();
    Double_t yPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetY();
    Double_t zPrimaryVertex = lAODEvent->GetPrimaryVertex()->GetZ();
    AliDebugClass(2, Form("Primary vertex: %f %f %f", xPrimaryVertex, yPrimaryVertex, zPrimaryVertex));
    
    // retrieve the V0 daughters
    AliAODTrack *pTrack = (AliAODTrack *) (Xi->GetDaughter(0));
    AliAODTrack *nTrack = (AliAODTrack *) (Xi->GetDaughter(1));
    AliAODTrack *bTrack = (AliAODTrack *) (Xi->GetDecayVertexXi()->GetDaughter(0));
    
    // check quality cuts
    UInt_t  filtermapP = 9999;
    UInt_t  filtermapN = 9999;
    UInt_t  filtermapB = 9999;
    filtermapP = pTrack->GetFilterMap();
    filtermapN = nTrack->GetFilterMap();
    filtermapB = bTrack->GetFilterMap();
    
    // filter like-sign V0
    if ( TMath::Abs( ((pTrack->Charge()) - (nTrack->Charge())) ) < 0.1) {
        AliDebugClass(2, "Failed like-sign V0 check");
        return kFALSE;
    }
    
    // check compatibility with expected species hypothesis
    Double_t mass = 0.0;
    
    if((fHypothesis==kXiMinus) || (fHypothesis==kXiPlusBar)) {
        mass = Xi->MassXi();
    }
    else if ((fHypothesis==kOmegaMinus) || (fHypothesis==kOmegaPlusBar)) {
        mass = Xi->MassOmega();
    }
    if ((TMath::Abs(mass - fMass)) > fMassTolerance) {
        AliDebugClass(2, Form("Cascade is not in the expected inv mass range  Mass: %d %f %f", fHypothesis, fMass, mass));
        return kFALSE;
    }
    AliDebugClass(2, Form("Mass: %d %f %f", fHypothesis, fMass, mass));
    
    // topological checks
    if (TMath::Abs(Xi->DcaV0ToPrimVertex()) > fV0MaxDCAVertex || TMath::Abs(Xi->DcaV0ToPrimVertex()) < fV0MinDCAVertex) {
        AliDebugClass(2, Form("Failed check on V0 DCA to primary vertex dca=%f maxdca=%f mindca=%f",TMath::Abs(Xi->DcaV0ToPrimVertex()),fV0MaxDCAVertex,fV0MinDCAVertex));
        return kFALSE;
    }
    
    if ( (Xi->DcaXiToPrimVertex() < fCascadeMinDCAVertex) || (Xi->DcaXiToPrimVertex() > fCascadeMaxDCAVertex) ) {
        AliDebugClass(2, Form("Failed check on Cascade DCA to primary vertex dca=%f maxdca=%f mindca=%f",TMath::Abs(Xi->DcaXiToPrimVertex()),fCascadeMaxDCAVertex,fCascadeMinDCAVertex));
        return kFALSE;
    }
    
    // next cut is effective (should it be in AODV0?)
    Double_t lPosXi[3];
    lPosXi[0] = Xi->DecayVertexXiX();
    lPosXi[1] = Xi->DecayVertexXiY();
    lPosXi[2] = Xi->DecayVertexXiZ();
    
    Double_t V0cospointangle = Xi->CosPointingAngle(lPosXi);
    Double_t Xicospointangle = Xi->CosPointingAngleXi(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
    if (TMath::Abs( V0cospointangle )  < fV0MinCosPointAngle) {
        AliDebugClass(2, "Failed check on V0 cosine of pointing angle");
        return kFALSE;
    }
    if (TMath::Abs( Xicospointangle )  < fCascadeMinCosPointAngle) {
        AliDebugClass(2, "Failed check on Cacade cosine of pointing angle");
        return kFALSE;
    }
    
    // next cut is effective (should it be in AODV0?)
    if (TMath::Abs(Xi->DcaV0Daughters()) > fV0MaxDaughtersDCA) {
        AliDebugClass(2, "Failed check on V0 DCA between daughters");
        return kFALSE;
    }
    if (TMath::Abs(Xi->DcaXiDaughters()) > fCascadeMaxDaughtersDCA) {
        AliDebugClass(2, "Failed check on Cascade DCA between daughters");
        return kFALSE;
    }
    
    if(fHypothesis==kXiMinus || fHypothesis==kXiPlusBar) {
        if (TMath::Abs(Xi->RapXi()) > fMaxRapidity) {
            AliDebugClass(2, "Failed check on Cascade rapidity");
            return kFALSE;
        }
    }
    else if(fHypothesis==kOmegaMinus || fHypothesis==kOmegaPlusBar) {
        if (TMath::Abs(Xi->RapOmega()) > fMaxRapidity) {
            AliDebugClass(2, "Failed check on Cascade rapidity");
            return kFALSE;
        }
    }
    
    Double_t pXi = sqrt(Xi->Ptot2Xi());
    Double_t pzXi = Xi->MomXiZ();
    Double_t etaXi = 0.5*TMath::Log((pXi+pzXi)/(pXi-pzXi+1.e-13));
    if (TMath::Abs(etaXi) > fMaxPseudorapidity){
        AliDebugClass(2, "Failed check on Cascade pseuorapidity");
        return kFALSE;
    }
    
    Double_t V0radius = Xi->RadiusV0();
    Double_t Xiradius = TMath::Sqrt(TMath::Power(lPosXi[0],2) + TMath::Power(lPosXi[1],2));
    
    if ( ( V0radius < fV0LowRadius ) || ( V0radius > fV0HighRadius ) ){
        AliDebugClass(2, "Failed V0 fiducial volume");
        return kFALSE;
    }
    if ( ( Xiradius < fCascadeLowRadius ) || ( Xiradius > fCascadeHighRadius ) ){
        AliDebugClass(2, "Failed Cascade fiducial volume");
        return kFALSE;
    }
    
    //-----------------------------------------------------------
    // check initialization of PID object
    AliPIDResponse *pid = fEvent->GetPIDResponse();
    if (!pid) {
        AliFatal("NULL PID response");
        return kFALSE;
    }
    Double_t posnsTPC   = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, fPID));
    Double_t posnsTPC2  = TMath::Abs(pid->NumberOfSigmasTPC(pTrack, fPID2));
    Double_t negnsTPC   = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, fPID));
    Double_t negnsTPC2  = TMath::Abs(pid->NumberOfSigmasTPC(nTrack, fPID2));
    Double_t bpinsTPC  = TMath::Abs(pid->NumberOfSigmasTPC(bTrack, fPID3));
    
    Double_t maxTPC = fPIDCutV0Proton;
    Double_t maxTPC2 = fPIDCutV0Pion;
    Double_t maxTPC3 =  fPIDCutBachelor;
    
    // applies the cut differently depending on the PID and the momentum
    
    if(fHypothesis==kXiMinus) {
        if (! ((posnsTPC <= maxTPC) && (negnsTPC2 <= maxTPC2) && (bpinsTPC <= maxTPC3)) ) {
            AliDebugClass(2, Form("Failed check on Cascade PID: ptrackPID: %f, nTrackPID %f, btrackPID %f",posnsTPC,negnsTPC2,bpinsTPC));
            return kFALSE;
        }
    }
    else if(fHypothesis==kXiPlusBar) {
        if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2) && (bpinsTPC <= maxTPC3)) ) {
            AliDebugClass(2, Form("Failed check on Cascade PID: ptrackPID: %f, nTrackPID %f, btrackPID %f",posnsTPC,negnsTPC2,bpinsTPC));
            return kFALSE;
        }
    }
    else if(fHypothesis==kOmegaMinus) {
        if (! ((posnsTPC <= maxTPC) && (negnsTPC2 <= maxTPC2) && (bpinsTPC <= maxTPC3)) ) {
            AliDebugClass(2, Form("Failed check on Cascade PID: ptrackPID: %f, nTrackPID %f, btrackPID %f",posnsTPC,negnsTPC2,bpinsTPC));
            return kFALSE;
        }
    }
    else if(fHypothesis==kOmegaPlusBar) {
        if(! ((negnsTPC <= maxTPC) && (posnsTPC2 <= maxTPC2) && (bpinsTPC <= maxTPC3)) ) {
            AliDebugClass(2, Form("Failed check on Cascade PID: ptrackPID: %f, nTrackPID %f, btrackPID %f",posnsTPC,negnsTPC2,bpinsTPC));
            return kFALSE;
        }
    }
    
    //---------------------------------------------------------------
    // if we reach this point, all checks were successful
    AliDebugClass(1, "Good AOD Cascade");
    AliDebugClass(1, Form("Mass: %d %f %f %d %d %d", fHypothesis, fMass, mass, filtermapP, filtermapN, filtermapB));
    return kTRUE;
}

//_________________________________________________________________________________________________
void AliRsnCutCascade::Print(const Option_t *) const
{
    //
    // Print information on this cut
    //
}
