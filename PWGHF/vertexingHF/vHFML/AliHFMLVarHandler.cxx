// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//*************************************************************************
// \class AliHFMLVarHandler
// \brief helper class to handle a tree and variables for ML analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFMLVarHandler.h"
#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliVertexingHFUtils.h"

/// \cond CLASSIMP
ClassImp(AliHFMLVarHandler);
/// \endcond

//________________________________________________________________
AliHFMLVarHandler::AliHFMLVarHandler() : AliHFMLVarHandler(kNsigmaPID)
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliHFMLVarHandler::AliHFMLVarHandler(int PIDopt) : TObject(),
                                                   fPidOpt(PIDopt)
{
    //
    // Standard constructor
    //
    for(unsigned int iProng = 0; iProng < knMaxProngs; iProng++) {
        for(unsigned int iDet = 0; iDet < knMaxDet4Pid; iDet++)
            fPIDrawVector[iProng][iDet] = -999.;
        for(unsigned int iDet = 0; iDet < knMaxDet4Pid + 1; iDet++) {
            for(unsigned int iHypo = 0; iHypo < knMaxHypo4Pid; iHypo++)
                fPIDNsigmaVector[iProng][iDet][iHypo] = -999.;  
        }
    }
}

//________________________________________________________________
AliHFMLVarHandler::~AliHFMLVarHandler()
{
    //
    // Destructor
    //
    if(fTreeVar) 
        delete fTreeVar;
    if(fPIDCombined)
        delete fPIDCombined;
}

//________________________________________________________________
void AliHFMLVarHandler::FillTree() {
    //if fill only signal and not signal/reflection candidate, do not store
    if(fFillOnlySignal && !(fCandType&kSignal) && !(fCandType&kRefl)) {
        fCandType = 0;
    }
    else {      
        fTreeVar->Fill(); 
        fCandType = 0;
    }
}

//________________________________________________________________
void AliHFMLVarHandler::SetCandidateType(bool issignal, bool isbkg, bool isprompt, bool isFD, bool isreflected) {  
    if(issignal) 
        fCandType |= kSignal;
    else 
        fCandType &= ~kSignal;
    if(isbkg) 
        fCandType |= kBkg;
    else 
        fCandType &= ~kBkg;
    if(isprompt) 
        fCandType |= kPrompt;
    else 
        fCandType &= ~kPrompt;
    if(isFD) 
        fCandType |= kFD;
    else 
        fCandType &= ~kFD;
    if(isreflected) 
        fCandType |= kRefl;
    else 
        fCandType &= ~kRefl;
}

//________________________________________________________________
void AliHFMLVarHandler::SetIsSignalWoQuark(bool isSignalWoQuark) {  
    if(isSignalWoQuark) 
        fCandType |= kSignalWoQuark;
    else 
        fCandType &= ~kSignalWoQuark;
}

//________________________________________________________________
void AliHFMLVarHandler::AddCommonDmesonVarBranches() {
    fTreeVar->Branch("cand_type", &fCandType);
    fTreeVar->Branch("inv_mass", &fInvMass);
    fTreeVar->Branch("pt_cand", &fPt);
    fTreeVar->Branch("d_len", &fDecayLength);
    fTreeVar->Branch("d_len_xy", &fDecayLengthXY);
    fTreeVar->Branch("norm_dl_xy", &fNormDecayLengthXY);
    fTreeVar->Branch("cos_p", &fCosP);
    fTreeVar->Branch("cos_p_xy", &fCosPXY);
    fTreeVar->Branch("imp_par_xy", &fImpParXY);
    fTreeVar->Branch("dca", &fDCA);
    if(fEnableBMotherPt)
        fTreeVar->Branch("pt_B", &fPtBMother);
} 

//________________________________________________________________
void AliHFMLVarHandler::AddSingleTrackBranches() {
    for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
        fTreeVar->Branch(Form("pt_prong%d", iProng), &fPtProng[iProng]);
}

//________________________________________________________________
void AliHFMLVarHandler::AddPidBranches(bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF) {
    if(fPidOpt > kBayesAndNsigmaCombPID) {
        AliWarning("Wrong PID setting!");
        return;
    }

    bool useHypo[knMaxHypo4Pid] = {usePionHypo, useKaonHypo, useProtonHypo};
    bool useDet[knMaxDet4Pid] = {useTPC, useTOF};
    TString partHypoName[knMaxHypo4Pid] = {"Pi", "K", "Pr"};
    TString detName[knMaxDet4Pid] = {"TPC", "TOF"};
    TString rawPidName[knMaxDet4Pid] = {"dEdxTPC", "ToF"};

    if(fPidOpt == kBayesPID || fPidOpt == kBayesAndNsigmaCombPID) {
        //initialise AliPIDCombined object for Bayesian PID
        if(!fPIDCombined)
            fPIDCombined = new AliPIDCombined();
        fPIDCombined->SetDefaultTPCPriors(); // only default priors for the moment
        fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC + AliPIDResponse::kDetTOF);
    }

    for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
        if(fPidOpt >= kNsigmaPID && fPidOpt != kNsigmaCombPID && fPidOpt != kBayesPID && fPidOpt != kBayesAndNsigmaCombPID) {
            for(unsigned int iDet = 0; iDet < knMaxDet4Pid; iDet++) {
                if(!useDet[iDet]) 
                    continue;
                for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                    if(!useHypo[iPartHypo]) 
                        continue;
                    fTreeVar->Branch(Form("nsig%s_%s_%d", detName[iDet].Data(), partHypoName[iPartHypo].Data(), iProng),
                                     &fPIDNsigmaVector[iProng][iDet][iPartHypo]);
                }
            }
        }
        if(fPidOpt == kNsigmaCombPID || fPidOpt == kNsigmaDetAndCombPID || fPidOpt == kBayesAndNsigmaCombPID) {
            for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                if(!useHypo[iPartHypo]) 
                    continue;
                fTreeVar->Branch(Form("nsigComb_%s_%d", partHypoName[iPartHypo].Data(), iProng),
                                 &fPIDNsigmaVector[iProng][kCombTPCTOF][iPartHypo]);
            }
        }
        if(fPidOpt == kRawPID || fPidOpt == kRawAndNsigmaPID) {
            for(unsigned int iDet = 0; iDet < knMaxDet4Pid; iDet++) {
                if(!useDet[iDet]) 
                    continue;
                fTreeVar->Branch(Form("%s_%d", rawPidName[iDet].Data(), iProng), &fPIDrawVector[iProng][iDet]);
            }
            if(useTPC) {
                fTreeVar->Branch(Form("pTPC_prong%d", iProng), &fTPCPProng[iProng]);
                fTreeVar->Branch(Form("nTPCclspid_prong%d", iProng), &fNTPCclsPidProng[iProng]);
            }
            if(useTOF) {
                fTreeVar->Branch(Form("pTOF_prong%d", iProng), &fTOFPProng[iProng]);
                fTreeVar->Branch(Form("trlen_prong%d", iProng), &fTrackIntegratedLengthProng[iProng]);
                fTreeVar->Branch(Form("start_time_res_prong%d", iProng), &fStartTimeResProng[iProng]);
            }
        }
        if(fPidOpt == kBayesPID || fPidOpt == kBayesAndNsigmaCombPID) {
            for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                if(!useHypo[iPartHypo]) 
                    continue;
                fTreeVar->Branch(Form("probBayes_%s_%d", partHypoName[iPartHypo].Data(), iProng),
                                 &fPIDNsigmaVector[iProng][kBayesTPCTOF][iPartHypo]);
        
            }
        }
    }
}

//________________________________________________________________
bool AliHFMLVarHandler::SetSingleTrackVars(AliAODTrack* prongtracks[]) {
    //Impact parameters of the prongs are defined as a species dependent variable because the prongs 
    //cannot be obtained in similar way for the different AliAODRecoDecay objects (AliAODTrack cannot
    //be used because of recomputation PV)
    for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
        if(!prongtracks[iProng]) {
            AliWarning("Prong track not found!");
            return false;
        }

    for(unsigned int iProng = 0; iProng < fNProngs; iProng++)
        fPtProng[iProng] = prongtracks[iProng]->Pt();

    return true;
}

//________________________________________________________________
bool AliHFMLVarHandler::SetPidVars(AliAODTrack* prongtracks[], AliAODPidHF* pidhf, bool usePionHypo, bool useKaonHypo, bool useProtonHypo, bool useTPC, bool useTOF) 
{
    if(!pidhf) 
        return false;
    for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
        if(!prongtracks[iProng]) {
            AliWarning("Prong track not found!");
            return false;
        }
    }

    //PID variables
    double sig[knMaxProngs][knMaxDet4Pid][knMaxHypo4Pid];
    double sigComb[knMaxProngs][knMaxHypo4Pid];
    double probBayes[knMaxProngs][knMaxHypo4Pid];
    double rawPID[knMaxProngs][knMaxDet4Pid];
    bool useHypo[knMaxHypo4Pid] = {usePionHypo, useKaonHypo, useProtonHypo};
    bool useDet[knMaxDet4Pid] = {useTPC, useTOF};
    AliPID::EParticleType parthypo[knMaxHypo4Pid] = {AliPID::kPion, AliPID::kKaon, AliPID::kProton};
    
    //compute PID variables for different options
    for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
        if(fPidOpt >= kNsigmaPID) {
            for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                if(useHypo[iPartHypo]) {
                    if(useTPC) {
                        double nSigmaTPC = -999.;
                        pidhf->GetnSigmaTPC(prongtracks[iProng], parthypo[iPartHypo], nSigmaTPC);
                        sig[iProng][kTPC][iPartHypo] = nSigmaTPC;
                    }
                    if(useTOF) {
                        double nSigmaTOF = -999.;
                        pidhf->GetnSigmaTOF(prongtracks[iProng], parthypo[iPartHypo], nSigmaTOF);
                        sig[iProng][kTOF][iPartHypo] = nSigmaTOF;
                    }
                    if((fPidOpt == kNsigmaCombPID || fPidOpt == kNsigmaDetAndCombPID || fPidOpt == kBayesAndNsigmaCombPID) && useTPC && useTOF) {
                        sigComb[iProng][iPartHypo] = AliVertexingHFUtils::CombineNsigmaTPCTOF(sig[iProng][kTPC][iPartHypo], 
                                                                                              sig[iProng][kTOF][iPartHypo]);
                    }
                }
            }
        }
        if(fPidOpt == kRawPID || fPidOpt == kRawAndNsigmaPID) {
            if(useTPC) 
                rawPID[iProng][kTPC] = prongtracks[iProng]->GetTPCsignal();
            if(useTOF) {
                if (!(prongtracks[iProng]->GetStatus() & AliESDtrack::kTOFout) || !(prongtracks[iProng]->GetStatus() & AliESDtrack::kTIME)) {
                    rawPID[iProng][kTOF] = -1.;
                }
                else {
                    float len = prongtracks[iProng]->GetIntegratedLength();
                    if (len < 350.f) {
                        rawPID[iProng][kTOF] = -1.;
                    }
                    else {
                        rawPID[iProng][kTOF] = prongtracks[iProng]->GetTOFsignal();
                        float time0 = pidhf->GetPidResponse()->GetTOFResponse().GetStartTime(prongtracks[iProng]->P());
                        rawPID[iProng][kTOF] -= time0;
                    }
                }
            }
        }
        if(fPidOpt == kBayesPID || fPidOpt == kBayesAndNsigmaCombPID) {
            double probTPCTOF[AliPID::kSPECIES] = {-1.};
            AliPIDResponse* pidrespo = pidhf->GetPidResponse(); //FIXME: TPC postcalibrations for Pb-Pb not implemented in AliPIDResponse
            unsigned int detUsed = fPIDCombined->ComputeProbabilities(prongtracks[iProng], pidrespo, probTPCTOF);
            if (detUsed == (unsigned int)fPIDCombined->GetDetectorMask()) { // Check that TPC-TOF worked
                for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                    if(useHypo[iPartHypo])
                        probBayes[iProng][iPartHypo] = probTPCTOF[parthypo[iPartHypo]];
                }
            }
            else { // if TOF information not available, try only TPC
                fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC);
                detUsed = fPIDCombined->ComputeProbabilities(prongtracks[iProng], pidrespo, probTPCTOF);
                if (detUsed == (unsigned int)fPIDCombined->GetDetectorMask()) { // Check that TPC-only worked. If not, then return -999. as probability
                    for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                        if(useHypo[iPartHypo])
                            probBayes[iProng][iPartHypo] = probTPCTOF[parthypo[iPartHypo]];
                    }
                }
                else {
                    for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                        if(useHypo[iPartHypo])
                            probBayes[iProng][iPartHypo] = -999.;
                    }
                }
                //Reset detector mask for PIDCombined object to TPC+TOF
                fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC + AliPIDResponse::kDetTOF);
            }
        }
    }

    //fill PID arrays for different options
    switch(fPidOpt) {
        case kRawPID:
        {
            for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
                for(int iDet = kTPC; iDet <= kTOF; iDet++) {
                    if(!useDet[iDet]) 
                        continue;
                    fPIDrawVector[iProng][iDet] = rawPID[iProng][iDet];
                }
                if(useTPC) {
                    fTPCPProng[iProng] = prongtracks[iProng]->GetTPCmomentum();
                    fNTPCclsPidProng[iProng] = prongtracks[iProng]->GetTPCsignalN();
                }
                if(useTOF) {
                    fTOFPProng[iProng] = GetTOFmomentum(prongtracks[iProng],pidhf);
                    fTrackIntegratedLengthProng[iProng] = prongtracks[iProng]->GetIntegratedLength();
                    fStartTimeResProng[iProng] = pidhf->GetPidResponse()->GetTOFResponse().GetStartTimeRes(prongtracks[iProng]->P());
                }
            }
            break;
        }
        case kNsigmaPID:
        {
            for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
                for(int iDet = kTPC; iDet <= kTOF; iDet++) {
                    if(!useDet[iDet]) 
                        continue;
                    for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                        if(!useHypo[iPartHypo]) 
                            continue;
                        fPIDNsigmaVector[iProng][iDet][iPartHypo] = sig[iProng][iDet][iPartHypo];
                    }
                }
            }
            break;
        }
        case kNsigmaCombPID:
        {
            for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
                for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                    if(!useHypo[iPartHypo]) 
                        continue;
                    fPIDNsigmaVector[iProng][kCombTPCTOF][iPartHypo] = sigComb[iProng][iPartHypo];
                }
            }
            break;
        }
        case kNsigmaDetAndCombPID:
        {
            for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
                for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                    if(!useHypo[iPartHypo]) 
                        continue;
                    fPIDNsigmaVector[iProng][kCombTPCTOF][iPartHypo] = sigComb[iProng][iPartHypo];
                    for(int iDet = kTPC; iDet <= kTOF; iDet++) {
                        if(!useDet[iDet]) 
                            continue;
                        fPIDNsigmaVector[iProng][iDet][iPartHypo] = sig[iProng][iDet][iPartHypo];
                    }
                }
            }
            break;
        }
        case kRawAndNsigmaPID:
        {
            for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
                for(int iDet = kTPC; iDet <= kTOF; iDet++) {
                    if(!useDet[iDet]) 
                        continue;
                    for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                        if(!useHypo[iPartHypo]) 
                            continue;
                        fPIDNsigmaVector[iProng][iDet][iPartHypo] = sig[iProng][iDet][iPartHypo];
                    }
                    fPIDrawVector[iProng][iDet] = rawPID[iProng][iDet];
                }
                if(useTPC) {
                    fTPCPProng[iProng] = prongtracks[iProng]->GetTPCmomentum();
                    fNTPCclsPidProng[iProng] = prongtracks[iProng]->GetTPCsignalN();
                }
                if(useTOF) {
                    fTOFPProng[iProng] = GetTOFmomentum(prongtracks[iProng],pidhf);
                    fTrackIntegratedLengthProng[iProng] = prongtracks[iProng]->GetIntegratedLength();
                    fStartTimeResProng[iProng] = pidhf->GetPidResponse()->GetTOFResponse().GetStartTimeRes(prongtracks[iProng]->P());
                }
            }
            break;
        }
        case kBayesPID:
        {
            for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
                for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                    if(!useHypo[iPartHypo]) 
                        continue;
                    fPIDNsigmaVector[iProng][kBayesTPCTOF][iPartHypo] = probBayes[iProng][iPartHypo];
                }
            }
            break;
        }
        case kBayesAndNsigmaCombPID:
        {
            for(unsigned int iProng = 0; iProng < fNProngs; iProng++) {
                for(unsigned int iPartHypo = 0; iPartHypo < knMaxHypo4Pid; iPartHypo++) {
                    if(!useHypo[iPartHypo]) 
                        continue;
                    fPIDNsigmaVector[iProng][kCombTPCTOF][iPartHypo] = sigComb[iProng][iPartHypo];
                    fPIDNsigmaVector[iProng][kBayesTPCTOF][iPartHypo] = probBayes[iProng][iPartHypo];
                }
            }
            break;
        }
        default:
            AliWarning("Wrong PID setting!");
            return false;
        break;
    }

    return true;
}

//________________________________________________________________
float AliHFMLVarHandler::GetTOFmomentum(AliAODTrack* track, AliAODPidHF* pidhf)
{
    AliPIDResponse *pidrespo = pidhf->GetPidResponse();
    float t_d = pidrespo->GetTOFResponse().GetExpectedSignal(track, AliPID::kTriton); //largest mass possible with Z=1
    float len = track->GetIntegratedLength();
    float beta_d = len / (t_d * kCSPEED);
    float mass = AliPID::ParticleMassZ(AliPID::kTriton); //largest mass possible with Z=1

    if(TMath::Abs(beta_d - 1.) < 1.e-12) 
        return track->GetTPCmomentum();
    else 
        return mass * beta_d / sqrt(1. - (beta_d * beta_d));
}
