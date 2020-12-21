#include "TChain.h"
#include "TH2F.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDHandler.h"
#include "AliAODMCHeader.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "AliAODMCParticle.h"
#include "AliGenHijingEventHeader.h"

#include "AliAnalysisTaskHFEEfficiency.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALGeometry.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliCentrality.h"
#include "AliOADBContainer.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

#include "AliGenPythiaEventHeader.h"

ClassImp(AliAnalysisTaskHFEEfficiency)
//________________________________________________________________________
AliAnalysisTaskHFEEfficiency::AliAnalysisTaskHFEEfficiency(const char *name)
: AliAnalysisTaskSE(name)
,fESD(0)
,fAOD(0)
,fVevent(0)
,fOutputList(0)
,fTrackCuts(new AliESDtrackCuts)
,fNoTrks(0)
,fCuts(0)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fRejectKinkMother(kFALSE)
,fNoEvents(0)
,fMCPdgmom(0)
,fMCElecPtPhoto(0)
,fMCElecPtIncl(0)
,fElecNos(0)
,fInclusiveElecPt(0)
,fPhotoElecPt(0)
,fTruePhotoElecPt(0)
,fMCtagPhotoElecPtAll(0)
,fMCtagTruePhotoElecPtAll(0)
,fMCtagGammaElecPtAll(0)
,fMCtagPi0ElecPtAll(0)
,fMCtagEtaElecPtAll(0)
,fSingleElecPt(0)
,fTrueSingleElecPt(0)
,fInvmassLS(0)
,fInvmassULS(0)
,fInvmassULSNHFE(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fTrackPtBefTrkCuts(0)
,fTrackPtAftTrkCuts(0)
,fTPCnsigma(0)
,fITSnsigma(0)
,fITSnsigmaElectron(0)
,fTOFnsigma(0)
,fTrkEovPBef(0)
,fdEdxBef(0)
,fPhotoElecCandiPdgMom(0)
,fInvmassULSNHFEAllPartn(0)
,fInvmassLSNHFEAllPartn(0)
,fCentralityPass(0)
,fCentralityNoPass(0)
,fPi0EtaSpectra(0)
,fHFEDenominator(0)
,fHFENumeratorTRKCUTS_Filt(0)
,fHFENumeratorTRKCUTS_ITSTPC(0)
//,fHFENumeratorTRKCUTS_ITStrkCut(0)
,fHFENumeratorTRKCUTS(0)
,fHFENumeratorTOF(0)
,fHFENumeratorITS(0)
,fHFENumerator(0)
,fInvMULS(0)
,fIncl(0)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(0)
,fkCentralityMethod(0)
,fVZEROA(0)
,fVZEROC(0)
,fAssoTPCCluster(0)
,fAssoITSRefit(0)
,fAssopTMin(0)
,fStackloop(0)
,fTPCnsigmaAft(0)
,fTPCnsigmaVSptAft(0)
,fTPCnsigmaAftITSTOF(0)
,fMassCut(0)
,fSparseMassULS(0)
,fSparseMassLS(0)
,fakepT(0)
,fakepTgraterFive(0)
,WeightsForEnhanced(0)
,fNoEventsStackHFE(0)
,fminITSnsigmaLowpT(-1)
,fmaxITSnsigmaLowpT(1)
,fminITSnsigmaHighpT(-2)
,fmaxITSnsigmaHighpT(2)
,fminTPCnsigmaLowpT(-1)
,fmaxTPCnsigmaLowpT(3)
,fminTPCnsigmaHighpT(0)
,fmaxTPCnsigmaHighpT(3)
,fminTOFnSigma(-2)
,fmaxTOFnSigma(2)
,fTPCS(0)
,fULS(0)
,fLS(0)
,Rconv_pT(0)
,tiltup(0)
,tiltdw(0)
,fcentral(0)
,fsemicentral(0)
,WeightsForHF(0)
,fmineta(-0.8)
,fmaxeta(0.8)
{
    
    fPID = new AliHFEpid("hfePid");
    
    //  fvalueTrkmat = new Double_t[9];
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //  DefineOutput(3, TTree::Class());
}
//________________________________________________________________________
AliAnalysisTaskHFEEfficiency::AliAnalysisTaskHFEEfficiency()
:AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskHFEEfficiency")
,fESD(0)
,fAOD(0)
,fVevent(0)
,fOutputList(0)
,fTrackCuts(new AliESDtrackCuts)
,fNoTrks(0)
,fCuts(0)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fRejectKinkMother(kFALSE)
,fNoEvents(0)
,fMCPdgmom(0)
,fMCElecPtPhoto(0)
,fMCElecPtIncl(0)
,fElecNos(0)
,fInclusiveElecPt(0)
,fPhotoElecPt(0)
,fTruePhotoElecPt(0)
,fMCtagPhotoElecPtAll(0)
,fMCtagTruePhotoElecPtAll(0)
,fMCtagGammaElecPtAll(0)
,fMCtagPi0ElecPtAll(0)
,fMCtagEtaElecPtAll(0)
,fSingleElecPt(0)
,fTrueSingleElecPt(0)
,fInvmassLS(0)
,fInvmassULS(0)
,fInvmassULSNHFE(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fTrackPtBefTrkCuts(0)
,fTrackPtAftTrkCuts(0)
,fTPCnsigma(0)
,fITSnsigma(0)
,fITSnsigmaElectron(0)
,fTOFnsigma(0)
,fTrkEovPBef(0)
,fdEdxBef(0)
,fPhotoElecCandiPdgMom(0)
,fInvmassULSNHFEAllPartn(0)
,fInvmassLSNHFEAllPartn(0)
,fCentralityPass(0)
,fCentralityNoPass(0)
,fPi0EtaSpectra(0)
,fHFEDenominator(0)
,fHFENumeratorTRKCUTS_Filt(0)
,fHFENumeratorTRKCUTS_ITSTPC(0)
//,fHFENumeratorTRKCUTS_ITStrkCut(0)
,fHFENumeratorTRKCUTS(0)
,fHFENumeratorTOF(0)
,fHFENumeratorITS(0)
,fHFENumerator(0)
,fInvMULS(0)
,fIncl(0)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(0)
,fkCentralityMethod(0)
,fVZEROA(0)
,fVZEROC(0)
,fAssoTPCCluster(0)
,fAssoITSRefit(0)
,fAssopTMin(0)
,fStackloop(0)
,fTPCnsigmaAft(0)
,fTPCnsigmaVSptAft(0)
,fTPCnsigmaAftITSTOF(0)
,fMassCut(0)
,fSparseMassULS(0)
,fSparseMassLS(0)
,fakepT(0)
,fakepTgraterFive(0)
,WeightsForEnhanced(0)
,fNoEventsStackHFE(0)
,fminITSnsigmaLowpT(-1)
,fmaxITSnsigmaLowpT(1)
,fminITSnsigmaHighpT(-2)
,fmaxITSnsigmaHighpT(2)
,fminTPCnsigmaLowpT(-1)
,fmaxTPCnsigmaLowpT(3)
,fminTPCnsigmaHighpT(0)
,fmaxTPCnsigmaHighpT(3)
,fminTOFnSigma(-2)
,fmaxTOFnSigma(2)
,fTPCS(0)
,fULS(0)
,fLS(0)
,Rconv_pT(0)
,tiltup(0)
,tiltdw(0)
,fcentral(0)
,fsemicentral(0)
,WeightsForHF(0)
,fmineta(-0.8)
,fmaxeta(0.8)
{
    
    fPID = new AliHFEpid("hfePid");
    //  fvalueTrkmat = new Double_t[9];
    // Constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //DefineOutput(3, TTree::Class());
}
//________________________________________________________________________
AliAnalysisTaskHFEEfficiency::~AliAnalysisTaskHFEEfficiency()
{
    
    
    delete fOutputList;
    delete fTrackCuts;
    //  delete fSparseTrkMat;
    //  delete [] fvalueTrkmat;
    delete fPID;
    delete fCFM;
    delete fPIDqa;
}
//________________________________________________________________________
void AliAnalysisTaskHFEEfficiency::UserExec(Option_t *)
{
    //Main loop
    //Called for each event
    // create pointer to event
    
    //    cout <<   "In UserExec()      " << endl;
    
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    
    if(!(fESD || fAOD)){
        printf("ERROR: fESD & fAOD not available\n");
        return;
    }
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVEvent not available\n");
        return;
    }
    
    AliKFParticle::SetField(fVevent->GetMagneticField());
    
    AliMCEvent *mcEvent;
    AliMCEventHandler *eventHandler;
    AliAODMCHeader *mcHeader;
    TClonesArray  *mcArray;
    //Int_t fNBG =0;
    // load MC array
    if(IsAODanalysis()){
        mcArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        if(!mcArray){
            AliError("Array of MC particles not found");
            return;
        }
        if(mcArray->GetEntries() < 1) return;
        
        
        mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!mcHeader) {
            AliError("Could not find MC Header in AOD");
            return;
        }
        //    TList *List = mcHeader->GetCocktailHeaders();
        //    AliGenHijingEventHeader* hijingH = dynamic_cast<AliGenHijingEventHeader*>(List->FindObject("Hijing"));
        //      if (!hijingH) cout << "no Genhijing header" <<endl;
        //      fNBG = hijingH->NProduced();
    }
    
    else{
        eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
        if (!eventHandler) {
            Printf("ERROR: Could not retrieve MC event handler");
            return;
        }
        
        mcEvent = eventHandler->MCEvent();
        if (!mcEvent) {
            Printf("ERROR: Could not retrieve MC event");
            return;
        }
    }
    
    if(!fCuts){
        AliError("HFE cuts not available");
        return;
    }
    
    if(!fPID->IsInitialized()){
        // Initialize PID with the given run number
        AliWarning("PID not initialised, get from Run no");
        if(IsAODanalysis())fPID->InitializePID(fAOD->GetRunNumber());
        else fPID->InitializePID(fESD->GetRunNumber());
    }
    
    // centrality selection
    Bool_t pass = kFALSE;
    CheckCentrality(fVevent,pass);
    if(!pass)return;
    
    
    
    double mcZvertex = mcHeader->GetVtxZ();//mcEvent->GetPrimaryVertex()->GetZ();
    
    //-------------AOD HFE spectrum--------------
    //    Bool_t MChijingPrim = kFALSE;
    if(TMath::Abs(mcZvertex) < 10){
        
        fNoEventsStackHFE->Fill(0);
        
        
        
        for(Int_t imcArrayL=0; imcArrayL< mcArray->GetEntries(); imcArrayL++){
            
            //
            //            // Fill Here after HFE track cuts
            Double_t ForHFEeffGen[4];
            Bool_t ShouldiFillit = kFALSE;
            //            FillNumerator(mcArray,track,ForHFEeffGen);
            //            fHFEDenominator->Fill(ForHFEeffGen);
            
            AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)mcArray->At(imcArrayL);
            if(AODMCtrack->Eta() < fmineta || AODMCtrack->Eta() > fmaxeta) continue;
            
            Int_t from_D_or_B;
            if(TMath::Abs(AODMCtrack->GetPdgCode())==11){
                
                Int_t fromhijing = kElse;
                Int_t Prim = GetPrimary(imcArrayL, mcArray);
                AliAODMCParticle *AODMCtrackPrim = (AliAODMCParticle*)mcArray->At(Prim);
                Int_t trkIndexPrim = AODMCtrackPrim->GetLabel();//gives index of the particle in original MCparticle
                if(IsFromBGEventAOD(trkIndexPrim)) fromhijing = kHijing;//check if the particle comes from hijing or from enhanced event
                
                Int_t fMCmom_HFE = AODMCtrack->GetMother();
                AliAODMCParticle *MCMother_HFE = (AliAODMCParticle*)mcArray->At(fMCmom_HFE);
                Int_t motherpdgforHFE = MCMother_HFE->GetPdgCode();
                
                if ( (int(TMath::Abs(motherpdgforHFE)/100.)%10) == 5 || (int(TMath::Abs(motherpdgforHFE)/1000.)%10) == 5 ) {
                    from_D_or_B = kelectronBeauty;
                    ForHFEeffGen[0] = AODMCtrack->Pt();
                    ForHFEeffGen[1] = fromhijing;
                    ForHFEeffGen[2] = from_D_or_B;
                    ForHFEeffGen[3] = AODMCtrack->Phi();
                    ShouldiFillit = kTRUE;
                }
                if ( (int(TMath::Abs(motherpdgforHFE)/100.)%10) == 4 || (int(TMath::Abs(motherpdgforHFE)/1000.)%10) == 4 ) {
                    from_D_or_B = kelectronCharm;
                    ForHFEeffGen[0] = AODMCtrack->Pt();
                    ForHFEeffGen[1] = fromhijing;
                    ForHFEeffGen[2] = from_D_or_B;
                    ForHFEeffGen[3] = AODMCtrack->Phi();
                    ShouldiFillit = kTRUE;
                }
                
                
            }//end is an electron
            if(ShouldiFillit){
                Double_t HFweights = 1;
                if(WeightsForHF){HFweights = GiveHFWeight(AODMCtrack->Pt());}
                fHFEDenominator->Fill(ForHFEeffGen,HFweights);
            }
            
        }//end loop stack HFE
    }//end vertex selection
    
    
    
    
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t pVtxZ = -999;
    pVtxZ = pVtx->GetZ();
    if(TMath::Abs(pVtxZ)>10) return;
    Int_t fNOtrks =  fVevent->GetNumberOfTracks();
    if(fNOtrks<2) return;
    fNoEvents->Fill(0);
    
    AliPIDResponse *pidResponse = fInputHandler->GetPIDResponse();
    if(!pidResponse){
        AliDebug(1, "Using default PID Response");
        pidResponse = AliHFEtools::GetDefaultPID(HasMCData(), fInputEvent->IsA() == AliAODEvent::Class());
    }
    
    fPID->SetPIDResponse(pidResponse);
    fCFM->SetRecEventInfo(fVevent);
    //  Double_t Bz = fInputEvent->GetMagneticField();
    
    Double_t qaweights[5];
    if(fStackloop){
        //-------------AOD pi0 and eta spectrum--------------
        //    Bool_t MChijingPrim = kFALSE;
        
        for(Int_t imcArrayL=0; imcArrayL< mcArray->GetEntries(); imcArrayL++){
            AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)mcArray->At(imcArrayL);
            //        Int_t PDGcode = TMath::Abs(AODMCtrack->GetPdgCode());
            if(AODMCtrack->Eta() < -1.2 || AODMCtrack->Eta() > 1.2) continue;
            
            
            // Pt
            qaweights[0] = AODMCtrack->Pt();
            // See from where it comes
            Int_t fromhijing = kElse;
            Int_t Prim = GetPrimary(imcArrayL, mcArray);
            AliAODMCParticle *AODMCtrackPrim = (AliAODMCParticle*)mcArray->At(Prim);
            Int_t trkIndexPrim = AODMCtrackPrim->GetLabel();//gives index of the particle in original MCparticle array
            if(IsFromBGEventAOD(trkIndexPrim)) fromhijing = kHijing;//check if the particle comes from hijing or from enhanced event
            qaweights[2]=fromhijing;
            
            // What pdg
            qaweights[1]=-1.;
            if (TMath::Abs(AODMCtrack->GetPdgCode())==111) qaweights[1]=0.2;
            if (TMath::Abs(AODMCtrack->GetPdgCode())==221) qaweights[1]=1.2;
            if (TMath::Abs(AODMCtrack->GetPdgCode())==22) qaweights[1]=2.2;
            
            // What type
            Int_t type = GetPi0EtaType(AODMCtrack,mcArray);
            qaweights[3]=type;
            
            // Fill
            if(qaweights[1]>0.) fPi0EtaSpectra->Fill(qaweights);
            
            
            
        }//end stack
    }//end if Bool stack
    
    
    // Look for kink mother for AOD
    Double_t *listofmotherkink =0;
    Int_t numberofvertices = 0, numberofmotherkink = 0;
    if(IsAODanalysis()){
        numberofvertices = fAOD->GetNumberOfVertices();
        listofmotherkink = new Double_t[numberofvertices];
        for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
            AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
            if(!aodvertex) continue;
            if(aodvertex->GetType()==AliAODVertex::kKink) {
                AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
                if(!mother) continue;
                Int_t idmother = mother->GetID();
                listofmotherkink[numberofmotherkink] = idmother;
                numberofmotherkink++;
            }
        }
    }
    
    //   cout <<   "After  kink      " << endl;
    
    Double_t incls[3];
    
    
    //------Reco stack-------------
    // Track loop
    for (Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) {
        AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
        if (!Vtrack) {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
        AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
        
        if(IsAODanalysis()) if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
        if(track->Eta() < fmineta || track->Eta() > fmaxeta) continue;
        
        // Fill Here after HFE track cuts
        Double_t ForHFEeffRecoFilt[4];
        Bool_t canIfilltrkFilt = FillNumerator(mcArray,track,ForHFEeffRecoFilt);
        if(canIfilltrkFilt){fHFENumeratorTRKCUTS_Filt->Fill(ForHFEeffRecoFilt);}
        
        fTrackPtBefTrkCuts->Fill(track->Pt());
        // RecKine: ITSTPC cuts
        if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
        
        // Reject kink mother
        if(fRejectKinkMother) { // Quick and dirty fix to reject both kink mothers and daughters
            if(IsAODanalysis()){
                Bool_t kinkmotherpass = kTRUE;
                for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
                    if(track->GetID() == listofmotherkink[kinkmother]) {
                        kinkmotherpass = kFALSE;
                        continue;
                    }
                }
                if(!kinkmotherpass) continue;
            }
            else{
                if(etrack->GetKinkIndex(0) != 0) continue;
            }
        }
        
        Double_t HFweightsReco = 1;
        if(WeightsForHF){HFweightsReco = GiveHFWeight(track->Pt());}
        
        // Fill Here after HFE track cuts
        Double_t ForHFEeffRecoITSTPCRef[4];
        Bool_t canIfilltrkITSTPC = FillNumerator(mcArray,track,ForHFEeffRecoITSTPCRef);
        if(canIfilltrkITSTPC){fHFENumeratorTRKCUTS_ITSTPC->Fill(ForHFEeffRecoITSTPCRef,HFweightsReco);}
        
        
        // RecPrim
        //      if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
        // HFEcuts: ITS layers cuts
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
        
        
        // Fill Here after HFE track cuts
        //          Double_t ForHFEeffRecoITStrkCut[4];
        //          Bool_t canIfilltrkITStrkCut = FillNumerator(mcArray,track,ForHFEeffRecoITStrkCut);
        //          if(canIfilltrkITStrkCut){fHFENumeratorTRKCUTS_ITStrkCut->Fill(ForHFEeffRecoITStrkCut);}
        
        
        // HFE cuts: TPC PID cleanup
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
        // Fill Here after HFE track cuts
        Double_t ForHFEeffRecoTrackCuts[4];
        Bool_t canIfilltrk = FillNumerator(mcArray,track,ForHFEeffRecoTrackCuts);
        if(canIfilltrk){fHFENumeratorTRKCUTS->Fill(ForHFEeffRecoTrackCuts,HFweightsReco);}
        
        
        fTrackPtAftTrkCuts->Fill(track->Pt());
        
        
        
        Double_t MCtrkPt=-999;
        Int_t trkLabel = TMath::Abs(track->GetLabel());
        // if(trkLabel <0) continue;
        
        
        AliAODMCParticle *MCtrack = (AliAODMCParticle*)mcArray->At(trkLabel);
        MCtrkPt = MCtrack->Pt();
        
        Double_t  p = -999, pt = -999, dEdx=-999, fTPCnSigma=0, fTOFnSigma=0, fITSnSigma=0;
        p = track->P();
        pt= track->Pt();
        dEdx = track->GetTPCsignal();
        
        if(track->GetTPCsignalN() < fTPCS) continue;
        fTPCnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(track, AliPID::kElectron) : 1000;
        fITSnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasITS(track, AliPID::kElectron) : 1000;
        fTOFnSigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTOF(track, AliPID::kElectron) : 1000;
        
        
        fdEdxBef->Fill(p,dEdx);
        fTPCnsigma->Fill(p,fTPCnSigma);
        fITSnsigma->Fill(p,fITSnSigma);
        fTOFnsigma->Fill(p,fTOFnSigma);
        
        
        
        Bool_t MCElectron=kFALSE;
        Int_t PrimElectron = GetPrimary(trkLabel, mcArray);
        AliAODMCParticle *AODMCElectron = (AliAODMCParticle*)mcArray->At(PrimElectron);
        Int_t trkIndexElectron = AODMCElectron->GetLabel();//gives index of the particle in original MCparticle array
        MCElectron = IsFromBGEventAOD(trkIndexElectron);//check if the particle comes from hijing or from enhanced event
        if(MCElectron){
            if(TMath::Abs(MCtrack->GetPdgCode()) == 11){
                fITSnsigmaElectron->Fill(p,fITSnSigma);
            }
        }
        
        //Electron id
        if(fTOFnSigma < fminTOFnSigma || fTOFnSigma > fmaxTOFnSigma) continue;
        
        // Fill Here after TOF cuts
        Double_t ForHFEeffRecoTOF[4];
        Bool_t canIfillTOF = FillNumerator(mcArray,track,ForHFEeffRecoTOF);
        if(canIfillTOF){fHFENumeratorTOF->Fill(ForHFEeffRecoTOF,HFweightsReco);}
        
        if( pt < 1.5){
            if(fITSnSigma < fminITSnsigmaLowpT || fITSnSigma > fmaxITSnsigmaLowpT)continue;
        }
        if( pt >= 1.5){
            if(fITSnSigma < fminITSnsigmaHighpT || fITSnSigma > fmaxITSnsigmaHighpT)continue;
        }//cuts on nsigma ITS
        fTPCnsigmaAftITSTOF->Fill(track->P(),fTPCnSigma);
        
        // Fill Here after ITS cuts
        Double_t ForHFEeffRecoITS[4];
        Bool_t canIfillITS = FillNumerator(mcArray,track,ForHFEeffRecoITS);
        if(canIfillITS){fHFENumeratorITS->Fill(ForHFEeffRecoITS,HFweightsReco);}
        
        if(pt < 1.5){
            if(fTPCnSigma < fminTPCnsigmaLowpT || fTPCnSigma > fmaxTPCnsigmaLowpT) continue;
        }//cuts on nsigma tpc
        if( pt >= 1.5){
            if(fTPCnSigma < fminTPCnsigmaHighpT || fTPCnSigma > fmaxTPCnsigmaHighpT) continue;
        }//cuts on nsigma tpc
        
        // Fill Here after TPC cuts
        Double_t ForHFEeffReco[4];
        Bool_t canIfill = FillNumerator(mcArray,track,ForHFEeffReco);
        if(canIfill){fHFENumerator->Fill(ForHFEeffReco,HFweightsReco);}
        
        
        
        fTPCnsigmaAft->Fill(track->P(),fTPCnSigma);
        fTPCnsigmaVSptAft->Fill(pt,fTPCnSigma);
        //=============FROM HERE efficiency NHFE===================================
        
        Int_t PhototrkLabel_forfake = TMath::Abs(track->GetLabel());
        //  if(PhototrkLabel_forfake <0) continue;
        AliAODMCParticle *MCPhotoEle_forfake = (AliAODMCParticle*)mcArray->At(PhototrkLabel_forfake);
        Int_t fMCmom_fake = MCPhotoEle_forfake->GetMother();
        AliAODMCParticle *MCMother_forfake = (AliAODMCParticle*)mcArray->At(fMCmom_fake);
        if(TMath::Abs(MCPhotoEle_forfake->GetPdgCode()) == 11 && MCMother_forfake->GetPdgCode()==22){
            Double_t fVx = MCPhotoEle_forfake->Xv();
            Double_t fVy = MCPhotoEle_forfake->Yv();
            Double_t Rconv = TMath::Sqrt(fVx*fVx+fVy*fVy);
            //   cout << "production X vertex  ==  " << fVx << endl;
            //   cout << "production Y vertex  ==  " << fVy << endl;
            
            //    cout << "production vertex  ==  " << Rconv << endl;
            // Double_t Rconv = MCPhotoEle_forfake->R(); // production vertex of particles.
            Rconv_pT->Fill(MCPhotoEle_forfake->Pt(),Rconv);
            
            fakepT->Fill(MCPhotoEle_forfake->Pt());
            if(Rconv > 5)fakepTgraterFive->Fill(MCPhotoEle_forfake->Pt());
        }//electron with gamma mother
        
        
        
        
        Bool_t fFlagMCPhotonicElec = kFALSE;
        
        fElecNos->Fill(1);
        
        //MC truth All photonic electron
        Int_t PhototrkLabel = TMath::Abs(track->GetLabel());
        //   if(PhototrkLabel <0) continue;
        
        AliAODMCParticle *MCPhotoEleAll = (AliAODMCParticle*)mcArray->At(PhototrkLabel);
        Int_t fMCmomAll = -999;
        if(TMath::Abs(MCPhotoEleAll->GetPdgCode())!=11) continue;
        
        fElecNos->Fill(2);
        fInclusiveElecPt->Fill(track->Pt());
        fMCmomAll = MCPhotoEleAll->GetMother();
        AliAODMCParticle *MCMotherAll = (AliAODMCParticle*)mcArray->At(fMCmomAll);
        
        
        if(MCMotherAll->GetPdgCode()==22 || MCMotherAll->GetPdgCode()==111 || MCMotherAll->GetPdgCode()==221)
        {
            fFlagMCPhotonicElec = kTRUE;
            fElecNos->Fill(3);
            fMCtagPhotoElecPtAll->Fill(track->Pt());
            fMCtagTruePhotoElecPtAll->Fill(MCPhotoEleAll->Pt());
            
            
            
            if(MCMotherAll->GetPdgCode()==22) fMCtagGammaElecPtAll->Fill(track->Pt());
            if(MCMotherAll->GetPdgCode()==111) fMCtagPi0ElecPtAll->Fill(track->Pt());
            if(MCMotherAll->GetPdgCode()==221) fMCtagEtaElecPtAll->Fill(track->Pt());
            
            // Pt
            incls[0] = track->Pt();
            
            //--------Weight calculation--------------
            //Check ele from hijing or enhanced
            Bool_t MChijingAll=kFALSE;
            Int_t PrimAll = GetPrimary(PhototrkLabel, mcArray);
            AliAODMCParticle *AODMCtrackPrimAll = (AliAODMCParticle*)mcArray->At(PrimAll);
            Int_t trkIndexPrimAll = AODMCtrackPrimAll->GetLabel();//gives index of the particle in original MCparticle array
            MChijingAll = IsFromBGEventAOD(trkIndexPrimAll);//check if the particle comes from hijing or from enhanced event
            if(MChijingAll) incls[1]=kHijing;
            else incls[1]=kElse;
            
            // Source
            if(MCMotherAll->GetPdgCode()==22) incls[2] = kGamma;
            if(MCMotherAll->GetPdgCode()==111) incls[2] = kPi0;
            if(MCMotherAll->GetPdgCode()==221) incls[2] = kEta;
            
            if(WeightsForEnhanced){//accept only Enhanced event using HIJING weights for MC closure test
                if(!MChijingAll) continue;
                
                
                // Weights
                Double_t weightsRaphaelle = -1.;
                Double_t ptmotherw = -1.;
                Int_t electronsource = GetElecSourceType(MCPhotoEleAll,mcArray,ptmotherw);
                if(fcentral){
                    if((electronsource==kPi0NoFeedDown) || (electronsource==kGPi0NoFeedDown)) {
                        if(!tiltup && !tiltdw)weightsRaphaelle = GetMCweightPi0(ptmotherw);
                        if(tiltup && !tiltdw)weightsRaphaelle = GetMCweightPi0tiltUp(ptmotherw);
                        if(!tiltup && tiltdw)weightsRaphaelle = GetMCweightPi0tiltDw(ptmotherw);
                        fIncl->Fill(incls,weightsRaphaelle);
                        SelectPhotonicElectronR(iTracks, track, incls[1], weightsRaphaelle, fMCmomAll, MCPhotoEleAll->GetPdgCode(), incls[2]);
                        // SelectPhotonicElectronR_ULSLS(iTracks, track, incls[1], weightsRaphaelle, MCPhotoEleAll->GetPdgCode(), incls[2]);
                    }
                    if((electronsource==kEtaNoFeedDown) || (electronsource==kGEtaNoFeedDown)) {
                        if(!tiltup && !tiltdw)weightsRaphaelle = GetMCweightEta(ptmotherw);
                        if(tiltup && !tiltdw)weightsRaphaelle = GetMCweightEtatiltUp(ptmotherw);
                        if(!tiltup && tiltdw)weightsRaphaelle = GetMCweightEtatiltDw(ptmotherw);
                        fIncl->Fill(incls,weightsRaphaelle);
                        SelectPhotonicElectronR(iTracks, track, incls[1], weightsRaphaelle, fMCmomAll, MCPhotoEleAll->GetPdgCode(), incls[2]);
                        //  SelectPhotonicElectronR_ULSLS(iTracks, track, incls[1], weightsRaphaelle, MCPhotoEleAll->GetPdgCode(), incls[2]);
                    }
                }
                if(fsemicentral){
                    if((electronsource==kPi0NoFeedDown) || (electronsource==kGPi0NoFeedDown)) {
                        if(!tiltup && !tiltdw)weightsRaphaelle = GetMCweightPi02040(ptmotherw);
                        if(tiltup && !tiltdw)weightsRaphaelle = GetMCweightPi0tiltUp2040(ptmotherw);
                        if(!tiltup && tiltdw)weightsRaphaelle = GetMCweightPi0tiltDw2040(ptmotherw);
                        fIncl->Fill(incls,weightsRaphaelle);
                        SelectPhotonicElectronR(iTracks, track, incls[1], weightsRaphaelle, fMCmomAll, MCPhotoEleAll->GetPdgCode(), incls[2]);
                        // SelectPhotonicElectronR_ULSLS(iTracks, track, incls[1], weightsRaphaelle, MCPhotoEleAll->GetPdgCode(), incls[2]);
                    }
                    if((electronsource==kEtaNoFeedDown) || (electronsource==kGEtaNoFeedDown)) {
                        if(!tiltup && !tiltdw)weightsRaphaelle = GetMCweightEta2040(ptmotherw);
                        if(tiltup && !tiltdw)weightsRaphaelle = GetMCweightEtatiltUp2040(ptmotherw);
                        if(!tiltup && tiltdw)weightsRaphaelle = GetMCweightEtatiltDw2040(ptmotherw);
                        fIncl->Fill(incls,weightsRaphaelle);
                        SelectPhotonicElectronR(iTracks, track, incls[1], weightsRaphaelle, fMCmomAll, MCPhotoEleAll->GetPdgCode(), incls[2]);
                        //  SelectPhotonicElectronR_ULSLS(iTracks, track, incls[1], weightsRaphaelle, MCPhotoEleAll->GetPdgCode(), incls[2]);
                    }
                }
                
            }//with whatever weights (or HIJING +Enh with cocktail weight or Enhanced only with HIJING weights)
            
            if(!WeightsForEnhanced){
                if(!MChijingAll) continue;
                // Weights == 1
                Double_t weightsRaphaelle = -1.;
                Double_t ptmotherw = -1.;
                Int_t electronsource = GetElecSourceType(MCPhotoEleAll,mcArray,ptmotherw);
                if((electronsource==kPi0NoFeedDown) || (electronsource==kGPi0NoFeedDown)) {
                    weightsRaphaelle = 1;
                    fIncl->Fill(incls,weightsRaphaelle);
                    SelectPhotonicElectronR(iTracks, track, incls[1], weightsRaphaelle, fMCmomAll, MCPhotoEleAll->GetPdgCode(), incls[2]);
                    //  SelectPhotonicElectronR_ULSLS(iTracks, track, incls[1], weightsRaphaelle, MCPhotoEleAll->GetPdgCode(), incls[2]);
                }
                if((electronsource==kEtaNoFeedDown) || (electronsource==kGEtaNoFeedDown)) {
                    weightsRaphaelle = 1;
                    fIncl->Fill(incls,weightsRaphaelle);
                    SelectPhotonicElectronR(iTracks, track, incls[1], weightsRaphaelle, fMCmomAll, MCPhotoEleAll->GetPdgCode(), incls[2]);
                    //  SelectPhotonicElectronR_ULSLS(iTracks, track, incls[1], weightsRaphaelle, MCPhotoEleAll->GetPdgCode(), incls[2]);
                }
                
            }
            
            
        }//mom
        
    } //track loop
    
    
    
    PostData(1, fOutputList);
}
//________________________________________________________________________
void AliAnalysisTaskHFEEfficiency::UserCreateOutputObjects()
{
    
    
    
    AliDebug(3, "Creating Output Objects");
    // Automatic determination of the analysis mode
    AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
        SetAODAnalysis();
    } else {
        SetESDAnalysis();
    }
    printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");
    
    //--------Initialize PID
    fPID->SetHasMCData(HasMCData());
    if(!fPID->GetNumberOfPIDdetectors())
    {
        fPID->AddDetector("ITS", 0);
        fPID->AddDetector("TOF", 1);
        fPID->AddDetector("TPC", 2);
    }
    
    fPID->SortDetectors();
    fPIDqa = new AliHFEpidQAmanager();
    fPIDqa->Initialize(fPID);
    
    //--------Initialize correction Framework and Cuts
    fCFM = new AliCFManager;
    const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
    fCFM->SetNStepParticle(kNcutSteps);
    for(Int_t istep = 0; istep < kNcutSteps; istep++)
        fCFM->SetParticleCutsList(istep, NULL);
    
    if(!fCuts){
        AliWarning("Cuts not available. Default cuts will be used");
        fCuts = new AliHFEcuts;
        fCuts->CreateStandardCuts();
    }
    
    if(IsAODanalysis()) fCuts->SetAOD();
    fCuts->Initialize(fCFM);
    
    fOutputList = new TList();
    fOutputList->SetOwner();
    fOutputList->Add(fPIDqa->MakeList("PIDQA"));
    
    fNoEvents = new TH1F("fNoEvents","",1,0,1) ;
    fOutputList->Add(fNoEvents);
    
    fMCPdgmom = new TH1F("fMCPdgmom", "MC pdg of electron mother; PDG code; count",1000,-0.5,999.5);
    fOutputList->Add(fMCPdgmom);
    
    fMCElecPtPhoto = new TH1F("fMCElecPtPhoto","MC photonic electron pt photonic; p_{T}; count",1000,0,50);
    fOutputList->Add(fMCElecPtPhoto);
    
    fMCElecPtIncl = new TH1F("fMCElecPtIncl","MC photonic electron pt inclusive; p_{T}; count",1000,0,50);
    fOutputList->Add(fMCElecPtIncl);
    
    fInclusiveElecPt = new TH1F("fInclusiveElecPt","Inclusive elec Pt; p_{T}; count",1000,0,50);
    fOutputList->Add(fInclusiveElecPt);
    
    fPhotoElecPt = new TH1F("fPhotoElecPt","Photonic elec Pt; p_{T}; count",1000,0,50);
    fOutputList->Add(fPhotoElecPt);
    
    fTruePhotoElecPt = new TH1F("fTruePhotoElecPt","True Photonic elec Pt; p_{T}; count",1000,0,50);
    fOutputList->Add(fTruePhotoElecPt);
    
    fSingleElecPt = new TH1F("fSingleElecPt","Single elec Pt; p_{T}; count",1000,0,50);
    fOutputList->Add(fSingleElecPt);
    
    fTrueSingleElecPt = new TH1F("fTrueSingleElecPt","True Single elec Pt; p_{T}; count",1000,0,50);
    fOutputList->Add(fTrueSingleElecPt);
    
    fElecNos = new TH1F("fElecNos","no of electrons",10,-0.5,9.5);
    fOutputList->Add(fElecNos);
    
    fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e); mass(GeV/c^2); counts;", 1000,0,1);
    fOutputList->Add(fInvmassLS);
    
    fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e); mass(GeV/c^2); counts;", 1000,0,1);
    fOutputList->Add(fInvmassULS);
    
    fInvmassULSNHFE = new TH1F("fInvmassULSNHFE", "Inv mass of NHFE (e,e) pair; mass(GeV/c^2); counts;", 1000,0,1);
    fOutputList->Add(fInvmassULSNHFE);
    
    fOpeningAngleLS = new TH1F("fOpeningAngleLS","Opening angle for LS pairs",100,0,1);
    fOutputList->Add(fOpeningAngleLS);
    
    fOpeningAngleULS = new TH1F("fOpeningAngleULS","Opening angle for ULS pairs",100,0,1);
    fOutputList->Add(fOpeningAngleULS);
    
    fTrackPtBefTrkCuts = new TH1F("fTrackPtBefTrkCuts","track pt before track cuts",1000,0,50);
    fOutputList->Add(fTrackPtBefTrkCuts);
    
    fTrackPtAftTrkCuts = new TH1F("fTrackPtAftTrkCuts","track pt after track cuts",1000,0,50);
    fOutputList->Add(fTrackPtAftTrkCuts);
    
    fTPCnsigma = new TH2F("fTPCnsigma", "TPC - n sigma",1000,0,50,400,-20,20);
    fOutputList->Add(fTPCnsigma);
    
    
    fITSnsigma = new TH2F("fITSnsigma", "ITS - n sigma",1000,0,50,400,-20,20);
    fOutputList->Add(fITSnsigma);
    
    fITSnsigmaElectron = new TH2F("fITSnsigmaElectron", "ITS - n sigma Elect",600,0,6,400,-20,20);
    fOutputList->Add(fITSnsigmaElectron);
    
    fTOFnsigma = new TH2F("fTOFnsigma", "TOF - n sigma",1000,0,50,400,-20,20);
    fOutputList->Add(fTOFnsigma);
    
    
    fTrkEovPBef = new TH2F("fTrkEovPBef","track E/p before HFE pid",1000,0,50,100,0,2);
    fOutputList->Add(fTrkEovPBef);
    
    fdEdxBef = new TH2F("fdEdxBef","track dEdx vs p before HFE pid",1000,0,50,150,0,150);
    fOutputList->Add(fdEdxBef);
    
    
    fCentralityPass = new TH1F("fCentralityPass", "Centrality Pass", 101, -1, 100);
    fOutputList->Add(fCentralityPass);
    
    fCentralityNoPass = new TH1F("fCentralityNoPass", "Centrality No Pass", 101, -1, 100);
    fOutputList->Add(fCentralityNoPass);
    
    fPhotoElecCandiPdgMom = new TH1F("fPhotoElecCandiPdgMom", "MC pdg of Photo elec candidate; PDG code; count",1000,-0.5,999.5);
    fOutputList->Add(fPhotoElecCandiPdgMom);
    
    fVZEROA = new TH1F("fVZEROA", "VZERO A Multiplicity", 1000, 0, 10000);
    fOutputList->Add(fVZEROA);
    
    fVZEROC = new TH1F("fVZEROC", "VZERO C Multiplicity", 1000, 0, 10000);
    fOutputList->Add(fVZEROC);
    
    fMCtagPhotoElecPtAll = new TH1F("fMCtagPhotoElecPtAll","All MC tagged Photonic elec Pt; p_{T}; count",1000,0,50);
    fMCtagPhotoElecPtAll->Sumw2();
    fOutputList->Add(fMCtagPhotoElecPtAll);
    
    fMCtagTruePhotoElecPtAll = new TH1F("fMCtagTruePhotoElecPtAll","All MC tagged True Photonic elec Pt; p_{T}; count",1000,0,50);
    fMCtagTruePhotoElecPtAll->Sumw2();
    fOutputList->Add(fMCtagTruePhotoElecPtAll);
    
    fMCtagGammaElecPtAll = new TH1F("fMCtagGammaElecPtAll","All MC tagged Gamma elec Pt; p_{T}; count",1000,0,50);
    fMCtagGammaElecPtAll->Sumw2();
    fOutputList->Add(fMCtagGammaElecPtAll);
    
    fMCtagPi0ElecPtAll = new TH1F("fMCtagPi0ElecPtAll","All MC tagged Pi0 elec Pt; p_{T}; count",1000,0,50);
    fMCtagPi0ElecPtAll->Sumw2();
    fOutputList->Add(fMCtagPi0ElecPtAll);
    
    fMCtagEtaElecPtAll = new TH1F("fMCtagEtaElecPtAll","All MC tagged Eta elec Pt; p_{T}; count",1000,0,50);
    fMCtagEtaElecPtAll->Sumw2();
    fOutputList->Add(fMCtagEtaElecPtAll);
    
    
    fTPCnsigmaAft = new TH2F("fTPCnsigmaAft", "TPC - n sigma after HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fTPCnsigmaAft);
    
    fTPCnsigmaVSptAft = new TH2F("fTPCnsigmaVSptAft", "TPC - n sigma after HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fTPCnsigmaVSptAft);
    
    fTPCnsigmaAftITSTOF = new TH2F("fTPCnsigmaAftITSTOF", "TPC - n sigma after HFE pid",600,0,6,400,-20,20);
    fOutputList->Add(fTPCnsigmaAftITSTOF);
    
    Rconv_pT = new TH2F("Rconv_pT", "Rconv_pT",500,0,50,100,0,30);
    fOutputList->Add(Rconv_pT);
    
    
    fakepT= new TH1F("fakepT", "fakepT",500,0,50);
    fOutputList->Add(fakepT);
    
    fakepTgraterFive= new TH1F("fakepTgraterFive", "fakepTgraterFive",500,0,50);
    fOutputList->Add(fakepTgraterFive);
    
    
    Int_t    binsmass[2] = { 100, 200};
    Double_t xminmass[2] = { 0.,  0};
    Double_t xmaxmass[2] = { 5., 1.};
    fSparseMassULS = new THnSparseF("fSparseMassULS", "pt:mass (GeV/c^{2})", 2, binsmass, xminmass, xmaxmass);
    fSparseMassULS->GetAxis(0)->SetTitle("pt (Gev/c)");
    fSparseMassULS->GetAxis(1)->SetTitle("mass");
    fOutputList->Add(fSparseMassULS);
    
    fSparseMassLS = new THnSparseF("fSparseMassLS", "pt:mass (GeV/c^{2})", 2, binsmass, xminmass, xmaxmass);
    fSparseMassLS->GetAxis(0)->SetTitle("pt (Gev/c)");
    fSparseMassLS->GetAxis(1)->SetTitle("mass");
    fOutputList->Add(fSparseMassLS);
    
    /*
     Int_t binsv1[9]={2000,2000,1000,200,1000,100, 2,2,2}; //trkPt, trkMCPt, p, dEdx, ClsE, E/p,fFlagPhotonicElec, fFlagMCPhotonicElec, fFlagMCPhotonicElecInvmass
     Double_t xminv1[9]={0,   0,  0,   0,    0,   0,  0,0,0};
     Double_t xmaxv1[]={100, 100,100, 200, 100, 2.5, 2,2,2};
     fSparseTrkMat = new THnSparseD ("Electrons","Electrons",9,binsv1,xminv1,xmaxv1);
     fOutputList->Add(fSparseTrkMat);
     */
    
    // QA weights
    Int_t nBinspdg = 3;
    Double_t minpdg = 0.;
    Double_t maxpdg = 3.;
    Double_t binLimpdg[nBinspdg+1];
    for(Int_t i=0; i<=nBinspdg; i++) binLimpdg[i]=(Double_t)minpdg + (maxpdg-minpdg)/nBinspdg*(Double_t)i ;
    
    Int_t nBinsg = 2;
    Double_t ming = 0.;
    Double_t maxg = 2.;
    Double_t binLimg[nBinsg+1];
    for(Int_t i=0; i<=nBinsg; i++) binLimg[i]=(Double_t)ming + (maxg-ming)/nBinsg*(Double_t)i ;
    
    Int_t nBinstype = 7;
    Double_t mintype = -1.;
    Double_t maxtype = 6.;
    Double_t binLimtype[nBinstype+1];
    for(Int_t i=0; i<=nBinstype; i++) binLimtype[i]=(Double_t)mintype + (maxtype-mintype)/nBinstype*(Double_t)i ;
    
    Int_t nBinspt = 500;
    Double_t minpt = 0.;
    Double_t maxpt = 50.;
    Double_t binLimpt[nBinspt+1];
    for(Int_t i=0; i<=nBinspt; i++) binLimpt[i]=(Double_t)minpt + (maxpt-minpt)/nBinspt*(Double_t)i ;
    
    const Int_t nDima=4;
    Int_t nBina[nDima] = {nBinspt,nBinspdg,nBinsg,nBinstype};
    fPi0EtaSpectra = new THnSparseF("Pi0EtaSpectra","Pi0EtaSpectra",nDima,nBina);
    fPi0EtaSpectra->SetBinEdges(0,binLimpt);
    fPi0EtaSpectra->SetBinEdges(1,binLimpdg);
    fPi0EtaSpectra->SetBinEdges(2,binLimg);
    fPi0EtaSpectra->SetBinEdges(3,binLimtype);
    fPi0EtaSpectra->Sumw2();
    fOutputList->Add(fPi0EtaSpectra);
    
    //    const Int_t BinsPtDefaultHF = 28;
    //    Double_t binLimPtDefaultHF[BinsPtDefaultHF+1] = {0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 4.25, 4.5, 4.75, 5., 6., 8., 10., 12., 14., 16., 18., 20.};
    
    const Int_t BinsPtDefaultHF = 12;
    Double_t binLimPtDefaultHF[BinsPtDefaultHF+1] = {0, 0.25, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2, 2.5, 3, 4, 5};
    
    
    Int_t nBinsPhi = 200;
    Double_t minPhi = 0;
    Double_t maxPhi = 6.5;
    Double_t binLimPhi[nBinsPhi+1];
    for(Int_t i=0; i<=nBinsPhi; i++) binLimPhi[i]=(Double_t)minPhi + (maxPhi-minPhi)/nBinsPhi*(Double_t)i ;
    
    
    const Int_t nDima3=4;
    Int_t nBina3[nDima3] = {BinsPtDefaultHF,nBinsg,nBinsg,nBinsPhi};
    fHFEDenominator = new THnSparseF("fHFEDenominator","fHFEDenominator;pt;hijing_Enh;beauty_or_Charm;",nDima3,nBina3);
    fHFEDenominator->SetBinEdges(0,binLimPtDefaultHF);
    fHFEDenominator->SetBinEdges(1,binLimg);
    fHFEDenominator->SetBinEdges(2,binLimg);
    fHFEDenominator->SetBinEdges(3,binLimPhi);
    fHFEDenominator->Sumw2();
    fOutputList->Add(fHFEDenominator);
    
    
    fHFENumeratorTRKCUTS_Filt = new THnSparseF("fHFENumeratorTRKCUTS_Filt","fHFENumeratorTRKCUTS_Filt;pt;hijing_Enh;beauty_or_Charm;",nDima3,nBina3);
    fHFENumeratorTRKCUTS_Filt->SetBinEdges(0,binLimPtDefaultHF);
    fHFENumeratorTRKCUTS_Filt->SetBinEdges(1,binLimg);
    fHFENumeratorTRKCUTS_Filt->SetBinEdges(2,binLimg);
    fHFENumeratorTRKCUTS_Filt->SetBinEdges(3,binLimPhi);
    fHFENumeratorTRKCUTS_Filt->Sumw2();
    fOutputList->Add(fHFENumeratorTRKCUTS_Filt);
    
    fHFENumeratorTRKCUTS_ITSTPC = new THnSparseF("fHFENumeratorTRKCUTS_ITSTPC","fHFENumeratorTRKCUTS_ITSTPC;pt;hijing_Enh;beauty_or_Charm;",nDima3,nBina3);
    fHFENumeratorTRKCUTS_ITSTPC->SetBinEdges(0,binLimPtDefaultHF);
    fHFENumeratorTRKCUTS_ITSTPC->SetBinEdges(1,binLimg);
    fHFENumeratorTRKCUTS_ITSTPC->SetBinEdges(2,binLimg);
    fHFENumeratorTRKCUTS_ITSTPC->SetBinEdges(3,binLimPhi);
    fHFENumeratorTRKCUTS_ITSTPC->Sumw2();
    fOutputList->Add(fHFENumeratorTRKCUTS_ITSTPC);
    
    
    //    fHFENumeratorTRKCUTS_ITStrkCut = new THnSparseF("fHFENumeratorTRKCUTS_ITStrkCut","fHFENumeratorTRKCUTS_ITStrkCut;pt;hijing_Enh;beauty_or_Charm;",nDima3,nBina3);
    //    fHFENumeratorTRKCUTS_ITStrkCut->SetBinEdges(0,binLimPtDefaultHF);
    //    fHFENumeratorTRKCUTS_ITStrkCut->SetBinEdges(1,binLimg);
    //    fHFENumeratorTRKCUTS_ITStrkCut->SetBinEdges(2,binLimg);
    //    fHFENumeratorTRKCUTS_ITStrkCut->SetBinEdges(3,binLimPhi);
    //    fHFENumeratorTRKCUTS_ITStrkCut->Sumw2();
    //    fOutputList->Add(fHFENumeratorTRKCUTS_ITStrkCut);
    //
    
    fHFENumeratorTRKCUTS = new THnSparseF("fHFENumeratorTRKCUTS","fHFENumeratorTRKCUTS;pt;hijing_Enh;beauty_or_Charm;",nDima3,nBina3);
    fHFENumeratorTRKCUTS->SetBinEdges(0,binLimPtDefaultHF);
    fHFENumeratorTRKCUTS->SetBinEdges(1,binLimg);
    fHFENumeratorTRKCUTS->SetBinEdges(2,binLimg);
    fHFENumeratorTRKCUTS->SetBinEdges(3,binLimPhi);
    fHFENumeratorTRKCUTS->Sumw2();
    fOutputList->Add(fHFENumeratorTRKCUTS);
    
    
    fHFENumeratorTOF = new THnSparseF("fHFENumeratorTOF","fHFENumeratorTOF;pt;hijing_Enh;beauty_or_Charm;",nDima3,nBina3);
    fHFENumeratorTOF->SetBinEdges(0,binLimPtDefaultHF);
    fHFENumeratorTOF->SetBinEdges(1,binLimg);
    fHFENumeratorTOF->SetBinEdges(2,binLimg);
    fHFENumeratorTOF->SetBinEdges(3,binLimPhi);
    fHFENumeratorTOF->Sumw2();
    fOutputList->Add(fHFENumeratorTOF);
    
    
    fHFENumeratorITS = new THnSparseF("fHFENumeratorITS","fHFENumeratorITS;pt;hijing_Enh;beauty_or_Charm;",nDima3,nBina3);
    fHFENumeratorITS->SetBinEdges(0,binLimPtDefaultHF);
    fHFENumeratorITS->SetBinEdges(1,binLimg);
    fHFENumeratorITS->SetBinEdges(2,binLimg);
    fHFENumeratorITS->SetBinEdges(3,binLimPhi);
    fHFENumeratorITS->Sumw2();
    fOutputList->Add(fHFENumeratorITS);
    
    fHFENumerator = new THnSparseF("fHFENumerator","fHFENumerator;pt;hijing_Enh;beauty_or_Charm;",nDima3,nBina3);
    fHFENumerator->SetBinEdges(0,binLimPtDefaultHF);
    fHFENumerator->SetBinEdges(1,binLimg);
    fHFENumerator->SetBinEdges(2,binLimg);
    fHFENumerator->SetBinEdges(3,binLimPhi);
    fHFENumerator->Sumw2();
    fOutputList->Add(fHFENumerator);
    
    
    
    Int_t nBinsInvMass = 100;
    Double_t minInvMass = 0.;
    Double_t maxInvMass = 1;
    Double_t binLimInvMass[nBinsInvMass+1];
    for(Int_t i=0; i<=nBinsInvMass; i++) binLimInvMass[i]=(Double_t)minInvMass + (maxInvMass-minInvMass)/nBinsInvMass*(Double_t)i ;
    
    Int_t nBinsSource = 6;
    Double_t minSource = 0.;
    Double_t maxSource = 6.;
    Double_t binLimSource[nBinsSource+1];
    for(Int_t i=0; i<=nBinsSource; i++) binLimSource[i]=(Double_t)minSource + (maxSource-minSource)/nBinsSource*(Double_t)i ;
    
    //    const Int_t kBinsPtDefault = 28;
    //    Double_t binLimPtDefault[kBinsPtDefault+1] = {0., 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 4.25, 4.5, 4.75, 5., 6., 8., 10., 12., 14., 16., 18., 20.};
    
    const Int_t kBinsPtDefault = 12;
    Double_t binLimPtDefault[kBinsPtDefault+1] = {0, 0.25, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2, 2.5, 3, 4, 5};
    
    const Int_t nDimb=3;
    Int_t nBinb[nDimb] = {kBinsPtDefault,nBinsg,nBinsSource};
    fIncl = new THnSparseF("Incl","Incl",nDimb,nBinb);
    fIncl->SetBinEdges(0,binLimPtDefault);
    fIncl->SetBinEdges(1,binLimg);
    fIncl->SetBinEdges(2,binLimSource);
    fIncl->Sumw2();
    fOutputList->Add(fIncl);
    
    
    const Int_t nDimc=4;
    Int_t nBinc[nDimc] = {kBinsPtDefault,nBinsg,nBinsSource,nBinsInvMass};
    fInvMULS = new THnSparseF("InvMULS","InvMULS",nDimc,nBinc);
    fInvMULS->SetBinEdges(0,binLimPtDefault);
    fInvMULS->SetBinEdges(1,binLimg);
    fInvMULS->SetBinEdges(2,binLimSource);
    fInvMULS->SetBinEdges(3,binLimInvMass);
    fInvMULS->Sumw2();
    fOutputList->Add(fInvMULS);
    
    fULS = new THnSparseF("fULS","fULS",nDimc,nBinc);
    fULS->SetBinEdges(0,binLimPtDefault);
    fULS->SetBinEdges(1,binLimg);
    fULS->SetBinEdges(2,binLimSource);
    fULS->SetBinEdges(3,binLimInvMass);
    fULS->Sumw2();
    fOutputList->Add(fULS);
    
    fLS = new THnSparseF("fLS","fLS",nDimc,nBinc);
    fLS->SetBinEdges(0,binLimPtDefault);
    fLS->SetBinEdges(1,binLimg);
    fLS->SetBinEdges(2,binLimSource);
    fLS->SetBinEdges(3,binLimInvMass);
    fLS->Sumw2();
    fOutputList->Add(fLS);
    
    
    //    fhHFStackBeauty = new TH1F("fhHFStackBeauty","fhHFStackBeauty",500,-5,5);
    //    fhHFStackBeauty->Sumw2();
    //    fOutputList->Add(fhHFStackBeauty);
    //    fhHFStackCharm = new TH1F("fhHFStackCharm","fhHFStackCharm",500,-5,5);
    //    fhHFStackCharm->Sumw2();
    //    fOutputList->Add(fhHFStackCharm);
    
    fNoEventsStackHFE = new TH1F("fNoEventsStackHFE","",1,0,1) ;
    fOutputList->Add(fNoEventsStackHFE);
    
    PostData(1,fOutputList);
}
//________________________________________________________________________
void AliAnalysisTaskHFEEfficiency::Terminate(Option_t *)
{
    // Info("Terminate");
    AliAnalysisTaskSE::Terminate();
}
//________________________________________________________________________
Bool_t AliAnalysisTaskHFEEfficiency::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
    // Check single track cuts for a given cut step
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    return kTRUE;
}
//_________________________________________
void AliAnalysisTaskHFEEfficiency::SelectPhotonicElectronR(Int_t itrack, AliVTrack *track, Int_t hijing, Double_t weight, Int_t motherindex, Int_t pdg, Int_t source)
{
    
    fTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fTrackCuts->SetRequireTPCRefit(kTRUE);
    fTrackCuts->SetRequireITSRefit(kTRUE);
    fTrackCuts->SetEtaRange(-0.9,0.9);
    fTrackCuts->SetRequireSigmaToVertex(kTRUE);
    fTrackCuts->SetMaxChi2PerClusterTPC(4);
    fTrackCuts->SetMinNClustersTPC(80);
    fTrackCuts->SetMaxDCAToVertexZ(3.2);
    fTrackCuts->SetMaxDCAToVertexXY(2.4);
    fTrackCuts->SetDCAToVertex2D(kTRUE);
    
    // load MC array
    AliMCEvent* mcEvent;
    AliMCEventHandler* eventHandler;
    AliAODMCHeader *mcHeader;
    TClonesArray* mcArray;
    
    if(IsAODanalysis()){
        mcArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        if(!mcArray){
            AliError("Array of MC particles not found");
            return;
        }
        
        mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!mcHeader) {
            AliError("Could not find MC Header in AOD");
            return;
        }
    }
    else{
        eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
        if (!eventHandler) {
            Printf("ERROR: Could not retrieve MC event handler");
            return;
        }
        
        mcEvent = eventHandler->MCEvent();
        if (!mcEvent) {
            Printf("ERROR: Could not retrieve MC event");
            return;
        }
        //AliStack * MCstack = mcEvent->Stack();
    }
    
    Double_t invms[4];
    for(Int_t jTracks = 0; jTracks<fVevent->GetNumberOfTracks(); jTracks++){
        AliVParticle* VtrackAsso = fVevent->GetTrack(jTracks);
        if (!VtrackAsso) {
            printf("ERROR: Could not receive track %d\n", jTracks);
            continue;
        }
        
        AliVTrack *trackAsso = dynamic_cast<AliVTrack*>(VtrackAsso);
        
        //track cuts applied
        Int_t pdgass = 0;
        if(IsAODanalysis()) {
            AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
            if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(atrackAsso->GetTPCNcls() < fAssoTPCCluster) continue;
            if(fAssoITSRefit){
                if(!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
            }
            if(!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
            AliAODMCParticle *MCass = (AliAODMCParticle*)mcArray->At(TMath::Abs(atrackAsso->GetLabel()));
            if(!(TMath::Abs(MCass->GetPdgCode())==11)) continue;
            pdgass=MCass->GetPdgCode();
            Int_t indexass = MCass->GetMother();
            if(TMath::Abs(indexass-motherindex)>0.8) continue;
            
        }
        else{
            AliESDtrack *etrackAsso = dynamic_cast<AliESDtrack*>(VtrackAsso);
            if(!fTrackCuts->AcceptTrack(etrackAsso)) continue;
        }
        
        if(jTracks==itrack) continue;
        
        Double_t  ptAsso=-999., nsigma=-999.0;
        Double_t mass=-999., width = -999;
        Bool_t fFlagULS=kFALSE;
        
        nsigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(trackAsso, AliPID::kElectron) : 1000;
        ptAsso = trackAsso->Pt();
        Int_t chargeAsso = trackAsso->Charge();
        Int_t charge = track->Charge();
        
        //---------------pt and track cuts-----------------
        if(ptAsso < fAssopTMin) continue;
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(nsigma < -3 || nsigma > 3) continue;
        
        
        //-------------------AliKFParticle-------------------
        Int_t PDGe1 = 11; Int_t PDGe2 = 11;
        if(charge>0) PDGe1 = -11;
        if(chargeAsso>0) PDGe2 = -11;
        
        if((pdgass*pdg)<0.) fFlagULS = kTRUE;
        
        AliKFParticle ge1(*track, PDGe1);
        AliKFParticle ge2(*trackAsso, PDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        recg.GetMass(mass,width);
        
        // Mother associated track
        invms[0]=track->Pt();
        invms[1]=hijing;
        invms[2]=source;
        invms[3]=mass;
        if(fFlagULS) fInvMULS->Fill(invms,weight);
        
        
    }
}
//___________________________________________________________________________________________________________________
void AliAnalysisTaskHFEEfficiency::SelectPhotonicElectronR_ULSLS(Int_t itrack, AliVTrack *track, Int_t hijing, Double_t weight, Int_t pdg, Int_t source)
{
    
    fTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fTrackCuts->SetRequireTPCRefit(kTRUE);
    fTrackCuts->SetRequireITSRefit(kTRUE);
    fTrackCuts->SetEtaRange(-0.9,0.9);
    fTrackCuts->SetRequireSigmaToVertex(kTRUE);
    fTrackCuts->SetMaxChi2PerClusterTPC(4);
    fTrackCuts->SetMinNClustersTPC(80);
    fTrackCuts->SetMaxDCAToVertexZ(3.2);
    fTrackCuts->SetMaxDCAToVertexXY(2.4);
    fTrackCuts->SetDCAToVertex2D(kTRUE);
    
    // load MC array
    AliMCEvent* mcEvent;
    AliMCEventHandler* eventHandler;
    AliAODMCHeader *mcHeader;
    TClonesArray* mcArray;
    
    if(IsAODanalysis()){
        mcArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        if(!mcArray){
            AliError("Array of MC particles not found");
            return;
        }
        
        mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!mcHeader) {
            AliError("Could not find MC Header in AOD");
            return;
        }
    }
    else{
        eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
        if (!eventHandler) {
            Printf("ERROR: Could not retrieve MC event handler");
            return;
        }
        
        mcEvent = eventHandler->MCEvent();
        if (!mcEvent) {
            Printf("ERROR: Could not retrieve MC event");
            return;
        }
        //AliStack * MCstack = mcEvent->Stack();
    }
    
    Double_t invmsULSLS[4];
    for(Int_t jTracks = 0; jTracks<fVevent->GetNumberOfTracks(); jTracks++){
        AliVParticle* VtrackAsso = fVevent->GetTrack(jTracks);
        if (!VtrackAsso) {
            printf("ERROR: Could not receive track %d\n", jTracks);
            continue;
        }
        
        AliVTrack *trackAsso = dynamic_cast<AliVTrack*>(VtrackAsso);
        
        //track cuts applied
        Int_t pdgass = 0;
        if(IsAODanalysis()) {
            AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
            if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(atrackAsso->GetTPCNcls() < fAssoTPCCluster) continue;
            if(fAssoITSRefit){
                if(!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
            }
            if(!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
            AliAODMCParticle *MCass = (AliAODMCParticle*)mcArray->At(TMath::Abs(atrackAsso->GetLabel()));
            if(!(TMath::Abs(MCass->GetPdgCode())==11)) continue;
            pdgass=MCass->GetPdgCode();
            // Int_t indexass = MCass->GetMother();
            // if(TMath::Abs(indexass-motherindex)>0.8) continue;
            
        }
        else{
            AliESDtrack *etrackAsso = dynamic_cast<AliESDtrack*>(VtrackAsso);
            if(!fTrackCuts->AcceptTrack(etrackAsso)) continue;
        }
        
        if(jTracks==itrack) continue;
        
        Double_t  ptAsso=-999., nsigma=-999.0;
        Double_t mass=-999., width = -999;
        Bool_t fFlagULS=kFALSE;
        Bool_t fFlagLS=kFALSE;
        
        nsigma = fPID->GetPIDResponse() ? fPID->GetPIDResponse()->NumberOfSigmasTPC(trackAsso, AliPID::kElectron) : 1000;
        ptAsso = trackAsso->Pt();
        Int_t chargeAsso = trackAsso->Charge();
        Int_t charge = track->Charge();
        
        //---------------pt and track cuts-----------------
        if(ptAsso < fAssopTMin) continue;
        if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
        if(nsigma < -3 || nsigma > 3) continue;
        
        
        //-------------------AliKFParticle-------------------
        Int_t PDGe1 = 11; Int_t PDGe2 = 11;
        if(charge>0) PDGe1 = -11;
        if(chargeAsso>0) PDGe2 = -11;
        
        if((pdgass*pdg)<0.) fFlagULS = kTRUE;
        if((pdgass*pdg)>0.) fFlagLS = kTRUE;
        
        AliKFParticle ge1(*track, PDGe1);
        AliKFParticle ge2(*trackAsso, PDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        recg.GetMass(mass,width);
        
        // Mother associated track
        invmsULSLS[0]=track->Pt();
        invmsULSLS[1]=hijing;
        invmsULSLS[2]=source;
        invmsULSLS[3]=mass;
        if(fFlagULS) fULS->Fill(invmsULSLS,weight);
        if(fFlagLS) fLS->Fill(invmsULSLS,weight);
        
        
    }
}

//_________________________________________
void AliAnalysisTaskHFEEfficiency::CheckCentrality(AliVEvent* event, Bool_t &centralitypass)
{
    // Check if event is within the set centrality range. Falls back to V0 centrality determination if no method is set
    if (!fkCentralityMethod) AliFatal("No centrality method set! FATAL ERROR!");
    fCentrality = event->GetCentrality()->GetCentralityPercentile(fkCentralityMethod);
    //  cout << "--------------Centrality evaluated-------------------------"<<endl;
    
    if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax))
    {
        fCentralityNoPass->Fill(fCentrality);
        //   cout << "--------------Fill no pass-------------------------"<<endl;
        centralitypass = kFALSE;
    }else
    {
        fCentralityPass->Fill(fCentrality);
        //   cout << "--------------Fill pass-------------------------"<<endl;
        centralitypass = kTRUE;
    }
    
}
//_____________________________________________________________________________
void AliAnalysisTaskHFEEfficiency::SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod)
{
    // Set a centrality range ]min, max] and define the method to use for centrality selection
    fCentralityMin = CentralityMin;
    fCentralityMax = CentralityMax;
    fkCentralityMethod = CentralityMethod;
}
//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskHFEEfficiency::PlotVZeroMultiplcities(const T* event) const
{
    // QA multiplicity plots
    fVZEROA->Fill(event->GetVZEROData()->GetMTotV0A());
    fVZEROC->Fill(event->GetVZEROData()->GetMTotV0C());
}
//___________________________________________________________________________
Bool_t AliAnalysisTaskHFEEfficiency::IsFromBGEventAOD(Int_t Index)
{
    //Check if the particle is from Hijing or Enhanced event
    AliAODMCHeader *mcHeader;
    Int_t fNBG =-1;
    
    mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
        AliError("Could not find MC Header in AOD");
        return (0);
    }
    
    TList *List = mcHeader->GetCocktailHeaders();
    AliGenHijingEventHeader* hijingH = dynamic_cast<AliGenHijingEventHeader*>(List->FindObject("Hijing"));
    if (!hijingH){
        AliError("no GenHijing header");
        return (0);
    }
    fNBG = hijingH->NProduced();
    
    return (Index < fNBG);
}
//_________________________________________

Int_t AliAnalysisTaskHFEEfficiency::GetPrimary(Int_t id, TClonesArray *mcArray){
    
    // Return number of primary that has generated track
    int current, parent;
    parent=id;
    while (1) {
        current=parent;
        AliAODMCParticle *Part = (AliAODMCParticle*)mcArray->At(current);
        parent=Part->GetMother();
        //  cout << "GetPartArr momid :"  << parent << endl;
        if(parent<0) return current;
    }
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GiveHFWeight(Double_t electronpt)
{
    double weight = 1.0;
    weight =    (5.38354e-01)/pow((exp(-(1.56041e-01)*electronpt-(3.74084e+00)*electronpt*electronpt)+(electronpt/(5.88171e-01))),(2.40354e+00)); //HFeWeights
    return weight;
}
//_________________________________________
Int_t AliAnalysisTaskHFEEfficiency::GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray){
    //
    // Return what type of pi0, eta it is
    //
    
    
    // IsPrimary
    Bool_t primMC = pi0eta->IsPrimary();
    if(!primMC) return kNoIsPrimary;
    
    // Mother
    Int_t motherlabel = pi0eta->GetMother();
    if(motherlabel<0) return kNoMother;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        //    if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113) return kLightMesons;
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323) return kLightMesons;
        
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
        return kNoFeedDown;
        
    }
}
//_________________________________________

Int_t AliAnalysisTaskHFEEfficiency::GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm){
    //
    // Return what type of gammax it is
    //
    
    
    
    // Mother
    Int_t motherlabel = electron->GetMother();
    if(motherlabel<0) return kNoMotherE;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        ptm=mother->Pt();
        if(motherpdg == 111) {
            Int_t typepi0eta = GetPi0EtaType(mother,mcArray);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) 	return kPi0NoFeedDown;
        }
        if(motherpdg == 221) {
            Int_t typepi0eta = GetPi0EtaType(mother,mcArray);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kEtaNoFeedDown;
        }
        if(motherpdg == 22) {
            Int_t gmotherlabel = mother->GetMother();
            if(gmotherlabel<0) return kDirectGamma;
            else {
                AliAODMCParticle *gmother = (AliAODMCParticle*)mcArray->At(gmotherlabel);
                ptm=gmother->Pt();
                Int_t gmotherpdg = TMath::Abs(gmother->GetPdgCode());
                if(gmotherpdg == 111) {
                    Int_t typepi0eta = GetPi0EtaType(gmother,mcArray);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                }
                if(gmotherpdg == 221) {
                    Int_t typepi0eta = GetPi0EtaType(gmother,mcArray);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                }
                if(gmotherpdg == 22) {
                    Int_t ggmotherlabel = gmother->GetMother();
                    if(ggmotherlabel<0) return kDirectGamma;
                    else {
                        AliAODMCParticle *ggmother = (AliAODMCParticle*)mcArray->At(ggmotherlabel);
                        ptm=ggmother->Pt();
                        Int_t ggmotherpdg = TMath::Abs(ggmother->GetPdgCode());
                        if(ggmotherpdg == 111) {
                            Int_t typepi0eta = GetPi0EtaType(ggmother,mcArray);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                        }
                        if(ggmotherpdg == 221) {
                            Int_t typepi0eta = GetPi0EtaType(ggmother,mcArray);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                        }
                    }
                }
            }
        }
    }
    
    return kOthersE;
    
}
//================================================================================================
//_________________________________________
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightPi0(Double_t mcPi0pT) // make the function in first iteration DONE (HERE 0-10%)
{
    double weight = 1.0;
    //  cout << "pi0 weigth in weigth cal : " << mcPi0pT <<endl;
    
    if(mcPi0pT <=5){
        weight =    (4.69681e-01)/pow((exp(-(1.23103e+00)*mcPi0pT-(-1.26331e-01)*mcPi0pT*mcPi0pT)+(mcPi0pT/(3.69697e+00))),(2.54998e+00)); //Met 1 From EMCal Task
    }
    if(mcPi0pT >5){
        weight =    (1.19426e+00)/pow((exp(-(-1.02710e+00)*mcPi0pT-(7.12190e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(9.80915e-01))),(5.34179e-01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal pi0 : " << weight << endl;
    return weight;
}
//_________________________________________
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightPi0tiltUp(Double_t mcPi0pT) // make the function in first iteration DONE (HERE 0-10%)
{
    double weight = 1.0;
    //  cout << "pi0 weigth in weigth cal : " << mcPi0pT <<endl;
    if(mcPi0pT <=5){
        weight =    (4.87817e-01)/pow((exp(-(1.22741e+00)*mcPi0pT-(-1.36457e-01)*mcPi0pT*mcPi0pT)+(mcPi0pT/(3.74946e+00))),(2.58713e+00)); //Met 1 From EMCal Task
    }
    if(mcPi0pT >5){
        weight =    (2.19686e+00)/pow((exp(-(-8.34637e-01)*mcPi0pT-(5.85892e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(1.13653e+00))),(8.18447e-01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal pi0 : " << weight << endl;
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightPi0tiltDw(Double_t mcPi0pT) // make the function in first iteration DONE (HERE 0-10%)
{
    double weight = 1.0;
    //  cout << "pi0 weigth in weigth cal : " << mcPi0pT <<endl;
    if(mcPi0pT <=5){
        weight =    (4.37575e-01)/pow((exp(-(1.25745e+00)*mcPi0pT-(-1.51273e-01)*mcPi0pT*mcPi0pT)+(mcPi0pT/(3.78738e+00))),(2.52279e+00)); //Met 1 From EMCal Task
    }
    if(mcPi0pT >5){
        weight =    (8.59089e-01)/pow((exp(-(-1.09749e+00)*mcPi0pT-(7.66834e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(1.30487e+00))),(4.36030e-01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal pi0 : " << weight << endl;
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightEta(Double_t mcEtapT) // make the function in first iteration DONE
{
    double weight = 1.0;
    // cout << "eta weigth in weigth cal : " << mcEtapT <<endl;
    
    if(mcEtapT <=5.5){
        weight =    (4.01348e-01)/pow((exp(-(9.64548e-01)*mcEtapT-(7.30011e-02)*mcEtapT*mcEtapT)+(mcEtapT/(2.92423e+00))),(2.32353e+00)); //Met 1 From EMCal Task
    }
    if(mcEtapT >5.5){
        weight =    (1.42169e-01)/pow((exp(-(-8.71455e+01)*mcEtapT-(7.94280e+00)*mcEtapT*mcEtapT)+(mcEtapT/(7.36452e+01))),(2.55991e-03)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal eta : " << weight << endl;
    
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightEtatiltUp(Double_t mcEtapT) // make the function in first iteration DONE
{
    double weight = 1.0;
    // cout << "eta weigth in weigth cal : " << mcEtapT <<endl;
    
    if(mcEtapT <=5.5){
        weight =    (4.01348e-01)/pow((exp(-(9.64548e-01)*mcEtapT-(7.30011e-02)*mcEtapT*mcEtapT)+(mcEtapT/(2.92423e+00))),(2.32353e+00)); //Met 1 From EMCal Task
    }
    if(mcEtapT >5.5){
        weight =    (1.42169e-01)/pow((exp(-(-8.71455e+01)*mcEtapT-(7.94280e+00)*mcEtapT*mcEtapT)+(mcEtapT/(7.36452e+01))),(2.55991e-03)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal eta : " << weight << endl;
    
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightEtatiltDw(Double_t mcEtapT) // make the function in first iteration DONE
{
    double weight = 1.0;
    // cout << "eta weigth in weigth cal : " << mcEtapT <<endl;
    
    if(mcEtapT <=5.5){
        weight =    (4.01348e-01)/pow((exp(-(9.64548e-01)*mcEtapT-(7.30011e-02)*mcEtapT*mcEtapT)+(mcEtapT/(2.92423e+00))),(2.32353e+00)); //Met 1 From EMCal Task
    }
    if(mcEtapT >5.5){
        weight =    (1.42169e-01)/pow((exp(-(-8.71455e+01)*mcEtapT-(7.94280e+00)*mcEtapT*mcEtapT)+(mcEtapT/(7.36452e+01))),(2.55991e-03)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal eta : " << weight << endl;
    
    return weight;
}
//_________________________________________

Bool_t AliAnalysisTaskHFEEfficiency::FillNumerator(TClonesArray *mcArray, AliVTrack *track, Double_t ForHFEReco[])
{
    //============Numerator HFE efficiency===============
    Int_t HFElabel = TMath::Abs(track->GetLabel());
    AliAODMCParticle *MCHFE = (AliAODMCParticle*)mcArray->At(HFElabel);
    
    Int_t from_D_or_B_Reco;
    Bool_t ShouldIfillevenHere = kFALSE;
    
    if(TMath::Abs(MCHFE->GetPdgCode())==11){
        
        Int_t MChijingHFE=kElse;
        Int_t PrimHFE = GetPrimary(HFElabel, mcArray);
        AliAODMCParticle *AODMCtrackPrimHFE = (AliAODMCParticle*)mcArray->At(PrimHFE);
        Int_t trkIndexPrimHFE= AODMCtrackPrimHFE->GetLabel();//gives index of the particle in original MCparticle array
        if(IsFromBGEventAOD(trkIndexPrimHFE)) MChijingHFE = kHijing;//check if the particle comes from hijing or from enhanced event
        
        
        Int_t fMCmom_HFEReco = MCHFE->GetMother();
        AliAODMCParticle *MCMother_HFEReco = (AliAODMCParticle*)mcArray->At(fMCmom_HFEReco);
        Int_t motherpdgforHFE_Reco = MCMother_HFEReco->GetPdgCode();
        
        if ( (int(TMath::Abs(motherpdgforHFE_Reco)/100.)%10) == 5 || (int(TMath::Abs(motherpdgforHFE_Reco)/1000.)%10) == 5 ){
            from_D_or_B_Reco = kelectronBeauty;
            ForHFEReco[0] = track->Pt();
            ForHFEReco[1] = MChijingHFE;
            ForHFEReco[2] = from_D_or_B_Reco;
            ForHFEReco[3] = track->Phi();
            ShouldIfillevenHere = kTRUE;
        }
        
        if ( (int(TMath::Abs(motherpdgforHFE_Reco)/100.)%10) == 4 || (int(TMath::Abs(motherpdgforHFE_Reco)/1000.)%10) == 4 ){
            from_D_or_B_Reco = kelectronCharm;
            ForHFEReco[0] = track->Pt();
            ForHFEReco[1] = MChijingHFE;
            ForHFEReco[2] = from_D_or_B_Reco;
            ForHFEReco[3] = track->Phi();
            ShouldIfillevenHere = kTRUE;
        }
    }
    
    if(ShouldIfillevenHere) return kTRUE;
    else return kFALSE;
    
}
//====================================================


//_____________________________________________________________________________
void AliAnalysisTaskHFEEfficiency::SetIDCuts(Double_t minTOFnSigma, Double_t maxTOFnSigma, Double_t minITSnsigmalowpt, Double_t maxITSnsigmalowpt,  Double_t minITSnsigmahighpt, Double_t maxITSnsigmahighpt, Double_t minTPCnsigmalowpt, Double_t maxTPCnsigmalowpt,  Double_t minTPCnsigmahighpt, Double_t maxTPCnsigmahighpt)
{
    //Set PID cuts
    fminTOFnSigma = minTOFnSigma;
    fmaxTOFnSigma = maxTOFnSigma;
    fminITSnsigmaLowpT = minITSnsigmalowpt;
    fmaxITSnsigmaLowpT = maxITSnsigmalowpt;
    fminITSnsigmaHighpT = minITSnsigmahighpt;
    fmaxITSnsigmaHighpT = maxITSnsigmahighpt;
    fminTPCnsigmaLowpT = minTPCnsigmalowpt;
    fmaxTPCnsigmaLowpT = maxTPCnsigmalowpt;
    fminTPCnsigmaHighpT = minTPCnsigmahighpt;
    fmaxTPCnsigmaHighpT = maxTPCnsigmahighpt;
    
}
//_____________________________________________________________________________

//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightPi02040(Double_t mcPi0pT) // make the function in first iteration DONE (HERE 20-40%)
{
    double weight = 1.0;
    //  cout << "pi0 weigth in weigth cal : " << mcPi0pT <<endl;
    
    if(mcPi0pT <=4){
        weight =    (8.52613e-01)/pow((exp(-(8.40269e-01)*mcPi0pT-(-4.34674e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(3.31267e+00))),(2.56504e+00)); //Met 1 From EMCal Task
    }
    if(mcPi0pT >4){
        weight =    (1.58170e+02)/pow((exp(-(-9.51650e-02)*mcPi0pT-(1.43834e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(1.01125e+01))),(1.29361e+01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal pi0 : " << weight << endl;
    return weight;
}
//_________________________________________
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightPi0tiltUp2040(Double_t mcPi0pT) // make the function in first iteration DONE (HERE 20-40%)
{
    double weight = 1.0;
    //  cout << "pi0 weigth in weigth cal : " << mcPi0pT <<endl;
    
    if(mcPi0pT <=4){
        weight =    (8.52613e-01)/pow((exp(-(8.40269e-01)*mcPi0pT-(-4.34674e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(3.31267e+00))),(2.56504e+00)); //Met 1 From EMCal Task
    }
    if(mcPi0pT >4){
        weight =    (1.58170e+02)/pow((exp(-(-9.51650e-02)*mcPi0pT-(1.43834e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(1.01125e+01))),(1.29361e+01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal pi0 : " << weight << endl;
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightPi0tiltDw2040(Double_t mcPi0pT) // make the function in first iteration DONE (HERE 20-40%)
{
    double weight = 1.0;
    //  cout << "pi0 weigth in weigth cal : " << mcPi0pT <<endl;
    
    if(mcPi0pT <=4){
        weight =    (8.52613e-01)/pow((exp(-(8.40269e-01)*mcPi0pT-(-4.34674e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(3.31267e+00))),(2.56504e+00)); //Met 1 From EMCal Task
    }
    if(mcPi0pT >4){
        weight =    (1.58170e+02)/pow((exp(-(-9.51650e-02)*mcPi0pT-(1.43834e-02)*mcPi0pT*mcPi0pT)+(mcPi0pT/(1.01125e+01))),(1.29361e+01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal pi0 : " << weight << endl;
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightEta2040(Double_t mcEtapT) // make the function in first iteration DONE
{
    double weight = 1.0;
    // cout << "eta weigth in weigth cal : " << mcEtapT <<endl;
    
    if(mcEtapT <=4){
        weight =    (8.38897e-01)/pow((exp(-(2.69050e-02)*mcEtapT-(8.25401e-01)*mcEtapT*mcEtapT)+(mcEtapT/(1.92761e+00))),(1.48011e+00)); //Met 1 From EMCal Task
    }
    if(mcEtapT >4){
        weight =    (8.83196e+01)/pow((exp(-(-1.24373e-02)*mcEtapT-(6.69093e-03)*mcEtapT*mcEtapT)+(mcEtapT/(1.57278e+01))),(3.14751e+01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal eta : " << weight << endl;
    
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightEtatiltUp2040(Double_t mcEtapT) // make the function in first iteration DONE
{
    double weight = 1.0;
    // cout << "eta weigth in weigth cal : " << mcEtapT <<endl;
    
    if(mcEtapT <=4){
        weight =    (8.38897e-01)/pow((exp(-(2.69050e-02)*mcEtapT-(8.25401e-01)*mcEtapT*mcEtapT)+(mcEtapT/(1.92761e+00))),(1.48011e+00)); //Met 1 From EMCal Task
    }
    if(mcEtapT >4){
        weight =    (8.83196e+01)/pow((exp(-(-1.24373e-02)*mcEtapT-(6.69093e-03)*mcEtapT*mcEtapT)+(mcEtapT/(1.57278e+01))),(3.14751e+01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal eta : " << weight << endl;
    
    return weight;
}
//_________________________________________
Double_t AliAnalysisTaskHFEEfficiency::GetMCweightEtatiltDw2040(Double_t mcEtapT) // make the function in first iteration DONE
{
    double weight = 1.0;
    // cout << "eta weigth in weigth cal : " << mcEtapT <<endl;
    
    if(mcEtapT <=4){
        weight =    (8.38897e-01)/pow((exp(-(2.69050e-02)*mcEtapT-(8.25401e-01)*mcEtapT*mcEtapT)+(mcEtapT/(1.92761e+00))),(1.48011e+00)); //Met 1 From EMCal Task
    }
    if(mcEtapT >4){
        weight =    (8.83196e+01)/pow((exp(-(-1.24373e-02)*mcEtapT-(6.69093e-03)*mcEtapT*mcEtapT)+(mcEtapT/(1.57278e+01))),(3.14751e+01)); //Met 1 From EMCal Task
    }
    //  cout << "weight in cal eta : " << weight << endl;
    
    return weight;
}
//_________________________________________

//---------------------------------------------------------------------------
void AliAnalysisTaskHFEEfficiency::SetEtaRange(Double_t etaminimum, Double_t etamaximum)
{
    //Set PID cuts
    fmineta = etaminimum;
    fmaxeta = etamaximum;
}
//_____________________________________________________________________________
