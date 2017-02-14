#ifndef AliAnalysisTaskDiJetCorrelationsAllb2b_H
#define AliAnalysisTaskDiJetCorrelationsAllb2b_H

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>
#include <THn.h>
#include <iosfwd>

#include "AliAnalysisTaskSE.h"
#include "AliEventPoolManager.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliAODTrack.h"

class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class AliAODEvent;
class AliVParticle;
class TObjArray;
class AliEventPoolManager;
class AliESDEvent;
class AliESDtrackCuts;
class AliMultiplicity;
class AliAODHandler;
class AliAODInputHandler;
class AliMCParticle;
class TDatabasePDG;
class TParticlePDG;
class AliMCEvent;
class AliStack;
class AliInputEventHandler;
class AliAnalysisTaskSE;


class AliAnalysisTaskDiJetCorrelationsAllb2b : public AliAnalysisTaskSE
{
    
public:
    
    AliAnalysisTaskDiJetCorrelationsAllb2b();
    AliAnalysisTaskDiJetCorrelationsAllb2b(const char *name);
    virtual ~AliAnalysisTaskDiJetCorrelationsAllb2b();
    
    virtual void     UserCreateOutputObjects();
    //virtual void   SetInputCorrection();
    //virtual void   Init();
    //virtual void   LocalInit() {Init();}
    virtual void    UserExec(Option_t *option);
    virtual void    Terminate(Option_t *);
    
    virtual void    SetSEorME(Bool_t flag) {fMixedEvent = flag; }
    virtual void    SetMESettings(Int_t mevt, Int_t mtrk, Int_t nminMix){fMEMaxPoolEvent = mevt, fMEMinTracks = mtrk,  fMEMinEventToMix = nminMix;}
    virtual void    SetSystem(Bool_t system) { fSetSystemValue = system; }
    virtual void    SetTrigger1PTValue(Double_t pTmin1, Double_t pTmax1){fTrigger1pTLowThr = pTmin1, fTrigger1pTHighThr = pTmax1;}
    virtual void    SetTrigger2PTValue(Double_t pTmin2, Double_t pTmax2){fTrigger2pTLowThr = pTmin2, fTrigger2pTHighThr = pTmax2;}
    virtual void    SetCentralityRange(Double_t minCent, Double_t maxCent){fMinCentrality = minCent, fMaxCentrality = maxCent;}
    virtual void    SetDataType(Bool_t DataOrPart){fRecoOrMontecarlo = DataOrPart;}
    virtual void    SetFilterBit(Bool_t filterbit){fSetFilterBit=filterbit;}//
    virtual void    SetFilterType(Int_t bittype){fbit= bittype;}
    virtual void    SetVarCentBin(Bool_t varCbin=kTRUE){fuseVarCentBins= varCbin;}
    virtual void    SetVarPtBin(Bool_t varPtbin=kTRUE){fuseVarPtBins= varPtbin;}
    
    
    void    SetCorr2plus1or1plus1(Bool_t twoplus1){ftwoplus1 = twoplus1;}
    void    SetDiJetAlphaAngle(Double_t alpha){fAlpha = alpha;}
    void    SetBkgSE(Bool_t BkgSE){fBkgSE = BkgSE;}
    void    SetBkgSEBothSide(Bool_t BkgSEbothside){fBkgSEBothSide = BkgSEbothside;}
    void    SetEfficiencyWeightMap(THnF* hEff){fThnEff = hEff;}
    void    SetResonanceCut(Bool_t resCut){fCutResonances = resCut;}
    void    SetConversionCut(Bool_t conversionCut){fCutConversions = conversionCut;}
    void    SetTwoTrackEfficiencyCut(Bool_t TTRcut){ftwoTrackEfficiencyCut = TTRcut;}
    void    SetEqualT1T2Demand(Bool_t T1T2eQ){fEqualT1T2 = T1T2eQ;}
    void    SetFinerBinsME(Bool_t MEfinebins){fFineBinsME = MEfinebins;}
    
    
private:
    AliAnalysisTaskDiJetCorrelationsAllb2b(const AliAnalysisTaskDiJetCorrelationsAllb2b &source);
    AliAnalysisTaskDiJetCorrelationsAllb2b& operator=(const AliAnalysisTaskDiJetCorrelationsAllb2b& source);
    
    
    void DefineHistoNames();
    //Double_t GetTrackbyTrackEffValue(AliAODTrack* track, Double_t CentrOrMult, TH3F *h);
    
    Bool_t ConversionResonanceCut(Double_t refmaxpT, Double_t phiMaxpT, Double_t etaMaxpT, Double_t Charge, AliAODTrack* AodTracks, TH2F*fControlConvResT, TH1F* fHistTCorrTrack);
    Bool_t TwoTrackEfficiencyCut(Double_t refmaxpT, Double_t phiMaxpT, Double_t etaMaxpT, Double_t Charge, AliAODTrack* AodTracks, Float_t bSigntmp);
    
    TObjArray* CloneAndReduceTrackList(TObjArray* tracks);
    
    Double_t GetTrackWeight(Double_t eta, Double_t pt, Double_t CentrOrMult, Double_t zVertex);
    
    
    //______________________________|  DPhi Star
    Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign){
        
        //calculate dphistar for two track efficiency cut!
        Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
        static const Double_t kPi = TMath::Pi();
        if (dphistar > kPi)
            dphistar = kPi*2 - dphistar;
        if (dphistar < -kPi)
            dphistar = -kPi*2 - dphistar;
        //if (dphistar > kPi)
        //  dphistar = kPi * 2 - dphistar;
        return dphistar;
    }
    
    //______________________________|  Mass Squared
    Float_t GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2){
        //calculate inv ass squared
        Float_t tantheta1 = 1e10;
        if (eta1 < -1e-10 || eta1 > 1e-10){
            Float_t expTmp = TMath::Exp(-eta1);
            tantheta1 = 2.0 * expTmp/(1.0 -expTmp*expTmp);
        }
        
        Float_t tantheta2 = 1e10;
        if (eta2 < -1e-10 || eta2 > 1e-10){
            Float_t expTmp = TMath::Exp(-eta2);
            tantheta2 = 2* expTmp/(1.0 - expTmp*expTmp);
        }
        
        Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
        Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
        Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * (TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2) ) );
        return mass2;
    }
    
    //______________________________|  Mass Cheap
    Float_t GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2){
        
        // calculate inv mass squared approximately
        Float_t tantheta1 = 1e10;
        if (eta1 < -1e-10 || eta1 > 1e-10){
            Float_t expTmp = 1.0-eta1+eta1*eta1/2-eta1*eta1*eta1/6+eta1*eta1*eta1*eta1/24;
            tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
        }
        
        Float_t tantheta2 = 1e10;
        if (eta2 < -1e-10 || eta2 > 1e-10){
            Float_t expTmp = 1.0-eta2+eta2*eta2/2-eta2*eta2*eta2/6+eta2*eta2*eta2*eta2/24;
            tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
        }
        
        Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
        Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
        
        // fold onto 0...pi
        Float_t deltaPhi = TMath::Abs(phi1 - phi2);
        while (deltaPhi > TMath::TwoPi())
            deltaPhi -= TMath::TwoPi();
        if (deltaPhi > TMath::Pi())
            deltaPhi = TMath::TwoPi() - deltaPhi;
        
        Float_t cosDeltaPhi = 0;
        if (deltaPhi < TMath::Pi()/3)
            cosDeltaPhi = 1.0 - deltaPhi*deltaPhi/2 + deltaPhi*deltaPhi*deltaPhi*deltaPhi/24;
        else if (deltaPhi < 2*TMath::Pi()/3)
            cosDeltaPhi = -(deltaPhi - TMath::Pi()/2) + 1.0/6 * TMath::Power((deltaPhi - TMath::Pi()/2), 3);
        else
            cosDeltaPhi = -1.0 + 1.0/2.0*(deltaPhi - TMath::Pi())*(deltaPhi - TMath::Pi()) - 1.0/24.0 * TMath::Power(deltaPhi - TMath::Pi(), 4);
        
        Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( cosDeltaPhi + 1.0 / tantheta1 / tantheta2 ) ) );
        //   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));
        return mass2;
    }
    
    //______________________________| Phi Settings
    Double_t AssignCorrectPhiRange(Double_t phi){
        if (phi > 1.5 * TMath::Pi())phi -= TMath::TwoPi();
        if (phi < -0.5 * TMath::Pi())phi += TMath::TwoPi();
        Double_t phiClone = 0.;
        phiClone = phi;
        return phiClone;
    }
    
    //______________________________| Mixed Event Pool
    Bool_t DefineMixedEventPoolPbPb(){
        Int_t  NofCentBins  = 12;
        Double_t MBins[]={0., 1., 2., 3., 4., 5., 7.5, 10., 20., 30., 40., 50., 100.1};
        Double_t * CentrORMultBins = MBins;
        
        Int_t NofZVrtxBins  = 10;
        Double_t ZBins[]={-10.0, -8.0, -6.0, -4.0, -2.0, 0, 2., 4., 6., 8., 10.};
        Double_t *ZVrtxBins = ZBins;
        
        fPoolMgr = new AliEventPoolManager(fMEMaxPoolEvent, fMEMinTracks, NofCentBins, CentrORMultBins, NofZVrtxBins, ZVrtxBins);
        fPoolMgr->SetTargetValues(fMEMinTracks, 0.1, 5);
        if(!fPoolMgr) return kFALSE;
        return kTRUE;
    }
    
    //______________________________| Mixed Event Pool
    Bool_t DefineMixedEventPoolpp(){
        Int_t  NofCentBins  = 2;
        Double_t MBins[]={0,35,250.};
        Double_t * CentrORMultBins = MBins;
        
        Int_t NofZVrtxBins  = 10;
        Double_t ZBins[]={-10.0, -8.0, -6.0, -4.0, -2., 0, 2.0, 4.0, 6.0, 8.0, 10.0};
        Double_t *ZVrtxBins = ZBins;
        
        fPoolMgr = new AliEventPoolManager(fMEMaxPoolEvent, fMEMinTracks, NofCentBins, CentrORMultBins, NofZVrtxBins, ZVrtxBins);
        fPoolMgr->SetTargetValues(fMEMinTracks, 0.1, 5);
        if(!fPoolMgr) return kFALSE;
        return kTRUE;
    }
    
    
    //______________________________| Mixed Event Pool finer bins
    Bool_t DefineMixedEventPoolPbPb_largebins(){
        Int_t  NofCentBins  = 7;
        Double_t MBins[]={0., 7.5, 10., 20., 30., 40., 50., 100.1};
        Double_t * CentrORMultBins = MBins;
        
        Int_t NofZVrtxBins  = 6;
        Double_t ZBins[]={-10.0, -8.0, -4.0,  0, 4.0, 8.0, 10.};
        Double_t *ZVrtxBins = ZBins;
        
        fPoolMgr = new AliEventPoolManager(fMEMaxPoolEvent, fMEMinTracks, NofCentBins, CentrORMultBins, NofZVrtxBins, ZVrtxBins);
        fPoolMgr->SetTargetValues(fMEMinTracks, 0.1, 5);
        if(!fPoolMgr) return kFALSE;
        return kTRUE;
    }
    
    //______________________________| Mixed Event Pool
    Bool_t DefineMixedEventPoolpp_largebins(){
        Int_t  NofCentBins  = 1;
        Double_t MBins[]={0,250.};
        Double_t * CentrORMultBins = MBins;
        
        Int_t NofZVrtxBins  = 10;
        Double_t ZBins[]={-10.0, -8.0, -6.0, -4.0, -2.0, 0, 2.0, 4.0, 6.0, 8.0, 10.0};
        Double_t *ZVrtxBins = ZBins;
        
        fPoolMgr = new AliEventPoolManager(fMEMaxPoolEvent, fMEMinTracks, NofCentBins, CentrORMultBins, NofZVrtxBins, ZVrtxBins);
        fPoolMgr->SetTargetValues(fMEMinTracks, 0.1, 5);
        if(!fPoolMgr) return kFALSE;
        return kTRUE;
    }

    
    //______________________________| Mixed Event Pool Process
    Bool_t ProcessMixedEventPool(){
        if(!fMixedEvent) return kFALSE;
        if(!fPool->IsReady()) return kFALSE;
        if(fPool->GetCurrentNEvents()<fMEMinEventToMix) return kFALSE;
        return kTRUE;
    }
    
    
    //______________________________| All Used Ojects
    Bool_t    ftwoplus1;
    Bool_t    fEqualT1T2;
    Bool_t    fFineBinsME;
    Bool_t    fSetSystemValue;
    Bool_t    fRecoOrMontecarlo;
    Bool_t    fReadMC;
    Bool_t    fSetFilterBit; //
    Int_t     fbit ;   //
    TClonesArray *farrayMC;//!
    
    Double_t  fCentrOrMult; // Multiplicity of Event for D eff
    Double_t  fMinCentrality; // Minimun Centrality Value
    Double_t  fMaxCentrality; // Maximum Centrality Value
    
    Double_t  fTrigger1pTLowThr;
    Double_t  fTrigger1pTHighThr;
    Double_t  fTrigger2pTLowThr;
    Double_t  fTrigger2pTHighThr;
    
    
    Bool_t    fCutResonances;
    Bool_t    fCutConversions;
    Bool_t    ftwoTrackEfficiencyCut;
    Bool_t    fuseVarCentBins;
    Bool_t    fuseVarPtBins;
    Double_t  fAlpha;
    Bool_t    fBkgSE;
    Bool_t    fBkgSEBothSide;
    
    TH1F     *fHistNEvents;
    TH1F     *fHistTrackCutsQA;
    TH1F     *f256Phi;
    TH1F     *f512Phi;
    TH1F     *f768Phi;
    TH1F     *fHistCent;
    TH1F     *fHistT1CorrTrack;
    TH1F     *fHistT2CorrTrack;
    TList    *fOutputQA;   //!QA Output list
    TList    *fOutputCorr;//!corr Output list
    
    THnF *fThnEff;// eff histogram
    //TH3F     *fThnEff_clone; //!
    
    AliEventPool *fPool; //! Pool for event mixing
    AliEventPoolManager  *fPoolMgr;         //! event pool manager
    Bool_t   fMixedEvent;		// enable event mixing (default: ON)
    Int_t  fMEMaxPoolEvent;
    Int_t  fMEMinTracks;
    Int_t  fMEMinEventToMix;
    
    TH1F *fHistQA[9];//! QA histos
    TH1F *fHistTrigDPhi;//!dphi between T1 and T2
    TH2F *fControlConvResT1;//!control res histo T1
    TH2F *fControlConvResT2;//!control res hist T2
    TH2F *fControlConvResMT1;//! control res histo T1 mixed
    TH2F *fControlConvResMT2;//!control res histo T2 mixed
    
    TH2F *fEtaPhi;//!eta-phi
    
    TH1F *fEffCheck;
    
    TH1F *fNoMixedEvents;//
    TH2F *fMixStatCentorMult; //no of events in pool vs cent/multplicity
    TH2F *fMixStatZvtx; //no of events in pool vs zvtx
    
    
    
    
    ClassDef(AliAnalysisTaskDiJetCorrelationsAllb2b, 2); // example of analysis
};


class AliDPhiBasicParticleDiJet : public AliVParticle
{
public:
    AliDPhiBasicParticleDiJet(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
    : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
    }
    ~AliDPhiBasicParticleDiJet() {}
    
    // kinematics
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pt() const { return fpT; }
    virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
    
    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }
    
    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
    
    
    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Short_t Charge()      const { return fCharge; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
    
    virtual Bool_t IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }
    
    virtual void SetPhi(Double_t phi) { fPhi = phi; }
    
private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCharge;   // charge
    
    ClassDef( AliDPhiBasicParticleDiJet, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};
#endif
