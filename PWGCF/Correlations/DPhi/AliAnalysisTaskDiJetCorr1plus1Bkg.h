#ifndef AliAnalysisTaskDiJetCorr1plus1Bkg_cxx
#define AliAnalysisTaskDiJetCorr1plus1Bkg_cxx

#include <THn.h>

class TH1F;
class TH2F;
class TH3F;
class THnSparse;

class AliESDEvent;
class AliESDtrackCuts;
class AliCentrality;
class AliMultiplicity;
class TList;
class TObjArray;
class AliEventPoolManager;

class AliAODEvent;
class AliAODTrack;
class AliAODHandler;
class AliAODInputHandler;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskDiJetCorr1plus1Bkg : public AliAnalysisTaskSE {
 public:
 AliAnalysisTaskDiJetCorr1plus1Bkg() : AliAnalysisTaskSE()
, fAOD(0)
, centrality(0)
, fOutputList(0)
, fAodTracks(0)
, fAodTracksT(0)
, fAodTracksA(0)
, vertex(0)
, vtxSPD(0)
, fCentrOrMult(0)
, fBit(0)
, fTriggerpTLowThr(0)
, fTriggerpTHighThr(0)
, fHistTrigDPhi(0)
, fHistTrigDPhiM(0)
, fHistTrigPEtaS(0)
, fHistTrigSEtaS(0) 
, fHistTrigPEtaM(0)
, fHistTrigSEtaM(0)
, fHistDeltaPhiT1T2BC(0)
, fHistDeltaPhiT1T2C1(0)
, fHistDeltaPhiT1T2AC(0)
, fHistDeltaEtaT1T2BC(0)
, fHistDeltaEtaT1T2AC(0)

    //, fHistEff(0)
, fEventCounter(0)
, fHistCent(0)
, f3DEffCor(0)
, fControlConvResT1(0)
, fControlConvResT2(0)
, fControlConvResMT1(0)
, fControlConvResMT2(0)
, fTHnCentZvtxDEtaDPhi1SE(0)
//, fTHnCentZvtxDEtaDPhi2SE(0)
, fTHnTrigCentZvtxpTtrig1SE(0)
//, fTHnTrigCentZvtxpTtrig2SE(0)
, fTHnCentZvtxDEtaDPhi1ME(0)
//, fTHnCentZvtxDEtaDPhi2ME(0)
, fTHnTrigCentZvtxpTtrig1ME(0)
, fThnEff(0)
, fEffCheck(0)
, fNoMixedEvents(0)
, fMixStatCentorMult(0)
, fMixStatZvtx(0)
//, fTHnTrigCentZvtxpTtrig2ME(0)
, fTrackArray(0)
, fMixingTracks(15000)//i ran on 25000 for without resonance correction
, fPoolMgr(0x0)
, useVarBins(kTRUE)
, fCutConversions(kTRUE)
, fCutResonances(kTRUE)
, twoTrackEfficiencyCut(kTRUE)
, fEtaOrdering(kTRUE)
, fSetSystemValue(kTRUE)
 
{
  for ( Int_t i = 0; i < 9; i++) 
    { 
      fHistQA[i] = NULL;
    }
 }
  
  AliAnalysisTaskDiJetCorr1plus1Bkg(const char *name);
  virtual ~AliAnalysisTaskDiJetCorr1plus1Bkg();
  
  virtual void   UserCreateOutputObjects();
  //virtual void   SetInputCorrection();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  virtual void	 SetEventMixing(Bool_t flag) { fFillMixed = flag; }
  virtual void   SetMixingTracks(Int_t tracks) {fMixingTracks = tracks;}
  virtual void   SetSystemType(Bool_t system) { fSetSystemValue = system;}
  virtual void   SetFilterBit(Int_t filterbit){fBit=filterbit;}//
  virtual void   SetTriggerpTValue(Double_t pTmin1, Double_t pTmax1){fTriggerpTLowThr = pTmin1, fTriggerpTHighThr = pTmax1;}
          void   SetEfficiencyWeightMap(THnF* hEff){fThnEff = hEff;}
          void    SetResonanceCut(Bool_t resCut){fCutResonances = resCut;}
          void    SetConversionCut(Bool_t conversionCut){fCutConversions = conversionCut;}
          void    SetTwoTrackEfficiencyCut(Bool_t TTRcut){twoTrackEfficiencyCut = TTRcut;}
  
 private:
  inline Float_t GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  inline Float_t GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  inline Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);
    
    Double_t GetTrackWeight(Double_t eta, Double_t pt, Double_t CentrOrMult, Double_t zVertex);
    
    
//______________________________| Mixed Event Pool PbPb
    Bool_t DefineMixedEventPoolPbPb(){
        

        Int_t poolsize = 500;  // Maximum number of events

        
        Int_t  NofCentBins  = 7;
        Double_t MBins[]={0.,7.5, 10., 20., 30., 40., 50., 100.1};
        Double_t * CentrORMultBins = MBins;
        
        Int_t NofZVrtxBins  = 10;
        Double_t ZBins[]={-10.0, -8.0, -6.0, -4.0, -2.0, 0, 2., 4., 6., 8., 10.};
        Double_t *ZVrtxBins = ZBins;
        
        fPoolMgr = new AliEventPoolManager(poolsize, fMixingTracks, NofCentBins, CentrORMultBins, NofZVrtxBins, ZVrtxBins);
        fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);
        if(!fPoolMgr) return kFALSE;
        return kTRUE;
    }
    
    //______________________________| Mixed Event Poolpp
    Bool_t DefineMixedEventPoolpp(){
        
        Int_t poolsize = 500;  // Maximum number of events
        Int_t  NofCentBins  = 1;
        Double_t MBins[]={0, 250.};
        Double_t * CentrORMultBins = MBins;
        
        Int_t NofZVrtxBins  = 10;
        Double_t ZBins[]={-10.0, -8.0, -6.0, -4.0, -2., 0, 2.0, 4.0, 6.0, 8.0, 10.0};
        Double_t *ZVrtxBins = ZBins;
        
        fPoolMgr = new AliEventPoolManager(poolsize, fMixingTracks, NofCentBins, CentrORMultBins, NofZVrtxBins, ZVrtxBins);
        fPoolMgr->SetTargetValues(fMixingTracks, 0.05, 5);
        if(!fPoolMgr) return kFALSE;
        return kTRUE;
    }

  
  AliAODEvent *fAOD;    //! ESD object
  AliCentrality *centrality;
  TList       *fOutputList; //! Output list
  AliAODTrack *fAodTracks;
  AliAODTrack *fAodTracksT;
  AliAODTrack *fAodTracksA;
  //AliAODtrackCuts *fTrackCuts;
  AliAODVertex    *vtxSPD;
  AliAODVertex *vertex;
    
  Double_t fCentrOrMult;
  Int_t fBit;
  Double_t fTriggerpTLowThr;
  Double_t  fTriggerpTHighThr;
  // TObjArray *fgoodTracks;
  TH1F *fEventCounter;//!
  TH1F *fHistCent;//!
  TH1F *fHistTrigDPhi;//!
  TH1F *fHistTrigDPhiM;//!
  TH1F *fHistTrigPEtaS;//!
  TH1F *fHistTrigSEtaS;//!
  TH1F * fHistTrigPEtaM;//!
  TH1F * fHistTrigSEtaM;//!
  //TH1D *fHistEff;
  TH1F *fHistQA[9];//!
  TH3F *f3DEffCor;//!
  TH2F *fControlConvResT1;//!
  TH2F *fControlConvResT2;//!
  TH2F *fControlConvResMT1;//!
  TH2F *fControlConvResMT2;//!
    
  TH3F * fHistDeltaPhiT1T2BC;//!
  TH2F * fHistDeltaPhiT1T2C1;//!
  TH2F * fHistDeltaPhiT1T2AC;//!
  TH3F * fHistDeltaEtaT1T2BC;//!
  TH2F * fHistDeltaEtaT1T2AC;//!
    
  THnSparse *fTHnCentZvtxDEtaDPhi1SE;
 // THnSparse *fTHnCentZvtxDEtaDPhi2SE;
  THnSparse *fTHnTrigCentZvtxpTtrig1SE;
 // THnSparse *fTHnTrigCentZvtxpTtrig2SE;
  THnSparse *fTHnCentZvtxDEtaDPhi1ME;
 // THnSparse *fTHnCentZvtxDEtaDPhi2ME;
  THnSparse *fTHnTrigCentZvtxpTtrig1ME;
 // THnSparse *fTHnTrigCentZvtxpTtrig2ME;
    
  THnF *fThnEff;
  TH1F *fEffCheck;
    
  TH1F *fNoMixedEvents;//
  TH2F *fMixStatCentorMult; //no of events in pool vs cent/multplicity
  TH2F *fMixStatZvtx; //no of events in pool vs zvtx
  
    


 
  //  TList *fTrackList; //! Output list
  TObjArray   *fTrackArray;
  AliEventPoolManager  *fPoolMgr;         //! event pool manager
  Bool_t   fFillMixed;		// enable event mixing (default: ON)
  Bool_t useVarBins;
  Int_t fMixingTracks;
  Bool_t fCutResonances;
  Bool_t fCutConversions;
  Bool_t twoTrackEfficiencyCut;
  Bool_t fEtaOrdering;
  Bool_t fSetSystemValue;
  // TList *fTrackListMixed;

  AliAnalysisTaskDiJetCorr1plus1Bkg(const AliAnalysisTaskDiJetCorr1plus1Bkg&); // not implemented
  AliAnalysisTaskDiJetCorr1plus1Bkg& operator=(const AliAnalysisTaskDiJetCorr1plus1Bkg&); // not implemented
  
  ClassDef(AliAnalysisTaskDiJetCorr1plus1Bkg, 1); // example of analysis
};

Float_t AliAnalysisTaskDiJetCorr1plus1Bkg::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign)
{
  //calculate dphistar for two track efficiency cut!
  Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);
  static const Double_t kPi = TMath::Pi();
  if (dphistar > kPi)
    dphistar = kPi*2 - dphistar;
  if (dphistar < -kPi)
    dphistar = -kPi*2 - dphistar;
  if (dphistar > kPi)
    dphistar = kPi * 2 - dphistar;
  return dphistar;
}

Float_t AliAnalysisTaskDiJetCorr1plus1Bkg::GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  //calculate inv ass squared
  Float_t tantheta1 = 1e10;

  if (eta1 < -1e-10 || eta1 > 1e-10)
    {
      Float_t expTmp = TMath::Exp(-eta1);
      tantheta1 = 2.0 * expTmp/(1.0 -expTmp*expTmp);
    }

  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
    {
      Float_t expTmp = TMath::Exp(-eta2);
      tantheta2 = 2* expTmp/(1.0 - expTmp*expTmp);
    }

  Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
  Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);

  Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * (TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2) ) );

  return mass2; 

}

Float_t AliAnalysisTaskDiJetCorr1plus1Bkg::GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
{
  // calculate inv mass squared approximately
  
  Float_t tantheta1 = 1e10;
  
  if (eta1 < -1e-10 || eta1 > 1e-10)
  {
    Float_t expTmp = 1.0-eta1+eta1*eta1/2-eta1*eta1*eta1/6+eta1*eta1*eta1*eta1/24;
    tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
  }
  
  Float_t tantheta2 = 1e10;
  if (eta2 < -1e-10 || eta2 > 1e-10)
  {
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


#endif
