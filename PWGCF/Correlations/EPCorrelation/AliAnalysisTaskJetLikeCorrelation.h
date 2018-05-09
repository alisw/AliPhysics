#ifndef ALIANALYSISTASKJETLIKECORRELATION_H
#define ALIANALYSISTASKJETLIKECORRELATION_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif


class TArrayI;
class TArrayD;
class TH2D;
class TH1D;
class TH3D;
class TList;
class TTree;
class TFile;
class THnSparse;
class AliAnalysisManager;
class AliExtractedTrack;
class AliExtractedEvent;
class AliAnalysisTaskFlowVectorCorrections;
class AliQnCorrectionsManager;
class AliAODEvent;
class AliMCEvent;
class AliAODMCHeader;

#include "AliEventPoolManager.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliQnCorrectionsManager.h"

#include "TObjArray.h"

#include <vector>
#include <string>

class AliExtractedTrack {
  public :
    AliExtractedTrack();
    AliExtractedTrack(const AliExtractedTrack &track);
    AliExtractedTrack &operator= (const AliExtractedTrack &track);
    AliExtractedTrack(AliExtractedTrack &&track);
    AliExtractedTrack &operator= (AliExtractedTrack &&track);
    virtual ~AliExtractedTrack();

    void SetPt(double pt) {fPt = pt;}
    void SetEta(double eta) {fEta = eta;}
    void SetPhi(double phi) {fPhi = phi;}
    void SetCharge(int charge) {fCharge = charge;}
    void SetFilterBit(unsigned long filterbit) { fFilterBit = filterbit; }
    void SetIsPrimary(bool isprimary) { kIsPrimary = isprimary; }

    double GetPt() {return fPt;}
    double GetPhi() {return fPhi;}
    double GetEta() {return fEta;}
    int GetCharge() {return fCharge;}
    unsigned long GetFilterBit() { return fFilterBit; }
    bool IsPrimary() { return kIsPrimary; }

  private :
    double fPt;
    double fEta;
    double fPhi;
    int fCharge;
    unsigned long fFilterBit;
    bool kIsPrimary;

  ClassDef(AliExtractedTrack, 4);
};

class AliExtractedEvent {
  public : 
    AliExtractedEvent();
    AliExtractedEvent(const AliExtractedEvent &event);
    AliExtractedEvent &operator= (const AliExtractedEvent &event);
    AliExtractedEvent(AliExtractedEvent &&event);
    AliExtractedEvent &operator= (AliExtractedEvent &&event);
    virtual ~AliExtractedEvent();

    void SetRunNumber(int runnumber) { fRunNumber = runnumber; }
    void SetEventID(unsigned long evid) { fEventID = evid; }
    void SetCentrality(double centrality) { fCentrality = centrality; }
    void SetZVertex (double zvtx) { fZVertex = zvtx; }
    void SetEventPlane(double eventplane) {fEventPlane = eventplane; }
    void SetEventPlaneV0A(double eventplane) {fEventPlaneV0A = eventplane; }
    void SetEventPlaneV0C(double eventplane) {fEventPlaneV0C = eventplane; }
    void SetEventPlaneTPC(double eventplane) {fEventPlaneTPC = eventplane; }
    void SetTracks (std::vector<AliExtractedTrack> tracks) { fTracks = tracks; }
    void AddTracks (AliExtractedTrack track) { fTracks.push_back(track); }
    void SetBSign (float bsign) { fBSign = bsign; }

    int GetRunNumber() { return fRunNumber; }
    unsigned long GetEventID() { return fEventID; }
    double GetCentrality() { return fCentrality; }
    double GetZVertex() { return fZVertex; }
    double GetEventPlane() { return fEventPlane; }
    double GetEventPlaneV0A() { return fEventPlaneV0A; }
    double GetEventPlaneV0C() { return fEventPlaneV0C; }
    double GetEventPlaneTPC() { return fEventPlaneTPC; }
    float GetBSign() { return fBSign; }
    std::vector<AliExtractedTrack> GetTracks() { return fTracks; }

  private : 
    int fRunNumber;
    unsigned long fEventID;
    float fBSign;
    double fCentrality;
    double fZVertex;
    std::vector<AliExtractedTrack> fTracks;
    double fEventPlane;
    double fEventPlaneV0A; 
    double fEventPlaneV0C; 
    double fEventPlaneTPC; 

    ClassDef(AliExtractedEvent, 4);
};

class AliAnalysisTaskJetLikeCorrelation : public AliAnalysisTaskSE {

  public : 
    
    enum ECollision {
      kpp,
      kpPb,
      kPbPb,
    };
    enum EPeriod {
      k10h,
      k11h,
      k11a,
      k15a,
      k15o,
      k17p,
    };
    enum EPlane {
      kIncl,
      kIn,
      kOut,
      kM1, 
      kM2,
      kM3,
      kM4,
      kPlane, 
    };
    AliAnalysisTaskJetLikeCorrelation();
    AliAnalysisTaskJetLikeCorrelation(const char *name);
//    AliAnalysisTaskJetLikeCorrelation(const AliAnalysisTaskJetLikeCorrelation &);
//    AliAnalysisTaskJetLikeCorrelation(const AliAnalysisTaskJetLikeCorrelation &&);
//    AliAnalysisTaskJetLikeCorrelation &operator=(AliAnalysisTaskJetLikeCorrelation &&);

    virtual ~AliAnalysisTaskJetLikeCorrelation();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    void ConnectInputHandler();


    int SetupMixing();

    void SetMCCorrection(bool mc) { fMCCorrection = mc; }
    void SetMCTruth(bool mc) { fMCTruth = mc; }
    void SetCentArray(TArrayD centarray) { fCentArray = centarray; }
//    void SetCentArray(vector<float> centarray) { fCentArray = centarray; }
    void SetZVertexArray(TArrayD zvertexarray) { fZVertexArray = zvertexarray; }
    void SetPttArray(TArrayD ptarray) { fPttArray = ptarray; }
    void SetPtaArray(TArrayD ptarray) { fPtaArray = ptarray; }
    void SetEtaCut(float etacut) { fEtaCut = etacut; }
    void SetPhiCut(float phicut) { fPhiCut = phicut; }
    void SetMinPtTrigCut(float ptcut) { fMinPtTrigCut = ptcut; }
    void SetMinimumPtABinForMerging(int ptamin) {fMinimumPtABinMerging = ptamin;}
    void SetTwoTrackEffCut(float twotrackcut) { fTwoTrackEffCut = twotrackcut; }
    void SetResonancesCut(float resocut) { fResonancesVCut = resocut; }
    void SetConversionsCut(float convcut) { fConversionsVCut = convcut; }
//    void SetEventMixingQueueSize(unsigned long size) { fEventMixingQueueSize = size; }
    void SetMinNumTrack(int minnumtrack) { fMinNumTrack = minnumtrack; }
    void SetTrackDepth(int trackdepth) { fTrackDepth = trackdepth; }
    void SetDebugOption(int debug) { fDebug = debug; } 
    void SetFilterBit(int filterbit) { fFilterBit = filterbit; }
    void SetCollision(ECollision col) { fCollision = col; }
    void SetMixingPoolSize(int size) { fMixingPoolSize = size; }
    void SetNumberOfPlanes(int nplane) { fNumberOfPlanes = nplane; }
    void UseSeparateMixingPool(bool pool) { bUseMixingPool = pool; }

    
    int GetCentBin(double);
    int GetZVertexBin(double);
    int GetPttBin(double);
    int GetPtaBin(double);
    double GetCentralityFromIP(double ip);
    int CheckInOut(double );
    int FillPt(double, int, int, double);

    float CalculatedPhiStar(float, float, float, float, float, float, float);
    int GetEventInformation(AliAODEvent *);
    int DoMixing(float , float , TObjArray *,float, int );
    
    void InitHistogramVectors();
    void InitHistograms();


  private :


    TList *fOutput[kPlane];      //!

    AliInputEventHandler *fInputHandler; //!
    AliInputEventHandler *fMCHandler; //!
    AliMCEvent *fMCEvent; //!
    AliAODMCHeader *fMCHeader; //!

    TH1D *fHistZVertex;     //!
    TH1D *fHistPt;     //!
    TH1D *fHistPhi;         //!
    TH1D *fHistCent;        //!
    TH1D *fHistEta;         //!
    TH1D *fHistPos;         //!
    TH1D *fHistNeg;         //!
    TH1D *fHistEventPlaneV0A;   //!
    TH1D *fHistEventPlaneV0C;   //!
    TH1D *fHistEventPlaneTPC;   //!

    AliEventPoolManager *fPoolMgr; //!
    TFile *fInputRoot; //!
    TTree *fInputTree; //!

    bool fMCCorrection; // Is MC or Data?
    bool fMCTruth; // MC Truth check
   
    bool bUseMixingPool; // Do we use separate mixing pool? 
    
    int fNumberOfPlanes; //
    int fMixingPoolSize;   //
    int fDebug;         //
    int fMinNumTrack;      //
    float fEtaCut;         //
    float fPhiCut;         //
    float fMinPtTrigCut;   //
    float fTwoTrackEffCut;         //
    int fFilterBit;         //
    int fTrackDepth;  // 
    ULong64_t fEventID; //!
    int fIsFirstEvent; //!
    int fMinimumPtABinMerging; //

    float fResonancesVCut;    //
    float fConversionsVCut;   //

    TArrayD fPttArray;         //
    TArrayD fPtaArray;         //
    TArrayD fCentArray;         //
    TArrayD fZVertexArray;         //

    std::vector< std::vector< std::vector< std::vector< std::vector<TH2D*> > > > >fHistdEtadPhiSame;          //!
    std::vector< std::vector< std::vector< std::vector< std::vector<TH2D*> > > > >fHistdEtadPhiMixed;         //!

    std::vector< std::vector< std::vector< std::vector< std::vector<TH2D*> > > > >fHistdEtadPhiSameMCCorrCont;          //!
    std::vector< std::vector< std::vector< std::vector< std::vector<TH2D*> > > > >fHistdEtadPhiSameMCCorrPrim;          //!


    TH2D *fHistNevtSame;         //!

    TH3D *fHistPtSame[kPlane];         //!

    THnSparse *fHistEtaSparse;   //!
    TProfile *fHistV2;             //! Only for V0A
    TProfile *fHistResolutionV2[3];   //!
    
    THnSparse *fHistContamination; //!
    AliExtractedEvent *fHighPtMixing; //!

    int fTreeStart; //!
    int fTreeEnd; //!
    std::vector< std::vector<int> > fPoolSize; //!

    double fCentPercentile; //!
    double fZVertex; //!
    double fEventPlane;  //!
    double fEventPlaneV0A;  //!
    double fEventPlaneV0C; //!
    double fEventPlaneTPC; //!
    double pi = TMath::Pi(); //
    
    ECollision fCollision; //
    EPeriod fPeriod; //!
    
    AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask; //! new Qn Vector framework
    AliQnCorrectionsManager *fFlowQnVectorMgr = NULL;  //! new ep

  protected :
    ULong64_t GetEventIdAsLong(AliVHeader *header ) {
      return ((ULong64_t)header->GetBunchCrossNumber()+
          (ULong64_t)header->GetOrbitNumber()*3564+
          (ULong64_t)header->GetPeriodNumber()*16777215*3564);
    }

    inline Float_t GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
    {
      // calculate inv mass squared
      // same can be achieved, but with more computing time with
      /*TLorentzVector photon, p1, p2;
        p1.SetPtEtaPhiM(triggerParticle->Pt(), triggerEta, triggerParticle->Phi(), 0.510e-3);
        p2.SetPtEtaPhiM(particle->Pt(), eta[j], particle->Phi(), 0.510e-3);
        photon = p1+p2;
        photon.M()*/

      Float_t tantheta1 = 1e10;

      if (eta1 < -1e-10 || eta1 > 1e-10)
      {
        Float_t expTmp = TMath::Exp(-eta1);
        tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
      }

      Float_t tantheta2 = 1e10;
      if (eta2 < -1e-10 || eta2 > 1e-10)
      {
        Float_t expTmp = TMath::Exp(-eta2);
        tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
      }

      Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
      Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);

      Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2 ) ) );

      //   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));

      return mass2;
    }

    inline Float_t GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
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

    ClassDef(AliAnalysisTaskJetLikeCorrelation, 3);
};

class MixedParticle : public TObject {
  public :
    MixedParticle(float pt, float eta, float phi, int charge) :
      fPt(pt), fEta(eta), fPhi(phi), fCharge(charge) {}
    virtual ~MixedParticle() {}

    float Eta() const { return fEta; }
    float Phi() const { return fPhi; }
    float Pt() const { return fPt; }
    int Charge() const { return fCharge; }
  private :
    MixedParticle(const MixedParticle&);
    MixedParticle& operator=(const MixedParticle&);

    float fPt;   //!
    float fEta;//!
    float fPhi;//!
    int fCharge;//!

    ClassDef(MixedParticle, 1);
};


#endif
