/// \class AliAnalysisTaskHypTritKf
/// \brief Hypertriton Analysis in two particle decay channel
///
/// Hypertriton candidates are identified using the KFParticle package.


#ifndef ALIANALYSISTASKHypTritKf_H
#define ALIANALYSISTASKHypTritKf_H

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDpid;
class AliESDtrackCuts;
class AliESDv0;
class AliESDVertex;
class AliESDInputHandler;
class AliESDtrack;

#include "AliTRDonlineTrackMatching.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliEventCuts.h"

#define HomogeneousField ///homogenous field in z direction
#include "AliITStrackerMI.h"
#include "AliITStrackMI.h"
#include "KFParticle.h"
#include "KFVertex.h"
#include "KFPVertex.h"
#include "KFPTrack.h"

class KFParticle;
class KFVertex;
class AliITStrackerMI;

class KFParticle3LH : public KFParticle {
  public:
    Bool_t CheckDaughter(KFParticle3LH daughter) {
      Float_t m[8], mV[36], D[3][3];
      if(KFParticleBase::GetMeasurement(daughter, m, mV, D)) return kTRUE;
      return kFALSE;    
    }
  
    Bool_t AddDaughters(KFParticle3LH daughter1, KFParticle3LH daughter2, int constructMethod) {
      if (constructMethod == 2) {
        KFParticle3LH KFCheck1;
        KFCheck1.AddDaughter(daughter1);
        if (!KFCheck1.CheckDaughter(daughter2)) return kFALSE;
        KFParticle3LH KFCheck2;
        KFCheck2.AddDaughter(daughter2);
        if (!KFCheck2.CheckDaughter(daughter1)) return kFALSE;
      } 
      this->SetConstructMethod(constructMethod);
      this->AddDaughter(daughter1);
      this->AddDaughter(daughter2);
      return kTRUE;
    } 
  
  Double_t GetParticleMass() {
    Double_t E  = this->GetE();
    Double_t Px = this->GetPx();
    Double_t Py = this->GetPy();
    Double_t Pz = this->GetPz();
    return TMath::Sqrt(E*E - Px*Px - Py*Py - Pz*Pz);
  }
};


class AliAnalysisTaskHypTritKf : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHypTritKf();
  AliAnalysisTaskHypTritKf(const char *name);
  virtual ~AliAnalysisTaskHypTritKf();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);


  void SetPeriod(Int_t period = 2015) {fPeriod = period;};
  void SetTriggerMask(UInt_t triggerMask = AliVEvent::kINT7) {fTriggerMask = triggerMask;};
  void SetBetheSplines(Bool_t betheSplines = kTRUE ) {fBetheSplines = betheSplines;};
  void SetParamsHe(Double_t params[6]) { for(Int_t i=0; i < 6; i++) fBetheParamsHe[i] = params[i];};
  void SetParamsT(Double_t params[6]) { for(Int_t i=0; i < 6; i++) fBetheParamsT[i] = params[i];};
  void SetUseExternalSplines(Bool_t extspl = kFALSE) {fUseExternalSplines = extspl;}; 
  
  
 private:
  AliESDInputHandler     *fInputHandler;        //!<! Input handler
  AliESDpid              *fPID;                 //!<! ESD pid
  AliESDEvent            *fESDevent;            //!<! ESD event
  TH2F                   *fHistdEdx;            //<   Histogram of Tpc dEdx for pid qa
  TH2F                   *fHistdEdxV0;          //<   Histogram of Tpc dEdx for pid qa
  TH1F                   *fHistNumEvents;       //<   Histogram of number of events
  TH1F                   *fHistTrigger;         //<   Histogram of trigger for all events 
  TH1F                   *fHistV0;              //<   Histogram of trigger for all V0s 
  TTree                  *fTree;                //<   Tree containing reduced events
  TTree                  *fTreeMCGen;           //<   Tree containing reduced MC events
  TList                  *fHistogramList;       //<   List of histograms
  TVector3                fPrimaryVertex;       //!<! Vector of primary vertex of collision
  TVector3                fSecondaryVertex;     //!<! Vector of decay vertex 
  Bool_t                  fMCtrue;              //< Flag for activating MC analysis (set automatically)
  AliEventCuts            fEventCuts;           //< 2015 event cuts as advised by PDG (AliEventCuts)
  UInt_t                  fTriggerMask;         //< Triggermask for event cuts
  Int_t                   fPeriod;              //< Data period for centrality selector
  Bool_t                  fBetheSplines;        //< Switch between built in BetheSplines and personal Fit
  Bool_t                  fUseExternalSplines;  //< Use Splines given in Add file
  Double_t                fBetheParamsHe[6];    //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Double_t                fBetheParamsT[6];     //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Int_t                   fYear;          
  
  Double_t                tPVx;									// primary Vertex position
  Double_t                tPVy;
  Double_t                tPVz;
  Double_t                tSVx;									// secondary Vertex position
  Double_t                tSVy;
  Double_t                tSVz;
      
  Double_t                tM;
  Double_t                tCt;
  Double_t                tPt; 
  Double_t                tP;
	Double_t                tPx;
	Double_t                tPy;
  Double_t                tPz;
  Double_t                tHePx;
  Double_t                tHePy;
  Double_t                tHePz;
	Double_t                tPiPx;
  Double_t                tPiPy;
  Double_t                tPiPz;												
  Double_t                tDca; 
  Double_t                tCosPA; 
  Double_t                tY; 
  Double_t                tPiDca; 
  Double_t                tHeDca; 
  Double_t                tHeDedxSigma; 
  Double_t                tTDedxSigma; 
  Double_t                tHeP; 
  Double_t                tPiP; 
  Double_t                tHeDedx; 
  Double_t                tPiDedx; 
  Double_t                tTrig; 
  Double_t                tHePt; 
  Double_t                tPiPt; 
  Double_t                tImpPiHe;
  Double_t                tChi2;
  Double_t                tNDF;
  Int_t                   tZ; 
  Double_t                tMc; 
  Double_t                tMagField;
  Int_t                   tRunnumber;
  Double_t                tKinkHe;
  Double_t                tTPCrefitHe;
  Double_t                tITSrefitHe;  
  Double_t                tnTPCclusterHe;
  Double_t                tnITSclusterHe;
  Double_t                tTPCchi2He;  
  Double_t                tKinkPi;
  Double_t                tTPCrefitPi;
  Double_t                tITSrefitPi;  
  Double_t                tnTPCclusterPi;
  Double_t                tnITSclusterPi;
  Double_t                tTPCchi2Pi;  
  
  Double_t                tMultV0M;      // multiplicity estimators
  Double_t                tMultOfV0M;      
  Double_t                tMultSPDTracklet;  
  Double_t                tMultSPDCluster;  
  Double_t                tMultRef05;      
  Double_t                tMultRef08;  

  Int_t                   tTRDvalid;    // has valid TRD track
  Int_t                   tTRDtrigHNU;  // HNU fired by track
  Int_t                   tTRDtrigHQU;  // HQU fired by track
  Int_t                   tTRDPid;
  Int_t                   tTRDnTracklets;
  Int_t                   tTRDPt;
  Int_t                   tTRDLayerMask;
  Double_t                tTRDSagitta;
  Int_t                   tTRDStack;
  Int_t                   tTRDSector;
  UInt_t                  tTRDPID0;
  UInt_t                  tTRDPID1;
  UInt_t                  tTRDPID2;
  UInt_t                  tTRDPID3;
  UInt_t                  tTRDPID4;
  UInt_t                  tTRDPID5;    
       
  void MCStackLoop(AliMCEvent* mcEvent);
  void CalculateV0(const AliESDtrack& trackN, const AliESDtrack& trackP, AliPID::EParticleType typeNeg, AliPID::EParticleType typePos, AliMCEvent* mcEvent);
  Bool_t TriggerSelection(AliMCEvent* mcEvent);
  Float_t GetInvPtDevFromBC(Int_t b, Int_t c);
  Double_t Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
  void SetBetheBlochParams(Int_t runNumber);
  Double_t TRDtrack(AliESDtrack* esdTrack);
  KFParticle3LH CreateKFParticle(AliExternalTrackParam& track, float Mass, Int_t Charge);
  KFVertex CreateKFVertex(const AliVVertex& vertex);
  AliAnalysisTaskHypTritKf(const AliAnalysisTaskHypTritKf&);
  AliAnalysisTaskHypTritKf &operator=(const AliAnalysisTaskHypTritKf&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskHypTritKf, 1);
  /// \endcond
};
#endif
