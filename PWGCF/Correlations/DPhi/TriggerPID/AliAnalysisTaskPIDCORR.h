#ifndef ALIANALYSISTASKPIDCORR_H
#define ALIANALYSISTASKPIDCORR_H

class TH1F;
class TH2F;
class TList;
class TObjArray;
class AliAODEvent;
class AliAODVertex;
class AliAODTrack;
class AliPIDResponse;
class AliPID;
class AliEventPoolManager;
class AliPIDCorrParticle;


#include <TObject.h> //PIDCorrParticle is a derived class from"TObject"

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskPIDCORR : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskPIDCORR();
    AliAnalysisTaskPIDCORR(const char *name);
    virtual ~AliAnalysisTaskPIDCORR();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    
 private:

 Bool_t SelectEvent(AliAODVertex *vertex);
 Int_t ClassifyTrack(AliAODTrack* track);

 void FillGlobalTracksArray();
 AliAODTrack* GetGlobalTrack(Int_t trackIndx);
 Double_t PhiRange(Double_t DPhi);
 Bool_t GetTrackStatus(AliAODTrack *track);
 Int_t GetTriggerPtBin(AliAODTrack *track);
 Bool_t TwoTrackEfficiency(AliAODTrack *Trig,AliAODTrack *Asso,Float_t TwoTrackEfficiencyCut,Float_t BSign);
 Bool_t TwoTrackEfficiencyBg(AliAODTrack *Trig,AliPIDCorrParticle *Asso,Float_t TwoTrackEfficiencyCut,Float_t BSign);
 Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t BSign);

 
    
AliAODEvent    *fAOD;
AliAODVertex   *fAODVertex;
AliPIDResponse *fPIDResponse;
TList          *fOutputList;
TH1F           *fHistPt;
TH2F           *fHistdEdx;
TH1F           *fHistNSigmaTPCPion[10];
TH2F           *fDihadronCorrelation[10];
TH2F           *fDihadronCorrelationPion[10];
TH2F           *fDihadronCorrelationProtonKaon[10];

TH2F           *fMixedEvent[10];
TH2F           *fMixedPion[10];
TH2F           *fMixedProtonKaon[10];

TH1F           *fTriggerPhiAll;
TH1F           *fTriggerPhiPion;
TH1F           *fTriggerPhiKaonProton;

TH1F           *fTriggerEtaAll;
TH1F           *fTriggerEtaPion;
TH1F           *fTriggerEtaKaonProton;

TH1F           *fAssoPhi;
TH1F           *fAssoEta;



TObjArray	 *fGlobalTracks;
TClonesArray     *fArrayMC;

TH1F             *fHistNSAll[3];
TH1F             *fHistNSPion[3];
TH1F		 *fHistNSProton[3];

TH1F             *fHistASAll[3];
TH1F             *fHistASPion[3];
TH1F		 *fHistASProton[3];

TH1F             *fHistBgAll[3];
TH1F             *fHistBgPion[3];
TH1F             *fHistBgProton[3];

TH1F             *fHistBulkAll[3];
TH1F             *fHistBulkPion[3];
TH1F		 *fHistBulkProton[3];
  
//Mixing functions
void DefineEventPool();
TObjArray *AcceptTracksforMixing(AliAODEvent *event);
AliEventPoolManager    *fPoolMgr; 

//Introduce Correction Factor
 Float_t GetEtaCorrectionFactorAsso(Double_t Eta);
 Float_t GetEtaCorrectionFactorTrigAll(Double_t Eta);
 Float_t GetEtaCorrectionFactorTrigPion(Double_t Eta);
 Float_t GetEtaCorrectionFactorTrigProton(Double_t Eta);
 Float_t EffEtaTrigPr;
 Float_t EffEtaTrigPi;
 Float_t EffEtaTrigAll;

//Identifying Associated
void IdentifyAssociated(Double_t ETA_trig,Double_t PHI_trig,Bool_t kPION,Bool_t kPROTON,AliAODTrack *track);
Int_t GetTOFPID(AliAODTrack *track);
Int_t GetTPCTOFPID(AliAODTrack *track);
Bool_t CheckTOF(AliVTrack *track);


    

	   // NEW HISTO to be declared here
    
    AliAnalysisTaskPIDCORR(const AliAnalysisTaskPIDCORR&); // not implemented
    AliAnalysisTaskPIDCORR& operator=(const AliAnalysisTaskPIDCORR&); // not implemented
    
    ClassDef(AliAnalysisTaskPIDCORR, 1); // example of analysis
};



class AliPIDCorrParticle : public TObject
{
  public:
    AliPIDCorrParticle(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
      : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
}
 
    virtual Double_t Pt() const { return fpT; }
    
    virtual Double_t Phi()        const { return fPhi; }
    
    
    virtual Double_t Eta()        const { return fEta; }
    
    
    virtual Short_t Charge()      const { return fCharge; }
    
    
    ~AliPIDCorrParticle() {}
private:
     AliPIDCorrParticle(const AliPIDCorrParticle&);  // not implemented
    AliPIDCorrParticle& operator=(const AliPIDCorrParticle&);  // not implemented
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCharge;   // charge
    
    ClassDef( AliPIDCorrParticle, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};



#endif

