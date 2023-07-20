#ifndef AliAnalysisTaskDowangppDCAfit_H
#define AliAnalysisTaskDowangppDCAfit_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TString.h>
#include <THnSparse.h>
#include <TClonesArray.h>

// AliRoot includes
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"
#include "AliAODpidUtil.h"
#include "AliAODInputHandler.h"
class AliAnalysisTaskDowangppDCAfit : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDowangppDCAfit();
  AliAnalysisTaskDowangppDCAfit(const char *name);

  virtual ~AliAnalysisTaskDowangppDCAfit();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *); 

  
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }     
  Double_t GetMinPt() { return fMinPt; }   

    virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
    virtual void  SetFilterbit(UInt_t filterbit){fFilterbit = filterbit;}
    virtual void  SetNoClus(Int_t noclus){fNoClus = noclus;}
    virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
    virtual void  SetMinPt(Double_t minPt){fMinPt = minPt;}
    virtual void  SetMaxPt(Double_t maxPt){fMaxPt = maxPt;}
    
    bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi, float TOFtime);
    bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime);


    bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP, float TOFtime);
    bool IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP, float nsigmaTPCD);
    bool IsDeuteronNSigma(float mom, float nsigmaTPCD, float nsigmaTOFD,float swithch_mom);
    bool IsTritonNSigma(float mom, float nsigmaTPCT, float nsigmaTOFT);
    bool IsHe3NSigma(float mom, float nsigmaTPCH, float nsigmaTOFH);


    bool IsKaonNSigmaReal(float mom, float nsigmaTPCK, float nsigmaTOFK, float TOFtime);

    bool RejectFakeP(double *NSigmaList, float mom);
    int ReCentrCode(Float_t v0Centr);

    int ReParticleLabel(Int_t pdgCode);
    Float_t GetTOFBeta(AliAODTrack *track);
    Float_t GetMass2TOF(Float_t beta, AliAODTrack *track);
    bool WiolaDCut(float mom, float nsigmaTPCD, float nsigmaTOFD);
    bool IsDeuteronTPCdEdx(float mom, float dEdx);
    bool WiolaRejectPion(float mom,float nsigmaTPCpi,float nsigmaTOFpi);
    double re_kstar(TLorentzVector &Particle1, TLorentzVector &Particle2);
 private:
    virtual Float_t GetVertex(AliAODEvent* aod) const;
    virtual void Analyze(AliAODEvent* aod,Float_t v0Centr);
    Short_t GetPidCode(Int_t pdgCode) const;
    int GetImostByPDGCode(Int_t pdgCode);
    float AODEventCut(AliAODEvent* fAOD);
    // dowang
    void GetDCADis(AliAODEvent* fEvent,AliAODTrack *tAodTrack,float *DCA_re,int fDCAglobalTrack);
    int GetImost(Short_t pidCode,int TrackCharge);


    AliAODEvent*      fAOD;            //! AOD object
	AliMCEvent	*mcEvent;			//! MCEvent
    

    TClonesArray*     fMCarray;        //! MC particle branch
    AliAODMCHeader*   fMCheader;       //! MC header
    AliAnalysisUtils* fAnalysisUtil;   ///< Event selection
    AliAODpidUtil  *fAODpidUtil;
    AliPIDResponse *fPIDResponse;

    //
    // Cuts and options
    //
    Double_t     fVtxCut;             // Vtx cut on z position in cm
    UInt_t       fFilterbit;          // filter bit
    Double_t     fEtaCut;             // Eta cut used to select particles
    Int_t        fNoClus;	          // No of TPC clusters
    Double_t     fMinPt;              // Min pt - for histogram limits
    Double_t     fMaxPt;              // Max pt - for histogram limits


    //
    // Output objects
    //
    TList*     fListOfObjects;     //! Output list of objects
    TH1F *     EventDis;// ddd

    // TH2F * DCAxy_noCut[8][3][5];
    // TH2F * DCAz_noCut[8][3][5];

    // TH2F * DCAxy_check_method[8][3][5];
    // TH2F * DCAz_check_method[8][3][5];

    // TH2F * DCAxy_allCut[8][3][5];
    // TH2F * DCAz_allCut[8][3][5];

    // TH2F * dEdxVspT[8][3][5];
    // TH2F * m2TOFVspT[8][3][5];
    // TH1F * pTdis[8][3][5];
    // TH1F * Etadis[8][3][5];


    // TH2F * dEdxVspT_noCut[8][3][5];
    // TH2F * dEdxVspT_allCut[8][3][5];

    // TH1F *fHistQA_Event;
    // TH1F *fHistQA_Track;
    // TH1F *fHistQA_pInEvent[3];

	// TH1F *PileUPTracl_QA[3];   

    // TH1F *MCTrackQA[8];

    TH2F * DCAxy_pDetailFrac[2][8][5];// ddd
    TH2F * MomSmearing[2][5];// ddd

    TH1F * pTdisReco[2][5];// ddd
    //TH1F * pTdisTrue[2][5];
    TH1F * pTdisRecoWhenPrimary[2][5];// ddd
    TH1F * pTdisRecoIDAsP[2][5];// ddd

    ClassDef(AliAnalysisTaskDowangppDCAfit, 1);    //Analysis task for high pt analysis
    
};

#endif



