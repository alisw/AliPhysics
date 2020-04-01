#ifndef _AliAnalysisTaskStronglyIntensiveCorrTree_H_
#define _AliAnalysisTaskStronglyIntensiveCorrTree_H_


////////////////////////////////////////////////////////////////////////
//
// Analysis class for 
//
// Look for correlations and strongly intensige quantity in eta 
//
// This class needs input AODs.
// The output is a list of analysis-specific containers.
//
//    Author:
//    Iwona Sputowska
// 
////////////////////////////////////////////////////////////////////////


class TAxis;
class TClonesArray;
class TList;
class TObjArray;

class AliAODEvent;
class AliAODHeader;
class AliEventPoolManager;

#include <TObject.h>
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
class AliAnalysisTaskStronglyIntensiveCorrTree : public AliAnalysisTaskSE { 
public: 
  AliAnalysisTaskStronglyIntensiveCorrTree();
  AliAnalysisTaskStronglyIntensiveCorrTree(const char *name);
  virtual ~AliAnalysisTaskStronglyIntensiveCorrTree();

  /*   virtual void NotifyRun(); */
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  
  
  void SelectCollisionCandidates(UInt_t triggerInfo)   {fTrigger = triggerInfo;}
  void SetTrackBit(Int_t trackFilter)                  {fTrackFilter = trackFilter;}
  void SetMCStatus(Bool_t isMC=kFALSE)                  {fIsMC = isMC;}
  void SetPtRange(Double_t ptMin, Double_t ptMax) {
    fPtMin = ptMin; fPtMax = ptMax;
  }
  void SetPhiRange(Double_t phiMin, Double_t phiMax) {
    fPhiMin = phiMin; fPhiMax = phiMax;
  }
  AliEventCuts 	        fEventCuts;                   
 

protected:
  TObjArray* GetAcceptedTracks(AliAODEvent*, TClonesArray*, Int_t);
  TObjArray* GetAcceptedTracks(TClonesArray*, Int_t);

private:
  TTree                *fTree;   //!
  TTree                *fTreeMC; //!
  TTree                *fTreePrim; //!
  Bool_t               fIsMC; //  MC flag            
  Int_t                fTrackFilter; // filter bit selection      
  UInt_t               fTrigger; // physics selection
  Double_t             fCentMin, fCentMax;  // centrality range
  Double_t             fPtMin, fPtMax;      // P_{T} range
  Double_t             fPhiMin, fPhiMax;    // #phi range
  Int_t                fSelectPrimaryMCParticles;     // 0: no, 1: yes, -1: only non-primary particles
  Int_t                fSelectPrimaryMCDataParticles; // 0: no, 1: yes, -1: only non-primary particles
  // Tree branches
  Int_t fRunNumber;  //!
  Float_t fVertexZ;  //!
  Float_t fMagneticField; //!
  Float_t fCent[4];  //!
  Short_t fNf[15];   //!
  Short_t fNb[15];   //!

  Short_t fNf_MC[15]; //!
  Short_t fNb_MC[15]; //!

  Short_t fNf_MCPrim[15]; //!
  Short_t fNb_MCPrim[15]; //!

  TList*               fOutputList; //! list of all Histo
  TH1F *               fHistEta ; //! eta Histogram
  TH1F *               fHistPt ;  //! pt Histogram
  TH1F *               fHistPhi ; //! phi Histogram
  TH1F *               fHistDcaX ;  //! dcaX Histogram
  TH1F *               fHistDcaY ;  //! dcaY Histogram
  TH1F *               fHistDcaZ ;  //! dcaZ Histogram
  TH2F *               fHist2D_EtaDcaX ; //! eta_dcaX Histogram
  TH2F *               fHist2D_EtaDcaY ; //! eta_dcaY Histogram
  TH2F *               fHist2D_EtaDcaZ ; //! eta_dcaZ Histogram
  TH2F *               fHist2D_EtaPt ;   //! eta_pt Histogram
  TH2F *               fHist2D_EtaPhi ;  //! eta_phi Histogram
  TH1F *               fHistEtaMC ; //! eta Histogram
  TH1F *               fHistPtMC ;  //! pt Histogram
  TH1F *               fHistPhiMC ; //! phi Histogram
  TH1F *               fHistDcaXMC ;  //! dcaX Histogram
  TH1F *               fHistDcaYMC ;  //! dcaY Histogram
  TH1F *               fHistDcaZMC ;  //! dcaZ Histogram
  TH1I *               fEventStatistics ;//! Events Stat
  TH2F *               fHist2D_EtaDcaXMC ; //! eta_dcaX Histogram
  TH2F *               fHist2D_EtaDcaYMC ; //! eta_dcaY Histogram
  TH2F *               fHist2D_EtaDcaZMC ; //! eta_dcaZ Histogram
  TH2F *               fHist2D_EtaPtMC ;   //! eta_pt Histogram
  TH2F *               fHist2D_EtaPhiMC ;  //! eta_phi Histogram
  TH1F *               fHistEtaPrim ; //! eta Histogram
  TH1F *               fHistPtPrim ;  //! pt Histogram
  TH1F *               fHistPhiPrim ; //! phi Histogram
  TH1F *               fHistDcaXPrim ;  //! dcaX Histogram
  TH1F *               fHistDcaYPrim ;  //! dcaY Histogram
  TH1F *               fHistDcaZPrim ;  //! dcaZ Histogram
  TH2F *               fHist2D_EtaDcaXPrim ; //! eta_dcaX Histogram
  TH2F *               fHist2D_EtaDcaYPrim ; //! eta_dcaY Histogram
  TH2F *               fHist2D_EtaDcaZPrim ; //! eta_dcaZ Histogram
  TH2F *               fHist2D_EtaPtPrim ;   //! eta_pt Histogram
  TH2F *               fHist2D_EtaPhiPrim ;  //! eta_phi Histogram
  
  
  
  
  AliAnalysisTaskStronglyIntensiveCorrTree(const AliAnalysisTaskStronglyIntensiveCorrTree&); // not implemented
  AliAnalysisTaskStronglyIntensiveCorrTree& operator=(const AliAnalysisTaskStronglyIntensiveCorrTree&); // not implemented
  
  ClassDef(AliAnalysisTaskStronglyIntensiveCorrTree, 1);
} ; 

class TrackInfoCorr : public TObject {
public:
   TrackInfoCorr(Double_t eta=0, Double_t phi=0,Double_t pt=0, Double_t dcaX=0, Double_t dcaY=0, Double_t dcaZ=0)
    : fEta(eta), fPhi(phi), fPt(pt), fDcaX(dcaX), fDcaY(dcaY), fDcaZ(dcaZ) {}
  virtual ~TrackInfoCorr() {}

  Double_t Eta()  const { return fEta; }
  Double_t Phi()  const { return fPhi; }
  Double_t Pt()  const { return fPt; }
  Double_t XAtDCA() const { return fDcaX; } 
  Double_t YAtDCA() const { return fDcaY; }
  Double_t ZAtDCA() const { return fDcaZ; }  

protected:
private:
  TrackInfoCorr(const TrackInfoCorr&);
  TrackInfoCorr& operator=(const TrackInfoCorr&);

  Double_t fEta;
  Double_t fPhi;
  Double_t fPt;
  Double_t fDcaX;
  Double_t fDcaY;
  Double_t fDcaZ;
  ClassDef(TrackInfoCorr, 1);
} ;

#endif // _AliAnalysisTaskStronglyIntensiveCorrTree_H_


