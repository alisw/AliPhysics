#ifndef AliAnalysisFBMultFluct_cxx
#define AliAnalysisFBMultFluct_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class TH2F;
class TList;
class TProfile;
class TTree;

class AliAODEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisFBMultFluct : public AliAnalysisTaskSE {
 public:
   AliAnalysisFBMultFluct() : AliAnalysisTaskSE(){}//, fAOD(0), fOutputList(0)
  //AliAnalysisFBMultFluct();//, fAOD(0), fOutputList(0)
    /*
    , fAnaHistList(0x0)
    , fHistVertexZ(0), fHistVertexXY(0)
    , fHistCentrality(0)
    , fHistZDCN1(0), fHistZDCN2(0) 
    , fHistPt(0) , fHistEta(0)//, fHistEtaAsymBin(0)
    , fHistZDCAsym(0)
    , fHistZDCAsymvsVz(0)
    , fHistLog(0) 
    */
  //  {}

  AliAnalysisFBMultFluct(const char *name);
  virtual ~AliAnalysisFBMultFluct() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  static const Int_t nslice = 10;
  static const Int_t NMaxTrack = 8000;
  static const Int_t NMaxTracklet = 8000;

 private:
  //  AliESDEvent *fESD;    //! ESD object
  AliAODEvent *fAOD;    //! AOD object
  TList       *fOutputList; //! Output list

  Float_t asymMin[nslice]; //
  Float_t asymMax[nslice]; //

  Int_t fRunNumber, fNumberOfTracklets,fNumberOfTracks; //
  // There are four centrlaity information for each event
  // [0]="VOM","VOMEq","TRK","ZEMvsZDC"
  Float_t fCentPercentile[4]; //
  Float_t fVertexX,fVertexY,fVertexZ; //
  Float_t fzdcn1en,fzdcn2en,fzdcp1en,fzdcp2en,fzdcEM1en,fzdcEM2en; //
  Float_t fTrackletEta[NMaxTracklet],fTrackletPhi[NMaxTracklet]; //
  Float_t fTrackPt[NMaxTrack],fTrackEta[NMaxTrack],fTrackPhi[NMaxTrack]; //
  Int_t fTrackCharge[NMaxTrack], fTrackFilterMap[NMaxTrack]; //
  Float_t fV0Mult[64]; //
  Double_t Qx,Qy,Q1x,Q1y,Q2x,Q2y; //

  TH1F        *fHistVertexZ; //
  TH2F        *fHistVertexXY; //
  TH1F        *fHistCentrality; //
  TH1F        *fHistZDCN1; //
  TH1F        *fHistZDCN2; //
  TH1F        *fHistTrackPt; //! Pt spectrum
  TH1F        *fHistTrackEta; //! Eta spectrum
  TH1F        *fHistTrackletEta; // Eta spectrum for Tracklets 
  TH1F        *fHistTrackletEtaAsymBin[nslice]; // Eta spectrum for diff Asym Bin 
  TH1F        *fHistNeventZDCAsym; // To hold number of eventss for each asymmetry bin
  TH1F        *fHistZDCAsym; // ZDC Asymmetry in data
  TH2F        *fHistZDCAsymvsVz; // ZDC Asymmetry v/s vertex Z position
  TH2F        *fHistZDCAsymvsCent; // ZDC Asymmetry v/s Event Centrality
  TH1F        *fHistTrackletMult; // Tracklet Multiplicity within |eta|<1.0 
  TProfile    *fHistZDCAsymvsTrackletAsym; //
  TH2F        *fHistTracklets_pos_neg; //
  TH1F        *fHistLog; //

  TH1F        *fV0SubeventAPlane, *fV0SubeventCPlane; //
  TTree       *fEventTree; // output Tree 

  Double_t avfzdcn1, avfzdcn2; //

  AliAnalysisFBMultFluct(const AliAnalysisFBMultFluct&); // not implemented
  AliAnalysisFBMultFluct& operator=(const AliAnalysisFBMultFluct&); // not implemented
  
  Double_t GetAsymmetry(Double_t, Double_t, Double_t, Double_t);

  ClassDef(AliAnalysisFBMultFluct, 1); // example of analysis
};

#endif
