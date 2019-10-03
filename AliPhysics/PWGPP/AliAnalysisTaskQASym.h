// Analysis task for basic QA exploiting symmetries
//of global, TPC, and ITS tracks

#ifndef ALIANALYSISTASKQASYM_H
#define ALIANALYSISTASKQASYM_H
 

class TH1F;
class TH2F;
class TH3F;
class TList;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"
#include "TFile.h"

class AliAnalysisTaskQASym : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskQASym(const char *name = "AliAnalysisTaskQASym");
  virtual ~AliAnalysisTaskQASym() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   SetTrackType(Int_t type) {fTrackType = type;}  // set the track type: global, its and tpc
  virtual void   SetStandAloneTrack(Bool_t standAlone = kFALSE) {fStandAlone = standAlone;}  //needed for ITS tracks

  virtual void   SetNChargedRange(Int_t low = 0, Int_t high=1.0*1e7)//set multiplicity region of analysis
                                                                    // in order to compare different trigger settings
  {
    fLow  = low;
    fHigh = high;
  } 
  
  virtual void   SetCuts(AliESDtrackCuts* cuts)
     {fCuts = cuts;}

  virtual void   SetFieldOn(Bool_t b = kTRUE){fFieldOn = b;} 

  
 private:

  Int_t       fTrackType;       // track type
  Bool_t      fStandAlone;      // needed for ITS tracks
  Int_t       fLow;             // low Ncharges cut
  Int_t       fHigh;            // high Ncharges cut
  Bool_t      fFieldOn;         // field status

  TList       *fHists;          // List of histos

  //old
  TH1F        *fHistRECpt;      // pt 
  TH1F        *fEta;            // eta

  TH1F        *fEtaWidth;       // eta of tracks in Nch range (fLow <= Ncharged <= fHigh)
  TH1F        *fPhiWidth;       // phi of tracks in Nch range (fLow <= Ncharged <= fHigh)
  TH1F        *fDcaWidth;       // dca of tracks in Nch range (fLow <= Ncharged <= fHigh)
  TH1F        *fPtWidth;        // pt of tracks in Nch range (fLow <= Ncharged <= fHigh)

  TH2F        *fEtaPhi;         // eta-phi
  TH1F        *fEtaPt;          // eta over pt 
  TH1F        *fQPt;            // charge over pt 
  TH1F        *fDca;            // distance of closest approach
  TH1F        *fDcaZ;            // distance of closest approach
  TH1F        *fqRec;           // reconstrcuted charge
  TH1F        *fSigmaPtHist;    // sigma_pT

   //positive und negative tracks
  TH1F        *fRecPtPos;      // pt of pos tracks
  TH1F        *fRecPtNeg;      // pt of neg tracks
  TH1F        *fRecPhiPos;     // phi of pos. tracks
  TH1F        *fRecPhiNeg;     // phi of neg. tracks
  TH1F        *fRecEtaPos;     // eta of neg. tracks
  TH1F        *fRecEtaNeg;     // eta of neg. tracks
  TH1F        *fRecEtaPtPos;   // eta over pt of neg. tracks
  TH1F        *fRecEtaPtNeg;   // eta over pt of neg. tracks
  TH1F        *fRecDcaPos;     // distance of closest approach of neg. tracks
  TH1F        *fRecDcaNeg;     // distance of closest of neg. tracks
  TH1F        *fRecDcaNegInv;  // invers dca of neg. tracks
  TH1F        *fRecDPos;       // impact parameter of neg. tracks
  TH1F        *fRecDNeg;       // impact parameter of neg. tracks

  // two sides of TPC -> Eta/Theta
  TH1F        *fRecQPtPosEta;   //charge/pT for pos. eta
  TH1F        *fRecQPtNegEta;   //charge/pT for neg. eta
  TH1F        *fRecPtPosEta;    //pT        for pos. eta
  TH1F        *fRecPtNegEta;    //pT        for neg. eta
  TH1F        *fRecPhiPosEta;   //phi       for pos. eta
  TH1F        *fRecPhiNegEta;   //phi       for neg. eta
  TH1F        *fRecDcaPosEta;   //dca       for pos. eta
  TH1F        *fRecDcaNegEta;   //dca       for neg. eta
  TH1F        *fRecDPosEta;     //d         for pos. eta
  TH1F        *fRecDNegEta;     //d         for neg. eta

  // 2D: all measures as function of z of first trackpoint 
  TH2F        *fRecPtPosVz;     //pt-zfirst of pos tracks
  TH2F        *fRecPtNegVz;     //pt-zfirst of neg tracks
  TH2F        *fRecEtaPosVz;    //eta-zfirst of pos tracks
  TH2F        *fRecEtaNegVz;    //eta-zfirst of neg tracks
  TH2F        *fRecPhiPosVz;    //phi-zfirst of pos tracks
  TH2F        *fRecPhiNegVz;    //phi-zfirst of neg tracks
  TH2F        *fSignedDcaPosVz; //dca-zfirst of pos tracks
  TH2F        *fSignedDcaNegVz; //dca-zfirst of neg tracks
  TH2F        *fRecQPtPosEtaVz; //charge/pT-zfirst of pos tracks
  TH2F        *fRecQPtNegEtaVz; //charge/pT-zfirst of neg tracks
  TH2F        *fRecEtaPtPosVz;  //eta/pT-zfirst of pos tracks
  TH2F        *fRecEtaPtNegVz;  //eta/pT-zfirst of neg tracks


  //high
  TH1F * fDeltaPhiAll;         // phiLeaingTracks-phiOthers
  TH2F * fDeltaPhiLeading;     // phiLeaingTracks-phiOthers vs. phiLeading
  TH1F * fDiffDcaD;            // d-dca

  //sim
  TH1F * fPhiRec;              //phi
  TH1F * fThetaRec;            //theta
  TH1F * fNumber;              //Number of tracks per event
  TH1F * fNumberAfterCut;      //Number of tracks per event after cuts
  TH1F * fVx;                  // x of first track point
  TH1F * fVy;                  // y of first track point
  TH1F * fVz;                  // z of first track point
  TH1F * fNVertexSPD;          //number of vertices SPD
  TH1F * fNVertexTracks;       //number of vertices of Tracks

  //new
  TH2F        *fRecDcaPosPhi;     //dca-phi for pos.
  TH2F        *fRecDcaNegPhi;     //dca-phi for neg.
  TH2F        *fRecPtPosPhi;      //pt-phi for pos.
  TH2F        *fRecPtNegPhi;      //pt-phi for neg.
  TH2F        *fRecEtaPosPhi;     //eta-phi for pos.
  TH2F        *fRecEtaNegPhi;     //eta-phi for neg.
  TH2F        *fRecQPtPhi;        //charge/pt-phi
  TH2F        *fRecEtaPtPosPhi;   //eta/pt-phi for neg.
  TH2F        *fRecEtaPtNegPhi;   //eta/pt-phi for pos.

  TH1F        *fRecPtPosEtaPos;   //pt for pos tracks and pos eta
  TH1F        *fRecPtNegEtaPos;   //pt for neg tracks and pos eta
  TH1F        *fRecPtPosEtaNeg;   //pt for pos tracks and neg eta
  TH1F        *fRecPtNegEtaNeg;   //pt for neg tracks and neg eta

  TH1F        *fRec1PtPosEtaPos;   //1/pt for pos tracks and pos eta
  TH1F        *fRec1PtNegEtaPos;   //1/pt for neg tracks and pos eta
  TH1F        *fRec1PtPosEtaNeg;   //1/pt for pos tracks and neg eta
  TH1F        *fRec1PtNegEtaNeg;   //1/pt for neg tracks and neg eta

  TH1F        *fRecPhiPosEtaPos;   //phi for pos tracks and pos eta
  TH1F        *fRecPhiNegEtaPos;   //phi for neg tracks and pos eta
  TH1F        *fRecPhiPosEtaNeg;   //phi for pos tracks and neg eta
  TH1F        *fRecPhiNegEtaNeg;   //phi for neg tracks and neg eta

  TH2F        *fRecDcaPosPhiEtaPos;  //dca-phi for pos tracks and pos eta
  TH2F        *fRecDcaNegPhiEtaPos;  //dca-phi for neg tracks and pos eta
  TH2F        *fRecDcaPosPhiEtaNeg;  //dca-phi for pos tracks and neg eta
  TH2F        *fRecDcaNegPhiEtaNeg;  //dca-phi for neg tracks and neg eta

  TH2F        *fRecDcaPosPtEtaPos;  //dca-pt for pos tracks and pos eta
  TH2F        *fRecDcaNegPtEtaPos;  //dca-pt for neg tracks and pos eta
  TH2F        *fRecDcaPosPtEtaNeg;  //dca-pt for pos tracks and neg eta
  TH2F        *fRecDcaNegPtEtaNeg;  //dca-pt for neg tracks and neg eta

  TH2F        *fRecPtPosPhiEtaPos;  //pt-phi for pos tracks and pos eta
  TH2F        *fRecPtNegPhiEtaPos;  //pt-phi for neg tracks and pos eta 
  TH2F        *fRecPtPosPhiEtaNeg;  //pt-phi for pos tracks and neg eta
  TH2F        *fRecPtNegPhiEtaNeg;  //pt-phi for neg tracks and neg eta

  //  TH3F        *fRecDcaPhiPtPosEtaPos; //dca-pt-phi for pos tracks and pos eta
  //  TH3F        *fRecDcaPhiPtNegEtaPos; //dca-pt-phi for neg tracks and pos eta
  //  TH3F        *fRecDcaPhiPtPosEtaNeg; //dca-pt-phi for pos tracks and neg eta
  //  TH3F        *fRecDcaPhiPtNegEtaNeg; //dca-pt-phi for neg tracks and neg eta

  TH2F        *fEtavPt;                 // eta vs pt 
  TH2F        *fPhivPt;                 // phi vs pt 

  TH2F        *fCompareTPCparam;        // TPC param

  TH1F        *fITSlayer;               // ITS layer
  TH2F        *fITSlayerEta;            // ITS layer vs eta
  TH2F        *fITSlayerPhi;            // ITS layer vs phi

  AliESDtrackCuts* fCuts;               // List of cuts

  // four different vertex types: primary, spd, tracks, tpc
  TH1F * fVertexX[4];             // x of vertex
  TH1F * fVertexY[4];             // y of vertex
  TH1F * fVertexZ[4];             // z of vertex

  // sectors of TPC 
  TH1F        *fRecPtTpcSector[18];     //pt for TPC sectors
  TH1F        *fRecEtaTpcSector[18];    //eta for TPC sectors
  TH1F        *fRecQPtTpcSector[18];    //charge/pt for TPC sectors
  TH1F        *fRecEtaPtTpcSector[18];  //eta/pt for TPC sectors
  TH1F        *fSignedDcaTpcSector[18]; //dca for TPC sectors


  // 7 different case of hit in ITS ladders
  TH1F        *fRecPtPosLadder[7];  //pt for pos tracks
  TH1F        *fRecPtNegLadder[7];  //pt for neg tracks
  TH1F        *fRecPhiPosLadder[7]; //phi for pos tracks
  TH1F        *fRecPhiNegLadder[7]; //phi for neg tracks
  TH1F        *fRecEtaPosLadder[7]; //eta for pos tracks
  TH1F        *fRecEtaNegLadder[7]; //eta for neg tracks
  TH1F        *fSignDcaPos[7];      //dca for pos tracks
  TH1F        *fSignDcaNeg[7];      //dca for neg tracks
  TH1F        *fSignDcaNegInv[7];   //-dca for neg tracks
  TH1F        *fPtSigmaPos[7];      //sigma_pT for pos tracks
  TH1F        *fPtSigmaNeg[7];      //sigma_pT for neg tracks
  TH1F        *fqPtRec[7];          // charge/pt 
  TH2F        *fDcaSigmaPos[7];     // dca - sigma_pT for pos tracks
  TH2F        *fDcaSigmaNeg[7];     // dca - sigma_pT for neg tracks

  TH1F        *fEtaBinPt[3][2];      // eta histogram for different pt bins and different charges
  TH1F        *fPhiBinPt[3][2];      // phi histogram for different pt bins and different charges
  TH1F        *fDcaBinPt[3][2];      // dca histogram for different pt bins and different charges
  TH2F        *fEtaPhiBinPt[3][2];   // eta-phi histogram for different pt bins and different charges

  
  


  AliAnalysisTaskQASym(const AliAnalysisTaskQASym&); // not implemented
  AliAnalysisTaskQASym& operator=(const AliAnalysisTaskQASym&); // not implemented
  
  ClassDef(AliAnalysisTaskQASym, 1); // Basic QA exploiting symmetries
};

#endif
