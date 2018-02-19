/*
 * AliFemtoDreamv0Cuts.h
 *
 *  Created on: Dec 12, 2017
 *      Author: gu74req
 */

#ifndef ALIFEMTODREAMV0CUTS_H_
#define ALIFEMTODREAMV0CUTS_H_
#include "Rtypes.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0MCHist.h"
#include "AliFemtoDreamv0Hist.h"
#include "AliFemtoDreamv0.h"
class AliFemtoDreamv0Cuts {
 public:
  AliFemtoDreamv0Cuts();
  virtual ~AliFemtoDreamv0Cuts();
  static AliFemtoDreamv0Cuts* LambdaCuts(bool isMC,bool CPAPlots,
                                         bool SplitContrib);
  //Setters for plots
  void SetIsMonteCarlo(bool isMC){fMCData=isMC;};
  bool GetIsMonteCarlo(){return fMCData;};
  void SetPlotCPADist(bool plot) {fCPAPlots=plot;};
  void SetPlotContrib(bool plot) {fContribSplitting=plot;};
  void SetAxisInvMassPlots(int nBins,double minMass,double maxMass) {
    fNumberXBins=nBins;fAxisMinMass=minMass;fAxisMaxMass=maxMass;
  }
  //Setters for the daughter track cuts
  void SetPosDaugterTrackCuts(AliFemtoDreamTrackCuts *cuts){fPosCuts=cuts;};
  void SetNegDaugterTrackCuts(AliFemtoDreamTrackCuts *cuts){fNegCuts=cuts;};
  //Setters for PDG Codes of the daughters+v0
  void SetPDGCodev0(int pdgCode) {fPDGv0=pdgCode;};
  int GetPDGv0() const {return fPDGv0;};
  void SetPDGCodePosDaug(int pdgCode) {fPDGDaugP=pdgCode;};
  int GetPDGPosDaug() const {return fPDGDaugP;};
  void SetPDGCodeNegDaug(int pdgCode) {fPDGDaugN=pdgCode;};
  int GetPDGNegDaug() const {return fPDGDaugN;};
  //Setters v0 cuts
  void SetCheckOnFlyStatus(bool val) {fOnFlyStatus=val;fCutOnFlyStatus=true;};
  void SetCutCharge(int charge) {fCutCharge=true;fCharge=charge;};
  void SetPtRange(double pmin,double pmax) {
    fpTmin=pmin;fpTmax=pmax;fCutPt=true;
  }
  void SetKaonRejection(double MassLow,double MassUp) {
    fKaonRejection=true;fKaonRejLow=MassLow;fKaonRejUp=MassUp;
  };
  void SetCutMaxDecayVtx(double maxDecayVtx)  {
    fCutDecayVtxXYZ=true;fMaxDecayVtxXYZ=maxDecayVtx;
  };
  void SetCutTransverseRadius(double minRadius,double maxRadius) {
    fMinTransRadius=minRadius;fMaxTransRadius=maxRadius;fCutTransRadius=true;
  }
  void SetCutDCADaugToPrimVtx(double minDCA) {
    fMinDCADaugToPrimVtx=minDCA;fCutMinDCADaugPrimVtx=true;
  };
  void SetCutDCADaugTov0Vtx(double maxDCA) {
    fMaxDCADaugToDecayVtx=maxDCA;fCutMaxDCADaugToDecayVtx=true;
  }
  void SetCutCPA(double cpa) {fMinCPA=cpa;fCutCPA=true;};
  void SetCutInvMass(double width){fInvMassCutWidth=width;fCutInvMass=true;};
  void Init();
  TList *GetQAHists() {return fHistList;};
//  TList *GetQAHistsPosDaug() {return fPosCuts->GetQAHists();};
//  TList *GetQAHistsNegDaug() {return fNegCuts->GetQAHists();};
  TList *GetMCQAHists() {return fMCHistList;};
//  TList *GetMCQAHistsPosDaug() {return fPosCuts->GetMCQAHists();};
//  TList *GetMCQAHistsNegDaug() {return fNegCuts->GetMCQAHists();};
  bool isSelected(AliFemtoDreamv0 *v0);
  TString ClassName(){return "v0Cuts";};
 private:
  bool RejectAsKaon(AliFemtoDreamv0 *v0);
  bool DaughtersPassCuts(AliFemtoDreamv0 *v0);
  bool MotherPassCuts(AliFemtoDreamv0 *v0);
  bool CPAandMassCuts(AliFemtoDreamv0 *v0);
  void BookQA(AliFemtoDreamv0 *v0);
  void BookMC(AliFemtoDreamv0 *v0);
  void BookTrackCuts();
  void FillMCContributions(AliFemtoDreamv0 *v0);
  double CalculateInvMass(AliFemtoDreamv0 *v0,int PDGPosDaug,int PDGNegDaug);
  TList *fHistList;                   //!
  TList *fMCHistList;                 //!
  AliFemtoDreamv0MCHist *fMCHist;     //!
  AliFemtoDreamv0Hist *fHist;         //!
  //Here the cut values for the Daughters are stored
  AliFemtoDreamTrackCuts *fPosCuts;   //
  AliFemtoDreamTrackCuts *fNegCuts;   //
  //These are all the cuts directly linked to the v0
  bool fMCData;                       //
  bool fCPAPlots;                     //
  bool fContribSplitting;             //
  bool fCutOnFlyStatus;               //
  bool fOnFlyStatus;                  //
  bool fCutCharge;                    //
  int fCharge;                        //
  bool fCutPt;                        //
  double fpTmin;                      //
  double fpTmax;                      //
  bool fKaonRejection;                //
  double fKaonRejLow;                 //
  double fKaonRejUp;                  //
  bool fCutDecayVtxXYZ;               //
  double fMaxDecayVtxXYZ;             //
  bool fCutTransRadius;               //
  double fMinTransRadius;             //
  double fMaxTransRadius;             //
  bool fCutMinDCADaugPrimVtx;         //
  double fMinDCADaugToPrimVtx;        //
  bool fCutMaxDCADaugToDecayVtx;      //
  double fMaxDCADaugToDecayVtx;       //
  bool fCutCPA;                       //
  double fMinCPA;                     //
  bool fCutInvMass;                   //
  double fInvMassCutWidth;            //
  //Range for the axis of the hists
  double fAxisMinMass;                //
  double fAxisMaxMass;                //
  int fNumberXBins;                   //
  //PDG Codes of the Mother and the Daughter needed for Inv Mass Calc. and
  //matching in the MC Sample
  int fPDGv0;
  int fPDGDaugP;
  int fPDGDaugN;
  ClassDef(AliFemtoDreamv0Cuts,1)
};

#endif /* ALIFEMTODREAMV0CUTS_H_ */
