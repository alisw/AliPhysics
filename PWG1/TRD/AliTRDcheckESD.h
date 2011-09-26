#ifndef ALITRDCHECKESD_H
#define ALITRDCHECKESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDcheckESD.h 27496 2008-07-22 08:35:45Z cblume $ */

/////////////////////////////////////////////////////
//
// Check basic detector results at ESD level
//
// Author
//   Alex Bercuci <A.Bercuci@gsi.de>
//
//////////////////////////////////////////////////////

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTaskSE.h"
#endif

class AliESDEvent;
class AliMCEvent;
class AliESDpid;
class AliCFContainer;
class AliAnalysisCuts;
class TH1;
class TH2;
class TH1F;
class TH1D;
class TH2F;
class TH3F;
class TH3;
class TObjArray;
class TGraph;
class TGraphErrors;

class AliTRDcheckESD : public AliAnalysisTaskSE {
public:
  enum ETRDcheckESDstatus {
     kMC        = BIT(0)  // use MC info
    ,kCollision = BIT(1)  // 
  };
  enum ETRDcheckESDhistos {
    kNCl = 1                // number of clusters per track
   ,kTRDstat                // TRD tracks status
   ,kTRDmom                 // TRD track momentum
   ,kPtRes                  // Pt resolution @ vertex for TRD
   ,kTPCVertex              // event vertex from TPC
   ,kEventVertex            // event vertex
   ,kNTracksAll             // ntracks - all
   ,kNTracksAcc             // ntracks - inside acc. and DCA cut
   ,kNTracksTPC             // additional cut on number of TPC clusters
   ,kDCAxy                  // transverse DCA 
   ,kDCAz                   // z - DCA
   ,kPt1                    // Pt distribution, eta and ptmin cuts
   ,kPt2                    // Pt distribution, cuts from kPt1 and DCA cuts
   ,kPt3pos                 // Pt distribution, cuts from kPt2 and cut on TPC clusters for positives (>100)
   ,kPt3neg                 // Pt distribution, cuts from kPt2 and cut on TPC clusters for negatives (>100)
   ,kPt4pos                 // Pt distribution, cuts from kPt3pos and at least one TRD tracklet
   ,kPt4neg                 // Pt distribution, cuts from kPt3neg and at least one TRD tracklet
   ,kTheta                  // distribution of theta for tracks passing the cuts from kPt4pos and kPt4neg
   ,kPhi                    // distribution of phi for tracks passing the cuts from kPt4pos and kPt4neg
   ,kNTPCCl                 // number of TPC clusters, cuts from kPt2
   ,kNTPCCl2                // number of TPC clusters, cuts from kPt2 + pt>1 GeV/c
   ,kTPCDedx                // TPC dE/dx, cuts from kPt3pos or kPt3neg
   ,kEtaPhi                 // (eta,phi) distrib. for tracks after the cuts from kPt3pos or kPt3neg
   ,kEtaNclsTPC             // (TPC_Ncls,eta) distrib. for tracks after the cuts from kPt3pos or kPt3neg
   ,kPhiNclsTPC             // (TPC_Ncls,phi) distrib. for tracks after the cuts from kPt3pos or kPt3neg
   ,kSPDMult                // SPD multiplicity
   ,kNTrackletsTRD          // (TRD tracklets per track, P) distribution, after cuts from kPt4pos or kPt4neg
   ,kNClsTrackTRD=kNTrackletsTRD+6    // (TRD clusters per track, P) distribution, after cuts from kPt4pos or kPt4neg
   ,kPHSlice=kNClsTrackTRD+6         // (slicePH,sliceNo) distribution, after cuts from kPt4pos or kPt4neg
   ,kPHSliceTPCpions=kPHSlice+6      // (slicePH,sliceNo) distribution for TPC pions, after cuts from kPt4pos or kPt4neg 
   ,kTPCdedxPions=kPHSliceTPCpions+6     // (TPC dedx,P) for selected TPC pions
   ,kPHSliceTPCelectrons=kTPCdedxPions+6    // (slicePH,sliceNo) distribution for TPC electrons, after cuts from kPt4pos or kPt4neg
   ,kTPCdedxElectrons=kPHSliceTPCelectrons+6   // (TPC dedx,P) for selected TPC electrons
   ,kQtotP=kTPCdedxElectrons+6              // (total Q from slices, momentum) distribution, after cuts from kPt4pos or kPt4neg
   ,kPropagXYvsP=kQtotP+6       // (X,Y,momentum) distribution after AliESDtrack::PropagateTo(r=300.)
   ,kPropagRZvsP            // (R,Z,momentum) distribution after AliESDtrack::PropagateTo(r=300.)
   ,kTPCRefTracksPos        // (eta,detector phi,Pt) distribution of reference TPC positive tracks (fulfill cuts from kPt3pos)
   ,kTPCRefTracksNeg=kTPCRefTracksPos+6    // (eta,detector phi,Pt) distribution of reference TPC negative tracks (fulfill cuts from kPt3neg)
   ,kTRDRefTracksPos=kTPCRefTracksNeg+6        // (eta,detector phi,Pt) distribution of reference TRD positive tracks (fulfill cuts from kPt4pos)
   ,kTRDRefTracksNeg=kTRDRefTracksPos+6        // (eta,detector phi,Pt) distribution of reference TRD negative tracks (fulfill cuts from kPt4neg)
   ,kTRDRefTracksPos4=kTRDRefTracksNeg+6        // (eta,detector phi,Pt) distribution of reference TRD positive tracks with 4 tracklets (fulfill cuts from kPt4pos)
   ,kTRDRefTracksNeg4=kTRDRefTracksPos4+6        // (eta,detector phi,Pt) distribution of reference TRD negative tracks with 4 tracklets (fulfill cuts from kPt4neg)
   ,kTRDRefTracksPos5=kTRDRefTracksNeg4+6
   ,kTRDRefTracksNeg5=kTRDRefTracksPos5+6
   ,kTRDRefTracksPos6=kTRDRefTracksNeg5+6
   ,kTRDRefTracksNeg6=kTRDRefTracksPos6+6
   ,kTRDEtaPhiAvNtrkl=kTRDRefTracksNeg6+6       // (eta, detector phi) profile of average number of tracklets
   ,kTRDEtaDeltaPhiAvNtrkl=kTRDEtaPhiAvNtrkl+6  // (eta, delta-phi) profile of average number of tracklets
                                         // delta-phi is the angle made by the track with the normal to the chamber entrance plane
   ,kTRDEtaPhiAvQtot=kTRDEtaDeltaPhiAvNtrkl+6      // (eta, detector phi) profile of total tracklet charge from slices			    
   ,kNhistos = kTRDEtaPhiAvQtot+36      // number of histograms
   ,kNrefs   = 4  // number of reference plots
  };
  enum ETrdCfVariables {
    kEventVtxZ=0,
    kEventMult,
    kEventBC,
    kTrackTOFdeltaBC,
    kTrackCharge,
    kTrackPhi,
    kTrackEta,
    kTrackPt,
    kTrackP,
    kTrackTrdTracklets,
    kTrackTrdClusters,
    kTrackQtot,
    kNTrdCfVariables=kTrackQtot+6
  };
  enum ETRDcheckESDbits {
    kTPCout = 1 // track left TPC
   ,kTRDin      // track reach TRD fiducial volume
   ,kTRDout     // track reconstructed in TRD
   ,kTRDpid     // PID calculated in TRD
   ,kTRDref     // track refitted in TRD
  };
  AliTRDcheckESD();
  AliTRDcheckESD(char* name);
  virtual ~AliTRDcheckESD();
  
  void          UserCreateOutputObjects();
  Bool_t        GetRefFigure(Int_t ifig);
  Int_t         GetNRefFigures() const  { return fNRefFigures; } 
  void          UserExec(Option_t *);

  void          SetRefTrackFilter(AliAnalysisCuts* const filter) {fReferenceTrackFilter = filter;}
  
  Bool_t        HasMC() const { return TESTBIT(fStatus, kMC);}
  Bool_t        IsCollision() const {return TESTBIT(fStatus, kCollision);}
  void          SetCollision(Bool_t set=kTRUE) {set ? SETBIT(fStatus, kCollision) : CLRBIT(fStatus, kCollision);}
  TObjArray*    Histos();
  AliCFContainer* GetCFContainer() {return fCfContainer;}
  Bool_t        Load(const Char_t *fn="AnalysisResults.root", const Char_t *dir="TRD_Performance", const Char_t *name=NULL);
  void          SetMC(Bool_t mc = kTRUE) { mc ? SETBIT(fStatus, kMC) : CLRBIT(fStatus, kMC);}
  Bool_t        PutTrendValue(const Char_t *name, Double_t val);
  void          Terminate(Option_t *);
  void          MakeSummary();
  void          MakeSummaryFromCF();

private:
  static const Float_t fgkxTPC; // end radial position of TPC
  static const Float_t fgkxTOF; // start radial position of TOF
  static const UChar_t fgkNgraph[kNrefs]; // number of graphs/ref plot

  void PlotTrackingSummary(Int_t centralityClass=1);     // 1 <= centralityClass <= 5; 0-all centrality classes together
  void PlotPidSummary(Int_t centralityClass=1);     // 1 <= centralityClass <= 5; 0-all centrality classes together
  void PlotCentSummary();  // centrality dependent plots

  void PlotTrackingSummaryFromCF(Int_t centralityClass=1);     // 1 <= centralityClass <= 5; 0-all centrality classes together
  void PlotPidSummaryFromCF(Int_t centralityClass=1);     // 1 <= centralityClass <= 5; 0-all centrality classes together
  void PlotCentSummaryFromCF();                            // centrality dependent plots
    
  AliTRDcheckESD(const AliTRDcheckESD&);
  AliTRDcheckESD& operator=(const AliTRDcheckESD&);
  Int_t         Pdg2Idx(Int_t pdg) const;
  void          Process(TH1 **h, TGraphErrors *g);
  void          Process2D(TH2 * const h, TGraphErrors **g);
  void          PrintStatus(ULong_t s);
  TH2F*         Proj3D(TH3* hist, TH2* accMap, Int_t binLow, Int_t binHigh, Float_t &entries);
  TH1D*         Proj2D(TH2* hist);
  TH1F*         EfficiencyTRD(TH3* tpc3D, TH3* trd3D, Bool_t useAcceptance=kTRUE);
  void          DrawTRDGrid();
  void          SetStyle(TH1* hist, Int_t lineStyle, Int_t lineColor, Int_t lineWidth,
                         Int_t markerStyle, Int_t markerColor, Int_t markerSize);
  void          CheckActiveSM(TH1D* phiProj, Bool_t activeSM[18]);
  
  Int_t            fStatus;            // bit mask for controlling the task
  Int_t            fNRefFigures;       // number of current ref plots
  AliESDEvent      *fESD;              //! ESD event
  AliMCEvent       *fMC;               //! MC event
  AliESDpid        *fESDpid;           //  ESD pid object 
  TObjArray        *fHistos;           //! QA histos
  TObjArray        *fResults;          // QA graphs
  static FILE      *fgFile;            //! trend file streamer
  
  AliCFContainer*  fCfContainer;       // CF container for computing efficiencies
  AliAnalysisCuts* fReferenceTrackFilter;     // reference track filter
  
  // Vertex selection
  static const Float_t fgkEvVertexZ;// cm
  static const Int_t   fgkEvVertexN;// cm
  // Track selection
  static const Float_t fgkTrkDCAxy; // cm
  static const Float_t fgkTrkDCAz;  // cm
  static const Int_t   fgkNclTPC;   // N clusters TPC
  static const Float_t fgkPt;       // min. pt
  static const Float_t fgkEta;      // eta range
  
  static const Float_t fgkQs;      // scale for the total charge

  ClassDef(AliTRDcheckESD, 7)          // user oriented TRD analysis based on ESD-MC data
};
#endif
