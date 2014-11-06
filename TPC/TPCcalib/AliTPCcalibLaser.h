#ifndef ALITPCCALIBLASER_H
#define ALITPCCALIBLASER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "TObject.h"
#include "TObjArray.h"
#include "TLinearFitter.h"
#include "AliTPCcalibBase.h"
#include "TH1.h"
#include "TH2F.h"
#include "THnSparse.h"


class AliExternalTrackParam;
//class AliESDtrack;
//class AliESDEvent;
//class AliESDfriend;
class AliVEvent;
class AliVTrack;
class AliVfriendEvent;
class TGraphErrors;
class TTree;
class TH2F;
class AliTPCLaserTrack;
class TCut;

class AliTPCcalibLaser:public AliTPCcalibBase {
public:
  AliTPCcalibLaser();
  AliTPCcalibLaser(const Text_t *name, const Text_t *title, Bool_t full=kTRUE);
  AliTPCcalibLaser(const AliTPCcalibLaser& laser);
  AliTPCcalibLaser & operator=(const AliTPCcalibLaser& calibLaser);
  virtual ~AliTPCcalibLaser();
  virtual void     Process(AliVEvent *event);
  Int_t   GetNtracks(){return fNtracks;}
  virtual void Analyze();
  static void        DumpLaser(const char *finput, Int_t run);
  static void        FitLaserClusters(Int_t run);
  virtual Long64_t Merge(TCollection *li);
  virtual void DumpMeanInfo(Int_t run=-1);
  static  void DumpScanInfo(TTree * tree, const char * cutUser="entries>300&&(gz2<0.15&&gphi2<0.1&&gp42<0.02&&abs(gp41)<0.03)");
  static  void DumpFitInfo(TTree * chainFit, Int_t id);
  static  TH1* GetLaserProjection(TH2F* his, Int_t laser){return his->ProjectionY("aaa",laser+1,laser+1);}
  //
  //
  virtual void DumpLaser(Int_t id);
  virtual void RefitLaserJW(Int_t id);
  void         FitDriftV();
  Bool_t       FitDriftV(Float_t minFraction);
  //
  void         MakeDistHisto(Int_t id);
  void         AddCut(Double_t xcut, Double_t ycut, Double_t ncl){fEdgeXcuts[fNcuts]=xcut; fEdgeYcuts[fNcuts]=ycut; fNClCuts[fNcuts]=ncl; fNcuts++;}

  Int_t  FindMirror(AliVTrack *track, AliTPCseed *seed);
  Bool_t AcceptLaser(Int_t id);
  Float_t GetDistance(AliExternalTrackParam *track, AliTPCLaserTrack *ltrp);
  void   MakeFitHistos();
  void   UpdateFitHistos();
  void   MergeFitHistos(AliTPCcalibLaser * add);
  //void     Process(AliESDtrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliVTrack *track, Int_t runNo=-1){AliTPCcalibBase::Process(track,runNo);}
  void     Process(AliTPCseed *track){return AliTPCcalibBase::Process(track);}
  //
  void SetBeamParameters(TVectorD& meanOffset, TVectorD& meanSlope,
			 TVectorD& sectorArray, Int_t option);
  void SetFixedDriftVConstant(Double_t aside0, Double_t aside1,
			      Double_t cside0, Double_t cside1)
    {fUseFixedDriftV = 1; fFixedFitAside0=aside0; fFixedFitAside1=aside1;
    fFixedFitCside0=cside0; fFixedFitCside1=cside1;}

  AliVEvent  * fEvent;             //! ESD event  - not OWNER
  AliVfriendEvent * fEventFriend;       //! ESD event  - not OWNER
  Int_t          fNtracks;         //! counter of associated laser tracks
  //
  TObjArray      fTracksMirror;    //! tracks with mirror information
  TObjArray      fTracks;       //! tracks with reconstructed information -
  //                               not owner ESD
  TObjArray      fTracksEsdParam;  //! tracks with reconstructed information - 
  //                               is owner ESD at mirror
  TObjArray      fTracksTPC;       //! tracks with reconstructed information - TPC
  Int_t          fCounter[336];    //! counter of usage
  Float_t        fClusterCounter[336]; //!couter of clusters in "sensitive are"
  Float_t        fClusterSatur[336];   //!couter of saturated clusters in "sensitive are"
  Bool_t         fFullCalib;            // do full calibrration
  Float_t        fFitZ[336];           //fitted z position
  //
  TObjArray      fDeltaZ;          //-> array of histograms of delta z for each track
  TObjArray      fDeltaP3;         //-> array of histograms of P3      for each track
  TObjArray      fDeltaP4;         //-> array of histograms of P4      for each track
  TObjArray      fDeltaPhi;        //-> array of histograms of delta Phi for each track
  TObjArray      fDeltaPhiP;       //-> array of histograms of delta Phi direction for each track
  TObjArray      fSignals;         //->Array of dedx signals

  //
  // Refit residuals histogram
  //
  THnSparseS     *fHisLaser;      //  N dim histogram of laser 
  THnSparseS     *fHisLaserPad;   //  N dim histogram of laser 
  THnSparseS     *fHisLaserTime;   //  N dim histogram of laser 
  //
  TH2F           *fHisNclIn;      //->Number of clusters inner
  TH2F           *fHisNclOut;     //->Number of clusters outer
  TH2F           *fHisNclIO;      //->Number of cluster inner outer
  TH2F           *fHisLclIn;      //->Level arm inner
  TH2F           *fHisLclOut;     //->Level arm outer
  TH2F           *fHisLclIO;      //->Level aram inner outer

  TH2F           *fHisdEdx;       //->dEdx histo
  TH2F           *fHisdZfit;      //->distance to the mirror after linear fit
  //
  //
  TH2F           *fHisChi2YIn1;      //->chi2 y inner - line
  TH2F           *fHisChi2YOut1;     //->chi2 y inner - line
  TH2F           *fHisChi2YIn2;      //->chi2 y inner - parabola
  TH2F           *fHisChi2YOut2;     //->chi2 y inner - parabola
  TH2F           *fHisChi2YIO1;      //->chi2 y IO    - common
  TH2F           *fHisChi2ZIn1;      //->chi2 z inner - line
  TH2F           *fHisChi2ZOut1;     //->chi2 z inner - line
  TH2F           *fHisChi2ZIn2;      //->chi2 z inner - parabola
  TH2F           *fHisChi2ZOut2;     //->chi2 z inner - parabola
  TH2F           *fHisChi2ZIO1;      //->chi2 z IO    - common
  //
  //
  TH2F           *fHisPy1vP0;     //-> delta y   P0outer-P0inner - line
  TH2F           *fHisPy2vP0;     //-> delta y   P0outer-P0inner - parabola
  TH2F           *fHisPy3vP0;     //-> delta y   P0outer-P0inner - common parabola
  TH2F           *fHisPy1vP1;     //-> delta ky  P1outer-P1inner - line
  TH2F           *fHisPy2vP1;     //-> delta ky  P1outer-P1inner - parabola
  TH2F           *fHisPy3vP1;     //-> delta ky  P1outer-P1inner - common parabola
  TH2F           *fHisPy2vP2In;   //-> Curv  P2inner - parabola
  TH2F           *fHisPy2vP2Out;  //-> Curv  P2outer - parabola
  TH2F           *fHisPy3vP2IO;   //-> Curv  P2outerinner - common parabola
  //
  //
  TH2F           *fHisPz1vP0;     //-> delta z   P0outer-P0inner - line
  TH2F           *fHisPz2vP0;     //-> delta z   P0outer-P0inner - parabola
  TH2F           *fHisPz3vP0;     //-> delta z   P0outer-P0inner - common parabola
  TH2F           *fHisPz1vP1;     //-> delta kz  P1outer-P1inner - line
  TH2F           *fHisPz2vP1;     //-> delta kz  P1outer-P1inner - parabola
  TH2F           *fHisPz3vP1;     //-> delta kz  P1outer-P1inner - common parabola
  TH2F           *fHisPz2vP2In;   //-> Curv  P2inner - parabola
  TH2F           *fHisPz2vP2Out;  //-> Curv  P2outer - parabola
  TH2F           *fHisPz3vP2IO;   //-> Curv  P2outerinner - common parabola
  //
  // Residual histograms
  //
  TObjArray      fDeltaYres;       //-> array of histograms of delta y residuals for each track
  TObjArray      fDeltaZres;       //-> array of histograms of delta z residuals for each track
  TObjArray      fDeltaYres2;       //-> array of histograms of delta y residuals for each track
  TObjArray      fDeltaZres2;       //-> array of histograms of delta z residuals for each track
  TObjArray      fDeltaYresAbs; //-> array of histograms of absolute delta y residuals for each track
  TH1F           *fHisYAbsErrors; //-> Number of errors (wrongly assigned tracks) per beam
  TObjArray      fDeltaZresAbs; //-> array of histograms of absolute delta z residuals for each track
  TH1F           *fHisZAbsErrors; //-> Number of errors (wrongly assigned tracks or missing drift velocity) per beam
  //  TObjArray      fDeltaYres3;       //-> array of histograms of delta y residuals for each track
  //TObjArray      fDeltaZres3;       //-> array of histograms of delta z residuals for each track

  //
  TVectorD*      fFitAside;        //! drift fit - A side
  TVectorD*      fFitCside;        //! drift fit - C- side
  TVectorD*      fFitACside;        //! drift fit - A+C- side
  //
  TVectorD       fEdgeXcuts;       //! cuts in local x direction; used in the refit of the laser tracks
  TVectorD       fEdgeYcuts;       //! cuts in local y direction; used in the refit of the laser tracks
  TVectorD       fNClCuts;         //! cuts on the number of clusters per tracklet; used in the refit of the laser tracks
  Int_t          fNcuts;           //! number of cuts
  TVectorD       fBeamSectorOuter;  //! sector map for beams in outer sector
  TVectorD       fBeamSectorInner;  //! sector map for beams in inner sector
  TVectorD       fBeamOffsetYOuter; //! absolute y beam offset in outer sector
  TVectorD       fBeamSlopeYOuter;  //! absolute y beam slope  in outer sector
  TVectorD       fBeamOffsetYInner; //! absolute y beam offset in inner sector
  TVectorD       fBeamSlopeYInner;  //! absolute y beam slope  in inner sector
  TVectorD       fBeamOffsetZOuter; //! absolute z beam offset in outer sectror 
  TVectorD       fBeamSlopeZOuter;  //! absolute z beam slope  in outer sector 
  TVectorD       fBeamOffsetZInner; //! absolute z beam offset in inner sectror 
  TVectorD       fBeamSlopeZInner;  //! absolute z beam slope  in inner sector
  Bool_t         fInverseSlopeZ;    //! invert slope in z - mismatch between database and lasers
  Int_t          fUseFixedDriftV;   // flag for fixed drift velocity for abs res
  Double_t       fFixedFitAside0;   // Fixed drift v constant 0 - A side
  Double_t       fFixedFitAside1;   // Fixed drift v constant 1 - A side
  Double_t       fFixedFitCside0;   // Fixed drift v constant 0 - C side
  Double_t       fFixedFitCside1;   // Fixed drift v constant 1 - C side
  //
private:
  ClassDef(AliTPCcalibLaser,6)
};





#endif
