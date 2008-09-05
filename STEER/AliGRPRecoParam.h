#ifndef ALIGRPRECOPARAM_H
#define ALIGRPRECOPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class with global reconstruction parameters                               //
// (initially, parameters for AliVertexerTracks)                             //
// Origin: andrea.dainese@lnl.infn.it                                        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliDetectorRecoParam.h"

class AliGRPRecoParam : public AliDetectorRecoParam
{
 public: 
  AliGRPRecoParam();
  virtual ~AliGRPRecoParam();

  static AliGRPRecoParam *GetLowFluxParam();// make reco parameters for low flux env.
  static AliGRPRecoParam *GetHighFluxParam();// make reco parameters for high flux env. 

  void  SetMostProbablePt(Double_t pt=0.350) { fMostProbablePt=pt; return; }
  Double_t GetMostProbablePt() const { return fMostProbablePt; }

  void  SetVertexerTracksCuts(Int_t mode,Int_t ncuts,Double_t cuts[10]);
  void  SetVertexerTracksCutsITS(Int_t ncuts,Double_t cuts[10])
    { SetVertexerTracksCuts(0,ncuts,cuts); return; }
  void  SetVertexerTracksCutsTPC(Int_t ncuts,Double_t cuts[10])
    { SetVertexerTracksCuts(1,ncuts,cuts); return; }
  Int_t GetVertexerTracksNCuts() const { return fVertexerTracksNCuts; }
  void  GetVertexerTracksCuts(Int_t mode,Double_t *cuts) const;
  void  GetVertexerTracksCutsITS(Double_t *cuts) const
    { GetVertexerTracksCuts(0,cuts); return; }
  void  GetVertexerTracksCutsTPC(Double_t *cuts) const
    { GetVertexerTracksCuts(1,cuts); return; }

  AliGRPRecoParam(const AliGRPRecoParam&);
  AliGRPRecoParam& operator=(const AliGRPRecoParam&);

 protected:
  //

  Double_t fMostProbablePt; // to be used for B=0 tracking
  Int_t    fVertexerTracksNCuts; // number of cuts for AliVertexerTracks
  // cuts for AliVertexerTracks: ITS mode
  Double_t fVertexerTracksITSdcacut; // general dca
  Double_t fVertexerTracksITSdcacutIter0; // dca in iteration 0
  Double_t fVertexerTracksITSmaxd0z0; // max d0z0
  Double_t fVertexerTracksITSminCls; // min clusters
  Double_t fVertexerTracksITSmintrks; // min tracks
  Double_t fVertexerTracksITSnsigma; // n sigma for d0 cut
  Double_t fVertexerTracksITSnindetfitter; // min det to try inversion
  Double_t fVertexerTracksITSmaxtgl; // max tgl 
  Double_t fVertexerTracksITSfidR; // fiducial radius
  Double_t fVertexerTracksITSfidZ; // fiducial z

  // cuts for AliVertexerTracks: TPC-only mode
  Double_t fVertexerTracksTPCdcacut; // general dca
  Double_t fVertexerTracksTPCdcacutIter0; // dca in iteration 0
  Double_t fVertexerTracksTPCmaxd0z0; // max d0z0
  Double_t fVertexerTracksTPCminCls; // min clusters
  Double_t fVertexerTracksTPCmintrks; // min tracks
  Double_t fVertexerTracksTPCnsigma; // n sigma for d0 cut
  Double_t fVertexerTracksTPCnindetfitter; // min det to try inversion
  Double_t fVertexerTracksTPCmaxtgl; // max tgl 
  Double_t fVertexerTracksTPCfidR; // fiducial radius
  Double_t fVertexerTracksTPCfidZ; // fiducial z

  ClassDef(AliGRPRecoParam,2) // global reco parameters
};

#endif
