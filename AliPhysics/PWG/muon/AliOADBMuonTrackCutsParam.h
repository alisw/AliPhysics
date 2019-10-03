#ifndef AliOADBMuonTrackCutsParam_H
#define AliOADBMuonTrackCutsParam_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     OADB interface for the Muon track cuts
//     Author: Diego Stocco
//    
// This class contains the parameters needed to configure the
// muon track cuts
// -------------------------------------------------------------------------

#include <TNamed.h>

class TVector3;

class AliOADBMuonTrackCutsParam : public TNamed {

 public :
  AliOADBMuonTrackCutsParam ();
  AliOADBMuonTrackCutsParam ( const char* name );
  virtual ~AliOADBMuonTrackCutsParam ();
  AliOADBMuonTrackCutsParam ( const AliOADBMuonTrackCutsParam& other ); 
  AliOADBMuonTrackCutsParam& operator = ( const AliOADBMuonTrackCutsParam& other );
  
  void SetMeanDCA ( Double_t xAtDca, Double_t yAtDca, Double_t zAtDca = 0. );
  TVector3 GetMeanDCA () const;
  
  void SetMeanPCorr ( Double_t pCorrThetaAbs23, Double_t pCorrThetaAbs310 ); 
  Double_t GetMeanPCorr23 ( ) const;
  Double_t GetMeanPCorr310 ( ) const;
  
  void SetSigmaPdca ( Double_t sigmaThetaAbs23, Double_t sigmaThetaAbs310 );
  Double_t GetSigmaPdca23 ( ) const;
  Double_t GetSigmaPdca310 ( ) const;
  
  void SetNSigmaPdca ( Double_t nSigmas );
  Double_t GetNSigmaPdca () const;
  
  void SetChi2NormCut ( Double_t chi2normCut );
  Double_t GetChi2NormCut () const;
  
  void SetRelPResolution ( Double_t relPResolution );
  Double_t GetRelPResolution () const;
  
  void SetSlopeResolution ( Double_t slopeResolution );
  Double_t GetSlopeResolution () const;
  
  void SetSharpPtCut ( Double_t valueApt, Double_t valueLpt, Double_t valueHpt );
  Double_t GetSharpPtCut ( Int_t trigPtCut, Bool_t warn = kTRUE ) const;
  
  void Print ( Option_t* option = "" ) const;

 private :
  Double_t fMeanDcaX;        ///< Average track DCA_x
  Double_t fMeanDcaY;        ///< Average track DCA_y
  Double_t fMeanDcaZ;        ///< Average track DCA_z
  Double_t fMeanPCorr23;     ///< Average momentum correction in 2-3 deg
  Double_t fMeanPCorr310;    ///< Average momentum correction in 3-10 deg
  Double_t fSigmaPdca23;     ///< Sigma_PxDCA in 2-3 deg
  Double_t fSigmaPdca310;    ///< Sigma_PxDCA in 3-10 deg
  Double_t fNSigmaPdcaCut;   ///< Cut value in units of sigma_PxDCA
  Double_t fChi2NormCut;     ///< Cut on the normalized chi2 of track
  Double_t fRelPResolution;  ///< Relative momentum resolution
  Double_t fSlopeResolution; ///< Slope resolution
  Double_t fSharpPtApt;      ///< Sharp tracker pt cut for Apt
  Double_t fSharpPtLpt;      ///< Sharp tracker pt cut for Lpt
  Double_t fSharpPtHpt;      ///< Sharp tracker pt cut for Hpt
  
  ClassDef(AliOADBMuonTrackCutsParam, 1);
};

#endif
