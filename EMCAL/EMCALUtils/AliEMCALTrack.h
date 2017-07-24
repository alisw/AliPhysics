#ifndef AliEMCALTRACK_H
#define AliEMCALTRACK_H

//========================================================================  
// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  
// See cxx source for full Copyright notice                                
//========================================================================  
///                       
/// \class AliEMCALTrack 
/// \ingroup EMCALUtils
/// \brief Matched track to cluster handling.
///
/// Implementation of a track to be used for EMCAL track matching, 
/// a track which is propagated to EMCAL and matches an EMCAL cluster.
///
/// This object is used to find track intersection with EMCAL surface
/// in order to find the most well matched EMCAL cluster to associate to it. 
///
/// NO Kalman-like parameter updating is done, but it allows track propagation 
/// and suitable energy loss correction even in an environment with a variable 
/// magnetic field, which is not well managed in the AliExternalTrackParam class.
///
/// ------------------------------------------------------------------------
/// \author A. Pulvirenti, <alberto.pulvirenti@ct.infn.it>
//=========================================================================

#include "AliExternalTrackParam.h"

class AliESDtrack;

class AliEMCALTrack : public AliExternalTrackParam
{
  
public:

	AliEMCALTrack();
	AliEMCALTrack(const AliESDtrack& t);
	AliEMCALTrack(const AliEMCALTrack& t);
	AliEMCALTrack& operator=(const AliEMCALTrack &t);
	
  Bool_t      PropagateTo(Double_t xr, Double_t d = 0.0, Double_t x0 = 21.82);
	Double_t    StraightPropagateTo(Double_t xk, Double_t &x, Double_t &y, Double_t &z);
	Bool_t      PropagateToGlobal(Double_t x, Double_t y, Double_t z, Double_t d = 0.0, Double_t x0 = 21.82);  
  
	Int_t       Compare(const TObject *o)        const ;
	Int_t       GetClusterIndex()                const { return fClusterIndex     ; }
	Double_t    GetClusterDist()                 const { return fClusterDist      ; }
	Double_t    GetMass()                        const { return fMass             ; }
	Int_t       GetSeedIndex()                   const { return fSeedIndex        ; }
	Int_t       GetSeedLabel()                   const { return fSeedLabel        ; }
	Bool_t      IsSortable()                     const { return kTRUE             ; }
  			
	void        SetClusterIndex(Int_t idx)             { fClusterIndex    = idx   ; }
	void        SetClusterDist(Double_t dist)          { fClusterDist     = dist  ; }
	void        SetMass(Double_t mass)                 { fMass            = mass  ; }
	void        SetSeedIndex(Int_t index)              { fSeedIndex       = index ; }
	void        SetSeedLabel(Int_t label)              { fSeedLabel       = label ; }
		
	static void SetUseOuterParams(Bool_t doit = kTRUE) { fgUseOuterParams = doit  ; }
	static void SetCorrectForEL  (Bool_t doit = kTRUE) { fgCorrectForEL   = doit  ; }
	static void SetSortByPt      (Bool_t doit = kTRUE) { fgSortByPt       = doit  ; }

protected:
	
	Bool_t      SimplePropagation(Double_t xk, Double_t d = 0.0, Double_t x0 = 0.0);
	
  static  Bool_t    fgUseOuterParams;    ///< Use outer parameters from AliESDtrack?
  static  Bool_t    fgCorrectForEL;      ///< Apply corrections for energy loss?
  static  Bool_t    fgSortByPt;          ///< true --> sort by pt, false --> sort by track-cluster distance
          Int_t     fClusterIndex;       ///< Index of matched cluster (if any)
          Double_t  fClusterDist;        ///< Distance between track propagation and matched cluster
          Double_t  fMass;               ///< Mass hypothesis (in GeV/c2)
          Int_t     fSeedIndex;          ///< Index of imported ESD track in its owner AliESD
          Int_t     fSeedLabel;          ///< GEANT label of imported ESD track

private:
	
  /// \cond CLASSIMP
  ClassDef(AliEMCALTrack, 0) ;
  /// \endcond

};

#endif //AliEMCALTRACK_H
