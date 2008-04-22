// @(#) $Id$
// Original: AliHLTConfMapFit.h,v 1.5 2004/07/05 09:03:11 loizides 

#ifndef ALIHLTTPCCONFMAPFIT_H
#define ALIHLTTPCCONFMAPFIT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCConfMapFit.h
    @author Anders Vestbo, maintained by Matthias Richter
    @date   
    @brief  Fit class for conformal mapping tracking.
*/

class AliHLTTPCConfMapTrack;
class AliHLTTPCVertex;

/** 
 * @class AliHLTTPCConfMapFit
 *
 * Fit class for conformal mapping tracking
 *
 * @ingroup alihlt_tpc
 */
class AliHLTTPCConfMapFit {

 public:
  /** default constructor */
  AliHLTTPCConfMapFit();
  /** constructor */
  AliHLTTPCConfMapFit (AliHLTTPCConfMapTrack *track,AliHLTTPCVertex *vertex);
  /** destructor */
  virtual ~AliHLTTPCConfMapFit();

  // helix fit
  Int_t FitHelix();
  Int_t FitCircle();
  Int_t FitLine();

  // straight line fit
  Int_t FitStraightLine();
  Int_t FitLineXY();
  Int_t FitLineSZ();
  
 private:
  /** copy constructor prohibited */
  AliHLTTPCConfMapFit(const AliHLTTPCConfMapFit&);
  /** assignment operator prohibited */
  AliHLTTPCConfMapFit& operator=(const AliHLTTPCConfMapFit&);

  AliHLTTPCConfMapTrack *fTrack; //!
  AliHLTTPCVertex *fVertex; //!
  
  ClassDef(AliHLTTPCConfMapFit,1) //Conformal mapping fit class
};

#endif // ALIHLTTPCCONFMAPFIT_H
