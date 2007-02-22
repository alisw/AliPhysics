// @(#) $Id$
// Original: AliHLTConfMapFit.h,v 1.5 2004/07/05 09:03:11 loizides 

#ifndef ALIHLTTPCCONFMAPFIT_H
#define ALIHLTTPCCONFMAPFIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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
 */
class AliHLTTPCConfMapFit {

 public:
  /** default constructor */
  AliHLTTPCConfMapFit();
  /** constructor */
  AliHLTTPCConfMapFit (AliHLTTPCConfMapTrack *track,AliHLTTPCVertex *vertex);
  /** not a valid copy constructor, defined according to effective C++ style */
  AliHLTTPCConfMapFit(const AliHLTTPCConfMapFit&);
  /** not a valid assignment op, but defined according to effective C++ style */
  AliHLTTPCConfMapFit& operator=(const AliHLTTPCConfMapFit&);
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
  AliHLTTPCConfMapTrack *fTrack; //!
  AliHLTTPCVertex *fVertex; //!
  
  ClassDef(AliHLTTPCConfMapFit,1) //Conformal mapping fit class
};

#endif // ALIHLTTPCCONFMAPFIT_H
