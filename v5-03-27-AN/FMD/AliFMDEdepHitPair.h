#ifndef ALIFMDEDEPHITPAIR_H
#define ALIFMDEDEPHITPAIR_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDEdepHitPair.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:50 2006
    @brief   Per strip map of energy deposited and number of hits 
    @ingroup FMD_sim
*/
//____________________________________________________________________
//                                                                          
// Contains a pair of energy deposited fEdep and number of hits  
// fN, fEdep is the summed energy deposition, and fN is the
// number of hits.  The map contains one such object or each strip.
// It is used to cache the data in the digitization classes
// AliFMDBaseDigitizer and so on. 
//
#ifndef ROOT_Rtypes
# include <Rtypes.h>
#endif 
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif

//____________________________________________________________________
/** @brief Cache of Energy deposited, hit information per strip.
    Contains a pair of energy deposited @c fEdep and 
    number of hits @c fN, @c fEdep is the summed energy deposition,
    and @c fN is the number of hits 
    @ingroup FMD_sim
*/
class AliFMDEdepHitPair 
{
public:
  Float_t  fEdep;   // summed energy deposition
  UShort_t fN;      // Number of hits
  UShort_t fNPrim;  // Number of primaries;
  TArrayI  fLabels; // Track labels.
  
  /** CTOR  */
  AliFMDEdepHitPair() : fEdep(0), fN(0), fNPrim(0), fLabels(0) {}
  /** DTOR */
  virtual ~AliFMDEdepHitPair() {}
  /** Assignment operator 
      @param o Object to assign from 
      @return Reference to this object */
  AliFMDEdepHitPair& operator=(const AliFMDEdepHitPair& o) 
  { 
    if (&o == this) return *this; 
    fEdep   = o.fEdep; 
    fN      = o.fN; 
    fNPrim  = o.fNPrim;
    fLabels = o.fLabels;
    return *this; 
  }
  /** Copy CTOR 
      @param o Object to copy from */
  AliFMDEdepHitPair(const AliFMDEdepHitPair& o) 
    : fEdep(o.fEdep), fN(o.fN), fNPrim(o.fNPrim), fLabels(o.fLabels)
  {}
  ClassDef(AliFMDEdepHitPair, 3)
};

#endif 
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//


