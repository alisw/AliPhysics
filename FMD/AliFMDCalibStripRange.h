#ifndef ALIFMDCALIBSTRIPRANGE_H
#define ALIFMDCALIBSTRIPRANGE_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//____________________________________________________________________
//                                                                          
// This class stores which strips are read-out.  
// In principle this can be set for each half-ring.   
// However, in real life, all the detectors will probably read out all
// strips, and dead areas can be handled off-line. 
// This information comes from DCS or the like.
//
/** @file    AliFMDCalibStripRange.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:32:14 2006
    @brief   Per digitizer card pulser calibration
    @ingroup FMD_base
*/
#include <iosfwd>
#ifndef ROOT_TObject
# include <TObject.h>
#endif
#ifndef ALIFMDUSHORTMAP_H
# include "AliFMDUShortMap.h"
#endif
#ifndef ROOT_TArrayI
# include <TArrayI.h>
#endif

//____________________________________________________________________
/** @brief Per digitizer card pulser calibration
    @ingroup FMD_base
*/
class AliFMDCalibStripRange : public TObject
{
public:
  /** CTOR */
  AliFMDCalibStripRange();
  /** Copy CTOR
      @param o Object to copy from  */
  AliFMDCalibStripRange(const AliFMDCalibStripRange& o);
  /** Assignment operator 
      @param o Object to assign from 
      @return Reference to assign from  */
  AliFMDCalibStripRange& operator=(const AliFMDCalibStripRange& o);
  /** Set sample for a DDL
      @param det   Detector #
      @param ring  Ring ID 
      @param sec   Sector # 
      @param str   Strip number (not used)
      @param min   Minimum strip (0-127)
      @param max   Maximum strip (0-127) */
  void Set(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, 
	   UShort_t min, UShort_t max);
  /** Get minimum strip read out (0-127)
      @param det  Detector #
      @param ring Ring ID 
      @param sec  Sector # 
      @param str  Strip number (not used)
      @return Minimum strip  */
  UShort_t Min(UShort_t det, Char_t ring, UShort_t sec, UShort_t str=0) const;
  /** Get maximum strip read out (0-127)
      @param det  Detector #
      @param ring Ring ID 
      @param sec  Sector # 
      @param str  Strip number (not used)
      @return Maximum strip  */
  UShort_t Max(UShort_t det, Char_t ring, UShort_t sec, UShort_t str=0) const;
  /**
     Dump stored strip ranges to file passed as ofstream
     @param outFile Outputfile
  */
  void WriteToFile(std::ostream &, Bool_t* detectors=0);
  /**
     Read information from file and set values
     @param inFile inputFile
   */
  void ReadFromFile(std::istream &);
  
  const AliFMDUShortMap& Ranges() const { return fRanges; }
protected:
  // TArrayI fRates; // Sample rates 
  AliFMDUShortMap fRanges; // Min max 
  ClassDef(AliFMDCalibStripRange,1); // Sample rates 
};

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//


