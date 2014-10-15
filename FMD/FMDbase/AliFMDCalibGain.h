#ifndef ALIFMDCALIBGAIN_H
#define ALIFMDCALIBGAIN_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDCalibGain.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 18:30:16 2006
    @brief   Per strip gain calibration 
*/
// Gain value and width for each strip in the FMD. 
// Foo 
// Bar 
//
#ifndef ALIFMDFLOATMAP_H
# include <AliFMDFloatMap.h>
#endif
#include <iosfwd>
class AliFMDBoolMap;

//____________________________________________________________________
/** @brief Gain value and width for each strip in the FMD 
    @ingroup FMD_base
*/
class AliFMDCalibGain : public TObject 
{
public:
  /** Constructor */
  AliFMDCalibGain();
  /** Destructor */
  ~AliFMDCalibGain() {}
  /** Copy constructor 
      @param o object to copy from */
  AliFMDCalibGain(const AliFMDCalibGain& o);
  /** Assignment operator 
      @param o Object to assign from */
  AliFMDCalibGain& operator=(const AliFMDCalibGain& o);
  /** Set the values for a strip. 
      @param det  Detector 
      @param ring Ring 
      @param sec  Sector 
      @param str  Strip
      @param val  Value of gain */
  void Set(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, Float_t val);
  /** Set the global threshold 
      @param thres Threshold */
  void Set(Float_t thres) { fThreshold = thres; }
  /** Get gain for a strip. 
      @param det  Detector 
      @param ring Ring 
      @param sec  Sector 
      @param str  Strip
      @return Gain for strip */
  Float_t Value(UShort_t det, Char_t ring, UShort_t sec, UShort_t str);

  /**
     Read information from file and set values
     @param inFile inputFile
   */
  Bool_t ReadFromFile(std::istream & inFile);

  /** @return threshold */
  Float_t Threshold() const { return fThreshold; }
  const AliFMDFloatMap& Values() const { return fValue; }
  /** 
   * Make a dead map based on the gain of the channels.  If the gain 
   * of a paraticular channel falls outside the range specified by @a
   * min and @a max, then the channel is marked as dead. 
   *
   * If the argument @a dead is non-null, then the map passed is
   * modified.  That is, channels marked as dead in the map will
   * remain marked.   Channels that meat the criterion (gain outside
   * the specified range) will in addition be marked as dead. 
   *
   * If the argument @a dead is null, then a new map is created and a
   * pointer to this will be returned. 
   * 
   * @param min Minimum value of gain for a channel before it is
   * @param max Maximum value of gain for a channel before it is
   * marked as dead. 
   * @param dead If non-null, then modify this map. 
   * 
   * @return A pointer to possibly newly allocated dead map. 
   */
  AliFMDBoolMap* MakeDeadMap(Float_t min, Float_t max, 
			     AliFMDBoolMap* dead=0) const;
private:
  AliFMDFloatMap fValue;       // Map
  Float_t        fThreshold;   // Global threshold
  ClassDef(AliFMDCalibGain, 1) // Gain data for the FMD 
};


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//


