#ifndef ALI_FMD_PREPROCESSOR_H
#define ALI_FMD_PREPRECESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
//___________________________________________________________________
// The class processes data points from DCS (via Amanada), and DAQ DA
// files (via FXS) to make calibration data for the FMD.
// More to come. 


#include "AliPreprocessor.h"
#ifndef ALIFMDUSHORTMAP_H
# include <AliFMDUShortMap.h>
#endif
#ifndef ALIFMDBOOLMAP_H
# include <AliFMDBoolMap.h>
#endif
typedef AliFMDUShortMap AliFMDCalibZeroSuppression;
typedef AliFMDBoolMap  AliFMDCalibDeadMap;
class AliFMDCalibPedestal;
class AliFMDCalibGain;
class AliFMDCalibSampleRate;
class AliFMDCalibStripRange;
class AliCDBEntry;
class TList;

//___________________________________________________________________
/** The class processes data points from DCS (via Amanada), and DAQ DA
    files (via FXS) to make calibration data for the FMD.

    Data points: 
       *  Nothing yet. 

    DAQ FXS file:
       * pedestals - a (ASCII) Comma Separated Values files with the
                     fields 
                          rcu	 DDL number 
                          board   FEC board number 
                          chip    ALTRO chip number on FEC
                          channel ALTRO channel number
                          strip   VA1 strip number
                          sample  Sample number
                          ped     Mean of ADC spectra
                          noise   Spread of ADC spectra
                          mu      Mean of Gaussian fit to ADC spectra
                          sigma   Variance of Gaussian fit to ADC spectra
                          chi2    Chi^2 per degrees of freedom of fit
       * Gains     - a (ASCII) Comma Separated Values files with the
                     fields 
                          rcu	 DDL number 
                          board   FEC board number 
                          chip    ALTRO chip number on FEC
                          channel ALTRO channel number
                          strip   VA1 strip number
                          gain    Slope of gain
                          error   Error on gain
                          chi2    Chi^2 per degrees of freedom of fit
*/
class AliFMDPreprocessor: public AliPreprocessor 
{
public:
  /** Constructor */
  AliFMDPreprocessor(): AliPreprocessor("FMD",0) {}
  /** Constructor 
      @param shuttle Shuttle */
  AliFMDPreprocessor(AliShuttleInterface* shuttle) ;
  
  
  /** Destructor */
  virtual ~AliFMDPreprocessor() {}
  /** Get an entry from OCDB */ 
  AliCDBEntry* GetFromCDB(const char* second, const char* third);
protected:
  /** 
   * Get the pedestal calibrations 
   *
   * @param list List of files 
   */
  AliFMDCalibPedestal*   GetPedestalCalibration(const TList* list);
  /** 
   * Get the gain calibrations 
   * 
   * @param list List of files 
   */
  AliFMDCalibGain*       GetGainCalibration(const TList*);
  /** 
   * Get the dead channels map based on the pedestal
   * an gain calibration objects. 
   */
  AliFMDCalibDeadMap*    GetDeadChannelMap(AliFMDCalibPedestal* pedcalib,
					   AliFMDCalibGain*     gaincalib);

  /** Get info calibrations. 
      @param files List of files. 
      @param s     On return, newly allocated object 
      @param r     On return, newly allocated object 
      @param z     On return, newly allocated object 
      @return kTRUE on success */
  Bool_t GetInfoCalibration(TList* files, 
			    AliFMDCalibSampleRate*&      s,
			    AliFMDCalibStripRange*&      r, 
			    AliFMDCalibZeroSuppression*& z);
  /** Convinience function 
      @param list On return, list of files. 
      @param system Alice system (DAQ, DCS, ...)
      @param id     File id
      @return kTRUE on success. */
  Bool_t GetAndCheckFileSources(TList*& list, Int_t system, const char* id);
  
  /** Entry method 
      @param dcsAliasMap Map of DCS data points */
  virtual UInt_t Process(TMap* dcsAliasMap);
private:
  ClassDef(AliFMDPreprocessor, 1)
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
