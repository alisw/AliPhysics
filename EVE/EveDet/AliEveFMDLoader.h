// -*- mode: C++ -*- 
//
/* Copyright(c) 2008, Christian Holm Christensen
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliEveFMDLoader.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 17:59:37 2006
    @brief   Declaration of AliEveFMDLoader singleton class 
*/
//____________________________________________________________________
// 
// Forward Multiplicity Detector based on Silicon wafers. 
// . 
// This class is the loader for the event display.
// 
#ifndef ALIEVEFMDLOADER_H
#define ALIEVEFMDLOADER_H
#include <TEveElement.h>
#include <TEveRGBAPalette.h>
#include <TClonesArray.h>
// Forward declarations
class TEveDigitSet;

/** @class AliEveFMDLoader 
    @brief Loader of FMD data for the EVE event display 
    @ingroup FMD_util

    This class is a singleton, meaning that there's only one instance
    of this.  This is done to speed up the processing by putting all
    things that are needed every time into the constructor.
*/
class AliEveFMDLoader : public TEveElementList
{
public:
  /** @{ 
      @name Loading and displaying data */
  /** Load and display hits */
  virtual void LoadHits();
  /** Load and display digits */
  virtual void LoadDigits();
  /** Load and display raw data digits */
  virtual void LoadRaw();
  /** Load and display ESD */
  virtual void LoadESD();
  /** @} */

  /** @{ 
      @name Hacks to optimise performance */
  /** Called when the element should be removed from the list. We
      overload this to allow clearing of signals.
      @param el Tree to remove from.
  */
  virtual Int_t RemoveFromListTrees(TEveElement* el);
  /** Called when the element should be removed from the list. We
      overload this to allow clearing of signals.
      @param el Parent to remove from.
  */  
  virtual void RemoveParent(TEveElement* el);
  /** @} */

  /** @{ 
      @name Singleton interface */
  /** Get the singleton instance.  If the instance has not been
      instantised yet, it will be after this call. */
  static AliEveFMDLoader* Instance();
  /** Destructor */
  virtual ~AliEveFMDLoader();
  /** @} */
protected:
  struct ModuleData 
  {
    Float_t  fScaledSum;
  };
  /** Constructor 
      @param name     Name of the folder. 
      @param useBoxes Whether to use boxes or Quads for the signals 
      @param old      Whether to enable reading old RCU data format */
  AliEveFMDLoader(const char* name="FMD", Bool_t useBoxes=true, 
		  Bool_t old=kTRUE);

  /** @{ 
      @name Service functions for loading data */
  /** Do the actual display of digits 
      @param type What to show 
      @param digits The digits */
  void          DoLoadDigits(const char* type, TClonesArray* digits);
  /** @} */

  /** @{ 
      @name Digit set management */
  /** Find a digit set corresponding to the passed parameters.  If it
      is not found, one is created 
      @param type   Type of data 
      @param d      Detector 
      @param r      Ring 
      @param s      Sector 
      @return a digit set */
  TEveDigitSet* FindDigitSet(const char* type,UShort_t d, Char_t r, UShort_t s);
  /** Make a digit set.  The type of digit set depends on the setting
      of fUseBoxDigits.  If this is true, we return a TEveBoxSet,
      otherwise a TEveQuadSet 
      @param name  Name of set. 
      @param nstr  Number of strips to make room for. 
      @return newly allocated digit set */
  virtual TEveDigitSet* MakeDigitSet(const char* name, UShort_t nstr);
  /** Clear signals of some type.   
      @param type Type of signals to clear 
      Type can be one of 
      - All    All signals 
      - Hits   Hits 
      - Digits Digits 
      - Raw    Raw 
      - ESD    ESD */
  virtual void ClearDigitSets(const char* type);
  /** @} */

  /** @{ 
      @name Adding signals to digit sets */
  /** Add a signal to a digit set
      @param type   Type of data 
      @param det    Detector 
      @param rng    Ring 
      @param sec    Sector 
      @param str    Strip
      @param signal Signal value 
      @param min    Minimum of this kind of signal 
      @param max    Maximum of this kind of signal 
      @param ref    Reference object  */
  void          AddSignal(const char* type, UShort_t det, Char_t ring, 
			  UShort_t sec, UShort_t str, 
			  Float_t signal, Float_t min, Float_t  max, 
			  TObject* ref=0);
  /** Add a signal to a digit set, with known (x,y,z) coordinates
      (this is for hits)
      @param type   Type of data 
      @param det    Detector 
      @param rng    Ring 
      @param sec    Sector 
      @param str    Strip
      @param x      X coordinate  
      @param y      Y coordinate  
      @param z      Z coordinate  
      @param signal Signal value 
      @param min    Minimum of this kind of signal 
      @param max    Maximum of this kind of signal 
      @param ref    Reference object */
  void          AddSignal(const char* type, UShort_t det, Char_t ring, 
			  UShort_t sec, UShort_t str, 
			  Double_t x, Double_t y, Double_t z, 
			  Float_t signal, Float_t min, Float_t  max, 
			  TObject* ref=0);
  /** Add a digit to a digit set. 
      @param signals Digit set. 
      @param x      X coordinate  
      @param y      Y coordinate  
      @param z      Z coordinate  
      @param w      strip pitch 
      @param scaled Scaled value 
      @param value  Signal value 
      @param ref    Reference object */
  virtual void AddDigit(TEveDigitSet* signals, 
			Double_t x, Double_t y, Double_t z, 
			Double_t w, Float_t scaled, Int_t value, 
			TObject* ref);
  /** @} */

  /** @{
      @name Various service functions */
  /** check if we shoul re-add ourselves to the current event node */
  virtual void CheckAdd();
  void SummarizeModule(TEveElement* module);
  void SummarizeModules();
  /** @} */


  /** @{ 
      @name Palettes */
  TEveRGBAPalette  fHitPalette;    // Palette for hits
  TEveRGBAPalette  fDigitPalette;  // Palette for ADC values
  TEveRGBAPalette  fMultPalette;   // Palette for multiplicity values
  /** @} */
  
  /** @{ 
      @name Settings */
  Bool_t           fUseBoxDigits;  // Whether to show as boxes or quads
  /** @} */

  /** @{ 
      @name Caches */
  TClonesArray     fHitCache;      // Cache of digits
  TClonesArray     fDigitCache;    // Cache of digits
  TClonesArray     fRawCache;      // Cache of raw
  /** @} */

  /** @{ 
      @name Singleton interface */
  static AliEveFMDLoader* fgInstance; // Singleton
  /** @} */

  ClassDef(AliEveFMDLoader,0)
};

  
  
#endif
//
// EOF
//
