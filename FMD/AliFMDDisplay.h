#ifndef AliFMDDISPLAY_H
#define AliFMDDISPLAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDDisplay.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:09 2006
    @brief   FMD Event display 
*/
//___________________________________________________________________
//
// The classes defined here, are utility classes for reading in data
// for the FMD.  They are  put in a seperate library to not polute the
// normal libraries.  The classes are intended to be used as base
// classes for customized class that do some sort of analysis on the
// various types of data produced by the FMD. 
//
#include "AliFMDInput.h"
#include <TObjArray.h>
#include <TTimer.h>
class TCanvas;
class TPad;
class TButton;
class TSlider;
class TH1;

//___________________________________________________________________
/** @class AliFMDDisplay 
    @brief Utility class to visualize FMD data in geometry. 
    @ingroup FMD_util
 */
class AliFMDDisplay : public AliFMDInput
{
public:
  /** 
   * Constructor
   *
   * @param onlyFMD Only show the FMD
   * @param gAliceFile galice file
   */
  AliFMDDisplay(Bool_t onlyFMD=kTRUE, const char* gAliceFile="galice.root");
  /** 
   * DTOR 
   */
  virtual ~AliFMDDisplay();
  /** 
   * Singleton access function
   *
   * @return Singleton object. 
   */
  static AliFMDDisplay* Instance();

  /** 
   * Continue to next event 
   */
  void  Continue() { fWait = kFALSE; }
  /** 
   * Run throug events as fast as possible 
   */ 
  void Start() { fContinous = kTRUE; fWait = kFALSE; }
  /** 
   * Pause the processing 
   */ 
  void Pause() { fContinous = kFALSE; fWait = kTRUE; }
  /** 
   * Zoom mode 
   */
  void  Zoom() { fZoomMode = kTRUE; }
  /** 
   * Pick mode 
   */
  void  Pick() { fZoomMode = kFALSE; }
  /** 
   * Redisplay the event 
   */ 
  virtual void Redisplay(); // *MENU*
  /** 
   * Break 
   */
  virtual void Break();
  /** 
   * Render in 3D 
   */
  virtual void Render();
  
  /** 
   * Change cut 
   */
  virtual void ChangeCut();
  /** 
   * Change cut 
   */
  virtual void ChangeFactor();
  /** 
   * Called when a mouse or similar event happens in the display. 
   *
   * @param event Event type
   * @param px    where the event happened in pixels along X
   * @param py    where the event happened in pixels along Y 
   */
  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
  /** 
   * Paint into canvas 
   *
   * @param option Not used 
   */
  virtual void  Paint(Option_t* option="") { (void)option; }

  /** 
   * Initialize
   *
   * @return  @c false on error 
   */
  virtual Bool_t Init();
  /** 
   * Called at beginning of an event 
   * 
   * @param event Event number
   * @return @c false on error  
   */
  virtual Bool_t Begin(Int_t event);
  /** 
   * Called at end of an event 
   *
   * @return @c false on error  
   */
  virtual Bool_t End();
  /** 
   * Visualize a hit
   *
   * @param hit Hit
   * @param p   Track
   * @return @c false on error  
   */
  virtual Bool_t ProcessHit(AliFMDHit* hit, TParticle* p);
  /** 
   * Visualize a digit
   *
   * @param digit Digit to draw
   * @return @c false on error  
   */
  virtual Bool_t ProcessDigit(AliFMDDigit* digit);
  /** 
   * Visualize a summable digit
   *
   * @param sdigit Summable digit to draw
   * @return @c false on error  
   */
  virtual Bool_t ProcessSDigit(AliFMDSDigit* sdigit);
  /** 
   * Visualize a raw digit
   *
   * @param digit Raw digit.
   * @return @c false on error  
   */
  virtual Bool_t ProcessRawDigit(AliFMDDigit* digit);
  /** 
   * Visualize a raw digit
   * 
   * @param digit Raw digit.
   * @return @c false on error  
   */
  virtual Bool_t ProcessRawCalibDigit(AliFMDDigit* digit);
  /** 
   * Visualize a reconstructed point.
   *
   * @param recpoint Reconstructed point
   * @return @c false on error  
   */
  virtual Bool_t ProcessRecPoint(AliFMDRecPoint* recpoint);
  /** 
   * Process ESD data for the FMD.  Users should overload this to
   * deal with ESD data. 
   *
   * @param d    Detector number (1-3)
   * @param r    Ring identifier ('I' or 'O')
   * @param s    Sector number (0-19, or 0-39)
   * @param t    Strip number (0-511, or 0-255)
   * @param eta  Psuedo-rapidity 
   * @param mult Psuedo-multiplicity 
   * 
   * @return  @c false on error  
   */
  virtual Bool_t ProcessESD(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
			    Float_t eta, Float_t mult);
  /** 
   * Look up a color index, based on the value @a x and the maximum
   * value of @a x
   *
   * @param x   Value 
   * @param max Maximum (for example 1023 for digits)
   * @return @c false on error  
   */
  virtual Int_t  LookupColor(Float_t x, Float_t min, Float_t max)  const;

  /** 
   * Set range of displayed values 
   */
  virtual void SetCut(Float_t l=0., Float_t h=1.); //*MENU*
  /** 
   * Set the noise factor
   * 
   * @param f  Noise factor 
   */
  virtual void SetFactor(Float_t f=1);
protected:
  /** 
   * Copy constructor 
   *
   * @param o Object to copy from  
   */
  AliFMDDisplay(const AliFMDDisplay& o) 
    : AliFMDInput(o),
      fWait(kFALSE),
      fMarkers(0),
      fHits(0),
      fCanvas(0),
      fPad(0),
      fButtons(0),
      fSlider(0),
      fFactor(0),
      fZoomMode(0),
      fX0(0),
      fY0(0),
      fX1(0),
      fY1(0),
      fXPixel(0),
      fYPixel(0),
      fOldXPixel(0),
      fOldYPixel(0),
      fLineDrawn(0), 
      fOnlyFMD(kTRUE),
      fSpec(0), 
      fSpecCut(0),
      fAux(0),
      fReturn(kFALSE),
      fContinous(kFALSE), 
      fTimeout(), 
      fInitialMin(0), 
      fInitialMax(1), 
      fInitialFactor(3./10)
  { } 
  /** 
   * Assignment operator 
   *
   * @return Reference to this object 
   */
  AliFMDDisplay& operator=(const AliFMDDisplay&) { return *this; } 
  /** 
   * Add a marker to the display
   *
   * @param x   X position
   * @param y   Y position
   * @param z   Z position
   * @param o   Object to refer to
   * @param s   Signal 
   * @param max Maximum of signal 
   */
  virtual void AddMarker(Float_t x, Float_t y, Float_t z, 
			 TObject* o, Float_t s, Float_t min, Float_t max);
  /** 
   * Add a marker to the display
   *
   * @param det Detector
   * @param rng Ring
   * @param sec Sector 
   * @param str Strip
   * @param o   Object to refer to
   * @param s   Signal 
   * @param max Maximum of signal 
   */
  virtual void AddMarker(UShort_t det, Char_t rng, UShort_t sec, UShort_t str, 
			 TObject* o, Float_t s, Float_t min, Float_t max);
  
  /** 
   * Show only the FMD detectors. 
   */
  void ShowOnlyFMD();
  /** 
   * Make base canvas 
   */ 
  virtual void MakeCanvas(const char** which);
  /** 
   * Make AUX canvas 
   * 
   */
  virtual void MakeAux();
  /** 
   * Draw AUX canvas 
   * 
   */
  virtual void DrawAux();
  /** 
   * Set the ROOT event loop to idle 
   * 
   */
  virtual void Idle();
  /** 
   * Called at end of event loop
   * 
   */
  virtual void AtEnd();
  /** 
   * Whether a point is inside 
   * 
   * @param v   Point
   * @param min Minimum
   * @param max Maximum
   * 
   * @return true if @a v is inside cut 
   */
  virtual Bool_t InsideCut(Float_t v, const Float_t& min, 
			 const Float_t& max) const;
  /** 
   * Get the ADC threshold
   * 
   * @param d Detector
   * @param r Ring 
   * @param s Sector 
   * @param t Strip 
   * 
   * @return The threshold 
   */
  virtual Double_t GetADCThreshold(UShort_t d, Char_t r, 
				   UShort_t s, UShort_t t) const;
  
  static AliFMDDisplay* fgInstance; // Static instance 
  Bool_t                fWait;      // Wait until user presses `Continue'
  TObjArray*            fMarkers;   // Cache of markers
  TObjArray*            fHits;      // Cache of `hits'
  TCanvas*              fCanvas;    // Canvas to draw in 
  TPad*                 fPad;       // View pad. 
  TObjArray             fButtons;   // Continue button
  TSlider*              fSlider;    // Cut slider
  TSlider*              fFactor;    // Factor slider
  Bool_t                fZoomMode;  // Whether we're in Zoom mode
  Float_t               fX0;        // X at lower left corner or range 
  Float_t               fY0;        // Y at lower left corner or range 
  Float_t               fX1;        // X at upper right corner or range 
  Float_t               fY1;        // Y at upper right corner or range 
  Int_t                 fXPixel;    // X pixel of mark
  Int_t                 fYPixel;    // Y pixel of mark
  Int_t                 fOldXPixel; // Old x pixel of mark
  Int_t                 fOldYPixel; // Old y pixel of mark
  Bool_t                fLineDrawn; // Whether we're drawing a box
  Bool_t                fOnlyFMD;   // Whether to only do FMD
  TH1*                  fSpec;      // Spectra
  TH1*                  fSpecCut;   // Cut spectra
  TCanvas*              fAux;       // Aux canvas.
  Bool_t                fReturn;    // Stop 
  Bool_t                fContinous; // Run continuous - do not stop
  TTimer                fTimeout;   // Timeout 
  Float_t               fInitialMin;// Initial minimum
  Float_t               fInitialMax;// Initial maximum
  Float_t               fInitialFactor;//Initial factor 

  struct Range_t { 
    UInt_t  fNbins;  // Number of bins
    Float_t fLow;    // Low cut
    Float_t fHigh;   // High cut
  };
  static const Range_t fgkEdepRange; // Energy loss range
  static const Range_t fgkAdcRange;  // ADC counts range 
  static const Range_t fgkMultRange; // Multiplicity range 

  ClassDef(AliFMDDisplay,0)  // FMD specialised event display
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
