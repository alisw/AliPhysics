// Author: Filimon Roukoutakis 
/******************************************************************************
  MOOD - Monitor Of On-line Data and Detector Debugger for ALICE Experiment
******************************************************************************/

#ifndef TMCAL_H
#define TMCAL_H

/******************************************************************************
 *                                                                            *
 * TMCAL                                                                      *
 *                                                                            *
 * CAL Module Class; base class for EMCAL and PHOS monitors                   *
 *                                                                            *
 * Author: Timo Alho (Jyvaskyla), original version                            *
 * [Consultant: D. Silvermyr (ORNL)]                                          *
 *                                                                            *
 ******************************************************************************/

#include "TMBaseModule.h"
#include "Mood.h"
#include "MoodOC.hh"
#include "TGSlider.h"

//This class will handle the pedestal data analysis
#include "AliRawReader.h"
#include "AliCaloRawStream.h"
#include "AliCaloCalibPedestal.h"

//TMCalConstants.hh includes all the defaults and GUI layout etc..
#include "TMCalConstants.hh" 

class TMCal: public TMBaseModule {
  
  RQ_OBJECT("TMCal")

 protected:

  virtual void ConstructGUI(void) {};
  virtual void ConstructMonitors(void) {};
  virtual void DestructMonitors(void) {};
  virtual void DestructGUI(void) {};
  virtual void InitMonitors(void) {};
  virtual void MonitorEvent(void);
  virtual void UpdateMonitors(void);
  virtual void ResetMonitors(void) {};
  virtual void PreMonitor(void) {};
  virtual void PostMonitor(void);
 
 private:
 
  AliCaloCalibPedestal::kDetType fDetType; //The detector type for this object
  int fNumModules;
 
  enum kTab {kPeds, kPeaks, kPedsDiff, kPeaksDiff, kPedsRatio, kPeaksRatio, 
	     kPedsSingle, kPeaksSingle, kPedsDiffSingle, kPeaksDiffSingle, kPedsRatioSingle, kPeaksRatioSingle};
 
  int	fEvents; //The number of events processed
  AliCaloCalibPedestal *fPedestals; //The object taking care of pedestal and pedestal-peak profile generation
  AliRawReader *fRawReader; //The raw data reader
  AliCaloRawStream *fCaloRawStream; //The object to sort this out to EMCAL/PHOS data
 
  /////////////////////////////
  //GUI components
 
  //Histocanvases. The last one, that is [fgkNumModules] is the single histo display canvas.
  //This isn't too neat or good programming style, but it does make handling the updates etc simpler and more unified -> less prone to errors
  TCanvas * fPedsCanvas[fgkMaxModules+1]; //!The pointers to the Ped canvases for the modules
  TCanvas * fPeaksCanvas[fgkMaxModules+1]; //!The pointers to the Peak canvases for the modules
  TCanvas * fPedsDiffCanvas[fgkMaxModules+1]; //!The pointers to the Ped difference canvases for the modules
  TCanvas * fPeaksDiffCanvas[fgkMaxModules+1]; //!The pointers to the Peak difference canvases for the modules
  TCanvas * fPedsRatioCanvas[fgkMaxModules+1]; //!The pointers to the Ped ratio canvases for the modules
  TCanvas * fPeaksRatioCanvas[fgkMaxModules+1]; //!The pointers to the Peak ratio anvases for the modules
  TCanvas * fDeadMapCanvas[fgkMaxModules+1];
 
  //Single histo canvases. These hold pointers to the single histo canvases. Note that a copy is in the main canvas arrays above
  TCanvas * fPedsCanvasSingle; //!
  TCanvas * fPeaksCanvasSingle; //!
  TCanvas * fPedsDiffCanvasSingle; //!
  TCanvas * fPeaksDiffCanvasSingle; //!
  TCanvas * fPedsRatioCanvasSingle; //!
  TCanvas * fPeaksRatioCanvasSingle; //!
  TCanvas * fDeadMapCanvasSingle; //!
  
  TGLabel * fDeadCountText;
  
  //The tab pointers
  TGCompositeFrame * fPedsTab;
  TGCompositeFrame * fPeaksTab;
  TGCompositeFrame * fPedsDiffTab;
  TGCompositeFrame * fPeaksDiffTab;
  TGCompositeFrame * fPedsRatioTab;
  TGCompositeFrame * fPeaksRatioTab;
  TGCompositeFrame * fDeadMapTab;
  
  //The vertical major frames for the multihisto displays in each of the tabs
  TGCompositeFrame * fPeaksDiffFrame;
  TGCompositeFrame * fPedsDiffFrame;
  TGCompositeFrame * fPeaksFrame;
  TGCompositeFrame * fPedsFrame;
  TGCompositeFrame * fPeaksRatioFrame;
  TGCompositeFrame * fPedsRatioFrame;
  TGCompositeFrame * fDeadMapFrame;
  
  //The same for the single histo displays
  TGCompositeFrame * fPeaksDiffFrameSingle;
  TGCompositeFrame * fPedsDiffFrameSingle;
  TGCompositeFrame * fPeaksFrameSingle;
  TGCompositeFrame * fPedsFrameSingle;
  TGCompositeFrame * fPeaksRatioFrameSingle;
  TGCompositeFrame * fPedsRatioFrameSingle;
  TGCompositeFrame * fDeadMapFrameSingle;
  
  //The buttons on the lower panel
  TGTextButton * fGainButton;
  TGTextButton * fSingleModeButton;
  TGNumberEntry * fRunNoEntry;
  
  //The range sliders
  TGVSlider * fPeaksRange;
  TGVSlider * fPedsRange;
  TGVSlider * fPeaksDiffRange;
  TGVSlider * fPedsDiffRange;
  TGVSlider * fPeaksRatioRange;
  TGVSlider * fPedsRatioRange;
  //The range sliders for single module
  TGVSlider * fPeaksRangeSingle;
  TGVSlider * fPedsRangeSingle;
  TGVSlider * fPeaksDiffRangeSingle;
  TGVSlider * fPedsDiffRangeSingle;
  TGVSlider * fPeaksRatioRangeSingle;
  TGVSlider * fPedsRatioRangeSingle;
  
  //GUI State
  Bool_t fLowGainMode; //Low or high gain to be displayed
  Bool_t fSingleModule; //If we are in a single module state
  int fVisibleModule; //The module which is visible, if we're in the single module mode
  Int_t fDeadMapPalette[AliCaloCalibPedestal::kNumDeadMapStates]; //The palette used for drawing the dead map. See TMCalConstants.hh for the color
  //definition (which are assigned in the constructor of this class)
  //The easiest way to define a string constant...
#define fgkDeadMapTabName "DeadMap"
  TPaveText * fDeadMapPaveText;

  double fMaxPed;
  double fMaxPeak;
  double fMaxPeakDiff;
  double fMaxPedDiff;
  double fMaxPedRatio;
  double fMaxPeakRatio;
  
  //////////////////////////////
  //GUI Creation helpers
  
  void DrawHistos();//Draw the histos				     
  TGCompositeFrame * MakeAllModulesTab(TGCompositeFrame * tab,
				       TCanvas **canvArray, //Array of pointers to the canvases to which the histos have been drawn. MUST be initialized.
				       TGVSlider **slider = 0, //The pointer to the pointer to the slider, which to put the slider to
				       kTab tabId = kPeds, //The tab Id which the slider should give as a parameter when connecting
				       double sliderMin = 0.0, double sliderMax = 0.0, double sliderPos = 0.0 //The range and position of the slider
				       );//The function to create a tab full of histos

  TGCompositeFrame * MakeSingleModuleTab(TGCompositeFrame * tab,
					 TCanvas **canvas,
					 TGVSlider **slider = 0, //The pointer to the pointer to the slider, which to put the slider to
					 kTab tabId = kPeds, //The tab Id which the slider should give as a parameter when connecting
					 double sliderMin = 0.0, double sliderMax = 0.0, double sliderPos = 0.0 //The range and position of the slider
					 );
  void ConstructFullGUI(); //Constructs the GUI with all modules visible
  void ConstructMenu(UInt_t w); //Constructs the GUI with all modules visible
  void InitRawReaderAndStream(unsigned int runNo); //Inits the Raw  reader and CaloRawStream
  void SetTitles();//Sets the title of the histos
  
 public:
 
  TMCal(const TGWindow *p, UInt_t w, UInt_t h, AliCaloCalibPedestal::kDetType detectorType);
  virtual ~TMCal();
 
  //GUI event handlers
  void OnGainButton();
  void OnModeButton(int button);
  void OnRange(kTab tab);
  void OnRunNoChange(Long_t val);
  void OnTabChange(Int_t Id);//Here we have to set the appropriate palette for the tab

  ClassDef(TMCal, 1); // Calorimeter Module Base Class
  
};

#endif
