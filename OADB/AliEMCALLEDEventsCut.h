#ifndef ALIEMCALLEDEventsCut_H
#define ALIEMCALLEDEventsCut_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliEMCALLEDEventsCut
/// \brief Class for event EMCAL LED events finding
///
///  EMCal LED events sometimes fire in accepted as physics events.
///  In those events LED system lights part or the full EMCal SMs
///  Here events are tagged as bad in such case, depending EMCal cell multiplicity or energy per SM 
///  Different prescription for LHC11a period and Run2 periods. 
///
/// More information on the problem can be found in [this presentation](https://indico.cern.ch/event/889785/contributions/3755598/attachments/1989660/3316682/Run2pp13TeV_EMCALDataFeaturesAndFixes_PWGGAMeeting.pdf).
/// Copy from AliPhysics/PWG/CaloTrackCorrBase/AliCaloTrackReader::RejectLEDEvents()
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS
//_________________________________________________________________________


class AliVEvent;
class TH2F;

#include "TObject.h" 
#include "TList.h" 

#include "AliEMCALGeometry.h"


class AliEMCALLEDEventsCut : public TObject {
  
public:
  
  AliEMCALLEDEventsCut() ; // ctor
  
  virtual         ~AliEMCALLEDEventsCut() {} // virtual dtor
  
  //
  // Main methods
  //
  
  void             InitControlHistograms();
  
  TList *          GetControlHistograms()  { return fHistogramContainer ; }
  
  Bool_t           IsEMCALLEDEvent(AliVEvent * event, Int_t runNumber) ; 

  void             PrintParam() const ;

  //
  // Settings
  //
  void             FillControlHistograms(Bool_t fill)      { fFillHisto             = fill   ; }
  
  void             SetLEDHighEnergyCutSM(Float_t e)        { fLEDHighEnergyCutSM    = e      ; }   
  
  void             SetLEDHighNCellsCutSM(Int_t   n)        { fLEDHighNCellsCutSM    = n      ; }   
  
  void             SetLEDLowEnergyCutSM3(Float_t e)        { fLEDLowEnergyCutSM3    = e      ; }   
  
  void             SetLEDLowNCellsCutSM3(Int_t   n)        { fLEDLowNCellsCutSM3    = n      ; }   
  
  void             SetLEDMinCellEnergy  (Float_t e)        { fLEDMinCellEnergy      = e      ; }   
  
  void             SetLEDMaxCellEnergy  (Float_t e)        { fLEDMaxCellEnergy      = e      ; }  
  
  void             SwitchOnLEDStripEventsRemoval()         { fRemoveLEDStripEvents  = kTRUE  ; }
  
  void             SwitchOffLEDStripEventsRemoval()        { fRemoveLEDStripEvents  = kFALSE ; }
  
 
  void             SetLEDStripHighEnergyCutSM(Float_t eFull, Float_t eThird)   
  { fLEDHighEnergyCutStrip[0] = eFull ; fLEDHighEnergyCutStrip[1] = eThird ; }   
  
  void             SetLEDStripHighNCellsCutSM(Int_t   nFull, Int_t   nThird) 
  { fLEDHighNCellsCutStrip[0] = nFull ; fLEDHighNCellsCutStrip[0] = nThird ; }   
  
  void             SetLEDStripLowEnergyCutSM3(Float_t e)   { fLEDLowEnergyCutSM3Strip    = e ; }   
  
  void             SetLEDStripLowNCellsCutSM3(Int_t   n)   { fLEDLowNCellsCutSM3Strip    = n ; }   
  
  void             SetLEDEventMaxNumberOfStrips(Int_t n)   { fLEDEventMaxNumberOfStrips  = n ; }   
  
private:
  
  AliEMCALGeometry * fEMCalGeom;                   //!< EMCal geometry pointer
  
  TList *          fHistogramContainer;            ///< output control histograms container
  
  Bool_t           fFillHisto;                     ///< Fill control histograms, set in AliEventCuts
  
  Float_t          fLEDHighEnergyCutSM;            ///<  SM is too active if energy above this value, likely LED event 
  
  Int_t            fLEDHighNCellsCutSM;            ///<  SM is too active if n cells above this value, likely LED event 
  
  Float_t          fLEDLowEnergyCutSM3;            ///<  SM3 low activity if energy below this value, check activity on other SM for LED event (Run2)
  
  Int_t            fLEDLowNCellsCutSM3;            ///<  SM3 low activity if n cells below this value, check activity on other SM LED event (Run2)
  
  Float_t          fLEDMinCellEnergy;              ///<  Count or sum cells energy above this value to determine if event had LEDs
  
  Float_t          fLEDMaxCellEnergy;              ///<  Count or sum cells energy below this value to determine if event had LEDs
  
  Int_t            fRemoveLEDStripEvents;          ///<  Remove events where an LED strip or more was wrongly firing - only EMCAL 
  
  Int_t            fLEDEventMaxNumberOfStrips;     ///<  Cut on events with a number of too active strips
  
  Float_t          fLEDHighEnergyCutStrip[2];      ///<  SM strip is too active if energy above this value, likely LED event. [0] Full SM, [1] 1/3 SM 
  
  Int_t            fLEDHighNCellsCutStrip[2];      ///<  SM strip is too active if n cells above this value, likely LED event. [0] Full SM, [1] 1/3 SM  
  
  Float_t          fLEDLowEnergyCutSM3Strip;       ///<  SM3 strip low activity if energy below this value, check activity on other SM for LED event (Run2)
  
  Int_t            fLEDLowNCellsCutSM3Strip;       ///<  SM3 strip low activity if n cells below this value, check activity on other SM LED event (Run2)
    
  //
  // QA check histograms
  //
  TH2F  *          fhEMCALNSumEnCellsPerSM;              //!<! Control histogram of LED events rejection
  
  TH2F  *          fhEMCALNSumEnCellsPerSMAfter;         //!<! Control histogram of LED events rejection, after cut
  
  TH2F  *          fhEMCALNSumEnCellsPerSMAfterStripCut; //!<! Control histogram of LED events rejection, after LED strip rejection
  
  TH2F  *          fhEMCALNSumEnCellsPerStrip;           //!<! Control histogram of LED events on strips rejection, after LED SM rejection
  
  TH2F  *          fhEMCALNSumEnCellsPerStripAfter;      //!<! Control histogram of LED events on strips rejection, after strip LED and SM rejection
  
  /// Copy constructor not implemented.
  AliEMCALLEDEventsCut(              const AliEMCALLEDEventsCut & r) ; 
  
  /// Assignment operator not implemented.
  AliEMCALLEDEventsCut & operator = (const AliEMCALLEDEventsCut & r) ; 
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALLEDEventsCut,1) ;
  /// \endcond 
  
} ;

#endif //ALIEMCALLEDEVENTSCUT_H



