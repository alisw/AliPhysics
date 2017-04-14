#ifndef ALIEMCALSDIGITIZER_H
#define ALIEMCALSDIGITIZER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
///  
/// \class AliEMCALSDigitizer
/// \ingroup EMCALbase
/// \brief EMCal summable digits maker  
///
/// This is a class that makes SDigits out of Hits
/// A Summable Digits is the sum of all hits originating 
/// from one in one tower of the EMCAL 
/// A threshold for assignment of the primary to SDigit is applied 
///
/// SDigits need to hold the energy sum of the hits, but AliEMCALDigit
/// can (should) only store amplitude.  Therefore, the SDigit energy is
/// "digitized" before being stored and must be "calibrated" back to an
/// energy before SDigits are summed to form true Digits
///
/// SDigits are written to TreeS, branch "EMCAL"
/// AliEMCALSDigitizer with all current parameters is written 
/// to TreeS branch "AliEMCALSDigitizer".
/// Both branches have the same title. If necessary one can produce 
/// another set of SDigits with different parameters. Two versions
/// can be distunguished using titles of the branches.
/// User case:
/// root [0] AliEMCALSDigitizer * s = new AliEMCALSDigitizer("galice.root")
/// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
/// root [1] s->Digitize()
///             // Makes SDigitis for all events stored in galice.root
/// root [2] s->SetPedestalParameter(0.001)
///             // One can change parameters of digitization
/// root [3] s->SetSDigitsBranch("Redestal 0.001")
///             // and write them into the new branch
/// root [4] s->Digitize("deb all tim")
///             // available parameters:
///             deb - print # of produced SDigitis
///             deb all  - print # and list of produced SDigits
///             tim - print benchmarking information
///
/// First versions based on : AliPHOSDigit
///
/// \author Sahal Yacoob (LBL)
/// \author Jenn Klay (LBNL)
///
////////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TNamed.h"
class TBrowser;

// --- AliRoot header files ---
#include "AliConfig.h"

class AliEMCALSDigitizer: public TNamed 
{

public:

  AliEMCALSDigitizer() ;          // ctor
  AliEMCALSDigitizer(const char * alirunFileName, const char * eventFolderName = AliConfig::GetDefaultEventFolderName()) ; 
  AliEMCALSDigitizer(const AliEMCALSDigitizer & sd) ;
  
  AliEMCALSDigitizer& operator = (const AliEMCALSDigitizer& source) ;
  
  Bool_t operator == (const AliEMCALSDigitizer & sd) const ;

  virtual ~AliEMCALSDigitizer(); // dtor

  Float_t       Digitize(Float_t energy) const; // convert energy in GeV to int amplitude
  
  Float_t       Calibrate(Float_t amp)   const; // opposite of Digitize()

  virtual void  Digitize(Option_t *option=""); 
  
  Int_t         GetSDigitsInRun() const { return fSDigitsInRun ; }  
  
  virtual void  Print(Option_t *option="") const;
  void          Print1(Option_t *option="all");  // *MENU*
  
  void          SetEventFolderName(TString name)            { fEventFolderName = name ; }
  void          SetEventRange(Int_t first=0, Int_t last=-1) { fFirstEvent=first ; fLastEvent=last ; }

  virtual void  Browse(TBrowser* b);

private:
  
  void     Init() ;
  
  void     InitParameters() ; 
  
  void     PrintSDigits(Option_t * option) ;
  
  void     Unload() const ;

private:
  
  Float_t fA ;                     ///<  Pedestal parameter
  Float_t fB ;                     ///<  Slope Digitizition parameters
  
  Float_t fECPrimThreshold ;       ///<  To store primary if EC Shower Elos > threshold
  
  Bool_t  fDefaultInit;            //!<! Says if the object was created by defaut ctor (only parameters are initialized)
  
  TString fEventFolderName;        ///<  Event folder name
  
  Bool_t  fInit ;                  //!<! Tells if initialisation went OK, will revent exec if not
  
  Int_t   fSDigitsInRun ;          //!<! Total number of sdigits in one run
  
  Int_t   fFirstEvent;             ///<  First event to process
  Int_t   fLastEvent;              ///<  Last  event to process
  
  Float_t fSampling;               ///< See AliEMCALGeometry
  
  /// Temporal array with hits
  TClonesArray* fHits;             //-> 
	
  /// \cond CLASSIMP
  ClassDef(AliEMCALSDigitizer,8) ;
  /// \endcond

};

#endif // AliEMCALSDIGITIZER_H

