#ifndef ALIPHOS_H
#define ALIPHOS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.70  2007/01/17 17:28:56  kharlov
 * Extract ALTRO sample generation to a separate class AliPHOSPulseGenerator
 *
 * Revision 1.69  2006/11/14 17:11:15  hristov
 * Removing inheritances from TAttLine, TAttMarker and AliRndm in AliModule. The copy constructor and assignment operators are moved to the private part of the class and not implemented. The corresponding changes are propagated to the detectors
 *
 * Revision 1.68  2006/08/11 12:36:25  cvetan
 * Update of the PHOS code needed in order to read and reconstruct the beam test raw data (i.e. without an existing galice.root)
 *
 * Revision 1.67  2006/04/07 08:42:00  hristov
 * Follow AliAlignObj framework and remove AliPHOSAlignData (Yu.Kharlov)
 *
 * Revision 1.66  2006/03/24 21:39:33  schutz
 * Modification needed to include PHOS in the global trigger framework
 *
 * Revision 1.65  2006/03/07 18:56:25  kharlov
 * CDB is passed via environment variable
 *
 * Revision 1.64  2005/11/03 13:09:19  hristov
 * Removing meaningless const declarations (linuxicc)
 *
 * Revision 1.63  2005/07/26 13:32:39  kharlov
 * Restoring raw data fit from version of 29-Aug-2004
 *
 * Revision 1.62  2005/07/06 10:10:32  hristov
 * Moving the functions used to initialize TF1 and TF2 to the pivate part of the class
 *
 * Revision 1.61  2005/05/28 12:10:07  schutz
 * Copy constructor is corrected (by T.P.)
 *
 */


//_________________________________________________________________________
//  Base Class for PHOS     
//                  
//*-- Author: Laurent Aphecetche & Yves Schutz (SUBATECH)


// --- ROOT system ---
class TString ; 
class TTask ;
class TFolder ;
class TTree ; 
class TRandom ; 

// --- AliRoot header files ---
#include "AliDetector.h" 
#include "AliLog.h"
#include "AliPHOSGeometry.h" 
#include "AliPHOSTrigger.h"

class AliPHOS : public AliDetector {

public:

  AliPHOS() ;
  AliPHOS(const char* name, const char* title="") ;  
  virtual ~AliPHOS() ; 
  virtual void   AddHit(Int_t, Int_t*, Float_t *) {
    // do not use this definition but the one below
    AliFatal(Form("do not use")) ;
    
  }
  virtual void   AddHit( Int_t shunt, Int_t primary, Int_t track, 
			 Int_t id, Float_t *hits ) = 0 ;   
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* manager) const;
  virtual void  CreateMaterials() ;            
  virtual void  Digits2Raw();
  virtual void  FinishRun() {;}
  virtual AliPHOSGeometry * GetGeometry() const 
  {return AliPHOSGeometry::GetInstance(GetTitle(),"") ;  }

  virtual void    Hits2SDigits();
  virtual Int_t   IsVersion(void) const = 0 ;  
  virtual void    Init();
  virtual AliTriggerDetector* CreateTriggerDetector() const 
    { return new AliPHOSTrigger(); }

  virtual AliLoader* MakeLoader(const char* topfoldername);
  virtual void    SetTreeAddress();   
  virtual const TString Version() const {return TString(" ") ; } 

 private:                                        
  AliPHOS(AliPHOS & phos);
  AliPHOS & operator = (const AliPHOS & /*rvalue*/);

  ClassDef(AliPHOS,6) // Photon Spectrometer Detector (base class)
} ;

#endif // ALIPHOS_H
