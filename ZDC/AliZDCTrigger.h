/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALIZDCTRIGGER_H
#define ALIZDCTRIGGER_H

/// \ingroup sim
/// \class AliZDCTrigger
/// \brief ZDC trigger class
///
/////////////////////////////////////////////////
///  ZDC Trigger Detector Class               //
/////////////////////////////////////////////////

#include "AliTriggerDetector.h"

class AliZDCTrigger : public AliTriggerDetector
{
 public:
   AliZDCTrigger();  // constructor
   virtual ~AliZDCTrigger(){}  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

   // Print method
/*   virtual void Print(Option_t *) const {
     printf("\t AliZDCTrigger: fZNMinCut = %1.0f, fZDCMinCut = %1.0f, fZEMMinCut= %1.0f \n"
     "fZDCLeftEMDCuts = [%1.0f, %1.0f], fZDCRightEMDCuts = [%1.0f, %1.0f], fZDCMBCut = %1.0f\n"
     "fZDCCentrCut = %1.0f, fZDCSemiCentrCut = %1.0f, fZEMCentrCut = %1.0f\n\n",
     fZNMinCut,fZDCMinCut,fZEMMinCut,fZDCLeftEMDCuts[0],fZDCLeftEMDCuts[1],
     fZDCRightEMDCuts[0],fZDCRightEMDCuts[1],fZDCMBCut,fZDCCentrCut,fZDCSemiCentrCut,
     fZEMCentrCut);
   }
*/
 
 protected:
   
   // Setters   
   void SetZNMinCut(Float_t ZNMinCut);
   void SetZDCLeftMinCut(Float_t ZDCLeftMinCut);
   void SetZDCRightMinCut(Float_t ZDCRightMinCut);
   void SetZEMMinCut(Float_t ZEMMinCut);
   void SetZDCLeftEMDCuts(Float_t *ZDCLeftEMDCuts);
   void SetZDCLeftEMDCuts(Float_t ZDCLeftEMDCutInf, Float_t ZDCLeftEMDCutSup);
   void SetZDCRightEMDCuts(Float_t *ZDCRightEMDCuts);
   void SetZDCRightEMDCuts(Float_t ZDCRightEMDCutInf, Float_t  ZDCRightEMDCutSup);
   void SetZDCLeftMBCut(Float_t ZDCLeftMBCut);
   void SetZDCRightMBCut(Float_t ZDCRightMBCut);
   void SetZDCLeftCentrCut(Float_t ZDCLeftCentrCuts);
   void SetZDCRightCentrCut(Float_t ZDCRightCentrCuts);
   void SetZDCLeftSemiCentrCut(Float_t ZDCLeftSemiCentrCut);
   void SetZDCRightSemiCentrCut(Float_t ZDCRightSemiCentrCut);
   void SetZEMCentrCut(Float_t ZEMCentrCut);

   // Data member
   Float_t fZNMinCut;   
   Float_t fZDCLeftMinCut;   
   Float_t fZDCRightMinCut;   
   Float_t fZEMMinCut;   
   Float_t fZDCLeftEMDCuts[2];
   Float_t fZDCRightEMDCuts[2];
   Float_t fZDCLeftMBCut;
   Float_t fZDCRightMBCut;
   Float_t fZDCLeftCentrCut;
   Float_t fZDCRightCentrCut;
   Float_t fZDCLeftSemiCentrCut;
   Float_t fZDCRightSemiCentrCut;
   Float_t fZEMCentrCut;   
    
  ClassDef(AliZDCTrigger,1)  // ZDC Trigger Detector class
};
#endif








