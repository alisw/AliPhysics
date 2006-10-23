#ifndef ALIMUONV1_H
#define ALIMUONV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONv1
/// \brief AliDetector class for MUON subsystem which implements
/// functions for simulation 
 
#include "AliMUON.h"

#include <TLorentzVector.h>

class TF1;
class TGeoCombiTrans;
class TString;
class TGeoHMatrix;

class AliMUONv1 : public AliMUON 
{
 public:
   AliMUONv1();
   AliMUONv1(const char *name, const char* title = "default");
   virtual  ~AliMUONv1();
   virtual void   CreateGeometry();
   virtual void   CreateMaterials();
   virtual void   AddAlignableVolumes() const;
   virtual void   Init();
   
                  /// Return the version number of this detector
   virtual Int_t  IsVersion() const {return 1;}
   
   virtual void   StepManager();

        /// Set option to take into account angle effect
   void SetAngleEffect(Bool_t Opt) 
     { fAngleEffect = Opt; }
        
        /// Set max step max in active gas 
   void SetStepMaxInActiveGas(Float_t StepMax)
     {fStepMaxInActiveGas = StepMax; }

 protected:
   AliMUONv1(const AliMUONv1& right);
   AliMUONv1&  operator = (const AliMUONv1& right);

   virtual Int_t  GetGeomModuleId(Int_t volId) const;
   TString CurrentVolumePath() const;	     

   Bool_t  fAngleEffect;           ///< Angle Effect along wires, Default is true
   Float_t fStepMaxInActiveGas;    ///< Step max in active gas default 0.6cm

   // StepManager 
   Float_t *  fStepSum;   //!< Sum of track steps per chamber
   Float_t *  fDestepSum; //!< Sum of energy deposits per chamber
  
   TLorentzVector fTrackMomentum; ///< Momentum of the particle entering in the active gas of chamber
   TLorentzVector fTrackPosition; ///< Position of the particle exiting the active gas of chamber
   TF1 *          fElossRatio;    ///< Ratio of particle mean eloss with respect MIP's 
   TF1 *          fAngleEffect10; ///< Angle effect in tracking chambers at theta =10 degres as a function of ElossRatio (Khalil BOUDJEMLINE sep 2003 Ph.D Thesis) (in micrometers)  
   TF1 *          fAngleEffectNorma;///< Angle effect: Normalisation form theta=10 degres to theta between 0 and 10 (Khalil BOUDJEMLINE sep 2003 Ph.D Thesis)
     
   ClassDef(AliMUONv1,5)  // MUON Detector class Version 1
};
#endif







