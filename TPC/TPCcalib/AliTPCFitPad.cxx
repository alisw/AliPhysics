/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////////////////////////////////
//                                                                        
//       === Class for fitting properties specific to pad regions ===
//
//    For each pad size region a separate TLinearFitter object is assigned.
//    Commonly used functions such as getting the center of a pad size
//    region or visualization functions are provided. Also, choosing the
//    segment and pad type is made easy due to different methods that
//    can calculate these informations from other coordinates.
//
////////////////////////////////////////////////////////////////////////////

#include "AliTPCFitPad.h"

ClassImp(AliTPCFitPad)




AliTPCFitPad::AliTPCFitPad(Int_t ndim, const char* formula, Option_t* opt) :
   AliTPCCalPadRegion("", ""),
   fNdim(ndim),
   fFormula(formula),
   fOpt(opt)
{
   //
   // Constructor. The parameters are used for generating new TLinearFitter
   // objects and are described in its documentation.
   //
}

AliTPCFitPad& AliTPCFitPad::operator=(const AliTPCFitPad& rhs)
{
  //
  // Assignment operator.
  //
  
  if (this != &rhs) {
    AliTPCCalPadRegion::operator=(rhs);
    fNdim = rhs.fNdim;
    fFormula = rhs.fFormula;
    fOpt = rhs.fOpt;
  }
  return *this;
}

AliTPCFitPad::AliTPCFitPad(const AliTPCFitPad& rhs):
  AliTPCCalPadRegion(rhs),
  fNdim(rhs.fNdim),
  fFormula(rhs.fFormula),
  fOpt(rhs.fOpt)

{
  //
  // Copy constructor
  //
}

AliTPCFitPad::~AliTPCFitPad() {
   //
   // Destructor.
   //

   Delete();
}



void AliTPCFitPad::Add(AliTPCFitPad* fit) {
   //
   // Adds another AliTPCFitPad object to this object. The formula should be the
   // same, though it won't be checked!
   //

   for (UInt_t iSegment = 0; iSegment < GetNSegments(); iSegment++) {
      for (UInt_t iPadType = 0; iPadType < GetNPadTypes(); iPadType++) {
         TLinearFitter* fitter = fit->GetFitterSimple(iSegment, iPadType);
         // parameter workaround == kTRUE because it is not possible to add another
         // TLinearFitter object to a "virgin" one. Thus a dummy data point is added
         // and cleared again immediately afterwards. Due to a broken TLinearFitter
         // copy constructor this is a necessary workaround.
         if (fitter) {
            cout << "TLinearFitter::Add called for " << iSegment << ", " << iPadType << endl;
            GetFitter(iSegment, iPadType, kTRUE)->Add(fitter);
         }
      }
   }
}

TLinearFitter* AliTPCFitPad::GetFitterSimple(UInt_t segment, UInt_t padType) {
   //
   // This method returns the fitter corresponding to segment and pad type.
   // In contrast to GetFitter() no fitter will be created, if it does
   // not exist, but a null pointer is returned.
   //
   
   return (TLinearFitter*)(GetObject(segment, padType));
}

TLinearFitter* AliTPCFitPad::GetFitter(UInt_t segment, UInt_t padType, Bool_t workaround) {
   //
   // This method returns the fitter corresponding
   // to segment and pad type.
   // If the fitter doesn't exist yet, it will be created on the fly
   // according to the parameters passed to the constructor.
   //
   // The workaround parameter should always be kFALSE. (It is only used by the Add method.)
   //

   TLinearFitter* fitter = GetFitterSimple(segment, padType);
   if (fitter == 0 || fitter->GetNumberTotalParameters() !=fNdim) {
     fitter = new TLinearFitter(fNdim, fFormula, fOpt);
     fitter->StoreData(kFALSE);
     SetObject(fitter, segment, padType);
     fitter = (TLinearFitter*)(GetObject(segment, padType));
     if (workaround) {
       Double_t x[1000];
       for (Int_t i = 0; i < fNdim; i++) x[i] = 3.141592;
       fitter->AddPoint(x, 31.41592);
       fitter->ClearPoints();
       //cout << "workaround called for " << segment << ", " << padType << endl;
     }
   }
   return fitter;
}

Int_t AliTPCFitPad::Evaluate(Bool_t robust, Double_t frac) {
   //
   // Evaluates all fitters. Returns 0 if successful, 1 in case of an error.
   // If the robust option is set to kTRUE a robust fit is performed with frac as
   // the minimal fraction of good points (see TLinearFitter::EvalRobust for details).
   // Beware: Robust fitting is much slower!
   //

   Int_t returnCode = 0;
   for (UInt_t iSegment = 0; iSegment < GetNSegments(); iSegment++) {
      for (UInt_t iPadType = 0; iPadType < GetNPadTypes(); iPadType++) {
         if (TLinearFitter* fitter = GetFitterSimple(iSegment, iPadType)) {
            Int_t status = robust ? fitter->EvalRobust(frac) : fitter->Eval();
            if (status != 0) {
               returnCode = 1;
               Error("Evaluate", "Error in evaluation of fitter in segment %d, pad region %d", iSegment, iPadType);
            }
         }
      }
   }
   return returnCode;
}
