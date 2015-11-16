/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2015                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
 
 
  

#include "AliQnCorrectionsQnVector.h"

ClassImp(AliQnCorrectionsQnVector)



//_______________________________________________________________________________
AliQnCorrectionsQnVector::AliQnCorrectionsQnVector() :
  //fQvecX(),
  //fQvecY(),
  fQvectorX(0x0),
  fQvectorY(0x0),
  fBin(-1),
  fQvectorNormalization(3),
  fSumW(0.0),
  fN(0),
  fEventPlaneStatus()
{
  //
  // Constructor
  //


  fQvectorX = TArrayF();
  fQvectorY = TArrayF();
  fEventPlaneStatus = TArrayC();
        

}


//_______________________________________________________________________________
AliQnCorrectionsQnVector::AliQnCorrectionsQnVector(Int_t nHarmonics1) :
  //fQvecX(),
  //fQvecY(),
  fQvectorX(0x0),
  fQvectorY(0x0),
  fBin(-1),
  fQvectorNormalization(3),
  fSumW(0.0),
  fN(0),
  fEventPlaneStatus()
{
  //
  // Constructor
  //

  fQvectorX = TArrayF(nHarmonics1);
  fQvectorY = TArrayF(nHarmonics1);
  fEventPlaneStatus = TArrayC(nHarmonics1);
        

}


//_______________________________________________________________________________
AliQnCorrectionsQnVector::~AliQnCorrectionsQnVector()
{
  //
  // De-Constructor
  //
}



//_______________________________________________________________________________
AliQnCorrectionsQnVector::AliQnCorrectionsQnVector(const AliQnCorrectionsQnVector &c) :
  TObject()
{
  //
  // Add a Q-vector
  //

  fSumW=c.SumOfWeights();
  fN=c.N();
  fBin=c.Bin();
  fQvectorNormalization=c.QvectorNormalization();
  fQvectorX=TArrayF(c.Qx());
  fQvectorY=TArrayF(c.Qy());
  fEventPlaneStatus=TArrayC(c.GetEventPlaneStatus());
  //fEventPlaneStatus.Set(c.Qx().GetSize());
  ////for(Int_t ih=1; ih<=c.Qx().GetSize(); ++ih){
  //  fQvecX[ih-1]=c.Qx(ih);
  //  fQvecY[ih-1]=c.Qy(ih);
  //  fEventPlaneStatus[ih-1]=c.GetEventPlaneStatus(ih);
  //}


}

////_______________________________________________________________________________
//void AliQnCorrectionsQnVector::SetEventPlaneStatus(Int_t harmonic, EventPlaneStatus status) { 
//  if(harmonic>0 && harmonic<=fgkEPMaxHarmonics) 
//    fEventPlaneStatus[harmonic-fgkEPMinHarmonics] |= (1<<status);
//}
//
////_______________________________________________________________________________
//void AliQnCorrectionsQnVector::UnsetEventPlaneStatus(Int_t harmonic, EventPlaneStatus status) { 
//  if(harmonic>0 && harmonic<=fgkEPMaxHarmonics) 
//    fEventPlaneStatus[harmonic-fgkEPMinHarmonics] |= (0<<status);
//}


