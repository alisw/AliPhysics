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
  fQvectorX(0x0),
  fQvectorY(0x0),
  fQnVectorStatus(),
  fBin(-1),
  fQvectorNormalization(3),
  fSumW(0.0),
  fN(0)
{
  //
  // Constructor
  //


  fQvectorX = TArrayF();
  fQvectorY = TArrayF();
  fQnVectorStatus = TArrayC();
        

}


//_______________________________________________________________________________
AliQnCorrectionsQnVector::AliQnCorrectionsQnVector(Int_t nHarmonics1) :
  fQvectorX(0x0),
  fQvectorY(0x0),
  fQnVectorStatus(),
  fBin(-1),
  fQvectorNormalization(3),
  fSumW(0.0),
  fN(0)
{
  //
  // Constructor
  //

  fQvectorX = TArrayF(nHarmonics1);
  fQvectorY = TArrayF(nHarmonics1);
  fQnVectorStatus = TArrayC(nHarmonics1);
        

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
  fQnVectorStatus=TArrayC(c.GetQnVectorStatus());

}


