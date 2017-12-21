/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliEmcalJetShapeProperties.h"

/**
 * Default constructor
 */
AliEmcalJetShapeProperties::AliEmcalJetShapeProperties():
  fJetShapeMassFirstDer(0),
  fJetShapeMassSecondDer(0),
  fJetShapeMassFirstSub(0),
  fJetShapeMassSecondSub(0),
  fGRNumerator(0),
  fGRDenominator(0),
  fGRNumeratorSub(0),
  fGRDenominatorSub(0),
  fJetShapeAngularityFirstDer(0),
  fJetShapeAngularitySecondDer(0),
  fJetShapeAngularityFirstSub(0),
  fJetShapeAngularitySecondSub(0),
  fJetShapepTDFirstDer(0),
  fJetShapepTDSecondDer(0),
  fJetShapepTDFirstSub(0),
  fJetShapepTDSecondSub(0),
  fJetShapeCircularityFirstDer(0),
  fJetShapeCircularitySecondDer(0),
  fJetShapeCircularityFirstSub(0),
  fJetShapeCircularitySecondSub(0),
  fJetShapeSigma2FirstDer(0),
  fJetShapeSigma2SecondDer(0),
  fJetShapeSigma2FirstSub(0),
  fJetShapeSigma2SecondSub(0),
  fJetShapeConstituentFirstDer(0),
  fJetShapeConstituentSecondDer(0),
  fJetShapeConstituentFirstSub(0),
  fJetShapeConstituentSecondSub(0),
  fJetShapeLeSubFirstDer(0),
  fJetShapeLeSubSecondDer(0),
  fJetShapeLeSubFirstSub(0),
  fJetShapeLeSubSecondSub(0),
  fJetShape1subjettinessktFirstDer(0),
  fJetShape1subjettinessktSecondDer(0),
  fJetShape1subjettinessktFirstSub(0),
  fJetShape1subjettinessktSecondSub(0),
  fJetShape2subjettinessktFirstDer(0),
  fJetShape2subjettinessktSecondDer(0),
  fJetShape2subjettinessktFirstSub(0),
  fJetShape2subjettinessktSecondSub(0),
  fJetShape3subjettinessktFirstDer(0),
  fJetShape3subjettinessktSecondDer(0),
  fJetShape3subjettinessktFirstSub(0),
  fJetShape3subjettinessktSecondSub(0),
  fJetShapeOpeningAnglektFirstDer(0),
  fJetShapeOpeningAnglektSecondDer(0),
  fJetShapeOpeningAnglektFirstSub(0),
  fJetShapeOpeningAnglektSecondSub(0),
  fJetShape1subjettinesscaFirstDer(0),
  fJetShape1subjettinesscaSecondDer(0),
  fJetShape1subjettinesscaFirstSub(0),
  fJetShape1subjettinesscaSecondSub(0),
  fJetShape2subjettinesscaFirstDer(0),
  fJetShape2subjettinesscaSecondDer(0),
  fJetShape2subjettinesscaFirstSub(0),
  fJetShape2subjettinesscaSecondSub(0),
  fJetShapeOpeningAnglecaFirstDer(0),
  fJetShapeOpeningAnglecaSecondDer(0),
  fJetShapeOpeningAnglecaFirstSub(0),
  fJetShapeOpeningAnglecaSecondSub(0),
  fJetShape1subjettinessakt02FirstDer(0),
  fJetShape1subjettinessakt02SecondDer(0),
  fJetShape1subjettinessakt02FirstSub(0),
  fJetShape1subjettinessakt02SecondSub(0),
  fJetShape2subjettinessakt02FirstDer(0),
  fJetShape2subjettinessakt02SecondDer(0),
  fJetShape2subjettinessakt02FirstSub(0),
  fJetShape2subjettinessakt02SecondSub(0),
  fJetShapeOpeningAngleakt02FirstDer(0),
  fJetShapeOpeningAngleakt02SecondDer(0),
  fJetShapeOpeningAngleakt02FirstSub(0),
  fJetShapeOpeningAngleakt02SecondSub(0),
  fJetShape1subjettinesscasdFirstDer(0),
  fJetShape1subjettinesscasdSecondDer(0),
  fJetShape1subjettinesscasdFirstSub(0),
  fJetShape1subjettinesscasdSecondSub(0),
  fJetShape2subjettinesscasdFirstDer(0),
  fJetShape2subjettinesscasdSecondDer(0),
  fJetShape2subjettinesscasdFirstSub(0),
  fJetShape2subjettinesscasdSecondSub(0),
  fJetShapeOpeningAnglecasdFirstDer(0),
  fJetShapeOpeningAnglecasdSecondDer(0),
  fJetShapeOpeningAnglecasdFirstSub(0),
  fJetShapeOpeningAnglecasdSecondSub(0),
  fSoftDropZg(0),
  fSoftDropdR(0),
  fSoftDropPtfrac(0),
  fSoftDropDropCount(0)
{
}

/**
 * Copy constructor
 *
 * @param[in] jet Const reference to copy the content from
 */
AliEmcalJetShapeProperties::AliEmcalJetShapeProperties(const AliEmcalJetShapeProperties &jet):
  fJetShapeMassFirstDer(jet.fJetShapeMassFirstDer),
  fJetShapeMassSecondDer(jet.fJetShapeMassSecondDer),
  fJetShapeMassFirstSub(jet.fJetShapeMassFirstSub),
  fJetShapeMassSecondSub(jet.fJetShapeMassSecondSub),
  fGRNumerator(jet.fGRNumerator),
  fGRDenominator(jet.fGRDenominator),
  fGRNumeratorSub(jet.fGRNumeratorSub),
  fGRDenominatorSub(jet.fGRDenominatorSub),
  fJetShapeAngularityFirstDer(jet.fJetShapeAngularityFirstDer),
  fJetShapeAngularitySecondDer(jet.fJetShapeAngularitySecondDer),
  fJetShapeAngularityFirstSub(jet.fJetShapeAngularityFirstSub),
  fJetShapeAngularitySecondSub(jet.fJetShapeAngularitySecondSub),
  fJetShapepTDFirstDer(jet.fJetShapepTDFirstDer),
  fJetShapepTDSecondDer(jet.fJetShapepTDSecondDer),
  fJetShapepTDFirstSub(jet.fJetShapepTDFirstSub),
  fJetShapepTDSecondSub(jet.fJetShapepTDSecondSub),
  fJetShapeCircularityFirstDer(jet.fJetShapeCircularityFirstDer),
  fJetShapeCircularitySecondDer(jet.fJetShapeCircularitySecondDer),
  fJetShapeCircularityFirstSub(jet.fJetShapeCircularityFirstSub),
  fJetShapeCircularitySecondSub(jet.fJetShapeCircularitySecondSub),
  fJetShapeSigma2FirstDer(jet.fJetShapeSigma2FirstDer),
  fJetShapeSigma2SecondDer(jet.fJetShapeSigma2SecondDer),
  fJetShapeSigma2FirstSub(jet.fJetShapeSigma2FirstSub),
  fJetShapeSigma2SecondSub(jet.fJetShapeSigma2SecondSub),
  fJetShapeConstituentFirstDer(jet.fJetShapeConstituentFirstDer),
  fJetShapeConstituentSecondDer(jet.fJetShapeConstituentSecondDer),
  fJetShapeConstituentFirstSub(jet.fJetShapeConstituentFirstSub),
  fJetShapeConstituentSecondSub(jet.fJetShapeConstituentSecondSub),
  fJetShapeLeSubFirstDer(jet.fJetShapeLeSubFirstDer),
  fJetShapeLeSubSecondDer(jet.fJetShapeLeSubSecondDer),
  fJetShapeLeSubFirstSub(jet.fJetShapeLeSubFirstSub),
  fJetShapeLeSubSecondSub(jet.fJetShapeLeSubSecondSub),
  fJetShape1subjettinessktFirstDer(jet.fJetShape1subjettinessktFirstDer),
  fJetShape1subjettinessktSecondDer(jet.fJetShape1subjettinessktSecondDer),
  fJetShape1subjettinessktFirstSub(jet.fJetShape1subjettinessktFirstSub),
  fJetShape1subjettinessktSecondSub(jet.fJetShape1subjettinessktSecondSub),
  fJetShape2subjettinessktFirstDer(jet.fJetShape2subjettinessktFirstDer),
  fJetShape2subjettinessktSecondDer(jet.fJetShape2subjettinessktSecondDer),
  fJetShape2subjettinessktFirstSub(jet.fJetShape2subjettinessktFirstSub),
  fJetShape2subjettinessktSecondSub(jet.fJetShape2subjettinessktSecondSub),
  fJetShape3subjettinessktFirstDer(jet.fJetShape3subjettinessktFirstDer),
  fJetShape3subjettinessktSecondDer(jet.fJetShape3subjettinessktSecondDer),
  fJetShape3subjettinessktFirstSub(jet.fJetShape3subjettinessktFirstSub),
  fJetShape3subjettinessktSecondSub(jet.fJetShape3subjettinessktSecondSub),
  fJetShapeOpeningAnglektFirstDer(jet.fJetShapeOpeningAnglektFirstDer),
  fJetShapeOpeningAnglektSecondDer(jet.fJetShapeOpeningAnglektSecondDer),
  fJetShapeOpeningAnglektFirstSub(jet.fJetShapeOpeningAnglektFirstSub),
  fJetShapeOpeningAnglektSecondSub(jet.fJetShapeOpeningAnglektSecondSub),
  fJetShape1subjettinesscaFirstDer(jet.fJetShape1subjettinesscaFirstDer),
  fJetShape1subjettinesscaSecondDer(jet.fJetShape1subjettinesscaSecondDer),
  fJetShape1subjettinesscaFirstSub(jet.fJetShape1subjettinesscaFirstSub),
  fJetShape1subjettinesscaSecondSub(jet.fJetShape1subjettinesscaSecondSub),
  fJetShape2subjettinesscaFirstDer(jet.fJetShape2subjettinesscaFirstDer),
  fJetShape2subjettinesscaSecondDer(jet.fJetShape2subjettinesscaSecondDer),
  fJetShape2subjettinesscaFirstSub(jet.fJetShape2subjettinesscaFirstSub),
  fJetShape2subjettinesscaSecondSub(jet.fJetShape2subjettinesscaSecondSub),
  fJetShapeOpeningAnglecaFirstDer(jet.fJetShapeOpeningAnglecaFirstDer),
  fJetShapeOpeningAnglecaSecondDer(jet.fJetShapeOpeningAnglecaSecondDer),
  fJetShapeOpeningAnglecaFirstSub(jet.fJetShapeOpeningAnglecaFirstSub),
  fJetShapeOpeningAnglecaSecondSub(jet.fJetShapeOpeningAnglecaSecondSub),
  fJetShape1subjettinessakt02FirstDer(jet.fJetShape1subjettinessakt02FirstDer),
  fJetShape1subjettinessakt02SecondDer(jet.fJetShape1subjettinessakt02SecondDer),
  fJetShape1subjettinessakt02FirstSub(jet.fJetShape1subjettinessakt02FirstSub),
  fJetShape1subjettinessakt02SecondSub(jet.fJetShape1subjettinessakt02SecondSub),
  fJetShape2subjettinessakt02FirstDer(jet.fJetShape2subjettinessakt02FirstDer),
  fJetShape2subjettinessakt02SecondDer(jet.fJetShape2subjettinessakt02SecondDer),
  fJetShape2subjettinessakt02FirstSub(jet.fJetShape2subjettinessakt02FirstSub),
  fJetShape2subjettinessakt02SecondSub(jet.fJetShape2subjettinessakt02SecondSub),
  fJetShapeOpeningAngleakt02FirstDer(jet.fJetShapeOpeningAngleakt02FirstDer),
  fJetShapeOpeningAngleakt02SecondDer(jet.fJetShapeOpeningAngleakt02SecondDer),
  fJetShapeOpeningAngleakt02FirstSub(jet.fJetShapeOpeningAngleakt02FirstSub),
  fJetShapeOpeningAngleakt02SecondSub(jet.fJetShapeOpeningAngleakt02SecondSub),
  fJetShape1subjettinesscasdFirstDer(jet.fJetShape1subjettinesscasdFirstDer),
  fJetShape1subjettinesscasdSecondDer(jet.fJetShape1subjettinesscasdSecondDer),
  fJetShape1subjettinesscasdFirstSub(jet.fJetShape1subjettinesscasdFirstSub),
  fJetShape1subjettinesscasdSecondSub(jet.fJetShape1subjettinesscasdSecondSub),
  fJetShape2subjettinesscasdFirstDer(jet.fJetShape2subjettinesscasdFirstDer),
  fJetShape2subjettinesscasdSecondDer(jet.fJetShape2subjettinesscasdSecondDer),
  fJetShape2subjettinesscasdFirstSub(jet.fJetShape2subjettinesscasdFirstSub),
  fJetShape2subjettinesscasdSecondSub(jet.fJetShape2subjettinesscasdSecondSub),
  fJetShapeOpeningAnglecasdFirstDer(jet.fJetShapeOpeningAnglecasdFirstDer),
  fJetShapeOpeningAnglecasdSecondDer(jet.fJetShapeOpeningAnglecasdSecondDer),
  fJetShapeOpeningAnglecasdFirstSub(jet.fJetShapeOpeningAnglecasdFirstSub),
  fJetShapeOpeningAnglecasdSecondSub(jet.fJetShapeOpeningAnglecasdSecondSub),
  fSoftDropZg(jet.fSoftDropZg),
  fSoftDropdR(jet.fSoftDropdR),
  fSoftDropPtfrac(jet.fSoftDropPtfrac),
  fSoftDropDropCount(jet.fSoftDropDropCount)

{
}

/**
 * Assignment operator
 *
 * @param[in] jet Const reference to copy the content from
 */
AliEmcalJetShapeProperties& AliEmcalJetShapeProperties::operator=(const AliEmcalJetShapeProperties &jet)
{
  fJetShapeMassFirstDer  = jet.fJetShapeMassFirstDer;
  fJetShapeMassSecondDer = jet.fJetShapeMassSecondDer;
  fJetShapeMassFirstSub  = jet.fJetShapeMassFirstSub;
  fJetShapeMassSecondSub = jet.fJetShapeMassSecondSub;
  fGRNumerator        = jet.fGRNumerator;
  fGRDenominator      = jet.fGRDenominator;
  fGRNumeratorSub     = jet.fGRNumeratorSub;
  fGRDenominatorSub   = jet.fGRDenominatorSub;
  fJetShapeAngularityFirstDer  = jet.fJetShapeAngularityFirstDer;
  fJetShapeAngularitySecondDer = jet.fJetShapeAngularitySecondDer;
  fJetShapeAngularityFirstSub  = jet.fJetShapeAngularityFirstSub;
  fJetShapeAngularitySecondSub = jet.fJetShapeAngularitySecondSub;
  fJetShapepTDFirstDer  = jet.fJetShapepTDFirstDer;
  fJetShapepTDSecondDer = jet.fJetShapepTDSecondDer;
  fJetShapepTDFirstSub  = jet.fJetShapepTDFirstSub;
  fJetShapepTDSecondSub = jet.fJetShapepTDSecondSub;
  fJetShapeCircularityFirstDer  = jet.fJetShapeCircularityFirstDer;
  fJetShapeCircularitySecondDer = jet.fJetShapeCircularitySecondDer;
  fJetShapeCircularityFirstSub  = jet.fJetShapeCircularityFirstSub;
  fJetShapeCircularitySecondSub = jet.fJetShapeCircularitySecondSub;
  fJetShapeSigma2FirstDer  = jet.fJetShapeSigma2FirstDer;
  fJetShapeSigma2SecondDer = jet.fJetShapeSigma2SecondDer;
  fJetShapeSigma2FirstSub  = jet.fJetShapeSigma2FirstSub;
  fJetShapeSigma2SecondSub = jet.fJetShapeSigma2SecondSub;
  fJetShapeConstituentFirstDer  = jet.fJetShapeConstituentFirstDer;
  fJetShapeConstituentSecondDer = jet.fJetShapeConstituentSecondDer;
  fJetShapeConstituentFirstSub  = jet.fJetShapeConstituentFirstSub;
  fJetShapeConstituentSecondSub = jet.fJetShapeConstituentSecondSub;
  fJetShapeLeSubFirstDer  = jet.fJetShapeLeSubFirstDer;
  fJetShapeLeSubSecondDer = jet.fJetShapeLeSubSecondDer;
  fJetShapeLeSubFirstSub  = jet.fJetShapeLeSubFirstSub;
  fJetShapeLeSubSecondSub = jet.fJetShapeLeSubSecondSub;
  fJetShape1subjettinessktFirstDer  = jet.fJetShape1subjettinessktFirstDer;
  fJetShape1subjettinessktSecondDer = jet.fJetShape1subjettinessktSecondDer;
  fJetShape1subjettinessktFirstSub  = jet.fJetShape1subjettinessktFirstSub;
  fJetShape1subjettinessktSecondSub = jet.fJetShape1subjettinessktSecondSub;
  fJetShape2subjettinessktFirstDer  = jet.fJetShape2subjettinessktFirstDer;
  fJetShape2subjettinessktSecondDer = jet.fJetShape2subjettinessktSecondDer;
  fJetShape2subjettinessktFirstSub  = jet.fJetShape2subjettinessktFirstSub;
  fJetShape2subjettinessktSecondSub = jet.fJetShape2subjettinessktSecondSub;
  fJetShape3subjettinessktFirstDer  = jet.fJetShape3subjettinessktFirstDer;
  fJetShape3subjettinessktSecondDer = jet.fJetShape3subjettinessktSecondDer;
  fJetShape3subjettinessktFirstSub  = jet.fJetShape3subjettinessktFirstSub;
  fJetShape3subjettinessktSecondSub = jet.fJetShape3subjettinessktSecondSub;
  fJetShapeOpeningAnglektFirstDer  = jet.fJetShapeOpeningAnglektFirstDer;
  fJetShapeOpeningAnglektSecondDer = jet.fJetShapeOpeningAnglektSecondDer;
  fJetShapeOpeningAnglektFirstSub  = jet.fJetShapeOpeningAnglektFirstSub;
  fJetShapeOpeningAnglektSecondSub = jet.fJetShapeOpeningAnglektSecondSub;
  fJetShape1subjettinesscaFirstDer  = jet.fJetShape1subjettinesscaFirstDer;
  fJetShape1subjettinesscaSecondDer = jet.fJetShape1subjettinesscaSecondDer;
  fJetShape1subjettinesscaFirstSub  = jet.fJetShape1subjettinesscaFirstSub;
  fJetShape1subjettinesscaSecondSub = jet.fJetShape1subjettinesscaSecondSub;
  fJetShape2subjettinesscaFirstDer  = jet.fJetShape2subjettinesscaFirstDer;
  fJetShape2subjettinesscaSecondDer = jet.fJetShape2subjettinesscaSecondDer;
  fJetShape2subjettinesscaFirstSub  = jet.fJetShape2subjettinesscaFirstSub;
  fJetShape2subjettinesscaSecondSub = jet.fJetShape2subjettinesscaSecondSub;
  fJetShapeOpeningAnglecaFirstDer  = jet.fJetShapeOpeningAnglecaFirstDer;
  fJetShapeOpeningAnglecaSecondDer = jet.fJetShapeOpeningAnglecaSecondDer;
  fJetShapeOpeningAnglecaFirstSub  = jet.fJetShapeOpeningAnglecaFirstSub;
  fJetShapeOpeningAnglecaSecondSub = jet.fJetShapeOpeningAnglecaSecondSub;
  fJetShape1subjettinessakt02FirstDer  = jet.fJetShape1subjettinessakt02FirstDer;
  fJetShape1subjettinessakt02SecondDer = jet.fJetShape1subjettinessakt02SecondDer;
  fJetShape1subjettinessakt02FirstSub  = jet.fJetShape1subjettinessakt02FirstSub;
  fJetShape1subjettinessakt02SecondSub = jet.fJetShape1subjettinessakt02SecondSub;
  fJetShape2subjettinessakt02FirstDer  = jet.fJetShape2subjettinessakt02FirstDer;
  fJetShape2subjettinessakt02SecondDer = jet.fJetShape2subjettinessakt02SecondDer;
  fJetShape2subjettinessakt02FirstSub  = jet.fJetShape2subjettinessakt02FirstSub;
  fJetShape2subjettinessakt02SecondSub = jet.fJetShape2subjettinessakt02SecondSub;
  fJetShapeOpeningAngleakt02FirstDer  = jet.fJetShapeOpeningAngleakt02FirstDer;
  fJetShapeOpeningAngleakt02SecondDer = jet.fJetShapeOpeningAngleakt02SecondDer;
  fJetShapeOpeningAngleakt02FirstSub  = jet.fJetShapeOpeningAngleakt02FirstSub;
  fJetShapeOpeningAngleakt02SecondSub = jet.fJetShapeOpeningAngleakt02SecondSub;
  fJetShape1subjettinesscasdFirstDer  = jet.fJetShape1subjettinesscasdFirstDer;
  fJetShape1subjettinesscasdSecondDer = jet.fJetShape1subjettinesscasdSecondDer;
  fJetShape1subjettinesscasdFirstSub  = jet.fJetShape1subjettinesscasdFirstSub;
  fJetShape1subjettinesscasdSecondSub = jet.fJetShape1subjettinesscasdSecondSub;
  fJetShape2subjettinesscasdFirstDer  = jet.fJetShape2subjettinesscasdFirstDer;
  fJetShape2subjettinesscasdSecondDer = jet.fJetShape2subjettinesscasdSecondDer;
  fJetShape2subjettinesscasdFirstSub  = jet.fJetShape2subjettinesscasdFirstSub;
  fJetShape2subjettinesscasdSecondSub = jet.fJetShape2subjettinesscasdSecondSub;
  fJetShapeOpeningAnglecasdFirstDer  = jet.fJetShapeOpeningAnglecasdFirstDer;
  fJetShapeOpeningAnglecasdSecondDer = jet.fJetShapeOpeningAnglecasdSecondDer;
  fJetShapeOpeningAnglecasdFirstSub  = jet.fJetShapeOpeningAnglecasdFirstSub;
  fJetShapeOpeningAnglecasdSecondSub = jet.fJetShapeOpeningAnglecasdSecondSub;
  //
  fSoftDropZg = jet.fSoftDropZg;
  fSoftDropdR = jet.fSoftDropdR;
  fSoftDropPtfrac = jet.fSoftDropPtfrac;
  fSoftDropDropCount = jet.fSoftDropDropCount;

  return *this;
}

/**
 * Print the list of the GR properties in the standard output
 */
void AliEmcalJetShapeProperties::PrintGR() const
{
  for(Int_t i = 0; i < fGRNumerator.GetSize(); i++) {
    Printf("num[%d] = %f", i, fGRNumerator.At(i));
  }
}
