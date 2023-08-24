////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoModelGausLCMSFreezeOutGeneratorGR - freeze-out                     ///
/// coordinates generator, generating a 3D gaussian ellipsoid in LCMS        ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifndef ALIFEMTOMODELGAUSLCMSFREEZEOUTGENERATORGR_H
#define ALIFEMTOMODELGAUSLCMSFREEZEOUTGENERATORGR_H

#include "AliFemtoModelFreezeOutGenerator.h"

#include "TRandom.h"

class AliFemtoModelGausLCMSFreezeOutGeneratorGR : public AliFemtoModelFreezeOutGenerator
{
public:
AliFemtoModelGausLCMSFreezeOutGeneratorGR();
AliFemtoModelGausLCMSFreezeOutGeneratorGR(const AliFemtoModelGausLCMSFreezeOutGeneratorGR &aModel);
virtual ~AliFemtoModelGausLCMSFreezeOutGeneratorGR();
AliFemtoModelGausLCMSFreezeOutGeneratorGR& operator=(const AliFemtoModelGausLCMSFreezeOutGeneratorGR &aModel);
virtual void GenerateFreezeOut(AliFemtoPair *aPair);

void SetSizeOut(Double_t aSizeOut);
void SetSizeSide(Double_t aSizeSide);
void SetSizeLong(Double_t aSizeLong);

void SetSelectPrimaryFromHidden(bool aUse);
Bool_t GetSelectPrimaryFromHidden();

Double_t GetSizeOut() const;
Double_t GetSizeSide() const;
Double_t GetSizeLong() const;

virtual AliFemtoModelFreezeOutGenerator* Clone() const;

protected:
Double_t fSizeOut;  // Size of the source in the out direction
Double_t fSizeSide; // Size of the source in the side direction
Double_t fSizeLong; // Size of the source in the long direction
Bool_t fSelectPrimary;    // If set to true, the existing hidden info is assumed
                          // to contain the particle creation point (in cm)
                          // and the model will try to guess whether the particle
                          // is primary based on that and assign creation point
                          // only for primary particles

private:
AliFemtoModelFreezeOutGenerator* GetGenerator() const;

#ifdef __ROOT__
/// \cond CLASSIMP
ClassDef(AliFemtoModelGausLCMSFreezeOutGeneratorGR, 1);
/// \endcond
#endif

};

#endif


