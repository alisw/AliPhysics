#ifndef ALI_MUON_ST1_SPECIAL_MOTIF_H
#define ALI_MUON_ST1_SPECIAL_MOTIF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup sim
/// \class AliMUONSt1SpecialMotif
/// \brief Helper class to encapsulate the distance between the daughter card
/// and the pad/kapton connector
///
/// Encapsulate the distance between the center of a given daughter card
/// and the pad/kapton connector.
///
/// Authors: David Guez, Ivana Hrivnacova, Marion MacCormick; IPN Orsay

#include <TVector2.h>

class AliMUONSt1SpecialMotif  
{
  public:
    AliMUONSt1SpecialMotif();
    AliMUONSt1SpecialMotif(const TVector2& delta,Double_t rotAngle=0.);
    AliMUONSt1SpecialMotif(const AliMUONSt1SpecialMotif& src);
    virtual ~AliMUONSt1SpecialMotif();
    AliMUONSt1SpecialMotif& operator=(const AliMUONSt1SpecialMotif& rhs);
             
	     /// Return offset
    TVector2 GetDelta()    const {return fDelta;}
    
             /// Return rotation angle in degrees (0 = vertical) 
    Double_t GetRotAngle() const {return fRotAngle;}

  private:
    TVector2  fDelta;   ///< offset of this motif
    Double_t  fRotAngle;///< rotation angle in degrees (0 = vertical) 
};

#endif //ALI_MUON_ST1_SPECIAL_MOTIF_H
