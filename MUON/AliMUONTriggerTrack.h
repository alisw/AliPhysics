#ifndef ALIMUONTRIGGERTRACK_H
#define ALIMUONTRIGGERTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTriggerTrack
/// \brief Reconstructed trigger track in ALICE dimuon spectrometer
/// \author Philippe Crochet

#include <TObject.h>

class AliMUONTrackReconstructor;

class AliMUONTriggerTrack : public TObject 
{
 public:
    AliMUONTriggerTrack(); // Constructor
    virtual ~AliMUONTriggerTrack(); // Destructor
    AliMUONTriggerTrack (const AliMUONTriggerTrack& AliMUONTriggerTrack); // copy constructor
    AliMUONTriggerTrack& operator=(const AliMUONTriggerTrack& AliMUONTriggerTrack); // assignment operator
    AliMUONTriggerTrack(Float_t x11, Float_t y11, Float_t thetax, Float_t thetay,
			Int_t iloTrg, Long_t theGTPattern, UShort_t hitsPatternInTrigCh=0); 
    
    // getters
    
    /// Return x position of fired Y strip in MC11
    Float_t GetX11()    const {return fx11;}
    /// Return y position of fired X strip in MC11
    Float_t GetY11()    const {return fy11;}
    /// Return track theta angle in X 
    Float_t GetThetax() const {return fthetax;}
    /// Return track theta angle in Y
    Float_t GetThetay() const {return fthetay;}
    /// Return local trigger number
    Int_t   GetLoTrgNum() const {return floTrgNum;}    

    // setters
    
    /// Set x position of fired Y strip in MC11
    void SetX11(Float_t x)     {fx11 = x;}
    /// Set y position of fired X strip in MC11
    void SetY11(Float_t y)     {fy11 = y;}
    /// Set track theta angle in X 
    void SetThetax(Float_t tx) {fthetax = tx;}
    /// Set track theta angle in Y
    void SetThetay(Float_t ty) {fthetay = ty;}
    /// Set local trigger number
    void SetLoTrgNum(Int_t loTrgNum) {floTrgNum = loTrgNum;}    

    /// Set Global trigger pattern  (do not work with static statement) 
    void SetGTPattern(UChar_t pat) {fGTPattern = pat;}    
    /// Return Global trigger pattern  (do not work with static statement) 
    UChar_t GetGTPattern() const {return fGTPattern;}

    /// set word telling which trigger chambers where hit by track
    UShort_t GetHitsPatternInTrigCh() const {return fHitsPatternInTrigCh;}
    /// set word telling which trigger chambers where hit by track
    void     SetHitsPatternInTrigCh(UShort_t hitsPatternInTrigCh) {fHitsPatternInTrigCh = hitsPatternInTrigCh;}
    
    virtual void Print(Option_t* opt="") const;
    
protected:
  private:
  Float_t fx11;    ///< x position of fired Y strip in MC11
  Float_t fy11;    ///< y position of fired X strip in MC11
  Float_t fthetax; ///< track theta angle in X   
  Float_t fthetay; ///< track theta angle in Y
  Int_t   floTrgNum; ///< local trigger number
  UChar_t fGTPattern; ///< Global trigger pattern  (do not work with static statement)
  UShort_t fHitsPatternInTrigCh; ///< Word containing info on the hits left in trigger chambers

  ClassDef(AliMUONTriggerTrack, 5) // Reconstructed trigger track in ALICE dimuon spectrometer
    };
	
#endif

