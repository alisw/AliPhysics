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
#include <TMatrixD.h>
#include <TMath.h>

class AliMUONTrackReconstructor;

class AliMUONTriggerTrack : public TObject 
{
 public:
    AliMUONTriggerTrack(); // Constructor
    virtual ~AliMUONTriggerTrack(); // Destructor
    AliMUONTriggerTrack (const AliMUONTriggerTrack& AliMUONTriggerTrack); // copy constructor
    AliMUONTriggerTrack& operator=(const AliMUONTriggerTrack& AliMUONTriggerTrack); // assignment operator
    AliMUONTriggerTrack(Float_t x11, Float_t y11, Float_t z11, Float_t z21,
			Float_t slopeX, Float_t slopeY,
			Int_t iloTrg, Long_t theGTPattern, Int_t ptCutLevel=1); 
    
    virtual void Clear(Option_t* opt = "");
    
    // getters
    
    /// Return x position of fired Y strip in MC11
    Float_t GetX11()    const {return fx11;}
    /// Return y position of fired X strip in MC11
    Float_t GetY11()    const {return fy11;}
    /// Return z position of fired X strip in MC11
    Float_t GetZ11()    const {return fz11;}
    /// Return z position of fired X strip in MC21
    Float_t GetZ21()    const {return fz21;}
    /// Return track theta angle in X 
    Float_t GetThetax() const {return TMath::ATan(fSlopeX);}
    /// Return track theta angle in Y
    Float_t GetThetay() const {return TMath::ATan(fSlopeY);}
    /// Return track slope in X 
    Float_t GetSlopeX() const {return fSlopeX;}
    /// Return track slope in Y
    Float_t GetSlopeY() const {return fSlopeY;}
    /// Return local trigger number
    Int_t   GetLoTrgNum() const {return floTrgNum;}    

    // setters
    
    /// Set x position of fired Y strip in MC11
    void SetX11(Float_t x)     {fx11 = x;}
    /// Set y position of fired X strip in MC11
    void SetY11(Float_t y)     {fy11 = y;}
    /// Set z position of fired X strip in MC11
    void SetZ11(Float_t z)     {fz11 = z;}
    /// Set z position of fired X strip in MC21
    void SetZ21(Float_t z)     {fz21 = z;}
    /// Set track slope in X 
    void SetSlopeX(Float_t slopeX) {fSlopeX = slopeX;}
    /// Set track slope in Y
    void SetSlopeY(Float_t slopeY) {fSlopeY = slopeY;}
    /// Set local trigger number
    void SetLoTrgNum(Int_t loTrgNum) {floTrgNum = loTrgNum;}    

    /// Set Global trigger pattern  (do not work with static statement) 
    void SetGTPattern(UChar_t pat) {fGTPattern = pat;}    
    /// Return Global trigger pattern  (do not work with static statement) 
    UChar_t GetGTPattern() const {return fGTPattern;}

    /// Set word telling which trigger chambers where hit by track
    void     SetHitsPatternInTrigCh(UShort_t hitsPatternInTrigCh) {fHitsPatternInTrigCh = hitsPatternInTrigCh;}
    /// Get word telling which trigger chambers where hit by track
    UShort_t GetHitsPatternInTrigCh() const {return fHitsPatternInTrigCh;}
  
    /// Set pt cut level
    void  SetPtCutLevel(Int_t ptCutLevel) { fPtCutLevel = ptCutLevel;}
    /// Get pt cut level
    Int_t GetPtCutLevel() const {return fPtCutLevel;}

    
    virtual void Print(Option_t* opt="FULL") const;

    const TMatrixD& GetCovariances() const;
    void  SetCovariances(const TMatrixD& covariances);
    void  SetCovariances(const Double_t matrix[3][3]);
  
    Bool_t Match(AliMUONTriggerTrack &track, Double_t sigmaCut) const;
    
protected:
  private:
  Float_t fx11;    ///< x position of fired Y strip in MC11
  Float_t fy11;    ///< y position of fired X strip in MC11
  Float_t fz11;    ///< z position of fired X strip in MC11
  Float_t fz21;    ///< z position of fired X strip in MC21
  Float_t fSlopeX; ///< track slope in X   
  Float_t fSlopeY; ///< track slope in Y
  Int_t   floTrgNum; ///< local trigger number
  UChar_t fGTPattern; ///< Global trigger pattern  (do not work with static statement)
  Int_t   fPtCutLevel;  ///< Trigger pt cut level as in ESDs (1->Apt; 2->Lpt; 3->Hpt)
  UShort_t fHitsPatternInTrigCh; ///< Word containing info on the hits left in trigger chambers
  mutable TMatrixD *fCovariances; ///< Covariance matrix of track parameters 

  ClassDef(AliMUONTriggerTrack, 7) // Reconstructed trigger track in ALICE dimuon spectrometer
    };
	
#endif

