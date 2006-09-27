#ifndef ALIMUONTRIGGERTRACK_H
#define ALIMUONTRIGGERTRACK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*$Id$*/
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONTriggerTrack
/// \brief Reconstructed trigger track in ALICE dimuon spectrometer

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
                        Long_t theGTPattern); 

    // getter
    Float_t GetX11()    const {return fx11;}
    Float_t GetY11()    const {return fy11;}
    Float_t GetThetax() const {return fthetax;}
    Float_t GetThetay() const {return fthetay;}    

    // setter
    void SetX11(Float_t x)     {fx11 = x;}
    void SetY11(Float_t y)     {fy11 = y;}
    void SetThetax(Float_t tx) {fthetax = tx;}
    void SetThetay(Float_t ty) {fthetay = ty;}    

    void SetGTPattern(Long_t pat) {fGTPattern = pat;}    
    Long_t GetGTPattern() const {return fGTPattern;}    
    
protected:
  private:
  Float_t fx11;    ///< x position of fired Y strip in MC11
  Float_t fy11;    ///< y position of fired X strip in MC11
  Float_t fthetax; ///< track theta angle in X   
  Float_t fthetay; ///< track theta angle in Y
  Long_t fGTPattern; ///< Global trigger pattern  (do not work with static statement) 

  ClassDef(AliMUONTriggerTrack, 3) // Reconstructed trigger track in ALICE dimuon spectrometer
    };
	
#endif

