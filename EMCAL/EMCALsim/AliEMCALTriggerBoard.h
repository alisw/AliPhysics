#ifndef ALIEMCALTRIGGERBOARD_H
#define ALIEMCALTRIGGERBOARD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//________________________________________________
/// \class AliEMCALTriggerBoard
/// \ingroup EMCALsim
/// \brief EMCal trigger board super class
///
/// EMCal trigger board super class
/// run the sliding window algorithm
///
/// \author: R. GUERNANE LPSC Grenoble CNRS/IN2P3
//________________________________________________

#include "TVector2.h"

class TClonesArray;

class AliEMCALTriggerBoard : public TObject 
{	
	
public:
           AliEMCALTriggerBoard();
           AliEMCALTriggerBoard(const TVector2& RegionSize);
  virtual ~AliEMCALTriggerBoard();
  
  virtual void SlidingWindow(Int_t Threshold);
  
  virtual void ZeroRegion();
  
  virtual void Scan();
  virtual void Reset();
  
  virtual void          PatchSize(TVector2& Size) const {Size =     *fPatchSize;}
  virtual TVector2*     PatchSize(              ) const {     return fPatchSize;}
  virtual void         RegionSize(TVector2& Size) const {Size =    *fRegionSize;}
  virtual TVector2*    RegionSize(              ) const {    return fRegionSize;}
  virtual void      SubRegionSize(TVector2& Size) const {Size = *fSubRegionSize;}
  virtual TVector2* SubRegionSize(              ) const { return fSubRegionSize;}
  
  virtual const TClonesArray& Patches() const {return *fPatches;}
  
  virtual void    SetRegionSize(const TVector2& Size) {    *fRegionSize = Size;}
  virtual void     SetPatchSize(const TVector2& Size) {     *fPatchSize = Size;}
  virtual void SetSubRegionSize(const TVector2& Size) { *fSubRegionSize = Size;}
  
  virtual void SetRegion(Int_t arr[][48]) {for (Int_t i = 0; i < fRegionSize->X(); i++) for (Int_t j = 0; j < fRegionSize->Y(); j++) fRegion[i][j] = arr[i][j];}
  virtual void SetRegion(Int_t arr[][64]) {for (Int_t i = 0; i < fRegionSize->X(); i++) for (Int_t j = 0; j < fRegionSize->Y(); j++) fRegion[i][j] = arr[i][j];}
  
  virtual Int_t** Region() {return fRegion;}
  virtual Int_t**    Map() {return    fMap;}
  virtual void       Map(Int_t arr[][48], const TVector2& Size) {for (Int_t i = 0; i < Size.X(); i++) for (Int_t j = 0; j < Size.Y(); j++) arr[i][j] = fMap[i][j];}
  virtual void       Map(Int_t arr[][64], const TVector2& Size) {for (Int_t i = 0; i < Size.X(); i++) for (Int_t j = 0; j < Size.Y(); j++) arr[i][j] = fMap[i][j];}

protected:
  
  Int_t**       fRegion;        //!<! add comment 
  Int_t**       fMap;           //!<! Map the position to digit index (which eq. to ADC channel)
  TVector2*     fRegionSize;    //!<! in FastOR unit
  TVector2*     fSubRegionSize; //!<! in FastOR unit
  TVector2*     fPatchSize;     //!<! in subregion unit
  TClonesArray* fPatches;       //!<! add comment
  
private:
  
  AliEMCALTriggerBoard           (const AliEMCALTriggerBoard& rhs); // NOT implemented
  AliEMCALTriggerBoard& operator=(const AliEMCALTriggerBoard& rhs); // NOT implemented
  
  /// \cond CLASSIMP
  ClassDef(AliEMCALTriggerBoard,1) ;
  /// \endcond

};
 
#endif //ALIEMCALTRIGGERBOARD_H
