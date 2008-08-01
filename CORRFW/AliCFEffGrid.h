#ifndef ALICFEFFGRID_H
#define ALICFEFFGRID_H

/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFEffGrid Class                                              //
// Class to handle efficiency grids                                   // 
//                                                                    //
//--------------------------------------------------------------------//

#include "AliCFGridSparse.h"
#include "AliCFContainer.h"
class TH1D;
class TH2D;
class TH3D;

class AliCFEffGrid : public AliCFGridSparse
{
 public:
  AliCFEffGrid();
  AliCFEffGrid(const Char_t* name,const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn, const Double_t  *binLimitsIn=0);
  AliCFEffGrid(const Char_t* name,const Char_t* title,const AliCFContainer &c);
  AliCFEffGrid(const AliCFEffGrid& eff);
  
  virtual ~AliCFEffGrid();
  AliCFEffGrid& operator=(const AliCFEffGrid& eff);
  virtual Double_t GetAverage() const ;
  virtual Double_t GetAverage(Double_t *varMin,Double_t *varMax) const ;
  virtual Int_t GetSelNumStep() const {return fSelNum;};
  virtual Int_t GetSelDenStep() const {return fSelDen;};
  virtual TH1D* Project( Int_t ivar) const;
  virtual TH2D* Project( Int_t ivar1, Int_t ivar2) const;
  virtual TH3D* Project( Int_t ivar1, Int_t ivar2,Int_t ivar3) const;

  //Efficiency calculation
  virtual void  CalculateEfficiency(Int_t istep1, Int_t istep2);
  virtual const AliCFVGrid*  GetNum() {return (AliCFVGrid*)fContainer->GetGrid(fSelNum);};
  virtual const AliCFVGrid*  GetDen() {return (AliCFVGrid*)fContainer->GetGrid(fSelDen);};
  virtual void  SetContainer(const AliCFContainer &c) {fContainer=&c;};

  //basic operations
  virtual void Copy(TObject& eff) const;
 
  
 private:
  const    AliCFContainer *fContainer; //pointer to the input AliContainer
  Int_t fSelNum; //numerator selection step
  Int_t fSelDen; //denominator selection step
  
  ClassDef(AliCFEffGrid,1);
};
    
#endif

