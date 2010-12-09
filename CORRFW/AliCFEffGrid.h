#ifndef ALICFEFFGRID_H
#define ALICFEFFGRID_H

/* $Id$ */

//-----------------------------------------------------------------//
//                                                                 //
// AliCFEffGrid Class                                              //
// Class to handle efficiency grids                                // 
//                                                                 //
//-----------------------------------------------------------------//

#include "AliCFGridSparse.h"
#include "AliCFContainer.h"
class TH1D;
class TH2D;
class TH3D;

class AliCFEffGrid : public AliCFGridSparse
{
 public:
  AliCFEffGrid();
  AliCFEffGrid(const Char_t* name,const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn);
  AliCFEffGrid(const Char_t* name,const Char_t* title, const AliCFContainer &c);
  virtual ~AliCFEffGrid();
  AliCFEffGrid(const AliCFEffGrid& eff);
  AliCFEffGrid& operator=(const AliCFEffGrid& eff);
  
  virtual Double_t GetAverage() const ;
  virtual Int_t GetSelNumStep() const {return fSelNum;};
  virtual Int_t GetSelDenStep() const {return fSelDen;};
  virtual TH1*  Project(Int_t ivar1, Int_t ivar2=-1,Int_t ivar3=-1) const;
  virtual AliCFGridSparse*  Project(Int_t, const Int_t*, const Double_t*, const Double_t*, Bool_t) const {AliWarning("should not be used"); return 0x0;}
  virtual AliCFEffGrid* MakeSlice(Int_t nVars, const Int_t* vars, const Double_t* varMin, const Double_t* varMax, Bool_t useBins=0) const;

  //Efficiency calculation
  virtual void  CalculateEfficiency(Int_t istep1, Int_t istep2, Option_t *option ="B" /*binomial*/);
  virtual AliCFGridSparse*  GetNum() const {return fContainer->GetGrid(fSelNum);};
  virtual AliCFGridSparse*  GetDen() const {return fContainer->GetGrid(fSelDen);};
  virtual void  SetContainer(const AliCFContainer &c) {fContainer=&c;};

 private:
  const AliCFContainer *fContainer; //pointer to the input AliContainer
  Int_t fSelNum;                    //numerator selection step
  Int_t fSelDen;                    //denominator selection step
  
  ClassDef(AliCFEffGrid,1);
};
    
#endif
