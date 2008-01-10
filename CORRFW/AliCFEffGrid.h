#ifndef ALICFEFFGRID_H
#define ALICFEFFGRID_H

/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFEffGrid Class                                              //
// Class to handle efficiency grids                                   // 
//                                                                    //
//--------------------------------------------------------------------//

#include "AliCFGrid.h"
#include "AliCFContainer.h"
#include <TNamed.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

class AliCFEffGrid : public AliCFGrid
{
 public:
  AliCFEffGrid();
  AliCFEffGrid(const Char_t* name,const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn, const Float_t  *binLimitsIn=0);
  AliCFEffGrid(const Char_t* name,const Char_t* title,const AliCFContainer &c);
  AliCFEffGrid(const AliCFEffGrid& eff);
  
  virtual ~AliCFEffGrid();
  AliCFEffGrid& operator=(const AliCFEffGrid& eff);
  virtual Float_t GetAverage() const ;
  virtual Float_t GetAverage(Float_t *varMin,Float_t *varMax) const ;
  virtual Int_t GetSelNumStep() const {return fSelNum;};
  virtual Int_t GetSelDenStep() const {return fSelDen;};
  virtual TH1F* Project( Int_t ivar) const;
  virtual TH2F* Project( Int_t ivar1, Int_t ivar2) const;
  virtual TH3F* Project( Int_t ivar1, Int_t ivar2,Int_t ivar3) const;

  //Efficiency calculation
  virtual void  CalculateEfficiency(Int_t istep1, Int_t istep2);
  virtual const AliCFGrid*  GetNum() {return fContainer->GetGrid(fSelNum);};
  virtual const AliCFGrid*  GetDen() {return fContainer->GetGrid(fSelDen);};
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

