#ifndef ALICFDATAGRID_H
#define ALICFDATAGRID_H

/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFDataGrid Class                                             //
// Class to handle observed data and correct them                     // 
//                                                                    //
//--------------------------------------------------------------------//

#include "AliCFGridSparse.h"
#include "AliCFEffGrid.h"
#include "AliCFContainer.h"

class AliCFDataGrid : public AliCFGridSparse
{
 public:
  AliCFDataGrid();
  AliCFDataGrid(const Char_t* name,const Char_t* title);
  AliCFDataGrid(const Char_t* name,const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn, const Double_t  *binLimitsIn=0);
  AliCFDataGrid(const Char_t* name,const Char_t* title,const AliCFContainer &c);
  AliCFDataGrid(const AliCFDataGrid& c);
  
  virtual ~AliCFDataGrid();
  AliCFDataGrid& operator=(const AliCFDataGrid& c);
  virtual Int_t GetSelDataStep() const {return fSelData;};

  // Methods for handling/correcting data 

  virtual void  SetMeasured(Int_t istep);
  virtual const AliCFVGrid*  GetData() {return (AliCFVGrid*)fContainer->GetGrid(fSelData);};
  virtual void  ApplyEffCorrection(const AliCFEffGrid &eff);
  virtual void  ApplyBGCorrection(const AliCFDataGrid &c);
  virtual void  SetContainer(const AliCFContainer &c) {fContainer=&c;};
  //basic operations
  virtual void   Copy(TObject& data) const;
 
  
 private:
  Int_t fSelData; //sel step of the observed data 
  const AliCFContainer *fContainer; //pointer to the input data Container
  ClassDef(AliCFDataGrid,2);
};
    
#endif

