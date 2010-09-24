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
  AliCFDataGrid(const Char_t* name,const Char_t* title, const AliCFContainer &c, Int_t step);      //create data grid from container
  AliCFDataGrid(const Char_t* name,const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn);  //create empty data grid to fill it yourself
  AliCFDataGrid(const AliCFDataGrid& c);
  AliCFDataGrid& operator=(const AliCFDataGrid& c);
  virtual ~AliCFDataGrid();

  virtual Int_t GetSelDataStep() const {return fSelData;};

  // Methods for handling/correcting data 

  virtual const AliCFGridSparse*  GetData() {return fContainer->GetGrid(fSelData);};
  virtual void  ApplyEffCorrection(const AliCFEffGrid &eff);
  virtual void  ApplyBGCorrection(const AliCFDataGrid &c);
  //basic operations
  virtual void   Copy(TObject& data) const;
 
  
 private:
  Int_t fSelData; //sel step of the observed data 
  const AliCFContainer *fContainer; //pointer to the input data Container
  ClassDef(AliCFDataGrid,2);
};
    
#endif
