#ifndef ALICFFRAME_H
#define ALICFFRAME_H

/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFFrame.cxx Class                                               //
// Class to handle input data for correction Framework                // 
//                                                                    //
//--------------------------------------------------------------------//

#include <TNamed.h>

class AliCFFrame : public TNamed
{
 public:
  AliCFFrame();
  AliCFFrame(const Char_t* name,const Char_t* title);
  AliCFFrame(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t* nBinIn, const Double_t  *binLimitsIn=0);
  AliCFFrame(const AliCFFrame & c);
  
  virtual ~AliCFFrame();
  AliCFFrame& operator=(const AliCFFrame& corr);
  virtual void  PrintBinLimits();
  virtual void  PrintNBins();
  virtual void  SetBinLimits(Int_t ivar, Double_t * array);
  virtual void  GetBinLimits(Int_t ivar, Double_t * array) const;
  virtual Int_t GetBinIndex(Int_t *ibin) const;
  virtual void  GetBinIndex(Int_t iel, Int_t *ibin) const;
  virtual Int_t GetBinIndex(Int_t ivar, Int_t ind) const;
  virtual Int_t GetNDim() const {return fNDim;};
  virtual Int_t GetNVar() const {return fNVar;};
  virtual Int_t GetNBins(Int_t ivar) const {return fNVarBins[ivar];};
  virtual Int_t GetNBinLimits() const {return fNVarBinLimits;};
  virtual Int_t   *GetNBins() const {return fNVarBins;};
  virtual Double_t *GetBinLimits() const {return fVarBinLimits;};
  virtual Double_t GetBinCenter(Int_t ivar,Int_t ibin) const;
  virtual Double_t GetBinSize(Int_t ivar,Int_t ibin) const;
  virtual void    GetBinCenters(Int_t *ibin, Double_t *binCenter) const;
  virtual void    GetBinSizes(Int_t *ibin, Double_t *binSizes) const;

  //basic operations

  virtual void Copy(TObject& c) const;
  virtual void Save(const Char_t *outfile) const;
  
 protected:
  Int_t    fNVar; //number of variables in the grid
  Int_t    fNDim; //Overall number of elements in the grid
  Int_t    fNVarBinLimits; //total number of bin limits
  Int_t    *fNVarBins; //[fNVar] size of the grid in each dimension (binning)
  Int_t    *fIndex;//[fNVar] current N-dim index on the grid
  Int_t    *fProduct;//[fNVar] current N-dim index on the grid
  Int_t    *fOffset;//[fNVar] current N-dim index on the grid
  Double_t  *fVarBinLimits;//[fNVarBinLimits] array defining the binLimits

  
  ClassDef(AliCFFrame,2);
};
    
#endif

