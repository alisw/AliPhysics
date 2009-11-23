#ifndef ALICFFRAME_H
#define ALICFFRAME_H

/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFFrame.cxx Class                                               //
// Class to handle input data for correction Framework                // 
//                                                                    //
//--------------------------------------------------------------------//

#include "TNamed.h"

class AliCFFrame : public TNamed
{
 public:
  AliCFFrame();
  AliCFFrame(const Char_t* name,const Char_t* title);
  virtual ~AliCFFrame() {} ;
  virtual void Copy(TObject& c) const {TNamed::Copy(c);}
  virtual void Save(const Char_t *outfile) const;
  virtual Int_t      GetNVar()                                                   const = 0 ; // number of variables
  virtual void       PrintBinLimits()                                            const = 0 ; // prints the bin limits for each variable
  virtual void       PrintNBins()                                                const = 0 ; // prints the number of bins for each variable
  virtual void       SetBinLimits(Int_t ivar, const Double_t * array)                  = 0 ; // sets the bin limits specified in array
  virtual void       GetBinLimits(Int_t ivar, Double_t * array)                  const = 0 ; // puts in array the bin limits for variable ivar
  virtual Double_t * GetBinLimits(Int_t ivar)                                    const = 0 ; // returns an array of bin limits for variable ivar
  virtual Long_t     GetNBinsTotal()                                             const = 0 ; // total number of bins 
  virtual Int_t      GetNBins(Int_t ivar)                                        const = 0 ; // number of bins for variable ivar
  virtual Int_t    * GetNBins()                                                  const = 0 ; // returns an array containing the bins for each variable
  virtual Float_t    GetBinCenter(Int_t ivar,Int_t ibin)                         const = 0 ; // the center of bin number ibin for variable ivar
  virtual Float_t    GetBinSize  (Int_t ivar,Int_t ibin)                         const = 0 ; // the   size of bin number ibin for variable ivar

  //  virtual void       Clear() = 0 ; // clear all the cells
  
  //virtual void       GetBinCenters(const Int_t *ibin, Float_t *binCenter)        const = 0 ; // 
  //virtual void       GetBinSizes  (const Int_t *ibin, Float_t *binSizes)         const = 0 ; //

  // probably not needed anymore
/*   virtual Int_t      GetBinIndex(const Int_t *ibin)                              const = 0 ; */
/*   virtual void       GetBinIndex(Int_t iel, const Int_t *ibin)                   const = 0 ; */
/*   virtual Int_t      GetBinIndex(Int_t ivar, Int_t ind)                          const = 0 ; */

  ClassDef(AliCFFrame,3);
};
    
#endif

