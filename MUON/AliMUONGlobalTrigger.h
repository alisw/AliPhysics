#ifndef ALIMUONGLOBALTRIGGER_H
#define ALIMUONGLOBALTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup base
/// \class AliMUONGlobalTrigger
/// \brief Global trigger object
//  Author Ph. Crochet

#include <TObject.h>

class AliMUONGlobalTrigger : public TObject {
 public:
  AliMUONGlobalTrigger();
  AliMUONGlobalTrigger(const AliMUONGlobalTrigger& rhs); // copy constructor
  virtual ~AliMUONGlobalTrigger();
  AliMUONGlobalTrigger& operator=(const  AliMUONGlobalTrigger& rhs);
        
	/// Return number of Single Low pt
  Int_t SingleLpt()  const {return fSingleLpt;} 
	/// Return number of Single High pt
  Int_t SingleHpt()  const {return fSingleHpt ;}    
	/// Return number of Unlike sign pair Low pt
  Int_t PairUnlikeLpt()  const {return fPairUnlikeLpt ;}   
	/// Return number of Unlike sign pair High pt
  Int_t PairUnlikeHpt()  const {return fPairUnlikeHpt ;}   
	/// Return number of Like sign pair Low pt
  Int_t PairLikeLpt()    const {return fPairLikeLpt ;}     
	/// Return number of Like sign pair High pt
  Int_t PairLikeHpt()    const {return fPairLikeHpt ;}     
  
  void  SetFromGlobalResponse(UShort_t globalResponse);
  UChar_t GetGlobalResponse() const;

  virtual void Print(Option_t* opt="") const;
  
private:
  Int_t fSingleLpt;      ///< Number of Single Low pt 
  Int_t fSingleHpt;      ///< Number of Single High pt 
  Int_t fPairUnlikeLpt;  ///< Number of Unlike sign pair Low pt
  Int_t fPairUnlikeHpt;  ///< Number of Unlike sign pair High pt
  Int_t fPairLikeLpt;    ///< Number of Like sign pair Low pt
  Int_t fPairLikeHpt;    ///< Number of Like sign pair High pt

 ClassDef(AliMUONGlobalTrigger,2)  // reconstructed Global Trigger object    
};
#endif






