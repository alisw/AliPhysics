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
  AliMUONGlobalTrigger(Int_t *singlePlus, Int_t *singleMinus,
		       Int_t *singleUndef, Int_t *pairUnlike, Int_t *pairLike);
  virtual ~AliMUONGlobalTrigger();
  AliMUONGlobalTrigger& operator=(const  AliMUONGlobalTrigger& rhs);
        
	/// Return number of Single Plus Low pt
  Int_t SinglePlusLpt()  const {return fSinglePlusLpt;} 
	/// Return number of Single Plus High pt
  Int_t SinglePlusHpt()  const {return fSinglePlusHpt ;}    
	/// Return number of Single Plus All pt
  Int_t SinglePlusApt()  const {return fSinglePlusApt ;}     
	/// Return number of Single Minus Low pt
  Int_t SingleMinusLpt() const {return fSingleMinusLpt ;}  
	/// Return number of Single Minus High pt
  Int_t SingleMinusHpt() const {return fSingleMinusHpt;}  
	/// Return number of Single Minus All pt
  Int_t SingleMinusApt() const {return fSingleMinusApt;}  
	/// Return number of Single Undefined Low pt
  Int_t SingleUndefLpt() const {return fSingleUndefLpt ;}  
	/// Return number of Single Undefined High pt
  Int_t SingleUndefHpt() const {return fSingleUndefHpt ;}   
	/// Return number of Single Undefined All pt
  Int_t SingleUndefApt() const {return fSingleUndefApt ;}  
	/// Return number of Unlike sign pair Low pt
  Int_t PairUnlikeLpt()  const {return fPairUnlikeLpt ;}   
	/// Return number of Unlike sign pair High pt
  Int_t PairUnlikeHpt()  const {return fPairUnlikeHpt ;}   
	/// Return number of Unlike sign pair All pt
  Int_t PairUnlikeApt()  const {return fPairUnlikeApt ;}   
	/// Return number of Like sign pair Low pt
  Int_t PairLikeLpt()    const {return fPairLikeLpt ;}     
	/// Return number of Like sign pair High pt
  Int_t PairLikeHpt()    const {return fPairLikeHpt ;}     
	/// Return number of Like sign pair All pt
  Int_t PairLikeApt()    const {return fPairLikeApt ;}     
  
  void  SetGlobalPattern(Int_t globalPattern);
  void  SetGlobalPattern(UShort_t globalResponse);
  void  SetFromGlobalResponse(UChar_t globalResponse);

  Int_t GetGlobalPattern() const;
  UChar_t GetGlobalResponse() const;

  virtual void Print(Option_t* opt="") const;
  
private:
  Int_t fSinglePlusLpt;  ///< Number of Single Plus Low pt 
  Int_t fSinglePlusHpt;  ///< Number of Single Plus High pt 
  Int_t fSinglePlusApt;  ///< Number of Single Plus All pt 
  Int_t fSingleMinusLpt; ///< Number of Single Minus Low pt
  Int_t fSingleMinusHpt; ///< Number of Single Minus High pt 
  Int_t fSingleMinusApt; ///< Number of Single Minus All pt
  Int_t fSingleUndefLpt; ///< Number of Single Undefined Low pt
  Int_t fSingleUndefHpt; ///< Number of Single Undefined High pt 
  Int_t fSingleUndefApt; ///< Number of Single Undefined All pt
  Int_t fPairUnlikeLpt;  ///< Number of Unlike sign pair Low pt
  Int_t fPairUnlikeHpt;  ///< Number of Unlike sign pair High pt
  Int_t fPairUnlikeApt;  ///< Number of Unlike sign pair All pt
  Int_t fPairLikeLpt;    ///< Number of Like sign pair Low pt
  Int_t fPairLikeHpt;    ///< Number of Like sign pair High pt
  Int_t fPairLikeApt;    ///< Number of Like sign pair All pt

 ClassDef(AliMUONGlobalTrigger,1)  // reconstructed Global Trigger object    
};
#endif






