#ifndef ALIRECINFOCUTS_H
#define ALIRECINFOCUTS_H

//------------------------------------------------------------------------------
// Class to keep selection cuts for reconstructed tracks. 
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

#include "AliESDtrackCuts.h"

class AliRecInfoCuts : public AliESDtrackCuts
{
public:
  AliRecInfoCuts(const Char_t* name ="AliRecInfoCuts", const Char_t *title ="");
  virtual ~AliRecInfoCuts() {;}
 
  // setters 
  void SetMinTPCsignalN(const Int_t min=0) 	 {fMinTPCsignalN = min;}
  void SetMaxAbsTanTheta(const Float_t max=1e99)  {fMaxAbsTanTheta = max;}

  // getters
  Int_t GetMinTPCsignalN()    const {return fMinTPCsignalN;}
  Float_t GetMaxAbsTanTheta() const {return fMaxAbsTanTheta;}

  // getters for selected AliESDtrackCuts data members
  Float_t GetPtMin()          const {return fPtMin;}
  Float_t GetPtMax()          const {return fPtMax;}
  Int_t GetMinNClustersTPC()  const {return fCutMinNClusterTPC;}

  // cuts init function
  void Init();

private:
  Int_t   fMinTPCsignalN;  // min. number of TPC hits
  Float_t fMaxAbsTanTheta; // max. absolute value of tan(theta)

  AliRecInfoCuts(const AliRecInfoCuts&); // not implemented
  AliRecInfoCuts& operator=(const AliRecInfoCuts&); // not implemented

  ClassDef(AliRecInfoCuts, 1)
};

#endif //ALIRECINFOCUTS_H
