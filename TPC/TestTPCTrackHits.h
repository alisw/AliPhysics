/// \class AliTPChitD
/// \brief Macro to compare TClonesArray hits with interpolated hits
/// \author MI

#ifndef TESTTPCTRACKHITS_H
#define TESTTPCTRACKHITS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

void ConvertHits(const char * benchmark="0", Bool_t debug=kFALSE);
void CompareHits(const char * benchmark="1", Bool_t debug=kFALSE);

class AliTPChitD : public AliTPChit {
public:
  AliTPChit * GetDelta() {return &fDelta;}
private:
  AliTPChit fDelta;     ///< delta of hit information
  ClassDef(AliTPChitD,1) 
};
#endif
