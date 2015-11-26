#ifndef ALI_GRP_DCS_H
#define ALI_GRP_DCS_H

//-------------------------------------------------------------------------
//                          Class AliGRPDCS
//   This class deals with the DCS related info of the GRP
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliGRPDCS                                     //
//                                                                      //
//           Implementation of the class that processes                 //
//           the DCS related fields of the GRP.                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "TString.h"

class AliGRPDCS: public TObject {
 public:
  AliGRPDCS();
  AliGRPDCS(TObjArray *dcsArray, UInt_t fStart, UInt_t fStop);
  AliGRPDCS(const AliGRPDCS& grpDcs);

  void SetTime(UInt_t fStart, UInt_t fStop) {fStartTime = fStart; fStopTime = fStop;}
  void SetObjArray(TObjArray *dcsSArray) {fDCSArray = dcsSArray;}
  TString ProcessDCS(Int_t iType);  
  
 private:
  UInt_t fStartTime, fStopTime; //start and stop time of the run (DAQ lb)
  TObjArray *fDCSArray; //TObjArray for a dcs data point
  
  TString ProcessInt();
  TString ProcessUInt();
  TString ProcessFloat();
  TString ProcessChar();
//  TString ProcessString();
  TString ProcessBoolean();
  
  AliGRPDCS & operator=(const AliGRPDCS & ) {return *this;}

  ClassDef(AliGRPDCS, 0);
};

#endif
