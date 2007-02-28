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

class AliGRPDCS: public TObject {
 public:
  AliGRPDCS();
  AliGRPDCS(TObjArray *dcsArray);
  AliGRPDCS(const AliGRPDCS& grpDcs);

  void SetObjArray(TObjArray *dcsSArray) {fDCSArray = dcsSArray;}
  const char *ProcessDCS(Int_t iType);  
  
 private:
  
  TObjArray *fDCSArray; //TObjArray for a dcs data point
  
  const char *ProcessInt();
  const char *ProcessUInt();
  const char *ProcessFloat();
  const char *ProcessString();
  const char *ProcessBoolean();
  AliGRPDCS & operator=(const AliGRPDCS & ) {return *this;}

  ClassDef(AliGRPDCS, 0);
};

#endif
