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

class TH1;

class AliGRPDCS: public TObject {
 public:
  AliGRPDCS();
  AliGRPDCS(TObjArray *dcsArray);
  AliGRPDCS(const AliGRPDCS& grpDcs);

  void SetObjArray(TObjArray *dcsSArray) {fDCSArray = dcsSArray;}
  const char *ProcessDCS(TH1 *h);  
  
 private:
  
  TObjArray *fDCSArray; //TObjArray for a dcs data point
  
  AliGRPDCS & operator=(const AliGRPDCS & ) {return *this;}

  ClassDef(AliGRPDCS, 0);
};

#endif
