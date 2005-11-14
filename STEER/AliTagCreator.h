#ifndef ALITAGCREATOR_H
#define ALITAGCREATOR_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliTagCreator
//   This is the AliTagCreator class for the tag creation (post process)
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        AliTagCreator                                 //
//                                                                      //
//           Implementation of the tag creation mechanism.              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//ROOT
#include <TObject.h>

class TFile;
class TGridResult;


class AliTagCreator : public TObject {

 public:
  AliTagCreator();
  AliTagCreator(const char *host, Int_t port, const char *username);
  AliTagCreator(const char *host, Int_t port, const char *username, const char *passwd);
  ~AliTagCreator(); 

  void SetSE(const char *se){fSE = se;}
  void SetStorage(Int_t storage);
  void SetGridPath(const char *gridpath){fgridpath = gridpath;}
  Bool_t ConnectToGrid(const char *host, Int_t port, const char *username);
  Bool_t ReadESDCollection(TGridResult *result);

 protected:
  TString fUser; //the username in AliEn
  TString fPasswd;   //the username's password
  TString fSE;   //the defined storage element
  TString fHost; //the defined AliEn host
  TString fgridpath;   //the alien location of the tag files
  Int_t fPort;  //the defined port for the host login
  Int_t fStorage;  //0:local - 1:grid
  //  TGridResult *fresult; //the results from the grid query

  void CreateTag(TFile* file, const char *guid, Int_t Counter);

  ClassDef(AliTagCreator,0)  
};

#endif

