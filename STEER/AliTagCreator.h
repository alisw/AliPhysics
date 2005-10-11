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
#include <TFile.h>
#include <TSystem.h>



class AliTagCreator : public TObject {

 protected:
  TString fUser; //the username in AliEn
  TString fPasswd;   //the username's password
  TString fSE;   //the defined storage element
  TString fCollectionFile; //the xml collection file
  TString fHost; //the defined AliEn host
  Int_t fPort;  //the defined port for the host login
  
  //void CreateTag(TAlienFile* file, const char *guid, Int_t Counter);
  void CreateTag(TFile* file, const char *guid, Int_t Counter);

 public:
  AliTagCreator();
  AliTagCreator(const char *host, Int_t port, const char *username);
  AliTagCreator(const char *host, Int_t port, const char *username, const char *passwd);
  ~AliTagCreator(); 

  void SetSE(const char *se);
  Bool_t ConnectToGrid(const char *host, Int_t port, const char *username);
  Bool_t ReadESDCollection(const char *CollectionFile);
  Bool_t StoreGridTagFile(const char *localpath, const char *gridpath);

  ClassDef(AliTagCreator,0)  
};

#endif

