//=================================================
// AliJBaseCard.h
// last modified FK 6.NOV 2009
//=================================================

#ifndef ALIJBASECARD_H
#define ALIJBASECARD_H

#include <TObject.h>

#include <iostream>
#include <fstream>

#include <TString.h>
#include <TVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TFile.h>
#include <TF1.h>
#include <vector>
#include <THashList.h>
#include <TNamed.h>

using namespace std;
#ifndef AliJMaxDimBuffer
#define AliJMaxDimBuffer
const int kMaxDimBuffer = 300;//max length of a line read to a buffe
#endif

class AliJBaseCard {

  //====   M e m b e r    F u n c t i o n s   ========

public:

  AliJBaseCard(); // constructor
  AliJBaseCard(const char *filename); // constructor
  AliJBaseCard& operator=(const AliJBaseCard& obj);

  virtual ~AliJBaseCard();

  void AddToKeyTable( TString key, int index ){
    TNamed * ko = new TNamed(key.Data(),"");
    ko->SetUniqueID(index);
    fKeyTable.Add(ko);// TODO check duplicate
  }

  float  Get(TString keyword, int VectorComponent=0); //get TVector component 
  TString  GetStr(TString keyword ); //get TVector component 
  TVector* GetVector( TString keyword );
  int    GetN(TString keyword);       //get TVector dimension
  void   PrintOut(); 
  void   WriteCard(TDirectory *file);  

  virtual void InitCard();
  void FinishCard();
  void ReadInputLine( const char* buffer );
  void ReadLine( const char * buffer );

protected:
  void    ReadInputCard();

  int     GetNwithIndex(int i){ return fValuesVector[i].GetNrows(); }

  unsigned int GetTVectorIndex(TString keyword, int tol=0);

  //====   D a t a    M e m b e r s  ========

  char fcardname[255];   // file name
  int  fnentry;                   //Number of lines in cfg file
  std::vector< TString >   fKeyWordVector;    //array of key words
  std::vector< TVector >   fValuesVector;     //array of float number confg parameter vectors 
  std::vector< TString >   fValueString;      // Storage of raw inut string for each item
  //std::map< TString, unsigned int > MapKeyWordToTVector;//mapping keywords to TVector
  THashList           fKeyTable;              // key map with hash algorithm

  //ClassDef(AliJBaseCard, 1); // EMCAL for jcorran
};

#endif






















