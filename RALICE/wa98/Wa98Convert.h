#ifndef Wa98Convert_h
#define Wa98Convert_h

// $Id$

#include "TObject.h"
#include "TChain.h"
#include "TFile.h"
#include "THbookFile.h"

#include "Wa98Event.h"

class Wa98Convert : public TObject
{
 public :
  Wa98Convert(TTree* tree=0);                                    // Constructor
  virtual ~Wa98Convert();                                        // Destructor
  void Loop(TTree* otree=0,Int_t nentries=-1,Int_t printfreq=1); // Perform the conversion

 protected :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   //Declaration of leaves types
   Int_t           Jrun;
   Int_t           Jevt;
   Int_t           Jdate;
   Int_t           Jtime;
   Int_t           Jevid;
   UInt_t          Jwscal;
   UInt_t          Itword;
   Float_t         Zdc;
   Float_t         Emir;
   Float_t         Emire;
   Float_t         Emirh;
   Float_t         Etm;
   Float_t         Etme;
   Float_t         Etmh;
   Int_t           Nmod;
   Int_t           Irowl[3000]; //[Nmod]
   UInt_t          Icoll[3000]; //[Nmod]
   Float_t         Adcl[3000];  //[Nmod]
   Int_t           Nclu;
   Int_t           Irowc[400];  //[Nclu]
   UInt_t          Icolc[400];  //[Nclu]
   Float_t         Adcc[400];   //[Nclu]
   Int_t           Ncluv;
   Int_t           Iadccv[1000];  //[Ncluv]
   Float_t         Thetacv[1000]; //[Ncluv]
   Float_t         Phicv[1000];   //[Ncluv]

 ClassDef(Wa98Convert,1) // Conversion Wa98 ntuple data into Wa98Event physics event structures.
};
#endif
