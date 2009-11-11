#ifndef ALIGLAUBERMC_H
#define ALIGLAUBERMC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////////////
//
//  Glauber MC
//
//  origin: PHOBOS experiment
//  alification: Mikolaj Krzewicki, Nikhef, mikolaj.krzewicki@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

#include "AliGlauberNucleus.h"

class TNamed;
class TObjArray;
class TNtuple;

class AliGlauberMC : public TNamed {
private:
   AliGlauberNucleus fANucleus;       //Nucleus A
   AliGlauberNucleus fBNucleus;       //Nucleus B
   Double_t     fXSect;          //Nucleon-nucleon cross section
   TObjArray*   fNucleonsA;      //Array of nucleons in nucleus A
   TObjArray*   fNucleonsB;      //Array of nucleons in nucleus B
   Int_t        fAN;             //Number of nucleons in nucleus A
   Int_t        fBN;             //Number of nucleons in nucleus B
   TNtuple*     fnt;             //Ntuple for results (created, but not deleted)
   Double_t     fMeanX2;         //<x^2> of wounded nucleons
   Double_t     fMeanY2;         //<y^2> of wounded nucleons
   Double_t     fMeanXY;         //<xy> of wounded nucleons
   Double_t     fMeanXParts;     //<x> of wounded nucleons
   Double_t     fMeanYParts;     //<x> of wounded nucleons
   Double_t     fMeanXSystem;    //<x> of all nucleons
   Double_t     fMeanYSystem;    //<x> of all nucleons  
   Double_t     fMeanX_A;        //<x> of nucleons in nucleus A
   Double_t     fMeanY_A;        //<x> of nucleons in nucleus A
   Double_t     fMeanX_B;        //<x> of nucleons in nucleus B
   Double_t     fMeanY_B;        //<x> of nucleons in nucleus B
   Double_t     fB_MC;           //Impact parameter (b)
   Int_t        fEvents;         //Number of events with at least one collision
   Int_t        fTotalEvents;    //All events within selected impact parameter range
   Double_t     fBMin;           //Minimum impact parameter to be generated
   Double_t     fBMax;           //Maximum impact parameter to be generated
   Int_t        fMaxNpartFound;  //Largest value of Npart obtained
   Int_t        fNpart;          //Number of wounded (participating) nucleons in current event
   Int_t        fNcoll;          //Number of binary collisions in current event
   Double_t     fSx2;            //Variance of x of wounded nucleons
   Double_t     fSy2;            //Variance of y of wounded nucleons
   Double_t     fSxy;            //Covariance of x and y of wounded nucleons

   Bool_t       CalcResults(Double_t bgen);
   Bool_t       CalcEvent(Double_t bgen);

public:
   AliGlauberMC(Option_t* NA = "Au", Option_t* NB = "Au", Double_t xsect = 42);
   virtual     ~AliGlauberMC() {Reset();}

   void         Draw(Option_t* option);
   Double_t     GetB()               const {return fB_MC;}
   Double_t     GetBMin()            const {return fBMin;}
   Double_t     GetBMax()            const {return fBMax;}
   Int_t        GetNcoll()           const {return fNcoll;}
   Int_t        GetNpart()           const {return fNpart;}
   Int_t        GetNpartFound()      const {return fMaxNpartFound;}
   TNtuple*     GetNtuple()          const {return fnt;}
   TObjArray   *GetNucleons();
   Double_t     GetTotXSect()        const;
   Double_t     GetTotXSectErr()     const;
   Bool_t       NextEvent(Double_t bgen=-1);
   void         Reset()                    {delete fnt; fnt=0; }
   void         Run(Int_t nevents);
   void         SetBmin(Double_t bmin)      {fBMin = bmin;}
   void         SetBmax(Double_t bmax)      {fBMax = bmax;}
   void         SetMinDistance(Double_t d)  {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}

   static void       PrintVersion()         {cout << "AliGlauberMC " << Version() << endl;}
   static const char *Version()             {return "v1.1";}
   static void       runAndSaveNucleons( Int_t n,                    
                                         Option_t *sysA="Au",           
                                         Option_t *sysB="Au",           
                                         Double_t signn=42,           
                                         Double_t mind=0.4,
                                         Bool_t verbose=0,
                                         const char *fname="glau_auau_nucleons.root" );
   static void       runAndSaveNtuple( Int_t n,
                                       Option_t *sysA="Au",
                                       Option_t *sysB="Au",
                                       Double_t signn=42,
                                       Double_t mind=0.4,
                                       const char *fname="glau_auau_ntuple.root");

   ClassDef(AliGlauberMC,1)
};

#endif
