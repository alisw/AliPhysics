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

#include <Riostream.h>
class TNamed;
class TObjArray;
class TNtuple;
class TArray;

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
   Double_t     fMeanYParts;     //<y> of wounded nucleons
   Double_t     fMeanXColl;      //<x> of binary collisions
   Double_t     fMeanYColl;      //<y> of binary collisions
   Double_t     fMeanX2Coll;     //<x^2> of binary collisions
   Double_t     fMeanY2Coll;     //<y^2> of binary collisions
   Double_t     fMeanXYColl;     //<xy> of binary collisions
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
   Double_t	    fdNdEtaParam[2];	   //Parameters: nnp, x
   Double_t     fdNdEtaGBWParam[3];  //Parameters: delta, lambda, snn
   Double_t     fdNdEtaNBDParam[3];       //Parameters:  k, nmean, beta
   Double_t     fdNdEtaTwoNBDParam[6];    //Parameters: k1, nmean1, k2, nmean2, alpha, beta
   Int_t        fMaxNpartFound;  //Largest value of Npart obtained
   Int_t        fNpart;          //Number of wounded (participating) nucleons in current event
   Int_t        fNcoll;          //Number of binary collisions in current event
   Double_t     fSx2;            //Variance of x of wounded nucleons
   Double_t     fSy2;            //Variance of y of wounded nucleons
   Double_t     fSxy;            //Covariance of x and y of wounded nucleons
   Double_t     fSx2Coll;            //Variance of x of binaruy collisions
   Double_t     fSy2Coll;            //Variance of y of binaruy collisions
   Double_t     fSxyColl;            //Covariance of x and y of binaruy collisions
   Double_t     fX;              //hard particle production fraction
   Double_t     fNpp;            //Multiplicity normalization
   Bool_t       CalcResults(Double_t bgen);

public:
   AliGlauberMC(Option_t* NA = "Pb", Option_t* NB = "Pb", Double_t xsect = 72);
   virtual     ~AliGlauberMC();
   AliGlauberMC(const AliGlauberMC& in);
   AliGlauberMC& operator=(const AliGlauberMC& in);
   void         Draw(Option_t* option);

   void         Run(Int_t nevents);
   Bool_t       NextEvent(Double_t bgen=-1);
   Bool_t       CalcEvent(Double_t bgen);

   //various ways to calculate multiplicity
   Double_t     GetdNdEta(    Double_t nnp=8.0,
                              Double_t x=0.13 );
   Double_t     GetdNdEtaGBW( Double_t delta=0.79,
                              Double_t lambda=0.288,
                              Double_t snn=30.25 );
   Double_t	GetdNdEtaNBD(     Int_t k=3,
                              Double_t nmean = 4,
                              Double_t beta = 0.13 );
   Double_t	GetdNdEtaTwoNBD(  Int_t k1=3,
                              Double_t nmean1=4,
                              Int_t k2=2,
                              Double_t nmean2=11,
                              Double_t alpha=0.4,
                              Double_t beta=0.13 );  

   Double_t     GetEccentricity()    const;
   Double_t     GetEccentricityColl()      const;
   Double_t     GetEccentricityPart()      const;
   Double_t     GetEccentricityPartColl()  const;
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
   void         Reset();
   static Double_t	NegativeBinomialDistribution(Int_t x, Int_t k, Double_t nmean);
   Int_t NegativeBinomialRandom(Int_t k, Double_t nmean);
   Int_t DoubleNegativeBinomialRandom(Int_t k1, Double_t nmean1, Int_t k2, Double_t nmean2, Double_t alpha);
   void   SetBmin(Double_t bmin)      {fBMin = bmin;}
   void   SetBmax(Double_t bmax)      {fBMax = bmax;}
   void   SetdNdEtaParam( Double_t nnp = 8., Double_t x = 0.13);
   void   SetdNdEtaGBWParam( Double_t delta = 0.79, Double_t lambda = 0.288, Double_t snn = 30.25);
   void   SetdNdEtaNBDParam(Double_t k=3, Double_t nmean=4, Double_t beta=0.13);
   void   SetdNdEtaTwoNBDParam(Double_t alpha = 0.4, Double_t k1 = 3, Double_t nmean1 = 4., Double_t k2 = 2., Double_t nmean2 = 11., Double_t beta=0.13);
   void   SetMinDistance(Double_t d)  {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}
   static void       PrintVersion()         {cout << "AliGlauberMC " << Version() << endl;}
   static const char *Version()             {return "v1.2";}
   static void       runAndSaveNtuple( Int_t n,
                                       Option_t *sysA="Au",
                                       Option_t *sysB="Au",
                                       Double_t signn=42,
                                       Double_t mind=0.4,
                                       const char *fname="glau_auau_ntuple.root");
   void runAndSaveNucleons( Int_t n,
                            Option_t *sysA,
                            Option_t *sysB,
                            Double_t signn,
                            Double_t mind,
                            Bool_t verbose,
                            const char *fname);
   
   ClassDef(AliGlauberMC,3)
};

#endif
