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
//  update:      You Zhou, Nikhef, yzhou@nikhef.nl :)
////////////////////////////////////////////////////////////////////////////////

#include "AliGlauberNucleus.h"
#include <Riostream.h>
#include <TNamed.h>

class TObjArray;
class TNtuple;

class AliGlauberMC : public TNamed {
public:
   enum EdNdEtaType { kSimple,
                      kNBD,
                      kNBDSV,
                      kTwoNBD,
                      kGBW,
                      kNone };

   AliGlauberMC(Option_t* NA = "Pb", Option_t* NB = "Pb", Double_t xsect = 64);
   virtual     ~AliGlauberMC();
   AliGlauberMC(const AliGlauberMC& in);
   AliGlauberMC& operator=(const AliGlauberMC& in);
   void         Draw(Option_t* option);

   void         Run(Int_t nevents);
   Bool_t       NextEvent(Double_t bgen=-1);
   Bool_t       CalcEvent(Double_t bgen);

   //various ways to calculate multiplicity
   Double_t GetdNdEta() const;
   Double_t GetdNdEtaSimple( const Double_t* param ) const;
   Double_t GetdNdEtaGBW(    const Double_t* param ) const;
   Double_t	GetdNdEtaNBD(    const Double_t* param ) const;
   Double_t	GetdNdEtaNBDSV(  const Double_t* param ) const;
   Double_t	GetdNdEtaTwoNBD( const Double_t* param ) const;

   Double_t     GetEccentricity()    const;
   Double_t     GetEccentricityColl()      const;
   Double_t     GetEccentricityPart()      const;
   Double_t     GetEpsilon2Part()      const;
   Double_t     GetEpsilon3Part()      const;
   Double_t     GetEpsilon4Part()      const;
   Double_t     GetEpsilon5Part()      const;
   Double_t     GetEpsilon2Coll()      const;
   Double_t     GetEpsilon3Coll()      const;
   Double_t     GetEpsilon4Coll()      const;
   Double_t     GetEpsilon5Coll()      const;
   Double_t     GetEccentricityPartColl()  const;
   Double_t     GetB()               const {return fBMC;}
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
   Int_t  NegativeBinomialRandom(Int_t k, Double_t nmean) const;
   Int_t  NegativeBinomialRandomSV(Double_t nbar, Double_t k) const;
   Int_t  DoubleNegativeBinomialRandom(Int_t k1, Double_t nmean1, Int_t k2, Double_t nmean2, Double_t alpha) const;
   Double_t* GetdNdEtaParam() {return fdNdEtaParam;}
   void      SetdNdEtaType(EdNdEtaType method) {fMultType=method;}
   void   SetBmin(Double_t bmin)      {fBMin = bmin;}
   void   SetBmax(Double_t bmax)      {fBMax = bmax;}
   void   SetMinDistance(Double_t d)  {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}
   void   SetDoPartProduction(Bool_t b) { fDoPartProd = b; }
   void   Setr(Double_t r)  {fANucleus.SetR(r); fBNucleus.SetR(r);}
   void   Seta(Double_t a)  {fANucleus.SetA(a); fBNucleus.SetA(a);}
   static void       PrintVersion()         {cout << "AliGlauberMC " << Version() << endl;}
   static const char *Version()             {return "v1.2";}
   static void       RunAndSaveNtuple( Int_t n,
                                       const Option_t *sysA="Pb",
                                       const Option_t *sysB="Pb",
                                       Double_t signn=64,
                                       Double_t mind=0.4,
				       Double_t r=6.62,
				       Double_t a=0.546,
                                       const char *fname="glau_pbpb_ntuple.root");
   void RunAndSaveNucleons( Int_t n,
                            const Option_t *sysA,
                            const Option_t *sysB,
                            Double_t signn,
                            Double_t mind,
                            Bool_t verbose,
                            const char *fname);
   
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
   Double_t     fMeanXA;        //<x> of nucleons in nucleus A
   Double_t     fMeanYA;        //<x> of nucleons in nucleus A
   Double_t     fMeanXB;        //<x> of nucleons in nucleus B
   Double_t     fMeanYB;        //<x> of nucleons in nucleus B 
   Double_t     fBMC;           //Impact parameter (b)
   Int_t        fEvents;         //Number of events with at least one collision
   Int_t        fTotalEvents;    //All events within selected impact parameter range
   Double_t     fBMin;           //Minimum impact parameter to be generated
   Double_t     fBMax;           //Maximum impact parameter to be generated
   Double_t	    fdNdEtaParam[10];//Parameters for multiplicity calculation: meaning depends on method selection
   EdNdEtaType fMultType;//mutliplicity method selection  
   Int_t        fMaxNpartFound;  //Largest value of Npart obtained
   Int_t        fNpart;          //Number of wounded (participating) nucleons in current event
   Int_t        fNcoll;          //Number of binary collisions in current event
   Double_t     fMeanr2;         //----------<r^2> of wounded nucleons
   Double_t     fMeanr3;         //----------<r^3> of wounded nucleons
   Double_t     fMeanr4;         //----------<r^4> of wounded nucleons
   Double_t     fMeanr5;         //----------<r^5> of wounded nucleons
   Double_t     fMeanr2Cos2Phi;   //------<r^2*cos2phi> of wounded nucleons
   Double_t     fMeanr2Sin2Phi;   //------<r^2*sin2phi> of wounded nucleons
   Double_t     fMeanr2Cos3Phi;   //------<r^2*cos3phi> of wounded nucleons
   Double_t     fMeanr2Sin3Phi;   //------<r^2*sin3phi> of wounded nucleons
   Double_t     fMeanr2Cos4Phi;   //------<r^2*cos4phi> of wounded nucleons
   Double_t     fMeanr2Sin4Phi;   //------<r^2*sin4phi> of wounded nucleons
   Double_t     fMeanr2Cos5Phi;   //------<r^2*cos5phi> of wounded nucleons
   Double_t     fMeanr2Sin5Phi;   //------<r^2*sin5phi> of wounded nucleons
   Double_t     fMeanr3Cos3Phi;   //------<r^3*cos3phi> of wounded nucleons
   Double_t     fMeanr3Sin3Phi;   //------<r^3*sin3phi> of wounded nucleons
   Double_t     fMeanr4Cos4Phi;   //------<r^4*cos4phi> of wounded nucleons
   Double_t     fMeanr4Sin4Phi;   //------<r^4*sin4phi> of wounded nucleons
   Double_t     fMeanr5Cos5Phi;   //------<r^5*cos5phi> of wounded nucleons
   Double_t     fMeanr5Sin5Phi;   //------<r^5*sin5phi> of wounded nucleons
   Double_t     fMeanr2Coll;         //----------<r^2> of wounded nucleons
   Double_t     fMeanr3Coll;         //----------<r^3> of wounded nucleons
   Double_t     fMeanr4Coll;         //----------<r^4> of wounded nucleons
   Double_t     fMeanr5Coll;         //----------<r^5> of wounded nucleons
   Double_t     fMeanr2Cos2PhiColl;   //------<r^2*cos2phi> 
   Double_t     fMeanr2Sin2PhiColl;   //------<r^2*sin2phi> 
   Double_t     fMeanr2Cos3PhiColl;   //------<r^2*cos3phi> 
   Double_t     fMeanr2Sin3PhiColl;   //------<r^2*sin3phi> 
   Double_t     fMeanr2Cos4PhiColl;   //------<r^2*cos4phi> 
   Double_t     fMeanr2Sin4PhiColl;   //------<r^2*sin4phi> 
   Double_t     fMeanr2Cos5PhiColl;   //------<r^2*cos5phi> 
   Double_t     fMeanr2Sin5PhiColl;   //------<r^2*sin5phi> 
   Double_t     fMeanr3Cos3PhiColl;   //------<r^3*cos3phi> 
   Double_t     fMeanr3Sin3PhiColl;   //------<r^3*sin3phi> 
   Double_t     fMeanr4Cos4PhiColl;   //------<r^4*cos4phi> 
   Double_t     fMeanr4Sin4PhiColl;   //------<r^4*sin4phi> 
   Double_t     fMeanr5Cos5PhiColl;   //------<r^5*cos5phi> 
   Double_t     fMeanr5Sin5PhiColl;   //------<r^5*sin5phi> 
   //Double_t     fPsi2;
   Double_t     fSx2;            //Variance of x of wounded nucleons
   Double_t     fSy2;            //Variance of y of wounded nucleons
   Double_t     fSxy;            //Covariance of x and y of wounded nucleons
   Double_t     fSx2Coll;            //Variance of x of binaruy collisions
   Double_t     fSy2Coll;            //Variance of y of binaruy collisions
   Double_t     fSxyColl;            //Covariance of x and y of binaruy collisions
   Double_t     fX;              //hard particle production fraction
   Double_t     fNpp;            //Multiplicity normalization
   Bool_t       fDoPartProd;     //=1 then particle production on
   Bool_t       CalcResults(Double_t bgen);

   ClassDef(AliGlauberMC,3)
};

#endif
