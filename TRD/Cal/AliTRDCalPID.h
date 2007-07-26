#ifndef ALITRDCALPID_H
#define ALITRDCALPID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Container for the distributions of dE/dx and the time bin of the          //
// max. cluster for electrons and pions                                      //
//                                                                           //
// Authors:                                                                  //
//   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>                    //
//   Alex Bercuci <A.Bercuci@gsi.de>                                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TH1;
class TObjArray;
class AliTRDCalPID : public TNamed {
public:
	enum {
		kNMom = 11,
		kNLength = 4
	};

  AliTRDCalPID();
  AliTRDCalPID(const Text_t *name, const Text_t *title);
  AliTRDCalPID(const AliTRDCalPID& pd);
  virtual        ~AliTRDCalPID();
  AliTRDCalPID&   operator=(const AliTRDCalPID &c);
  virtual void    Copy(TObject &c) const;
 
           Bool_t   LoadLQReferences(Char_t* refFile);
           Bool_t   LoadNNReferences(Char_t* /*refFile*/) {return kTRUE;}
  //void         SetMeanChargeRatio(Double_t ratio)     { fMeanChargeRatio = ratio;  }

  //Double_t     GetMeanChargeRatio() const             { return fMeanChargeRatio;   }
  static   Double_t GetMomentum(Int_t ip)            { return (ip<0 || ip>=kNMom)    ? -1. : fTrackMomentum[ip];  }
  static   Double_t GetLength(Int_t il)              { return (il<0 || il>=kNLength) ? -1. : fTrackSegLength[il]; }
  //Double_t   GetMean(Int_t iType, Int_t ip) const;
  //Double_t   GetNormalization(Int_t iType, Int_t ip) const;

           TH1*     GetHistogram(Int_t iType, Int_t ip/*, Int_t il*/) const;
           TH1*     GetHistogramT(Int_t iType, Int_t ip) const;
           Double_t GetProbability(Int_t spec, Float_t mom, Float_t *dedx, Float_t length) const;
           Double_t GetProbabilityT(Int_t spec, Double_t mom, Int_t timbin) const;

 protected:

           void     Init();
  inline  Int_t     GetHistID(Int_t part, Int_t mom/*, Int_t length=0*/) const { return part*kNMom + mom; }
           void     CleanUp();
 
 public:
  static  Char_t   *fpartName[5];     //! Names of particle species
  static  Char_t   *fpartSymb[5];     //! Symbols of particle species
  
 protected:
	static Float_t   fTrackMomentum[kNMom];     // Track momenta for which response functions are available
        static Float_t   fTrackSegLength[kNLength]; // Track segment lengths for which response functions are available
         Double_t  fMeanChargeRatio;  //  Ratio of mean charge from real Det. to prob. dist.
         TObjArray *fHistdEdx;        //  Prob. of dEdx for 5 particles and for several momenta
         TObjArray *fHistTimeBin;     //  Prob. of max time bin for 5 particles and for several momenta

	
  ClassDef(AliTRDCalPID, 1)           //   The TRD PID response container class

};

#endif

