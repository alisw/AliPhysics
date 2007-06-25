#ifndef ALITRDCALPIDLQ_H
#define ALITRDCALPIDLQ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Container for the distributions of dE/dx and the time bin of the          //
// max. cluster for electrons and pions                                      //
//                                                                           //
// Author:                                                                   //
//   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class TH1;
class TObjArray;
class AliTRDCalPIDLQRef;

class AliTRDCalPIDLQ : public TNamed {

 friend class AliTRDCalPIDLQRef;

 public:

  AliTRDCalPIDLQ();
  AliTRDCalPIDLQ(const Text_t *name, const Text_t *title);
  AliTRDCalPIDLQ(const AliTRDCalPIDLQ& pd);
  virtual        ~AliTRDCalPIDLQ();
  AliTRDCalPIDLQ &operator=(const AliTRDCalPIDLQ &c);

  virtual void Copy(TObject &c) const;

  Bool_t       ReadReferences(Char_t *responseFile);
  void         SetMeanChargeRatio(Double_t ratio)     { fMeanChargeRatio = ratio;  }

  Double_t     GetMeanChargeRatio() const             { return fMeanChargeRatio;   }
  Double_t     GetMomentum(Int_t ip) const            { return (ip<0 || ip>=fNMom)    ? -1. : fTrackMomentum[ip];  }
  Double_t     GetLength(Int_t il) const              { return (il<0 || il>=fNLength) ? -1. : fTrackSegLength[il]; }
  //Double_t   GetMean(Int_t iType, Int_t ip) const;
  //Double_t   GetNormalization(Int_t iType, Int_t ip) const;

  TH1*         GetHistogram(Int_t iType, Int_t ip/*, Int_t il*/) const;
  TH1*         GetHistogramT(Int_t iType, Int_t ip) const;

  Double_t     GetProbability(Int_t spec, Float_t mom, Float_t *dedx, Float_t length) const;
  Double_t     GetProbabilityT(Int_t spec, Double_t mom, Int_t timbin) const;
  //Int_t      GetNbins() const                       { return fNbins;   }
  //Double_t   GetBinSize() const                     { return fBinSize; }
  Bool_t       WriteReferences(Char_t *File="TRDPIDHistograms.root", Char_t *dir =".");

 protected:

  void         Init();
  void         CleanUp();

  inline Int_t GetHistID(Int_t part, Int_t mom/*, Int_t length=0*/) const { return part*fNMom + mom; }
  void	       SaveMaxTimeBin(const Int_t mom, const char *fn);
				
 protected:

  static  Char_t    *fpartName[5];     //! Names of particle species
  static  Char_t    *fpartSymb[5];     //! Symbols of particle species  
  Int_t              fNMom;            //  Number of momenta
  Int_t              fNLength;         //  Number of track segment length intervals

  Double_t          *fTrackMomentum;   //[fNMom]    Track momenta for which response functions are available
  Double_t          *fTrackSegLength;  //[fNLength] Track segment lengths for which response functions are available
  Int_t              fNTimeBins;       //  Number of time bins

  Double_t           fMeanChargeRatio; //  Ratio of mean charge from real Det. to prob. dist.
  Int_t              fNbins;           //  Number of energy bins
  Double_t           fBinSize;         //  Size of energy bin

  TObjArray         *fHistdEdx;        //  Prob. of dEdx for 5 particles and for several momenta
  TObjArray         *fHistTimeBin;     //  Prob. of max time bin for 5 particles and for several momenta

 private:

  TH1 	            *h1MaxTB[2];      //!
  AliTRDCalPIDLQRef *fEstimator;      //!
	
  ClassDef(AliTRDCalPIDLQ, 2)         //   The TRD PID response container class

};

#endif

