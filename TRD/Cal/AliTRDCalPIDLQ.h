#ifndef ALITRDCALPIDLQ_H
#define ALITRDCALPIDLQ_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*-----------------------------------------------------------------
   Class for dE/dx and Time Bin of Max. Cluster for Electrons and 
   pions in TRD. 
   It is instantiated in class AliTRDpidESD for particle identification
   in TRD
   Prashant Shukla <shukla@pi0.physi.uni-heidelberg.de>
   -----------------------------------------------------------------*/

#include <TNamed.h>
#include <AliPID.h>

class TH1F;
class TObjArray;

class AliTRDCalPIDLQ : public TNamed {

  public:

    AliTRDCalPIDLQ(); 
    AliTRDCalPIDLQ(const Text_t *name, const Text_t *title);
    AliTRDCalPIDLQ(const AliTRDCalPIDLQ& pd);
    virtual ~AliTRDCalPIDLQ();              

    AliTRDCalPIDLQ &operator=(const AliTRDCalPIDLQ &c);
    virtual void      Copy(TObject &c) const;

            Bool_t    ReadData(Char_t *responseFile);     
            void      SetMeanChargeRatio(Double_t ratio)    { fMeanChargeRatio = ratio; }  

            Double_t  GetMeanChargeRatio() const            { return fMeanChargeRatio;  } 
            Double_t  GetMomentum(Int_t ip) const           { return fTrackMomentum[ip];}
            Double_t  GetMean(Int_t iType, Int_t ip) const;        
            Double_t  GetNormalization(Int_t iType, Int_t ip) const;

            TH1F     *GetHistogram(Int_t iType, Int_t ip) const;
            TH1F     *GetHistogramT(Int_t iType, Int_t ip) const;

            Double_t  GetProbability(Int_t iType, Double_t mom, Double_t dedx) const;
            Double_t  GetProbabilityT(Int_t iType, Double_t mom, Int_t timbin) const;
            Int_t     GetNbins() const                      { return fNbins;   }        
            Double_t  GetBinSize() const                    { return fBinSize; } 

  protected:

            void  Init();      
            void  CleanUp();   
    inline  Int_t GetHistID(Int_t part, Int_t mom) const    { return part*fNMom + mom; }
    
    static Char_t *fpartName[AliPID::kSPECIES]; //! Names of particle species
    
    Int_t      fNMom;            // Number of momenta  
    Double_t  *fTrackMomentum;   //[fNMom] Track momenta for which response functions are available
    Double_t   fMeanChargeRatio; // Ratio of mean charge from real Det. to prob. dist.

    Int_t      fNbins;           // Number of energy bins
    Double_t   fBinSize;         // Size of energy bin
    
    TObjArray *fHistdEdx;        // Prob. of dEdx for 5 particles and for several momenta
    TObjArray *fHistTimeBin;     // Prob. of max time bin for 5 particles and for several momenta
    
    ClassDef(AliTRDCalPIDLQ, 1)  // The TRD PID response container class

};

#endif

