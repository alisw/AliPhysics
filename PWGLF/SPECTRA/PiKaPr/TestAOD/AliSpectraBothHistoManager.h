#ifndef ALISPECTRABOTHHISTOMANAGER_H
#define ALISPECTRABOTHHISTOMANAGER_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraBothHistoManager
//
//
//
//
// Authors: Michele Floris, CERN, Philip Versteeg, UU, Redmer Bertens, UU
//-------------------------------------------------------------------------

class AliAODEvent;
class TH1F;
class TH2F;
class TH1;
class TH2;

//class TList;
#include "TNamed.h"
#include "TList.h"
#include "HistogramsBoth.h" // Change this file if you want to add an histogram
#include "HistogramNamesBoth.h" // generate this automatically running createNames.py 

namespace AliSpectraNameSpaceBoth
{

   
   enum BothParticleSpecies_t
   {
     kSpPion,
     kSpKaon,
     kSpProton,
     kNSpecies,
     kSpUndefined,
   }; // Particle species used in plotting

   //extern const char * kParticleSpecies[];

   enum BothHistoType_t
   {
     kHistPtGenTruePrimary,
     kHistPtRecSigma,
     kHistPtRecTrue,
     kHistPtRecTruePrimary,
     kHistPtRecPrimary,	
     kHistPtRecSigmaPrimary,
     kHistPtRecSigmaSecondaryMaterial,
     kHistPtRecSigmaSecondaryWeakDecay,
     kHistNSigTPC,
     kHistNSigTOF,
     kHistNSigTPCTOF,
     kNHistoTypes
   }; // Types of histos

   enum BothCharge_t
   {
       kChPos = 0,
       kChNeg,
       kNCharge
   };

}

using namespace AliSpectraNameSpaceBoth;

class AliSpectraBothHistoManager : public TNamed
{
public:
   AliSpectraBothHistoManager() :  TNamed(), fOutputList(0), fNRebin(0) ,fIncludecorrectlyidentifiedinMCtemplates(kFALSE){}
  AliSpectraBothHistoManager(const char *name,Int_t nrebin,Bool_t pidqa=kTRUE);
   virtual  ~AliSpectraBothHistoManager() ;


   TH2F*   BookPtGenHistogram(const char * name);
   TH2F*   BookPtGenAllChHistogram(const char * name);
   TH2F*   BookPtRecHistogram(const char * name);
   TH2F*   BookPtRecAllChHistogram(const char * name);
   TH2F*   BookPIDHistogram(const char * name);
   TH2F*   BookNSigHistogram(const char * name);
   TH2F*   BookGenMulvsRawMulHistogram(const char * name);
   TH2F*   BookDoubleCountsHistogram(const char * name);
   TH1F*   BookEventStatHist(); // due to zip bug in merging it is good to store in one list everything
   TH2F*   BookNSigTOFmissmatchHistogram(const char * name);	 

   TH1F*   GetPtHistogram1D(const char * name,Double_t minDCA,Double_t maxDCA);
   TH1F*   GetDCAHistogram1D(const char * name,Double_t minPt,Double_t maxPt);
   TH2*     GetHistogram(UInt_t id)      {      return (TH2*) fOutputList->At(id);   }
   //   TH1*     GetHistogram(BothHistoType_t histoType, BothParticleSpecies_t particleType, UInt_t charge);
   TH1*     GetHistogram1D(UInt_t histoType, UInt_t particleType, UInt_t charge);
   TH2*     GetHistogram2D(UInt_t histoType, UInt_t particleType, UInt_t charge);
   TH2*     GetPtHistogram(UInt_t id)    {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetPtHistogram(const char * name)   {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetPtHistogramByName(UInt_t id)     {      return (TH2*) fOutputList->FindObject(kHistNameBoth[id]); }  // Use this if you want to read a file saved with a different histo list   
   TH2*     GetPIDHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetPIDHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetPIDHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistNameBoth[id]);  }// Use this if you want to read a file saved with a different histo list
   TH2*     GetNSigHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetNSigHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetNSigHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistNameBoth[id]);  }// Use this if you want to read a file saved with a different histo list
   TH2*     GetqVecHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetqVecHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetqVecHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistNameBoth[id]);  }// Use this if you want to read a file saved with a different histo list
   TH2*     GetGenMulvsRawMulHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetGenMulvsRawMulHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistNameBoth[id]);  }// Use this if you want to read a file saved with a different histo list
   TH1F*    GetEventStatHist()  {      return (TH1F*) fOutputList->FindObject("EventStatHisto");  }	
   //TH1F*   GetTH1F(UInt_t id)            {      return (TH1F*) GetPtHistogram(id);   }
   //TH2F*   GetTH2F(UInt_t id)            {      return (TH2F*) GetPIDHistogram(id);   }

  TList * GetOutputList() {return fOutputList;}
  void    SetNRebin(Int_t nreb){fNRebin=nreb;}
  Int_t   GetNRebin() {return fNRebin;}
  void SetIncludecorrectlyidentifiedinMCtemplates(Bool_t flag=kFALSE){fIncludecorrectlyidentifiedinMCtemplates=flag;}
  Bool_t GetIncludecorrectlyidentifiedinMCtemplates() {return fIncludecorrectlyidentifiedinMCtemplates;}

  Long64_t Merge(TCollection* list);


private:
   TList     *fOutputList;  // List of Pt Histo's
   Int_t      fNRebin; //rebin of histos
   Bool_t      fIncludecorrectlyidentifiedinMCtemplates; // if set to true secondary templates are only filed after checking MC PID	
   AliSpectraBothHistoManager(const AliSpectraBothHistoManager&);
   AliSpectraBothHistoManager& operator=(const AliSpectraBothHistoManager&);

   ClassDef(AliSpectraBothHistoManager, 3);

};
#endif

