
#ifndef ALISPECTRAAODHISTOMANAGER_H
#define ALISPECTRAAODHISTOMANAGER_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraAODHistoManager
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
#include "Histograms.h" // Change this file if you want to add an histogram
#include "HistogramNames.h"

namespace AliSpectraNameSpace
{

   
   enum AODParticleSpecies_t
   {
     kSpPion,
     kSpKaon,
     kSpProton,
     kNSpecies,
     kSpUndefined,
   }; // Particle species used in plotting

   const char * kParticleSpecies[] =
   {
       "PionPlus",
       "KaonPlus",
       "ProtonPlus",
       "PionMinus",
       "KaonMinus",
       "ProtonMinus",
   };

   enum AODHistoType_t
   {
     kHistPtGenTruePrimary,
     kHistPtRecSigma,
     kHistPtRecTrue,
     kHistPtRecTruePrimary,
     kHistPtRecSigmaPrimary,
     kHistPtRecSigmaSecondaryMaterial,
     kHistPtRecSigmaSecondaryWeakDecay,
     kHistNSigTPC,
     kHistNSigTOF,
     kHistNSigTPCTOF,
     kNHistoTypes
   }; // Types of histos

   enum AODCharge_t
   {
       kChPos = 0,
       kChNeg,
       kNCharge
   };

}

using namespace AliSpectraNameSpace;

class AliSpectraAODHistoManager : public TNamed
{
public:
   AliSpectraAODHistoManager() :  TNamed(), fOutputList(0) {}
   AliSpectraAODHistoManager(const char *name);
   virtual  ~AliSpectraAODHistoManager() {}


   TH2F*   BookPtGenHistogram(const char * name);
   TH2F*   BookPtRecHistogram(const char * name);
   TH2F*   BookPIDHistogram(const char * name);
   TH2F*   BookNSigHistogram(const char * name);
   TH2F*   BookqVecHistogram(const char * name);
   TH1F*   GetPtHistogram1D(const char * name,Double_t minDCA,Double_t maxDCA);
   TH1F*   GetDCAHistogram1D(const char * name,Double_t minPt,Double_t maxPt);
   TH2*     GetHistogram(UInt_t id)      {      return (TH2*) fOutputList->At(id);   }
  //   TH1*     GetHistogram(AODHistoType_t histoType, AODParticleSpecies_t particleType, UInt_t charge);
   TH1*     GetHistogram1D(UInt_t histoType, UInt_t particleType, UInt_t charge);
   TH2*     GetHistogram2D(UInt_t histoType, UInt_t particleType, UInt_t charge);
   TH2*     GetPtHistogram(UInt_t id)    {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetPtHistogram(const char * name)   {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetPtHistogramByName(UInt_t id)     {      return (TH2*) fOutputList->FindObject(kHistName[id]); }  // Use this if you want to read a file saved with a different histo list   
   TH2*     GetPIDHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetPIDHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetPIDHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistName[id]);  }// Use this if you want to read a file saved with a different histo list
   TH2*     GetNSigHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetNSigHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetNSigHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistName[id]);  }// Use this if you want to read a file saved with a different histo list
   TH2*     GetqVecHistogram(UInt_t id)   {      return (TH2*) fOutputList->At(id);   }
   TH2*     GetqVecHistogram(const char * name)  {      return (TH2*) fOutputList->FindObject(name);   }
   TH2*     GetqVecHistogramByName(UInt_t id)    {      return (TH2*) fOutputList->FindObject(kHistName[id]);  }// Use this if you want to read a file saved with a different histo list
  
   //TH1F*   GetTH1F(UInt_t id)            {      return (TH1F*) GetPtHistogram(id);   }
   //TH2F*   GetTH2F(UInt_t id)            {      return (TH2F*) GetPIDHistogram(id);   }

  TList * GetOutputList() {return fOutputList;}

  Long64_t Merge(TCollection* list);


private:
   TList     *fOutputList;  // List of Pt Histo's

   AliSpectraAODHistoManager(const AliSpectraAODHistoManager&);
   AliSpectraAODHistoManager& operator=(const AliSpectraAODHistoManager&);

   ClassDef(AliSpectraAODHistoManager, 1);

};
#endif

