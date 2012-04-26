
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
class TList;
#include "TNamed.h"
#include "Histograms.h" // Change this file if you want to add an histogram
#include "HistogramNames.h"

namespace AliSpectraNameSpace
{

   
   enum AODParticleSpecies_t
   {
     kProton = 0,
     kKaon,
     kPion,
     kNSpecies = kPion,
   }; // Particle species used in plotting

   const char * kParticleSpecies[] =
   {
       "ProtonPlus",
       "ProtonMinus",
       "KaonPlus",
       "KaonMinus",
       "PionPlus",
       "PionMinus",
   };
   enum AODHistoType_t
   {
       kHistSpectraRec = 0,
       kHistSpectraGen,
       kHNHistoTypes = kHistSpectraGen,
   }; // Types of histos

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

