#ifndef ALISPECTRAAODEVENTCUTS_H
#define ALISPECTRAAODEVENTCUTS_H

/*  See cxx source for full Copyright notice */

//-------------------------------------------------------------------------
//                      AliSpectraAODEventCuts
//
//
//
//
// Authors: Michele Floris, CERN, Philip Versteeg, UU, Redmer Bertens, UU
//-------------------------------------------------------------------------

class AliAODEvent;
class AliSpectraAODHistoManager;

#include "TNamed.h"

class AliSpectraAODEventCuts : public TNamed
{
public:
   enum { kAcceptedEvents = 0, kVtxRange, kVtxCentral, kVtxNoEvent, kNVtxCuts};

   // Constructors
   AliSpectraAODEventCuts() : TNamed(), fAOD(0), fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fHistoCuts(0) {}
   AliSpectraAODEventCuts(const char *name);
   virtual  ~AliSpectraAODEventCuts() {}

   // Methods
   Bool_t IsSelected(AliAODEvent * aod);
   Bool_t CheckVtxRange();
   Bool_t CheckCentralityCut();
   void  SetCentralityCutMin(Float_t cut)  { fCentralityCutMin = cut; }
   void  SetCentralityCutMax(Float_t cut)  { fCentralityCutMax = cut; }

   
   TH1I * GetHistoCuts()         {  return fHistoCuts; }
   Float_t  GetCentralityMin()  const {  return fCentralityCutMin; }
   Float_t  GetCentralityMax()  const {  return fCentralityCutMax; }
   void   PrintCuts();
   Float_t  NumberOfEvents()     { return fHistoCuts->GetBinContent(kAcceptedEvents+1); }

  Long64_t Merge(TCollection* list);


private:

   AliAODEvent     *fAOD;              //! AOD event
   Bool_t          fIsSelected;        // True if cuts are selected
   Float_t         fCentralityCutMin;     // minimum centrality percentile
  Float_t         fCentralityCutMax;     // maximum centrality percentile
   TH1I            *fHistoCuts;        // Cuts statistics
   AliSpectraAODEventCuts(const AliSpectraAODEventCuts&);
   AliSpectraAODEventCuts& operator=(const AliSpectraAODEventCuts&);

   ClassDef(AliSpectraAODEventCuts, 2);

};
#endif

