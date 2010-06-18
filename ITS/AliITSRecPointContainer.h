#ifndef ALIITSRECPOINTCONTAINER_H
#define ALIITSRECPOINTCONTAINER_H
/* Copyright(c) 2009-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//           Container class for ITS rec points                       //
////////////////////////////////////////////////////////////////////////

#include <TClonesArray.h>
#include <TString.h>

class AliITSRecoParam;

class AliITSRecPointContainer : public TObject {

 public:

  virtual ~AliITSRecPointContainer();  //Destructor

  Bool_t IsSPDActive() const {return fDet.Contains("SPD");}
  Bool_t IsSDDActive() const {return fDet.Contains("SDD");}
  Bool_t IsSSDActive() const {return fDet.Contains("SSD");}
  Bool_t IsITSComplete() const {return fDet.Contains("ALL");}
  Bool_t GetStatusOK() const {return fStatusOK;}
  Int_t GetNumberOfModules() const {return fActualSize; }

  static AliITSRecPointContainer* Instance(const AliITSRecoParam *ptr=NULL);
  void PrepareToRead(){if(fNextEvent<0){fNextEvent=0;} else {++fNextEvent;}}
  TClonesArray* FetchClusters(Int_t mod, TTree* tR);
  TClonesArray* FetchClusters(Int_t mod, TTree* tR,Int_t cureve);
  TClonesArray* UncheckedGetClusters(Int_t mod) const {return fArray[mod];}

  // In the following two methods: 1<=lay<=6  (i.e. layers numbered from 1) 
  UInt_t GetNClustersInLayer(Int_t lay, TTree* tR, Int_t eventN=-1);
  UInt_t GetNClustersInLayerFast(Int_t lay) const;
  void FullReset(){fCurrentEve=-1000; Reset();}
  void ResetSPD(); // clears only SPD parts - see implementation for usage
  void ResetSDD(); // clears only SPD parts - see implementation for usage
  void ResetSSD(); // clears only SPD parts - see implementation for usage

 private:
  // methods
  AliITSRecPointContainer(const AliITSRecoParam* krp=NULL);   // Default constructor
  AliITSRecPointContainer(const AliITSRecPointContainer& rec);
  AliITSRecPointContainer& operator=(const AliITSRecPointContainer &source);
  Bool_t CheckBoundaries(Int_t i)const { return (i>=0 && i<fgkNModules);}
  void CookEntries();
  void Reset();
  void ClearClus(Int_t first, Int_t lastpp){ // clears clusters for modules 
    // ranging from first to lastpp-1 included
    for(Int_t i=first;i<lastpp;i++)(fArray[i])->Clear();
  }
  //Data members
  static AliITSRecPointContainer* fgInstance; //! AliITSRecPointContainer 
                                              //  singleton
  static const Int_t fgkNModules=2198;  //! total number of ITS modules

  Int_t fSPDNModules; //! number of SPD modules
  Int_t fSDDNModules; //! number of SDD modules
  Int_t fSSDNModules; //! number of SDD modules
  TClonesArray* fArray[fgkNModules];  //! container - 1 TClonesArray per module
  Int_t fCurrentEve; //! event number
  Int_t fNextEvent; //! number of the next event to be read; used only when
                    //! the run loader is not available. It is just a counter.
  Int_t fActualSize; //! actual number of ITS modules in TTree R 
  TString fDet; //! ITS subdetectors active for the current run 
  Bool_t fStatusOK; //! kFALSE is RP branch is absent or if there are anomalies
                    //! in the number of active modules
  UInt_t fNClusters[6]; //! Number of clusters per layer

  ClassDef(AliITSRecPointContainer,0)
};

#endif
