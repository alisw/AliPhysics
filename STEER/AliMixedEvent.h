// -*- mode: C++ -*- 
#ifndef ALIMIXEDEVENT_H
#define ALIMIXEDEVENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliMixedEvent
// VEvent which is the container of several VEvents 
// Use Case: Event Mixing     
// Origin: Andreas Morsch, CERN, Andreas.Morsch@cern.ch 
//-------------------------------------------------------------------------

#include "AliVEvent.h"
#include <TList.h>
#include <TObjArray.h>
class AliCentrality;
class AliEventplane;
class AliVVertex;

class AliMixedEvent : public AliVEvent {

public:
    
    AliMixedEvent();
    virtual ~AliMixedEvent() ;
    AliMixedEvent(const AliMixedEvent& vEvnt); 
    AliMixedEvent& operator=(const AliMixedEvent& vEvnt);
    
    // Services from VEvent
    virtual void AddObject(TObject* /*obj*/) {;}
    virtual TObject* FindListObject(const char* /*name*/) const {return 0;}
    virtual TList* GetList() const {return 0;}
    virtual void CreateStdContent() {;}
    virtual void GetStdContent() {;}
    virtual void ReadFromTree(TTree* /*tree*/, Option_t* /*opt*/) {;} 
    virtual void WriteToTree(TTree* /*tree*/) const {;}
    virtual void SetStdNames() {;} 
    virtual void Print(Option_t * /*option*/) const  {;} 
    // Specific Services
    virtual void AddEvent(AliVEvent* evt);
    virtual void SetPrimaryVertex(AliVVertex* newVertex) {fMeanVertex = newVertex;}
    virtual void Reset();
    virtual void Init();
    const AliVVertex* GetEventVertex(Int_t i) const;
    const AliVVertex* GetVertexOfEvent(Int_t i) const;
    Int_t        GetNumberOfEvents() {return fNEvents;}      
  
    // Header
    virtual AliVHeader* GetHeader() const {return 0;}

    // Delegated methods for fESDRun or AODHeader
  
    virtual void     SetRunNumber(Int_t /*n*/)         {;} 
    virtual void     SetPeriodNumber(UInt_t /*n*/)     {;} 
    virtual void     SetMagneticField(Double_t /*mf*/) {;} 
    
    virtual Int_t    GetRunNumber()     const  {return -999 ;} 
    virtual UInt_t   GetPeriodNumber()  const  {return    0 ;} 
    virtual Double_t GetMagneticField() const;
    
    
    virtual Double_t GetDiamondX() const {return -999.;}
    virtual Double_t GetDiamondY() const {return -999.;}
    virtual void     GetDiamondCovXY(Float_t cov[3]) const
	{cov[0]=-999.; return;}
  
  // Delegated methods for fHeader
    virtual void      SetOrbitNumber(UInt_t /*n*/)        {;} 
    virtual void      SetBunchCrossNumber(UShort_t /*n*/) {;}
    virtual void      SetEventType(UInt_t /*eventType*/)  {;}
    virtual void      SetTriggerMask(ULong64_t /*n*/)     {;}
    virtual void      SetTriggerCluster(UChar_t /*n*/)    {;}

    virtual UInt_t    GetOrbitNumber()      const  {return    0;}
    virtual UShort_t  GetBunchCrossNumber() const  {return    0;}
    virtual UInt_t    GetEventType()        const  {return    0;}
    virtual ULong64_t GetTriggerMask()      const  {return    0;}
    virtual UChar_t   GetTriggerCluster()   const  {return    0;}
    
    virtual Double_t  GetZDCN1Energy()            const {return -999.;}
    virtual Double_t  GetZDCP1Energy()            const {return -999.;}
    virtual Double_t  GetZDCN2Energy()            const {return -999.;}
    virtual Double_t  GetZDCP2Energy()            const {return -999.;}
    virtual Double_t  GetZDCEMEnergy(Int_t /*i*/) const {return -999.;}
 
  // Tracks
    virtual AliVParticle *GetTrack(Int_t i)    const;
    virtual Int_t        GetNumberOfTracks()   const {return fNumberOfTracks;}
    virtual Int_t        GetNumberOfV0s()      const {return -999;}
    virtual Int_t        GetNumberOfCascades() const {return -999;}
  
    // Calo Clusters and cells
  virtual AliVCluster *GetCaloCluster(Int_t i)     const;
  virtual Int_t        GetNumberOfCaloClusters()   const {return fNumberOfCaloClusters;}
  virtual AliVCaloCells *GetPHOSCells()            const {return fPHOSCells;}
  virtual AliVCaloCells *GetEMCALCells()           const {return fEMCALCells;}
  virtual Int_t        GetNumberOfPHOSCells()      const {return fNumberOfPHOSCells;}
  virtual Int_t        GetNumberOfEMCALCells()     const {return fNumberOfEMCALCells;}
  virtual Int_t*       GetPHOSCellsCumul()         const {return fNPHOSCellsCumul;} 
  virtual Int_t*       GetEMCALCellsCumul()        const {return fNEMCALCellsCumul;} 
  virtual Int_t        GetNCaloClustersCumul(Int_t iev) const {return fNCaloClustersCumul[iev];}


  virtual Int_t        EventIndex(Int_t itrack) const;
  virtual Int_t        EventIndexForCaloCluster(Int_t iclu) const;
  virtual Int_t        EventIndexForPHOSCell(Int_t icell) const;
  virtual Int_t        EventIndexForEMCALCell(Int_t icell) const;

  virtual AliCentrality* GetCentrality() {return 0;}
  virtual AliEventplane* GetEventplane() {return 0;}
  // Primary vertex
    virtual const AliVVertex   *GetPrimaryVertex() const {return fMeanVertex;}
    virtual Bool_t ComputeVtx(const TObjArray *vertices, Double_t *pos,Double_t *sig,Int_t *nContributors);
  // VZERO
  virtual AliVVZERO *GetVZEROData() const {return 0;}
private:
  TList   fEventList;            //! List of Events
  Int_t   fNEvents;              //! Number of Events 
  Int_t   fNumberOfTracks;       //! Total number of tracks
  Int_t   fNumberOfCaloClusters; //! Total number of calo clusters
  Int_t   fNumberOfPHOSCells;    //! Total number of PHOS Cells
  Int_t   fNumberOfEMCALCells;   //! Total number of EMCAL Cells
  Int_t*  fNTracksCumul;         //! Cumulant
  Int_t*  fNCaloClustersCumul;   //! Cumulant
  Int_t*  fNPHOSCellsCumul;      //! Cumulant
  Int_t*  fNEMCALCellsCumul;     //! Cumulant
  AliVCaloCells* fPHOSCells;     //! array ofPHOS cells  
  AliVCaloCells* fEMCALCells;    //! array ofPHOS cells  
  AliVVertex* fMeanVertex;    //! Mean vertex
    
    ClassDef(AliMixedEvent, 0)  // Container for mixed events
};
#endif 

