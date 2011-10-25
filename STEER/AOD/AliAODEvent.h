#ifndef AliAODEvent_H
#define AliAODEvent_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD base class
//     Author: Markus Oldenburg, CERN
//-------------------------------------------------------------------------

#include <TBuffer.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TTree.h>
#include <TNamed.h>

#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliAODHeader.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODTracklets.h"
#include "AliAODJet.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODPmdCluster.h"
#include "AliAODFmdCluster.h"
#include "AliAODDimuon.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#ifdef MFT_UPGRADE
#include "AliAODMFT.h"
#endif

class TTree;
class TFolder;
class AliCentrality;
class AliEventplane;

class AliAODEvent : public AliVEvent {

 public :
  enum AODListIndex_t {kAODHeader,
		       kAODTracks,
		       kAODVertices,
		       kAODv0,
		       kAODcascade,
		       kAODTracklets,
		       kAODJets,
		       kAODEmcalCells,
		       kAODPhosCells,
		       kAODCaloClusters,
		       kAODFmdClusters,
		       kAODPmdClusters,
		       kAODDimuons,
		       kAODVZERO,
		       kAODZDC,
		       kAODListN
	           #ifdef MFT_UPGRADE
	           ,kAODVZERO
			   #endif
  };

  AliAODEvent();
  virtual ~AliAODEvent();

  AliAODEvent(const AliAODEvent& aodevent);
  AliAODEvent& operator=(const AliAODEvent& aodevent);

  void          AddObject(TObject *obj);
  void          RemoveObject(TObject *obj);
  TObject      *FindListObject(const char *objName) const;
  TList        *GetList()                const { return fAODObjects; }
  void          SetConnected(Bool_t conn=kTRUE) {fConnected=conn;}
  Bool_t        GetConnected() const {return fConnected;}

  // -- Header
  AliAODHeader *GetHeader()              const { return fHeader; }
  void          AddHeader(const AliAODHeader* hdx)
    {
	delete fHeader; fHeader = new AliAODHeader(*hdx);
	(fAODObjects->FirstLink())->SetObject(fHeader);
    }

  // setters and getters for header information
  void     SetRunNumber(Int_t n) {if (fHeader) fHeader->SetRunNumber(n);}
  void     SetPeriodNumber(UInt_t n){if (fHeader) fHeader->SetPeriodNumber(n);}
  void     SetOrbitNumber(UInt_t n) {if (fHeader) fHeader->SetOrbitNumber(n);}
  void     SetBunchCrossNumber(UShort_t n) {if (fHeader) fHeader->SetBunchCrossNumber(n);}
  void     SetMagneticField(Double_t mf){if (fHeader) fHeader->SetMagneticField(mf);}
  void     SetMuonMagFieldScale(Double_t mf){if (fHeader) fHeader->SetMuonMagFieldScale(mf);}
  void     SetDiamond(Float_t xy[2],Float_t cov[3]){if (fHeader) fHeader->SetDiamond(xy,cov);}
  void     SetDiamondZ(Float_t z, Float_t sig2z){if (fHeader) fHeader->SetDiamondZ(z,sig2z);}
  Int_t    GetRunNumber() const {return fHeader ? fHeader->GetRunNumber() : -999;}
  UInt_t   GetPeriodNumber() const {return fHeader ? fHeader->GetPeriodNumber() : 0;}
  UInt_t   GetOrbitNumber() const {return fHeader ? fHeader->GetOrbitNumber() : 0;}
  UShort_t GetBunchCrossNumber() const {return fHeader ? fHeader->GetBunchCrossNumber() : 0;}
  Double_t GetMagneticField() const {return fHeader ? fHeader->GetMagneticField() : -999.;}
  Double_t GetMuonMagFieldScale() const {return fHeader ? fHeader->GetMuonMagFieldScale() : -999.;}
  Double_t GetDiamondX() const {return fHeader ? fHeader->GetDiamondX() : -999.;}
  Double_t GetDiamondY() const {return fHeader ? fHeader->GetDiamondY() : -999.;}
  Double_t GetDiamondZ() const {return fHeader ? fHeader->GetDiamondZ() : -999.;}
  void     GetDiamondCovXY(Float_t cov[3]) const {cov[0]=-999.; if(fHeader) fHeader->GetDiamondCovXY(cov);}
  Double_t GetSigma2DiamondX() const {return fHeader ? fHeader->GetSigma2DiamondX() : -999.;}
  Double_t GetSigma2DiamondY() const {return fHeader ? fHeader->GetSigma2DiamondY() : -999.;}
  Double_t GetSigma2DiamondZ() const {return fHeader ? fHeader->GetSigma2DiamondZ() : -999.;}
  
  void      SetEventType(UInt_t eventType){fHeader->SetEventType(eventType);}
  void      SetTriggerMask(ULong64_t n) {fHeader->SetTriggerMask(n);}
  void      SetTriggerCluster(UChar_t n) {fHeader->SetTriggerCluster(n);}

  UInt_t    GetEventType()          const { return fHeader ? fHeader->GetEventType() : 0;}
  ULong64_t GetTriggerMask()        const { return fHeader ? fHeader->GetTriggerMask() : 0;}
  UChar_t   GetTriggerCluster()     const { return fHeader ? fHeader->GetTriggerCluster() : 0;}
  TString   GetFiredTriggerClasses()const { return fHeader->GetFiredTriggerClasses();};
  Double_t  GetZDCN1Energy()        const { return fHeader ? fHeader->GetZDCN1Energy() : -999.; }
  Double_t  GetZDCP1Energy()        const { return fHeader ? fHeader->GetZDCP1Energy() : -999.; }
  Double_t  GetZDCN2Energy()        const { return fHeader ? fHeader->GetZDCN2Energy() : -999.; }
  Double_t  GetZDCP2Energy()        const { return fHeader ? fHeader->GetZDCP2Energy() : -999.; }
  Double_t  GetZDCEMEnergy(Int_t i) const { return fHeader ? fHeader->GetZDCEMEnergy(i) : -999.; }

  // -- Tracks
  TClonesArray *GetTracks()              const { return fTracks; }
  Int_t         GetNTracks()             const { return fTracks? fTracks->GetEntriesFast() : 0; }
  Int_t         GetNumberOfTracks()      const { return GetNTracks(); }
  AliAODTrack  *GetTrack(Int_t nTrack)   const { return fTracks ? (AliAODTrack*)fTracks->UncheckedAt(nTrack):0; }
  Int_t         AddTrack(const AliAODTrack* trk)
    {new((*fTracks)[fTracks->GetEntriesFast()]) AliAODTrack(*trk); return fTracks->GetEntriesFast()-1;}
  Int_t         GetMuonTracks(TRefArray *muonTracks) const;
  Int_t         GetNumberOfMuonTracks() const;

  // -- Vertex
  TClonesArray *GetVertices()            const { return fVertices; }
  Int_t         GetNumberOfVertices()    const { return fVertices?fVertices->GetEntriesFast():0; }
  AliAODVertex *GetVertex(Int_t nVertex) const { return fVertices?(AliAODVertex*)fVertices->At(nVertex):0; }
  Int_t         AddVertex(const AliAODVertex* vtx)
  {new((*fVertices)[fVertices->GetEntriesFast()]) AliAODVertex(*vtx); return fVertices->GetEntriesFast()-1;}
  
  // primary vertex
  virtual AliAODVertex *GetPrimaryVertex() const { return GetVertex(0); }
  virtual AliAODVertex *GetPrimaryVertexSPD() const;

  // -- Pileup vertices 
  Int_t         GetNumberOfPileupVerticesTracks()   const;
  Int_t         GetNumberOfPileupVerticesSPD()    const;
  virtual AliAODVertex *GetPileupVertexSPD(Int_t iV=0) const;
  virtual AliAODVertex *GetPileupVertexTracks(Int_t iV=0) const;
  virtual Bool_t  IsPileupFromSPD(Int_t minContributors=3, Double_t minZdist=0.8, Double_t nSigmaZdist=3., Double_t nSigmaDiamXY=2., Double_t nSigmaDiamZ=5.) const;
  virtual Bool_t IsPileupFromSPDInMultBins() const;


  // V0
  TClonesArray *GetV0s()                 const { return fV0s; }
  Int_t         GetNumberOfV0s()         const { return fV0s->GetEntriesFast(); }
  AliAODv0     *GetV0(Int_t nV0)         const { return (AliAODv0*)fV0s->UncheckedAt(nV0); }
  Int_t         AddV0(const AliAODv0* v0)
  {new((*fV0s)[fV0s->GetEntriesFast()]) AliAODv0(*v0); return fV0s->GetEntriesFast()-1;}

  // Cascades
  TClonesArray  *GetCascades()            const { return fCascades; }
  Int_t          GetNumberOfCascades()    const { return fCascades->GetEntriesFast(); }
  AliAODcascade *GetCascade(Int_t nCasc)  const { return (AliAODcascade*)fCascades->UncheckedAt(nCasc); }
  Int_t          AddCascade(const AliAODcascade* cascade)
  {new((*fCascades)[fCascades->GetEntriesFast()]) AliAODcascade(*cascade); return fCascades->GetEntriesFast()-1;}

  // -- EMCAL and PHOS Cluster
  TClonesArray *GetCaloClusters()          const { return fCaloClusters; }
  Int_t         GetNumberOfCaloClusters()  const { return fCaloClusters->GetEntriesFast(); }
  AliAODCaloCluster *GetCaloCluster(Int_t nCluster) const { return (AliAODCaloCluster*)fCaloClusters->UncheckedAt(nCluster); }
  Int_t         AddCaloCluster(const AliAODCaloCluster* clus)
  {new((*fCaloClusters)[fCaloClusters->GetEntriesFast()]) AliAODCaloCluster(*clus); return fCaloClusters->GetEntriesFast()-1;}

  Int_t GetEMCALClusters(TRefArray *clusters) const;
  Int_t GetPHOSClusters(TRefArray *clusters) const;


  // -- FMD Cluster
  TClonesArray *GetFmdClusters()        const { return fFmdClusters; }
  Int_t         GetNFmdClusters()       const { return fFmdClusters->GetEntriesFast(); }
  AliAODFmdCluster *GetFmdCluster(Int_t nCluster) const { return (AliAODFmdCluster*)fFmdClusters->UncheckedAt(nCluster); }
  Int_t         AddFmdCluster(const AliAODFmdCluster* clus)
  {new((*fFmdClusters)[fFmdClusters->GetEntriesFast()]) AliAODFmdCluster(*clus); return fFmdClusters->GetEntriesFast()-1;}

  // -- PMD Cluster
  TClonesArray *GetPmdClusters()        const { return fPmdClusters; }
  Int_t         GetNPmdClusters()       const { return fPmdClusters->GetEntriesFast(); }
  AliAODPmdCluster *GetPmdCluster(Int_t nCluster) const { return (AliAODPmdCluster*)fPmdClusters->UncheckedAt(nCluster); }
  Int_t         AddPmdCluster(const AliAODPmdCluster* clus)
  {new((*fPmdClusters)[fPmdClusters->GetEntriesFast()]) AliAODPmdCluster(*clus); return fPmdClusters->GetEntriesFast()-1;}

  // -- Jet
  TClonesArray *GetJets()            const { return fJets; }
  Int_t         GetNJets()           const { return fJets?fJets->GetEntriesFast():0; }
  AliAODJet    *GetJet(Int_t nJet) const { return fJets?(AliAODJet*)fJets->UncheckedAt(nJet):0; }
  Int_t         AddJet(const AliAODJet* vtx)
    {new((*fJets)[fJets->GetEntriesFast()]) AliAODJet(*vtx); return fJets->GetEntriesFast()-1;}

  // -- Tracklets
  AliAODTracklets *GetTracklets() const { return fTracklets; }

  // -- Calorimeter Cells
  AliAODCaloCells *GetEMCALCells() const { return fEmcalCells; }
  AliAODCaloCells *GetPHOSCells() const { return fPhosCells; }
  const TGeoHMatrix* GetPHOSMatrix(Int_t /*i*/) const { return NULL; }
  const TGeoHMatrix* GetEMCALMatrix(Int_t /*i*/)const { return NULL; }


  // -- Dimuons
  TClonesArray *GetDimuons()              const { return fDimuons; }
  Int_t         GetNDimuons()             const { return fDimuons->GetEntriesFast(); }
  Int_t         GetNumberOfDimuons()      const { return GetNDimuons(); }
  AliAODDimuon *GetDimuon(Int_t nDimu)    const { return (AliAODDimuon*)fDimuons->UncheckedAt(nDimu); }
  Int_t         AddDimuon(const AliAODDimuon* dimu)
    {new((*fDimuons)[fDimuons->GetEntriesFast()]) AliAODDimuon(*dimu); return fDimuons->GetEntriesFast()-1;}
  
  // -- Services
  void    CreateStdContent();
  void    SetStdNames();
  void    GetStdContent();
  void    CreateStdFolders();
  void    ResetStd(Int_t trkArrSize = 0, 
		   Int_t vtxArrSize = 0, 
		   Int_t v0ArrSize = 0, 
		   Int_t cascadeArrSize = 0,
		   Int_t jetSize = 0, 
		   Int_t caloClusSize = 0, 
		   Int_t fmdClusSize = 0, 
		   Int_t pmdClusSize = 0,
		   Int_t dimuonArrsize =0
		   );
  void    ClearStd();
  void    Reset() {ClearStd();} 
  void    ReadFromTree(TTree *tree, Option_t* opt = "");
  void    WriteToTree(TTree* tree) const {tree->Branch(fAODObjects);}

  void  Print(Option_t *option="") const;
  void  MakeEntriesReferencable();
  static void AssignIDtoCollection(const TCollection* col);
  
    //Following needed only for mixed event
  virtual Int_t        EventIndex(Int_t)       const {return 0;}
  virtual Int_t        EventIndexForCaloCluster(Int_t) const {return 0;}
  virtual Int_t        EventIndexForPHOSCell(Int_t)    const {return 0;}
  virtual Int_t        EventIndexForEMCALCell(Int_t)   const {return 0;} 
  AliCentrality*       GetCentrality() {return fHeader->GetCentralityP();} 
  AliEventplane*       GetEventplane() {return fHeader->GetEventplaneP();}

  // VZERO 
  AliAODVZERO *GetVZEROData() const { return fAODVZERO; }
  virtual const Float_t* GetVZEROEqFactors() const {return fHeader?fHeader->GetVZEROEqFactors():0x0;}
  virtual Float_t        GetVZEROEqMultiplicity(Int_t i) const;
  void           SetVZEROEqFactors(const Float_t *factors) const {
    if(fHeader && factors)
      fHeader->SetVZEROEqFactors(factors);}

  //ZDC
  AliAODZDC   *GetZDCData() const { return fAODZDC; }

#ifdef MFT_UPGRADE
  // MFT 
  AliAODMFT *GetMFTData() const { return fAODMFT; }
#endif

  private :

  TList   *fAODObjects; //  list of AODObjects
  TFolder *fAODFolder;  //  folder structure of branches
  Bool_t   fConnected;  //! flag if leaves are alreday connected 
  // standard content
  AliAODHeader    *fHeader;       //! event information
  TClonesArray    *fTracks;       //! charged tracks
  TClonesArray    *fVertices;     //! vertices
  TClonesArray    *fV0s;          //! V0s
  TClonesArray    *fCascades;     //! Cascades
  AliAODTracklets *fTracklets;    //! SPD tracklets
  TClonesArray    *fJets;         //! jets
  AliAODCaloCells *fEmcalCells;   //! EMCAL calorimenter cells
  AliAODCaloCells *fPhosCells;    //! PHOS calorimenter cells
  TClonesArray    *fCaloClusters; //! calorimeter clusters
  TClonesArray    *fFmdClusters;  //! FMDclusters
  TClonesArray    *fPmdClusters;  //! PMDclusters
  TClonesArray    *fDimuons;      //! dimuons
  AliAODVZERO     *fAODVZERO;     //! VZERO AOD
  AliAODZDC       *fAODZDC;       //! ZDC AOD
#ifdef MFT_UPGRADE
  AliAODMFT       *fAODMFT;       //! VZERO AOD
#endif
  
  static const char* fAODListName[kAODListN]; //!

  ClassDef(AliAODEvent,88);
};

#endif
