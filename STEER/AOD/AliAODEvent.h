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
#include "AliVHeader.h"
#include "AliAODHeader.h"
#include "AliNanoAODHeader.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODTracklets.h"
#include "AliAODJet.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloTrigger.h"
#include "AliAODPmdCluster.h"
#include "AliAODFmdCluster.h"
#include "AliAODDimuon.h"
#include "AliAODTZERO.h"
#include "AliAODVZERO.h"
#include "AliAODHMPIDrings.h"
#include "AliAODZDC.h"
#include "AliAODAD.h"
#include "AliAODTrdTrack.h"

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
	           kAODEMCALTrigger,
	           kAODPHOSTrigger,
		       kAODFmdClusters,
		       kAODPmdClusters,
                       kAODHMPIDrings,
		       kAODDimuons,
		       kAODTZERO,
		       kAODVZERO,
		       kAODZDC,
		       kAODAD,
		       kTOFHeader,                       
		       kAODTrdTracks,
		       kAODListN
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
  Bool_t        AreTracksConnected() const {return fTracksConnected;}

  // -- Header
  AliVHeader    *GetHeader()              const { return fHeader; }
  void          AddHeader(const AliVHeader* hdx)
    {
      delete fHeader; 
      if(dynamic_cast<const AliAODHeader*>(hdx)) {
	fHeader = new AliAODHeader(*(const AliAODHeader*)hdx);
      } else if (dynamic_cast<const AliNanoAODHeader*>(hdx)) {
        fHeader = new AliNanoAODHeader(*(const AliNanoAODHeader*)hdx);
      }
      else {
        AliError(Form("Unknown header type %s", hdx->ClassName()));
      }
        (fAODObjects->FirstLink())->SetObject(fHeader);
    }

  virtual  Bool_t InitMagneticField() const {return fHeader ? fHeader->InitMagneticField() : kFALSE;}

  void   SetDAQAttributes(UInt_t attributes) {if (fHeader) fHeader->SetDAQAttributes(attributes);}
  UInt_t GetDAQAttributes() const            {return fHeader ? fHeader->GetDAQAttributes() : 0;}
  Bool_t IsIncompleteDAQ();


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
  Int_t     GetNumberOfESDTracks()  const { return fHeader ? fHeader->GetNumberOfESDTracks() : 0; }
  Int_t     GetNumberOfITSClusters(Int_t lr) const {return fHeader ? (int)fHeader->GetNumberOfITSClusters(lr) : 0;}
  void SetTOFHeader(const AliTOFHeader * tofEventTime);
  const AliTOFHeader *GetTOFHeader() const {return fTOFHeader;}
  Float_t GetEventTimeSpread() const {if (fTOFHeader) return fTOFHeader->GetT0spread(); else return 0.;}
  Float_t GetTOFTimeResolution() const {if (fTOFHeader) return fTOFHeader->GetTOFResolution(); else return 0.;}
  Float_t GetT0spread(Int_t i) const {return fHeader->GetT0spread(i);}


  // -- Tracks
  TClonesArray *GetTracks()              const { return fTracks; }
  void          ConnectTracks();
  Int_t         GetNumberOfTracks()      const { return fTracks? fTracks->GetEntriesFast() : 0; }
  AliVTrack    *GetTrack(Int_t nTrack)   const { return fTracks ? (AliVTrack*)fTracks->UncheckedAt(nTrack):0; }
  Int_t         AddTrack(const AliAODTrack* trk);
  Int_t         GetMuonTracks(TRefArray *muonTracks) const;
  Int_t         GetNumberOfMuonTracks() const;
  Int_t         GetMuonGlobalTracks(TRefArray *muonGlobalTracks) const;    // AU
  Int_t         GetNumberOfMuonGlobalTracks() const;                       // AU

  // -- Vertex
  TClonesArray *GetVertices()            const { return fVertices; }
  Int_t         GetNumberOfVertices()    const { return fVertices?fVertices->GetEntriesFast():0; }
  AliAODVertex *GetVertex(Int_t nVertex) const { return fVertices?(AliAODVertex*)fVertices->At(nVertex):0; }
  Int_t         AddVertex(const AliAODVertex* vtx)
  {new((*fVertices)[fVertices->GetEntriesFast()]) AliAODVertex(*vtx); return fVertices->GetEntriesFast()-1;}
  
  // primary vertex
  using AliVEvent::GetPrimaryVertex;
  using AliVEvent::GetPrimaryVertexSPD;
  using AliVEvent::GetPrimaryVertexTPC;
  virtual AliAODVertex *GetPrimaryVertex() const { return GetVertex(0); }
  virtual AliAODVertex *GetPrimaryVertexSPD() const;
  virtual AliAODVertex *GetVertex() const { return GetPrimaryVertexSPD(); }
  virtual AliAODVertex *GetPrimaryVertexTPC() const;

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
  using AliVEvent::GetV0;
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
  Int_t         GetNumberOfCaloClusters()  const { return fCaloClusters?fCaloClusters->GetEntriesFast():0; }
  AliAODCaloCluster *GetCaloCluster(Int_t nCluster) const { return fCaloClusters?(AliAODCaloCluster*)fCaloClusters->UncheckedAt(nCluster):0x0; }
  Int_t         AddCaloCluster(const AliAODCaloCluster* clus)
  {new((*fCaloClusters)[fCaloClusters->GetEntriesFast()]) AliAODCaloCluster(*clus); return fCaloClusters->GetEntriesFast()-1;}
  AliAODCaloTrigger *GetCaloTrigger(TString calo) const 
  { 	  
     if (calo.Contains("EMCAL")) return fEMCALTrigger;
     else
     return fPHOSTrigger; 
  }	
	
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

  // -- HMPID objects 
  TClonesArray *GetHMPIDrings()       const {return fHMPIDrings; } 
  Int_t         GetNHMPIDrings()      const;
  AliAODHMPIDrings *GetHMPIDring(Int_t nRings) const;
  Int_t         AddHMPIDrings(const  AliAODHMPIDrings* ring) 
  {new((*fHMPIDrings)[fHMPIDrings->GetEntriesFast()]) AliAODHMPIDrings(*ring); return fHMPIDrings->GetEntriesFast()-1;}
  
  AliAODHMPIDrings *GetHMPIDringForTrackID(Int_t trackID) const;
  
  
  // -- Jet
  TClonesArray *GetJets()            const { return fJets; }
  Int_t         GetNJets()           const { return fJets?fJets->GetEntriesFast():0; }
  AliAODJet    *GetJet(Int_t nJet) const { return fJets?(AliAODJet*)fJets->UncheckedAt(nJet):0; }
  Int_t         AddJet(const AliAODJet* vtx)
    {new((*fJets)[fJets->GetEntriesFast()]) AliAODJet(*vtx); return fJets->GetEntriesFast()-1;}

  // -- Tracklets
  AliAODTracklets *GetTracklets() const { return fTracklets; }  
  AliAODTracklets *GetMultiplicity() const {return GetTracklets();}
  // -- Calorimeter Cells
  AliAODCaloCells *GetEMCALCells() const { return fEmcalCells; }
  AliAODCaloCells *GetPHOSCells() const { return fPhosCells; }
  const TGeoHMatrix* GetPHOSMatrix(Int_t /*i*/) const { return NULL; }
  const TGeoHMatrix* GetEMCALMatrix(Int_t /*i*/)const { return NULL; }


  // -- Dimuons (\deprecated)
  TClonesArray *GetDimuons() const;
  Int_t         GetNDimuons() const;
  Int_t         GetNumberOfDimuons() const;
  AliAODDimuon *GetDimuon(Int_t nDimu) const;
  Int_t         AddDimuon(const AliAODDimuon* dimu);
  
  // // -- TRD
  Int_t GetNumberOfTrdTracks() const { return fTrdTracks ? fTrdTracks->GetEntriesFast() : 0; }
  AliAODTrdTrack* GetTrdTrack(Int_t i) const {
    return (AliAODTrdTrack *) (fTrdTracks ? fTrdTracks->At(i) : 0x0);
  }
  AliAODTrdTrack& AddTrdTrack(const AliVTrdTrack *track);

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
                   Int_t hmpidRingsSize = 0,
		   Int_t dimuonArrsize =0,
		   Int_t nTrdTracks = 0
		   );
  void    ClearStd();
  void    Reset(); 
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

  // TZERO 
  AliAODTZERO *GetTZEROData() const { return fAODTZERO; }
  Double32_t GetT0TOF(Int_t icase) const { return fAODTZERO?fAODTZERO->GetT0TOF(icase):999999;}
  const Double32_t * GetT0TOF() const { return fAODTZERO?fAODTZERO->GetT0TOF():0x0;}
 
  // VZERO 
  AliAODVZERO *GetVZEROData() const { return fAODVZERO; }
  virtual const Float_t* GetVZEROEqFactors() const {return fHeader?fHeader->GetVZEROEqFactors():0x0;}
  virtual Float_t        GetVZEROEqMultiplicity(Int_t i) const;
  virtual void   SetVZEROEqFactors(Float_t factors[64]) const {
    if(fHeader)
      fHeader->SetVZEROEqFactors(factors);}

  //ZDC
  AliAODZDC   *GetZDCData() const { return fAODZDC; }
  
  //AD
  AliAODAD   *GetADData() const { return fAODAD; }

  virtual AliVEvent::EDataLayoutType GetDataLayoutType() const;

  private :

  TList   *fAODObjects; //  list of AODObjects
  TFolder *fAODFolder;  //  folder structure of branches
  Bool_t   fConnected;  //! flag if leaves are alreday connected 
  Bool_t   fTracksConnected;      //! flag if tracks have already pointer to event set
  // standard content
  AliVAODHeader   *fHeader;       //! event information
  TClonesArray    *fTracks;       //! charged tracks
  TClonesArray    *fVertices;     //! vertices
  TClonesArray    *fV0s;          //! V0s
  TClonesArray    *fCascades;     //! Cascades
  AliAODTracklets *fTracklets;    //! SPD tracklets
  TClonesArray    *fJets;         //! jets
  AliAODCaloCells *fEmcalCells;   //! EMCAL calorimenter cells
  AliAODCaloCells *fPhosCells;    //! PHOS calorimenter cells
  TClonesArray    *fCaloClusters; //! calorimeter clusters
  AliAODCaloTrigger *fEMCALTrigger; //! EMCAL Trigger information
  AliAODCaloTrigger *fPHOSTrigger;  //! PHOS Trigger information
  TClonesArray    *fFmdClusters;  //! FMDclusters
  TClonesArray    *fPmdClusters;  //! PMDclusters
  TClonesArray    *fHMPIDrings;   //! HMPID signals
  TClonesArray    *fDimuons;      //! dimuons
  AliAODTZERO     *fAODTZERO;     //! TZERO AOD
  AliAODVZERO     *fAODVZERO;     //! VZERO AOD
  AliAODZDC       *fAODZDC;       //! ZDC AOD
  AliAODAD        *fAODAD;        //! AD AOD
  AliTOFHeader    *fTOFHeader;  //! event times (and sigmas) as estimated by TOF
			     //  combinatorial algorithm.
                             //  It contains also TOF time resolution
                             //  and T0spread as written in OCDB
  TClonesArray    *fTrdTracks;    //! TRD AOD tracks (triggered)
  
  static const char* fAODListName[kAODListN]; //!

  ClassDef(AliAODEvent,94);
};

#endif
