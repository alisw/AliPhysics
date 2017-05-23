// Class for extended track information
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 

#ifndef ALIREDUCEDTRACKINFO_H
#define ALIREDUCEDTRACKINFO_H

#include <TMath.h>
#include "AliReducedBaseTrack.h"

//_____________________________________________________________________
class AliReducedTrackInfo : public AliReducedBaseTrack {

  friend class AliAnalysisTaskReducedTreeMaker;  // friend analysis task which fills the object
  
 public:
  AliReducedTrackInfo();
  virtual ~AliReducedTrackInfo();

  // getters
  UShort_t TrackId()                     const {return fTrackId;}
  ULong_t  Status()                      const {return fStatus;}
  Bool_t   CheckTrackStatus(UInt_t flag) const {return (flag<8*sizeof(ULong_t) ? (fStatus&(1<<flag)) : kFALSE);}
  Float_t  PxTPC()                       const {return fTPCPt*TMath::Cos(fTPCPhi);}
  Float_t  PyTPC()                       const {return fTPCPt*TMath::Sin(fTPCPhi);}
  Float_t  PzTPC()                       const {return fTPCPt*TMath::SinH(fTPCEta);}
  Float_t  PTPC()                        const {return fTPCPt*TMath::CosH(fTPCEta);};
  Float_t  PhiTPC()                      const {return fTPCPhi;}
  Float_t  PtTPC()                       const {return fTPCPt;}
  Float_t  EtaTPC()                      const {return fTPCEta;}
  Float_t  ThetaTPC()                    const {return TMath::ACos(TMath::TanH(fTPCEta));}
  Float_t  DCAxyTPC()                    const {return fTPCDCA[0];}
  Float_t  DCAzTPC()                     const {return fTPCDCA[1];}
  Float_t  Pin()                         const {return fMomentumInner;}
  Float_t  DCAxy()                       const {return fDCA[0];}
  Float_t  DCAz()                        const {return fDCA[1];}
  Float_t  TrackLength()                 const {return fTrackLength;}
  Float_t  MassForTracking()        const {return fMassForTracking;}
  Float_t  Chi2TPCConstrainedVsGlobal()  const {return fChi2TPCConstrainedVsGlobal;}
  Float_t  HelixX()                      const {return fHelixCenter[0];}
  Float_t  HelixY()                      const {return fHelixCenter[1];}
  Float_t  HelixR()                      const {return fHelixRadius;}
  
  UShort_t ITSncls()                const;
  UChar_t  ITSclusterMap()          const {return fITSclusterMap;}
  UChar_t  ITSSharedClusterMap()   const {return fITSSharedClusterMap;}
  Bool_t   ITSLayerHit(Int_t layer) const {return (layer>=0 && layer<6 ? (fITSclusterMap&(1<<layer)) : kFALSE);};
  Bool_t   ITSClusterIsShared(Int_t layer) const {return (layer>=0 && layer<6 ? (fITSSharedClusterMap&(1<<layer)) : kFALSE);};
  UShort_t     ITSnSharedCls() const;
  Float_t  ITSsignal()              const {return fITSsignal;}
  Float_t  ITSnSig(Int_t specie)    const {return (specie>=0 && specie<=3 ? fITSnSig[specie] : -999.);}
  Float_t  ITSchi2()                const {return fITSchi2;}
  
  UChar_t TPCncls()                        const {return fTPCNcls;}
  UChar_t TPCFindableNcls()                const {return fTPCNclsF;}
  UChar_t TPCCrossedRows()                 const {return fTPCCrossedRows;}
  UChar_t TPCnclsShared()                  const {return fTPCNclsShared;}
  UChar_t TPCClusterMap()                  const {return fTPCClusterMap;}
  Int_t   TPCClusterMapBitsFired()         const;
  Bool_t  TPCClusterMapBitFired(Int_t bit) const {return (bit>=0 && bit<8 ? (fTPCClusterMap&(1<<bit)) : kFALSE);};
  Float_t TPCsignal()                      const {return fTPCsignal;}
  UChar_t TPCsignalN()                     const {return fTPCsignalN;}
  Float_t TPCnSig(Int_t specie)            const {return (specie>=0 && specie<=3 ? fTPCnSig[specie] : -999.);}
  Float_t TPCchi2()                        const {return fTPCchi2;}
  Float_t TPCActiveLength()         const {return fTPCActiveLength;}
  Float_t TPCGeomLength()       const {return fTPCGeomLength;}
  
  Float_t  TOFbeta()             const {return fTOFbeta;}
  Float_t  TOFtime()             const {return fTOFtime;}
  Float_t  TOFdx()               const {return fTOFdx;}
  Float_t  TOFdz()               const {return fTOFdz;}
  Float_t  TOFmismatchProbab()   const {return fTOFmismatchProbab;}
  Float_t  TOFchi2()             const {return fTOFchi2;}
  Float_t  TOFnSig(Int_t specie) const {return (specie>=0 && specie<=3 ? fTOFnSig[specie] : -999.);}
  Short_t  TOFdeltaBC()          const {return fTOFdeltaBC;}
  
  Int_t    TRDntracklets(Int_t type)  const {return (type==0 || type==1 ? fTRDntracklets[type] : -1);}
  Float_t  TRDpid(Int_t specie)       const {return (specie>=0 && specie<=1 ? fTRDpid[specie] : -999.);}
  Float_t  TRDpidLQ1D(Int_t specie)   const {return (specie>=0 && specie<=1 ? fTRDpid[specie] : -999.);}
  Float_t  TRDpidLQ2D(Int_t specie)   const {return (specie>=0 && specie<=1 ? fTRDpidLQ2D[specie] : -999.);}
  
  Int_t    CaloClusterId() const {return fCaloClusterId;}
  
  Float_t TrackParam(Int_t iPar = 0) {return (iPar>=0 && iPar<6 ? fTrackParam[iPar] : 0.0);}
  Float_t CovMatrix(Int_t iCov = 0) {return (iCov>=0 && iCov<21 ? fCovMatrix[iCov] : 0.0);}
  
  Float_t MCmom(Int_t dim) {return (dim>=0 && dim<3 ? fMCMom[dim] : 0.0);}
  Float_t PtMC() {return TMath::Sqrt(fMCMom[0]*fMCMom[0]+fMCMom[1]*fMCMom[1]);}
  Float_t PMC()   const {return TMath::Sqrt(fMCMom[0]*fMCMom[0]+fMCMom[1]*fMCMom[1]+fMCMom[2]*fMCMom[2]);}
  Float_t PhiMC() const;
  Float_t EtaMC() const;
  Float_t RapidityMC(Float_t massAssumption) const;
  Float_t ThetaMC() const;
  Float_t EnergyMC(Float_t massAssumption) const {return TMath::Sqrt(massAssumption*massAssumption+PMC()*PMC());}
  Float_t MCFreezeout(Int_t dim) {return (dim>=0 && dim<3 ? fMCFreezeout[dim] : 0.0);}
  Int_t MCLabel(Int_t history=0) {return (history>=0 && history<4 ? fMCLabels[history] : -9999);}
  Int_t MCPdg(Int_t history=0) {return (history>=0 && history<4 ? fMCPdg[history] : -9999);}
  Short_t MCGeneratorIndex() {return fMCGeneratorIndex;}
  

     
 protected:
  UShort_t fTrackId;            // track id 
  ULong_t fStatus;              // tracking status
  Float_t fTPCPhi;              // inner param phi
  Float_t fTPCPt;               // inner param pt  
  Float_t fTPCEta;              // inner param eta 
  Float_t fMomentumInner;       // inner param momentum (only the magnitude)
  Float_t fDCA[2];              // DCA xy,z
  Float_t fTPCDCA[2];           // TPConly DCA xy,z
  Float_t fTrackLength;         // track length
  Float_t fMassForTracking;    // mass hypothesis used for tracking
  Float_t fChi2TPCConstrainedVsGlobal;   // AliESDtrack::GetChi2TPCConstrainedVsGlobal()
  Float_t fHelixCenter[2];      // Helix Center x,y
  Float_t fHelixRadius;         // Radius of the Helix
  
  // ITS
  UChar_t  fITSclusterMap;      // ITS cluster map
  UChar_t  fITSSharedClusterMap;  // ITS shared cluster map
  Float_t  fITSsignal;          // ITS signal
  Float_t  fITSnSig[4];         // 0-electron; 1-pion; 2-kaon; 3-proton
  Float_t  fITSchi2;            // ITS chi2 / cls
  
  // TPC
  UChar_t fTPCNcls;            // TPC ncls                          
  UChar_t fTPCCrossedRows;     // TPC crossed rows                  
  UChar_t fTPCNclsF;           // TPC findable ncls                 
  UChar_t fTPCNclsShared;      // TPC number of shared clusters
  UChar_t fTPCClusterMap;      // TPC cluster distribution map
  Float_t fTPCsignal;          // TPC de/dx
  UChar_t fTPCsignalN;         // TPC no clusters de/dx
  Float_t fTPCnSig[4];         // 0-electron; 1-pion; 2-kaon; 3-proton
  Float_t fTPCchi2;            // TPC chi2 / cls
  Float_t fTPCActiveLength;   // track length in active parts of the TPC
  Float_t fTPCGeomLength;   // geometrical track length in the TPC
    
  // TOF
  Float_t fTOFbeta;             // TOF pid info
  Float_t fTOFtime;             // TOF flight time
  Float_t fTOFdx;               // TOF matching dx
  Float_t fTOFdz;               // TOF matching dz
  Float_t fTOFmismatchProbab;   // TOF mismatch probability
  Float_t fTOFchi2;             // TOF chi2
  Float_t fTOFnSig[4];          // TOF n-sigma deviation from expected signal
  Short_t fTOFdeltaBC;          // BC(event) - BC(track) estimated by TOF

  // TRD
  UChar_t fTRDntracklets[2];       // 0 - AliESDtrack::GetTRDntracklets(); 1 - AliESDtrack::GetTRDntrackletsPID()   TODO: use only 1 char
  Float_t fTRDpid[2];              // TRD pid 1D likelihoods, [0]-electron , [1]- pion
  Float_t fTRDpidLQ2D[2];          // TRD pid 2D likelihoods, [0]-electron , [1]- pion
  
  // EMCAL/PHOS
  Int_t  fCaloClusterId;          // ID for the calorimeter cluster (if any)
  
  // Track parameters stored at the primary vertex
  Float_t fTrackParam[6];     // parameters: x, y, z, px, py, pz
  Float_t fCovMatrix[21];     // covariance matrix for the track parameter
  
  
  // Monte-Carlo truth information
  Float_t fMCMom[3];             // MC truth 3-momentum information in cartezian coordinates
  Float_t fMCFreezeout[3];    // MC truth 3-position information in cartezian coordinates
  Int_t    fMCLabels[4];           // MC label for: [0] - the current track, [1] - mother, [2] - grand mother, [3] - grand grand mother 
  Int_t    fMCPdg[4];                // MC PDG code for: [0] - the current track, [1] - mother, [2] - grand mother, [3] - grand grand mother 
  Short_t fMCGeneratorIndex;    // generator index (used for cocktail generators ?)
  

  
          
  AliReducedTrackInfo(const AliReducedTrackInfo &c);
  AliReducedTrackInfo& operator= (const AliReducedTrackInfo &c);

  ClassDef(AliReducedTrackInfo, 4);
};

//_______________________________________________________________________________
inline Float_t AliReducedTrackInfo::PhiMC() const {
   //
   // Return the azimuthal angle of this particle
   //
   Float_t phi=TMath::ATan2(fMCMom[1],fMCMom[0]); 
   if(phi>=0.0) 
      return phi;
   else 
      return (TMath::TwoPi()+phi);
}

//_______________________________________________________________________________
inline Float_t AliReducedTrackInfo::ThetaMC() const {
   //
   // Return the polar angle for this particle
   //
   Float_t p=PMC(); 
   if(p>=1.0e-6) 
      return TMath::ACos(fMCMom[2]/p);
   else 
      return 0.0;
}

//_______________________________________________________________________________
inline Float_t AliReducedTrackInfo::EtaMC() const {
   //
   // Return the pseudorapidity of this particle
   //
   Float_t eta = TMath::Tan(0.5*ThetaMC());
   if(eta>1.0e-6) 
      return -1.0*TMath::Log(eta);
   else 
      return 0.0;
}

//_______________________________________________________________________________
inline Float_t AliReducedTrackInfo::RapidityMC(Float_t massAssumption) const {
   //
   // Return the rapidity of this particle using a massAssumption
   //
   Float_t e = EnergyMC(massAssumption);
   Float_t factor = e-fMCMom[2];
   if(TMath::Abs(factor)<1.0e-6) return 0.0;
   factor = (e+fMCMom[2])/factor;
   if(factor<1.0e-6) return 0.0;
   return 0.5*TMath::Log(factor);
}


//_______________________________________________________________________________
inline UShort_t AliReducedTrackInfo::ITSncls() const
{
  //
  // ITS number of clusters from the cluster map
  //
  UShort_t ncls=0;
  for(Int_t i=0; i<6; ++i) ncls += (ITSLayerHit(i) ? 1 : 0);
  return ncls;
}

//_______________________________________________________________________________
inline UShort_t AliReducedTrackInfo::ITSnSharedCls() const
{
   //
   // ITS number of shared clusters from the shared cluster map
   //
   UShort_t ncls=0;
   for(Int_t i=0; i<6; ++i) ncls += (ITSClusterIsShared(i) ? 1 : 0);
   return ncls;
}


//_______________________________________________________________________________
inline Int_t AliReducedTrackInfo::TPCClusterMapBitsFired()  const
{
  //
  // Count the number of bits fired in the TPC cluster map
  //
  Int_t nbits=0;
  for(Int_t i=0; i<8; ++i) nbits += (TPCClusterMapBitFired(i) ? 1 : 0);
  return nbits;
}

#endif
