//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 17 18:28:08 2006 by ROOT version 5.11/06
// from TTree esdTree/Tree with ESD objects
// found on file: root://lxgate06.cern.ch//alice/cern.ch/user/a/aliprod/psaiz/prod2006/output_pp/1/002/root_archive.zip#AliESDs.root
//////////////////////////////////////////////////////////

#ifndef esdAna_h
#define esdAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1F.h>
   const Int_t kMaxfTracks = 62;
   const Int_t kMaxfHLTConfMapTracks = 37;
   const Int_t kMaxfHLTHoughTracks = 1;
   const Int_t kMaxfMuonTracks = 1;
   const Int_t kMaxfPmdTracks = 1;
   const Int_t kMaxfTrdTracks = 1;
   const Int_t kMaxfV0s = 1;
   const Int_t kMaxfCascades = 1;
   const Int_t kMaxfKinks = 1;
   const Int_t kMaxfV0MIs = 1;
   const Int_t kMaxfCaloClusters = 45;

class esdAna : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leave types
// AliESD          *fESD;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           fEventNumber;
   Int_t           fRunNumber;
   ULong64_t       fTriggerMask;
   UChar_t         fTriggerCluster;
   Int_t           fRecoVersion;
   Float_t         fMagneticField;
   Float_t         fZDCN1Energy;
   Float_t         fZDCP1Energy;
   Float_t         fZDCN2Energy;
   Float_t         fZDCP2Energy;
   Float_t         fZDCEMEnergy;
   Int_t           fZDCParticipants;
   Float_t         fT0zVertex;
   UInt_t          fSPDVertex_fUniqueID;
   UInt_t          fSPDVertex_fBits;
   TString         fSPDVertex_fName;
   TString         fSPDVertex_fTitle;
   Double_t        fSPDVertex_fPosition[3];
   Double_t        fSPDVertex_fSigma;
   Int_t           fSPDVertex_fNContributors;
   Double_t        fSPDVertex_fCovXX;
   Double_t        fSPDVertex_fCovXY;
   Double_t        fSPDVertex_fCovYY;
   Double_t        fSPDVertex_fCovXZ;
   Double_t        fSPDVertex_fCovYZ;
   Double_t        fSPDVertex_fCovZZ;
   Double_t        fSPDVertex_fSNR[3];
   Double_t        fSPDVertex_fChi2;
   Double_t        fSPDVertex_fTruePos[3];
   UInt_t          fPrimaryVertex_fUniqueID;
   UInt_t          fPrimaryVertex_fBits;
   TString         fPrimaryVertex_fName;
   TString         fPrimaryVertex_fTitle;
   Double_t        fPrimaryVertex_fPosition[3];
   Double_t        fPrimaryVertex_fSigma;
   Int_t           fPrimaryVertex_fNContributors;
   Double_t        fPrimaryVertex_fCovXX;
   Double_t        fPrimaryVertex_fCovXY;
   Double_t        fPrimaryVertex_fCovYY;
   Double_t        fPrimaryVertex_fCovXZ;
   Double_t        fPrimaryVertex_fCovYZ;
   Double_t        fPrimaryVertex_fCovZZ;
   Double_t        fPrimaryVertex_fSNR[3];
   Double_t        fPrimaryVertex_fChi2;
   Double_t        fPrimaryVertex_fTruePos[3];
   UInt_t          fSPDMult_fUniqueID;
   UInt_t          fSPDMult_fBits;
   Int_t           fSPDMult_fNtracks;
   Float_t         fSPDMult_fTh[25];   //[fSPDMult.fNtracks]
   Float_t         fSPDMult_fPhi[25];   //[fSPDMult.fNtracks]
   Float_t         fSPDMult_fDeltPhi[25];   //[fSPDMult.fNtracks]
   Float_t         fT0timeStart;
   Float_t         fT0time[24];
   Float_t         fT0amplitude[24];
   Int_t           fTracks_;
   UInt_t          fTracks_fUniqueID[kMaxfTracks];   //[fTracks_]
   UInt_t          fTracks_fBits[kMaxfTracks];   //[fTracks_]
   Double_t        fTracks_fX[kMaxfTracks];   //[fTracks_]
   Double_t        fTracks_fAlpha[kMaxfTracks];   //[fTracks_]
   Double_t        fTracks_fP[kMaxfTracks][5];   //[fTracks_]
   Double_t        fTracks_fC[kMaxfTracks][15];   //[fTracks_]
   ULong_t         fTracks_fFlags[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fLabel[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fID[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTrackLength[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fD[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fZ[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fCdd[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fCdz[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fCzz[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTrackTime[kMaxfTracks][5];   //[fTracks_]
   Float_t         fTracks_fR[kMaxfTracks][5];   //[fTracks_]
   Int_t           fTracks_fStopVertex[kMaxfTracks];   //[fTracks_]
   Double_t        fTracks_fCchi2[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fITSchi2[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fITSncls[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fITSsignal[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fITSr[kMaxfTracks][5];   //[fTracks_]
   Int_t           fTracks_fITSLabel[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fITSFakeRatio[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTPCchi2[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fTPCncls[kMaxfTracks];   //[fTracks_]
   UShort_t        fTracks_fTPCnclsF[kMaxfTracks];   //[fTracks_]
   UInt_t          fTracks_fTPCClusterMap_fUniqueID[kMaxfTracks];   //[fTracks_]
   UInt_t          fTracks_fTPCClusterMap_fBits[kMaxfTracks];   //[fTracks_]
   UInt_t          fTracks_fTPCClusterMap_fNbits[kMaxfTracks];   //[fTracks_]
   UInt_t          fTracks_fTPCClusterMap_fNbytes[kMaxfTracks];   //[fTracks_]
   UChar_t        *fTracks_fTPCClusterMap_fAllBits[kMaxfTracks];   //[fTracks_fTPCClusterMap_fNbytes]
   Float_t         fTracks_fTPCsignal[kMaxfTracks];   //[fTracks_]
   UShort_t        fTracks_fTPCsignalN[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTPCsignalS[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTPCr[kMaxfTracks][5];   //[fTracks_]
   Int_t           fTracks_fTPCLabel[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTPCPoints[kMaxfTracks][4];   //[fTracks_]
   Int_t           fTracks_fKinkIndexes[kMaxfTracks][3];   //[fTracks_]
   Int_t           fTracks_fV0Indexes[kMaxfTracks][3];   //[fTracks_]
   Float_t         fTracks_fTRDchi2[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fTRDncls[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fTRDncls0[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTRDsignal[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTRDsignals[kMaxfTracks][6][3];   //[fTracks_]
   Int_t           fTracks_fTRDTimBin[kMaxfTracks][6];   //[fTracks_]
   Float_t         fTracks_fTRDr[kMaxfTracks][5];   //[fTracks_]
   Int_t           fTracks_fTRDLabel[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTRDQuality[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTRDBudget[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTOFchi2[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fTOFindex[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fTOFCalChannel[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTOFsignal[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTOFsignalToT[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fTOFr[kMaxfTracks][5];   //[fTracks_]
   Int_t           fTracks_fTOFLabel[kMaxfTracks][3];   //[fTracks_]
   Float_t         fTracks_fHMPIDchi2[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fHMPIDncls[kMaxfTracks];   //[fTracks_]
   Int_t           fTracks_fHMPIDindex[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fHMPIDsignal[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fHMPIDr[kMaxfTracks][5];   //[fTracks_]
   Float_t         fTracks_fHMPIDtheta[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fHMPIDphi[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fHMPIDdx[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fHMPIDdy[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fHMPIDmipX[kMaxfTracks];   //[fTracks_]
   Float_t         fTracks_fHMPIDmipY[kMaxfTracks];   //[fTracks_]
   Int_t           fHLTConfMapTracks_;
   UInt_t          fHLTConfMapTracks_fUniqueID[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   UInt_t          fHLTConfMapTracks_fBits[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   UShort_t        fHLTConfMapTracks_fNHits[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Int_t           fHLTConfMapTracks_fMCid[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   UShort_t        fHLTConfMapTracks_fWeight[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Bool_t          fHLTConfMapTracks_fFromMainVertex[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Int_t           fHLTConfMapTracks_fRowRange[kMaxfHLTConfMapTracks][2];   //[fHLTConfMapTracks_]
   UShort_t        fHLTConfMapTracks_fSector[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fFirstPoint[kMaxfHLTConfMapTracks][3];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fLastPoint[kMaxfHLTConfMapTracks][3];   //[fHLTConfMapTracks_]
   Int_t           fHLTConfMapTracks_fQ[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fTanl[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fPsi[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fPt[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fPterr[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fPsierr[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fTanlerr[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fBinX[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fBinY[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fSizeX[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fSizeY[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Float_t         fHLTConfMapTracks_fPID[kMaxfHLTConfMapTracks];   //[fHLTConfMapTracks_]
   Int_t           fHLTHoughTracks_;
   UInt_t          fHLTHoughTracks_fUniqueID[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   UInt_t          fHLTHoughTracks_fBits[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   UShort_t        fHLTHoughTracks_fNHits[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Int_t           fHLTHoughTracks_fMCid[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   UShort_t        fHLTHoughTracks_fWeight[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Bool_t          fHLTHoughTracks_fFromMainVertex[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Int_t           fHLTHoughTracks_fRowRange[kMaxfHLTHoughTracks][2];   //[fHLTHoughTracks_]
   UShort_t        fHLTHoughTracks_fSector[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fFirstPoint[kMaxfHLTHoughTracks][3];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fLastPoint[kMaxfHLTHoughTracks][3];   //[fHLTHoughTracks_]
   Int_t           fHLTHoughTracks_fQ[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fTanl[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fPsi[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fPt[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fPterr[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fPsierr[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fTanlerr[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fBinX[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fBinY[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fSizeX[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fSizeY[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Float_t         fHLTHoughTracks_fPID[kMaxfHLTHoughTracks];   //[fHLTHoughTracks_]
   Int_t           fMuonTracks_;
   UInt_t          fMuonTracks_fUniqueID[kMaxfMuonTracks];   //[fMuonTracks_]
   UInt_t          fMuonTracks_fBits[kMaxfMuonTracks];   //[fMuonTracks_]
   Double_t        fMuonTracks_fInverseBendingMomentum[kMaxfMuonTracks];   //[fMuonTracks_]
   Double_t        fMuonTracks_fThetaX[kMaxfMuonTracks];   //[fMuonTracks_]
   Double_t        fMuonTracks_fThetaY[kMaxfMuonTracks];   //[fMuonTracks_]
   Double_t        fMuonTracks_fZ[kMaxfMuonTracks];   //[fMuonTracks_]
   Double_t        fMuonTracks_fBendingCoor[kMaxfMuonTracks];   //[fMuonTracks_]
   Double_t        fMuonTracks_fNonBendingCoor[kMaxfMuonTracks];   //[fMuonTracks_]
   Double_t        fMuonTracks_fChi2[kMaxfMuonTracks];   //[fMuonTracks_]
   UInt_t          fMuonTracks_fNHit[kMaxfMuonTracks];   //[fMuonTracks_]
   Bool_t          fMuonTracks_fMatchTrigger[kMaxfMuonTracks];   //[fMuonTracks_]
   Double_t        fMuonTracks_fChi2MatchTrigger[kMaxfMuonTracks];   //[fMuonTracks_]
   Int_t           fPmdTracks_;
   UInt_t          fPmdTracks_fUniqueID[kMaxfPmdTracks];   //[fPmdTracks_]
   UInt_t          fPmdTracks_fBits[kMaxfPmdTracks];   //[fPmdTracks_]
   Int_t           fPmdTracks_fDet[kMaxfPmdTracks];   //[fPmdTracks_]
   Float_t         fPmdTracks_fX[kMaxfPmdTracks];   //[fPmdTracks_]
   Float_t         fPmdTracks_fY[kMaxfPmdTracks];   //[fPmdTracks_]
   Float_t         fPmdTracks_fZ[kMaxfPmdTracks];   //[fPmdTracks_]
   Float_t         fPmdTracks_fCluADC[kMaxfPmdTracks];   //[fPmdTracks_]
   Float_t         fPmdTracks_fNcell[kMaxfPmdTracks];   //[fPmdTracks_]
   Float_t         fPmdTracks_fCluPID[kMaxfPmdTracks];   //[fPmdTracks_]
   Int_t           fTrdTracks_;
   UInt_t          fTrdTracks_fUniqueID[kMaxfTrdTracks];   //[fTrdTracks_]
   UInt_t          fTrdTracks_fBits[kMaxfTrdTracks];   //[fTrdTracks_]
   Float_t         fTrdTracks_fYproj[kMaxfTrdTracks];   //[fTrdTracks_]
   Float_t         fTrdTracks_fZproj[kMaxfTrdTracks];   //[fTrdTracks_]
   Float_t         fTrdTracks_fSlope[kMaxfTrdTracks];   //[fTrdTracks_]
   Int_t           fTrdTracks_fDetector[kMaxfTrdTracks];   //[fTrdTracks_]
   Int_t           fTrdTracks_fNtracklets[kMaxfTrdTracks];   //[fTrdTracks_]
   Int_t           fTrdTracks_fNplanes[kMaxfTrdTracks];   //[fTrdTracks_]
   Int_t           fTrdTracks_fNclusters[kMaxfTrdTracks];   //[fTrdTracks_]
   Float_t         fTrdTracks_fPt[kMaxfTrdTracks];   //[fTrdTracks_]
   Float_t         fTrdTracks_fPhi[kMaxfTrdTracks];   //[fTrdTracks_]
   Float_t         fTrdTracks_fEta[kMaxfTrdTracks];   //[fTrdTracks_]
   Int_t           fTrdTracks_fLabel[kMaxfTrdTracks];   //[fTrdTracks_]
   Float_t         fTrdTracks_fPID[kMaxfTrdTracks];   //[fTrdTracks_]
   Bool_t          fTrdTracks_fIsElectron[kMaxfTrdTracks];   //[fTrdTracks_]
   Int_t           fV0s_;
   UInt_t          fV0s_fUniqueID[kMaxfV0s];   //[fV0s_]
   UInt_t          fV0s_fBits[kMaxfV0s];   //[fV0s_]
   Int_t           fV0s_fPdgCode[kMaxfV0s];   //[fV0s_]
   Double_t        fV0s_fEffMass[kMaxfV0s];   //[fV0s_]
   Double_t        fV0s_fDcaDaughters[kMaxfV0s];   //[fV0s_]
   Double_t        fV0s_fChi2[kMaxfV0s];   //[fV0s_]
   Double_t        fV0s_fPos[kMaxfV0s][3];   //[fV0s_]
   Double_t        fV0s_fPosCov[kMaxfV0s][6];   //[fV0s_]
   Int_t           fV0s_fNidx[kMaxfV0s];   //[fV0s_]
   Double_t        fV0s_fNmom[kMaxfV0s][3];   //[fV0s_]
   Double_t        fV0s_fNmomCov[kMaxfV0s][6];   //[fV0s_]
   Int_t           fV0s_fPidx[kMaxfV0s];   //[fV0s_]
   Double_t        fV0s_fPmom[kMaxfV0s][3];   //[fV0s_]
   Double_t        fV0s_fPmomCov[kMaxfV0s][6];   //[fV0s_]
   Int_t           fCascades_;
   UInt_t          fCascades_fUniqueID[kMaxfCascades];   //[fCascades_]
   UInt_t          fCascades_fBits[kMaxfCascades];   //[fCascades_]
   Int_t           fCascades_fPdgCode[kMaxfCascades];   //[fCascades_]
   Double_t        fCascades_fEffMass[kMaxfCascades];   //[fCascades_]
   Double_t        fCascades_fChi2[kMaxfCascades];   //[fCascades_]
   Double_t        fCascades_fPos[kMaxfCascades][3];   //[fCascades_]
   Double_t        fCascades_fPosCov[kMaxfCascades][6];   //[fCascades_]
   Int_t           fCascades_fV0idx[kMaxfCascades][2];   //[fCascades_]
   Double_t        fCascades_fV0mom[kMaxfCascades][2][3];   //[fCascades_]
   Double_t        fCascades_fV0momCov[kMaxfCascades][6];   //[fCascades_]
   Int_t           fCascades_fBachIdx[kMaxfCascades];   //[fCascades_]
   Double_t        fCascades_fBachMom[kMaxfCascades][3];   //[fCascades_]
   Double_t        fCascades_fBachMomCov[kMaxfCascades][6];   //[fCascades_]
   Int_t           fKinks_;
   UInt_t          fKinks_fUniqueID[kMaxfKinks];   //[fKinks_]
   UInt_t          fKinks_fBits[kMaxfKinks];   //[fKinks_]
   Int_t           fKinks_fID[kMaxfKinks];   //[fKinks_]
   UInt_t          fKinks_fParamDaughter_fUniqueID[kMaxfKinks];   //[fKinks_]
   UInt_t          fKinks_fParamDaughter_fBits[kMaxfKinks];   //[fKinks_]
   Double_t        fKinks_fParamDaughter_fX[kMaxfKinks];   //[fKinks_]
   Double_t        fKinks_fParamDaughter_fAlpha[kMaxfKinks];   //[fKinks_]
   Double_t        fKinks_fParamDaughter_fP[kMaxfKinks][5];   //[fKinks_]
   Double_t        fKinks_fParamDaughter_fC[kMaxfKinks][15];   //[fKinks_]
   UInt_t          fKinks_fParamMother_fUniqueID[kMaxfKinks];   //[fKinks_]
   UInt_t          fKinks_fParamMother_fBits[kMaxfKinks];   //[fKinks_]
   Double_t        fKinks_fParamMother_fX[kMaxfKinks];   //[fKinks_]
   Double_t        fKinks_fParamMother_fAlpha[kMaxfKinks];   //[fKinks_]
   Double_t        fKinks_fParamMother_fP[kMaxfKinks][5];   //[fKinks_]
   Double_t        fKinks_fParamMother_fC[kMaxfKinks][15];   //[fKinks_]
   Double_t        fKinks_fDist1[kMaxfKinks];   //[fKinks_]
   Double_t        fKinks_fDist2[kMaxfKinks];   //[fKinks_]
   Double_t        fKinks_fPdr[kMaxfKinks][3];   //[fKinks_]
   Double_t        fKinks_fXr[kMaxfKinks][3];   //[fKinks_]
   Double_t        fKinks_fPm[kMaxfKinks][3];   //[fKinks_]
   Double_t        fKinks_fAngle[kMaxfKinks][3];   //[fKinks_]
   Double_t        fKinks_fRr[kMaxfKinks];   //[fKinks_]
   Int_t           fKinks_fLab[kMaxfKinks][2];   //[fKinks_]
   Int_t           fKinks_fIndex[kMaxfKinks][2];   //[fKinks_]
   Char_t          fKinks_fStatus[kMaxfKinks][12];   //[fKinks_]
   Float_t         fKinks_fTPCdensity[kMaxfKinks][2][2];   //[fKinks_]
   Float_t         fKinks_fTPCdensity2[kMaxfKinks][2][2];   //[fKinks_]
   Float_t         fKinks_fShapeFactor[kMaxfKinks];   //[fKinks_]
   Int_t           fKinks_fRow0[kMaxfKinks];   //[fKinks_]
   Int_t           fKinks_fMultiple[kMaxfKinks][2];   //[fKinks_]
   Int_t           fKinks_fTPCncls[kMaxfKinks][2];   //[fKinks_]
   Int_t           fV0MIs_;
   UInt_t          fV0MIs_fUniqueID[kMaxfV0MIs];   //[fV0MIs_]
   UInt_t          fV0MIs_fBits[kMaxfV0MIs];   //[fV0MIs_]
   Int_t           fV0MIs_fPdgCode[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fEffMass[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fDcaDaughters[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fChi2[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fPos[kMaxfV0MIs][3];   //[fV0MIs_]
   Double_t        fV0MIs_fPosCov[kMaxfV0MIs][6];   //[fV0MIs_]
   Int_t           fV0MIs_fNidx[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fNmom[kMaxfV0MIs][3];   //[fV0MIs_]
   Double_t        fV0MIs_fNmomCov[kMaxfV0MIs][6];   //[fV0MIs_]
   Int_t           fV0MIs_fPidx[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fPmom[kMaxfV0MIs][3];   //[fV0MIs_]
   Double_t        fV0MIs_fPmomCov[kMaxfV0MIs][6];   //[fV0MIs_]
   UInt_t          fV0MIs_fParamP_fUniqueID[kMaxfV0MIs];   //[fV0MIs_]
   UInt_t          fV0MIs_fParamP_fBits[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fParamP_fX[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fParamP_fAlpha[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fParamP_fP[kMaxfV0MIs][5];   //[fV0MIs_]
   Double_t        fV0MIs_fParamP_fC[kMaxfV0MIs][15];   //[fV0MIs_]
   UInt_t          fV0MIs_fParamM_fUniqueID[kMaxfV0MIs];   //[fV0MIs_]
   UInt_t          fV0MIs_fParamM_fBits[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fParamM_fX[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fParamM_fAlpha[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fParamM_fP[kMaxfV0MIs][5];   //[fV0MIs_]
   Double_t        fV0MIs_fParamM_fC[kMaxfV0MIs][15];   //[fV0MIs_]
   Float_t         fV0MIs_fRP[kMaxfV0MIs][5];   //[fV0MIs_]
   Float_t         fV0MIs_fRM[kMaxfV0MIs][5];   //[fV0MIs_]
   Int_t           fV0MIs_fID[kMaxfV0MIs];   //[fV0MIs_]
   Int_t           fV0MIs_fLab[kMaxfV0MIs][2];   //[fV0MIs_]
   Int_t           fV0MIs_fIndex[kMaxfV0MIs][2];   //[fV0MIs_]
   Float_t         fV0MIs_fNormDCAPrim[kMaxfV0MIs][2];   //[fV0MIs_]
   Double_t        fV0MIs_fDist1[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fDist2[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fPP[kMaxfV0MIs][3];   //[fV0MIs_]
   Double_t        fV0MIs_fPM[kMaxfV0MIs][3];   //[fV0MIs_]
   Double_t        fV0MIs_fXr[kMaxfV0MIs][3];   //[fV0MIs_]
   Double_t        fV0MIs_fAngle[kMaxfV0MIs][3];   //[fV0MIs_]
   Double_t        fV0MIs_fRr[kMaxfV0MIs];   //[fV0MIs_]
   Int_t           fV0MIs_fStatus[kMaxfV0MIs];   //[fV0MIs_]
   Int_t           fV0MIs_fRow0[kMaxfV0MIs];   //[fV0MIs_]
   Int_t           fV0MIs_fOrder[kMaxfV0MIs][3];   //[fV0MIs_]
   Double_t        fV0MIs_fDistNorm[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fDistSigma[kMaxfV0MIs];   //[fV0MIs_]
   Float_t         fV0MIs_fCausality[kMaxfV0MIs][4];   //[fV0MIs_]
   Float_t         fV0MIs_fChi2Before[kMaxfV0MIs];   //[fV0MIs_]
   Float_t         fV0MIs_fNBefore[kMaxfV0MIs];   //[fV0MIs_]
   Float_t         fV0MIs_fChi2After[kMaxfV0MIs];   //[fV0MIs_]
   Float_t         fV0MIs_fNAfter[kMaxfV0MIs];   //[fV0MIs_]
   Float_t         fV0MIs_fPointAngleFi[kMaxfV0MIs];   //[fV0MIs_]
   Float_t         fV0MIs_fPointAngleTh[kMaxfV0MIs];   //[fV0MIs_]
   Double_t        fV0MIs_fPointAngle[kMaxfV0MIs];   //[fV0MIs_]
   Int_t           fCaloClusters_;
   UInt_t          fCaloClusters_fUniqueID[kMaxfCaloClusters];   //[fCaloClusters_]
   UInt_t          fCaloClusters_fBits[kMaxfCaloClusters];   //[fCaloClusters_]
   Int_t           fCaloClusters_fID[kMaxfCaloClusters];   //[fCaloClusters_]
   Int_t           fCaloClusters_fClusterType[kMaxfCaloClusters];   //[fCaloClusters_]
   Bool_t          fCaloClusters_fEMCALCluster[kMaxfCaloClusters];   //[fCaloClusters_]
   Bool_t          fCaloClusters_fPHOSCluster[kMaxfCaloClusters];   //[fCaloClusters_]
   Float_t         fCaloClusters_fGlobalPos[kMaxfCaloClusters][3];   //[fCaloClusters_]
   Float_t         fCaloClusters_fEnergy[kMaxfCaloClusters];   //[fCaloClusters_]
   Float_t         fCaloClusters_fDispersion[kMaxfCaloClusters];   //[fCaloClusters_]
   Float_t         fCaloClusters_fChi2[kMaxfCaloClusters];   //[fCaloClusters_]
   Float_t         fCaloClusters_fPID[kMaxfCaloClusters][10];   //[fCaloClusters_]
   Int_t           fCaloClusters_fPrimaryIndex[kMaxfCaloClusters];   //[fCaloClusters_]
   Float_t         fCaloClusters_fM20[kMaxfCaloClusters];   //[fCaloClusters_]
   Float_t         fCaloClusters_fM02[kMaxfCaloClusters];   //[fCaloClusters_]
   Float_t         fCaloClusters_fM11[kMaxfCaloClusters];   //[fCaloClusters_]
   UShort_t        fCaloClusters_fNExMax[kMaxfCaloClusters];   //[fCaloClusters_]
   Float_t         fCaloClusters_fEmcCpvDistance[kMaxfCaloClusters];   //[fCaloClusters_]
   Int_t           fCaloClusters_fNumberOfDigits[kMaxfCaloClusters];   //[fCaloClusters_]
   UShort_t       *fCaloClusters_fDigitAmplitude[kMaxfCaloClusters];   //[fCaloClusters_fNumberOfDigits]
   UShort_t       *fCaloClusters_fDigitTime[kMaxfCaloClusters];   //[fCaloClusters_fNumberOfDigits]
   UShort_t       *fCaloClusters_fDigitIndex[kMaxfCaloClusters];   //[fCaloClusters_fNumberOfDigits]
   Int_t           fEMCALClusters;
   Int_t           fFirstEMCALCluster;
   Int_t           fPHOSClusters;
   Int_t           fFirstPHOSCluster;

   // List of branches
   TBranch        *b_ESD_fUniqueID;   //!
   TBranch        *b_ESD_fBits;   //!
   TBranch        *b_ESD_fEventNumber;   //!
   TBranch        *b_ESD_fRunNumber;   //!
   TBranch        *b_ESD_fTriggerMask;   //!
   TBranch        *b_ESD_fTriggerCluster;   //!
   TBranch        *b_ESD_fRecoVersion;   //!
   TBranch        *b_ESD_fMagneticField;   //!
   TBranch        *b_ESD_fZDCN1Energy;   //!
   TBranch        *b_ESD_fZDCP1Energy;   //!
   TBranch        *b_ESD_fZDCN2Energy;   //!
   TBranch        *b_ESD_fZDCP2Energy;   //!
   TBranch        *b_ESD_fZDCEMEnergy;   //!
   TBranch        *b_ESD_fZDCParticipants;   //!
   TBranch        *b_ESD_fT0zVertex;   //!
   TBranch        *b_ESD_fSPDVertex_fUniqueID;   //!
   TBranch        *b_ESD_fSPDVertex_fBits;   //!
   TBranch        *b_ESD_fSPDVertex_fName;   //!
   TBranch        *b_ESD_fSPDVertex_fTitle;   //!
   TBranch        *b_ESD_fSPDVertex_fPosition;   //!
   TBranch        *b_ESD_fSPDVertex_fSigma;   //!
   TBranch        *b_ESD_fSPDVertex_fNContributors;   //!
   TBranch        *b_ESD_fSPDVertex_fCovXX;   //!
   TBranch        *b_ESD_fSPDVertex_fCovXY;   //!
   TBranch        *b_ESD_fSPDVertex_fCovYY;   //!
   TBranch        *b_ESD_fSPDVertex_fCovXZ;   //!
   TBranch        *b_ESD_fSPDVertex_fCovYZ;   //!
   TBranch        *b_ESD_fSPDVertex_fCovZZ;   //!
   TBranch        *b_ESD_fSPDVertex_fSNR;   //!
   TBranch        *b_ESD_fSPDVertex_fChi2;   //!
   TBranch        *b_ESD_fSPDVertex_fTruePos;   //!
   TBranch        *b_ESD_fPrimaryVertex_fUniqueID;   //!
   TBranch        *b_ESD_fPrimaryVertex_fBits;   //!
   TBranch        *b_ESD_fPrimaryVertex_fName;   //!
   TBranch        *b_ESD_fPrimaryVertex_fTitle;   //!
   TBranch        *b_ESD_fPrimaryVertex_fPosition;   //!
   TBranch        *b_ESD_fPrimaryVertex_fSigma;   //!
   TBranch        *b_ESD_fPrimaryVertex_fNContributors;   //!
   TBranch        *b_ESD_fPrimaryVertex_fCovXX;   //!
   TBranch        *b_ESD_fPrimaryVertex_fCovXY;   //!
   TBranch        *b_ESD_fPrimaryVertex_fCovYY;   //!
   TBranch        *b_ESD_fPrimaryVertex_fCovXZ;   //!
   TBranch        *b_ESD_fPrimaryVertex_fCovYZ;   //!
   TBranch        *b_ESD_fPrimaryVertex_fCovZZ;   //!
   TBranch        *b_ESD_fPrimaryVertex_fSNR;   //!
   TBranch        *b_ESD_fPrimaryVertex_fChi2;   //!
   TBranch        *b_ESD_fPrimaryVertex_fTruePos;   //!
   TBranch        *b_ESD_fSPDMult_fUniqueID;   //!
   TBranch        *b_ESD_fSPDMult_fBits;   //!
   TBranch        *b_ESD_fSPDMult_fNtracks;   //!
   TBranch        *b_fSPDMult_fTh;   //!
   TBranch        *b_fSPDMult_fPhi;   //!
   TBranch        *b_fSPDMult_fDeltPhi;   //!
   TBranch        *b_ESD_fT0timeStart;   //!
   TBranch        *b_ESD_fT0time;   //!
   TBranch        *b_ESD_fT0amplitude;   //!
   TBranch        *b_ESD_fTracks_;   //!
   TBranch        *b_fTracks_fUniqueID;   //!
   TBranch        *b_fTracks_fBits;   //!
   TBranch        *b_fTracks_fX;   //!
   TBranch        *b_fTracks_fAlpha;   //!
   TBranch        *b_fTracks_fP;   //!
   TBranch        *b_fTracks_fC;   //!
   TBranch        *b_fTracks_fFlags;   //!
   TBranch        *b_fTracks_fLabel;   //!
   TBranch        *b_fTracks_fID;   //!
   TBranch        *b_fTracks_fTrackLength;   //!
   TBranch        *b_fTracks_fD;   //!
   TBranch        *b_fTracks_fZ;   //!
   TBranch        *b_fTracks_fCdd;   //!
   TBranch        *b_fTracks_fCdz;   //!
   TBranch        *b_fTracks_fCzz;   //!
   TBranch        *b_fTracks_fTrackTime;   //!
   TBranch        *b_fTracks_fR;   //!
   TBranch        *b_fTracks_fStopVertex;   //!
   TBranch        *b_fTracks_fCchi2;   //!
   TBranch        *b_fTracks_fITSchi2;   //!
   TBranch        *b_fTracks_fITSncls;   //!
   TBranch        *b_fTracks_fITSsignal;   //!
   TBranch        *b_fTracks_fITSr;   //!
   TBranch        *b_fTracks_fITSLabel;   //!
   TBranch        *b_fTracks_fITSFakeRatio;   //!
   TBranch        *b_fTracks_fTPCchi2;   //!
   TBranch        *b_fTracks_fTPCncls;   //!
   TBranch        *b_fTracks_fTPCnclsF;   //!
   TBranch        *b_fTracks_fTPCClusterMap_fUniqueID;   //!
   TBranch        *b_fTracks_fTPCClusterMap_fBits;   //!
   TBranch        *b_fTracks_fTPCClusterMap_fNbits;   //!
   TBranch        *b_fTracks_fTPCClusterMap_fNbytes;   //!
   TBranch        *b_fTracks_fTPCClusterMap_fAllBits;   //!
   TBranch        *b_fTracks_fTPCsignal;   //!
   TBranch        *b_fTracks_fTPCsignalN;   //!
   TBranch        *b_fTracks_fTPCsignalS;   //!
   TBranch        *b_fTracks_fTPCr;   //!
   TBranch        *b_fTracks_fTPCLabel;   //!
   TBranch        *b_fTracks_fTPCPoints;   //!
   TBranch        *b_fTracks_fKinkIndexes;   //!
   TBranch        *b_fTracks_fV0Indexes;   //!
   TBranch        *b_fTracks_fTRDchi2;   //!
   TBranch        *b_fTracks_fTRDncls;   //!
   TBranch        *b_fTracks_fTRDncls0;   //!
   TBranch        *b_fTracks_fTRDsignal;   //!
   TBranch        *b_fTracks_fTRDsignals;   //!
   TBranch        *b_fTracks_fTRDTimBin;   //!
   TBranch        *b_fTracks_fTRDr;   //!
   TBranch        *b_fTracks_fTRDLabel;   //!
   TBranch        *b_fTracks_fTRDQuality;   //!
   TBranch        *b_fTracks_fTRDBudget;   //!
   TBranch        *b_fTracks_fTOFchi2;   //!
   TBranch        *b_fTracks_fTOFindex;   //!
   TBranch        *b_fTracks_fTOFCalChannel;   //!
   TBranch        *b_fTracks_fTOFsignal;   //!
   TBranch        *b_fTracks_fTOFsignalToT;   //!
   TBranch        *b_fTracks_fTOFr;   //!
   TBranch        *b_fTracks_fTOFLabel;   //!
   TBranch        *b_fTracks_fHMPIDchi2;   //!
   TBranch        *b_fTracks_fHMPIDncls;   //!
   TBranch        *b_fTracks_fHMPIDindex;   //!
   TBranch        *b_fTracks_fHMPIDsignal;   //!
   TBranch        *b_fTracks_fHMPIDr;   //!
   TBranch        *b_fTracks_fHMPIDtheta;   //!
   TBranch        *b_fTracks_fHMPIDphi;   //!
   TBranch        *b_fTracks_fHMPIDdx;   //!
   TBranch        *b_fTracks_fHMPIDdy;   //!
   TBranch        *b_fTracks_fHMPIDmipX;   //!
   TBranch        *b_fTracks_fHMPIDmipY;   //!
   TBranch        *b_ESD_fHLTConfMapTracks_;   //!
   TBranch        *b_fHLTConfMapTracks_fUniqueID;   //!
   TBranch        *b_fHLTConfMapTracks_fBits;   //!
   TBranch        *b_fHLTConfMapTracks_fNHits;   //!
   TBranch        *b_fHLTConfMapTracks_fMCid;   //!
   TBranch        *b_fHLTConfMapTracks_fWeight;   //!
   TBranch        *b_fHLTConfMapTracks_fFromMainVertex;   //!
   TBranch        *b_fHLTConfMapTracks_fRowRange;   //!
   TBranch        *b_fHLTConfMapTracks_fSector;   //!
   TBranch        *b_fHLTConfMapTracks_fFirstPoint;   //!
   TBranch        *b_fHLTConfMapTracks_fLastPoint;   //!
   TBranch        *b_fHLTConfMapTracks_fQ;   //!
   TBranch        *b_fHLTConfMapTracks_fTanl;   //!
   TBranch        *b_fHLTConfMapTracks_fPsi;   //!
   TBranch        *b_fHLTConfMapTracks_fPt;   //!
   TBranch        *b_fHLTConfMapTracks_fPterr;   //!
   TBranch        *b_fHLTConfMapTracks_fPsierr;   //!
   TBranch        *b_fHLTConfMapTracks_fTanlerr;   //!
   TBranch        *b_fHLTConfMapTracks_fBinX;   //!
   TBranch        *b_fHLTConfMapTracks_fBinY;   //!
   TBranch        *b_fHLTConfMapTracks_fSizeX;   //!
   TBranch        *b_fHLTConfMapTracks_fSizeY;   //!
   TBranch        *b_fHLTConfMapTracks_fPID;   //!
   TBranch        *b_ESD_fHLTHoughTracks_;   //!
   TBranch        *b_fHLTHoughTracks_fUniqueID;   //!
   TBranch        *b_fHLTHoughTracks_fBits;   //!
   TBranch        *b_fHLTHoughTracks_fNHits;   //!
   TBranch        *b_fHLTHoughTracks_fMCid;   //!
   TBranch        *b_fHLTHoughTracks_fWeight;   //!
   TBranch        *b_fHLTHoughTracks_fFromMainVertex;   //!
   TBranch        *b_fHLTHoughTracks_fRowRange;   //!
   TBranch        *b_fHLTHoughTracks_fSector;   //!
   TBranch        *b_fHLTHoughTracks_fFirstPoint;   //!
   TBranch        *b_fHLTHoughTracks_fLastPoint;   //!
   TBranch        *b_fHLTHoughTracks_fQ;   //!
   TBranch        *b_fHLTHoughTracks_fTanl;   //!
   TBranch        *b_fHLTHoughTracks_fPsi;   //!
   TBranch        *b_fHLTHoughTracks_fPt;   //!
   TBranch        *b_fHLTHoughTracks_fPterr;   //!
   TBranch        *b_fHLTHoughTracks_fPsierr;   //!
   TBranch        *b_fHLTHoughTracks_fTanlerr;   //!
   TBranch        *b_fHLTHoughTracks_fBinX;   //!
   TBranch        *b_fHLTHoughTracks_fBinY;   //!
   TBranch        *b_fHLTHoughTracks_fSizeX;   //!
   TBranch        *b_fHLTHoughTracks_fSizeY;   //!
   TBranch        *b_fHLTHoughTracks_fPID;   //!
   TBranch        *b_ESD_fMuonTracks_;   //!
   TBranch        *b_fMuonTracks_fUniqueID;   //!
   TBranch        *b_fMuonTracks_fBits;   //!
   TBranch        *b_fMuonTracks_fInverseBendingMomentum;   //!
   TBranch        *b_fMuonTracks_fThetaX;   //!
   TBranch        *b_fMuonTracks_fThetaY;   //!
   TBranch        *b_fMuonTracks_fZ;   //!
   TBranch        *b_fMuonTracks_fBendingCoor;   //!
   TBranch        *b_fMuonTracks_fNonBendingCoor;   //!
   TBranch        *b_fMuonTracks_fChi2;   //!
   TBranch        *b_fMuonTracks_fNHit;   //!
   TBranch        *b_fMuonTracks_fMatchTrigger;   //!
   TBranch        *b_fMuonTracks_fChi2MatchTrigger;   //!
   TBranch        *b_ESD_fPmdTracks_;   //!
   TBranch        *b_fPmdTracks_fUniqueID;   //!
   TBranch        *b_fPmdTracks_fBits;   //!
   TBranch        *b_fPmdTracks_fDet;   //!
   TBranch        *b_fPmdTracks_fX;   //!
   TBranch        *b_fPmdTracks_fY;   //!
   TBranch        *b_fPmdTracks_fZ;   //!
   TBranch        *b_fPmdTracks_fCluADC;   //!
   TBranch        *b_fPmdTracks_fNcell;   //!
   TBranch        *b_fPmdTracks_fCluPID;   //!
   TBranch        *b_ESD_fTrdTracks_;   //!
   TBranch        *b_fTrdTracks_fUniqueID;   //!
   TBranch        *b_fTrdTracks_fBits;   //!
   TBranch        *b_fTrdTracks_fYproj;   //!
   TBranch        *b_fTrdTracks_fZproj;   //!
   TBranch        *b_fTrdTracks_fSlope;   //!
   TBranch        *b_fTrdTracks_fDetector;   //!
   TBranch        *b_fTrdTracks_fNtracklets;   //!
   TBranch        *b_fTrdTracks_fNplanes;   //!
   TBranch        *b_fTrdTracks_fNclusters;   //!
   TBranch        *b_fTrdTracks_fPt;   //!
   TBranch        *b_fTrdTracks_fPhi;   //!
   TBranch        *b_fTrdTracks_fEta;   //!
   TBranch        *b_fTrdTracks_fLabel;   //!
   TBranch        *b_fTrdTracks_fPID;   //!
   TBranch        *b_fTrdTracks_fIsElectron;   //!
   TBranch        *b_ESD_fV0s_;   //!
   TBranch        *b_fV0s_fUniqueID;   //!
   TBranch        *b_fV0s_fBits;   //!
   TBranch        *b_fV0s_fPdgCode;   //!
   TBranch        *b_fV0s_fEffMass;   //!
   TBranch        *b_fV0s_fDcaDaughters;   //!
   TBranch        *b_fV0s_fChi2;   //!
   TBranch        *b_fV0s_fPos;   //!
   TBranch        *b_fV0s_fPosCov;   //!
   TBranch        *b_fV0s_fNidx;   //!
   TBranch        *b_fV0s_fNmom;   //!
   TBranch        *b_fV0s_fNmomCov;   //!
   TBranch        *b_fV0s_fPidx;   //!
   TBranch        *b_fV0s_fPmom;   //!
   TBranch        *b_fV0s_fPmomCov;   //!
   TBranch        *b_ESD_fCascades_;   //!
   TBranch        *b_fCascades_fUniqueID;   //!
   TBranch        *b_fCascades_fBits;   //!
   TBranch        *b_fCascades_fPdgCode;   //!
   TBranch        *b_fCascades_fEffMass;   //!
   TBranch        *b_fCascades_fChi2;   //!
   TBranch        *b_fCascades_fPos;   //!
   TBranch        *b_fCascades_fPosCov;   //!
   TBranch        *b_fCascades_fV0idx;   //!
   TBranch        *b_fCascades_fV0mom;   //!
   TBranch        *b_fCascades_fV0momCov;   //!
   TBranch        *b_fCascades_fBachIdx;   //!
   TBranch        *b_fCascades_fBachMom;   //!
   TBranch        *b_fCascades_fBachMomCov;   //!
   TBranch        *b_ESD_fKinks_;   //!
   TBranch        *b_fKinks_fUniqueID;   //!
   TBranch        *b_fKinks_fBits;   //!
   TBranch        *b_fKinks_fID;   //!
   TBranch        *b_fKinks_fParamDaughter_fUniqueID;   //!
   TBranch        *b_fKinks_fParamDaughter_fBits;   //!
   TBranch        *b_fKinks_fParamDaughter_fX;   //!
   TBranch        *b_fKinks_fParamDaughter_fAlpha;   //!
   TBranch        *b_fKinks_fParamDaughter_fP;   //!
   TBranch        *b_fKinks_fParamDaughter_fC;   //!
   TBranch        *b_fKinks_fParamMother_fUniqueID;   //!
   TBranch        *b_fKinks_fParamMother_fBits;   //!
   TBranch        *b_fKinks_fParamMother_fX;   //!
   TBranch        *b_fKinks_fParamMother_fAlpha;   //!
   TBranch        *b_fKinks_fParamMother_fP;   //!
   TBranch        *b_fKinks_fParamMother_fC;   //!
   TBranch        *b_fKinks_fDist1;   //!
   TBranch        *b_fKinks_fDist2;   //!
   TBranch        *b_fKinks_fPdr;   //!
   TBranch        *b_fKinks_fXr;   //!
   TBranch        *b_fKinks_fPm;   //!
   TBranch        *b_fKinks_fAngle;   //!
   TBranch        *b_fKinks_fRr;   //!
   TBranch        *b_fKinks_fLab;   //!
   TBranch        *b_fKinks_fIndex;   //!
   TBranch        *b_fKinks_fStatus;   //!
   TBranch        *b_fKinks_fTPCdensity;   //!
   TBranch        *b_fKinks_fTPCdensity2;   //!
   TBranch        *b_fKinks_fShapeFactor;   //!
   TBranch        *b_fKinks_fRow0;   //!
   TBranch        *b_fKinks_fMultiple;   //!
   TBranch        *b_fKinks_fTPCncls;   //!
   TBranch        *b_ESD_fV0MIs_;   //!
   TBranch        *b_fV0MIs_fUniqueID;   //!
   TBranch        *b_fV0MIs_fBits;   //!
   TBranch        *b_fV0MIs_fPdgCode;   //!
   TBranch        *b_fV0MIs_fEffMass;   //!
   TBranch        *b_fV0MIs_fDcaDaughters;   //!
   TBranch        *b_fV0MIs_fChi2;   //!
   TBranch        *b_fV0MIs_fPos;   //!
   TBranch        *b_fV0MIs_fPosCov;   //!
   TBranch        *b_fV0MIs_fNidx;   //!
   TBranch        *b_fV0MIs_fNmom;   //!
   TBranch        *b_fV0MIs_fNmomCov;   //!
   TBranch        *b_fV0MIs_fPidx;   //!
   TBranch        *b_fV0MIs_fPmom;   //!
   TBranch        *b_fV0MIs_fPmomCov;   //!
   TBranch        *b_fV0MIs_fParamP_fUniqueID;   //!
   TBranch        *b_fV0MIs_fParamP_fBits;   //!
   TBranch        *b_fV0MIs_fParamP_fX;   //!
   TBranch        *b_fV0MIs_fParamP_fAlpha;   //!
   TBranch        *b_fV0MIs_fParamP_fP;   //!
   TBranch        *b_fV0MIs_fParamP_fC;   //!
   TBranch        *b_fV0MIs_fParamM_fUniqueID;   //!
   TBranch        *b_fV0MIs_fParamM_fBits;   //!
   TBranch        *b_fV0MIs_fParamM_fX;   //!
   TBranch        *b_fV0MIs_fParamM_fAlpha;   //!
   TBranch        *b_fV0MIs_fParamM_fP;   //!
   TBranch        *b_fV0MIs_fParamM_fC;   //!
   TBranch        *b_fV0MIs_fRP;   //!
   TBranch        *b_fV0MIs_fRM;   //!
   TBranch        *b_fV0MIs_fID;   //!
   TBranch        *b_fV0MIs_fLab;   //!
   TBranch        *b_fV0MIs_fIndex;   //!
   TBranch        *b_fV0MIs_fNormDCAPrim;   //!
   TBranch        *b_fV0MIs_fDist1;   //!
   TBranch        *b_fV0MIs_fDist2;   //!
   TBranch        *b_fV0MIs_fPP;   //!
   TBranch        *b_fV0MIs_fPM;   //!
   TBranch        *b_fV0MIs_fXr;   //!
   TBranch        *b_fV0MIs_fAngle;   //!
   TBranch        *b_fV0MIs_fRr;   //!
   TBranch        *b_fV0MIs_fStatus;   //!
   TBranch        *b_fV0MIs_fRow0;   //!
   TBranch        *b_fV0MIs_fOrder;   //!
   TBranch        *b_fV0MIs_fDistNorm;   //!
   TBranch        *b_fV0MIs_fDistSigma;   //!
   TBranch        *b_fV0MIs_fCausality;   //!
   TBranch        *b_fV0MIs_fChi2Before;   //!
   TBranch        *b_fV0MIs_fNBefore;   //!
   TBranch        *b_fV0MIs_fChi2After;   //!
   TBranch        *b_fV0MIs_fNAfter;   //!
   TBranch        *b_fV0MIs_fPointAngleFi;   //!
   TBranch        *b_fV0MIs_fPointAngleTh;   //!
   TBranch        *b_fV0MIs_fPointAngle;   //!
   TBranch        *b_ESD_fCaloClusters_;   //!
   TBranch        *b_fCaloClusters_fUniqueID;   //!
   TBranch        *b_fCaloClusters_fBits;   //!
   TBranch        *b_fCaloClusters_fID;   //!
   TBranch        *b_fCaloClusters_fClusterType;   //!
   TBranch        *b_fCaloClusters_fEMCALCluster;   //!
   TBranch        *b_fCaloClusters_fPHOSCluster;   //!
   TBranch        *b_fCaloClusters_fGlobalPos;   //!
   TBranch        *b_fCaloClusters_fEnergy;   //!
   TBranch        *b_fCaloClusters_fDispersion;   //!
   TBranch        *b_fCaloClusters_fChi2;   //!
   TBranch        *b_fCaloClusters_fPID;   //!
   TBranch        *b_fCaloClusters_fPrimaryIndex;   //!
   TBranch        *b_fCaloClusters_fM20;   //!
   TBranch        *b_fCaloClusters_fM02;   //!
   TBranch        *b_fCaloClusters_fM11;   //!
   TBranch        *b_fCaloClusters_fNExMax;   //!
   TBranch        *b_fCaloClusters_fEmcCpvDistance;   //!
   TBranch        *b_fCaloClusters_fNumberOfDigits;   //!
   TBranch        *b_fCaloClusters_fDigitAmplitude;   //!
   TBranch        *b_fCaloClusters_fDigitTime;   //!
   TBranch        *b_fCaloClusters_fDigitIndex;   //!
   TBranch        *b_ESD_fEMCALClusters;   //!
   TBranch        *b_ESD_fFirstEMCALCluster;   //!
   TBranch        *b_ESD_fPHOSClusters;   //!
   TBranch        *b_ESD_fFirstPHOSCluster;   //!

   esdAna(TTree *tree=0) { }
   virtual ~esdAna() { }
   virtual Int_t   Version() const { return 1; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   TH1F* h1;
   TH1F* h2;
   TH1F* h3;
   TFile* hfile;
   ClassDef(esdAna,0);
};

#endif

#ifdef esdAna_cxx
void esdAna::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses of the tree
   // will be set. It is normaly not necessary to make changes to the
   // generated code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running with PROOF.

   // Set branch addresses
   if (tree == 0) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("fUniqueID",&fUniqueID);
   fChain->SetBranchAddress("fBits",&fBits);
   fChain->SetBranchAddress("fEventNumber",&fEventNumber);
   fChain->SetBranchAddress("fRunNumber",&fRunNumber);
   fChain->SetBranchAddress("fTriggerMask",&fTriggerMask);
   fChain->SetBranchAddress("fTriggerCluster",&fTriggerCluster);
   fChain->SetBranchAddress("fRecoVersion",&fRecoVersion);
   fChain->SetBranchAddress("fMagneticField",&fMagneticField);
   fChain->SetBranchAddress("fZDCN1Energy",&fZDCN1Energy);
   fChain->SetBranchAddress("fZDCP1Energy",&fZDCP1Energy);
   fChain->SetBranchAddress("fZDCN2Energy",&fZDCN2Energy);
   fChain->SetBranchAddress("fZDCP2Energy",&fZDCP2Energy);
   fChain->SetBranchAddress("fZDCEMEnergy",&fZDCEMEnergy);
   fChain->SetBranchAddress("fZDCParticipants",&fZDCParticipants);
   fChain->SetBranchAddress("fT0zVertex",&fT0zVertex);
   fChain->SetBranchAddress("fSPDVertex.fUniqueID",&fSPDVertex_fUniqueID);
   fChain->SetBranchAddress("fSPDVertex.fBits",&fSPDVertex_fBits);
   fChain->SetBranchAddress("fSPDVertex.fName",&fSPDVertex_fName);
   fChain->SetBranchAddress("fSPDVertex.fTitle",&fSPDVertex_fTitle);
   fChain->SetBranchAddress("fSPDVertex.fPosition[3]",fSPDVertex_fPosition);
   fChain->SetBranchAddress("fSPDVertex.fSigma",&fSPDVertex_fSigma);
   fChain->SetBranchAddress("fSPDVertex.fNContributors",&fSPDVertex_fNContributors);
   fChain->SetBranchAddress("fSPDVertex.fCovXX",&fSPDVertex_fCovXX);
   fChain->SetBranchAddress("fSPDVertex.fCovXY",&fSPDVertex_fCovXY);
   fChain->SetBranchAddress("fSPDVertex.fCovYY",&fSPDVertex_fCovYY);
   fChain->SetBranchAddress("fSPDVertex.fCovXZ",&fSPDVertex_fCovXZ);
   fChain->SetBranchAddress("fSPDVertex.fCovYZ",&fSPDVertex_fCovYZ);
   fChain->SetBranchAddress("fSPDVertex.fCovZZ",&fSPDVertex_fCovZZ);
   fChain->SetBranchAddress("fSPDVertex.fSNR[3]",fSPDVertex_fSNR);
   fChain->SetBranchAddress("fSPDVertex.fChi2",&fSPDVertex_fChi2);
   fChain->SetBranchAddress("fSPDVertex.fTruePos[3]",fSPDVertex_fTruePos);
   fChain->SetBranchAddress("fPrimaryVertex.fUniqueID",&fPrimaryVertex_fUniqueID);
   fChain->SetBranchAddress("fPrimaryVertex.fBits",&fPrimaryVertex_fBits);
   fChain->SetBranchAddress("fPrimaryVertex.fName",&fPrimaryVertex_fName);
   fChain->SetBranchAddress("fPrimaryVertex.fTitle",&fPrimaryVertex_fTitle);
   fChain->SetBranchAddress("fPrimaryVertex.fPosition[3]",fPrimaryVertex_fPosition);
   fChain->SetBranchAddress("fPrimaryVertex.fSigma",&fPrimaryVertex_fSigma);
   fChain->SetBranchAddress("fPrimaryVertex.fNContributors",&fPrimaryVertex_fNContributors);
   fChain->SetBranchAddress("fPrimaryVertex.fCovXX",&fPrimaryVertex_fCovXX);
   fChain->SetBranchAddress("fPrimaryVertex.fCovXY",&fPrimaryVertex_fCovXY);
   fChain->SetBranchAddress("fPrimaryVertex.fCovYY",&fPrimaryVertex_fCovYY);
   fChain->SetBranchAddress("fPrimaryVertex.fCovXZ",&fPrimaryVertex_fCovXZ);
   fChain->SetBranchAddress("fPrimaryVertex.fCovYZ",&fPrimaryVertex_fCovYZ);
   fChain->SetBranchAddress("fPrimaryVertex.fCovZZ",&fPrimaryVertex_fCovZZ);
   fChain->SetBranchAddress("fPrimaryVertex.fSNR[3]",fPrimaryVertex_fSNR);
   fChain->SetBranchAddress("fPrimaryVertex.fChi2",&fPrimaryVertex_fChi2);
   fChain->SetBranchAddress("fPrimaryVertex.fTruePos[3]",fPrimaryVertex_fTruePos);
   fChain->SetBranchAddress("fSPDMult.fUniqueID",&fSPDMult_fUniqueID);
   fChain->SetBranchAddress("fSPDMult.fBits",&fSPDMult_fBits);
   fChain->SetBranchAddress("fSPDMult.fNtracks",&fSPDMult_fNtracks);
   fChain->SetBranchAddress("fSPDMult.fTh",fSPDMult_fTh);
   fChain->SetBranchAddress("fSPDMult.fPhi",fSPDMult_fPhi);
   fChain->SetBranchAddress("fSPDMult.fDeltPhi",fSPDMult_fDeltPhi);
   fChain->SetBranchAddress("fT0timeStart",&fT0timeStart);
   fChain->SetBranchAddress("fT0time[24]",fT0time);
   fChain->SetBranchAddress("fT0amplitude[24]",fT0amplitude);
   fChain->SetBranchAddress("fTracks",&fTracks_);
   fChain->SetBranchAddress("fTracks.fUniqueID",fTracks_fUniqueID);
   fChain->SetBranchAddress("fTracks.fBits",fTracks_fBits);
   fChain->SetBranchAddress("fTracks.fX",fTracks_fX);
   fChain->SetBranchAddress("fTracks.fAlpha",fTracks_fAlpha);
   fChain->SetBranchAddress("fTracks.fP[5]",fTracks_fP);
   fChain->SetBranchAddress("fTracks.fC[15]",fTracks_fC);
   fChain->SetBranchAddress("fTracks.fFlags",fTracks_fFlags);
   fChain->SetBranchAddress("fTracks.fLabel",fTracks_fLabel);
   fChain->SetBranchAddress("fTracks.fID",fTracks_fID);
   fChain->SetBranchAddress("fTracks.fTrackLength",fTracks_fTrackLength);
   fChain->SetBranchAddress("fTracks.fD",fTracks_fD);
   fChain->SetBranchAddress("fTracks.fZ",fTracks_fZ);
   fChain->SetBranchAddress("fTracks.fCdd",fTracks_fCdd);
   fChain->SetBranchAddress("fTracks.fCdz",fTracks_fCdz);
   fChain->SetBranchAddress("fTracks.fCzz",fTracks_fCzz);
   fChain->SetBranchAddress("fTracks.fTrackTime[5]",fTracks_fTrackTime);
   fChain->SetBranchAddress("fTracks.fR[5]",fTracks_fR);
   fChain->SetBranchAddress("fTracks.fStopVertex",fTracks_fStopVertex);
   fChain->SetBranchAddress("fTracks.fCchi2",fTracks_fCchi2);
   fChain->SetBranchAddress("fTracks.fITSchi2",fTracks_fITSchi2);
   fChain->SetBranchAddress("fTracks.fITSncls",fTracks_fITSncls);
   fChain->SetBranchAddress("fTracks.fITSsignal",fTracks_fITSsignal);
   fChain->SetBranchAddress("fTracks.fITSr[5]",fTracks_fITSr);
   fChain->SetBranchAddress("fTracks.fITSLabel",fTracks_fITSLabel);
   fChain->SetBranchAddress("fTracks.fITSFakeRatio",fTracks_fITSFakeRatio);
   fChain->SetBranchAddress("fTracks.fTPCchi2",fTracks_fTPCchi2);
   fChain->SetBranchAddress("fTracks.fTPCncls",fTracks_fTPCncls);
   fChain->SetBranchAddress("fTracks.fTPCnclsF",fTracks_fTPCnclsF);
   fChain->SetBranchAddress("fTracks.fTPCClusterMap.fUniqueID",fTracks_fTPCClusterMap_fUniqueID);
   fChain->SetBranchAddress("fTracks.fTPCClusterMap.fBits",fTracks_fTPCClusterMap_fBits);
   fChain->SetBranchAddress("fTracks.fTPCClusterMap.fNbits",fTracks_fTPCClusterMap_fNbits);
   fChain->SetBranchAddress("fTracks.fTPCClusterMap.fNbytes",fTracks_fTPCClusterMap_fNbytes);
   fChain->SetBranchAddress("fTracks.fTPCClusterMap.fAllBits",fTracks_fTPCClusterMap_fAllBits);
   fChain->SetBranchAddress("fTracks.fTPCsignal",fTracks_fTPCsignal);
   fChain->SetBranchAddress("fTracks.fTPCsignalN",fTracks_fTPCsignalN);
   fChain->SetBranchAddress("fTracks.fTPCsignalS",fTracks_fTPCsignalS);
   fChain->SetBranchAddress("fTracks.fTPCr[5]",fTracks_fTPCr);
   fChain->SetBranchAddress("fTracks.fTPCLabel",fTracks_fTPCLabel);
   fChain->SetBranchAddress("fTracks.fTPCPoints[4]",fTracks_fTPCPoints);
   fChain->SetBranchAddress("fTracks.fKinkIndexes[3]",fTracks_fKinkIndexes);
   fChain->SetBranchAddress("fTracks.fV0Indexes[3]",fTracks_fV0Indexes);
   fChain->SetBranchAddress("fTracks.fTRDchi2",fTracks_fTRDchi2);
   fChain->SetBranchAddress("fTracks.fTRDncls",fTracks_fTRDncls);
   fChain->SetBranchAddress("fTracks.fTRDncls0",fTracks_fTRDncls0);
   fChain->SetBranchAddress("fTracks.fTRDsignal",fTracks_fTRDsignal);
   fChain->SetBranchAddress("fTracks.fTRDsignals[6][3]",fTracks_fTRDsignals);
   fChain->SetBranchAddress("fTracks.fTRDTimBin[6]",fTracks_fTRDTimBin);
   fChain->SetBranchAddress("fTracks.fTRDr[5]",fTracks_fTRDr);
   fChain->SetBranchAddress("fTracks.fTRDLabel",fTracks_fTRDLabel);
   fChain->SetBranchAddress("fTracks.fTRDQuality",fTracks_fTRDQuality);
   fChain->SetBranchAddress("fTracks.fTRDBudget",fTracks_fTRDBudget);
   fChain->SetBranchAddress("fTracks.fTOFchi2",fTracks_fTOFchi2);
   fChain->SetBranchAddress("fTracks.fTOFindex",fTracks_fTOFindex);
   fChain->SetBranchAddress("fTracks.fTOFCalChannel",fTracks_fTOFCalChannel);
   fChain->SetBranchAddress("fTracks.fTOFsignal",fTracks_fTOFsignal);
   fChain->SetBranchAddress("fTracks.fTOFsignalToT",fTracks_fTOFsignalToT);
   fChain->SetBranchAddress("fTracks.fTOFr[5]",fTracks_fTOFr);
   fChain->SetBranchAddress("fTracks.fTOFLabel[3]",fTracks_fTOFLabel);
   fChain->SetBranchAddress("fTracks.fHMPIDchi2",fTracks_fHMPIDchi2);
   fChain->SetBranchAddress("fTracks.fHMPIDncls",fTracks_fHMPIDncls);
   fChain->SetBranchAddress("fTracks.fHMPIDindex",fTracks_fHMPIDindex);
   fChain->SetBranchAddress("fTracks.fHMPIDsignal",fTracks_fHMPIDsignal);
   fChain->SetBranchAddress("fTracks.fHMPIDr[5]",fTracks_fHMPIDr);
   fChain->SetBranchAddress("fTracks.fHMPIDtheta",fTracks_fHMPIDtheta);
   fChain->SetBranchAddress("fTracks.fHMPIDphi",fTracks_fHMPIDphi);
   fChain->SetBranchAddress("fTracks.fHMPIDdx",fTracks_fHMPIDdx);
   fChain->SetBranchAddress("fTracks.fHMPIDdy",fTracks_fHMPIDdy);
   fChain->SetBranchAddress("fTracks.fHMPIDmipX",fTracks_fHMPIDmipX);
   fChain->SetBranchAddress("fTracks.fHMPIDmipY",fTracks_fHMPIDmipY);
   fChain->SetBranchAddress("fHLTConfMapTracks",&fHLTConfMapTracks_);
   fChain->SetBranchAddress("fHLTConfMapTracks.fUniqueID",fHLTConfMapTracks_fUniqueID);
   fChain->SetBranchAddress("fHLTConfMapTracks.fBits",fHLTConfMapTracks_fBits);
   fChain->SetBranchAddress("fHLTConfMapTracks.fNHits",fHLTConfMapTracks_fNHits);
   fChain->SetBranchAddress("fHLTConfMapTracks.fMCid",fHLTConfMapTracks_fMCid);
   fChain->SetBranchAddress("fHLTConfMapTracks.fWeight",fHLTConfMapTracks_fWeight);
   fChain->SetBranchAddress("fHLTConfMapTracks.fFromMainVertex",fHLTConfMapTracks_fFromMainVertex);
   fChain->SetBranchAddress("fHLTConfMapTracks.fRowRange[2]",fHLTConfMapTracks_fRowRange);
   fChain->SetBranchAddress("fHLTConfMapTracks.fSector",fHLTConfMapTracks_fSector);
   fChain->SetBranchAddress("fHLTConfMapTracks.fFirstPoint[3]",fHLTConfMapTracks_fFirstPoint);
   fChain->SetBranchAddress("fHLTConfMapTracks.fLastPoint[3]",fHLTConfMapTracks_fLastPoint);
   fChain->SetBranchAddress("fHLTConfMapTracks.fQ",fHLTConfMapTracks_fQ);
   fChain->SetBranchAddress("fHLTConfMapTracks.fTanl",fHLTConfMapTracks_fTanl);
   fChain->SetBranchAddress("fHLTConfMapTracks.fPsi",fHLTConfMapTracks_fPsi);
   fChain->SetBranchAddress("fHLTConfMapTracks.fPt",fHLTConfMapTracks_fPt);
   fChain->SetBranchAddress("fHLTConfMapTracks.fPterr",fHLTConfMapTracks_fPterr);
   fChain->SetBranchAddress("fHLTConfMapTracks.fPsierr",fHLTConfMapTracks_fPsierr);
   fChain->SetBranchAddress("fHLTConfMapTracks.fTanlerr",fHLTConfMapTracks_fTanlerr);
   fChain->SetBranchAddress("fHLTConfMapTracks.fBinX",fHLTConfMapTracks_fBinX);
   fChain->SetBranchAddress("fHLTConfMapTracks.fBinY",fHLTConfMapTracks_fBinY);
   fChain->SetBranchAddress("fHLTConfMapTracks.fSizeX",fHLTConfMapTracks_fSizeX);
   fChain->SetBranchAddress("fHLTConfMapTracks.fSizeY",fHLTConfMapTracks_fSizeY);
   fChain->SetBranchAddress("fHLTConfMapTracks.fPID",fHLTConfMapTracks_fPID);
   fChain->SetBranchAddress("fHLTHoughTracks",&fHLTHoughTracks_);
   fChain->SetBranchAddress("fHLTHoughTracks.fUniqueID",&fHLTHoughTracks_fUniqueID);
   fChain->SetBranchAddress("fHLTHoughTracks.fBits",&fHLTHoughTracks_fBits);
   fChain->SetBranchAddress("fHLTHoughTracks.fNHits",&fHLTHoughTracks_fNHits);
   fChain->SetBranchAddress("fHLTHoughTracks.fMCid",&fHLTHoughTracks_fMCid);
   fChain->SetBranchAddress("fHLTHoughTracks.fWeight",&fHLTHoughTracks_fWeight);
   fChain->SetBranchAddress("fHLTHoughTracks.fFromMainVertex",&fHLTHoughTracks_fFromMainVertex);
   fChain->SetBranchAddress("fHLTHoughTracks.fRowRange[2]",&fHLTHoughTracks_fRowRange);
   fChain->SetBranchAddress("fHLTHoughTracks.fSector",&fHLTHoughTracks_fSector);
   fChain->SetBranchAddress("fHLTHoughTracks.fFirstPoint[3]",&fHLTHoughTracks_fFirstPoint);
   fChain->SetBranchAddress("fHLTHoughTracks.fLastPoint[3]",&fHLTHoughTracks_fLastPoint);
   fChain->SetBranchAddress("fHLTHoughTracks.fQ",&fHLTHoughTracks_fQ);
   fChain->SetBranchAddress("fHLTHoughTracks.fTanl",&fHLTHoughTracks_fTanl);
   fChain->SetBranchAddress("fHLTHoughTracks.fPsi",&fHLTHoughTracks_fPsi);
   fChain->SetBranchAddress("fHLTHoughTracks.fPt",&fHLTHoughTracks_fPt);
   fChain->SetBranchAddress("fHLTHoughTracks.fPterr",&fHLTHoughTracks_fPterr);
   fChain->SetBranchAddress("fHLTHoughTracks.fPsierr",&fHLTHoughTracks_fPsierr);
   fChain->SetBranchAddress("fHLTHoughTracks.fTanlerr",&fHLTHoughTracks_fTanlerr);
   fChain->SetBranchAddress("fHLTHoughTracks.fBinX",&fHLTHoughTracks_fBinX);
   fChain->SetBranchAddress("fHLTHoughTracks.fBinY",&fHLTHoughTracks_fBinY);
   fChain->SetBranchAddress("fHLTHoughTracks.fSizeX",&fHLTHoughTracks_fSizeX);
   fChain->SetBranchAddress("fHLTHoughTracks.fSizeY",&fHLTHoughTracks_fSizeY);
   fChain->SetBranchAddress("fHLTHoughTracks.fPID",&fHLTHoughTracks_fPID);
   fChain->SetBranchAddress("fMuonTracks",&fMuonTracks_);
   fChain->SetBranchAddress("fMuonTracks.fUniqueID",&fMuonTracks_fUniqueID);
   fChain->SetBranchAddress("fMuonTracks.fBits",&fMuonTracks_fBits);
   fChain->SetBranchAddress("fMuonTracks.fInverseBendingMomentum",&fMuonTracks_fInverseBendingMomentum);
   fChain->SetBranchAddress("fMuonTracks.fThetaX",&fMuonTracks_fThetaX);
   fChain->SetBranchAddress("fMuonTracks.fThetaY",&fMuonTracks_fThetaY);
   fChain->SetBranchAddress("fMuonTracks.fZ",&fMuonTracks_fZ);
   fChain->SetBranchAddress("fMuonTracks.fBendingCoor",&fMuonTracks_fBendingCoor);
   fChain->SetBranchAddress("fMuonTracks.fNonBendingCoor",&fMuonTracks_fNonBendingCoor);
   fChain->SetBranchAddress("fMuonTracks.fChi2",&fMuonTracks_fChi2);
   fChain->SetBranchAddress("fMuonTracks.fNHit",&fMuonTracks_fNHit);
   fChain->SetBranchAddress("fMuonTracks.fMatchTrigger",&fMuonTracks_fMatchTrigger);
   fChain->SetBranchAddress("fMuonTracks.fChi2MatchTrigger",&fMuonTracks_fChi2MatchTrigger);
   fChain->SetBranchAddress("fPmdTracks",&fPmdTracks_);
   fChain->SetBranchAddress("fPmdTracks.fUniqueID",&fPmdTracks_fUniqueID);
   fChain->SetBranchAddress("fPmdTracks.fBits",&fPmdTracks_fBits);
   fChain->SetBranchAddress("fPmdTracks.fDet",&fPmdTracks_fDet);
   fChain->SetBranchAddress("fPmdTracks.fX",&fPmdTracks_fX);
   fChain->SetBranchAddress("fPmdTracks.fY",&fPmdTracks_fY);
   fChain->SetBranchAddress("fPmdTracks.fZ",&fPmdTracks_fZ);
   fChain->SetBranchAddress("fPmdTracks.fCluADC",&fPmdTracks_fCluADC);
   fChain->SetBranchAddress("fPmdTracks.fNcell",&fPmdTracks_fNcell);
   fChain->SetBranchAddress("fPmdTracks.fCluPID",&fPmdTracks_fCluPID);
   fChain->SetBranchAddress("fTrdTracks",&fTrdTracks_);
   fChain->SetBranchAddress("fTrdTracks.fUniqueID",&fTrdTracks_fUniqueID);
   fChain->SetBranchAddress("fTrdTracks.fBits",&fTrdTracks_fBits);
   fChain->SetBranchAddress("fTrdTracks.fYproj",&fTrdTracks_fYproj);
   fChain->SetBranchAddress("fTrdTracks.fZproj",&fTrdTracks_fZproj);
   fChain->SetBranchAddress("fTrdTracks.fSlope",&fTrdTracks_fSlope);
   fChain->SetBranchAddress("fTrdTracks.fDetector",&fTrdTracks_fDetector);
   fChain->SetBranchAddress("fTrdTracks.fNtracklets",&fTrdTracks_fNtracklets);
   fChain->SetBranchAddress("fTrdTracks.fNplanes",&fTrdTracks_fNplanes);
   fChain->SetBranchAddress("fTrdTracks.fNclusters",&fTrdTracks_fNclusters);
   fChain->SetBranchAddress("fTrdTracks.fPt",&fTrdTracks_fPt);
   fChain->SetBranchAddress("fTrdTracks.fPhi",&fTrdTracks_fPhi);
   fChain->SetBranchAddress("fTrdTracks.fEta",&fTrdTracks_fEta);
   fChain->SetBranchAddress("fTrdTracks.fLabel",&fTrdTracks_fLabel);
   fChain->SetBranchAddress("fTrdTracks.fPID",&fTrdTracks_fPID);
   fChain->SetBranchAddress("fTrdTracks.fIsElectron",&fTrdTracks_fIsElectron);
   fChain->SetBranchAddress("fV0s",&fV0s_);
   fChain->SetBranchAddress("fV0s.fUniqueID",fV0s_fUniqueID);
   fChain->SetBranchAddress("fV0s.fBits",fV0s_fBits);
   fChain->SetBranchAddress("fV0s.fPdgCode",fV0s_fPdgCode);
   fChain->SetBranchAddress("fV0s.fEffMass",fV0s_fEffMass);
   fChain->SetBranchAddress("fV0s.fDcaDaughters",fV0s_fDcaDaughters);
   fChain->SetBranchAddress("fV0s.fChi2",fV0s_fChi2);
   fChain->SetBranchAddress("fV0s.fPos[3]",fV0s_fPos);
   fChain->SetBranchAddress("fV0s.fPosCov[6]",fV0s_fPosCov);
   fChain->SetBranchAddress("fV0s.fNidx",fV0s_fNidx);
   fChain->SetBranchAddress("fV0s.fNmom[3]",fV0s_fNmom);
   fChain->SetBranchAddress("fV0s.fNmomCov[6]",fV0s_fNmomCov);
   fChain->SetBranchAddress("fV0s.fPidx",fV0s_fPidx);
   fChain->SetBranchAddress("fV0s.fPmom[3]",fV0s_fPmom);
   fChain->SetBranchAddress("fV0s.fPmomCov[6]",fV0s_fPmomCov);
   fChain->SetBranchAddress("fCascades",&fCascades_);
   fChain->SetBranchAddress("fCascades.fUniqueID",&fCascades_fUniqueID);
   fChain->SetBranchAddress("fCascades.fBits",&fCascades_fBits);
   fChain->SetBranchAddress("fCascades.fPdgCode",&fCascades_fPdgCode);
   fChain->SetBranchAddress("fCascades.fEffMass",&fCascades_fEffMass);
   fChain->SetBranchAddress("fCascades.fChi2",&fCascades_fChi2);
   fChain->SetBranchAddress("fCascades.fPos[3]",&fCascades_fPos);
   fChain->SetBranchAddress("fCascades.fPosCov[6]",&fCascades_fPosCov);
   fChain->SetBranchAddress("fCascades.fV0idx[2]",&fCascades_fV0idx);
   fChain->SetBranchAddress("fCascades.fV0mom[2][3]",&fCascades_fV0mom);
   fChain->SetBranchAddress("fCascades.fV0momCov[6]",&fCascades_fV0momCov);
   fChain->SetBranchAddress("fCascades.fBachIdx",&fCascades_fBachIdx);
   fChain->SetBranchAddress("fCascades.fBachMom[3]",&fCascades_fBachMom);
   fChain->SetBranchAddress("fCascades.fBachMomCov[6]",&fCascades_fBachMomCov);
   fChain->SetBranchAddress("fKinks",&fKinks_);
   fChain->SetBranchAddress("fKinks.fUniqueID",fKinks_fUniqueID);
   fChain->SetBranchAddress("fKinks.fBits",fKinks_fBits);
   fChain->SetBranchAddress("fKinks.fID",fKinks_fID);
   fChain->SetBranchAddress("fKinks.fParamDaughter.fUniqueID",fKinks_fParamDaughter_fUniqueID);
   fChain->SetBranchAddress("fKinks.fParamDaughter.fBits",fKinks_fParamDaughter_fBits);
   fChain->SetBranchAddress("fKinks.fParamDaughter.fX",fKinks_fParamDaughter_fX);
   fChain->SetBranchAddress("fKinks.fParamDaughter.fAlpha",fKinks_fParamDaughter_fAlpha);
   fChain->SetBranchAddress("fKinks.fParamDaughter.fP[5]",fKinks_fParamDaughter_fP);
   fChain->SetBranchAddress("fKinks.fParamDaughter.fC[15]",fKinks_fParamDaughter_fC);
   fChain->SetBranchAddress("fKinks.fParamMother.fUniqueID",fKinks_fParamMother_fUniqueID);
   fChain->SetBranchAddress("fKinks.fParamMother.fBits",fKinks_fParamMother_fBits);
   fChain->SetBranchAddress("fKinks.fParamMother.fX",fKinks_fParamMother_fX);
   fChain->SetBranchAddress("fKinks.fParamMother.fAlpha",fKinks_fParamMother_fAlpha);
   fChain->SetBranchAddress("fKinks.fParamMother.fP[5]",fKinks_fParamMother_fP);
   fChain->SetBranchAddress("fKinks.fParamMother.fC[15]",fKinks_fParamMother_fC);
   fChain->SetBranchAddress("fKinks.fDist1",fKinks_fDist1);
   fChain->SetBranchAddress("fKinks.fDist2",fKinks_fDist2);
   fChain->SetBranchAddress("fKinks.fPdr[3]",fKinks_fPdr);
   fChain->SetBranchAddress("fKinks.fXr[3]",fKinks_fXr);
   fChain->SetBranchAddress("fKinks.fPm[3]",fKinks_fPm);
   fChain->SetBranchAddress("fKinks.fAngle[3]",fKinks_fAngle);
   fChain->SetBranchAddress("fKinks.fRr",fKinks_fRr);
   fChain->SetBranchAddress("fKinks.fLab[2]",fKinks_fLab);
   fChain->SetBranchAddress("fKinks.fIndex[2]",fKinks_fIndex);
   fChain->SetBranchAddress("fKinks.fStatus[12]",fKinks_fStatus);
   fChain->SetBranchAddress("fKinks.fTPCdensity[2][2]",fKinks_fTPCdensity);
   fChain->SetBranchAddress("fKinks.fTPCdensity2[2][2]",fKinks_fTPCdensity2);
   fChain->SetBranchAddress("fKinks.fShapeFactor",fKinks_fShapeFactor);
   fChain->SetBranchAddress("fKinks.fRow0",fKinks_fRow0);
   fChain->SetBranchAddress("fKinks.fMultiple[2]",fKinks_fMultiple);
   fChain->SetBranchAddress("fKinks.fTPCncls[2]",fKinks_fTPCncls);
   fChain->SetBranchAddress("fV0MIs",&fV0MIs_);
   fChain->SetBranchAddress("fV0MIs.fUniqueID",fV0MIs_fUniqueID);
   fChain->SetBranchAddress("fV0MIs.fBits",fV0MIs_fBits);
   fChain->SetBranchAddress("fV0MIs.fPdgCode",fV0MIs_fPdgCode);
   fChain->SetBranchAddress("fV0MIs.fEffMass",fV0MIs_fEffMass);
   fChain->SetBranchAddress("fV0MIs.fDcaDaughters",fV0MIs_fDcaDaughters);
   fChain->SetBranchAddress("fV0MIs.fChi2",fV0MIs_fChi2);
   fChain->SetBranchAddress("fV0MIs.fPos[3]",fV0MIs_fPos);
   fChain->SetBranchAddress("fV0MIs.fPosCov[6]",fV0MIs_fPosCov);
   fChain->SetBranchAddress("fV0MIs.fNidx",fV0MIs_fNidx);
   fChain->SetBranchAddress("fV0MIs.fNmom[3]",fV0MIs_fNmom);
   fChain->SetBranchAddress("fV0MIs.fNmomCov[6]",fV0MIs_fNmomCov);
   fChain->SetBranchAddress("fV0MIs.fPidx",fV0MIs_fPidx);
   fChain->SetBranchAddress("fV0MIs.fPmom[3]",fV0MIs_fPmom);
   fChain->SetBranchAddress("fV0MIs.fPmomCov[6]",fV0MIs_fPmomCov);
   fChain->SetBranchAddress("fV0MIs.fParamP.fUniqueID",fV0MIs_fParamP_fUniqueID);
   fChain->SetBranchAddress("fV0MIs.fParamP.fBits",fV0MIs_fParamP_fBits);
   fChain->SetBranchAddress("fV0MIs.fParamP.fX",fV0MIs_fParamP_fX);
   fChain->SetBranchAddress("fV0MIs.fParamP.fAlpha",fV0MIs_fParamP_fAlpha);
   fChain->SetBranchAddress("fV0MIs.fParamP.fP[5]",fV0MIs_fParamP_fP);
   fChain->SetBranchAddress("fV0MIs.fParamP.fC[15]",fV0MIs_fParamP_fC);
   fChain->SetBranchAddress("fV0MIs.fParamM.fUniqueID",fV0MIs_fParamM_fUniqueID);
   fChain->SetBranchAddress("fV0MIs.fParamM.fBits",fV0MIs_fParamM_fBits);
   fChain->SetBranchAddress("fV0MIs.fParamM.fX",fV0MIs_fParamM_fX);
   fChain->SetBranchAddress("fV0MIs.fParamM.fAlpha",fV0MIs_fParamM_fAlpha);
   fChain->SetBranchAddress("fV0MIs.fParamM.fP[5]",fV0MIs_fParamM_fP);
   fChain->SetBranchAddress("fV0MIs.fParamM.fC[15]",fV0MIs_fParamM_fC);
   fChain->SetBranchAddress("fV0MIs.fRP[5]",fV0MIs_fRP);
   fChain->SetBranchAddress("fV0MIs.fRM[5]",fV0MIs_fRM);
   fChain->SetBranchAddress("fV0MIs.fID",fV0MIs_fID);
   fChain->SetBranchAddress("fV0MIs.fLab[2]",fV0MIs_fLab);
   fChain->SetBranchAddress("fV0MIs.fIndex[2]",fV0MIs_fIndex);
   fChain->SetBranchAddress("fV0MIs.fNormDCAPrim[2]",fV0MIs_fNormDCAPrim);
   fChain->SetBranchAddress("fV0MIs.fDist1",fV0MIs_fDist1);
   fChain->SetBranchAddress("fV0MIs.fDist2",fV0MIs_fDist2);
   fChain->SetBranchAddress("fV0MIs.fPP[3]",fV0MIs_fPP);
   fChain->SetBranchAddress("fV0MIs.fPM[3]",fV0MIs_fPM);
   fChain->SetBranchAddress("fV0MIs.fXr[3]",fV0MIs_fXr);
   fChain->SetBranchAddress("fV0MIs.fAngle[3]",fV0MIs_fAngle);
   fChain->SetBranchAddress("fV0MIs.fRr",fV0MIs_fRr);
   fChain->SetBranchAddress("fV0MIs.fStatus",fV0MIs_fStatus);
   fChain->SetBranchAddress("fV0MIs.fRow0",fV0MIs_fRow0);
   fChain->SetBranchAddress("fV0MIs.fOrder[3]",fV0MIs_fOrder);
   fChain->SetBranchAddress("fV0MIs.fDistNorm",fV0MIs_fDistNorm);
   fChain->SetBranchAddress("fV0MIs.fDistSigma",fV0MIs_fDistSigma);
   fChain->SetBranchAddress("fV0MIs.fCausality[4]",fV0MIs_fCausality);
   fChain->SetBranchAddress("fV0MIs.fChi2Before",fV0MIs_fChi2Before);
   fChain->SetBranchAddress("fV0MIs.fNBefore",fV0MIs_fNBefore);
   fChain->SetBranchAddress("fV0MIs.fChi2After",fV0MIs_fChi2After);
   fChain->SetBranchAddress("fV0MIs.fNAfter",fV0MIs_fNAfter);
   fChain->SetBranchAddress("fV0MIs.fPointAngleFi",fV0MIs_fPointAngleFi);
   fChain->SetBranchAddress("fV0MIs.fPointAngleTh",fV0MIs_fPointAngleTh);
   fChain->SetBranchAddress("fV0MIs.fPointAngle",fV0MIs_fPointAngle);
   fChain->SetBranchAddress("fCaloClusters",&fCaloClusters_);
   fChain->SetBranchAddress("fCaloClusters.fUniqueID",fCaloClusters_fUniqueID);
   fChain->SetBranchAddress("fCaloClusters.fBits",fCaloClusters_fBits);
   fChain->SetBranchAddress("fCaloClusters.fID",fCaloClusters_fID);
   fChain->SetBranchAddress("fCaloClusters.fClusterType",fCaloClusters_fClusterType);
   fChain->SetBranchAddress("fCaloClusters.fEMCALCluster",fCaloClusters_fEMCALCluster);
   fChain->SetBranchAddress("fCaloClusters.fPHOSCluster",fCaloClusters_fPHOSCluster);
   fChain->SetBranchAddress("fCaloClusters.fGlobalPos[3]",fCaloClusters_fGlobalPos);
   fChain->SetBranchAddress("fCaloClusters.fEnergy",fCaloClusters_fEnergy);
   fChain->SetBranchAddress("fCaloClusters.fDispersion",fCaloClusters_fDispersion);
   fChain->SetBranchAddress("fCaloClusters.fChi2",fCaloClusters_fChi2);
   fChain->SetBranchAddress("fCaloClusters.fPID[10]",fCaloClusters_fPID);
   fChain->SetBranchAddress("fCaloClusters.fPrimaryIndex",fCaloClusters_fPrimaryIndex);
   fChain->SetBranchAddress("fCaloClusters.fM20",fCaloClusters_fM20);
   fChain->SetBranchAddress("fCaloClusters.fM02",fCaloClusters_fM02);
   fChain->SetBranchAddress("fCaloClusters.fM11",fCaloClusters_fM11);
   fChain->SetBranchAddress("fCaloClusters.fNExMax",fCaloClusters_fNExMax);
   fChain->SetBranchAddress("fCaloClusters.fEmcCpvDistance",fCaloClusters_fEmcCpvDistance);
   fChain->SetBranchAddress("fCaloClusters.fNumberOfDigits",fCaloClusters_fNumberOfDigits);
   fChain->SetBranchAddress("fCaloClusters.fDigitAmplitude",fCaloClusters_fDigitAmplitude);
   fChain->SetBranchAddress("fCaloClusters.fDigitTime",fCaloClusters_fDigitTime);
   fChain->SetBranchAddress("fCaloClusters.fDigitIndex",fCaloClusters_fDigitIndex);
   fChain->SetBranchAddress("fEMCALClusters",&fEMCALClusters);
   fChain->SetBranchAddress("fFirstEMCALCluster",&fFirstEMCALCluster);
   fChain->SetBranchAddress("fPHOSClusters",&fPHOSClusters);
   fChain->SetBranchAddress("fFirstPHOSCluster",&fFirstPHOSCluster);
}

Bool_t esdAna::Notify()
{
  static int nfiles=0;
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. Typically here the branch pointers
   // will be retrieved. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed.
  nfiles++;
  if (fChain) {
    if (fChain->GetCurrentFile()) {
      Info("Notify","Opening File: %s [%d]", fChain->GetCurrentFile()->GetName(),nfiles);
    }
  }

   // Get branch pointers
   b_ESD_fUniqueID = fChain->GetBranch("fUniqueID");
   b_ESD_fBits = fChain->GetBranch("fBits");
   b_ESD_fEventNumber = fChain->GetBranch("fEventNumber");
   b_ESD_fRunNumber = fChain->GetBranch("fRunNumber");
   b_ESD_fTriggerMask = fChain->GetBranch("fTriggerMask");
   b_ESD_fTriggerCluster = fChain->GetBranch("fTriggerCluster");
   b_ESD_fRecoVersion = fChain->GetBranch("fRecoVersion");
   b_ESD_fMagneticField = fChain->GetBranch("fMagneticField");
   b_ESD_fZDCN1Energy = fChain->GetBranch("fZDCN1Energy");
   b_ESD_fZDCP1Energy = fChain->GetBranch("fZDCP1Energy");
   b_ESD_fZDCN2Energy = fChain->GetBranch("fZDCN2Energy");
   b_ESD_fZDCP2Energy = fChain->GetBranch("fZDCP2Energy");
   b_ESD_fZDCEMEnergy = fChain->GetBranch("fZDCEMEnergy");
   b_ESD_fZDCParticipants = fChain->GetBranch("fZDCParticipants");
   b_ESD_fT0zVertex = fChain->GetBranch("fT0zVertex");
   b_ESD_fSPDVertex_fUniqueID = fChain->GetBranch("fSPDVertex.fUniqueID");
   b_ESD_fSPDVertex_fBits = fChain->GetBranch("fSPDVertex.fBits");
   b_ESD_fSPDVertex_fName = fChain->GetBranch("fSPDVertex.fName");
   b_ESD_fSPDVertex_fTitle = fChain->GetBranch("fSPDVertex.fTitle");
   b_ESD_fSPDVertex_fPosition = fChain->GetBranch("fSPDVertex.fPosition[3]");
   b_ESD_fSPDVertex_fSigma = fChain->GetBranch("fSPDVertex.fSigma");
   b_ESD_fSPDVertex_fNContributors = fChain->GetBranch("fSPDVertex.fNContributors");
   b_ESD_fSPDVertex_fCovXX = fChain->GetBranch("fSPDVertex.fCovXX");
   b_ESD_fSPDVertex_fCovXY = fChain->GetBranch("fSPDVertex.fCovXY");
   b_ESD_fSPDVertex_fCovYY = fChain->GetBranch("fSPDVertex.fCovYY");
   b_ESD_fSPDVertex_fCovXZ = fChain->GetBranch("fSPDVertex.fCovXZ");
   b_ESD_fSPDVertex_fCovYZ = fChain->GetBranch("fSPDVertex.fCovYZ");
   b_ESD_fSPDVertex_fCovZZ = fChain->GetBranch("fSPDVertex.fCovZZ");
   b_ESD_fSPDVertex_fSNR = fChain->GetBranch("fSPDVertex.fSNR[3]");
   b_ESD_fSPDVertex_fChi2 = fChain->GetBranch("fSPDVertex.fChi2");
   b_ESD_fSPDVertex_fTruePos = fChain->GetBranch("fSPDVertex.fTruePos[3]");
   b_ESD_fPrimaryVertex_fUniqueID = fChain->GetBranch("fPrimaryVertex.fUniqueID");
   b_ESD_fPrimaryVertex_fBits = fChain->GetBranch("fPrimaryVertex.fBits");
   b_ESD_fPrimaryVertex_fName = fChain->GetBranch("fPrimaryVertex.fName");
   b_ESD_fPrimaryVertex_fTitle = fChain->GetBranch("fPrimaryVertex.fTitle");
   b_ESD_fPrimaryVertex_fPosition = fChain->GetBranch("fPrimaryVertex.fPosition[3]");
   b_ESD_fPrimaryVertex_fSigma = fChain->GetBranch("fPrimaryVertex.fSigma");
   b_ESD_fPrimaryVertex_fNContributors = fChain->GetBranch("fPrimaryVertex.fNContributors");
   b_ESD_fPrimaryVertex_fCovXX = fChain->GetBranch("fPrimaryVertex.fCovXX");
   b_ESD_fPrimaryVertex_fCovXY = fChain->GetBranch("fPrimaryVertex.fCovXY");
   b_ESD_fPrimaryVertex_fCovYY = fChain->GetBranch("fPrimaryVertex.fCovYY");
   b_ESD_fPrimaryVertex_fCovXZ = fChain->GetBranch("fPrimaryVertex.fCovXZ");
   b_ESD_fPrimaryVertex_fCovYZ = fChain->GetBranch("fPrimaryVertex.fCovYZ");
   b_ESD_fPrimaryVertex_fCovZZ = fChain->GetBranch("fPrimaryVertex.fCovZZ");
   b_ESD_fPrimaryVertex_fSNR = fChain->GetBranch("fPrimaryVertex.fSNR[3]");
   b_ESD_fPrimaryVertex_fChi2 = fChain->GetBranch("fPrimaryVertex.fChi2");
   b_ESD_fPrimaryVertex_fTruePos = fChain->GetBranch("fPrimaryVertex.fTruePos[3]");
   b_ESD_fSPDMult_fUniqueID = fChain->GetBranch("fSPDMult.fUniqueID");
   b_ESD_fSPDMult_fBits = fChain->GetBranch("fSPDMult.fBits");
   b_ESD_fSPDMult_fNtracks = fChain->GetBranch("fSPDMult.fNtracks");
   b_fSPDMult_fTh = fChain->GetBranch("fSPDMult.fTh");
   b_fSPDMult_fPhi = fChain->GetBranch("fSPDMult.fPhi");
   b_fSPDMult_fDeltPhi = fChain->GetBranch("fSPDMult.fDeltPhi");
   b_ESD_fT0timeStart = fChain->GetBranch("fT0timeStart");
   b_ESD_fT0time = fChain->GetBranch("fT0time[24]");
   b_ESD_fT0amplitude = fChain->GetBranch("fT0amplitude[24]");
   b_ESD_fTracks_ = fChain->GetBranch("fTracks");
   b_fTracks_fUniqueID = fChain->GetBranch("fTracks.fUniqueID");
   b_fTracks_fBits = fChain->GetBranch("fTracks.fBits");
   b_fTracks_fX = fChain->GetBranch("fTracks.fX");
   b_fTracks_fAlpha = fChain->GetBranch("fTracks.fAlpha");
   b_fTracks_fP = fChain->GetBranch("fTracks.fP[5]");
   b_fTracks_fC = fChain->GetBranch("fTracks.fC[15]");
   b_fTracks_fFlags = fChain->GetBranch("fTracks.fFlags");
   b_fTracks_fLabel = fChain->GetBranch("fTracks.fLabel");
   b_fTracks_fID = fChain->GetBranch("fTracks.fID");
   b_fTracks_fTrackLength = fChain->GetBranch("fTracks.fTrackLength");
   b_fTracks_fD = fChain->GetBranch("fTracks.fD");
   b_fTracks_fZ = fChain->GetBranch("fTracks.fZ");
   b_fTracks_fCdd = fChain->GetBranch("fTracks.fCdd");
   b_fTracks_fCdz = fChain->GetBranch("fTracks.fCdz");
   b_fTracks_fCzz = fChain->GetBranch("fTracks.fCzz");
   b_fTracks_fTrackTime = fChain->GetBranch("fTracks.fTrackTime[5]");
   b_fTracks_fR = fChain->GetBranch("fTracks.fR[5]");
   b_fTracks_fStopVertex = fChain->GetBranch("fTracks.fStopVertex");
   b_fTracks_fCchi2 = fChain->GetBranch("fTracks.fCchi2");
   b_fTracks_fITSchi2 = fChain->GetBranch("fTracks.fITSchi2");
   b_fTracks_fITSncls = fChain->GetBranch("fTracks.fITSncls");
   b_fTracks_fITSsignal = fChain->GetBranch("fTracks.fITSsignal");
   b_fTracks_fITSr = fChain->GetBranch("fTracks.fITSr[5]");
   b_fTracks_fITSLabel = fChain->GetBranch("fTracks.fITSLabel");
   b_fTracks_fITSFakeRatio = fChain->GetBranch("fTracks.fITSFakeRatio");
   b_fTracks_fTPCchi2 = fChain->GetBranch("fTracks.fTPCchi2");
   b_fTracks_fTPCncls = fChain->GetBranch("fTracks.fTPCncls");
   b_fTracks_fTPCnclsF = fChain->GetBranch("fTracks.fTPCnclsF");
   b_fTracks_fTPCClusterMap_fUniqueID = fChain->GetBranch("fTracks.fTPCClusterMap.fUniqueID");
   b_fTracks_fTPCClusterMap_fBits = fChain->GetBranch("fTracks.fTPCClusterMap.fBits");
   b_fTracks_fTPCClusterMap_fNbits = fChain->GetBranch("fTracks.fTPCClusterMap.fNbits");
   b_fTracks_fTPCClusterMap_fNbytes = fChain->GetBranch("fTracks.fTPCClusterMap.fNbytes");
   b_fTracks_fTPCClusterMap_fAllBits = fChain->GetBranch("fTracks.fTPCClusterMap.fAllBits");
   b_fTracks_fTPCsignal = fChain->GetBranch("fTracks.fTPCsignal");
   b_fTracks_fTPCsignalN = fChain->GetBranch("fTracks.fTPCsignalN");
   b_fTracks_fTPCsignalS = fChain->GetBranch("fTracks.fTPCsignalS");
   b_fTracks_fTPCr = fChain->GetBranch("fTracks.fTPCr[5]");
   b_fTracks_fTPCLabel = fChain->GetBranch("fTracks.fTPCLabel");
   b_fTracks_fTPCPoints = fChain->GetBranch("fTracks.fTPCPoints[4]");
   b_fTracks_fKinkIndexes = fChain->GetBranch("fTracks.fKinkIndexes[3]");
   b_fTracks_fV0Indexes = fChain->GetBranch("fTracks.fV0Indexes[3]");
   b_fTracks_fTRDchi2 = fChain->GetBranch("fTracks.fTRDchi2");
   b_fTracks_fTRDncls = fChain->GetBranch("fTracks.fTRDncls");
   b_fTracks_fTRDncls0 = fChain->GetBranch("fTracks.fTRDncls0");
   b_fTracks_fTRDsignal = fChain->GetBranch("fTracks.fTRDsignal");
   b_fTracks_fTRDsignals = fChain->GetBranch("fTracks.fTRDsignals[6][3]");
   b_fTracks_fTRDTimBin = fChain->GetBranch("fTracks.fTRDTimBin[6]");
   b_fTracks_fTRDr = fChain->GetBranch("fTracks.fTRDr[5]");
   b_fTracks_fTRDLabel = fChain->GetBranch("fTracks.fTRDLabel");
   b_fTracks_fTRDQuality = fChain->GetBranch("fTracks.fTRDQuality");
   b_fTracks_fTRDBudget = fChain->GetBranch("fTracks.fTRDBudget");
   b_fTracks_fTOFchi2 = fChain->GetBranch("fTracks.fTOFchi2");
   b_fTracks_fTOFindex = fChain->GetBranch("fTracks.fTOFindex");
   b_fTracks_fTOFCalChannel = fChain->GetBranch("fTracks.fTOFCalChannel");
   b_fTracks_fTOFsignal = fChain->GetBranch("fTracks.fTOFsignal");
   b_fTracks_fTOFsignalToT = fChain->GetBranch("fTracks.fTOFsignalToT");
   b_fTracks_fTOFr = fChain->GetBranch("fTracks.fTOFr[5]");
   b_fTracks_fTOFLabel = fChain->GetBranch("fTracks.fTOFLabel[3]");
   b_fTracks_fHMPIDchi2 = fChain->GetBranch("fTracks.fHMPIDchi2");
   b_fTracks_fHMPIDncls = fChain->GetBranch("fTracks.fHMPIDncls");
   b_fTracks_fHMPIDindex = fChain->GetBranch("fTracks.fHMPIDindex");
   b_fTracks_fHMPIDsignal = fChain->GetBranch("fTracks.fHMPIDsignal");
   b_fTracks_fHMPIDr = fChain->GetBranch("fTracks.fHMPIDr[5]");
   b_fTracks_fHMPIDtheta = fChain->GetBranch("fTracks.fHMPIDtheta");
   b_fTracks_fHMPIDphi = fChain->GetBranch("fTracks.fHMPIDphi");
   b_fTracks_fHMPIDdx = fChain->GetBranch("fTracks.fHMPIDdx");
   b_fTracks_fHMPIDdy = fChain->GetBranch("fTracks.fHMPIDdy");
   b_fTracks_fHMPIDmipX = fChain->GetBranch("fTracks.fHMPIDmipX");
   b_fTracks_fHMPIDmipY = fChain->GetBranch("fTracks.fHMPIDmipY");
   b_ESD_fHLTConfMapTracks_ = fChain->GetBranch("fHLTConfMapTracks");
   b_fHLTConfMapTracks_fUniqueID = fChain->GetBranch("fHLTConfMapTracks.fUniqueID");
   b_fHLTConfMapTracks_fBits = fChain->GetBranch("fHLTConfMapTracks.fBits");
   b_fHLTConfMapTracks_fNHits = fChain->GetBranch("fHLTConfMapTracks.fNHits");
   b_fHLTConfMapTracks_fMCid = fChain->GetBranch("fHLTConfMapTracks.fMCid");
   b_fHLTConfMapTracks_fWeight = fChain->GetBranch("fHLTConfMapTracks.fWeight");
   b_fHLTConfMapTracks_fFromMainVertex = fChain->GetBranch("fHLTConfMapTracks.fFromMainVertex");
   b_fHLTConfMapTracks_fRowRange = fChain->GetBranch("fHLTConfMapTracks.fRowRange[2]");
   b_fHLTConfMapTracks_fSector = fChain->GetBranch("fHLTConfMapTracks.fSector");
   b_fHLTConfMapTracks_fFirstPoint = fChain->GetBranch("fHLTConfMapTracks.fFirstPoint[3]");
   b_fHLTConfMapTracks_fLastPoint = fChain->GetBranch("fHLTConfMapTracks.fLastPoint[3]");
   b_fHLTConfMapTracks_fQ = fChain->GetBranch("fHLTConfMapTracks.fQ");
   b_fHLTConfMapTracks_fTanl = fChain->GetBranch("fHLTConfMapTracks.fTanl");
   b_fHLTConfMapTracks_fPsi = fChain->GetBranch("fHLTConfMapTracks.fPsi");
   b_fHLTConfMapTracks_fPt = fChain->GetBranch("fHLTConfMapTracks.fPt");
   b_fHLTConfMapTracks_fPterr = fChain->GetBranch("fHLTConfMapTracks.fPterr");
   b_fHLTConfMapTracks_fPsierr = fChain->GetBranch("fHLTConfMapTracks.fPsierr");
   b_fHLTConfMapTracks_fTanlerr = fChain->GetBranch("fHLTConfMapTracks.fTanlerr");
   b_fHLTConfMapTracks_fBinX = fChain->GetBranch("fHLTConfMapTracks.fBinX");
   b_fHLTConfMapTracks_fBinY = fChain->GetBranch("fHLTConfMapTracks.fBinY");
   b_fHLTConfMapTracks_fSizeX = fChain->GetBranch("fHLTConfMapTracks.fSizeX");
   b_fHLTConfMapTracks_fSizeY = fChain->GetBranch("fHLTConfMapTracks.fSizeY");
   b_fHLTConfMapTracks_fPID = fChain->GetBranch("fHLTConfMapTracks.fPID");
   b_ESD_fHLTHoughTracks_ = fChain->GetBranch("fHLTHoughTracks");
   b_fHLTHoughTracks_fUniqueID = fChain->GetBranch("fHLTHoughTracks.fUniqueID");
   b_fHLTHoughTracks_fBits = fChain->GetBranch("fHLTHoughTracks.fBits");
   b_fHLTHoughTracks_fNHits = fChain->GetBranch("fHLTHoughTracks.fNHits");
   b_fHLTHoughTracks_fMCid = fChain->GetBranch("fHLTHoughTracks.fMCid");
   b_fHLTHoughTracks_fWeight = fChain->GetBranch("fHLTHoughTracks.fWeight");
   b_fHLTHoughTracks_fFromMainVertex = fChain->GetBranch("fHLTHoughTracks.fFromMainVertex");
   b_fHLTHoughTracks_fRowRange = fChain->GetBranch("fHLTHoughTracks.fRowRange[2]");
   b_fHLTHoughTracks_fSector = fChain->GetBranch("fHLTHoughTracks.fSector");
   b_fHLTHoughTracks_fFirstPoint = fChain->GetBranch("fHLTHoughTracks.fFirstPoint[3]");
   b_fHLTHoughTracks_fLastPoint = fChain->GetBranch("fHLTHoughTracks.fLastPoint[3]");
   b_fHLTHoughTracks_fQ = fChain->GetBranch("fHLTHoughTracks.fQ");
   b_fHLTHoughTracks_fTanl = fChain->GetBranch("fHLTHoughTracks.fTanl");
   b_fHLTHoughTracks_fPsi = fChain->GetBranch("fHLTHoughTracks.fPsi");
   b_fHLTHoughTracks_fPt = fChain->GetBranch("fHLTHoughTracks.fPt");
   b_fHLTHoughTracks_fPterr = fChain->GetBranch("fHLTHoughTracks.fPterr");
   b_fHLTHoughTracks_fPsierr = fChain->GetBranch("fHLTHoughTracks.fPsierr");
   b_fHLTHoughTracks_fTanlerr = fChain->GetBranch("fHLTHoughTracks.fTanlerr");
   b_fHLTHoughTracks_fBinX = fChain->GetBranch("fHLTHoughTracks.fBinX");
   b_fHLTHoughTracks_fBinY = fChain->GetBranch("fHLTHoughTracks.fBinY");
   b_fHLTHoughTracks_fSizeX = fChain->GetBranch("fHLTHoughTracks.fSizeX");
   b_fHLTHoughTracks_fSizeY = fChain->GetBranch("fHLTHoughTracks.fSizeY");
   b_fHLTHoughTracks_fPID = fChain->GetBranch("fHLTHoughTracks.fPID");
   b_ESD_fMuonTracks_ = fChain->GetBranch("fMuonTracks");
   b_fMuonTracks_fUniqueID = fChain->GetBranch("fMuonTracks.fUniqueID");
   b_fMuonTracks_fBits = fChain->GetBranch("fMuonTracks.fBits");
   b_fMuonTracks_fInverseBendingMomentum = fChain->GetBranch("fMuonTracks.fInverseBendingMomentum");
   b_fMuonTracks_fThetaX = fChain->GetBranch("fMuonTracks.fThetaX");
   b_fMuonTracks_fThetaY = fChain->GetBranch("fMuonTracks.fThetaY");
   b_fMuonTracks_fZ = fChain->GetBranch("fMuonTracks.fZ");
   b_fMuonTracks_fBendingCoor = fChain->GetBranch("fMuonTracks.fBendingCoor");
   b_fMuonTracks_fNonBendingCoor = fChain->GetBranch("fMuonTracks.fNonBendingCoor");
   b_fMuonTracks_fChi2 = fChain->GetBranch("fMuonTracks.fChi2");
   b_fMuonTracks_fNHit = fChain->GetBranch("fMuonTracks.fNHit");
   b_fMuonTracks_fMatchTrigger = fChain->GetBranch("fMuonTracks.fMatchTrigger");
   b_fMuonTracks_fChi2MatchTrigger = fChain->GetBranch("fMuonTracks.fChi2MatchTrigger");
   b_ESD_fPmdTracks_ = fChain->GetBranch("fPmdTracks");
   b_fPmdTracks_fUniqueID = fChain->GetBranch("fPmdTracks.fUniqueID");
   b_fPmdTracks_fBits = fChain->GetBranch("fPmdTracks.fBits");
   b_fPmdTracks_fDet = fChain->GetBranch("fPmdTracks.fDet");
   b_fPmdTracks_fX = fChain->GetBranch("fPmdTracks.fX");
   b_fPmdTracks_fY = fChain->GetBranch("fPmdTracks.fY");
   b_fPmdTracks_fZ = fChain->GetBranch("fPmdTracks.fZ");
   b_fPmdTracks_fCluADC = fChain->GetBranch("fPmdTracks.fCluADC");
   b_fPmdTracks_fNcell = fChain->GetBranch("fPmdTracks.fNcell");
   b_fPmdTracks_fCluPID = fChain->GetBranch("fPmdTracks.fCluPID");
   b_ESD_fTrdTracks_ = fChain->GetBranch("fTrdTracks");
   b_fTrdTracks_fUniqueID = fChain->GetBranch("fTrdTracks.fUniqueID");
   b_fTrdTracks_fBits = fChain->GetBranch("fTrdTracks.fBits");
   b_fTrdTracks_fYproj = fChain->GetBranch("fTrdTracks.fYproj");
   b_fTrdTracks_fZproj = fChain->GetBranch("fTrdTracks.fZproj");
   b_fTrdTracks_fSlope = fChain->GetBranch("fTrdTracks.fSlope");
   b_fTrdTracks_fDetector = fChain->GetBranch("fTrdTracks.fDetector");
   b_fTrdTracks_fNtracklets = fChain->GetBranch("fTrdTracks.fNtracklets");
   b_fTrdTracks_fNplanes = fChain->GetBranch("fTrdTracks.fNplanes");
   b_fTrdTracks_fNclusters = fChain->GetBranch("fTrdTracks.fNclusters");
   b_fTrdTracks_fPt = fChain->GetBranch("fTrdTracks.fPt");
   b_fTrdTracks_fPhi = fChain->GetBranch("fTrdTracks.fPhi");
   b_fTrdTracks_fEta = fChain->GetBranch("fTrdTracks.fEta");
   b_fTrdTracks_fLabel = fChain->GetBranch("fTrdTracks.fLabel");
   b_fTrdTracks_fPID = fChain->GetBranch("fTrdTracks.fPID");
   b_fTrdTracks_fIsElectron = fChain->GetBranch("fTrdTracks.fIsElectron");
   b_ESD_fV0s_ = fChain->GetBranch("fV0s");
   b_fV0s_fUniqueID = fChain->GetBranch("fV0s.fUniqueID");
   b_fV0s_fBits = fChain->GetBranch("fV0s.fBits");
   b_fV0s_fPdgCode = fChain->GetBranch("fV0s.fPdgCode");
   b_fV0s_fEffMass = fChain->GetBranch("fV0s.fEffMass");
   b_fV0s_fDcaDaughters = fChain->GetBranch("fV0s.fDcaDaughters");
   b_fV0s_fChi2 = fChain->GetBranch("fV0s.fChi2");
   b_fV0s_fPos = fChain->GetBranch("fV0s.fPos[3]");
   b_fV0s_fPosCov = fChain->GetBranch("fV0s.fPosCov[6]");
   b_fV0s_fNidx = fChain->GetBranch("fV0s.fNidx");
   b_fV0s_fNmom = fChain->GetBranch("fV0s.fNmom[3]");
   b_fV0s_fNmomCov = fChain->GetBranch("fV0s.fNmomCov[6]");
   b_fV0s_fPidx = fChain->GetBranch("fV0s.fPidx");
   b_fV0s_fPmom = fChain->GetBranch("fV0s.fPmom[3]");
   b_fV0s_fPmomCov = fChain->GetBranch("fV0s.fPmomCov[6]");
   b_ESD_fCascades_ = fChain->GetBranch("fCascades");
   b_fCascades_fUniqueID = fChain->GetBranch("fCascades.fUniqueID");
   b_fCascades_fBits = fChain->GetBranch("fCascades.fBits");
   b_fCascades_fPdgCode = fChain->GetBranch("fCascades.fPdgCode");
   b_fCascades_fEffMass = fChain->GetBranch("fCascades.fEffMass");
   b_fCascades_fChi2 = fChain->GetBranch("fCascades.fChi2");
   b_fCascades_fPos = fChain->GetBranch("fCascades.fPos[3]");
   b_fCascades_fPosCov = fChain->GetBranch("fCascades.fPosCov[6]");
   b_fCascades_fV0idx = fChain->GetBranch("fCascades.fV0idx[2]");
   b_fCascades_fV0mom = fChain->GetBranch("fCascades.fV0mom[2][3]");
   b_fCascades_fV0momCov = fChain->GetBranch("fCascades.fV0momCov[6]");
   b_fCascades_fBachIdx = fChain->GetBranch("fCascades.fBachIdx");
   b_fCascades_fBachMom = fChain->GetBranch("fCascades.fBachMom[3]");
   b_fCascades_fBachMomCov = fChain->GetBranch("fCascades.fBachMomCov[6]");
   b_ESD_fKinks_ = fChain->GetBranch("fKinks");
   b_fKinks_fUniqueID = fChain->GetBranch("fKinks.fUniqueID");
   b_fKinks_fBits = fChain->GetBranch("fKinks.fBits");
   b_fKinks_fID = fChain->GetBranch("fKinks.fID");
   b_fKinks_fParamDaughter_fUniqueID = fChain->GetBranch("fKinks.fParamDaughter.fUniqueID");
   b_fKinks_fParamDaughter_fBits = fChain->GetBranch("fKinks.fParamDaughter.fBits");
   b_fKinks_fParamDaughter_fX = fChain->GetBranch("fKinks.fParamDaughter.fX");
   b_fKinks_fParamDaughter_fAlpha = fChain->GetBranch("fKinks.fParamDaughter.fAlpha");
   b_fKinks_fParamDaughter_fP = fChain->GetBranch("fKinks.fParamDaughter.fP[5]");
   b_fKinks_fParamDaughter_fC = fChain->GetBranch("fKinks.fParamDaughter.fC[15]");
   b_fKinks_fParamMother_fUniqueID = fChain->GetBranch("fKinks.fParamMother.fUniqueID");
   b_fKinks_fParamMother_fBits = fChain->GetBranch("fKinks.fParamMother.fBits");
   b_fKinks_fParamMother_fX = fChain->GetBranch("fKinks.fParamMother.fX");
   b_fKinks_fParamMother_fAlpha = fChain->GetBranch("fKinks.fParamMother.fAlpha");
   b_fKinks_fParamMother_fP = fChain->GetBranch("fKinks.fParamMother.fP[5]");
   b_fKinks_fParamMother_fC = fChain->GetBranch("fKinks.fParamMother.fC[15]");
   b_fKinks_fDist1 = fChain->GetBranch("fKinks.fDist1");
   b_fKinks_fDist2 = fChain->GetBranch("fKinks.fDist2");
   b_fKinks_fPdr = fChain->GetBranch("fKinks.fPdr[3]");
   b_fKinks_fXr = fChain->GetBranch("fKinks.fXr[3]");
   b_fKinks_fPm = fChain->GetBranch("fKinks.fPm[3]");
   b_fKinks_fAngle = fChain->GetBranch("fKinks.fAngle[3]");
   b_fKinks_fRr = fChain->GetBranch("fKinks.fRr");
   b_fKinks_fLab = fChain->GetBranch("fKinks.fLab[2]");
   b_fKinks_fIndex = fChain->GetBranch("fKinks.fIndex[2]");
   b_fKinks_fStatus = fChain->GetBranch("fKinks.fStatus[12]");
   b_fKinks_fTPCdensity = fChain->GetBranch("fKinks.fTPCdensity[2][2]");
   b_fKinks_fTPCdensity2 = fChain->GetBranch("fKinks.fTPCdensity2[2][2]");
   b_fKinks_fShapeFactor = fChain->GetBranch("fKinks.fShapeFactor");
   b_fKinks_fRow0 = fChain->GetBranch("fKinks.fRow0");
   b_fKinks_fMultiple = fChain->GetBranch("fKinks.fMultiple[2]");
   b_fKinks_fTPCncls = fChain->GetBranch("fKinks.fTPCncls[2]");
   b_ESD_fV0MIs_ = fChain->GetBranch("fV0MIs");
   b_fV0MIs_fUniqueID = fChain->GetBranch("fV0MIs.fUniqueID");
   b_fV0MIs_fBits = fChain->GetBranch("fV0MIs.fBits");
   b_fV0MIs_fPdgCode = fChain->GetBranch("fV0MIs.fPdgCode");
   b_fV0MIs_fEffMass = fChain->GetBranch("fV0MIs.fEffMass");
   b_fV0MIs_fDcaDaughters = fChain->GetBranch("fV0MIs.fDcaDaughters");
   b_fV0MIs_fChi2 = fChain->GetBranch("fV0MIs.fChi2");
   b_fV0MIs_fPos = fChain->GetBranch("fV0MIs.fPos[3]");
   b_fV0MIs_fPosCov = fChain->GetBranch("fV0MIs.fPosCov[6]");
   b_fV0MIs_fNidx = fChain->GetBranch("fV0MIs.fNidx");
   b_fV0MIs_fNmom = fChain->GetBranch("fV0MIs.fNmom[3]");
   b_fV0MIs_fNmomCov = fChain->GetBranch("fV0MIs.fNmomCov[6]");
   b_fV0MIs_fPidx = fChain->GetBranch("fV0MIs.fPidx");
   b_fV0MIs_fPmom = fChain->GetBranch("fV0MIs.fPmom[3]");
   b_fV0MIs_fPmomCov = fChain->GetBranch("fV0MIs.fPmomCov[6]");
   b_fV0MIs_fParamP_fUniqueID = fChain->GetBranch("fV0MIs.fParamP.fUniqueID");
   b_fV0MIs_fParamP_fBits = fChain->GetBranch("fV0MIs.fParamP.fBits");
   b_fV0MIs_fParamP_fX = fChain->GetBranch("fV0MIs.fParamP.fX");
   b_fV0MIs_fParamP_fAlpha = fChain->GetBranch("fV0MIs.fParamP.fAlpha");
   b_fV0MIs_fParamP_fP = fChain->GetBranch("fV0MIs.fParamP.fP[5]");
   b_fV0MIs_fParamP_fC = fChain->GetBranch("fV0MIs.fParamP.fC[15]");
   b_fV0MIs_fParamM_fUniqueID = fChain->GetBranch("fV0MIs.fParamM.fUniqueID");
   b_fV0MIs_fParamM_fBits = fChain->GetBranch("fV0MIs.fParamM.fBits");
   b_fV0MIs_fParamM_fX = fChain->GetBranch("fV0MIs.fParamM.fX");
   b_fV0MIs_fParamM_fAlpha = fChain->GetBranch("fV0MIs.fParamM.fAlpha");
   b_fV0MIs_fParamM_fP = fChain->GetBranch("fV0MIs.fParamM.fP[5]");
   b_fV0MIs_fParamM_fC = fChain->GetBranch("fV0MIs.fParamM.fC[15]");
   b_fV0MIs_fRP = fChain->GetBranch("fV0MIs.fRP[5]");
   b_fV0MIs_fRM = fChain->GetBranch("fV0MIs.fRM[5]");
   b_fV0MIs_fID = fChain->GetBranch("fV0MIs.fID");
   b_fV0MIs_fLab = fChain->GetBranch("fV0MIs.fLab[2]");
   b_fV0MIs_fIndex = fChain->GetBranch("fV0MIs.fIndex[2]");
   b_fV0MIs_fNormDCAPrim = fChain->GetBranch("fV0MIs.fNormDCAPrim[2]");
   b_fV0MIs_fDist1 = fChain->GetBranch("fV0MIs.fDist1");
   b_fV0MIs_fDist2 = fChain->GetBranch("fV0MIs.fDist2");
   b_fV0MIs_fPP = fChain->GetBranch("fV0MIs.fPP[3]");
   b_fV0MIs_fPM = fChain->GetBranch("fV0MIs.fPM[3]");
   b_fV0MIs_fXr = fChain->GetBranch("fV0MIs.fXr[3]");
   b_fV0MIs_fAngle = fChain->GetBranch("fV0MIs.fAngle[3]");
   b_fV0MIs_fRr = fChain->GetBranch("fV0MIs.fRr");
   b_fV0MIs_fStatus = fChain->GetBranch("fV0MIs.fStatus");
   b_fV0MIs_fRow0 = fChain->GetBranch("fV0MIs.fRow0");
   b_fV0MIs_fOrder = fChain->GetBranch("fV0MIs.fOrder[3]");
   b_fV0MIs_fDistNorm = fChain->GetBranch("fV0MIs.fDistNorm");
   b_fV0MIs_fDistSigma = fChain->GetBranch("fV0MIs.fDistSigma");
   b_fV0MIs_fCausality = fChain->GetBranch("fV0MIs.fCausality[4]");
   b_fV0MIs_fChi2Before = fChain->GetBranch("fV0MIs.fChi2Before");
   b_fV0MIs_fNBefore = fChain->GetBranch("fV0MIs.fNBefore");
   b_fV0MIs_fChi2After = fChain->GetBranch("fV0MIs.fChi2After");
   b_fV0MIs_fNAfter = fChain->GetBranch("fV0MIs.fNAfter");
   b_fV0MIs_fPointAngleFi = fChain->GetBranch("fV0MIs.fPointAngleFi");
   b_fV0MIs_fPointAngleTh = fChain->GetBranch("fV0MIs.fPointAngleTh");
   b_fV0MIs_fPointAngle = fChain->GetBranch("fV0MIs.fPointAngle");
   b_ESD_fCaloClusters_ = fChain->GetBranch("fCaloClusters");
   b_fCaloClusters_fUniqueID = fChain->GetBranch("fCaloClusters.fUniqueID");
   b_fCaloClusters_fBits = fChain->GetBranch("fCaloClusters.fBits");
   b_fCaloClusters_fID = fChain->GetBranch("fCaloClusters.fID");
   b_fCaloClusters_fClusterType = fChain->GetBranch("fCaloClusters.fClusterType");
   b_fCaloClusters_fEMCALCluster = fChain->GetBranch("fCaloClusters.fEMCALCluster");
   b_fCaloClusters_fPHOSCluster = fChain->GetBranch("fCaloClusters.fPHOSCluster");
   b_fCaloClusters_fGlobalPos = fChain->GetBranch("fCaloClusters.fGlobalPos[3]");
   b_fCaloClusters_fEnergy = fChain->GetBranch("fCaloClusters.fEnergy");
   b_fCaloClusters_fDispersion = fChain->GetBranch("fCaloClusters.fDispersion");
   b_fCaloClusters_fChi2 = fChain->GetBranch("fCaloClusters.fChi2");
   b_fCaloClusters_fPID = fChain->GetBranch("fCaloClusters.fPID[10]");
   b_fCaloClusters_fPrimaryIndex = fChain->GetBranch("fCaloClusters.fPrimaryIndex");
   b_fCaloClusters_fM20 = fChain->GetBranch("fCaloClusters.fM20");
   b_fCaloClusters_fM02 = fChain->GetBranch("fCaloClusters.fM02");
   b_fCaloClusters_fM11 = fChain->GetBranch("fCaloClusters.fM11");
   b_fCaloClusters_fNExMax = fChain->GetBranch("fCaloClusters.fNExMax");
   b_fCaloClusters_fEmcCpvDistance = fChain->GetBranch("fCaloClusters.fEmcCpvDistance");
   b_fCaloClusters_fNumberOfDigits = fChain->GetBranch("fCaloClusters.fNumberOfDigits");
   b_fCaloClusters_fDigitAmplitude = fChain->GetBranch("fCaloClusters.fDigitAmplitude");
   b_fCaloClusters_fDigitTime = fChain->GetBranch("fCaloClusters.fDigitTime");
   b_fCaloClusters_fDigitIndex = fChain->GetBranch("fCaloClusters.fDigitIndex");
   b_ESD_fEMCALClusters = fChain->GetBranch("fEMCALClusters");
   b_ESD_fFirstEMCALCluster = fChain->GetBranch("fFirstEMCALCluster");
   b_ESD_fPHOSClusters = fChain->GetBranch("fPHOSClusters");
   b_ESD_fFirstPHOSCluster = fChain->GetBranch("fFirstPHOSCluster");

   return kTRUE;
}

#endif // #ifdef esdAna_cxx
