/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* mailto: svallero@to.infn.it, ycorrale@cern.ch */

/* Class defining containers for the HF b-jets analysis */

//--Root--
#include <TClonesArray.h>

//--AliRoot--
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliEmcalJet.h"

//--AliHFJetsClass--
#include "AliHFJetsContainerVertex.h"

//_____________________________________________________________________________________

ClassImp(AliHFJetsContainerVertex)

//_____________________________________________________________________________________
AliHFJetsContainerVertex::AliHFJetsContainerVertex() : AliHFJetsContainer("", kTRUE),
fType(kJetVtxSim)
{
  fTagger = new AliHFJetsTaggingVertex("tagger");
  // dummy
}

//_____________________________________________________________________________________
AliHFJetsContainerVertex::AliHFJetsContainerVertex(const char *name, EContTYPE contType) :
AliHFJetsContainer(name, kTRUE),
fType(kJetVtxSim)
{
  // Constructor
  fTagger = new AliHFJetsTaggingVertex("tagger");
  CreateContainerVertex(contType);
}

//_____________________________________________________________________________________
AliHFJetsContainerVertex::AliHFJetsContainerVertex(const AliHFJetsContainerVertex &c) :
AliHFJetsContainer("", kTRUE)
{
  // AliHFJetsContainerVertex copy constructor

  ((AliHFJetsContainerVertex &)c).Copy(*this);
}

//_____________________________________________________________________________________
AliHFJetsContainerVertex::~AliHFJetsContainerVertex() {

  // Destructor
  if (fTagger) {
    delete fTagger; fTagger = NULL;
  }
}

//_____________________________________________________________________________________
AliHFJetsContainerVertex &AliHFJetsContainerVertex::operator=(const AliHFJetsContainerVertex &c)
{
  // assigment operator
  if (this != &c)
    ((AliHFJetsContainerVertex &) c).Copy(*this);

  return *this;
}

//_____________________________________________________________________________________
void AliHFJetsContainerVertex::FillStepJetVtxSim(CFSteps                step,
                                                 const Int_t            nSVtx,
                                                 Double_t               evtCent,          // Event multiplicity percentile (based on ZNA in pA)
                                                 Double_t               jetPt_wBkgRej,    // Jet's pt after background subtraction
                                                 vctr_pair_dbl_int     &arrVtxDisp,       // Vector of pair with SV_vtx and SV_sigma
                                                 const TClonesArray    *arrVtxHF,
                                                 const AliAODVertex    *primVtx,
                                                 const AliEmcalJet     *jet,
                                                 const TClonesArray    *mcPart,
                                                 const Double_t        *partonnat,
                                                 const Double_t        *partpt,
                                                 Double_t               wght)
{
  if (fType != kJetVtxSim) {
    AliErrorF(MSGERROR("This method is available only for container type kJetVtxSim: \
                       you are trying to fill %s!"), fgStrContType(fType));
  }

  const Int_t kNvar = 23;
  Double_t point[kNvar] = {
    evtCent,                            /* 0 */
    jet->Pt(),                          /* 1 */
    jet->Eta(),                         /* 2 */
    jet->Phi()-TMath::Pi(),             /* 3 */
    jet->Area(),                        /* 4 */
    jetPt_wBkgRej,                      /* 5 */
    static_cast<Double_t>(nSVtx),       /* 6 */
    -1, -1, -1, -1, -1, -1, -1, -1,     /* 7-14 */   // mass1, sigVtx1, Lxy1, SLxy1, SIP1, SIP2, SIP3, nV0Trks
    partonnat[0],                       /* 15 */
    partonnat[1],                       /* 16 */
    partpt[0],                          /* 17 */
    partpt[1],                          /* 18 */
    -99,-99,-99,                        /* 19-21 */   // nReal1, nFromB1, nFromPrD1
    -99                                 /* 22 */      // sumMass3MostDispl
  };
  if (step == kCFStepEventSelected) {

    TArrayD *apoint = new TArrayD(kNvar, point);
    FillStep(step, apoint, wght);
    
    return;
  }

  if (arrVtxHF && primVtx) {
    
    Double_t xyz[3], vtxVect[3], jetP[3];
    Double_t xyzPrim[3];
    
    primVtx->GetXYZ(xyzPrim);
    jet->PxPyPz(jetP);
    
    Int_t    *idxLxy = new Int_t[nSVtx];
    Int_t    *nTrkV0 = new Int_t[nSVtx];
    
    Double_t *decLenXY        = new Double_t[nSVtx];
    Double_t *errdecLenXY     = new Double_t[nSVtx];
    Double_t *sigdecLenXY     = new Double_t[nSVtx];
    Double_t *invMasses       = new Double_t[nSVtx];
    Double_t *nRealVtx        = new Double_t[nSVtx];
    Double_t *nFromBVtx       = new Double_t[nSVtx];
    Double_t *nFromPromptDVtx = new Double_t[nSVtx];
    Double_t *sigmavertex     = new Double_t[nSVtx];
    
    Float_t *ipR1 = new Float_t[nSVtx];
    Float_t *ipR2 = new Float_t[nSVtx];
    Float_t *ipR3 = new Float_t[nSVtx];
    
    Double_t *sigipR1 = new Double_t[nSVtx];
    Double_t *sigipR2 = new Double_t[nSVtx];
    Double_t *sigipR3 = new Double_t[nSVtx];
    
    Double_t xMC,yMC;
    Double_t signLxy;
    
    Int_t nvtxMC = 0;
    
    for (Int_t jj = 0; jj < nSVtx; ++jj) {
      
      xMC = -99999.;
      yMC = -99999.;
      
      AliAODVertex *svtx = (AliAODVertex *)arrVtxHF->UncheckedAt(jj);
      invMasses[jj] = fTagger->GetVertexInvariantMass(svtx);
      nRealVtx[jj]  = -1;
      nFromBVtx[jj] = -1;
      
      svtx->GetXYZ(xyz);
      vtxVect[0] = xyz[0]-xyzPrim[0];
      vtxVect[1] = xyz[1]-xyzPrim[1];
      vtxVect[2] = xyz[2]-xyzPrim[2];
      
      signLxy = vtxVect[0] * jetP[0] + vtxVect[1] * jetP[1];
      
      sigmavertex[jj] = arrVtxDisp[jj].first;
      nTrkV0[jj]      = arrVtxDisp[jj].second;
      
      decLenXY[jj] = primVtx->DistanceXYToVertex(svtx);
      
      if (signLxy < 0.) decLenXY[jj] *= -1.;
      
      errdecLenXY[jj] = primVtx->ErrorDistanceXYToVertex(svtx);
      sigdecLenXY[jj] = decLenXY[jj]/errdecLenXY[jj];
      
      Double_t dca0aod[2] = {-999., -999.};
      Double_t dca1aod[2] = {-999., -999.};
      Double_t dca2aod[2] = {-999., -999.};
      Double_t cov0aod[3] = {-999., -999., -999.};
      Double_t cov1aod[3] = {-999., -999., -999.};
      Double_t cov2aod[3] = {-999., -999., -999.};
      
      AliAODTrack *track0 = (AliAODTrack *)svtx->GetDaughter(0);
      AliAODTrack *track1 = (AliAODTrack *)svtx->GetDaughter(1);
      AliAODTrack *track2 = (AliAODTrack *)svtx->GetDaughter(2);

      track0->PropagateToDCA(primVtx, 0., 10000., dca0aod, cov0aod);
      track1->PropagateToDCA(primVtx, 0., 10000., dca1aod, cov1aod);
      track2->PropagateToDCA(primVtx, 0., 10000., dca2aod, cov2aod);
      
      ipR1[jj] = dca0aod[0];
      ipR2[jj] = dca1aod[0];
      ipR3[jj] = dca2aod[0];
      
      if (cov0aod[0]) sigipR1[jj] = dca0aod[0]/TMath::Sqrt(cov0aod[0]);;
      if (cov1aod[0]) sigipR2[jj] = dca1aod[0]/TMath::Sqrt(cov1aod[0]);
      if (cov2aod[0]) sigipR3[jj] = dca2aod[0]/TMath::Sqrt(cov2aod[0]);
      
      if (mcPart) {
        
        Int_t nfromBandD   = 0;
        Int_t nfromD       = 0;
        Int_t nfromPromptD = 0;
        fTagger->GetNTracksFromCommonVertex(svtx, mcPart, nvtxMC, xMC, yMC,
                                            nfromBandD, nfromD, nfromPromptD);
        
        nRealVtx[jj]        = nvtxMC;
        nFromBVtx[jj]       = nfromBandD;
        nFromPromptDVtx[jj] = nfromPromptD;
      }
    }
    
    TMath::Sort(nSVtx, decLenXY, idxLxy);
    if (nSVtx > 0) {
      
      point[7]  = invMasses[idxLxy[0]];
      point[8]  = sigmavertex[idxLxy[0]];
      point[9]  = decLenXY[idxLxy[0]];
      point[10] = sigdecLenXY[idxLxy[0]];
      point[11] = sigipR1[idxLxy[0]];
      point[12] = sigipR2[idxLxy[0]];
      point[13] = sigipR3[idxLxy[0]];
      point[14] = nTrkV0[idxLxy[0]];
      
      point[19]  = nRealVtx[idxLxy[0]];
      point[20]  = nFromBVtx[idxLxy[0]];
      point[21]  = nFromPromptDVtx[idxLxy[0]];
    }
    // Calculate sum of inv masses of the 3 more displaced vertices
    for (Int_t ivtx = 0; ivtx < 3; ++ivtx) {
      if (nSVtx > ivtx) point[22] += invMasses[idxLxy[ivtx]];
    }
    
    delete[] idxLxy;
    delete[] nTrkV0;
    
    delete[] decLenXY;
    delete[] invMasses;
    delete[] errdecLenXY;
    delete[] sigdecLenXY;
    delete[] nRealVtx;
    delete[] nFromBVtx;
    delete[] nFromPromptDVtx;
    delete[] sigmavertex;
    
    delete[] ipR1;
    delete[] ipR2;
    delete[] ipR3;
    
    delete[] sigipR1;
    delete[] sigipR2;
    delete[] sigipR3;
    
  } // end if (vertices && primVtx)
  
  TArrayD *apoint = new TArrayD(kNvar, point);
  FillStep(step, apoint, wght);
  
  return;
}

//_____________________________________________________________________________________
void AliHFJetsContainerVertex::FillStepJetVtxData(CFSteps               step,
                                                  const Int_t           nSVtx,
                                                  Double_t              evtCent,
                                                  Double_t              jetPt_wBkgRej,
                                                  vctr_pair_dbl_int    &arrVtxDisp,
                                                  const TClonesArray   *arrVtxHF,
                                                  const AliAODVertex   *primVtx,
                                                  const AliEmcalJet    *jet,
                                                  Double_t              wght)
{
  if (fType != kJetVtxData) {
    AliErrorF(MSGERROR("This method is available only for container type kJetVtxData: \
                       you are trying to fill %s!"), fgStrContType(fType));
  }
  
  const Int_t kNvar = 15;
  Double_t point[kNvar] = {
    evtCent,                            /* 0 */
    jet->Pt(),                          /* 1 */
    jet->Eta(),                         /* 2 */
    jet->Phi()-TMath::Pi(),             /* 3 */
    jet->Area(),                        /* 4 */
    jetPt_wBkgRej,                      /* 5 */
    static_cast<Double_t>(nSVtx),       /* 6 */
    -1, -1, -1, -1, -1, -1, -1, -1,     /* 7-14 */   // mass1, sigVtx1, Lxy1, SLxy1, SIP1, SIP2, SIP3, nV0Trks
  };
  
  if (arrVtxHF && primVtx) {
    
    Double_t xyz[3], vtxVect[3], jetP[3];
    Double_t xyzPrim[3];
    
    primVtx->GetXYZ(xyzPrim);
    jet->PxPyPz(jetP);
    
    Int_t    *idxLxy = new Int_t[nSVtx];
    Int_t    *nTrkV0 = new Int_t[nSVtx];
    
    Double_t *decLenXY        = new Double_t[nSVtx];
    Double_t *errdecLenXY     = new Double_t[nSVtx];
    Double_t *sigdecLenXY     = new Double_t[nSVtx];
    Double_t *invMasses       = new Double_t[nSVtx];
    Double_t *sigmavertex     = new Double_t[nSVtx];
    
    Float_t *ipR1 = new Float_t[nSVtx];
    Float_t *ipR2 = new Float_t[nSVtx];
    Float_t *ipR3 = new Float_t[nSVtx];
    
    Double_t *sigipR1 = new Double_t[nSVtx];
    Double_t *sigipR2 = new Double_t[nSVtx];
    Double_t *sigipR3 = new Double_t[nSVtx];
    
    Double_t *sumOfSqs    = new Double_t[nSVtx];
    Double_t *sigsumOfSqs = new Double_t[nSVtx];
    
    Double_t vtxP[3], signLxy;
    
    for (Int_t jj = 0; jj < nSVtx; ++jj) {
      
      AliAODVertex *svtx = (AliAODVertex *)arrVtxHF->UncheckedAt(jj);
      invMasses[jj] = fTagger->GetVertexInvariantMass(svtx);
      
      svtx->GetXYZ(xyz);
      vtxVect[0] = xyz[0]-xyzPrim[0];
      vtxVect[1] = xyz[1]-xyzPrim[1];
      vtxVect[2] = xyz[2]-xyzPrim[2];
      
      signLxy = vtxVect[0] * jetP[0] + vtxVect[1] * jetP[1];
      
      sigmavertex[jj] = arrVtxDisp[jj].first;
      nTrkV0[jj]      = arrVtxDisp[jj].second;
      
      if (nTrkV0[jj] < 0) {
        AliWarning(MSGWARNING("Error in V0s"));
      }
      
      decLenXY[jj] = primVtx->DistanceXYToVertex(svtx);
      
      if (signLxy < 0.) decLenXY[jj] *= -1.;
      
      errdecLenXY[jj] = primVtx->ErrorDistanceXYToVertex(svtx);
      sigdecLenXY[jj] = decLenXY[jj]/errdecLenXY[jj];
      
      Double_t dca0aod[2] = {-999., -999.};
      Double_t dca1aod[2] = {-999., -999.};
      Double_t dca2aod[2] = {-999., -999.};
      Double_t cov0aod[3] = {-999., -999., -999.};
      Double_t cov1aod[3] = {-999., -999., -999.};
      Double_t cov2aod[3] = {-999., -999., -999.};
      
      AliAODTrack *track0 = (AliAODTrack *)svtx->GetDaughter(0);
      AliAODTrack *track1 = (AliAODTrack *)svtx->GetDaughter(1);
      AliAODTrack *track2 = (AliAODTrack *)svtx->GetDaughter(2);
      
      track0->PropagateToDCA(primVtx, 0., 10000., dca0aod, cov0aod);
      track1->PropagateToDCA(primVtx, 0., 10000., dca1aod, cov1aod);
      track2->PropagateToDCA(primVtx, 0., 10000., dca2aod, cov2aod);
      
      ipR1[jj] = dca0aod[0];
      ipR2[jj] = dca1aod[0];
      ipR3[jj] = dca2aod[0];
      
      if (cov0aod[0]) sigipR1[jj] = dca0aod[0]/TMath::Sqrt(cov0aod[0]);;
      if (cov1aod[0]) sigipR2[jj] = dca1aod[0]/TMath::Sqrt(cov1aod[0]);
      if (cov2aod[0]) sigipR3[jj] = dca2aod[0]/TMath::Sqrt(cov2aod[0]);
      
      fTagger->GetVtxPxy(svtx, vtxP);
    }
    
    TMath::Sort(nSVtx, decLenXY, idxLxy);
    if (nSVtx > 0) {
      
      point[7]  = invMasses[idxLxy[0]];
      point[8]  = sigmavertex[idxLxy[0]];
      point[9]  = decLenXY[idxLxy[0]];
      point[10] = sigdecLenXY[idxLxy[0]];
      point[11] = sigipR1[idxLxy[0]];
      point[12] = sigipR2[idxLxy[0]];
      point[13] = sigipR3[idxLxy[0]];
      point[14] = nTrkV0[idxLxy[0]];
    }
    // for now we only take the most dispaced vertex
    
    // if(nvtx>1){
    //   point[7]=decLenXY[indexLxy[1]];
    //   point[10]=invMasses[indexLxy[1]];
    //   point[13]=nRealVtx[indexLxy[1]];
    //   point[16]=nFromBVtx[indexLxy[1]];
    //   point[19]=nFromPromptDVtx[indexLxy[1]];
    // 	 point[25]=sigmavertex[indexLxy[1]];
    // 	 point[28]=sumOfSqs[indexLxy[1]];
    // 	 // point[31]=errdecLenXY[indexLxy[0]];
    //   }
    
    // if(nvtx>2){
    
    //   point[8]=decLenXY[indexLxy[2]];
    //   point[11]=invMasses[indexLxy[2]];
    //   point[14]=nRealVtx[indexLxy[2]];
    //   point[17]=nFromBVtx[indexLxy[2]];
    //   point[20]=nFromPromptDVtx[indexLxy[2]];
    // 	 point[26]=sigmavertex[indexLxy[2]];
    // 	 point[29]=sumOfSqs[indexLxy[2]];
    // 	 // point[32]=errdecLenXY[indexLxy[0]];
    //   }
    
    // Calculate sum of inv masses of the 3 more displaced vertices
    // for(Int_t ivtx=0;ivtx<3;ivtx++){
    //   if(nvtx>ivtx) point[5]+=invMasses[indexLxy[ivtx]];
    //   }
    
    
    // Calculate sum of inv masses of the 3 more displaced vertices
    //    for (Int_t ivtx = 0; ivtx < 3; ++ivtx) {
    //      if (nSVtx > ivtx) point[22] += invMasses[idxLxy[ivtx]];
    //    }
    
    delete[] idxLxy;
    delete[] nTrkV0;
    
    delete[] decLenXY;
    delete[] invMasses;
    delete[] errdecLenXY;
    delete[] sigdecLenXY;
    delete[] sigmavertex;
    
    delete[] ipR1;
    delete[] ipR2;
    delete[] ipR3;
    
    delete[] sigipR1;
    delete[] sigipR2;
    delete[] sigipR3;
    
    delete[] sumOfSqs;
    delete[] sigsumOfSqs;
    
  } // end if (vertices && primVtx)
  
  TArrayD *apoint = new TArrayD(kNvar, point);
  FillStep(step, apoint, wght);
  
  return;
}

//_____________________________________________________________________________________
void AliHFJetsContainerVertex::FillStepQaVtx(CFSteps                    step,
                                             const Int_t                nSVtx,
                                             Double_t                   evtCent,
                                             const AliAODVertex        *primVtx,
                                             const AliEmcalJet         *jet,
                                             const TClonesArray        *arrVtxHF,
                                             const TClonesArray        *mcPart,
                                             vctr_pair_dbl_int         &arrVtxDisp,
                                             Double_t                  *p,
                                             Double_t                   wght)
{
  
  if (fType != kQaVtx) {
    AliErrorF(MSGERROR("This method is available only for container type kQaVtx: \
                       you are trying to fill %s!"), fgStrContType(fType));
  }
  
  Double_t xyz[3],vtxVect[3],jetP[3];
  Double_t xyzPrim[3];
  Double_t cosTheta;
  
  primVtx->GetXYZ(xyzPrim);
  jet->PxPyPz(jetP);
  
  Int_t *indexLxy = new Int_t[nSVtx];
  const Int_t kNvar = 20;
  Double_t pointVtxProp[kNvar] = {
    evtCent*1.,                         /* 1 */
    jet->Pt(),                          /* 2 */
    jet->Eta(),                         /* 3 */
    jet->Phi()-TMath::Pi(),             /* 4 */
    jet->GetNumberOfTracks()*1.,        /* 5 */
    0.,0.,0.,0.,0,0,0,0,0,0,0,0,0,     /* 6-18 */
    p[0],                               /* 19 */
    p[1],                               /* 20 */
  };
  
  Double_t *decLengthXY     = new Double_t[nSVtx];
  Double_t *invMasses       = new Double_t[nSVtx];
  Double_t *nRealVtx        = new Double_t[nSVtx];
  Double_t *nFromBVtx       = new Double_t[nSVtx];
  Double_t *nFromPromptDVtx = new Double_t[nSVtx];
  Double_t xMC,yMC;
  Double_t vtxP[3],vtxPt,signLxy;
  Int_t nvtxMC=0;
  
  for (Int_t jj = 0; jj < nSVtx; jj++) {
    xMC=-99999.;
    yMC=-99999.;
    AliAODVertex *vtx=(AliAODVertex *)arrVtxHF->UncheckedAt(jj);
    
    Double_t chi2ndf=vtx->GetChi2perNDF();
    invMasses[jj]=fTagger->GetVertexInvariantMass(vtx);
    Double_t sigvert = arrVtxDisp[jj].first;
    
    nRealVtx[jj]=-1;
    nFromBVtx[jj]=-1;
    
    vtx->GetXYZ(xyz);
    vtxVect[0]=xyz[0]-xyzPrim[0];
    vtxVect[1]=xyz[1]-xyzPrim[1];
    vtxVect[2]=xyz[2]-xyzPrim[2];
    signLxy=vtxVect[0]*jetP[0]+vtxVect[1]*jetP[1];
    
    Double_t absJetPt=TMath::Sqrt(jetP[0]*jetP[0]+jetP[1]*jetP[1]);
    Double_t absVtxVect=TMath::Sqrt(vtxVect[0]*vtxVect[0]+vtxVect[1]*vtxVect[1]);
    cosTheta=signLxy/(absJetPt*absVtxVect);//angle between jet and Lxy
    
    decLengthXY[jj]=TMath::Sqrt((xyz[0]-xyzPrim[0])*(xyz[0]-xyzPrim[0])+(xyz[1]-xyzPrim[1])*(xyz[1]-xyzPrim[1]));
    if(signLxy<0.){
      decLengthXY[jj]*=-1.;
    }
    
    fTagger->GetVtxPxy(vtx,vtxP);
    vtxPt=TMath::Sqrt(vtxP[0]*vtxP[0]+vtxP[1]*vtxP[1]);
    pointVtxProp[5]=vtxPt;
    pointVtxProp[6]=invMasses[jj];
    pointVtxProp[7]=decLengthXY[jj];
    pointVtxProp[8]=chi2ndf;
    
    if(mcPart){
      Int_t nfromBandD=0,nfromD=0,nfromPromptD=0;
      fTagger->GetNTracksFromCommonVertex(vtx,mcPart,nvtxMC,xMC,yMC,nfromBandD,nfromD,nfromPromptD);
      pointVtxProp[9]=nvtxMC;
      pointVtxProp[10]=xyz[0]-xMC;
      pointVtxProp[11]=xyz[1]-yMC;
      pointVtxProp[16]=nfromBandD;
      pointVtxProp[17]=nfromD;
      
      nRealVtx[jj]=nvtxMC;
      nFromBVtx[jj]=nfromBandD;
      nFromPromptDVtx[jj]=nfromPromptD;
    }
    
    pointVtxProp[12]=sigvert;
    pointVtxProp[13]=cosTheta*TMath::Abs(decLengthXY[jj]);
    pointVtxProp[14]=TMath::Sqrt(decLengthXY[jj]*decLengthXY[jj]-pointVtxProp[14]*pointVtxProp[14]);
    pointVtxProp[15]=cosTheta;
    TArrayD *apointVtxProp = new TArrayD(kNvar, pointVtxProp);
    FillStep(step, apointVtxProp, wght);
  }
  
  delete[] indexLxy;
  delete[] decLengthXY;
  delete[] invMasses;
  delete[] nRealVtx;
  delete[] nFromBVtx;
  delete[] nFromPromptDVtx;
  
  return;
}

//_____________________________________________________________________________________
void AliHFJetsContainerVertex::Copy(TObject &c) const
{
  // copy function
  AliHFJetsContainerVertex &target = (AliHFJetsContainerVertex &) c;
  if (fType)
    target.fType = fType;
  
}

//_____________________________________________________________________________________
void AliHFJetsContainerVertex::CreateContainerVertex(EContTYPE contType) {

  TString stype, vars;
  
  fType = contType;
  
  switch (fType) {
    
    case kJetVtxSim :
      
      stype.Form("Jet vertices properties simulation.");
      // Relevant variables
      vars.Form("JetArea;JetPt_wBkgRej;nRecoVtx;mass1;sigVtx1;Lxy1;SLxy1;SIP1;SIP2;SIP3;nV0Trks;partBP;partBH;ptBP;ptBH;nReal1;nFromB1;nFromPrD1;SumMass3MostDispl");
      
      break;
    
    case kJetVtxData :
      
      stype.Form("Jet vertices properties data.");
      // Relevant variables
      vars.Form("JetArea;JetPt_wBkgRej;nRecoVtx;mass1;sigVtx1;Lxy1;SLxy1;SIP1;SIP2;SIP3;nV0Trks");

      break;
  
    case kQaVtx :
      
      stype.Form("QA secondary vertices properties.");
      // Relevant variables
      vars.Form("nTrk;ptVtx;mass;Lxy;chi2/ndf;nRealVtx;deltaX;deltaY;sigVtx;LxyJet;LxyJetPerp;cosTheta;nFromB;nFromBD;partBP;partBH");
      
      break;
    
    default:
      AliError(MSGERROR("Not a valid container type!"));
      return;
  }

  AliInfoF(MSGINFO("Container type set to %s: %s"), fgStrContType(fType), stype.Data());
  
  // Get binning for each variable
  TObjArray  *arr = vars.Tokenize(";");
  TIter next(arr);
  
  const Int_t nvars = arr->GetEntriesFast();
  Int_t nbins[nvars];       // number of bins for each variable
  
  const char *axistitle[nvars]; // axis title for each variable
  const char *axisname[nvars]; // axis name for each variable
  
  Double_t *binning[nvars]; // array of bins for each variable

  TObjString *objstr;

  Int_t i = 0;
  while ((objstr=(TObjString *)next())) {
    
    binning[i] = new Double_t[1000];
    
    GetBinningVertex(objstr->GetString(), nbins[i], binning[i], axistitle[i]);
    
    axisname[i] = objstr->GetName();
    i++;
  }
  
  AliHFJetsContainer::CreateCustomContainer(nvars, axisname, nbins, binning, axistitle);

  // Delete arrays
  for (Int_t k = 0; k < nvars; ++k)
    delete [] binning[k];
  
  delete arr;
}

//_____________________________________________________________________________________
void AliHFJetsContainerVertex::GetBinningVertex(TString     var,
                                               Int_t      &nBins,
                                               Double_t   *bins,
                                               const char *&axistitle)
{
  // Assigns variable-specific bin edges
  // (you can define array of bins "by hand" to allow for non uniform bin width)
  
  Float_t binmin=0., binmax=0.;    
  
  if (var.EqualTo("JetArea")) {
    
    axistitle="Charged area";
    nBins = 50; binmin= 0.; binmax= 1.;
  }
  else if (var.EqualTo("JetPt_wBkgRej")) {
    
    axistitle = "p_{T,jet} (GeV/c)";
    nBins = 100; binmin= 0.; binmax= 100.;
  }
  else if (var.EqualTo("nRecoVtx")) {

    axistitle="N_{vtx,reco}";
    nBins = 21; binmin= -0.5; binmax= 20.5;
  }
  else if ( var.EqualTo("nV0Trks") ) {
    
    axistitle = "nV0 tracks";
    nBins = 8; binmin= -1.5; binmax= 6.5;
  }
  else if (var.EqualTo("nTrk")) {
    
      axistitle="N_{trk}";
      nBins = 20; binmin= 0.99; binmax= 20.99;
  }
  else if (var.EqualTo("nEle")) {
    
      axistitle="N_{ele}";
      nBins = 5; binmin= 0.; binmax= 4.99;
  }
  else if (var.EqualTo("partDP")) {
    
      axistitle="ID_{parton} DP";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("partBP")) {
    
      axistitle="ID_{parton} BP";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("partBH")) {
    
      axistitle="ID_{parton} BH";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("ptDP")) {
    
      axistitle="p_{T,part} DP (GeV/c)";
      nBins = 60; binmin= 0.; binmax= 60.;
  }
  else if (var.EqualTo("ptBP")) {
    
      axistitle="p_{T,part} BP (GeV/c)";
      nBins = 60; binmin= 0.; binmax= 60.;
  }
  else if (var.EqualTo("ptBH")) {
    
      axistitle="p_{T,part} BH (GeV/c)";
      nBins = 60; binmin= 0.; binmax= 60.;
  }
  else if (var.EqualTo("ptVtx")) {
    
      axistitle="p_{T,vtx} (GeV/c)";
      nBins = 20; binmin= 0.; binmax= 20.;
  }
  else if (var.EqualTo("mass")) {
    
      axistitle="mass (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 5.;
  }
  else if (var.EqualTo("mass1")) {
    
      axistitle="mass_{1} (GeV/c^{2})";
      nBins = 40; binmin= 0.; binmax= 5.;
  }
  else if (var.EqualTo("mass2")) {
    
      axistitle="mass_{2} (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 5.;
  }
  else if (var.EqualTo("mass3")) {
    
      axistitle="mass_{3} (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 5.;
  }
  else if (var.EqualTo("SumMass3MostDispl")) {
    
      axistitle="#Sigma mass (GeV/c^{2})";
      nBins = 20; binmin= 0.; binmax= 20.;
  }
  else if (var.EqualTo("Lxy")) {
    
      axistitle="L_{xy} (cm)";
      nBins = 200; binmin= -1.; binmax= 1.;
  }
  else if (var.EqualTo("Lxy1")) {
    
      axistitle="L_{xy,1} (cm)";
      nBins = 50; binmin= 0.; binmax= 1.;
  }
  else if (var.EqualTo("Lxy2")) {
    
      axistitle="L_{xy,2} (cm)";
      nBins = 50; binmin= 0.; binmax= 0.5;
  }
  else if (var.EqualTo("Lxy3")) {
    
      axistitle="L_{xy,3} (cm)";
      nBins = 50; binmin= 0.; binmax= 0.5;
  }
  else if (var.EqualTo("chi2/ndf")){
      axistitle="#chi^{2}/NDF";
      nBins = 10; binmin= 0.; binmax= 10.;
  }
  else if (var.EqualTo("nRealVtx")) {
    
      axistitle="N_{vtx,real}";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("deltaX")) {
    
      axistitle="#Delta x (cm)";
      nBins = 15; binmin= -0.03; binmax= 0.03;
  }
  else if (var.EqualTo("deltaY")) {
    
      axistitle="#Delta y (cm)";
      nBins = 15; binmin= -0.03; binmax= 0.03;
  }
  else if (var.EqualTo("sigVtx")) {
    
      axistitle="#sigma_{vtx}";
      nBins = 20; binmin= 0.; binmax= 0.1;
  }
  else if (var.EqualTo("LxyJet")) {
    
      axistitle="#L_{xy,jet} (cm)";
      nBins = 100; binmin= -1.; binmax= 1.;
  }
  else if (var.EqualTo("LxyJetPerp")) {
    
      axistitle="#L_{xy,jet} #perp (cm)";
      nBins = 30; binmin= 0.; binmax= 0.2;
  }
  else if (var.EqualTo("cosTheta")) {
    
      axistitle="cos#theta";
      nBins = 50; binmin= -1.; binmax= 1.;
  }
  else if (var.EqualTo("nFromB")) {
    
      axistitle="N (from B)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nFromB1")) {
    
      axistitle="N (from B1)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nFromB2")) {
    
      axistitle="N (from B2)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nFromB3")) {
    
      axistitle="N (from B3)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nFromBD")) {
    
      axistitle="N (from D from B)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nFromPrD1")) {
    
      axistitle="N (from prompt D1)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nFromPrD2")) {
    
      axistitle="N (from prompt D2)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nFromPrD3")) {
    
      axistitle="N (from prompt D3)";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nReal1")) {
    
      axistitle="N_{real,1}";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nReal2")) {
    
      axistitle="N_{real,2}";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("nReal3")) {
    
      axistitle="N_{real,3}";
      nBins = 5; binmin= -0.5; binmax= 4.5;
  }
  else if (var.EqualTo("sigVtx1")) {
    
    axistitle="#sigma_{vertex,1}";
    nBins = 20; binmin= 0.; binmax= 0.1;
  }
  else if (var.EqualTo("sigvtx2")) {
    
    axistitle="#sigma_{vertex,2}";
    nBins = 20; binmin= 0.; binmax= 0.1;
  }
  else if (var.EqualTo("sigvtx3")) {
    
    axistitle="#sigma_{vertex,3}";
    nBins = 20; binmin= 0.; binmax= 0.1;
  }
  else if (var.EqualTo("SLxy1")) {
    
    axistitle="L_{xy}/#sigma_{L_{xy}}";
    nBins = 100; binmin= 0.; binmax= 100.;
  }
  else if (var.EqualTo("SIP1")) {
    
    axistitle="IP_{xy}/#sigma_{IP_{xy}}";
    nBins = 40; binmin= -20.; binmax= 20.;
  }
  else if (var.EqualTo("SIP2")) {
    
    axistitle="IP_{xy}/#sigma_{IP_{xy}}";
    nBins = 40; binmin= -20.; binmax= 20.;
  }
  else if (var.EqualTo("SIP3")) {
    
    axistitle = "IP_{xy}/#sigma_{IP_{xy}}";
    nBins = 40; binmin= -20.; binmax= 20.;
  }
  else {
    
    AliErrorF("Variable %s not defined!", var.Data());
  }
  
  // Define regular binning
  Double_t binwidth = (binmax - binmin) / (1. * nBins);
  for (Int_t j = 0; j <= nBins; ++j)
    bins[j] = binmin + j * binwidth;
}
