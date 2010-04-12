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
//
// Secondary vertexing construction Class
//  Construct secondary vertex from Beauty hadron with electron and
//  hadrons, then apply selection criteria
//
// Authors:
//   MinJung Kweon <minjung@physi.uni-heidelberg.de>
//

#include <TH2F.h>
#include <TIterator.h>
#include <TParticle.h>

#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliVTrack.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include "AliHFEsecVtx.h"
#include <AliKFParticle.h>
#include <AliKFVertex.h>
#include <AliLog.h>
#include <AliStack.h>
#include <AliAODMCParticle.h>
#include "AliHFEpairs.h"
#include "AliHFEsecVtxs.h"
#include "AliHFEtrackFilter.h"

ClassImp(AliHFEsecVtx)

//_______________________________________________________________________________________________
AliHFEsecVtx::AliHFEsecVtx():
  fFilter(0x0)
  ,fESD1(0x0)
  ,fAOD1(0x0)
  ,fStack(0x0)
  ,fUseMCPID(kFALSE)
  ,fkSourceLabel()
  ,fNparents(0)
  ,fParentSelect()
  ,fPtRng()
  ,fDcaCut()
  ,fNoOfHFEpairs(0)
  ,fNoOfHFEsecvtxs(0)
  ,fHFEpairs(0x0)
  ,fHFEsecvtxs(0x0)
  ,fMCArray(0x0)
  ,fPVx(0)
  ,fPVy(0)
  ,fCosPhi(-1)
  ,fSignedLxy(-1)
  ,fKFchi2(-1)
  ,fInvmass(-1)
  ,fInvmassSigma(-1)
  ,fKFip(0)
  ,fPairQA(0x0)
  ,fSecvtxQA(0x0)
  ,fSecVtxList(0x0)
{ 
  //
  // Default constructor
  //

  Init();
}

//_______________________________________________________________________________________________
AliHFEsecVtx::AliHFEsecVtx(const AliHFEsecVtx &p):
  TObject(p)
  ,fFilter(0x0)
  ,fESD1(0x0)
  ,fAOD1(0x0)
  ,fStack(0x0)
  ,fUseMCPID(p.fUseMCPID)
  ,fkSourceLabel()
  ,fNparents(p.fNparents)
  ,fParentSelect()
  ,fPtRng()
  ,fDcaCut()
  ,fNoOfHFEpairs(p.fNoOfHFEpairs)
  ,fNoOfHFEsecvtxs(p.fNoOfHFEsecvtxs)
  ,fHFEpairs(0x0)
  ,fHFEsecvtxs(0x0)
  ,fMCArray(0x0)
  ,fPVx(p.fPVx)
  ,fPVy(p.fPVy)
  ,fCosPhi(p.fCosPhi)
  ,fSignedLxy(p.fSignedLxy)
  ,fKFchi2(p.fKFchi2)
  ,fInvmass(p.fInvmass)
  ,fInvmassSigma(p.fInvmassSigma)
  ,fKFip(p.fKFip)
  ,fPairQA(0x0)
  ,fSecvtxQA(0x0)
  ,fSecVtxList(0x0)
{ 
  //
  // Copy constructor
  //
  fFilter = new AliHFEtrackFilter(*p.fFilter);
}

//_______________________________________________________________________________________________
AliHFEsecVtx&
AliHFEsecVtx::operator=(const AliHFEsecVtx &)
{
  //
  // Assignment operator
  //

  AliInfo("Not yet implemented.");
  return *this;
}

//_______________________________________________________________________________________________
AliHFEsecVtx::~AliHFEsecVtx()
{
  //
  // Destructor
  //

  //cout << "Analysis Done." << endl;
  delete fFilter;
}

//__________________________________________
void AliHFEsecVtx::Init()
{
  //
  // set pdg code and index
  //

  fNparents = 7;

  fParentSelect[0][0] =  411;
  fParentSelect[0][1] =  421;
  fParentSelect[0][2] =  431;
  fParentSelect[0][3] = 4122;
  fParentSelect[0][4] = 4132;
  fParentSelect[0][5] = 4232;
  fParentSelect[0][6] = 4332;

  fParentSelect[1][0] =  511;
  fParentSelect[1][1] =  521;
  fParentSelect[1][2] =  531;
  fParentSelect[1][3] = 5122;
  fParentSelect[1][4] = 5132;
  fParentSelect[1][5] = 5232;
  fParentSelect[1][6] = 5332;

  // momentum ranges to apply pt dependent cuts
  fPtRng[0] = 1.0;
  fPtRng[1] = 2.0;
  fPtRng[2] = 2.5;
  fPtRng[3] = 3.0;
  fPtRng[4] = 5.0;
  fPtRng[5] = 12.0;
  fPtRng[6] = 20.0;

  // momentum dependent DCA cut values (preliminary)
  fDcaCut[0] = 0.04;
  fDcaCut[1] = 0.03;
  fDcaCut[2] = 0.02;
  fDcaCut[3] = 0.015;
  fDcaCut[4] = 0.01;
  fDcaCut[5] = 0.005;

  fFilter = new AliHFEtrackFilter("SecVtxFilter");
  fFilter->MakeCutStepRecKineITSTPC();
  fFilter->MakeCutStepPrimary();
} 

void AliHFEsecVtx::Process(AliVTrack *signalTrack){ 
  if(signalTrack->Pt() < 1.0) return;
  AliESDtrack *track = dynamic_cast<AliESDtrack *>(signalTrack);
  InitHFEpairs();
  InitHFEsecvtxs();
  AliESDtrack *htrack = 0x0; 
  fFilter->Flush();
  fFilter->FilterTracks(fESD1);
  TObjArray *candidates = fFilter->GetFilteredTracks();
  TIterator *trackIter = candidates->MakeIterator();
  while((htrack = dynamic_cast<AliESDtrack *>(trackIter->Next()))){
    if(track->GetID() == htrack->GetID()) continue; // since it is for tagging single electron, don't need additional condition
    if (htrack->Pt()<1.0) continue;
    PairAnalysis(track, htrack, htrack->GetID()); // e-h pairing
  }
  delete trackIter;
  /*for(int ip=0; ip<fSecVtx->HFEpairs()->GetEntriesFast(); ip++){
      if(HasMCData()){
        AliHFEpairs *pair = (AliHFEpairs*) (fSecVtx->HFEpairs()->UncheckedAt(ip));
        if(!(pair->GetPairCode()>1. && pair->GetPairCode()<4.))  // apply various cuts
        fSecVtx->HFEpairs()->RemoveAt(ip);
      }
    }*/
  HFEpairs()->Compress();
  RunSECVTX(track); // secondary vertexing with e,h1,h2,.. tracks
  for(int ip=0; ip<HFEsecvtxs()->GetEntriesFast(); ip++){
    AliHFEsecVtxs *secvtx=0x0;
    secvtx = (AliHFEsecVtxs*) (HFEsecvtxs()->UncheckedAt(ip));
    // here you apply cuts, then if it doesn't pass the cut, remove it from the fSecVtx->HFEsecvtxs() 
  }
  DeleteHFEpairs();
  DeleteHFEsecvtxs();
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::CreateHistograms(TList *qaList)
{ 
  //
  // create histograms
  //

  fSecVtxList = new TList;
  fSecVtxList->SetName("SecVtx");

  MakeContainer();
  /*
  fkSourceLabel[kAll]="all";
  fkSourceLabel[kDirectCharm]="directCharm";
  fkSourceLabel[kDirectBeauty]="directBeauty";
  fkSourceLabel[kBeautyCharm]="beauty2charm";
  fkSourceLabel[kGamma]="gamma";
  fkSourceLabel[kPi0]="pi0";
  fkSourceLabel[kElse]="others";
  fkSourceLabel[kBeautyGamma]="beauty22gamma";
  fkSourceLabel[kBeautyPi0]="beauty22pi0";
  fkSourceLabel[kBeautyElse]="beauty22others";

  TString hname;
  TString hnopt="secvtx_";
  for (Int_t isource = 0; isource < 10; isource++ ){
  }
  */

  qaList->AddLast(fSecVtxList);
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::GetPrimaryCondition()
{
  //
  // get primary characteristics and set
  //

  if (fESD1) {
    AliKFVertex primVtxCopy(*(fESD1->GetPrimaryVertex()));
    if( primVtxCopy.GetNDF() <1 ) return;
    fPVx = primVtxCopy.GetX();
    fPVx = primVtxCopy.GetY();
  }
  else if(fAOD1) {
    AliKFVertex primVtxCopy(*(fAOD1->GetPrimaryVertex()));
    if( primVtxCopy.GetNDF() <1 ) return;
    fPVx = primVtxCopy.GetX();
    fPVx = primVtxCopy.GetY();
  }
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::PairAnalysis(AliVTrack* track1, AliVTrack* track2, Int_t index2)
{
  //
  // calculate e-h pair characteristics and tag pair 
  //

  //later consider the below
  Float_t dca1[2]={-999.,-999.}, dca2[2]={-999.,-999.};
  Float_t cov1[3]={-999.,-999.,-999.}, cov2[3]={-999.,-999.,-999.};

  if (IsAODanalysis()){
    AliESDtrack esdTrk1(track1);
    AliESDtrack esdTrk2(track2);
    esdTrk1.PropagateToDCA(fAOD1->GetPrimaryVertex(),0,10000,(Double_t*)dca1,(Double_t*)cov1);
    esdTrk2.PropagateToDCA(fAOD1->GetPrimaryVertex(),0,10000,(Double_t*)dca2,(Double_t*)cov2);
  }
  else {
    ((AliESDtrack*)track1)->GetImpactParameters(dca1,cov1);
    ((AliESDtrack*)track2)->GetImpactParameters(dca2,cov2);
  }

  // apply pt dependent dca cut on hadrons
  for(int ibin=0; ibin<6; ibin++){
    if((track2->Pt()>fPtRng[ibin] && track2->Pt()<fPtRng[ibin+1]) && TMath::Abs(dca2[0])<fDcaCut[ibin]) return;
  }

  // get KF particle input pid
  Int_t pdg1 = GetPDG(track1);
  Int_t pdg2 = GetPDG(track2);
        
  if(pdg1==-1 || pdg2==-1) {
    //printf("out if considered pid range \n");
    return;
  }

  // create KF particle of pair
  if(IsAODanalysis()) AliKFParticle::SetField(fAOD1->GetMagneticField());
  else AliKFParticle::SetField(fESD1->GetMagneticField()); 
  AliKFParticle kfTrack1(*track1, pdg1);
  AliKFParticle kfTrack2(*track2, pdg2);

  AliKFParticle kfSecondary(kfTrack1,kfTrack2);

  //secondary vertex point from kf particle
  Double_t kfx = kfSecondary.GetX();
  Double_t kfy = kfSecondary.GetY();
  //Double_t kfz = kfSecondary.GetZ();
        
  //momentum at the decay point from kf particle
  Double_t kfpx = kfSecondary.GetPx();
  Double_t kfpy = kfSecondary.GetPy();
  //Double_t kfpz = kfSecondary.GetPz();

  Double_t dx = kfx-fPVx;
  Double_t dy = kfy-fPVy;

  // discriminating variables 
  // invariant mass of the KF particle
  Double_t invmass = -1;
  Double_t invmassSigma = -1;
  kfSecondary.GetMass(invmass,invmassSigma);

  // chi2 of the KF particle
  Double_t kfchi2 = -1;
  if(kfSecondary.GetNDF()>0) kfchi2 = TMath::Sqrt(TMath::Abs(kfSecondary.GetChi2()/kfSecondary.GetNDF()));

  // opening angle between two particles in XY plane
  Double_t phi = kfTrack1.GetAngleXY(kfTrack2);
  Double_t cosphi = TMath::Cos(phi);

  // projection of kf vertex vector to the kf momentum direction 
  Double_t signedLxy=-999.;
  if((dx*kfpx+dy*kfpy)>0) signedLxy = TMath::Sqrt(dx*dx+dy*dy);  
  if((dx*kfpx+dy*kfpy)<0) signedLxy = -1*TMath::Sqrt(dx*dx+dy*dy);  
  //[the other way to think about]
  //Double_t psqr = kfpx*kfpx+kfpy*kfpy;
  //if(psqr>0) signedLxy=(dx*kfpx+dy*kfpy)/TMath::Sqrt(psqr);  

  // DCA from primary to e-h KF particle (impact parameter of KF particle)
  Double_t vtx[2]={fPVx, fPVy}; 
  Double_t kfip = kfSecondary.GetDistanceFromVertexXY(vtx);

  Int_t paircode = -1;
  if (HasMCData()) paircode = GetPairCode(track1,track2); 

  AliHFEpairs hfepair; 
  hfepair.SetTrkLabel(index2);
  hfepair.SetInvmass(invmass);
  hfepair.SetKFChi2(kfchi2);
  hfepair.SetOpenangle(phi);
  hfepair.SetCosOpenangle(cosphi);
  hfepair.SetSignedLxy(signedLxy);
  hfepair.SetKFIP(kfip);
  hfepair.SetPairCode(paircode);
  AddHFEpairToArray(&hfepair);
  fNoOfHFEpairs++; 

  // fill into container for later QA
  Double_t dataE[6];
  dataE[0]=invmass;
  dataE[1]=kfchi2;
  dataE[2]=phi;
  dataE[3]=signedLxy;
  dataE[4]=kfip;
  dataE[5]=paircode;
  /*
  dataE[6]=TMath::Abs(dca1[0]);
  dataE[7]=TMath::Abs(dca2[0]);
  //if(cov1[0]>0) dataE[6]=Double_t(dca1[0]/cov1[0]);
  //if(cov2[0]>0) dataE[7]=Double_t(dca2[0]/cov2[0]);
  dataE[8]=track1->Pt();
  dataE[9]=track2->Pt();
  */
  fPairQA->Fill(dataE);

}

//_______________________________________________________________________________________________
void AliHFEsecVtx::RunSECVTX(AliVTrack *track)
{
  //
  // run secondary vertexing algorithm and do tagging
  //

  //printf("number of considered pairs= %d\n",HFEpairs()->GetEntriesFast());
  FindSECVTXCandid(track);         
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::FindSECVTXCandid(AliVTrack *track) 
{
  //
  // find secondary vertex candidate and store them into container
  //

  AliVTrack *htrack[20];
  Int_t htracklabel[20];
  Double_t vtxchi2cut=3.; // testing cut 
  Double_t dataE[6]={-999.,-999.,-999.,-999.,-1.,0};  
  if (HFEpairs()->GetEntriesFast()>20){
    AliDebug(3, "number of paired hadron is over maximum(20)");
    return; 
  }
        
  // get paired track objects  
  AliHFEpairs *pair=0x0;
  if (IsAODanalysis()){
    for (int ip=0; ip<HFEpairs()->GetEntriesFast(); ip++){
       pair = (AliHFEpairs*) (HFEpairs()->UncheckedAt(ip));
       htracklabel[ip] = pair->GetTrkLabel();
       htrack[ip] = fAOD1->GetTrack(pair->GetTrkLabel());
    }
  }
  else{
    for (int ip=0; ip<HFEpairs()->GetEntriesFast(); ip++){
       pair = (AliHFEpairs*) (HFEpairs()->UncheckedAt(ip));
       htracklabel[ip] = pair->GetTrkLabel();
       htrack[ip] = fESD1->GetTrack(pair->GetTrkLabel());
    }
  }
  // in case there is only one paired track with the electron, put pair characteristics into secvtx container
  // for the moment, I only apply pair vertex chi2 cut
  if (HFEpairs()->GetEntriesFast() == 1){
    if (pair->GetKFChi2()<vtxchi2cut) { // you can also put single track cut
      AliHFEsecVtxs hfesecvtx;
      hfesecvtx.SetTrkLabel1(pair->GetTrkLabel());
      hfesecvtx.SetTrkLabel2(-999);
      hfesecvtx.SetInvmass(pair->GetInvmass());
      hfesecvtx.SetKFChi2(pair->GetKFChi2());
      hfesecvtx.SetSignedLxy(pair->GetSignedLxy());
      hfesecvtx.SetKFIP(pair->GetKFIP());
      AddHFEsecvtxToArray(&hfesecvtx);
      fNoOfHFEsecvtxs++; 

      dataE[0]=pair->GetInvmass();
      dataE[1]=pair->GetKFChi2();
      dataE[2]=pair->GetSignedLxy();
      dataE[3]=pair->GetKFIP();
      if(HasMCData()) dataE[4]=GetElectronSource(TMath::Abs(track->GetLabel()));
      dataE[5]=2;
      fSecvtxQA->Fill(dataE);
    }
    return;
  }

  // in case there are multiple paired track with the electron, calculate secvtx characteristics
  // put the secvtx characteristics into container if it passes cuts
  for (int i=0; i<HFEpairs()->GetEntriesFast()-1; i++){
     for (int j=i+1; j<HFEpairs()->GetEntriesFast(); j++){
        CalcSECVTXProperty(track, htrack[i], htrack[j]);
        if (fKFchi2<vtxchi2cut) {
          AliHFEsecVtxs hfesecvtx;
          hfesecvtx.SetTrkLabel1(htracklabel[i]);
          hfesecvtx.SetTrkLabel2(htracklabel[j]);
          hfesecvtx.SetKFChi2(fKFchi2);
          hfesecvtx.SetInvmass(fInvmass);
          hfesecvtx.SetSignedLxy(fSignedLxy);
          hfesecvtx.SetKFIP(fKFip);
          AddHFEsecvtxToArray(&hfesecvtx);
          fNoOfHFEsecvtxs++; 

          dataE[0]=fInvmass;
          dataE[1]=fKFchi2;
          dataE[2]=fSignedLxy;
          dataE[3]=fKFip;
          if(HasMCData()) dataE[4]=GetElectronSource(TMath::Abs(track->GetLabel()));
          dataE[5]=3;
          fSecvtxQA->Fill(dataE);
        }
     }
  }
}

//_______________________________________________________________________________________________
void AliHFEsecVtx::CalcSECVTXProperty(AliVTrack* track1, AliVTrack* track2, AliVTrack* track3)
{
  //
  // calculate secondary vertex properties
  //

  // get KF particle input pid
  Int_t pdg1 = GetPDG(track1);
  Int_t pdg2 = GetPDG(track2);
  Int_t pdg3 = GetPDG(track3);

  if(pdg1==-1 || pdg2==-1 || pdg3==-1) {
    //printf("out if considered pid range \n");
    return;
  }

  // create KF particle of pair
  if(IsAODanalysis()) AliKFParticle::SetField(fAOD1->GetMagneticField());
  else AliKFParticle::SetField(fESD1->GetMagneticField());
  AliKFParticle kfTrack1(*track1, pdg1);
  AliKFParticle kfTrack2(*track2, pdg2);
  AliKFParticle kfTrack3(*track3, pdg3);

  AliKFParticle kfSecondary(kfTrack1,kfTrack2,kfTrack3);
        
  //secondary vertex point from kf particle
  Double_t kfx = kfSecondary.GetX();
  Double_t kfy = kfSecondary.GetY();
  //Double_t kfz = kfSecondary.GetZ();

  //momentum at the decay point from kf particle
  Double_t kfpx = kfSecondary.GetPx();
  Double_t kfpy = kfSecondary.GetPy();
  //Double_t kfpz = kfSecondary.GetPz();

  Double_t dx = kfx-fPVx;
  Double_t dy = kfy-fPVy;

  // discriminating variables ----------------------------------------------------------

  if(kfSecondary.GetNDF()>0) fKFchi2 = TMath::Sqrt(TMath::Abs(kfSecondary.GetChi2()/kfSecondary.GetNDF())); 

  // invariant mass of the KF particle
  kfSecondary.GetMass(fInvmass,fInvmassSigma);

  // DCA from primary to e-h KF particle (impact parameter of KF particle)
  Double_t vtx[2]={fPVx, fPVy};
  fKFip = kfSecondary.GetDistanceFromVertexXY(vtx);

  if((dx*kfpx+dy*kfpy)>0) fSignedLxy= TMath::Sqrt(dx*dx+dy*dy);
  if((dx*kfpx+dy*kfpy)<0) fSignedLxy= -1*TMath::Sqrt(dx*dx+dy*dy);
  //[the other way to think about] - projection of kf vertex vector to the kf momentum direction
  //Double_t psqr = kfpx*kfpx+kfpy*kfpy;
  //if(psqr>0) fSignedLxy=(dx*kfpx+dy*kfpy)/TMath::Sqrt(psqr);  
}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetMCPID(AliESDtrack *track) 
{
  //      
  // return mc pid
  //      

  Int_t label = TMath::Abs(track->GetLabel());
  TParticle* mcpart = fStack->Particle(label);
  if ( !mcpart ) return 0;
  Int_t pdgCode = mcpart->GetPdgCode();

  return pdgCode;
}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetPairOriginESD(AliESDtrack* trk1, AliESDtrack* trk2)
{
  //
  // return pdg code of the origin(source) of the pair 
  // 
  //
  // ---*---*---*-----ancester A----- track1
  //                        |____*______ 
  //                             |______ track2
  // => if they originated from same ancester, 
  //    then return "the absolute value of pdg code of ancester A"
  //
  // ---*---*---B hadron-----ancester A----- track1
  //                               |____*______ 
  //                                    |______ track2
  // => if they originated from same ancester, and this ancester originally comes from B hadrons
  //    then return -1*"the absolute value of pdg code of ancester A"
  //
  // caution : it can also return parton pdg code if it originated from same string or gluon spliting. 
  //           

  if (trk1->GetLabel()<0 || trk2->GetLabel()<0) return 0;
  TParticle* part1 = fStack->Particle(trk1->GetLabel());
  TParticle* part2 = fStack->Particle(trk2->GetLabel());
  TParticle* part2cp = part2;
  if (!(part1) || !(part2)) return 0;

  Int_t srcpdg = 0;

  //if the two tracks' mother's label is same, get the mother info
  //in case of charm, check if it originated from beauty
  for (Int_t i=0; i<10; i++){ //check up to 10 ancesters
     Int_t label1 = part1->GetFirstMother(); 
     if (label1 < 0) return 0;

     for (Int_t j=0; j<10; j++){ //check up to 10 ancesters
        Int_t label2 = part2->GetFirstMother(); 
        if (label2 < 0) break; 

        if (label1 == label2){ //check if two tracks are originated from same mother
          TParticle* commonmom = fStack->Particle(label2); 
          srcpdg = abs(commonmom->GetPdgCode()); 

          //check ancester to see if it is originally from beauty 
          for (Int_t k=0; k<10; k++){ //check up to 10 ancesters
             Int_t ancesterlabel = commonmom->GetFirstMother();
             if (ancesterlabel < 0) return srcpdg; // if there is no more commonancester, return commonmom's pdg  

             TParticle* commonancester = fStack->Particle(ancesterlabel);
             Int_t ancesterpdg = abs(commonancester->GetPdgCode());

             for (Int_t l=0; l<fNparents; l++){
                if (abs(ancesterpdg)==fParentSelect[1][l]){
                  srcpdg = -1*srcpdg; //multiply -1 for hadron from bottom
                  return srcpdg;
                }
             }
             commonmom = commonancester;
          }
        }
        part2 = fStack->Particle(label2); //if their mother is different, go to earlier generation of 2nd particle
        if (!(part2)) break;
     }
     part1 = fStack->Particle(label1); //if their mother is different, go to earlier generation of 1st particle
     part2 = part2cp;
     if (!(part1)) return 0;
  }

  return srcpdg; 
}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetPairOriginAOD(AliAODTrack* trk1, AliAODTrack* trk2)
{

  //
  // return pdg code of the origin(source) of the pair 
  // 
  //
  // ---*---*---*-----ancester A----- track1
  //                        |____*______ 
  //                             |______ track2
  // => if they originated from same ancester, 
  //    then return "the absolute value of pdg code of ancester A"
  //
  // ---*---*---B hadron-----ancester A----- track1
  //                               |____*______ 
  //                                    |______ track2
  // => if they originated from same ancester, and this ancester originally comes from B hadrons
  //    then return -1*"the absolute value of pdg code of ancester A"
  //
  // caution : it can also return parton pdg code if it originated from same string or gluon spliting. 
  //           

  if (trk1->GetLabel()<0 || trk2->GetLabel()<0) return 0;
  AliAODMCParticle *part1 = (AliAODMCParticle*)fMCArray->At(trk1->GetLabel());
  AliAODMCParticle *part2 = (AliAODMCParticle*)fMCArray->At(trk2->GetLabel());
  AliAODMCParticle *part2cp = part2;
  if (!(part1) || !(part2)) return 0;

  Int_t srcpdg = 0;

  //if the two tracks' mother's label is same, get the mother info
  //in case of charm, check if it originated from beauty
  for (Int_t i=0; i<10; i++){ //check up to 10 ancesters
     Int_t label1 = part1->GetMother(); 
     if (label1 < 0) return 0;

     for (Int_t j=0; j<10; j++){ //check up to 10 ancesters
        Int_t label2 = part2->GetMother(); 
        if (label2 < 0) return 0; 

        if (label1 == label2){ //check if two tracks are originated from same mother
          AliAODMCParticle *commonmom = (AliAODMCParticle*)fMCArray->At(label1);
          srcpdg = abs(commonmom->GetPdgCode()); 

          //check ancester to see if it is originally from beauty 
          for (Int_t k=0; k<10; k++){ //check up to 10 ancesters
             Int_t ancesterlabel = commonmom->GetMother();
             if (ancesterlabel < 0) return srcpdg; // if there is no more commonancester, return commonmom's pdg  

             AliAODMCParticle *commonancester = (AliAODMCParticle*)fMCArray->At(ancesterlabel);
             Int_t ancesterpdg = abs(commonancester->GetPdgCode());

             for (Int_t l=0; l<fNparents; l++){
                if (abs(ancesterpdg)==fParentSelect[1][l]){
                  srcpdg = -1*srcpdg; //multiply -1 for charm from bottom
                  return srcpdg;
                }
             }
             commonmom = commonancester;
          }
        }
        part2 = (AliAODMCParticle*)fMCArray->At(label2); //if their mother is different, go to earlier generation of 2nd particle
        if (!(part2)) break;
     }
     part1 = (AliAODMCParticle*)fMCArray->At(label1); //if their mother is different, go to earlier generation of 2nd particle
     part2 = part2cp;
     if (!(part1)) return 0;
  }

  return srcpdg; 
}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetPairCode(const AliVTrack* const trk1, const AliVTrack* const trk2)
{
  //           
  // return pair code which is predefinded as:
  //  kDirectCharm, kDirectBeauty, kBeautyCharm, kGamma, kPi0, 
  //  kElse, kBeautyGamma, kBeautyPi0, kBeautyElse
  //           

  Int_t srcpdg = -1;
  Int_t srccode = kElse;

  if(IsAODanalysis()) srcpdg = GetPairOriginAOD((AliAODTrack*)trk1,(AliAODTrack*)trk2);
  else srcpdg = GetPairOriginESD((AliESDtrack*)trk1,(AliESDtrack*)trk2);

  if (srcpdg < 0) srccode = kBeautyElse;
  for (Int_t i=0; i<fNparents; i++){
    if (abs(srcpdg)==fParentSelect[0][i]) {
      if (srcpdg>0) srccode = kDirectCharm;
      if (srcpdg<0) srccode = kBeautyCharm;
    }
    if (abs(srcpdg)==fParentSelect[1][i]) {
      if (srcpdg>0) srccode = kDirectBeauty;
      if (srcpdg<0)  return kElse;
    }
  }
  if (srcpdg == 22) srccode = kGamma;
  if (srcpdg == -22) srccode = kBeautyGamma;
  if (srcpdg == 111) srccode = kPi0;
  if (srcpdg == -111) srccode = kBeautyPi0;

  return srccode;
}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetElectronSource(Int_t iTrack) 
{
  //
  // return decay electron's origin 
  //

  if (iTrack < 0) {
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  }

  TParticle* mcpart = fStack->Particle(iTrack);

//  if ( abs(mcpart->GetPdgCode()) != 11 ) return -1; // check if it is electron !

  Int_t iLabel = mcpart->GetFirstMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  }

  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  TParticle *partMother = fStack->Particle(iLabel);
  Int_t maPdgcode = partMother->GetPdgCode();

  // if the mother is charmed hadron  
  if ( int(abs(maPdgcode)/100.) == kCharm || int(abs(maPdgcode)/1000.) == kCharm ) {

    for (Int_t i=0; i<fNparents; i++){
      if (abs(maPdgcode)==fParentSelect[0][i]){
        isFinalOpenCharm = kTRUE;
      }
    }
    if (!isFinalOpenCharm) return -1;

    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<100; i++){ // check back to the 100 generation older

      Int_t jLabel = partMother->GetFirstMother();
      if (jLabel == -1){
        origin = kDirectCharm;
        return origin;
      }
      if (jLabel < 0){ // safety protection
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }

      // if there is an ancester
      TParticle* grandMa = fStack->Particle(jLabel);
      Int_t grandMaPDG = grandMa->GetPdgCode();

      for (Int_t j=0; j<fNparents; j++){
        if (abs(grandMaPDG)==fParentSelect[1][j]){
          origin = kBeautyCharm;
          return origin;
        }
      }

      partMother = grandMa;
    } // end of iteration 
  } // end of if
  else if ( int(abs(maPdgcode)/100.) == kBeauty || int(abs(maPdgcode)/1000.) == kBeauty ) {
    for (Int_t i=0; i<fNparents; i++){
      if (abs(maPdgcode)==fParentSelect[1][i]){
        origin = kDirectBeauty;
        return origin;
      }
    }
  } // end of if

  //============ gamma ================
  else if ( abs(maPdgcode) == 22 ) {
    origin = kGamma;

    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<100; i++){ // check back to the 100 generation older

      Int_t jLabel = partMother->GetFirstMother();
      if (jLabel == -1){
        origin = kGamma;
        return origin;
      }
      if (jLabel < 0){ // safety protection
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }

      // if there is an ancester
      TParticle* grandMa = fStack->Particle(jLabel);
      Int_t grandMaPDG = grandMa->GetPdgCode();

      for (Int_t j=0; j<fNparents; j++){
        if (abs(grandMaPDG)==fParentSelect[1][j]){
          origin = kBeautyGamma;
          return origin;
        }
      }

      partMother = grandMa;
    } // end of iteration 

    return origin;
  } // end of if

  //============ pi0 ================
  else if ( abs(maPdgcode) == 111 ) {
    origin = kPi0;

    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<100; i++){ // check back to the 100 generation older

      Int_t jLabel = partMother->GetFirstMother();
      if (jLabel == -1){
        origin = kPi0;
        return origin;
      }
      if (jLabel < 0){ // safety protection
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }

      // if there is an ancester
      TParticle* grandMa = fStack->Particle(jLabel);
      Int_t grandMaPDG = grandMa->GetPdgCode();

      for (Int_t j=0; j<fNparents; j++){
        if (abs(grandMaPDG)==fParentSelect[1][j]){
          origin = kBeautyPi0;
          return origin;
        }
      }

      partMother = grandMa;
    } // end of iteration 

    return origin;
  } // end of if

  else {
    origin = kElse;
    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<100; i++){ // check back to the 100 generation older

      Int_t jLabel = partMother->GetFirstMother();
      if (jLabel == -1){
        origin = kElse;
        return origin;
      }
      if (jLabel < 0){ // safety protection
        AliDebug(1, "Stack label is negative, return\n");
        return -1;
      }

      // if there is an ancester
      TParticle* grandMa = fStack->Particle(jLabel);
      Int_t grandMaPDG = grandMa->GetPdgCode();

      for (Int_t j=0; j<fNparents; j++){
        if (abs(grandMaPDG)==fParentSelect[1][j]){
          origin = kBeautyElse;
          return origin;
        }
      }

      partMother = grandMa;
    } // end of iteration 
  }

  return origin;
}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetPDG(AliVTrack *track)
{
  //
  // get KF particle input pdg for mass hypothesis
  //

  Int_t pdgCode=-1;

  if (fUseMCPID && HasMCData()){
    pdgCode = GetMCPDG(track);
  }
  else if(fESD1){
    Int_t pid=0;
    Double_t prob;
    GetESDPID((AliESDtrack*)track, pid, prob);
    switch(pid){
      case 0:  pdgCode = 11; break;
      case 1:  pdgCode = 13; break;
      case 2:  pdgCode = 211; break;
      case 3:  pdgCode = 321; break;
      case 4:  pdgCode = 2212; break;
      default: pdgCode = -1;
    }
  }
  else if(fAOD1){
    Int_t pid = ((AliAODTrack*)track)->GetMostProbablePID();
    switch(pid){
      case 0:  pdgCode = 11; break;
      case 1:  pdgCode = 13; break;
      case 2:  pdgCode = 211; break;
      case 3:  pdgCode = 321; break;
      case 4:  pdgCode = 2212; break;
      default: pdgCode = -1;
    }
  }

  return pdgCode;
}

//_______________________________________________________________________________________________
Int_t AliHFEsecVtx::GetMCPDG(AliVTrack *track)
{
  //
  // return mc pdg code
  //

  Int_t label = TMath::Abs(track->GetLabel());
  Int_t pdgCode; 

  if (IsAODanalysis()) {
    AliAODMCParticle *mcpart = (AliAODMCParticle*)fMCArray->At(label);
    if ( !mcpart ) return 0;
      pdgCode = mcpart->GetPdgCode();
  }
  else {
    TParticle* mcpart = fStack->Particle(label);
    if ( !mcpart ) return 0;
      pdgCode = mcpart->GetPdgCode();
  }

    return pdgCode;
}

//______________________________________________________________________________
void AliHFEsecVtx::GetESDPID(AliESDtrack *track, Int_t &recpid, Double_t &recprob) 
{
  //
  // calculate likehood for esd pid
  // 

  recpid = -1;
  recprob = -1;

  Int_t ipid=-1;
  Double_t max=0.;

  Double_t probability[5];

  // get probability of the diffenrent particle types
  track->GetESDpid(probability);

  // find most probable particle in ESD pid
  // 0:Electron - 1:Muon - 2:Pion - 3:Kaon - 4:Proton
  ipid = TMath::LocMax(5,probability);
  max = TMath::MaxElement(5,probability);

  recpid = ipid;
  recprob = max;
}

//_____________________________________________________________________________
void AliHFEsecVtx::AddHFEpairToArray(const AliHFEpairs* const pair)
{
  //
  // Add a HFE pair to the array
  //

  Int_t n = HFEpairs()->GetEntriesFast();
  if(n!=fNoOfHFEpairs)AliError(Form("fNoOfHFEpairs != HFEpairs()->GetEntriesFast %i != %i \n", fNoOfHFEpairs, n));
  new((*HFEpairs())[n]) AliHFEpairs(*pair);
}

//_____________________________________________________________________________
TClonesArray *AliHFEsecVtx::HFEpairs()
{
  //
  // Returns the list of HFE pairs
  //

  if (!fHFEpairs) {
      fHFEpairs = new TClonesArray("AliHFEpairs", 200);
  }
  return fHFEpairs;
}

//_____________________________________________________________________________
void AliHFEsecVtx::DeleteHFEpairs()
{
  //
  // Delete the list of HFE pairs
  //

  if (fHFEpairs) {
    fHFEpairs->Delete();
    //delete fHFEpairs;
  }
}

//_____________________________________________________________________________
void AliHFEsecVtx::InitHFEpairs()
{
  //
  // Initialization should be done before make all possible pairs for a given electron candidate
  //

  fNoOfHFEpairs = 0;
}

//_____________________________________________________________________________
void AliHFEsecVtx::AddHFEsecvtxToArray(const AliHFEsecVtxs* const secvtx)
{
  //
  // Add a HFE secondary vertex to the array
  //

  Int_t n = HFEsecvtxs()->GetEntriesFast();
  if(n!=fNoOfHFEsecvtxs)AliError(Form("fNoOfHFEsecvtxs != HFEsecvtxs()->GetEntriesFast %i != %i \n", fNoOfHFEsecvtxs, n));
  new((*HFEsecvtxs())[n]) AliHFEsecVtxs(*secvtx);
}

//_____________________________________________________________________________
TClonesArray *AliHFEsecVtx::HFEsecvtxs()
{
  //
  // Returns the list of HFE secvtx
  //

  if (!fHFEsecvtxs) {
      fHFEsecvtxs = new TClonesArray("AliHFEsecVtxs", 200);
  }
  return fHFEsecvtxs;
}

//_____________________________________________________________________________
void AliHFEsecVtx::DeleteHFEsecvtxs()
{
  //
  // Delete the list of HFE pairs
  //

  if (fHFEsecvtxs) {
    fHFEsecvtxs->Delete();
    //delete fHFEsecvtx;
  }
}

//_____________________________________________________________________________
void AliHFEsecVtx::InitHFEsecvtxs()
{
  //
  // Initialization should be done 
  //

  fNoOfHFEsecvtxs = 0;
}

//_______________________________________________________________________________________________
Bool_t AliHFEsecVtx::SingleTrackCut(AliESDtrack* track) const
{
  //if (track->Pt() < 1.0) return kFALSE;
  //if (TMath::Abs(track->Eta()) > 0.9) return kFALSE;
  //if (!(track->GetStatus() & AliESDtrack::kITSrefit)) return kFALSE;
  //if (!(track->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
  if (!(TESTBIT(track->GetITSClusterMap(),0))) return kFALSE; // ask hit on the first pixel layer
  //if (!(TESTBIT(track->GetITSClusterMap(),0) | TESTBIT(track->GetITSClusterMap(),1))) return kFALSE;

/*
  Float_t dcaR=-1;
  Float_t dcaZ=-1;
  track->GetImpactParameters(dcaR,dcaZ);
  if (TMath::Abs(TMath::Sqrt(dcaR*dcaR + dcaZ*dcaZ)) < 0.005) return kFALSE;
  if (TMath::Abs(TMath::Sqrt(dcaR*dcaR + dcaZ*dcaZ)) > 0.3) return kFALSE;
*/
  return kTRUE;
}
//____________________________________________________________
void AliHFEsecVtx::MakeContainer(){

  //
  // make container
  //

  const Int_t nDimPair=6;
  Int_t nBinPair[nDimPair] = {200, 500, 314, 2000, 2000, 11};
  //Int_t nBinPair[nDimPair] = {200, 500, 314, 2000, 2000, 11, 1000, 1000, 60, 60};
  const Double_t kInvmassmin = 0., kInvmassmax = 20.;
  const Double_t kKFChi2min = 0, kKFChi2max= 50;
  const Double_t kOpenanglemin = 0, kOpenanglemax = 3.14;
  const Double_t kSignedLxymin = -10, kSignedLxymax= 10;
  const Double_t kKFIPmin = -10, kKFIPmax= 10;
  const Double_t kPairCodemin = -1, kPairCodemax= 10;
  //const Double_t kDCAsigmin = 0, kDCAsigmax= 5;
  //const Double_t kPtmin = 0, kPtmax= 30;

  Double_t* binEdgesPair[nDimPair];
  for(Int_t ivar = 0; ivar < nDimPair; ivar++)
    binEdgesPair[ivar] = new Double_t[nBinPair[ivar] + 1];

  for(Int_t i=0; i<=nBinPair[0]; i++) binEdgesPair[0][i]=(Double_t)kInvmassmin + (kInvmassmax - kInvmassmin)/nBinPair[0]*(Double_t)i;
  for(Int_t i=0; i<=nBinPair[1]; i++) binEdgesPair[1][i]=(Double_t)kKFChi2min + (kKFChi2max - kKFChi2min)/nBinPair[1]*(Double_t)i;
  for(Int_t i=0; i<=nBinPair[2]; i++) binEdgesPair[2][i]=(Double_t)kOpenanglemin + (kOpenanglemax - kOpenanglemin)/nBinPair[2]*(Double_t)i;
  for(Int_t i=0; i<=nBinPair[3]; i++) binEdgesPair[3][i]=(Double_t)kSignedLxymin + (kSignedLxymax - kSignedLxymin)/nBinPair[3]*(Double_t)i;
  for(Int_t i=0; i<=nBinPair[4]; i++) binEdgesPair[4][i]=(Double_t)kKFIPmin + (kKFIPmax - kKFIPmin)/nBinPair[4]*(Double_t)i;
  for(Int_t i=0; i<=nBinPair[5]; i++) binEdgesPair[5][i]=(Double_t)kPairCodemin + (kPairCodemax - kPairCodemin)/nBinPair[5]*(Double_t)i;
  /*for(Int_t i=0; i<=nBinPair[6]; i++) binEdgesPair[6][i]=(Double_t)kDCAsigmin + (kDCAsigmax - kDCAsigmin)/nBinPair[6]*(Double_t)i;
  for(Int_t i=0; i<=nBinPair[7]; i++) binEdgesPair[7][i]=binEdgesPair[6][i];
  for(Int_t i=0; i<=nBinPair[8]; i++) binEdgesPair[8][i]=(Double_t)kPtmin + (kPtmax - kPtmin)/nBinPair[8]*(Double_t)i;
  for(Int_t i=0; i<=nBinPair[9]; i++) binEdgesPair[9][i]=binEdgesPair[8][i];*/

  fPairQA = new THnSparseF("pairQA", "QA for Pair; invmass[GeV/c^2]; KF chi2; opening angle; signed Lxy; KF ip; pair code", nDimPair, nBinPair);
  //fPairQA = new THnSparseF("pairQA", "QA for Pair; invmass[GeV/c^2]; KF chi2; opening angle; signed Lxy; KF ip; pair code; dca sig trk1; dca sig trk2; pt trk1; pt trk2 ", nDimPair, nBinPair);
  for(Int_t idim = 0; idim < nDimPair; idim++){
    fPairQA->SetBinEdges(idim, binEdgesPair[idim]);
  }

  fSecVtxList->AddAt(fPairQA,0);

  const Int_t nDimSecvtx=6;
  Double_t* binEdgesSecvtx[nDimSecvtx];
  Int_t nBinSecvtx[nDimSecvtx] = {200, 500, 2000, 2000, 11, 3};
  const Double_t kNtrksmin = 0, kNtrksmax= 3;
  for(Int_t ivar = 0; ivar < nDimSecvtx; ivar++)
    binEdgesSecvtx[ivar] = new Double_t[nBinSecvtx[ivar] + 1];

  for(Int_t i=0; i<=nBinSecvtx[0]; i++) binEdgesSecvtx[0][i]=binEdgesPair[0][i];
  for(Int_t i=0; i<=nBinSecvtx[1]; i++) binEdgesSecvtx[1][i]=binEdgesPair[1][i];
  for(Int_t i=0; i<=nBinSecvtx[2]; i++) binEdgesSecvtx[2][i]=binEdgesPair[3][i];
  for(Int_t i=0; i<=nBinSecvtx[3]; i++) binEdgesSecvtx[3][i]=binEdgesPair[4][i];
  for(Int_t i=0; i<=nBinSecvtx[4]; i++) binEdgesSecvtx[4][i]=binEdgesPair[5][i];
  for(Int_t i=0; i<=nBinSecvtx[5]; i++) binEdgesSecvtx[5][i]=(Double_t)kNtrksmin + (kNtrksmax - kNtrksmin)/nBinSecvtx[5]*(Double_t)i;

  fSecvtxQA = new THnSparseF("secvtxQA", "QA for Secvtx; invmass[GeV/c^2]; KF chi2; signed Lxy; KF ip; pair code; n tracks ", nDimSecvtx, nBinSecvtx);
  for(Int_t idim = 0; idim < nDimSecvtx; idim++){
    fSecvtxQA->SetBinEdges(idim, binEdgesSecvtx[idim]);
  }

  fSecVtxList->AddAt(fSecvtxQA,1);
  for(Int_t ivar = 0; ivar < nDimPair; ivar++)
    delete binEdgesPair[ivar];
  for(Int_t ivar = 0; ivar < nDimSecvtx; ivar++)
    delete binEdgesSecvtx[ivar];
}
