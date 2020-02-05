/*************************************************************************
* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
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

#include <TObjArray.h>
#include <TClonesArray.h>
#include <TGraph.h>
#include <TFile.h>
#include <TList.h>
#include <TNtuple.h>

#include "AliVertex.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliExternalTrackParam.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliNeutralTrackParam.h"
#include "AliAnalysisTaskSEImproveITS3.h"
#include "AliAnalysisVertexingHF.h"

//
// Implementation of the "hybrid-approach" for ITS upgrade studies.
// The tastk smears the track parameters according to estimations
// from single-track upgrade studies. Afterwards it recalculates
// the parameters of the reconstructed decays.
//
// WARNING: This will affect all tasks in a train after this one
// (which is typically desired, though).
//

AliAnalysisTaskSEImproveITS3::AliAnalysisTaskSEImproveITS3() 
  :AliAnalysisTaskSE(),
   fD0ZResPCur  (0),
   fD0ZResKCur  (0),
   fD0ZResPiCur (0),
   fD0ZResDCur (0),
   fD0RPResPCur (0),
   fD0RPResKCur (0),
   fD0RPResPiCur(0),
   fD0RPResDCur(0),
   fPt1ResPCur  (0),
   fPt1ResKCur  (0),
   fPt1ResPiCur (0),
   fPt1ResDCur (0),
   fD0ZResPUpg  (0),
   fD0ZResKUpg  (0),
   fD0ZResPiUpg (0),
   fD0ZResDUpg (0),
   fD0RPResPUpg (0),
   fD0RPResKUpg (0),
   fD0RPResPiUpg(0),
   fD0RPResDUpg(0),
   fPt1ResPUpg  (0),
   fPt1ResKUpg  (0),
   fPt1ResPiUpg (0),
   fPt1ResDUpg (0),
/*   fD0ZResPCurSA  (0),
   fD0ZResKCurSA  (0),
   fD0ZResPiCurSA (0),
   fD0RPResPCurSA (0),
   fD0RPResKCurSA (0),
   fD0RPResPiCurSA(0),
   fPt1ResPCurSA  (0),
   fPt1ResKCurSA  (0),
   fPt1ResPiCurSA (0),
   fD0ZResPUpgSA  (0),
   fD0ZResKUpgSA  (0),
   fD0ZResPiUpgSA (0),
   fD0RPResPUpgSA (0),
   fD0RPResKUpgSA (0),
   fD0RPResPiUpgSA(0),
   fPt1ResPUpgSA  (0),
   fPt1ResKUpgSA  (0),
   fPt1ResPiUpgSA (0),
*/ fRunInVertexing(kFALSE),
   fImproveTracks(kTRUE),
   fUpdateSTCovMatrix(kTRUE),
   fUpdateSecVertCovMat(kTRUE),
   fDebugOutput (0),
   fDebugNtuple (0),
   fDebugVars   (0), 
   fNDebug      (0),
   fOnlyProcessFilledCand(kFALSE)
{
  //
  // Default constructor.
  //
}

AliAnalysisTaskSEImproveITS3::AliAnalysisTaskSEImproveITS3(const char *name,
                           const char *resfileCurURI,
                           const char *resfileUpgURI,
                           Bool_t isRunInVertexing,
			                     Int_t ndebug)
  :AliAnalysisTaskSE(name),
   fD0ZResPCur  (0),
   fD0ZResKCur  (0),
   fD0ZResPiCur (0),
   fD0ZResDCur (0),
   fD0RPResPCur (0),
   fD0RPResKCur (0),
   fD0RPResPiCur(0),
   fD0RPResDCur(0),
   fPt1ResPCur  (0),
   fPt1ResKCur  (0),
   fPt1ResPiCur (0),
   fPt1ResDCur (0),
   fD0ZResPUpg  (0),
   fD0ZResKUpg  (0),
   fD0ZResPiUpg (0),
   fD0ZResDUpg (0),
   fD0RPResPUpg (0),
   fD0RPResKUpg (0),
   fD0RPResPiUpg(0),
   fD0RPResDUpg(0),
   fPt1ResPUpg  (0),
   fPt1ResKUpg  (0),
   fPt1ResPiUpg (0),
   fPt1ResDUpg (0),
/*   fD0ZResPCurSA  (0),
   fD0ZResKCurSA  (0),
   fD0ZResPiCurSA (0),
   fD0RPResPCurSA (0),
   fD0RPResKCurSA (0),
   fD0RPResPiCurSA(0),
   fPt1ResPCurSA  (0),
   fPt1ResKCurSA  (0),
   fPt1ResPiCurSA (0),
   fD0ZResPUpgSA  (0),
   fD0ZResKUpgSA  (0),
   fD0ZResPiUpgSA (0),
   fD0RPResPUpgSA (0),
   fD0RPResKUpgSA (0),
   fD0RPResPiUpgSA(0),
   fPt1ResPUpgSA  (0),
   fPt1ResKUpgSA  (0),
   fPt1ResPiUpgSA (0),
  */ fRunInVertexing(isRunInVertexing),
   fImproveTracks(kTRUE),
   fUpdateSTCovMatrix(kTRUE),
   fUpdateSecVertCovMat(kTRUE),
   fDebugOutput (0),
   fDebugNtuple (0),
   fDebugVars   (0),
   fNDebug      (ndebug),
   fOnlyProcessFilledCand(kFALSE)
{
  //
  // Constructor to be used to create the task.
  // The the URIs specify the resolution files to be used. 
  // They are expected to contain TGraphs with the resolutions
  // for the current and the upgraded ITS (see code for details).
  // One may also specify for how many tracks debug information
  // is written to the output.
  //
  TFile *resfileCur=TFile::Open(resfileCurURI);
  fD0RPResPCur  = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("D0RPResP" )))->Clone("fD0RPResPCur"));
  fD0RPResKCur  = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("D0RPResK" )))->Clone("fD0RPResKCur"));
  fD0RPResPiCur = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("D0RPResPi" )))->Clone("fD0RPResPiCur"));
  fD0RPResDCur = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("D0RPResD" )))->Clone("fD0RPResDCur"));
  fD0ZResPCur   = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("D0ZResP" )))->Clone("fD0ZResPCur"));
  fD0ZResKCur   = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("D0ZResK" )))->Clone("fD0ZResKCur"));
  fD0ZResPiCur  = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("D0ZResPi" )))->Clone("fD0ZResPiCur"));
  fD0ZResDCur  = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("D0ZResD" )))->Clone("fD0ZResDCur"));
  fPt1ResPCur   = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("Pt1ResP" )))->Clone("fPt1ResPCur"));
  fPt1ResKCur   = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("Pt1ResK" )))->Clone("fPt1ResKCur"));
  fPt1ResPiCur  = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("Pt1ResPi" )))->Clone("fPt1ResPiCur"));
  fPt1ResDCur  = static_cast<TGraph*>((static_cast<TGraph*>(resfileCur->Get("Pt1ResD" )))->Clone("fPt1ResDCur"));


/*  fD0RPResPCurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResPSA" )));
  fD0RPResKCurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResKSA" )));
  fD0RPResPiCurSA=new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResPiSA")));
  fD0ZResPCurSA  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResPSA"  )));
  fD0ZResKCurSA  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResKSA"  )));
  fD0ZResPiCurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResPiSA" )));
  fPt1ResPCurSA  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResPSA"  )));
  fPt1ResKCurSA  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResKSA"  )));
  fPt1ResPiCurSA =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResPiSA" )));
  fD0RPResPCurSA ->SetName("D0RPResPCurSA" );
  fD0RPResKCurSA ->SetName("D0RPResKCurSA" );
  fD0RPResPiCurSA->SetName("D0RPResPiCurSA");
  fD0ZResPCurSA  ->SetName("D0ZResPCurSA"  ); 
  fD0ZResKCurSA  ->SetName("D0ZResKCurSA"  );
  fD0ZResPiCurSA ->SetName("D0ZResPiCurSA" );
  fPt1ResPCurSA  ->SetName("Pt1ResPCurSA"  );
  fPt1ResKCurSA  ->SetName("Pt1ResKCurSA"  );
  fPt1ResPiCurSA ->SetName("Pt1ResPiCurSA" );
*/
  delete resfileCur;

  TFile *resfileUpg=TFile::Open(resfileUpgURI );
  fD0RPResPUpg  = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("D0RPResP" )))->Clone("D0RPResPUpg"));
  fD0RPResKUpg  = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("D0RPResK" )))->Clone("fD0RPResKUpg"));
  fD0RPResPiUpg = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("D0RPResPi" )))->Clone("fD0RPResPiUpg"));
  fD0RPResDUpg = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("D0RPResD" )))->Clone("fD0RPResDUpg"));
  fD0ZResPUpg   = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("D0ZResP" )))->Clone("fD0ZResPUpg"));
  fD0ZResKUpg   = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("D0ZResK" )))->Clone("fD0ZResKUpg"));
  fD0ZResPiUpg  = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("D0ZResPi" )))->Clone("fD0ZResPiUpg"));
  fD0ZResDUpg  = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("D0ZResD" )))->Clone("fD0ZResDUpg"));
  fPt1ResPUpg   = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("Pt1ResP" )))->Clone("fPt1ResPUpg"));
  fPt1ResKUpg   = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("Pt1ResK" )))->Clone("fPt1ResKUpg"));
  fPt1ResPiUpg  = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("Pt1ResPi" )))->Clone("fPt1ResPiUpg"));
  fPt1ResDUpg  = static_cast<TGraph*>((static_cast<TGraph*>(resfileUpg->Get("Pt1ResD" )))->Clone("fPt1ResDUpg"));
  
  /*fD0RPResPUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResPSA" )));
  fD0RPResKUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResKSA" )));
  fD0RPResPiUpgSA=new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResPiSA")));
  fD0ZResPUpgSA  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResPSA"  )));
  fD0ZResKUpgSA  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResKSA"  )));
  fD0ZResPiUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResPiSA" )));
  fPt1ResPUpgSA  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResPSA"  )));
  fPt1ResKUpgSA  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResKSA"  )));
  fPt1ResPiUpgSA =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResPiSA" )));
  fD0RPResPUpgSA ->SetName("D0RPResPUpgSA" );
  fD0RPResKUpgSA ->SetName("D0RPResKUpgSA" );
  fD0RPResPiUpgSA->SetName("D0RPResPiUpgSA");
  fD0ZResPUpgSA  ->SetName("D0ZResPUpgSA"  );
  fD0ZResKUpgSA  ->SetName("D0ZResKUpgSA"  );
  fD0ZResPiUpgSA ->SetName("D0ZResPiUpgSA" );
  fPt1ResPUpgSA  ->SetName("Pt1ResPUpgSA"  );
  fPt1ResKUpgSA  ->SetName("Pt1ResKUpgSA"  );
  fPt1ResPiUpgSA ->SetName("Pt1ResPiUpgSA" );
*/
  delete resfileUpg;

  DefineOutput(1,TList::Class());
}

AliAnalysisTaskSEImproveITS3::~AliAnalysisTaskSEImproveITS3() {
  //
  // Destructor.
  //
  delete fDebugOutput;
}

void AliAnalysisTaskSEImproveITS3::UserCreateOutputObjects() {
  //
  // Creation of user output objects.
  //
  fDebugOutput=new TList();
  fDebugOutput->SetOwner();
  fDebugNtuple=new TNtuple("fDebugNtuple","Smearing","pdg:ptmc:d0rpo:d0zo:pt1o:sd0rpo:sd0zo:spt1o:d0rpn:d0zn:pt1n:sd0rpn:sd0zn:spt1n:d0rpmc:d0zmc:pt1mc");
  fDebugVars=new Float_t[fDebugNtuple->GetNvar()];
  
  fDebugOutput->Add(fDebugNtuple );

  fDebugOutput->Add(fD0RPResPCur );
  fDebugOutput->Add(fD0RPResKCur );
  fDebugOutput->Add(fD0RPResPiCur);
  fDebugOutput->Add(fD0RPResDCur);
  fDebugOutput->Add(fD0ZResPCur  ); 
  fDebugOutput->Add(fD0ZResKCur  );
  fDebugOutput->Add(fD0ZResPiCur );
  fDebugOutput->Add(fD0ZResDCur );
  fDebugOutput->Add(fPt1ResPCur  );
  fDebugOutput->Add(fPt1ResKCur  );
  fDebugOutput->Add(fPt1ResPiCur );
  fDebugOutput->Add(fPt1ResDCur );
  fDebugOutput->Add(fD0RPResPUpg );
  fDebugOutput->Add(fD0RPResKUpg );
  fDebugOutput->Add(fD0RPResPiUpg);
  fDebugOutput->Add(fD0RPResDUpg);
  fDebugOutput->Add(fD0ZResPUpg  );
  fDebugOutput->Add(fD0ZResKUpg  );
  fDebugOutput->Add(fD0ZResPiUpg );
  fDebugOutput->Add(fD0ZResDUpg );
  fDebugOutput->Add(fPt1ResPUpg  );
  fDebugOutput->Add(fPt1ResKUpg  );
  fDebugOutput->Add(fPt1ResPiUpg );
  fDebugOutput->Add(fPt1ResDUpg );

  PostData(1,fDebugOutput);
}

void AliAnalysisTaskSEImproveITS3::UserExec(Option_t*) {
  //
  // The event loop
  //
  AliAODEvent *ev=0x0;
  if(!fRunInVertexing) {
    ev=dynamic_cast<AliAODEvent*>(InputEvent());
  } else {
    if(AODEvent() && IsStandardAOD()) ev = dynamic_cast<AliAODEvent*> (AODEvent());
  }  
  if(!ev) return;
  Double_t bz=ev->GetMagneticField();

  // Smear all tracks
  TClonesArray *mcs=static_cast<TClonesArray*>(ev->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
  if (!mcs) return;
  
  // first loop on candidates to fill them in case of reduced AODs
  // this is done to have the same behaviour of the improver with full (pp, p-Pb) and recuced (Pb-Pb) candidates
  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
  
  // D0->Kpi
  TClonesArray *array2Prong=static_cast<TClonesArray*>(ev->GetList()->FindObject("D0toKpi"));
  if (array2Prong && !fOnlyProcessFilledCand) {
    Int_t ncand = array2Prong->GetEntriesFast();
    for (Int_t icand=0;icand<ncand;++icand) {
      AliAODRecoDecayHF2Prong *decay=static_cast<AliAODRecoDecayHF2Prong*>(array2Prong->At(icand));
      vHF->FillRecoCand(ev,(AliAODRecoDecayHF2Prong*)decay);
    }
  }
  // Dstar->Kpipi
  TClonesArray *arrayCascade=static_cast<TClonesArray*>(ev->GetList()->FindObject("Dstar"));
  if (arrayCascade && !fOnlyProcessFilledCand) {
    Int_t ncand = arrayCascade->GetEntriesFast();
    for (Int_t icand=0;icand<ncand;++icand) {
      AliAODRecoCascadeHF *decayDstar=static_cast<AliAODRecoCascadeHF*>(arrayCascade->At(icand));
      vHF->FillRecoCasc(ev,((AliAODRecoCascadeHF*)decayDstar),kTRUE);
    }
  }
  // Three prong
  TClonesArray *array3Prong=static_cast<TClonesArray*>(ev->GetList()->FindObject("Charm3Prong"));
  if (array3Prong && !fOnlyProcessFilledCand) {
    Int_t ncand = array3Prong->GetEntriesFast();
    for (Int_t icand=0;icand<ncand;++icand) {
      AliAODRecoDecayHF3Prong *decay=static_cast<AliAODRecoDecayHF3Prong*>(array3Prong->At(icand));
      vHF->FillRecoCand(ev,(AliAODRecoDecayHF3Prong*)decay);
    }
  }
  
  if (fImproveTracks) {
    for(Int_t itrack=0;itrack<ev->GetNumberOfTracks();++itrack) {
      AliAODTrack * trk = dynamic_cast<AliAODTrack*>(ev->GetTrack(itrack));
      if(!trk) AliFatal("Not a standard AOD");
      SmearTrack(trk,mcs);
    }
  }

  // TODO: recalculated primary vertex
  AliVVertex *primaryVertex=ev->GetPrimaryVertex();

  // Recalculate all candidates
  // D0->Kpi
  if (array2Prong) {
    Int_t ncand = array2Prong->GetEntriesFast();
    for (Int_t icand=0;icand<ncand;++icand) {
      AliAODRecoDecayHF2Prong *decay=static_cast<AliAODRecoDecayHF2Prong*>(array2Prong->At(icand));
      if(fOnlyProcessFilledCand && decay->GetIsFilled()!=2) continue;

      if(!vHF->FillRecoCand(ev,(AliAODRecoDecayHF2Prong*)decay))continue;

      // recalculate vertices
      AliVVertex *oldSecondaryVertex=decay->GetSecondaryVtx();

      AliExternalTrackParam et1; et1.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(0)));
      AliExternalTrackParam et2; et2.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(1)));

      TObjArray ta12;
     
      ta12.Add(&et1); ta12.Add(&et2); 
      AliESDVertex *v12 =RecalculateVertex(oldSecondaryVertex,&ta12 ,bz);
     

      // update secondary vertex
      Double_t pos[3];
      Double_t covpos[6];
      v12->GetXYZ(pos);
      v12->GetCovMatrix(covpos);
      decay->GetSecondaryVtx()->SetPosition(pos[0],pos[1],pos[2]);
      if(fUpdateSecVertCovMat) decay->GetSecondaryVtx()->SetCovMatrix(covpos);
      decay->GetSecondaryVtx()->SetChi2perNDF(v12->GetChi2toNDF()); 
     
      // update d0 
      Double_t d0z0[2],covd0z0[3];
      Double_t d0[2],d0err[2];
      et1.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d0[0]=d0z0[0];
      d0err[0] = TMath::Sqrt(covd0z0[0]);
      et2.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d0[1]=d0z0[0];
      d0err[1] = TMath::Sqrt(covd0z0[0]);   
      decay->Setd0Prongs(2,d0);
      decay->Setd0errProngs(2,d0err);
      // 


      Double_t xdummy=0.,ydummy=0.;
      Double_t dca;
      dca=et1.GetDCA(&et2,bz,xdummy,ydummy);
      decay->SetDCA(dca);

      
      Double_t px[2],py[2],pz[2];
      for (Int_t i=0;i<2;++i) {
        AliExternalTrackParam t;
        t.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(i)));
        t.PropagateToDCA(v12,bz,100.,d0z0,covd0z0);
        px[i]=t.Px();
        py[i]=t.Py();
        pz[i]=t.Pz();
      }
      decay->SetPxPyPzProngs(2,px,py,pz);
      delete v12;
    }
  }


  // Dstar->Kpipi
  if (arrayCascade) {
    Int_t ncand = arrayCascade->GetEntriesFast();
    for (Int_t icand=0;icand<ncand;++icand) {
      AliAODRecoCascadeHF *decayDstar=static_cast<AliAODRecoCascadeHF*>(arrayCascade->At(icand));
      if(fOnlyProcessFilledCand && decayDstar->GetIsFilled()!=2) continue;

      if(!vHF->FillRecoCasc(ev,((AliAODRecoCascadeHF*)decayDstar),kTRUE))continue;
      //Get D0 from D*
      AliAODRecoDecayHF2Prong* decay=(AliAODRecoDecayHF2Prong*)decayDstar->Get2Prong();
      
      // recalculate vertices
      //AliVVertex *oldSecondaryVertex=decay->GetSecondaryVtx();

      //soft pion
      AliExternalTrackParam et3; et3.CopyFromVTrack(static_cast<AliAODTrack*>(decayDstar->GetBachelor()));
       
      //track D0
      AliNeutralTrackParam *trackD0 = new AliNeutralTrackParam(decay);

      //!!!!TODO: covariance matrix
      
      // update d0 
      Double_t d0z0[2],covd0z0[3];
      Double_t d01[2],d01err[2];

      //the D*
      et3.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d01[0]=d0z0[0];
      d01err[0] = TMath::Sqrt(covd0z0[0]); 
      trackD0->PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d01[1]=d0z0[0];
      d01err[1] = TMath::Sqrt(covd0z0[0]);  
      decayDstar->Setd0Prongs(2,d01);
      decayDstar->Setd0errProngs(2,d01err);
        
      // delete v12;
      delete trackD0; trackD0=NULL;

       // a run for D*
      Double_t px1[2],py1[2],pz1[2];
      for (Int_t i=0;i<2;++i) {
      	const AliAODTrack *t1=static_cast<AliAODTrack*>(decayDstar->GetDaughter(i));
      	px1[i]=t1->Px();
      	py1[i]=t1->Py();
      	pz1[i]=t1->Pz();
      }
      decayDstar->SetPxPyPzProngs(2,px1,py1,pz1);
      
     }
  }


  // Three prong
  if (array3Prong) {
    Int_t ncand = array3Prong->GetEntriesFast();
    for (Int_t icand=0;icand<ncand;++icand) {
      AliAODRecoDecayHF3Prong *decay=static_cast<AliAODRecoDecayHF3Prong*>(array3Prong->At(icand));
      if(fOnlyProcessFilledCand && decay->GetIsFilled()!=2) continue;

      if(!vHF->FillRecoCand(ev,(AliAODRecoDecayHF3Prong*)decay))continue;
      
      // recalculate vertices
      AliVVertex *oldSecondaryVertex=decay->GetSecondaryVtx();
      AliExternalTrackParam et1; et1.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(0)));
      AliExternalTrackParam et2; et2.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(1)));
      AliExternalTrackParam et3; et3.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(2)));
      TObjArray ta123,ta12,ta23;
      ta123.Add(&et1);ta123.Add(&et2);ta123.Add(&et3);
      ta12. Add(&et1);ta12 .Add(&et2);
                      ta23 .Add(&et2);ta23 .Add(&et3);
      AliESDVertex *v123=RecalculateVertex(oldSecondaryVertex,&ta123,bz);
      AliESDVertex *v12 =RecalculateVertex(oldSecondaryVertex,&ta12 ,bz);
      AliESDVertex *v23 =RecalculateVertex(oldSecondaryVertex,&ta23 ,bz);

      // update secondary vertex
      Double_t pos[3];
      Double_t covpos[6];
      v123->GetXYZ(pos);
      v123->GetCovMatrix(covpos);
      decay->GetSecondaryVtx()->SetPosition(pos[0],pos[1],pos[2]);
      if(fUpdateSecVertCovMat) decay->GetSecondaryVtx()->SetCovMatrix(covpos);
      decay->GetSecondaryVtx()->SetChi2perNDF(v123->GetChi2toNDF()); 

      // update d0 for all progs
      Double_t d0z0[2],covd0z0[3];
      Double_t d0[3],d0err[3];
      et1.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d0[0]=d0z0[0];
      d0err[0] = TMath::Sqrt(covd0z0[0]);
      et2.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d0[1]=d0z0[0];
      d0err[1] = TMath::Sqrt(covd0z0[0]);
      et3.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d0[2]=d0z0[0];
      d0err[2] = TMath::Sqrt(covd0z0[0]);
      decay->Setd0Prongs   (3,d0   );
      decay->Setd0errProngs(3,d0err);
      // TODO: setter missing

      // update dca for prong combinations
      Double_t xdummy=0.,ydummy=0.;
      Double_t dca[3];
      dca[0]=et1.GetDCA(&et2,bz,xdummy,ydummy);
      dca[1]=et3.GetDCA(&et2,bz,xdummy,ydummy);
      dca[2]=et1.GetDCA(&et3,bz,xdummy,ydummy);
      decay->SetDCAs(3,dca);
      
      //update sigmavertex = dispersion
      Float_t sigmaV=v123->GetDispersion();
      decay->SetSigmaVert(sigmaV);
      // update dist12 and dist23
      primaryVertex->GetXYZ(pos);
      decay->SetDist12toPrim(TMath::Sqrt((v12->GetX()-pos[0])*(v12->GetX()-pos[0])
                                        +(v12->GetY()-pos[1])*(v12->GetY()-pos[1])
                                        +(v12->GetZ()-pos[2])*(v12->GetZ()-pos[2])));
      decay->SetDist23toPrim(TMath::Sqrt((v23->GetX()-pos[0])*(v23->GetX()-pos[0])
                                        +(v23->GetY()-pos[1])*(v23->GetY()-pos[1])
                                        +(v23->GetZ()-pos[2])*(v23->GetZ()-pos[2])));
 

      Double_t px[3],py[3],pz[3];
      for (Int_t i=0;i<3;++i) {
        AliExternalTrackParam t;
        t.CopyFromVTrack(static_cast<AliAODTrack*>(decay->GetDaughter(i)));
        t.PropagateToDCA(v123,bz,100.,d0z0,covd0z0);
        px[i]=t.Px();
        py[i]=t.Py();
        pz[i]=t.Pz();
      }
      decay->SetPxPyPzProngs(3,px,py,pz);

      delete v123;delete v12;delete v23;
    }
  }
  delete vHF;
}

void AliAnalysisTaskSEImproveITS3::SmearTrack(AliAODTrack *track,const TClonesArray *mcs) {
  // Early exit, if this track has nothing in common with the ITS
  if (!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)))
    return;

  // Check if the track was already "improved" (this is done with a trick using layer 7 (ie the 8th))
  if (TESTBIT(track->GetITSClusterMap(),7)) return;
  //


  // Get reconstructed track parameters
  AliExternalTrackParam et; et.CopyFromVTrack(track);
  Double_t *param=const_cast<Double_t*>(et.GetParameter());
  // Get covariance
  Double_t *covar=const_cast<Double_t*>(et.GetCovariance());
    

  // Get MC info
  Int_t imc=track->GetLabel();
  if (imc<=0) return;
  const AliAODMCParticle *mc=static_cast<AliAODMCParticle*>(mcs->At(imc));
  Double_t mcx[3];
  Double_t mcp[3];
  Double_t mccv[36]={0.};
  Short_t  mcc;
  mc->XvYvZv(mcx);
  mc->PxPyPz(mcp);
  mcc=mc->Charge();
  AliExternalTrackParam mct(mcx,mcp,mccv,mcc);
  const Double_t *parammc=mct.GetParameter();
//TODO:  const Double_t *covermc=mct.GetCovariance();
  AliVertex vtx(mcx,1.,1);

  // Correct reference points and frames according to MC
  // TODO: B-Field correct?
  // TODO: failing propagation....
  et.PropagateToDCA(&vtx,track->GetBz(),10.);
  et.Rotate(mct.GetAlpha());

  // Select appropriate smearing functions
  Double_t ptmc=TMath::Abs(mc->Pt());
  Double_t sd0rpn=0.;
  Double_t sd0zn =0.;
  Double_t spt1n =0.;
  Double_t sd0rpo=0.;
  Double_t sd0zo =0.;
  Double_t spt1o =0.;
  switch (mc->GetPdgCode()) {
  case 2212: case -2212:
    sd0rpo=EvalGraph(ptmc,fD0RPResPCur/*,fD0RPResPCurSA*/);
    sd0zo =EvalGraph(ptmc,fD0ZResPCur/*,fD0ZResPCurSA*/);
    spt1o =EvalGraph(ptmc,fPt1ResPCur/*,fPt1ResPCurSA*/);
    sd0rpn=EvalGraph(ptmc,fD0RPResPUpg/*,fD0RPResPUpgSA*/);
    sd0zn =EvalGraph(ptmc,fD0ZResPUpg/*,fD0ZResPUpgSA*/);
    spt1n =EvalGraph(ptmc,fPt1ResPUpg/*,fPt1ResPUpgSA*/);
    break;
  case 321: case -321:
    sd0rpo=EvalGraph(ptmc,fD0RPResKCur/*,fD0RPResKCurSA*/);
    sd0zo =EvalGraph(ptmc,fD0ZResKCur/*,fD0ZResKCurSA*/);
    spt1o =EvalGraph(ptmc,fPt1ResKCur/*,fPt1ResKCurSA*/);
    sd0rpn=EvalGraph(ptmc,fD0RPResKUpg/*,fD0RPResKUpgSA*/);
    sd0zn =EvalGraph(ptmc,fD0ZResKUpg/*,fD0ZResKUpgSA*/);
    spt1n =EvalGraph(ptmc,fPt1ResKUpg/*,fPt1ResKUpgSA*/);
    break;
  case 211: case -211:
    sd0rpo=EvalGraph(ptmc,fD0RPResPiCur/*,fD0RPResPiCurSA*/);
    sd0zo =EvalGraph(ptmc,fD0ZResPiCur/*,fD0ZResPiCurSA*/);
    spt1o =EvalGraph(ptmc,fPt1ResPiCur/*,fPt1ResPiCurSA*/);
    sd0rpn=EvalGraph(ptmc,fD0RPResPiUpg/*,fD0RPResPiUpgSA*/);
    sd0zn =EvalGraph(ptmc,fD0ZResPiUpg/*,fD0ZResPiUpgSA*/);
    spt1n =EvalGraph(ptmc,fPt1ResPiUpg/*,fPt1ResPiUpgSA*/);
    break;
  case 1000010020: case -1000010020:
    sd0rpo=EvalGraph(ptmc,fD0RPResDCur/*,fD0RPResDCurSA*/);
    sd0zo =EvalGraph(ptmc,fD0ZResDCur/*,fD0ZResDCurSA*/);
    spt1o =EvalGraph(ptmc,fPt1ResDCur/*,fPt1ResDCurSA*/);
    sd0rpn=EvalGraph(ptmc,fD0RPResDUpg/*,fD0RPResDUpgSA*/);
    sd0zn =EvalGraph(ptmc,fD0ZResDUpg/*,fD0ZResDUpgSA*/);
    spt1n =EvalGraph(ptmc,fPt1ResDUpg/*,fPt1ResDUpgSA*/);
    break;
  default:
    return;
  }

  // Use the same units (i.e. cm and GeV/c)! TODO: pt!
  sd0rpo*=1.e-4;
  sd0zo *=1.e-4;
  sd0rpn*=1.e-4;
  sd0zn *=1.e-4;

  // Apply the smearing
  Double_t d0zo  =param  [1];
  Double_t d0zmc =parammc[1];
  Double_t d0rpo =param  [0];
  Double_t d0rpmc=parammc[0];
  Double_t pt1o  =param  [4];
  Double_t pt1mc =parammc[4];
  Double_t dd0zo =d0zo-d0zmc;
  Double_t dd0zn =dd0zo *(sd0zo >0. ? (sd0zn /sd0zo ) : 1.);
  Double_t d0zn  =d0zmc+dd0zn;
  Double_t dd0rpo=d0rpo-d0rpmc;
  Double_t dd0rpn=dd0rpo*(sd0rpo>0. ? (sd0rpn/sd0rpo) : 1.);
  Double_t d0rpn =d0rpmc+dd0rpn;
  Double_t dpt1o =pt1o-pt1mc;
  Double_t dpt1n =dpt1o *(spt1o >0. ? (spt1n /spt1o ) : 1.);
  Double_t pt1n  =pt1mc+dpt1n;
  param[0]=d0rpn;
  param[1]=d0zn ;
  param[4]=pt1n ;
  /*
  Double_t d0zoinsigma = 0.;
  if(covar[0] > 0.) d0zoinsigma = d0zo/TMath::Sqrt(covar[2]);
  Double_t d0rpoinsigma = 0.;
  if(covar[2] > 0.) d0rpoinsigma = d0rpo/TMath::Sqrt(covar[0]);
  */
  //update the covariance matix
  if(fUpdateSTCovMatrix){
        if(sd0rpo>0.)            covar[0]*=(sd0rpn/sd0rpo)*(sd0rpn/sd0rpo);//yy
        if(sd0zo>0. && sd0rpo>0.)covar[1]*=(sd0rpn/sd0rpo)*(sd0zn/sd0zo);//yz
        if(sd0zo>0.)             covar[2]*=(sd0zn/sd0zo)*(sd0zn/sd0zo);//zz
        if(sd0rpo>0.)            covar[3]*=(sd0rpn/sd0rpo);//yl
        if(sd0zo>0.)             covar[4]*=(sd0zn/sd0zo);//zl
        if(sd0rpo>0.)            covar[6]*=(sd0rpn/sd0rpo);//ysenT
        if(sd0zo>0.)             covar[7]*=(sd0zn/sd0zo);//zsenT
        if(sd0rpo>0. && spt1o>0.)covar[10]*=(sd0rpn/sd0rpo)*(spt1n/spt1o);//ypt
        if(sd0zo>0. && spt1o>0.) covar[11]*=(sd0zn/sd0zo)*(spt1n/spt1o);//zpt
        if(spt1o>0.)             covar[12]*=(spt1n/spt1o);//sinPhipt
        if(spt1o>0.)             covar[13]*=(spt1n/spt1o);//tanTpt
        if(spt1o>0.)             covar[14]*=(spt1n/spt1o)*(spt1n/spt1o);//ptpt
  }
  /*
  Double_t d0zninsigma = 0.;
  if(covar[0] > 0.) d0zninsigma = d0zn/TMath::Sqrt(covar[2]);
  Double_t d0rpninsigma = 0.;
  if(covar[2] > 0.) d0rpninsigma = d0rpn/TMath::Sqrt(covar[0]);
   */
    
  // Copy the smeared parameters to the AOD track
  Double_t x[3];
  Double_t p[3];
  et.GetXYZ(x);
  et.GetPxPyPz(p);
  Double_t cv[21];
  et.GetCovarianceXYZPxPyPz(cv);
    
  track->SetPosition(x,kFALSE);
  track->SetP(p,kTRUE);
  track->SetCovMatrix(cv);
    


  // Mark the track as "improved" with a trick (this is done with a trick using layer 7 (ie the 8th))
  UChar_t itsClusterMap = track->GetITSClusterMap();
  SETBIT(itsClusterMap,7);
  track->SetITSClusterMap(itsClusterMap);
  //

  // write out debug infos
  if (fDebugNtuple->GetEntries()<fNDebug) {
    Int_t idbg=0;
    fDebugVars[idbg++]=mc->GetPdgCode();
    fDebugVars[idbg++]=ptmc  ;
    fDebugVars[idbg++]=d0rpo ;
    fDebugVars[idbg++]=d0zo  ;
    fDebugVars[idbg++]=pt1o  ;
    fDebugVars[idbg++]=sd0rpo;
    fDebugVars[idbg++]=sd0zo ;
    fDebugVars[idbg++]=spt1o ;
    fDebugVars[idbg++]=d0rpn ;
    fDebugVars[idbg++]=d0zn  ;
    fDebugVars[idbg++]=pt1n  ;
    fDebugVars[idbg++]=sd0rpn;
    fDebugVars[idbg++]=sd0zn ;
    fDebugVars[idbg++]=spt1n ;
    fDebugVars[idbg++]=d0rpmc;
    fDebugVars[idbg++]=d0zmc ;
    fDebugVars[idbg++]=pt1mc ;
    fDebugNtuple->Fill(fDebugVars);
    PostData(1,fDebugOutput);
  }
}

AliESDVertex* AliAnalysisTaskSEImproveITS3::RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField) {
  //
  // Helper function to recalculate a vertex.
  //

  static UShort_t ids[]={1,2,3}; //TODO: unsave...
  AliVertexerTracks vertexer(bField);
  vertexer.SetVtxStart(old->GetX(),old->GetY(),old->GetZ());
  AliESDVertex *vertex=vertexer.VertexForSelectedTracks(tracks,ids);
  return vertex;
}

Double_t AliAnalysisTaskSEImproveITS3::EvalGraph(Double_t x,const TGraph *graph/*,const TGraph *graphSA*/) const {
  //
  // Evaluates a TGraph without linear extrapolation. Instead the last
  // valid point of the graph is used when out of range.
  // The function assumes an ascending order of X.
  //
  // TODO: find a pretty solution for this:
  Int_t    n   =graph->GetN();
  Double_t xmin=graph->GetX()[0  ];
  Double_t xmax=graph->GetX()[n-1];
  if (x<xmin) {
   /* if(!graphSA)*/ return graph->Eval(xmin);
  //  Double_t xminSA=graphSA->GetX()[0];
  //  if(x<xminSA) return graphSA->Eval(xminSA);
  //  return graphSA->Eval(x);
  }
  if (x>xmax) return graph->Eval(xmax);
  return graph->Eval(x);
}

ClassImp(AliAnalysisTaskSEImproveITS3);

