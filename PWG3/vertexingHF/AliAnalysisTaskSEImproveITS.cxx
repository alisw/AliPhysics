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
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisTaskSEImproveITS.h"

//
// Implementation of the "hybrid-approach" for ITS upgrade studies.
// The tastk smears the track parameters according to estimations
// from single-track upgrade studies. Afterwards it recalculates
// the parameters of the reconstructed decays.
//
// WARNING: This will affect all tasks in a train after this one
// (which is typically desired, though).
//

AliAnalysisTaskSEImproveITS::AliAnalysisTaskSEImproveITS() 
  :AliAnalysisTaskSE(),
   fD0ZResPCur  (0),
   fD0ZResKCur  (0),
   fD0ZResPiCur (0),
   fD0RPResPCur (0),
   fD0RPResKCur (0),
   fD0RPResPiCur(0),
   fPt1ResPCur  (0),
   fPt1ResKCur  (0),
   fPt1ResPiCur (0),
   fD0ZResPUpg  (0),
   fD0ZResKUpg  (0),
   fD0ZResPiUpg (0),
   fD0RPResPUpg (0),
   fD0RPResKUpg (0),
   fD0RPResPiUpg(0),
   fPt1ResPUpg  (0),
   fPt1ResKUpg  (0),
   fPt1ResPiUpg (0),
   fDebugOutput (0),
   fDebugNtuple (0),
   fDebugVars   (0), 
   fNDebug      (0)
{
  //
  // Default constructor.
  //
}

AliAnalysisTaskSEImproveITS::AliAnalysisTaskSEImproveITS(const char *name,
                           const char *resfileCurURI,
                           const char *resfileUpgURI,
                           Int_t ndebug)
  :AliAnalysisTaskSE(name),
   fD0ZResPCur  (0),
   fD0ZResKCur  (0),
   fD0ZResPiCur (0),
   fD0RPResPCur (0),
   fD0RPResKCur (0),
   fD0RPResPiCur(0),
   fPt1ResPCur  (0),
   fPt1ResKCur  (0),
   fPt1ResPiCur (0),
   fD0ZResPUpg  (0),
   fD0ZResKUpg  (0),
   fD0ZResPiUpg (0),
   fD0RPResPUpg (0),
   fD0RPResKUpg (0),
   fD0RPResPiUpg(0),
   fPt1ResPUpg  (0),
   fPt1ResKUpg  (0),
   fPt1ResPiUpg (0),
   fDebugOutput (0),
   fDebugNtuple (0),
   fDebugVars   (0),
   fNDebug      (ndebug)
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
  fD0RPResPCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResP" )));
  fD0RPResKCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResK" )));
  fD0RPResPiCur=new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResPi")));
  fD0ZResPCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResP"  )));
  fD0ZResKCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResK"  )));
  fD0ZResPiCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResPi" )));
  fPt1ResPCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResP"  )));
  fPt1ResKCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResK"  )));
  fPt1ResPiCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResPi" )));
  fD0RPResPCur ->SetName("D0RPResPCur" );
  fD0RPResKCur ->SetName("D0RPResKCur" );
  fD0RPResPiCur->SetName("D0RPResPiCur");
  fD0ZResPCur  ->SetName("D0ZResPCur"  ); 
  fD0ZResKCur  ->SetName("D0ZResKCur"  );
  fD0ZResPiCur ->SetName("D0ZResPiCur" );
  fPt1ResPCur  ->SetName("Pt1ResPCur"  );
  fPt1ResKCur  ->SetName("Pt1ResKCur"  );
  fPt1ResPiCur ->SetName("Pt1ResPiCur" );
  delete resfileCur;
  TFile *resfileUpg=TFile::Open(resfileUpgURI );
  fD0RPResPUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResP" )));
  fD0RPResKUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResK" )));
  fD0RPResPiUpg=new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResPi")));
  fD0ZResPUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResP"  )));
  fD0ZResKUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResK"  )));
  fD0ZResPiUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResPi" )));
  fPt1ResPUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResP"  )));
  fPt1ResKUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResK"  )));
  fPt1ResPiUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResPi" )));
  fD0RPResPUpg ->SetName("D0RPResPUpg" );
  fD0RPResKUpg ->SetName("D0RPResKUpg" );
  fD0RPResPiUpg->SetName("D0RPResPiUpg");
  fD0ZResPUpg  ->SetName("D0ZResPUpg"  );
  fD0ZResKUpg  ->SetName("D0ZResKUpg"  );
  fD0ZResPiUpg ->SetName("D0ZResPiUpg" );
  fPt1ResPUpg  ->SetName("Pt1ResPUpg"  );
  fPt1ResKUpg  ->SetName("Pt1ResKUpg"  );
  fPt1ResPiUpg ->SetName("Pt1ResPiUpg" );
  delete resfileUpg;

  DefineOutput(1,TNtuple::Class());
}

AliAnalysisTaskSEImproveITS::~AliAnalysisTaskSEImproveITS() {
  //
  // Destructor.
  //
  delete fDebugOutput;
}

void AliAnalysisTaskSEImproveITS::UserCreateOutputObjects() {
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
  fDebugOutput->Add(fD0ZResPCur  ); 
  fDebugOutput->Add(fD0ZResKCur  );
  fDebugOutput->Add(fD0ZResPiCur );
  fDebugOutput->Add(fPt1ResPCur  );
  fDebugOutput->Add(fPt1ResKCur  );
  fDebugOutput->Add(fPt1ResPiCur );
  fDebugOutput->Add(fD0RPResPUpg );
  fDebugOutput->Add(fD0RPResKUpg );
  fDebugOutput->Add(fD0RPResPiUpg);
  fDebugOutput->Add(fD0ZResPUpg  );
  fDebugOutput->Add(fD0ZResKUpg  );
  fDebugOutput->Add(fD0ZResPiUpg );
  fDebugOutput->Add(fPt1ResPUpg  );
  fDebugOutput->Add(fPt1ResKUpg  );
  fDebugOutput->Add(fPt1ResPiUpg );

  PostData(1,fDebugOutput);
}

void AliAnalysisTaskSEImproveITS::UserExec(Option_t*) {
  //
  // The event loop
  //
  AliAODEvent *ev=dynamic_cast<AliAODEvent*>(InputEvent());
  if(!ev) return;
  Double_t bz=ev->GetMagneticField();

  // Smear all tracks
  TClonesArray *mcs=static_cast<TClonesArray*>(ev->GetList()->FindObject(AliAODMCParticle::StdBranchName()));
  if (!mcs) return;
  for (Int_t itrack=0;itrack<ev->GetNumberOfTracks();++itrack)
    SmearTrack(ev->GetTrack(itrack),mcs);

  // TODO: recalculated primary vertex
  AliVVertex *primaryVertex=ev->GetPrimaryVertex();

  // Recalculate all candidates
  TClonesArray *array3Prong=static_cast<TClonesArray*>(ev->GetList()->FindObject("Charm3Prong"));
  if (array3Prong) {
    for (Int_t icand=0;icand<array3Prong->GetEntries();++icand) {
      AliAODRecoDecayHF3Prong *decay=static_cast<AliAODRecoDecayHF3Prong*>(array3Prong->At(icand));

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
      v123->GetXYZ(pos);
      decay->GetSecondaryVtx()->SetPosition(pos[0],pos[1],pos[2]);
      decay->GetSecondaryVtx()->SetChi2perNDF(v123->GetChi2toNDF()); 
      //TODO: covariance matrix

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

      // update dist12 and dist23
      primaryVertex->GetXYZ(pos);
      decay->SetDist12toPrim(TMath::Sqrt((v12->GetX()-pos[0])*(v12->GetX()-pos[0])
                                        +(v12->GetY()-pos[1])*(v12->GetY()-pos[1])
                                        +(v12->GetZ()-pos[2])*(v12->GetZ()-pos[2])));
      decay->SetDist23toPrim(TMath::Sqrt((v23->GetX()-pos[0])*(v23->GetX()-pos[0])
                                        +(v23->GetY()-pos[1])*(v23->GetY()-pos[1])
                                        +(v23->GetZ()-pos[2])*(v23->GetZ()-pos[2])));
 
      delete v123;delete v12;delete v23;

      Double_t px[3],py[3],pz[3];
      for (Int_t i=0;i<3;++i) {
        const AliAODTrack *t=static_cast<AliAODTrack*>(decay->GetDaughter(i));
        px[i]=t->Px();
        py[i]=t->Py();
        pz[i]=t->Pz();
      }
      decay->SetPxPyPzProngs(3,px,py,pz);
    }
  }
}

void AliAnalysisTaskSEImproveITS::SmearTrack(AliAODTrack *track,const TClonesArray *mcs) {
  // Early exit, if this track has nothing in common with the ITS
  if (!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)))
    return;

  // Get reconstructed track parameters
  AliExternalTrackParam et; et.CopyFromVTrack(track);
  Double_t *param=const_cast<Double_t*>(et.GetParameter());
//TODO:  Double_t *covar=const_cast<Double_t*>(et.GetCovariance());

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
    sd0rpo=EvalGraph(fD0RPResPCur,ptmc);
    sd0zo =EvalGraph(fD0ZResPCur ,ptmc);
    spt1o =EvalGraph(fPt1ResPCur ,ptmc);
    sd0rpn=EvalGraph(fD0RPResPUpg,ptmc);
    sd0zn =EvalGraph(fD0ZResPUpg ,ptmc);
    spt1n =EvalGraph(fPt1ResPUpg ,ptmc);
    break;
  case 321: case -321:
    sd0rpo=EvalGraph(fD0RPResKCur,ptmc); 
    sd0zo =EvalGraph(fD0ZResKCur ,ptmc);
    spt1o =EvalGraph(fPt1ResKCur ,ptmc);
    sd0rpn=EvalGraph(fD0RPResKUpg,ptmc);
    sd0zn =EvalGraph(fD0ZResKUpg ,ptmc);
    spt1n =EvalGraph(fPt1ResKUpg ,ptmc);
    break;
  case 211: case -211:
    sd0rpo=EvalGraph(fD0RPResPiCur,ptmc); 
    sd0zo =EvalGraph(fD0ZResPiCur ,ptmc);
    spt1o =EvalGraph(fPt1ResPiCur ,ptmc);
    sd0rpn=EvalGraph(fD0RPResPiUpg,ptmc);
    sd0zn =EvalGraph(fD0ZResPiUpg ,ptmc);
    spt1n =EvalGraph(fPt1ResPiUpg ,ptmc);
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

  // Copy the smeared parameters to the AOD track
  Double_t x[3];
  Double_t p[3];
  et.GetXYZ(x);
  et.GetPxPyPz(p);
  track->SetPosition(x,kFALSE);
  track->SetP(p,kTRUE);

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

AliESDVertex* AliAnalysisTaskSEImproveITS::RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField) {
  //
  // Helper function to recalculate a vertex.
  //

  static UShort_t ids[]={1,2,3}; //TODO: unsave...
  AliVertexerTracks vertexer(bField);
  vertexer.SetVtxStart(old->GetX(),old->GetY(),old->GetZ());
  AliESDVertex *vertex=vertexer.VertexForSelectedTracks(tracks,ids);
  return vertex;
}

Double_t AliAnalysisTaskSEImproveITS::EvalGraph(const TGraph *graph,Double_t x) const {
  //
  // Evaluates a TGraph without linear extrapolation. Instead the last
  // valid point of the graph is used when out of range.
  // The function assumes an ascending order of X.
  //

  // TODO: find a pretty solution for this:
  Int_t    n   =graph->GetN();
  Double_t xmin=graph->GetX()[0  ];
  Double_t xmax=graph->GetX()[n-1];
  if (x<xmin) return graph->Eval(xmin);
  if (x>xmax) return graph->Eval(xmax);
  return graph->Eval(x);
}

ClassImp(AliAnalysisTaskSEImproveITS);

