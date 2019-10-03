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
// Get improved dca info. by ITS upgrade implemented by the ALICE Heavy Flavour Electron Group
//
// Authors:
//   MinJung Kweon <minjung@physi.uni-heidelberg.de>
//
#include <TBits.h>
#include <TClass.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TString.h>
#include <TMath.h>
#include <TGraph.h>
#include <TFile.h>

#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVParticle.h"
#include "AliVertexerTracks.h"
#include "AliVVertex.h"
#include "AliExternalTrackParam.h"
#include "AliMCEvent.h"

#include "AliHFEsmearDCA.h"

ClassImp(AliHFEsmearDCA)


//______________________________________________________
AliHFEsmearDCA::AliHFEsmearDCA(const Char_t */*name*/, const char *resfileCurURI, const char *resfileUpgURI, const Char_t */*title*/):
   fEvent(NULL),
   fMCEvent(NULL),
   fD0ZResPCur  (0),
   fD0ZResKCur  (0),
   fD0ZResPiCur (0),
   fD0ZResECur  (0),
   fD0RPResPCur (0),
   fD0RPResKCur (0),
   fD0RPResPiCur(0),
   fD0RPResECur (0),
   fPt1ResPCur  (0),
   fPt1ResKCur  (0),
   fPt1ResPiCur (0),
   fPt1ResECur  (0),
   fD0ZResPUpg  (0),
   fD0ZResKUpg  (0),
   fD0ZResPiUpg (0),
   fD0ZResEUpg  (0),
   fD0RPResPUpg (0),
   fD0RPResKUpg (0),
   fD0RPResPiUpg(0),
   fD0RPResEUpg (0),
   fPt1ResPUpg  (0),
   fPt1ResKUpg  (0),
   fPt1ResPiUpg (0),
   fPt1ResEUpg  (0)
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
  fD0RPResECur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0RPResE")));
  fD0ZResPCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResP"  )));
  fD0ZResKCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResK"  )));
  fD0ZResPiCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResPi" )));
  fD0ZResECur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("D0ZResE")));
  fPt1ResPCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResP"  )));
  fPt1ResKCur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResK"  )));
  fPt1ResPiCur =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResPi" )));
  fPt1ResECur  =new TGraph(*static_cast<TGraph*>(resfileCur->Get("Pt1ResE" )));
  fD0RPResPCur ->SetName("D0RPResPCur" );
  fD0RPResKCur ->SetName("D0RPResKCur" );
  fD0RPResPiCur->SetName("D0RPResPiCur");
  fD0RPResECur ->SetName("D0RPResECur");
  fD0ZResPCur  ->SetName("D0ZResPCur"  ); 
  fD0ZResKCur  ->SetName("D0ZResKCur"  );
  fD0ZResPiCur ->SetName("D0ZResPiCur" );
  fD0ZResECur  ->SetName("D0ZResECur" );
  fPt1ResPCur  ->SetName("Pt1ResPCur"  );
  fPt1ResKCur  ->SetName("Pt1ResKCur"  );
  fPt1ResPiCur ->SetName("Pt1ResPiCur" );
  fPt1ResECur  ->SetName("Pt1ResECur" );
  delete resfileCur;
  TFile *resfileUpg=TFile::Open(resfileUpgURI );
  fD0RPResPUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResP" )));
  fD0RPResKUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResK" )));
  fD0RPResPiUpg=new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResPi")));
  fD0RPResEUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0RPResE")));
  fD0ZResPUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResP"  )));
  fD0ZResKUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResK"  )));
  fD0ZResPiUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResPi" )));
  fD0ZResEUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("D0ZResE" )));
  fPt1ResPUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResP"  )));
  fPt1ResKUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResK"  )));
  fPt1ResPiUpg =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResPi" )));
  fPt1ResEUpg  =new TGraph(*static_cast<TGraph*>(resfileUpg->Get("Pt1ResE" )));
  fD0RPResPUpg ->SetName("D0RPResPUpg" );
  fD0RPResKUpg ->SetName("D0RPResKUpg" );
  fD0RPResPiUpg->SetName("D0RPResPiUpg");
  fD0RPResEUpg ->SetName("D0RPResEUpg");
  fD0ZResPUpg  ->SetName("D0ZResPUpg"  );
  fD0ZResKUpg  ->SetName("D0ZResKUpg"  );
  fD0ZResPiUpg ->SetName("D0ZResPiUpg" );
  fD0ZResEUpg  ->SetName("D0ZResEUpg" );
  fPt1ResPUpg  ->SetName("Pt1ResPUpg"  );
  fPt1ResKUpg  ->SetName("Pt1ResKUpg"  );
  fPt1ResPiUpg ->SetName("Pt1ResPiUpg" );
  fPt1ResEUpg  ->SetName("Pt1ResEUpg" );
  delete resfileUpg;

}

//______________________________________________________
AliHFEsmearDCA::AliHFEsmearDCA(const AliHFEsmearDCA &c):
   TObject(c),
   fEvent(c.fEvent),
   fMCEvent(c.fMCEvent),
   fD0ZResPCur  (c.fD0ZResPCur),
   fD0ZResKCur  (c.fD0ZResKCur),
   fD0ZResPiCur (c.fD0ZResPiCur),
   fD0ZResECur  (c.fD0ZResECur),
   fD0RPResPCur (c.fD0RPResPCur),
   fD0RPResKCur (c.fD0RPResKCur),
   fD0RPResPiCur(c.fD0RPResPiCur),
   fD0RPResECur (c.fD0RPResECur),
   fPt1ResPCur  (c.fPt1ResPCur),
   fPt1ResKCur  (c.fPt1ResKCur),
   fPt1ResPiCur (c.fPt1ResPiCur),
   fPt1ResECur  (c.fPt1ResECur),
   fD0ZResPUpg  (c.fD0ZResPUpg),
   fD0ZResKUpg  (c.fD0ZResKUpg),
   fD0ZResPiUpg (c.fD0ZResPiUpg),
   fD0ZResEUpg  (c.fD0ZResEUpg),
   fD0RPResPUpg (c.fD0RPResPUpg),
   fD0RPResKUpg (c.fD0RPResKUpg),
   fD0RPResPiUpg(c.fD0RPResPiUpg),
   fD0RPResEUpg (c.fD0RPResEUpg),
   fPt1ResPUpg  (c.fPt1ResPUpg),
   fPt1ResKUpg  (c.fPt1ResKUpg),
   fPt1ResPiUpg (c.fPt1ResPiUpg),
   fPt1ResEUpg  (c.fPt1ResEUpg)
{
  //
  // Copy constructor
  // Performs a deep copy
  //
}

//______________________________________________________
AliHFEsmearDCA::~AliHFEsmearDCA(){
  //
  // Destructor
  //
}

//______________________________________________________
void AliHFEsmearDCA::SetRecEventInfo(const TObject *event){
  //
  // Set Virtual event an make a copy
  //
  if (!event) {
    AliError("Pointer to AliVEvent !");
    return;
  }
  TString className(event->ClassName());
  if (! (className.CompareTo("AliESDEvent")==0 || className.CompareTo("AliAODEvent")==0)) {
    AliError("argument must point to an AliESDEvent or AliAODEvent !");
    return ;
  }
  fEvent = (AliVEvent*) event;

}

//______________________________________________________
void AliHFEsmearDCA::GetImproveITSImpactParameters(AliVTrack *track, Double_t &dcaxyn, Double_t &dcaxyo, Double_t &dcaxySign, Double_t &dcaxySigo, Double_t &dcazn, Double_t &dcazo, Double_t &dcazSign, Double_t &dcazSigo){
  //
  // Get HFE impact parameter (with recalculated primary vertex)
  //
  if(!fEvent){
    AliDebug(1, "No Input event available\n");
    return;
  }
  const Double_t kBeampiperadius=3.;
  const Double_t kBeampiperadiusnew=1.9;
  TString type = track->IsA()->GetName();
  Double_t dcan[2]={-999.,-999.};
  Double_t covn[3]={-999.,-999.,-999.};
  Double_t dcao[2]={-999.,-999.};
  Double_t covo[3]={-999.,-999.,-999.};

  // Get reconstructed track parameters

  Double_t xyz[3],pxpypz[3],cv[21];
  track->GetXYZ(xyz);
  pxpypz[0]=track->Px();
  pxpypz[1]=track->Py(); 
  pxpypz[2]=track->Pz(); 
  track->GetCovarianceXYZPxPyPz(cv);
  Short_t sign = (Short_t)track->Charge();

  //AliExternalTrackParam et(track); //it doesn't work, i don't know why yet! so, I use the line below
  AliExternalTrackParam et(xyz,pxpypz,cv,sign);
  Double_t *param=const_cast<Double_t*>(et.GetParameter());

  const AliMCParticle *mc=static_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(track->GetLabel())));
  Double_t mcx[3];
  Double_t mcp[3];
  Double_t mccv[36]={0.};
  Short_t  mcc;
  mc->XvYvZv(mcx);
  mc->PxPyPz(mcp);
  mcc=mc->Charge();
  AliExternalTrackParam mct(mcx,mcp,mccv,mcc);
  const Double_t *parammc=mct.GetParameter();
  const AliVVertex *tmpvtx = fMCEvent->GetPrimaryVertex();

  AliVertex vtx(mcx,1.,1);
  AliVertex *tmpvtx2 = (AliVertex *)fMCEvent->GetPrimaryVertex();

  mct.PropagateToDCA(tmpvtx2,track->GetBz(),10.); //mjtest
  // Correct reference points and frames according to MC
  // TODO: B-Field correct?
  // TODO: failing propagation....
  et.PropagateToDCA(tmpvtx2,track->GetBz(),10.);
  et.Rotate(mct.GetAlpha());


  // Select appropriate smearing functions
  Double_t ptmc=TMath::Abs(mc->Pt());
  Double_t sd0rpn=0.;
  Double_t sd0zn =0.;
  Double_t spt1n =0.;
  Double_t sd0rpo=0.;
  Double_t sd0zo =0.;
  Double_t spt1o =0.;

  switch (mc->Particle()->GetPdgCode()) {
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
  case 11: case -11:
    sd0rpo=EvalGraph(fD0RPResECur,ptmc); 
    sd0zo =EvalGraph(fD0ZResECur ,ptmc);
    spt1o =EvalGraph(fPt1ResECur ,ptmc);
    sd0rpn=EvalGraph(fD0RPResEUpg,ptmc);
    sd0zn =EvalGraph(fD0ZResEUpg ,ptmc);
    spt1n =EvalGraph(fPt1ResEUpg ,ptmc);
    break;
  default:
    return;
  }

  //mj special to check resolution smeared by 10 %
  sd0rpn = sd0rpo * 1.1;
  sd0zn = sd0zo * 1.1;
  spt1n = spt1o;

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

  //printf("d0rpo= %lf d0rpmc= %lf dd0rpo= %lf dd0rpn= %lf\n",d0rpo,d0rpmc,dd0rpo,dd0rpn);
  // Copy the smeared parameters to the AOD track
  Double_t x[3];
  Double_t p[3];
  et.GetXYZ(x);
  et.GetPxPyPz(p);
//  printf("x[0]= %lf x[1]= %lf x[2]= %lf \n",x[0],x[1],x[2]);
//  track->SetPosition(x,kFALSE);
//  track->SetP(p,kTRUE);

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fEvent);
  if(!esdevent) return;
  //const AliVVertex *vtxESDSkip = esdevent->GetPrimaryVertex();
   // recalculate primary vertex for peri. and pp
   AliVertexerTracks vertexer(fEvent->GetMagneticField());
   vertexer.SetITSMode();
   vertexer.SetMinClusters(4);
   Int_t skipped[2];
   skipped[0] = track->GetID();
   vertexer.SetSkipTracks(1,skipped);
   vertexer.SetConstraintOn();
 //----diamond constraint-----------------------------
   Float_t diamondcovxy[3];
   esdevent->GetDiamondCovXY(diamondcovxy);
   Double_t pos[3]={esdevent->GetDiamondX(),esdevent->GetDiamondY(),0.};
   Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
   AliESDVertex *diamond = new AliESDVertex(pos,cov,1.,1);
   vertexer.SetVtxStart(diamond);
   delete diamond; diamond=NULL;
 //----------------------------------------------------
   const AliVVertex *vtxESDSkip = (const AliVVertex *) vertexer.FindPrimaryVertex(fEvent);   
   //vtxESDSkip = vertexer.FindPrimaryVertex(fEvent);

   // Getting the DCA
   // Propagation always done on a working copy to not disturb the track params of the original track
   AliESDtrack *esdtrack = NULL;
   if(!TString(track->IsA()->GetName()).CompareTo("AliESDtrack")){
     // Case ESD track: take copy constructor
     AliESDtrack *tmptrack = dynamic_cast<AliESDtrack *>(track);
     if(tmptrack) esdtrack = new AliESDtrack(*tmptrack);
   } else {
    // Case AOD track: take different constructor
    esdtrack = new AliESDtrack(track);
   }
   if(esdtrack){
    et.PropagateToDCA(tmpvtx, fEvent->GetMagneticField(), kBeampiperadiusnew, dcan, covn);
    //et.PropagateToDCA(vtxESDSkip, fEvent->GetMagneticField(), kBeampiperadiusnew, dcan, covn);
    esdtrack->PropagateToDCA(vtxESDSkip, fEvent->GetMagneticField(), kBeampiperadius, dcao, covo);
    dcaxyn = dcan[0];
    dcaxyo = dcao[0];
    dcazn = dcan[1];
    dcazo = dcao[1];
    Double_t resfactor=1., resfactorz=1.;
    if(sd0rpo) resfactor = sd0rpn/sd0rpo;
    if(sd0zo) resfactorz = sd0zn/sd0zo;
    if(covo[0]) dcaxySign = dcan[0]/(TMath::Sqrt(covo[0])*resfactor);
    if(covo[0]) dcaxySigo = dcao[0]/TMath::Sqrt(covo[0]);
    if(covo[2]) dcazSign = dcan[1]/(TMath::Sqrt(covo[2])*resfactorz);
    if(covo[2]) dcazSigo = dcao[1]/TMath::Sqrt(covo[2]);
    //if(dcaxyo) printf("pt= %lf resolusion xy= %lf dcaxyo= %lf dcaxyn=  %lf ratio= %lf\n",pt1mc,resfactor,dcaxyo,dcaxyn,dcaxyn/dcaxyo); 

    if(!covo[0]){
     dcaxySigo = -49.;
     dcaxySign = -49.;
    }
    if(!covo[2]){
     dcazSigo = -49.;
     dcazSign = -49.;
    }
   }
   delete esdtrack;
   //delete vtxESDSkip;

}

//______________________________________________________
Double_t AliHFEsmearDCA::EvalGraph(const TGraph *graph,Double_t x) const {
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
