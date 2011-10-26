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
// Static function member which can be used in standalone cases
// especially as utils for AliTPCCosmicTrackfit
//
// detailed description can be found inside individual function
//
// grep "exitreport" in output log to check abnormal termination
//
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch

//

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TTreeStream.h>
#include <TVector3.h>

#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliTPCseed.h"
#include "AliTrackerBase.h"
#include "AliTrackPointArray.h"

#include "AliTPCCosmicUtils.h"

Double_t AliTPCCosmicUtils::Point2LineDist(const TVector3 p0, const TVector3 l1, const TVector3 l2)
{
  //
  //return distance of p0 to line (l2-l1)
  //

  const TVector3 va = p0 - l1;
  const TVector3 vb = l2 - l1;

  const TVector3 dd = va.Cross(vb);

  return dd.Mag()/vb.Mag();
}

Double_t AliTPCCosmicUtils::AngleInRange(Double_t phi)
{
  //
  //get the phi value (in rad) in -pi ~ pi, so that fabs() works naiively
  //
  const Double_t pmin = -TMath::Pi();
  const Double_t pmax = pmin+TMath::TwoPi();

  while(phi<pmin)
    phi+= TMath::TwoPi();

  while(phi>=pmax)
    phi-=  TMath::TwoPi();

  return phi;
}

Bool_t AliTPCCosmicUtils::RotateSafe(AliExternalTrackParam *trackPar, const Double_t rawalpha)
{
  //1. in AliExternalTrackParam::GetXYZ
  //r[0]=fX; r[1]=fP[0]; r[2]=fP[1];
  //return Local2GlobalPosition(r,fAlpha);
  //in AliVParticle::Local2GlobalMomentum 
  //Double_t cs=TMath::Cos(alpha), sn=TMath::Sin(alpha);
  //Double_t r=TMath::Sqrt((1. - p[1])*(1. + p[1]));
  //p[0]=pt*(r*cs - p[1]*sn); p[1]=pt*(p[1]*cs + r*sn); p[2]=pt*p[2];
  //which is cos(phi_local=Snp+fAlpha) > 0 always.  
  //2. also in Bool_t AliExternalTrackParam::Rotate(Double_t alpha) 
  // Double_t sf=fP2, cf=TMath::Sqrt((1.- fP2)*(1.+fP2));
  //in AliExternalTrackParam::Set
  //mom.RotateZ(-fAlpha);
  //fP[2] = TMath::Sin(mom.Phi());
  //since only sin is used for mom.Phi(), that is assuming cos(mom.Phi())>0

  //the range-changing in AliExternalTrackParam::Rotate not precise enough
  const Double_t aa = AngleInRange(rawalpha);

  const Double_t a0 = trackPar->GetAlpha();
  const Double_t p2 = trackPar->GetParameter()[2];

  //copied from AliExternalTrackParam::Rotate
  const Double_t ca=TMath::Cos(aa-a0), sa=TMath::Sin(aa-a0);
  const Double_t sf=p2, cf=TMath::Sqrt((1.- p2)*(1.+p2)); 

  if((cf*ca+sf*sa)<0) {
    return kFALSE;
  }

  return trackPar->Rotate(aa);

  /*
  Double_t xyz[3], pxpypz[3];
  trackPar->GetXYZ(xyz);
  trackPar->GetPxPyPz(pxpypz);
  const TVector3 pos0(xyz);
  const TVector3 mom0(pxpypz);
  TVector3 pos1, mom1, dpos, dmom;
  const Double_t eps = 1e-6;
  AliExternalTrackParam tmppar = (*trackPar);

  if(!tmppar.Rotate(aa)){
    return 0;
  }

  tmppar.GetXYZ(xyz);
  tmppar.GetPxPyPz(pxpypz);
  pos1.SetXYZ(xyz[0], xyz[1], xyz[2]);
  mom1.SetXYZ(pxpypz[0], pxpypz[1], pxpypz[2]);
  dpos = pos1-pos0;
  dmom = mom1-mom0;
  if(dpos.Mag()<eps && dmom.Mag()<eps){
    (*trackPar)=tmppar;
    return 1;
  }
  else 
    return 0;
  */

  /*
  //flip
  tmppar = (*trackPar);
  if(!tmppar.Rotate(aa+TMath::Pi())){
    return 0;
  }

  tmppar.GetXYZ(xyz);
  tmppar.GetPxPyPz(pxpypz);
  pos1.SetXYZ(xyz[0], xyz[1], xyz[2]);
  mom1.SetXYZ(pxpypz[0], pxpypz[1], pxpypz[2]);
  dpos = pos1-pos0;
  dmom = mom1-mom0;
  if(dpos.Mag()<eps && dmom.Mag()<eps){
    (*trackPar)=tmppar;
    return -1;
  }

  return 0;
  */
}

void AliTPCCosmicUtils::PrintTrackParam(const Int_t id, const AliExternalTrackParam * trackpar, const char *tag)
{
  //
  //print out TrackParam
  //
  Double_t xyz[3]={-999,-999,-999};
  trackpar->GetXYZ(xyz);
  Double_t pxpypz[3]={-999,-999,-999};
  trackpar->GetPxPyPz(pxpypz);

  printf("PrintTrackPar %2s %3d : %6.1f %5.1f : %6.6f %6.6f %6.6f :  %6.6f  : %8.3f %8.3f %8.3f : (%2d) %8.3f %11.6f %.2f\n", tag, id, trackpar->GetX(), sqrt(pow(xyz[0],2)+pow(xyz[1],2)), xyz[0], xyz[1], xyz[2], trackpar->Phi(), pxpypz[0], pxpypz[1], pxpypz[2], trackpar->Charge(), trackpar->Pt(), trackpar->P(), trackpar->GetAlpha()*TMath::RadToDeg());
}

void AliTPCCosmicUtils::DrawTracks(AliESDtrack *esdtrks[], const TString tag, const TString outputformat)
{
  //
  //draw esdtracks
  //
  const AliTPCseed *seeds[]={GetTPCseed(esdtrks[0]), GetTPCseed(esdtrks[1])};
  DrawSeeds(seeds, tag, outputformat);
}


void AliTPCCosmicUtils::DrawSeeds(const AliTPCseed * seeds[], const TString tag, const TString outputformat)
{
  //
  //draw seed and output to file
  //

  if(!seeds[0] || !seeds[1])
    return;

  TGraph *grsyx[]={new TGraph, new TGraph};
  TGraph *grsyz[]={new TGraph, new TGraph};

  for(Int_t itrk=0; itrk<2; itrk++){
    if(!seeds[itrk])
      continue;

    Int_t ncl = 0;
    for(Int_t irow=0; irow<fgkNRow; irow++){
      const AliTPCclusterMI * cls = seeds[itrk]->GetClusterPointer(irow);
      if(!cls)
        continue;

      Float_t xyz[3]={-999,-999,-999};
      cls->GetGlobalXYZ(xyz);
      //printf("Test %f %f %f \n", xyz[0], xyz[1], xyz[2]);
      grsyx[itrk]->SetPoint(ncl, xyz[0], xyz[1]);
      grsyz[itrk]->SetPoint(ncl, xyz[2], xyz[1]);
      ncl++;
    }
  }

  grsyx[0]->SetTitle(tag);
  grsyx[0]->SetMaximum(250);
  grsyx[0]->SetMinimum(-250);
  grsyx[0]->GetXaxis()->SetLimits(-250,250);
  grsyx[0]->GetXaxis()->SetTitle("X (cm)");
  grsyx[0]->GetYaxis()->SetTitle("Y (cm)");

  grsyx[0]->SetMarkerStyle(20);
  grsyx[0]->SetMarkerColor(kRed);
  grsyx[1]->SetMarkerStyle(22);
  grsyx[1]->SetMarkerColor(kBlue);

  grsyz[0]->SetTitle(tag);
  grsyz[0]->SetMaximum(250);
  grsyz[0]->SetMinimum(-250);
  grsyz[0]->GetXaxis()->SetLimits(-250,250);
  grsyz[0]->GetXaxis()->SetTitle("Z (cm)");
  grsyz[0]->GetYaxis()->SetTitle("Y (cm)");

  grsyz[0]->SetMarkerStyle(20);
  grsyz[0]->SetMarkerColor(kRed);
  grsyz[1]->SetMarkerStyle(22);
  grsyz[1]->SetMarkerColor(kBlue);

  TCanvas *cc=new TCanvas("cc","",1200,600);
  cc->Divide(2);
  cc->cd(1);
  for(int itrk=0; itrk<2; itrk++){
    grsyx[itrk]->SetMarkerSize(1);
    grsyx[itrk]->Draw(itrk?"lp same":"alp");
  }
  cc->cd(2);
  for(int itrk=0; itrk<2; itrk++){
    grsyz[itrk]->SetMarkerSize(1);
    grsyz[itrk]->Draw(itrk?"lp same":"alp");
  }

  gErrorIgnoreLevel = 1001;
  cc->Print(Form("drawTrack%s.%s", tag.Data(), outputformat.Data()));

  for(Int_t ii=0; ii<2; ii++){
    delete grsyx[ii];
    delete grsyz[ii];
  }
  delete cc;
}

//================================================================================================
AliTPCseed * AliTPCCosmicUtils::GetTPCseed(const AliESDtrack *esdtrack)
{
  //
  //Get TPC seeds from ESDfriendTrack
  //

  AliESDfriendTrack * friendtrk = (AliESDfriendTrack *) esdtrack->GetFriendTrack();
  if(!friendtrk){
    printf("exitreport AliTPCCosmicUtils::GetTPCseed no friend track!\n");
    exit(1);
  }
    
  TObject *calibObject=0x0;
  AliTPCseed *tseed = 0x0;
  for(Int_t l=0; (calibObject=friendtrk->GetCalibObject(l)); l++) {
    if( (tseed=dynamic_cast<AliTPCseed*>(calibObject)) )
      break;
  }
  return tseed;
}

//================================================================================================
AliExternalTrackParam *AliTPCCosmicUtils::MakeSeed(const AliTPCseed *tseed)
{
  //
  //make seed for propagation of TrackParam, using np = 3 outer clusters (separated by deltancls clusters) in TPCseed
  //

  const Int_t rowstart = fgkNRow-1;
  const Int_t rowstop = 0;  
  const Int_t drow = -1;

  //---
  const Int_t np = 3;
  AliTrackPoint *tpos[np];
  for(Int_t ii=0; ii<np; ii++)
    tpos[ii] = 0x0;

  const Float_t cov[6]={0,0,0, 0.01*0.01 ,0, 0.01*0.01};

  //---
  Int_t npos = 0;
  Int_t icl = 0;
  const Int_t deltancls = fgkNclsMin/2-1;
  Int_t oldcl = -deltancls;

  for(Int_t irow=rowstart; drow*irow<=drow*rowstop; irow+=drow){
    AliTPCclusterMI *cls = tseed->GetClusterPointer(irow);
    if(!cls) {
      continue;
    }

    if( icl == (oldcl+deltancls) ){
      Float_t txyz[3];
      cls->GetGlobalXYZ(txyz);
      tpos[npos++] = new AliTrackPoint(txyz, cov, 0);
      //printf("------ %d %f %f %f\n", npos, txyz[0], txyz[1], txyz[2]);

      oldcl = icl;
      if(npos==np) break;
    }
    icl++;
  }
  if(npos!=np){
    for(Int_t ii=0; ii<npos; ii++)
      delete tpos[ii];

    return 0x0;
  }

  AliExternalTrackParam * trackparam = AliTrackerBase::MakeSeed(*(tpos[0]), *(tpos[1]), *(tpos[2]));
  if(!trackparam || trackparam->Pt()==0){
    for(Int_t ii=0; ii<npos; ii++)
      delete tpos[ii];
    delete trackparam;

    return 0x0;
  }

  //------

  Double_t sxyz[3]={-999,-999,-999}, spxpypz[3]={-999,-999,-999};
  trackparam->GetXYZ(sxyz);
  trackparam->GetPxPyPz(spxpypz);
  Double_t scov[21];
  Int_t sign = trackparam->Charge();
  
  //reset covariance matrix -- necessary, otherwise bad fitting result: trackparam->GetCovariance()[ii] has problem: some are strange values
  for(Int_t ii=0; ii<21; ii++) {
    scov[ii]=0;
  }
  trackparam->Set(sxyz, spxpypz, scov, sign);

  for(Int_t ii=0; ii<np; ii++)
    delete tpos[ii];

  return trackparam;
}

//================================================================================================
void AliTPCCosmicUtils::IniCov(AliExternalTrackParam *trackPar, const Double_t ncl)
{
  //
  //initialize covariance matrix
  //

  const Double_t ksigma=5.;
  Double_t acov[16];
  for (Int_t i=0;i<15;i++)
    acov[i]=0;

  acov[0]=ksigma*ksigma;
  acov[2]=ksigma*ksigma;
  acov[5]=ksigma*ksigma;
  acov[9]=ksigma*ksigma;
  acov[14]=0.2*0.2;

  acov[5] = ksigma*ksigma/(ncl*ncl);
  acov[9] = ksigma*ksigma/(ncl*ncl);

  const Double_t resetcov = 4; 

  trackPar->ResetCovariance(resetcov);
  //the following helps a lot!!
  trackPar->AddCovariance(acov);
}

void AliTPCCosmicUtils::SingleFit(AliExternalTrackParam * trackInOld, AliExternalTrackParam * trackOutOld, const AliTPCseed *tseed, const Bool_t kinward,  const Int_t rowstartshift, const Int_t rowstep, const Double_t xmin, const Double_t xmax, Int_t &nfit, Int_t &nmiss, Double_t &pchi2, Double_t &lfit, TTreeSRedirector *debugstreamer)
{
  //
  //fit single track
  //
  //kinward is the true geometry of the track. Incomming track: 1; outgoing track: 0
  const Double_t inde  = kinward? -1 :  1;
  const Double_t outde = kinward?  1 : -1;
  
  //PrintTrackParam(9000, trackOutOld);

  AliExternalTrackParam trackOut = *trackOutOld;
  AliExternalTrackParam trackIn;
  Int_t ksite = -999;

  //nmiss is from the 2 FitKernel of the last iteration
  //nfit, pchi2 is from the last FitKernel of the last iteration

  Int_t rowouter = fgkNRow-1-rowstartshift;

  //so that when reversed, the same rows are read! important!!!
  Int_t rowinner = rowouter - rowouter/rowstep * rowstep;

  TVector3 gposStart;
  TVector3 gposStop;
  TVector3 dpos;
  lfit = 0;

  for(Int_t ii=0; ii<fgkNiter; ii++){
    nmiss = 0;

    gposStart.SetXYZ(-999,-999,-999);
    ksite = -999;
    trackIn  = trackOut;
    FitKernel(&trackIn,  tseed, rowouter, rowinner, -rowstep, xmin, xmax, inde , ksite, nfit, nmiss, pchi2, gposStart, gposStop, 0x0, kTRUE); 
    //PrintTrackParam(9010+ii, &trackIn);                

    dpos = gposStart-gposStop;
    lfit += dpos.Pt();

    //---------------------------

    nfit = 0;
    pchi2 = 0;                                        

    gposStart.SetXYZ(-999,-999,-999);
    ksite = -999;
    trackOut = trackIn;
    FitKernel(&trackOut, tseed, rowinner, rowouter,  rowstep, xmin, xmax, outde, ksite, nfit, nmiss, pchi2, gposStart, gposStop, (ii==fgkNiter-1 ? debugstreamer : 0x0), kTRUE);
    //PrintTrackParam(90020+ii, &trackOut);

    dpos = gposStart-gposStop;
    lfit += dpos.Pt();
  }

  lfit /= 2.*fgkNiter;

  if(trackInOld)
    (*trackInOld)  = trackIn;

  (*trackOutOld) = trackOut;
}

void AliTPCCosmicUtils::CombinedFit(AliExternalTrackParam *trackPars[],  const AliTPCseed *seeds[],  const Int_t rowstartshift, const Int_t rowstep, const Double_t xmin, const Double_t xmax, Int_t &nfit, Int_t &nmiss, Double_t &pchi2, Double_t &lfit, Double_t &vtxD, Double_t &vtxZ, TTreeSRedirector *debugstreamer)
{
  //
  //combined propagation
  //

  //seen from 1/pt_fit - 1/pt_gen, for differnt de:
  //u+d+: combUp good, combLow bad
  //u+d-: combUp good, combLow good
  //u-d-: combUp bad,  combLow good
  //u-d+: combUp bad,  combLow bad
  const Double_t upde  = 1;
  const Double_t lowde = -1;

  AliExternalTrackParam trackLow = *(trackPars[1]);
  AliExternalTrackParam trackUp;

  for(Int_t ii=0; ii<fgkNiter; ii++){
    lfit = 0;
    vtxD = 0;
    vtxZ = 0;
    Double_t tmpl = -999, tmpd = -999, tmpz = -999;

    //lower->upper
    trackUp = trackLow;
    SubCombined(&trackUp,  seeds, 1, 0, rowstartshift, rowstep, xmin, xmax, upde,  nfit, nmiss, pchi2, tmpl, tmpd, tmpz);
    lfit += tmpl;
    vtxD += tmpd;
    vtxZ += tmpz;
    
    //upper->lower
    trackLow = trackUp;    
    SubCombined(&trackLow, seeds, 0, 1, rowstartshift, rowstep, xmin, xmax, lowde, nfit, nmiss, pchi2, tmpl, tmpd, tmpz, (ii==fgkNiter-1? debugstreamer : 0x0));
    lfit += tmpl;
    vtxD += tmpd;
    vtxZ += tmpz;
  }

  *(trackPars[0]) = trackUp;
  *(trackPars[1]) = trackLow;

  //only last iteration used
  lfit /= 2;
  vtxD /= 2;
  vtxZ /= 2;
}

void AliTPCCosmicUtils::SubCombined(AliExternalTrackParam *trackPar, const AliTPCseed *seeds[], const Int_t tk0, const Int_t tk1, const Int_t rowstartshift, const Int_t rowstep, const Double_t xmin, const Double_t xmax, const Double_t eloss, Int_t &nfit, Int_t &nmiss, Double_t &pchi2, Double_t &lfit, Double_t &vtxD, Double_t &vtxZ, TTreeSRedirector *debugstreamer)
{
  //
  //sub-routine for combined propagation
  //

  IniCov(trackPar, seeds[0]->GetNumberOfClusters()+seeds[1]->GetNumberOfClusters());

  Int_t dtk=1;
  if(tk0>tk1)
    dtk=-1;

  //reset counters
  Int_t ksite = -999;
  nfit = 0;
  nmiss = 0;
  pchi2 = 0;
  vtxD = -1e10;
  vtxZ = -1e10;

  //always nrow -> 1 -> nrow
  Int_t rowstart= fgkNRow-1-rowstartshift;

  //so that when reversed, the same rows are read! important!!!
  Int_t rowstop = rowstart - rowstart/rowstep * rowstep;
  Int_t drow = -rowstep;
    
  TVector3 gposStart(-999,-999,-999);
  TVector3 gposStop(-999,-999,-999);

  for(Int_t itrk=tk0; dtk*itrk<=dtk*tk1; itrk+=dtk){
    if(itrk==tk1){
      Int_t tmprow = rowstart;
      rowstart = rowstop;
      rowstop = tmprow;
      drow *= -1;
    }
   
    FitKernel(trackPar, seeds[itrk], rowstart, rowstop, drow, xmin, xmax, eloss, ksite, nfit, nmiss, pchi2, gposStart, gposStop, debugstreamer, kFALSE);

    //get the impact parameters at the end of the first propagation (X=0)
    if(itrk==tk0){
      AliExternalTrackParam vertex(*trackPar);
      const Double_t maxStep = 1;
      const Bool_t rotateTo = kFALSE;
      const Double_t maxSnp = 0.8;
      if(AliTrackerBase::PropagateTrackToBxByBz(&vertex, 0, fgkMass, maxStep, rotateTo, maxSnp, eloss)){
        vtxD = TMath::Abs(vertex.GetParameter()[0]);
        vtxZ = vertex.GetParameter()[1];
      }
    }
  }

  TVector3 dpos = gposStart-gposStop;
  lfit = dpos.Pt();
}

void AliTPCCosmicUtils::FitKernel(AliExternalTrackParam *trackPar, const AliTPCseed *tseed, const Int_t rowstart, const Int_t rowstop, const Int_t drow, const Double_t xmin, const Double_t xmax, const Double_t eloss, Int_t &ksite, Int_t &nfit, Int_t &nmiss, Double_t &pchi2, TVector3 &gposStart, TVector3 &gposStop, TTreeSRedirector *debugstreamer, const Bool_t kinicov)
{
  //
  //routine for propagation
  //

  //only for SingleFit-->
  if(kinicov)
    IniCov(trackPar, tseed->GetNumberOfClusters());
  //<--

  Int_t checkstop = -999;

  for(Int_t irow=rowstart; drow*irow<=drow*rowstop; irow+=drow){
    checkstop = irow;
    //printf("test irow %d\n", irow);

    AliTPCclusterMI *cl=tseed->GetClusterPointer(irow);
    if (!cl) continue;
    if (cl->GetX()< fgkXMin) continue;

    //cut on cluster X (i.e. r) to specify leverarm
    if(cl->GetX()< xmin){
      continue;
    }
    if(cl->GetX()> xmax){
      continue;
    }
    //if propagation not successful, trackPar is not changed
    AliExternalTrackParam tmppar(*trackPar);

    //--------------------------- to cure large chi2 in z, 2011-05-11. Thandks to M. Ivanov!  ---------------------------
    const Int_t tmpsite = ((cl->GetDetector()%36)<18);
    if(tmpsite!=ksite && ksite!=-999){
      Double_t *clscov=(Double_t*)tmppar.GetCovariance();
      clscov[2]+=3*3;
    }
    ksite = tmpsite;

    //--------------------------- rotate cluster position to trackPar position ---------------------------
    //DO NOT rotate trackPar, there is "bug" in trackparam::rotate!! stay in the initial one defined as alpha = ATan2(posy0,posx0)
    Float_t gxyz[3];
    cl->GetGlobalXYZ(gxyz);

    const TVector3 gptmp(gxyz[0], gxyz[1], gxyz[2]);

    TVector3 ptogo(gxyz);

    //printf("test %d : %f %f %f\n", irow, gxyz[0], gxyz[1], gxyz[2]);

    const Double_t trkalpha = tmppar.GetAlpha();
    ptogo.RotateZ(-trkalpha);
    const Double_t xTogo = ptogo.X();
    const Double_t yTogo = ptogo.Y();
    const Double_t zTogo = ptogo.Z();
      
    //--------------------------- propagate ---------------------------
    //AliTrackerBase::PropagateTrackToBxByBz Double_t step = dir*TMath::Min(TMath::Abs(xToGo-xpos), maxStep); 
    const Double_t maxStep = 1;
    const Bool_t rotateTo = kFALSE;
    const Double_t maxSnp = 0.8;

    //eloss critical only for combine fit!!!
    const Bool_t kpro = AliTrackerBase::PropagateTrackToBxByBz(&tmppar, xTogo, fgkMass, maxStep, rotateTo, maxSnp, eloss);
    if(!kpro){
      nmiss++;
      continue;
    }

    //--------------------------- set error ---------------------------
    //necessary, help at mostly low p, and also high p
    Double_t cov[3]={0.01, 0., 0.01}; 
    AliTPCseed::GetError(cl, &tmppar, cov[0],cov[2]);
    cov[0]*=cov[0];
    cov[2]*=cov[2];

    if(debugstreamer){
      TVector3 tmpgxyz(gxyz);
      (*debugstreamer)<<"ErrParam"<<
        "Cl.="<<cl<<
        "T.="<<&tmppar<<
        "covy="<<cov[0]<<
        "covz="<<cov[2]<<
        "clpos.="<<&ptogo<<
        "gpos.="<<&tmpgxyz<<
        "\n";
    }

    //--------------------------- get chi2 and THEN update ---------------------------
    Double_t yz[2]={yTogo, zTogo};

    pchi2 += tmppar.GetPredictedChi2(yz, cov);
      
    tmppar.Update(yz,cov);  
      
    //--------------------------- save change ---------------------------
    (*trackPar) = tmppar;

    nfit++;

    gposStop = gptmp;
    if(gposStart.X()<-998)
      gposStart = gptmp;
  }

  //to make sure rowstart and rowstop are the actual ending
  if(checkstop != rowstop){
    printf("exitreport AliTPCCosmicUtils::FitKernel wrong rowstart, stop, drow!! %d %d %d\n", rowstart, rowstop, drow);
    exit(1);
  }
}
