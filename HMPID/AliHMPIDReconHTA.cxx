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

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliHMPIDReconHTA                                                         //
//                                                                      //
// HMPID class to perfom pattern recognition based on Hough transfrom    //
// for single chamber                                                   //
//////////////////////////////////////////////////////////////////////////

#include "AliHMPIDReconHTA.h"//class header
#include "AliHMPIDCluster.h" //CkovHiddenTrk()
#include "AliHMPIDRecon.h"   //FunMinPhot()
#include <TFile.h>           //Database()
#include <TMinuit.h>         //FitFree()
#include <TClonesArray.h>    //CkovHiddenTrk()
#include <AliESDtrack.h>     //CkovHiddenTrk()
#include <TH2F.h>            //InitDatabase()
#include <TGraph.h>          //ShapeModel()
#include <TSpline.h>         //ShapeModel()
#include "TStopwatch.h"      //

Int_t AliHMPIDReconHTA::fgDB[500][150]={75000*0};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDReconHTA::AliHMPIDReconHTA():
  TTask("RichRec","RichPat"),
  fMipX(-999),
  fMipY(-999),
  fMipQ(-999),
  fRadX(-999),
  fRadY(-999),
  fIdxMip(0),
  fNClu(0),
  fXClu(0),
  fYClu(0),
  fClCk(0),
  fThTrkIn(-999),
  fPhTrkIn(-999),
  fThTrkFit(-999),
  fPhTrkFit(-999),
  fCkovFit(-999),
  fCkovSig2(0),
  fParam(AliHMPIDParam::Instance())
{
//..
//hidden algorithm
//..
  fParam->SetRefIdx(fParam->MeanIdxRad()); // initialization of ref index to a default one
  InitDatabase();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
AliHMPIDReconHTA::~AliHMPIDReconHTA()
{
  DeleteVars();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconHTA::InitVars(Int_t n)
{
//..
//Init some variables
//..
  fXClu = new Double_t[n];
  fYClu = new Double_t[n];
  fClCk = new Bool_t[n];
  for(Int_t i=0;i<n;i++) fClCk[i] = kTRUE;
//
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconHTA::DeleteVars()const
{
//..
//Delete variables
//..
  if(fXClu) delete fXClu;
  if(fYClu) delete fYClu;
  if(fClCk) delete fClCk;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDReconHTA::CkovHiddenTrk(AliESDtrack *pTrk,TClonesArray *pCluLst,Int_t index, Double_t nmean)
{
// Pattern recognition method without any infos from tracking:HTA (Hidden Track Algorithm)...
// The method finds in the chmuber the cluster with the highest charge
// compatibile with a MIP, then the strategy is applied
// Arguments:  pTrk     - pointer to ESD track
//             pCluLs   - list of clusters for a given chamber 
//             pNmean   - pointer to ref. index
//             pQthre   - pointer to qthre
//   Returns:           - 0=ok,1=not fitted 
  
  AliHMPIDParam *pParam = AliHMPIDParam::Instance(); 

  if(!CluPreFilter(pCluLst)) return kFALSE;

  Int_t nCh=0;
  Int_t sizeClu=0;
  
  fNClu = pCluLst->GetEntriesFast();
    
  for (Int_t iClu=0;iClu<fNClu;iClu++){                                                       //clusters loop
    AliHMPIDCluster *pClu=(AliHMPIDCluster*)pCluLst->UncheckedAt(iClu);                       //get pointer to current cluster    
    fXClu[iClu] = pClu->X();fYClu[iClu] = pClu->Y();                                          //store x,y for fitting procedure
    fClCk[iClu] = kTRUE;                                                                      //all cluster are accepted at this stage to be reconstructed
    
    if(iClu == index) {

      fMipX = pClu->X();
      fMipY = pClu->Y();
      fMipQ = pClu->Q();
      sizeClu = pClu->Size();
      nCh = pClu->Ch();
      fClCk[index] = kFALSE;
      fIdxMip = index;
 //    Printf(" n. %d x %f y %f Q %f",iClu,pClu->X(),pClu->Y(),pClu->Q());
    }
  }//clusters loop
  
  pParam->SetRefIdx(nmean);
  
  //
  Float_t xra,yra,th,ph; pTrk->GetHMPIDtrk(xra,yra,th,ph);
//  Printf(" simulated phi %6.2f ",ph*TMath::RadToDeg());
  //
  
  if(!DoRecHiddenTrk()) {
    pTrk->SetHMPIDsignal(pParam->kNoPhotAccept);
    return kFALSE;
  }                                                                           //Do track and ring reconstruction,if problems returns 1
//  Printf("    fitted phi %6.2f ",fPhTrkFit*TMath::RadToDeg());
  
  pTrk->SetHMPIDtrk(fRadX,fRadY,fThTrkFit,fPhTrkFit);                                        //store track intersection info
  pTrk->SetHMPIDmip(fMipX,fMipY,(Int_t)fMipQ,fNClu);                                         //store mip info 
  pTrk->SetHMPIDcluIdx(nCh,fIdxMip+1000*sizeClu);                                            //set cham number, index of cluster + cluster size
  pTrk->SetHMPIDsignal(fCkovFit);                                                            //find best Theta ckov for ring i.e. track
  pTrk->SetHMPIDchi2(fCkovSig2);                                                             //errors squared
//    Printf(" n clusters tot %i accepted %i",pCluLst->GetEntriesFast(),fNClu);
//    Printf("CkovHiddenTrk: thetaC %f th %f ph %f",fCkovFit,fThTrkFit,fPhTrkFit);
  
  return kTRUE;
  
}//CkovHiddenTrk()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDReconHTA::DoRecHiddenTrk()
{
// Pattern recognition method without any infos from tracking...
// First a preclustering filter to avoid part of the noise
// Then only ellipsed-rings are fitted (no possibility, 
// for the moment, to reconstruct very inclined tracks)
// Finally a fitting with (th,ph) free, starting by very close values
// previously evaluated.
// Arguments:   none
//   Returns:   none
  Double_t thTrkRec,phiTrkRec,thetaCRec;
  
  if(!FindShape(thTrkRec,phiTrkRec,thetaCRec)) {
//    Printf("FindShape failed...!");
    return kFALSE;
  }

  if(!FitFree(thTrkRec,phiTrkRec)) {
//    Printf("FitFree failed...!");
    return kFALSE;
  }
   
  return kTRUE;
}//DoRecHiddenTrk()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDReconHTA::CluPreFilter(TClonesArray *pCluLst)
{
// Pre-filter of bkg clusters
// Arguments:    pSluLst  -  List of the clusters for a given chamber
//   Returns:    status   -  TRUE if filtering leaves enough photons, FALSE if not
//
  Int_t nClusTot = pCluLst->GetEntriesFast();
  if(nClusTot<4||nClusTot>100) {
    return kFALSE; 
  } else { 
    InitVars(nClusTot);
    return kTRUE;
  }
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDReconHTA::FindShape(Double_t &thTrkRec,Double_t &phiTrkRec,Double_t &thetaCRec)
{
// Finds the estimates for phi and theta of the track and the ThetaCerenkov
// by using a database of the shapes of the rings
// Arguments:   none  
//   Returns:   thTrkRec  - estimate of theta track
//              phiTrkRec - estimate of phi   track
//              thetaCRec - estimate of ThetaCerenkov
//              status    - TRUE if a good solution is found, FALSE if not

  Double_t *phiphot = new Double_t[fNClu];  
  Double_t *dist    = new Double_t[fNClu];  
  Int_t    *indphi  = new    Int_t[fNClu];  

  Bool_t status;
    
// Sort in phi angle...
//  Printf(" mipX %f mipy %f",fMipX,fMipY);
  for(Int_t i=0;i<fNClu;i++) {
    if(!fClCk[i]) {
      phiphot[i] = 999.;
      dist[i] = 999.;
      continue;
    }
    phiphot[i] = (TMath::ATan2(fMipY-fYClu[i],fMipX-fXClu[i])+TMath::Pi())*TMath::RadToDeg();
    dist[i]=TMath::Sqrt((fMipX-fXClu[i])*(fMipX-fXClu[i])+(fMipY-fYClu[i])*(fMipY-fYClu[i]));
//    Printf(" n.%3i  phiphot %f dist %f check %i",i,phiphot[i],dist[i],fClCk[i]);
  }
  
  TMath::Sort(fNClu,phiphot,indphi,kFALSE);
  
// Purify with a truncated mean;
  Int_t np=0;
  Double_t dMean  = 0;
  Double_t dMean2 = 0;
  for(Int_t i=0;i<fNClu;i++) {
    if(!fClCk[indphi[i]]) continue;                                                  // Check if a good photon candidate or not
    dMean +=dist[indphi[i]];
    dMean2+=dist[indphi[i]]*dist[indphi[i]];
    np++;
  }
  
  dMean  /=(Double_t)np;
  dMean2 /=(Double_t)np;
  Double_t rms = TMath::Sqrt(dMean2 - dMean*dMean);
  
  for(Int_t i=0;i<fNClu;i++) {
    if(!fClCk[indphi[i]]) continue;                                                  // Check if a good photon candidate or not
    if(TMath::Abs(dMean-dist[indphi[i]]) > 2*rms) {
      fClCk[indphi[i]] = kFALSE;
      continue;
    }
  }

//
//  purify vectors for good photon candidates
//
  Int_t npeff=0;
  Double_t *phiphotP = new Double_t[fNClu+1];  
  Double_t *distP    = new Double_t[fNClu+1];  
  for(Int_t i=0;i<fNClu;i++) {
    if(!fClCk[indphi[i]]) continue;                                                  // Check if a good photon candidate or not
    phiphotP[npeff] = phiphot[indphi[i]];
    distP[npeff]    = dist[indphi[i]];
//    Printf("n. %2i phi %f dist %f",npeff,phiphotP[npeff],distP[npeff]);
    npeff++;
  }
  
  delete [] phiphot;
  delete [] dist;
  delete [] indphi;

  if(npeff<3) {
    Printf("FindShape failed: no enough photons = %i...",npeff);
    delete [] phiphotP;
    delete [] distP;
    return kFALSE;
  }

//  for(Int_t i=0;i<npeff;i++) {Printf(" n. %d phiphot %f dist %f",i,phiphotP[i],distP[i]);}
  
  Double_t xA,xB;
  status = kFALSE;
  
  if(!ShapeModel(npeff,phiphotP,distP,xA,xB,phiTrkRec)) {/*Printf("ShapeModel failed            ");*/ return kFALSE;}
    
//  if(xA > 50 || xB > 15)                                {Printf("xA and xB failed out of range"); return kFALSE;}

  Int_t binxDB,binyDB;
  FindBinDB(xA,xB,binxDB,binyDB);
  if(binxDB<0 || binyDB<0)                              {/*Printf("bin < 0 ! failed             ");*/ return kFALSE;}
  
  Int_t compactDB = fgDB[binxDB][binyDB];
  
  if(compactDB<0)                                       {/*Printf("compact< 0! failed           ");*/ return kFALSE;} 
  //
  //
  thetaCRec =  (Double_t)(compactDB%1000);
  thTrkRec  =  (Double_t)(compactDB/1000);

  thTrkRec  *= TMath::DegToRad(); 
  thetaCRec *= TMath::DegToRad();

  status = kTRUE;

  delete [] phiphotP;
  delete [] distP;
  
  return status;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDReconHTA::ShapeModel(Int_t np,Double_t *phiphot,Double_t *dist,Double_t &xA,Double_t &xB,Double_t &phiStart)
{
// Find a Spline curve to define dist. vs. phi angle
// in order to estimate the phi of the track
// Arguments:   np     - # points corresponding to # photon candidates
//             dist    - distance of each photon from MIP
//             phiphot - phi of the photon in the DRS
//   Returns:  xA      - min. distance from MIP
//             xB      - dist. from mip perpedicular to the major axis 
//             phiStart- estimate of the track phi

  TGraph *phigr = new TGraph(np,phiphot,dist);
  TSpline3 *sphi = new TSpline3("sphi",phigr);
  if(!sphi) {Printf("Spline not created!Bye.");return kFALSE;}

  Int_t locMin = TMath::LocMin(np,dist);
  Int_t locMax = TMath::LocMax(np,dist);
  
  Double_t minX = phiphot[locMin];
//  Double_t minY =    dist[locMin];
  Double_t maxX = phiphot[locMax];
//  Double_t maxY =    dist[locMax];
  
  Int_t ip[3] = {-1,0,1};
  if(locMin==0   ) {ip[0]= 0;ip[1]= 1;ip[2]= 2;}
  if(locMin==np-1) {ip[0]=-2;ip[1]=-1;ip[2]= 0;}
  
  Double_t minXf = VertParab(phiphot[locMin+ip[0]],dist[locMin+ip[0]], 
                             phiphot[locMin+ip[1]],dist[locMin+ip[1]],
                             phiphot[locMin+ip[2]],dist[locMin+ip[2]]);
  if(minXf< phiphot[locMin+ip[0]] || minXf > phiphot[locMin+ip[2]]) minXf = minX;
  
  ip[0]=-1;ip[1]=0;ip[2]=1;
  if(locMax==0   ) {ip[0]= 0;ip[1]= 1;ip[2]= 2;}
  if(locMax==np-1) {ip[0]=-2;ip[1]=-1;ip[2]= 0;}
  
  Double_t maxXf = VertParab(phiphot[locMax+ip[0]],dist[locMax+ip[0]], 
                             phiphot[locMax+ip[1]],dist[locMax+ip[1]],
                             phiphot[locMax+ip[2]],dist[locMax+ip[2]]);
  if(maxXf< phiphot[locMax+ip[0]] || maxXf > phiphot[locMax+ip[2]]) maxXf = maxX;
  
//  
  if(TMath::Abs(maxXf-minXf)>30) {
    xA = sphi->Eval(minXf);
    if(xA < 0) xA = dist[sphi->FindX(xA)];
    xB = sphi->Eval(minXf-90);
    if(xB < 0) xB = dist[sphi->FindX(xB)];
    phiStart = minXf-180;  //open ring or acceptance effect...so believe to min phi angle!
  } else {
    phiStart = 0.5*(maxXf-180+minXf);
    xA = sphi->Eval(phiStart+180);
    if(xA < 0) xA = dist[sphi->FindX(xA)];
    xB = sphi->Eval(phiStart+90);
    if(xB < 0) xB = dist[sphi->FindX(xB)];
  }
  //
//  Printf("ShapeModel: phiStart %f xA %f xB %f",phiStart,xA,xB);
  
  
  phiStart*=TMath::DegToRad();
  //----
  Double_t phitest = FindSimmPhi();   
  phiStart = phitest;
  //---
//  Printf("phiStart %f phiTest %f",phiStart*TMath::RadToDeg(),phitest*TMath::RadToDeg());
 
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t AliHMPIDReconHTA::VertParab(Double_t x1,Double_t y1,Double_t x2, Double_t y2, Double_t x3, Double_t y3)const
{
// It uses parabola from 3 points to evaluate the x-coord of the parab 
// Arguments:    xi,yi - points
//   Returns:    x-coord of the vertex 
  
  Double_t a = ((x1-x3)*(y1-y2)-(x1-x2)*(y1-y3))/((x1*x1-x2*x2)*(x1-x3)-(x1*x1-x3*x3)*(x1-x2));
  Double_t b = (y1-y2 - a*(x1*x1-x2*x2))/(x1-x2);
//  Double_t c = y1 - a*x1*x1-b*x1;
  return -0.5*b/a;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t AliHMPIDReconHTA::FitFree(Double_t thTrkRec,Double_t phiTrkRec)
{
// Fit performed by minimizing RMS/sqrt(n) of the
// photons reconstructed. First phi is fixed and theta
// is fouond, then (th,ph) of the track
// as free parameters
// Arguments:    PhiRec phi of the track
//   Returns:    none
  
  TMinuit *pMinuit = new TMinuit(2);
  pMinuit->mncler();
  gMinuit->SetObjectFit((TObject*)this);  gMinuit->SetFCN(AliHMPIDReconHTA::FunMinPhot);  //set fit function
  Double_t aArg=-1,parStep,parLow,parHigh;     Int_t iErrFlg;                 //tmp vars for TMinuit
  Double_t d1,d2,d3;
  TString sName;
  Double_t th,ph;
  
  pMinuit->mnexcm("SET PRI",&aArg,1,iErrFlg);                                          //suspend all printout from TMinuit
  pMinuit->mnexcm("SET NOW",&aArg,0,iErrFlg);

  if(thTrkRec==0) thTrkRec = 3.*TMath::DegToRad();    // not to start from the edge...
  
  pMinuit->mnparm(0," thTrk  ",thTrkRec ,parStep=0.01,parLow=0,parHigh=TMath::PiOver4(),iErrFlg);
  pMinuit->mnparm(1," phiTrk ",phiTrkRec,parStep=0.01,parLow=0,parHigh=TMath::TwoPi(),iErrFlg);
  
  pMinuit->FixParameter(1);
  pMinuit->mnexcm("SIMPLEX" ,&aArg,0,iErrFlg);   
  pMinuit->mnexcm("MIGRAD"  ,&aArg,0,iErrFlg);
  pMinuit->Release(1);  
  pMinuit->mnexcm("MIGRAD"  ,&aArg,0,iErrFlg);
  
  pMinuit->mnpout(0,sName,th,d1,d2,d3,iErrFlg);
  pMinuit->mnpout(1,sName,ph,d1,d2,d3,iErrFlg);   
  
  Double_t f,par[2];
  Double_t *grad=0x0;
  par[0] = th;par[1] = ph;
  pMinuit->Eval(2,grad,f,par,3);

  SetTrkFit(th,ph);
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconHTA::FunMinPhot(Int_t &/* */,Double_t* /* */,Double_t &f,Double_t *par,Int_t iflag)
{
// Minimization function to find best track and thetaC parameters
// Arguments:    f = function value to minimize
//             par = list of parameter to find
//           iflag = flag status. See Minuit instructions
//   Returns:    none
//
// Note: it is necessary to call an instance of AlihMPIDParam. Not possible to use fParam
// because of the static instantiation of the function in Minuit
  
  AliHMPIDParam *pParam=AliHMPIDParam::Instance();
  AliHMPIDReconHTA *pRecHTA=(AliHMPIDReconHTA*)gMinuit->GetObjectFit();
  AliHMPIDRecon pRec;
  Double_t sizeCh = 0.5*pParam->RadThick()+pParam->WinThick()+pParam->GapThick();
  Double_t thTrk = par[0]; 
  Double_t phTrk = par[1];
  Double_t xrad = pRecHTA->MipX() - sizeCh*TMath::Tan(thTrk)*TMath::Cos(phTrk);
  Double_t yrad = pRecHTA->MipY() - sizeCh*TMath::Tan(thTrk)*TMath::Sin(phTrk);
  pRecHTA->SetRadXY(xrad,yrad);
  pRec.SetTrack(xrad,yrad,thTrk,phTrk);

  Double_t meanCkov =0;
  Double_t meanCkov2=0;
  Double_t thetaCer,phiCer;
  Int_t nClAcc = 0;
  Int_t nClTot=pRecHTA->NClu();
    
  for(Int_t i=0;i<nClTot;i++) {
    if(!(pRecHTA->ClCk(i))) continue;
    pRec.FindPhotCkov(pRecHTA->XClu(i),pRecHTA->YClu(i),thetaCer,phiCer);  
    meanCkov  += thetaCer;
    meanCkov2 += thetaCer*thetaCer;
    nClAcc++;
  }
  if(nClAcc==0) {f=999;return;}
  meanCkov  /=(Double_t)nClAcc;
  meanCkov2 /=(Double_t)nClAcc;
  Double_t rms = TMath::Sqrt(TMath::Abs(meanCkov2 - meanCkov*meanCkov));
  f = rms/TMath::Sqrt((Double_t)nClAcc);
  
  if(iflag==3) {
/*
    Printf("FunMinPhot before: photons candidates %i used %i",nClTot,nClAcc);
    nClAcc = 0;
    Double_t meanCkov1=0;
    Double_t meanCkov2=0;
    for(Int_t i=0;i<nClTot;i++) {
      if(!(pRec->ClCk(i))) continue;
      pRec->FindPhotCkov(pRec->XClu(i),pRec->YClu(i),thetaCer,phiCer);  
      if(TMath::Abs(thetaCer-meanCkov)<2*rms) {
        meanCkov1 += thetaCer;
        meanCkov2 += thetaCer*thetaCer;
        nClAcc++;
      } else pRec->SetClCk(i,kFALSE);
    }
    meanCkov1/=nClAcc;
    Double_t rms2 = (meanCkov2 - meanCkov*meanCkov*nClAcc)/nClAcc;
    Printf("FunMinPhot after: photons candidates %i used %i thetaCer %f",nClTot,nClAcc,meanCkov1);
    pRec->SetCkovFit(meanCkov1);
    pRec->SetCkovSig2(rms2);
    pRec->SetNClu(nClAcc);
*/
//    Printf("FunMinPhot: reconstructed theta Cerenkov %f with %d photons",meanCkov,nClAcc);
    pRecHTA->SetCkovFit(meanCkov);
    pRecHTA->SetCkovSig2(rms*rms);
    pRecHTA->SetNClu(nClAcc);
  }
}//FunMinPhot()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconHTA::InitDatabase()
{
// Construction a database of ring shapes on fly
//   Arguments: none
//   Returns  : none
//  N.B. fgDB is the distance with x-min from MIP
//                                 y-dist from the ring of the MIP perpendicular to major axis
//        The content is the packed info of track theta and thetaC in degrees
//                        thetaC+1000*thTrk
//
//  TFile *pout = new TFile("./database.root","recreate");

  static Bool_t isDone = kFALSE;
  
  TStopwatch timer;
  
  if(isDone) {
   return;
  }
  
  if(!isDone) {
    timer.Start();
  }
 
  AliInfo(Form("database HTA is being built.Please, wait..."));
//  
  Double_t x[3]={0,0,0},y[3];

  AliHMPIDRecon rec;

  if(!fParam) fParam=AliHMPIDParam::Instance();
  Double_t thetaMax = TMath::ACos(1./fParam->MeanIdxRad());
  Double_t thTrkMax = 1./TMath::ASin(fParam->MeanIdxRad());    

  Int_t nstepx = 1000;
  Int_t nstepy = 1000;

  //
  Double_t xrad = 0;
  Double_t yrad = 0;
  Double_t phTrk = 0;

  for(Int_t i=0;i<nstepx;i++) {     //loop on thetaC
    for(Int_t j=0;j<nstepy;j++) {   //loop on theta particle
      Double_t thetaC = thetaMax/nstepx*((Double_t)i+0.5);
      Double_t thTrk  = thTrkMax/nstepy*((Double_t)j+0.5);
      //
      //mip position
      //
      Double_t sizeCh = 0.5*fParam->RadThick()+fParam->WinThick()+fParam->GapThick();
      Double_t xmip = xrad + sizeCh*TMath::Tan(thTrk)*TMath::Cos(phTrk);
      Double_t ymip = yrad + sizeCh*TMath::Tan(thTrk)*TMath::Sin(phTrk);

      Double_t dist1,dist2;
      //
      //first point at phi=0
      //
      rec.SetTrack(xrad,yrad,thTrk,phTrk);
      TVector2 pos;
      pos=rec.TracePhot(thetaC,0);

      if(pos.X()==-999) {
        dist1 = 0;            //open ring...only the distance btw mip and point at 180 will be considered
      } else {
        x[0] = pos.X(); y[0] = pos.Y();
        dist1   = TMath::Sqrt((x[0]-xmip)*(x[0]-xmip)+(y[0]-ymip)*(y[0]-ymip));
      }
      //
      //second point at phi=180
      //
      rec.SetTrack(xrad,yrad,thTrk,phTrk);
      pos=rec.TracePhot(thetaC,TMath::Pi());

      if(pos.X()==-999) {Printf("it should not happens!Bye");return;}
      x[1] = pos.X(); y[1] = pos.Y();
      if((x[1]-xmip)*(x[0]-xmip)>0) continue; // to avoid circles out mips (for very low ThetaC)
      dist2   = TMath::Sqrt((x[1]-xmip)*(x[1]-xmip)+(y[1]-ymip)*(y[1]-ymip));

//      Double_t distA = dist1+dist2;
      Double_t distA = dist2;     // only the minimum: problem of acceptance
      //
      //second point at phi=90
      //
      rec.SetTrack(xrad,yrad,thTrk,phTrk);
      pos=rec.TracePhot(thetaC,TMath::PiOver2());

      if(pos.X()==-999) continue;
      x[2] = pos.X(); y[2] = pos.Y();
      Double_t distB   = TMath::Sqrt((x[2]-xmip)*(x[2]-xmip)+(y[2]-ymip)*(y[2]-ymip));
// compact the infos...      
      Int_t compact = (Int_t)(thetaC*TMath::RadToDeg())+1000*(Int_t)(thTrk*TMath::RadToDeg());
      Int_t binxDB,binyDB;
      FindBinDB(distA,distB,binxDB,binyDB);
      if(fgDB[binxDB][binyDB]==0) fgDB[binxDB][binyDB] = compact;
    }
  }

  FillZeroChan();

  if(!isDone) {
    timer.Stop();
    Double_t nSecs = timer.CpuTime();  
    AliInfo(Form("database HTA successfully open in %3.1f sec.(CPU). Reconstruction is started.",nSecs));
    isDone = kTRUE;
  }
  
//  pout->Write();
//  pout->Close();
  
}//InitDatabase()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconHTA::FillZeroChan()const
{
  //If fills eventually channel without entries
  //inthe histo "database" jyst interpolating the neighboring cells
  // Arguments: histogram pointer of the database
  //   Returns: none
  //
  const Int_t nxDB = 500;
  const Int_t nyDB = 150;

  for(Int_t i = 0;i<nxDB;i++) {
    for(Int_t j = 0;j<nyDB;j++) {
      if(fgDB[i][j] == 0) {
        fgDB[i][j] = -1;
        Int_t nXmin = i-1; Int_t nXmax=i+1;
        Int_t nYmin = j-1; Int_t nYmax=j+1;
        Int_t nc = 0;
        Double_t meanC =0;
        Double_t meanTrk =0;
        for(Int_t ix=nXmin;ix<=nXmax;ix++) {
          if(ix<0||ix>=nxDB) continue;
          for(Int_t iy=nYmin;iy<=nYmax;iy++) {
            if(iy<0||iy>=nyDB) continue;
            meanC  +=  (Int_t)(fgDB[ix][iy]%1000);
            meanTrk+=  (Int_t)(fgDB[ix][iy]/1000);
            nc++;
          }
          meanC/=nc; meanTrk/=nc;
          Int_t compact = (Int_t)meanC+1000*(Int_t)meanTrk;
          if(compact>0)fgDB[i][j] = compact;
        }
      }
    }
  }
}//FillZeroChan()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Int_t AliHMPIDReconHTA::r2(Double_t *coef, Double_t &x1, Double_t &x2)
{
  //2nd deg. equation
  //solution
  // Arguments: coef 2 1 0: ax^2+bx+c=0
  //   Returns: n. of solutions
  //            x1= 1st sol
  //            x2= 2nd sol
  Double_t a,b,c;
  a = coef[2];
  b = coef[1];
  c = coef[0];
  Double_t delta = b*b-4*a*c;
  if(delta<0) {return 0;}
  if(delta==0) {
    x1=x2=-b/(2*a);
    return 1;
  }
  // delta>0
  x1 = (-b+TMath::Sqrt(delta))/(2*a);
  x2 = (-b-TMath::Sqrt(delta))/(2*a);
  return 2;
}//r2()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t AliHMPIDReconHTA::FindSimmPhi() 
{     
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// RESTITUISCE  IN  OUTPUT  IL  VALORE  FINALE  DELL'ANGOLO  RICOSTRUITO 
// CON  I  DUE  METODI:
// - metodo dei minimi quadrati con le distanze effettive;................................(PER RING CHIUSI)
// - metodo della determin della pendenza individuando la distanza minima mip-fotone;.....(PER RING APERTI)
  
  Float_t coeff1ord=0;     Float_t coeff2ord=0;     Float_t coeff0ord=0;  
  Float_t xrotsumm =0;     Float_t yrotsumm =0;     Float_t xx =0;           
  Float_t yy =0;           Float_t xy =0;           Float_t yx =0; 
  Double_t xmin=0;
  Double_t ymin=0;

  Int_t np=0;    
    
  Double_t distMin = 999.;
  
  for(Int_t i=0;i<fNClu;i++) {
    if(!fClCk[i]) continue;
    np++;
    xrotsumm+=fXClu[i];         // summ xi
    yrotsumm+=fYClu[i];         // summ yi
    xx+=fXClu[i]*fXClu[i];      // summ xixi     
    yy+=fYClu[i]*fYClu[i];      // summ yiyi
    xy+=fXClu[i]*fYClu[i];      // summ yixi     
    Double_t dist2= (fXClu[i]-fMipX)*(fXClu[i]-fMipX)+(fYClu[i]-fMipY)*(fYClu[i]-fMipY);
    if(dist2<distMin) {
      distMin = dist2;
      xmin = fXClu[i];
      ymin = fYClu[i];
    }
  }

  Double_t AngM=0;
  if(ymin <  fMipY && xmin >  fMipX)  {AngM  =  TMath::ATan((ymin-fMipY)/(xmin-fMipX))*TMath::RadToDeg()+180;}
  if(ymin >  fMipY && xmin <  fMipX)  {AngM  =  TMath::ATan((ymin-fMipY)/(xmin-fMipX))*TMath::RadToDeg()+360;}
  if(ymin >  fMipY && xmin >  fMipX)  {AngM  =  TMath::ATan((ymin-fMipY)/(xmin-fMipX))*TMath::RadToDeg()+180;}
  if(ymin <  fMipY && xmin <  fMipX)  {AngM  =  TMath::ATan((ymin-fMipY)/(xmin-fMipX))*TMath::RadToDeg();}
  if(ymin == fMipY && xmin >  fMipX)  {AngM  =  TMath::ATan((ymin-fMipY)/(xmin-fMipX))*TMath::RadToDeg()+180;}
  if(ymin == fMipY && xmin <  fMipX)  {AngM  =  TMath::ATan((ymin-fMipY)/(xmin-fMipX))*TMath::RadToDeg();}
  if(ymin <  fMipY && xmin == fMipX)  {AngM  =  90;}
  if(ymin >  fMipY && xmin == fMipX)  {AngM  =  270;}
  
  //_____calc. met min quadr using effective distance _________________________________________________
  
  coeff2ord= xy-xrotsumm*yrotsumm/np;    
  coeff1ord= yrotsumm*yrotsumm/np - xrotsumm*xrotsumm/np - yy + xx;
  coeff0ord= xrotsumm*yrotsumm/np -  yx;
  Double_t m1=0, m2=0;  Double_t n1=0, n2=0;
                            // c           // b         // a
  Double_t coeff[3]={coeff0ord,coeff1ord,coeff2ord};    
  
  r2(coeff,m1,m2);
  
  n1=(yrotsumm-m1*xrotsumm)/np;                         
  n2=(yrotsumm-m2*xrotsumm)/np;
  // 2 solutions.................
  
  Double_t PhiTrk1= TMath::ATan(m1);                              
  //Double_t PhiTrk2= TMath::ATan(m2);
  
  // negative angles solved...
  
  Double_t PhiTrk1Positive=0;
  
  if(PhiTrk1<0)    PhiTrk1Positive= PhiTrk1 + 180*TMath::DegToRad();
  if(PhiTrk1>=0)   PhiTrk1Positive= PhiTrk1;
  
  Double_t d1 =(1/(m1*m1+1))*(yy+m1*m1*xx+np*n1*n1-2*m1*xy-2*n1*yrotsumm+2*m1*n1*xrotsumm);
  Double_t d2 =(1/(m2*m2+1))*(yy+m2*m2*xx+np*n2*n2-2*m2*xy-2*n2*yrotsumm+2*m2*n2*xrotsumm);
  
  Double_t mMin;
  if(d1 > d2) mMin = m2; else mMin = m1;
  
  //Double_t PhiTrk = TMath::ATan(mMin)*TMath::RadToDeg();
  Double_t PhiTrkPositive=0;
  // 
  if(ymin <  fMipY && xmin >  fMipX)  {PhiTrkPositive  =  TMath::ATan(mMin)*TMath::RadToDeg()+180;}
  if(ymin >  fMipY && xmin <  fMipX)  {PhiTrkPositive  =  TMath::ATan(mMin)*TMath::RadToDeg()+360;}
  if(ymin >  fMipY && xmin >  fMipX)  {PhiTrkPositive  =  TMath::ATan(mMin)*TMath::RadToDeg()+180;}
  if(ymin <  fMipY && xmin <  fMipX)  {PhiTrkPositive  =  TMath::ATan(mMin)*TMath::RadToDeg();}
  if(ymin == fMipY && xmin >  fMipX)  {PhiTrkPositive  =  TMath::ATan(mMin)*TMath::RadToDeg()+180;}
  if(ymin == fMipY && xmin <  fMipX)  {PhiTrkPositive  =  TMath::ATan(mMin)*TMath::RadToDeg();}
  if(ymin <  fMipY && xmin == fMipX)  {PhiTrkPositive  =   90;}
  if(ymin >  fMipY && xmin == fMipX)  {PhiTrkPositive  =  270;}
  
  //  ------------------------- choose the best-----------------------
  
  
  Double_t PhiTrkFinal=0;
  if( AngM-40 <=  PhiTrkPositive  &&  AngM+40 >= PhiTrkPositive)   PhiTrkFinal = PhiTrkPositive;  else PhiTrkFinal = AngM;
  
  return PhiTrkFinal*TMath::DegToRad();
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void AliHMPIDReconHTA::FindBinDB(Double_t x,Double_t y,Int_t &binX,Int_t &binY)
{
  const Int_t nxDB = 500;
  const Int_t nyDB = 150;
  const Double_t xlowDB =  0;   
  const Double_t xhigDB = 50;   
  const Double_t ylowDB =  0;
  const Double_t yhigDB = 15;

  binX = -1;
  binY = -1;
  if(x<xlowDB && x>xhigDB &&
     y<ylowDB && y>yhigDB)    return;
  binX = Int_t((x-xlowDB)/(xhigDB-xlowDB)*nxDB);
  binY = Int_t((y-ylowDB)/(yhigDB-ylowDB)*nyDB);
  if(binX>=nxDB || binY>=nyDB) {
    binX = -1;
    binY = -1;
  }
  
}
