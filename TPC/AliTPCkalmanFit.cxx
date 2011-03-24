/*
  marian.ivanov@cern.ch 
  
  AliTPCkalmanFit:

  Kalman filter(s) for fitting of the tracks together with calibration/transformation 
  parameters.

  Correction/Transformation are currently described by set (TObjArray) of primitive 
  correction/transformatio - AliTPCTransformation.  Currently we assume that transformation 
  comute (in first order). AliTPCTransformation describe general non linear transformation.   
  
  Current calibration parameters and covariance stored (fCalibParam, fCalibCovar).
  
  Currenly only linear track model implemented.
  Fits to be implemented:
   0. Plane fitting (Laser CE)
   1. Primary vertex fitting.
   2. Propagation in magnetic field and fit of planes


  
  How to use it - see  AliTPCkalmanFit::Test function

  Simple test (see AliTPCkalmanFit::Test)  
  AliTPCTransformation::BuildBasicFormulas();
  AliTPCkalmanFit *kalmanFit0 = AliTPCkalmanFit::Test(2000);
  TFile f("kalmanfitTest.root");
 
  Transformation visualization:
  Transforamtion can be visualized using TFn (TF1,TF2 ...) and using tree->Draw()
  e.g:
  kalmanFit0->SetInstance(kalmanFit0);   // 
  kalmanFit0->InitTransformation();      //
  //
  TF2 fxRZdz("fxRZdz","sign(y)*1000*(AliTPCkalmanFit::SGetTPCDeltaXYZ(0,-1,0,x,0,y)-AliTPCkalmanFit::SGetTPCDeltaXYZ(0,-1,0,x,0,y-1))",85,245,-250,250);
  fxRZdz->Draw("");

  TF2 fxRZ("fxRZ","sign(y)*10*(AliTPCkalmanFit::SGetTPCDeltaXYZ(0,-1,x,0,y))",85,245,-250,250);
  fxRZ->Draw("");



*/

#include "TRandom.h"
#include "TMath.h"
#include "TBits.h"
#include "TFormula.h"
#include "TF1.h"
#include "TLinearFitter.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TAxis.h"


#include "TTreeStream.h"
#include "AliTrackPointArray.h"
#include "AliLog.h"
#include "AliTPCTransformation.h"
#include "AliTPCkalmanFit.h"

ClassImp(AliTPCkalmanFit)

AliTPCkalmanFit* AliTPCkalmanFit::fgInstance = 0;

AliTPCkalmanFit::AliTPCkalmanFit():
  TNamed(),
  fCalibration(0),
  fCalibParam(0),
  fCalibCovar(0),
  fLinearParam(0),
  fLinearCovar(0),
  fLastTimeStamp(-1),
  fCurrentAlpha(0),  //! current rotation frame
  fCA(0),            //! cosine of current angle
  fSA(0)            //! sinus of current angle    
{
  //
  // Default constructor
  //  
  for (Int_t ihis=0; ihis<12; ihis++){
    fLinearTrackDelta[ihis]=0;
    fLinearTrackPull[ihis]=0;
  }
}

void AliTPCkalmanFit::InitTransformation(){
  //
  // Initialize pointers to the transforamtion functions
  //
  //
  Int_t ncalibs = fCalibration->GetEntries();
  for (Int_t icalib=0;icalib<ncalibs; icalib++){
    AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
    transform->Init();
  }
}

void AliTPCkalmanFit::Add(const AliTPCkalmanFit * kalman){
  //
  //
  //
  Update(kalman);
  for (Int_t i=0;i<12;i++){
    if (fLinearTrackDelta[i] && kalman->fLinearTrackDelta[i]){
      fLinearTrackDelta[i]->Add(kalman->fLinearTrackDelta[i]);
    }
    if (fLinearTrackPull[i] && kalman->fLinearTrackPull[i]){
      fLinearTrackPull[i]->Add(kalman->fLinearTrackPull[i]);
    }
  }
  
}


void AliTPCkalmanFit::Init(){
  //
  // Initialize parameter vector and covariance matrix
  // To be called after initialization of all of the transformations
  //
  //
  Int_t ncalibs = fCalibration->GetEntries();
  fCalibParam = new TMatrixD(ncalibs,1);
  fCalibCovar = new TMatrixD(ncalibs,ncalibs);
  for (Int_t icalib=0;icalib<ncalibs; icalib++){
    AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
    (*fCalibParam)(icalib,0) = transform->GetParam();
    for (Int_t jcalib=0;jcalib<ncalibs; jcalib++){
      if (icalib!=jcalib) (*fCalibCovar)(icalib,jcalib)= 0;
      if (icalib==jcalib) (*fCalibCovar)(icalib,jcalib) = transform->GetSigma()*transform->GetSigma();    
    }
  }
  //
  // Build QA histograms
  //
  Double_t mpi = TMath::Pi();
  //
  //                    delta  alpha   y0     z0    ky   kz 
  Int_t binsQA[6]    = {300,   36*4,   30,    25,   20,  20};
  Double_t xminQA[6] = {-0.5,  -mpi,  -120, -250,   -1,  -1.};
  Double_t xmaxQA[6] = { 0.5,   mpi,   120,  250,    1,   1.};  
  TString  axisName[6]={"#Delta",
			"#alpha",
			"y_{0}",
			"z_{0}",
			"tan(#phi)",
			"tan(theta)"};
  //
  TString deltaName[12]={"#Delta_{y}(cm)",
			 "100*#Delta_{#phi}(cm)",
			 "100^{2}dy0^{2}/dx0^{2}(cm)",
			 "100^{2}dy1^{2}/dx1^{2}(cm)",
			 "#Delta_{z}(cm)",
			 "100*#Delta_{#theta}(cm)",
			 "100^{2}*dz0^{2}/dx0^{2}(cm)",
			 "100^{2}*dz1^{2}/dx1^{2}(cm)",
			 "RMSy_{0} (cm)",
			 "RMSy_{1} (cm)",
			 "RMSz_{0} (cm)",
			 "RMSz_{1} (cm)"};

  TString pullName[12]={"#Delta_{y}(unit)",
			"100*#Delta_{#phi}(unit)",
			"100^{2}dy0^{2}/dx0^{2}(unit)",
			"100^{2}dy1^{2}/dx1^{2}(unit)",
			"#Delta_{z}(unit)",
			"100*#Delta_{#theta}(unit)",
			"100^{2}*dz0^{2}/dx0^{2}(unit)",
			"100^{2}*dz1^{2}/dx1^{2}(unit)"
			"RMSy_{0} (cm)",
			"RMSy_{1} (cm)",
			"RMSz_{0} (cm)",
			"RMSz_{1} (cm)"};
  //
  //
  //
  for (Int_t ihis=0; ihis<12; ihis++){
    fLinearTrackDelta[ihis]=0;
    fLinearTrackPull[ihis]=0;
    xminQA[0]=-0.5; xmaxQA[0] = 0.5; 
    fLinearTrackDelta[ihis]  = new THnSparseS(deltaName[ihis],deltaName[ihis], 6, binsQA,xminQA, xmaxQA);
    xminQA[0]=-10; xmaxQA[0]  = 10;     
    fLinearTrackPull[ihis]   = new THnSparseS(pullName[ihis],pullName[ihis],   6, binsQA,xminQA, xmaxQA);
    for (Int_t iaxis=1; iaxis<6; iaxis++){
      fLinearTrackDelta[ihis]->GetAxis(iaxis)->SetName(axisName[iaxis]);
      fLinearTrackDelta[ihis]->GetAxis(iaxis)->SetTitle(axisName[iaxis]);
      fLinearTrackPull[ihis]->GetAxis(iaxis)->SetName(axisName[iaxis]);
      fLinearTrackPull[ihis]->GetAxis(iaxis)->SetTitle(axisName[iaxis]);
    }
    fLinearTrackDelta[ihis]->GetAxis(0)->SetName(deltaName[ihis]);
    fLinearTrackDelta[ihis]->GetAxis(0)->SetTitle(deltaName[ihis]);
    fLinearTrackPull[ihis]->GetAxis(0)->SetName(deltaName[ihis]);
    fLinearTrackPull[ihis]->GetAxis(0)->SetTitle(deltaName[ihis]);
  }

  
}

void AliTPCkalmanFit::SetStatus(const char * mask, Bool_t setOn, Bool_t isOr){
  //
  // 0. To activate all transforamtion call SetStatus(0,kTRUE)
  // 1. To disable everything               SetStatus(0,kFALSE)
  // 2. To activate/desactivate             SetStatus("xxx",kTRUE/kFALSE,kFALSE)
  Int_t ncalibs = fCalibration->GetEntries();
  if (mask==0) {
    for (Int_t i=0; i<ncalibs;i++){
      AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(i);
      transform->SetActive(setOn); 
    }
  }
  else{
    for (Int_t i=0; i<ncalibs;i++){
      AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(i);
      TString  strName(transform->GetName());
      if (strName.Contains(mask)){
	if (!isOr) transform->SetActive( transform->IsActive() && setOn); 
	if (isOr)  transform->SetActive( transform->IsActive() || setOn); 
      }
    }
  }
}


void AliTPCkalmanFit::Update(const AliTPCkalmanFit * kalman){
  //
  // Update Kalman filter
  //
  Int_t ncalibs = fCalibration->GetEntries();
  TMatrixD vecXk=*fCalibParam;       // X vector
  TMatrixD covXk=*fCalibCovar;       // X covariance
  TMatrixD &vecZk = *(kalman->fCalibParam);
  TMatrixD &measR = *(kalman->fCalibCovar);

  TMatrixD matHk(ncalibs,ncalibs);   // vector to mesurement
  TMatrixD vecYk(ncalibs,1);         // Innovation or measurement residual
  TMatrixD matHkT(ncalibs,ncalibs);  // helper matrix Hk transpose
  TMatrixD matSk(ncalibs,ncalibs);   // Innovation (or residual) covariance
  TMatrixD matKk(ncalibs,ncalibs);   // Optimal Kalman gain
  TMatrixD covXk2(ncalibs,ncalibs);  // helper matrix
  TMatrixD covXk3(ncalibs,ncalibs);  // helper matrix
  //
  for (Int_t i=0;i<ncalibs;i++){
    for (Int_t j=0;j<ncalibs;j++) matHk(i,j)=0;
    matHk(i,i)=1;
  }
   vecYk = vecZk-matHk*vecXk;               // Innovation or measurement residual
   matHkT=matHk.T(); matHk.T();
   matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
   matSk.Invert();
   matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
   vecXk += matKk*vecYk;                    //  updated vector
   covXk2= (matHk-(matKk*matHk));
   covXk3 =  covXk2*covXk;
   covXk = covXk3;
   (*fCalibParam) = vecXk;
   (*fCalibCovar) = covXk;
}





void AliTPCkalmanFit::FitTrackLinear(AliTrackPointArray &points, TTreeSRedirector *debug,  Float_t scalingRMSY, Float_t scalingRMSZ){
  //
  //

  if (!fCalibParam) {
    AliError("Kalman Fit not initialized");
    return;
  }
  //
  // 1. Make initial track hypothesis
  //
  TLinearFitter lfitY(2,"pol1");
  TLinearFitter lfitZ(2,"pol1");
  TVectorD vecZ(2);
  TVectorD vecY(2);
  //
  lfitY.StoreData(kTRUE);
  lfitZ.StoreData(kTRUE);
  lfitY.ClearPoints();
  lfitZ.ClearPoints();  
  Int_t npoints = points.GetNPoints();
  if (npoints<2) return;
  const Double_t kFac=npoints*npoints*100;
  const Double_t kOff0=0.01*0.01;
  const Double_t kOff1=kOff0/(250.*250.);
  //
  // 0. Make seed
  //
  //  AliTrackPointArray pointsc(points);
  //ApplyCalibration(&pointsc, -1.);
  //
  // 1.a Choosing reference rotation alpha
  //
  fCurrentAlpha = TMath::ATan2(points.GetY()[npoints-1]-points.GetY()[0], points.GetX()[npoints-1]-points.GetX()[0]);  
  fCA = TMath::Cos(fCurrentAlpha);
  fSA = TMath::Sin(fCurrentAlpha);
  //
  // 1.b Fit the track in the rotated frame - MakeSeed 
  //
  for (Int_t ipoint=0; ipoint<npoints-1; ipoint++){
    Double_t rx =   fCA*points.GetX()[ipoint]+fSA*points.GetY()[ipoint];
    Double_t ry =  -fSA*points.GetX()[ipoint]+fCA*points.GetY()[ipoint];
    Double_t rz =  points.GetZ()[ipoint];
    lfitY.AddPoint(&rx,ry,1);
    lfitZ.AddPoint(&rx,rz,1);
  }
  lfitY.Eval();
  lfitZ.Eval();
  lfitY.GetParameters(vecY);
  lfitZ.GetParameters(vecZ);
  //
  // 2. Set initial parameters and the covariance matrix
  //
  Int_t ncalibs = fCalibration->GetEntries();
  if (!fLinearParam) {
    fLinearParam = new TMatrixD(ncalibs+4,1);
    fLinearCovar = new TMatrixD(ncalibs+4,ncalibs+4);    
  }
  //
  for (Int_t i=0; i<ncalibs+4;i++){
    (*fLinearParam)(i,0)=0;
    if (i<ncalibs) (*fLinearParam)(i,0) = (*fCalibParam)(i,0);
    for (Int_t j=0; j<ncalibs+4;j++){
      (*fLinearCovar)(i,j)=0;
      if (i<ncalibs&&j<ncalibs) (*fLinearCovar)(i,j) = (*fCalibCovar)(i,j);
    }
  }
  Double_t chi2Y = lfitY.GetChisquare()/lfitY.GetNpoints();
  Double_t chi2Z = lfitZ.GetChisquare()/lfitZ.GetNpoints();
  (*fLinearParam)(ncalibs+0,0) =  lfitY.GetParameter(0);
  (*fLinearCovar)(ncalibs+0,ncalibs+0)= lfitY.GetCovarianceMatrixElement(0,0)*chi2Y*kFac+kOff0;
  (*fLinearParam)(ncalibs+1,0) =  lfitY.GetParameter(1);
  (*fLinearCovar)(ncalibs+1,ncalibs+1)= lfitY.GetCovarianceMatrixElement(1,1)*chi2Y*kFac+kOff1;
  //
  (*fLinearParam)(ncalibs+2,0) =  lfitZ.GetParameter(0);
  (*fLinearCovar)(ncalibs+2,ncalibs+2)= lfitZ.GetCovarianceMatrixElement(0,0)*chi2Z*kFac+kOff0;
  (*fLinearParam)(ncalibs+3,0) =  lfitZ.GetParameter(1);
  (*fLinearCovar)(ncalibs+3,ncalibs+3)= lfitZ.GetCovarianceMatrixElement(1,1)*chi2Z*kFac+kOff1;
  //
  // Fit thetrack together with correction
  //
  AliTrackPoint point;
  for (Int_t ipoint=0; ipoint<npoints-1; ipoint++){
    //
    if (!points.GetPoint(point,ipoint)) continue;
    Double_t erry2 = chi2Y;
    Double_t errz2 = chi2Z;
    //set error - it is hack
    Float_t *cov = (Float_t*) point.GetCov();
    cov[1] = (erry2+kOff0)*scalingRMSY;
    cov[2] = (errz2+kOff0)*scalingRMSZ;
    UpdateLinear(point,debug);
    if (!points.GetPoint(point,npoints-1-ipoint)) continue;
    //set error - it is hack
    cov = (Float_t*) point.GetCov();
    cov[1] = (erry2+kOff0)*scalingRMSY;
    cov[2] = (errz2+kOff0)*scalingRMSZ;
    UpdateLinear(point,debug);
  }
  //
  // save current param and covariance
  for (Int_t i=0; i<ncalibs;i++){
    AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(i);
    transform->SetParam( (*fLinearParam)(i,0));
    (*fCalibParam)(i,0) = (*fLinearParam)(i,0);
    for (Int_t j=0; j<ncalibs;j++){
      (*fCalibCovar)(i,j) = (*fLinearCovar)(i,j);
    }
  }
  if (debug){ // dump debug info
    (*debug)<<"fitLinear"<<
      "vecY.="<<&vecY<<
      "vecZ.="<<&vecZ<<
      "chi2Y="<<chi2Y<<
      "chi2Z="<<chi2Z<<
      "lP.="<<fLinearParam<<
      "cP.="<<fCalibParam<<
      "lC.="<<fLinearCovar<<
      "cC.="<<fCalibCovar<<
      "\n";
  }
}

void AliTPCkalmanFit::DumpTrackLinear(AliTrackPointArray &points, TTreeSRedirector *debug){
  //
  // Dump the track parameters before and after current calibration
  //
  // Track is divided to two parts - 
  // X mean is defined as middle point 
  //

  if (!fCalibParam) {
    AliError("Kalman Fit not initialized");
    return;
  }
  //
  // 
  //
  TLinearFitter *fitters[16];
  TVectorD      *params[16];
  TVectorD      *errs[16];
  TVectorD      chi2N(16);
  Int_t npoints = points.GetNPoints();
  AliTrackPointArray pointsTrans(points);
  ApplyCalibration(&pointsTrans,-1.);
  for (Int_t ifit=0; ifit<8;ifit++){
    fitters[ifit]  = new TLinearFitter(2,"pol1");
    params[ifit]   = new TVectorD(2);
    fitters[ifit+8]= new TLinearFitter(3,"pol2");
    params[ifit+8] = new TVectorD(3);
    errs[ifit]     = new TVectorD(2);
    errs[ifit+8]   = new TVectorD(3);
  }
  //
  // calculate mean x point  and corrdinate frame
  //
  fCurrentAlpha = TMath::ATan2(points.GetY()[npoints-1]-points.GetY()[0], points.GetX()[npoints-1]-points.GetX()[0]);  
  fCA = TMath::Cos(fCurrentAlpha);
  fSA = TMath::Sin(fCurrentAlpha);
  Double_t xmean=0, sum=0;
  for (Int_t ipoint=0; ipoint<npoints-1; ipoint++){
    Double_t rx  =   fCA*points.GetX()[ipoint]+fSA*points.GetY()[ipoint];
    xmean+=rx;
    sum++;
  }
  xmean/=sum;
  //
  for (Int_t ipoint=0; ipoint<npoints-1; ipoint++){
    Double_t rx  =   fCA*points.GetX()[ipoint]+fSA*points.GetY()[ipoint];
    Double_t ry  =  -fSA*points.GetX()[ipoint]+fCA*points.GetY()[ipoint];
    Double_t rz  =  points.GetZ()[ipoint];
    //
    Double_t rxT =   fCA*pointsTrans.GetX()[ipoint]+fSA*pointsTrans.GetY()[ipoint];
    Double_t ryT =  -fSA*pointsTrans.GetX()[ipoint]+fCA*pointsTrans.GetY()[ipoint];
    Double_t rzT =  pointsTrans.GetZ()[ipoint];
    if (rx>xmean){
      fitters[0]->AddPoint(&rxT, ryT,1);
      fitters[2]->AddPoint(&rxT, rzT,1);
      fitters[4]->AddPoint(&rx, ry,1);
      fitters[6]->AddPoint(&rx, rz,1);
      fitters[8]->AddPoint(&rxT, ryT,1);
      fitters[10]->AddPoint(&rxT, rzT,1);
      fitters[12]->AddPoint(&rx, ry,1);
      fitters[14]->AddPoint(&rx, rz,1);
    }else{
      fitters[1]->AddPoint(&rxT, ryT,1);
      fitters[3]->AddPoint(&rxT, rzT,1);
      fitters[5]->AddPoint(&rx, ry,1);
      fitters[7]->AddPoint(&rx, rz,1);
      fitters[9]->AddPoint(&rxT, ryT,1);
      fitters[11]->AddPoint(&rxT, rzT,1);
      fitters[13]->AddPoint(&rx, ry,1);
      fitters[15]->AddPoint(&rx, rz,1);
    }
  }
  for (Int_t ifit=0;ifit<16;ifit++){
    fitters[ifit]->Eval();
    fitters[ifit]->GetParameters(*(params[ifit]));
    fitters[ifit]->GetErrors(*(errs[ifit]));
    chi2N[ifit] = TMath::Sqrt(fitters[ifit]->GetChisquare()/(fitters[ifit]->GetNpoints()-1));
    (*(errs[ifit]))[0]*=chi2N[ifit];
    (*(errs[ifit]))[1]*=chi2N[ifit];
    if (ifit>=8) (*(errs[ifit]))[2]*=chi2N[ifit];  // second derivative
  }
  if (debug){
    (*debug)<<"dumpLinear"<<
      "alpha="<<fCurrentAlpha<<
      "xmean="<<xmean<<
      "y0T.="<<params[0]<<
      "y1T.="<<params[1]<<
      "z0T.="<<params[2]<<
      "z1T.="<<params[3]<<
      "y0O.="<<params[4]<<
      "y1O.="<<params[5]<<
      "z0O.="<<params[6]<<
      "z1O.="<<params[7]<<
      "y0T2.="<<params[8]<<
      "y1T2.="<<params[9]<<
      "z0T2.="<<params[10]<<
      "z1T2.="<<params[11]<<
      "y0O2.="<<params[12]<<
      "y1O2.="<<params[13]<<
      "z0O2.="<<params[14]<<
      "z1O2.="<<params[15]<<
      "y0TErr.="<<errs[0]<<
      "y1TErr.="<<errs[1]<<
      "z0TErr.="<<errs[2]<<
      "z1TErr.="<<errs[3]<<
      "y0OErr.="<<errs[4]<<
      "y1OErr.="<<errs[5]<<
      "z0OErr.="<<errs[6]<<
      "z1OErr.="<<errs[7]<<
      "y0T2Err.="<<errs[8]<<
      "y1T2Err.="<<errs[9]<<
      "z0T2Err.="<<errs[10]<<
      "z1T2Err.="<<errs[11]<<
      "y0O2Err.="<<errs[12]<<
      "y1O2Err.="<<errs[13]<<
      "z0O2Err.="<<errs[14]<<
      "z1O2Err.="<<errs[15]<<
      "chi2.="<<&chi2N<<
      "\n";
  }
  //
  //                    delta  alpha   y0     z0    ky   kz 
  Double_t x[6]={0,0,0,0,0,0};
  x[1]=fCurrentAlpha;
  x[2]=(*params[0])[0];
  x[3]=(*params[2])[0];
  x[4]=(*params[0])[1];
  x[5]=(*params[2])[1];
  //
  //Delta y
  //
  x[0]= (*params[1])[0]-(*params[0])[0];
  fLinearTrackDelta[0]->Fill(x);
  x[0]/=TMath::Sqrt((*errs[1])[0]*(*errs[1])[0]+(*errs[0])[0]*(*errs[0])[0]);
  fLinearTrackPull[0]->Fill(x);
  //
  //Delta ky
  //
  x[0]= 100*((*params[1])[1]-(*params[0])[1]);
  fLinearTrackDelta[1]->Fill(x);
  x[0]/=100*TMath::Sqrt((*errs[1])[1]*(*errs[1])[1]+(*errs[0])[1]*(*errs[0])[1]);
  fLinearTrackPull[1]->Fill(x);
  //
  // ky2_0
  //
  x[0]= 100*100*((*params[8])[2]);
  fLinearTrackDelta[2]->Fill(x);
  x[0]/=100*100*TMath::Sqrt((*errs[8])[2]*(*errs[8])[2]);
  fLinearTrackPull[2]->Fill(x);
  //
  // ky2_1
  //
  x[0]= 100*100*((*params[9])[2]);
  fLinearTrackDelta[3]->Fill(x);
  x[0]/=100*100*TMath::Sqrt((*errs[9])[2]*(*errs[9])[2]);
  fLinearTrackPull[3]->Fill(x);
  //
  //
  //Delta z
  //
  x[0]= (*params[3])[0]-(*params[2])[0];
  fLinearTrackDelta[4]->Fill(x);
  x[0]/=TMath::Sqrt((*errs[3])[0]*(*errs[3])[0]+(*errs[2])[0]*(*errs[2])[0]);
  fLinearTrackPull[4]->Fill(x);
  //
  //Delta kz
  //
  x[0]= 100*((*params[3])[1]-(*params[2])[1]);
  fLinearTrackDelta[5]->Fill(x);
  x[0]/=100*TMath::Sqrt((*errs[3])[1]*(*errs[3])[1]+(*errs[2])[1]*(*errs[2])[1]);
  fLinearTrackPull[5]->Fill(x);
  //
  // kz2_0
  //
  x[0]= 100*100*((*params[10])[2]);
  fLinearTrackDelta[6]->Fill(x);
  x[0]/=100*100*TMath::Sqrt((*errs[10])[2]*(*errs[10])[2]);
  fLinearTrackPull[6]->Fill(x);
  //
  // kz2_1
  //
  x[0]= 100*100*((*params[11])[2]);
  fLinearTrackDelta[7]->Fill(x);
  x[0]/=100*100*TMath::Sqrt((*errs[11])[2]*(*errs[11])[2]);
  fLinearTrackPull[7]->Fill(x);

  //
  // rms of track
  //
  x[0]= chi2N[0];
  fLinearTrackDelta[8]->Fill(x);
  x[0]= chi2N[1];
  fLinearTrackDelta[9]->Fill(x);
  x[0]= chi2N[2];
  fLinearTrackDelta[10]->Fill(x);
  x[0]= chi2N[3];
  fLinearTrackDelta[11]->Fill(x);
  //


  for (Int_t ifit=0; ifit<8;ifit++){
    delete fitters[ifit];
    delete params[ifit];
    delete fitters[ifit+8];
    delete params[ifit+8];
    delete errs[ifit];
    delete errs[ifit+8];
  }
}


void AliTPCkalmanFit::Propagate(TTreeSRedirector */*debug*/){
  //
  // Propagate the Kalman
  //
}

void AliTPCkalmanFit::AddCovariance(const char * varName, Double_t sigma){
  //
  //
  //
  if (fCalibCovar) return;
  if (!fCalibration) return;
  if (!fCalibration->FindObject(varName)) return;
  Int_t ncalibs = fCalibration->GetEntries();
  TString strVar(varName);
  for (Int_t icalib=0;icalib<ncalibs; icalib++){
    AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
    if (strVar.CompareTo(transform->GetName())==0){
      (*fCalibCovar)(icalib,icalib)+=sigma*sigma;
    }
  }
}


void AliTPCkalmanFit::PropagateTime(Int_t time){
  //
  // Propagate the calibration in time
  // - Increase covariance matrix
  //
  if (!fCalibCovar) return;
  Int_t ncalibs = fCalibration->GetEntries();
  Double_t deltaT = 0;
  if (fLastTimeStamp>0) deltaT = (fLastTimeStamp-time)/(60*60.);  // delta T in hours
  fLastTimeStamp = time;  
  for (Int_t icalib=0;icalib<ncalibs; icalib++){
    AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
    if ((*fCalibCovar)(icalib,icalib)<transform->GetSigmaMax()*transform->GetSigmaMax())
      (*fCalibCovar)(icalib,icalib)+=  transform->GetSigma2Time()*TMath::Abs(deltaT);
  }
}


void  AliTPCkalmanFit::UpdateLinear(AliTrackPoint &point, TTreeSRedirector *debug){
  //
  // Update Kalman
  //
  //
  Int_t ncalibs = fCalibration->GetEntries();
  Int_t kNmeas  = 2; 
  Int_t nelem   = ncalibs+4;
  TMatrixD &vecXk=*fLinearParam;     // X vector
  TMatrixD &covXk=*fLinearCovar;     // X covariance
  //
  TMatrixD matHk(kNmeas,nelem);     // vector to mesurement
  TMatrixD vecYk(kNmeas,1);         // Innovation or measurement residual
  TMatrixD vecZk(kNmeas,1);         // Innovation or measurement residual
  TMatrixD measR(kNmeas,kNmeas);
  TMatrixD matHkT(nelem,kNmeas);    // helper matrix Hk transpose
  TMatrixD matSk(kNmeas,kNmeas);    // Innovation (or residual) covariance
  TMatrixD matKk(nelem,kNmeas);     // Optimal Kalman gain
  TMatrixD covXk2(nelem,nelem);     // helper matrix
  TMatrixD covXk3(nelem,nelem);     // helper matrix
  TMatrixD mat1(nelem,nelem);
  //
  // reset matHk
  for (Int_t iel=0;iel<nelem;iel++){
    for (Int_t ip=0;ip<kNmeas;ip++) {
      matHk(ip,iel)=0; 
    }
  }
  for (Int_t iel0=0;iel0<nelem;iel0++)
    for (Int_t iel1=0;iel1<nelem;iel1++){
      if (iel0!=iel1) mat1(iel0,iel1)=0;
      if (iel0==iel1) mat1(iel0,iel1)=1;
    }
  //
  // fill matrix Hk
  //
  Int_t volId = point.GetVolumeID();
  Double_t gxyz[3]={point.GetX(), point.GetY(),point.GetZ()};
  Double_t rxyz[3]={fCA*gxyz[0]+fSA*gxyz[1],-fSA*gxyz[0]+fCA*gxyz[1] ,point.GetZ()};
  //
  Double_t dxdydz[3]={0,0,0};
  Double_t rdxdydz[3]={0,0,0};
  
  for (Int_t icalib=0;icalib<ncalibs;icalib++){
    AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
    dxdydz[0] = transform->GetDeltaXYZ(0,volId, 1, gxyz[0], gxyz[1], gxyz[2]); 
    dxdydz[1] = transform->GetDeltaXYZ(1,volId, 1, gxyz[0], gxyz[1], gxyz[2]);  
    dxdydz[2] = transform->GetDeltaXYZ(2,volId, 1, gxyz[0], gxyz[1], gxyz[2]);
    rdxdydz[0] =  fCA*dxdydz[0]+fSA*dxdydz[1]; 
    rdxdydz[1] = -fSA*dxdydz[0]+fCA*dxdydz[1]; 
    rdxdydz[2] =  dxdydz[2]; 
    //
    matHk(0,icalib)= rdxdydz[1]-rdxdydz[0]*(*fLinearParam)(ncalibs+1,0);  // shift y + shift x * angleY
    matHk(1,icalib)= rdxdydz[2]-rdxdydz[0]*(*fLinearParam)(ncalibs+3,0);  // shift z + shift x * angleZ
  }
  matHk(0,ncalibs+0)=1;
  matHk(0,ncalibs+1)=rxyz[0];
  matHk(1,ncalibs+2)=1;
  matHk(1,ncalibs+3)=rxyz[0];
  //
  //
  //
  vecZk(0,0) =  rxyz[1];
  vecZk(1,0) =  rxyz[2];
  measR(0,0) = point.GetCov()[1]; measR(0,1)=0;
  measR(1,1) = point.GetCov()[2]; measR(1,0)=0;
  //
  vecYk = vecZk-matHk*vecXk;               // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;    // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;            //  Optimal Kalman gain
  //
  covXk2= (mat1-(matKk*matHk));
  covXk3 =  covXk2*covXk;          
  //
  CheckCovariance(covXk3,100);
  vecXk += matKk*vecYk;                    //  updated vector 
  covXk = covXk3;  
  if (debug){ // dump debug info
    (*debug)<<"updateLinear"<<
      //
      "point.="<<&point<<
      "vecYk.="<<&vecYk<<
      "vecZk.="<<&vecZk<<
      "measR.="<<&measR<<
      "matHk.="<<&matHk<<
      "matHkT.="<<&matHkT<<
      "matSk.="<<&matSk<<
      "matKk.="<<&matKk<<
      "covXk2.="<<&covXk2<<
      "covXk.="<<&covXk<<
      "vecXk.="<<&vecXk<<
      "\n";
  }  
}


AliTrackPointArray * AliTPCkalmanFit::SortPoints(AliTrackPointArray &points){
  //
  //Creates the array  - points sorted according radius - neccessay for kalman fit
  // 
  //
  // 0. choose the frame - rotation angle
  //
  Int_t npoints = points.GetNPoints();
  if (npoints<1) return 0;
  Double_t currentAlpha = TMath::ATan2(points.GetY()[npoints-1]-points.GetY()[0], points.GetX()[npoints-1]-points.GetX()[0]);  
  Double_t ca = TMath::Cos(currentAlpha);
  Double_t sa = TMath::Sin(currentAlpha);
  //
  // 1. sort the points
  //
  Double_t *rxvector = new Double_t[npoints];
  Int_t    *indexes  = new Int_t[npoints];
  for (Int_t ipoint=0; ipoint<npoints-1; ipoint++){
    rxvector[ipoint]=ca*points.GetX()[ipoint]+sa*points.GetY()[ipoint];
  }
  TMath::Sort(npoints, rxvector,indexes,kFALSE);
  AliTrackPoint point;
  AliTrackPointArray *pointsSorted= new AliTrackPointArray(npoints);
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    if (!points.GetPoint(point,indexes[ipoint])) continue;
    pointsSorted->AddPoint(ipoint,&point);
  }
  delete [] rxvector;
  delete [] indexes;
  return pointsSorted;
}



AliTrackPointArray * AliTPCkalmanFit::MakePointArrayLinear(Double_t alpha, Double_t y0, Double_t z0, Double_t ky, Double_t kz, Double_t err){
  //
  //
  //
  AliTrackPointArray array(500);
  Float_t cov[10];  // dummy covariance
  Int_t npoints=0;
  for (Int_t i=0;i<6;i++) cov[i]=0.001;
  for (Int_t i=0;i<500;i++){    
    AliTrackPoint point(0, 0, 0, cov, 0,0,0);
    array.AddPoint(npoints, &point);
    npoints++;
  }
  npoints=0;
  for (Float_t ir = -245; ir<245; ir++){
    //
    //
    if (TMath::Abs(ir)<80) continue;
    Double_t ca = TMath::Cos(alpha);
    Double_t sa = TMath::Sin(alpha);
    Double_t lx = ir;
    Double_t ly = y0+lx*ky+gRandom->Gaus(0,err);
    Double_t lz = z0+lx*kz+gRandom->Gaus(0,err);
    Double_t gx = lx*ca-ly*sa;
    Double_t gy = lx*sa+ly*ca;
    Double_t gz = lz;
    Double_t galpha = TMath::ATan2(gy,gx);
    Int_t isec = TMath::Nint((9*galpha/TMath::Pi()+0.5));
    if (isec<0) isec+=18;
    if (gz<0) isec+=18;
    if (ir>130) isec+=36;
    //
    AliTrackPoint point(gx, gy, gz, cov, isec,0,0);
    array.AddPoint(npoints, &point);
    npoints++;
  }
  AliTrackPointArray *parray = new AliTrackPointArray(npoints);
  AliTrackPoint point;
  for (Int_t i=0;i<npoints;i++){
    array.GetPoint(point,i);
    parray->AddPoint(i,&point);
  }

  return parray;
}

void AliTPCkalmanFit::AddCalibration(AliTPCTransformation * calib){
  //
  // Add calibration
  //
  if (!fCalibration) fCalibration = new TObjArray(2000);
  fCalibration->AddLast(calib);
}

Int_t AliTPCkalmanFit::GetTransformationIndex(const char * trName){
  //
  //
  //
  if (!fCalibration) return -1;
  if (!fCalibration->FindObject(trName)) return -1;
  //
  Int_t ncalibs = fCalibration->GetEntries();
  TString strVar(trName);
  for (Int_t icalib=0;icalib<ncalibs; icalib++){
    AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
    if (strVar.CompareTo(transform->GetName())==0){
      return icalib;
    }
  }
  return -1;
}



void AliTPCkalmanFit::ApplyCalibration(AliTrackPointArray *array, Double_t csign){
  //
  //
  //
  if (!fCalibration) return;
  Int_t ncalibs = fCalibration->GetEntries();
  if (ncalibs==0) return;
  Int_t npoints = array->GetNPoints();
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Int_t volId = array->GetVolumeID()[ipoint];
    Double_t xyz[3]={array->GetX()[ipoint], array->GetY()[ipoint],array->GetZ()[ipoint]};
    Double_t dxdydz[3]={0,0,0};
    for (Int_t icalib=0; icalib<ncalibs; icalib++){
      AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
      dxdydz[0] += transform->GetDeltaXYZ(0,volId, transform->GetParam(), xyz[0], xyz[1], xyz[2]); 
      dxdydz[1] += transform->GetDeltaXYZ(1,volId, transform->GetParam(), xyz[0], xyz[1], xyz[2]);  
      dxdydz[2] += transform->GetDeltaXYZ(2,volId, transform->GetParam(), xyz[0], xyz[1], xyz[2]);  
    }
    ((Float_t*)array->GetX())[ipoint]+=csign*dxdydz[0];
    ((Float_t*)array->GetY())[ipoint]+=csign*dxdydz[1];
    ((Float_t*)array->GetZ())[ipoint]+=csign*dxdydz[2];
  }
}

Bool_t AliTPCkalmanFit::DumpCorelation(Double_t threshold,  const char *mask0, const char *mask1){
  //
  //
  //
  TMatrixD &mat = *fCalibCovar;
  Int_t nrow= mat.GetNrows();
  for (Int_t irow=0; irow<nrow; irow++){
    AliTPCTransformation * trans0 = GetTransformation(irow);
    TString  strName0(trans0->GetName());
    if (mask0){
      if (!strName0.Contains(mask0)) continue;
    }
    for (Int_t icol=irow+1; icol<nrow; icol++){
      AliTPCTransformation * trans1 = GetTransformation(icol);
      TString  strName1(trans1->GetName());
      if (mask1){
	if (!strName1.Contains(mask1)) continue;
      }	
      //
      Double_t diag = TMath::Sqrt(TMath::Abs(mat(irow,irow)*mat(icol,icol))); 
      if (diag<=0){
	printf("Negative covariance\t%d\t%d\t%f\n",irow,icol, mat(irow,icol));
	continue;
      }
      Double_t corr0 = mat(irow,icol)/diag;
      if (TMath::Abs(corr0)>threshold){
	printf("%d\t%d\t%s\t%s\t%f\t%f\t%f\n", irow,icol, trans0->GetName(), trans1->GetName(),
	       TMath::Sqrt(mat(irow,irow)), TMath::Sqrt(mat(icol,icol)), corr0);
      }
    }
  }
  return (nrow>0);   
}

Bool_t AliTPCkalmanFit::DumpCalib(const char *mask, Float_t correlationCut){
  //
  // Print calibration entries - name, value, error
  //
  TMatrixD &mat = *fCalibCovar;
  Int_t nrow= mat.GetNrows();
  TString  strMask(mask);
  TVectorD vecCorrSum(nrow);
  for (Int_t irow=0; irow<nrow; irow++){
    vecCorrSum[irow]=0;
    for (Int_t icol=0; icol<nrow; icol++){
      if (icol==irow) continue;
      Double_t diag = TMath::Sqrt(TMath::Abs(mat(irow,irow)*mat(icol,icol))); 
      Double_t corr0 = mat(irow,icol)/diag;
      vecCorrSum[irow]+=TMath::Abs(corr0);
    }
    vecCorrSum[irow]*=0.5;
  }

  for (Int_t irow=0; irow<nrow; irow++){
    AliTPCTransformation * trans0 = GetTransformation(irow);
    TString  strName(trans0->GetName());
    if (mask){
      if (!strName.Contains(mask)) continue;
    }
    if (vecCorrSum[irow]<correlationCut) continue;
    printf("%d\t%s\t%f\t%f\t%f\n", 
	   irow, 
	   trans0->GetName(),
	   (*fCalibParam)(irow,0),
	   TMath::Sqrt(mat(irow,irow)),  vecCorrSum[irow]);          
  }
  return (nrow>0);
}


Bool_t  AliTPCkalmanFit::CheckCovariance(TMatrixD &mat, Float_t /*maxEl*/){
  //
  // Check consistency of covariance matrix
  // + symetrize coavariance matrix
  Bool_t isOK=kTRUE;
  Int_t nrow= mat.GetNrows();
  for (Int_t irow=0; irow<nrow; irow++){
    if (mat(irow,irow)<=0){
      printf("Negative covariance\t%d\t%f\n",irow,mat(irow,irow));
      isOK=kFALSE;
    }
//     if (mat(irow,irow)>maxEl){
//       printf("Too big  covariance\t%d\t%f\n",irow,mat(irow,irow));
//       isOK=kFALSE;
//     }
    
    for (Int_t icol=0; icol<nrow; icol++){
      //      if (mat(irow,irow)
      Double_t diag = TMath::Sqrt(TMath::Abs(mat(irow,irow)*mat(icol,icol))); 
      if (diag<=0){
	printf("Negative covariance\t%d\t%d\t%f\n",irow,icol, mat(irow,icol));
	isOK=kFALSE;
	continue;
      }
      Double_t cov0 = mat(irow,icol)/diag;
      Double_t cov1 = mat(icol,irow)/diag;
      if (TMath::Abs(cov0)>1 || TMath::Abs(cov1)>1 ){
	printf("Covariance Problem %d\t%d\t%f\t%f\n",irow,icol, cov0, cov1);
	isOK=kFALSE;
      }
      if (TMath::Abs(cov0-cov1)>0.0000001){
	printf("Asymetry problem %d\t%d\t%f\t%f\n",irow,icol, cov0, cov1);
	isOK=kFALSE;
      }
      //
      // symetrize the covariance matrix
      Double_t mean = (mat(irow,icol)+ mat(icol,irow))*0.5;
      mat(irow,icol)=mean;
      mat(icol,irow)=mean;
    }
  }
  return isOK;
}


AliTPCkalmanFit *  AliTPCkalmanFit::Test(Int_t ntracks){
  //
  // This is test example
  //

  //
  // 1. Setup transformation
  //


  TVectorD fpar(10);
  AliTPCTransformation * transformation=0;
  AliTPCkalmanFit * kalmanFit0 =  new AliTPCkalmanFit;
  AliTPCkalmanFit * kalmanFit2 =  new AliTPCkalmanFit;
  //
  // Radial scaling
  //
  for (Int_t iside=0; iside<=1; iside++)
    for (Int_t ipar0=0; ipar0<3; ipar0++)
      for (Int_t ipar1=0; ipar1<3; ipar1++){
	fpar[0]=ipar0; 
	fpar[1]=ipar1;
	if (ipar0+ipar1==0) continue;
	Double_t param = (gRandom->Rndm()-0.5)*0.5;  // generate random parameters
	char tname[100];
	snprintf(tname,100,"tscalingR%d%dSide%d",ipar0,ipar1,iside);
	transformation = new AliTPCTransformation(tname,AliTPCTransformation::BitsSide(iside),"TPCscalingRPol",0,0,1);
	transformation->SetParams(0,5*0.25,0,&fpar);
	kalmanFit0->AddCalibration(transformation);
	transformation = new AliTPCTransformation(tname,AliTPCTransformation::BitsSide(iside),"TPCscalingRPol",0,0,1);
	transformation->SetParams(param,0.25,0,&fpar);
	kalmanFit2->AddCalibration(transformation);
      }


  //
  // 2. Init transformation
  //
  kalmanFit2->Init();    
  kalmanFit0->Init();
  
  //
  // 3. Run kalman filter
  //
  Int_t ncalibs = kalmanFit0->fCalibration->GetEntries();
  TVectorD err(ncalibs);
  TTreeSRedirector *pcstream = new TTreeSRedirector("kalmanfitTest.root");
  for (Int_t i=0;i<ntracks;i++){
    if (i%100==0) printf("%d\n",i);
     Double_t alpha = gRandom->Rndm()*TMath::TwoPi();
     Double_t   y0  = (gRandom->Rndm()-0.5)*180; 
     Double_t   z0  = (gRandom->Rndm()-0.5)*250*2; 
     Double_t   ky  = (gRandom->Rndm()-0.5)*1; 
     Double_t   kz  = (gRandom->Rndm()-0.5)*1;
     //generate random TPC track
     AliTrackPointArray * array  = AliTPCkalmanFit::MakePointArrayLinear(alpha,y0,z0, ky, kz,0.04);
     AliTrackPointArray * arrayB = new AliTrackPointArray(*array);  // backup
     kalmanFit2->ApplyCalibration(array,1.);  // misalign ideal track
     for (Int_t icalib=0; icalib<ncalibs; icalib++){
       err[icalib] = TMath::Sqrt((*kalmanFit0->fCalibCovar)(icalib,icalib));
     }
     //
     (*pcstream)<<"dump0"<<
       "alpha="<<alpha<<
       "y0="<<y0<<
       "z0="<<z0<<
       "ky="<<ky<<
       "lz="<<kz<<
       "p.="<<array<<
       "pB.="<<arrayB<<
       "cparam.="<<kalmanFit0->fCalibParam<<
       "ccovar.="<<kalmanFit0->fCalibCovar<<
       "err.="<<&err<<
       "gparam.="<<kalmanFit2->fCalibParam<<
       "gcovar.="<<kalmanFit2->fCalibCovar<<
       "\n";
     if (i%20==0) {
       kalmanFit0->FitTrackLinear(*array,pcstream); // fit track - dump intermediate results
     }else{
       kalmanFit0->FitTrackLinear(*array,0);        // fit track + calibration
     }
     kalmanFit0->DumpTrackLinear(*array,pcstream);    // dump track residuals to the tree + fill histograms
  }
  pcstream->GetFile()->cd();
  kalmanFit0->Write("kalmanFit");
  delete pcstream;
  return kalmanFit0;
}


// Double_t AliTPCkalmanFit::GetTPCDeltaXYZ(Int_t coord, Int_t volID, Double_t x, Double_t y, Double_t z){
//   //
//   // function for visualization purposes
//   //
//   if (!fCalibration) return 0;
//   Int_t ncalibs = fCalibration->GetEntries();
//   if (ncalibs==0) return 0;
//   Double_t dxdydz[3]={0,0,0};
//   //
//   if (volID<0){
//     Double_t alpha       = TMath::ATan2(y,x);
//     Double_t r           = TMath::Sqrt(y*y+x*x);
//     volID                = TMath::Nint(9*alpha/TMath::Pi()-0.5);
//     if (volID<0) volID+=18;
//     if (z<0) volID+=18;
//     if (r>120) volID+=36;
//   }
//   for (Int_t icalib=0; icalib<ncalibs; icalib++){
//     AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
//     Double_t param = (*fCalibParam)(icalib,0);
//     dxdydz[coord] += transform->GetDeltaXYZ(coord,volID, param, x, y,z); 
//   }
//   return dxdydz[coord];
// }

// Double_t AliTPCkalmanFit::SGetTPCDeltaXYZ(Int_t coord, Int_t volID, Double_t x, Double_t y, Double_t z){
//   //
//   //
//   //
//   if (AliTPCkalmanFit::fgInstance==0) return 0;
//   return AliTPCkalmanFit::fgInstance->GetTPCDeltaXYZ(coord, volID,x,y,z);
// }





Double_t AliTPCkalmanFit::GetTPCDeltaXYZ(Int_t coord, Int_t volID, Int_t icoordsys, Double_t x, Double_t y, Double_t z){
  //
  // function for visualization purposes
  //
  // coord - coordinate for output
  //       - 0 -X
  //         1 -Y
  //         2 -Z
  //         3 -R
  //         4 -RPhi
  //         5 -Z
  //
  //icoordsys - type of coordinate system for input
  //         - 0  - x,y,z
  //         - 1  - r,phi,z
  //
  if (!fCalibration) return 0;
  Int_t ncalibs = fCalibration->GetEntries();
  if (ncalibs==0) return 0;
  Double_t xyz[3]={0,0,0};
  Double_t dxdydz[6]={0,0,0,0,0,0};
  Double_t alpha=0;
  Double_t r=0;
  if(icoordsys==0){alpha=TMath::ATan2(y,x); r =TMath::Sqrt(y*y+x*x);}
  if(icoordsys==1){alpha=y; r=x;}
  Double_t ca    = TMath::Cos(alpha);
  Double_t sa    = TMath::Sin(alpha);
  if(icoordsys==0){xyz[0]=x; xyz[1]=y; xyz[2]=z;}
  if(icoordsys==1){xyz[0]=x*ca; xyz[1]=x*sa; xyz[2]=z;}
  //
  if (volID<0){
    // Double_t alpha       = TMath::ATan2(y,x);
    //Double_t r           = TMath::Sqrt(y*y+x*x);
    volID                = TMath::Nint(9*alpha/TMath::Pi()-0.5);
    if (volID<0) volID+=18;
    if (z<0) volID+=18;
    if (r>120) volID+=36;
  }
  for (Int_t icalib=0; icalib<ncalibs; icalib++){
    AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
    Double_t param = (*fCalibParam)(icalib,0);
    for (Int_t icoord=0;icoord<6;icoord++){
      dxdydz[icoord] += transform->GetDeltaXYZ(icoord,volID, param, xyz[0],xyz[1],xyz[2]);
    } 
  }

  return dxdydz[coord];
}


Double_t AliTPCkalmanFit::SGetTPCDeltaXYZ(Int_t coord, Int_t volID, Int_t icoordsys, Double_t x, Double_t y, Double_t z){
  //
  //
  //
  if (AliTPCkalmanFit::fgInstance==0) return 0;
  return AliTPCkalmanFit::fgInstance->GetTPCDeltaXYZ(coord, volID, icoordsys,x,y,z);
}





Double_t AliTPCkalmanFit::GetTPCtransXYZ(Int_t coord, Int_t volID, Int_t calibID, Int_t icoordsys, Double_t x, Double_t y, Double_t z){

  Int_t ncalibs = fCalibration->GetEntries();
  if (calibID>=ncalibs) return 0;
  //Int_t volID=-1;
  //Double_t xyz[3]={x,y,z};
  Double_t r=0;
  Double_t alpha=0;
  if(icoordsys==0){r=TMath::Sqrt(x*x+y*y); alpha = TMath::ATan2(y,x);}
  if(icoordsys==1){r=x; alpha = y;}
  Double_t ca    = TMath::Cos(alpha);
  Double_t sa    = TMath::Sin(alpha);
  Double_t xyz[3]={0,0,0};
  if(icoordsys==0){xyz[0]=x;xyz[1]=y;xyz[2]=z;} 
  if(icoordsys==1){xyz[0]=x*ca; xyz[1]=x*sa; xyz[2]=z;}
  //xyz[3]=param; xyz[4]=volID;

  if (volID<0){
    //Double_t alpha       = TMath::ATan2(xyz[1],xyz[0]);
    //Double_t r           = TMath::Sqrt(xyz[1]*xyz[1]+xyz[0]*xyz[0]);
    volID                = TMath::Nint(9*alpha/TMath::Pi()-0.5);
    if (volID<0) volID+=18;
    if (xyz[2]<0) volID+=18;
    if (r>120) volID+=36;
  }
  AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(calibID);
  //transform->SetInstance(transform);
  Double_t param = (*fCalibParam)(calibID,0);
  Double_t delta = (Double_t)transform->GetDeltaXYZ(coord,volID, param, xyz[0],xyz[1],xyz[2]);
      
  return delta;
}

Double_t AliTPCkalmanFit::SGetTPCtransXYZ(Int_t coord, Int_t volID, Int_t calibID, Int_t icoordsys, Double_t x, Double_t y, Double_t z){
  //
  //
  //
  if (AliTPCkalmanFit::fgInstance==0) return 0;
  return AliTPCkalmanFit::fgInstance->GetTPCtransXYZ(coord, volID, calibID,icoordsys,x,y,z);
}


void AliTPCkalmanFit::MakeTreeTrans(TTreeSRedirector *debug, const char *treeName){
  //
  // Make the Tree before and after current calibration
  //
  if (!fCalibParam) {
    AliError("Kalman Fit not initialized");
    return;
  }
  //
  //
  //
  const Int_t ncalibs = fCalibration->GetEntries();
  TMatrixD dxdydz(ncalibs,5);
  Double_t * adx    = new Double_t[ncalibs];
  Double_t * ady    = new Double_t[ncalibs];
  Double_t * adz    = new Double_t[ncalibs];
  Double_t * adr    = new Double_t[ncalibs];
  Double_t * adrphi = new Double_t[ncalibs];

  Double_t x[3];
  for (x[0]=-250.;x[0]<=250.;x[0]+=10.){
    for (x[1]=-250.;x[1]<=250.;x[1]+=10.){
      for (x[2]=-250.;x[2]<=250.;x[2]+=20.) {
	Double_t r=TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
	if (r<20) continue;
	if (r>260) continue;
	//Double_t z = x[2];
	Double_t phi=TMath::ATan2(x[1],x[0]);
	Double_t ca=TMath::Cos(phi);
	Double_t sa=TMath::Sin(phi);
	Double_t dx=0;
	Double_t dy=0;
	Double_t dz=0;
	Double_t dr=0;
	Double_t rdphi=0;

	Int_t volID= TMath::Nint(9*phi/TMath::Pi()-0.5);
        if (volID<0) volID+=18;
        if (x[2]<0) volID+=18; //C-side
        if (r>120) volID+=36; //outer
        Double_t volalpha=(volID+0.5)*TMath::Pi()/9;
        Double_t cva=TMath::Cos(volalpha);
        Double_t sva=TMath::Sin(volalpha);

        Double_t lx=x[0]*cva+x[1]*sva;
        Double_t ly=-x[0]*sva+x[1]*cva;


	for(Int_t icalib=0;icalib<ncalibs;icalib++){
	  for(Int_t icoord=0;icoord<5;icoord++){
	    dxdydz(icalib,icoord)= GetTPCtransXYZ(icoord, -1, icalib, 0, x[0], x[1], x[2]);
	  }
	}
	dx=GetTPCDeltaXYZ(0, -1, 0, x[0], x[1], x[2]);
	dy=GetTPCDeltaXYZ(1, -1, 0, x[0], x[1], x[2]);
	dz=GetTPCDeltaXYZ(2, -1, 0, x[0], x[1], x[2]);
	dr=GetTPCDeltaXYZ(3, -1, 0, x[0], x[1], x[2]);
	rdphi=GetTPCDeltaXYZ(4, -1, 0, x[0], x[1], x[2]);


	if(debug){
	  TTreeStream &cstream=
	    (*debug)<<treeName<<
	    "x="<<x[0]<<
	    "y="<<x[1]<<
	    "z="<<x[2]<<
	    "r="<<r<<
	    "ca="<<ca<<
	    "sa="<<sa<<
            "lx="<<lx<<
            "ly="<<ly<<
            "sector="<<volID<<
 	    "phi="<<phi<<
	    "dx="<<dx<<
	    "dy="<<dy<<
	    "dz="<<dz<<
	    "dr="<<dr<<
	    "rdphi="<<rdphi<<
	    "dxdydz.="<<&dxdydz;
	    for (Int_t icalib=0;icalib<ncalibs;icalib++){
	      AliTPCTransformation * transform = (AliTPCTransformation *)fCalibration->At(icalib);
	      char tname[1000];
	      //
	      snprintf(tname,1000,"dx%s=",transform->GetName());
	      adx[icalib] =dxdydz(icalib,0); 
	      cstream<<tname<<adx[icalib];
	      snprintf(tname,1000,"dy%s=",transform->GetName());
	      ady[icalib] =dxdydz(icalib,1); 
	      cstream<<tname<<ady[icalib];
	      snprintf(tname,1000,"dz%s=",transform->GetName());
	      adz[icalib] =dxdydz(icalib,2); 
	      cstream<<tname<<adz[icalib];
	      //
	      snprintf(tname,1000,"dr%s=",transform->GetName());
	      adr[icalib] =dxdydz(icalib,3); 
	      cstream<<tname<<adr[icalib];
	      snprintf(tname,1000,"rdphi%s=",transform->GetName());
	      adrphi[icalib] =dxdydz(icalib,4); 
	      cstream<<tname<<adrphi[icalib];
	    }
	    cstream<<"\n";
	}
      }
    }
    Printf("x0=%f finished",x[0]);
  }
  
}
