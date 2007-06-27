
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH3F.h>
//
#include <TPDGCode.h>
#include <TStyle.h>
#include "TLinearFitter.h"
#include "TMatrixD.h"
#include "TTreeStream.h"
#include "TF1.h"



#include "AliMagF.h"
#include "AliTracker.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h" 
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "AliTPCROC.h"


#include "AliTPCParamSR.h"
#include "AliTPCClusterParam.h"
#include "AliTrackPointArray.h"
#include "TCint.h"
#include "AliTPCcalibTracks.h"

ClassImp(AliTPCcalibTracks)

AliTPCParam  param;



AliTPCcalibTracks::AliTPCcalibTracks() : 
   TNamed(),
   fHclus(0),
   fFileNo(0)     
 {
   G__SetCatchException(0);     
   param.Update();
   TFile fparam("/u/miranov/TPCClusterParam.root");
   fClusterParam  =  (AliTPCClusterParam *) fparam.Get("Param"); 
   if (fClusterParam){
     //fClusterParam->SetInstance(fClusterParam);
   }else{
     printf("Cluster Param not found\n");
   } 
   fDebugStream = new TTreeSRedirector("TPCSelectorDebug.root");
 }   


void    AliTPCcalibTracks::ProcessTrack(AliTPCseed * seed){
}


Int_t   AliTPCcalibTracks::GetBin(Float_t q, Int_t pad){
  //
  // calculate bins for given q and pad type 
  // used in TObjArray
  //
  Int_t   res  = TMath::Max(TMath::Nint((TMath::Sqrt(q)-3.)),0);  
  res*=3;
  res+=pad;
  return res;
}

Int_t   AliTPCcalibTracks::GetBin(Int_t iq, Int_t pad){
  //
  // calculate bins for given iq and pad type 
  // used in TObjArray
  //
  return iq*3+pad;;
}

Float_t AliTPCcalibTracks::GetQ(Int_t bin){
  Int_t bin0 = bin/3;
  bin0+=3;
  return bin0*bin0;
}

 





void AliTPCcalibTracks::ProofSlaveBegin(TList * output)
{
  // Called on PROOF - fill output list
  //fChain = tree;
  //Init(tree);

  char chname[1000];
  TProfile * prof1=0;
  TH1F     * his1 =0;
  fHclus = new TH1I("hclus","Number of clusters",100,0,200);
  output->AddLast(fHclus);


  //
  // Amplitude  - sector -row histograms 
  //
  fArrayAmpRow = new TObjArray(72);
  fArrayAmp    = new TObjArray(72);
  for (Int_t i=0; i<36; i++){   
    sprintf(chname,"Amp_row_Sector%d",i);
    prof1 = new TProfile(chname,chname,63,0,64);
    prof1->SetXTitle("Pad row");
    prof1->SetYTitle("Mean Max amplitude");
    fArrayAmpRow->AddAt(prof1,i);
    output->AddLast(prof1);
    sprintf(chname,"Amp_row_Sector%d",i+36);
    prof1 = new TProfile(chname,chname,96,0,97);
    prof1->SetXTitle("Pad row");
    prof1->SetYTitle("Mean Max  amplitude");
    fArrayAmpRow->AddAt(prof1,i+36);
    output->AddLast(prof1);
    //
    // amplitude
    sprintf(chname,"Amp_Sector%d",i);
    his1 = new TH1F(chname,chname,250,0,500);
    his1->SetXTitle("Max Amplitude (ADC)");
    fArrayAmp->AddAt(his1,i);
    output->AddLast(his1);
    sprintf(chname,"Amp_Sector%d",i+36);
    his1 = new TH1F(chname,chname,200,0,600);
    his1->SetXTitle("Max Amplitude (ADC)");
    fArrayAmp->AddAt(his1,i+36);
    output->AddLast(his1);
    //
  }

  fDeltaY = new TH1F("DeltaY","DeltaY",100,-1,1);
  fDeltaZ = new TH1F("DeltaZ","DeltaZ",100,-1,1);
  output->AddLast(fDeltaY);
  output->AddLast(fDeltaZ);

  fResolY = new TObjArray(3);
  fResolZ = new TObjArray(3);
  fRMSY   = new TObjArray(3);
  fRMSZ   = new TObjArray(3);
  TH3F * his3D;
  //
  his3D = new TH3F("Resol Y0","Resol Y0", 5,20,250, 4, 0,1., 50, -1,1);
  fResolY->AddAt(his3D,0);	
  output->AddLast(his3D);
  his3D = new TH3F("Resol Y1","Resol Y1", 5,20,250, 4, 0,1., 50, -1,1);
  fResolY->AddAt(his3D,1);
  output->AddLast(his3D);
  his3D = new TH3F("Resol Y2","Resol Y2", 5,20,250, 4, 0,0.8, 50, -1,1);
  fResolY->AddAt(his3D,2);
  output->AddLast(his3D);
  //
  his3D = new TH3F("Resol Z0","Resol Z0", 5,20,250, 4, 0,1, 50, -1,1);
  fResolZ->AddAt(his3D,0);
  output->AddLast(his3D);
  his3D = new TH3F("Resol Z1","Resol Z1", 5,20,250, 4, 0,1, 50, -1,1);
  fResolZ->AddAt(his3D,1);
  output->AddLast(his3D);
  his3D = new TH3F("Resol Z2","Resol Z2", 5,20,250, 4, 0,1, 50, -1,1);
  fResolZ->AddAt(his3D,2);
  output->AddLast(his3D);
  //
  his3D = new TH3F("RMS Y0","RMS Y0", 5,20,250, 4, 0,1., 50, 0,0.8);
  fRMSY->AddAt(his3D,0);
  output->AddLast(his3D);
  his3D = new TH3F("RMS Y1","RMS Y1", 5,20,250, 4, 0,1., 50, 0,0.8);
  fRMSY->AddAt(his3D,1);
  output->AddLast(his3D);
  his3D = new TH3F("RMS Y2","RMS Y2", 5,20,250, 4, 0,0.8, 50, 0,0.8);
  fRMSY->AddAt(his3D,2);
  output->AddLast(his3D);
  //
  his3D = new TH3F("RMS Z0","RMS Z0", 5,20,250, 4, 0,1, 50, 0,0.8);
  fRMSZ->AddAt(his3D,0);
  output->AddLast(his3D);
  his3D = new TH3F("RMS Z1","RMS Z1", 5,20,250, 4, 0,1, 50, 0,0.8);
  fRMSZ->AddAt(his3D,1);
  output->AddLast(his3D);
  his3D = new TH3F("RMS Z2","RMS Z2", 5,20,250, 4, 0,1, 50, 0,0.8);
  fRMSZ->AddAt(his3D,2);
  output->AddLast(his3D);
  //
  fArrayQDY = new TObjArray(300);
  fArrayQDZ = new TObjArray(300);
  fArrayQRMSY = new TObjArray(300);
  fArrayQRMSZ = new TObjArray(300);
  for (Int_t iq=0; iq<10; iq++){
    for (Int_t ipad=0; ipad<3; ipad++){
      Int_t   bin   = GetBin(iq,ipad);
      Float_t qmean = GetQ(bin);
      char name[200];
      sprintf(name,"ResolY Pad%d Qmiddle%f",ipad, qmean);
      his3D = new TH3F(name, name, 20,10,250, 20, 0,1.5, 50, -1,1);
      fArrayQDY->AddAt(his3D,bin);
      output->AddLast(his3D);
      sprintf(name,"ResolZ Pad%d Qmiddle%f",ipad, qmean);
      his3D = new TH3F(name, name, 20,10,250, 20, 0,1.5, 50, -1,1);
      fArrayQDZ->AddAt(his3D,bin);
      output->AddLast(his3D);
      //
      sprintf(name,"RMSY Pad%d Qmiddle%f",ipad, qmean);
      his3D = new TH3F(name, name, 20,10,250, 20, 0,1.5, 50, 0,1);
      fArrayQRMSY->AddAt(his3D,bin);
      output->AddLast(his3D);
      sprintf(name,"RMSZ Pad%d Qmiddle%f",ipad, qmean);
      his3D = new TH3F(name, name, 20,10,250, 20, 0,1.5, 50, 0,1);
      fArrayQRMSZ->AddAt(his3D,bin);
      output->AddLast(his3D);
    }
  }  
}




Float_t AliTPCcalibTracks::TPCBetheBloch(Float_t bg)
{
 //
 // Bethe-Bloch energy loss formula
 //
 const Double_t kp1=0.76176e-1;
 const Double_t kp2=10.632;
 const Double_t kp3=0.13279e-4;
 const Double_t kp4=1.8631;
 const Double_t kp5=1.9479;
 Double_t dbg = (Double_t) bg;
 Double_t beta = dbg/TMath::Sqrt(1.+dbg*dbg);
 Double_t aa = TMath::Power(beta,kp4);
 Double_t bb = TMath::Power(1./dbg,kp5);
 bb=TMath::Log(kp3+bb);
 return ((Float_t)((kp2-aa-bb)*kp1/aa));
}


Bool_t AliTPCcalibTracks::AcceptTrack(AliTPCseed * track){
  //
  //
  //
  const Int_t   kMinClusters  = 20;
  const Float_t kMinRatio     = 0.4;
  const Float_t kMax1pt       = 0.5;
  const Float_t kEdgeYXCutNoise    = 0.13;
  const Float_t kEdgeThetaCutNoise = 0.018;
  //
  // edge induced noise tracks - NEXT RELEASE will be removed during tracking
  if (TMath::Abs(track->GetY()/track->GetX())> kEdgeYXCutNoise)
    if (TMath::Abs(track->GetTgl())<kEdgeThetaCutNoise) return kFALSE;
  
  //
  if (track->GetNumberOfClusters()<kMinClusters) return kFALSE;
  Float_t ratio = track->GetNumberOfClusters()/(track->GetNFoundable()+1.);
  if (ratio<kMinRatio) return kFALSE;
  Float_t mpt = track->Get1Pt();
  if (TMath::Abs(mpt)>kMax1pt) return kFALSE;
  //if (TMath::Abs(track->GetZ())>240.) return kFALSE;
  //if (TMath::Abs(track->GetZ())<10.) return kFALSE;
  //if (TMath::Abs(track->GetTgl())>0.03) return kFALSE;
  
  return kTRUE;
}

void AliTPCcalibTracks::FillHistoCluster(AliTPCseed * track){
  //
  //
  //
  const Int_t kFirstLargePad = 127;
  const Float_t kLargePadSize = 1.5;
  for (Int_t irow=0; irow<159; irow++){
    AliTPCclusterMI * cluster = track->GetClusterPointer(irow);
    if (!cluster) continue;
    Int_t sector = cluster->GetDetector();
    if (cluster->GetQ()<=0) continue;
    Float_t max = cluster->GetMax();
    printf ("irow, kFirstLargePad = %d, %d \n",irow,kFirstLargePad);
    if ( irow >= kFirstLargePad) {
      max /= kLargePadSize;
    }
    TProfile *profAmpRow =  (TProfile*)fArrayAmpRow->At(sector);
    profAmpRow->Fill(cluster->GetRow(), max);
  }  
}

void  AliTPCcalibTracks::FillResolutionHistoLocal(AliTPCseed * track){
  //
  // fill resolution histograms - localy - trcklet in the neighborhood
  //
  const Int_t   kDelta  = 10;      // delta rows to fit
  const Float_t kMinRatio = 0.75;  // minimal ratio
  const Float_t kCutChi2 = 6.;     // cut chi2 - left right  - kink removal
  const Float_t kErrorFraction = 0.5;  // use only clusters with small intrpolation error - for error param
  const Int_t   kFirstLargePad = 127;
  const Float_t kLargePadSize = 1.5;
  static TLinearFitter fitterY2(3,"pol2");
  static TLinearFitter fitterZ2(3,"pol2");
  static TLinearFitter fitterY0(2,"pol1");
  static TLinearFitter fitterZ0(2,"pol1");
  static TLinearFitter fitterY1(2,"pol1");
  static TLinearFitter fitterZ1(2,"pol1");
  TVectorD      paramY0(2);
  TVectorD      paramY1(2);
  TVectorD      paramY2(3);
  TVectorD      paramZ0(2);
  TVectorD      paramZ1(2);
  TVectorD      paramZ2(3);
  TMatrixD    matrixY0(2,2);
  TMatrixD    matrixZ0(2,2);
  TMatrixD    matrixY1(2,2);
  TMatrixD    matrixZ1(2,2);
  //
  // estimate mean error
  //
  Int_t nTrackletsAll  = 0;
  Int_t nClusters      = 0;
  Float_t csigmaY       = 0;
  Float_t csigmaZ       = 0;
  Int_t sectorG       = -1;
  for (Int_t irow=0; irow<159; irow++){
    AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
    if (!cluster0) continue;
    Int_t sector = cluster0->GetDetector();
    if (sector!=sectorG){
      nClusters=0;
      fitterY2.ClearPoints();
      fitterZ2.ClearPoints();
      sectorG=sector;
    }else{
      nClusters++;
      Double_t x = cluster0->GetX();
      fitterY2.AddPoint(&x,cluster0->GetY(),1);
      fitterZ2.AddPoint(&x,cluster0->GetZ(),1);
      //
      if (nClusters>=kDelta+3){
	fitterY2.Eval();
	fitterZ2.Eval();
	nTrackletsAll++;
	csigmaY+=fitterY2.GetChisquare()/(nClusters-3.);
	csigmaZ+=fitterZ2.GetChisquare()/(nClusters-3.);
	nClusters=-1;
	fitterY2.ClearPoints();
	fitterZ2.ClearPoints();
      }
    }
  }
  csigmaY = TMath::Sqrt(csigmaY/nTrackletsAll);
  csigmaZ = TMath::Sqrt(csigmaZ/nTrackletsAll);
  //
  //
  //
  for (Int_t irow=0; irow<159; irow++){
    Int_t nclFound     = 0;
    Int_t nclFoundable = 0;
    AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
    if (!cluster0) continue;
    Int_t sector = cluster0->GetDetector();
    Float_t xref = cluster0->GetX();
    //
    // check the neighborhood occupancy - (Delta ray - noise removal)
    //    
    for (Int_t idelta= -kDelta; idelta<=kDelta; idelta++){
      if (idelta==0) continue;
      if (idelta+irow<0) continue;
      if (idelta+irow>159) continue;
      AliTPCclusterMI * clusterD = track->GetClusterPointer(irow);
      if ( clusterD && clusterD->GetDetector()!= sector) continue;
      if (clusterD->GetType()<0) continue;      
      nclFoundable++;
      if (clusterD) nclFound++;
    }
    if (nclFound<kDelta*kMinRatio) continue;    
    if (Float_t(nclFound)/Float_t(nclFoundable)<kMinRatio) continue;
    //
    // Make Fit
    //
    fitterY2.ClearPoints();
    fitterZ2.ClearPoints();
    fitterY0.ClearPoints();
    fitterZ0.ClearPoints();
    fitterY1.ClearPoints();
    fitterZ1.ClearPoints();
    //
    nclFound=0;
    Int_t ncl0=0;
    Int_t ncl1=0;
    for (Int_t idelta=-kDelta; idelta<=kDelta; idelta++){
      if (idelta==0) continue;
      if (idelta+irow<0) continue;
      if (idelta+irow>159) continue;
      AliTPCclusterMI * cluster = track->GetClusterPointer(irow+idelta);
      if (!cluster) continue;
      if (cluster->GetType()<0) continue;
      if (cluster->GetDetector()!=sector) continue;
      Double_t x = cluster->GetX()-xref;
      nclFound++;
      if (idelta<0){
	ncl0++;
	fitterY0.AddPoint(&x, cluster->GetY(),csigmaY);
	fitterZ0.AddPoint(&x, cluster->GetZ(),csigmaZ);
      }
      if (idelta>0){
	ncl1++;
	fitterY1.AddPoint(&x, cluster->GetY(),csigmaY);
	fitterZ1.AddPoint(&x, cluster->GetZ(),csigmaZ);
      }
      fitterY2.AddPoint(&x, cluster->GetY(),csigmaY);
      fitterZ2.AddPoint(&x, cluster->GetZ(),csigmaZ);
    }
    if (nclFound<kDelta*kMinRatio) continue;    
    fitterY2.Eval();
    fitterZ2.Eval();
    Double_t chi2 = (fitterY2.GetChisquare()+fitterZ2.GetChisquare())/(2.*nclFound-6.);
    if (chi2>kCutChi2) continue;
    //
    //
    // REMOVE KINK
    //
    if (ncl0>4){
      fitterY0.Eval();
      fitterZ0.Eval();
    }
    if (ncl1>4){
      fitterY1.Eval();
      fitterZ1.Eval();
    }
    //
    //
    if (ncl0>4&&ncl1>4){
      fitterY0.GetCovarianceMatrix(matrixY0);
      fitterY1.GetCovarianceMatrix(matrixY1);
      fitterZ0.GetCovarianceMatrix(matrixZ0);
      fitterZ1.GetCovarianceMatrix(matrixZ1);
      fitterY1.GetParameters(paramY1);
      fitterZ1.GetParameters(paramZ1);
      fitterY0.GetParameters(paramY0);
      fitterZ0.GetParameters(paramZ0);
      paramY0-= paramY1;
      paramZ0-= paramZ1;
      matrixY0+= matrixY1;
      matrixZ0+= matrixZ1;
      Double_t chi2 =0;
      TMatrixD difY(2,1,paramY0.GetMatrixArray());
      TMatrixD difYT(1,2,paramY0.GetMatrixArray());
      matrixY0.Invert();
      TMatrixD mulY(matrixY0,TMatrixD::kMult,difY);
      TMatrixD chi2Y(difYT,TMatrixD::kMult,mulY);
      chi2+=chi2Y(0,0);
      TMatrixD difZ(2,1,paramZ0.GetMatrixArray());
      TMatrixD difZT(1,2,paramZ0.GetMatrixArray());
      matrixZ0.Invert();
      TMatrixD mulZ(matrixZ0,TMatrixD::kMult,difZ);
      TMatrixD chi2Z(difZT,TMatrixD::kMult,mulZ);
      chi2+= chi2Z(0,0);      
      if (chi2*0.25>kCutChi2) continue;
    } 
    Double_t paramY[4], paramZ[4];
    paramY[0] = fitterY2.GetParameter(0);
    paramY[1] = fitterY2.GetParameter(1);
    paramY[2] = fitterY2.GetParameter(2);
    paramZ[0] = fitterZ2.GetParameter(0);
    paramZ[1] = fitterZ2.GetParameter(1);
    paramZ[2] = fitterZ2.GetParameter(2);    
    //
    //
    //
    Double_t tracky = paramY[0];
    Double_t trackz = paramZ[0];
    Float_t  deltay = tracky-cluster0->GetY();
    Float_t  deltaz = trackz-cluster0->GetZ();
    Float_t  angley = paramY[1]-paramY[0]/xref;
    Float_t  anglez = paramZ[1];
    //
    //
    Float_t max = cluster0->GetMax();
    TProfile *profAmpRow =  (TProfile*)fArrayAmpRow->At(sector);
    if ( irow >= kFirstLargePad) max /= kLargePadSize;
    profAmpRow->Fill(cluster0->GetRow(), max);
    TH1F *hisAmp =  (TH1F*)fArrayAmp->At(sector);
    hisAmp->Fill(max);
    //
    //
    Int_t  ipad=0;
    if (cluster0->GetDetector()>=36) {
      ipad=1;
      if (cluster0->GetRow()>63) ipad=2;
    }
    //
    //
    TH3F * his3=0;
    his3 = (TH3F*)fRMSY->At(ipad);
    if (his3) his3->Fill(250-TMath::Abs(cluster0->GetZ()),TMath::Abs(angley),TMath::Sqrt(cluster0->GetSigmaY2()));
    his3 = (TH3F*)fRMSZ->At(ipad);
    if (his3) his3->Fill(250-TMath::Abs(cluster0->GetZ()),TMath::Abs(anglez),TMath::Sqrt(cluster0->GetSigmaZ2()));
      his3 = (TH3F*)fArrayQRMSY->At(GetBin(cluster0->GetMax(),ipad));
      if (his3) his3->Fill(250-TMath::Abs(cluster0->GetZ()),TMath::Abs(angley), TMath::Sqrt(cluster0->GetSigmaY2()));
      //
      his3 = (TH3F*)fArrayQRMSZ->At(GetBin(cluster0->GetMax(),ipad));
      if (his3) his3->Fill(250-TMath::Abs(cluster0->GetZ()),TMath::Abs(anglez),TMath::Sqrt(cluster0->GetSigmaZ2()));

    //
    // Fill resolution histograms
    //
    Bool_t useForResol= kTRUE;
    if (fitterY2.GetParError(0)>kErrorFraction*csigmaY) useForResol=kFALSE;

    if (useForResol){
      fDeltaY->Fill(deltay);
      fDeltaZ->Fill(deltaz);
      his3 = (TH3F*)fResolY->At(ipad);
      if (his3) his3->Fill(250-TMath::Abs(cluster0->GetZ()),TMath::Abs(angley), deltay);
      his3 = (TH3F*)fResolZ->At(ipad);
      if (his3) his3->Fill(250-TMath::Abs(cluster0->GetZ()),TMath::Abs(anglez), deltaz);
      his3 = (TH3F*)fArrayQDY->At(GetBin(cluster0->GetMax(),ipad));
      if (his3) his3->Fill(250-TMath::Abs(cluster0->GetZ()),TMath::Abs(angley), deltay);
      //
      his3 = (TH3F*)fArrayQDZ->At(GetBin(cluster0->GetMax(),ipad));
      if (his3) his3->Fill(250-TMath::Abs(cluster0->GetZ()),TMath::Abs(anglez), deltaz);
    }
    //
    //
    //
    if (useForResol&&nclFound>2*kMinRatio*kDelta){
      //
      // fill resolution trees
      //
      static TLinearFitter fitY0(3,"pol2");
      static TLinearFitter fitZ0(3,"pol2");
      static TLinearFitter fitY2(5,"hyp4");
      static TLinearFitter fitZ2(5,"hyp4");
      static TLinearFitter fitY2Q(5,"hyp4");
      static TLinearFitter fitZ2Q(5,"hyp4");
      static TLinearFitter fitY2S(5,"hyp4");
      static TLinearFitter fitZ2S(5,"hyp4");
      fitY0.ClearPoints();
      fitZ0.ClearPoints();
      fitY2.ClearPoints();
      fitZ2.ClearPoints();
      fitY2Q.ClearPoints();
      fitZ2Q.ClearPoints();
      fitY2S.ClearPoints();
      fitZ2S.ClearPoints();
      
      for (Int_t idelta=-kDelta; idelta<=kDelta; idelta++){
	if (idelta==0) continue;
	if (idelta+irow<0) continue;
	if (idelta+irow>159) continue;
	AliTPCclusterMI * cluster = track->GetClusterPointer(irow+idelta);
	if (!cluster) continue;
	if (cluster->GetType()<0) continue;
	if (cluster->GetDetector()!=sector) continue;
	Double_t x = cluster->GetX()-xref;
	Double_t sigmaY0 = fClusterParam->GetError0Par(0,ipad,(250.0-TMath::Abs(cluster->GetZ())),TMath::Abs(angley));
	Double_t sigmaZ0 = fClusterParam->GetError0Par(1,ipad,(250.0-TMath::Abs(cluster->GetZ())),TMath::Abs(anglez));
	//
	Double_t sigmaYQ = fClusterParam->GetErrorQPar(0,ipad,(250.0-TMath::Abs(cluster->GetZ())),
						       TMath::Abs(angley), TMath::Abs(cluster->GetMax()));
	Double_t sigmaZQ = fClusterParam->GetErrorQPar(1,ipad,(250.0-TMath::Abs(cluster->GetZ())),
						       TMath::Abs(anglez),TMath::Abs(cluster->GetMax()));
	Double_t sigmaYS = fClusterParam->GetErrorQParScaled(0,ipad,(250.0-TMath::Abs(cluster->GetZ())),
						       TMath::Abs(angley), TMath::Abs(cluster->GetMax()));
	Double_t sigmaZS = fClusterParam->GetErrorQParScaled(1,ipad,(250.0-TMath::Abs(cluster->GetZ())),
						       TMath::Abs(anglez),TMath::Abs(cluster->GetMax()));
	Float_t rmsYFactor = fClusterParam->GetShapeFactor(0,ipad,(250.0-TMath::Abs(cluster->GetZ())),
							   TMath::Abs(anglez), TMath::Abs(cluster->GetMax()),
							   TMath::Sqrt(cluster0->GetSigmaY2()),0);
	Float_t rmsZFactor = fClusterParam->GetShapeFactor(0,ipad,(250.0-TMath::Abs(cluster->GetZ())),
							   TMath::Abs(anglez), TMath::Abs(cluster->GetMax()),
							   TMath::Sqrt(cluster0->GetSigmaZ2()),0);
	sigmaYS  = TMath::Sqrt(sigmaYS*sigmaYS+rmsYFactor*rmsYFactor/12.);
	sigmaZS  = TMath::Sqrt(sigmaZS*sigmaZS+rmsZFactor*rmsZFactor/12.+rmsYFactor*rmsYFactor/24.);
	//
	if (kDelta!=0){
	  fitY0.AddPoint(&x, cluster->GetY(), sigmaY0);
	  fitZ0.AddPoint(&x, cluster->GetZ(), sigmaZ0);
	}
	Double_t xxx[4];
	xxx[0] =  ((idelta+irow)%2==0)? 1:0;
	xxx[1] =  x;
	xxx[2] =  ((idelta+irow)%2==0)? x:0;
	xxx[3] =  x*x;	
	fitY2.AddPoint(xxx, cluster->GetY(), sigmaY0);
	fitY2Q.AddPoint(xxx, cluster->GetY(), sigmaYQ);
	fitY2S.AddPoint(xxx, cluster->GetY(), sigmaYS);
	fitZ2.AddPoint(xxx, cluster->GetZ(), sigmaZ0);
	fitZ2Q.AddPoint(xxx, cluster->GetZ(), sigmaZQ);
	fitZ2S.AddPoint(xxx, cluster->GetZ(), sigmaZS);
	//
      }
      //
      fitY0.Eval();
      fitZ0.Eval();
      fitY2.Eval();
      fitZ2.Eval();
      fitY2Q.Eval();
      fitZ2Q.Eval();
      fitY2S.Eval();
      fitZ2S.Eval();
      Float_t chi2Y0 = fitY0.GetChisquare()/(nclFound-3.);
      Float_t chi2Z0 = fitZ0.GetChisquare()/(nclFound-3.);
      Float_t chi2Y2 = fitY2.GetChisquare()/(nclFound-5.);
      Float_t chi2Z2 = fitZ2.GetChisquare()/(nclFound-5.);
      Float_t chi2Y2Q = fitY2Q.GetChisquare()/(nclFound-5.);
      Float_t chi2Z2Q = fitZ2Q.GetChisquare()/(nclFound-5.);
      Float_t chi2Y2S = fitY2S.GetChisquare()/(nclFound-5.);
      Float_t chi2Z2S = fitZ2S.GetChisquare()/(nclFound-5.);
      //
      static  TVectorD    parY0(3);
      static  TMatrixD    matY0(3,3);
      static  TVectorD    parZ0(3);
      static  TMatrixD    matZ0(3,3);
      fitY0.GetParameters(parY0);
      fitY0.GetCovarianceMatrix(matY0);
      fitZ0.GetParameters(parZ0);
      fitZ0.GetCovarianceMatrix(matZ0);
      //
      static  TVectorD    parY2(5);
      static  TMatrixD    matY2(5,5);
      static  TVectorD    parZ2(5);
      static  TMatrixD    matZ2(5,5);
      fitY2.GetParameters(parY2);
      fitY2.GetCovarianceMatrix(matY2);
      fitZ2.GetParameters(parZ2);
      fitZ2.GetCovarianceMatrix(matZ2);
      //
      static  TVectorD    parY2Q(5);
      static  TMatrixD    matY2Q(5,5);
      static  TVectorD    parZ2Q(5);
      static  TMatrixD    matZ2Q(5,5);
      fitY2Q.GetParameters(parY2Q);
      fitY2Q.GetCovarianceMatrix(matY2Q);
      fitZ2Q.GetParameters(parZ2Q);
      fitZ2Q.GetCovarianceMatrix(matZ2Q);
      static  TVectorD    parY2S(5);
      static  TMatrixD    matY2S(5,5);
      static  TVectorD    parZ2S(5);
      static  TMatrixD    matZ2S(5,5);
      fitY2S.GetParameters(parY2S);
      fitY2S.GetCovarianceMatrix(matY2S);
      fitZ2S.GetParameters(parZ2S);
      fitZ2S.GetCovarianceMatrix(matZ2S);
      Float_t sigmaY0   = TMath::Sqrt(matY0(0,0));
      Float_t sigmaZ0   = TMath::Sqrt(matZ0(0,0));
      Float_t sigmaDY0  = TMath::Sqrt(matY0(1,1));
      Float_t sigmaDZ0  = TMath::Sqrt(matZ0(1,1));
      Float_t sigmaY2   = TMath::Sqrt(matY2(1,1));
      Float_t sigmaZ2   = TMath::Sqrt(matZ2(1,1));
      Float_t sigmaDY2  = TMath::Sqrt(matY2(3,3));
      Float_t sigmaDZ2  = TMath::Sqrt(matZ2(3,3));
      Float_t sigmaY2Q  = TMath::Sqrt(matY2Q(1,1));
      Float_t sigmaZ2Q  = TMath::Sqrt(matZ2Q(1,1));
      Float_t sigmaDY2Q = TMath::Sqrt(matY2Q(3,3));
      Float_t sigmaDZ2Q = TMath::Sqrt(matZ2Q(3,3));
      Float_t sigmaY2S  = TMath::Sqrt(matY2S(1,1));
      Float_t sigmaZ2S  = TMath::Sqrt(matZ2S(1,1));
      Float_t sigmaDY2S = TMath::Sqrt(matY2S(3,3));
      Float_t sigmaDZ2S = TMath::Sqrt(matZ2S(3,3));
      //
      // Error parameters
      //
      Float_t csigmaY0 = fClusterParam->GetError0Par(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),TMath::Abs(angley));
      Float_t csigmaZ0 = fClusterParam->GetError0Par(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),TMath::Abs(anglez));
      //
      Float_t csigmaYQ = fClusterParam->GetErrorQPar(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						     TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
      Float_t csigmaZQ = fClusterParam->GetErrorQPar(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						       TMath::Abs(anglez),TMath::Abs(cluster0->GetMax()));
      Float_t csigmaYS = fClusterParam->GetErrorQParScaled(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						     TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
      Float_t csigmaZS = fClusterParam->GetErrorQParScaled(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						       TMath::Abs(anglez),TMath::Abs(cluster0->GetMax()));
      //
      // RMS parameters
      //
      Float_t meanRMSY = 0;
      Float_t meanRMSZ = 0;
      Int_t   nclRMS=0;
      for (Int_t idelta=-2; idelta<=2; idelta++){
	if (idelta+irow<0) continue;
	if (idelta+irow>159) continue;
	AliTPCclusterMI * cluster = track->GetClusterPointer(irow+idelta);
	if (!cluster) continue;
	meanRMSY += TMath::Sqrt(cluster->GetSigmaY2());
	meanRMSZ += TMath::Sqrt(cluster->GetSigmaZ2());
	nclRMS++;
      }
      meanRMSY /= nclRMS; 
      meanRMSZ /= nclRMS; 

      Float_t rmsY     = TMath::Sqrt(cluster0->GetSigmaY2());  
      Float_t rmsZ     = TMath::Sqrt(cluster0->GetSigmaZ2());
      Float_t rmsYT    = fClusterParam->GetRMSQ(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
      Float_t rmsZT    = fClusterParam->GetRMSQ(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
      Float_t rmsYT0    = fClusterParam->GetRMS0(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						 TMath::Abs(angley));
      Float_t rmsZT0    = fClusterParam->GetRMS0(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						 TMath::Abs(anglez));
      Float_t rmsYSigma = fClusterParam->GetRMSSigma(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						     TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
      Float_t rmsZSigma = fClusterParam->GetRMSSigma(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
						     TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()));
      Float_t rmsYFactor = fClusterParam->GetShapeFactor(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
							 TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()),
							 rmsY,meanRMSY);
      Float_t rmsZFactor = fClusterParam->GetShapeFactor(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
							 TMath::Abs(anglez), TMath::Abs(cluster0->GetMax()),
							 rmsZ,meanRMSZ);
      //
      // cluster debug
      //
      (*fDebugStream)<<"ResolCl"<<	
	"Sector="<<sector<<
	"Cl.="<<cluster0<<
	"CSigmaY0="<<csigmaY0<<   // cluster errorY
	"CSigmaYQ="<<csigmaYQ<<   // cluster errorY - q scaled
	"CSigmaYS="<<csigmaYS<<   // cluster errorY - q scaled
	"CSigmaZ0="<<csigmaZ0<<   // 
	"CSigmaZQ="<<csigmaZQ<<
	"CSigmaZS="<<csigmaZS<<
	"shapeYF="<<rmsYFactor<<
	"shapeZF="<<rmsZFactor<<
	"rmsY="<<rmsY<<
	"rmsZ="<<rmsZ<<
	"rmsYM="<<meanRMSY<<
	"rmsZM="<<meanRMSZ<<
	"rmsYT="<<rmsYT<<
	"rmsZT="<<rmsZT<<
	"rmsYT0="<<rmsYT0<<
	"rmsZT0="<<rmsZT0<<
	"rmsYS="<<rmsYSigma<<  
	"rmsZS="<<rmsZSigma<<
	"IPad="<<ipad<<
	"Ncl="<<nclFound<<	
	"PY0.="<<&parY0<<
	"PZ0.="<<&parZ0<<
	"SigmaY0="<<sigmaY0<< 
	"SigmaZ0="<<sigmaZ0<< 
	"\n";
      //
      // tracklet dubug
      //
      (*fDebugStream)<<"ResolTr"<<	
	"IPad="<<ipad<<
	"Sector="<<sector<<
	"Ncl="<<nclFound<<	
	"chi2Y0="<<chi2Y0<<
	"chi2Z0="<<chi2Z0<<
	"chi2Y2="<<chi2Y2<<
	"chi2Z2="<<chi2Z2<<
	"chi2Y2Q="<<chi2Y2Q<<
	"chi2Z2Q="<<chi2Z2Q<<
	"chi2Y2S="<<chi2Y2S<<
	"chi2Z2S="<<chi2Z2S<<
	"PY0.="<<&parY0<<
	"PZ0.="<<&parZ0<<
	"PY2.="<<&parY2<<
	"PZ2.="<<&parZ2<<
	"PY2Q.="<<&parY2Q<<
	"PZ2Q.="<<&parZ2Q<<
	"PY2S.="<<&parY2S<<
	"PZ2S.="<<&parZ2S<<
	"SigmaY0="<<sigmaY0<< 
	"SigmaZ0="<<sigmaZ0<< 
	"SigmaDY0="<<sigmaDY0<< 
	"SigmaDZ0="<<sigmaDZ0<< 
	"SigmaY2="<<sigmaY2<< 
	"SigmaZ2="<<sigmaZ2<< 
	"SigmaDY2="<<sigmaDY2<< 
	"SigmaDZ2="<<sigmaDZ2<< 
	"SigmaY2Q="<<sigmaY2Q<< 
	"SigmaZ2Q="<<sigmaZ2Q<< 
	"SigmaDY2Q="<<sigmaDY2Q<< 
	"SigmaDZ2Q="<<sigmaDZ2Q<< 
	"SigmaY2S="<<sigmaY2S<< 
	"SigmaZ2S="<<sigmaZ2S<< 
	"SigmaDY2S="<<sigmaDY2S<< 
	"SigmaDZ2S="<<sigmaDZ2S<< 
	"\n";
    }
  }    
}


void  AliTPCcalibTracks::AlignUpDown(AliTPCseed * track, AliESDtrack * esdTrack){
  //
  // Make simple parabolic fit
  //
  const Int_t kMinClusters = 60;
  const Int_t kMinClustersSector =15;
  const Float_t kSigmaCut = 6;
  const Float_t kMaxTan = TMath::Tan(TMath::Pi()*10./180.);
  const Float_t kDeadZone = 6.;
  const Float_t kMinZ     = 15;
  if (track->GetNumberOfClusters()<kMinClusters) return;
  if (TMath::Abs(track->GetZ())<kMinZ) return;
  //
  Int_t nclUp   = 0;
  Int_t nclDown = 0;
  Int_t rSector =-1;
  Float_t refX  = (param.GetInnerRadiusUp()+param.GetOuterRadiusLow())*0.5;
  for (Int_t irow=0; irow<159; irow++){
    AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
    if (!cluster0) continue;
    Int_t sector = cluster0->GetDetector();
    if (rSector<0) rSector=sector%36;
    if (sector%36 != rSector) continue;
    if (  ((TMath::Abs(cluster0->GetY())-kDeadZone)/cluster0->GetX())>kMaxTan) continue;  //remove edge clusters
    if (sector>35) nclUp++;
    if (sector<36) nclDown++;
  }
  if (nclUp<kMinClustersSector) return;
  if (nclDown<kMinClustersSector) return;
  //
  //
  TLinearFitter fitterY(5,"hyp4");  //fitter with common 2 nd derivation
  TLinearFitter fitterZ(5,"hyp4");
  //
  TLinearFitter fitterY0(3,"pol2"); 
  TLinearFitter fitterZ0(3,"pol2");
  TLinearFitter fitterY1(3,"pol2");
  TLinearFitter fitterZ1(3,"pol2");
  //
  Float_t msigmay =1;
  Float_t msigmaz =1;
  Float_t param0[3];
  Float_t param1[3];
  Float_t angley=0;
  Float_t anglez=0;
  //
  for (Int_t iter=0; iter<3; iter++){
    nclUp  = 0;
    nclDown= 0;
    for (Int_t irow=0; irow<159; irow++){
      AliTPCclusterMI * cluster0 = track->GetClusterPointer(irow);
      if (!cluster0) continue;
      Int_t sector = cluster0->GetDetector();
      if (sector%36 != rSector) continue;
      Double_t y = cluster0->GetY();
      Double_t z = cluster0->GetZ();
      //remove edge clusters
      if ( (iter==0) && ((TMath::Abs(cluster0->GetY())-kDeadZone)/cluster0->GetX())>kMaxTan ) continue;  
      if (iter>0){
	Float_t tx = cluster0->GetX()-refX;
	Float_t ty = 0;
	if (sector<36){
	  ty = param0[0]+param0[1]*tx+param0[2]*tx*tx;
	}else{
	  ty = param1[0]+param1[1]*tx+param1[2]*tx*tx;	  
	}
	if (((TMath::Abs(ty)-kDeadZone)/cluster0->GetX())>kMaxTan) continue;
	if (TMath::Abs(ty-y)>kSigmaCut*(msigmay+0.2)) continue;
      }
      Int_t  ipad=0;
      if (cluster0->GetDetector()>=36) {
	ipad=1;
	if (cluster0->GetRow()>63) ipad=2;
      }
      //
      Float_t sigmaY =msigmay;
      Float_t sigmaZ =msigmay;      
      if (iter==2){
	sigmaY = fClusterParam->GetErrorQParScaled(0,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
							   TMath::Abs(angley), TMath::Abs(cluster0->GetMax()));
	sigmaZ = fClusterParam->GetErrorQParScaled(1,ipad,(250.0-TMath::Abs(cluster0->GetZ())),
							   TMath::Abs(anglez),TMath::Abs(cluster0->GetMax()));
      }
      Double_t deltaX = cluster0->GetX()-refX;
      Double_t x[5];
      x[0] = (ipad==0) ? 0:1;
      x[1] = deltaX;
      x[2] = (ipad==0) ? 0:deltaX;
      x[3] = deltaX*deltaX;
      if (ipad<2){
	fitterY.AddPoint(x,y,sigmaY);
	fitterZ.AddPoint(x,z,sigmaZ);
      }
      if (ipad==0){
	nclDown++;
	fitterY0.AddPoint(&deltaX,y,sigmaY);
	fitterZ0.AddPoint(&deltaX,z,sigmaZ);
      }
      if (ipad==1){
	nclUp++;
	fitterY1.AddPoint(&deltaX,y,sigmaY);
	fitterZ1.AddPoint(&deltaX,z,sigmaZ);
      }
    }
    if (nclUp<kMinClustersSector) continue;
    if (nclDown<kMinClustersSector) continue;
    fitterY.Eval();
    fitterZ.Eval();
    fitterY0.Eval();
    fitterZ0.Eval();
    fitterY1.Eval();
    fitterZ1.Eval();
    param0[0] = fitterY0.GetParameter(0);
    param0[1] = fitterY0.GetParameter(1);
    param0[2] = fitterY0.GetParameter(2);
    param1[0] = fitterY1.GetParameter(0);
    param1[1] = fitterY1.GetParameter(1);
    param1[2] = fitterY1.GetParameter(2);
    //
    angley = fitterY.GetParameter(2);
    anglez = fitterZ.GetParameter(2);
    //
    TVectorD    parY(5);
    TMatrixD    matY(5,5);
    TVectorD    parZ(5);
    TMatrixD    matZ(5,5);
    Double_t    chi2Y= fitterY.GetChisquare()/(nclUp+nclDown); 
    Double_t    chi2Z= fitterZ.GetChisquare()/(nclUp+nclDown); 
    fitterY.GetParameters(parY);
    fitterY.GetCovarianceMatrix(matY);
    fitterZ.GetParameters(parZ);
    fitterZ.GetCovarianceMatrix(matZ); 
    if (iter==0) {
      msigmay = msigmay*TMath::Sqrt(chi2Y);
      msigmaz = msigmaz*TMath::Sqrt(chi2Z);
    }
    Float_t sigmaY  = TMath::Sqrt(matY(1,1)*chi2Y);
    Float_t sigmaDY = TMath::Sqrt(matY(3,3)*chi2Y);
    Float_t sigmaDDY = TMath::Sqrt(matY(4,4)*chi2Y);
    Float_t sigmaZ  = TMath::Sqrt(matZ(1,1)*chi2Z);
    Float_t sigmaDZ = TMath::Sqrt(matZ(3,3)*chi2Z);
    Float_t sigmaDDZ = TMath::Sqrt(matZ(4,4)*chi2Z);
    //
    TVectorD    parY0(3);
    TMatrixD    matY0(3,3);
    TVectorD    parZ0(3);
    TMatrixD    matZ0(3,3);
    Double_t    chi2Y0= fitterY0.GetChisquare()/(nclDown); 
    Double_t    chi2Z0= fitterZ0.GetChisquare()/(nclDown); 
    fitterY0.GetParameters(parY0);
    fitterY0.GetCovarianceMatrix(matY0);
    fitterZ0.GetParameters(parZ0);
    fitterZ0.GetCovarianceMatrix(matZ0); 
    Float_t sigmaY0  = TMath::Sqrt(matY0(0,0)*chi2Y0);
    Float_t sigmaDY0 = TMath::Sqrt(matY0(1,1)*chi2Y0);
    Float_t sigmaDDY0 = TMath::Sqrt(matY0(2,2)*chi2Y0);
    Float_t sigmaZ0  = TMath::Sqrt(matZ0(0,0)*chi2Z0);
    Float_t sigmaDZ0 = TMath::Sqrt(matZ0(1,1)*chi2Z0);
    Float_t sigmaDDZ0 = TMath::Sqrt(matZ0(2,2)*chi2Z0);
    //
    TVectorD    parY1(3);
    TMatrixD    matY1(3,3);
    TVectorD    parZ1(3);
    TMatrixD    matZ1(3,3);
    Double_t    chi2Y1= fitterY1.GetChisquare()/(nclUp); 
    Double_t    chi2Z1= fitterZ1.GetChisquare()/(nclUp); 
    fitterY1.GetParameters(parY1);
    fitterY1.GetCovarianceMatrix(matY1);
    fitterZ1.GetParameters(parZ1);
    fitterZ1.GetCovarianceMatrix(matZ1); 
    Float_t sigmaY1  = TMath::Sqrt(matY1(0,0)*chi2Y1);
    Float_t sigmaDY1 = TMath::Sqrt(matY1(1,1)*chi2Y1);
    Float_t sigmaDDY1 = TMath::Sqrt(matY1(2,2)*chi2Y1);
    Float_t sigmaZ1  = TMath::Sqrt(matZ1(0,0)*chi2Z1);
    Float_t sigmaDZ1 = TMath::Sqrt(matZ1(1,1)*chi2Z1);
    Float_t sigmaDDZ1 = TMath::Sqrt(matZ1(2,2)*chi2Z1);
    const AliESDfriendTrack * ftrack = esdTrack->GetFriendTrack();
    AliTrackPointArray *points = (AliTrackPointArray*)ftrack->GetTrackPointArray();

    if (iter>0) (*fDebugStream)<<"Align"<<
      "track.="<<track<<
      "Iter="<<iter<<
      "xref="<<refX<<
      "Points="<<points<<
      "Sector="<<rSector<<
      "nclUp="<<nclUp<<
      "nclDown="<<nclDown<<
      "angley="<<angley<<
      "anglez="<<anglez<<
      //
      "chi2Y="<<chi2Y<<
      "chi2Z="<<chi2Z<<
      "parY.="<<&parY<<
      "parZ.="<<&parZ<<
      "matY.="<<&matY<<
      "matZ.="<<&matZ<<
      "sigmaY="<<sigmaY<<
      "sigmaZ="<<sigmaZ<<
      "sigmaDY="<<sigmaDY<<
      "sigmaDZ="<<sigmaDZ<<
      "sigmaDDY="<<sigmaDDY<<
      "sigmaDDZ="<<sigmaDDZ<<
      //
      "chi2Y0="<<chi2Y0<<
      "chi2Z0="<<chi2Z0<<
      "parY0.="<<&parY0<<
      "parZ0.="<<&parZ0<<
      "matY0.="<<&matY0<<
      "matZ0.="<<&matZ0<<
      "sigmaY0="<<sigmaY0<<
      "sigmaZ0="<<sigmaZ0<<
      "sigmaDY0="<<sigmaDY0<<
      "sigmaDZ0="<<sigmaDZ0<<
      "sigmaDDY0="<<sigmaDDY0<<
      "sigmaDDZ0="<<sigmaDDZ0<<
      //
      "chi2Y1="<<chi2Y1<<
      "chi2Z1="<<chi2Z1<<
      "parY1.="<<&parY1<<
      "parZ1.="<<&parZ1<<
      "matY1.="<<&matY1<<
      "matZ1.="<<&matZ1<<
      "sigmaY1="<<sigmaY1<<
      "sigmaZ1="<<sigmaZ1<<
      "sigmaDY1="<<sigmaDY1<<
      "sigmaDZ1="<<sigmaDZ1<<
      "sigmaDDY1="<<sigmaDDY1<<
      "sigmaDDZ1="<<sigmaDDZ1<<
      "\n";
  }
}
