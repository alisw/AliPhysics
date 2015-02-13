/// \file makeBBFit.C
/// \brief Macro to make the MC based dEdx fits, and study the influence of thre track selection to the PID.
/// 
/// Motivation
/// ----------
/// In the ALICE Geant3 MC the input Bethe-Bloch parameterization of the primary ionization can be 
/// parameterized by user defined formula.
/// 
/// In detector Input \f$ dE/dx(BG)_{in} \f$ (more exact \f$ dN_{prim}/dx \f$) is transformed to the output reconstructed \f$ dE/dx_{rec} \f$. 
/// While original input function is just function of particle \beta\gama, random variable, reconstructed 
/// dEdx estimate, is influenced by detection processes and is sensitive to other aspects namily 
/// diffusion, track inclination angle (\f$\phi\f$,\f$\theta\f$), gain ...
/// 
/// In the following we will calibrate transform function:
/// \f[
/// \frac{dE}{dx}_{BB}= f_{tr}\left(\frac{dE}{dx}_{rec}, \phi, \eta\right)
/// \f]
/// 
///  Example code for author:
/// ~~~{.cpp} 
/// .x $HOME/rootlogon.C
/// .x $HOME/NimStyle.C
/// .L $ALICE_ROOT/TPC/macros/data2011/makeBBFit.C+
/// Init();
/// MakeESDCutsPID();
/// makeBBfit(10000000);
/// ~~~

#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TStatToolkit.h"
#include "TGraphErrors.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TLatex.h"
//
#include "TTreeStream.h"
#include "AliTPCParam.h"
#include "AliTPCcalibBase.h"

//
// Global variables
//

TCut cutGeom, cutNcr, cutNcl;
TCut cutFiducial;
Int_t pdgs[4]={11,211,321,2212};
Double_t massPDG[4]={0};
const char *partName[4]={"Electron","Pion","Kaon","Proton"};
const Int_t colors[5]={kGreen,kBlack,kRed,kBlue,5};
const Int_t markers[5]={24,20,21,25,25};
const Float_t markerSize[5]={1.2,1.2,1.0,1.2,1};
TTree * treeFit=0;

void Init();            // init (ALICE) BB parameterization

//
void FitTransferFunctionTheta(Bool_t isMax, Int_t dType, Int_t tgl, Bool_t skipKaon, Double_t maxP, TTreeSRedirector *pcstream, Bool_t doDraw);
void FitTransferFunctionAll(Bool_t isMax, Int_t dType, Bool_t skipKaon, Double_t maxP, TTreeSRedirector *pcstream, Bool_t doDraw);

void Init(){
  /// 1.) Register default ALICE parametrization

  AliTPCParam param;
  param.RegisterBBParam(param.GetBetheBlochParamAlice(),1);
  for (Int_t ihis=0; ihis<4; ihis++)     massPDG [ihis]= TDatabasePDG::Instance()->GetParticle(pdgs[ihis])->Mass();
}


void MakeESDCutsPID(Double_t fgeom=1, Double_t fcr=0.85, Double_t fcl=0.7 ){
  /// Cuts to be used in the sequenece
  /// 1.)  cutGeom - stronger cut on geometry  - perfect agreement data-MC
  ///              - Default length factor 1
  /// 2.)  cutNcr  - relativally stong cut on the number of crossed rows to gaurantee tracking performance
  ///                relativally good agrement in the MC except of the week decay propability
  ///              - Default length factor 0.85
  /// 3.)  cutNcl  - very week cut on the numebr of clusters
  ///              - Default length factor 0.7
  ///
  /// Combined cuts should be combination:
  /// cutGeom+cutNcr+cutNcl
  ///  Default  factors 1, 0.85 and 0.7 should be suited for most of analysis
  ///  Modification can be expected for the Jet analysis and further tuning of the parameters for the
  ///     relativistic PID analysis

  cutGeom=TString::Format("esdTrack.GetLengthInActiveZone(0,3,236, -5 ,0,0)> %f*(130-5*abs(esdTrack.fP[4]))",fgeom); 
  //  Geomtrical cut in fiducial volume - proper description  in the MC
  //  GetLengthInActiveZone(Int_t mode, Double_t deltaY, Double_t deltaZ, Double_t bz, Double_t exbPhi = 0, TTreeSRedirector* pcstream = 0) const
  //  mode   ==0     - inner param used
  //  deltaY ==3     - cut on edges of sector 
  //                   a.) to avoid dead zone - bias for tracking
  //                   b.) zone with lower Q qualit
  //                   c.) non homogenous Q sample - preferable IROC is skipped -bais in the dEdx
  //  deltaZ ==236   - 
  //  bz = 5 KG      - proper sign should be used ohterwise wrong calculation
  cutNcr=TString::Format("esdTrack.GetTPCClusterInfo(3,1,0,159,1)>%f*(130-5*abs(esdTrack.fP[4]))",fcr);
  //
  // Cut on the number of crossed raws
  // GetTPCClusterInfo(Int_t nNeighbours = 3, Int_t type = 0, Int_t row0 = 0, Int_t row1 = 159, Int_t bitType = 0) 
  // pad row is decalared as crossed if some clusters was found in small neighborhood +- nNeighbours
  //    nNeighbours =3  - +-3 padrows used to define the crossed rows
  //    type = 0        - return number of rows (type==1 returns fraction of clusters)
  //    row0-row1       - upper and lower part of integration
  //  
  cutNcl=TString::Format("esdTrack.fTPCncls>%f*(130-5*abs(esdTrack.fP[4]))",fcl);
  //
  // cut un the number of clusters
  //
  cutFiducial="abs(esdTrack.fP[3])<1&&abs(esdTrack.fD)<1";
}

void makeBBfit(Int_t ntracks=100000){
  /// Make dEdx  fits of indiviadual particle species in bins:
  ///  a.)  momenta
  ///  b.)  tan(theta) -tgl
  ///  c.)  per detector segmen
  ///  d.)  Qmax/Qtot

  TFile * ff = TFile::Open("Filtered.root");
  TTreeSRedirector *pcstream = new TTreeSRedirector("dedxFit.root","update");
  TTree * treeMC = (TTree*)ff->Get("highPt");
  //  treeMC->SetCacheSize(6000000000);
  //
  // 1.) Register default ALICE parametrization
  //
  AliTPCParam param;
  param.RegisterBBParam(param.GetBetheBlochParamAlice(),1);
  TH3D * hisdEdx3D[4]={0};
  //
  for (Int_t ihis=0; ihis<4; ihis++) {
    massPDG [ihis]= TDatabasePDG::Instance()->GetParticle(pdgs[ihis])->Mass();
    treeMC->SetAlias(TString::Format("cut%s",partName[ihis]), TString::Format("abs(particle.fPdgCode)==%d",pdgs[ihis]));
    treeMC->SetAlias(TString::Format("dEdxp%s",partName[ihis]), TString::Format("AliTPCParam::BetheBlochAleph(esdTrack.fIp.P()/%f,1)",massPDG[ihis]));
  }
  //
  // 2.) Fill dEdx histograms if not done before
  //

  if (pcstream->GetFile()->Get("RatioP_QMax0Pion3D")==0){    
    //
    TStopwatch timer;
    for (Int_t iDetType=0; iDetType<9; iDetType++){
      for (Int_t ihis=0; ihis<4; ihis++){
	printf("%d\t%d\n",iDetType,ihis);
	timer.Print();
	TString detType="All";
	TString dedx   ="esdTrack.fTPCsignal";
	if (iDetType>0){
	  detType=TString::Format("Q%s%d",(iDetType-1)/4>0?"Max":"Tot",(iDetType-1)%4);
	  if (iDetType<5) dedx=TString::Format("esdTrack.fTPCdEdxInfo.GetSignalTot(%d)",(iDetType-1)%4);
	  if (iDetType>5) dedx=TString::Format("esdTrack.fTPCdEdxInfo.GetSignalMax(%d)",(iDetType-1)%4);
	}  

	TCut cutPDG=TString::Format("abs(particle.fPdgCode)==%d",pdgs[ihis]).Data();
	TString hname= TString::Format("RatioP_%s%s3D",detType.Data(),partName[ihis]);	
	TString query= TString::Format("%s/AliTPCParam::BetheBlochAleph(esdTrack.fIp.P()/%f,1):abs(esdTrack.fP[3]):esdTrack.fIp.P()>>%s",dedx.Data(), massPDG[ihis],hname.Data());
	hisdEdx3D[ihis]  = new TH3D(hname.Data(),hname.Data(), 50, 0.2,25, 10,0,1, 100,20,80);
	AliTPCcalibBase::BinLogX(hisdEdx3D[ihis]->GetXaxis());
	treeMC->Draw(query,cutFiducial+cutGeom+cutNcr+cutNcl+cutPDG,"goff",ntracks);
	
	hisdEdx3D[ihis]->GetXaxis()->SetTitle("p (GeV/c)");
	hisdEdx3D[ihis]->GetYaxis()->SetTitle("dEdx/dEdx_{BB} (a.u.)");
	pcstream->GetFile()->cd();
	hisdEdx3D[ihis]->Write(hname.Data());
      }
    }
  }
  delete pcstream;
  //
  // Fit histograms
  //
  pcstream = new TTreeSRedirector("dedxFit.root","update");
  TF1 fg("fg","gaus");
  //
  for (Int_t iDetType=0; iDetType<9; iDetType++){
    TString detType="All";
    if (iDetType>0){
      detType=TString::Format("Q%s%d",(iDetType-1)/4>0?"Max":"Tot",(iDetType-1)%4);
    }    
    for (Int_t ihis=0; ihis<4; ihis++){ 
      TString hname= TString::Format("RatioP_%s%s3D",detType.Data(),partName[ihis]);	
      hisdEdx3D[ihis] = (TH3D*)pcstream->GetFile()->Get(hname.Data()); 
      Int_t nbinsP =  hisdEdx3D[0]->GetXaxis()->GetNbins();
      Int_t nbinsTgl =  hisdEdx3D[0]->GetYaxis()->GetNbins();
      //
      for (Int_t ibinP=2; ibinP<nbinsP; ibinP++){
	for (Int_t ibinTgl=2; ibinTgl<nbinsTgl; ibinTgl++){
	  //
	  Double_t pCenter =  hisdEdx3D[0]->GetXaxis()->GetBinCenter(ibinP);
	  Double_t tglCenter =  hisdEdx3D[0]->GetYaxis()->GetBinCenter(ibinTgl);
	  TH1D * hisProj =hisdEdx3D[ihis]->ProjectionZ("xxx", ibinP-1,ibinP+1, ibinTgl-1,ibinTgl+1);
	  Double_t entries = hisProj->GetEntries();
	  if (entries<10) continue;
	  Double_t mean = hisProj->GetMean();
	  Double_t rms = hisProj->GetRMS();
	  hisProj->Fit(&fg,"","");
	  TVectorD vecFit(3, fg.GetParameters());
	  TVectorD vecFitErr(3, fg.GetParErrors());
	  Double_t chi2=fg.GetChisquare();
	  Double_t mass    = massPDG[ihis];
	  Double_t dEdxExp =  AliTPCParam::BetheBlochAleph(pCenter/mass,1);
	  Bool_t isAll=(iDetType==0);
	  Bool_t isMax=(iDetType-1)/4>0;
	  Bool_t isTot=(iDetType-1)/4==0;
	  Int_t dType=(iDetType-1)%4;
	  Double_t dEdx = mean*dEdxExp;
	  //if (
	  (*pcstream)<<"fitdEdxG"<<
	    // dEdx type
	    "iDet="<<iDetType<<      // detector internal number
	    "isAll="<<isAll<<        // full TPC used?
	    "isMax="<<isMax<<        // Qmax charge used ?
	    "isTot="<<isTot<<        // Qtot charge used?
	    "dType="<<dType<<        // Detector region 0-IROC, 1- OROCmedium, 2- OROClong, 3- OROCboth
	    "pType="<<ihis<<         // particle type
	    //
	    "entries="<<entries<<      // entries in histogram
	    "ibinTgl="<<ibinTgl<<      //
	    "ibinP="<<ibinP<<          // 
	    "pCenter="<<pCenter<<      // momentum of center bin
	    "tglCenter="<<tglCenter<<  // tangent lambda
	    //
	    "mass="<<mass<<             // particle mass
	    "dEdxExp="<<dEdxExp<<       // mean expected dEdx in bin 
	    "dEdx="<<dEdx<<             // mean measured dEdx in bin
	    "mean="<<mean<<             // mean measured/expected
	    "rms="<<rms<<               // 
	    "chi2="<<chi2<<             // chi2 of the gausian fit
	    "vecFit.="<<&vecFit<<       // gaus fit param
	    "vecFitErr.="<<&vecFitErr<< // gaus fit error
	    "\n";
	}
      }
    }
  }
  delete pcstream;
  //
  //
  //
}


void FitTransferFunctionScanAll(){
  /// Make fit of transfer functions - parabolic in theta

  TTreeSRedirector * pcstream = new TTreeSRedirector("dedxFit.root","update");
  treeFit=(TTree*)pcstream->GetFile()->Get("fitdEdxG"); 
  treeFit = (TTree*)pcstream->GetFile()->Get("fitdEdxG");
  treeFit->SetMarkerStyle(25);
  treeFit->SetCacheSize(200000000);
  for (Int_t isMax=0; isMax<=1; isMax++){
    for (Int_t dType=0; dType<4; dType++){      
      for (Int_t skipKaon=0; skipKaon<=1; skipKaon++){
	for (Float_t maxP=3; maxP<21; maxP+=3){
	  printf("%d\t%d\t%d\t%f\n",isMax,dType,skipKaon,maxP);
	  FitTransferFunctionAll(isMax,dType,skipKaon,maxP,pcstream,1);
	}
      }
    }
  }
  delete pcstream;
}


void FitTransferFunctionScanTheta(){
  /// Make fit of transfer functions - bin by bin in Theta

  TTreeSRedirector * pcstream = new TTreeSRedirector("dedxFit.root","update");
  treeFit=(TTree*)pcstream->GetFile()->Get("fitdEdxG"); 
  treeFit = (TTree*)pcstream->GetFile()->Get("fitdEdxG");
  treeFit->SetMarkerStyle(25);
  treeFit->SetCacheSize(200000000);
  for (Int_t isMax=0; isMax<=1; isMax++){
    for (Int_t dType=0; dType<4; dType++){
      
      for (Int_t tgl=2; tgl<10; tgl++){
	for (Int_t skipKaon=0; skipKaon<=1; skipKaon++){
	  for (Float_t maxP=3; maxP<21; maxP+=3){
	    printf("%d\t%d\t%d\t%d\t%f\n",isMax,dType,tgl,skipKaon,maxP);
	    FitTransferFunctionTheta(isMax,dType,tgl,skipKaon,maxP,pcstream,1);
	  }
	}
      }
    }
  }
  delete pcstream;
}


void FitTransferFunctionTheta(Bool_t isMax, Int_t dType, Int_t tgl, Bool_t skipKaon, Double_t maxP, TTreeSRedirector *pcstream, Bool_t doDraw){
  /// dEdx_{BB}= f_{tr}(dEdx_{rec}, #phi,#eta)
  ///
  /// Fit the parematers of transfer function
  /// 4 models of transfer function used:
  ///
  /// Example usage:
  ///    FitTransferFunctionTheta(1,3,8,1,5,0)
  ///    isMax=1; dType=1; tgl=5; skipKaon=0; Double_t maxP=10

  if (!pcstream)  pcstream = new TTreeSRedirector("dedxFit.root","update");
  //
  if (!treeFit){
    treeFit = (TTree*)pcstream->GetFile()->Get("fitdEdxG");
    treeFit->SetMarkerStyle(25);
    treeFit->SetCacheSize(100000000);
  }
  const char *chfitType[4] = { "dEdx_{BB}=k(#theta)dEdx_{rec}", 
			       "dEdx_{BB}=k(#theta)dEdx_{rec}+#Delta(#theta)", 
			       "dEdx_{BB}=k(#theta,#phi)dEdx_{rec}+#Delta(#theta,#phi)",
			       "dEdx_{BB}=k(#theta,#phi,1/Q)dEdx_{rec}+#Delta(#theta,#phi,1/Q)",
  };
  //
  treeFit->SetAlias("isOK","entries>100&&chi2<80");
  treeFit->SetAlias("dEdxM","(vecFit.fElements[1]*dEdxExp)/50.");
  treeFit->SetAlias("dEdxMErr","(vecFitErr.fElements[1]*dEdxExp)/50.");
  //
  //
  Int_t  npointsMax=30000000;
  TStatToolkit toolkit;
  Double_t chi20=0,chi21=0,chi22=0,chi23=0;
  Int_t    npoints=0;
  TVectorD param0,param1,param2,param3;
  TMatrixD covar0,covar1,covar2,covar3;
  TString fstringFast0="";
  fstringFast0+="dEdxM++";
  //
  TString fstringFast="";
  fstringFast+="dEdxM++";                  // 1
  fstringFast+="(1/pCenter)++";            // 2
  fstringFast+="(1/pCenter)^2++";          // 3
  fstringFast+="(1/pCenter)*dEdxM++";      // 4
  fstringFast+="((1/pCenter)^2)*dEdxM++";  // 5
  //
  //
  TCut cutUse=TString::Format("isOK&&isMax==%d&&dType==%d&&ibinTgl==%d",isMax,dType,tgl).Data();
  TCut cutFit=cutUse+TString::Format("pCenter<%f",maxP).Data();
  if (skipKaon) cutFit+="pType!=3";
  TCut cutDraw =  "(ibinP%2)==0";
  TString *strDeltaFit0 = TStatToolkit::FitPlane(treeFit,"dEdxExp:dEdxMErr", fstringFast0.Data(),cutFit, chi20,npoints,param0,covar0,-1,0, npointsMax, 1);
  TString *strDeltaFit1 = TStatToolkit::FitPlane(treeFit,"dEdxExp:dEdxMErr", fstringFast0.Data(),cutFit, chi21,npoints,param1,covar1,-1,0, npointsMax, 0);
  TString *strDeltaFit2 = TStatToolkit::FitPlane(treeFit,"dEdxExp:dEdxMErr", fstringFast.Data(),cutFit, chi22,npoints,param2,covar2,-1,0, npointsMax, 0);
  //
  fstringFast+="(1/pCenter)/dEdxM++";       //6
  fstringFast+="((1/pCenter)^2)/dEdxM++";   //7
  TString *strDeltaFit3 = TStatToolkit::FitPlane(treeFit,"dEdxExp:dEdxMErr", fstringFast.Data(),cutFit, chi23,npoints,param3,covar3,-1,0, npointsMax, 0);

  treeFit->SetAlias("fitdEdx0",strDeltaFit0->Data());
  treeFit->SetAlias("fitdEdx1",strDeltaFit1->Data());
  treeFit->SetAlias("fitdEdx2",strDeltaFit2->Data());
  treeFit->SetAlias("fitdEdx3",strDeltaFit3->Data());
  
  strDeltaFit0->Tokenize("++")->Print();
  strDeltaFit1->Tokenize("++")->Print();
  strDeltaFit2->Tokenize("++")->Print();
  strDeltaFit3->Tokenize("++")->Print();
  //
  (*pcstream)<<"fitTheta"<<
    "isMax="<<isMax<<          // switch is Qmax/Qtot used
    "dType="<<dType<<          // detector Type
    "tgl="<<tgl<<              // tgl number 
    "skipKaon="<<skipKaon<<    // Was kaon dEdx used in the calibration?
    "maxP="<<maxP<<            // Maximal p used in the fit 
    "npoints="<<npoints<<      // number of points for fit
    // model 0
    "chi20="<<chi20<<          // chi2
    "param0.="<<&param0<<      // parameters
    "covar0.="<<&covar0<<      // covariance
    // model 1
    "chi21="<<chi21<<          // chi2
    "param1.="<<&param1<<      // parameters
    "covar1.="<<&covar1<<      // covariance
    // model 2
    "chi22="<<chi22<<          // chi2
    "param2.="<<&param2<<      // parameters
    "covar2.="<<&covar2<<      // covariance
    // model 3
    "chi23="<<chi23<<          // chi2
    "param3.="<<&param3<<      // parameters
    "covar3.="<<&covar3<<      // covariance
    "\n";
  //
  if (!doDraw) return;
  TGraphErrors * graphs[4] ={0};
  //  TCanvas *canvasdEdxFit = new TCanvas(TString::Format("canvasdEdxFit%d",fitType),"canvasdEdxFit",800,800);
  TCanvas *canvasdEdxFit = new TCanvas("canvasdEdxFit","canvasdEdxFit",800,800);
  canvasdEdxFit->SetLeftMargin(0.15);
  canvasdEdxFit->SetRightMargin(0.1);
  canvasdEdxFit->SetBottomMargin(0.15);
  canvasdEdxFit->Divide(2,2,0,0);
 
  for (Int_t fitType=0; fitType<4; fitType++){
    canvasdEdxFit->cd(fitType+1);
    TLegend * legendPart = new TLegend(0.16+0.1*(fitType%2==0),0.16-0.14*(fitType<2),0.89,0.4,TString::Format("%s",chfitType[fitType]));
    legendPart->SetTextSize(0.05);
    legendPart->SetBorderSize(0);
    for (Int_t ihis=3; ihis>=0;ihis--){
      gPad->SetLogx(1);
      TString expr= TString::Format("100*(dEdxExp-fitdEdx%d)/dEdxExp:pCenter:100*dEdxMErr",fitType);
      graphs[ihis]= TStatToolkit::MakeGraphErrors(treeFit,expr,cutUse+cutDraw+TString::Format("pType==%d",ihis).Data(), markers[ihis],colors[ihis],markerSize[ihis]);
      graphs[ihis]->SetMinimum(-5);
      graphs[ihis]->SetMaximum(5);
      graphs[ihis]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
      graphs[ihis]->GetXaxis()->SetTitleSize(0.06);
      graphs[ihis]->GetYaxis()->SetTitleSize(0.06);      
      graphs[ihis]->GetYaxis()->SetTitle("#Delta(d#it{E}dx)/d#it{E}dx_{BB}(#beta#gamma) (%)");
      if (ihis==3) graphs[ihis]->Draw("alp");
      graphs[ihis]->Draw("lp");
      legendPart->AddEntry(graphs[ihis],partName[ihis],"p");
    }
    legendPart->Draw();
  }
  TLatex  latexDraw;
  latexDraw.SetTextSize(0.045);

  const char *chDType[4]={"IROC","OROC medium","OROC long","OROC"};
  {
    if (isMax) latexDraw.DrawLatex(1,4,"Q_{Max}");
    if (!isMax) latexDraw.DrawLatex(1,4,"Q_{tot}");
    latexDraw.DrawLatex(1,3.5,chDType[dType%4]);
    latexDraw.DrawLatex(1,3,TString::Format("#Theta=%1.1f",tgl/10.));
    latexDraw.DrawLatex(1,2.5,TString::Format("Fit: p_{t}<%1.1f (GeV/c)",maxP));
    latexDraw.DrawLatex(1,2.,TString::Format("Fit: Skip Kaon=%d",skipKaon));
  }
  //
  //
  canvasdEdxFit->SaveAs(TString::Format("fitTransfer_Max%d_Det_%dTheta%d_SkipKaon%d_MaxP%1.0f.pdf",isMax,dType,tgl, skipKaon,maxP));
  canvasdEdxFit->SaveAs(TString::Format("fitTransfer_Max%d_Det_%dTheta%d_SkipKaon%d_MaxP%1.0f.png",isMax,dType,tgl, skipKaon,maxP));

}
//
//
//
void FitTransferFunctionAll(Bool_t isMax, Int_t dType, Bool_t skipKaon, Double_t maxP, TTreeSRedirector *pcstream, Bool_t doDraw){
  /// dEdx_{BB}= f_{tr}(dEdx_{rec}, #phi,#eta)
  ///
  /// Example usage:
  ///    FitTransferFunctionAll(1,1,1,15,0,1)
  ///    isMax=1; dType=1; skipKaon=0; Double_t maxP=10

  if (!pcstream)  pcstream = new TTreeSRedirector("dedxFit.root","update");
  //
  if (!treeFit){
    treeFit = (TTree*)pcstream->GetFile()->Get("fitdEdxG");
    treeFit->SetMarkerStyle(25);
    treeFit->SetCacheSize(100000000);
  }
  const char *chfitType[4] = { "dEdx_{BB}=k(#theta)dEdx_{rec}", 
			       "dEdx_{BB}=k(#theta)dEdx_{rec}+#Delta(#theta)", 
			       "dEdx_{BB}=k(#theta,#phi)dEdx_{rec}+#Delta(#theta,#phi)",
			       "dEdx_{BB}=k(#theta,#phi,1/Q)dEdx_{rec}+#Delta(#theta,#phi,1/Q)",
  };
  //
  treeFit->SetAlias("isOK","entries>100&&chi2<80");
  treeFit->SetAlias("dEdxM","(vecFit.fElements[1]*dEdxExp)/50.");
  treeFit->SetAlias("dEdxMErr","(vecFitErr.fElements[1]*dEdxExp)/50.");
  //
  //
  Int_t  npointsMax=30000000;
  TStatToolkit toolkit;
  Double_t chi20=0,chi21=0,chi22=0,chi23=0;
  Int_t    npoints=0;
  TVectorD param0,param1,param2,param3;
  TMatrixD covar0,covar1,covar2,covar3;
  //
  //
  //
  TString fstringFast0="";
  fstringFast0+="dEdxM++";
  fstringFast0+="(tglCenter-0.5)++";
  fstringFast0+="(tglCenter-0.5)*dEdxM++";
  //
  TString fstringFast1="";
  fstringFast1+="dEdxM++";
  fstringFast1+="(tglCenter-0.5)++";
  fstringFast1+="(tglCenter-0.5)*dEdxM++";
  fstringFast1+="(tglCenter-0.5)^2++";
  fstringFast1+="((tglCenter-0.5)^2)*dEdxM++";
  //
  TString fstringFast2="";
  fstringFast2+="dEdxM++";                  // 1
  fstringFast2+="(1/pCenter)++";            // 2
  fstringFast2+="(1/pCenter)^2++";          // 3
  fstringFast2+="(1/pCenter)*dEdxM++";      // 4
  fstringFast2+="((1/pCenter)^2)*dEdxM++";  // 5
  //
  fstringFast2+="(tglCenter-0.5)*dEdxM++";                  // 8
  fstringFast2+="(tglCenter-0.5)*(1/pCenter)++";            // 9
  fstringFast2+="(tglCenter-0.5)*(1/pCenter)^2++";          // 10
  fstringFast2+="(tglCenter-0.5)*(1/pCenter)*dEdxM++";      // 11
  fstringFast2+="(tglCenter-0.5)*((1/pCenter)^2)*dEdxM++";  // 12
  //
  //
  fstringFast2+="((tglCenter-0.5)^2)*dEdxM++";                  // 15
  fstringFast2+="((tglCenter-0.5)^2)*(1/pCenter)++";            // 16
  fstringFast2+="((tglCenter-0.5)^2)*(1/pCenter)^2++";          // 17
  fstringFast2+="((tglCenter-0.5)^2)*(1/pCenter)*dEdxM++";      // 18
  fstringFast2+="((tglCenter-0.5)^2)*((1/pCenter)^2)*dEdxM++";  // 19
  TString fstringFast3=fstringFast2;
  //
  fstringFast3+="(1/pCenter)/dEdxM++";       // 6
  fstringFast3+="((1/pCenter)^2)/dEdxM++";   // 7
  fstringFast3+="(tglCenter-0.5)*(1/pCenter)/dEdxM++";       // 13
  fstringFast3+="(tglCenter-0.5)*((1/pCenter)^2)/dEdxM++";   // 14
  fstringFast3+="((tglCenter-0.5)^2)*(1/pCenter)/dEdxM++";       // 20
  fstringFast3+="((tglCenter-0.5)^2)*((1/pCenter)^2)/dEdxM++";   // 21
  //
  //
  //
  TCut cutUse=TString::Format("isOK&&isMax==%d&&dType==%d",isMax,dType).Data();
  TCut cutFit=cutUse+TString::Format("pCenter<%f",maxP).Data();
  if (skipKaon) cutFit+="pType!=3";
  TCut cutDraw =  "(ibinP%2)==0";
  TString *strDeltaFit0 = TStatToolkit::FitPlane(treeFit,"dEdxExp:dEdxMErr", fstringFast0.Data(),cutFit, chi20,npoints,param0,covar0,-1,0, npointsMax, 1);
  TString *strDeltaFit1 = TStatToolkit::FitPlane(treeFit,"dEdxExp:dEdxMErr", fstringFast1.Data(),cutFit, chi21,npoints,param1,covar1,-1,0, npointsMax, 0);
  TString *strDeltaFit2 = TStatToolkit::FitPlane(treeFit,"dEdxExp:dEdxMErr", fstringFast2.Data(),cutFit, chi22,npoints,param2,covar2,-1,0, npointsMax, 0);
  //
  TString *strDeltaFit3 = TStatToolkit::FitPlane(treeFit,"dEdxExp:dEdxMErr", fstringFast3.Data(),cutFit, chi23,npoints,param3,covar3,-1,0, npointsMax, 0);

  treeFit->SetAlias("fitdEdx0",strDeltaFit0->Data());
  treeFit->SetAlias("fitdEdx1",strDeltaFit1->Data());
  treeFit->SetAlias("fitdEdx2",strDeltaFit2->Data());
  treeFit->SetAlias("fitdEdx3",strDeltaFit3->Data());
  
  strDeltaFit0->Tokenize("++")->Print();
  strDeltaFit1->Tokenize("++")->Print();
  strDeltaFit2->Tokenize("++")->Print();
  strDeltaFit3->Tokenize("++")->Print();
  //
  (*pcstream)<<"fitAll"<<
    "isMax="<<isMax<<          // switch is Qmax/Qtot used
    "dType="<<dType<<          // detector Type
    "skipKaon="<<skipKaon<<    // Was kaon dEdx used in the calibration?
    "maxP="<<maxP<<            // Maximal p used in the fit 
    "npoints="<<npoints<<      // number of points for fit
    // model 0
    "chi20="<<chi20<<          // chi2
    "param0.="<<&param0<<      // parameters
    "covar0.="<<&covar0<<      // covariance
    // model 1
    "chi21="<<chi21<<          // chi2
    "param1.="<<&param1<<      // parameters
    "covar1.="<<&covar1<<      // covariance
    // model 2
    "chi22="<<chi22<<          // chi2
    "param2.="<<&param2<<      // parameters
    "covar2.="<<&covar2<<      // covariance
    // model 3
    "chi23="<<chi23<<          // chi2
    "param3.="<<&param3<<      // parameters
    "covar3.="<<&covar3<<      // covariance
    "\n";
  //
  if (!doDraw) return;
  TGraphErrors * graphs[4] ={0};
  //  TCanvas *canvasdEdxFit = new TCanvas(TString::Format("canvasdEdxFit%d",fitType),"canvasdEdxFit",800,800);
  TCanvas *canvasdEdxFit = new TCanvas("canvasdEdxFit","canvasdEdxFit",800,800);
  canvasdEdxFit->SetLeftMargin(0.15);
  canvasdEdxFit->SetRightMargin(0.1);
  canvasdEdxFit->SetBottomMargin(0.15);
  canvasdEdxFit->Divide(2,2,0,0);
 
  for (Int_t fitType=0; fitType<4; fitType++){
    canvasdEdxFit->cd(fitType+1);
    TLegend * legendPart = new TLegend(0.16+0.1*(fitType%2==0),0.16-0.14*(fitType<2),0.89,0.4,TString::Format("%s",chfitType[fitType]));
    legendPart->SetTextSize(0.05);
    legendPart->SetBorderSize(0);
    for (Int_t ihis=3; ihis>=0;ihis--){
      gPad->SetLogx(1);
      TString expr= TString::Format("100*(dEdxExp-fitdEdx%d)/dEdxExp:pCenter:100*dEdxMErr",fitType);
      graphs[ihis]= TStatToolkit::MakeGraphErrors(treeFit,expr,cutUse+cutDraw+TString::Format("pType==%d",ihis).Data(), markers[ihis],colors[ihis],markerSize[ihis]);
      graphs[ihis]->SetMinimum(-5);
      graphs[ihis]->SetMaximum(5);
      graphs[ihis]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
      graphs[ihis]->GetXaxis()->SetTitleSize(0.06);
      graphs[ihis]->GetYaxis()->SetTitleSize(0.06);      
      graphs[ihis]->GetYaxis()->SetTitle("#Delta(d#it{E}dx)/d#it{E}dx_{BB}(#beta#gamma) (%)");
      if (ihis==3) graphs[ihis]->Draw("alp");
      graphs[ihis]->Draw("lp");
      legendPart->AddEntry(graphs[ihis],partName[ihis],"p");
    }
    legendPart->Draw();
  }
  TLatex  latexDraw;
  latexDraw.SetTextSize(0.045);

  const char *chDType[4]={"IROC","OROC medium","OROC long","OROC"};
  {
    if (isMax) latexDraw.DrawLatex(1,4,"Q_{Max}");
    if (!isMax) latexDraw.DrawLatex(1,4,"Q_{tot}");
    latexDraw.DrawLatex(1,3.5,chDType[dType%4]);
    latexDraw.DrawLatex(1,2.5,TString::Format("Fit: p_{t}<%1.1f (GeV/c)",maxP));
    latexDraw.DrawLatex(1,2.,TString::Format("Fit: Skip Kaon=%d",skipKaon));
  }
  //
  //
  canvasdEdxFit->SaveAs(TString::Format("fitTransfer_Max%d_Det_%dSkipKaon%d_MaxP%1.0f.pdf",isMax,dType, skipKaon,maxP));
  canvasdEdxFit->SaveAs(TString::Format("fitTransfer_Max%d_Det_%dSkipKaon%d_MaxP%1.0f.png",isMax,dType, skipKaon,maxP));

}


void DrawFit(){
  ///

  TTreeSRedirector * pcstream = new TTreeSRedirector("dedxFit.root","update");
  TTree * treeTheta=(TTree*)pcstream->GetFile()->Get("fitTheta"); 
  treeTheta->SetMarkerStyle(25);
  TTree * treeAll=(TTree*)pcstream->GetFile()->Get("fitAll"); 
  treeAll->SetMarkerStyle(25);
  //
  //
  //
  //treeFit->Draw("vecFit.fElements[2]/vecFit.fElements[1]*pow(dEdxExp*sqrt(1+tglCenter**2),0.25):dEdxExp:tglCenter","isTot&&dType==1&&entries>500","colz",1000000);

}
