#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include "TROOT.h"
#include "TSelector.h"
#include "TGeoManager.h"
#include "TArrayI.h"
#include "TText.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "AliTRDgeometry.h"
#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliAlignObjParams.h"
#include "AliAlignmentTracks.h"
#include "AliTrackFitterStraight.h"
#include "AliTrackFitterRieman.h"
#include "AliTrackFitterKalman.h"
#include "AliTrackResidualsLinear.h"
#include "AliTrackResidualsFast.h"
#include "AliTrackResidualsChi2.h"
#include "AliTRDalignment.h"
#endif

//#include "/misc/misko/root/peakfit.C"

AliAlignmentTracks    *alt;
AliTrackFitter        *trf;
AliTrackResiduals     *res;
TNtuple               *ntresi;      // residuals after minimization
TNtuple               *ntalpa;      // alignment params
TNtuple               *nteva;       // results of minimization
double aml[100][40][6];             // results from align_module_layers [iter][layer][phi,z,r,rphi,rz,rr]
enum {kStraight, kRieman, kKalman}; // track fitters
enum {kLinear, kFast, kChi2};       // minimizer

//=============================================================================
TNtuple *makeNtuple(AliTrackResiduals *res, char *name, int flag) {

  // make and fill ntuple containing residuals so they can be visualized
  // flag=1 means including transformation by AliTrackResiduals::fAlignObj
  // use flag 0 to see the status before minimization
  // use flag 1 to see the status after minimization

  TNtuple *nt = new TNtuple(name,name,"x0:y0:z0:x1:y1:z1:v:dy",320000);
  AliTrackPoint p0,p1; // 0 - reference; 1 - actual point
  
  for (Int_t itrack = 0; itrack < res->GetNFilledTracks(); itrack++) {
    AliTrackPointArray *volarray;
    AliTrackPointArray *traarray;
    if (!res->GetTrackPointArrays(itrack, volarray, traarray)) continue;
    for (Int_t ipoint = 0; ipoint < volarray->GetNPoints(); ipoint++) {
      traarray->GetPoint(p0,ipoint);
      volarray->GetPoint(p1,ipoint);
      if (flag) res->GetAlignObj()->Transform(p1);
      AliTrackPoint q0,q1; // in tracking system
      Double_t alpha = p1.GetAngle();
      q0=p0.Rotate(alpha);
      q1=p1.Rotate(alpha);
      UShort_t volid = p1.GetVolumeID();
      double x0 = q0.GetX();
      double y0 = q0.GetY();
      double z0 = q0.GetZ();
      double x1 = q1.GetX();
      double y1 = q1.GetY();
      double z1 = q1.GetZ();
      double dy = y1-y0; 
      AliGeomManager::ELayerID layer = AliGeomManager::VolUIDToLayer(volid);
      dy += 0.035*(z1-z0)*pow(-1.0,layer); // pad tilt
      nt->Fill(x0,y0,z0,x1,y1,z1,volid,dy);
    }
  }
  return nt;
}
//=============================================================================
void filterTracks(AliTrackResiduals *res, double level) {

  // find the peak in y-residuals; remove outlier tracks

  // first, find the peak and prepare the track quality cut 

  char *hinam = "filter_tracks_tmp";
  TH1D *hist = new TH1D(hinam,hinam,250,-5,5);
  TNtuple *nt = makeNtuple(res,"kuku",0);
  nt->Project(hinam,"dy");
  delete nt;

  if (hist->GetEntries()<100) return;

  // find peak in the histogram 
  Int_t bpeak=(Int_t) hist->GetMaximumBin();
  Int_t nbins=(Int_t) hist->GetXaxis()->GetNbins();
  Int_t ypeak=(Int_t) hist->GetBinContent(bpeak);
  if (bpeak==0 || bpeak==nbins) return;

  // find the edges of the distribution
  Int_t bleft,brigt;
  {for (bleft=bpeak; bleft>0; bleft--) {
    if (hist->GetBinContent(bleft)<level*ypeak) break;
  }}
  {for (brigt=bpeak; brigt<nbins; brigt++) {
    if (hist->GetBinContent(brigt)<level*ypeak) break;
  }}
  bleft=bpeak+(bleft-bpeak);
  brigt=bpeak+(brigt-bpeak);
  if (bleft<0) bleft=0;
  if (brigt>nbins) brigt=nbins;
  Axis_t xleft=hist->GetBinCenter(bleft);
  Axis_t xrigt=hist->GetBinCenter(brigt);
  delete hist;
  
  // second, make copy of original tracks

  Int_t ntracks = res->GetNFilledTracks();
  if (ntracks<1) return;
  AliTrackPointArray **tmpvol = new AliTrackPointArray*[ntracks];
  AliTrackPointArray **tmptra = new AliTrackPointArray*[ntracks];
  for (Int_t itrack = 0; itrack < ntracks; itrack++){
    AliTrackPointArray *volarray;
    AliTrackPointArray *traarray;
    if (!res->GetTrackPointArrays(itrack, volarray, traarray)) {printf("cannot access track point array"); return;}
    tmpvol[itrack] = new AliTrackPointArray(*volarray);
    tmptra[itrack] = new AliTrackPointArray(*traarray);
  }

  // third, fire all employees and hire only the good ones

  res->SetNTracks(ntracks);
  {for (Int_t itrack = 0; itrack < ntracks; itrack++) {
    AliTrackPointArray *volarray = tmpvol[itrack];
    AliTrackPointArray *traarray = tmptra[itrack];
    int good_track=1;
    {for (Int_t ipoint = 0; ipoint < volarray->GetNPoints(); ipoint++) {
      AliTrackPoint p0,p1;
      traarray->GetPoint(p0,ipoint);
      volarray->GetPoint(p1,ipoint);
      //      res->GetAlignObj()->Transform(p1);
      AliTrackPoint q0,q1; // in tracking system
      Double_t alpha = p1.GetAngle();
      q0=p0.Rotate(alpha);
      q1=p1.Rotate(alpha);
      double y0 = q0.GetY();
      double z0 = q0.GetZ();
      double y1 = q1.GetY();
      double z1 = q1.GetZ();
      UShort_t volid = p1.GetVolumeID();
      AliGeomManager::ELayerID layer = AliGeomManager::VolUIDToLayer(volid);
      double dy = y1-y0-0.035*(z1-z0)*pow(-1.0,layer); 
      if (dy<xleft || dy>xrigt) good_track=0;
    }}
    if (good_track) res->AddTrackPointArrays(volarray, traarray);
  }}
}
//=============================================================================
void drawWithTitles(TTree *tr, char *what, char *xtitle, char *ytitle) {
  tr->Draw(what);
  tr->GetHistogram()->GetXaxis()->SetTitle(xtitle);
  tr->GetHistogram()->GetYaxis()->SetTitle(ytitle);
  tr->GetHistogram()->DrawCopy();
}
//=============================================================================
TCanvas *show(TNtuple *nt) {
  printf("%d track points\n",(int) nt->GetEntries());
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleOffset(3.5,"y");
  gStyle->SetTitleOffset(2.5,"x");
  gStyle->SetStatX(0.99);
  gStyle->SetStatY(0.99);
  gStyle->SetStatW(0.3);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.10);
  TCanvas *c = new TCanvas();
  c->Divide(3,2);
  c->cd(1); drawWithTitles(nt, "y1-y0:y0>>hi(60,-60,60,40,-2,2)",   "y_{track} (cm)", "y-y_{track} (cm)");
  c->cd(2); drawWithTitles(nt, "y1-y0:x0>>hi(60,280,400,40,-2,2)",  "x_{track} (cm)", "y-y_{track} (cm)");
  c->cd(3); drawWithTitles(nt, "y1-y0:z0>>hi(70,-350,350,40,-2,2)", "z_{track} (cm)", "y-y_{track} (cm)");
  c->cd(4); drawWithTitles(nt, "y1-y0:z1-z0>>hi(50,-10,10,40,-2,2)", "z-z_{track} (cm)", "y-y_{track} (cm)");
  c->cd(5); drawWithTitles(nt, "dy>>dyhist(150,-2,4)", "dy +- 0.035*dz", ""); // same binning as in FilterTracks
  //  double par[10]; 
  //  peakfit(nt->GetHistogram(),0.1,0.1,par);
  nt->GetHistogram()->Fit("gaus");
  ((TH1F*) gDirectory->Get("dyhist"))->DrawCopy();
  return c;
}
//=============================================================================
void init_alt(const char *trapo, int fitter_flag, int res_flag) {

  // prepare the AliAlignmentTracks object
  // read space points from file trapo

  alt = new AliAlignmentTracks();
  alt->SetPointsFilename(trapo);
  //alt->Misalign("kuku.root","TRDAlignObjs");
  
  if (fitter_flag == kStraight) {
    trf = new AliTrackFitterStraight();
  } else if (fitter_flag == kRieman) {
    trf = new AliTrackFitterRieman();
    ((AliTrackFitterRieman *) trf)->SetMaxDelta(6);
  } else if (fitter_flag == kKalman) {
    trf = new AliTrackFitterKalman();
    ((AliTrackFitterKalman *) trf)->SetMaxChi2(1000);
  }
  trf->SetMinNPoints(30); // if working with clusters
  trf->SetMinNPoints(4);  // if working with tracklets
  alt->SetTrackFitter(trf);

  if (res_flag == kLinear) {
    res = new AliTrackResidualsLinear();
    ((AliTrackResidualsLinear *)res)->SetRobust(0.7);
  } else if (res_flag == kFast) {
    res = new AliTrackResidualsFast();
  } else if (res_flag == kChi2) {
    res = new AliTrackResidualsChi2();
    res->FixParameter(0,0);
    res->FixParameter(1,0);
    res->FixParameter(2,0);
    res->FixParameter(3,0);
    res->FixParameter(4,0);
    // only rotation around z allowed
  }
  res->SetMinNPoints(10); // if working with clusters
  res->SetMinNPoints(1);  // if working with tracklets
  alt->SetMinimizer(res);  

  if (nteva) delete nteva;
  nteva = new TNtuple("eva","eva","chamber:module:layer:iteration:p0:p1:p2:p3:p4:p5:fp0:fp1:fp2");
}
//=============================================================================
void misa_alt(UShort_t volid, Double_t *mis) {

  // Misalign volid in alt by (local) pa[6].
  // Actually, both fAlignObjs and fMisalignObjs are applied in alt.
  // The proper way to introduce misalignment would be via fMisalignObjs. 
  // Because of lazyness I am using fAlignObj (they exist already). 
  // As a consequence, one should not expect the resulting transformations 
  // to be the inverse of the introduced misalignment, but to be zero.

  AliAlignObjParams *al = (AliAlignObjParams *) alt->GetAlignObj(volid);
  al->SetLocalPars(mis[0],mis[1],mis[2],mis[3],mis[4],mis[5]);
}
//=============================================================================
int align_volume(TArrayI volIds, TArrayI volIdRefs, int iterations, 
		 double level, double *pa, int showflag) {

  // fit tracks in volumes volIdRefs
  // determine residuals in volIds and save on AliTrackResiduals res
  // find alignment parameters such as to minimize the residuals
  // repeat the procedure iteration times
  // return resulting params in pa
  // if level>0 then suppress tracks from the tails of the distribution
  // tails being defined as the places where the distribution falls below 
  // level (relative to the peak height)
  // showflag = 1 show the residuals before minimization
  // showflag = 2 show the residuals after minimization

  // align

  printf("Aligning volume");
  if (volIds.GetSize()>1) printf("s");
  for (int i=0; i<volIds.GetSize(); i++) printf(" %d (%s)",volIds.At(i),AliGeomManager::SymName(volIds.At(i)));
  printf(" to reference volume");
  if (volIdRefs.GetSize()>1) printf("s");
  for (int i=0; i<volIdRefs.GetSize(); i++) printf(" %d (%s)",volIdRefs.At(i),AliGeomManager::SymName(volIdRefs.At(i)));
  printf("\n");

  int resu = alt->AlignVolumes(&volIds, &volIdRefs, AliGeomManager::kTPC1, AliGeomManager::kTOF,iterations);
  if (!resu) return 0;

  // repeat taking only peak

  if (level) {
    filterTracks(res,level);
    //    alt->AlignVolumes(&volIds, &volIdRefs, AliGeomManager::kTPC1, AliGeomManager::kTOF,iterations); 
    // this would reload res... just do
    AliAlignObjParams oldal(*res->GetAlignObj());   // old alignment object
    res->Minimize();
    AliAlignObjParams newal(*res->GetAlignObj());   // new alignment object
    AliAlignObjParams delta(oldal.Inverse());
    delta *= newal;                                 // new/old
    oldal.GetPars(pa,pa+3); printf("oldal    %10.3f\n",pa[1]);
    newal.GetPars(pa,pa+3); printf("newal    %10.3f\n",pa[1]);
    delta.GetPars(pa,pa+3); printf("delta    %10.3f\n",pa[1]);
    if (alt->GetUpdate()) for (Int_t i = 0; i < volIds.GetSize(); i++) {
      AliAlignObj *alignObj = alt->GetAlignObj(volIds.At(i));
      *alignObj *= delta;
    }
  }

  printf("minimizer: %12d tracks\n",res->GetNFilledTracks());
  printf("minimizer: %12d degrees of freedom\n",res->GetNdf());
  printf("minimizer: %12.2f chi2\n",res->GetChi2());
  printf("minimizer: %12.2f normalized chi2\n",res->GetChi2()/res->GetNdf());

  {for (int i=0; i<volIds.GetSize(); i++) {
    //  TGeoHMatrix m;
    //  res->GetAlignObj()->GetMatrix(m);
    //  alt->GetAlignObj(volIds.At(i))->SetMatrix(m); // delta matrix rather than final matrix, for debugging
    alt->GetAlignObj(volIds.At(i))->GetLocalPars(pa,pa+3);
    printf("updated alignment object for %d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
	   volIds.At(i),pa[0],pa[1],pa[2],pa[3],pa[4],pa[5]);
  }}

  // show result

  if (showflag) {
    if (ntresi) delete ntresi;
    ntresi = makeNtuple(res,"kuku",showflag-1);
    TCanvas *c = show(ntresi);

    c->cd(6);
    TText te;
    te.SetTextFont(43);
    te.SetTextColor(1);
    te.SetTextSizePixels(11);
    double xpos = 0.1;
    double ypos = 1.0;
    te.SetTextAlign(13);
    char buf[2000];
    
    te.DrawTextNDC(xpos,ypos-=0.04,res->GetName());
    te.DrawTextNDC(xpos,ypos-=0.07,"Aligning volumes");
    for (int i=0; i<volIds.GetSize(); i++) {
      sprintf(buf,"%d (%s) ",volIds.At(i),AliGeomManager::SymName(volIds.At(i)));
      te.DrawTextNDC(xpos+0.05,ypos-=0.04,buf);
      if (volIds.GetSize()>4 && i==1) { // cut it short
	te.DrawTextNDC(xpos+0.05,ypos-=0.04,"...etc...");
	i=volIds.GetSize()-2;
      }
    }
    te.DrawTextNDC(xpos,ypos-=0.04,"to reference volumes");
    for (int i=0; i<volIdRefs.GetSize(); i++) {
      sprintf(buf,"%d (%s) ",volIdRefs.At(i),AliGeomManager::SymName(volIdRefs.At(i)));
      te.DrawTextNDC(xpos+0.05,ypos-=0.04,buf);
    }
    
    te.DrawTextNDC(xpos,ypos-=0.07,"Result");
    double ypos1 = ypos;
    te.DrawTextNDC(xpos+0.05,ypos1-=0.04,"shift in phi");
    te.DrawTextNDC(xpos+0.05,ypos1-=0.04,"shift in z");
    te.DrawTextNDC(xpos+0.05,ypos1-=0.04,"shift in r");
    te.DrawTextNDC(xpos+0.05,ypos1-=0.04,"tilt in phi");
    te.DrawTextNDC(xpos+0.05,ypos1-=0.04,"tilt in z");
    te.DrawTextNDC(xpos+0.05,ypos1-=0.04,"tilt in r");
    te.SetTextAlign(33);    
    sprintf(buf,"%.3f",pa[0]); te.DrawTextNDC(xpos+0.45,ypos-=0.04,buf);
    sprintf(buf,"%.3f",pa[1]); te.DrawTextNDC(xpos+0.45,ypos-=0.04,buf);
    sprintf(buf,"%.3f",pa[2]); te.DrawTextNDC(xpos+0.45,ypos-=0.04,buf);
    sprintf(buf,"%.4f",pa[3]); te.DrawTextNDC(xpos+0.45,ypos-=0.04,buf);
    sprintf(buf,"%.4f",pa[4]); te.DrawTextNDC(xpos+0.45,ypos-=0.04,buf);
    sprintf(buf,"%.4f",pa[5]); te.DrawTextNDC(xpos+0.45,ypos-=0.04,buf);
    te.SetTextAlign(13);    
    if (showflag==1) te.DrawTextNDC(xpos,ypos-=0.07,"Residuals shown before alignment");
    if (showflag==2) te.DrawTextNDC(xpos,ypos-=0.07,"Residuals shown after alignment");
    c->Print("c1.ps");
  }

  return res->GetNFilledTracks();
}
//=============================================================================
int align_volume(int module, int layer0, int layer1, int reflayer0, int reflayer1, 
		 int iterations, double level, double *pa, int showflag) {

  // align the chambers (module,layer0-layer1) to the layers reflayer0-reflayer1 
  // (but excluding the layers to be aligned)

  // prepare volume arrays
  int n = 0;
  int temp[1000]={};
  int sec = module/5;
  for (int i=reflayer0; i<=reflayer1; i++) {
    if (i>=layer0 && i<=layer1) continue;
    //if (i!=reflayer0 && i!=reflayer1) continue; // take only first and last ref layer
    if (i>=AliGeomManager::kTPC1 && i<=AliGeomManager::kTPC2) {
      // both z-halves of TPC
      temp[n++] = AliGeomManager::LayerToVolUID(i,sec);
      temp[n++] = AliGeomManager::LayerToVolUID(i,sec+18);
    }
    if (i>=AliGeomManager::kTRD1 && i<=AliGeomManager::kTRD6) {
      temp[n++] = AliGeomManager::LayerToVolUID(i,module);
    }
  }
  TArrayI volIdRefs(n);
  for (int i=0; i<n; i++) volIdRefs.AddAt(temp[i],i);

  TArrayI volIds(layer1-layer0+1);  
  
  for (int i=layer0; i<=layer1; i++)  volIds.AddAt(AliGeomManager::LayerToVolUID(i,module),i-layer0);

  return align_volume(volIds,volIdRefs,iterations,level,pa,showflag); 
}
//=============================================================================
int align_sm(int sec, int iterations, double level, double *pa, int showflag) {

  // align the whole supermodule to the TPC
  // could be unified with the align_volume above

  // prepare volume arrays
  TArrayI volIdRefs(4);
  volIdRefs.AddAt(AliGeomManager::LayerToVolUID(7,sec),0);
  volIdRefs.AddAt(AliGeomManager::LayerToVolUID(8,sec),1);
  volIdRefs.AddAt(AliGeomManager::LayerToVolUID(7,sec+18),2);
  volIdRefs.AddAt(AliGeomManager::LayerToVolUID(8,sec+18),3);

  TArrayI volIds(30);
  int n = 0;
  for (int layer=9; layer<=14; layer++) for (int module=5*sec; module<5*sec+5; module++) 
    volIds.AddAt(AliGeomManager::LayerToVolUID(layer,module),n++);

  return align_volume(volIds,volIdRefs,iterations,level,pa,showflag); 
}
//=============================================================================
void align_layer_modules(int layer, int reflayer0, int reflayer1, 
			 int iterations, double level) {

// Loop through modules (TRD chambers) within given layer. 
// Align each module using other layers.
// Store result on an ntuple. 

  double a[6]={};
  double b[6]={};
  ntalpa = new TNtuple("ntalpa","ntalpa","a0:a1:a2:a3:a4:a5:b0:b1:b2:b3:b4:b5");
  Int_t maxmodule = AliGeomManager::LayerSize(layer);
  {for (int i=0; i<maxmodule; i++) {
      printf("*************  module %d  of layer %d *************\n",i,layer);
      int ntra = align_volume(i,layer,layer,reflayer0,reflayer1,iterations,level,b,0); 
      if (ntra) ntalpa->Fill(a[0],a[1],a[2],a[3],a[4],a[5],b[0],b[1],b[2],b[3],b[4],b[5]);
  }}
}
//=============================================================================
void align_module_layers(int module, int iterations, double level) {

// Loop through the 4 inner layers (TRD chambers) within given stack. 
// Find the alignment needed to match each chamber to the other 5.
// Apply all 4 alignments. 
// Repeat this procedure iteration times. 

  int flayer = AliGeomManager::kTRD1;
  int llayer = AliGeomManager::kTRD6;
  char buf[1000];
  double par[10]; 
  TH1F *dy=0;
  for (int i=0; i<iterations; i++) {
    // go through the layers once
    {for (int layer=flayer+1; layer<llayer; layer++) {
      //if (layer==11) continue;
      align_volume(module, layer, layer, flayer, llayer, 1, level, aml[i][layer], 0);
      int sector = module/5;
      int stack = module%5;
      int det = AliTRDgeometry::GetDetector(layer-AliGeomManager::kTRD1,stack,sector);
      nteva->Fill(det,module,layer,i,
		  aml[i][layer][0],aml[i][layer][1],aml[i][layer][2],
		  aml[i][layer][3],aml[i][layer][4],aml[i][layer][5],
		  0,0,0);
      UShort_t voluid = AliGeomManager::LayerToVolUID(layer,module);
      alt->GetAlignObj(voluid)->SetLocalPars(0,0,0,0,0,0);
      if (ntresi) delete ntresi;
      ntresi = makeNtuple(res,"kuku",0);
      sprintf(buf,"dy>>hidy%02d_%d_befo(100,-2,2)",i,layer); 
      ntresi->Draw(buf);
      sprintf(buf,"hidy%02d_%d_befo",i,layer); 
      dy=(TH1F*) gDirectory->Get(buf);
      //      peakfit(dy,0.1,0.1,par);
      nteva->Fill(det,module,layer,i-0.2,0,0,0,0,0,0,
		  par[0],par[1],par[2]);

      if (ntresi) delete ntresi;
      ntresi = makeNtuple(res,"kuku",1);
      sprintf(buf,"dy>>hidy%02d_%d_afte(100,-2,2)",i,layer); 
      ntresi->Draw(buf);
      sprintf(buf,"hidy%02d_%d_afte",i,layer); 
      dy=(TH1F*) gDirectory->Get(buf);
      //      peakfit(dy,0.1,0.1,par);
      dy->DrawCopy();
      nteva->Fill(det,module,layer,i+0.2,0,0,0,0,0,0,
		  par[0],par[1],par[2]);
    }}
    // update all 6 alignment objects
    {for (int layer=flayer+1; layer<llayer; layer++) {
      UShort_t voluid = AliGeomManager::LayerToVolUID(layer,module);
      alt->GetAlignObj(voluid)->SetLocalPars(
					     aml[i][layer][0],
					     aml[i][layer][1],
					     aml[i][layer][2],
					     aml[i][layer][3],
					     aml[i][layer][4],
					     aml[i][layer][5]);
    }}
  }

}
//=============================================================================
void align_module_layers_plot(int iterations, double ymax) {
  // show histograms produced by align_module_layers

  int flayer = AliGeomManager::kTRD1;
  int llayer = AliGeomManager::kTRD6;
  int ncol = 8;
  if (2*iterations>ncol) ncol=2*iterations;
  TCanvas *c = new TCanvas();
  /*
  char *xdes[]={"befo","afte","befo","afte","befo","afte","befo","afte","befo","afte",
		"befo","afte","befo","afte","befo","afte","befo","afte","befo","afte"};
  char *ydes[]={"plane 5","plane 4","plane 3","plane 2","plane 1","plane 0"};
  multi(ncol,6, -2,2,0,ymax, 0,0, 4,5, "dy (cm)", "", xdes, ydes, "");
  */
  c->Divide(ncol,6);
  char buf[1000];
  for (int layer=flayer; layer<=llayer; layer++) {
    int col=0;
    for (int i=0; i<iterations; i++) {
      int row = AliGeomManager::kTRD6-layer;
      sprintf(buf,"hidy%02d_%d_befo",i,layer); 
      TH1F *hidy = (TH1F*) gDirectory->Get(buf);
      if (!hidy) break;
      hidy->SetMaximum(ymax);
      //printf("drawing %s %ld in pad %d\n",buf,hidy,row*ncol+i+1);
      c->cd(row*ncol+col+1); hidy->Draw("same"); col++;
      sprintf(buf,"hidy%02d_%d_afte",i,layer); 
      hidy = (TH1F*) gDirectory->Get(buf);
      hidy->SetMaximum(ymax);
      //printf("drawing %s %ld in pad %d\n",buf,hidy,row*ncol+i+1);
      c->cd(row*ncol+col+1); hidy->Draw("same"); col++;
      c->Modified();
      c->Update();
    }    
  }
}
//=============================================================================
void col(TNtuple *nt, int i) {
  nt->SetLineColor(i);
  nt->SetMarkerColor(i);
}
//=============================================================================
void align_module_layers_summary_plot() {

  //  didi(1);
  //  fosi(16);
  TCanvas *c = new TCanvas("c","c");
  c->Divide(2,2);
  TH2D *du0 = new TH2D("du0","du0",20,-0.5,9.5,100,-0.2,0.2);
  TH2D *du1 = new TH2D("du1","du1",20,-0.5,9.5,100,-0.2,0.2);
  TH2D *du2 = new TH2D("du2","du2",20,-0.5,9.5,100,0,0.5);
  du0->SetXTitle("iteration");
  du1->SetXTitle("iteration");
  du2->SetXTitle("iteration");
  du0->SetYTitle("suggested shift in phi (cm)");
  du1->SetYTitle("residuals peak position (cm)");
  du2->SetYTitle("residuals peak width (cm)");
  
  c->cd(1); du0->Draw();
  col(nteva,1); nteva->Draw("p0:iteration","fp0==0 && layer==10","same*l");
  col(nteva,2); nteva->Draw("p0:iteration","fp0==0 && layer==11","same*l");
  col(nteva,3); nteva->Draw("p0:iteration","fp0==0 && layer==12","same*l");
  col(nteva,4); nteva->Draw("p0:iteration","fp0==0 && layer==13","same*l");
  c->cd(3); du1->Draw();
  col(nteva,1); nteva->Draw("fp1:iteration","p0==0 && layer==10","same*l");
  col(nteva,2); nteva->Draw("fp1:iteration","p0==0 && layer==11","same*l");
  col(nteva,3); nteva->Draw("fp1:iteration","p0==0 && layer==12","same*l");
  col(nteva,4); nteva->Draw("fp1:iteration","p0==0 && layer==13","same*l");
  c->cd(4); du2->Draw();
  col(nteva,1); nteva->Draw("fp2:iteration","p0==0 && layer==10","same*l");
  col(nteva,2); nteva->Draw("fp2:iteration","p0==0 && layer==11","same*l");
  col(nteva,3); nteva->Draw("fp2:iteration","p0==0 && layer==12","same*l");
  col(nteva,4); nteva->Draw("fp2:iteration","p0==0 && layer==13","same*l");
  c->Print("amlsp.gif");
}
//=============================================================================
void align_module_layers_print(int iterations) {

  // print the result of align_module_layers

  int flayer = AliGeomManager::kTRD1;
  int llayer = AliGeomManager::kTRD6;
  for (int i=0; i<iterations; i++) {
    printf("\niteration %d\n",i);
    for (int layer=flayer; layer<=llayer; layer++) {
      printf("layer %2d ",layer);
      for (int j=0; j<6; j++) printf("%10.3f ",aml[i][layer][j]);
      printf("\n");
    }
  }
}
//=============================================================================
void align_totpc(char *inpfil, int what, int fitter_flag, int minim_flag, char *outfil) {

  // loop through detectors and align 
  // track points will be read from inpfil??.root
  // what=0 ... detectors
  // what=1 ... stacks
  // what=2 ... supermodules
  // result will be stored on outfil.dat and outfil.stat

  AliTRDalignment ch; // misalignment in term of chambers
  AliTRDalignment sm; // misalignment in terms of supermodule (applies if what==3)
  int smpat[18] = {1,1,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1};
  int ntracks[540] = {};
  enum {kDetector,kStack,kSupermodule};
  for (int sec=0; sec<18; sec++) {
    if (!smpat[sec]) continue;
    TString filnam = Form("%s%02d.root",inpfil,sec);
    init_alt(filnam.Data(), fitter_flag, minim_flag);
    double pa[6];
    int ntra=0;
    if (what==kSupermodule) {  // align supermodels
      ntra = align_sm(sec, 1,0.1,pa,1);
      // use the mean between the two central chambers for the whole supermodule
      double pa0[6];
      double pa1[6];
      alt->GetAlignObj(AliGeomManager::LayerToVolUID(11,sec*5+2))->GetLocalPars(pa0,pa0+3);
      alt->GetAlignObj(AliGeomManager::LayerToVolUID(12,sec*5+2))->GetLocalPars(pa1,pa1+3);
      for (int i=0; i<6; i++) pa[i] = (pa0[i]+pa1[i])/2.0;
      sm.SetSm(sec,pa);
    }
    for (int st=0; st<5; st++) {
      int module = sec*5+st;
      if (what==kStack) ntra = align_volume(module, 9,14, 7,8, 1,0.1,pa,2); // align stacks
      for (int lay=0; lay<6; lay++) {
	if (sec==8 && st==4 && lay==3) continue; 
	if (what==kDetector) ntra = align_volume(module, 9+lay,9+lay, 7,8, 1,0.1,pa,0);
	int det = AliTRDgeometry::GetDetector(lay,st,sec);
	ntracks[det] = ntra;
	alt->GetAlignObj(AliGeomManager::LayerToVolUID(9+lay,module))->GetLocalPars(pa,pa+3); // alternative way
	ch.SetCh(det,pa);
	//	getchar();
      }
    }
  }
  if (what==kSupermodule) sm.WriteAscii(Form("%ssm.dat",outfil));
  ch.WriteAscii(Form("%s.dat",outfil));
  FILE *fp=fopen(Form("%s.stat",outfil),"w");
  for (int i=0; i<540; i++) fprintf(fp,"%3d %12d\n",i,ntracks[i]);
  fclose(fp);
}
//=============================================================================
