/// \file TestAnalisys.C
///
/// ~~~{.cxx}
/// .L AliGenInfo.C+
/// .L TestAnalisys.C+
/// AddChains(868);    // AddChains(runNumber);
///  Select();          // make default selection of data
///  MakePictures("pic868");
/// ~~~

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TCut.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TEventList.h"

#include "AliTrackReference.h"
#include "AliTPCParam.h"
#include "AliDetector.h"
#include "AliStack.h" 
#include "AliGenInfo.h"






void Select();                // make default selection 
void SelectLaser();           // make default selection  for laser tracks

void AddChains(Int_t run);    // add all the trees with selected run number to the chain
void MakePictures(char *dirname);           // make default pictures

void PRFYZ(TCut cut0, TCut cut1,  char * description);  
void PRFZZ(TCut cut0, TCut cut1,  char * description);  
void P5Z(TCut cut0, TCut cut1,  char * description);  
void P3Z(TCut cut0, TCut cut1,  char * description);  
void ResYZ(TCut cut0, TCut cut1,  char * description);
void SysYX(TCut cut0,  char * description);
void SysZX(TCut cut0,  char * description);

TProfile * ProfileMaxRow(TCut cut0, char *name, Int_t max);
TProfile * ProfileMaxPhi(TCut cut0, char *name, Int_t max);
TProfile * ProfileMaxZ(TCut cut0, char *name, Int_t max);
TProfile * ProfileQRow(TCut cut0, char *name, Int_t max);
TProfile * ProfileQPhi(TCut cut0, char *name, Int_t max);
TProfile * ProfileQZ(TCut cut0, char *name, Int_t max);
TCanvas *  NoiseSector(TCut cut0,  char * description, Int_t maxrow, Int_t maxpad);

//
// global variables
//
TChain chaincl("Tracks","Tracks");   // tpc tracks and clusters
TChain chaincl2("Tracks","Tracks");   // tpc tracks and clusters
TChain chainSignal("SignalB","SignalB");   // signals over threshold 50

TChain chainFit("Fit","Fit");        // fitted signals with fit parameters
TChain chainPed("Fit","Fit");        // fitted pedestal with noise
TString runDesc="Run ";              // run descriptor
//
AliComparisonDraw comp;
AliComparisonDraw compCl;
AliComparisonDraw compF;
AliComparisonDraw compP;
//
//
// selection of data for analysis
//
TEventList * listTracks    = new TEventList("listTracks","listTracks");
TEventList * listFitS      = new TEventList("listFitS","listFitS");
TEventList * listFitPed    = new TEventList("listFitPed","listFitPed");




void MakePictures(char *dirname){
  /// Define Uli Style

  gROOT->SetStyle("Plain");
  gStyle->SetFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetStatColor(10);
  
  gStyle->SetPalette(1,0);
  gStyle->SetNumberContours(50);
  //
  const Int_t kMinCl = 200;
  char chshell[100];
  sprintf(chshell,"mkdir %s", dirname);
  gSystem->Exec(chshell); 
  sprintf(chshell,"cd %s", dirname);
  gSystem->Exec(chshell);
  //
  //
  //
  TCanvas c(dirname,dirname);
  for (Int_t isector=0; isector<36; isector++){
    char chcut1[100];
    char chcut2[100];
    char chdesc[100];
    TProfile * prof;
    sprintf(chshell,"Cl.fX>0&&Cl.fDetector==%d",isector);
    Int_t ncl = comp.fTree->Draw("Cl.fY",chshell);
    if (ncl<kMinCl) continue;
    printf("MakePictures sector %d\n",isector);
    //
    // charge pictures
    //
    //
    // charge row
    //
    c.cd();
    sprintf(chdesc,"%s Sector %d IROC",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d", isector);
    prof = ProfileQRow(chcut1, chdesc, 70);
    sprintf(chshell,"%s/qrow_sec%dIROC.eps", dirname,isector);
    prof->Draw();
    c.Update();
    c.Print(chshell);
    prof = ProfileMaxRow(chcut1, chdesc, 70);
    sprintf(chshell,"%s/maxrow_sec%dIROC.eps", dirname,isector);
    prof->Draw();
    c.Update();
    c.Print(chshell);
    //
    sprintf(chdesc,"%s Sector %d OROC",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d", isector+36);
    prof = ProfileQRow(chcut1, chdesc, 100);
    sprintf(chshell,"%s/qrow_sec%dOROC.eps", dirname,isector);
    prof->Draw();
    c.Update();
    c.Print(chshell);
    prof = ProfileMaxRow(chcut1, chdesc, 100);
    sprintf(chshell,"%s/maxrow_sec%dOROC.eps", dirname,isector);
    prof->Draw();
    c.Update();
    c.Print(chshell);
    //
    // charge phi
    //
    sprintf(chdesc,"%s Sector %d IROC",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d", isector);
    prof = ProfileMaxPhi(chcut1, chdesc,20);
    sprintf(chshell,"%s/qphi_sec%dIROC.eps", dirname,isector);
    prof->Draw();
    c.Update();
    c.Print(chshell);
    //
    sprintf(chdesc,"%s Sector %d OROC",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d", isector+36);
    prof = ProfileMaxPhi(chcut1, chdesc,20);
    sprintf(chshell,"%s/qphi_sec%dOROC.eps", dirname,isector);
    prof->Draw();
    c.Update();
    c.Print(chshell);
    //
    //   charge z
    //
    c.cd();
    sprintf(chdesc,"%s Sector %d IROC",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d", isector);
    prof = ProfileQZ(chcut1, chdesc,20);
    sprintf(chshell,"%s/qz_sec%dIROC.eps", dirname,isector);
    //    prof->Draw();
    c.Update();
    c.Print(chshell);
    //
    c.cd();
    sprintf(chdesc,"%s Sector %d OROC",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d", isector+36);
    prof = ProfileQZ(chcut1, chdesc,20);
    sprintf(chshell,"%s/qz_sec%dOROC.eps", dirname,isector);
    //prof->Draw();
    c.Update();
    c.Print(chshell);
    //
    // Picture noise
    //
    sprintf(chdesc,"%s Sector %d IROC",runDesc.Data(), isector);
    sprintf(chcut1,"Sector==%d", isector);
    TCanvas *cnoise = NoiseSector(chcut1, chdesc,70,70);
    sprintf(chshell,"%s/noise_sec%dIROC.eps", dirname,isector);
    cnoise->Print(chshell);
    sprintf(chdesc,"%s Sector %d OROC",runDesc.Data(), isector);
    sprintf(chcut1,"Sector==%d", isector+36);
    cnoise = NoiseSector(chcut1, chdesc,70,70);
    sprintf(chshell,"%s/noise_sec%dOROC.eps", dirname,isector);
    cnoise->Print(chshell);
    //
    // systematic
    //
    c.cd();
    sprintf(chdesc,"%s Sector %d",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d||Cl.fDetector==%d", isector, isector+36);
    SysYX(chcut1,chdesc);
    sprintf(chshell,"%s/deltayx_sec%d.eps", dirname,isector);
    c.Print(chshell);
    c.cd();
    SysZX(chcut1,chdesc);
    sprintf(chshell,"%s/deltazx_sec%d.eps", dirname,isector);
    c.Print(chshell);
    
    //
    // picture prf
    //  
    if (ncl<500) continue;  //not enough statistic
    //
    sprintf(chdesc,"%s Sector %d",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d", isector);
    sprintf(chcut2,"Cl.fDetector==%d", isector+36);
    c.cd();
    PRFYZ(chcut1, chcut2,chdesc);
    sprintf(chshell,"%s/prfyz_sec%d.eps", dirname,isector);
    c.Print(chshell);
    sprintf(chcut1,"Sector==%d", isector);
    sprintf(chcut2,"Sector==%d", isector+36);
    PRFZZ(chcut1, chcut2,chdesc);
    sprintf(chshell,"%s/prfzz_sec%d.eps", dirname,isector);
    c.Print(chshell);
    //
    // y resolution
    //
    sprintf(chdesc,"%s Sector %d",runDesc.Data(), isector);
    sprintf(chcut1,"Cl.fDetector==%d", isector);
    sprintf(chcut2,"Cl.fDetector==%d", isector+36);
    c.cd();
    ResYZ(chcut1, chcut2,chdesc);
    sprintf(chshell,"%s/resyz_sec%d.eps", dirname,isector);
    c.Print(chshell);
    //
  }
}






void AddChains(Int_t run){
  /// add files to the chains + check consistency

  ifstream in0;
  ifstream in1;
  ifstream in2;
  ifstream in3;
  ifstream in4;
  TString sfile;
  char strcl[100];
  runDesc+=run;
  // TPC tracks
  //
  sprintf(strcl,"ls  *%d*/TPCtracks.root > files.txt", run);
  gSystem->Exec(strcl);
  in0.open("files.txt");
  for (;in0>>sfile;){
    if (sfile.Length()==0) break;
    printf("%s\n",sfile.Data());
    TFile f(sfile.Data());
    TTree * tree = (TTree*)f.Get("Tracks");
    if (tree){      
      f.Close();
      chaincl.Add(sfile.Data());
    }
  }
  //
  // Fitted signals
  sprintf(strcl,"ls  *%d*/FitSignal.root > files.txt", run);
  gSystem->Exec(strcl);
  in1.open("files.txt");
  for (;in1>>sfile;){
    if (sfile.Length()==0) break;
    printf("%s\n",sfile.Data()); 
    TFile f(sfile.Data());
    TTree * tree =(TTree*)f.Get("Fit");
    if (tree){      
      f.Close();
      chainFit.Add(sfile.Data());
    }
  }
  //
  // Fitted pedestal
  sprintf(strcl,"ls  *%d*/TPCsignal.root > files.txt", run);
  gSystem->Exec(strcl);
  in2.open("files.txt");
  for (;in2>>sfile;){
    if (sfile.Length()==0) break;
    printf("%s\n",sfile.Data());
    TFile f(sfile.Data());
    TTree * tree =(TTree*)f.Get("Fit");
    if (tree){      
      f.Close();
      chainPed.Add(sfile.Data());
    }
    //    chainPed.Add(sfile.Data());
  }
  //
  // Random signals
  sprintf(strcl,"ls  *%d*/TPCsignal.root > files.txt", run);
  gSystem->Exec(strcl);
  in4.open("files.txt");
  for (;in4>>sfile;){
    if (sfile.Length()==0) break;
    printf("%s\n",sfile.Data());
    TFile f(sfile.Data());
    TTree * tree =(TTree*)f.Get("SignalB");
    if (tree){      
      f.Close();
      chainSignal.Add(sfile.Data());
    }
    //    chainPed.Add(sfile.Data());
  }
  //
  // Rec points trees
  //
  printf("\n IMPORT REC points");
  sprintf(strcl,"ls  *%d*/*RecPoints* > files.txt", run);
  gSystem->Exec(strcl);
  in3.open("files.txt");
  for (;in3>>sfile;){
    if (sfile.Length()==0) break;
    printf("%s\n",sfile.Data());    
    TFile fcl(sfile.Data());
    char tname[100];
    sprintf(tname,"%s/%s/TreeR",sfile.Data(),fcl.GetListOfKeys()->At(0)->GetName());
    chaincl2.Add(tname);
    //    chainPed.Add(sfile.Data());
  }

  comp.fTree = &chaincl;
  compF.fTree = &chainFit;
  compP.fTree = &chainPed;
}

void Select(){
  /// base cut on the tracks

  comp.fTree->Draw(">>listTracks","Etrack.fTPCncls>30&&abs(Etrack.fIp.fP[4])<1");
  comp.fTree->SetEventList(listTracks);
  //
  compF.fTree->Draw(">>listFitS","p2>0&&p2<5&&p1<900&&p0<10000&&p4<1&&p4>0&&p5<p3&&chi2<150");
  compF.fTree->SetEventList(listFitS);
}

void SelectLaser(){
  /// base cut on the tracks

  comp.fTree->Draw(">>listTracks","Etrack.fTPCncls>20&&abs(Etrack.fIp.fP[4])<1&&abs(Etrack.fIp.fP[3])<0.01");
  comp.fTree->SetEventList(listTracks);
  //
  compF.fTree->Draw(">>listFitS","p2>0&&p2<5&&p1<900&&p0<10000&&p4<1&&p4>0&&p5<p3&&chi2<150");
  compF.fTree->SetEventList(listFitS);
  //
  // make default aliases
  //
  //  laser z beam
  comp.fTree->SetAlias("lz0","abs(Etrack.fIp.fP[1]-20)<5");
  comp.fTree->SetAlias("lz1","abs(Etrack.fIp.fP[1]-70)<20");
  comp.fTree->SetAlias("lz2","abs(Etrack.fIp.fP[1]-150)<20");
  comp.fTree->SetAlias("lz3","abs(Etrack.fIp.fP[1]-210)<20");

}



void PRFYZ(TCut cut0, TCut cut1,  char * description){
  /// plot Pad response function as funtion of drift z

  TF1 * f1 = new TF1("fdiff","sqrt([0]*[0]+(250-x)*[1]*[1])");
  f1->SetParameter(1,0.2);
  f1->SetParameter(0,0.2);
  comp.DrawXY("abs(Cl.fZ)","sqrt(Cl.fSigmaY2)","abs(Track.fTrackPoints.GetAngleY())<0.05","Track.fTrackPoints.fTX>0"+cut0,5,10,240,-0,1);
  TH1F * prfInnerY = (TH1F*)comp.fMean->Clone();

  comp.DrawXY("abs(Cl.fZ)","sqrt(Cl.fSigmaY2)","abs(Track.fTrackPoints.GetAngleY())<0.05","Track.fTrackPoints.fTX>0"+cut1,5,10,240,-0,1);
  TH1F * prfOuterY = (TH1F*)comp.fMean->Clone();
  //
  //
  prfOuterY->SetMinimum(0);
  prfOuterY->SetMarkerStyle(23);
  prfInnerY->SetMarkerStyle(24);
  prfOuterY->SetXTitle("Z position (cm)");
  prfOuterY->SetYTitle("PRF width (cm)");
  char chouter[100];
  char chinner[100];
  prfOuterY->Fit(f1);
  sprintf(chouter,"Outer sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfInnerY->Fit(f1);
  sprintf(chinner,"Inner sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfOuterY->Draw();
  prfInnerY->Draw("same");
  TString desc = description;
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.35, desc+"\nTPC cluster shape Fit: #sigma = #sqrt{p_{0}^{2}+(z_{d}-z)p_{1}^{2}}");
  legend->SetBorderSize(1);
  legend->AddEntry(prfOuterY,chouter);
  legend->AddEntry(prfInnerY,chinner);
  legend->Draw();
}



void PRFZZ(TCut cut0, TCut cut1,  char * description){
  TF1 * f1 = new TF1("fdiff","sqrt([0]*[0]+x*[1]*[1])");
  f1->SetParameter(1,0.2);
  f1->SetParameter(0,0.2);
  compF.DrawXY("p1*0.285","p2*0.285","p2>0",cut0,8,20,250,-0,2);
  TH1F * prfInnerY = (TH1F*)compF.fMean->Clone();
  compF.DrawXY("p1*0.285","p2*0.285","p2>0",cut1,8,20,250,-0,2);
  TH1F * prfOuterY = (TH1F*)compF.fMean->Clone();
  //
  //
  prfOuterY->SetMinimum(0);
  prfOuterY->SetMarkerStyle(23);
  prfInnerY->SetMarkerStyle(24);
  prfOuterY->SetXTitle("Drift length(cm)");
  prfOuterY->SetYTitle("Z Sigma (cm)");
  char chouter[100];
  char chinner[100];
  prfOuterY->Fit(f1);
  sprintf(chouter,"Outer sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfInnerY->Fit(f1);
  sprintf(chinner,"Inner sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfOuterY->Draw();
  prfInnerY->Draw("same");
  TString desc = description;
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.35, desc+"TPC signal shape Fit: #sigma = #sqrt{p_{0}^{2}+(z)p_{1}^{2}}");
  legend->SetBorderSize(1);
  legend->AddEntry(prfOuterY,chouter);
  legend->AddEntry(prfInnerY,chinner);
  legend->Draw();
}


void ResYZ(TCut cut0, TCut cut1,  char * description){
  /// resolution in y coordinate as function of z

  TF1 * f1 = new TF1("fdiff","sqrt([0]*[0]+(250-x)*[1]*[1])");
  f1->SetParameter(1,0.2);
  f1->SetParameter(0,0.2);
  comp.DrawXY("abs(Cl.fZ)","Track.fTrackPoints.GetY()-Cl.GetY()","abs(Track.fTrackPoints.GetAngleY())<0.05","Track.fTrackPoints.fTX>0"+cut0,5,10,240,-0.5,0.5);
  TH1F * prfInnerY = (TH1F*)comp.fRes->Clone();

  comp.DrawXY("abs(Cl.fZ)","Track.fTrackPoints.GetY()-Cl.GetY()","abs(Track.fTrackPoints.GetAngleY())<0.05","Track.fTrackPoints.fTX>0"+cut1,5,10,240,-0.5,0.5);
  TH1F * prfOuterY = (TH1F*)comp.fRes->Clone();
  //
  //
  prfOuterY->SetMinimum(0);
  prfOuterY->SetMaximum(0.15);  
  prfOuterY->SetMarkerStyle(23);
  prfInnerY->SetMarkerStyle(24);
  prfOuterY->SetXTitle("Z position (cm)");
  prfOuterY->SetYTitle("Y resolution (cm)");
  char chouter[100];
  char chinner[100];
  prfOuterY->Fit(f1);
  sprintf(chouter,"Outer sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfInnerY->Fit(f1);
  sprintf(chinner,"Inner sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfOuterY->Draw();
  prfInnerY->Draw("same");
  TString desc = description;
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.35, desc+"TPC cluster resolution: #sigma = #sqrt{p_{0}^{2}+(z_{d}-z)p_{1}^{2}}");
  legend->SetBorderSize(1);
  legend->AddEntry(prfOuterY,chouter);
  legend->AddEntry(prfInnerY,chinner);
  legend->Draw();
}

void SysYX(TCut cut0,  char * description){
  ///

  TProfile * profA = new TProfile("profY","profY",70,89,250);
  comp.fTree->Draw("Cl.fY-Track.fTrackPoints.GetY():Track.fTrackPoints.GetX()>>profY","abs(Cl.fY-Track.fTrackPoints.GetY())<1&&Track.fTrackPoints.fTX>10"+cut0,"prof");
  profA->SetXTitle("Local X (cm)");
  profA->SetYTitle("Mean #Delta Y (cm)");  
  TLegend *legend = new TLegend(0.55,0.25,0.85,0.30, description);
  legend->Draw();
}

void SysZX(TCut cut0,  char * description){
  ///

  TProfile * profA = new TProfile("profZ","profZ",70,89,250);
  comp.fTree->Draw("abs(Cl.fZ)-abs(Track.fTrackPoints.GetZ()):Track.fTrackPoints.GetX()>>profZ","abs(abs(Cl.fZ)-abs(Track.fTrackPoints.GetZ()))<1&&Track.fTrackPoints.fTX>10"+cut0,"prof");
  profA->SetXTitle("Local X (cm)");
  profA->SetYTitle("Mean #Delta Z (cm)");  
  TLegend *legend = new TLegend(0.55,0.25,0.85,0.30, description);
  legend->Draw(); 
}

TProfile * ProfileMaxRow(TCut cut0, char *name, Int_t max){ 
  /// make profile histrogram of amplitudes

  TProfile *profA = new TProfile(name,name,max,0,max-1);
  char expr[100];
  sprintf(expr,"Cl.fMax:Cl.fRow>>%s",name);
  comp.fTree->Draw(expr,"abs(Cl.fZ)>0&&Cl.fMax<500"+cut0,"prof");
  profA->SetXTitle("Pad Row");
  profA->SetYTitle("Amplitude at maxima (ADC)");
  return profA;
}

TProfile * ProfileMaxPhi(TCut cut0, char *name, Int_t max){ 
  /// make profile histrogram of amplitudes

  TProfile *profA = new TProfile(name,name,max,-0.14,0.14);
  char expr[100];
  sprintf(expr,"Cl.fMax:Cl.fY/Cl.fX>>%s",name);
  comp.fTree->Draw(expr,"abs(Cl.fZ)>0&&Cl.fMax<500"+cut0,"prof");
  profA->SetXTitle("Local #phi(rad)");
  profA->SetYTitle("Amplitude at maxima (ADC)");
  return profA; 
}

TProfile * ProfileQRow(TCut cut0, char *name, Int_t max){ 
  /// make profile histrogram of amplitudes

  TProfile *profA = new TProfile(name,name,max,0,max-1);
  char expr[100];
  sprintf(expr,"Cl.fQ:Cl.fRow>>%s",name);
  comp.fTree->Draw(expr,"abs(Cl.fZ)>0&&Cl.fMax<500"+cut0,"prof");
  profA->SetXTitle("Pad Row");
  profA->SetYTitle("Total charge(ADC)");
  return profA;
}

TProfile * ProfileQPhi(TCut cut0, char *name, Int_t max){ 
  /// make profile histrogram of amplitudes

  TProfile *profA = new TProfile(name,name,max,-0.14,0.14);
  char expr[100];
  sprintf(expr,"Cl.fQ:Cl.fY/Cl.fX>>%s",name);
  comp.fTree->Draw(expr,"abs(Cl.fZ)>0&&Cl.fMax<500"+cut0,"prof");
  profA->SetXTitle("Local #phi(rad)");
  profA->SetYTitle("Total charge (ADC)");
  return profA;
}

TProfile * ProfileQZ(TCut cut0, char *name, Int_t max){ 
  /// make profile histrogram of amplitudes

  TF1 * fline = new TF1("fline","[0]+[1]*[0]*(250-x)");
  TF1 * f1 = new TF1("f1","[0]*exp(-[1]*(250-x))");
  TProfile *profA = new TProfile(name,name,max,0,250);
  char expr[100];
  sprintf(expr,"Cl.fQ:abs(Cl.fZ)>>%s",name);
  comp.fTree->Draw(expr,"abs(Cl.fZ)>0&&Cl.fMax<500"+cut0,"prof");
  profA->SetXTitle("Z position (cm)"); 
  profA->SetYTitle("Amplitude (ADC)");
  char chc[100];
  profA->Fit(fline);
  f1->SetParameter(0,fline->GetParameter(0));
  f1->SetParameter(1,fline->GetParameter(1));
  profA->Fit(f1);
  sprintf(chc,"Exponential fit params: p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  printf("%s",chc);
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.25, chc);
  legend->Draw();
  return profA;
}

TProfile * ProfileMaxZ(TCut cut0, char *name, Int_t max){ 
  /// make profile histrogram of amplitudes

  TF1 * f1 = new TF1("f1","[0]+[1]*[0]*(250-x)");
  TProfile *profA = new TProfile(name,name,max,0,250);
  char expr[100];
  sprintf(expr,"Cl.fMax:abs(Cl.fZ)>>%s",name);
  comp.fTree->Draw(expr,"abs(Cl.fZ)>0&&Cl.fMax<500"+cut0,"prof");
  profA->SetXTitle("Z position (cm)"); 
  profA->SetYTitle("Amplitude at maxima (ADC)");
  char chc[100];
  profA->Fit(f1);
  sprintf(chc,"p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.35, chc);
  legend->Draw();
  return profA;
}


void P3Z(TCut cut0, TCut cut1,  char * description){
  /// first exponenent as function of z drift

  TF1 * f1 = new TF1("fdiff","[0]+[1]/[0]*x");
  f1->SetParameter(1,0.2);   
  f1->SetParameter(0,0.2);
  compF.DrawXY("p1*0.285","p3","Max>250",cut0,5,20,250,-0,2);
  TH1F * prfInnerY = (TH1F*)compF.fMean->Clone();
  compF.DrawXY("p1*0.285","p3","Max>250",cut1,5,20,250,-0,2);
  TH1F * prfOuterY = (TH1F*)compF.fMean->Clone();
  //
  //
  prfOuterY->SetMinimum(0);
  prfOuterY->SetMaximum(1);
  prfOuterY->SetMarkerStyle(23);
  prfInnerY->SetMarkerStyle(24);
  prfOuterY->SetXTitle("Drift length (cm)");
  prfOuterY->SetYTitle("Lambda 0 (Time Bin)");
  char chouter[100];
  char chinner[100];
  prfOuterY->Fit(f1);
  sprintf(chouter,"Outer sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfInnerY->Fit(f1);
  sprintf(chinner,"Inner sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfOuterY->Draw();
  prfInnerY->Draw("same");
  TString desc = description;
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.35, desc+"TPC cluster shape fit Parameter Lambda0 - P3:");
  legend->SetBorderSize(1);
  legend->AddEntry(prfOuterY,chouter);
  legend->AddEntry(prfInnerY,chinner);
  legend->Draw();
}


void P5Z(TCut cut0, TCut cut1,  char * description){
  /// second exponenent as function of z drift

  TF1 * f1 = new TF1("fdiff","[0]+[1]/[0]*x");
  f1->SetParameter(1,0.2);
  f1->SetParameter(0,0.2);
  compF.DrawXY("p1*0.285","p5","Max>250",cut0,5,20,250,-0,0.2);
  TH1F * prfInnerY = (TH1F*)compF.fMean->Clone();
  compF.DrawXY("p1*0.285","p5","Max>250",cut1,5,20,250,-0,0.2);
  TH1F * prfOuterY = (TH1F*)compF.fMean->Clone();
  //
  //
  prfOuterY->SetMinimum(0);
  prfOuterY->SetMaximum(0.15);
  prfOuterY->SetMarkerStyle(23);
  prfInnerY->SetMarkerStyle(24);
  prfOuterY->SetXTitle("Drift length (Time Bin)");
  prfOuterY->SetYTitle("Lambda1 (Time bin)");
  char chouter[100];
  char chinner[100];
  prfOuterY->Fit(f1);
  sprintf(chouter,"Outer sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfInnerY->Fit(f1);
  sprintf(chinner,"Inner sector : p_{0} = %f  p_{1} = %f",f1->GetParameter(0),f1->GetParameter(1));
  prfOuterY->Draw();
  prfInnerY->Draw("same");
   TString desc = description;
  TLegend *legend = new TLegend(0.25,0.12,0.85,0.35, desc+"TPC cluster shape fit Parameter Lambda1 - P5");
  legend->SetBorderSize(1);
  legend->AddEntry(prfOuterY,chouter);
  legend->AddEntry(prfInnerY,chinner);
  legend->Draw();
}

TCanvas *  NoiseSector(TCut cut0,  char * description, Int_t maxrow, Int_t maxpad){
  /// draw plots of the noise

  TCanvas * c = new TCanvas;
  c->Divide(2,1);
  c->Draw();
  c->cd(1);
  compP.fTree->Draw("GSigma","GSigma<5"+cut0);
  c->cd(2);
  Float_t rand  = gRandom->Rndm();
  char name[100];
  sprintf(name,"prof%f",rand);
  TProfile2D * prof= new TProfile2D(name,name,maxrow, 0, maxrow-1, 2*maxpad,-maxpad,maxpad);
  char expr[100];
  sprintf(expr,"GSigma:RPad:Row>>%s",name);
  prof->SetXTitle("Pad row");
  prof->SetYTitle("Pad number");
  compP.fTree->Draw(expr,cut0,"profcolz");
  c->cd(1);
  TString desc = description;
  TLegend *legend = new TLegend(0.25,0.30,0.85,0.85, desc+"Noise map");
  legend->Draw();
  
  return c;
}


