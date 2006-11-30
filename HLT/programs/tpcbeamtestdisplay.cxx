// $Id$

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include "AliLevel3.h"
#include "AliHLTTransform.h"
#include "AliHLTRawDataFileHandler.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTClustFinderNew.h"
#include "AliHLTConfMapper.h"
#include "AliHLTVertex.h"
#include "AliHLTMemHandler.h"
#include "AliHLTLogging.h"
#include "AliHLTTrackArray.h"
#include "AliHLTTrack.h"
#include <TCanvas.h>
#include <TH2.h>
#include <TView.h>
#include <TApplication.h>
#include <TColor.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TFormula.h>
#include <TF1.h>
#include <TLine.h>
#include <TAttLine.h>
#if __GNUC__== 3
using namespace std;
#else
#include <stream.h>
#include <string.h>
#include <stdlib.h>
#endif

//This program does the clusterfinding and the tracking.
int main(Int_t argc,Char_t **argv){

  //Display all clusters.
    
  Char_t cfile[1024];
  Char_t path[1024];  
  Char_t fname[1024];  

  if(argc<3){
    cout<<"Usage: tpcbeamtestdisplay filename path_to_files"<<endl;
    return -1;
  }
  if (argc>2) {
    sprintf(cfile,"%s",argv[1]);
    sprintf(path,"%s",argv[2]);
  }

  //displaying the tracks:
  UInt_t fNcl;
  AliHLTMemHandler *clusterfile = new AliHLTMemHandler();
  sprintf(fname,"%s/%s-points_0_-1.raw",path,cfile);
  cout<<"file: "<<fname<<endl;
  if(!clusterfile->SetBinaryInput(fname)){
    cout<<"file: "<<fname<<" can not be set as binary input!"<<endl;
    LOG(AliHLTLog::kError,"AliHLTEvaluation::Setup","File Open")
      <<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
    delete clusterfile;
    clusterfile = 0; 
  }
  clusterfile->SetBinaryInput(fname);
  AliHLTSpacePointData *fClusters = new AliHLTSpacePointData();
  memset(fClusters,0,sizeof(AliHLTSpacePointData*));
  fClusters = (AliHLTSpacePointData*)clusterfile->Allocate();
  clusterfile->Binary2Memory(fNcl,fClusters);
  clusterfile->CloseBinaryInput();
  TApplication theApp("App", &argc, argv);
  TCanvas *c1 = new TCanvas("c1","",900,700);
  c1->cd();

  c1->Clear();
  c1->SetFillColor(1);
  
  AliHLTSpacePointData *points = fClusters;
  TH2F *clusters = new TH2F("clusters","",65,-1,64,77,17,94);
  if(!points)
    cout<<"no points"<<endl;
  Int_t npoints = fNcl;
  //  cout<<"number of points: "<<npoints<<endl;
  for(Int_t i=0; i<npoints; i++){
    float x_tmp = points[i].fX;
    float y_tmp = points[i].fY;
    //    cout<<"x: "<<x_tmp<<" y: "<<y_tmp<<endl;
    clusters->Fill(x_tmp,y_tmp,1); 
  }
  TH2F *nothing = new TH2F(cfile,"",65,-1,64,77,17,94);
//  nothing->SetName("Test");
  //Adding a regression method to calculate the line parameters.
  float *x=new float[npoints];
  float *y=new float[npoints];
  float xyParamA=0;
  float xyParamB=0;
  
  Int_t numofXYPoints=npoints;

  float sumX=0;
  float sumY=0;

  float sumXtimesY=0;
  float sumXsquare=0;
  float sumYsquare=0;


  for(Int_t i=0; i<npoints; i++){
    x[i] = points[i].fX;
    y[i] = points[i].fY;
    sumX+=x[i];
    sumY+=y[i];
    sumXtimesY+=x[i]*y[i];
    sumXsquare+=x[i]*x[i];
    sumYsquare+=y[i]*y[i];
  }
  


  if(numofXYPoints*sumXsquare-sumX*sumX!=0)
    xyParamB=(numofXYPoints*sumXtimesY-sumX*sumY)/(numofXYPoints*sumXsquare-sumX*sumX);
  else
    cout<<"Divident is zero calculating the xParamB, numofXYPoints*sumXsquare-sumX*sumX=0"<<endl;
  if(numofXYPoints!=0)
    xyParamA=(sumY-xyParamB*sumX)/numofXYPoints;
  else
    cout<<"Divident is zero calculating the xParamA, numofXYPoints=0"<<endl;
  
  //cout<<"y= a + bx : "<<"y= "<<xyParamA<<" + "<<xyParamB<<" x"<<endl;

  TF1 *line = new TF1("line","[0]+[1]*x",0,63);
  line->SetParameters(xyParamA,xyParamB);
  TLine *leftline= new TLine(0,(0-AliHLTTransform::GetNPads(0)/2+55),0,AliHLTTransform::GetNPads(0)-(AliHLTTransform::GetNPads(0)/2)+55);
  TLine *rightline=new TLine(63,(0-AliHLTTransform::GetNPads(63)/2+55),63,AliHLTTransform::GetNPads(63)-(AliHLTTransform::GetNPads(63)/2)+55);
  TLine *upperline=new TLine(0,(AliHLTTransform::GetNPads(0)-AliHLTTransform::GetNPads(0)/2+55),63,AliHLTTransform::GetNPads(63)-(AliHLTTransform::GetNPads(63)/2)+55);
  TLine *underline=new TLine(0,(0-AliHLTTransform::GetNPads(0)/2+55),63,0-(AliHLTTransform::GetNPads(63)/2)+55);
  
  nothing->Draw("colz");
  nothing->SetXTitle("Padrows");
  nothing->SetYTitle("Pads");
  clusters->SetMarkerStyle(22);
  clusters->Draw("same");

  line->Draw("same");
  line->SetLineColor(2);
  leftline->Draw("same");
  leftline->SetLineWidth(2);
  rightline->Draw("same");
  rightline->SetLineWidth(2);
  upperline->Draw("same");
  upperline->SetLineWidth(2);
  underline->Draw("same");
  underline->SetLineWidth(2);
  

  c1->SetFillColor(10);
  c1->Update();
  c1->Draw();
  while(1){

  }
}
