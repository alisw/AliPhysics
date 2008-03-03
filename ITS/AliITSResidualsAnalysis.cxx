/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//   Implementation of the tracks residuals analysis class
//   It provides an access to the track space points
//   written along the esd tracks. The class enables
//   the user to plug any track fitter (deriving from
//   AliTrackFitter class) and minimization fo the
//   track residual sums (deriving from the AliTrackResiduals).
//-----------------------------------------------------------------

#include <Riostream.h> 
//#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TArrayI.h>
#include <TClonesArray.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include "TGeoMatrix.h"
#include "TGeoManager.h"
#include "TGeoPhysicalNode.h"
#include "TMatrixDSymEigen.h"
#include "TString.h"

#include "AliAlignmentTracks.h"
#include "AliTrackPointArray.h"
#include "AliAlignObjParams.h"
#include "AliTrackResiduals.h"
#include "AliTrackFitter.h"
#include "AliTrackFitterKalman.h"
#include "AliTrackFitterRieman.h"
#include "AliTrackResiduals.h"
#include "AliTrackResidualsChi2.h"
#include "AliTrackResidualsFast.h"
#include "AliLog.h"

#include "AliITSResidualsAnalysis.h"


// Structure of the RealignmentAnalysisLayer.C
  typedef struct {
    Int_t nentries;
    Float_t rms;
    Float_t meanFit;
    Float_t errmeanFit;
    Float_t sigmaFit;
  }  histProperties_t;

ClassImp(AliITSResidualsAnalysis)
  
//____________________________________________________________________________
  AliITSResidualsAnalysis::AliITSResidualsAnalysis():
    AliAlignmentTracks(),
    fnHist(0),
    fnPhi(0),
    fnZ(0),
    fvolidsToBin(0),
    fLastVolVolid(0), 
    fCoordToBinTable(0),
    fVolResHistRPHI(0),
    fResHistZ(0),
    fPullHistRPHI(0), 
    fPullHistZ(0), 
    fTrackDirPhi(0),
    fTrackDirLambda(0),
    fTrackDirLambda2(0),
    fTrackDirAlpha(0),
    fTrackDirPhiAll(0),
    fTrackDirLambdaAll(0),
    fTrackDirLambda2All(0),
    fTrackDirAlphaAll(0),
    fTrackDir(0), 
    fTrackDirAll(0), 
    fTrackDir2All(0),
    fTrackDirXZAll(0), 
    fResHistGlob(0),  
    fhistCorrVol(0),
    fVolNTracks(0),
    fhEmpty(0),
    fhistVolNptsUsed(0),
    fhistVolUsed(0),
    fSigmaVolZ(0),
    fsingleLayer(0),
    fWriteHist(0),
    fpTrackVolIDs(0),
    fVolVolids(0),
    fVolUsed(0),         
    fRealignObjFileIsOpen(kFALSE),
    fClonesArray(0),
    fAliTrackPoints("AliTrackPoints.root"),
    fGeom("geometry.root")

{

  //
  // Defaults
  //

 
}

//____________________________________________________________________________
AliITSResidualsAnalysis::AliITSResidualsAnalysis(const TString aliTrackPoints,  
						 const TString geom):
  AliAlignmentTracks(),
  fnHist(0),
  fnPhi(0),
  fnZ(0),
  fvolidsToBin(0),
  fLastVolVolid(0), 
  fCoordToBinTable(0),
  fVolResHistRPHI(0),
  fResHistZ(0),
  fPullHistRPHI(0), 
  fPullHistZ(0), 
  fTrackDirPhi(0),
  fTrackDirLambda(0),
  fTrackDirLambda2(0),
  fTrackDirAlpha(0),
  fTrackDirPhiAll(0),
  fTrackDirLambdaAll(0),
  fTrackDirLambda2All(0),
  fTrackDirAlphaAll(0),
  fTrackDir(0),
  fTrackDirAll(0),
  fTrackDir2All(0),
  fTrackDirXZAll(0),
  fResHistGlob(0),  
  fhistCorrVol(0),
  fVolNTracks(0),
  fhEmpty(0),
  fhistVolNptsUsed(0),
  fhistVolUsed(0),
  fSigmaVolZ(0),
  fsingleLayer(0),
  fWriteHist(0),
  fpTrackVolIDs(0),
  fVolVolids(0),
  fVolUsed(0),
  fRealignObjFileIsOpen(kFALSE),
  fClonesArray(0),
  fAliTrackPoints(aliTrackPoints),
  fGeom(geom)
{ 
  //
  // Constructor (alitrackpoints)
  //


}

//____________________________________________________________________________
AliITSResidualsAnalysis::AliITSResidualsAnalysis(const TArrayI *volIDs):
  AliAlignmentTracks(),
  fnHist(0),
  fnPhi(0),
  fnZ(0),
  fvolidsToBin(0),
  fLastVolVolid(0), 
  fCoordToBinTable(0),
  fVolResHistRPHI(0),
  fResHistZ(0),
  fPullHistRPHI(0), 
  fPullHistZ(0), 
  fTrackDirPhi(0),
  fTrackDirLambda(0),
  fTrackDirLambda2(0),
  fTrackDirAlpha(0),
  fTrackDirPhiAll(0),
  fTrackDirLambdaAll(0),
  fTrackDirLambda2All(0),
  fTrackDirAlphaAll(0),
  fTrackDir(0),
  fTrackDirAll(0),
  fTrackDir2All(0),
  fTrackDirXZAll(0),
  fResHistGlob(0),  
  fhistCorrVol(0),
  fVolNTracks(0),
  fhEmpty(0),
  fhistVolNptsUsed(0),
  fhistVolUsed(0),
  fSigmaVolZ(0),
  fsingleLayer(0),
  fWriteHist(0),
  fpTrackVolIDs(0),
  fVolVolids(0),
  fVolUsed(0),
  fRealignObjFileIsOpen(kFALSE),
  fClonesArray(0),
  fAliTrackPoints("AliTrackPoints.root"),
  fGeom("geometry.root")

{ 
  //
  // Original Constructor
  //

  InitHistograms(volIDs);

}

//____________________________________________________________________________
AliITSResidualsAnalysis::AliITSResidualsAnalysis(TArrayI *volIDs,AliTrackPointArray **tracksClustArray,AliTrackPointArray **tracksFitPointsArray):
  AliAlignmentTracks(),
  fnHist(0),
  fnPhi(0),
  fnZ(0),
  fvolidsToBin(0),
  fLastVolVolid(0), 
  fCoordToBinTable(0),
  fVolResHistRPHI(0),
  fResHistZ(0),
  fPullHistRPHI(0), 
  fPullHistZ(0),
  fTrackDirPhi(0),
  fTrackDirLambda(0),
  fTrackDirLambda2(0),
  fTrackDirAlpha(0),
  fTrackDirPhiAll(0),
  fTrackDirLambdaAll(0),
  fTrackDirLambda2All(0),
  fTrackDirAlphaAll(0),
  fTrackDir(0), 
  fTrackDirAll(0), 
  fTrackDir2All(0), 
  fTrackDirXZAll(0), 
  fResHistGlob(0),  
  fhistCorrVol(0),
  fVolNTracks(0),
  fhEmpty(0),
  fhistVolNptsUsed(0),
  fhistVolUsed(0),
  fSigmaVolZ(0),
  fsingleLayer(0),
  fWriteHist(0),
  fpTrackVolIDs(0),
  fVolVolids(0),
  fVolUsed(0),
  fRealignObjFileIsOpen(kFALSE),
  fClonesArray(0),
  fAliTrackPoints("AliTrackPoints.root"),
  fGeom("geometry.root")
  
{ 
  //
  // Detailed Constructor (deprecated)
  //


  TString histnameRPHI="HistRPHI_volID_",aux;
  TString histnameZ="HistZ_volID_";
  TString histnameGlob="HistGlob_volID_";
  TString histnameCorrVol="HistCorrVol_volID";
  TString histnamePullRPHI="HistPullRPHI_volID_";
  TString histnamePullZ="HistPullZ_volID_";
  fnHist=volIDs->GetSize();
  fVolResHistRPHI=new TH1F*[fnHist];
  fResHistGlob=new TH1F*[fnHist];
  fResHistZ=new TH1F*[fnHist];
  fPullHistRPHI=new TH1F*[fnHist];
  fPullHistZ=new TH1F*[fnHist];
  fhistCorrVol=new TH2F*[fnHist];
  Float_t **binningZPhi=CheckSingleLayer(volIDs);
  fvolidsToBin=new Int_t*[fnPhi*fnZ];
  Float_t *binningphi=binningZPhi[0];
  Float_t *binningz=binningZPhi[1];
  Bool_t binning=SetBinning(volIDs,binningphi,binningz);
  if(binning){
    fVolNTracks=new TH2F("fVolNTracks","Hist with N tracks passing through a given module==(r,phi) zone",fnPhi,binningphi,fnZ,binningz);
    fhistVolNptsUsed=new TH2F("fhistVolNptsUsed","Hist with N points used for given module==(r,phi) ",fnPhi,binningphi,fnZ,binningz);
    fhistVolUsed=new TH2F("fhistVolUsed","Hist with N modules used for a given module==(r,phi) zone",fnPhi,binningphi,fnZ,binningz);
    fSigmaVolZ=new TH2F("fSigmaVolZ","Hist with Sigma of Residual distribution for each module",fnPhi,binningphi,fnZ,binningz);
    fhEmpty=new TH2F("fhEmpty","Hist for getting binning",fnPhi,binningphi,fnZ,binningz);
    fVolNTracks->SetXTitle("Volume #phi");
    fVolNTracks->SetYTitle("Volume z ");
    fhistVolNptsUsed->SetXTitle("Volume #phi");
    fhistVolNptsUsed->SetYTitle("Volume z ");
    fhistVolUsed->SetXTitle("Volume #phi");
    fhistVolUsed->SetYTitle("Volume z ");
    fSigmaVolZ->SetXTitle("Volume #phi");
    fSigmaVolZ->SetYTitle("Volume z ");
  }
  else{
    fVolNTracks=new TH2F("fVolNTracks","Hist with N tracks passing through a given module==(r,phi) zone",50,-3.2,3.2,100,-80,80);
    fhistVolNptsUsed=new TH2F("fhistVolNptsUsed","Hist with N points used for given module==(r,phi) ",50,-3.2,3.2,100,-80,80);
    fhistVolUsed=new TH2F("fhistVolUsed","Hist with N modules used for a given module==(r,phi) zone",50,-3.2,3.2,100,-80,80);
    fSigmaVolZ=new TH2F("fSigmaVolZ","Hist with Sigma of Residual distribution for each module",50,-3.2,3.2,100,-80,80);
    fhEmpty=new TH2F("fhEmpty","Hist for getting binning",50,-3.2,3.2,100,-80,80);
    fVolNTracks->SetXTitle("Volume #phi");
    fVolNTracks->SetYTitle("Volume z ");
    fhistVolNptsUsed->SetXTitle("Volume #phi");
    fhistVolNptsUsed->SetYTitle("Volume z ");
    fhistVolUsed->SetXTitle("Volume #phi");
    fhistVolUsed->SetYTitle("Volume z ");
    fSigmaVolZ->SetXTitle("Volume #phi");
    fSigmaVolZ->SetYTitle("Volume z ");
  }
  
  fpTrackVolIDs=new TArrayI(fnHist);
  fVolUsed=new TArrayI*[fnHist];
  fVolVolids=new TArrayI*[fnHist]; 
  fLastVolVolid=new Int_t[fnHist];
 
  for (Int_t nhist=0;nhist<fnHist;nhist++){
    fpTrackVolIDs->AddAt(volIDs->At(nhist),nhist);   
    aux=histnameRPHI;
    aux+=volIDs->At(nhist);
    fVolResHistRPHI[nhist]=new TH1F("namehist","histname",200,-0.02,0.02);   
    fVolResHistRPHI[nhist]->SetName(aux.Data());
    fVolResHistRPHI[nhist]->SetTitle(aux.Data());
 
    aux=histnameZ;
    aux+=volIDs->At(nhist);
    fResHistZ[nhist]=new TH1F("histname","histname",400,-0.08,0.08);   
    fResHistZ[nhist]->SetName(aux.Data()); 
    fResHistZ[nhist]->SetTitle(aux.Data()); 

    aux=histnamePullRPHI;
    aux+=volIDs->At(nhist);
    fPullHistRPHI[nhist]=new TH1F("histname","histname",100,-7.,7.);   
    fPullHistRPHI[nhist]->SetName(aux.Data()); 
    fPullHistRPHI[nhist]->SetTitle(aux.Data()); 
    
    aux=histnamePullZ;
    aux+=volIDs->At(nhist);
    fPullHistZ[nhist]=new TH1F("histname","histname",100,-7.,7.);   
    fPullHistZ[nhist]->SetName(aux.Data()); 
    fPullHistZ[nhist]->SetTitle(aux.Data()); 
    
    aux=histnameGlob;
    aux+=volIDs->At(nhist);
    fResHistGlob[nhist]=new TH1F("histname","histname",400,-0.08,0.08);   
    fResHistGlob[nhist]->SetName(aux.Data()); 
    fResHistGlob[nhist]->SetTitle(aux.Data());    
    
    aux=histnameCorrVol;
    aux+=volIDs->At(nhist);
    fhistCorrVol[nhist]=new TH2F("histname","histname",50,-3.2,3.2,100,-80,80);   
    fhistCorrVol[nhist]->SetName(aux.Data()); 
    fhistCorrVol[nhist]->SetTitle(aux.Data()); 
    fhistCorrVol[nhist]->SetXTitle("Volume #varphi");
    fhistCorrVol[nhist]->SetYTitle("Volume z ");

    fVolVolids[nhist]=new TArrayI(1000);
    fVolUsed[nhist]=new TArrayI(1000);  
    fLastVolVolid[nhist]=0;   
    FillResHists(tracksClustArray[nhist],tracksFitPointsArray[nhist]);
  } 
  fWriteHist=kFALSE;
  DrawHists();

  SetFileNameTrackPoints("");  // Filename with the AliTrackPoints
  SetFileNameGeometry(""); // Filename with the Geometry


}

//____________________________________________________________________________
AliITSResidualsAnalysis::AliITSResidualsAnalysis(const AliITSResidualsAnalysis& obj):
    AliAlignmentTracks(),
    fnHist(0),
    fnPhi(0),
    fnZ(0),
    fvolidsToBin(0),
    fLastVolVolid(0), 
    fCoordToBinTable(0),
    fVolResHistRPHI(0),
    fResHistZ(0),
    fPullHistRPHI(0), 
    fPullHistZ(0), 
    fTrackDirPhi(0),
    fTrackDirLambda(0),
    fTrackDirLambda2(0),
    fTrackDirAlpha(0),
    fTrackDirPhiAll(0),
    fTrackDirLambdaAll(0),
    fTrackDirLambda2All(0),
    fTrackDirAlphaAll(0),
    fTrackDir(0), 
    fTrackDirAll(0), 
    fTrackDir2All(0),
    fTrackDirXZAll(0), 
    fResHistGlob(0),  
    fhistCorrVol(0),
    fVolNTracks(0),
    fhEmpty(0),
    fhistVolNptsUsed(0),
    fhistVolUsed(0),
    fSigmaVolZ(0),
    fsingleLayer(0),
    fWriteHist(0),
    fpTrackVolIDs(0),
    fVolVolids(0),
    fVolUsed(0),         
    fRealignObjFileIsOpen(kFALSE),
    fClonesArray(0),
    fAliTrackPoints("AliTrackPoints.root"),
    fGeom("geometry.root")

{
  // copy constructor. This is not allowed.

  AliFatal("Copy constructor not allowed\n");
 
}

//____________________________________________________________________________
AliITSResidualsAnalysis& AliITSResidualsAnalysis::operator = (const AliITSResidualsAnalysis& obj) {
  // assignment operator. This is not allowed
  AliFatal("Assignment operator not allowed\n");
  return *this;
}

//____________________________________________________________________________
AliITSResidualsAnalysis::~AliITSResidualsAnalysis()
{
  //
  //  Destructor
  //

  if(fvolidsToBin)        delete[] fvolidsToBin;         
  if(fLastVolVolid)       delete[] fLastVolVolid;        
  if(fCoordToBinTable)    delete[] fCoordToBinTable;
  if(fVolResHistRPHI)     delete fVolResHistRPHI; 
  if(fResHistZ)           delete fResHistZ;         
  if(fPullHistRPHI)       delete fPullHistRPHI;      
  if(fPullHistZ)          delete fPullHistZ;        
  if(fTrackDirPhi)        delete fTrackDirPhi;         
  if(fTrackDirLambda)     delete fTrackDirLambda;
  if(fTrackDirLambda2)    delete fTrackDirLambda2;     
  if(fTrackDirAlpha)      delete fTrackDirAlpha;   
  if(fTrackDirPhiAll)     delete fTrackDirPhiAll;     
  if(fTrackDirLambdaAll)  delete fTrackDirLambdaAll;   
  if(fTrackDirLambda2All) delete fTrackDirLambda2All;  
  if(fTrackDirAlphaAll)   delete fTrackDirAlphaAll;      
  if(fTrackDir)           delete fTrackDir;           
  if(fTrackDirAll)        delete fTrackDirAll;          
  if(fTrackDir2All)       delete fTrackDir2All;          
  if(fTrackDirXZAll)      delete fTrackDirXZAll;       
  if(fResHistGlob)        delete fResHistGlob;        
  if(fhistCorrVol)        delete fhistCorrVol;        
  if(fVolNTracks)         delete fVolNTracks;      
  if(fhEmpty)             delete fhEmpty;          
  if(fhistVolNptsUsed)    delete fhistVolNptsUsed;  
  if(fhistVolUsed)        delete fhistVolUsed;     
  if(fSigmaVolZ)          delete fSigmaVolZ;       
  if(fpTrackVolIDs)       delete fpTrackVolIDs;
  if(fVolVolids)          delete fVolVolids;
  if(fVolUsed)            delete fVolUsed;    
  if(fClonesArray)        delete fClonesArray;  

  SetFileNameTrackPoints(""); 
  SetFileNameGeometry(""); 

}

//____________________________________________________________________________
void AliITSResidualsAnalysis::InitHistograms(const TArrayI *volIDs)
{
  //
  // Method that sets and creates the required hisstrograms
  // with the correct binning (it dos not fill them)
  //

  TString histnameRPHI="HistRPHI_volID_",aux;
  TString histnameZ="HistZ_volID_";
  TString histnameGlob="HistGlob_volID_";
  TString histnameCorrVol="HistCorrVol_volID";
  TString histnamePullRPHI="HistPullRPHI_volID_";
  TString histnamePullZ="HistPullZ_volID_";

  TString histnameDirPhi="HistTrackDirPhi_volID_";
  TString histnameDirLambda="HistTrackDirLambda_volID_";
  TString histnameDirLambda2="HistTrackDirLambda2_volID_";
  TString histnameDirAlpha="HistTrackDirAlpha_volID_";
  TString histnameDir="HistTrackDir_volID_";


  fnHist=volIDs->GetSize();
  fVolResHistRPHI=new TH1F*[fnHist];
  fResHistGlob=new TH1F*[fnHist];
  fResHistZ=new TH1F*[fnHist];
  fPullHistRPHI=new TH1F*[fnHist];
  fPullHistZ=new TH1F*[fnHist];
  fhistCorrVol=new TH2F*[fnHist];
 

  fTrackDirPhi=new TH1F*[fnHist];
  fTrackDirLambda=new TH1F*[fnHist];
  fTrackDirLambda2=new TH1F*[fnHist];
  fTrackDirAlpha=new TH1F*[fnHist];

	
  fTrackDirPhiAll=new TH1F("fTrackDirPhiAll","fTrackDirPhiAll",100,-180,180);
  fTrackDirLambdaAll=new TH1F("fTrackDirLambdaAll","fTrackDirLambdaAll",100,-180,180);
  fTrackDirLambda2All=new TH1F("fTrackDirLambda2All","fTrackDirLambda2All",100,0,180);
  fTrackDirAlphaAll=new TH1F("fTrackDirAlphaAll","fTrackDirAlphaAll",100,-180,180);

  fTrackDirAll=new TH2F("fTrackDirAll","Hist with trakcs directions",100,-180,180,100,-180,180);
  fTrackDir2All=new TH2F("fTrackDir2All","Hist with trakcs directions",100,-180,180,100,-180,180);
 fTrackDirXZAll=new TH2F("fTrackDirXZAll","Hist with trakcs directions from XZ ",100,-3.,3.,100,-3.,3.);

  fTrackDir=new TH2F*[fnHist];

  Float_t **binningZPhi=CheckSingleLayer(volIDs);
  fvolidsToBin=new Int_t*[fnPhi*fnZ];

  Float_t *binningphi=binningZPhi[0];
  Float_t *binningz=binningZPhi[1];
  Bool_t binning=SetBinning(volIDs,binningphi,binningz);

  if(binning){
    fVolNTracks=new TH2F("fVolNTracks","Hist with N tracks passing through a given module==(r,phi) zone",fnPhi,binningphi,fnZ,binningz);
    fhistVolNptsUsed=new TH2F("fhistVolNptsUsed","Hist with N points used for given module==(r,phi) ",fnPhi,binningphi,fnZ,binningz);
    fhistVolUsed=new TH2F("fhistVolUsed","Hist with N modules used for a given module==(r,phi) zone",fnPhi,binningphi,fnZ,binningz);
    fSigmaVolZ=new TH2F("fSigmaVolZ","Hist with Sigma of Residual distribution for each module",fnPhi,binningphi,fnZ,binningz);
    fhEmpty=new TH2F("fhEmpty","Hist for getting binning",fnPhi,binningphi,fnZ,binningz);
    fVolNTracks->SetXTitle("Volume #phi");
    fVolNTracks->SetYTitle("Volume z ");
    fhistVolNptsUsed->SetXTitle("Volume #phi");
    fhistVolNptsUsed->SetYTitle("Volume z ");
    fhistVolUsed->SetXTitle("Volume #phi");
    fhistVolUsed->SetYTitle("Volume z ");
    fSigmaVolZ->SetXTitle("Volume #phi");
    fSigmaVolZ->SetYTitle("Volume z ");
  }
  else{
    fVolNTracks=new TH2F("fVolNTracks","Hist with N tracks passing through a given module==(r,phi) zone",50,-3.2,3.2,100,-80,80);
    fhistVolNptsUsed=new TH2F("fhistVolNptsUsed","Hist with N points used for given module==(r,phi) ",50,-3.2,3.2,100,-80,80);
    fhistVolUsed=new TH2F("fhistVolUsed","Hist with N modules used for a given module==(r,phi) zone",50,-3.2,3.2,100,-80,80);
    fSigmaVolZ=new TH2F("fSigmaVolZ","Hist with Sigma of Residual distribution for each module",50,-3.2,3.2,100,-80,80);
    fhEmpty=new TH2F("fhEmpty","Hist for getting binning",50,-3.2,3.2,100,-80,80);
    fVolNTracks->SetXTitle("Volume #phi");
    fVolNTracks->SetYTitle("Volume z ");
    fhistVolNptsUsed->SetXTitle("Volume #phi");
    fhistVolNptsUsed->SetYTitle("Volume z ");
    fhistVolUsed->SetXTitle("Volume #phi");
    fhistVolUsed->SetYTitle("Volume z ");
    fSigmaVolZ->SetXTitle("Volume #phi");
    fSigmaVolZ->SetYTitle("Volume z ");
  }
  fpTrackVolIDs=new TArrayI(fnHist);
  fVolUsed=new TArrayI*[fnHist];
  fVolVolids=new TArrayI*[fnHist]; 
  fLastVolVolid=new Int_t[fnHist];

  for (Int_t nhist=0;nhist<fnHist;nhist++){
    fpTrackVolIDs->AddAt(volIDs->At(nhist),nhist);   
    aux=histnameRPHI;
    aux+=volIDs->At(nhist);
    fVolResHistRPHI[nhist]=new TH1F("histname","histname",200,-0.02,0.02);   
    fVolResHistRPHI[nhist]->SetName(aux.Data()); 
    fVolResHistRPHI[nhist]->SetTitle(aux.Data()); 
    
    aux=histnameZ;
    aux+=volIDs->At(nhist);
    fResHistZ[nhist]=new TH1F("histname","histname",400,-0.08,0.08);   
    fResHistZ[nhist]->SetName(aux.Data()); 
    fResHistZ[nhist]->SetTitle(aux.Data()); 

    aux=histnamePullRPHI;
    aux+=volIDs->At(nhist);
    fPullHistRPHI[nhist]=new TH1F("histname","histname",100,-7.,7.);   
    fPullHistRPHI[nhist]->SetName(aux.Data()); 
    fPullHistRPHI[nhist]->SetTitle(aux.Data()); 
    
    aux=histnamePullZ;
    aux+=volIDs->At(nhist);
    fPullHistZ[nhist]=new TH1F("histname","histname",100,-7.,7.);   
    fPullHistZ[nhist]->SetName(aux.Data()); 
    fPullHistZ[nhist]->SetTitle(aux.Data()); 

    aux=histnameDirPhi;
    aux+=volIDs->At(nhist);
    fTrackDirPhi[nhist]=new TH1F("histname","histname",100,-180,180);
    fTrackDirPhi[nhist]->SetName(aux.Data()); 
    fTrackDirPhi[nhist]->SetTitle(aux.Data()); 

    aux=histnameDirLambda;
    aux+=volIDs->At(nhist);
    fTrackDirLambda[nhist]=new TH1F("histname","histname",100,0,180);
    fTrackDirLambda[nhist]->SetName(aux.Data()); 
    fTrackDirLambda[nhist]->SetTitle(aux.Data()); 

    aux=histnameDirLambda2;
    aux+=volIDs->At(nhist);
    fTrackDirLambda2[nhist]=new TH1F("histname","histname",100,0,180);
    fTrackDirLambda2[nhist]->SetName(aux.Data()); 
    fTrackDirLambda2[nhist]->SetTitle(aux.Data()); 
    
    aux=histnameDirAlpha;
    aux+=volIDs->At(nhist);
    fTrackDirAlpha[nhist]=new TH1F("histname","histname",100,-180,180);
    fTrackDirAlpha[nhist]->SetName(aux.Data()); 
    fTrackDirAlpha[nhist]->SetTitle(aux.Data()); 

    aux=histnameDir;
    aux+=volIDs->At(nhist);
    fTrackDir[nhist]=new TH2F("histname","histname",100,-90.,90.,100,-180.,180.);    
    fTrackDir[nhist]->SetName(aux.Data()); 
    fTrackDir[nhist]->SetTitle(aux.Data()); 

    aux=histnameGlob;
    aux+=volIDs->At(nhist);
    fResHistGlob[nhist]=new TH1F("histname","histname",400,-0.08,0.08);   
    fResHistGlob[nhist]->SetName(aux.Data());
    fResHistGlob[nhist]->SetTitle(aux.Data());

    aux=histnameCorrVol;
    aux+=volIDs->At(nhist);
    fhistCorrVol[nhist]=new TH2F("histname","histname",50,-3.2,3.2,100,-80,80);   

  
    fhistCorrVol[nhist]->SetName(aux.Data());
    fhistCorrVol[nhist]->SetTitle(aux.Data()); 
    fhistCorrVol[nhist]->SetXTitle("Volume #varphi");
    fhistCorrVol[nhist]->SetYTitle("Volume z ");
    fVolVolids[nhist]=new TArrayI(100);
    fVolUsed[nhist]=new TArrayI(1000);  
    fLastVolVolid[nhist]=0;
 
  }
  fWriteHist=kFALSE;

  return;
}

//____________________________________________________________________________
void AliITSResidualsAnalysis::ListVolUsed(TTree *pointsTree,TArrayI ***arrayIndex,Int_t **lastIndex)
{
  //
  // This is copied from AliAlignmentClass::LoadPoints() method
  //

  Int_t volIDalignable,volIDpoint,iModule; 
  AliTrackPoint p;
  AliTrackPointArray* array = 0;
  pointsTree->SetBranchAddress("SP", &array);
  
  
  for(Int_t ivol=0;ivol<fnHist;ivol++){
    Int_t lastused=0;
    volIDalignable=fpTrackVolIDs->At(ivol);
    AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer((UShort_t)volIDalignable,iModule);
    
    Int_t nArraysId = lastIndex[iLayer-AliGeomManager::kFirstLayer][iModule];
    printf("volume %d (Layer %d, Modulo %d) , numero di elementi per volume %d \n",volIDalignable,iLayer,iModule,nArraysId);
    TArrayI *index = arrayIndex[iLayer-AliGeomManager::kFirstLayer][iModule];
    for (Int_t iArrayId = 0;iArrayId < nArraysId; iArrayId++) {

      // Get tree entry
      Int_t entry = (*index)[iArrayId];

      pointsTree->GetEvent(entry);
      if (!array) {
	AliWarning("Wrong space point array index!");
	continue;
      }
      
      // Get the space-point array
      Int_t modnum,nPoints = array->GetNPoints();
  
      for (Int_t iPoint = 0; iPoint < nPoints; iPoint++) {
	array->GetPoint(p,iPoint);
	
	AliGeomManager::ELayerID layer = AliGeomManager::VolUIDToLayer(p.GetVolumeID(),modnum);
	// check if the layer id is valid
	if ((layer < AliGeomManager::kFirstLayer) ||
	    (layer >= AliGeomManager::kLastLayer)) {
	  AliError(Form("Layer index is invalid: %d (%d -> %d) !",
			layer,AliGeomManager::kFirstLayer,AliGeomManager::kLastLayer-1));
	  continue;
	}
	if ((modnum >= AliGeomManager::LayerSize(layer)) ||
	    (modnum < 0)) {
	  AliError(Form("Module number inside layer %d is invalid: %d (0 -> %d)",
			layer,modnum,AliGeomManager::LayerSize(layer)));
	  continue;
	}
	if (layer > AliGeomManager::kSSD2) continue; // ITS only
	
	volIDpoint=(Int_t)p.GetVolumeID();
	if (volIDpoint==volIDalignable)continue;
	Int_t size = fVolVolids[ivol]->GetSize();
	// If needed allocate new size
	if (fLastVolVolid[ivol]>=size){// Warning: fLAST[NHIST] is useless
	  fVolVolids[ivol]->Set(size + 1000);
	}
     	fVolVolids[ivol]->AddAt(volIDpoint,fLastVolVolid[ivol]);
	fLastVolVolid[ivol]++;
	Bool_t usedVol=kFALSE;
	for(Int_t used=0;used<lastused;used++){
	  if(fVolUsed[ivol]->At(used)==volIDpoint){
	    usedVol=kTRUE;
	    break;
	  }
	}
	if (!usedVol){
	  size = fVolUsed[ivol]->GetSize();
	  // If needed allocate new size
	  if (lastused>= size){
	    fVolUsed[ivol]->Set(size + 1000);
	  }
	  fVolUsed[ivol]->AddAt(volIDpoint,lastused);
	  lastused++;
	}
	
	FillVolumeCorrelationHists(ivol,volIDalignable,volIDpoint,usedVol);	
      }
    }
  }
  fWriteHist=kTRUE;
  return;
}

//____________________________________________________________________________
void AliITSResidualsAnalysis::FillVolumeCorrelationHists(Int_t ivol,Int_t volIDalignable,Int_t volIDpoint,Bool_t usedVol) const
{
  //
  // Fill the histograms with the correlations between volumes
  //
  
  if(!gGeoManager)AliGeomManager::LoadGeometry(GetFileNameGeometry());
  Double_t *transGlobal,radius,phi;
  const char *symname,*volpath;
  TGeoPNEntry *pne;
  TGeoPhysicalNode *pn;
  TGeoHMatrix *globMatrix;  
  
  symname = AliGeomManager::SymName(volIDalignable);
  pne = gGeoManager->GetAlignableEntry(symname);
  volpath=pne->GetTitle();
  pn=gGeoManager->MakePhysicalNode(volpath);
  globMatrix=pn->GetMatrix();
  
  transGlobal=globMatrix->GetTranslation();
  radius=TMath::Sqrt(transGlobal[0]*transGlobal[0]+transGlobal[1]*transGlobal[1]);
  phi=TMath::ATan2(transGlobal[1],transGlobal[0]);
  fhistVolNptsUsed->Fill(phi,transGlobal[2]);
  if(!usedVol)fhistVolUsed->Fill(phi,transGlobal[2]);

  symname = AliGeomManager::SymName(volIDpoint);
  pne = gGeoManager->GetAlignableEntry(symname);
  volpath=pne->GetTitle();
  pn=gGeoManager->MakePhysicalNode(volpath);
  globMatrix=pn->GetMatrix();
  
  transGlobal=globMatrix->GetTranslation();
  radius=TMath::Sqrt(transGlobal[0]*transGlobal[0]+transGlobal[1]*transGlobal[1]);
  phi=TMath::ATan2(transGlobal[1],transGlobal[0]);
  
  fhistCorrVol[ivol]->Fill(phi,transGlobal[2]);

  return;
}
  
   
//____________________________________________________________________________
void AliITSResidualsAnalysis::FillResHists(AliTrackPointArray *points,AliTrackPointArray *pTrack) const
{
  //
  // Method that fills the histograms with the residuals
  //
  
  Int_t volIDpoint;  
  Float_t xyz[3],xyz2[3];
  const Float_t *cov,*cov2;
  Float_t resRPHI,resGlob,resZ;
  Double_t pullz,pullrphi,sign;
  Double_t phi,lambda,lambda2,alpha,xovery,zovery;
  AliTrackPoint p,pTr;
  for(Int_t ipoint=0;ipoint<points->GetNPoints();ipoint++){
    points->GetPoint(p,ipoint);
    volIDpoint=(Int_t)p.GetVolumeID();
    p.GetXYZ(xyz);
    cov=p.GetCov();
    pTrack->GetPoint(pTr,ipoint);
    GetTrackDirClusterCov(&pTr,phi,lambda,lambda2,alpha,xovery,zovery);
    pTr.GetXYZ(xyz2);
    cov2=pTr.GetCov();
    resRPHI=TMath::Sqrt((xyz2[0]-xyz[0])*(xyz2[0]-xyz[0])+(xyz2[1]-xyz[1])*(xyz2[1]-xyz[1]));
    //resRPHI is always positive value
    sign=TMath::ATan2(xyz2[1],xyz2[0])-TMath::ATan2(xyz[1],xyz[0]);
    if(sign!=0.){
      sign=sign/TMath::Abs(sign);
      resRPHI=resRPHI*sign;
      pullrphi=sign*resRPHI*resRPHI/TMath::Sqrt((xyz2[0]-xyz[0])*(xyz2[0]-xyz[0])*(cov2[0]/100000000.+cov[0])+(xyz2[1]-xyz[1])*(xyz2[1]-xyz[1])*(cov2[3]/100000000.+cov[3]));
    }
    else{
      pullrphi=0.;
      resRPHI=0.;
    }
    
    resZ=xyz2[2]-xyz[2];
    pullz=resZ/(TMath::Sqrt(cov2[5])/10000.);
    resGlob=TMath::Sqrt((xyz2[0]-xyz[0])*(xyz2[0]-xyz[0])+(xyz2[1]-xyz[1])*(xyz2[1]-xyz[1])+(xyz2[2]-xyz[2])*(xyz2[2]-xyz[2]));
    for(Int_t ivolIDs=0;ivolIDs<fpTrackVolIDs->GetSize();ivolIDs++){
      if(volIDpoint==fpTrackVolIDs->At(ivolIDs)){
	fVolResHistRPHI[ivolIDs]->Fill(resRPHI);
	fResHistZ[ivolIDs]->Fill(resZ);
	fResHistGlob[ivolIDs]->Fill(resGlob);


	fTrackDirPhi[ivolIDs]->Fill(phi);
	fTrackDirLambda[ivolIDs]->Fill(lambda);
	fTrackDirLambda2[ivolIDs]->Fill(lambda2);
	fTrackDirAlpha[ivolIDs]->Fill(alpha);
	
	fTrackDirPhiAll->Fill(phi);
	fTrackDirLambdaAll->Fill(lambda);
	fTrackDirLambda2All->Fill(lambda2);
	fTrackDirAlphaAll->Fill(alpha);

	fTrackDirAll->Fill(lambda,alpha);
	fTrackDir2All->Fill(lambda2,phi);
	fTrackDirXZAll->Fill(xovery,zovery);
	fTrackDir[ivolIDs]->Fill(lambda,alpha);

	fPullHistRPHI[ivolIDs]->Fill(pullrphi);
	fPullHistZ[ivolIDs]->Fill(pullz);
	
	if(fsingleLayer){
	  Int_t binz,binphi;
	  Float_t globalPhi,globalZ;
	  if(kTRUE||(fvolidsToBin[ivolIDs][0]!=volIDpoint)){
	    binphi=GetBinPhiZ((Int_t)volIDpoint,&binz);
	  }
	  else{//this in the case of alignment of one entire layer (fnHIst=layersize) may reduce iterations: remind of that fsingleLayer->fnHista<layerSize
	    binphi=fvolidsToBin[ivolIDs][1];
	    binz=fvolidsToBin[ivolIDs][2];
	  }
	  globalPhi=fCoordToBinTable[binphi][binz][0];
	  globalZ=fCoordToBinTable[binphi][binz][1];
	  
	  fVolNTracks->Fill(globalPhi,globalZ);
	}
	else fVolNTracks->Fill(TMath::ATan2(xyz[1],xyz[0]),xyz[2]);
      }
    }
  }
}


//____________________________________________________________________________
Bool_t AliITSResidualsAnalysis::AnalyzeHists(Int_t minNpoints) const
{
  //  
  // Saves the histograms into a tree and saves the tree into a file
  //

  TString outname = "ResidualsAnalysisTree.root";
  TFile *hFile=new TFile(outname.Data(),"RECREATE","The Files containing the TREE with Align. Vol. hists nd Prop.");
  TTree *analysisTree=new TTree("analysisTree","Tree whith residuals analysis data for alignable volumes");
  static histProperties_t histRPHIprop,histZprop,histGlobprop;
  Int_t volID;

  TF1 *gauss;
  TH1F *histRPHI,*histZ,*histGlob,*histPullRPHI,*histPullZ,*hTrackDirPhi,*hTrackDirLambda,*hTrackDirLambda2,*hTrackDirAlpha;

  TH2F *histCorrVol,*hTrackDir;

  histRPHI=new TH1F();
  histZ=new TH1F();
  histPullRPHI=new TH1F();
  histPullZ=new TH1F();
  hTrackDirPhi=new TH1F();
  hTrackDirLambda=new TH1F();
  hTrackDirLambda2=new TH1F();
  hTrackDirAlpha=new TH1F();
  hTrackDir=new TH2F();
  histGlob=new TH1F();
  histCorrVol=new TH2F();
  Float_t globalPhi,globalZ;
  Double_t rms;
  Int_t nHistAnalyzed=0,entries;
  analysisTree->Branch("volID",&volID,"volID/I");
  if(fsingleLayer){
    analysisTree->Branch("globalPhi",&globalPhi,"globalPhi/F");
    analysisTree->Branch("globalZ",&globalZ,"globalZ/F");
  }
  analysisTree->Branch("histRPHI","TH1F",&histRPHI,128000,0);
  analysisTree->Branch("histPullRPHI","TH1F",&histPullRPHI,128000,0);
  
  analysisTree->Branch("histRPHIprop",&histRPHIprop,"nentries/I:rms/F:meanFit:errmeanFit:sigmaFit");
  analysisTree->Branch("histZ","TH1F",&histZ,128000,0);
  analysisTree->Branch("histPullZ","TH1F",&histPullZ,128000,0);
  analysisTree->Branch("hTrackDirPhi","TH1F",&hTrackDirPhi,128000,0);
  analysisTree->Branch("hTrackDirLambda","TH1F",&hTrackDirLambda,128000,0);
  analysisTree->Branch("hTrackDirLambda2","TH1F",&hTrackDirLambda2,128000,0);
  analysisTree->Branch("hTrackDirAlpha","TH1F",&hTrackDirAlpha,128000,0);
  analysisTree->Branch("hTrackDir","TH2F",&hTrackDir,128000,0);

  analysisTree->Branch("histZprop",&histZprop,"nentries/I:rms/F:meanFit:errmeanFit:sigmaFit");
  analysisTree->Branch("histGlob","TH1F",&histGlob,128000,0);
  analysisTree->Branch("histGlobprop",&histGlobprop,"nentries/I:rms/F:meanFit:errmeanFit:sigmaFit");
  if(fWriteHist){
    analysisTree->Branch("histCorrVol","TH2F",&histCorrVol,128000,0);
  }
  
  for(Int_t j=0;j<fnHist;j++){
    volID=fpTrackVolIDs->At(j);
    if((entries=(fResHistGlob[j]->GetEntries())>=minNpoints)||fsingleLayer){
      nHistAnalyzed++;
      //histRPHI
      histRPHI=fVolResHistRPHI[j];
      histPullRPHI=fPullHistRPHI[j];
      histRPHIprop.nentries=(Int_t)fVolResHistRPHI[j]->GetEntries();
      rms=fVolResHistRPHI[j]->GetRMS();
      histRPHIprop.rms=rms;
      gauss=new TF1("gauss","gaus",-3*rms,3*rms);
      fVolResHistRPHI[j]->Fit("gauss","RN");
      histRPHIprop.meanFit=gauss->GetParameter(1);
      histRPHIprop.errmeanFit=gauss->GetParError(1);
      histRPHIprop.sigmaFit=gauss->GetParameter(2);     
      //histZ
      histZ=fResHistZ[j];
      histPullZ=fPullHistZ[j];
      histZprop.nentries=(Int_t)fResHistZ[j]->GetEntries();
      rms=fResHistZ[j]->GetRMS();
      histZprop.rms=rms;
      gauss=new TF1("gauss","gaus",-3*rms,3*rms);
      fResHistZ[j]->Fit("gauss","RN");
      histZprop.meanFit=gauss->GetParameter(1);
      histZprop.errmeanFit=gauss->GetParError(1);
      histZprop.sigmaFit=gauss->GetParameter(2);
      //histGlob
      histGlob=fResHistGlob[j];
      histGlobprop.nentries=(Int_t)fResHistGlob[j]->GetEntries();
      rms=fResHistGlob[j]->GetRMS();
      histGlobprop.rms=rms;
      gauss=new TF1("gauss","gaus",-3*rms,3*rms);
      fResHistGlob[j]->Fit("gauss","RN");
      histGlobprop.meanFit=gauss->GetParameter(1);
      histGlobprop.errmeanFit=gauss->GetParError(1);
      histGlobprop.sigmaFit=gauss->GetParameter(2);

      //histTrackDir
      hTrackDirPhi=fTrackDirPhi[j];
      hTrackDirLambda=fTrackDirLambda[j];
      hTrackDirLambda2=fTrackDirLambda2[j];
      hTrackDirAlpha=fTrackDirAlpha[j];
      hTrackDir=fTrackDir[j];
      
      if(fsingleLayer){
	Int_t binz,binphi;
	if (fvolidsToBin[j][0]!=volID)binphi=GetBinPhiZ((Int_t)volID,&binz);
	else{
	  binphi=fvolidsToBin[j][1];
	  binz=fvolidsToBin[j][2];
	}
	globalPhi=fCoordToBinTable[binphi][binz][0];
	globalZ=fCoordToBinTable[binphi][binz][1];
	
	
	histCorrVol=fhistCorrVol[j];
	fSigmaVolZ->SetBinContent(binphi+1,binz+1,histRPHIprop.sigmaFit);//+1 is for underflow
      }
      analysisTree->Fill();
    }
    else continue;
    
  }
  if(nHistAnalyzed>0){ 
    analysisTree->Print();
    fVolNTracks->Write();
    hFile->Write();
    fhEmpty->Write();
    if(fWriteHist){
      fhistVolUsed->Draw();
      fSigmaVolZ->Draw();
      fSigmaVolZ->Write();
      fhistVolUsed->Write();
      fTrackDirPhiAll->Write();
      fTrackDirLambdaAll->Write();
      fTrackDirLambda2All->Write();
      fTrackDirAlphaAll->Write();
      fTrackDirAll->Write();
      fTrackDir2All->Write();
      fTrackDirXZAll->Write();
      fhistVolNptsUsed->Write();
      hFile->Close();
    }
    return kTRUE;
  }
  else {
    delete analysisTree;
    delete hFile;
    return kFALSE;}
}


//____________________________________________________________________________
void AliITSResidualsAnalysis::DrawHists() const
{
  //
  // Draws the histograms of the residuals and of the number of tracks
  //

  TString cname;
  for(Int_t canv=0;canv<fnHist;canv++){
    cname="canv Residuals";
    cname+=canv;
    TCanvas *c=new TCanvas(cname.Data(),cname.Data(),700,700);
    c->Divide(3,1);
    c->cd(1);
    fVolResHistRPHI[canv]->Draw();
    c->cd(2);
    fResHistZ[canv]->Draw();
    c->cd(3);
    fResHistGlob[canv]->Draw();
  }
  cname="canv NVolTracks";
  
  TCanvas *c2=new TCanvas(cname.Data(),cname.Data(),700,700);
  c2->cd();
  fVolNTracks->Draw();  
  
  return;
}

//____________________________________________________________________________
Float_t** AliITSResidualsAnalysis::CheckSingleLayer(const TArrayI *volids)
{
  //
  // Checks if volumes array is a single (ITS) layer or not
  //
  
  Float_t **binningzphi=new Float_t*[2];
  Int_t iModule;
  AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer((UShort_t)volids->At(0),iModule);
  //Check that one single Layer is going to be aligned
  for(Int_t nvol=0;nvol<volids->GetSize();nvol++){
    if(iLayer != AliGeomManager::VolUIDToLayer((UShort_t)volids->At(nvol),iModule)){
      printf("Wrong Layer! \n %d , %d , %d ,%d \n",(Int_t)AliGeomManager::VolUIDToLayer((UShort_t)volids->At(nvol),iModule),nvol,volids->GetSize(),iModule);
      fsingleLayer=kFALSE;
      return 0;
    }
  }
  
  //Bool_t used=kFALSE;
  switch (iLayer) {
  case AliGeomManager::kSPD1:{
    fnPhi=kPhiSPD1;//kPhiSPD1;
    fnZ=kZSPD1;//nZSPD1;
    binningzphi[0]=new Float_t[kPhiSPD1+1];
    binningzphi[1]=new Float_t[kZSPD1+1];
    fCoordToBinTable=new Double_t**[kPhiSPD1];
    for(Int_t j=0;j<fnPhi;j++){
      fCoordToBinTable[j]=new Double_t*[kZSPD1];
    }
  }; break;
  case AliGeomManager::kSPD2:{
    fnPhi=kPhiSPD2;//kPhiSPD1;
    fnZ=kZSPD2;//nZSPD1;
    binningzphi[0]=new Float_t[kPhiSPD2+1];
    binningzphi[1]=new Float_t[kZSPD2+1];
    fCoordToBinTable=new Double_t**[kPhiSPD2];
    for(Int_t j=0;j<fnPhi;j++){
      fCoordToBinTable[j]=new Double_t*[kZSPD2];
    }

  }; break; case AliGeomManager::kSDD1:{
    fnPhi=kPhiSDD1;//kPhiSPD1;
    fnZ=kZSDD1;//nZSPD1;
    binningzphi[0]=new Float_t[kPhiSDD1+1];
    binningzphi[1]=new Float_t[kZSDD1+1];
    fCoordToBinTable=new Double_t**[kPhiSDD1];
    for(Int_t j=0;j<fnPhi;j++){
      fCoordToBinTable[j]=new Double_t*[kZSDD1];
    }
  }; break; case AliGeomManager::kSDD2:{
    fnPhi=kPhiSDD2;//kPhiSPD1;
    fnZ=kZSDD2;//nZSPD1;
    binningzphi[0]=new Float_t[kPhiSDD2+1];
    binningzphi[1]=new Float_t[kZSDD2+1];
    fCoordToBinTable=new Double_t**[kPhiSDD2];
    for(Int_t j=0;j<fnPhi;j++){
      fCoordToBinTable[j]=new Double_t*[kZSDD2];
    }
  }; break; case AliGeomManager::kSSD1:{
    fnPhi=kPhiSSD1;//kPhiSPD1;
    fnZ=kZSSD1;//nZSPD1;
    binningzphi[0]=new Float_t[kPhiSSD1+1];
    binningzphi[1]=new Float_t[kZSSD1+1];
    fCoordToBinTable=new Double_t**[kPhiSSD1];
    for(Int_t j=0;j<fnPhi;j++){
      fCoordToBinTable[j]=new Double_t*[kZSSD1];
    }
  }; break; case AliGeomManager::kSSD2:{
    fnPhi=kPhiSSD2;//kPhiSPD1;
    fnZ=kZSSD2;//nZSPD1;
    binningzphi[0]=new Float_t[kPhiSSD2+1];
    binningzphi[1]=new Float_t[kZSSD2+1];
    fCoordToBinTable=new Double_t**[kPhiSSD2];
    for(Int_t j=0;j<fnPhi;j++){
      fCoordToBinTable[j]=new Double_t*[kZSSD2];
    }
  }; break;
  
  default:{
    printf("Wrong Layer Label! \n");    
    fsingleLayer=kFALSE;
    return 0x0;
  } 
  }
  fsingleLayer=kTRUE;
  return binningzphi;
}

//____________________________________________________________________________
Bool_t AliITSResidualsAnalysis::SetBinning(const TArrayI *volids,Float_t *phiBin,Float_t *zBin)
{
  //
  // Sets the correct binning
  //
  
  if(!fsingleLayer)return kFALSE;
  const char *volpath,*symname;
  Int_t iModule;
  Int_t *orderArrayPhi,*orderArrayZ;
  UShort_t volID;
  Double_t *phiArray,*zArray,*transGlobal,*phiArrayOrdered,*zArrayOrdered; 
  Double_t lastPhimin=-10;
  Double_t lastZmin=-99;
  Int_t ***orderPhiZ;
  TGeoPNEntry *pne;
  TGeoPhysicalNode *pn;
  TGeoHMatrix *globMatrix;
  
  Bool_t used=kFALSE;
  
  orderPhiZ=new Int_t**[fnPhi];
  phiArray=new Double_t[fnPhi];//phiBin[nModulesPhi+1];
  zArray=new Double_t[fnZ];//zBin[nModulesZ+1];
  phiArrayOrdered=new Double_t[fnPhi];
  zArrayOrdered=new Double_t[fnZ];
  orderArrayPhi=new Int_t[fnPhi];
  orderArrayZ=new Int_t[fnZ];
  for(Int_t k=0;k<fnZ;k++){
    orderArrayZ[k]=0;
    zArray[k]=-1000;
  }
  for(Int_t k=0;k<fnPhi;k++){
    orderArrayPhi[k]=0;
    phiArray[k]=-10;
    orderPhiZ[k]=new Int_t*[fnZ];
  }
  
  
  AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer((UShort_t)volids->At(0),iModule);  
  
  Int_t lastPhi=0,lastZ=0;
  for(iModule=0;iModule<AliGeomManager::LayerSize(iLayer);iModule++){
    fvolidsToBin[iModule]=new Int_t[3];
    volID=AliGeomManager::LayerToVolUID(iLayer,iModule);
    fvolidsToBin[iModule][0]=volID;
    symname = AliGeomManager::SymName(volID);
    pne = gGeoManager->GetAlignableEntry(symname);
    volpath=pne->GetTitle();
    pn=gGeoManager->MakePhysicalNode(volpath);
    globMatrix=pn->GetMatrix();
    transGlobal=globMatrix->GetTranslation();
    
    for(Int_t j=0;j<lastPhi;j++){
      used=kFALSE;
      if(TMath::Abs(phiArray[j]-TMath::ATan2(transGlobal[1],transGlobal[0]))<2*TMath::Pi()/(10*fnPhi)){//10 is a safety factor but....
	fvolidsToBin[iModule][1]=j;
	used=kTRUE;
	break;
      }
    }
    if(!used){
      phiArray[lastPhi]=TMath::ATan2(transGlobal[1],transGlobal[0]);
      fvolidsToBin[iModule][1]=lastPhi;
      if(phiArray[lastPhi]<lastPhimin)lastPhimin=phiArray[lastPhi];
      lastPhi++;
      if(lastPhi>fnPhi){
	printf("Wrong Phi! \n");
	return kFALSE;}
    }
    
    for(Int_t j=0;j<lastZ;j++){
      used=kFALSE;
      if(TMath::Abs(zArray[j]-transGlobal[2])<0.1){
	fvolidsToBin[iModule][2]=j;
	used=kTRUE;
	break;
      }
    }
    if(!used){
      fvolidsToBin[iModule][2]=lastZ;
      zArray[lastZ]=transGlobal[2];
      if(zArray[lastZ]<lastZmin)lastZmin=zArray[lastZ];
      lastZ++;
      if(lastZ>fnZ){
	printf("Wrong Z! \n");
	return kFALSE;}
    }
  }
  
  
  //ORDERING THE ARRAY OF PHI AND Z VALUES
  for(Int_t order=0;order<fnPhi;order++){
    for(Int_t j=0;j<fnPhi;j++){
      if((j!=order)&&(phiArray[order]>phiArray[j]))orderArrayPhi[order]++;
    }
  }
  
  for(Int_t order=0;order<fnPhi;order++){
    for(Int_t j=0;j<fnPhi;j++){
      if(orderArrayPhi[j]==order){
	phiArrayOrdered[order]=phiArray[j];
	break;
	}
    }
  }
  
  
  for(Int_t order=0;order<fnZ;order++){
    for(Int_t j=0;j<fnZ;j++){
      if((j!=order)&&(zArray[order]>zArray[j]))orderArrayZ[order]++;
    }
  }
  
  
  for(Int_t order=0;order<fnZ;order++){
    for(Int_t j=0;j<fnZ;j++){
      if(orderArrayZ[j]==order){
	zArrayOrdered[order]=zArray[j];
	break;
      }
    }
  }

  
  //Filling the  fCoordToBinTable
  for(Int_t j=0;j<fnPhi;j++){
    for(Int_t i=0;i<fnZ;i++){
      orderPhiZ[j][i]=new Int_t[2];
      orderPhiZ[j][i][0]=orderArrayPhi[j];
      orderPhiZ[j][i][1]=orderArrayZ[i];
    }
  }
  
  
  for(Int_t j=0;j<fnPhi;j++){
    for(Int_t i=0;i<fnZ;i++){
      fCoordToBinTable[j][i]=new Double_t[2];
      fCoordToBinTable[j][i][0]=phiArrayOrdered[j];
      fCoordToBinTable[j][i][1]=zArrayOrdered[i];
      printf("ecco (binphi,binz)= %d, %d e (phi,z)=%f,%f \n",j,i,fCoordToBinTable[j][i][0],fCoordToBinTable[j][i][1]);
    }
  }
  Int_t istar,jstar;
  for(iModule=0;iModule<fnPhi*fnZ;iModule++){
    istar=fvolidsToBin[iModule][1];
    jstar=fvolidsToBin[iModule][2];
    fvolidsToBin[iModule][1]=orderPhiZ[istar][jstar][0];
    fvolidsToBin[iModule][2]=orderPhiZ[istar][jstar][1];
  }
  
    
  //now constructing the binning
  for(Int_t iModPhi=0;iModPhi<fnPhi-1;iModPhi++){
    phiBin[iModPhi+1]=(Float_t)phiArrayOrdered[iModPhi]+0.5*(phiArrayOrdered[iModPhi+1]-phiArrayOrdered[iModPhi]);
  }

  phiBin[0]=(Float_t)phiArrayOrdered[0]-(phiArrayOrdered[1]-phiArrayOrdered[0])/2.;

  phiBin[fnPhi]=(Float_t)phiArrayOrdered[fnPhi-1]+(phiArrayOrdered[fnPhi-1]-phiArrayOrdered[fnPhi-2])/2.;
  for(Int_t iModPhi=0;iModPhi<fnPhi+1;iModPhi++){
    printf("Mean Phi mod %d su %d:  %f \n",iModPhi+1,fnPhi,phiBin[iModPhi]);
  }
  
  for(Int_t iModZ=0;iModZ<fnZ-1;iModZ++){
    zBin[iModZ+1]=(Float_t)zArrayOrdered[iModZ]+0.5*(zArrayOrdered[iModZ+1]-zArrayOrdered[iModZ]);
  }
  zBin[0]=(Float_t)zArrayOrdered[0]-0.5*(zArrayOrdered[1]-zArrayOrdered[0]);
  zBin[fnZ]=(Float_t)zArrayOrdered[fnZ-1]+0.5*(zArrayOrdered[1]-zArrayOrdered[0]);
  
  
  for(Int_t iModPhi=0;iModPhi<fnZ+1;iModPhi++){
     printf("Mean Z mod %d su %d:  %f \n",iModPhi+1,fnZ,zBin[iModPhi]);
  }
  return kTRUE;
}


//____________________________________________________________________________
Int_t AliITSResidualsAnalysis::GetBinPhiZ(const Int_t volID,Int_t *binz) const
{
  //
  // Returns the correct Phi-Z bin
  //

  if (!fsingleLayer){
    printf("No Single Layer reAlignment! \n");
    return 100;
  }
  
  for(Int_t j=0;j<fnPhi*fnZ;j++){
    if(j==fnZ*fnPhi){
      printf("Wrong volume UID introduced! fnHist: %d  volID: %d iter: %d \n",fnHist,volID,j);
      return 100;
    }
    if(fvolidsToBin[j][0]==volID){
      
      *binz=fvolidsToBin[j][2];
      return fvolidsToBin[j][1];
    }
  }

  return 100;

}

//____________________________________________________________________________
TArrayI* AliITSResidualsAnalysis::GetSingleLayerVolids(Int_t layer) const
{
  //
  // Translates the layer number into a Volumes Array
  //

  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry());

  if(layer<1 || layer>6){
    printf("WRONG LAYER SET! \n");
    return 0;
  }
  Int_t iModule,size;
  UShort_t volid;
  size = AliGeomManager::LayerSize(layer);
  TArrayI *volIDs = new TArrayI(size);
  for(iModule=0;iModule<size;iModule++){
    volid = AliGeomManager::LayerToVolUID(layer,iModule);
    volIDs->AddAt(volid,iModule);
  }

  return volIDs;
  
}

//____________________________________________________________________________
void AliITSResidualsAnalysis::GetTrackDirClusterCov(AliTrackPoint *point,Double_t &phi,Double_t &lambda,Double_t &lambda2,Double_t &alpha,Double_t &xovery,Double_t &zovery) const
{
  //
  // ...
  //
  
  TMatrixDSym cov(3);
  const Float_t *covvector=point->GetCov();
  cov(0,0)=covvector[0];
  cov(1,0)=cov(0,1)=covvector[1];
  cov(2,0)=cov(0,2)=covvector[2];
  cov(1,1)=covvector[3];
  cov(1,2)=cov(2,1)=covvector[4];
  cov(2,2)=covvector[5];
  
  Double_t determinant=cov.Determinant();
  if(determinant!=0.){
    TMatrixD vect(3,3);
    TVectorD eigenvalues(3);
    const TMatrixDSymEigen keigen(cov);
    eigenvalues=keigen.GetEigenValues();
    vect=keigen.GetEigenVectors();
    Double_t mainvect[3];
    mainvect[0]=vect(0,0);
    mainvect[1]=vect(1,0);
    mainvect[2]=vect(2,0);
    if(mainvect[1]!=0.){
      xovery=mainvect[0]/mainvect[1];
      zovery=mainvect[2]/mainvect[1];
    }
    else {
      xovery=9999.;
      zovery=9999.;
    }
    if(mainvect[1]<0.){
      mainvect[0]=-1.*mainvect[0];
      mainvect[1]=-1.*mainvect[1];
      mainvect[2]=-1.*mainvect[2];
    }
    lambda2=TMath::ATan2(TMath::Sqrt(mainvect[0]*mainvect[0]+mainvect[2]*mainvect[2]),mainvect[1])*TMath::RadToDeg();
    lambda=TMath::ATan2(mainvect[2],TMath::Sqrt(mainvect[0]*mainvect[0]+mainvect[1]*mainvect[1]))*TMath::RadToDeg();
    phi=TMath::ATan2(mainvect[0],mainvect[2])*TMath::RadToDeg();
    alpha=TMath::ATan2(mainvect[1],mainvect[0])*TMath::RadToDeg();
  }
  else printf("determinant =0!, skip this point \n");
  
  return;
}

//____________________________________________________________________________
void AliITSResidualsAnalysis::CalculateResiduals(const TArrayI *volids, 
      const TArrayI *volidsfit,
      AliGeomManager::ELayerID layerRangeMin,
      AliGeomManager::ELayerID layerRangeMax,
      Int_t iterations,
      Bool_t draw)
{
  // CalculateResiduals for a set of detector volumes.
  // Tracks are fitted only within
  // the range defined by the user
  // (by layerRangeMin and layerRangeMax)
  // or within the set of volidsfit
  // Repeat the procedure 'iterations' times


  Int_t nVolIds = volids->GetSize();
  if (nVolIds == 0) {
    AliError("Volume IDs array is empty!");
    return;
  }

  // Load only the tracks with at least one
  // space point in the set of volume (volids)

  //AliAlignmentTracks::SetPointsFilename(GetFileNameTrackPoints()); 
  AliAlignmentTracks::BuildIndex();

  cout<<" Index Built!"<<endl;

  if(draw) ListVolUsed(fPointsTree,fArrayIndex,fLastIndex);
  
  AliTrackPointArray **points;

  // Start the iterations
  while (iterations > 0) {
    Int_t nArrays = LoadPoints(volids, points);
    if (nArrays == 0) return;
    
    AliTrackResiduals *minimizer = CreateMinimizer();
    minimizer->SetNTracks(nArrays);
    minimizer->InitAlignObj();
    AliTrackFitter *fitter = CreateFitter();
    
    for (Int_t iArray = 0; iArray < nArrays; iArray++) {
      if (!points[iArray]) continue;

     
      fitter->SetTrackPointArray(points[iArray],kFALSE);
      if (fitter->Fit(volids,volidsfit,layerRangeMin,layerRangeMax) == kFALSE) continue;
      AliTrackPointArray *pVolId,*pTrack;


      fitter->GetTrackResiduals(pVolId,pTrack);
      if(draw) FillResHists(pVolId,pTrack); // WARNING!

      minimizer->AddTrackPointArrays(pVolId,pTrack);
      
    }

    if(minimizer->GetNFilledTracks()<=1){
      printf("No good tracks found: could not find parameter for volume %d (and following in volids)\n",volids->At(0));
      UnloadPoints(nArrays, points);
      return;
    }

    minimizer->Minimize();
   
    // Update the alignment object(s)
    Int_t last=0;
 
    if(fRealignObjFileIsOpen)last=fClonesArray->GetLast(); 
    
    
    if (fDoUpdate) for (Int_t iVolId = 0; iVolId < nVolIds; iVolId++) {
      UShort_t volid = (*volids)[iVolId];
      Int_t iModule;
      AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer(volid,iModule);
      AliAlignObj *alignObj = fAlignObjs[iLayer-AliGeomManager::kFirstLayer][iModule];      
      *alignObj *= *minimizer->GetAlignObj();
      
      if(iterations==1)alignObj->Print("");
      if(iterations==1&&fRealignObjFileIsOpen){
	TClonesArray &alo=*fClonesArray;
	new(alo[last+1+(Int_t)iVolId])AliAlignObjParams(*alignObj);
      }
      
      
    }


    UnloadPoints(nArrays, points);
    
    iterations--;


    if(draw && iterations==0) AnalyzeHists(30);

  }

  return;

}


//______________________________________________________________________________
void AliITSResidualsAnalysis::ProcessPoints(TString minimizer,
      Int_t fit,
      AliGeomManager::ELayerID iLayerToAlign,
      AliGeomManager::ELayerID iLayerToExclude,
      TString misalignmentFile)
{
  //
  // This function process the AliTrackPoints (into residuals)
  //
 
  SetPointsFilename(GetFileNameTrackPoints());
  AliTrackFitter *fitter;

  if(fit==1){
    fitter = new AliTrackFitterKalman();
  }else fitter = new AliTrackFitterRieman();

  fitter->SetMinNPoints(4);
  SetTrackFitter(fitter);


  AliTrackResiduals *res;
  
  if(minimizer=="minuit"){
    res = new AliTrackResidualsChi2();
  }else if(minimizer=="minuitnorot"){
    res = new AliTrackResidualsChi2();
    res->FixParameter(3);
    res->FixParameter(4);
    res->FixParameter(5);
  }else if(minimizer=="fast"){
    res = new AliTrackResidualsFast();
  }else {
    printf("Trying to set a non existing minimizer! \n");
    return;
  }

  res->SetMinNPoints(1);
  SetMinimizer(res);
  Bool_t draw = kTRUE;

  if(misalignmentFile=="")printf("NO FAKE MISALIGNMENT\n");
  else {
    Bool_t misal=Misalign(misalignmentFile,"ITSAlignObjs");
    if(!misal)return;
  }
  
  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry());
  
  TArrayI volIDsFit(2200);
  Int_t iLayer,j=0;
  for (iLayer=(Int_t)AliGeomManager::kSPD1;iLayer<(Int_t)AliGeomManager::kTPC1;iLayer++){
    for (Int_t iModule=0;iModule<AliGeomManager::LayerSize(iLayer);iModule++){
      UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iModule);

      if((iLayer!=iLayerToAlign)&&(iLayer!=iLayerToExclude))volIDsFit.AddAt(volid,j);
      
      j++;
    }
  }
  
  Int_t size=AliGeomManager::LayerSize(iLayerToAlign);
  TArrayI volIDs(size);
  
  j=0;
  for (Int_t iModule=0;iModule<AliGeomManager::LayerSize(iLayerToAlign);iModule++){

    UShort_t volid = AliGeomManager::LayerToVolUID(iLayerToAlign,iModule);
    volIDs.AddAt(volid,j);
    j++;
  }
  
    if(j==0){printf("j=0 \n");return;}

    CalculateResiduals(&volIDs,&volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSSD2,1,draw);

  
    return;


}
 
//______________________________________________________________________________
void AliITSResidualsAnalysis::ExtractResiduals(Int_t layer,
					       Int_t minEnt,
					       TString filename) const
					   
{

  //
  // Function that saves the residuals into a Entuple
  //

  TString title,strminEnt=" (Npts > ";
  histProperties_t histRPHIprop,histZprop;

  // Labels for the plots
  strminEnt+=minEnt;
  strminEnt.Append(")");
  
  // name of the output file
  title="resPlot_MA_layer";
  title+=layer;
  title.Append(".root");
  
  // Load INfiles, OUTfiles and TTrees and labels them
  TFile *f1=TFile::Open(filename.Data());
  TFile &f2=*f1;
  TFile *outfile2=new TFile(title.Data(),"RECREATE");
  TFile &outfile=*outfile2;
  TTree *tRealign2=(TTree*)f2.Get("analysisTree"); // TTree with the DATA
  TTree &tRealign=*tRealign2;


  // Setting variables
  Int_t nEntries;
  Int_t *volid;
  Float_t z,phi;
  TH2F *hVolCorrBranch;
  TH1F *hResRPhi;
  TH1F *hResZ;
  
  TString layer2="";
  layer2+=layer;


  TH2F *hEmpty=(TH2F*)f2.Get("fhEmpty"); 
  hEmpty->SetName("hEmpty");


  // Creates histos using the "hEmpty" template (binning!)
  TH2F *hSigmaMeanRPHI=new TH2F();
  TH2F *hSigmaRPHI=new TH2F();
  TH2F *hSigmaMeanZ=new TH2F();
  hEmpty->Copy(*hSigmaMeanRPHI);
  hSigmaMeanRPHI->SetName("hSigmaMeanRPHI");
  hSigmaMeanRPHI->GetZaxis()->SetRangeUser(0.,200);
  hEmpty->Copy(*hSigmaRPHI);
  hSigmaRPHI->SetName("hSigmaRPHI");
  hSigmaRPHI->GetZaxis()->SetRangeUser(0.,200);
  hEmpty->Copy(*hSigmaMeanZ);
  hSigmaMeanZ->SetName("hSigmaMeanZ");
  hSigmaMeanZ->GetZaxis()->SetRangeUser(0.,400);


  // Branching of the DATA TTree
  tRealign.SetBranchAddress("globalPhi",&phi);
  tRealign.SetBranchAddress("globalZ",&z);
  tRealign.SetBranchAddress("histZ",&hResZ);
  tRealign.SetBranchAddress("histRPHI",&hResRPhi);
  tRealign.SetBranchAddress("volID",&volid);
  tRealign.SetBranchAddress("histCorrVol",&hVolCorrBranch);
  tRealign.SetBranchAddress("histRPHIprop",&histRPHIprop);  
  tRealign.SetBranchAddress("histZprop",&histZprop);  

  TNtuple *ntMonA = new TNtuple("ntMonA","Residuals","layer:volID:phi:z:nentries:meanFitRPHI:meanFitZ:RMS_RPHI");
  nEntries=tRealign.GetEntries();
  printf("entries: %d\n",nEntries);
  Float_t volidfill = 0;

  for(Int_t j=0;j<nEntries;j++){ // LOOP OVER THE ENTRIES

    printf(" Loading Event %d \n",j);

    tRealign.GetEvent(j);

    // Saving data in an entuple -> To be turned into a Tree
    ntMonA->Fill((Float_t)layer,
		volidfill, // WRONG! To be corrected!
		(Float_t)phi,
		(Float_t)z,
		10000*(Float_t)histRPHIprop.nentries,
		10000*(Float_t)histRPHIprop.meanFit,
		10000*(Float_t)histZprop.meanFit,
		10000*(Float_t)histRPHIprop.rms);

  } // END LOOP OVER ENTRIES
  

  //write histograms
  outfile.cd();//return to top directory
  ntMonA->Write();
  hEmpty->Write();

  delete  tRealign2;
  f2.Close();

  return;

}

//______________________________________________________________________________
Int_t AliITSResidualsAnalysis::PlotResiduals(Int_t layer,TString filename) const
{
  //
  // Function that plots the residual distributions
  //
  filename+=layer;
  filename.Append(".root");
  TFile *f1 = TFile::Open(filename.Data());
  if(!f1) return 1;

  TH2F *hEmpty=(TH2F*)f1->Get("hEmpty"); 
  TNtuple *ntData=(TNtuple*)f1->Get("ntMonA"); 
  if(!ntData) return 2;


  TH2F *hMeanZ = new TH2F("hMeanZ","Hist for getting banged",40,-20,20,30,-15,15);


  Int_t nn=4;
  Float_t x[4],y[4],ex[4],ey[4],yZ[4],eyZ[4];
  x[0]=10.5;
  x[1]=3.5;
  x[2]=-3.5;
  x[3]=-10.5;

  // Declaring Histos
  TH2F *hStatGlob = new TH2F();
  TH2F *hMeanGlob = new TH2F();

  TH1F **hMeanPHI;
  TH1F **hMeanPHIz;
  TH1F *hGlobPhi = new TH1F("hGlobPhi","hGlobPhi",41,-(TMath::Pi())-(TMath::Pi()/40),(TMath::Pi())+(TMath::Pi()/40));
  //TH1F *hGlobPhi = new TH1F("hGlobPhi","hGlobPhi",40,-(TMath::Pi()),(TMath::Pi()));

  hMeanPHI = new TH1F*[40]; //watch out!
  hMeanPHIz = new TH1F*[40];

  TString title;

  for(Int_t bPhi = 0; bPhi<40; bPhi++){
    title="hMeanPHI";
    title+=bPhi;
    hMeanPHI[bPhi]=new TH1F(title.Data(),title.Data(),300,-150,150);
    title="hMeanPHIz";
    title+=bPhi;
    hMeanPHIz[bPhi]=new TH1F(title.Data(),title.Data(),300,-150,150);
  }

  // Setting the binning of the histos
  hEmpty->Copy(*hStatGlob);
  hStatGlob->SetName("hStatGlob");
  hEmpty->Copy(*hMeanGlob);
  hMeanGlob->SetName("hMeanGlob");

  Int_t entries;
  Float_t volID,phi,z,nentries,meanFitRPHI,meanFitZ,rms;
  entries = (Int_t)ntData->GetEntries();

  // Branching ...
  //ntData->SetBranchAddress("layer",&layernt);
  ntData->SetBranchAddress("volID",&volID);
  ntData->SetBranchAddress("phi",&phi);
  ntData->SetBranchAddress("z",&z);
  ntData->SetBranchAddress("nentries",&nentries);
  ntData->SetBranchAddress("meanFitRPHI",&meanFitRPHI);
  ntData->SetBranchAddress("meanFitZ",&meanFitZ);
  ntData->SetBranchAddress("RMS_RPHI",&rms);

  Int_t nbytes = 0;
  Int_t bin,bin2,ban;
  Double_t c1,c2,c3,c4;
  Double_t m1,m2,m3,m4;
  Double_t n1,n2,n3,n4;
  c1=1e-10;
  c2=1e-10;
  c3=1e-10;
  c4=1e-10;

  TH1F *hMeanFit1 = new TH1F("hMeanFit1","lol",1000,-500,500);
  TH1F *hMeanFit2 = new TH1F("hMeanFit2","lol",1000,-500,500);
  TH1F *hMeanFit3 = new TH1F("hMeanFit3","lol",1000,-500,500);
  TH1F *hMeanFit4 = new TH1F("hMeanFit4","lol",1000,-500,500);

  TH1F *hMeanFitZ1 = new TH1F("hMeanFitZ1","lol",1000,-500,500);
  TH1F *hMeanFitZ2 = new TH1F("hMeanFitZ2","lol",1000,-500,500);
  TH1F *hMeanFitZ3 = new TH1F("hMeanFitZ3","lol",1000,-500,500);
  TH1F *hMeanFitZ4 = new TH1F("hMeanFitZ4","lol",1000,-500,500);

  for(Int_t j=0;j<entries;j++){

    nbytes += ntData->GetEvent(j);

    // Set binning
    bin=hStatGlob->FindBin(phi,z);
    bin2=hMeanZ->FindBin(meanFitRPHI,z);

    // Global Histograms
    hStatGlob->AddBinContent(bin,nentries);
    hMeanGlob->AddBinContent(bin,meanFitRPHI);
    hMeanZ->AddBinContent(bin2,1);

    bin=hGlobPhi->FindBin(phi);
    bin2=hMeanPHI[bin-2]->FindBin(meanFitRPHI);
    hMeanPHI[bin-2]->AddBinContent(bin2,1);
    bin2=hMeanPHIz[bin-2]->FindBin(meanFitZ);
    hMeanPHIz[bin-2]->AddBinContent(bin2,1);


    if(z<12 && z>9) {
      c1+=nentries;
      m1+=(meanFitRPHI*nentries);
      n1+=rms*nentries;
      ban=hMeanFit1->FindBin(meanFitRPHI);
      //hMeanFit1->AddBinContent(ban,1);
      hMeanFit1->AddBinContent(ban,1);
      ban=hMeanFitZ1->FindBin(meanFitZ);
      hMeanFitZ1->AddBinContent(ban,1);
    } else if(z<5 && z>2){
      c2+=nentries;
      m2+=(meanFitRPHI*nentries);
      n2+=rms*nentries;
      ban=hMeanFit2->FindBin(meanFitRPHI);
      //hMeanFit2->AddBinContent(ban,1);
      hMeanFit2->AddBinContent(ban,1);
      ban=hMeanFitZ2->FindBin(meanFitZ);
      hMeanFitZ2->AddBinContent(ban,1);
    } else if(z<-2 && z>-5){
      c3+=nentries;
      m3+=(meanFitRPHI*nentries);
      n3+=rms*nentries;
      ban=hMeanFit3->FindBin(meanFitRPHI);
      //hMeanFit3->AddBinContent(ban,1);
      hMeanFit3->AddBinContent(ban,1);
      ban=hMeanFitZ3->FindBin(meanFitZ);
      hMeanFitZ3->AddBinContent(ban,1);
    } else if(z<-9 && z>-12){
      c4+=nentries;
      m4+=(meanFitRPHI*nentries);
      n4+=rms*nentries;
      ban=hMeanFit4->FindBin(meanFitRPHI);
      //hMeanFit4->AddBinContent(ban,1);
      hMeanFit4->AddBinContent(ban,1);
      ban=hMeanFitZ4->FindBin(meanFitZ);
      hMeanFitZ4->AddBinContent(ban,1);
    }

  }

  ex[0]=3;
  ex[1]=3;
  ex[2]=3;
  ex[3]=3;

  y[0]=hMeanFit1->GetMean();
  y[1]=hMeanFit2->GetMean();
  y[2]=hMeanFit3->GetMean();
  y[3]=hMeanFit4->GetMean();
  
  ey[0]=hMeanFit1->GetRMS();
  ey[1]=hMeanFit2->GetRMS();
  ey[2]=hMeanFit3->GetRMS();
  ey[3]=hMeanFit4->GetRMS();
  
  yZ[0]=hMeanFitZ1->GetMean();
  yZ[1]=hMeanFitZ2->GetMean();
  yZ[2]=hMeanFitZ3->GetMean();
  yZ[3]=hMeanFitZ4->GetMean();
  
  eyZ[0]=hMeanFitZ1->GetRMS();
  eyZ[1]=hMeanFitZ2->GetRMS();
  eyZ[2]=hMeanFitZ3->GetRMS();
  eyZ[3]=hMeanFitZ4->GetRMS();

  TGraphErrors *gZres = new TGraphErrors(nn,x,y,ex,ey);
  TGraphErrors *gZresZ = new TGraphErrors(nn,x,yZ,ex,eyZ);

  TCanvas *cc1 = new TCanvas("cc1","Title1",1);
  cc1->cd();
  hStatGlob->DrawCopy("LEGO2");
  
  Double_t xx[40],yy[40],exx[40],eyy[40];

  for(Int_t bp = 0; bp<40;bp++){
    xx[bp]=(bp*(TMath::Pi()/20))-TMath::Pi();
    if(TMath::Abs(hMeanPHI[bp]->GetMean())<1e-6) continue;
    yy[bp]=hMeanPHI[bp]->GetMean();
    eyy[bp]=hMeanPHI[bp]->GetRMS();
    exx[bp]=(TMath::Pi())/41;
  }
  TGraphErrors *gPHIres = new TGraphErrors(40,xx,yy,exx,eyy);
  //gPHIres->Fit("pol1","","same",-3,3);

  Double_t xxz[40],yyz[40],exxz[40],eyyz[40];

  for(Int_t bp = 0; bp<40;bp++){
    xxz[bp]=(bp*(TMath::Pi()/20))-TMath::Pi();
    if(TMath::Abs(hMeanPHIz[bp]->GetMean())<1e-6) continue;
    yyz[bp]=hMeanPHIz[bp]->GetMean();
    eyyz[bp]=hMeanPHIz[bp]->GetRMS();
    exxz[bp]=(TMath::Pi())/41;
  }
  TGraphErrors *gPHIresZ = new TGraphErrors(40,xxz,yyz,exxz,eyyz);

  TCanvas *cc4 = new TCanvas("cc4","meanRes VS Z",1);
  cc4->Divide(1,2);
  cc4->cd(1);
  gZres->Draw("AP");
  cc4->cd(2);
  gZresZ->Draw("AP");
  
  TCanvas *cc3 = new TCanvas("cc3","Title3",1);
  cc3->Divide(2,2);
  cc3->cd(1);
  hMeanFitZ1->DrawCopy();
  cc3->cd(2);
  hMeanFitZ2->DrawCopy();
  cc3->cd(3);
  hMeanFitZ3->DrawCopy();
  cc3->cd(4);
  hMeanFitZ4->DrawCopy();
  
  TCanvas *cc6 = new TCanvas("cc6","meanRes(RPHI) VS PHI",1);
  cc6->Divide(1,2);
  cc6->cd(1);
  gPHIres->Draw("AP");

  cc6->cd(2);
  gPHIresZ->Draw("AP");
  
  TCanvas *cc7 = new TCanvas("cc7","Title7",1);
  cc7->Divide(2,2);
  cc7->cd(1);
  hMeanPHI[1]->DrawCopy();
  cc7->cd(2);
  hMeanPHI[2]->DrawCopy();
  cc7->cd(3);
  hMeanPHI[29]->DrawCopy();
  cc7->cd(4);
  hMeanPHI[30]->DrawCopy();



  f1->Close();

  return 0;
}

