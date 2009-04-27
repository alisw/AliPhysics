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
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixD.h"
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
#include "AliITSgeomTGeo.h"

#include "AliITSResidualsAnalysis.h"

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
    fResHistX(0),
    fResHistXLocsddL(0),
    fResHistXLocsddR(0),
    fHistCoordGlobY(0),
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
    fGeom("geometry.root"),
    fUseGausFit(kFALSE)
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
  fResHistX(0),
  fResHistXLocsddL(0),
  fResHistXLocsddR(0),
  fHistCoordGlobY(0),
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
  fGeom(geom),
  fUseGausFit(kFALSE)
{ 
  //
  // Standard Constructor (alitrackpoints)
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
  fResHistX(0),
  fResHistXLocsddL(0),
  fResHistXLocsddR(0),
  fHistCoordGlobY(0),
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
  fGeom("geometry.root"),
  fUseGausFit(kFALSE)
{ 
  //
  // Original Constructor
  //
  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry());
  InitHistograms(volIDs);

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
  if(fHistCoordGlobY)     delete[] fHistCoordGlobY;
  if(fVolResHistRPHI)     delete fVolResHistRPHI; 
  if(fResHistZ)           delete fResHistZ;
  if(fResHistX){         
    for(Int_t i=0; i<fnHist; i++) delete fResHistX[i];
    delete [] fResHistX;
  }
  if(fResHistXLocsddL){         
    for(Int_t i=0; i<fnHist; i++) delete fResHistXLocsddL[i];
    delete [] fResHistXLocsddL;
  }
  if(fResHistXLocsddR){         
    for(Int_t i=0; i<fnHist; i++) delete fResHistXLocsddR[i];
    delete [] fResHistXLocsddR;
  }
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


  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry());
 
  TString histnameRPHI="HistRPHI_volID_",aux;
  TString histnameZ="HistZ_volID_";
  TString histnameX="HistX_volID_";
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
  fResHistX=new TH1F*[fnHist];
  fResHistXLocsddL=new TH1F*[fnHist];
  fResHistXLocsddR=new TH1F*[fnHist];
  fHistCoordGlobY=new TH1F*[fnHist];
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

  Bool_t binning;
  Float_t **binningZPhi;
  Float_t *binningz;
  Float_t *binningphi;

  binningZPhi=CheckSingleLayer(volIDs); 
  fvolidsToBin=new Int_t*[fnPhi*fnZ];
  binningphi=binningZPhi[0];
  binningz=binningZPhi[1];
  binning=SetBinning(volIDs,binningphi,binningz);
    
  if(binning){ //ONLY FOR A SINGLE LAYER!
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
  } else{
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
    fVolResHistRPHI[nhist]=new TH1F("histname","histname",20000,-5.0,5.0);   
    fVolResHistRPHI[nhist]->SetName(aux.Data()); 
    fVolResHistRPHI[nhist]->SetTitle(aux.Data()); 
    
    aux=histnameZ;
    aux+=volIDs->At(nhist);
    fResHistZ[nhist]=new TH1F("histname","histname",20000,-5.0,5.0);   
    fResHistZ[nhist]->SetName(aux.Data()); 
    fResHistZ[nhist]->SetTitle(aux.Data()); 

    aux=histnameX;
    aux+=volIDs->At(nhist);
    fResHistX[nhist]=new TH1F("histname","histname",20000,-5.0,5.0);   
    fResHistX[nhist]->SetName(aux.Data()); 
    fResHistX[nhist]->SetTitle(aux.Data()); 

    aux=histnameX;
    aux+=volIDs->At(nhist);
    aux.Append("LocalSDDLeft");
    fResHistXLocsddL[nhist]=new TH1F("histname","histname",20000,-5.0,5.0);   
    fResHistXLocsddL[nhist]->SetName(aux.Data()); 
    fResHistXLocsddL[nhist]->SetTitle(aux.Data()); 
    aux=histnameX;

    aux+=volIDs->At(nhist);
    aux.Append("LocalSDDRight");
    fResHistXLocsddR[nhist]=new TH1F("histname","histname",20000,-5.0,5.0);   
    fResHistXLocsddR[nhist]->SetName(aux.Data()); 
    fResHistXLocsddR[nhist]->SetTitle(aux.Data()); 

    aux="fHistCoordGlobY";
    fHistCoordGlobY[nhist]=new TH1F("histname","histname",24000,-30.,30.);   
    fHistCoordGlobY[nhist]->SetName(aux.Data()); 
    fHistCoordGlobY[nhist]->SetTitle(aux.Data()); 

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
  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry());
  Int_t volIDalignable,volIDpoint,iModule; 
  AliTrackPoint p;
  AliTrackPointArray* array = 0;
  pointsTree->SetBranchAddress("SP", &array);
  
  
  for(Int_t ivol=0;ivol<fnHist;ivol++){
    Int_t lastused=0;
    volIDalignable=fpTrackVolIDs->At(ivol);
    AliGeomManager::ELayerID iLayer = AliGeomManager::VolUIDToLayer((UShort_t)volIDalignable,iModule);
    
    Int_t nArraysId = lastIndex[iLayer-AliGeomManager::kFirstLayer][iModule];
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
	if (volIDpoint==volIDalignable) continue;
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
	
      }// end loop
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
  

  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry());
  Double_t translGlobal[3];
  //  Double_t radius;
  Double_t phi;
  //  const char *symname,*volpath;
  /*  TGeoPNEntry *pne;
      TGeoPhysicalNode *pn;
      TGeoHMatrix *globMatrix;  
  
  
      symname = AliGeomManager::SymName(volIDalignable);
      pne = gGeoManager->GetAlignableEntry(symname);
      volpath=pne->GetTitle();
      pn=gGeoManager->MakePhysicalNode(volpath);
      globMatrix=pn->GetMatrix();
  */
  
  AliGeomManager::GetOrigTranslation(volIDalignable,translGlobal);
  //  radius=TMath::Sqrt(transGlobal[0]*transGlobal[0]+transGlobal[1]*transGlobal[1]);
  phi=TMath::ATan2(translGlobal[1],translGlobal[0]);
  fhistVolNptsUsed->Fill(phi,translGlobal[2]);
  if(!usedVol){
    fhistVolUsed->Fill(phi,translGlobal[2]);
  }

  /*  symname = AliGeomManager::SymName(volIDpoint);
      pne = gGeoManager->GetAlignableEntry(symname);
      volpath=pne->GetTitle();
      pn=gGeoManager->MakePhysicalNode(volpath);
      globMatrix=pn->GetMatrix();
      transGlobal=globMatrix->GetTranslation();
  */
  AliGeomManager::GetOrigTranslation(volIDpoint,translGlobal);
  //  radius=TMath::Sqrt(transGlobal[0]*transGlobal[0]+transGlobal[1]*transGlobal[1]);
  phi=TMath::ATan2(translGlobal[1],translGlobal[0]);

  fhistCorrVol[ivol]->Fill(phi,translGlobal[2]);

  return;
}
  
//____________________________________________________________________________
void AliITSResidualsAnalysis::FillResidualsH(AliTrackPointArray *points,
					     AliTrackPointArray *pTrack) const
{
  //
  // Method that fills the histograms with the residuals
  //
  
  Int_t volIDpoint;  
  Float_t xyz[3],xyz2[3];
  Double_t xyzD[3],xyz2D[3];
  Double_t loc[3],loc2[3];

  Float_t resRPHI,resGlob,resZ,resX;

  Double_t pullrphi,sign,phi;
  AliTrackPoint p,pTr;

  for(Int_t ipoint=0;ipoint<points->GetNPoints();ipoint++){

    //pTrack->GetPoint(pTr,ipoint);
    points->GetPoint(p,ipoint);
    volIDpoint=(Int_t)p.GetVolumeID();
    p.GetXYZ(xyz);

    pTrack->GetPoint(pTr,ipoint);
    pTr.GetXYZ(xyz2);

    for(Int_t i=0;i<3;i++){
      xyzD[i]=xyz[i];
      xyz2D[i]=xyz2[i];
    }

    phi = TMath::ATan2(xyz[1],xyz[0]);//<-watch out: phi of the pPoints!
 
    resZ=xyz2[2]-xyz[2];
    resX=xyz2[0]-xyz[0];

    resRPHI=TMath::Sqrt((xyz2[0]-xyz[0])*(xyz2[0]-xyz[0])+(xyz2[1]-xyz[1])*(xyz2[1]-xyz[1]));

    sign=TMath::ATan2(xyz2[1],xyz2[0])-TMath::ATan2(xyz[1],xyz[0]);
    if(sign!=0.){
      sign=sign/TMath::Abs(sign);
      resRPHI=resRPHI*sign;

    }
    else{
      pullrphi=0.;
      resRPHI=0.;
    }
    
    resGlob=TMath::Sqrt((xyz2[0]-xyz[0])*(xyz2[0]-xyz[0])+(xyz2[1]-xyz[1])*(xyz2[1]-xyz[1])+(xyz2[2]-xyz[2])*(xyz2[2]-xyz[2]));

    for(Int_t ivolIDs=0;ivolIDs<fpTrackVolIDs->GetSize();ivolIDs++){
      if(volIDpoint==fpTrackVolIDs->At(ivolIDs)){

	fVolResHistRPHI[ivolIDs]->Fill(resRPHI);
	fResHistZ[ivolIDs]->Fill(resZ);
	fResHistX[ivolIDs]->Fill(resX);
	fHistCoordGlobY[ivolIDs]->Fill(xyz[1]); 

	Int_t modIndex = -1; // SDD Section
	if(AliGeomManager::VolUIDToLayer(volIDpoint)==3) modIndex=volIDpoint-6144+240;
	if(AliGeomManager::VolUIDToLayer(volIDpoint)==4) modIndex=volIDpoint-8192+240+84;
	if(modIndex>0){
	  AliITSgeomTGeo::GlobalToLocal(modIndex,xyzD,loc); // error here!?
	  AliITSgeomTGeo::GlobalToLocal(modIndex,xyz2D,loc2);
	  Float_t rexloc=loc2[0]-loc[0];
	  //cout<<"Residual: "<<volIDpoint<<" "<<loc[0]<<" -> "<<rexloc<<endl;
	  if(loc[0]<0){
	    fResHistXLocsddR[ivolIDs]->Fill(rexloc);
	  }else{
	    fResHistXLocsddL[ivolIDs]->Fill(rexloc);
	  }
	}
	fResHistGlob[ivolIDs]->Fill(resGlob);

	fTrackDirPhiAll->Fill(phi);
	fTrackDirPhi[ivolIDs]->Fill(phi);

	if(fsingleLayer){
	  Int_t binz,binphi;
	  Float_t globalPhi,globalZ;
	  if(kTRUE||(fvolidsToBin[ivolIDs][0]!=volIDpoint)){
	    binphi=GetBinPhiZ((Int_t)volIDpoint,&binz);
	  }
	  else{
	    // This in the case of alignment of one entire layer 
	    // (fnHIst=layersize) may reduce iterations: 
	    // remind of that fsingleLayer->fnHista<layerSize
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
Bool_t AliITSResidualsAnalysis::SaveHists(Int_t minNpoints, TString outname) const
{
  //  
  // Saves the histograms into a tree and saves the tree into a file
  //


  // Output file
  TFile *hFile=new TFile(outname.Data(),"RECREATE","File containing the Residuals Tree");

  // TTree with the residuals
  TTree *analysisTree=new TTree("analysisTree","Tree with the residuals");

  // Declares Variables to be stored into the TTree
  TF1 *gauss=new TF1("gauss","gaus",-10.,10.);
  Int_t volID,entries,nHistAnalyzed=0;
  Double_t meanResRPHI,meanResZ,meanResX,rmsResRPHI,rmsResZ,rmsResX,coordVol[3],x,y,z;
  TH1F *histRPHI = new TH1F();
  TH1F *histZ = new TH1F();
  TH1F *histX = new TH1F();
  TH1F *histXLocsddL = new TH1F();
  TH1F *histXLocsddR = new TH1F();
  TH1F *histCoordGlobY = new TH1F();
  // Note: 0 = RPHI, 1 = Z


  // Branching the TTree
  analysisTree->Branch("volID",&volID,"volID/I");
  analysisTree->Branch("x",&x,"x/D");
  analysisTree->Branch("y",&y,"y/D");
  analysisTree->Branch("z",&z,"z/D");
  analysisTree->Branch("meanResRPHI",&meanResRPHI,"meanResRPHI/D");
  analysisTree->Branch("meanResZ",&meanResZ,"meanResZ/D");
  analysisTree->Branch("meanResX",&meanResX,"meanResX/D");
  analysisTree->Branch("rmsResRPHI",&rmsResRPHI,"rmsResRPHI/D");
  analysisTree->Branch("rmsResZ",&rmsResZ,"rmsResZ/D");

  analysisTree->Branch("histRPHI","TH1F",&histRPHI,128000,0);
  analysisTree->Branch("histZ","TH1F",&histZ,128000,0);
  analysisTree->Branch("histX","TH1F",&histX,128000,0);
  analysisTree->Branch("histXLocsddL","TH1F",&histXLocsddL,128000,0);
  analysisTree->Branch("histXLocsddR","TH1F",&histXLocsddR,128000,0);
  analysisTree->Branch("histCoordGlobY","TH1F",&histCoordGlobY,128000,0);

  Int_t blimps=0;

  for(Int_t j=0;j<fnHist;j++){

    volID=fpTrackVolIDs->At(j);
    AliGeomManager::GetTranslation(volID,coordVol);
    x=coordVol[0];
    y=coordVol[1];
    z=coordVol[2];
    
    entries=(Int_t)(fResHistGlob[j]->GetEntries());
    blimps+=entries;

    if(entries>=minNpoints){
      nHistAnalyzed++;

      // Entries
      //entries=(Int_t)fVolResHistRPHI[j]->GetEntries();

      // Filling the RPHI
      histRPHI=fVolResHistRPHI[j];
      rmsResRPHI=fVolResHistRPHI[j]->GetRMS();
      if(fUseGausFit){
        // Fit (for average)
	gauss->SetRange(-3*rmsResRPHI,3*rmsResRPHI);
	fVolResHistRPHI[j]->Fit("gauss","QRN");
	meanResRPHI=gauss->GetParameter(1);
      }else{
	meanResRPHI=fVolResHistRPHI[j]->GetMean();
      }

      // Filling the Z
      histZ=fResHistZ[j];
      rmsResZ=fResHistZ[j]->GetRMS();
      if(fUseGausFit){
        // Fit (for average)
	gauss->SetRange(-3*rmsResZ,3*rmsResZ);
	fResHistZ[j]->Fit("gauss","QRN");
	meanResZ=gauss->GetParameter(1);
      }else{
	meanResZ=fResHistZ[j]->GetMean();
      }

      // Filling the X
      histX=fResHistX[j];
      rmsResX=fResHistX[j]->GetRMS();
      if(fUseGausFit){
        // Fit (for average)
	gauss->SetRange(-3*rmsResX,3*rmsResX);
	fResHistX[j]->Fit("gauss","QRN");
	meanResX=gauss->GetParameter(1);
      }else{
	meanResX=fResHistX[j]->GetMean();
      }

      histXLocsddL=fResHistXLocsddL[j];
      histXLocsddR=fResHistXLocsddR[j];
      histCoordGlobY=fHistCoordGlobY[j];

      analysisTree->Fill();
    }else{

      // Entries
      //entries=(Int_t)fVolResHistRPHI[j]->GetEntries();

      // Filling the RPHI
      histRPHI=fVolResHistRPHI[j];
      rmsResRPHI=-1.0;
      meanResRPHI=0.0;

      // Filling the Z
      histZ=fResHistZ[j];
      rmsResZ=-1.0;
      meanResZ=0.0;

      // Filling the X
      histX=fResHistX[j];
      rmsResX=-1.0;
      meanResX=0.0;
      histXLocsddL=fResHistXLocsddL[j];
      histXLocsddR=fResHistXLocsddR[j];
      histCoordGlobY=fHistCoordGlobY[j];
 
      analysisTree->Fill();

    }

  }
  delete gauss;

  cout<<"-> Modules Analyzed: "<<nHistAnalyzed<<endl;
  cout<<"   With "<<blimps<<" events"<<endl;

  if(blimps>0){ 
    hFile->cd();
    analysisTree->Write();
    fVolNTracks->Write();
    fhEmpty->Write();
    if(fWriteHist){
      //TCanvas *color = new TCanvas("color","fhistVolUsed",800,600);
      //fhistVolUsed->DrawCopy("COLZ");
      fSigmaVolZ->Write();
      fhistVolUsed->Write();
      /*      fTrackDirPhiAll->Write();
	      fTrackDirLambdaAll->Write();
	      fTrackDirLambda2All->Write();
	      fTrackDirAlphaAll->Write();
	      fTrackDirAll->Write();
	      fTrackDir2All->Write();
	      fTrackDirXZAll->Write();
	      hFile->Close();*/
      fhistVolNptsUsed->Write();
      hFile->mkdir("CorrVol");
      hFile->cd("CorrVol");
      for(Int_t corr=0;corr<fnHist;corr++)fhistCorrVol[corr]->Write();
    }
    hFile->cd();
    //    fhistVolNptsUsed->Write();
    hFile->Close();
    return kTRUE;
  }else {
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
      return binningzphi;
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
    //////////////////////
    fnPhi=1;//kPhiSPD1;
    fnZ=1;//nZSPD1;
    binningzphi[0]=new Float_t[1];
    binningzphi[1]=new Float_t[1];
    fCoordToBinTable=new Double_t**[1];
    for(Int_t j=0;j<fnPhi;j++){
      fCoordToBinTable[j]=new Double_t*[1];
    }
    return binningzphi;
    /////////////////////
    // return 0x0;
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
  //  const char *volpath,*symname;
  Int_t iModule;
  Int_t *orderArrayPhi,*orderArrayZ;
  UShort_t volID;
  Double_t *phiArray,*zArray,*phiArrayOrdered,*zArrayOrdered; 
  Double_t translGlobal[3];
  Double_t lastPhimin=-10;
  Double_t lastZmin=-99;
  Int_t ***orderPhiZ;
  /*  TGeoPNEntry *pne;
      TGeoPhysicalNode *pn;
      TGeoHMatrix *globMatrix;*/
  
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
    /*    symname = AliGeomManager::SymName(volID);
	  pne = gGeoManager->GetAlignableEntry(symname);
	  volpath=pne->GetTitle();
	  pn=gGeoManager->MakePhysicalNode(volpath);
	  globMatrix=pn->GetMatrix();
	  translGlobal=globMatrix->GetTranslation();
	  
    */
    AliGeomManager::GetOrigTranslation(volID,translGlobal);
    
    for(Int_t j=0;j<lastPhi;j++){
      used=kFALSE;
      if(TMath::Abs(phiArray[j]-TMath::ATan2(translGlobal[1],translGlobal[0]))<2*TMath::Pi()/(10*fnPhi)){//10 is a safety factor but....
	fvolidsToBin[iModule][1]=j;
	used=kTRUE;
	break;
      }
    }
    if(!used){
      phiArray[lastPhi]=TMath::ATan2(translGlobal[1],translGlobal[0]);
      fvolidsToBin[iModule][1]=lastPhi;
      if(phiArray[lastPhi]<lastPhimin)lastPhimin=phiArray[lastPhi];
      lastPhi++;
      if(lastPhi>fnPhi){
	printf("Wrong Phi! \n");
	return kFALSE;}
    }
    
    for(Int_t j=0;j<lastZ;j++){
      used=kFALSE;
      if(TMath::Abs(zArray[j]-translGlobal[2])<0.1){
	fvolidsToBin[iModule][2]=j;
	used=kTRUE;
	break;
      }
    }
    if(!used){
      fvolidsToBin[iModule][2]=lastZ;
      zArray[lastZ]=translGlobal[2];
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
      printf("Now (binphi,binz)= %d, %d e (phi,z)=%f,%f \n",j,i,fCoordToBinTable[j][i][0],fCoordToBinTable[j][i][1]);
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

//___________________________________________________________________________
Int_t AliITSResidualsAnalysis::WhichSector(Int_t module) const
{
  //
  // This method returns the number of the SPD Sector
  // to which belongs the module (Sectors 0-9)
   
  //--->cSect = 0 <---
  if(module==2048
     || module==2049
     || module==2050
     || module==2051
     || module==2052
     || module==2053
     || module==2054
     || module==2055
     || module==4096
     || module==4097
     || module==4098
     || module==4099
     || module==4100
     || module==4101
     || module==4102
     || module==4103
     || module==4104
     || module==4105
     || module==4106
     || module==4107
     || module==4108
     || module==4109
     || module==4110
     || module==4111) return 0;
     
  //--->cSect = 1 <---    
  if(module==2056
     || module==2057
     || module==2058
     || module==2059
     || module==2060
     || module==2061
     || module==2062
     || module==2063
     || module==4112
     || module==4113
     || module==4114
     || module==4115
     || module==4116
     || module==4117
     || module==4118
     || module==4119
     || module==4120
     || module==4121
     || module==4122
     || module==4123
     || module==4124
     || module==4125
     || module==4126
     || module==4127) return 1;

  //--->cSect = 2 <---
  if(module==2064
     || module==2065
     || module==2066
     || module==2067
     || module==2068
     || module==2069
     || module==2070
     || module==2071
     || module==4128
     || module==4129
     || module==4130
     || module==4131
     || module==4132
     || module==4133
     || module==4134
     || module==4135
     || module==4136
     || module==4137
     || module==4138
     || module==4139
     || module==4140
     || module==4141
     || module==4142
     || module==4143) return 2;

  //--->cSect = 3 <---
  if(module==2072
     || module==2073
     || module==2074
     || module==2075
     || module==2076
     || module==2077
     || module==2078
     || module==2079
     || module==4144
     || module==4145
     || module==4146
     || module==4147
     || module==4148
     || module==4149
     || module==4150
     || module==4151
     || module==4152
     || module==4153
     || module==4154
     || module==4155
     || module==4156
     || module==4157
     || module==4158
     || module==4159) return 3;

  //--->cSect = 4 <---
  if(module==2080
     || module==2081
     || module==2082
     || module==2083
     || module==2084
     || module==2085
     || module==2086
     || module==2087
     || module==4160
     || module==4161
     || module==4162
     || module==4163
     || module==4164
     || module==4165
     || module==4166
     || module==4167
     || module==4168
     || module==4169
     || module==4170
     || module==4171
     || module==4172
     || module==4173
     || module==4174
     || module==4175) return 4;
  
  //--->cSect = 5 <---
  if(module==2088
     || module==2089
     || module==2090
     || module==2091
     || module==2092
     || module==2093
     || module==2094
     || module==2095
     || module==4176
     || module==4177
     || module==4178
     || module==4179
     || module==4180
     || module==4181
     || module==4182
     || module==4183
     || module==4184
     || module==4185
     || module==4186
     || module==4187
     || module==4188
     || module==4189
     || module==4190
     || module==4191) return 5;

  //--->cSect = 6 <---
  if(module==2096
     || module==2097
     || module==2098
     || module==2099
     || module==2100
     || module==2101
     || module==2102
     || module==2103
     || module==4192
     || module==4193
     || module==4194
     || module==4195
     || module==4196
     || module==4197
     || module==4198
     || module==4199
     || module==4200
     || module==4201
     || module==4202
     || module==4203
     || module==4204
     || module==4205
     || module==4206
     || module==4207) return 6;

  //--->cSect = 7 <---
  if(module==2104
     || module==2105
     || module==2106
     || module==2107
     || module==2108
     || module==2109
     || module==2110
     || module==2111
     || module==4208
     || module==4209
     || module==4210
     || module==4211
     || module==4212
     || module==4213
     || module==4214
     || module==4215
     || module==4216
     || module==4217
     || module==4218
     || module==4219
     || module==4220
     || module==4221
     || module==4222
     || module==4223) return 7;

  //--->cSect = 8 <---
  if(module==2112
     || module==2113
     || module==2114
     || module==2115
     || module==2116
     || module==2117
     || module==2118
     || module==2119
     || module==4224
     || module==4225
     || module==4226
     || module==4227
     || module==4228
     || module==4229
     || module==4230
     || module==4231
     || module==4232
     || module==4233
     || module==4234
     || module==4235
     || module==4236
     || module==4237
     || module==4238
     || module==4239) return 8;

  //--->cSect = 9 <---
  if(module==2120
     || module==2121
     || module==2122
     || module==2123
     || module==2124
     || module==2125
     || module==2126
     || module==2127
     || module==4240
     || module==4241
     || module==4242
     || module==4243
     || module==4244
     || module==4245
     || module==4246
     || module==4247
     || module==4248
     || module==4249
     || module==4250
     || module==4251
     || module==4252
     || module==4253
     || module==4254
     || module==4255) return 9;

  //printf("Module not belonging to SPD, sorry!");
  return -1;

}

//____________________________________________________________________________
TArrayI* AliITSResidualsAnalysis::GetSPDSectorsVolids(Int_t sectors[10]) const
{
  //
  // This method gets the volID Array for the chosen sectors.
  // You have to pass an array with a 1 for each selected sector.
  // i.e. sectors[10] = {1,1,0,0,0,0,0,0,1,0} -> Sector 0, 1, 9 selected.
  //

  Int_t nSect=0;
  Int_t iModule=0;

  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry());

  for(Int_t co=0;co<10;co++){ //counts the number of sectors chosen
    if(sectors[co]==1) nSect++;
  }
  
  if(nSect<1){ //if no sector chosen -> exit
    Printf("Error! No Sector/s Selected!");
    return 0x0;
  }

  TArrayI *volIDs = new TArrayI(nSect*24);
  
    if(sectors[0]==1){ //--->cSect = 0 <---
      volIDs->AddAt(2048,iModule); iModule++;
      volIDs->AddAt(2049,iModule); iModule++;
      volIDs->AddAt(2050,iModule); iModule++;
      volIDs->AddAt(2051,iModule); iModule++;
      volIDs->AddAt(2052,iModule); iModule++;
      volIDs->AddAt(2053,iModule); iModule++;
      volIDs->AddAt(2054,iModule); iModule++;
      volIDs->AddAt(2055,iModule); iModule++;
      volIDs->AddAt(4096,iModule); iModule++;
      volIDs->AddAt(4097,iModule); iModule++;
      volIDs->AddAt(4098,iModule); iModule++;
      volIDs->AddAt(4099,iModule); iModule++;
      volIDs->AddAt(4100,iModule); iModule++;
      volIDs->AddAt(4101,iModule); iModule++;
      volIDs->AddAt(4102,iModule); iModule++;
      volIDs->AddAt(4103,iModule); iModule++;
      volIDs->AddAt(4104,iModule); iModule++;
      volIDs->AddAt(4105,iModule); iModule++;
      volIDs->AddAt(4106,iModule); iModule++;
      volIDs->AddAt(4107,iModule); iModule++;
      volIDs->AddAt(4108,iModule); iModule++;
      volIDs->AddAt(4109,iModule); iModule++;
      volIDs->AddAt(4110,iModule); iModule++;
      volIDs->AddAt(4111,iModule); iModule++;
    }
    if(sectors[1]==1){ //--->cSect = 1 <//---
      volIDs->AddAt(2056,iModule); iModule++;
      volIDs->AddAt(2057,iModule); iModule++;
      volIDs->AddAt(2058,iModule); iModule++;
      volIDs->AddAt(2059,iModule); iModule++;
      volIDs->AddAt(2060,iModule); iModule++;
      volIDs->AddAt(2061,iModule); iModule++;
      volIDs->AddAt(2062,iModule); iModule++;
      volIDs->AddAt(2063,iModule); iModule++;
      volIDs->AddAt(4112,iModule); iModule++;
      volIDs->AddAt(4113,iModule); iModule++;
      volIDs->AddAt(4114,iModule); iModule++;
      volIDs->AddAt(4115,iModule); iModule++;
      volIDs->AddAt(4116,iModule); iModule++;
      volIDs->AddAt(4117,iModule); iModule++;
      volIDs->AddAt(4118,iModule); iModule++;
      volIDs->AddAt(4119,iModule); iModule++;
      volIDs->AddAt(4120,iModule); iModule++;
      volIDs->AddAt(4121,iModule); iModule++;
      volIDs->AddAt(4122,iModule); iModule++;
      volIDs->AddAt(4123,iModule); iModule++;
      volIDs->AddAt(4124,iModule); iModule++;
      volIDs->AddAt(4125,iModule); iModule++;
      volIDs->AddAt(4126,iModule); iModule++;
      volIDs->AddAt(4127,iModule); iModule++;
    }
    if(sectors[2]==1){//--->cSect = 2 <//---
      volIDs->AddAt(2064,iModule); iModule++;
      volIDs->AddAt(2065,iModule); iModule++;
      volIDs->AddAt(2066,iModule); iModule++;
      volIDs->AddAt(2067,iModule); iModule++;
      volIDs->AddAt(2068,iModule); iModule++;
      volIDs->AddAt(2069,iModule); iModule++;
      volIDs->AddAt(2070,iModule); iModule++;
      volIDs->AddAt(2071,iModule); iModule++;
      volIDs->AddAt(4128,iModule); iModule++;
      volIDs->AddAt(4129,iModule); iModule++;
      volIDs->AddAt(4130,iModule); iModule++;
      volIDs->AddAt(4131,iModule); iModule++;
      volIDs->AddAt(4132,iModule); iModule++;
      volIDs->AddAt(4133,iModule); iModule++;
      volIDs->AddAt(4134,iModule); iModule++;
      volIDs->AddAt(4135,iModule); iModule++;
      volIDs->AddAt(4136,iModule); iModule++;
      volIDs->AddAt(4137,iModule); iModule++;
      volIDs->AddAt(4138,iModule); iModule++;
      volIDs->AddAt(4139,iModule); iModule++;
      volIDs->AddAt(4140,iModule); iModule++;
      volIDs->AddAt(4141,iModule); iModule++;
      volIDs->AddAt(4142,iModule); iModule++;
      volIDs->AddAt(4143,iModule); iModule++;
    }
    if(sectors[3]==1){//--->cSect = 3 <//---
      volIDs->AddAt(2072,iModule); iModule++;
      volIDs->AddAt(2073,iModule); iModule++;
      volIDs->AddAt(2074,iModule); iModule++;
      volIDs->AddAt(2075,iModule); iModule++;
      volIDs->AddAt(2076,iModule); iModule++;
      volIDs->AddAt(2077,iModule); iModule++;
      volIDs->AddAt(2078,iModule); iModule++;
      volIDs->AddAt(2079,iModule); iModule++;
      volIDs->AddAt(4144,iModule); iModule++;
      volIDs->AddAt(4145,iModule); iModule++;
      volIDs->AddAt(4146,iModule); iModule++;
      volIDs->AddAt(4147,iModule); iModule++;
      volIDs->AddAt(4148,iModule); iModule++;
      volIDs->AddAt(4149,iModule); iModule++;
      volIDs->AddAt(4150,iModule); iModule++;
      volIDs->AddAt(4151,iModule); iModule++;
      volIDs->AddAt(4152,iModule); iModule++;
      volIDs->AddAt(4153,iModule); iModule++;
      volIDs->AddAt(4154,iModule); iModule++;
      volIDs->AddAt(4155,iModule); iModule++;
      volIDs->AddAt(4156,iModule); iModule++;
      volIDs->AddAt(4157,iModule); iModule++;
      volIDs->AddAt(4158,iModule); iModule++;
      volIDs->AddAt(4159,iModule); iModule++;
    }
    if(sectors[4]==1){//--->cSect = 4 <//---
      volIDs->AddAt(2080,iModule); iModule++;
      volIDs->AddAt(2081,iModule); iModule++;
      volIDs->AddAt(2082,iModule); iModule++;
      volIDs->AddAt(2083,iModule); iModule++;
      volIDs->AddAt(2084,iModule); iModule++;
      volIDs->AddAt(2085,iModule); iModule++;
      volIDs->AddAt(2086,iModule); iModule++;
      volIDs->AddAt(2087,iModule); iModule++;
      volIDs->AddAt(4160,iModule); iModule++;
      volIDs->AddAt(4161,iModule); iModule++;
      volIDs->AddAt(4162,iModule); iModule++;
      volIDs->AddAt(4163,iModule); iModule++;
      volIDs->AddAt(4164,iModule); iModule++;
      volIDs->AddAt(4165,iModule); iModule++;
      volIDs->AddAt(4166,iModule); iModule++;
      volIDs->AddAt(4167,iModule); iModule++;
      volIDs->AddAt(4168,iModule); iModule++;
      volIDs->AddAt(4169,iModule); iModule++;
      volIDs->AddAt(4170,iModule); iModule++;
      volIDs->AddAt(4171,iModule); iModule++;
      volIDs->AddAt(4172,iModule); iModule++;
      volIDs->AddAt(4173,iModule); iModule++;
      volIDs->AddAt(4174,iModule); iModule++;
      volIDs->AddAt(4175,iModule); iModule++;
    }
    if(sectors[5]==1){//--->cSect = 5 <//---
      volIDs->AddAt(2088,iModule); iModule++;
      volIDs->AddAt(2089,iModule); iModule++;
      volIDs->AddAt(2090,iModule); iModule++;
      volIDs->AddAt(2091,iModule); iModule++;
      volIDs->AddAt(2092,iModule); iModule++;
      volIDs->AddAt(2093,iModule); iModule++;
      volIDs->AddAt(2094,iModule); iModule++;
      volIDs->AddAt(2095,iModule); iModule++;
      volIDs->AddAt(4176,iModule); iModule++;
      volIDs->AddAt(4177,iModule); iModule++;
      volIDs->AddAt(4178,iModule); iModule++;
      volIDs->AddAt(4179,iModule); iModule++;
      volIDs->AddAt(4180,iModule); iModule++;
      volIDs->AddAt(4181,iModule); iModule++;
      volIDs->AddAt(4182,iModule); iModule++;
      volIDs->AddAt(4183,iModule); iModule++;
      volIDs->AddAt(4184,iModule); iModule++;
      volIDs->AddAt(4185,iModule); iModule++;
      volIDs->AddAt(4186,iModule); iModule++;
      volIDs->AddAt(4187,iModule); iModule++;
      volIDs->AddAt(4188,iModule); iModule++;
      volIDs->AddAt(4189,iModule); iModule++;
      volIDs->AddAt(4190,iModule); iModule++;
      volIDs->AddAt(4191,iModule); iModule++;
    }
    if(sectors[6]==1){//--->cSect = 6 <//---
      volIDs->AddAt(2096,iModule); iModule++;
      volIDs->AddAt(2097,iModule); iModule++;
      volIDs->AddAt(2098,iModule); iModule++;
      volIDs->AddAt(2099,iModule); iModule++;
      volIDs->AddAt(2100,iModule); iModule++;
      volIDs->AddAt(2101,iModule); iModule++;
      volIDs->AddAt(2102,iModule); iModule++;
      volIDs->AddAt(2103,iModule); iModule++;
      volIDs->AddAt(4192,iModule); iModule++;
      volIDs->AddAt(4193,iModule); iModule++;
      volIDs->AddAt(4194,iModule); iModule++;
      volIDs->AddAt(4195,iModule); iModule++;
      volIDs->AddAt(4196,iModule); iModule++;
      volIDs->AddAt(4197,iModule); iModule++;
      volIDs->AddAt(4198,iModule); iModule++;
      volIDs->AddAt(4199,iModule); iModule++;
      volIDs->AddAt(4200,iModule); iModule++;
      volIDs->AddAt(4201,iModule); iModule++;
      volIDs->AddAt(4202,iModule); iModule++;
      volIDs->AddAt(4203,iModule); iModule++;
      volIDs->AddAt(4204,iModule); iModule++;
      volIDs->AddAt(4205,iModule); iModule++;
      volIDs->AddAt(4206,iModule); iModule++;
      volIDs->AddAt(4207,iModule); iModule++;
    }
     if(sectors[7]==1){ //--->cSect = 7 <//---
       volIDs->AddAt(2104,iModule); iModule++;
       volIDs->AddAt(2105,iModule); iModule++;
       volIDs->AddAt(2106,iModule); iModule++;
       volIDs->AddAt(2107,iModule); iModule++;
       volIDs->AddAt(2108,iModule); iModule++;
       volIDs->AddAt(2109,iModule); iModule++;
       volIDs->AddAt(2110,iModule); iModule++;
       volIDs->AddAt(2111,iModule); iModule++;
       volIDs->AddAt(4208,iModule); iModule++;
       volIDs->AddAt(4209,iModule); iModule++;
       volIDs->AddAt(4210,iModule); iModule++;
       volIDs->AddAt(4211,iModule); iModule++;
       volIDs->AddAt(4212,iModule); iModule++;
       volIDs->AddAt(4213,iModule); iModule++;
       volIDs->AddAt(4214,iModule); iModule++;
       volIDs->AddAt(4215,iModule); iModule++;
       volIDs->AddAt(4216,iModule); iModule++;
       volIDs->AddAt(4217,iModule); iModule++;
       volIDs->AddAt(4218,iModule); iModule++;
       volIDs->AddAt(4219,iModule); iModule++;
       volIDs->AddAt(4220,iModule); iModule++;
       volIDs->AddAt(4221,iModule); iModule++;
       volIDs->AddAt(4222,iModule); iModule++;
       volIDs->AddAt(4223,iModule); iModule++;
     }
     if(sectors[8]==1){//--->cSect = 8 <//---
       volIDs->AddAt(2112,iModule); iModule++;
       volIDs->AddAt(2113,iModule); iModule++;
       volIDs->AddAt(2114,iModule); iModule++;
       volIDs->AddAt(2115,iModule); iModule++;
       volIDs->AddAt(2116,iModule); iModule++;
       volIDs->AddAt(2117,iModule); iModule++;
       volIDs->AddAt(2118,iModule); iModule++;
       volIDs->AddAt(2119,iModule); iModule++;
       volIDs->AddAt(4224,iModule); iModule++;
       volIDs->AddAt(4225,iModule); iModule++;
       volIDs->AddAt(4226,iModule); iModule++;
       volIDs->AddAt(4227,iModule); iModule++;
       volIDs->AddAt(4228,iModule); iModule++;
       volIDs->AddAt(4229,iModule); iModule++;
       volIDs->AddAt(4230,iModule); iModule++;
       volIDs->AddAt(4231,iModule); iModule++;
       volIDs->AddAt(4232,iModule); iModule++;
       volIDs->AddAt(4233,iModule); iModule++;
       volIDs->AddAt(4234,iModule); iModule++;
       volIDs->AddAt(4235,iModule); iModule++;
       volIDs->AddAt(4236,iModule); iModule++;
       volIDs->AddAt(4237,iModule); iModule++;
       volIDs->AddAt(4238,iModule); iModule++;
       volIDs->AddAt(4239,iModule); iModule++;
     }
     if(sectors[9]==1){//--->cSect = 9 <//---
       volIDs->AddAt(2120,iModule); iModule++;
       volIDs->AddAt(2121,iModule); iModule++;
       volIDs->AddAt(2122,iModule); iModule++;
       volIDs->AddAt(2123,iModule); iModule++;
       volIDs->AddAt(2124,iModule); iModule++;
       volIDs->AddAt(2125,iModule); iModule++;
       volIDs->AddAt(2126,iModule); iModule++;
       volIDs->AddAt(2127,iModule); iModule++;
       volIDs->AddAt(4240,iModule); iModule++;
       volIDs->AddAt(4241,iModule); iModule++;
       volIDs->AddAt(4242,iModule); iModule++;
       volIDs->AddAt(4243,iModule); iModule++;
       volIDs->AddAt(4244,iModule); iModule++;
       volIDs->AddAt(4245,iModule); iModule++;
       volIDs->AddAt(4246,iModule); iModule++;
       volIDs->AddAt(4247,iModule); iModule++;
       volIDs->AddAt(4248,iModule); iModule++;
       volIDs->AddAt(4249,iModule); iModule++;
       volIDs->AddAt(4250,iModule); iModule++;
       volIDs->AddAt(4251,iModule); iModule++;
       volIDs->AddAt(4252,iModule); iModule++;
       volIDs->AddAt(4253,iModule); iModule++;
       volIDs->AddAt(4254,iModule); iModule++;
       volIDs->AddAt(4255,iModule); iModule++;
     }

  return volIDs;

}

//____________________________________________________________________________
TArrayI* AliITSResidualsAnalysis::GetITSLayersVolids(Int_t layers[6]) const
{
  //
  // This method gets the volID Array for the chosen layers.
  // You have to pass an array with a 1 for each selected layer.
  // i.e. layers[6] = {1,1,0,0,1,1} -> SPD + SSD
  //

  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry());

  Int_t size=0,last=0;

  // evaluates the size of the array
  for(Int_t i=0;i<6;i++) if(layers[i]==1) size+=AliGeomManager::LayerSize(i+1);

  if(size==0){
    printf("Error: no layer selected");
    return 0x0;
  }

  TArrayI *volids = new TArrayI(size);

  // fills the volId array only for the chosen layers
  for(Int_t ilayer=1;ilayer<7;ilayer++){
    
    if(layers[ilayer-1]!=1) continue;
    
    for(Int_t imod=0;imod<AliGeomManager::LayerSize(ilayer);imod++){
      volids->AddAt(AliGeomManager::LayerToVolUID(ilayer,imod),last);
      last++;
    }
  }
  
  return volids;

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
						 TString outname)
{
  // CalculateResiduals for a set of detector volumes.
  // Tracks are fitted only within
  // the range defined by the user
  // (by layerRangeMin and layerRangeMax)
  // or within the set of volidsfit
  // Repeat the procedure 'iterations' times

  Int_t nVolIds = volids->GetSize();
  if (nVolIds == 0) { AliError("Volume IDs array is empty!"); return; }

  // Load only the tracks with at least one
  // space point in the set of volume (volids)


  //AliAlignmentTracks::SetPointsFilename(GetFileNameTrackPoints()); 
  AliAlignmentTracks::BuildIndex();

  ListVolUsed(fPointsTree,fArrayIndex,fLastIndex);  
  AliTrackPointArray **points;  
  
  Int_t pointsDim;
  LoadPoints(volids, points,pointsDim);

  Int_t nArrays = fPointsTree->GetEntries();

  if (nArrays == 0){ AliError("Points array is empty!"); return; }
  AliTrackFitter *fitter = CreateFitter();

  Int_t ecount=0;
  Int_t totcount=0;

  for (Int_t iArray = 0; iArray < nArrays; iArray++){
 
    cout<<"Investigating "<<iArray<<"/"<<nArrays<<endl;
    
    if (!points[iArray]){
      cout<<" Skipping: "<<iArray<<endl;
      continue;
    }
         
    fitter->SetTrackPointArray(points[iArray],kTRUE); // Watch out, problems
                                                      // when few sectors
                               
    totcount++;

    // *** FITTING ***
    if(fitter->Fit(volids,volidsfit,layerRangeMin,layerRangeMax) == kFALSE){ 
      ecount++;
      cout<<"->BAD: "<<iArray<<endl;
      continue;
    } //else cout<<"->GOOD: "<<iArray<<endl;

    AliTrackPointArray *pVolId,*pTrack;

    fitter->GetTrackResiduals(pVolId,pTrack);
    FillResidualsH(pVolId,pTrack);

  }
  
  cout<<"   -> nVolIds: "<<nVolIds<<endl;
  cout<<"   -> Non-Fitted tracks: "<<ecount<<"/"<<totcount<<endl; 
  
  UnloadPoints(pointsDim,points);
  SaveHists(3,outname);

  
  return;
  
}


//______________________________________________________________________________
void AliITSResidualsAnalysis::ProcessVolumes(Int_t fit,
					     TArrayI *volIDs,
					     TArrayI *volIDsFit,
					     TString misalignmentFile,
					     TString outname,
					     Int_t minPoints)
{
  //
  // This function process the AliTrackPoints and volID (into residuals) 
  //

  // setting up geometry and the trackpoints file
  if(!gGeoManager) AliGeomManager::LoadGeometry(GetFileNameGeometry()); 

  SetPointsFilename(GetFileNameTrackPoints());

  // creating some tools
  AliTrackFitter *fitter;
  if(fit==1){
    fitter = new AliTrackFitterKalman();
  }else fitter = new AliTrackFitterRieman();

  fitter->SetMinNPoints(minPoints);

  SetTrackFitter(fitter);

  if(misalignmentFile=="")printf("NO FAKE MISALIGNMENT\n");
  else {
    Bool_t misal=Misalign(misalignmentFile,"ITSAlignObjs");
    if(!misal){ 
      printf("PROBLEM WITH FAKE MISALIGNMENT!");
      return;
    }
  }

  CalculateResiduals(volIDs,volIDsFit,AliGeomManager::kSPD1,AliGeomManager::kSSD2,outname);

    return;

}
