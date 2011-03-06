/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliITSAlignMille
/// Alignment class for the ALICE ITS detector
///
/// ITS specific alignment class which interface to AliMillepede.   
/// For each track ProcessTrack calculates the local and global derivatives
/// at each hit and fill the corresponding local equations. Provide methods for
/// fixing or constraining detection elements for best results. 
///
/// \author M. Lunardon (thanks to J. Castillo)
//-----------------------------------------------------------------------------

#include <TF1.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TGraph.h>
#include <TMath.h>
#include <TGraphErrors.h>

#include "AliITSAlignMilleModule.h"
#include "AliITSAlignMille.h"
#include "AliITSAlignMilleData.h"
#include "AliITSgeomTGeo.h"
#include "AliGeomManager.h"
#include "AliMillepede.h"
#include "AliTrackPointArray.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"
#include <TSystem.h>
#include "AliTrackFitterRieman.h"

/// \cond CLASSIMP
ClassImp(AliITSAlignMille)
/// \endcond
  
Int_t AliITSAlignMille::fgNDetElem = ITSMILLENDETELEM;
Int_t AliITSAlignMille::fgNParCh = ITSMILLENPARCH;

AliITSAlignMille::AliITSAlignMille(const Char_t *configFilename, Bool_t initmille) 
  : TObject(),
    fMillepede(0),
    fStartFac(16.), 
    fResCutInitial(100.), 
    fResCut(100.),
    fNGlobal(ITSMILLENDETELEM*ITSMILLENPARCH),
    fNLocal(4),
    fNStdDev(ITSMILLENSTDEV),
    fIsMilleInit(kFALSE),
    fParSigTranslations(0.0100),
    fParSigRotations(0.1),
    fTrack(NULL),
    fCluster(),
    fGlobalDerivatives(NULL),
    fSigmaXfactor(1.0),
    fSigmaZfactor(1.0),
    fTempAlignObj(NULL),
    fDerivativeXLoc(0),
    fDerivativeZLoc(0),
    fMinNPtsPerTrack(3),
    fInitTrackParamsMeth(1),
    fProcessedPoints(NULL),
    fTotBadLocEqPoints(0),
    fRieman(NULL),
    fRequirePoints(kFALSE),
    fTempExcludedModule(-1),
    fGeometryFileName("geometry.root"),
    fPreAlignmentFileName(""),
    fGeoManager(0),
    fCurrentModuleIndex(0),
    fCurrentModuleInternalIndex(0),
    fCurrentSensVolIndex(0),
    fNModules(0),
    fUseLocalShifts(kTRUE),
    fUseSuperModules(kFALSE),
    fUsePreAlignment(kFALSE),
    fUseSortedTracks(kTRUE),
    fBOn(kFALSE),
    fBField(0.0),
    fNSuperModules(0),
    fCurrentModuleHMatrix(NULL),
    fIsConfigured(kFALSE),
    fBug(0)
{
  /// main constructor that takes input from configuration file
  
  fMillepede = new AliMillepede();
  fGlobalDerivatives = new Double_t[fNGlobal];

  for (Int_t i=0; i<ITSMILLENDETELEM*2; i++) {
    fPreAlignQF[i]=-1;
    fSensVolSigmaXfactor[i]=1.0;
    fSensVolSigmaZfactor[i]=1.0;
  }

  for (Int_t i=0; i<6; i++) {
    fNReqLayUp[i]=0;
    fNReqLayDown[i]=0;
    fNReqLay[i]=0;
  }
  for (Int_t i=0; i<3; i++) {
    fNReqDetUp[i]=0;
    fNReqDetDown[i]=0;
    fNReqDet[i]=0;
  }

  Int_t lc=LoadConfig(configFilename);
  if (lc) {
    AliInfo(Form("Error %d loading configuration from %s",lc,configFilename));
  }
  else {    
    fIsConfigured=kTRUE;
    if (initmille && fNGlobal<20000) {
      AliInfo(Form("Initializing Millepede with %d gpar, %d lpar and %d stddev ...",fNGlobal, fNLocal, fNStdDev));
      Init(fNGlobal, fNLocal, fNStdDev);      
      ResetLocalEquation();    
      AliInfo("Parameters initialized to zero");
    }
    else {
      AliInfo("Millepede has not been initialized ... ");
    }
  }
  
  if (fNModules) {
    fProcessedPoints=new Int_t[fNModules];
    for (Int_t ipp=0; ipp<fNModules; ipp++) fProcessedPoints[ipp]=0;
  }
}

AliITSAlignMille::~AliITSAlignMille() {
  /// Destructor
  if (fMillepede) delete fMillepede;
  delete [] fGlobalDerivatives;
  for (int i=0; i<fNModules; i++) delete fMilleModule[i];
  for (int i=0; i<fNSuperModules; i++) delete fSuperModule[i];
  if (fNModules) delete [] fProcessedPoints;
  if (fRieman) delete fRieman;
}

///////////////////////////////////////////////////////////////////////
Int_t AliITSAlignMille::LoadConfig(const Char_t *cfile) {
  /// return 0 if success
  ///        1 if error in module index or voluid
  
  FILE *pfc=fopen(cfile,"r");
  if (!pfc) return -1;
  
  Char_t st[200],st2[200];
  Char_t tmp[100];
  Int_t idx,itx,ity,itz,ith,ips,iph;
  Float_t f1,f2;
  UShort_t voluid;
  Int_t nmod=0;

  while (fgets(st,200,pfc)) {

    // skip comments
    for (int i=0; i<int(strlen(st)); i++) {
      if (st[i]=='#') st[i]=0;
    }

    if (strstr(st,"GEOMETRY_FILE")) {
      memset(tmp,0,100*sizeof(char));
      memset(st2,0,200*sizeof(char));
      sscanf(st,"%99s %199s",tmp,st2);
      if (gSystem->AccessPathName(st2)) {
	AliInfo("*** WARNING! *** geometry file not found! ");
	fclose(pfc);
	return -1;
      }  
      fGeometryFileName=st2;
      InitGeometry();
    }

    if (strstr(st,"PREALIGNMENT_FILE")) {
      memset(tmp,0,100*sizeof(char));
      memset(st2,0,200*sizeof(char));
      sscanf(st,"%99s %199s",tmp,st2);
      if (gSystem->AccessPathName(st2)) {
	AliInfo("*** WARNING! *** prealignment file not found! ");
	fclose(pfc);
	return -1;
      }  
      fPreAlignmentFileName=st2;
      itx=ApplyToGeometry();
      if (itx) {
	AliInfo(Form("*** WARNING! *** error %d reading prealignment file! ",itx));
	fclose(pfc);
	return -6;
      }
    }

    if (strstr(st,"SUPERMODULE_FILE")) {
      memset(tmp,0,100*sizeof(char));
      memset(st2,0,200*sizeof(char));
      sscanf(st,"%99s %199s",tmp,st2);
      if (gSystem->AccessPathName(st2)) {
	AliInfo("*** WARNING! *** supermodule file not found! ");
	fclose(pfc);
	return -1;
      }  
      if (LoadSuperModuleFile(st2)) {fclose(pfc); return -1;}
    }

    if (strstr(st,"SET_B_FIELD")) {
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %f",tmp,&f1);
      if (f1>0) {
	fBField = f1;
	fBOn = kTRUE;
	fNLocal = 5; // helices
	fRieman = new AliTrackFitterRieman();
      }  
      else {
	fBField = 0.0;
	fBOn = kFALSE;
	fNLocal = 4;
      }
    }

    if (strstr(st,"SET_PARSIG_TRA")) {
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %f",tmp,&f1);
      fParSigTranslations=f1;
    }

    if (strstr(st,"SET_PARSIG_ROT")) {
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %f",tmp,&f1);
      fParSigRotations=f1;
    }

    if (strstr(st,"SET_NSTDDEV")) {
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %d",tmp,&idx);
      fNStdDev=idx;
    }

    if (strstr(st,"SET_RESCUT_INIT")) {
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %f",tmp,&f1);
      fResCutInitial=f1;
    }

    if (strstr(st,"SET_RESCUT_OTHER")) {
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %f",tmp,&f1);
      fResCut=f1;
    }

    if (strstr(st,"SET_LOCALSIGMAFACTOR")) {
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %f %f",tmp,&f1,&f2);
      if (f1>0 && f2>0) {
	fSigmaXfactor=f1;
	fSigmaZfactor=f2;
      }
    }

    if (strstr(st,"SET_STARTFAC")) {
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %f",tmp,&f1);
      fStartFac=f1;
    }

    if (strstr(st,"REQUIRE_POINT")) {
      // syntax:   REQUIRE_POINT where ndet updw nreqpts
      //    where = LAYER or DETECTOR
      //    ndet = detector number: 1-6 for LAYER and 1-3 for DETECTOR (SPD=1, SDD=2, SSD=3)
      //    updw = 1 for Y>0, -1 for Y<0, 0 if not specified
      //    nreqpts = minimum number of points of that type
      memset(tmp,0,100*sizeof(char));
      memset(st2,0,200*sizeof(char));
      sscanf(st,"%99s %199s %d %d %d",tmp,st2,&itx,&ity,&itz);
      itx--;
      if (strstr(st2,"LAYER")) {
	if (itx<0 || itx>5) {fclose(pfc); return -7;}
	if (ity>0) fNReqLayUp[itx]=itz;
	else if (ity<0) fNReqLayDown[itx]=itz;
	else fNReqLay[itx]=itz;
	fRequirePoints=kTRUE;
      }
      else if (strstr(st2,"DETECTOR")) { // DETECTOR
	if (itx<0 || itx>2) {fclose(pfc); return -7;}
	if (ity>0) fNReqDetUp[itx]=itz;
	else if (ity<0) fNReqDetDown[itx]=itz;
	else fNReqDet[itx]=itz;	
	fRequirePoints=kTRUE;
      }
    }
    

    if (strstr(st,"SET_LOCAL_SHIFTS")) { // only enabled mode...
      fUseLocalShifts = kTRUE;
    }

    if (strstr(st,"MODULE_INDEX")) { // works only for sensitive modules
      f1=0; f2=0;
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %d %d %d %d %d %d %d %f %f",tmp,&idx,&itx,&ity,&itz,&iph,&ith,&ips,&f1,&f2);
      if (idx<0 || idx>2197) {fclose(pfc); return 1;} // bad index
      voluid=GetModuleVolumeID(idx);
      if (!voluid || voluid>14300) {fclose(pfc); return 1;} // bad index
      fModuleIndex[nmod]=idx;
      fModuleVolumeID[nmod]=voluid;
      fFreeParam[nmod][0]=itx;
      fFreeParam[nmod][1]=ity;
      fFreeParam[nmod][2]=itz;
      fFreeParam[nmod][3]=iph;
      fFreeParam[nmod][4]=ith;
      fFreeParam[nmod][5]=ips;
      fMilleModule[nmod] = new AliITSAlignMilleModule(voluid);
      if (f1>0) fSensVolSigmaXfactor[idx]=f1;
      if (f2>0) fSensVolSigmaZfactor[idx]=f2;
      nmod++;
    }
   
    if (strstr(st,"MODULE_VOLUID")) {
      f1=0; f2=0;
      memset(tmp,0,100*sizeof(char));
      sscanf(st,"%99s %d %d %d %d %d %d %d %f %f",tmp,&idx,&itx,&ity,&itz,&iph,&ith,&ips,&f1,&f2);
      voluid=UShort_t(idx);
      if (voluid>14335 && fUseSuperModules) { // custom supermodule
	int ism=-1;
	for (int j=0; j<fNSuperModules; j++) {
	  if (voluid==fSuperModule[j]->GetVolumeID()) ism=j;
	}
	if (ism<0) {fclose(pfc); return -1;} // bad volid
	fModuleIndex[nmod]=fSuperModule[ism]->GetIndex();
	fModuleVolumeID[nmod]=voluid;
	fFreeParam[nmod][0]=itx;
	fFreeParam[nmod][1]=ity;
	fFreeParam[nmod][2]=itz;
	fFreeParam[nmod][3]=iph;
	fFreeParam[nmod][4]=ith;
	fFreeParam[nmod][5]=ips;
	fMilleModule[nmod] = new AliITSAlignMilleModule(*fSuperModule[ism]);
	if (f1>0) {
	  for (int kk=0; kk<fMilleModule[nmod]->GetNSensitiveVolumes(); kk++) {
	    idx=AliITSAlignMilleModule::GetIndexFromVolumeID(fMilleModule[nmod]->GetSensitiveVolumeVolumeID()[kk]);
	    if (idx>=0) fSensVolSigmaXfactor[idx]=f1;
	  }
	}
	if (f2>0) {
	  for (int kk=0; kk<fMilleModule[nmod]->GetNSensitiveVolumes(); kk++) {
	    idx=AliITSAlignMilleModule::GetIndexFromVolumeID(fMilleModule[nmod]->GetSensitiveVolumeVolumeID()[kk]);
	    if (idx>=0) fSensVolSigmaZfactor[idx]=f2;
	  }
	}	nmod++;
      }
      else { // sensitive volume
	idx=GetModuleIndex(voluid);
	if (idx<0 || idx>2197) {fclose(pfc); return 1;} // bad index
	fModuleIndex[nmod]=idx;
	fModuleVolumeID[nmod]=voluid;
	fFreeParam[nmod][0]=itx;
	fFreeParam[nmod][1]=ity;
	fFreeParam[nmod][2]=itz;
	fFreeParam[nmod][3]=iph;
	fFreeParam[nmod][4]=ith;
	fFreeParam[nmod][5]=ips;
	fMilleModule[nmod] = new AliITSAlignMilleModule(voluid);
	if (f1>0) fSensVolSigmaXfactor[idx]=f1;
	if (f2>0) fSensVolSigmaZfactor[idx]=f2;
	nmod++;
      }
    }
    //----------

  } // end while

  fNModules = nmod;
  fNGlobal = fNModules*fgNParCh;
 
  fclose(pfc);
  return 0;
}

void AliITSAlignMille::SetRequiredPoint(Char_t* where, Int_t ndet, Int_t updw, Int_t nreqpts) 
{
  // set minimum number of points in specific detector or layer
  // where = LAYER or DETECTOR
  // ndet = detector number: 1-6 for LAYER and 1-3 for DETECTOR (SPD=1, SDD=2, SSD=3)
  // updw = 1 for Y>0, -1 for Y<0, 0 if not specified
  // nreqpts = minimum number of points of that type
  ndet--;
  if (strstr(where,"LAYER")) {
    if (ndet<0 || ndet>5) return;
    if (updw>0) fNReqLayUp[ndet]=nreqpts;
    else if (updw<0) fNReqLayDown[ndet]=nreqpts;
    else fNReqLay[ndet]=nreqpts;
    fRequirePoints=kTRUE;
  }
  else if (strstr(where,"DETECTOR")) {
    if (ndet<0 || ndet>2) return;
    if (updw>0) fNReqDetUp[ndet]=nreqpts;
    else if (updw<0) fNReqDetDown[ndet]=nreqpts;
    else fNReqDet[ndet]=nreqpts;	
    fRequirePoints=kTRUE;
  }
}

Int_t AliITSAlignMille::GetModuleIndex(const Char_t *symname) {
  /// index from symname
  if (!symname) return -1;
  for (Int_t i=0; i<2198; i++) {
    if (!strcmp(symname,AliITSgeomTGeo::GetSymName(i))) return i;
  }
  return -1;
}

Int_t AliITSAlignMille::GetModuleIndex(UShort_t voluid) {
  /// index from volume ID
  AliGeomManager::ELayerID lay = AliGeomManager::VolUIDToLayer(voluid);
  if (lay<1|| lay>6) return -1;
  Int_t idx=Int_t(voluid)-2048*lay;
  if (idx>=AliGeomManager::LayerSize(lay)) return -1;
  for (Int_t ilay=1; ilay<lay; ilay++) 
    idx += AliGeomManager::LayerSize(ilay);
  return idx;
}

UShort_t AliITSAlignMille::GetModuleVolumeID(const Char_t *symname) {
  /// volume ID from symname
  /// works for sensitive volumes only
  if (!symname) return 0;

  for (UShort_t voluid=2000; voluid<13300; voluid++) {
    Int_t modId;
    AliGeomManager::ELayerID layerId = AliGeomManager::VolUIDToLayer(voluid,modId);
    if (layerId>0 && layerId<7 && modId>=0 && modId<AliGeomManager::LayerSize(layerId)) {
      if (!strcmp(symname,AliGeomManager::SymName(layerId,modId))) return voluid;
    }
  }

  return 0;
}

UShort_t AliITSAlignMille::GetModuleVolumeID(Int_t index) {
  /// volume ID from index
  if (index<0) return 0;
  if (index<2198)
    return GetModuleVolumeID(AliITSgeomTGeo::GetSymName(index));
  else {
    for (int i=0; i<fNSuperModules; i++) {
      if (fSuperModule[i]->GetIndex()==index) return fSuperModule[i]->GetVolumeID();
    }
  }
  return 0;
}

void AliITSAlignMille::InitGeometry() {
  /// initialize geometry
  AliGeomManager::LoadGeometry(fGeometryFileName.Data());
  fGeoManager = AliGeomManager::GetGeometry();
  if (!fGeoManager) {
    AliInfo("Couldn't initialize geometry");
    return;
  }
  // temporary align object, just use the rotation...
  fTempAlignObj=new AliAlignObjParams;
}

void AliITSAlignMille::Init(Int_t nGlobal,  /* number of global paramers */
			   Int_t nLocal,   /* number of local parameters */
			   Int_t nStdDev   /* std dev cut */ )
{
  /// Initialization of AliMillepede. Fix parameters, define constraints ...
  fMillepede->InitMille(nGlobal,nLocal,nStdDev,fResCut,fResCutInitial);
  fIsMilleInit = kTRUE;
  
  /// Fix non free parameters
  for (Int_t i=0; i<fNModules; i++) {
    for (Int_t j=0; j<ITSMILLENPARCH; j++) {
      if (!fFreeParam[i][j]) FixParameter(i*ITSMILLENPARCH+j,0.0);
      else {
	// pepopepo: da verificare il settaggio delle sigma, ma forse va bene...
	Double_t parsig=0;
	if (j<3) parsig=fParSigTranslations; // translations (0.0100 cm)
	else parsig=fParSigRotations; // rotations (1/10 deg)
	FixParameter(i*ITSMILLENPARCH+j,parsig);
      }
    }    
  }
  
  
//   // Set iterations
  if (fStartFac>1) fMillepede->SetIterations(fStartFac);          
}


void AliITSAlignMille::AddConstraint(Double_t *par, Double_t value) {
  /// Constrain equation defined by par to value
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->SetGlobalConstraint(par, value);
  AliInfo("Adding constraint");
}

void AliITSAlignMille::InitGlobalParameters(Double_t *par) {
  /// Initialize global parameters with par array
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->SetGlobalParameters(par);
  AliInfo("Init Global Parameters");
}
 
void AliITSAlignMille::FixParameter(Int_t iPar, Double_t value) {
  /// Parameter iPar is encourage to vary in [-value;value]. 
  /// If value == 0, parameter is fixed
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->SetParSigma(iPar, value);
  if (value==0) AliInfo(Form("Parameter %i Fixed", iPar));
}

void AliITSAlignMille::ResetLocalEquation()
{
  /// Reset the derivative vectors
  for(int i=0; i<fNLocal; i++) {
    fLocalDerivatives[i] = 0.0;
  }
  for(int i=0; i<fNGlobal; i++) {
    fGlobalDerivatives[i] = 0.0;
  }
}

// newpep
Int_t AliITSAlignMille::ApplyToGeometry() {
  /// apply starting realignment to ideal geometry
  if(!AliGeomManager::GetGeometry()) return -1;

  TFile *pref = new TFile(fPreAlignmentFileName.Data());
  if (!pref->IsOpen()) return -2;
  TClonesArray *prea=(TClonesArray*)pref->Get("ITSAlignObjs");
  if (!prea) return -3;  
  Int_t nprea=prea->GetEntriesFast();
  AliInfo(Form("Array of input misalignments with %d entries",nprea));

  AliGeomManager::ApplyAlignObjsToGeom(*prea); // apply all levels of objs

  // set prealignment factor if defined...
  for (int ix=0; ix<nprea; ix++) {
    AliAlignObjParams *preo=(AliAlignObjParams*) prea->UncheckedAt(ix);
    Int_t index=AliITSAlignMilleModule::GetIndexFromVolumeID(preo->GetVolUID());
    if (index>=0) {
      fPreAlignQF[index] = (int) preo->GetUniqueID();
      //printf("index=%d   QF=%d\n",index,preo->GetUniqueID());
    }
    //if (!preo->ApplyToGeometry()) return -4;
  }
  pref->Close();
  delete pref;

  fUsePreAlignment = kTRUE;
  return 0;
}
// endnewpep

Int_t AliITSAlignMille::GetPreAlignmentQualityFactor(Int_t index) const {
  /// works for sensitive volumes
  if (!fUsePreAlignment || index<0 || index>2197) return -1;
  return fPreAlignQF[index];
}

AliTrackPointArray *AliITSAlignMille::PrepareTrack(const AliTrackPointArray *atp) {
  /// create a new AliTrackPointArray keeping only defined modules
  /// move points according to a given prealignment, if any
  /// sort alitrackpoints w.r.t. global Y direction, if selected

  AliTrackPointArray *atps=NULL;
  Int_t idx[20];
  Int_t npts=atp->GetNPoints();

  /// checks if AliTrackPoints belong to defined modules
  Int_t ngoodpts=0;
  Int_t intidx[20];
  
  for (int j=0; j<npts; j++) {
    intidx[j] = IsContained(atp->GetVolumeID()[j]);
    if (intidx[j]>=0) ngoodpts++;
  }
  AliDebug(3,Form("Number of points in defined modules: %d",ngoodpts));

  // reject track if not enough points are left
  if (ngoodpts<fMinNPtsPerTrack) {
    AliInfo("Track with not enough points!");
    return NULL;
  }

  AliTrackPoint p;
  // check points in specific places
  if (fRequirePoints) {
    Int_t nlayup[6],nlaydown[6],nlay[6];
    Int_t ndetup[3],ndetdown[3],ndet[3];
    for (Int_t j=0; j<6; j++) {nlayup[j]=0; nlaydown[j]=0; nlay[j]=0;}
    for (Int_t j=0; j<3; j++) {ndetup[j]=0; ndetdown[j]=0; ndet[j]=0;}
    
    for (int i=0; i<npts; i++) {
      // skip not defined points
      if (intidx[i]<0) continue;
      Float_t xx=atp->GetX()[i];
      Float_t yy=atp->GetY()[i];
      Float_t r=TMath::Sqrt(xx*xx + yy*yy);
      int lay=-1;
      if (r<5) lay=0;
      else if (r>5 && r<10) lay=1;
      else if (r>10 && r<18) lay=2;
      else if (r>18 && r<30) lay=3;
      else if (r>30 && r<40) lay=4;
      else if (r>40) lay=5;
      if (lay<0) continue;
      int det=lay/2;
      //printf("Point %d - x=%f  y=%f  R=%f  lay=%d  det=%d\n",i,xx,yy,r,lay,det);

      if (yy>=0.0) { // UP point
	nlayup[lay]++;
	nlay[lay]++;
	ndetup[det]++;
	ndet[det]++;
      }
      else {
	nlaydown[lay]++;
	nlay[lay]++;
	ndetdown[det]++;
	ndet[det]++;
      }
    }
    
    // checks minimum values
    Bool_t isok=kTRUE;
    for (Int_t j=0; j<6; j++) {
      if (nlayup[j]<fNReqLayUp[j]) isok=kFALSE; 
      if (nlaydown[j]<fNReqLayDown[j]) isok=kFALSE; 
      if (nlay[j]<fNReqLay[j]) isok=kFALSE; 
    }
    for (Int_t j=0; j<3; j++) {
      if (ndetup[j]<fNReqDetUp[j]) isok=kFALSE; 
      if (ndetdown[j]<fNReqDetDown[j]) isok=kFALSE; 
      if (ndet[j]<fNReqDet[j]) isok=kFALSE; 
    }
    if (!isok) {
      AliDebug(2,Form("Track does not meet all location point requirements!"));
      return NULL;
    }
  }
  
  // build a new track with (sorted) (prealigned) good points
  atps=new AliTrackPointArray(ngoodpts);

  for (int i=0; i<npts; i++) idx[i]=i;
  // sort track if required
  if (fUseSortedTracks) TMath::Sort(npts,atp->GetY(),idx); // sort descending...

  Int_t npto=0;
  for (int i=0; i<npts; i++) {
    // skip not defined points
    if (intidx[idx[i]]<0) continue;
    atp->GetPoint(p,idx[i]);

    // prealign point if required
    // get IDEAL matrix
    TGeoHMatrix *svOrigMatrix = fMilleModule[intidx[idx[i]]]->GetSensitiveVolumeOrigGlobalMatrix(p.GetVolumeID());
    // get back real local coordinates: use OriginalGlobalMatrix because AliTrackPoints were written
    // with idel geometry  
    Double_t pg[3],pl[3];
    pg[0]=p.GetX();
    pg[1]=p.GetY();
    pg[2]=p.GetZ();
    AliDebug(3,Form("Global coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]));
    svOrigMatrix->MasterToLocal(pg,pl);

    AliDebug(3,Form("Local coordinates of measured point : X=%f  Y=%f  Z=%f \n",pl[0],pl[1],pl[2]));

    // update covariance matrix
    TGeoHMatrix hcov;
    Double_t hcovel[9];
    hcovel[0]=double(p.GetCov()[0]);
    hcovel[1]=double(p.GetCov()[1]);
    hcovel[2]=double(p.GetCov()[2]);
    hcovel[3]=double(p.GetCov()[1]);
    hcovel[4]=double(p.GetCov()[3]);
    hcovel[5]=double(p.GetCov()[4]);
    hcovel[6]=double(p.GetCov()[2]);
    hcovel[7]=double(p.GetCov()[4]);
    hcovel[8]=double(p.GetCov()[5]);
    hcov.SetRotation(hcovel);
    // now rotate in local system
    hcov.Multiply(svOrigMatrix);
    hcov.MultiplyLeft(&svOrigMatrix->Inverse());
    // now hcov is LOCAL COVARIANCE MATRIX


    // pepopepo
    if (fBug==1) {
      // correzione bug LAYER 5  SSD temporanea..
      int ssdidx=AliITSAlignMilleModule::GetIndexFromVolumeID(p.GetVolumeID());
      if (ssdidx>=500 && ssdidx<1248) {
	int ladder=(ssdidx-500)%22;
      if (ladder==18) p.SetVolumeID(AliITSAlignMilleModule::GetVolumeIDFromIndex(ssdidx+1));
      if (ladder==19) p.SetVolumeID(AliITSAlignMilleModule::GetVolumeIDFromIndex(ssdidx-1));
      }
    }

    /// get (evenctually prealigned) matrix of sens. vol.
    TGeoHMatrix *svMatrix = fMilleModule[intidx[idx[i]]]->GetSensitiveVolumeMatrix(p.GetVolumeID());
    // modify global coordinates according with pre-aligment
    svMatrix->LocalToMaster(pl,pg);
    // now rotate in local system
    hcov.Multiply(&svMatrix->Inverse());
    hcov.MultiplyLeft(svMatrix);
    // hcov is back in GLOBAL RF
    Float_t pcov[6];
    pcov[0]=hcov.GetRotationMatrix()[0];
    pcov[1]=hcov.GetRotationMatrix()[1];
    pcov[2]=hcov.GetRotationMatrix()[2];
    pcov[3]=hcov.GetRotationMatrix()[4];
    pcov[4]=hcov.GetRotationMatrix()[5];
    pcov[5]=hcov.GetRotationMatrix()[8];

    p.SetXYZ(pg[0],pg[1],pg[2],pcov);
    AliDebug(3,Form("New global coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]));
    atps->AddPoint(npto,&p);
    AliDebug(2,Form("Adding point[%d] = ( %f , %f , %f )     volid = %d",npto,atps->GetX()[npto],atps->GetY()[npto],atps->GetZ()[npto],atps->GetVolumeID()[npto] ));

    npto++;
  }

  return atps;
}



AliTrackPointArray *AliITSAlignMille::SortTrack(const AliTrackPointArray *atp) {
  /// sort alitrackpoints w.r.t. global Y direction
  AliTrackPointArray *atps=NULL;
  Int_t idx[20];
  Int_t npts=atp->GetNPoints();
  AliTrackPoint p;
  atps=new AliTrackPointArray(npts);

  TMath::Sort(npts,atp->GetY(),idx);

  for (int i=0; i<npts; i++) {
    atp->GetPoint(p,idx[i]);
    atps->AddPoint(i,&p);
    AliDebug(2,Form("Point[%d] = ( %f , %f , %f )     volid = %d",i,atps->GetX()[i],atps->GetY()[i],atps->GetZ()[i],atps->GetVolumeID()[i] ));
  }
  return atps;
}


Int_t AliITSAlignMille::InitModuleParams() {
  /// initialize geometry parameters for a given detector
  /// for current cluster (fCluster)
  /// fGlobalInitParam[] is set as:
  ///    [tx,ty,tz,psi,theta,phi]
  ///
  /// return 0 if success

  if (!fGeoManager) {
    AliInfo("ITS geometry not initialized!");
    return -1;
  }

  // now 'voluid' is the volumeID of a SENSITIVE VOLUME (coming from a cluster)

  // set the internal index (index in module list)
  UShort_t voluid=fCluster.GetVolumeID();
  Int_t k=fNModules-1;
  while (k>=0 && !(fMilleModule[k]->IsIn(voluid)) ) k--; 
  if (k<0) return -3;    
  fCurrentModuleInternalIndex=k; // the internal index of the SUPERMODULE

  fCurrentModuleIndex=fMilleModule[k]->GetIndex(); // index of the SUPERMODULE

  // set the index
  Int_t index = GetModuleIndex(voluid);
  if (index<0) return -2;
  fCurrentSensVolIndex = index; // the index of the SENSITIVE VOLUME

  fModuleInitParam[0] = 0.0;
  fModuleInitParam[1] = 0.0;
  fModuleInitParam[2] = 0.0;
  fModuleInitParam[3] = 0.0; // psi   (X)
  fModuleInitParam[4] = 0.0; // theta (Y)
  fModuleInitParam[5] = 0.0; // phi   (Z)
  
  fCurrentModuleHMatrix = fMilleModule[fCurrentModuleInternalIndex]->GetMatrix();

  for (int ii=0; ii<3; ii++)
    fCurrentModuleTranslation[ii]=fCurrentModuleHMatrix->GetTranslation()[ii];

  /// get (evenctually prealigned) matrix of sens. vol.
  TGeoHMatrix *svMatrix = fMilleModule[fCurrentModuleInternalIndex]->GetSensitiveVolumeMatrix(voluid);
  
  fMeasGlo[0] = fCluster.GetX();
  fMeasGlo[1] = fCluster.GetY();
  fMeasGlo[2] = fCluster.GetZ(); 
  svMatrix->MasterToLocal(fMeasGlo,fMeasLoc);  
  AliDebug(2,Form("Local coordinates of measured point : X=%f  Y=%f  Z=%f \n",fMeasLoc[0] ,fMeasLoc[1] ,fMeasLoc[2] ));
  
  // set stdev from cluster
  TGeoHMatrix hcov;
  Double_t hcovel[9];
  hcovel[0]=double(fCluster.GetCov()[0]);
  hcovel[1]=double(fCluster.GetCov()[1]);
  hcovel[2]=double(fCluster.GetCov()[2]);
  hcovel[3]=double(fCluster.GetCov()[1]);
  hcovel[4]=double(fCluster.GetCov()[3]);
  hcovel[5]=double(fCluster.GetCov()[4]);
  hcovel[6]=double(fCluster.GetCov()[2]);
  hcovel[7]=double(fCluster.GetCov()[4]);
  hcovel[8]=double(fCluster.GetCov()[5]);
  hcov.SetRotation(hcovel);
  // now rotate in local system
  hcov.Multiply(svMatrix);
  hcov.MultiplyLeft(&svMatrix->Inverse());

  // set local sigmas
  fSigmaLoc[0] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[0]));
  fSigmaLoc[1] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[4]));
  fSigmaLoc[2] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[8]));

  // set minimum value for SigmaLoc to 10 micron 
  if (fSigmaLoc[0]<0.0010) fSigmaLoc[0]=0.0010;
  if (fSigmaLoc[2]<0.0010) fSigmaLoc[2]=0.0010;

  // multiply local sigmas by global factor
  fSigmaLoc[0] *= fSigmaXfactor;
  fSigmaLoc[2] *= fSigmaZfactor;

  // multiply local sigmas by individual factor
  fSigmaLoc[0] *= fSensVolSigmaXfactor[index];
  fSigmaLoc[2] *= fSensVolSigmaZfactor[index];

  AliDebug(2,Form("Setting StDev from CovMat : fSigmaLocX=%g  fSigmaLocY=%g fSigmaLocZ=%g \n",fSigmaLoc[0] ,fSigmaLoc[1] ,fSigmaLoc[2] ));
   
  return 0;
}

void AliITSAlignMille::SetCurrentModule(Int_t index) {
  /// set as current the SuperModule that contains the 'index' sens.vol.
  if (index<0 || index>2197) {
    AliInfo("index does not correspond to a sensitive volume!");
    return;
  }
  UShort_t voluid=GetModuleVolumeID(index);
  Int_t k=IsContained(voluid);
  if (k>=0){
    fCluster.SetVolumeID(voluid);
    fCluster.SetXYZ(0,0,0);
    InitModuleParams();
  }
  else
    AliInfo(Form("module %d not defined\n",index));    
}

void AliITSAlignMille::SetCurrentSensitiveModule(Int_t index) {
  /// set as current the SuperModule that contains the 'index' sens.vol.
  if (index<0 || index>2197) {
    AliInfo("index does not correspond to a sensitive volume!");
    return;
  }
  UShort_t voluid=AliITSAlignMilleModule::GetVolumeIDFromIndex(index);
  Int_t k=IsDefined(voluid);
  if (k>=0){
    fCluster.SetVolumeID(voluid);
    fCluster.SetXYZ(0,0,0);
    InitModuleParams();
  }
  else
    AliInfo(Form("module %d not defined\n",index));    
}

void AliITSAlignMille::Print(Option_t*) const 
{
  /// print infos
  printf("*** AliMillepede for ITS ***\n");
  printf("    number of defined super modules: %d\n",fNModules);
  
  if (fGeoManager)
    printf("    geometry loaded from %s\n",fGeometryFileName.Data());
  else
    printf("    geometry not loaded\n");
  
  if (fUseSuperModules) 
    printf("    using custom supermodules ( %d defined )\n",fNSuperModules);
  else
    printf("    custom supermodules not used\n");    

  if (fUsePreAlignment) 
    printf("    using prealignment from %s \n",fPreAlignmentFileName.Data());
  else
    printf("    prealignment not used\n");    

  if (fBOn) 
    printf("    B Field set to %f T - using Riemann's helices\n",fBField);
  else
    printf("    B Field OFF - using straight lines \n");

  if (fUseLocalShifts) 
    printf("    Alignment shifts will be computed in LOCAL RS\n");
  else
    printf("    Alignment shifts will be computed in GLOBAL RS\n");    

  if (fRequirePoints) printf("    Required points in tracks:\n");
  for (Int_t i=0; i<6; i++) {
    if (fNReqLayUp[i]>0) printf("        Layer %d : %d points with Y>0\n",i+1,fNReqLayUp[i]);
    if (fNReqLayDown[i]>0) printf("        Layer %d : %d points with Y<0\n",i+1,fNReqLayDown[i]);
    if (fNReqLay[i]>0) printf("        Layer %d : %d points \n",i+1,fNReqLay[i]);
  }
  for (Int_t i=0; i<3; i++) {
    if (fNReqDetUp[i]>0) printf("        Detector %d : %d points with Y>0\n",i+1,fNReqDetUp[i]);
    if (fNReqDetDown[i]>0) printf("        Detector %d : %d points with Y<0\n",i+1,fNReqDetDown[i]);
    if (fNReqDet[i]>0) printf("        Detector %d : %d points \n",i+1,fNReqDet[i]);
  }
  
  printf("\n    Millepede configuration parameters:\n");
  printf("        init parsig for translations  : %.4f\n",fParSigTranslations);
  printf("        init parsig for rotations     : %.4f\n",fParSigRotations);
  printf("        init value for chi2 cut       : %.4f\n",fStartFac);
  printf("        first iteration cut value     : %.4f\n",fResCutInitial);
  printf("        other iterations cut value    : %.4f\n",fResCut);
  printf("        number of stddev for chi2 cut : %d\n",fNStdDev);
  printf("        mult. fact. for local sigmas  : %.4f %.4f\n",fSigmaXfactor, fSigmaZfactor);

  printf("List of defined modules:\n");
  printf("  intidx\tindex\tvoluid\tname\n");
  for (int i=0; i<fNModules; i++)
    printf("  %d\t%d\t%d\t%s\n",i,fModuleIndex[i],fModuleVolumeID[i],fMilleModule[i]->GetName());
   
}

AliITSAlignMilleModule  *AliITSAlignMille::GetMilleModule(UShort_t voluid) const
{
  // return pointer to a define supermodule
  // return NULL if error
  Int_t i=IsDefined(voluid);
  if (i<0) return NULL;
  return fMilleModule[i];
}

AliITSAlignMilleModule *AliITSAlignMille::GetCurrentModule() const
{
  if (fNModules) return fMilleModule[fCurrentModuleInternalIndex];
  return NULL;
}

void AliITSAlignMille::PrintCurrentModuleInfo() 
{
  ///
  Int_t k=fCurrentModuleInternalIndex;
  if (k<0 || k>=fNModules) return;
  fMilleModule[k]->Print("");
}

Bool_t AliITSAlignMille::InitRiemanFit() {
  // Initialize Riemann Fitter for current track
  // return kFALSE if error

  if (!fBOn) return kFALSE;

  Int_t npts=0;
  AliTrackPoint ap;
  npts = fTrack->GetNPoints();
  AliDebug(3,Form("Fitting track with %d points",npts));

  fRieman->Reset();
  fRieman->SetTrackPointArray(fTrack);

  TArrayI ai(npts);
  for (Int_t ipt=0; ipt<npts; ipt++) ai[ipt]=fTrack->GetVolumeID()[ipt];
  
  // fit track with 5 params in his own tracking-rotated reference system
  // xc = -p[1]/p[0];
  // yc = 1/p[0];
  // R  = sqrt( x0*x0 + y0*y0 - y0*p[2]);
  if (!fRieman->Fit(&ai,NULL,(AliGeomManager::ELayerID)1,(AliGeomManager::ELayerID)6)) {
    return kFALSE;
  }

  for (int i=0; i<5; i++)
    fLocalInitParam[i] = fRieman->GetParam()[i];
  
  return kTRUE;
}


void AliITSAlignMille::InitTrackParams(int meth) {
  /// initialize local parameters with different methods
  /// for current track (fTrack)
  
  Int_t npts=0;
  TF1 *f1=NULL;
  TGraph *g=NULL;
  Float_t sigmax[20],sigmay[20],sigmaz[20];
  AliTrackPoint ap;
  TGraphErrors *ge=NULL;

  switch (meth) {
  case 1:   // simple linear interpolation
    // get local starting parameters (to be substituted by ESD track parms)
    // local parms (fLocalInitParam[]) are:
    //      [0] = global x coord. of straight line intersection at y=0 plane
    //      [1] = global z coord. of straight line intersection at y=0 plane
    //      [2] = px/py  
    //      [3] = pz/py
    
    // test #1: linear fit in x(y) and z(y)
    npts = fTrack->GetNPoints();
    AliDebug(3,Form("*** initializing track with %d points ***",npts));

    f1=new TF1("f1","[0]+x*[1]",-50,50);

    g=new TGraph(npts,fTrack->GetY(),fTrack->GetX());
    g->Fit(f1,"RNQ");
    fLocalInitParam[0] = f1->GetParameter(0);
    fLocalInitParam[2] = f1->GetParameter(1);
    AliDebug(2,Form("X = p0gx + ugx*Y : p0gx = %f +- %f    ugx = %f +- %f\n",fLocalInitParam[0],f1->GetParError(0),fLocalInitParam[2],f1->GetParError(1)));
    delete g; g=NULL;

    g=new TGraph(npts,fTrack->GetY(),fTrack->GetZ());
    g->Fit(f1,"RNQ");
    fLocalInitParam[1] = f1->GetParameter(0);
    fLocalInitParam[3] = f1->GetParameter(1);
    AliDebug(2,Form("Z = p0gz + ugz*Y : p0gz=%f  ugz=%f\n",fLocalInitParam[1],fLocalInitParam[3]));
    delete g;
    delete f1;

    break;
    
  case 2:   // simple linear interpolation weighted using sigmas
    // get local starting parameters (to be substituted by ESD track parms)
    // local parms (fLocalInitParam[]) are:
    //      [0] = global x coord. of straight line intersection at y=0 plane
    //      [1] = global z coord. of straight line intersection at y=0 plane
    //      [2] = px/py  
    //      [3] = pz/py
    
    // test #1: linear fit in x(y) and z(y)
    npts = fTrack->GetNPoints();
    AliDebug(3,Form("*** initializing track with %d points ***",npts));

    for (Int_t isig=0; isig<npts; isig++) {
      fTrack->GetPoint(ap,isig);
      sigmax[isig]=ap.GetCov()[0]; 
      if (sigmax[isig]<1.0e-07) sigmax[isig]=1.0e-07; // minimum sigma=3 mu
      sigmax[isig]=TMath::Sqrt(sigmax[isig]);

      sigmay[isig]=ap.GetCov()[3]; 
      if (sigmay[isig]<1.0e-07) sigmay[isig]=1.0e-07; // minimum sigma=3 mu
      sigmay[isig]=TMath::Sqrt(sigmay[isig]);

      sigmaz[isig]=ap.GetCov()[5]; 
      if (sigmaz[isig]<1.0e-07) sigmaz[isig]=1.0e-07; // minimum sigma=3 mu
      sigmaz[isig]=TMath::Sqrt(sigmaz[isig]);      

      if (fTempExcludedModule>=0 && fTempExcludedModule==AliITSAlignMilleModule::GetIndexFromVolumeID(ap.GetVolumeID())) {
	sigmax[isig] *= 1000.;
	sigmay[isig] *= 1000.;
	sigmaz[isig] *= 1000.;
	fTempExcludedModule=-1;
      }
    }

    f1=new TF1("f1","[0]+x*[1]",-50,50);

    ge=new TGraphErrors(npts,fTrack->GetY(),fTrack->GetX(),sigmay,sigmax);
    ge->Fit(f1,"RNQ");
    fLocalInitParam[0] = f1->GetParameter(0);
    fLocalInitParam[2] = f1->GetParameter(1);
    AliDebug(2,Form("X = p0gx + ugx*Y : p0gx = %f +- %f    ugx = %f +- %f\n",fLocalInitParam[0],f1->GetParError(0),fLocalInitParam[2],f1->GetParError(1)));
    delete ge; ge=NULL;
    
    ge=new TGraphErrors(npts,fTrack->GetY(),fTrack->GetZ(),sigmay,sigmaz);
    ge->Fit(f1,"RNQ");
    fLocalInitParam[1] = f1->GetParameter(0);
    fLocalInitParam[3] = f1->GetParameter(1);
    AliDebug(2,Form("Z = p0gz + ugz*Y : p0gz=%f  ugz=%f\n",fLocalInitParam[1],fLocalInitParam[3]));
    delete ge;
    delete f1;
    
    break;
    
  }
}

Int_t AliITSAlignMille::IsDefined(UShort_t voluid) const
{
  // checks if supermodule 'voluid' is defined and return the internal index
  // return -1 if error
  Int_t k=fNModules-1;
  while (k>=0 && !(voluid==fModuleVolumeID[k]) ) k--;  
  if (k<0) return -1; 
  return k;
}

Int_t AliITSAlignMille::IsContained(UShort_t voluid) const
{
  // checks if the sensitive module 'voluid' is contained inside a supermodule and return the internal index of the last identified supermodule
  // return -1 if error
  if (AliITSAlignMilleModule::GetIndexFromVolumeID(voluid)<0) return -1;
  Int_t k=fNModules-1;
  while (k>=0 && !(fMilleModule[k]->IsIn(voluid)) ) k--;  
  if (k<0) return -1; 
  return k;
}

Bool_t AliITSAlignMille::CheckVolumeID(UShort_t voluid) const 
{
  /// check if a sensitive volume is contained inside one of the defined supermodules
  Int_t k=fNModules-1;
  while (k>=0 && !(fMilleModule[k]->IsIn(voluid)) ) k--;  
  if (k>=0) return kTRUE;
  return kFALSE;
}

Int_t AliITSAlignMille::CheckCurrentTrack() {
  /// checks if AliTrackPoints belongs to defined modules
  /// return number of good poins
  /// return 0 if not enough points

  Int_t npts = fTrack->GetNPoints();
  Int_t ngoodpts=0;
  // debug points
  for (int j=0; j<npts; j++) {
    UShort_t voluid = fTrack->GetVolumeID()[j];    
    if (CheckVolumeID(voluid)) {
      ngoodpts++;
    }
  }

  if (ngoodpts<fMinNPtsPerTrack) return 0;

  return ngoodpts;
}

Int_t AliITSAlignMille::ProcessTrack(AliTrackPointArray *track) {
  /// Process track; Loop over hits and set local equations
  /// here 'track' is a AliTrackPointArray
  /// return 0 if success;
  
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return -1;
  }

  Int_t npts = track->GetNPoints();
  AliDebug(2,Form("*** Input track with %d points ***",npts));

  // preprocessing of the input track: keep only points in defined volumes,
  // move points if prealignment is set, sort by Yglo if required
  
  fTrack=PrepareTrack(track);
  if (!fTrack) return -1;

  npts = fTrack->GetNPoints();
  AliDebug(2,Form("*** Processing prepared track with %d points ***",npts));
  
  if (!fBOn) { // straight lines  
    // set local starting parameters (to be substituted by ESD track parms)
    // local parms (fLocalInitParam[]) are:
    //      [0] = global x coord. of straight line intersection at y=0 plane
    //      [1] = global z coord. of straight line intersection at y=0 plane
    //      [2] = px/py  
    //      [3] = pz/py
    InitTrackParams(fInitTrackParamsMeth);  
  } 
  else {
    // local parms (fLocalInitParam[]) are the Riemann Fitter params
    if (!InitRiemanFit()) {
      AliInfo("Riemann fit failed! skipping this track...");
      delete fTrack;
      fTrack=NULL;
      return -5;
    }
  }
  
  Int_t nloceq=0;
  AliITSAlignMilleData *md = new AliITSAlignMilleData[npts];
  
  for (Int_t ipt=0; ipt<npts; ipt++) {
    fTrack->GetPoint(fCluster,ipt);
    AliDebug(2,Form("\n--- processing point %d --- \n",ipt));    

    // set geometry parameters for the the current module
    if (InitModuleParams()) continue;
    AliDebug(2,Form("    VolID=%d  Index=%d  InternalIdx=%d  symname=%s\n", track->GetVolumeID()[ipt], fCurrentModuleIndex ,fCurrentModuleInternalIndex, AliGeomManager::SymName(track->GetVolumeID()[ipt]) ));
    AliDebug(2,Form("    Preprocessed Point = ( %f , %f , %f ) \n",fCluster.GetX(),fCluster.GetY(),fCluster.GetZ()));
    
    if (!AddLocalEquation(md[nloceq])) {
      nloceq++;    
      fProcessedPoints[fCurrentModuleInternalIndex]++;
    }
    else {
      fTotBadLocEqPoints++;
    }
    
  } // end loop over points
  
  delete fTrack;
  fTrack=NULL;

  // not enough good points!
  if (nloceq<fMinNPtsPerTrack) {
    delete [] md;      
    return -1;
  }
  
  // finally send local equations to millepede
  SetLocalEquations(md,nloceq);

  delete [] md;
  
  return 0;
}

Int_t AliITSAlignMille::CalcIntersectionPoint(const Double_t *lpar, const Double_t *gpar) {
  /// calculate intersection point of track with current module in local coordinates
  /// according with a given set of parameters (local(4/5) and global(6))
  /// and fill fPintLoc/Glo
  ///    local are:   pgx0, pgz0, ugx, ugz   OR   riemann fitters pars
  ///    global are:  tx,ty,tz,psi,theta,phi (Raff's delta angles in deg.)
  /// return 0 if success
  
  AliDebug(3,Form("lpar = %g %g %g %g %g\ngpar= %g %g %g %g %g %g\n",lpar[0],lpar[1],lpar[2],lpar[3],lpar[4],gpar[0],gpar[1],gpar[2],gpar[3],gpar[4],gpar[5]));
  AliDebug(3,Form("deltalpar = %g %g %g %g %g\n",lpar[0]-fLocalInitParam[0],lpar[1]-fLocalInitParam[1],lpar[2]-fLocalInitParam[2],lpar[3]-fLocalInitParam[3],lpar[4]-fLocalInitParam[4]));

  
  // prepare the TGeoHMatrix
  TGeoHMatrix *fTempHMat = fMilleModule[fCurrentModuleInternalIndex]->GetSensitiveVolumeModifiedMatrix(fCluster.GetVolumeID(),gpar);
  if (!fTempHMat) return -1;
  
  Double_t v0g[3]; // vector with straight line direction in global coord.
  Double_t p0g[3]; // point of the straight line (glo)
  
  if (fBOn) { // B FIELD!
    AliTrackPoint prf; 
    for (int ip=0; ip<5; ip++)
      fRieman->SetParam(ip,lpar[ip]);

    if (!fRieman->GetPCA(fCluster,prf))  {
      AliInfo(Form("error in GetPCA for point %d",fCluster.GetVolumeID()));
      return -3;
    }
    // now determine straight line passing tangent to fit curve at prf
    // ugx = dX/dY_glo = DeltaX/DeltaY_glo
    // mo' P1=(X,Y,Z)_glo_prf
    //       => (x,y,Z)_trk_prf ruotando di alpha...
    Double_t alpha=fRieman->GetAlpha();
    Double_t x1g = prf.GetX();
    Double_t y1g = prf.GetY();
    Double_t z1g = prf.GetZ();
    Double_t x1t =  x1g*TMath::Cos(alpha) + y1g*TMath::Sin(alpha);
    Double_t y1t = -x1g*TMath::Sin(alpha) + y1g*TMath::Cos(alpha);
    Double_t z1t =  z1g;    

    Double_t x2t = x1t+1.0;
    Double_t y2t = y1t+fRieman->GetDYat(x1t);
    Double_t z2t = z1t+fRieman->GetDZat(x1t);
    Double_t x2g =  x2t*TMath::Cos(alpha) - y2t*TMath::Sin(alpha);
    Double_t y2g =  x2t*TMath::Sin(alpha) + y2t*TMath::Cos(alpha);
    Double_t z2g =  z2t;  

    AliDebug(3,Form("Riemann frame:  fAlpha = %f  =  %f  ",alpha,alpha*180./TMath::Pi()));
    AliDebug(3,Form("   prf_glo=( %f , %f , %f )  prf_rf=( %f , %f , %f )\n", x1g,y1g,z1g, x1t,y1t,z1t));
    AliDebug(3,Form("   mov_glo=( %f , %f , %f )      rf=( %f , %f , %f )\n",x2g,y2g,z2g, x2t,y2t,z2t));
        
    if (TMath::Abs(y2g-y1g)<1e-15) {
      AliInfo("DeltaY=0! Cannot proceed...");
      return -1;
    }
    // ugx,1,ugz
    v0g[0] = (x2g-x1g)/(y2g-y1g);
    v0g[1]=1.0;
    v0g[2] = (z2g-z1g)/(y2g-y1g);
    
    // point: just keep prf
    p0g[0]=x1g;
    p0g[1]=y1g;
    p0g[2]=z1g;
  }  
  else { // staight line
    // vector of initial straight line direction in glob. coord
    v0g[0]=lpar[2];
    v0g[1]=1.0;
    v0g[2]=lpar[3];
    
    // intercept in yg=0 plane in glob coord
    p0g[0]=lpar[0];
    p0g[1]=0.0;
    p0g[2]=lpar[1];
  }
  AliDebug(3,Form("Line vector: ( %f , %f , %f )  point:( %f , %f , %f )\n",v0g[0],v0g[1],v0g[2],p0g[0],p0g[1],p0g[2]));
  
  // same in local coord.
  Double_t p0l[3],v0l[3];
  fTempHMat->MasterToLocalVect(v0g,v0l);
  fTempHMat->MasterToLocal(p0g,p0l);
  
  if (TMath::Abs(v0l[1])<1e-15) {
    AliInfo("Track Y direction in local frame is zero! Cannot proceed...");
    return -1;
  }
  
  // local intersection point
  fPintLoc[0] = p0l[0] - (v0l[0]/v0l[1])*p0l[1];
  fPintLoc[1] = 0;
  fPintLoc[2] = p0l[2] - (v0l[2]/v0l[1])*p0l[1];
  
  // global intersection point
  fTempHMat->LocalToMaster(fPintLoc,fPintGlo);
  AliDebug(3,Form("Intesect. point: L( %f , %f , %f )  G( %f , %f , %f )\n",fPintLoc[0],fPintLoc[1],fPintLoc[2],fPintGlo[0],fPintGlo[1],fPintGlo[2]));
  
  return 0;
}

Int_t AliITSAlignMille::CalcDerivatives(Int_t paridx, Bool_t islpar) {
  /// calculate numerically (ROOT's style) the derivatives for
  /// local X intersection  and local Z intersection
  /// parlist: local  (islpar=kTRUE)  pgx0, pgz0, ugx0, ugz0  OR riemann's params
  ///          global (islpar=kFALSE) tx, ty, tz, psi, theta, phi (Raf's angles in deg)
  /// return 0 if success
  
  // copy initial parameters
  Double_t lpar[ITSMILLENLOCAL];
  Double_t gpar[ITSMILLENPARCH];
  for (Int_t i=0; i<ITSMILLENLOCAL; i++) lpar[i]=fLocalInitParam[i];
  for (Int_t i=0; i<ITSMILLENPARCH; i++) gpar[i]=fModuleInitParam[i];

  // trial with fixed dpar...
  Double_t dpar=0.0;

  if (islpar) { // track parameters
    //dpar=fLocalInitParam[paridx]*0.001;
    // set minimum dpar
    if (!fBOn) {
      if (paridx<2) dpar=1.0e-4; // translations
      else dpar=1.0e-6; // direction
    }
    else { // B Field
      // pepo: proviamo con 1/1000, poi evenctually 1/100...
      Double_t dfrac=0.01;
      switch(paridx) {
      case 0:
	// RMS cosmics: 1e-4
	dpar = TMath::Max(1.0e-6,TMath::Abs(fLocalInitParam[paridx]*dfrac)); 
	break;
      case 1: 
	// RMS cosmics: 0.2
	dpar = TMath::Max(0.002,TMath::Abs(fLocalInitParam[paridx]*dfrac)); 
	break;
      case 2: 
	// RMS cosmics: 9
	dpar = TMath::Max(0.09,TMath::Abs(fLocalInitParam[paridx]*dfrac)); 
	break;
      case 3: 
	// RMS cosmics: 7
	dpar = TMath::Max(0.07,TMath::Abs(fLocalInitParam[paridx]*dfrac)); 
	break;
      case 4: 
	// RMS cosmics: 0.3
	dpar = TMath::Max(0.003,TMath::Abs(fLocalInitParam[paridx]*dfrac)); 
	break;
      }
    }
  }
  else { // alignment global parameters
    //dpar=fModuleInitParam[paridx]*0.001;
    if (paridx<3) dpar=1.0e-4; // translations
    else dpar=1.0e-2; // angles    
  }

  AliDebug(3,Form("+++ using dpar=%g",dpar));
  
  // calculate derivative ROOT's like:
  //  using f(x+h),f(x-h),f(x+h/2),f(x-h2)...
  Double_t pintl1[3]; // f(x-h)
  Double_t pintl2[3]; // f(x-h/2)
  Double_t pintl3[3]; // f(x+h/2)
  Double_t pintl4[3]; // f(x+h)
    
  // first values
  if (islpar) lpar[paridx] -= dpar;
  else gpar[paridx] -= dpar;
  if (CalcIntersectionPoint(lpar, gpar)) return -2;
  for (Int_t i=0; i<3; i++) pintl1[i]=fPintLoc[i];

  // second values
  if (islpar) lpar[paridx] += dpar/2;
  else gpar[paridx] += dpar/2;
  if (CalcIntersectionPoint(lpar, gpar)) return -2;
  for (Int_t i=0; i<3; i++) pintl2[i]=fPintLoc[i];

  // third values
  if (islpar) lpar[paridx] += dpar;
  else gpar[paridx] += dpar;
  if (CalcIntersectionPoint(lpar, gpar)) return -2;
  for (Int_t i=0; i<3; i++) pintl3[i]=fPintLoc[i];

  // fourth values
  if (islpar) lpar[paridx] += dpar/2;
  else gpar[paridx] += dpar/2;
  if (CalcIntersectionPoint(lpar, gpar)) return -2;
  for (Int_t i=0; i<3; i++) pintl4[i]=fPintLoc[i];

  Double_t h2 = 1./(2.*dpar);
  Double_t d0 = pintl4[0]-pintl1[0];
  Double_t d2 = 2.*(pintl3[0]-pintl2[0]);
  fDerivativeXLoc = h2*(4*d2 - d0)/3.;
  if (TMath::Abs(fDerivativeXLoc) < 1.0e-9) fDerivativeXLoc=0.0;

  d0 = pintl4[2]-pintl1[2];
  d2 = 2.*(pintl3[2]-pintl2[2]);
  fDerivativeZLoc = h2*(4*d2 - d0)/3.;
  if (TMath::Abs(fDerivativeZLoc) < 1.0e-9) fDerivativeZLoc=0.0;

  AliDebug(3,Form("\n+++ derivatives +++ \n"));
  AliDebug(3,Form("+++ dXLoc/dpar = %g +++\n",fDerivativeXLoc));
  AliDebug(3,Form("+++ dZLoc/dpar = %g +++\n\n",fDerivativeZLoc));
  
  return 0;
}


Int_t AliITSAlignMille::AddLocalEquation(AliITSAlignMilleData &m) {
  /// Define local equation for current cluster in X and Z coor.
  /// and store them to memory
  /// return 0 if success
  
  // store first interaction point
  if (CalcIntersectionPoint(fLocalInitParam, fModuleInitParam)) return -4;  
  for (Int_t i=0; i<3; i++) fPintLoc0[i]=fPintLoc[i];
  AliDebug(2,Form("Intesect. point: L( %f , %f , %f )",fPintLoc[0],fPintLoc[1],fPintLoc[2]));
  
  // calculate local derivatives numerically
  Double_t dXdL[ITSMILLENLOCAL],dZdL[ITSMILLENLOCAL];
  for (Int_t i=0; i<fNLocal; i++) {
    if (CalcDerivatives(i,kTRUE)) return -1;
    dXdL[i]=fDerivativeXLoc;
    dZdL[i]=fDerivativeZLoc;
  }

  Double_t dXdG[ITSMILLENPARCH],dZdG[ITSMILLENPARCH];
  for (Int_t i=0; i<ITSMILLENPARCH; i++) {
    if (CalcDerivatives(i,kFALSE)) return -1;
    dXdG[i]=fDerivativeXLoc;
    dZdG[i]=fDerivativeZLoc;
  }

  AliDebug(2,Form("\n***************\n"));
  for (Int_t i=0; i<fNLocal; i++)
    AliDebug(2,Form("Local parameter %d - dXdpar = %g  - dZdpar = %g\n",i,dXdL[i],dZdL[i]));
  for (Int_t i=0; i<ITSMILLENPARCH; i++)
    AliDebug(2,Form("Global parameter %d - dXdpar = %g  - dZdpar = %g\n",i,dXdG[i],dZdG[i]));
  AliDebug(2,Form("\n***************\n"));

  // check if at least 1 local and 1 global derivs. are not null
  Double_t nonzero=0.0;
  for (Int_t i=0; i<fNLocal; i++) nonzero += TMath::Abs(dXdL[i]);
  if (nonzero==0.0) {
    AliInfo("Discarding local equations for this point beacuse of zero local X derivatives!");
    return -2;
  }
  nonzero=0.0;
  for (Int_t i=0; i<ITSMILLENPARCH; i++) nonzero += TMath::Abs(dXdG[i]);
  if (nonzero==0.0) {
    AliInfo("Discarding local equations for this point beacuse of zero global X derivatives!");
    return -2;
  }
  nonzero=0.0;
  for (Int_t i=0; i<fNLocal; i++) nonzero += TMath::Abs(dZdL[i]);
  if (nonzero==0.0) {
    AliInfo("Discarding local equations for this point beacuse of zero local Z derivatives!");
    return -2;
  }
  nonzero=0.0;
  for (Int_t i=0; i<ITSMILLENPARCH; i++) nonzero += TMath::Abs(dZdG[i]);
  if (nonzero==0.0) {
    AliInfo("Discarding local equations for this point beacuse of zero global Z derivatives!");
    return -2;
  }

  // ok, can copy to m

  AliDebug(2,Form("Adding local equation X with fMeas=%.6f  and fSigma=%.6f",(fMeasLoc[0]-fPintLoc0[0]), fSigmaLoc[0]));
  // set equation for Xloc coordinate
  for (Int_t i=0; i<fNLocal; i++) {
    m.GetIdxlocX()[i]=i;
    m.GetDerlocX()[i]=dXdL[i];
  }
  for (Int_t i=0; i<ITSMILLENPARCH; i++) {
    m.GetIdxgloX()[i]=fCurrentModuleInternalIndex*ITSMILLENPARCH+i;
    m.GetDergloX()[i]=dXdG[i];    
  }
  m.SetMeasX(fMeasLoc[0]-fPintLoc0[0]);
  m.SetSigmaX(fSigmaLoc[0]);
  
  AliDebug(2,Form("Adding local equation Z with fMeas=%.6f  and fSigma=%.6f",(fMeasLoc[2]-fPintLoc0[2]), fSigmaLoc[2]));
  // set equation for Zloc coordinate
  for (Int_t i=0; i<fNLocal; i++) {
    m.GetIdxlocZ()[i]=i;
    m.GetDerlocZ()[i]=dZdL[i];
  }
  for (Int_t i=0; i<ITSMILLENPARCH; i++) {
    m.GetIdxgloZ()[i]=fCurrentModuleInternalIndex*ITSMILLENPARCH+i;
    m.GetDergloZ()[i]=dZdG[i];    
  }
  m.SetMeasZ(fMeasLoc[2]-fPintLoc0[2]);
  m.SetSigmaZ(fSigmaLoc[2]);
 
  return 0;
}


void AliITSAlignMille::SetLocalEquations(const AliITSAlignMilleData *m, Int_t neq) {
  /// Set local equations with data stored in m
  /// return 0 if success
  
  for (Int_t j=0; j<neq; j++) {
    
    AliDebug(2,Form("setting local equation X with fMeas=%.6f  and fSigma=%.6f",m[j].GetMeasX(), m[j].GetSigmaX()));
    // set equation for Xloc coordinate
    for (Int_t i=0; i<fNLocal; i++) 
      SetLocalDerivative( m[j].GetIdxlocX()[i], m[j].GetDerlocX()[i] );
    
    for (Int_t i=0; i<ITSMILLENPARCH; i++)
      SetGlobalDerivative( m[j].GetIdxgloX()[i], m[j].GetDergloX()[i] );
    
    fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, m[j].GetMeasX(), m[j].GetSigmaX());  
    
    
    AliDebug(2,Form("setting local equation Z with fMeas=%.6f  and fSigma=%.6f",m[j].GetMeasZ(), m[j].GetSigmaZ()));
    // set equation for Zloc coordinate
    for (Int_t i=0; i<fNLocal; i++) 
      SetLocalDerivative( m[j].GetIdxlocZ()[i], m[j].GetDerlocZ()[i] );
    
    for (Int_t i=0; i<ITSMILLENPARCH; i++)
      SetGlobalDerivative( m[j].GetIdxgloZ()[i], m[j].GetDergloZ()[i] );
    
    fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, m[j].GetMeasZ(), m[j].GetSigmaZ());  
  }
}


void AliITSAlignMille::LocalFit(Int_t iTrack, Double_t *lTrackParam, Int_t lSingleFit) {
  /// Call local fit for this track
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  Int_t iRes = fMillepede->LocalFit(iTrack,lTrackParam,lSingleFit);
  AliDebug(2,Form("iRes = %d",iRes));
  //if (iRes && !lSingleFit) {
  if (!lSingleFit) { // Ruben Shahoyan's bug fix
    fMillepede->SetNLocalEquations(fMillepede->GetNLocalEquations()+1);
  }
}

void AliITSAlignMille::GlobalFit(Double_t *parameters,Double_t *errors,Double_t *pulls) {
  /// Call global fit; Global parameters are stored in parameters
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->GlobalFit(parameters,errors,pulls);
  AliInfo("Done fitting global parameters!");
}

Double_t AliITSAlignMille::GetParError(Int_t iPar) {
  /// Get error of parameter iPar
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return 0;
  }
  Double_t lErr = fMillepede->GetParError(iPar);
  return lErr;
}

void AliITSAlignMille::PrintGlobalParameters() {
  /// Print global parameters
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->PrintGlobalParameters();
}

// //_________________________________________________________________________
Int_t AliITSAlignMille::LoadSuperModuleFile(const Char_t *sfile)
{ 
  // load definitions of supermodules from a root file
  // return 0 if success

  TFile *smf=TFile::Open(sfile);
  if (!smf->IsOpen()) {
    AliInfo(Form("Cannot open supermodule file %s",sfile));
    return -1;
  }

  TClonesArray *sma=(TClonesArray*)smf->Get("ITSMilleSuperModules");
  if (!sma) {
    AliInfo(Form("Cannot find ITSMilleSuperModules array in file"));
    return -2;  
  }  
  Int_t nsma=sma->GetEntriesFast();
  AliInfo(Form("Array of SuperModules with %d entries\n",nsma));
  
  Char_t st[2048];
  char symname[250];
  UShort_t volid;
  TGeoHMatrix m;

  for (Int_t i=0; i<nsma; i++) {
    AliAlignObjParams *a = (AliAlignObjParams*)sma->UncheckedAt(i);
    volid=a->GetVolUID();
    strncpy(st,a->GetSymName(),TMath::Min(sizeof(st),strlen(a->GetSymName())+1));
    a->GetMatrix(m);
    memset(symname,0,250*sizeof(char));
    sscanf(st,"%249s",symname);
    // decode module list
    char *stp=strstr(st,"ModuleList:");
    if (!stp) return -3;
    stp += 11;
    int idx[2200];
    char spp[200]; int jp=0;
    char cl[20];
    strncpy(st,stp,TMath::Min(sizeof(st),strlen(stp)+1));
    int l=strlen(st);
    int j=0;
    int n=0;

    while (j<=l) {
      if (st[j]==9 || st[j]==32 || st[j]==10 || st[j]==0) {
	spp[jp]=0;
	jp=0;
	if (strlen(spp)) {
	  int k=strcspn(spp,"-");
	  if (k<int(strlen(spp))) { // c'e' il -
	    strncpy(cl,&(spp[k+1]), TMath::Min(sizeof(cl),strlen(&spp[k+1])+1));
	    spp[k]=0;
	    int ifrom=atoi(spp); int ito=atoi(cl);
	    for (int b=ifrom; b<=ito; b++) {
	      idx[n]=b;
	      n++;
	    }
	  }
	  else { // numerillo singolo
	    idx[n]=atoi(spp);
	    n++;
	  }
	}
      }
      else {
	spp[jp]=st[j];
	jp++;
      }
      j++;
    }
    UShort_t volidsv[2198];
    for (j=0;j<n;j++) {
      volidsv[j]=AliITSAlignMilleModule::GetVolumeIDFromIndex(idx[j]);
      if (!volidsv[j]) {
	AliInfo(Form("Index %d not valid (range 0->2197)",idx[j]));
	return -5;
      }
    }
    Int_t smindex=int(2198+volid-14336); // virtual index
    fSuperModule[fNSuperModules]=new AliITSAlignMilleModule(smindex,volid,symname,&m,n,volidsv);

    //-------------
    fNSuperModules++;
  }

  smf->Close();

  fUseSuperModules=1;
  return 0;
}

