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
/// Alignment class fro the ALICE ITS detector
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
#include <TGeoMatrix.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TVirtualFitter.h>

#include "AliITSAlignMille2.h"
#include "AliITSgeomTGeo.h"
#include "AliGeomManager.h"
#include "AliMillePede2.h"
#include "AliTrackPointArray.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"
#include "TSystem.h"  // come si fa?
#include "AliTrackFitterRieman.h"

/// \cond CLASSIMP
ClassImp(AliITSAlignMille2)
/// \endcond

AliITSAlignMille2* AliITSAlignMille2::fgInstance = 0;  
Int_t              AliITSAlignMille2::fgInstanceID = 0;

AliITSAlignMille2::AliITSAlignMille2(const Char_t *configFilename, Bool_t initmille) 
: TObject(),
  fMillepede(0),
  fStartFac(16.), 
  fResCutInitial(100.), 
  fResCut(100.),
  fNGlobal(0),
  fNLocal(4),
  fNStdDev(3),
  fIsMilleInit(kFALSE),
  fSensorsIn(kFALSE),
  fParSigTranslations(0.0100),
  fParSigRotations(0.1),
//
  fCurrentModule(0),
  fTrack(0),
  fCluster(),
  fGlobalDerivatives(0), 
//
  fMinNPtsPerTrack(3),
  fInitTrackParamsMeth(1),
  fTotBadLocEqPoints(0),
  fRieman(0),
  //
  fUseGlobalDelta(kFALSE),
  fRequirePoints(kFALSE),
  fTempExcludedModule(-1),
  //
  fGeometryFileName("geometry.root"),
  fPreAlignmentFileName(""),
  fGeoManager(0),
  fIsConfigured(kFALSE),
  fPreAlignQF(0),
//
  fPrealignment(0),
  fMilleModule(2),
  fSuperModule(2),
  fNModules(0),
  fNSuperModules(0),
  fUsePreAlignment(kFALSE),
  fUseSortedTracks(kTRUE),
  fBOn(kFALSE),
  fBField(0.0),
  fBug(0)
{
  /// main constructor that takes input from configuration file
  //  
  fMillepede = new AliMillePede2();
  for (int i=3;i--;) fSigmaFactor[i] = 1.0;
  //
  // new RS
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
  //
  Int_t lc=LoadConfig(configFilename);
  if (lc) {
    AliInfo(Form("Error %d loading configuration from %s",lc,configFilename));
  }
  else {    
    fIsConfigured=kTRUE;
    if (initmille) {
      AliInfo(Form("Initializing Millepede with %d gpar, %d lpar and %d stddev ...",fNGlobal, fNLocal, fNStdDev));
      Init(fNGlobal, fNLocal, fNStdDev);      
      ResetLocalEquation();    
      AliInfo("Parameters initialized to zero");
    }
    else {
      AliInfo("Millepede has not been initialized ... ");
    }
  }
  //
  fgInstance = this;
  fgInstanceID++;
  //
}

AliITSAlignMille2::~AliITSAlignMille2() {
  /// Destructor
  if (fMillepede)         delete fMillepede;
  if (fGlobalDerivatives) delete[] fGlobalDerivatives;
  if (fRieman)            delete fRieman;
  if (fPrealignment)      delete fPrealignment;
  fMilleModule.Delete();
  fSuperModule.Delete();
  if (--fgInstanceID==0) fgInstance = 0;
}

///////////////////////////////////////////////////////////////////////

Int_t AliITSAlignMille2::LoadConfig(const Char_t *cfile) {
  /// return 0 if success
  ///        1 if error in module index or voluid
  
  FILE *pfc=fopen(cfile,"r");
  if (!pfc) return -1;
  
  Char_t st[200],st2[200];
  Char_t tmp[100];
  Int_t idx,itx,ity,itz,ith,ips,iph;
  Float_t f1,f2,f3;
  UShort_t voluid;
  Int_t nmod=0;
  //
  while (fgets(st,200,pfc)) {
    //
    for (int i=0; i<int(strlen(st)); i++) if (st[i]=='#') st[i]=0; // skip comments
    //
    if (strstr(st,"GEOMETRY_FILE")) {
      sscanf(st,"%s %s",tmp,st2);
      if (gSystem->AccessPathName(st2)) { AliInfo("*** WARNING! *** geometry file not found! "); return -1;}  
      fGeometryFileName=st2;
      InitGeometry();
    }
    //
    if (strstr(st,"PREALIGNMENT_FILE")) {
      sscanf(st,"%s %s",tmp,st2);
      if (gSystem->AccessPathName(st2)) { AliInfo("*** WARNING! *** prealignment file not found! "); return -1;}  
      fPreAlignmentFileName=st2;
      itx=ApplyToGeometry();
      if (itx) { AliInfo(Form("*** WARNING! *** error %d reading prealignment file! ",itx)); return -6;}
    }
    //
    if (strstr(st,"SUPERMODULE_FILE")) {
      sscanf(st,"%s %s",tmp,st2);
      if (gSystem->AccessPathName(st2)) { AliInfo("*** WARNING! *** supermodule file not found! "); return -1;}  
      if (LoadSuperModuleFile(st2)) return -1;
    }
    //
    if (strstr(st,"SET_B_FIELD")) {
      sscanf(st,"%s %f",tmp,&f1);
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
    //
    if (strstr(st,"SET_MINPNT_TRA")) {
      sscanf(st,"%s %d",tmp,&idx);
      fMinNPtsPerTrack=idx;
    }
    //
    if (strstr(st,"SET_PARSIG_TRA")) {
      sscanf(st,"%s %f",tmp,&f1);
      fParSigTranslations=f1;
    }
    //
    if (strstr(st,"SET_PARSIG_ROT")) {
      sscanf(st,"%s %f",tmp,&f1);
      fParSigRotations=f1;
    }
    //
    if (strstr(st,"SET_NSTDDEV")) {
      sscanf(st,"%s %d",tmp,&idx);
      fNStdDev=idx;
    }
    //
    if (strstr(st,"SET_RESCUT_INIT")) {
      sscanf(st,"%s %f",tmp,&f1);
      fResCutInitial=f1;
    }
    //
    if (strstr(st,"SET_RESCUT_OTHER")) {
      sscanf(st,"%s %f",tmp,&f1);
      fResCut=f1;
    }
    //
    if (strstr(st,"SET_LOCALSIGMAFACTOR")) {
      f1=f2=f3=0;
      sscanf(st,"%s %f %f %f",tmp,&f1,&f2,&f3);
      if (f1>0) fSigmaFactor[0] = f1;
      if (f2>0) fSigmaFactor[1] = f2; else fSigmaFactor[1]=fSigmaFactor[0];
      if (f3>0) fSigmaFactor[2] = f3; else fSigmaFactor[2]=fSigmaFactor[1];
    }
    //
    if (strstr(st,"SET_STARTFAC")) {
      sscanf(st,"%s %f",tmp,&f1);
      fStartFac=f1;
    }
    //
    // >> RS
    if (strstr(st,"REQUIRE_POINT")) {
      // syntax:   REQUIRE_POINT where ndet updw nreqpts
      //    where = LAYER or DETECTOR
      //    ndet = detector number: 1-6 for LAYER and 1-3 for DETECTOR (SPD=1, SDD=2, SSD=3)
      //    updw = 1 for Y>0, -1 for Y<0, 0 if not specified
      //    nreqpts = minimum number of points of that type
      sscanf(st,"%s %s %d %d %d",tmp,st2,&itx,&ity,&itz);
      itx--;
      if (strstr(st2,"LAYER")) {
	if (itx<0 || itx>5) return -7;
	if (ity>0) fNReqLayUp[itx]=itz;
	else if (ity<0) fNReqLayDown[itx]=itz;
	else fNReqLay[itx]=itz;
	fRequirePoints=kTRUE;
      }
      else if (strstr(st2,"DETECTOR")) { // DETECTOR
	if (itx<0 || itx>2) return -7;
	if (ity>0) fNReqDetUp[itx]=itz;
	else if (ity<0) fNReqDetDown[itx]=itz;
	else fNReqDet[itx]=itz;	
	fRequirePoints=kTRUE;
      }
    }
    // << RS
    
    //
    if (strstr(st,"MODULE_INDEX") || strstr(st,"MODULE_VOLUID")) { 
      f1=f2=f3=0;
      sscanf(st,"%s %d %d %d %d %d %d %d %f %f %f",tmp,&idx,&itx,&ity,&itz,&iph,&ith,&ips,&f1,&f2,&f3);
      //
      if (idx<=kMaxITSSensID) voluid=GetModuleVolumeID(idx);
      else voluid = UShort_t(idx);
      //
      if (voluid>=kMinITSSupeModuleID) { // custom supermodule
	int ism=-1;
	for (int j=0; j<fNSuperModules; j++) if (voluid==GetSuperModule(j)->GetVolumeID()) ism=j;
	if (ism<0) return -1; // bad volid
	fMilleModule.AddAtAndExpand(new AliITSAlignMille2Module(*GetSuperModule(ism)),nmod);
	// >> RS
// 	if (f1>0) {
// 	  for (int kk=0; kk<GetMilleModule(nmod)->GetNSensitiveVolumes(); kk++) {
// 	    idx=AliITSAlignMille2Module::GetIndexFromVolumeID(GetMilleModule(nmod)->GetSensitiveVolumeVolumeID()[kk]);
// 	    if (idx>=0) fSensVolSigmaXfactor[idx]=f1;
// 	  }
// 	}
// 	if (f2>0) {
// 	  for (int kk=0; kk<GetMilleModule(nmod)->GetNSensitiveVolumes(); kk++) {
// 	    idx=AliITSAlignMille2Module::GetIndexFromVolumeID(GetMilleModule(nmod)->GetSensitiveVolumeVolumeID()[kk]);
// 	    if (idx>=0) fSensVolSigmaZfactor[idx]=f2;
// 	  }
// 	}
	// << RS
      }
      else if (idx<=kMaxITSSensVID) {
	fMilleModule.AddAtAndExpand(new AliITSAlignMille2Module(voluid),nmod);
	AliITSAlignMille2Module* md = (AliITSAlignMille2Module*) fMilleModule[nmod];
	fSensorsIn = kTRUE;
	md->SetSensorsProvided();
      }
      else return -1;  // bad volid
      //
      AliITSAlignMille2Module* mod = GetMilleModule(nmod);
      mod->SetFreeDOF(kDOFTX,itx!=0);
      mod->SetFreeDOF(kDOFTY,ity!=0);
      mod->SetFreeDOF(kDOFTZ,itz!=0);
      mod->SetFreeDOF(kDOFPH,iph!=0);
      mod->SetFreeDOF(kDOFTH,ith!=0);
      mod->SetFreeDOF(kDOFPS,ips!=0);
      //
      mod->SetUniqueID(nmod);
      if (f1>0) mod->SetSigmaXFactor(f1);
      if (f2>0) mod->SetSigmaYFactor(f2); else mod->SetSigmaYFactor(mod->GetSigmaXFactor());
      if (f3>0) mod->SetSigmaZFactor(f3); else mod->SetSigmaZFactor(mod->GetSigmaYFactor());
      nmod++;
    }
    //
  } // end while
  //
  fNModules = nmod;
  fNGlobal = fNModules*kNParCh;
  //
  fclose(pfc);
  //
  // set parent/child relationship for modules to align
  printf("Setting parent/child relationships\n");
  //
  for (int ipar=0;ipar<nmod;ipar++) {
    AliITSAlignMille2Module* parent = GetMilleModule(ipar);
    if (parent->GetIndex()<=kMaxITSSensID) continue; // sensor cannot be a parent
    //
    for (int icld=0;icld<nmod;icld++) {
      if (icld==ipar) continue;
      AliITSAlignMille2Module* child = GetMilleModule(icld);
      if (!child->BelongsTo(parent)) continue;
      //
      AliITSAlignMille2Module* parOld = child->GetParent();
      if (parOld && parOld->GetNSensitiveVolumes()<parent->GetNSensitiveVolumes()) continue; // parOld is closer
      child->SetParent(parent);
    }
    //
  }
  //
  // reorder the modules in such a way that parents come first
  for (int icld=0;icld<nmod;icld++) {
    AliITSAlignMille2Module* child  = GetMilleModule(icld);
    AliITSAlignMille2Module* parent; 
    while ( (parent=child->GetParent()) &&  (parent->GetUniqueID()<child->GetUniqueID()) ) {
      // swap
      fMilleModule[icld] = parent;
      fMilleModule[parent->GetUniqueID()] = child;
      child->SetUniqueID(parent->GetUniqueID());
      parent->SetUniqueID(icld);
      child = parent;
    }
    //
  } 
  //
  // go over the child->parent chain and mark modules with explicitly provided sensors
  for (int icld=nmod;icld--;) {
    AliITSAlignMille2Module* child = GetMilleModule(icld);
    AliITSAlignMille2Module* parent = child->GetParent();
    if (!parent) continue;
    parent->SetSensorsProvided( child->AreSensorsProvided() );
    if (!parent->AreSensorsProvided()) continue;
    // suppress unused sensors
    for (int isn=0;isn<parent->GetNSensitiveVolumes();isn++) {
      int snVID = parent->GetSensVolVolumeID(isn);
      // check if this sensor is explicitly requested
      for (int imd=nmod;imd--;) if (GetMilleModule(imd)->GetVolumeID() == snVID) {snVID = -1; break;}
      //
      if (snVID==-1) continue; // found this sensor, do nothing
      //
      // otherwise, remove this sensor from the module list
      AliInfo(Form("Removing sensor %d from %s",snVID,parent->GetName()));
      parent->DelSensitiveVolume(isn--);
    }
  }
  //
  fGlobalDerivatives = new Double_t[fNGlobal];
  memset(fGlobalDerivatives,0,fNGlobal*sizeof(Double_t));
  //
  return 0;
}


void AliITSAlignMille2::SetRequiredPoint(Char_t* where, Int_t ndet, Int_t updw, Int_t nreqpts) 
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


Int_t AliITSAlignMille2::GetModuleIndex(const Char_t *symname) {
  /// index from symname
  if (!symname) return -1;
  for (Int_t i=0;i<=kMaxITSSensID; i++) {
    if (!strcmp(symname,AliITSgeomTGeo::GetSymName(i))) return i;
  }
  return -1;
}

Int_t AliITSAlignMille2::GetModuleIndex(UShort_t voluid) {
  /// index from volume ID
  AliGeomManager::ELayerID lay = AliGeomManager::VolUIDToLayer(voluid);
  if (lay<1|| lay>6) return -1;
  Int_t idx=Int_t(voluid)-2048*lay;
  if (idx>=AliGeomManager::LayerSize(lay)) return -1;
  for (Int_t ilay=1; ilay<lay; ilay++) 
    idx += AliGeomManager::LayerSize(ilay);
  return idx;
}

UShort_t AliITSAlignMille2::GetModuleVolumeID(const Char_t *symname) {
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

UShort_t AliITSAlignMille2::GetModuleVolumeID(Int_t index) {
  /// volume ID from index
  if (index<0) return 0;
  if (index<2198)
    return GetModuleVolumeID(AliITSgeomTGeo::GetSymName(index));
  else {
    for (int i=0; i<fNSuperModules; i++) {
      if (GetSuperModule(i)->GetIndex()==index) return GetSuperModule(i)->GetVolumeID();
    }
  }
  return 0;
}

void AliITSAlignMille2::InitGeometry() {
  /// initialize geometry
  AliGeomManager::LoadGeometry(fGeometryFileName.Data());
  fGeoManager = AliGeomManager::GetGeometry();
  if (!fGeoManager) {
    AliInfo("Couldn't initialize geometry");
    return;
  }
}

void AliITSAlignMille2::Init(Int_t nGlobal,  /* number of global paramers */
			   Int_t nLocal,   /* number of local parameters */
			   Int_t nStdDev   /* std dev cut */ )
{
  /// Initialization of AliMillepede. Fix parameters, define constraints ...
  fMillepede->InitMille(nGlobal,nLocal,nStdDev,fResCut,fResCutInitial);
  fIsMilleInit = kTRUE;
  
  /// Fix non free parameters
  for (Int_t i=0; i<fNModules; i++) {
    for (Int_t j=0; j<kNParCh; j++) {
      if (!GetMilleModule(i)->IsFreeDOF(j)) FixParameter(i*kNParCh+j,0.0);
      else {
	// pepopepo: da verificare il settaggio delle sigma, ma forse va bene...
	Double_t parsig=0;
	if (j<3) parsig = fParSigTranslations; // translations (0.0100 cm)
	else     parsig = fParSigRotations; // rotations (1/10 deg)
	FixParameter(i*kNParCh+j,parsig);
      }
    }    
  }
  //
  // Set iterations
  if (fStartFac>1) fMillepede->SetIterations(fStartFac);          
}


void AliITSAlignMille2::AddConstraint(Double_t *par, Double_t value) {
  /// Constrain equation defined by par to value
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->SetGlobalConstraint(par, value);
  AliInfo("Adding constraint");
}

void AliITSAlignMille2::InitGlobalParameters(Double_t *par) {
  /// Initialize global parameters with par array
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->SetGlobalParameters(par);
  AliInfo("Init Global Parameters");
}
 
void AliITSAlignMille2::FixParameter(Int_t iPar, Double_t value) {
  /// Parameter iPar is encourage to vary in [-value;value]. 
  /// If value == 0, parameter is fixed
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->SetParSigma(iPar, value);
  if (value==0) AliInfo(Form("Parameter %i Fixed", iPar));
}

void AliITSAlignMille2::ResetLocalEquation()
{
  /// Reset the derivative vectors
  for(int i=fNLocal;i--;)  fLocalDerivatives[i] = 0.0;
  for(int i=fNGlobal;i--;) fGlobalDerivatives[i] = 0.0;
}

Int_t AliITSAlignMille2::ApplyToGeometry() 
{
  // apply starting realignment to ideal geometry
  //
  if (!fGeoManager) return -1; 
  TFile *pref = new TFile(fPreAlignmentFileName.Data());
  if (!pref->IsOpen()) return -2;
  fPrealignment = (TClonesArray*)pref->Get("ITSAlignObjs");
  if (!fPrealignment) return -3;  
  Int_t nprea = fPrealignment->GetEntriesFast();
  AliInfo(Form("Array of input misalignments with %d entries",nprea));
  //
  for (int ix=0; ix<nprea; ix++) {
    AliAlignObjParams *preo=(AliAlignObjParams*) fPrealignment->At(ix);
    Int_t index=AliITSAlignMille2Module::GetIndexFromVolumeID(preo->GetVolUID());
    if (index>=0) {
      if (index>=fPreAlignQF.GetSize()) fPreAlignQF.Set(index+10);
      fPreAlignQF[index] = (int) preo->GetUniqueID()+1;
    }
    //TString nms = preo->GetSymName();
    //if (!nms.Contains("Ladder")) continue; //RRR
    if (!preo->ApplyToGeometry()) return -4;
  }
  //
  pref->Close();
  delete pref;
  //
  fUsePreAlignment = kTRUE;
  return 0;
}

Int_t AliITSAlignMille2::GetPreAlignmentQualityFactor(Int_t index) const
{
  if (!fUsePreAlignment || index<0 || index>=fPreAlignQF.GetSize()) return -1;
  return fPreAlignQF[index]-1;
}

AliTrackPointArray *AliITSAlignMille2::PrepareTrack(const AliTrackPointArray *atp) {
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
  // >> RS
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
  
  // << RS
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
    TGeoHMatrix *svOrigMatrix = GetMilleModule(intidx[idx[i]])->GetSensitiveVolumeOrigGlobalMatrix(p.GetVolumeID());
    // get back real local coordinates: use OriginalGlobalMatrix because AliTrackPoints were written
    // with idel geometry  
    Double_t pg[3],pl[3];
    pg[0]=p.GetX();
    pg[1]=p.GetY();
    pg[2]=p.GetZ();
    //    printf("Global coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]);
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

    // >> RS
    if (fBug==1) {
      // correzione bug LAYER 5  SSD temporanea..
      int ssdidx=AliITSAlignMille2Module::GetIndexFromVolumeID(p.GetVolumeID());
      if (ssdidx>=500 && ssdidx<1248) {
	int ladder=(ssdidx-500)%22;
	if (ladder==18) p.SetVolumeID(AliITSAlignMille2Module::GetVolumeIDFromIndex(ssdidx+1));
	if (ladder==19) p.SetVolumeID(AliITSAlignMille2Module::GetVolumeIDFromIndex(ssdidx-1));
      }
    }
    
    // << RS

    /// get (evenctually prealigned) matrix of sens. vol.
    TGeoHMatrix *svMatrix = GetMilleModule(intidx[idx[i]])->GetSensitiveVolumeMatrix(p.GetVolumeID());
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
    //    printf("New Gl coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]);
    AliDebug(3,Form("New global coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]));
    atps->AddPoint(npto,&p);
    AliDebug(2,Form("Adding point[%d] = ( %f , %f , %f )     volid = %d",npto,atps->GetX()[npto],atps->GetY()[npto],atps->GetZ()[npto],atps->GetVolumeID()[npto] ));

    npto++;
  }

  return atps;
}



AliTrackPointArray *AliITSAlignMille2::SortTrack(const AliTrackPointArray *atp) {
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


Int_t AliITSAlignMille2::InitModuleParams() {
  /// initialize geometry parameters for a given detector
  /// for current cluster (fCluster)
  /// fGlobalInitParam[] is set as:
  ///    [tx,ty,tz,psi,theta,phi]
  ///    (old was [tx,ty,tz,theta,psi,phi] ROOT's angles...)
  /// *** At the moment: using Raffalele's angles definition ***
  ///
  /// return 0 if success

  if (!fGeoManager) {
    AliInfo("ITS geometry not initialized!");
    return -1;
  }

  // now 'voluid' is the volumeID of a SENSITIVE VOLUME (coming from a cluster)

  // set the internal index (index in module list)
  UShort_t voluid=fCluster.GetVolumeID();
  //
  // IT IS VERY IMPORTANT: start from the end of the list, where the childs are located !!!
  Int_t k=fNModules-1;
  fCurrentModule = 0;
  while (k>=0 && ! (fCurrentModule=GetMilleModule(k))->IsIn(voluid)) {
     // VERY IMPORTANT: if the sensors were explicitly provided, don't look in the supermodules  
    if (fSensorsIn && fCurrentModule->GetVolumeID() > kMaxITSSensVID) {k=-1; break;} 
    k--; 
  }
  if (k<0) return -3;
  //  fCurrentModule = GetMilleModule(k);
  //
  fModuleInitParam[0] = 0.0;
  fModuleInitParam[1] = 0.0;
  fModuleInitParam[2] = 0.0;
  fModuleInitParam[3] = 0.0; // psi   (X)
  fModuleInitParam[4] = 0.0; // theta (Y)
  fModuleInitParam[5] = 0.0; // phi   (Z)
  
  /// get (evenctually prealigned) matrix of sens. vol.
  TGeoHMatrix *svMatrix = fCurrentModule->GetSensitiveVolumeMatrix(voluid);
  
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
  fSigmaLoc[1] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[4])); // RS
  fSigmaLoc[2] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[8]));

  // set minimum value for SigmaLoc to 10 micron 
  if (fSigmaLoc[0]<0.0010) fSigmaLoc[0]=0.0010;
  if (fSigmaLoc[2]<0.0010) fSigmaLoc[2]=0.0010;

  // multiply local sigmas by global and module specific factor 
  for (int i=3;i--;) fSigmaLoc[i] *= fSigmaFactor[i]*fCurrentModule->GetSigmaFactor(i);
  //
  AliDebug(2,Form("Setting StDev from CovMat : fSigmaLocX=%g  fSigmaLocY=%g fSigmaLocZ=%g \n",fSigmaLoc[0] ,fSigmaLoc[1] ,fSigmaLoc[2] ));
   
  return 0;
}

void AliITSAlignMille2::Print(Option_t*) const 
{
  ///
  printf("*** AliMillepede for ITS ***\n");
  printf("    Number of defined super modules: %d\n",fNModules);
  printf("    Obtained parameters refer to %s Deltas\n",fUseGlobalDelta ? "GLOBAL":"LOCAL");
  //
  if (fGeoManager)
    printf("    geometry loaded from %s\n",fGeometryFileName.Data());
  else
    printf("    geometry not loaded\n");
  //  
  if (fUsePreAlignment) 
    printf("    using prealignment from %s \n",fPreAlignmentFileName.Data());
  else
    printf("    prealignment not used\n");    
  //
  //
  if (fBOn) 
    printf("    B Field set to %f T - using Riemann's helices\n",fBField);
  else
    printf("    B Field OFF - using straight lines \n");
  //
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
  //  
  printf("\n    Millepede configuration parameters:\n");
  printf("        init parsig for translations  : %.4f\n",fParSigTranslations);
  printf("        init parsig for rotations     : %.4f\n",fParSigRotations);
  printf("        init value for chi2 cut       : %.4f\n",fStartFac);
  printf("        first iteration cut value     : %.4f\n",fResCutInitial);
  printf("        other iterations cut value    : %.4f\n",fResCut);
  printf("        number of stddev for chi2 cut : %d\n",fNStdDev);
  printf("        mult. fact. for local sigmas  : %.4f %.4f %.4f\n",fSigmaFactor[0],fSigmaFactor[1],fSigmaFactor[2]);

  printf("List of defined modules:\n");
  printf("  intidx\tindex\tvoluid\tname\n");
  for (int i=0; i<fNModules; i++) {
    AliITSAlignMille2Module* md = GetMilleModule(i); 
    printf("  %d\t%d\t%d\t%s\n",i,md->GetIndex(),md->GetVolumeID(),md->GetName());
  }
}

AliITSAlignMille2Module  *AliITSAlignMille2::GetMilleModuleByVID(UShort_t voluid) const
{
  // return pointer to a define supermodule
  // return NULL if error
  Int_t i=IsDefined(voluid);
  if (i<0) return NULL;
  return GetMilleModule(i);
}

Bool_t AliITSAlignMille2::InitRiemanFit() 
{
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

Bool_t fullErr2D = kTRUE;

void trackFit2D(Int_t &, Double_t *, double &chi2, double *par, int)
{
  const double kTiny = 1.e-14;
  chi2 = 0;
  static AliTrackPoint pnt;
  //
  enum {kAX,kAZ,kBX,kBZ};
  enum {kXX=0,kXY=1,kXZ=2,kYX=kXY,kYY=3,kYZ=4,kZX=kXZ,kZY=kYZ,kZZ=5};
  //
  AliITSAlignMille2* alig = AliITSAlignMille2::GetInstance();
  AliTrackPointArray* track = alig->GetCurrentTrack();
  //
  int npts = track->GetNPoints();
  for (int ip=0;ip<npts;ip++) {
    track->GetPoint(pnt,ip);
    const float *cov = pnt.GetCov();
    double y  = pnt.GetY();
    double dx = pnt.GetX() - (par[kAX]+y*par[kBX]);
    double dz = pnt.GetZ() - (par[kAZ]+y*par[kBZ]);
    double xxe = cov[kXX];
    double zze = cov[kZZ];
    double xze = cov[kXZ];
    //
    if (fullErr2D) {
      xxe += par[kBX]*par[kBX]*cov[kYY]-2.*par[kBX]*cov[kXY];
      zze += par[kBZ]*par[kBZ]*cov[kYY]-2.*par[kBZ]*cov[kZY];
      xze += par[kBX]*par[kBZ]*cov[kYY]-cov[kYZ]*par[kBZ]-cov[kXY]*par[kBX];
    }
    //
    double det = xxe*zze - xze*xze;
    if (det<kTiny) {
      printf("Negative diag. error (det=%+e) |sxx:%+e szz:%+e sxz:%+e| bx:%+e bz:%+e|\n"
	     "Discarding correlation term\n",det,xxe,zze,xze,par[kBX],par[kBZ]);
      xxe = cov[kXX];
      zze = cov[kZZ];
      xze = cov[kXZ];
      fullErr2D = kFALSE;
    }
    double xxeI = zze/det;
    double zzeI = xxe/det;
    double xzeI =-xze/det;
    //
    chi2 += dx*dx*xxeI + dz*dz*zzeI + 2.*dx*dz*xzeI;
    // 
    //    printf("%d | %+e %+e %+e %+e %+e -> %+e\n",ip,dx,dz,xxeI,zzeI,xzeI,  chi2);
  }
  //
}

void AliITSAlignMille2::InitTrackParams(int meth) 
{
  /// initialize local parameters with different methods
  /// for current track (fTrack)
  Int_t npts=0;
  AliTrackPoint ap;
  double sX=0,sXY=0,sZ=0,sZY=0,sY=0,sYY=0,det=0;
  // simple linear interpolation
  // get local starting parameters (to be substituted by ESD track parms)
  // local parms (fLocalInitParam[]) are:
  //      [0] = global x coord. of straight line intersection at y=0 plane
  //      [1] = global z coord. of straight line intersection at y=0 plane
  //      [2] = px/py  
  //      [3] = pz/py
  // test #1: linear fit in x(y) and z(y)
  npts = fTrack->GetNPoints();
  AliDebug(3,Form("*** initializing track with %d points ***",npts));
  for (int i=npts;i--;) {
    sY  += fTrack->GetY()[i];
    sYY += fTrack->GetY()[i]*fTrack->GetY()[i];
    sX  += fTrack->GetX()[i];
    sXY += fTrack->GetX()[i]*fTrack->GetY()[i];
    sZ  += fTrack->GetZ()[i];
    sZY += fTrack->GetZ()[i]*fTrack->GetY()[i];
  }
  det = sYY*npts-sY*sY;
  if (det==0) det = 1E-20;
  fLocalInitParam[0] = (sX*sYY-sY*sXY)/det;
  fLocalInitParam[2] = (sXY*npts-sY*sX)/det;
  //
  fLocalInitParam[1] = (sZ*sYY-sY*sZY)/det;
  fLocalInitParam[3] = (sZY*npts-sY*sZ)/det;
  //AliDebug(2,Form("X = p0gx + ugx*Y : p0gx = %f +- %f    ugx = %f +- %f\n",fLocalInitParam[0],f1->GetParError(0),fLocalInitParam[2],f1->GetParError(1)));
  AliDebug(2,Form("X = p0gx + ugx*Y : p0gx = %f    ugx = %f\n",fLocalInitParam[0],fLocalInitParam[2]));
  //
  if (meth==1) return;
  //
  // perform full fit accounting for cov.matrix
  static TVirtualFitter *minuit = 0;
  static Double_t step[5]   = {1E-3,1E-3,1E-4,1E-4,1E-5};
  static Double_t arglist[10];
  //
  if (!minuit) {
    minuit = TVirtualFitter::Fitter(0,4);
    minuit->SetFCN(trackFit2D);
    arglist[0] = 1;
    minuit->ExecuteCommand("SET ERR",arglist, 1);
    //
    arglist[0] = -1;
    minuit->ExecuteCommand("SET PRINT",arglist,1);
    //
  }
  //
  minuit->SetParameter(0, "ax",   fLocalInitParam[0], step[0], 0,0);
  minuit->SetParameter(1, "az",   fLocalInitParam[1], step[1], 0,0);
  minuit->SetParameter(2, "bx",   fLocalInitParam[2], step[2], 0,0);
  minuit->SetParameter(3, "bz",   fLocalInitParam[3], step[3], 0,0);
  //
  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.001; // tolerance 
  fullErr2D = kTRUE;
  minuit->ExecuteCommand("MIGRAD",arglist,2);
  //
  for (int i=0;i<4;i++) fLocalInitParam[i] = minuit->GetParameter(i);
  for (int i=0;i<4;i++) for (int j=0;j<4;j++) fLocalInitParEr[i][j] = minuit->GetCovarianceMatrixElement(i,j);
  //
}


Int_t AliITSAlignMille2::IsDefined(UShort_t voluid) const
{
  // checks if supermodule 'voluid' is defined and return the internal index
  // return -1 if error
  Int_t k = fNModules-1;
  while ( k>=0 && !(voluid==GetMilleModule(k)->GetVolumeID()) ) k--;  
  if (k<0) return -1; 
  return k;
}

Int_t AliITSAlignMille2::IsContained(UShort_t voluid) const
{
  // checks if the sensitive module 'voluid' is contained inside a supermodule and return the internal index of the last identified supermodule
  // return -1 if error
  if (AliITSAlignMille2Module::GetIndexFromVolumeID(voluid)<0) return -1;
  Int_t k=fNModules-1;
  while (k>=0 && !(GetMilleModule(k)->IsIn(voluid)) ) k--;  
  if (k<0) return -1; 
  return k;
}

Bool_t AliITSAlignMille2::CheckVolumeID(UShort_t voluid) const 
{
  /// check if a sensitive volume is contained inside one of the defined supermodules
  Int_t k=fNModules-1;
  while (k>=0 && !(GetMilleModule(k)->IsIn(voluid)) ) k--;  
  if (k>=0) return kTRUE;
  return kFALSE;
}

Int_t AliITSAlignMille2::CheckCurrentTrack() {
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

Int_t AliITSAlignMille2::ProcessTrack(const AliTrackPointArray *track) {
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
  static Mille2Data md[100];
  
  for (Int_t ipt=0; ipt<npts; ipt++) {
    fTrack->GetPoint(fCluster,ipt);
    AliDebug(2,Form("\n--- processing point %d --- \n",ipt));    

    // set geometry parameters for the the current module
    if (InitModuleParams()) continue;
    AliDebug(2,Form("    VolID=%d  Index=%d  InternalIdx=%d  symname=%s\n", 
		    track->GetVolumeID()[ipt], fCurrentModule->GetIndex(),
		    fCurrentModule->GetUniqueID(), AliGeomManager::SymName(track->GetVolumeID()[ipt]) ));
    AliDebug(2,Form("    Preprocessed Point = ( %f , %f , %f ) \n",fCluster.GetX(),fCluster.GetY(),fCluster.GetZ()));
    
    if (!AddLocalEquation(md[nloceq])) nloceq++;    
    else fTotBadLocEqPoints++;
  } // end loop over points
  //
  delete fTrack;
  fTrack=NULL;
  // not enough good points!
  if (nloceq<fMinNPtsPerTrack) return -1;
  //
  // finally send local equations to millepede
  SetLocalEquations(md,nloceq);
  fMillepede->SaveRecordData(); // RRR
  //
  return 0;
}

Int_t AliITSAlignMille2::CalcIntersectionPoint(Double_t *lpar, Double_t *gpar) {
  /// calculate track intersection point in local coordinates
  /// according with a given set of parameters (local(4) and global(6))
  /// and fill fPintLoc/Glo
  ///    local are:   pgx0, pgz0, ugx, ugz   OR   riemann fitters pars
  ///    global are:  tx,ty,tz,psi,theta,phi (Raff's delta angles in deg.)
  /// return 0 if success
  
  AliDebug(3,Form("lpar = %g %g %g %g %g\ngpar= %g %g %g %g %g %g\n",lpar[0],lpar[1],lpar[2],lpar[3],lpar[4],gpar[0],gpar[1],gpar[2],gpar[3],gpar[4],gpar[5]));
  AliDebug(3,Form("deltalpar = %g %g %g %g %g\n",lpar[0]-fLocalInitParam[0],lpar[1]-fLocalInitParam[1],lpar[2]-fLocalInitParam[2],lpar[3]-fLocalInitParam[3],lpar[4]-fLocalInitParam[4]));

  
  // prepare the TGeoHMatrix
  TGeoHMatrix *fTempHMat = fCurrentModule->GetSensitiveVolumeModifiedMatrix(fCluster.GetVolumeID(),gpar,
									    !fUseGlobalDelta);
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

Int_t AliITSAlignMille2::CalcDerivatives(Int_t paridx, Bool_t islpar) {
  /// calculate numerically (ROOT's style) the derivatives for
  /// local X intersection  and local Z intersection
  /// parlist: local  (islpar=kTRUE)  pgx0, pgz0, ugx0, ugz0  OR riemann's params
  ///          global (islpar=kFALSE) tx, ty, tz, psi, theta, phi (Raf's angles in deg)
  /// return 0 if success
  
  // copy initial parameters
  Double_t lpar[ITSMILLE2_NLOCAL];
  Double_t gpar[ITSMILLE2_NPARCH];
  for (Int_t i=0; i<ITSMILLE2_NLOCAL; i++) lpar[i]=fLocalInitParam[i];
  for (Int_t i=0; i<ITSMILLE2_NPARCH; i++) gpar[i]=fModuleInitParam[i];

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
  fDerivativeLoc[0] = h2*(4*d2 - d0)/3.;
  if (TMath::Abs(fDerivativeLoc[0]) < 1.0e-9) fDerivativeLoc[0] = 0.0;

  d0 = pintl4[2]-pintl1[2];
  d2 = 2.*(pintl3[2]-pintl2[2]);
  fDerivativeLoc[2] = h2*(4*d2 - d0)/3.;
  if (TMath::Abs(fDerivativeLoc[2]) < 1.0e-9) fDerivativeLoc[2]=0.0;

  AliDebug(3,Form("\n+++ derivatives +++ \n"));
  AliDebug(3,Form("+++ dXLoc/dpar = %g +++\n",fDerivativeLoc[0]));
  AliDebug(3,Form("+++ dZLoc/dpar = %g +++\n\n",fDerivativeLoc[0]));
  
  return 0;
}


Int_t AliITSAlignMille2::AddLocalEquation(Mille2Data &m) {
  /// Define local equation for current cluster in X and Z coor.
  /// and store them to memory
  /// return 0 if success
  int nLev = 0; 
  // store first interaction point
  if (CalcIntersectionPoint(fLocalInitParam, fModuleInitParam)) return -4;  
  for (Int_t i=0; i<3; i++) fPintLoc0[i]=fPintLoc[i];
  AliDebug(2,Form("Intesect. point: L( %f , %f , %f )",fPintLoc[0],fPintLoc[1],fPintLoc[2]));
  
  // calculate local derivatives numerically
  Bool_t zeroX = kTRUE;
  Bool_t zeroZ = kTRUE;
  //
  for (Int_t i=0; i<fNLocal; i++) {
    if (CalcDerivatives(i,kTRUE)) return -1;
    m.derlocX[i] = fDerivativeLoc[0];
    m.derlocZ[i] = fDerivativeLoc[2];
    if (zeroX) zeroX = fDerivativeLoc[0]==0;
    if (zeroZ) zeroZ = fDerivativeLoc[2]==0;
  }
  //  for (Int_t i=0; i<fNLocal; i++) AliDebug(2,Form("Local parameter %d - dXdpar = %g  - dZdpar = %g\n",i,dXdL[i],dZdL[i]));
  if (zeroX) {AliInfo("Aborting: zero local X derivatives!"); return -2;}
  if (zeroZ) {AliInfo("Aborting: zero local Z derivatives!"); return -2;}
  //
  //
  AliITSAlignMille2Module* endModule = fCurrentModule;
  //
  do {
    if (nLev==0 || !fUseGlobalDelta) zeroX = zeroZ = kTRUE;
    int shiftL = nLev*ITSMILLE2_NPARCH;
    for (Int_t i=0; i<ITSMILLE2_NPARCH; i++) {
      if (nLev==0 || !fUseGlobalDelta) {
	if (CalcDerivatives(i,kFALSE)) return -1;
	m.dergloX[shiftL + i] = fDerivativeLoc[0];
	m.dergloZ[shiftL + i] = fDerivativeLoc[2];
	if (zeroX) zeroX = fDerivativeLoc[0]==0;
	if (zeroZ) zeroZ = fDerivativeLoc[2]==0;      
      } 
      else { // for the global delta the deravites of different levels are the same
	m.dergloX[shiftL + i] = m.dergloX[shiftL + i - ITSMILLE2_NPARCH];
	m.dergloZ[shiftL + i] = m.dergloZ[shiftL + i - ITSMILLE2_NPARCH];
      }
    }
    //  for (Int_t i=0; i<ITSMILLE2_NPARCH; i++) AliDebug(2,Form("Global parameter %d - dXdpar = %g  - dZdpar = %g\n",i,dXdG[i],dZdG[i]));
    if (zeroX) {AliInfo("Aborting: zero global X derivatives!");return -2;}
    if (zeroZ) {AliInfo("Aborting: zero global Z derivatives!");return -2;}
    // set equation for Xloc coordinate
    m.moduleIDX[nLev] = fCurrentModule->GetUniqueID();
    nLev++;
    //
  } while( (fCurrentModule=fCurrentModule->GetParent()) );
  //
  // ok, can copy to m
  AliDebug(2,Form("Adding local equation X with fMeas=%.6f  and fSigma=%.6f",(fMeasLoc[0]-fPintLoc0[0]), fSigmaLoc[0]));
  m.measX = fMeasLoc[0]-fPintLoc0[0];
  m.sigmaX = fSigmaLoc[0];
  //
  AliDebug(2,Form("Adding local equation Z with fMeas=%.6f  and fSigma=%.6f",(fMeasLoc[2]-fPintLoc0[2]), fSigmaLoc[2]));
  m.measZ = fMeasLoc[2]-fPintLoc0[2];
  m.sigmaZ = fSigmaLoc[2];
  //
  m.levFilled = nLev;
  fCurrentModule = endModule;
  //
  return 0;
}

void AliITSAlignMille2::SetLocalEquations(const Mille2Data *marr, Int_t neq) {
  /// Set local equations with data stored in m
  /// return 0 if success
  //
  for (Int_t j=0; j<neq; j++) {
    //
    const Mille2Data &m = marr[j];
    //
    // set equation for Xloc coordinate
    AliDebug(2,Form("setting local equation X with fMeas=%.6f  and fSigma=%.6f",m.measX, m.sigmaX));
    for (int i=fNLocal; i--;) SetLocalDerivative( i, m.derlocX[i] );
    for (int il=m.levFilled;il--;) {
      GetMilleModule(m.moduleIDX[il])->IncNProcessedPoints();
      int hlev = m.moduleIDX[il]*ITSMILLE2_NPARCH;       // id of the supermodule 
      int llev = il*ITSMILLE2_NPARCH;
      for (int i=ITSMILLE2_NPARCH; i--;) SetGlobalDerivative( hlev+i, m.dergloX[llev+i] );
    }
    //
    fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, m.measX, m.sigmaX);  
    //
    // set equation for Zloc coordinate
    AliDebug(2,Form("setting local equation Z with fMeas=%.6f  and fSigma=%.6f",m.measZ, m.sigmaZ));
    for (int i=fNLocal; i--;) SetLocalDerivative( i, m.derlocZ[i] );
    for (int il=m.levFilled;il--;) {
      int hlev = m.moduleIDX[il]*ITSMILLE2_NPARCH;       // id of the supermodule 
      int llev = il*ITSMILLE2_NPARCH;
      for (int i=ITSMILLE2_NPARCH; i--;) SetGlobalDerivative(hlev+i, m.dergloZ[llev+i] );
    }
    fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, m.measZ, m.sigmaZ);  
  }
}

Int_t AliITSAlignMille2::GlobalFit(Double_t *parameters,Double_t *errors,Double_t *pulls) {
  /// Call global fit; Global parameters are stored in parameters
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return 0;
  }
  int res = fMillepede->GlobalFit(parameters,errors,pulls);
  AliInfo(Form("%s fitting global parameters!",res ? "Done":"Failed"));
  return res;
}

Double_t AliITSAlignMille2::GetParError(Int_t iPar) {
  /// Get error of parameter iPar
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return 0;
  }
  Double_t lErr = fMillepede->GetParError(iPar);
  return lErr;
}

void AliITSAlignMille2::PrintGlobalParameters() {
  /// Print global parameters
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->PrintGlobalParameters();
}

// //_________________________________________________________________________
Int_t AliITSAlignMille2::LoadSuperModuleFile(const Char_t *sfile)
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
  //
  Char_t st[250];
  char symname[150];
  UShort_t volid;
  TGeoHMatrix m;
  //
  for (Int_t i=0; i<nsma; i++) {
    AliAlignObjParams *a = (AliAlignObjParams*)sma->UncheckedAt(i);
    volid=a->GetVolUID();
    strcpy(st,a->GetSymName());
    a->GetMatrix(m);
    //
    sscanf(st,"%s",symname);
    // decode module list
    char *stp=strstr(st,"ModuleList:");
    if (!stp) return -3;
    stp += 11;
    int idx[2200];
    char spp[200]; int jp=0;
    char cl[20];
    strcpy(st,stp);
    int l=strlen(st);
    int j=0;
    int n=0;
    //
    while (j<=l) {
      if (st[j]==9 || st[j]==32 || st[j]==10 || st[j]==0) {
	spp[jp]=0;
	jp=0;
	if (strlen(spp)) {
	  int k=strcspn(spp,"-");
	  if (k<int(strlen(spp))) { // c'e' il -
	    strcpy(cl,&(spp[k+1]));
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
      volidsv[j]=AliITSAlignMille2Module::GetVolumeIDFromIndex(idx[j]);
      if (!volidsv[j]) {
	AliInfo(Form("Index %d not valid (range 0->%d)",idx[j],kMaxITSSensID));
	return -5;
      }
    }
    Int_t smindex=int(2198+volid-14336); // virtual index
    fSuperModule.AddAtAndExpand(new AliITSAlignMille2Module(smindex,volid,symname,&m,n,volidsv),fNSuperModules);
    //
    fNSuperModules++;
  }

  smf->Close();
  //
  return 0;
}

//_________________________________________________________________________
void AliITSAlignMille2::ConstrainModuleSubUnits(Int_t idm, Double_t val, UInt_t pattern)
{
  // require that sum of modifications for the childs of this module is = val, i.e.
  // the internal corrections moves the module as a whole by fixed value (0 by default).
  // pattern is the bit pattern for the parameters to constrain
  //
  static TObjArray childs;
  childs.Clear();
  //
  // build list of childs for this module
  int nChilds = 0;
  AliITSAlignMille2Module* parent = GetMilleModule(idm);
  if (!parent) return;
  //
  for (int i=fNModules;i--;) {
    AliITSAlignMille2Module* child = GetMilleModule(i);
    if (child->GetParent() == parent) childs.AddAtAndExpand(child,nChilds++);
  }
  if (nChilds<1) return;
  //
  int npc = 0;
  for (int ip=0;ip<kNParCh;ip++) {
    if ( !((pattern>>ip)&0x1) /*|| !parent->IsFreeDOF(ip)*/) continue;
    ResetLocalEquation();
    for (int ich=nChilds;ich--;) fGlobalDerivatives[childs[ich]->GetUniqueID()*kNParCh+ip] = 1.0;
    AddConstraint(fGlobalDerivatives,val);
    npc++;
  }
  //
  AliInfo(Form("Constrained %d params for %d submodules of module #%d: %s",npc,nChilds,idm,parent->GetName()));
  //
}

//_________________________________________________________________________
void AliITSAlignMille2::PostConstrainModuleSubUnitsMedian(Int_t idm, Double_t val, UInt_t pattern)
{
  // require that median of the modifications for the childs of this module is = val, i.e.
  // the internal corrections moves the module as a whole by fixed value (0 by default) 
  // module the outliers.
  // pattern is the bit pattern for the parameters to constrain
  // The difference between the mean and the median will be transfered to the parent
  //
  static TObjArray childs;
  childs.Clear();
  //
  // build list of childs for this module
  int nChilds = 0;
  AliITSAlignMille2Module* parent = GetMilleModule(idm);
  if (!parent) return;
  //
  for (int i=fNModules;i--;) {
    AliITSAlignMille2Module* child = GetMilleModule(i);
    if (child->GetParent() == parent) childs.AddAtAndExpand(child,nChilds++);
  }
  if (nChilds<1) return;
  //
  int npc = 0;
  double *deltas = fMillepede->GetDeltaPars();
  double *tmpArr = new double[nChilds]; 
  //
  for (int ip=0;ip<kNParCh;ip++) {
    if ( !((pattern>>ip)&0x1) /*|| !parent->IsFreeDOF(ip)*/) continue;
    // compute the median of the deltas
    for (int ich=nChilds;ich--;) tmpArr[ich] = deltas[childs[ich]->GetUniqueID()*kNParCh+ip];
    for (int ic0=0;ic0<nChilds;ic0++) // order the deltas 
      for (int ic1=ic0+1;ic1<nChilds;ic1++) 
	if (tmpArr[ic0]>tmpArr[ic1]) {double tv=tmpArr[ic0]; tmpArr[ic0]=tmpArr[ic1]; tmpArr[ic1]=tv;}
    //
    int kmed = nChilds/2;
    double median = (tmpArr[kmed]+tmpArr[nChilds-kmed-1])/2.;
    //
    for (int ich=nChilds;ich--;) deltas[childs[ich]->GetUniqueID()*kNParCh+ip] -= median - val;
    deltas[parent->GetUniqueID()*kNParCh+ip] += median - val;
    npc++;
  }
  delete[] tmpArr;  
  //
  AliInfo(Form("Applied median constraint to %d params for %d submodules of module #%d: %s",npc,nChilds,idm,parent->GetName()));
  //
}

//_________________________________________________________________________
void AliITSAlignMille2::ConstrainOrphans(Double_t val, UInt_t pattern)
{
  // require that median of the modifications for the supermodules which have no parents is = val, i.e.
  // the corrections moves the whole setup by fixed value (0 by default) modulo the outliers.
  // pattern is the bit pattern for the parameters to constrain
  //
  static TObjArray modSet;
  modSet.Clear();
  //
  // build list of childs for this module
  int nModules = 0;
  int *nFree = new int[kNParCh];
  for (int i=0;i<kNParCh;i++) nFree[i] = 0;
  //
  for (int i=fNModules;i--;) {
    AliITSAlignMille2Module* module = GetMilleModule(i);
    if (module->GetParent()) continue;   // skip this    
    for (int ip=0;ip<kNParCh;ip++) if ( ((pattern>>ip)&0x1) && module->IsFreeDOF(ip)) nFree[ip]++;
    modSet.AddAtAndExpand(module,nModules++);
  }
  if (nModules<1) return;
  //
  int npc = 0;
  for (int ip=0;ip<kNParCh;ip++) {
    if (nFree[ip]<1) continue; // nothing to do
    ResetLocalEquation();
    for (int ich=nModules;ich--;) {
      AliITSAlignMille2Module* module = (AliITSAlignMille2Module*) modSet[ich];
      if ( !((pattern>>ip)&0x1) || !module->IsFreeDOF(ip)) continue;
      fGlobalDerivatives[module->GetUniqueID()*kNParCh+ip] = 1.0;
    }
    AddConstraint(fGlobalDerivatives,val);
    npc++;
  }
  //
  delete[] nFree;
  AliInfo(Form("Constrained %d params for %d orphan modules",npc,nModules));
  //
}
//_________________________________________________________________________
void AliITSAlignMille2::PostConstrainOrphansMedian(Double_t val, UInt_t pattern)
{
  // require that sum of modifications for the supermodules which have no parents is = val, i.e.
  // the corrections moves the whole setup by fixed value (0 by default).
  // pattern is the bit pattern for the parameters to constrain
  //
  static TObjArray modSet;
  modSet.Clear();
  //
  // build list of childs for this module
  int nModules = 0;
  int *nFree = new int[kNParCh];
  for (int i=0;i<kNParCh;i++) nFree[i] = 0;
  //
  for (int i=fNModules;i--;) {
    AliITSAlignMille2Module* module = GetMilleModule(i);
    if (module->GetParent()) continue;   // skip this    
    for (int ip=0;ip<kNParCh;ip++) if ( ((pattern>>ip)&0x1) && module->IsFreeDOF(ip)) nFree[ip]++;
    modSet.AddAtAndExpand(module,nModules++);
  }
  if (nModules<1) return;
  //
  int npc = 0;
  double *deltas = fMillepede->GetDeltaPars();
  double *tmpArr = new double[nModules]; 
  //
  for (int ip=0;ip<kNParCh;ip++) {
    if (nFree[ip]<1) continue; // nothing to do
    // compute the median of the deltas
    for (int ich=nModules;ich--;) tmpArr[ich] = deltas[modSet[ich]->GetUniqueID()*kNParCh+ip];
    for (int ic0=0;ic0<nModules;ic0++) // order the deltas 
      for (int ic1=ic0+1;ic1<nModules;ic1++) 
	if (tmpArr[ic0]>tmpArr[ic1])  {double tv=tmpArr[ic0]; tmpArr[ic0]=tmpArr[ic1]; tmpArr[ic1]=tv;};
    //
    int kmed = nModules/2;
    double median = (tmpArr[kmed]+tmpArr[nModules-kmed-1])/2.;
    //    
    for (int ich=nModules;ich--;) {
      AliITSAlignMille2Module* module = (AliITSAlignMille2Module*) modSet[ich];
      if ( !((pattern>>ip)&0x1) || !module->IsFreeDOF(ip)) continue;
      deltas[module->GetUniqueID()*kNParCh+ip] -= median - val;
    }
    npc++;
  }
  //
  delete[] nFree;
  delete[] tmpArr;
  AliInfo(Form("Applied median constraint to %d params for %d orphan modules",npc,nModules));
  //
}

//_________________________________________________________________________
void AliITSAlignMille2::ConstrainLinComb(const Int_t *vidLst, const Float_t *wghLst, Int_t nmd, Double_t val, UInt_t pattern)
{
  // require that the linear combinations of the nmd modules (refered by their volume ID) from the 
  // modList with the coefficients wghLst adds up to val.
  // pattern is the bit pattern for the parameters to constrain.
  //
  static TObjArray modSet;
  modSet.Clear();
  //
  // build list of childs for this module
  int nModules = 0;
  int *nFree = new int[kNParCh];
  for (int i=0;i<kNParCh;i++) nFree[i] = 0;
  //
  for (int imd=nmd;imd--;) {
    UShort_t vid = (UShort_t)vidLst[imd];
    for (int i=fNModules;i--;) {
      AliITSAlignMille2Module* module = GetMilleModule(i);
      if (module->GetVolumeID() == vid) {
	modSet.AddAtAndExpand(module,nModules++);
	for (int ip=0;ip<kNParCh;ip++) if ( ((pattern>>ip)&0x1) && module->IsFreeDOF(ip)) nFree[ip]++;
	break;
      }
    }
  }
  //
  if (nModules != nmd) {
    AliInfo(Form("Error: constraint for %d modules requested but %d are found",nmd,nModules));
    delete[] nFree;
    return;
  }
  int npc = 0;
  for (int ip=0;ip<kNParCh;ip++) {
    if (nFree[ip]<1) continue; // nothing to do
    ResetLocalEquation();
    for (int ich=nModules;ich--;) {
      AliITSAlignMille2Module* module = (AliITSAlignMille2Module*) modSet[ich];
      if ( !((pattern>>ip)&0x1) || !module->IsFreeDOF(ip)) continue;
      fGlobalDerivatives[module->GetUniqueID()*kNParCh+ip] = wghLst[ich];
    }
    AddConstraint(fGlobalDerivatives,val);
    npc++;
  }
  //
  delete[] nFree;
  AliInfo(Form("Constrained %d params for linerar combination of %d  modules",npc,nModules));
  //
}
