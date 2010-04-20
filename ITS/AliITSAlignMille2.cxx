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
//
//  Interface to AliMillePede2 alignment class for the ALICE ITS detector
// 
//  ITS specific alignment class which interface to AliMillepede.   
//  For each track ProcessTrack calculates the local and global derivatives
//  at each hit and fill the corresponding local equations. Provide methods for
//  fixing or constraning detection elements for best results. 
// 
//  author M. Lunardon (thanks to J. Castillo), ruben.shahoyan@cern.ch
//-----------------------------------------------------------------------------

#include <TFile.h>
#include <TGrid.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <TGeoManager.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TCollection.h>
#include <TGeoPhysicalNode.h>
#include <TMap.h>
#include <TObjString.h>
#include <TString.h>
#include "AliITSAlignMille2.h"
#include "AliITSgeomTGeo.h"
#include "AliGeomManager.h"
#include "AliMillePede2.h"
#include "AliTrackPointArray.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"
#include "AliTrackFitterRieman.h"
#include "AliITSAlignMille2Constraint.h"
#include "AliITSAlignMille2ConstrArray.h"
#include "AliITSresponseSDD.h"
#include "AliITSTPArrayFit.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSDriftSpeedArraySDD.h"
#include "AliESDVertex.h"

ClassImp(AliITSAlignMille2)

const Char_t* AliITSAlignMille2::fgkRecKeys[] = {
  "OCDB_PATH",
  "OCDB_SPECIFIC",
  "GEOMETRY_FILE",
  "SUPERMODULE_FILE",
  "CONSTRAINTS_REFERENCE_FILE",
  "PREALIGNMENT_FILE",
  "PRECALIBSDD_FILE",
  "PREVDRIFTSDD_FILE",
  "INITCALBSDD_FILE",
  "INITVDRIFTSDD_FILE",
  "INITDELTA_FILE",
  "SET_GLOBAL_DELTAS",
  "CONSTRAINT_LOCAL",
  "MODULE_VOLUID",
  "MODULE_INDEX",
  "SET_PSEUDO_PARENTS",
  "SET_TRACK_FIT_METHOD",
  "SET_MINPNT_TRA",
  "SET_NSTDDEV",
  "SET_RESCUT_INIT",
  "SET_RESCUT_OTHER",
  "SET_LOCALSIGMAFACTOR",
  "SET_STARTFAC",
  "SET_FINALFAC",
  "SET_B_FIELD",
  "SET_SPARSE_MATRIX",
  "REQUIRE_POINT",
  "CONSTRAINT_ORPHANS",
  "CONSTRAINT_SUBUNITS",
  "APPLY_CONSTRAINT",
  "SET_EXTRA_CLUSTERS_MODE",
  "SET_USE_TPAFITTER",
  "SET_USE_LOCAL_YERROR",
  "SET_MIN_POINTS_PER_MODULE",
  "SET_USE_SDDVDCORRMULT",
  "SET_WEIGHT_PT",
  "SET_USE_DIAMOND",
  "SET_SAME_SDDT0"
};

const Char_t AliITSAlignMille2::fgkXYZ[] = "XYZ";

//========================================================================================================

AliITSAlignMille2* AliITSAlignMille2::fgInstance = 0;  
Int_t              AliITSAlignMille2::fgInstanceID = 0;

//________________________________________________________________________________________________________
AliITSAlignMille2::AliITSAlignMille2(const Char_t *configFilename,TList *userInfo ) 
: TObject(),
  fMillepede(0),
  fStartFac(16.), 
  fFinalFac(1.), 
  fResCutInitial(100.), 
  fResCut(100.),
  fNGlobal(0),
  fNLocal(4),
  fNStdDev(3),
  fIsMilleInit(kFALSE),
  fAllowPseudoParents(kFALSE),
  //
  fTPAFitter(0),
  fCurrentModule(0),
  fTrack(0),
  fTrackBuff(0),
  fCluster(),
  fCurrentSensID(-1),
  fClusLoc(12*3),
  fClusGlo(12*3),
  fClusSigLoc(12*3),
  fGlobalDerivatives(0),
  fMeasLoc(0),
  fMeasGlo(0),
  fSigmaLoc(0),
  fConstrPT(-1),
  fConstrPTErr(-1),
  fConstrCharge(0),
  //
  fMinNPtsPerTrack(3),
  fIniTrackParamsMeth(1),
  fTotBadLocEqPoints(0),
  fRieman(0),
  //
  fConstraints(0),
  fCacheMatrixOrig(kMaxITSSensID+1),
  fCacheMatrixCurr(kMaxITSSensID+1),
  //
  fUseGlobalDelta(kFALSE),
  fTempExcludedModule(-1),
  //
  fIniUserInfo(userInfo),
  fIniDeltaPath(""),
  fIniSDDRespPath(""),
  fPreCalSDDRespPath(""),
  fIniSDDVDriftPath(""),
  fPreSDDVDriftPath(""),
  fGeometryPath(""),
  fPreDeltaPath(""),
  fConstrRefPath(""),
  fDiamondPath(""),
  fGeoManager(0),
  fIsConfigured(kFALSE),
  fPreAlignQF(0),
//
  fIniRespSDD(0),
  fPreRespSDD(0),
  fIniVDriftSDD(0),
  fPreVDriftSDD(0),
  fSegmentationSDD(0),
  fPrealignment(0),
  fConstrRef(0),
  fMilleModule(2),
  fSuperModule(2),
  fNModules(0),
  fNSuperModules(0),
  fUsePreAlignment(kFALSE),
  fUseLocalYErr(kFALSE),
  fBOn(kFALSE),
  fBField(0.0),
  fDataType(kCosmics),
  fMinPntPerSens(0),
  fBug(0),
  fMilleVersion(2),
  fExtraClustersMode(0),
  fTrackWeight(1),
  fWeightPt(0.),
  fIsSDDVDriftMult(kFALSE),
  fDiamond(),
  fDiamondI(),
  fUseDiamond(kFALSE),
  fDiamondPointID(-1),
  fDiamondModID(-1)
{
  /// main constructor that takes input from configuration file
  for (int i=3;i--;) fSigmaFactor[i] = 1.0;
  //
  // new RS
  for (int i=0;i<3;i++) {
    
  }
  for (int itp=0;itp<kNDataType;itp++) {
    fRequirePoints[itp] = kFALSE;
    for (Int_t i=0; i<6; i++) {
      fNReqLayUp[itp][i]=0;
      fNReqLayDown[itp][i]=0;
      fNReqLay[itp][i]=0;
    }
    for (Int_t i=0; i<3; i++) {
      fNReqDetUp[itp][i]=0;
      fNReqDetDown[itp][i]=0;
      fNReqDet[itp][i]=0;
    }
  }
  //
  //  if (ProcessUserInfo(userInfo)) exit(1);
  //
  fDiamond.SetVolumeID(kVtxSensVID);
  fDiamondI.SetVolumeID(kVtxSensVID);
  float xyzd[3] = {0,0,0};
  float covd[6] = {1,0,0,1,0,1e4}; 
  fDiamond.SetXYZ(xyzd,covd); // dummy diamond
  covd[5] = 1e-4;
  fDiamondI.SetXYZ(xyzd,covd);
  //
  Int_t lc=LoadConfig(configFilename);
  if (lc) {
    AliError(Form("Error %d loading configuration from %s",lc,configFilename));
    exit(1);
  }
  //
  fMillepede = new AliMillePede2();  
  fgInstance = this;
  fgInstanceID++;
  //
}

//________________________________________________________________________________________________________
AliITSAlignMille2::~AliITSAlignMille2()
{
  /// Destructor
  delete fMillepede;
  delete[] fGlobalDerivatives;
  delete fRieman;
  delete fPrealignment;
  delete fConstrRef;
  delete fPreRespSDD;
  delete fIniRespSDD;
  delete fSegmentationSDD;
  delete fIniVDriftSDD;
  delete fPreVDriftSDD;
  delete fTPAFitter;
  fCacheMatrixOrig.Delete();
  fCacheMatrixCurr.Delete();
  fTrackBuff.Delete();
  fConstraints.Delete();
  fMilleModule.Delete();
  fSuperModule.Delete();
  if (--fgInstanceID==0) fgInstance = 0;
}

///////////////////////////////////////////////////////////////////////

//________________________________________________________________________________________________________
TObjArray* AliITSAlignMille2::GetConfigRecord(FILE* stream, TString& recTitle, TString& recOpt, Bool_t rew)
{
  // read new record from config file
  TString record;
  static TObjArray* recElems = 0;
  if (recElems) {delete recElems; recElems = 0;}
  recOpt = "";
  //
  TString keyws = recTitle;
  if (!keyws.IsNull()) {
    keyws.ToUpper();
    //    keyws += " ";
  }
  while (record.Gets(stream)) {
    int cmt=record.Index("#"); 
    if (cmt>=0) record.Remove(cmt);  // skip comment
    record.ReplaceAll("\t"," ");
    record.ReplaceAll("\r"," ");
    record.Remove(TString::kBoth,' '); 
    if (record.IsNull()) continue;      // nothing to decode 
    if (!keyws.IsNull() && !record.BeginsWith(keyws.Data())) continue; // specific record was requested
    //
    recElems = record.Tokenize(" ");
    recTitle = recElems->At(0)->GetName();
    recTitle.ToUpper();
    recOpt = recElems->GetLast()>0 ? recElems->At(1)->GetName() : "";
    break;
  }
  if (rew || !recElems) rewind(stream);
  return recElems;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::CheckConfigRecords(FILE* stream)
{  
  TString record,recTitle;
  int lineCnt = 0;
  rewind(stream);
  while (record.Gets(stream)) {
    int cmt=record.Index("#"); 
    lineCnt++;
    if (cmt>=0) record.Remove(cmt);  // skip comment
    record.ReplaceAll("\t"," ");
    record.ReplaceAll("\r"," ");
    record.Remove(TString::kBoth,' ');
    if (record.IsNull()) continue;   // nothing to decode  
    // extract keyword
    int spc = record.Index(" ");
    if (spc>0) recTitle = record(0,spc);
    else     recTitle = record;
    recTitle.ToUpper();
    Bool_t strOK = kFALSE;
    for (int ik=kNKeyWords;ik--;) if (recTitle == fgkRecKeys[ik]) {strOK = kTRUE; break;}
    if (strOK) continue;
    //
    AliError(Form("Unknown keyword %s at line %d",
		  recTitle.Data(),lineCnt));
    return -1;
    //
  }
  //
  rewind(stream);
  return 0;
}


//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::LoadConfig(const Char_t *cfile)
{  
  // return 0 if success
  //        1 if error in module index or voluid
  //
  AliInfo(Form("Loading MillePede2 configuration from %s",cfile));
  AliCDBManager::Instance()->SetCacheFlag(kFALSE);
  FILE *pfc=fopen(cfile,"r");
  if (!pfc) return -1;
  //
  TString record,recTitle,recOpt,recExt;
  Int_t nrecElems,irec;
  TObjArray *recArr=0;
  //
  fNModules = 0;
  Bool_t stopped = kFALSE;
  //
  if (CheckConfigRecords(pfc)<0) return -1;
  //
  while(1) { 
    //
    // ============= 1: we read some important records in predefined order ================
    //  
    recTitle = fgkRecKeys[kOCDBDefaultPath];
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) && !recOpt.IsNull() ) {
      AliInfo(Form("Configuration sets OCDB default storage to %s",recOpt.Data()));
      AliCDBManager::Instance()->SetDefaultStorage( gSystem->ExpandPathName(recOpt.Data()) );
      TObjString* objStr = (TObjString*)AliCDBManager::Instance()->GetStorageMap()->GetValue("default");
      if (!objStr) {stopped = kTRUE; break;}
      objStr->SetUniqueID(1); // mark as user set
    }
    //
    if (fIniUserInfo && ProcessUserInfo(fIniUserInfo)) { AliError("Failed to process intial User Info"); stopped = kTRUE; break;}
    //  
    recTitle = fgkRecKeys[kGeomFile];
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) ) fGeometryPath = gSystem->ExpandPathName(recOpt.Data()); 
    if ( InitGeometry() ) { AliError("Failed to find/load Geometry"); stopped = kTRUE; break;}
    //
    // Do we use new TrackPointArray fitter ?
    recTitle = fgkRecKeys[kTPAFitter];
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) ) fTPAFitter = new AliITSTPArrayFit(kNLocal);
    //
    recTitle = fgkRecKeys[kSuperModileFile];
    if ( !GetConfigRecord(pfc,recTitle,recOpt,1) || 
	 recOpt.IsNull()                         || 
	 gSystem->ExpandPathName(recOpt)         ||
	 gSystem->AccessPathName(recOpt.Data())  ||
	 LoadSuperModuleFile(recOpt.Data()))
      { AliError("Failed to find/load SuperModules"); stopped = kTRUE; break;}
    //
    recTitle = fgkRecKeys[kConstrRefFile];      // LOCAL_CONSTRAINTS are defined wrt these deltas
    if ( (recArr = GetConfigRecord(pfc,recTitle,recOpt,1)) ) {
      if (recOpt.IsNull() || recOpt=="IDEAL") SetConstraintWrtRef( "IDEAL" );
      else {
	for (int i=2;i<=recArr->GetLast();i++) {recOpt += " "; recOpt += recArr->At(i)->GetName();} // in case of OCDB string
	if ( SetConstraintWrtRef(recOpt.Data()) )
	  { AliError("Failed to load reference deltas for local constraints"); stopped = kTRUE; break;}
      }
    }
    //
    //	 
    recTitle = fgkRecKeys[kInitDeltaFile];
    if ( (recArr = GetConfigRecord(pfc,recTitle,recOpt,1))  && !recOpt.IsNull() ) {
      for (int i=2;i<=recArr->GetLast();i++) {recOpt += " "; recOpt += recArr->At(i)->GetName();} // in case of OCDB string
      fIniDeltaPath = recOpt;
      gSystem->ExpandPathName(fIniDeltaPath);
      AliInfo(Form("Configuration sets Production Deltas to %s",fIniDeltaPath.Data()));
    }
    //
    // if initial deltas were provided, load them, apply to geometry and store are "original" matrices
    if (CacheMatricesOrig()) {stopped = kTRUE; break;}
    //	 
    recTitle = fgkRecKeys[kPreDeltaFile];
    if ( (recArr = GetConfigRecord(pfc,recTitle,recOpt,1)) ) {
      if (!recOpt.IsNull()) {
	for (int i=2;i<=recArr->GetLast();i++) {recOpt += " "; recOpt += recArr->At(i)->GetName();} // in case of OCDB string
	fPreDeltaPath = recOpt;
	gSystem->ExpandPathName(fPreDeltaPath);
      }
      else if (!fIniDeltaPath.IsNull()) {
	AliInfo("PreAlignment Deltas keyword is present but empty, will set to Init Deltas");
	fPreDeltaPath = fIniDeltaPath;	
      }
      AliInfo(Form("Configuration sets PreAlignment Deltas to %s",fPreDeltaPath.Data()));
    }
    if (LoadDeltas(fPreDeltaPath,fPrealignment)) {stopped = kTRUE; break;}
    if (fPrealignment && ApplyToGeometry()) {stopped = kTRUE; break;}
    //
    recTitle = fgkRecKeys[ kInitCalSDDFile ];
    if ( (recArr = GetConfigRecord(pfc,recTitle,recOpt,1)) && !recOpt.IsNull()) {
      for (int i=2;i<=recArr->GetLast();i++) {recOpt += " "; recOpt += recArr->At(i)->GetName();} // in case of OCDB string
      fIniSDDRespPath = recOpt;
      gSystem->ExpandPathName(fIniSDDRespPath);
      AliInfo(Form("Configuration sets Production SDD Response to %s",fIniSDDRespPath.Data()));
    }
    if (LoadSDDResponse(fIniSDDRespPath, fIniRespSDD) ) {stopped = kTRUE; break;}
    //
    recTitle = fgkRecKeys[kPreCalSDDFile];
    if ( (recArr = GetConfigRecord(pfc,recTitle,recOpt,1)) ) {
      if (!recOpt.IsNull()) {
	for (int i=2;i<=recArr->GetLast();i++) {recOpt += " "; recOpt += recArr->At(i)->GetName();} // in case of OCDB string
	fPreCalSDDRespPath = recOpt;
	gSystem->ExpandPathName(fPreCalSDDRespPath);
      }
      else if (!fIniSDDRespPath.IsNull()) {
	AliInfo("PreCalibration SDD response keyword is present but empty, will set to Init SDD repsonse");
	fPreCalSDDRespPath = fIniSDDRespPath;	
      }
      AliInfo(Form("Configuration sets PreCalibration SDD Response to %s",fPreCalSDDRespPath.Data()));
    }
    //
    if (LoadSDDResponse(fPreCalSDDRespPath, fPreRespSDD) ) {stopped = kTRUE; break;}
    //
    //
    recTitle = fgkRecKeys[ kInitVDriftSDDFile ];
    if ( (recArr = GetConfigRecord(pfc,recTitle,recOpt,1)) && !recOpt.IsNull()) {
      for (int i=2;i<=recArr->GetLast();i++) {recOpt += " "; recOpt += recArr->At(i)->GetName();} // in case of OCDB string
      fIniSDDVDriftPath = recOpt;
      gSystem->ExpandPathName(fIniSDDVDriftPath);
      AliInfo(Form("Configuration sets Production SDD VDrift to %s",fIniSDDVDriftPath.Data()));
    }
    if (LoadSDDVDrift(fIniSDDVDriftPath, fIniVDriftSDD) ) {stopped = kTRUE; break;}
    //
    recTitle = fgkRecKeys[ kPreVDriftSDDFile ];
    if ( (recArr = GetConfigRecord(pfc,recTitle,recOpt,1)) && !recOpt.IsNull()) {
      for (int i=2;i<=recArr->GetLast();i++) {recOpt += " "; recOpt += recArr->At(i)->GetName();} // in case of OCDB string
      fPreSDDVDriftPath = recOpt;
      gSystem->ExpandPathName(fPreSDDVDriftPath);
      AliInfo(Form("Configuration sets PreCalibration SDD VDrift to %s",fPreSDDVDriftPath.Data()));
      if (LoadSDDVDrift(fPreSDDVDriftPath, fPreVDriftSDD) ) {stopped = kTRUE; break;}
    }
    //
    recTitle = fgkRecKeys[ kGlobalDeltas ];
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) ) SetUseGlobalDelta(kTRUE);
    //
    recTitle = fgkRecKeys[ kUseDiamond ];
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) ) {
      if (!GetUseGlobalDelta()) {
	AliError("Diamond constraint is supported only for Global Frame mode");
	stopped = kTRUE; 
	break;
      }
      fUseDiamond = kTRUE;
      if (!recOpt.IsNull()) {
	fDiamondPath = recOpt;
	gSystem->ExpandPathName(fDiamondPath);
	AliInfo(Form("Configuration sets Diamond constraint to %s",fDiamondPath.Data()));
      }
    }
    // =========== 2: see if there are local gaussian constraints defined =====================
    //            Note that they should be loaded before the modules declaration
    //
    recTitle = fgkRecKeys[ kConstrLocal ];
    while( (recArr=GetConfigRecord(pfc,recTitle,recOpt,0)) ) {
      nrecElems = recArr->GetLast()+1;
      if (recOpt.IsFloat()) {stopped = kTRUE; break;} // wrong name
      if (GetConstraint(recOpt.Data())) {
	AliError(Form("Existing constraint %s repeated",recOpt.Data()));
	stopped = kTRUE; break;
      }
      recExt = recArr->At(2)->GetName();
      if (!recExt.IsFloat()) {stopped = kTRUE; break;}
      double val = recExt.Atof();      
      recExt = recArr->At(3)->GetName();
      if (!recExt.IsFloat()) {stopped = kTRUE; break;}
      double err = recExt.Atof();      
      int nwgh = nrecElems - 4;
      double *wgh = new double[nwgh];
      for (nwgh=0,irec=4;irec<nrecElems;irec++) {
	recExt = recArr->At(irec)->GetName();
	if (!recExt.IsFloat()) {stopped = kTRUE; break;}
	wgh[nwgh++] = recExt.Atof();
      }
      if (stopped) {delete[] wgh; break;}
      //
      ConstrainLocal(recOpt.Data(),wgh,nwgh,val,err);
      delete[] wgh;
      //
    } // end while for loop over local constraints
    if (stopped) break;
    //
    // =========== 3: now read modules to align ===================================
    //
    rewind(pfc);
    // create fixed modules
    for (int j=0; j<fNSuperModules; j++) {
      AliITSAlignMille2Module* proto = GetSuperModule(j);
      if (!proto->IsAlignable()) continue;
      AliITSAlignMille2Module* mod = new AliITSAlignMille2Module(*proto);
      // the matrix might be updated in case some prealignment was applied, check 
      TGeoHMatrix* mup = AliGeomManager::GetMatrix(mod->GetName());
      if (mup) *(mod->GetMatrix()) = *mup;
      fMilleModule.AddAtAndExpand(mod,fNModules);
      mod->SetGeomParamsGlobal(fUseGlobalDelta);
      mod->SetUniqueID(fNModules++);
      mod->SetNotInConf(kTRUE);
    }
    CreateVertexModule();
    //
    while( (recArr=GetConfigRecord(pfc,recTitle="",recOpt,0)) ) {
      if (!(recTitle==fgkRecKeys[ kModVolID ] || recTitle==fgkRecKeys[ kModIndex ])) continue;
      // Expected format: MODULE id tolX tolY tolZ tolPsi tolTh tolPhi [[sigX sigY sigZ]  extra params]
      // where tol* is the tolerance (sigma) for given DOF. 0 means fixed
      // sig* is the scaling parameters for the errors of the clusters of this module
      // extra params are defined for specific modules, e.g. t0 and vdrift corrections of SDD
      //
      nrecElems = recArr->GetLast()+1;
      if (nrecElems<2 || !recOpt.IsDigit()) {stopped = kTRUE; break;}
      int idx = recOpt.Atoi(); 
      UShort_t voluid =  (idx<=kMaxITSSensID) ? GetModuleVolumeID(idx) : idx;
      AliITSAlignMille2Module* mod = 0;
      //
      if (voluid>=kMinITSSupeModuleID) { // custom supermodule
	mod = GetMilleModuleByVID(voluid);
	if (!mod) { // need to create
	  for (int j=0; j<fNSuperModules; j++) {
	    if (voluid==GetSuperModule(j)->GetVolumeID()) {
	      mod = new AliITSAlignMille2Module(*GetSuperModule(j));
	      // the matrix might be updated in case some prealignment was applied, check 
	      TGeoHMatrix* mup = AliGeomManager::GetMatrix(mod->GetName());
	      if (mup) *(mod->GetMatrix()) = *mup;
	      fMilleModule.AddAtAndExpand(mod,fNModules);
	      mod->SetGeomParamsGlobal(fUseGlobalDelta);
	      mod->SetUniqueID(fNModules++);
	      break;
	    }	
	  }
	}
	mod->SetNotInConf(kFALSE);
      }
      else if (idx<=kMaxITSSensVID) {
	mod = new AliITSAlignMille2Module(voluid);
	fMilleModule.AddAtAndExpand(mod,fNModules);
	mod->SetGeomParamsGlobal(fUseGlobalDelta);
	mod->SetUniqueID(fNModules++);
      }
      if (!mod) {stopped = kTRUE; break;}  // bad volid
      //
      // geometry variation settings
      for (int i=0;i<AliITSAlignMille2Module::kMaxParGeom;i++) {
	irec = i+2;
	if (irec >= nrecElems) break;
	recExt = recArr->At(irec)->GetName();
	if (!recExt.IsFloat()) {stopped = kTRUE; break;}
	mod->SetFreeDOF(i, recExt.Atof() );	
      }
      if (stopped) break;
      //
      // scaling factors for cluster errors
      // first set default ones
      for (int i=0;i<3;i++) mod->SetSigmaFactor(i, fSigmaFactor[i]);	
      for (int i=0;i<3;i++) {
	irec = i+8;
	if (irec >= nrecElems) break;
	recExt = recArr->At(irec)->GetName();
	if (!recExt.IsFloat()) {stopped = kTRUE; break;}
	mod->SetSigmaFactor(i, recExt.Atof() );	
      }     
      if (stopped) break;
      //
      // now comes special detectors treatment
      if (mod->IsSDD()) {
	double vl = 0;
	if (nrecElems>11) {
	  recExt = recArr->At(11)->GetName();
	  if (recExt.IsFloat()) vl = recExt.Atof();
	  else {stopped = kTRUE; break;}
	  irec = 11;
	}
	mod->SetFreeDOF(AliITSAlignMille2Module::kDOFT0,vl);
	//
	Bool_t cstLR = kFALSE;
	for (int lr=0;lr<2;lr++) { // left right side vdrift corrections
	  vl = 0;
	  if (nrecElems>12+lr) {
	    recExt = recArr->At(12+lr)->GetName();
	    if (recExt.IsFloat()) vl = recExt.Atof();
	    else {stopped = kTRUE; break;}
	    irec = 12+lr;
	  }
	  mod->SetFreeDOF(lr==0 ? AliITSAlignMille2Module::kDOFDVL : AliITSAlignMille2Module::kDOFDVR,vl);
	  if (lr==1 && vl>=10) cstLR = kTRUE;  // the right side should be constrained to left one 
	}
	if (cstLR) mod->SetVDriftLRSame();
      }
      //
      mod->EvaluateDOF();
      //
      // now check if there are local constraints on this module
      for (++irec;irec<nrecElems;irec++) {
	recExt = recArr->At(irec)->GetName();
	if (recExt.IsFloat()) {stopped=kTRUE;break;}
	AliITSAlignMille2ConstrArray* cstr = (AliITSAlignMille2ConstrArray*)GetConstraint(recExt.Data());
	if (!cstr) {
	  AliInfo(Form("No Local constraint %s was declared",recExt.Data())); 
	  stopped=kTRUE; 
	  break;
	}
	cstr->AddModule(mod);
      }
      if (stopped) break;
    } // end while for loop over modules
    if (stopped) break;
    //
    if (fNModules==0) {AliError("Failed to find any MODULE"); stopped = kTRUE; break;}  
    BuildHierarchy();  // preprocess loaded modules
    //
    // =========== 4: the rest may come in arbitrary order =======================================
    rewind(pfc);
    while ( (recArr=GetConfigRecord(pfc,recTitle="",recOpt,0))!=0 ) {
      //
      nrecElems = recArr->GetLast()+1;
      //
      // some simple flags -----------------------------------------------------------------------
      //
      if      (recTitle == fgkRecKeys[ kPseudoParents ])  SetAllowPseudoParents(kTRUE);
      //
      // some optional parameters ----------------------------------------------------------------
      else if (recTitle == fgkRecKeys[ kTrackFitMethod ]) {
	if (recOpt.IsNull() || !recOpt.IsDigit() ) {stopped = kTRUE; break;}
	SetInitTrackParamsMeth(recOpt.Atoi());
      }
      //
      else if (recTitle == fgkRecKeys[ kMinPntTrack ]) {
	if (recOpt.IsNull() || !recOpt.IsDigit() ) {stopped = kTRUE; break;}
	fMinNPtsPerTrack = recOpt.Atoi();
      }
      //
      else if (recTitle == fgkRecKeys[ kNStDev ]) {
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fNStdDev = (Int_t)recOpt.Atof();
      }
      //
      else if (recTitle == fgkRecKeys[ kResCutInit  ]) {
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fResCutInitial = recOpt.Atof();
      }
      //
      else if (recTitle == fgkRecKeys[ kResCutOther ]) {
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fResCut = recOpt.Atof();
      }
      //
      else if (recTitle == fgkRecKeys[ kLocalSigFactor ]) { //-------------------------
	for (irec=0;irec<3;irec++) if (nrecElems>irec+1) {
	    fSigmaFactor[irec] = ((TObjString*)recArr->At(irec+1))->GetString().Atof();
	    if (fSigmaFactor[irec]<=0.) stopped = kTRUE;
	  }
	if (stopped) break; 
      }
      //
      else if (recTitle == fgkRecKeys[ kStartFactor ]) {        //-------------------------
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fStartFac = recOpt.Atof();
      }
      //
      else if (recTitle == fgkRecKeys[ kFinalFactor ]) {        //-------------------------
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fFinalFac = recOpt.Atof();
      }
      //
      // pepo2708909
      else if (recTitle == fgkRecKeys[ kExtraClustersMode ]) {        //-------------------------
	if (recOpt.IsNull() || !recOpt.IsDigit() ) {stopped = kTRUE; break;}
	fExtraClustersMode = recOpt.Atoi();
      }
      // endpepo270809
      //
      else if (recTitle == fgkRecKeys[ kBField ]) {         //-------------------------
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	SetBField( recOpt.Atof() );
      }
      //
      else if (recTitle == fgkRecKeys[ kSDDVDCorrMult ]) {         //-------------------------
	SetSDDVDCorrMult( recOpt.IsNull() || (recOpt.IsFloat() && (recOpt.Atof())>-0.5) ); 
      }
      //
      else if (recTitle == fgkRecKeys[ kWeightPt ]) {         //-------------------------
	double wgh = 1;
	if (!recOpt.IsNull()) {
	  if (!recOpt.IsFloat()) {stopped = kTRUE; break;}
	  else wgh = recOpt.Atof();
	}
	SetWeightPt(wgh);
      }
      //
      else if (recTitle == fgkRecKeys[ kSparseMatrix ]) {   // matrix solver type
	//
	AliMillePede2::SetGlobalMatSparse(kTRUE);
	if (recOpt.IsNull()) continue;
	// solver type and settings
	if      (recOpt == "MINRES") AliMillePede2::SetIterSolverType( AliMinResSolve::kSolMinRes );
	else if (recOpt == "FGMRES") AliMillePede2::SetIterSolverType( AliMinResSolve::kSolFGMRes );
	else {stopped = kTRUE; break;}
	//
	if (nrecElems>=3) { // preconditioner type
	  recExt = recArr->At(2)->GetName();
	  if (!recExt.IsDigit()) {stopped = kTRUE; break;}
	  AliMillePede2::SetMinResPrecondType( recExt.Atoi() );
	}
	//
	if (nrecElems>=4) { // tolerance
	  recExt = recArr->At(3)->GetName();
	  if (!recExt.IsFloat()) {stopped = kTRUE; break;}
	  AliMillePede2::SetMinResTol( recExt.Atof() );
	}
	//
	if (nrecElems>=5) { // maxIter
	  recExt = recArr->At(4)->GetName();
	  if (!recExt.IsDigit()) {stopped = kTRUE; break;}
	  AliMillePede2::SetMinResMaxIter( recExt.Atoi() );
	}	
      }
      //
      else if (recTitle == fgkRecKeys[ kRequirePoint ]) {       //-------------------------
	// syntax:   REQUIRE_POINT where ndet updw nreqpts
	//    where = LAYER or DETECTOR
	//    ndet = detector number: 1-6 for LAYER and 1-3 for DETECTOR (SPD=1, SDD=2, SSD=3)
	//    updw = 1 for Y>0, -1 for Y<0, 0 if not specified
	//    nreqpts = minimum number of points of that type
	if (nrecElems>=5) {
	  recOpt.ToUpper();
	  int lr = ((TObjString*)recArr->At(2))->GetString().Atoi() - 1;
	  int hb = ((TObjString*)recArr->At(3))->GetString().Atoi();
	  int np = ((TObjString*)recArr->At(4))->GetString().Atoi();
	  //
	  int rtp = -1; // use for run type
	  if (nrecElems>5) {
	    TString tpstr = ((TObjString*)recArr->At(5))->GetString();
	    if ( tpstr.Contains("cosmics",TString::kIgnoreCase) ) rtp = kCosmics;
	    else if ( tpstr.Contains("collision",TString::kIgnoreCase) ) rtp = kCollision;
	    else {stopped = kTRUE; break;}
	  }
	  //
	  int tpmn= rtp<0 ? 0 : rtp;
	  int tpmx= rtp<0 ? kNDataType-1 : rtp;
	  for (int itp=tpmn;itp<=tpmx;itp++) {
	    fRequirePoints[itp]=kTRUE;
	    if (recOpt == "LAYER") {
	      if (lr<0 || lr>5) {stopped = kTRUE; break;}
	      if (hb>0) fNReqLayUp[itp][lr]=np;
	      else if (hb<0) fNReqLayDown[itp][lr]=np;
	      else fNReqLay[itp][lr]=np;
	    }
	    else if (recOpt == "DETECTOR") {
	      if (lr<0 || lr>2) {stopped = kTRUE; break;}
	      if (hb>0) fNReqDetUp[itp][lr]=np;
	      else if (hb<0) fNReqDetDown[itp][lr]=np;
	      else fNReqDet[itp][lr]=np;
	    }
	    else {stopped = kTRUE; break;}
	  }
	  if (stopped) break;
	}
	else {stopped = kTRUE; break;}
      }
      //
      // global constraints on the subunits/orphans 
      else if (recTitle == fgkRecKeys[ kConstrOrphans ]) {    //------------------------
	// expect CONSTRAINT_ORPHANS MEAN/MEDIAN Value parID0 ... parID1 ...
	if (nrecElems<4) {stopped = kTRUE; break;}
	recExt = recArr->At(2)->GetName();
	if (!recExt.IsFloat()) {stopped = kTRUE; break;}
	double val = recExt.Atof();
	UInt_t pattern = 0;
	for (irec=3;irec<nrecElems;irec++) { // read params to constraint
	  recExt = recArr->At(irec)->GetName();
	  if (!recExt.IsDigit()) {stopped = kTRUE; break;}
	  pattern |= 0x1 << recExt.Atoi();
	}
	if (stopped) break;
	if      (recOpt == "MEAN")   ConstrainOrphansMean(val,pattern);
	else if (recOpt == "MEDIAN") ConstrainOrphansMedian(val,pattern);
	else {stopped = kTRUE; break;}
      }
      //
      else if (recTitle == fgkRecKeys[ kConstrSubunits ]) {    //------------------------
	// expect CONSTRAINT_SUBUNITS MEAN/MEDIAN Value parID0 ... parID1 ... VolID1 ... VolIDn - VolIDm
	if (nrecElems<5) {stopped = kTRUE; break;}
	recExt = recArr->At(2)->GetName();
	if (!recExt.IsFloat()) {stopped = kTRUE; break;}
	double val = recExt.Atof();
	UInt_t pattern = 0;
	for (irec=3;irec<nrecElems;irec++) { // read params to constraint
	  recExt = recArr->At(irec)->GetName();
	  if (!recExt.IsDigit()) {stopped = kTRUE; break;}
	  int parid = recExt.Atoi();
	  if (parid<kMaxITSSensID) pattern |= 0x1 << recExt.Atoi();
	  else break;           // list of params is over 
	}
	if (stopped) break;
	//
	Bool_t meanC;
	if      (recOpt == "MEAN")   meanC = kTRUE;
	else if (recOpt == "MEDIAN") meanC = kFALSE;
	else    {stopped = kTRUE; break;}
	//
	int curID = -1;
	int rangeStart = -1;
	for (;irec<nrecElems;irec++) { // read modules to apply this constraint
	  recExt = recArr->At(irec)->GetName();
	  if (recExt == "-") {rangeStart = curID; continue;}  // range is requested
	  else if (!recExt.IsDigit()) {stopped = kTRUE; break;}
	  else curID = recExt.Atoi();
	  //
	  if (curID<=kMaxITSSensID) curID = GetModuleVolumeID(curID);
	  // this was a range start or single 
	  int start;
	  if (rangeStart>=0) {start = rangeStart+1; rangeStart=-1;} // continue the range
	  else start = curID;  // create constraint either for single module (or 1st in the range)
	  for (int id=start;id<=curID;id++) {
	    int id0 = IsVIDDefined(id);
	    if (id0<0) {AliDebug(3,Form("Undefined module %d requested in the SubUnits constraint, skipping",id)); continue;}
	    if (meanC) ConstrainModuleSubUnitsMean(id0,val,pattern);
	    else       ConstrainModuleSubUnitsMedian(id0,val,pattern);
	  }
	}
	if (rangeStart>=0) stopped = kTRUE; // unfinished range
	if (stopped) break;
      } 
      // 
      // association of modules with local constraints
      else if (recTitle == fgkRecKeys[ kApplyConstr ]) {            //------------------------
	// expect APPLY_CONSTRAINT NAME [NAME1...] [VolID1 ... VolIDn - VolIDm]
	if (nrecElems<3) {stopped = kTRUE; break;}
	int nmID0=-1,nmID1=-1;
	for (irec=1;irec<nrecElems;irec++) { // find the range of constraint names
	  recExt = recArr->At(irec)->GetName();
	  if (recExt.IsFloat()) break;
	  // check if such a constraint was declared
	  if (!GetConstraint(recExt.Data())) {
	    AliInfo(Form("No Local constraint %s was declared",recExt.Data())); 
	    stopped=kTRUE; 
	    break;
	  }
	  if (nmID0<0) nmID0 = irec;
	  nmID1 = irec;
	}
	if (stopped) break;
	//
	if (irec>=nrecElems) {stopped = kTRUE; break;} // no modules provided
	//
	// now read the list of modules to constrain
	int curID = -1;
	int rangeStart = -1;
	for (;irec<nrecElems;irec++) { // read modules to apply this constraint
	  recExt = recArr->At(irec)->GetName();
	  if (recExt == "-") {rangeStart = curID; continue;}  // range is requested
	  else if (!recExt.IsDigit()) {stopped = kTRUE; break;}
	  else curID = recExt.Atoi();
	  //
	  if (curID<=kMaxITSSensID) curID = GetModuleVolumeID(curID);
	  //
	  // this was a range start or single 
	  int start;
	  if (rangeStart>=0) {start = rangeStart+1; rangeStart=-1;} // continue the range
	  else start = curID;  // create constraint either for single module (or 1st in the range)
	  for (int id=start;id<=curID;id++) {
	    AliITSAlignMille2Module *md = GetMilleModuleByVID(id);
	    if (!md) {AliDebug(3,Form("Undefined module %d requested in the Local constraint, skipping",id)); continue;}
	    for (int nmid=nmID0;nmid<=nmID1;nmid++) 
	      ((AliITSAlignMille2ConstrArray*)GetConstraint(recArr->At(nmid)->GetName()))->AddModule(md);
	  }
	}
	if (rangeStart>=0) stopped = kTRUE; // unfinished range
	if (stopped) break;
      }
      //
      // request of the same T0 for group of SDD modules
      else if (recTitle == fgkRecKeys[ kSameSDDT0 ]) {            //------------------------
	// expect SET_SAME_SDDT0 [SensID1 ... SensIDn - SensIDm]
	if (nrecElems<3) {stopped = kTRUE; break;}
	//
	// now read the list of modules to constrain
	int curID = -1;
	int rangeStart = -1;
	AliITSAlignMille2ConstrArray *cstrT0 = new AliITSAlignMille2ConstrArray("SDDT0",0,0,0,0);
	int naddM = 0;
	cstrT0->SetPattern(BIT(AliITSAlignMille2Module::kDOFT0));
	for (irec=1;irec<nrecElems;irec++) { // read modules to apply this constraint
	  recExt = recArr->At(irec)->GetName();
	  if (recExt == "-") {rangeStart = curID; continue;}  // range is requested
	  else if (!recExt.IsDigit()) {stopped = kTRUE; break;}
	  else curID = recExt.Atoi();
	  //
	  if (curID<kSDDoffsID || curID>=kSDDoffsID+kNSDDmod) {stopped = kTRUE; break;}
	  //
	  // this was a range start or single 
	  int start;
	  if (rangeStart>=0) {start = rangeStart+1; rangeStart=-1;} // continue the range
	  else start = curID;  // create constraint either for single module (or 1st in the range)
	  for (int id=start;id<=curID;id++) {
	    int vid = AliITSAlignMille2Module::GetVolumeIDFromIndex(id);
	    if (vid<=1) {AliDebug(3,Form("Undefined module index %d requested in the SAME_SDDT0 constraint, skipping",id)); continue;}
	    AliITSAlignMille2Module *md = GetMilleModuleByVID(vid);
	    if (!md) {AliDebug(3,Form("Undefined module %d requested in the Local constraint, skipping",id)); continue;}
	    cstrT0->AddModule(md,kFALSE);
	    naddM++;
	  }	  
	}
	if (rangeStart>=0) stopped = kTRUE; // unfinished range
	if (stopped) break;
	if (naddM<2) delete cstrT0;
	else {
	  cstrT0->SetConstraintID(GetNConstraints());
	  fConstraints.Add(cstrT0);
	}
      }
      //
      // Do we use new local Y errors?
      else if (recTitle == fgkRecKeys[ kUseLocalYErr ]) {
	// expect SET_TPAFITTER 
	fUseLocalYErr = kTRUE;
      }
      //
      else if (recTitle == fgkRecKeys[ kMinPointsSens ]) {         //-------------------------
	if (recOpt.IsNull() || !recOpt.IsDigit() ) {stopped = kTRUE; break;}
	SetMinPointsPerSensor( recOpt.Atoi() );
      }
      //
      else if (recTitle == fgkRecKeys[ kOCDBSpecificPath ]) {         //-------------------------
	if (recOpt.IsNull() || nrecElems<3 ) {stopped = kTRUE; break;}
	AliCDBManager::Instance()->SetSpecificStorage(recOpt.Data(), gSystem->ExpandPathName(recArr->At(2)->GetName()));
	AliInfo(Form("Configuration sets OCDB specific storage %s to %s",recOpt.Data(),recArr->At(2)->GetName()));
	TObjString *pths = (TObjString*)AliCDBManager::Instance()->GetStorageMap()->GetValue(recOpt.Data());
	if (!pths) { stopped = kTRUE; break; }
	pths->SetUniqueID(1); // mark as set by user
      }
      //
      else continue; // already processed record
      //
    } // end of while loop 4 over the various params 
    //
    break;
  } // end of while(1) loop 
  //
  fclose(pfc);
  if (!fDiamondPath.IsNull() && IsDiamondUsed() && LoadDiamond(fDiamondPath) ) stopped = kTRUE;
  if (stopped) {
    AliError(Form("Failed on record %s %s ...\n",recTitle.Data(),recOpt.Data()));
    return -1;
  }
  //
  if (CacheMatricesCurr()) return -1;
  SetUseLocalYErrors(fUseLocalYErr); // YErr used only with TPAFitter 
  fSegmentationSDD = new AliITSsegmentationSDD();
  //
  fIsConfigured = kTRUE;
  return 0;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::BuildHierarchy()
{
  // build the hieararhy of the modules to align
  //
  if (!GetUseGlobalDelta() && PseudoParentsAllowed()) {
    AliInfo("PseudoParents mode is allowed only when the deltas are global\n"
	    "Since Deltas are local, switching to NoPseudoParents");
    SetAllowPseudoParents(kFALSE);
  }
  // set parent/child relationship for modules to align
  AliInfo("Setting parent/child relationships\n");
  //
  // 1) child -> parent reference
  for (int ipar=0;ipar<fNModules;ipar++) {
    AliITSAlignMille2Module* parent = GetMilleModule(ipar);
    if (parent->IsSensor()) continue; // sensor cannot be a parent
    //
    for (int icld=0;icld<fNModules;icld++) {
      if (icld==ipar) continue;
      AliITSAlignMille2Module* child = GetMilleModule(icld);
      if (!child->BelongsTo(parent)) continue;
      // child cannot have more sensors than the parent
      if (child->GetNSensitiveVolumes() > parent->GetNSensitiveVolumes()) continue;
      //
      AliITSAlignMille2Module* parOld = child->GetParent();
      // is this parent candidate closer than the old parent ? 
      if (parOld && parOld->GetNSensitiveVolumes()<parent->GetNSensitiveVolumes()) continue; // parOld is closer
      child->SetParent(parent);
    }
    //
  }
  //
  // add parent -> children reference
  for (int icld=0;icld<fNModules;icld++) {
    AliITSAlignMille2Module* child = GetMilleModule(icld);
    AliITSAlignMille2Module* parent = child->GetParent();
    if (parent) parent->AddChild(child);
  }  
  //
  // reorder the modules in such a way that parents come first
  for (int icld=0;icld<fNModules;icld++) {
    AliITSAlignMille2Module* child  = GetMilleModule(icld);
    AliITSAlignMille2Module* parent; 
    while ( (parent=child->GetParent()) &&  (parent->GetUniqueID()>child->GetUniqueID()) ) {
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
  // Go over the child->parent chain and mark modules with explicitly provided sensors.
  // If the sensors of the unit are explicitly declared, all undeclared sensors are 
  // suppresed in this unit.
  for (int icld=fNModules;icld--;) {
    AliITSAlignMille2Module* child = GetMilleModule(icld);
    AliITSAlignMille2Module* parent = child->GetParent();
    if (!parent) continue;
    //
    // check if this parent was already processed
    if (!parent->AreSensorsProvided()) {
      parent->DelSensitiveVolumes();
      parent->SetSensorsProvided(kTRUE);
    }
    // reattach sensors to parent
    for (int isc=child->GetNSensitiveVolumes();isc--;) {
      UShort_t senVID = child->GetSensVolVolumeID(isc);
      if (!parent->IsIn(senVID)) parent->AddSensitiveVolume(senVID);
    }
  }
  //
}

// pepo
//________________________________________________________________________________________________________
void AliITSAlignMille2::SetCurrentModule(Int_t id)
{
  // set the current supermodule
  // new meaning
  if (fMilleVersion>=2) {
    fCurrentModule = GetMilleModule(id);
    return;
  }
  // old meaning
  if (fMilleVersion<=1) {
    Int_t index=id;
    /// set as current the SuperModule that contains the 'index' sens.vol.
    if (index<0 || index>2197) {
      AliInfo("index does not correspond to a sensitive volume!");
      return;
    }
    UShort_t voluid=AliITSAlignMille2Module::GetVolumeIDFromIndex(index);
    Int_t k=IsContained(voluid);
    if (k>=0){
      fCurrentSensID = index;
      fCluster.SetVolumeID(voluid);
      fCluster.SetXYZ(0,0,0);
      InitModuleParams();
    }
    else
      AliInfo(Form("module %d not defined\n",index));    
  }
}
// endpepo
//________________________________________________________________________________________________________
void AliITSAlignMille2::SetRequiredPoint(Char_t* where, Int_t ndet, Int_t updw, Int_t nreqpts,Int_t runtype) 
{
  // set minimum number of points in specific detector or layer
  // where = LAYER or DETECTOR
  // ndet = detector number: 1-6 for LAYER and 1-3 for DETECTOR (SPD=1, SDD=2, SSD=3)
  // updw = 1 for Y>0, -1 for Y<0, 0 if not specified
  // nreqpts = minimum number of points of that type
  ndet--;
  int tpmn= runtype<0 ? 0 : runtype;
  int tpmx= runtype<0 ? kNDataType-1 : runtype;
  //
  for (int itp=tpmn;itp<=tpmx;itp++) {
    fRequirePoints[itp]=kTRUE;
    if (strstr(where,"LAYER")) {
      if (ndet<0 || ndet>5) return;
      if (updw>0) fNReqLayUp[itp][ndet]=nreqpts;
      else if (updw<0) fNReqLayDown[itp][ndet]=nreqpts;
      else fNReqLay[itp][ndet]=nreqpts;
    }
    else if (strstr(where,"DETECTOR")) {
      if (ndet<0 || ndet>2) return;
      if (updw>0) fNReqDetUp[itp][ndet]=nreqpts;
      else if (updw<0) fNReqDetDown[itp][ndet]=nreqpts;
      else fNReqDet[itp][ndet]=nreqpts;	
    }
  }
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::GetModuleIndex(const Char_t *symname) 
{
  /// index from symname
  if (!symname) return -1;
  for (Int_t i=0;i<=kMaxITSSensID; i++) {
    if (!strcmp(symname,AliITSgeomTGeo::GetSymName(i))) return i;
  }
  return -1;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::GetModuleIndex(UShort_t voluid) 
{
  /// index from volume ID
  AliGeomManager::ELayerID lay = AliGeomManager::VolUIDToLayer(voluid);
  if (lay<1|| lay>6) return -1;
  Int_t idx=Int_t(voluid)-2048*lay;
  if (idx>=AliGeomManager::LayerSize(lay)) return -1;
  for (Int_t ilay=1; ilay<lay; ilay++) 
    idx += AliGeomManager::LayerSize(ilay);
  return idx;
}

//________________________________________________________________________________________________________
UShort_t AliITSAlignMille2::GetModuleVolumeID(const Char_t *symname) 
{
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

//________________________________________________________________________________________________________
UShort_t AliITSAlignMille2::GetModuleVolumeID(Int_t index) 
{
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

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::InitGeometry() 
{
  /// initialize geometry
  AliInfo("Loading initial geometry");
  if (!fGeometryPath.IsNull() && gSystem->AccessPathName(fGeometryPath.Data()) ) {
    AliError(Form("Explicitly provided geometry file %s is not accessible",fGeometryPath.Data()));
    return -1;
  }
  //
  AliGeomManager::LoadGeometry(fGeometryPath.Data());
  fGeoManager = AliGeomManager::GetGeometry();
  if (!fGeoManager) {
    AliInfo("Couldn't initialize geometry");
    return -1;
  }
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::SetConstraintWrtRef(const char* reffname) 
{
  // Load the global deltas from this file. The local gaussian constraints on some modules 
  // will be defined with respect to the deltas from this reference file, converted to local
  // delta format. Note: conversion to local format requires reloading the geometry!
  //
  AliInfo(Form("Loading reference deltas for local constraints from %s",reffname));
  if (!fGeoManager) return -1; 
  fConstrRefPath = reffname;
  if (fConstrRefPath == "IDEAL") { // the reference is the ideal geometry, just create dummy reference array
    fConstrRef = new TClonesArray("AliAlignObjParams",1);
    return 0;
  }
  if (LoadDeltas(fConstrRefPath,fConstrRef)) return -1;
  //
  // we need ideal geometry to convert global deltas to local ones
  if (fUsePreAlignment) {
    AliError("The call of SetConstraintWrtRef must be done before application of the prealignment");
    return -1;
  }
  //
  AliInfo("Converting global reference deltas to local ones");
  Int_t nprea = fConstrRef->GetEntriesFast();
  for (int ix=0; ix<nprea; ix++) {
    AliAlignObjParams *preo=(AliAlignObjParams*) fConstrRef->At(ix);
    if (!preo->ApplyToGeometry()) return -1;
  }
  //
  // now convert the global reference deltas to local ones
  for (int i=fConstrRef->GetEntriesFast();i--;) {
    AliAlignObjParams *preo = (AliAlignObjParams*)fConstrRef->At(i);
    TGeoHMatrix * mupd = AliGeomManager::GetMatrix(preo->GetSymName());
    if (!mupd) {  // this is not alignable entry, need to look in the supermodules
      for (int im=fNSuperModules;im--;) {
	AliITSAlignMille2Module* mod = GetSuperModule(im);
	if ( strcmp(mod->GetName(), preo->GetSymName()) ) continue;
	mupd = mod->GetMatrix();
	break;
      }
      if (!mupd) {
	AliError(Form("Failed to find the volume for reference %s",preo->GetSymName()));
	return -1;
      }
    } 
    TGeoHMatrix preMat;
    preo->GetMatrix(preMat);                     //  Delta_Glob
    TGeoHMatrix tmpMat    = *mupd;               //  Delta_Glob * Delta_Glob_Par * M
    preMat.MultiplyLeft( &tmpMat.Inverse() );    //  M^-1 * Delta_Glob_Par^-1 = (Delta_Glob_Par * M)^-1
    tmpMat.MultiplyLeft( &preMat );              //  (Delta_Glob_Par * M)^-1 * Delta_Glob * Delta_Glob_Par * M = Delta_loc
    preo->SetMatrix(tmpMat);     // local corrections 
  }
  //
  // we need to reload the geometry spoiled by this reference deltas...
  delete fGeoManager;
  AliInfo("Reloading initial geometry");
  return InitGeometry();
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::Init()
{
  // perform global initialization
  //
  if (fIsMilleInit) {
    AliInfo("Millepede has been already initialized!");
    return;
  }
  // range constraints in such a way that the childs are constrained before their parents
  // orphan constraints come last
  for (int ic=0;ic<GetNConstraints();ic++) {
    for (int ic1=ic+1;ic1<GetNConstraints();ic1++) {
      AliITSAlignMille2Constraint *cst0 = GetConstraint(ic);
      AliITSAlignMille2Constraint *cst1 = GetConstraint(ic1);
      if (cst0->GetModuleID()<cst1->GetModuleID()) {
	// swap
	fConstraints[ic] = cst1;
	fConstraints[ic1] = cst0;
      }
    }
  }
  //
  if (!GetUseGlobalDelta()) {
    AliInfo("ATTENTION: The parameters are defined in the local frame, no check for degeneracy will be done");
    for (int imd=fNModules;imd--;) {
      AliITSAlignMille2Module* mod = GetMilleModule(imd);
      int npar = mod->GetNParTot();
      // the parameter may have max 1 free instance, otherwise the equations are underdefined
      for (int ipar=0;ipar<npar;ipar++) {
	if (!mod->IsFreeDOF(ipar)) continue;
	mod->SetParOffset(ipar,fNGlobal++);
      }
    }
  }
  else {
    // init millepede, decide which parameters are to be fitted explicitly
    for (int imd=fNModules;imd--;) {
      AliITSAlignMille2Module* mod = GetMilleModule(imd);
      int npar = mod->GetNParTot();
      // the parameter may have max 1 free instance, otherwise the equations are underdefined
      for (int ipar=0;ipar<npar;ipar++) {
	if (!mod->IsFreeDOF(ipar)) continue;  // fixed
	//
	int nFreeInstances = 0;
	//
	AliITSAlignMille2Module* parent = mod;
	Bool_t cstMeanMed=kFALSE,cstGauss=kFALSE;
	//
	Bool_t addToFit = kFALSE;	
	// the parameter may be ommitted from explicit fit (if PseudoParentsAllowed is true) if
	// 1) it is not explicitly constrained or its does not participate in Gaussian constraint
	// 2) the same applies to all of its parents
	// 3) it has at least 1 unconstrained direct child
	while(parent) {
	  if (!parent->IsFreeDOF(ipar)) {parent = parent->GetParent(); continue;}
	  nFreeInstances++;
	  if (IsParModConstrained(parent,ipar, cstMeanMed, cstGauss)) nFreeInstances--;
	  if (cstGauss) addToFit = kTRUE;
	  parent = parent->GetParent();
	}
	if (nFreeInstances>1) {
	  AliError(Form("Parameter#%d of module %s\nhas %d free instances in the "
			"unconstrained parents\nSystem is undefined",ipar,mod->GetName(),nFreeInstances));
	  exit(1);
	}
	//
	// i) Are PseudoParents allowed?
	if (!PseudoParentsAllowed()) addToFit = kTRUE;
	// ii) check if this module has no child with such a free parameter. Since the order of this check 
	// goes from child to parent, by this moment such a parameter must have been already added
	else if (!IsParModFamilyVaried(mod,ipar))  addToFit = kTRUE;  // no varied children at all
	else if (!IsParFamilyFree(mod,ipar,1))     addToFit = kTRUE;  // no unconstrained direct children
	// otherwise the value of this parameter can be extracted from simple contraint and the values of 
	// the relevant parameters of its children the fit is done. Hence it is not included
	if (!addToFit) continue;
	//
	// shall add this parameter to explicit fit
	//	printf("Adding %s %d -> %d\n",mod->GetName(), ipar, fNGlobal);
	mod->SetParOffset(ipar,fNGlobal++);
      }
    }
  }
  //
  AliInfo(Form("Initializing Millepede with %d gpar, %d lpar and %d stddev ...",fNGlobal, kNLocal, fNStdDev));
  fGlobalDerivatives = new Double_t[fNGlobal];
  memset(fGlobalDerivatives,0,fNGlobal*sizeof(Double_t));
  //
  fMillepede->InitMille(fNGlobal,kNLocal,fNStdDev,fResCut,fResCutInitial);
  fMillepede->SetMinPntValid(fMinPntPerSens);
  fIsMilleInit = kTRUE;
  //
  ResetLocalEquation();    
  AliInfo("Parameters initialized to zero");
  //
  /// Fix non free parameters
  for (Int_t i=0; i<fNModules; i++) {
    AliITSAlignMille2Module* mod = GetMilleModule(i);
    for (Int_t j=0; j<mod->GetNParTot(); j++) {
      if (mod->GetParOffset(j)<0) continue; // not varied
      FixParameter(mod->GetParOffset(j),mod->GetParConstraint(j));
      fMillepede->SetParamGrID(i, mod->GetParOffset(j));
    }
  }
  //
  // Set iterations
  if (fStartFac>1) fMillepede->SetIterations(fStartFac);    
  if (fFinalFac>1) fMillepede->SetChi2CutRef(fFinalFac);    
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::AddConstraint(Double_t *par, Double_t value, Double_t sigma) 
{
  /// Constrain equation defined by par to value
  if (!fIsMilleInit) Init();
  fMillepede->SetGlobalConstraint(par, value, sigma);
  AliInfo("Adding constraint");
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::InitGlobalParameters(Double_t *par) 
{
  /// Initialize global parameters with par array
  if (!fIsMilleInit) Init();
  fMillepede->SetGlobalParameters(par);
  AliInfo("Init Global Parameters");
}

//________________________________________________________________________________________________________ 
void AliITSAlignMille2::FixParameter(Int_t iPar, Double_t value) 
{
  /// Parameter iPar is encourage to vary in [-value;value]. 
  /// If value == 0, parameter is fixed
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->SetParSigma(iPar, value);
  if (IsZero(value)) AliInfo(Form("Parameter %i Fixed", iPar));
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ResetLocalEquation()
{
  /// Reset the derivative vectors
  for(int i=kNLocal;i--;)  fLocalDerivatives[i] = 0.0;
  memset(fGlobalDerivatives, 0, fNGlobal*sizeof(double) );
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::ApplyToGeometry() 
{
  // apply prealignment to ideal geometry
  Int_t nprea = fPrealignment->GetEntriesFast();
  AliInfo(Form("Array of prealignment deltas: %d entries",nprea));
  //
  for (int ix=0; ix<nprea; ix++) {
    AliAlignObjParams *preo=(AliAlignObjParams*) fPrealignment->At(ix);
    Int_t index=AliITSAlignMille2Module::GetIndexFromVolumeID(preo->GetVolUID());
    if (index>=0) {
      if (index>=fPreAlignQF.GetSize()) fPreAlignQF.Set(index+10);
      fPreAlignQF[index] = (int) preo->GetUniqueID()+1;
    }
    if (!preo->ApplyToGeometry()) {
      AliError(Form("Failed on ApplyToGeometry at %s",preo->GetSymName()));
      return -1;
    }
  }
  //
  fUsePreAlignment = kTRUE;
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::GetPreAlignmentQualityFactor(Int_t index) const
{
  // quality factors from prealignment
  if (!fUsePreAlignment || index<0 || index>=fPreAlignQF.GetSize()) return -1;
  return fPreAlignQF[index]-1;
}

//________________________________________________________________________________________________________
AliTrackPointArray *AliITSAlignMille2::PrepareTrack(const AliTrackPointArray *atp) 
{
  /// create a new AliTrackPointArray keeping only defined modules
  /// move points according to a given prealignment, if any
  /// sort alitrackpoints w.r.t. global Y direction, if selected
  const Double_t kRad2L[6] = {5*5,10*10,18*18,30*30,40*40,60*60};
  const Float_t kSensSigY2[6] = {200e-4*200e-4/12, 200e-4*200e-4/12, 
				 300e-4*300e-4/12, 300e-4*300e-4/12, 
				 300e-4*300e-4/12, 300e-4*300e-4/12}; // thickness^2/12
  //
  fTrack = NULL;
  Int_t   idx[20];
  Short_t lrID[20];
  Int_t npts=atp->GetNPoints();
  if (npts<fMinNPtsPerTrack) return NULL;
  TGeoHMatrix hcov;
  //
  /// checks if AliTrackPoints belong to defined modules
  Int_t ngoodpts=0;
  Int_t intidx[20];
  for (int j=0; j<npts; j++) {
    intidx[j] = GetRequestedModID(atp->GetVolumeID()[j]);
    if (intidx[j]<0) continue;
    ngoodpts++;
    Float_t xx=atp->GetX()[j];
    Float_t yy=atp->GetY()[j];
    Float_t r=xx*xx + yy*yy;
    int lay;
    for (lay=0;lay<6;lay++) if (r<kRad2L[lay]) break;
    if (lay>5) continue;
    lrID[j] = lay;
  }
  //
  AliDebug(3,Form("Number of points in defined modules: %d out of %d",ngoodpts,npts));

  // pepo270809
  Int_t nextra=0;
  // extra clusters selection mode  
  if (fExtraClustersMode) {
    // 1 = keep one cluster, remove randomly the extra
    // 2 = keep one cluster, remove the internal one
    // 10 = keep tracks only if at least one extra is present
    
    int iextra1[20],iextra2[20],layovl[20];
    // extra clusters mapping
    for (Int_t ipt=0; ipt<npts; ipt++) {
      if (intidx[ipt]<0) continue; // looks only defined modules...
      float p1x=atp->GetX()[ipt];
      float p1y=atp->GetY()[ipt];
      float p1z=atp->GetZ()[ipt];
      int lay1=int(AliGeomManager::VolUIDToLayer(atp->GetVolumeID()[ipt]));
      float r1 = p1x*p1x + p1y*p1y;
      UShort_t volid1=atp->GetVolumeID()[ipt];
      
      for (int ik=ipt+1; ik<npts; ik++) {
	if (intidx[ik]<0) continue;
	// compare point ipt with next ones
	int lay2=int(AliGeomManager::VolUIDToLayer(atp->GetVolumeID()[ik]));
	// check if same layer
	if (lay2 != lay1) continue;
	UShort_t volid2=atp->GetVolumeID()[ik];
	// check if different module
	if (volid1 == volid2) continue;

	float p2x=atp->GetX()[ik];
	float p2y=atp->GetY()[ik];
	float p2z=atp->GetZ()[ik];
	float r2 = p2x*p2x + p2y*p2y;	
	float dr= (p1x-p2x)*(p1x-p2x) + (p1y-p2y)*(p1y-p2y) + (p1z-p2z)*(p1z-p2z);
	
	// looks for pairs with dr<1 cm, same layer but different module
	if (dr<1.0) {
	  // extra1 is the one with smaller radius in rphi plane
	  if (r1<r2) {
	    iextra1[nextra]=ipt;
	    iextra2[nextra]=ik;
	  }
	  else {
	    iextra1[nextra]=ik;
	    iextra2[nextra]=ipt;
	  }
	  layovl[nextra]=lay1;	  
	  nextra++;
	}
      }
    } // end overlaps mapping
    
    // mode=1: keep only one clusters and remove the other randomly
    if (fExtraClustersMode==1 && nextra) {
      for (int ie=0; ie<nextra; ie++) {
	if (gRandom->Rndm()<0.5) 
	  intidx[iextra1[ie]]=-1;
	else
	  intidx[iextra2[ie]]=-1;	  
      }
    }

    // mode=2: keep only one clusters and remove the other...
    if (fExtraClustersMode==2 && nextra) {
      for (int ie=0; ie<nextra; ie++) {
	if (layovl[ie]==1) intidx[iextra2[ie]]=-1;
	else if (layovl[ie]==2) intidx[iextra1[ie]]=-1;
	else intidx[iextra1[ie]]=-1;	  
      }
    }

    // mode=10: reject track if no overlaps are present
    if (fExtraClustersMode==10 && nextra==0) {
      AliInfo("Track with no extra clusters: rejected!");
      return NULL;
    }
    
    // recalculate ngoodpts
    ngoodpts=0;
    for (int i=0; i<npts; i++) {
      if (intidx[i]>=0) ngoodpts++;
    }
  }
  // endpepo270809

  // reject track if not enough points are left
  if (ngoodpts<fMinNPtsPerTrack) {
    AliDebug(2,"Track with not enough points!");
    return NULL;
  }
  // >> RS
  AliTrackPoint p;
  // check points in specific places
  if (fRequirePoints[fDataType]) {
    Int_t nlayup[6],nlaydown[6],nlay[6];
    Int_t ndetup[3],ndetdown[3],ndet[3];
    for (Int_t j=0; j<6; j++) {nlayup[j]=0; nlaydown[j]=0; nlay[j]=0;}
    for (Int_t j=0; j<3; j++) {ndetup[j]=0; ndetdown[j]=0; ndet[j]=0;}
    
    for (int i=0; i<npts; i++) {
      // skip not defined points
      if (intidx[i]<0) continue;
      //      
      Float_t yy=atp->GetY()[i];
      int lay = lrID[i];
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
    //
    // checks minimum values
    Bool_t isok=kTRUE;
    for (Int_t j=0; j<6; j++) {
      if (nlayup[j]<fNReqLayUp[fDataType][j]) isok=kFALSE; 
      if (nlaydown[j]<fNReqLayDown[fDataType][j]) isok=kFALSE; 
      if (nlay[j]<fNReqLay[fDataType][j]) isok=kFALSE; 
    }
    for (Int_t j=0; j<3; j++) {
      if (ndetup[j]<fNReqDetUp[fDataType][j]) isok=kFALSE; 
      if (ndetdown[j]<fNReqDetDown[fDataType][j]) isok=kFALSE; 
      if (ndet[j]<fNReqDet[fDataType][j]) isok=kFALSE; 
    }
    if (!isok) {
      AliDebug(2,Form("Track does not meet all location point requirements!"));
      return NULL;
    }
  }
  // build a new track with (sorted) (prealigned) good points
  // pepo200709
  //fTrack = (AliTrackPointArray*)fTrackBuff[ngoodpts-fMinNPtsPerTrack];
  Int_t addVertex = IsTypeCollision()&&IsDiamondUsed() ? 1 : 0;
  fTrack = (AliTrackPointArray*)fTrackBuff[ngoodpts + addVertex ];
  if (!fTrack) {
    fTrack = new AliTrackPointArray(ngoodpts + addVertex);
    //    fTrackBuff.AddAtAndExpand(fTrack,ngoodpts-fMinNPtsPerTrack);
    fTrackBuff.AddAtAndExpand(fTrack,ngoodpts + addVertex);
  }  
  //  fTrack = new AliTrackPointArray(ngoodpts);
  // endpepo200709
  //
  //
  for (int i=0; i<npts; i++) idx[i]=i;
  // sort track if required
  TMath::Sort(npts,atp->GetY(),idx); // sort descending...
  //
  Int_t npto=0;
  if (fClusLoc.GetSize()<3*npts)    fClusLoc.Set(3*npts);
  if (fClusGlo.GetSize()<3*npts)    fClusGlo.Set(3*npts);
  if (fClusSigLoc.GetSize()<3*npts) fClusSigLoc.Set(3*npts);
  //
  for (int i=0; i<npts; i++) {
    // skip not defined points
    if (intidx[idx[i]]<0) continue;
    atp->GetPoint(p,idx[i]);
    int sid = AliITSAlignMille2Module::GetIndexFromVolumeID(p.GetVolumeID());
    //
    // prealign point if required
    // get matrix used to produce the digits
    AliITSAlignMille2Module *mod = GetMilleModule(intidx[idx[i]]);
    TGeoHMatrix *svOrigMatrix = GetSensorOrigMatrixSID(sid); //mod->GetSensitiveVolumeOrigGlobalMatrix(p.GetVolumeID());
    // get back real local coordinate
    fMeasLoc  = fClusLoc.GetArray() + npto*3;
    fMeasGlo  = fClusGlo.GetArray() + npto*3;
    fSigmaLoc = fClusSigLoc.GetArray() + npto*3;
    fMeasGlo[0]=p.GetX();
    fMeasGlo[1]=p.GetY();
    fMeasGlo[2]=p.GetZ();
    AliDebug(3,Form("Global coordinates of measured point : X=%+f  Y=%+f  Z=%+f \n",fMeasGlo[0],fMeasGlo[1],fMeasGlo[2]));
    svOrigMatrix->MasterToLocal(fMeasGlo,fMeasLoc);
    AliDebug(3,Form("Local coordinates of measured point : X=%+f  Y=%+f  Z=%+f \n",fMeasLoc[0],fMeasLoc[1],fMeasLoc[2]));
    //
    if (p.GetDriftTime()>0) ProcessSDDPointInfo(&p,sid, npto);     // for SDD points extract vdrift
    //
    // update covariance matrix
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
    //
    if (AliLog::GetGlobalDebugLevel()>=2) {
      AliInfo("Original Global Cov Matrix");
      printf("%+.4e %+.4e %+.4e\n%+.4e %+.4e\n%+.4e\n",hcovel[0],hcovel[1],hcovel[2],hcovel[4],hcovel[5],hcovel[8]);
    } 
    //
    // now rotate in local system
    hcov.Multiply(svOrigMatrix);
    hcov.MultiplyLeft(&svOrigMatrix->Inverse());
    // now hcov is LOCAL COVARIANCE MATRIX
    // apply sigma scaling
    Double_t *hcovscl = hcov.GetRotationMatrix();
    if (AliLog::GetGlobalDebugLevel()>=2) {
      AliInfo("Original Local Cov Matrix");
      printf("%+.4e %+.4e %+.4e\n%+.4e %+.4e\n%+.4e\n",hcovscl[0],hcovscl[1],hcovscl[2],hcovscl[4],hcovscl[5],hcovscl[8]);
    } 
    hcovscl[4] = fUseLocalYErr ? kSensSigY2[lrID[idx[i]]] : 1E-8; // error due to the sensor thickness
    //
    for (int ir=3;ir--;) for (int ic=3;ic--;) {
	if (ir==ic) {	  
	  if ( IsZero(hcovscl[ir*3+ic],1e-8) ) hcovscl[ir*3+ic] = 1E-8;
	  else hcovscl[ir*3+ic] *= mod->GetSigmaFactor(ir)*mod->GetSigmaFactor(ic); //RRR
	  fSigmaLoc[ir] = TMath::Sqrt(hcovscl[ir*3+ic]);
	}
	else hcovscl[ir*3+ic]  = 0;
      }
    //
    if (AliLog::GetGlobalDebugLevel()>=2) {
      AliInfo("Modified Local Cov Matrix");
      printf("%+.4e %+.4e %+.4e\n%+.4e %+.4e\n%+.4e\n",hcovscl[0],hcovscl[1],hcovscl[2],hcovscl[4],hcovscl[5],hcovscl[8]);
    } 
    //
    if (fBug==1) {
      // correzione bug LAYER 5  SSD temporanea..
      int ssdidx=AliITSAlignMille2Module::GetIndexFromVolumeID(p.GetVolumeID());
      if (ssdidx>=500 && ssdidx<1248) {
	int ladder=(ssdidx-500)%22;
	if (ladder==18) p.SetVolumeID(AliITSAlignMille2Module::GetVolumeIDFromIndex(ssdidx+1));
	if (ladder==19) p.SetVolumeID(AliITSAlignMille2Module::GetVolumeIDFromIndex(ssdidx-1));
      }
    }
    /// get (evenctually prealigned) matrix of sens. vol.
    TGeoHMatrix *svMatrix = GetSensorCurrMatrixSID(sid);    //mod->GetSensitiveVolumeMatrix(p.GetVolumeID());
    // modify global coordinates according with pre-aligment
    svMatrix->LocalToMaster(fMeasLoc,fMeasGlo);
    // now rotate in local system
    hcov.Multiply(&svMatrix->Inverse());
    hcov.MultiplyLeft(svMatrix);         // hcov is back in GLOBAL RF
    // cure once more
    for (int ir=3;ir--;) for (int ic=3;ic--;) if (IsZero(hcovscl[ir*3+ic])) hcovscl[ir*3+ic] = 0.;
    //    printf("\nErrMatGlob: after\n"); hcov.Print(""); //RRR
    //
    if (AliLog::GetGlobalDebugLevel()>=2) {
      AliInfo("Modified Global Cov Matrix");
      printf("%+.4e %+.4e %+.4e\n%+.4e %+.4e\n%+.4e\n",hcovscl[0],hcovscl[1],hcovscl[2],hcovscl[4],hcovscl[5],hcovscl[8]);
    } 
    //
    Float_t pcov[6];
    pcov[0]=hcovscl[0];
    pcov[1]=hcovscl[1];
    pcov[2]=hcovscl[2];
    pcov[3]=hcovscl[4];
    pcov[4]=hcovscl[5];
    pcov[5]=hcovscl[8];
    //
    p.SetXYZ(fMeasGlo[0],fMeasGlo[1],fMeasGlo[2],pcov);
    //    printf("New Gl coordinates of measured point : X=%f  Y=%f  Z=%f \n",fMeasGlo[0],fMeasGlo[1],fMeasGlo[2]);
    AliDebug(3,Form("New global coordinates of measured point : X=%+f  Y=%+f  Z=%+f \n",fMeasGlo[0],fMeasGlo[1],fMeasGlo[2]));
    fTrack->AddPoint(npto,&p);
    AliDebug(2,Form("Adding point[%d] = ( %+f , %+f , %+f )     volid = %d",npto,fTrack->GetX()[npto],
		    fTrack->GetY()[npto],fTrack->GetZ()[npto],fTrack->GetVolumeID()[npto] ));
    //    printf("Adding %d %d %f\n",npto, p.GetVolumeID(), p.GetY()); 
    npto++;
  }
  //
  fDiamondPointID = -1;
  if (addVertex) {
    fTrack->AddPoint(npto,&fDiamond);
    fMeasLoc  = fClusLoc.GetArray() + npto*3;
    fMeasGlo  = fClusGlo.GetArray() + npto*3;
    fSigmaLoc = fClusSigLoc.GetArray() + npto*3;
    fMeasLoc[0] = fMeasGlo[0] = fDiamond.GetX();
    fMeasLoc[1] = fMeasGlo[1] = fDiamond.GetY();
    fMeasLoc[2] = fMeasGlo[2] = fDiamond.GetZ();
    fSigmaLoc[0] = fDiamond.GetCov()[0];
    fSigmaLoc[1] = fDiamond.GetCov()[3];
    fSigmaLoc[2] = fDiamond.GetCov()[5];
    fDiamondPointID = npto++;
  }
  //
  return fTrack;
}

//________________________________________________________________________________________________________
AliTrackPointArray *AliITSAlignMille2::SortTrack(const AliTrackPointArray *atp) 
{
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
    AliDebug(2,Form("Point[%d] = ( %+f , %+f , %+f )     volid = %d",i,atps->GetX()[i],atps->GetY()[i],atps->GetZ()[i],atps->GetVolumeID()[i] ));
  }
  return atps;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::GetCurrentLayer() const 
{
  // get current layer id
  if (!fGeoManager) {
    AliInfo("ITS geometry not initialized!");
    return -1;
  }
  return (Int_t)AliGeomManager::VolUIDToLayer(fCluster.GetVolumeID());
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::InitModuleParams() 
{
  /// initialize geometry parameters for a given detector
  /// for current cluster (fCluster)
  /// fGlobalInitParam[] is set as:
  ///    [tx,ty,tz,psi,theta,phi]
  ///    (old was [tx,ty,tz,theta,psi,phi] ROOT's angles...)
  /// *** At the moment: using Raffalele's angles definition ***
  ///
  /// return 0 if success
  /// If module is found but has no parameters to vary, return 1

  if (!fGeoManager) {
    AliInfo("ITS geometry not initialized!");
    return -1;
  }

  // now 'voluid' is the volumeID of a SENSITIVE VOLUME (coming from a cluster)

  // set the internal index (index in module list)
  UShort_t voluid=fCluster.GetVolumeID();
  fCurrentSensID = AliITSAlignMille2Module::GetIndexFromVolumeID(voluid);
  //
  if (fCurrentSensID==-1) { // this is a special "vertex" module
    fCurrentModule = GetMilleModuleByVID(voluid);
    fCurrentSensID = fCurrentModule->GetIndex();

  }
  else {
    //
    // IT IS VERY IMPORTANT: start from the end of the list, where the childs are located !!!
    Int_t k=fNModules-1;
    fCurrentModule = 0;
    // VERY IMPORTANT: if the sensors were explicitly provided, don't look in the supermodules  
    while (k>=0 && ! (fCurrentModule=GetMilleModule(k))->IsIn(voluid)) k--;
    if (k<0) return -3;
  }
  //
  for (int i=AliITSAlignMille2Module::kMaxParTot;i--;) fModuleInitParam[i] = 0.0;
  //
  int clID = fCluster.GetUniqueID()-1;
  if (clID<0) { // external cluster
    fMeasGlo  = &fExtClusterPar[0];
    fMeasLoc  = &fExtClusterPar[3];
    fSigmaLoc = &fExtClusterPar[6];
    fExtClusterPar[0] = fCluster.GetX();
    fExtClusterPar[1] = fCluster.GetY();
    fExtClusterPar[2] = fCluster.GetZ();
    //
    TGeoHMatrix *svMatrix = fCurrentModule->GetSensitiveVolumeMatrix(voluid);
    svMatrix->MasterToLocal(fMeasGlo,fMeasLoc);  
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
    if (fSigmaLoc[0]<0.0010) fSigmaLoc[0]=0.0010;
    if (fSigmaLoc[2]<0.0010) fSigmaLoc[2]=0.0010;
    //
  }
  else {
    int offs = 3*clID;
    fMeasGlo  = fClusGlo.GetArray()  + offs;
    fMeasLoc  = fClusLoc.GetArray()  + offs;
    fSigmaLoc = fClusSigLoc.GetArray() + offs;
  }
  //
  // set minimum value for SigmaLoc to 10 micron 
  if (fSigmaLoc[0]<0.0010) fSigmaLoc[0]=0.0010;
  if (fSigmaLoc[2]<0.0010) fSigmaLoc[2]=0.0010;
  //
  AliDebug(2,Form("Local coordinates of measured point : X=%+f  Y=%+f  Z=%+f \n",fMeasLoc[0] ,fMeasLoc[1] ,fMeasLoc[2] ));
  AliDebug(2,Form("Setting StDev from CovMat : fSigmaLocX=%g  fSigmaLocY=%g fSigmaLocZ=%g \n",fSigmaLoc[0] ,fSigmaLoc[1] ,fSigmaLoc[2] ));
  //   
  return 0;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::Print(Option_t*) const 
{
  // print current status 
  printf("*** AliMillepede for ITS ***\n");
  printf("    Number of defined super modules: %d\n",fNModules);
  printf("    Obtained parameters refer to %s Deltas\n",fUseGlobalDelta ? "GLOBAL":"LOCAL");
  //
  if (fGeoManager)
    printf("    geometry loaded from %s\n",fGeometryPath.Data());
  else
    printf("    geometry not loaded\n");
  //  
  if (fUsePreAlignment) 
    printf("    using prealignment from %s \n",fPreDeltaPath.Data());
  else
    printf("    prealignment not used\n");    
  //
  //
  if (fBOn) 
    printf("    B Field set to %+f T - using helices\n",fBField);
  else
    printf("    B Field OFF - using straight lines \n");
  //
  if (fTPAFitter)
    printf("    Using AliITSTPArrayFit class for track fitting\n");
  else 
    printf("    Using StraightLine/Riemann fitter for track fitting\n");
  //
  printf("Using local Y error due to the sensor thickness: %s\n",(fUseLocalYErr && fTPAFitter) ? "ON":"OFF");
  //
  for (int itp=0;itp<kNDataType;itp++) {
    if (fRequirePoints[itp]) printf("    Required points in %s tracks:\n",itp==kCosmics? "cosmics" : "collisions");
    for (Int_t i=0; i<6; i++) {
      if (fNReqLayUp[itp][i]>0) printf("        Layer %d : %d points with Y>0\n",i+1,fNReqLayUp[itp][i]);
      if (fNReqLayDown[itp][i]>0) printf("        Layer %d : %d points with Y<0\n",i+1,fNReqLayDown[itp][i]);
      if (fNReqLay[itp][i]>0) printf("        Layer %d : %d points \n",i+1,fNReqLay[itp][i]);
    }
    for (Int_t i=0; i<3; i++) {
      if (fNReqDetUp[itp][i]>0) printf("        Detector %d : %d points with Y>0\n",i+1,fNReqDetUp[itp][i]);
      if (fNReqDetDown[itp][i]>0) printf("        Detector %d : %d points with Y<0\n",i+1,fNReqDetDown[itp][i]);
      if (fNReqDet[itp][i]>0) printf("        Detector %d : %d points \n",i+1,fNReqDet[itp][i]);
    }
  }
  printf("        SDD VDrift correction         : %s",fIsSDDVDriftMult ? "Mult":"Add");
  printf("        Weight acc. to pT in power    : %f",fWeightPt);
  //  
  printf("\n    Millepede configuration parameters:\n");
  printf("        init factor for chi2 cut      : %.4f\n",fStartFac);
  printf("        final factor for chi2 cut     : %.4f\n",fFinalFac);
  printf("        first iteration cut value     : %.4f\n",fResCutInitial);
  printf("        other iterations cut value    : %.4f\n",fResCut);
  printf("        number of stddev for chi2 cut : %d\n",fNStdDev);
  printf("        def.scaling for local sigmas  : %.4f %.4f %.4f\n",fSigmaFactor[0],fSigmaFactor[1],fSigmaFactor[2]);
  printf("        min.tracks per module         : %d\n",fMinPntPerSens);
  //
  printf("List of defined modules:\n");
  printf("  intidx\tindex\tvoluid\tname\n");
  for (int i=0; i<fNModules; i++) {
    AliITSAlignMille2Module* md = GetMilleModule(i); 
    printf("  %d\t%d\t%d\t%s\n",i,md->GetIndex(),md->GetVolumeID(),md->GetName());
  }
}

//________________________________________________________________________________________________________
AliITSAlignMille2Module  *AliITSAlignMille2::GetMilleModuleByVID(UShort_t voluid) const
{
  // return pointer to a defined supermodule
  // return NULL if error
  Int_t i=IsVIDDefined(voluid);
  if (i<0) return NULL;
  return GetMilleModule(i);
}

//________________________________________________________________________________________________________
AliITSAlignMille2Module  *AliITSAlignMille2::GetMilleModuleBySymName(const Char_t* symname) const
{
  // return pointer to a defined supermodule
  // return NULL if error
  UShort_t vid = AliITSAlignMille2Module::GetVolumeIDFromSymname(symname);
  if (vid>0) return GetMilleModuleByVID(vid);
  else {    // this is not alignable module, need to look within defined supermodules
    int i = IsSymDefined(symname);
    if (i>=0) return  GetMilleModule(i);
  }
  return 0;
}

//________________________________________________________________________________________________________
AliITSAlignMille2Module  *AliITSAlignMille2::GetMilleModuleIfContained(const Char_t* symname) const
{
  // return pointer to a defined/contained supermodule
  // return NULL otherwise
  int i = IsSymContained(symname);
  return i<0 ? 0 : GetMilleModule(i);
}

//________________________________________________________________________________________________________
AliAlignObjParams* AliITSAlignMille2::GetPrealignedObject(const Char_t* symname) const
{
  // get delta from prealignment for given volume
  if (!fPrealignment) return 0;
  for (int ipre=fPrealignment->GetEntriesFast();ipre--;) { // was the corresponding object prealigned?
    AliAlignObjParams* preob = (AliAlignObjParams*)fPrealignment->At(ipre);
    if (!strcmp(preob->GetSymName(),symname)) return preob;
  }
  return 0;
}

//________________________________________________________________________________________________________
AliAlignObjParams* AliITSAlignMille2::GetConstrRefObject(const Char_t* symname) const
{
  // get delta with respect to which the constraint is declared
  if (!fConstrRef) return 0;
  for (int ipre=fConstrRef->GetEntriesFast();ipre--;) { // was the corresponding object prealigned?
    AliAlignObjParams* preob = (AliAlignObjParams*)fConstrRef->At(ipre);
    if (!strcmp(preob->GetSymName(),symname)) return preob;
  }
  return 0;
}

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2::InitRiemanFit() 
{
  // Initialize Riemann Fitter for current track
  // return kFALSE if error

  if (!fBOn) return kFALSE;

  Int_t npts=0;
  AliTrackPoint ap;
  npts = fTrack->GetNPoints();
  AliDebug(3,Form("Fitting track with %d points",npts));
  if (!fRieman) fRieman = new AliTrackFitterRieman();
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

//________________________________________________________________________________________________________
void trackFit2D(Int_t &, Double_t *, double &chi2, double *par, int flag)
{
  // local function for minuit
  const double kTiny = 1.e-14;
  chi2 = 0;
  static AliTrackPoint pnt;
  static Bool_t fullErr2D;
  //
  if (flag==1) fullErr2D = kFALSE;//kTRUE;
  //  fullErr2D = kTRUE;
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

//________________________________________________________________________________________________________
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
  if (IsZero(det)) det = 1E-16;
  fLocalInitParam[0] = (sX*sYY-sY*sXY)/det;
  fLocalInitParam[2] = (sXY*npts-sY*sX)/det;
  //
  fLocalInitParam[1] = (sZ*sYY-sY*sZY)/det;
  fLocalInitParam[3] = (sZY*npts-sY*sZ)/det;
  // pepo200709
  fLocalInitParam[4] = 0.0;
  // endpepo200709

  AliDebug(2,Form("X = p0gx + ugx*Y : p0gx = %+f    ugx = %+f\n",fLocalInitParam[0],fLocalInitParam[2]));
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
  minuit->ExecuteCommand("MIGRAD",arglist,2);
  //
  for (int i=0;i<4;i++) fLocalInitParam[i] = minuit->GetParameter(i);
  for (int i=0;i<4;i++) for (int j=0;j<4;j++) fLocalInitParEr[i][j] = minuit->GetCovarianceMatrixElement(i,j);
  /*
  double amin,edm,errdef;
  int nvpar,nparx;
  minuit->GetStats(amin,edm,errdef,nvpar,nparx);
  amin /= (2*npts - 4);
  printf("Mchi2: %+e\n",amin);
  */
  //
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::IsSymDefined(const Char_t* symname) const
{
  // checks if supermodule with this symname is defined and return the internal index
  // return -1 if not.
  for (int k=fNModules;k--;) if (!strcmp(symname,GetMilleModule(k)->GetName())) return k;
  return -1; 
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::IsSymContained(const Char_t* symname) const
{
  // checks if module with this symname is defined and return the internal index
  // return -1 if not.
  UShort_t vid = AliITSAlignMille2Module::GetVolumeIDFromSymname(symname);
  if (vid>0) return IsVIDContained(vid);
  // only sensors have real vid, but maybe we have a supermodule with fake vid? 
  // IMPORTANT: always start from the end to start from the sensors
  return IsSymDefined(symname);
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::IsVIDDefined(UShort_t voluid) const
{
  // checks if supermodule 'voluid' is defined and return the internal index
  // return -1 if not.
  for (int k=fNModules;k--;) if (voluid==GetMilleModule(k)->GetVolumeID()) return k;
  return -1; 
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::IsVIDContained(UShort_t voluid) const
{
  // checks if the sensitive module 'voluid' is contained inside a supermodule 
  // and return the internal index of the last identified supermodule
  // return -1 if error
  // IMPORTANT: always start from the end to start from the sensors
  if (AliITSAlignMille2Module::GetIndexFromVolumeID(voluid)<0) return -1;
  for (int k=fNModules;k--;) if (GetMilleModule(k)->IsIn(voluid)) return k;
  return -1; 
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::GetRequestedModID(UShort_t voluid) const
{
  // checks if the sensitive module 'voluid' is contained inside a supermodule 
  // and return the internal index of the last identified supermodule
  // return -1 if error
  // IMPORTANT: always start from the end to start from the sensors
  if (AliITSAlignMille2Module::GetIndexFromVolumeID(voluid)<0) return -1;
  int k;
  for (k=fNModules;k--;) if (GetMilleModule(k)->IsIn(voluid)) break;
  if (k<0) return -1;
  AliITSAlignMille2Module* md = GetMilleModule(k);
  while (md && md->IsNotInConf()) md = md->GetParent();
  return md ? md->GetUniqueID() : -1; 
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::CheckCurrentTrack() 
{
  /// checks if AliTrackPoints belongs to defined modules
  /// return number of good poins
  /// return 0 if not enough points

  Int_t npts = fTrack->GetNPoints();
  Int_t ngoodpts=0;
  // debug points
  for (int j=0; j<npts; j++) if (IsVIDContained(fTrack->GetVolumeID()[j])>=0) ngoodpts++;
  //
  if (ngoodpts<fMinNPtsPerTrack) return 0;

  return ngoodpts;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::ProcessTrack(const AliTrackPointArray *track, Double_t wgh) 
{
  /// Process track; Loop over hits and set local equations
  /// here 'track' is a AliTrackPointArray
  /// return 0 if success;
  //
  if (!fIsMilleInit) Init();
  //
  Int_t npts = track->GetNPoints();
  AliDebug(2,Form("*** Input track with %d points ***",npts));

  // preprocessing of the input track: keep only points in defined volumes,
  // move points if prealignment is set, sort by Yglo if required
  fTrackWeight = wgh;
  fTrack=PrepareTrack(track);
  if (!fTrack) {
    RemoveHelixFitConstraint();
    return -1;
  }
  npts = fTrack->GetNPoints();
  if (npts>kMaxPoints) {
    AliError(Form("Compiled with kMaxPoints=%d, current track has %d points",kMaxPoints,npts));
  }
  AliDebug(2,Form("*** Processing prepared track with %d points ***",npts));
  //
  npts = FitTrack();
  if (npts<0) return npts;
  //
  //  printf("Params: "); for (int i=0;i<fNLocal;i++) printf("%+.2e ",fLocalInitParam[i]); printf("\n");//RRR
  Int_t nloceq=0;
  Int_t ngloeq=0;
  static Mille2Data md[kMaxPoints];
  //
  for (Int_t ipt=0; ipt<npts; ipt++) {
    fTrack->GetPoint(fCluster,ipt);
    fCluster.SetUniqueID(ipt+1);
    AliDebug(2,Form("\n--- processing point %d --- \n",ipt));    

    // set geometry parameters for the the current module
    if (InitModuleParams()) continue;
    AliDebug(2,Form("    VolID=%d  Index=%d  InternalIdx=%d  symname=%s\n", 
		    track->GetVolumeID()[ipt], fCurrentModule->GetIndex(),
		    fCurrentModule->GetUniqueID(), AliGeomManager::SymName(track->GetVolumeID()[ipt]) ));
    AliDebug(2,Form("    Preprocessed Point = ( %+f , %+f , %+f ) \n",fCluster.GetX(),fCluster.GetY(),fCluster.GetZ()));
    int res = fTPAFitter ? AddLocalEquationTPA(md[nloceq]) : AddLocalEquation(md[nloceq]);
    if (res<0) {fTotBadLocEqPoints++; nloceq = 0; break;}
    else if (res==0) nloceq++;
    else {nloceq++; ngloeq++;}
  } // end loop over points
  //
  fTrack=NULL;
  // not enough good points?
  if (nloceq<fMinNPtsPerTrack || ngloeq<1) return -1;
  //
  // finally send local equations to millepede
  SetLocalEquations(md,nloceq);
  fMillepede->SaveRecordData(); // RRR
  //
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::FitTrack() 
{
  // Fit the track with selected constraints
  //
  const Double_t kfDiamondTolerance = 0.1;  //diamond tolerance on top of the MS error
  if (!fTrack) return -1;
  int npts = fTrack->GetNPoints();
  //
  if (fTPAFitter) {  // use dediacted fitter
    //
    // if the diamond point is attached, for the moment don't include it in the fit
    fTPAFitter->AttachPoints(fTrack,0, fDiamondPointID>0 ? fDiamondPointID-1 : npts-1); 
    fTPAFitter->SetBz(fBField);
    fTPAFitter->SetTypeCosmics(IsTypeCosmics());
    if (fIniTrackParamsMeth==1) fTPAFitter->SetIgnoreCov();
    //
    double chi2;
    double chi2f = 0;
    double dca2err;
    double dca2 = 0.;
    Bool_t fitIsDone = kFALSE;
    if (fDiamondPointID>0) { // vertex constraint was added, check if the track looks like prompt
      chi2f = fTPAFitter->Fit(fConstrCharge,fConstrPT,fConstrPTErr);
      if ( chi2f<0 || (chi2f>fNStdDev*fStartFac && fTPAFitter->GetNIterations()>=fTPAFitter->GetMaxIterations()) ) { //RRR
	AliInfo(Form("Track fit failed on checking if it is prompt! skipping this track... Chi2:%+e",chi2f));
	fTPAFitter->Reset();
	//      fTrack = NULL;
	return -5;
      }
      double xyzRes[3];
      fTPAFitter->GetResiduals(xyzRes,&fDiamondI,kTRUE);
      dca2 = xyzRes[0]*xyzRes[0] + xyzRes[1]*xyzRes[1];
      double pT = IsFieldON() ? fTPAFitter->GetPt() : 0.45;
      if (pT<0.1) pT = 0.1;
      dca2err = kfDiamondTolerance + 0.02/pT;
      if (dca2>dca2err*dca2err) { // this is secondary
	int* clst = (int*) fTrack->GetClusterType();
	clst[fDiamondPointID] = -1;;
	fDiamondPointID = -1; 
	fitIsDone = kTRUE;
	npts--;
      }
      else fTPAFitter->SetFirstLast(0,fDiamondPointID); // fit with diamond
    }
    //    fTPAFitter->SetParAxis(1);
    if (!fitIsDone) chi2 = fTPAFitter->Fit(fConstrCharge,fConstrPT,fConstrPTErr);
    //
    RemoveHelixFitConstraint();  // suppress eventual constraints to not affect fit of the next track
    //
    if ( !fitIsDone && (chi2<0 || (chi2>fNStdDev*fStartFac && fTPAFitter->GetNIterations()>=fTPAFitter->GetMaxIterations())) ) { //RRR
      AliInfo(Form("Track fit failed! skipping this track... Chi2:%+e",chi2));
      if (fDiamondPointID>0) AliInfo(Form("VertexFree fit gave Chi2:%+e with residual %+e",chi2f,TMath::Sqrt(dca2)));
      /*
	fTrack->Print("");
	fTPAFitter->FitHelixCrude();
	fTPAFitter->SetFitDone();
	fTPAFitter->Print();
      */
      fTPAFitter->Reset();
      //      fTrack = NULL;
      return -5;
    }
    fNLocal = fTPAFitter->IsFieldON() ? 5:4; // Attention: the fitter might have decided to work in line mode
    npts  = fTPAFitter->GetLast() - fTPAFitter->GetFirst() + 1; // actual number of points
    /*
      double *pr = fTPAFitter->GetParams();
      printf("FtPar: %+.5e  %+.5e  %+.5e  %+.5e | chi2:%.3e\n",pr[2],pr[0],pr[3],pr[1],chi2); // RRR
    */
  }
  else {
    //
    if (!fBOn) { // straight lines  
      // set local starting parameters (to be substituted by ESD track parms)
      // local parms (fLocalInitParam[]) are:
      //      [0] = global x coord. of straight line intersection at y=0 plane
      //      [1] = global z coord. of straight line intersection at y=0 plane
      //      [2] = px/py  
      //      [3] = pz/py
      InitTrackParams(fIniTrackParamsMeth); 
      /*
      double *pr = fLocalInitParam;
      printf("FtPar: %+.5e  %+.5e  %+.5e  %+.5e |\n",pr[0],pr[1],pr[2],pr[3]); // RRR
      */
    } 
    else {
      // local parms (fLocalInitParam[]) are the Riemann Fitter params
      if (!InitRiemanFit()) {
	AliInfo("Riemann fit failed! skipping this track...");
	fTrack=NULL;
	return -5;
      }
    }
  }
  return npts;
  //
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::CalcIntersectionPoint(Double_t *lpar, Double_t *gpar) 
{
  /// calculate track intersection point in local coordinates
  /// according with a given set of parameters (local(4) and global(6))
  /// and fill fPintLoc/Glo
  ///    local are:   pgx0, pgz0, ugx, ugz   OR   riemann fitters pars
  ///    global are:  tx,ty,tz,psi,theta,phi (Raff's delta angles in deg.)
  /// return 0 if success
  
  AliDebug(3,Form("lpar = %g %g %g %g %g\ngpar= %g %g %g %g %g %g\n",lpar[0],lpar[1],lpar[2],lpar[3],lpar[4],gpar[0],gpar[1],gpar[2],gpar[3],gpar[4],gpar[5]));
  AliDebug(3,Form("deltalpar = %g %g %g %g %g\n",lpar[0]-fLocalInitParam[0],lpar[1]-fLocalInitParam[1],lpar[2]-fLocalInitParam[2],lpar[3]-fLocalInitParam[3],lpar[4]-fLocalInitParam[4]));

  
  // prepare the TGeoHMatrix
  TGeoHMatrix *tempHMat = fCurrentModule->GetSensitiveVolumeModifiedMatrix(fCluster.GetVolumeID(),gpar,
									   !fUseGlobalDelta);
  if (!tempHMat) return -1;
  
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

    AliDebug(3,Form("Riemann frame:  fAlpha = %+f  =  %+f  ",alpha,alpha*180./TMath::Pi()));
    AliDebug(3,Form("   prf_glo=( %+f , %+f , %+f )  prf_rf=( %+f , %+f , %+f )\n", x1g,y1g,z1g, x1t,y1t,z1t));
    AliDebug(3,Form("   mov_glo=( %+f , %+f , %+f )      rf=( %+f , %+f , %+f )\n",x2g,y2g,z2g, x2t,y2t,z2t));
        
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
  AliDebug(3,Form("Line vector: ( %+f , %+f , %+f )  point:( %+f , %+f , %+f )\n",v0g[0],v0g[1],v0g[2],p0g[0],p0g[1],p0g[2]));
  
  // same in local coord.
  Double_t p0l[3],v0l[3];
  tempHMat->MasterToLocalVect(v0g,v0l);
  tempHMat->MasterToLocal(p0g,p0l);
  
  if (TMath::Abs(v0l[1])<1e-15) {
    AliInfo("Track Y direction in local frame is zero! Cannot proceed...");
    return -1;
  }
  
  // local intersection point
  fPintLoc[0] = p0l[0] - (v0l[0]/v0l[1])*p0l[1];
  fPintLoc[1] = 0;
  fPintLoc[2] = p0l[2] - (v0l[2]/v0l[1])*p0l[1];
  
  // global intersection point
  tempHMat->LocalToMaster(fPintLoc,fPintGlo);
  AliDebug(3,Form("Intesect. point: L( %+f , %+f , %+f )  G( %+f , %+f , %+f )\n",fPintLoc[0],fPintLoc[1],fPintLoc[2],fPintGlo[0],fPintGlo[1],fPintGlo[2]));
  
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::CalcDerivatives(Int_t paridx, Bool_t islpar) 
{
  /// calculate numerically (ROOT's style) the derivatives for
  /// local X intersection  and local Z intersection
  /// parlist: local  (islpar=kTRUE)  pgx0, pgz0, ugx0, ugz0  OR riemann's params
  ///          global (islpar=kFALSE) tx, ty, tz, psi, theta, phi (Raf's angles in deg)
  /// return 0 if success
  
  // copy initial parameters
  Double_t lpar[kNLocal];
  Double_t gpar[kNParCh];
  Double_t *derivative;
  for (Int_t i=0; i<kNLocal; i++) lpar[i]=fLocalInitParam[i];
  for (Int_t i=0; i<kNParCh; i++) gpar[i]=fModuleInitParam[i];

  // trial with fixed dpar...
  Double_t dpar = 0.;

  if (islpar) { // track parameters
    //dpar=fLocalInitParam[paridx]*0.001;
    // set minimum dpar
    derivative = fDerivativeLoc[paridx];
    if (!fBOn) {
      if (paridx<3) dpar=1.0e-4; // translations
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
    derivative = fDerivativeGlo[paridx];
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
  derivative[0] = h2*(4*d2 - d0)/3.;
  if (TMath::Abs(derivative[0]) < 1.0e-9) derivative[0] = 0.0;

  d0 = pintl4[2]-pintl1[2];
  d2 = 2.*(pintl3[2]-pintl2[2]);
  derivative[2] = h2*(4*d2 - d0)/3.;
  if (TMath::Abs(derivative[2]) < 1.0e-9) derivative[2]=0.0;

  AliDebug(3,Form("\n+++ derivatives +++ \n"));
  AliDebug(3,Form("+++ dXLoc/dpar = %g +++\n",derivative[0]));
  AliDebug(3,Form("+++ dZLoc/dpar = %g +++\n\n",derivative[2]));
  
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::AddLocalEquation(Mille2Data &m) 
{
  /// Define local equation for current cluster in X and Z coor.
  /// and store them to memory
  /// return -1 in case of failure to build some equation
  ///         0 if no free global parameters were found but local eq is built
  ///         1 if both local and global eqs are built
  //
  // store first intersection point
  if (CalcIntersectionPoint(fLocalInitParam, fModuleInitParam)) return -1;  
  for (Int_t i=0; i<3; i++) fPintLoc0[i]=fPintLoc[i];

  AliDebug(2,Form("Intersect. point: L( %+f , %+f , %+f )",fPintLoc[0],fPintLoc[1],fPintLoc[2]));
  
  // calculate local derivatives numerically
  Bool_t zeroX = kTRUE;
  Bool_t zeroZ = kTRUE;
  //
  for (Int_t i=0; i<fNLocal; i++) {
    if (CalcDerivatives(i,kTRUE)) return -1;
    m.fDerLoc[i][kX] = fDerivativeLoc[i][0];
    m.fDerLoc[i][kZ] = fDerivativeLoc[i][2];
    if (zeroX) zeroX = IsZero(fDerivativeLoc[i][0]);
    if (zeroZ) zeroZ = IsZero(fDerivativeLoc[i][2]);
  }
  //  for (Int_t i=0; i<fNLocal; i++) AliDebug(2,Form("Local parameter %d - dXdpar = %g  - dZdpar = %g\n",i,dXdL[i],dZdL[i]));
  //
  if (zeroX) {AliInfo("Skipping: zero local X derivatives!"); return -1;}
  if (zeroZ) {AliInfo("Skipping: zero local Z derivatives!"); return -1;}
  //
  int status = 0;
  int ifill = 0;
  //
  AliITSAlignMille2Module* endModule = fCurrentModule;
  //
  zeroX = zeroZ = kTRUE;
  Bool_t dfDone[kNParCh];
  for (int i=kNParCh;i--;) dfDone[i] = kFALSE;
  m.fNModFilled = 0;
  // 
  // special block for SDD derivatives
  Double_t jacobian[kNParChGeom];
  Int_t nmodTested = 0;
  //
  do {
    if (fCurrentModule->GetNParFree()==0) continue;
    nmodTested++;
    for (Int_t i=0; i<kNParChGeom; i++) {   // common for all sensors: derivatives over geom params 
      //
      if (!fUseGlobalDelta) dfDone[i] = kFALSE; // for global deltas the derivatives at diff. levels are different
      if (fCurrentModule->GetParOffset(i)<0) continue; // this parameter is not explicitly fitted
      if (!dfDone[i]) { 
	if (CalcDerivatives(i,kFALSE)) return -1; 
	else {
	  dfDone[i] = kTRUE;
	  if (zeroX) zeroX = IsZero(fDerivativeGlo[i][0]);
	  if (zeroZ) zeroZ = IsZero(fDerivativeGlo[i][2]);
	}
      }
      //
      m.fDerGlo[ifill][kX] = fDerivativeGlo[i][0];
      m.fDerGlo[ifill][kZ] = fDerivativeGlo[i][2];
      m.fParMilleID[ifill++] = fCurrentModule->GetParOffset(i);
    }
    //
    // specific for special sensors
    Int_t sddLR = -1;
    if ( fCurrentModule->IsSDD() && 
	 (fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFT0)>=0  ||
	  //	  fCurrentModule->GetParOffset(sddLR = fMeasLoc[kX]>0 ?
	  fCurrentModule->GetParOffset(sddLR = GetVDriftSDD()>0 ? 
				       AliITSAlignMille2Module::kDOFDVL : AliITSAlignMille2Module::kDOFDVR)>=0)
	 ) {
      //
      // assume for sensor local xloc = xloc0 + V0*dT0+dV*(T-T0)
      // where V0 and T are the nominal drift velocity, time and time0
      // and the dT0 and dV are the corrections:
      // dX/dT0 = dX/dxloc * dxloc/dT0 = dX/dxloc * V0
      // dX/dV  = dX/dxloc * dxloc/dV =  dX/dxloc * (T-T0)
      // IMPORTANT: the geom derivatives are over the SENSOR LOCAL parameters
      //
      if (!dfDone[AliITSAlignMille2Module::kDOFT0] ||  !dfDone[sddLR]) {
	//
	double dXdxlocsens=0., dZdxlocsens=0.;
	//
	// if the current module is the sensor itself and we work with local params, then 
	// we can directly take dX/dxloc_sens dZ/dxloc_sens
	if (!fUseGlobalDelta && fCurrentModule->GetVolumeID()==fCluster.GetVolumeID()) {
	  if (!dfDone[AliITSAlignMille2Module::kDOFTX]) {
	    CalcDerivatives(AliITSAlignMille2Module::kDOFTX,kFALSE); 
	    dfDone[AliITSAlignMille2Module::kDOFTX] = kTRUE;
	  }
	  dXdxlocsens = fDerivativeGlo[AliITSAlignMille2Module::kDOFTX][0];
	  dZdxlocsens = fDerivativeGlo[AliITSAlignMille2Module::kDOFTX][2];
	}
	//
	else { // need to perform some transformations
	  // fetch the jacobian of the transformation from the sensors local frame to the frame
	  // where the parameters are defined:
	  // Global: dX/dxloc_sens = dX/dxgl*dxgl/dxloc_sens + ...dX/dphigl*dphigl/dxloc_sens
	  if (fUseGlobalDelta) fCurrentModule->CalcDerivGloLoc(fCluster.GetVolumeID(),
							       AliITSAlignMille2Module::kDOFTX, jacobian);
	  // Local:  dX/dxloc_sens = dX/dxcurr*dxcurr/dxloc_sens +..+dX/dphicurr * dphicurr/dxloc_sens 
	  else                 fCurrentModule->CalcDerivCurLoc(fCluster.GetVolumeID(),
							       AliITSAlignMille2Module::kDOFTX, jacobian);
	  //
	  for (int j=0;j<kNParChGeom;j++) {
	    // need global derivative even if the j-th param is locked
	    if (!dfDone[j]) {CalcDerivatives(j,kFALSE); dfDone[j] = kTRUE;}
	    dXdxlocsens += fDerivativeGlo[j][0] * jacobian[j];
	    dZdxlocsens += fDerivativeGlo[j][2] * jacobian[j];
	  }
	}
	//
	if (zeroX) zeroX = IsZero(dXdxlocsens);
	if (zeroZ) zeroZ = IsZero(dZdxlocsens);
	//
	double vdrift = GetVDriftSDD();
	double tdrift = GetTDriftSDD();
	//
	fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][0] = dXdxlocsens*vdrift;
	fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][2] = dZdxlocsens*vdrift;
	dfDone[AliITSAlignMille2Module::kDOFT0] = kTRUE;
	//
	double mltCorr = fIsSDDVDriftMult ? TMath::Abs(vdrift) : 1;
	fDerivativeGlo[sddLR][0] = -dXdxlocsens*mltCorr*TMath::Sign(tdrift,vdrift);
	fDerivativeGlo[sddLR][2] = -dZdxlocsens*mltCorr*TMath::Sign(tdrift,vdrift);
	dfDone[sddLR] = kTRUE;
	//
      }
      //
      if (fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFT0)>=0) {
	m.fDerGlo[ifill][kX] = fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][0];
	m.fDerGlo[ifill][kZ] = fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][2];
	m.fParMilleID[ifill++] = fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFT0);      
      }
      //
      if (fCurrentModule->GetParOffset(sddLR)>=0) {
	m.fDerGlo[ifill][kX] = fDerivativeGlo[sddLR][0];
	m.fDerGlo[ifill][kZ] = fDerivativeGlo[sddLR][2];
	m.fParMilleID[ifill++] = fCurrentModule->GetParOffset(sddLR);      
      }
    }
    //
    m.fModuleID[m.fNModFilled++] = fCurrentModule->GetUniqueID();
  } while( (fCurrentModule=fCurrentModule->GetParent()) );
  //
  if (nmodTested>0 && zeroX) {AliInfo("Skipping: zero global X derivatives!");return -1;}
  if (nmodTested>0 && zeroZ) {AliInfo("Skipping: zero global Z derivatives!");return -1;}
  //
  // ok, can copy to m
  AliDebug(2,Form("Adding local equation X with fMeas=%.6f  and fSigma=%.6f",(fMeasLoc[0]-fPintLoc0[0]), fSigmaLoc[0]));
  m.fMeas[kX] = fMeasLoc[0]-fPintLoc0[0];
  m.fSigma[kX] = fSigmaLoc[0];
  //
  AliDebug(2,Form("Adding local equation Z with fMeas=%.6f  and fSigma=%.6f",(fMeasLoc[2]-fPintLoc0[2]), fSigmaLoc[2]));
  m.fMeas[kZ] = fMeasLoc[2]-fPintLoc0[2];
  m.fSigma[kZ] = fSigmaLoc[2];
  //
  m.fNGlobFilled = ifill;
  fCurrentModule = endModule;
  //
  status += Int_t(!zeroX && !zeroZ); // 0 - only locals, 1 locals + globals
  return status;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::AddLocalEquationTPA(Mille2Data &m) 
{
  /// Define local equation for current cluster in X Y and Z coor.
  /// and store them to memory
  /// return -1 in case of failure to build some equation
  ///         0 if no free global parameters were found but local eq is built
  ///         1 if both local and global eqs are built
  //
  int curpoint = fCluster.GetUniqueID()-1;
  TGeoHMatrix *tempHMat = GetSensorCurrMatrixSID(fCurrentSensID);// fCurrentModule->GetSensitiveVolumeMatrix(fCluster.GetVolumeID());
  //
  fTPAFitter->GetDResDParams(&fDerivativeLoc[0][0], curpoint);    // resid. derivatives over the track parameters 
  for (Int_t i=fNLocal; i--;) tempHMat->MasterToLocalVect(fDerivativeLoc[i],m.fDerLoc[i]); 
  //
  int status = 0;
  // derivatives over the global parameters ---------------------------------------->>>
  Double_t dGL[3];     // derivative of global position vs local X (for SDD)
  Double_t dRdP[3][3]; // derivative of local residuals vs local position
  Double_t dPdG[AliITSAlignMille2Module::kMaxParGeom][3]; // derivatives of local position vs global params
  fTPAFitter->GetDResDPos(&fDerivativeGlo[0][0], curpoint);
  for (int i=3;i--;) tempHMat->MasterToLocalVect(fDerivativeGlo[i],dRdP[i]);
  //
  UInt_t ifill=0, dfDone = 0;
  m.fNModFilled = 0;
  // 
  AliITSAlignMille2Module* endModule = fCurrentModule;
  //
  do {
    if (fCurrentModule->GetNParFree()==0) continue;
    status = 1;
    if (!fUseGlobalDelta) dfDone = 0; // for local deltas the derivatives at diff. levels are different
    Bool_t jacobOK = kFALSE;
    //
    for (Int_t i=0; i<kNParChGeom; i++) {              // common for all sensors: derivatives over geom params
      if (fCurrentModule->GetParOffset(i)<0) continue; // this parameter is not explicitly fitted
      //
      if (!TestWordBit(dfDone,i)) {                    // need to calculate new derivative
	if (!jacobOK) {
	  if (fCurrentSensID!=kVtxSensID) fCurrentModule->CalcDerivDPosDPar(fCluster.GetVolumeID(),fMeasLoc,&dPdG[0][0]); 
	  else for (int ip=AliITSAlignMille2Module::kMaxParGeom;ip--;) for (int jp=3;jp--;) dPdG[ip][jp] = (ip==jp) ? 1:0;	  
	  jacobOK = kTRUE;
	}	
	// dRes_j/dGlo_i = \sum_{k=1:3}  dRes_j/dPos_k * dPos_k/dGlo_i
	fDerivativeGlo[i][kX] = dRdP[kX][kX]*dPdG[i][kX] + dRdP[kY][kX]*dPdG[i][kY] + dRdP[kZ][kX]*dPdG[i][kZ];
	fDerivativeGlo[i][kY] = dRdP[kX][kY]*dPdG[i][kX] + dRdP[kY][kY]*dPdG[i][kY] + dRdP[kZ][kY]*dPdG[i][kZ];
	fDerivativeGlo[i][kZ] = dRdP[kX][kZ]*dPdG[i][kX] + dRdP[kY][kZ]*dPdG[i][kY] + dRdP[kZ][kZ]*dPdG[i][kZ];
	SetWordBit(dfDone,i);
      }
      //
      m.fDerGlo[ifill][kX] = fDerivativeGlo[i][kX];
      m.fDerGlo[ifill][kY] = fDerivativeGlo[i][kY];
      m.fDerGlo[ifill][kZ] = fDerivativeGlo[i][kZ];
      m.fParMilleID[ifill++] = fCurrentModule->GetParOffset(i);
      //
    }
    //
    if ( fCurrentModule->IsSDD() ) {     // specific for SDD
      //
      // assume for sensor local xloc = xloc0 + V0*dT0+dV*(T-T0)
      // where V0 and T are the nominal drift velocity, time and time0
      // and the dT0 and dV are the corrections:
      // drloc_i/dT0 = sum_j drloc_i/dMeasGlo_j * dMeasGlo_j/dT0 = 
      //             = sum_j drloc_i/dMeasGlo_j sum_k dMeasGlo_j/dMeasLoc_k * dMeasLoc_k/dT0
      //             = sum_j drloc_i/dMeasGlo_j dMeasGlo_j/dMeasLoc_X * V0
      //
      // drloc_i/dV0 = sum_j drloc_i/dMeasGlo_j * dMeasGlo_j/dV0 = 
      //             = sum_j drloc_i/dMeasGlo_j sum_k dMeasGlo_j/dMeasLoc_k * dMeasLoc_k/dV0
      //             = sum_j drloc_i/dMeasGlo_j dMeasGlo_j/dMeasLoc_X * T0

      // IMPORTANT: the geom derivatives are over the SENSOR LOCAL parameters
      //
      Bool_t jacOK = kFALSE;
      //Int_t sddLR = fMeasLoc[kX]>0 ? AliITSAlignMille2Module::kDOFDVL : AliITSAlignMille2Module::kDOFDVR;
      Int_t sddLR = GetVDriftSDD()>0 ? AliITSAlignMille2Module::kDOFDVL : AliITSAlignMille2Module::kDOFDVR;
      if (fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFT0)>=0) {
	if (!TestWordBit(dfDone, AliITSAlignMille2Module::kDOFT0)) {
	  double vdrift = GetVDriftSDD();
	  JacobianPosGloLoc(kX,dGL);
	  jacOK = kTRUE;
	  fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][kX] = 
	    vdrift*(dRdP[kX][kX]*dGL[kX] + dRdP[kY][kX]*dGL[kY] + dRdP[kZ][kX]*dGL[kZ]);
	  fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][kY] = 
	    vdrift*(dRdP[kX][kY]*dGL[kX] + dRdP[kY][kY]*dGL[kY] + dRdP[kZ][kY]*dGL[kZ]);
	  fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][kZ] = 
	    vdrift*(dRdP[kX][kZ]*dGL[kX] + dRdP[kY][kZ]*dGL[kY] + dRdP[kZ][kZ]*dGL[kZ]);
	  //
	  SetWordBit(dfDone, AliITSAlignMille2Module::kDOFT0);
	}
	m.fDerGlo[ifill][kX] = fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][kX];
	m.fDerGlo[ifill][kY] = fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][kY];
	m.fDerGlo[ifill][kZ] = fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][kZ];
	m.fParMilleID[ifill++] = fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFT0);      
      }
      //
      if (fCurrentModule->GetParOffset(sddLR)>=0) {
	if (!TestWordBit(dfDone, sddLR)) {
	  double tdrift = TMath::Sign(GetTDriftSDD(), GetVDriftSDD());
	  double vdrift = fIsSDDVDriftMult ? TMath::Abs(GetVDriftSDD()) : 1;
	  if (!jacOK) JacobianPosGloLoc(kX,dGL);
	  fDerivativeGlo[sddLR][kX] = 
	    -tdrift*vdrift*(dRdP[kX][kX]*dGL[kX] + dRdP[kY][kX]*dGL[kY] + dRdP[kZ][kX]*dGL[kZ]);
	  fDerivativeGlo[sddLR][kY] = 
	    -tdrift*vdrift*(dRdP[kX][kY]*dGL[kX] + dRdP[kY][kY]*dGL[kY] + dRdP[kZ][kY]*dGL[kZ]);
	  fDerivativeGlo[sddLR][kZ] = 
	    -tdrift*vdrift*(dRdP[kX][kZ]*dGL[kX] + dRdP[kY][kZ]*dGL[kY] + dRdP[kZ][kZ]*dGL[kZ]);
	  SetWordBit(dfDone, sddLR);
	}
	m.fDerGlo[ifill][kX] = fDerivativeGlo[sddLR][kX];
	m.fDerGlo[ifill][kY] = fDerivativeGlo[sddLR][kY];
	m.fDerGlo[ifill][kZ] = fDerivativeGlo[sddLR][kZ];
	m.fParMilleID[ifill++] = fCurrentModule->GetParOffset(sddLR);      
      }
    }
    //
    m.fModuleID[m.fNModFilled++] = fCurrentModule->GetUniqueID();
  } while( (fCurrentModule=fCurrentModule->GetParent()) );
  //
  // store first local residuals
  fTPAFitter->GetResiduals(fPintLoc , curpoint);       // lab residuals
  for (int i=3;i--;) fPintLoc[i] = -fPintLoc[i];
  tempHMat->MasterToLocalVect(fPintLoc,m.fMeas);       // local residuals 
  m.fSigma[kX] = fSigmaLoc[kX];
  m.fSigma[kY] = fSigmaLoc[kY];
  m.fSigma[kZ] = fSigmaLoc[kZ];
  //
  m.fNGlobFilled = ifill;
  fCurrentModule = endModule;
  //
  return status;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::SetLocalEquations(const Mille2Data *marr, Int_t neq) 
{
  /// Set local equations with data stored in m
  /// return 0 if success
  //
  for (Int_t j=0; j<neq; j++) {
    //
    const Mille2Data &m = marr[j];
    //
    Bool_t filled = kFALSE;
    for (int ic=3;ic--;) {
      // for the diamond point (if any) the Y residual is accounted
      if (ic==kY && !fUseLocalYErr && !(m.fModuleID[0]==fDiamondModID)) continue;
      AliDebug(2,Form("setting local equation %c with fMeas=%.6f  and fSigma=%.6f",fgkXYZ[ic],m.fMeas[ic], m.fSigma[ic]));      
      Int_t nzero = 0;
      for (int i=fNLocal; i--;) nzero += SetLocalDerivative(i,m.fDerLoc[i][ic] );
      if (nzero==fNLocal) { 
	AliInfo(Form("Skipping %c residual due to the zero derivatives!",fgkXYZ[ic])); 
	continue; 
      }
      for (int i=m.fNGlobFilled;i--;) SetGlobalDerivative( m.fParMilleID[i] , m.fDerGlo[i][ic] );
      fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, m.fMeas[ic], m.fSigma[ic]);  
      filled = kTRUE;
      //
    }
    //
    if (filled) for (int i=m.fNModFilled;i--;) GetMilleModule(m.fModuleID[i])->IncNProcessedPoints();
  }
  //
  double wgh = 1;
  if (GetWeightPt() && fTPAFitter) {
    wgh = fTPAFitter->GetPt();
    if (wgh>10) wgh = 10.;
    if (wgh<0) wgh = fTPAFitter->IsTypeCosmics() ? 7 : 0.5;
    if (GetWeightPt()>0) wgh = TMath::Power(wgh,GetWeightPt());
  }
  fMillepede->SetRecordWeight(wgh*fTrackWeight);
  //
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::GlobalFit()
{
  /// Call global fit; Global parameters are stored in parameters
  if (!fIsMilleInit) Init();
  //
  ApplyPreConstraints();
  int res = fMillepede->GlobalFit();
  AliInfo(Form("%s fitting global parameters!",res ? "Done":"Failed"));
  if (res) {
    // fetch the parameters
    for (int imd=fNModules;imd--;) {
      AliITSAlignMille2Module* mod = GetMilleModule(imd);
      int nprocp = 0;
      for (int ip=mod->GetNParTot();ip--;) {
	int idp = mod->GetParOffset(ip);
	if (idp<0) continue;    // was not in the explicit fit
	mod->SetParVal(ip,fMillepede->GetFinalParam(idp));
	mod->SetParErr(ip,fMillepede->GetFinalError(idp));
	int np = fMillepede->GetProcessedPoints(idp);
	if (TMath::Abs(np)>TMath::Abs(nprocp)) nprocp = np;
      }
      if (!mod->GetNProcessedPoints()) mod->SetNProcessedPoints(nprocp);
    }

  }
  ApplyPostConstraints();
  return res;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::PrintGlobalParameters() 
{
  /// Print global parameters
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  fMillepede->PrintGlobalParameters();
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::LoadSuperModuleFile(const Char_t *sfile)
{ 
  // load definitions of supermodules from a root file
  // return 0 if success
  AliInfo(Form("Loading SuperModule definitions from %s",sfile));
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
  // pepo200709
  Char_t st[2048];
  char symname[250];
  // end pepo200709

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
    //
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
    //
    fSuperModule.AddAtAndExpand(new AliITSAlignMille2Module(smindex,volid,symname,&m,n,volidsv),fNSuperModules);
    //
    fNSuperModules++;
  }
  //
  smf->Close();
  //
  return 0;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainModuleSubUnitsMean(Int_t idm, Double_t val, UInt_t pattern)
{
  // require that sum of modifications for the childs of this module is = val, i.e.
  // the internal corrections moves the module as a whole by fixed value (0 by default).
  // pattern is the bit pattern for the parameters to constrain
  //
  if (fIsMilleInit) {
    AliInfo("Millepede has been already initialized: no constrain may be added!");
    return;
  }
  if (!GetMilleModule(idm)->GetNChildren()) return;
  TString nm = "cstrSUMean";
  nm += GetNConstraints();
  AliITSAlignMille2Constraint *cstr = new AliITSAlignMille2Constraint(nm.Data(),AliITSAlignMille2Constraint::kTypeMean,
								      idm,val,pattern);
  cstr->SetConstraintID(GetNConstraints());
  fConstraints.Add(cstr);
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainModuleSubUnitsMedian(Int_t idm, Double_t val, UInt_t pattern)
{
  // require that median of the modifications for the childs of this module is = val, i.e.
  // the internal corrections moves the module as a whole by fixed value (0 by default) 
  // module the outliers.
  // pattern is the bit pattern for the parameters to constrain
  // The difference between the mean and the median will be transfered to the parent
  if (fIsMilleInit) {
    AliInfo("Millepede has been already initialized: no constrain may be added!");
    return;
  }
  if (!GetMilleModule(idm)->GetNChildren()) return;
  TString nm = "cstrSUMed";
  nm += GetNConstraints();
  AliITSAlignMille2Constraint *cstr = new AliITSAlignMille2Constraint(nm.Data(),AliITSAlignMille2Constraint::kTypeMedian,
								      idm,val,pattern);
  cstr->SetConstraintID(GetNConstraints());
  fConstraints.Add(cstr);
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainOrphansMean(Double_t val, UInt_t pattern)
{
  // require that median of the modifications for the supermodules which have no parents is = val, i.e.
  // the corrections moves the whole setup by fixed value (0 by default) modulo the outliers.
  // pattern is the bit pattern for the parameters to constrain
  //
  if (fIsMilleInit) {
    AliInfo("Millepede has been already initialized: no constrain may be added!");
    return;
  }
  TString nm = "cstrOMean";
  nm += GetNConstraints();
  AliITSAlignMille2Constraint *cstr = new AliITSAlignMille2Constraint(nm.Data(),AliITSAlignMille2Constraint::kTypeMean,
								      -1,val,pattern);
  cstr->SetConstraintID(GetNConstraints());
  fConstraints.Add(cstr);
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainOrphansMedian(Double_t val, UInt_t pattern)
{
  // require that median of the modifications for the supermodules which have no parents is = val, i.e.
  // the corrections moves the whole setup by fixed value (0 by default) modulo the outliers.
  // pattern is the bit pattern for the parameters to constrain
  //
  if (fIsMilleInit) {
    AliInfo("Millepede has been already initialized: no constrain may be added!");
    return;
  }
  TString nm = "cstrOMed";
  nm += GetNConstraints();
  AliITSAlignMille2Constraint *cstr = new AliITSAlignMille2Constraint(nm.Data(),AliITSAlignMille2Constraint::kTypeMedian,
								      -1,val,pattern);
  cstr->SetConstraintID(GetNConstraints());
  fConstraints.Add(cstr);
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainLocal(const Char_t* name,Double_t *parcf,Int_t npar,Double_t val,Double_t err)
{
  // apply constraint on parameters in the local frame
  if (fIsMilleInit) {
    AliInfo("Millepede has been already initialized: no constrain may be added!");
    return;
  }
  AliITSAlignMille2ConstrArray *cstr = new AliITSAlignMille2ConstrArray(name,parcf,npar,val,err);
  cstr->SetConstraintID(GetNConstraints());
  fConstraints.Add(cstr);
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ApplyGaussianConstraint(const AliITSAlignMille2ConstrArray* cstr)
{
  // apply the constraint on the local corrections of a list of modules
  int nmod = cstr->GetNModules();
  double jacobian[AliITSAlignMille2Module::kMaxParGeom][AliITSAlignMille2Module::kMaxParGeom];
  //
  // check if this not special SDDT0 constraint
  if (cstr->GetPattern()==BIT(AliITSAlignMille2Module::kDOFT0)) {
    for (int i=0;i<cstr->GetNModules()-1;i++) {
      AliITSAlignMille2Module *mdI = GetMilleModule(cstr->GetModuleID(i));
      if (!mdI->IsFreeDOF(AliITSAlignMille2Module::kDOFT0)) continue;
      for (int j=i+1;j<cstr->GetNModules();j++) {
	AliITSAlignMille2Module *mdJ = GetMilleModule(cstr->GetModuleID(j));
	if (!mdJ->IsFreeDOF(AliITSAlignMille2Module::kDOFT0)) continue;
	//
	ResetLocalEquation();
	fGlobalDerivatives[mdI->GetParOffset(AliITSAlignMille2Module::kDOFT0)] = 1;
	fGlobalDerivatives[mdJ->GetParOffset(AliITSAlignMille2Module::kDOFT0)] =-1;
	AddConstraint(fGlobalDerivatives, 0, 1.E-6);
      }
    }
    return;
  }

  for (int imd=nmod;imd--;) {
    int modID = cstr->GetModuleID(imd);
    AliITSAlignMille2Module* mod = GetMilleModule(modID);
    ResetLocalEquation();
    int nadded = 0;
    double value = cstr->GetValue();
    double sigma = cstr->GetError();
    //
    // in case the reference (survey) deltas were imposed for Gaussian constraints
    // already accumulated corrections: they must be subtracted from the constraint value.
    if (IsConstraintWrtRef()) {
      //
      Double_t precal[AliITSAlignMille2Module::kMaxParTot];
      Double_t refcal[AliITSAlignMille2Module::kMaxParTot];
      for (int ip=AliITSAlignMille2Module::kMaxParTot;ip--;) {precal[ip]=0; refcal[ip] = 0.;}
      //
      // check if there was a reference delta provided for this module
      AliAlignObjParams* parref = GetConstrRefObject(mod->GetName());
      if (parref) parref->GetPars(refcal, refcal+3);    // found reference delta
      //
      // extract already applied local corrections for this module
      if (fPrealignment) {
	//
	AliAlignObjParams *preo = GetPrealignedObject(mod->GetName());
	if (preo) {
	  TGeoHMatrix preMat,tmpMat = *mod->GetMatrix(); //  Delta_Glob * Delta_Glob_Par * M
	  preo->GetMatrix(preMat);                       //  Delta_Glob
	  preMat.MultiplyLeft( &tmpMat.Inverse() );      //  M^-1 * Delta_Glob_Par^-1 = (Delta_Glob_Par * M)^-1
	  tmpMat.MultiplyLeft( &preMat );                //  (Delta_Glob_Par * M)^-1 * Delta_Glob * Delta_Glob_Par * M = Delta_loc
	  AliAlignObjParams algob;
	  algob.SetMatrix(tmpMat);
	  algob.GetPars(precal,precal+3); // local corrections for geometry
	}
      }
      //
      // subtract the contribution to constraint from precalibration 
      for (int ipar=cstr->GetNCoeffs();ipar--;) value += (refcal[ipar]-precal[ipar])*cstr->GetCoeff(ipar);
      //
    } 
    //    
    if (fUseGlobalDelta) mod->CalcDerivLocGlo(&jacobian[0][0]);
    //
    for (int ipar=cstr->GetNCoeffs();ipar--;) {
      double coef = cstr->GetCoeff(ipar);
      if (IsZero(coef)) continue;
      //
      if (!fUseGlobalDelta || ipar>= AliITSAlignMille2Module::kMaxParGeom) { // 
	// we are working with local params or if the given param is not related to geometry, 
	// apply the constraint directly
	int parPos = mod->GetParOffset(ipar);
	if (parPos<0) continue; // not in the fit
	fGlobalDerivatives[parPos] += coef;
	nadded++;
      }
      else { // we are working with global params, while the constraint is on local ones -> jacobian
	for (int jpar=AliITSAlignMille2Module::kMaxParGeom;jpar--;) {
	  int parPos = mod->GetParOffset(jpar);
	  if (parPos<0) continue;
	  fGlobalDerivatives[parPos] += coef*jacobian[ipar][jpar];
	  nadded++;
	}
      }      
    }
    if (nadded) AddConstraint(fGlobalDerivatives, value, sigma);
  }
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ApplyPreConstraints()
{
  // apply constriants which cannot be imposed after the fit
  int nconstr = GetNConstraints();
  for (int i=0;i<nconstr;i++) {
    AliITSAlignMille2Constraint* cstr = GetConstraint(i);
    //
    if (cstr->GetType() == AliITSAlignMille2ConstrArray::kTypeGaussian) {
      ApplyGaussianConstraint( (AliITSAlignMille2ConstrArray*)cstr);
      continue;
    } 
    //
    if (cstr->GetType() == AliITSAlignMille2Constraint::kTypeMedian) continue; // post type constraint
    //
    if (!fUseGlobalDelta) continue; // mean/med constraints must be applied to global deltas
    // apply constraint on the mean's before the fit
    int imd = cstr->GetModuleID();
    if (imd>=0) {
      AliITSAlignMille2Module* mod = GetMilleModule(imd);
      UInt_t pattern = 0;
      for (int ipar=mod->GetNParTot();ipar--;) {
	if (!cstr->IncludesParam(ipar)) continue;
	if (mod->GetParOffset(ipar)<0) continue; // parameter is not in the explicit fit -> post constraint
	pattern |= 0x1<<ipar;
	cstr->SetApplied(ipar);
      }
      ConstrainModuleSubUnits(imd,cstr->GetValue(),pattern);
      //
    }
    else if (!PseudoParentsAllowed()) {
      ConstrainOrphans(cstr->GetValue(),(UInt_t)cstr->GetPattern());
      cstr->SetApplied(-1);
    }
  }
  //
  // do we need to tie the SDD left/right VDrift corrections
  for (int imd=0;imd<fNModules;imd++) {
    AliITSAlignMille2Module* mod = GetMilleModule(imd);
    if (mod->IsSDD() && mod->IsVDriftLRSame()) TieSDDVDriftsLR(mod);
  }
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ApplyPostConstraints()
{
  // apply constraints which can be imposed after the fit
  int nconstr = GetNConstraints();
  Bool_t convGlo      = kFALSE;
  // check if there is something to do
  int ntodo = 0;
  for (int i=0;i<nconstr;i++) {
    AliITSAlignMille2Constraint* cstr = GetConstraint(i);
    if (cstr->GetType() == AliITSAlignMille2ConstrArray::kTypeGaussian) continue;
    if (cstr->GetRemainingPattern() == 0) continue;
    ntodo++;
  }
  if (!ntodo) return;
  //
  if (!fUseGlobalDelta) { // need to convert to global params
    ConvertParamsToGlobal();
    convGlo = kTRUE;
  }
  //
  for (int i=0;i<nconstr;i++) {
    AliITSAlignMille2Constraint* cstr = GetConstraint(i);
    if (cstr->GetType() == AliITSAlignMille2ConstrArray::kTypeGaussian) continue;
    //
    int imd = cstr->GetModuleID();
    //
    if (imd>=0) {
      AliITSAlignMille2Module* mod = GetMilleModule(imd);
      UInt_t pattern = 0;
      for (int ipar=mod->GetNParTot();ipar--;) {
	if (cstr->IsApplied(ipar))      continue;
	if (!cstr->IncludesParam(ipar)) continue;
	if (!mod->IsFreeDOF(ipar))      continue; // parameter is fixed, will not apply constraint
	pattern |= 0x1<<ipar;
	cstr->SetApplied(ipar);
      }
      if (pattern) PostConstrainModuleSubUnits(cstr->GetType(),cstr->GetModuleID(),cstr->GetValue(),pattern);
      //
    }
    else if (PseudoParentsAllowed()) {
      UInt_t pattern = (UInt_t)cstr->GetRemainingPattern();
      PostConstrainOrphans(cstr->GetType(),cstr->GetValue(),pattern);
      cstr->SetApplied(-1);
    }
  }
  // if there was a conversion, rewind it
  if (convGlo) ConvertParamsToLocal();
  // 
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainModuleSubUnits(Int_t idm, Double_t val, UInt_t pattern)
{
  // require that sum of modifications for the childs of this module is = val, i.e.
  // the internal corrections moves the module as a whole by fixed value (0 by default).
  // pattern is the bit pattern for the parameters to constrain
  //
  //
  AliITSAlignMille2Module* mod = GetMilleModule(idm);
  //
  for (int ip=0;ip<kNParCh;ip++) {
    if ( !((pattern>>ip)&0x1) /*|| !parent->IsFreeDOF(ip)*/) continue;
    ResetLocalEquation();
    int nadd = 0;
    for (int ich=mod->GetNChildren();ich--;) {
      int idpar = ((AliITSAlignMille2Module*)mod->GetChild(ich))->GetParOffset(ip);
      if (idpar<0) continue;
      fGlobalDerivatives[idpar] = 1.0;
      nadd++;
    }
    //
    if (nadd>0) {
      AddConstraint(fGlobalDerivatives,val);
      AliInfo(Form("Constrained param %d for %d submodules of module #%d: %s",ip,nadd,idm,mod->GetName()));
    }
  }
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainOrphans(Double_t val, UInt_t pattern)
{
  // require that median of the modifications for the supermodules which have no parents is = val, i.e.
  // the corrections moves the whole setup by fixed value (0 by default) modulo the outliers.
  // pattern is the bit pattern for the parameters to constrain
  //
  for (int ip=0;ip<kNParCh;ip++) {
    //
    if ( !((pattern>>ip)&0x1) ) continue;
    ResetLocalEquation();
    int nadd = 0;
    for (int imd=fNModules;imd--;) {
      AliITSAlignMille2Module* mod = GetMilleModule(imd);
      if (mod->GetParent()) continue; // this is not an orphan
      int idpar = mod->GetParOffset(ip);
      if (idpar<0) continue;
      fGlobalDerivatives[idpar] = 1.0;
      nadd++;
    }
    if (nadd>0) {
      AddConstraint(fGlobalDerivatives,val);
      AliInfo(Form("Constrained param %d for %d orphan modules",ip,nadd));
    }
  }
  //
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::PostConstrainModuleSubUnits(Int_t type,Int_t idm, Double_t val, UInt_t pattern)
{
  // require that median or mean of the modifications for the childs of this module is = val, i.e.
  // the internal corrections moves the module as a whole by fixed value (0 by default) 
  // module the outliers.
  // pattern is the bit pattern for the parameters to constrain
  // The difference between the mean and the median will be transfered to the parent
  //
  AliITSAlignMille2Module* parent = GetMilleModule(idm);
  int nc = parent->GetNChildren();
  //
  double *tmpArr = new double[nc]; 
  //
  for (int ip=0;ip<kNParCh;ip++) {
    int npc = 0;
    if ( !((pattern>>ip)&0x1) || !parent->IsFreeDOF(ip)) continue;
    // compute the mean and median of the deltas
    int nfree = 0;
    for (int ich=nc;ich--;) {
      AliITSAlignMille2Module* child = parent->GetChild(ich);
      //      if (!child->IsFreeDOF(ip)) continue; 
      tmpArr[nfree++] = child->GetParVal(ip);
    }
    double median=0,mean=0;
    for (int ic0=0;ic0<nfree;ic0++) {// order the deltas 
      mean += tmpArr[ic0];
      for (int ic1=ic0+1;ic1<nfree;ic1++) 
	if (tmpArr[ic0]>tmpArr[ic1]) {double tv=tmpArr[ic0]; tmpArr[ic0]=tmpArr[ic1]; tmpArr[ic1]=tv;}
    }
    //
    int kmed = nfree/2;
    median = (tmpArr[kmed]+tmpArr[nfree-kmed-1])/2.;
    if (nfree>0) mean /= nfree;
    //
    double shift = val - (type==AliITSAlignMille2Constraint::kTypeMean ? mean : median);
    //
    for (int ich=nc;ich--;) {
      AliITSAlignMille2Module* child = parent->GetChild(ich);
      //    if (!child->IsFreeDOF(ip)) continue; 
      child->SetParVal(ip, child->GetParVal(ip) + shift);
      npc++;
    }
    //
    parent->SetParVal(ip, parent->GetParVal(ip) - shift);
    AliInfo(Form("%s constraint: added %+f shift to param[%d] of %d children of module %d: %s",
		 type==AliITSAlignMille2Constraint::kTypeMean ? "MEAN" : "MEDIAN",shift,
		 ip,npc,idm,parent->GetName()));
  }
  delete[] tmpArr;  
  //
  //
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::PostConstrainOrphans(Int_t type,Double_t val, UInt_t pattern)
{
  // require that median or mean of modifications for the supermodules which have no parents is = val, i.e.
  // the corrections moves the whole setup by fixed value (0 by default).
  // pattern is the bit pattern for the parameters to constrain
  //
  int nc = fNModules;
  //
  int norph = 0;
  for (int ich=nc;ich--;) if (!GetMilleModule(ich)->GetParent()) norph ++;
  if (!norph) return;
  double *tmpArr = new double[norph]; 
  //
  for (int ip=0;ip<kNParCh;ip++) {
    int npc = 0;
    if ( !((pattern>>ip)&0x1)) continue;
    // compute the mean and median of the deltas
    int nfree = 0;
    for (int ich=nc;ich--;) {
      AliITSAlignMille2Module* child = GetMilleModule(ich);
      //      if (child->GetParent() || !child->IsFreeDOF(ip)) continue; 
      if (child->GetParent()) continue; 
      tmpArr[nfree++] = child->GetParVal(ip);
    }
    double median=0,mean=0;
    for (int ic0=0;ic0<nfree;ic0++) {// order the deltas 
      mean += tmpArr[ic0];
      for (int ic1=ic0+1;ic1<nfree;ic1++) 
	if (tmpArr[ic0]>tmpArr[ic1]) {double tv=tmpArr[ic0]; tmpArr[ic0]=tmpArr[ic1]; tmpArr[ic1]=tv;}
    }
    //
    int kmed = nfree/2;
    median = (tmpArr[kmed]+tmpArr[nfree-kmed-1])/2.;
    if (nfree>0) mean /= nfree;
    //
    double shift = val - (type==AliITSAlignMille2Constraint::kTypeMean ? mean : median);
    //
    for (int ich=nc;ich--;) {
      AliITSAlignMille2Module* child = GetMilleModule(ich);
      //      if (child->GetParent() || !child->IsFreeDOF(ip)) continue; 
      if (child->GetParent()) continue; 
      child->SetParVal(ip, child->GetParVal(ip) + shift);
      npc++;
    }
    //
    AliInfo(Form("%s constraint: added %+f shift to param[%d] of %d orphan modules",
		 type==AliITSAlignMille2Constraint::kTypeMean ? "MEAN" : "MEDIAN",shift,
		 ip,npc));
  }
  delete[] tmpArr;  
  //
}

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2::IsParModConstrained(const AliITSAlignMille2Module* mod,Int_t par, Bool_t &meanmed, Bool_t &gaussian) const
{
  // check if par of the module participates in some constraint, and set the flag for their types
  meanmed = gaussian = kFALSE;
  //
  if ( mod->IsParConstrained(par) ) gaussian = kTRUE;     // direct constraint on this param
  //
  for (int icstr=GetNConstraints();icstr--;) {
    AliITSAlignMille2Constraint* cstr = GetConstraint(icstr);
    //
    if (!cstr->IncludesModPar(mod,par)) continue;
    if (cstr->GetType()==AliITSAlignMille2ConstrArray::kTypeGaussian) gaussian = kTRUE;
    else meanmed = kTRUE;
    //
    if (meanmed && gaussian) break; // no sense to check further
  }
  //
  return meanmed||gaussian;
}

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2::IsParModFamilyVaried(const AliITSAlignMille2Module* mod,Int_t par,Int_t depth) const
{
  // check if parameter par is varied for this module or its children up to the level depth
  if (depth<0) return kFALSE;
  if (mod->GetParOffset(par)>=0) return kTRUE;
  for (int icld=mod->GetNChildren();icld--;) {
    AliITSAlignMille2Module* child = mod->GetChild(icld);
    if (IsParModFamilyVaried(child, par, depth-1)) return kTRUE;
  }
  return kFALSE;
  //
}

/*
//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2::IsParFamilyFree(AliITSAlignMille2Module* mod,Int_t par,Int_t depth) const
{
  // check if parameter par is varied and is not subjected to gaussian constraint for the children up to the level depth
  if (depth<0) return kTRUE;
  for (int icld=mod->GetNChildren();icld--;) {
    AliITSAlignMille2Module* child = mod->GetChild(icld);
    //if (child->GetParOffset(par)<0) continue;                  // fixed
    Bool_t cstMM=kFALSE,cstGS=kFALSE;
    // does this child have gaussian constraint ?
    if (!IsParModConstrained(child,par,cstMM,cstGS) || !cstGS ) return kTRUE;
    // check its children
    if (!IsParFamilyFree(child,par,depth-1)) return kTRUE;
  }
  return kFALSE;
  //
}
*/

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2::IsParFamilyFree(const AliITSAlignMille2Module* mod,Int_t par,Int_t depth) const
{
  // check if parameter par is varied and is not subjected to gaussian constraint for the children up to the level depth
  if (depth<0) return kFALSE;
  for (int icld=mod->GetNChildren();icld--;) {
    AliITSAlignMille2Module* child = mod->GetChild(icld);
    //if (child->GetParOffset(par)<0) continue;                  // fixed
    Bool_t cstMM=kFALSE,cstGS=kFALSE;
    // does this child have gaussian constraint ?
    if (!IsParModConstrained(child,par,cstMM,cstGS) || !cstGS ) return kTRUE;
    // check its children
    if (IsParFamilyFree(child,par,depth-1)) return kTRUE;
  }
  return kFALSE;
  //
}

//________________________________________________________________________________________________________
Double_t AliITSAlignMille2::GetTDriftSDD() const 
{
  // obtain drift time corrected for t0
  double t = fCluster.GetDriftTime();
  return t - fDriftTime0[ fCluster.GetUniqueID()-1 ];
}

//________________________________________________________________________________________________________
Double_t AliITSAlignMille2::GetVDriftSDD() const 
{
  // obtain corrected drift speed
  return fDriftSpeed[ fCluster.GetUniqueID()-1 ];
}

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2::FixedOrphans() const
{
  // are there fixed modules with no parent (normally in such a case 
  // the constraints on the orphans should not be applied
  if (!IsConfigured()) {
    AliInfo("Still not configured");
    return kFALSE;
  }
  for (int i=0;i<fNModules;i++) {
    AliITSAlignMille2Module* md = GetMilleModule(i);
    if (md->GetParent()==0 && md->GetNParFree()==0) return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConvertParamsToGlobal()
{
  // convert params in local frame to global one
  double pars[AliITSAlignMille2Module::kMaxParGeom];
  for (int imd=fNModules;imd--;) {
    AliITSAlignMille2Module* mod = GetMilleModule(imd);
    if (mod->GeomParamsGlobal()) continue;
    mod->GetGeomParamsGlo(pars);
    mod->SetParVals(pars,AliITSAlignMille2Module::kMaxParGeom);
    mod->SetGeomParamsGlobal(kTRUE);
  }
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConvertParamsToLocal()
{
  // convert params in global frame to local one
  double pars[AliITSAlignMille2Module::kMaxParGeom];
  for (int imd=fNModules;imd--;) {
    AliITSAlignMille2Module* mod = GetMilleModule(imd);
    if (!mod->GeomParamsGlobal()) continue;
    mod->GetGeomParamsLoc(pars);
    mod->SetParVals(pars,AliITSAlignMille2Module::kMaxParGeom);
    mod->SetGeomParamsGlobal(kFALSE);
  }
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::SetBField(Double_t b)
{
  // set Bz value
  if (IsZero(b,1e-5)) {
    fBField = 0.0;
    fBOn = kFALSE;
    fNLocal = 4;
  }
  else {
    fBField = b;
    fBOn = kTRUE;
    fNLocal = 5; // helices
  }
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::ProcessUserInfo(TList* userInfo)
{
  // extract calibration information used for TrackPointArray creation from run info
  //
  if (!userInfo) { AliInfo("No UserInfo is provided"); return 0;}
  //
  TMap *cdbMap=0;
  TList* cdbList=0;
  TObjString *objStr,*objStr1,*keyStr;
  TString cdbStr;
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetCacheFlag(kFALSE);
  //
  int run = userInfo->GetUniqueID();
  AliInfo(Form("UserInfo corresponds to run#%d",run));
  cdbMap  = (TMap*)userInfo->FindObject("cdbMap");
  const TMap *curMap = man->GetStorageMap();
  if (!cdbMap) {AliInfo("No CDB Map found in UserInfo");}
  else {
    if ((objStr=(TObjString*)cdbMap->GetValue("default"))) { // first set default CDB path    
      if ((objStr1=(TObjString*)curMap->GetValue("default")) && objStr1->GetUniqueID()) {
	AliInfo(Form("OCDB default path from UserInfo: %s is overriden by user setting %s",objStr->GetName(),objStr1->GetName()));
      }
      else {
	cdbStr = objStr->GetString();
	man->UnsetDefaultStorage();
	if (man->GetRaw()) man->SetRaw(kFALSE);
	if (cdbStr.BeginsWith("raw://")) cdbStr = "raw://";
	AliInfo(Form("Default CDB Storage from UserInfo: %s",cdbStr.Data()));
	man->SetDefaultStorage( cdbStr.Data() ); // this may be overriden later by configuration file
      }
    }
    if (man->GetRaw() && run>0) man->SetRun(run);
    //    
    // set specific paths relevant for alignment
    TIter itMap(cdbMap);
    while( (keyStr=(TObjString*)itMap.Next()) ) {
      TString keyS = keyStr->GetString();
      if ( keyS == "default" ) continue;
      //
      TObjString* curPath = (TObjString*)curMap->GetValue(keyStr->GetName());
      if (curPath && curPath->GetUniqueID()) {
	AliInfo(Form("Storage for %s from UserInfo\n is overriden by user setting %s",keyS.Data(),curPath->GetName()));
	continue;
      }
      man->SetSpecificStorage( keyS.Data(), cdbMap->GetValue(keyS)->GetName() );
    }
  }
  //
  cdbList = (TList*)userInfo->FindObject("cdbList");  
  if (!cdbList) {AliInfo("No CDB List found in UserInfo");}
  else {
    // Deltas used for TrackPointArray production
    TIter itList(cdbList);
    ResetBit(kSameInitDeltasBit);
    while( (objStr=(TObjString*)itList.Next()) )
      if (objStr->GetString().Contains("ITS/Align/Data")) {
	TString newpath = objStr->GetString();
	AliInfo(Form("Production Misalignment from UserInfo: %s",newpath.Data()));
	if (newpath !=  fIniDeltaPath) fIniDeltaPath = newpath;
	else {
	  AliInfo("Production Misalignment is the same as already loaded");
	  SetBit(kSameInitDeltasBit);
	}
	break;
      }
    // SDD response (time0 and drift speed correction) used for TrackPointArray production
    itList.Reset();
    ResetBit(kSameInitSDDRespBit);
    while( (objStr=(TObjString*)itList.Next()) )
      if (objStr->GetString().Contains("ITS/Calib/RespSDD")) {
	TString newpath = objStr->GetString();
	AliInfo(Form("Production SDD Response from UserInfo: %s",newpath.Data()));
	if (newpath != fIniSDDRespPath) fIniSDDRespPath = newpath; 
	else {
	  AliInfo("Production SDD Response is the same as already loaded");
	  SetBit(kSameInitSDDRespBit);
	}
	break;
      }
    //
    // SDD vdrift used for TrackPointArray production
    itList.Reset();
    ResetBit(kSameInitSDDVDriftBit);
    while( (objStr=(TObjString*)itList.Next()) )
      if (objStr->GetString().Contains("ITS/Calib/DriftSpeedSDD")){
	TString newpath = objStr->GetString();
	AliInfo(Form("Production SDD VDrift from UserInfo: %s",newpath.Data()));
	if (newpath != fIniSDDVDriftPath) fIniSDDVDriftPath = newpath; 
	else {
	  AliInfo("Production SDD VDrift is the same as already loaded");
	  SetBit(kSameInitSDDVDriftBit);
	}
	break;
      }
    // Diamond constraint    
    itList.Reset();
    ResetBit(kSameDiamondBit);
    while( (objStr=(TObjString*)itList.Next()) )
      if (objStr->GetString().Contains("GRP/Calib/MeanVertexSPD")){
	TString newpath = objStr->GetString();
	AliInfo(Form("Diamond constraint from UserInfo: %s",newpath.Data()));
	if (newpath != fDiamondPath) fDiamondPath = newpath; 
	else {
	  AliInfo("Production Diamond Constraint is the same as already loaded");
	  SetBit(kSameDiamondBit);
	}
	break;
      }
    //
  }  
  //
  TList *bzlst = (TList*)userInfo->FindObject("BzkGauss");
  if (bzlst && bzlst->At(0)) {
    objStr = (TObjString*)bzlst->At(0);
    SetBField( objStr->GetString().Atof() );
    AliInfo(Form("Magentic field from UserInfo: %+.2e",GetBField()));
  }
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::LoadSDDResponse(TString& path, AliITSresponseSDD *&resp)
{
  if (path.IsNull()) return 0;
  AliInfo(Form("Loading SDD response from %s",path.Data()));
  //
  AliCDBEntry *entry = 0;
  delete resp;
  resp = 0;
  while(1) {
    if (path.BeginsWith("path: ")) { // must load from OCDB
      entry = GetCDBEntry(path.Data());
      if (!entry) break;
      resp = (AliITSresponseSDD*) entry->GetObject();
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      //      AliCDBManager::Instance()->UnloadFromCache(cdbId->GetPath()); // don't want cached object, read new copy
      //      delete cdbId;
      //      delete entry;
      break;
    }
    //
    if (gSystem->AccessPathName(path.Data())) break;
    TFile* precf = TFile::Open(path.Data());
    if (precf->FindKey("AliITSresponseSDD")) resp = (AliITSresponseSDD*)precf->Get("AliITSresponseSDD");
    else if (precf->FindKey("AliCDBEntry") && (entry=(AliCDBEntry*)precf->Get("AliCDBEntry"))) {
      resp = (AliITSresponseSDD*) entry->GetObject();
      if (resp && resp->InheritsFrom(AliITSresponseSDD::Class())) entry->SetObject(NULL);
      else resp = 0;
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      delete entry;
    }
    //
    precf->Close();
    delete precf;
    break;
  } 
  //
  if (!resp) {AliError(Form("Failed to load SDD response from %s",path.Data())); return -1;}
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::LoadSDDVDrift(TString& path, TObjArray *&arr)
{
  if (path.IsNull()) return 0;
  AliInfo(Form("Loading SDD VDrift from %s",path.Data()));
  //
  AliCDBEntry *entry = 0;
  delete arr;
  arr = 0;
  while(1) {
    if (path.BeginsWith("path: ")) { // must load from OCDB
      entry = GetCDBEntry(path.Data());
      if (!entry) break;
      arr = (TObjArray*) entry->GetObject();
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      //      AliCDBManager::Instance()->UnloadFromCache(cdbId->GetPath()); // don't want cached object, read new copy
      //      delete cdbId;
      //      delete entry;
      break;
    }
    //
    if (gSystem->AccessPathName(path.Data())) break;
    TFile* precf = TFile::Open(path.Data());
    if (precf->FindKey("TObjArray")) arr = (TObjArray*)precf->Get("TObjArray");
    else if (precf->FindKey("AliCDBEntry") && (entry=(AliCDBEntry*)precf->Get("AliCDBEntry"))) {
      arr = (TObjArray*) entry->GetObject();
      if (arr && arr->InheritsFrom(TObjArray::Class())) entry->SetObject(NULL);
      else arr = 0;
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      delete entry;
    }
    //
    precf->Close();
    delete precf;
    break;
  } 
  //
  if (!arr) {AliError(Form("Failed to load SDD vdrift from %s",path.Data())); return -1;}
  arr->SetOwner(kTRUE);
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::LoadDiamond(TString& path)
{
  if (path.IsNull()) return 0;
  AliInfo(Form("Loading Diamond Constraint from %s",path.Data()));
  //
  AliCDBEntry *entry = 0;
  AliESDVertex *vtx = 0;
  while(1) {
    if (path.BeginsWith("path: ")) { // must load from OCDB
      entry = GetCDBEntry(path.Data());
      if (!entry) break;
      vtx = (AliESDVertex*) entry->GetObject();
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      //      AliCDBManager::Instance()->UnloadFromCache(cdbId->GetPath()); // don't want cached object, read new copy
      //      delete cdbId;
      //      delete entry;
      break;
    }
    //
    if (gSystem->AccessPathName(path.Data())) break;
    TFile* precf = TFile::Open(path.Data());
    if (precf->FindKey("AliESDVertex")) vtx = (AliESDVertex*)precf->Get("AliESDVertex");
    else if (precf->FindKey("AliCDBEntry") && (entry=(AliCDBEntry*)precf->Get("AliCDBEntry"))) {
      vtx = (AliESDVertex*) entry->GetObject();
      if (vtx && vtx->InheritsFrom(AliESDVertex::Class())) entry->SetObject(NULL);
      else vtx = 0;
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      delete entry;
    }
    //
    precf->Close();
    delete precf;
    break;
  } 
  //
  if (!vtx) {AliError(Form("Failed to load Diamond constraint from %s",path.Data())); return -1;}
  //
  double cmat[6];
  float cmatF[6];
  vtx->GetCovMatrix(cmat);
  AliITSAlignMille2Module* diamMod = GetMilleModuleByVID(kVtxSensVID);
  if (diamMod) {
    cmat[0] *= diamMod->GetSigmaXFactor()*diamMod->GetSigmaXFactor();
    cmat[2] *= diamMod->GetSigmaYFactor()*diamMod->GetSigmaYFactor();
    cmat[5] *= diamMod->GetSigmaZFactor()*diamMod->GetSigmaZFactor();
    cmat[1] *= diamMod->GetSigmaXFactor()*diamMod->GetSigmaYFactor();
    cmat[3] *= diamMod->GetSigmaXFactor()*diamMod->GetSigmaZFactor();
    cmat[4] *= diamMod->GetSigmaYFactor()*diamMod->GetSigmaZFactor();
  }
  cmatF[0] = cmat[0]; // xx
  cmatF[1] = cmat[1]; // xy
  cmatF[2] = cmat[3]; // xz
  cmatF[3] = cmat[2]; // yy
  cmatF[4] = cmat[4]; // yz
  cmatF[5] = cmat[5]; // zz
  fDiamond.SetXYZ(vtx->GetX(),vtx->GetY(),vtx->GetZ(), cmatF);
  //
  Double_t t0 = cmatF[3]*cmatF[5] - cmatF[4]*cmatF[4];
  Double_t t1 = cmatF[1]*cmatF[5] - cmatF[2]*cmatF[4];
  Double_t t2 = cmatF[1]*cmatF[4] - cmatF[2]*cmatF[3];
  Double_t det = cmatF[0]*t0 - cmatF[1]*t1 + cmatF[2]*t2;
  float cmatFI[6];
  if (TMath::Abs(det)<1e-36) {
    AliError("Diamond constraint cov.matrix is singular");
    vtx->Print();
    exit(1);
  }
  cmatFI[0] =  t0/det;
  cmatFI[1] = -t1/det;
  cmatFI[2] =  t2/det;
  cmatFI[3] =  (cmatF[0]*cmatF[5] - cmatF[2]*cmatF[2])/det;
  cmatFI[4] =  (cmatF[1]*cmatF[2] - cmatF[0]*cmatF[4])/det;
  cmatFI[5] =  (cmatF[0]*cmatF[3] - cmatF[1]*cmatF[1])/det;
  fDiamondI.SetXYZ(vtx->GetX(),vtx->GetY(),vtx->GetZ(), cmatFI);
  AliInfo("Will use following Diamond Constraint (errors inverted):");
  fDiamondI.Print("");
  delete vtx;
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::LoadDeltas(TString& path, TClonesArray *&arr)
{
  if (path.IsNull()) return 0;
  AliInfo(Form("Loading Alignment Deltas from %s",path.Data()));
  //
  AliCDBEntry *entry = 0;
  delete arr;
  arr = 0;
  while(1) {
    if (path.BeginsWith("path: ")) { // must load from OCDB
      entry = GetCDBEntry(path.Data());
      if (!entry) break;
      arr = (TClonesArray*) entry->GetObject();
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      //      AliCDBManager::Instance()->UnloadFromCache(cdbId->GetPath()); // don't want cached object, read new copy
      //      delete cdbId;
      //      delete entry;
      break;
    }
    //
    if (gSystem->AccessPathName(path.Data())) break;
    TFile* precf = TFile::Open(path.Data());
    if (precf->FindKey("ITSAlignObjs")) arr = (TClonesArray*)precf->Get("ITSAlignObjs");
    else if (precf->FindKey("AliCDBEntry") && (entry=(AliCDBEntry*)precf->Get("AliCDBEntry"))) {
      arr = (TClonesArray*) entry->GetObject();
      if (arr && arr->InheritsFrom(TClonesArray::Class())) entry->SetObject(NULL);
      else arr = 0;
      entry->SetObject(NULL);
      entry->SetOwner(kTRUE);
      delete entry;
    }
    precf->Close();
    delete precf;
    break;
  } 
  //
  if (!arr) {AliError(Form("Failed to load Deltas from %s",path.Data())); return -1;}
  //
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::CacheMatricesCurr()
{
  // build arrays for the fast access to sensor matrices from their sensor ID
  //
  TGeoHMatrix mdel;
  AliInfo("Building sensors current matrices cache");
  //
  fCacheMatrixCurr.Delete();
  for (int idx=0;idx<=kMaxITSSensID;idx++) {
    int volID = AliITSAlignMille2Module::GetVolumeIDFromIndex(idx);
    TGeoHMatrix *mcurr = new TGeoHMatrix();
    AliITSAlignMille2Module::SensVolMatrix(volID, mcurr);
    fCacheMatrixCurr.AddAtAndExpand(mcurr,idx);
    //
  }
  //
  TGeoHMatrix *mcurr = new TGeoHMatrix();
  fCacheMatrixCurr.AddAtAndExpand(mcurr,kVtxSensID); // special unit matrix for diamond constraint
  //
  fCacheMatrixCurr.SetOwner(kTRUE);
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::CacheMatricesOrig()
{
  // build arrays for the fast access to sensor original matrices (used for production)
  //
  TGeoHMatrix mdel;
  AliInfo("Building sensors original matrices cache");
  //
  fCacheMatrixOrig.Delete();
  if (!fIniDeltaPath.IsNull()) {
    TClonesArray* prealSav = fPrealignment;
    fPrealignment = 0;
    if (LoadDeltas(fIniDeltaPath,fPrealignment) || ApplyToGeometry()) 
      { AliError("Failed to load/apply initial deltas used to produce points"); return -1;}
    delete fPrealignment; 
    fPrealignment = prealSav; 
  }
  //
  for (int idx=0;idx<=kMaxITSSensID;idx++) {
    int volID = AliITSAlignMille2Module::GetVolumeIDFromIndex(idx);
    TGeoHMatrix *morig = new TGeoHMatrix();
    AliITSAlignMille2Module::SensVolMatrix(volID,morig);
    fCacheMatrixOrig.AddAtAndExpand(morig,idx);
  }
  //
  //
  TGeoHMatrix *mcurr = new TGeoHMatrix();
  fCacheMatrixOrig.AddAtAndExpand(mcurr,kVtxSensID); // special unit matrix for diamond constraint
  //
  fCacheMatrixOrig.SetOwner(kTRUE);

  fUsePreAlignment = 0; 
  InitGeometry();
  //
  return 0;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::RemoveHelixFitConstraint()
{
  // suppress constraint
  fConstrCharge = 0;
  fConstrPT = fConstrPTErr = -1;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainHelixFitPT(Int_t q,Double_t pt,Double_t pterr)
{
  // constrain q and pT of the helical fit of the track (should be set before process.track)
  //
  fConstrCharge = q==0 ? q:TMath::Sign(1,q);
  fConstrPT = pt;
  fConstrPTErr = pterr;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ConstrainHelixFitCurv(Int_t q,Double_t crv,Double_t crverr)
{
  // constrain charge and curvature of the helical fit of the track (should be set before process.track)
  //
  const double kCQConv = 0.299792458e-3;// R = PT/Bz/fgkCQConv with GeV,kGauss,cm
  
  fConstrCharge = q==0 ? q:TMath::Sign(1,q);
  if (crv<0 || IsZero(crv)) {
    fConstrPT    = -1;
    fConstrPTErr = -1;
  }
  else {
    fConstrPT    = 1./crv*fBField*kCQConv;
    fConstrPTErr = crverr>1e-10 ? fConstrPT/crv*crverr : 0.;
  }
}

//________________________________________________________________________________________________________
TClonesArray* AliITSAlignMille2::CreateDeltas()
{
  // Create \Deltas for every explicitly or implicitly (via non-alignable volumes) varied
  // or prealigned module.
  // If the module has inded J in the hierarchy of alignable volumes (0 - the top, most 
  // coarse level), then its Delta is expressed via MP2 \deltas (in global frame) and 
  // prealignment \DeltaP's as:
  // \Delta_J = Y X Y^-1
  // where X = \delta_J * \DeltaP_J
  // Y = Prod_{K=0,J-1} \delta_K
  // Note that \delta_L accounts not only for its own correction but also of all non-alignable
  // modules in the hierarchy chain from L up to the closest alignable: 
  // while (parent && !parent->IsAlignable()) {
  //   \delta_L->MultiplyLeft( \delta_parent ); 
  //   parent = parent->GetParent();
  // }
  //  
  Bool_t convLoc = kFALSE;
  if (!GetUseGlobalDelta()) {
    ConvertParamsToGlobal();
    convLoc = kTRUE;
  }
  //
  AliAlignObjParams tempAlignObj;
  TGeoHMatrix tempMatX,tempMatY,tempMat1;
  //
  TClonesArray *array = new TClonesArray("AliAlignObjParams",10);
  TClonesArray &alobj = *array;
  int idx = 0;
  //
  TGeoManager* geoManager = AliGeomManager::GetGeometry();  
  int nalgtot = geoManager->GetNAlignable();
  //
  for (int ialg=0;ialg<nalgtot;ialg++) {             // loop over all alignable entries
    //
    const char* algname = geoManager->GetAlignableEntry(ialg)->GetName();
    //
    AliITSAlignMille2Module* md     = GetMilleModuleBySymName(algname); // explicitly varied?
    AliITSAlignMille2Module* parent = md ? md->GetParent(): GetMilleModuleIfContained(algname);
    AliAlignObjParams*       preob  = GetPrealignedObject(algname);  // was it prealigned ?
    //
    if (!preob && !md && (!parent || parent->IsAlignable())) continue; // noting to do 
    //
    // create matrix X (see comment) ------------------------------------------------->>>
    // start from unity matrix
    tempMatX.Clear();
    if (preob) {   // account prealigngment
      preob->GetMatrix(tempMat1);
      tempMatX.MultiplyLeft(&tempMat1);
    }
    //
    if (md) {
      tempAlignObj.SetTranslation( md->GetParVal(0),md->GetParVal(1),md->GetParVal(2));
      tempAlignObj.SetRotation(    md->GetParVal(3),md->GetParVal(4),md->GetParVal(5));
      tempAlignObj.GetMatrix(tempMat1);
      tempMatX.MultiplyLeft(&tempMat1);  // acount correction to varied module
    }
    //
    // the corrections to all non-alignable modules from current on 
    // till first alignable should add up to its matrix
    while (parent && !parent->IsAlignable()) {
      tempAlignObj.SetTranslation( parent->GetParVal(0),parent->GetParVal(1),parent->GetParVal(2));
      tempAlignObj.SetRotation(    parent->GetParVal(3),parent->GetParVal(4),parent->GetParVal(5));
      tempAlignObj.GetMatrix(tempMat1);
      tempMatX.MultiplyLeft(&tempMat1);  // add matrix of non-alignable module
      parent = parent->GetParent();
    } 
    // create matrix X (see comment) ------------------------------------------------<<<
    //
    // create matrix Y (see comment) ------------------------------------------------>>>
    // start from unity matrix
    tempMatY.Clear(); 
    while ( parent ) {
      tempAlignObj.SetTranslation( parent->GetParVal(0),parent->GetParVal(1),parent->GetParVal(2));
      tempAlignObj.SetRotation(    parent->GetParVal(3),parent->GetParVal(4),parent->GetParVal(5));
      tempAlignObj.GetMatrix(tempMat1);
      tempMatY.MultiplyLeft(&tempMat1); 
      parent = parent->GetParent();
    }
    // create matrix Y (see comment) ------------------------------------------------<<<
    //
    tempMatX.MultiplyLeft(&tempMatY);
    tempMatX.Multiply(&tempMatY.Inverse());
    //
    if (tempMatX.IsIdentity()) continue; // do not store dummy matrices
    UShort_t vid = AliITSAlignMille2Module::GetVolumeIDFromSymname(algname);
    new(alobj[idx++]) AliAlignObjParams(algname,vid,tempMatX,kTRUE);
    //
  }
  //
  if (convLoc) ConvertParamsToLocal();
  //
  return array;
  //
}

//_______________________________________________________________________________________
AliITSresponseSDD* AliITSAlignMille2::CreateSDDResponse()
{
  // create object with SDD repsonse (t0 and vdrift corrections) accounting for 
  // eventual precalibration
  //
  // if there was a precalibration provided, copy it to new arrray
  AliITSresponseSDD *precal = GetSDDPrecalResp();
  if (!precal)       precal = GetSDDInitResp();
  Bool_t isPreCalMult = precal&&precal->IsVDCorrMult() ? kTRUE : kFALSE; 
  AliITSresponseSDD *calibSDD = new AliITSresponseSDD();
  calibSDD->SetVDCorrMult(fIsSDDVDriftMult);
  //
  // copy initial values to the new object
  if (precal) {
    calibSDD->SetTimeOffset(precal->GetTimeOffset());
    calibSDD->SetADC2keV(precal->GetADC2keV());
    calibSDD->SetChargevsTime(precal->GetChargevsTime());
    for (int ind=kSDDoffsID;ind<kSDDoffsID+kNSDDmod;ind++) {
      calibSDD->SetModuleTimeZero(ind, precal->GetTimeZero(ind));
      calibSDD->SetDeltaVDrift(ind, precal->GetDeltaVDrift(ind),kFALSE);
      calibSDD->SetDeltaVDrift(ind, precal->GetDeltaVDrift(ind),kTRUE);
      calibSDD->SetADCtokeV(ind,precal->GetADCtokeV(ind));
    }
  }
  else for (int ind=kSDDoffsID;ind<kSDDoffsID+kNSDDmod;ind++) calibSDD->SetModuleTimeZero(ind,0);
  //
  Bool_t save = kFALSE;
  for (int imd=GetNModules();imd--;) {
    AliITSAlignMille2Module* md = GetMilleModule(imd);
    if (!md->IsSDD()) continue;
    if (md->IsFreeDOF(AliITSAlignMille2Module::kDOFT0)  ||
	md->IsFreeDOF(AliITSAlignMille2Module::kDOFDVL) ||
	md->IsFreeDOF(AliITSAlignMille2Module::kDOFDVR)) save = kTRUE;
	//
    for (int is=0;is<md->GetNSensitiveVolumes();is++) {
      int ind  = md->GetSensVolIndex(is);
      float t0  = calibSDD->GetTimeZero(ind)      + md->GetParVal(AliITSAlignMille2Module::kDOFT0);
      double dvL = md->GetParVal(AliITSAlignMille2Module::kDOFDVL);
      double dvR = md->GetParVal(AliITSAlignMille2Module::kDOFDVR);
      if (!calibSDD->IsVDCorrMult()) { // save as additive correction
	dvL *= 1e4;
	dvR *= 1e4;
	//
	double conv = 1;
	if (isPreCalMult) conv = 6.4; // convert multiplicative precal correction to additive
	dvL += calibSDD->GetDeltaVDrift(ind,kFALSE)*conv;
	dvR += calibSDD->GetDeltaVDrift(ind,kTRUE)*conv;
      }
      else { // save as multipicative correction
	double conv = 1;
	if (!isPreCalMult) conv = 1./6.4; // convert additive precal correction to multiplicative
	dvL += calibSDD->GetDeltaVDrift(ind,kFALSE)*conv;
	dvR += calibSDD->GetDeltaVDrift(ind,kTRUE)*conv;
      }
      //
      calibSDD->SetModuleTimeZero(ind, t0);
      calibSDD->SetDeltaVDrift(ind, dvL, kFALSE); // left  side correction
      calibSDD->SetDeltaVDrift(ind, dvR, kTRUE); // right side correction
    }
  }
  //
  if (!save) {
    AliInfo("No free parameters for SDD calibration, nothing to save");
    delete calibSDD;
    calibSDD = 0;
  }
  //
  return calibSDD;  
}

//_______________________________________________________________________________________
Int_t AliITSAlignMille2::ReloadInitCalib(TList *userInfo)
{
  // Use provided UserInfo to
  // load the initial calib parameters (geometry, SDD response...)
  // Can be used if set of data was processed with different calibration
  //
  if (!userInfo) {
    AliInfo("Reloading of the Calibration parameters was called with empty userInfo");
    return 1;
  }
  if (ProcessUserInfo(userInfo)) {
    AliInfo("Error in processing user info");
    userInfo->Print();
    exit(1);
  }
  return ReloadInitCalib();
}

//_______________________________________________________________________________________
Int_t AliITSAlignMille2::ReloadInitCalib()
{
  // Load the initial calib parameters (geometry, SDD response...)
  // Can be used if set of data was processed with different calibration
  //
  // 1st cache original matrices
  if (!TestBit(kSameInitDeltasBit)) { // need to reload geometry
    if (InitGeometry()) {
      AliInfo("Failed to re-load ideal geometry");
      exit(1);
    }
    if (CacheMatricesOrig()) {
      AliInfo("Failed to cache new initial geometry");
      exit(1);
    }
    //
    // then reload the prealignment geometry
    if (LoadDeltas(fPreDeltaPath,fPrealignment)) {
      AliInfo(Form("Failed to reload the prealigned geometry %s",fPreDeltaPath.Data()));
      exit(1);
    }
    //
    if (fPrealignment && ApplyToGeometry()) {
      AliInfo(Form("Failed re-apply prealigned geometry %s",fPreDeltaPath.Data()));
      exit(1);
    }
    //
    // usually no need to re-cache the prealignment geometry, it was not changed
    if (fCacheMatrixCurr.GetEntriesFast() != fCacheMatrixOrig.GetEntriesFast()) {
      //      CacheMatricesCurr();
      AliInfo(Form("Failed to cache the prealigned geometry %s",fPreDeltaPath.Data()));
      exit(1);
    }
  }
  else ResetBit(kSameInitDeltasBit);
  //
  // reload initial SDD response
  if (!TestBit(kSameInitSDDRespBit)) {
    if (LoadSDDResponse(fIniSDDRespPath, fIniRespSDD) ) {
      AliInfo(Form("Failed to load new SDD response %s",fIniSDDRespPath.Data()));
      exit(1);
    }
  }
  else ResetBit(kSameInitSDDRespBit);
  //
  // reload initial SDD vdrift
  if (!TestBit(kSameInitSDDVDriftBit)) {
    if (LoadSDDVDrift(fIniSDDVDriftPath, fIniVDriftSDD) ) {
      AliInfo(Form("Failed to load new SDD VDrift %s",fIniSDDVDriftPath.Data()));
      exit(1);
    }
  }
  else ResetBit(kSameInitSDDRespBit);
  //
  // reload diamond info
  if (!TestBit(kSameDiamondBit)) {
    if (LoadDiamond(fDiamondPath) ) {
      AliInfo(Form("Failed to load new Diamond constraint %s",fDiamondPath.Data()));
      exit(1);
    }
  }
  else ResetBit(kSameInitSDDRespBit);
  //
  

  return 0;
}

//_______________________________________________________________________________________
void AliITSAlignMille2::JacobianPosGloLoc(int locid,double* jacobian)
{
  // calculate the locid row of the jacobian for transformation of the local coordinate to global at current point
  TGeoHMatrix* mat = GetSensorCurrMatrixSID(fCurrentSensID);
  const Double_t dpar = 1e-2;
  double sav = fMeasLoc[locid];
  fMeasLoc[locid] += dpar;
  mat->LocalToMaster(fMeasLoc,jacobian);
  fMeasLoc[locid] = sav; // recover original value
  for (int i=3;i--;) jacobian[i] = (jacobian[i]-fMeasGlo[i])/dpar; // the transformation is linear!!!
}

//_______________________________________________________________________________________
void AliITSAlignMille2::TieSDDVDriftsLR(AliITSAlignMille2Module* mod)
{
  // impose equality of Left/Right sides VDrift correction for SDD
  ResetLocalEquation();
  if ( (mod->IsFreeDOF(AliITSAlignMille2Module::kDOFDVL) + mod->IsFreeDOF(AliITSAlignMille2Module::kDOFDVR))==1) {
    AliError("Left/Right VDrift equality is requested for SDD module with only one side VDrift free");
    mod->Print();
    return;
  }
  SetGlobalDerivative(mod->GetParOffset(AliITSAlignMille2Module::kDOFDVL),  1.);
  SetGlobalDerivative(mod->GetParOffset(AliITSAlignMille2Module::kDOFDVR), -1.);
  AddConstraint(fGlobalDerivatives, 0, 1e-12);
  //
}

//_______________________________________________________________________________________
void AliITSAlignMille2::ProcessSDDPointInfo(const AliTrackPoint* pnt,Int_t sID, Int_t pntID)
{
  // extract the drift information from SDD track point
  //
  fDriftTime0[pntID] = fIniRespSDD ? fIniRespSDD->GetTimeZero(sID) : 0.;
  double tdif = pnt->GetDriftTime() - fDriftTime0[pntID];
  if (tdif<0) tdif = 1;
  //
  // VDrift extraction
  double vdrift = 0;
  Bool_t sddSide = kFALSE;
  int sID0 = 2*(sID-kSDDoffsID);
  double zanode = -999;
  //
  if (fIniVDriftSDD) { // SDD VDrift object is provided, use the vdrift from it
    AliITSDriftSpeedArraySDD* drarr;
    double vdR,vdL,xlR,xlL;
    // sometimes xlocal on right side is negative due to the wrong calibration, need to test both hypothesis
    double xlabs = TMath::Abs(fMeasLoc[kX]); 
    drarr  = (AliITSDriftSpeedArraySDD*)fIniVDriftSDD->At(sID0); // left side, xloc>0
    zanode = fSegmentationSDD->GetAnodeFromLocal(xlabs,fMeasLoc[kZ]);
    vdL    = drarr->GetDriftSpeed(0, zanode);
    if (fIniRespSDD) {
      double corr = fIniRespSDD->GetDeltaVDrift(sID, kFALSE);
      if (fIniRespSDD->IsVDCorrMult()) vdL *= (1+corr);
      else vdL += corr;
    }
    xlL    = (fSegmentationSDD->Dx() - vdL*tdif)*1e-4;
    //
    drarr  = (AliITSDriftSpeedArraySDD*)fIniVDriftSDD->At(sID0+1); // right side, xloc<0
    zanode = fSegmentationSDD->GetAnodeFromLocal(-xlabs,fMeasLoc[kZ]) - 256;
    vdR    = drarr->GetDriftSpeed(0, zanode);
    if (fIniRespSDD) {
      double corr = fIniRespSDD->GetDeltaVDrift(sID, kTRUE);
      if (fIniRespSDD->IsVDCorrMult()) vdR *= (1+corr);
      else vdR += corr;
    }
    xlR    = -(fSegmentationSDD->Dx() - vdR*tdif)*1e-4;
    //
    if (TMath::Abs(xlL-fMeasLoc[kX])<TMath::Abs(xlR-fMeasLoc[kX])) {
      sddSide = 0; // left side
      vdrift  = vdL*1e-4;
    }
    else {         // right side
      sddSide = 1;
      vdrift  = vdR*1e-4;
    }
    //
  }
  else { // try to determine the vdrift from the xloc
    vdrift = (fSegmentationSDD->Dx()*1e-4 - TMath::Abs(fMeasLoc[kX]))/tdif;
    sddSide = fMeasLoc[kX]<0; // 0 = left (xloc>0) ; 1 = right (xloc<1)
  }
  //
  if (fPreVDriftSDD) { // use imposed vdrift as a starting point
    zanode = fSegmentationSDD->GetAnodeFromLocal(0.5-sddSide,fMeasLoc[kZ]);
    if (sddSide) zanode -= 256;
    vdrift = ((AliITSDriftSpeedArraySDD*)fPreVDriftSDD->At(sID0+sddSide))->GetDriftSpeed(0, zanode)*1e-4;
  }
  //
  if (vdrift<0) vdrift = 0;
  // at this point we have vdrift and t0 used to create the original point.
  // see if precalibration was provided
  if (fPreRespSDD) {
    float t0Upd = fPreRespSDD->GetTimeZero(sID);
    double corr = fPreRespSDD->GetDeltaVDrift(sID, sddSide);
    if (fPreRespSDD->IsVDCorrMult()) vdrift *= 1+corr; // right side (xloc<0) may have different correction
    else                             vdrift += corr*1e-4;
    tdif    = pnt->GetDriftTime() - t0Upd;
    // correct Xlocal
    fMeasLoc[0] = fSegmentationSDD->Dx()*1e-4 - vdrift*tdif;
    if (sddSide) fMeasLoc[0] = -fMeasLoc[0];
    fDriftTime0[pntID] =  t0Upd;
  }
  // TEMPORARY CORRECTION (if provided) --------------<<<
  fDriftSpeed[pntID] = sddSide ? -vdrift : vdrift;
  //
  //  printf("#%d: t:%+e x:%+e v:%+e: side:%d\n",pntID,fDriftTime0[pntID],fMeasLoc[0],fDriftSpeed[pntID],sddSide);
}

//_______________________________________________________________________________________
AliITSAlignMille2Module* AliITSAlignMille2::CreateVertexModule()
{
  // creates dummy module for vertex constraint
  TGeoHMatrix mt;
  AliITSAlignMille2Module* mod = new AliITSAlignMille2Module(kVtxSensID,kVtxSensVID,"VTX",&mt,0,0);
  fMilleModule.AddAtAndExpand(mod,fNModules);
  mod->SetGeomParamsGlobal(fUseGlobalDelta);
  fDiamondModID = fNModules;
  mod->SetUniqueID(fNModules++);
  mod->SetNotInConf(kTRUE);
  return mod;
  //
}

//_______________________________________________________________________________________
AliCDBEntry* AliITSAlignMille2::GetCDBEntry(const char* path)
{
  // return object from the OCDB
  AliCDBEntry *entry = 0;
  AliInfo(Form("Loading object %s",path));
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBId* cdbId = AliCDBId::MakeFromString(path);
  if (!cdbId) {
    AliError("Failed to create cdbId");
    return 0;
  }
  //
  AliCDBStorage* stor = man->GetDefaultStorage();
  if (!stor && !man->GetRaw()) man->SetDefaultStorage("raw://");
  if (man->GetRaw()) man->SetRun(cdbId->GetFirstRun());
  if (stor) {
    TString tp = stor->GetType();
    if (tp.Contains("alien",TString::kIgnoreCase) && !gGrid) TGrid::Connect("alien:"); 
  } 
  entry = man->Get( *cdbId );
  man->ClearCache();
  //
  delete cdbId;
  return entry;
  //
}

