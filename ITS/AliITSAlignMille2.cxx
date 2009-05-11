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
//  fixing or constraining detection elements for best results. 
// 
//  author M. Lunardon (thanks to J. Castillo), ruben.shahoyan@cern.ch
//-----------------------------------------------------------------------------

#include <TF1.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TGraph.h>
#include <TGeoMatrix.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TVirtualFitter.h>
#include <TGeoManager.h>

#include "AliITSAlignMille2.h"
#include "AliITSgeomTGeo.h"
#include "AliGeomManager.h"
#include "AliMillePede2.h"
#include "AliTrackPointArray.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"
#include "TSystem.h"  // come si fa?
#include "AliTrackFitterRieman.h"


ClassImp(AliITSAlignMille2)


//========================================================================================================

AliITSAlignMille2* AliITSAlignMille2::fgInstance = 0;  
Int_t              AliITSAlignMille2::fgInstanceID = 0;

//________________________________________________________________________________________________________
AliITSAlignMille2::AliITSAlignMille2(const Char_t *configFilename  ) 
: TObject(),
  fMillepede(0),
  fStartFac(16.), 
  fResCutInitial(100.), 
  fResCut(100.),
  fNGlobal(0),
  fNLocal(4),
  fNStdDev(3),
  fIsMilleInit(kFALSE),
  fAllowPseudoParents(kFALSE),
  //
  fCurrentModule(0),
  fTrack(0),
  fTrackBuff(0),
  fCluster(),
  fGlobalDerivatives(0), 
  //
  fMinNPtsPerTrack(3),
  fInitTrackParamsMeth(1),
  fTotBadLocEqPoints(0),
  fRieman(0),
  //
  fConstraints(0),
  //
  fUseGlobalDelta(kFALSE),
  fRequirePoints(kFALSE),
  fTempExcludedModule(-1),
  //
  fGeometryFileName("geometry.root"),
  fPreAlignmentFileName(""),
  fConstrRefFileName(""),
  fGeoManager(0),
  fIsConfigured(kFALSE),
  fPreAlignQF(0),
//
  fCorrectSDD(0),
  fInitialRecSDD(0),
  fPrealignment(0),
  fConstrRef(0),
  fMilleModule(2),
  fSuperModule(2),
  fNModules(0),
  fNSuperModules(0),
  fUsePreAlignment(kFALSE),
  fBOn(kFALSE),
  fBField(0.0),
  fBug(0),
  fMilleVersion(2)
{
  /// main constructor that takes input from configuration file
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
  if (fMillepede)         delete fMillepede;            fMillepede = 0;
  if (fGlobalDerivatives) delete[] fGlobalDerivatives;  fGlobalDerivatives = 0;
  if (fRieman)            delete fRieman;               fRieman = 0;
  if (fPrealignment)      delete fPrealignment;         fPrealignment = 0;
  if (fConstrRef)         delete fConstrRef;            fConstrRef = 0;
  if (fCorrectSDD)        delete fCorrectSDD;           fCorrectSDD = 0;
  if (fInitialRecSDD)     delete fInitialRecSDD;        fInitialRecSDD = 0;
  fTrackBuff.Delete();
  fConstraints.Delete();
  fMilleModule.Delete();
  fSuperModule.Delete();
  if (--fgInstanceID==0) fgInstance = 0;
}

///////////////////////////////////////////////////////////////////////
TObjArray* AliITSAlignMille2::GetConfigRecord(FILE* stream, TString& recTitle, TString& recOpt, Bool_t rew)
{
  TString record;
  static TObjArray* recElems = 0;
  if (recElems) {delete recElems; recElems = 0;}
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
Int_t AliITSAlignMille2::LoadConfig(const Char_t *cfile)
{  /// return 0 if success
  ///        1 if error in module index or voluid
  //
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
  while(1) { 
    //
    // ============= 1: we read some obligatory records in predefined order ================
    //  
    recTitle = "GEOMETRY_FILE";
    if ( !GetConfigRecord(pfc,recTitle,recOpt,1) || 
	 (fGeometryFileName=recOpt).IsNull()     || 
	 gSystem->AccessPathName(recOpt.Data())  ||
	 InitGeometry()	)
      { AliError("Failed to find/load Geometry"); stopped = kTRUE; break;}
    //
    recTitle = "SUPERMODULE_FILE";
    if ( !GetConfigRecord(pfc,recTitle,recOpt,1) || 
	 recOpt.IsNull()                         || 
	 gSystem->AccessPathName(recOpt.Data())  ||
	 LoadSuperModuleFile(recOpt.Data()))
      { AliError("Failed to find/load SuperModules"); stopped = kTRUE; break;}
    //
    recTitle = "CONSTRAINTS_REFERENCE_FILE";      // LOCAL_CONSTRAINTS are defined wrt these deltas
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) ) {
      if (recOpt.IsNull() || recOpt=="IDEAL") SetConstraintWrtRef( "IDEAL" );
      else if (gSystem->AccessPathName(recOpt.Data()) || SetConstraintWrtRef(recOpt.Data()) )
	{ AliError("Failed to load reference deltas for local constraints"); stopped = kTRUE; break;}
    }
    //	 
    recTitle = "PREALIGNMENT_FILE";
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) )
      if ( (fPreAlignmentFileName=recOpt).IsNull() || 
	   gSystem->AccessPathName(recOpt.Data())   ||
	   ApplyToGeometry()) 
	{ AliError(Form("Failed to load Prealignment file %s",recOpt.Data())); stopped = kTRUE; break;}
    //
    recTitle = "PRECALIBSDD_FILE";
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) ) {
      if ( recOpt.IsNull() || gSystem->AccessPathName(recOpt.Data()) ) {stopped = kTRUE; break;}
      AliInfo(Form("Using %s for SDD precalibration",recOpt.Data()));
      TFile* precfi = TFile::Open(recOpt.Data());
      if (!precfi->IsOpen()) {stopped = kTRUE; break;}
      fCorrectSDD = (AliITSresponseSDD*)precfi->Get("AliITSresponseSDD");
      precfi->Close();
      delete precfi;
      if (!fCorrectSDD) {AliError("Precalibration SDD object is not found"); stopped = kTRUE; break;}
    }
    //
    recTitle = "INITCALBSDD_FILE";
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) ) {
      if ( recOpt.IsNull() || gSystem->AccessPathName(recOpt.Data()) ) {stopped = kTRUE; break;}
      AliInfo(Form("Using %s as SDD calibration used in TrackPoints",recOpt.Data()));
      TFile* precf = TFile::Open(recOpt.Data());
      if (!precf->IsOpen()) {stopped = kTRUE; break;}
      fInitialRecSDD = (AliITSresponseSDD*)precf->Get("AliITSresponseSDD");
      precf->Close();
      delete precf;
      if (!fInitialRecSDD) {AliError("Initial Calibration SDD object is not found"); stopped = kTRUE; break;}
    }
    //
    recTitle = "SET_GLOBAL_DELTAS";
    if ( GetConfigRecord(pfc,recTitle,recOpt,1) ) SetUseGlobalDelta(kTRUE);
    //
    // =========== 2: see if there are local gaussian constraints defined =====================
    //            Note that they should be loaded before the modules declaration
    //
    recTitle = "CONSTRAINT_LOCAL";
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
    while( (recArr=GetConfigRecord(pfc,recTitle="",recOpt,0)) ) {
      if (!(recTitle=="MODULE_VOLUID" || recTitle=="MODULE_INDEX")) continue;
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
	for (int j=0; j<fNSuperModules; j++) {
	  if (voluid==GetSuperModule(j)->GetVolumeID()) {
	    mod = new AliITSAlignMille2Module(*GetSuperModule(j));
	    // the matrix might be updated in case some prealignment was applied, check 
	    TGeoHMatrix* mup = AliGeomManager::GetMatrix(mod->GetName());
	    if (mup) *(mod->GetMatrix()) = *mup;
	    fMilleModule.AddAtAndExpand(mod,fNModules);
	    break;
	  }	
	}
      }
      else if (idx<=kMaxITSSensVID) {
	mod = new AliITSAlignMille2Module(voluid);
	fMilleModule.AddAtAndExpand(mod,fNModules);
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
      mod->SetGeomParamsGlobal(fUseGlobalDelta);
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
	vl = 0;
	if (nrecElems>12) {
	  recExt = recArr->At(12)->GetName();
	  if (recExt.IsFloat()) vl = recExt.Atof();
	  else {stopped = kTRUE; break;}
	  irec = 12;
	}
	mod->SetFreeDOF(AliITSAlignMille2Module::kDOFDV,vl);
      }
      //
      mod->SetUniqueID(fNModules);
      mod->EvaluateDOF();
      fNModules++;
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
      if      (recTitle == "SET_PSEUDO_PARENTS")  SetAllowPseudoParents(kTRUE);
      //
      // some optional parameters ----------------------------------------------------------------
      else if (recTitle == "SET_TRACK_FIT_METHOD") {
	if (recOpt.IsNull() || !recOpt.IsDigit() ) {stopped = kTRUE; break;}
	SetInitTrackParamsMeth(recOpt.Atoi());
      }
      //
      else if (recTitle == "SET_MINPNT_TRA") {
	if (recOpt.IsNull() || !recOpt.IsDigit() ) {stopped = kTRUE; break;}
	fMinNPtsPerTrack = recOpt.Atoi();
      }
      //
      else if (recTitle == "SET_NSTDDEV") {
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fNStdDev = (Int_t)recOpt.Atof();
      }
      //
      else if (recTitle == "SET_RESCUT_INIT") {
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fResCutInitial = recOpt.Atof();
      }
      //
      else if (recTitle == "SET_RESCUT_OTHER") {
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fResCut = recOpt.Atof();
      }
      //
      else if (recTitle == "SET_LOCALSIGMAFACTOR") { //-------------------------
	for (irec=0;irec<3;irec++) if (nrecElems>irec+1) {
	    fSigmaFactor[irec] = ((TObjString*)recArr->At(irec+1))->GetString().Atof();
	    if (fSigmaFactor[irec]<=0.) stopped = kTRUE;
	  }
	if (stopped) break; 
      }
      //
      else if (recTitle == "SET_STARTFAC") {        //-------------------------
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fStartFac = recOpt.Atof();
      }
      //
      else if (recTitle == "SET_B_FIELD") {         //-------------------------
	if (recOpt.IsNull() || !recOpt.IsFloat() ) {stopped = kTRUE; break;}
	fBField = recOpt.Atof();
	if (fBField>0) {
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
      else if (recTitle == "SET_SPARSE_MATRIX") {   // matrix solver type
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
      else if (recTitle == "REQUIRE_POINT") {       //-------------------------
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
	  fRequirePoints = kTRUE;
	  if (recOpt == "LAYER") {
	    if (lr<0 || lr>5) {stopped = kTRUE; break;}
	    if (hb>0) fNReqLayUp[lr] = np;
	    else if (hb<0) fNReqLayDown[lr] = np;
	    else fNReqLay[lr] = np;
	  }
	  else if (recOpt == "DETECTOR") {
	    if (lr<0 || lr>2) {stopped = kTRUE; break;}
	    if (hb>0) fNReqDetUp[lr] = np;
	    else if (hb<0) fNReqDetDown[lr] = np;
	    else fNReqDet[lr] = np;
	  }
	  else {stopped = kTRUE; break;}
	}
	else {stopped = kTRUE; break;}
      }
      //
      // global constraints on the subunits/orphans 
      else if (recTitle == "CONSTRAINT_ORPHANS") {    //------------------------
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
      else if (recTitle == "CONSTRAINT_SUBUNITS") {    //------------------------
	// expect ONSTRAINT_SUBUNITS MEAN/MEDIAN Value parID0 ... parID1 ... VolID1 ... VolIDn - VolIDm
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
      else if (recTitle == "APPLY_CONSTRAINT") {            //------------------------
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
      else continue; // already processed record
      //
    } // end of while loop 4 over the various params 
    //
    break;
  } // end of while(1) loop 
  //
  fclose(pfc);
  if (stopped) {
    AliError(Form("Failed on record %s %s ...\n",recTitle.Data(),recOpt.Data()));
    return -1;
  }
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
  /// set the current supermodule
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
  AliGeomManager::LoadGeometry(fGeometryFileName.Data());
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
  fConstrRefFileName = reffname;
  if (fConstrRefFileName == "IDEAL") { // the reference is the ideal geometry, just create dummy reference array
    fConstrRef = new TClonesArray("AliAlignObjParams",1);
    return 0;
  }
  TFile *pref = TFile::Open(fConstrRefFileName.Data());
  if (!pref->IsOpen()) return -2;   
  fConstrRef = (TClonesArray*)pref->Get("ITSAlignObjs");
  pref->Close();
  delete pref;
  if (!fConstrRef) {
    AliError(Form("Did not find reference prealignment deltas in %s",reffname));
    return -1;
  }
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
  AliGeomManager::LoadGeometry(fGeometryFileName.Data());
  fGeoManager = AliGeomManager::GetGeometry();
  return 0;
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::Init()
{
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
	Bool_t AddToFit = kFALSE;	
	// the parameter may be ommitted from explicit fit (if PseudoParentsAllowed is true) if
	// 1) it is not explicitly constrained or its does not participate in Gaussian constraint
	// 2) the same applies to all of its parents
	// 3) it has at least 1 unconstrained direct child
	while(parent) {
	  if (!parent->IsFreeDOF(ipar)) {parent = parent->GetParent(); continue;}
	  nFreeInstances++;
	  if (IsParModConstrained(parent,ipar, cstMeanMed, cstGauss)) nFreeInstances--;
	  if (cstGauss) AddToFit = kTRUE;
	  parent = parent->GetParent();
	}
	if (nFreeInstances>1) {
	  AliError(Form("Parameter#%d of module %s\nhas %d free instances in the "
			"unconstrained parents\nSystem is undefined",ipar,mod->GetName(),nFreeInstances));
	  exit(1);
	}
	//
	// i) Are PseudoParents allowed?
	if (!PseudoParentsAllowed()) AddToFit = kTRUE;
	// ii) check if this module has no child with such a free parameter. Since the order of this check 
	// goes from child to parent, by this moment such a parameter must have been already added
	else if (!IsParModFamilyVaried(mod,ipar))  AddToFit = kTRUE;  // no varied children at all
	else if (!IsParFamilyFree(mod,ipar,1))     AddToFit = kTRUE;  // no unconstrained direct children
	// otherwise the value of this parameter can be extracted from simple contraint and the values of 
	// the relevant parameters of its children the fit is done. Hence it is not included
	if (!AddToFit) continue;
	//
	// shall add this parameter to explicit fit
	//	printf("Adding %s %d -> %d\n",mod->GetName(), ipar, fNGlobal);
	mod->SetParOffset(ipar,fNGlobal++);
      }
    }
  }
  //
  AliInfo(Form("Initializing Millepede with %d gpar, %d lpar and %d stddev ...",fNGlobal, fNLocal, fNStdDev));
  fGlobalDerivatives = new Double_t[fNGlobal];
  memset(fGlobalDerivatives,0,fNGlobal*sizeof(Double_t));
  //
  fMillepede->InitMille(fNGlobal,fNLocal,fNStdDev,fResCut,fResCutInitial);
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
  if (value==0) AliInfo(Form("Parameter %i Fixed", iPar));
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ResetLocalEquation()
{
  /// Reset the derivative vectors
  for(int i=fNLocal;i--;)  fLocalDerivatives[i] = 0.0;
  memset(fGlobalDerivatives, 0, fNGlobal*sizeof(double) );
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::ApplyToGeometry() 
{
  // apply starting realignment to ideal geometry
  AliInfo(Form("Using %s for prealignment",fPreAlignmentFileName.Data()));
  if (!fGeoManager) return -1; 
  TFile *pref = TFile::Open(fPreAlignmentFileName.Data());
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
    //printf("Applying#%4d %s\n",ix,preo->GetSymName());
    if (!preo->ApplyToGeometry()) return -4;
  }
  //
  pref->Close();
  delete pref;
  //
  fUsePreAlignment = kTRUE;
  return 0;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::GetPreAlignmentQualityFactor(Int_t index) const
{
  if (!fUsePreAlignment || index<0 || index>=fPreAlignQF.GetSize()) return -1;
  return fPreAlignQF[index]-1;
}

//________________________________________________________________________________________________________
AliTrackPointArray *AliITSAlignMille2::PrepareTrack(const AliTrackPointArray *atp) 
{
  /// create a new AliTrackPointArray keeping only defined modules
  /// move points according to a given prealignment, if any
  /// sort alitrackpoints w.r.t. global Y direction, if selected
  const double kTiny = 1E-12;
  //
  AliTrackPointArray *atps=NULL;
  Int_t idx[20];
  Int_t npts=atp->GetNPoints();

  /// checks if AliTrackPoints belong to defined modules
  Int_t ngoodpts=0;
  Int_t intidx[20];
  
  for (int j=0; j<npts; j++) {
    intidx[j] = IsVIDContained(atp->GetVolumeID()[j]);
    if (intidx[j]>=0) ngoodpts++;
  }
  AliDebug(3,Form("Number of points in defined modules: %d out of %d",ngoodpts,npts));

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
  // build a new track with (sorted) (prealigned) good points
  atps = (AliTrackPointArray*)fTrackBuff[ngoodpts-fMinNPtsPerTrack];
  if (!atps) {
    atps = new AliTrackPointArray(ngoodpts);
    fTrackBuff.AddAtAndExpand(atps,ngoodpts-fMinNPtsPerTrack);
  }
  //
  //
  for (int i=0; i<npts; i++) idx[i]=i;
  // sort track if required
  TMath::Sort(npts,atp->GetY(),idx); // sort descending...
  //
  Int_t npto=0;
  for (int i=0; i<npts; i++) {
    // skip not defined points
    if (intidx[idx[i]]<0) continue;
    atp->GetPoint(p,idx[i]);

    // prealign point if required
    // get IDEAL matrix
    AliITSAlignMille2Module *mod = GetMilleModule(intidx[idx[i]]);
    TGeoHMatrix *svOrigMatrix = mod->GetSensitiveVolumeOrigGlobalMatrix(p.GetVolumeID());
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
    //
    // this is a temporary code to extract the drift speed used for given point
    if (p.GetDriftTime()>0) { // RRR
      // calculate the drift speed
      int sid = AliITSAlignMille2Module::GetIndexFromVolumeID(p.GetVolumeID());// - kSDDoffsID;
      fDriftTime0[npto] = fInitialRecSDD ? fInitialRecSDD->GetTimeZero(sid) : 0.;
      /*
      AliGeomManager::ELayerID lay = AliGeomManager::VolUIDToLayer(p.GetVolumeID());
      if      (lay==3) fDriftTime0[npto] = pg[2]<0 ? 169.5 : 140.1;
      else if (lay==4) fDriftTime0[npto] = pg[2]<0 ? 158.3 : 139.0;
      else {
	AliError(Form("Strange layer %d for moduleID %d",lay,p.GetVolumeID()));
	exit(1);
      }
      */
      double tdif = p.GetDriftTime() - fDriftTime0[npto];
      if (tdif<=0) tdif = 1;
      double vdrift = (3.5085-TMath::Abs(pl[0]))/tdif;
      if (vdrift<0) vdrift = 0;
      //
      // TEMPORARY CORRECTION (if provided) -------------->>>
      if (fCorrectSDD) {
	float t0Upd = fCorrectSDD->GetTimeZero(sid);
	vdrift += fCorrectSDD->GetDeltaVDrift(sid);
	tdif    = p.GetDriftTime() - t0Upd;
	// correct Xlocal
	pl[0] = TMath::Sign(3.5085 - vdrift*tdif,pl[0]);
	fDriftTime0[npto] =  t0Upd;
      }
      // TEMPORARY CORRECTION (if provided) --------------<<<
      fDriftSpeed[npto] = TMath::Sign(vdrift,pl[0]);
      //
      /*
      printf("%d  %+6.2f %+6.2f %+6.2f  %+5.2f %+5.2f %+5.2f  %+6.1f  %+6.1f %+f %+f\n",
	     p.GetVolumeID(),pg[0],pg[1],pg[2],pl[0],pl[1],pl[2],p.GetDriftTime(), fDriftTime0[npto], fDriftSpeed[npto],tdif);
      */
    }

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
    //    printf("\nErrMatGlob: before\n"); hcov.Print(""); //RRR
    hcov.Multiply(svOrigMatrix);
    hcov.MultiplyLeft(&svOrigMatrix->Inverse());
    // now hcov is LOCAL COVARIANCE MATRIX
    // apply sigma scaling
    //    printf("\nErrMatLoc: before\n"); hcov.Print(""); //RRR
    Double_t *hcovscl = hcov.GetRotationMatrix(); 
    //    for (int ir=3;ir--;) for (int ic=3;ic--;) hcovscl[ir*3+ic] *= mod->GetSigmaFactor(ir)*mod->GetSigmaFactor(ic); //RRR
    // RS TEMPORARY: nullify non-diagonal elements and sigY
    hcovscl[5] = 0;
    for (int ir=3;ir--;) for (int ic=3;ic--;) {
	if (ir==ic) {
	  if (TMath::Abs(hcovscl[ir*3+ic])<kTiny) hcovscl[ir*3+ic] = 0.;
	  else hcovscl[ir*3+ic] *= mod->GetSigmaFactor(ir)*mod->GetSigmaFactor(ic); //RRR
	}
	else hcovscl[ir*3+ic]  = 0;
      }
    //
    //    printf("\nErrMatLoc: after\n"); hcov.Print(""); //RRR
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
    TGeoHMatrix *svMatrix = mod->GetSensitiveVolumeMatrix(p.GetVolumeID());
    // modify global coordinates according with pre-aligment
    svMatrix->LocalToMaster(pl,pg);
    // now rotate in local system
    hcov.Multiply(&svMatrix->Inverse());
    hcov.MultiplyLeft(svMatrix);
    // hcov is back in GLOBAL RF
    // cure once more
    for (int ir=3;ir--;) for (int ic=3;ic--;) if (TMath::Abs(hcovscl[ir*3+ic])<kTiny) hcovscl[ir*3+ic] = 0.;
    //    printf("\nErrMatGlob: after\n"); hcov.Print(""); //RRR
    //
    Float_t pcov[6];
    pcov[0]=hcovscl[0];
    pcov[1]=hcovscl[1];
    pcov[2]=hcovscl[2];
    pcov[3]=hcovscl[4];
    pcov[4]=hcovscl[5];
    pcov[5]=hcovscl[8];

    p.SetXYZ(pg[0],pg[1],pg[2],pcov);
    //    printf("New Gl coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]);
    AliDebug(3,Form("New global coordinates of measured point : X=%f  Y=%f  Z=%f \n",pg[0],pg[1],pg[2]));
    atps->AddPoint(npto,&p);
    AliDebug(2,Form("Adding point[%d] = ( %f , %f , %f )     volid = %d",npto,atps->GetX()[npto],
		    atps->GetY()[npto],atps->GetZ()[npto],atps->GetVolumeID()[npto] ));
    //    printf("Adding %d %d %f\n",npto, p.GetVolumeID(), p.GetY()); 
    npto++;
  }

  return atps;
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
    AliDebug(2,Form("Point[%d] = ( %f , %f , %f )     volid = %d",i,atps->GetX()[i],atps->GetY()[i],atps->GetZ()[i],atps->GetVolumeID()[i] ));
  }
  return atps;
}

//________________________________________________________________________________________________________
Int_t AliITSAlignMille2::GetCurrentLayer() const 
{
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
  //
  // IT IS VERY IMPORTANT: start from the end of the list, where the childs are located !!!
  Int_t k=fNModules-1;
  fCurrentModule = 0;
  // VERY IMPORTANT: if the sensors were explicitly provided, don't look in the supermodules  
  while (k>=0 && ! (fCurrentModule=GetMilleModule(k))->IsIn(voluid)) k--;
  if (k<0) return -3;
  //
  /*
  // Check if the module has free params. If not, go over the parents
  AliITSAlignMille2Module* mdtmp = fCurrentModule;
  while (mdtmp && mdtmp->GetNParFree()==0) mdtmp = mdtmp->GetParent();
  if (!mdtmp) return 1; // nothing to vary here
  fCurrentModule = mdtmp;
  */
  //
  fModuleInitParam[0] = 0.0;
  fModuleInitParam[1] = 0.0;
  fModuleInitParam[2] = 0.0;
  fModuleInitParam[3] = 0.0; // psi   (X)
  fModuleInitParam[4] = 0.0; // theta (Y)
  fModuleInitParam[5] = 0.0; // phi   (Z)
  fModuleInitParam[6] = 0.0;
  fModuleInitParam[7] = 0.0;
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
  //
  // set local sigmas
  fSigmaLoc[0] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[0]));
  fSigmaLoc[1] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[4])); // RS
  fSigmaLoc[2] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[8]));

  // set minimum value for SigmaLoc to 10 micron 
  if (fSigmaLoc[0]<0.0010) fSigmaLoc[0]=0.0010;
  if (fSigmaLoc[2]<0.0010) fSigmaLoc[2]=0.0010;
  //
  /* RRR the rescaling is moved to PrepareTrack
  // multiply local sigmas by global and module specific factor 
  for (int i=3;i--;) fSigmaLoc[i] *= fSigmaFactor[i]*fCurrentModule->GetSigmaFactor(i);
  //
  */
  AliDebug(2,Form("Setting StDev from CovMat : fSigmaLocX=%g  fSigmaLocY=%g fSigmaLocZ=%g \n",fSigmaLoc[0] ,fSigmaLoc[1] ,fSigmaLoc[2] ));
   
  return 0;
}

//________________________________________________________________________________________________________
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
  printf("        init value for chi2 cut       : %.4f\n",fStartFac);
  printf("        first iteration cut value     : %.4f\n",fResCutInitial);
  printf("        other iterations cut value    : %.4f\n",fResCut);
  printf("        number of stddev for chi2 cut : %d\n",fNStdDev);
  printf("        def.scaling for local sigmas  : %.4f %.4f %.4f\n",fSigmaFactor[0],fSigmaFactor[1],fSigmaFactor[2]);

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
  if (det==0) det = 1E-20;
  fLocalInitParam[0] = (sX*sYY-sY*sXY)/det;
  fLocalInitParam[2] = (sXY*npts-sY*sX)/det;
  //
  fLocalInitParam[1] = (sZ*sYY-sY*sZY)/det;
  fLocalInitParam[3] = (sZY*npts-sY*sZ)/det;
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
  fullErr2D = kFALSE;//kTRUE;
  minuit->ExecuteCommand("MIGRAD",arglist,2);
  //
  for (int i=0;i<4;i++) fLocalInitParam[i] = minuit->GetParameter(i);
  for (int i=0;i<4;i++) for (int j=0;j<4;j++) fLocalInitParEr[i][j] = minuit->GetCovarianceMatrixElement(i,j);
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
Int_t AliITSAlignMille2::ProcessTrack(const AliTrackPointArray *track) 
{
  /// Process track; Loop over hits and set local equations
  /// here 'track' is a AliTrackPointArray
  /// return 0 if success;
  
  if (!fIsMilleInit) Init();
  //
  Int_t npts = track->GetNPoints();
  AliDebug(2,Form("*** Input track with %d points ***",npts));

  // preprocessing of the input track: keep only points in defined volumes,
  // move points if prealignment is set, sort by Yglo if required
  
  fTrack=PrepareTrack(track);
  if (!fTrack) return -1;

  npts = fTrack->GetNPoints();
  if (npts>kMaxPoints) {
    AliError(Form("Compiled with kMaxPoints=%d, current track has %d points",kMaxPoints,npts));
  }
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
      fTrack=NULL;
      return -5;
    }
  }
  
  Int_t nloceq=0;
  Int_t ngloeq=0;
  static Mille2Data md[kMaxPoints];
  //
  for (Int_t ipt=0; ipt<npts; ipt++) {
    fTrack->GetPoint(fCluster,ipt);
    fCluster.SetUniqueID(ipt);
    AliDebug(2,Form("\n--- processing point %d --- \n",ipt));    

    // set geometry parameters for the the current module
    if (InitModuleParams()) continue;
    AliDebug(2,Form("    VolID=%d  Index=%d  InternalIdx=%d  symname=%s\n", 
		    track->GetVolumeID()[ipt], fCurrentModule->GetIndex(),
		    fCurrentModule->GetUniqueID(), AliGeomManager::SymName(track->GetVolumeID()[ipt]) ));
    AliDebug(2,Form("    Preprocessed Point = ( %f , %f , %f ) \n",fCluster.GetX(),fCluster.GetY(),fCluster.GetZ()));
    int res = AddLocalEquation(md[nloceq]);
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
  AliDebug(3,Form("Intesect. point: L( %f , %f , %f )  G( %f , %f , %f )\n",fPintLoc[0],fPintLoc[1],fPintLoc[2],fPintGlo[0],fPintGlo[1],fPintGlo[2]));
  
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
  // store first interaction point
  if (CalcIntersectionPoint(fLocalInitParam, fModuleInitParam)) return -1;  
  for (Int_t i=0; i<3; i++) fPintLoc0[i]=fPintLoc[i];
  AliDebug(2,Form("Intesect. point: L( %f , %f , %f )",fPintLoc[0],fPintLoc[1],fPintLoc[2]));
  
  // calculate local derivatives numerically
  Bool_t zeroX = kTRUE;
  Bool_t zeroZ = kTRUE;
  //
  for (Int_t i=0; i<fNLocal; i++) {
    if (CalcDerivatives(i,kTRUE)) return -1;
    m.derlocX[i] = fDerivativeLoc[i][0];
    m.derlocZ[i] = fDerivativeLoc[i][2];
    if (zeroX) zeroX = fDerivativeLoc[i][0]==0;
    if (zeroZ) zeroZ = fDerivativeLoc[i][2]==0;
  }
  //  for (Int_t i=0; i<fNLocal; i++) AliDebug(2,Form("Local parameter %d - dXdpar = %g  - dZdpar = %g\n",i,dXdL[i],dZdL[i]));
  //
  if (zeroX) {AliInfo("Aborting: zero local X derivatives!"); return -1;}
  if (zeroZ) {AliInfo("Aborting: zero local Z derivatives!"); return -1;}
  //
  int ifill = 0;
  //
  AliITSAlignMille2Module* endModule = fCurrentModule;
  //
  zeroX = zeroZ = kTRUE;
  Bool_t dfDone[kNParCh];
  for (int i=kNParCh;i--;) dfDone[i] = kFALSE;
  m.nModFilled = 0;
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
	  if (zeroX) zeroX = fDerivativeGlo[i][0]==0;
	  if (zeroZ) zeroZ = fDerivativeGlo[i][2]==0;
	}
      }
      //
      m.dergloX[ifill] = fDerivativeGlo[i][0];
      m.dergloZ[ifill] = fDerivativeGlo[i][2];
      m.parMilleID[ifill++] = fCurrentModule->GetParOffset(i);
    }
    //
    // specific for special sensors
    if ( fCurrentModule->IsSDD() && 
	 (fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFT0)>=0 ||
	  fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFDV)>=0) ) {
      //
      // assume for sensor local xloc = xloc0 + V0*dT0+dV*(T-T0)
      // where V0 and T are the nominal drift velocity, time and time0
      // and the dT0 and dV are the corrections:
      // dX/dT0 = dX/dxloc * dxloc/dT0 = dX/dxloc * V0
      // dX/dV  = dX/dxloc * dxloc/dV =  dX/dxloc * (T-T0)
      // IMPORTANT: the geom derivatives are over the SENSOR LOCAL parameters
      //
      if (!dfDone[AliITSAlignMille2Module::kDOFT0] || !dfDone[AliITSAlignMille2Module::kDOFDV]) {
	//
	double dXdxlocsens=0., dZdxlocsens=0.;
	//
	// if the current module is the sensor itself and we work with local params, then 
	// we can directly take dX/dxloc_sens dZ/dxloc_sens
	if (!fUseGlobalDelta && fCurrentModule->GetVolumeID()==fCluster.GetVolumeID()) {
	  if (dfDone[AliITSAlignMille2Module::kDOFTX]) {
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
	if (zeroX) zeroX = dXdxlocsens == 0;
	if (zeroZ) zeroZ = dZdxlocsens == 0;
	//
	double vdrift = GetVDriftSDD();
	double tdrift = GetTDriftSDD();
	//
	fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][0] = dXdxlocsens*vdrift;
	fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][2] = dZdxlocsens*vdrift;
	dfDone[AliITSAlignMille2Module::kDOFT0] = kTRUE;
	//
	fDerivativeGlo[AliITSAlignMille2Module::kDOFDV][0] = -dXdxlocsens*TMath::Sign(tdrift,vdrift);
	fDerivativeGlo[AliITSAlignMille2Module::kDOFDV][2] = -dZdxlocsens*TMath::Sign(tdrift,vdrift);
	dfDone[AliITSAlignMille2Module::kDOFDV] = kTRUE;
	//
      }
      //
      if (fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFT0)>=0) {
	m.dergloX[ifill] = fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][0];
	m.dergloZ[ifill] = fDerivativeGlo[AliITSAlignMille2Module::kDOFT0][2];
	m.parMilleID[ifill++] = fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFT0);      
      }
      //
      if (fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFDV)>=0) {
	m.dergloX[ifill] = fDerivativeGlo[AliITSAlignMille2Module::kDOFDV][0];
	m.dergloZ[ifill] = fDerivativeGlo[AliITSAlignMille2Module::kDOFDV][2];
	m.parMilleID[ifill++] = fCurrentModule->GetParOffset(AliITSAlignMille2Module::kDOFDV);      
      }
    }
    //
    m.moduleID[m.nModFilled++] = fCurrentModule->GetUniqueID();
  } while( (fCurrentModule=fCurrentModule->GetParent()) );
  //
  if (nmodTested>0 && zeroX) {AliInfo("Aborting: zero global X derivatives!");return -1;}
  if (nmodTested>0 && zeroZ) {AliInfo("Aborting: zero global Z derivatives!");return -1;}
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
  m.nGlobFilled = ifill;
  fCurrentModule = endModule;
  //
  return Int_t(!zeroX && !zeroZ);
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
    // set equation for Xloc coordinate
    AliDebug(2,Form("setting local equation X with fMeas=%.6f  and fSigma=%.6f",m.measX, m.sigmaX));
    for (int i=fNLocal; i--;) SetLocalDerivative( i, m.derlocX[i] );
    for (int i=m.nGlobFilled;i--;) SetGlobalDerivative( m.parMilleID[i] , m.dergloX[i] );
    fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, m.measX, m.sigmaX);  
    //
    // set equation for Zloc coordinate
    AliDebug(2,Form("setting local equation Z with fMeas=%.6f  and fSigma=%.6f",m.measZ, m.sigmaZ));
    for (int i=fNLocal; i--;) SetLocalDerivative( i, m.derlocZ[i] );
    for (int i=m.nGlobFilled;i--;) SetGlobalDerivative( m.parMilleID[i] , m.dergloZ[i] );
    fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, m.measZ, m.sigmaZ);  
    //
    for (int i=m.nModFilled;i--;) GetMilleModule(m.moduleID[i])->IncNProcessedPoints();
    //
  }
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
  if (fIsMilleInit) {
    AliInfo("Millepede has been already initialized: no constrain may be added!");
    return;
  }
  AliITSAlignMille2ConstrArray *cstr = new AliITSAlignMille2ConstrArray(name,parcf,npar,val,err);
  cstr->SetConstraintID(GetNConstraints());
  fConstraints.Add(cstr);
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ApplyGaussianConstraint(AliITSAlignMille2ConstrArray* cstr)
{
  // apply the constraint on the local corrections of a list of modules
  int nmod = cstr->GetNModules();
  double jacobian[AliITSAlignMille2Module::kMaxParGeom][AliITSAlignMille2Module::kMaxParGeom];
  //
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
      if (coef==0) continue;
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
}

//________________________________________________________________________________________________________
void AliITSAlignMille2::ApplyPostConstraints()
{
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
    AliInfo(Form("%s constraint: added %f shift to param[%d] of %d children of module %d: %s",
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
    AliInfo(Form("%s constraint: added %f shift to param[%d] of %d orphan modules",
		 type==AliITSAlignMille2Constraint::kTypeMean ? "MEAN" : "MEDIAN",shift,
		 ip,npc));
  }
  delete[] tmpArr;  
  //
}

//________________________________________________________________________________________________________
Bool_t AliITSAlignMille2::IsParModConstrained(AliITSAlignMille2Module* mod,Int_t par, Bool_t &meanmed, Bool_t &gaussian) const
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
Bool_t AliITSAlignMille2::IsParModFamilyVaried(AliITSAlignMille2Module* mod,Int_t par,Int_t depth) const
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
Bool_t AliITSAlignMille2::IsParFamilyFree(AliITSAlignMille2Module* mod,Int_t par,Int_t depth) const
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
  double t = fCluster.GetDriftTime();
  return t - fDriftTime0[ fCluster.GetUniqueID() ];
}

//________________________________________________________________________________________________________
Double_t AliITSAlignMille2::GetVDriftSDD() const 
{
  return fDriftSpeed[ fCluster.GetUniqueID() ];
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
  double pars[AliITSAlignMille2Module::kMaxParGeom];
  for (int imd=fNModules;imd--;) {
    AliITSAlignMille2Module* mod = GetMilleModule(imd);
    if (!mod->GeomParamsGlobal()) continue;
    mod->GetGeomParamsLoc(pars);
    mod->SetParVals(pars,AliITSAlignMille2Module::kMaxParGeom);
    mod->SetGeomParamsGlobal(kFALSE);
  }
}


