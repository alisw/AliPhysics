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
#include <TGraph.h>
#include <TGeoMatrix.h>
#include <TMath.h>

#include "AliITSAlignMille.h"
#include "AliITSgeomTGeo.h"
#include "AliGeomManager.h"
#include "AliMillepede.h"
#include "AliTrackPointArray.h"
#include "AliAlignObjParams.h"
#include "AliLog.h"
#include "TSystem.h"  // come si fa?

/// \cond CLASSIMP
ClassImp(AliITSAlignMille)
/// \endcond
  
Int_t AliITSAlignMille::fgNDetElem = ITSMILLE_NDETELEM;
Int_t AliITSAlignMille::fgNParCh = ITSMILLE_NPARCH;

AliITSAlignMille::AliITSAlignMille(const Char_t *configFilename, Bool_t initmille) 
  : TObject(),
    fMillepede(0),
    fStartFac(16.), 
    fResCutInitial(100.), 
    fResCut(100.),
    fNGlobal(ITSMILLE_NDETELEM*ITSMILLE_NPARCH),
    fNLocal(ITSMILLE_NLOCAL),
    fNStdDev(ITSMILLE_NSTDEV),
    fIsMilleInit(kFALSE),
    fParSigTranslations(0.0100),
    fParSigRotations(0.1),
    fTrack(NULL),
    fCluster(),
    fGlobalDerivatives(NULL),
    fTempHMat(NULL),
    fTempAlignObj(NULL),
    fDerivativeXLoc(0),
    fDerivativeZLoc(0),
    fDeltaPar(0),
    fMinNPtsPerTrack(3),
    fGeometryFileName("geometry.root"),
    fGeoManager(0),
    fCurrentModuleIndex(0),
    fCurrentModuleInternalIndex(0),
    fNModules(0),
    fUseLocalShifts(kTRUE),
    fCurrentModuleHMatrix(NULL)
{
  /// main constructor that takes input from configuration file
  
  fMillepede = new AliMillepede();
  fGlobalDerivatives = new Double_t[fNGlobal];
  fTempHMat = new TGeoHMatrix;
  fCurrentModuleHMatrix = new TGeoHMatrix;
  
  Int_t lc=LoadConfig(configFilename);
  if (lc) {
    AliInfo(Form("Error %d loading configuration from %s",lc,configFilename));
  }
  else {    
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
  
  fDeltaPar=0.0; // not used at the moment - to be checked later...
  
}

AliITSAlignMille::~AliITSAlignMille() {
  /// Destructor
  if (fMillepede) delete fMillepede;
  delete [] fGlobalDerivatives;
  delete fCurrentModuleHMatrix;
  delete fTempHMat;
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
  Float_t f1;
  UShort_t voluid;
  Int_t nmod=0;

  while (fgets(st,200,pfc)) {

    // skip comments
    for (int i=0; i<int(strlen(st)); i++) {
      if (st[i]=='#') st[i]=0;
    }

    if (strstr(st,"GEOMETRY_FILE")) {
      sscanf(st,"%s %s",tmp,st2);
      if (gSystem->AccessPathName(st2)) {
	AliInfo("*** WARNING! *** geometry file not found! ");
	return -1;
      }  
      fGeometryFileName=st2;
      InitGeometry();
    }

    if (strstr(st,"SET_PARSIG_TRA")) {
      sscanf(st,"%s %f",tmp,&f1);
      fParSigTranslations=f1;
    }

    if (strstr(st,"SET_PARSIG_ROT")) {
      sscanf(st,"%s %f",tmp,&f1);
      fParSigRotations=f1;
    }

    if (strstr(st,"SET_NSTDDEV")) {
      sscanf(st,"%s %d",tmp,&idx);
      fNStdDev=idx;
    }

    if (strstr(st,"SET_RESCUT_INIT")) {
      sscanf(st,"%s %f",tmp,&f1);
      fResCutInitial=f1;
    }

    if (strstr(st,"SET_RESCUT_OTHER")) {
      sscanf(st,"%s %f",tmp,&f1);
      fResCut=f1;
    }

    if (strstr(st,"SET_STARTFAC")) {
      sscanf(st,"%s %f",tmp,&f1);
      fStartFac=f1;
    }

    if (strstr(st,"SET_LOCAL_SHIFTS")) { // only enabled mode...
      fUseLocalShifts = kTRUE;
    }

    if (strstr(st,"MODULE_INDEX")) {
      sscanf(st,"%s %d %d %d %d %d %d %d",tmp,&idx,&itx,&ity,&itz,&iph,&ith,&ips);
      voluid=GetModuleVolumeID(idx);
      if (!voluid) return 1; // bad index
      fModuleIndex[nmod]=idx;
      fModuleVolumeID[nmod]=voluid;
      fFreeParam[nmod][0]=itx;
      fFreeParam[nmod][1]=ity;
      fFreeParam[nmod][2]=itz;
      fFreeParam[nmod][3]=iph;
      fFreeParam[nmod][4]=ith;
      fFreeParam[nmod][5]=ips;
      nmod++;
    }
   
    if (strstr(st,"MODULE_VOLUID")) {
      // to be implemented
    }

  } // end while

  fNModules = nmod;
  fNGlobal = fNModules*fgNParCh;
 
  fclose(pfc);
  return 0;
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
  if (index<0 || index>2197) return 0;
  return GetModuleVolumeID(AliITSgeomTGeo::GetSymName(index));
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
  fTempAlignObj=new AliAlignObjParams(AliITSgeomTGeo::GetSymName(7),2055,0,0,0,0,0,0,kFALSE);
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
    for (Int_t j=0; j<ITSMILLE_NPARCH; j++) {
      if (!fFreeParam[i][j]) FixParameter(i*ITSMILLE_NPARCH+j,0.0);
      else {
	// pepopepo: da sistemare il settaggio delle sigma...
	Double_t parsig=0;
	if (j<3) parsig=fParSigTranslations; // translations (0.0100 cm)
	else parsig=fParSigRotations; // rotations (1/10 deg)
	FixParameter(i*ITSMILLE_NPARCH+j,parsig);
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

Int_t AliITSAlignMille::InitModuleParams() {
  /// initialize geometry parameters for a given detector
  /// for current cluster (fCluster)
  /// fGlobalInitParam[] is set as:
  ///    [tx,ty,tz,psi,theta,phi]
  ///    (old was [tx,ty,tz,theta,psi,phi] ROOT's angles...)
  /// *** At the moment: using Raffalele's angles definition ***
  ///
  /// Num of Dets: ITSMILLE_NDETELEM = fgNDetElem
  /// Num of Par : ITSMILLE_NPARCH = fgNParCh
  /// return 0 if success

  if (!fGeoManager) {
    AliInfo("ITS geometry not initialized!");
    return -1;
  }

  // set the internal index (index in module list)
  UShort_t voluid=fCluster.GetVolumeID();
  Int_t k=fNModules-1;
  while (k>=0 && !(voluid==fModuleVolumeID[k]) ) k--;  
  if (k<0) return -3;    
  fCurrentModuleInternalIndex=k;

  // set the index
  Int_t index = GetModuleIndex(AliGeomManager::SymName(voluid));
  if (index<0) return -2;
  fCurrentModuleIndex = index;

  fModuleInitParam[0] = 0.0;
  fModuleInitParam[1] = 0.0;
  fModuleInitParam[2] = 0.0;
  fModuleInitParam[3] = 0.0; // psi   (X)
  fModuleInitParam[4] = 0.0; // theta (Y)
  fModuleInitParam[5] = 0.0; // phi   (Z)

  /// get global (corrected) matrix  
  //  if (!AliITSgeomTGeo::GetOrigMatrix(index,*fCurrentModuleHMatrix)) return -3;
  Double_t rott[9];
  if (!AliITSgeomTGeo::GetRotation(index,rott)) return -3;
  fCurrentModuleHMatrix->SetRotation(rott);
  Double_t oLoc[3]={0,0,0};
  if (!AliITSgeomTGeo::LocalToGlobal(index,oLoc,fCurrentModuleTranslation)) return -4;
  fCurrentModuleHMatrix->SetTranslation(fCurrentModuleTranslation);

  /// get back local coordinates
  fMeasGlo[0] = fCluster.GetX();
  fMeasGlo[1] = fCluster.GetY();
  fMeasGlo[2] = fCluster.GetZ();
  fCurrentModuleHMatrix->MasterToLocal(fMeasGlo,fMeasLoc);

  // set stdev from cluster
  TGeoHMatrix hcov;
  Double_t hcovel[9];
  hcovel[0]=double(fCluster.GetCov()[0]);
  hcovel[1]=double(fCluster.GetCov()[1]);
  hcovel[2]=double(fCluster.GetCov()[3]);
  hcovel[3]=double(fCluster.GetCov()[1]);
  hcovel[4]=double(fCluster.GetCov()[2]);
  hcovel[5]=double(fCluster.GetCov()[4]);
  hcovel[6]=double(fCluster.GetCov()[3]);
  hcovel[7]=double(fCluster.GetCov()[4]);
  hcovel[8]=double(fCluster.GetCov()[5]);
  hcov.SetRotation(hcovel);
  // now rotate in local system
  hcov.MultiplyLeft(&fCurrentModuleHMatrix->Inverse());
  hcov.Multiply(fCurrentModuleHMatrix);

  // per i ruotati c'e' delle sigmaY che compaiono... prob
  // e' un problema di troncamento
  fSigmaLoc[0] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[0]));
  fSigmaLoc[1] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[4]));
  fSigmaLoc[2] = TMath::Sqrt(TMath::Abs(hcov.GetRotationMatrix()[8]));

    AliDebug(2,Form("Setting StDev from CovMat : fSigmaLocX=%f  fSigmaLocY=%f fSigmaLocZ=%f \n",fSigmaLoc[0] ,fSigmaLoc[1] ,fSigmaLoc[2] ));
   
  return 0;
}

void AliITSAlignMille::SetCurrentModule(Int_t index) {
  ///
  UShort_t voluid=GetModuleVolumeID(index);
  if (voluid) {
    fCluster.SetVolumeID(voluid);
    fCluster.SetXYZ(0,0,0);
    InitModuleParams();
  }
}

void AliITSAlignMille::Print() {
  ///
  printf("*** AliMillepede for ITS ***\n");
  printf("    number of defined modules: %d\n",fNModules);
  if (fGeoManager)
    printf("    geometry loaded from %s\n",fGeometryFileName.Data());
  else
    printf("    geometry not loaded\n");
  if (fUseLocalShifts) 
    printf("    Alignment shifts will be computed in LOCAL RS\n");
  else
    printf("    Alignment shifts will be computed in GLOBAL RS\n");    
  printf("    Millepede configuration parameters:\n");
  printf("       init parsig for translations  : %.4f\n",fParSigTranslations);
  printf("       init parsig for rotations     : %.4f\n",fParSigRotations);
  printf("       init value for chi2 cut       : %.4f\n",fStartFac);
  printf("       first iteration cut value     : %.4f\n",fResCutInitial);
  printf("       other iterations cut value    : %.4f\n",fResCut);
  printf("       number of stddev for chi2 cut : %d\n",fNStdDev);

  printf("List of defined modules:\n");
  printf("  intidx\tindex\tvoluid\tsymname\n");
  for (int i=0; i<fNModules; i++)
    printf("  %d\t%d\t%d\t%s\n",i,fModuleIndex[i],fModuleVolumeID[i],AliITSgeomTGeo::GetSymName(fModuleIndex[i]));
}

void AliITSAlignMille::PrintCurrentModuleInfo() {
  ///
  if (fCurrentModuleIndex<0 || fCurrentModuleIndex>2197) return;
  UShort_t voluid=fModuleVolumeID[fCurrentModuleInternalIndex];
  printf("Current module: index=%d   voluid=%d\n",fCurrentModuleIndex,voluid);
  printf("                symname:%s\n",AliGeomManager::SymName(voluid));
  printf("  TGeoHMatrix: \n");
  fCurrentModuleHMatrix->Print();
}


void AliITSAlignMille::InitTrackParams(int meth) {
  /// initialize local parameters with different methods
  /// for current track (fTrack)
  
  switch (meth) {
  case 1:   // simple linear interpolation
    // get local starting parameters (to be substituted by ESD track parms)
    // local parms (fLocalInitParam[]) are:
    //      [0] = global x coord. of straight line intersection at y=0 plane
    //      [1] = global z coord. of straight line intersection at y=0 plane
    //      [2] = px/py  
    //      [3] = pz/py
    
    // test #1: linear fit in x(y) and z(y)
    Int_t npts = fTrack->GetNPoints();

    TF1 *f1=new TF1("f1","[0]+x*[1]",-50,50);

    TGraph *g=new TGraph(npts,fTrack->GetY(),fTrack->GetX());
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
  }

}
Bool_t AliITSAlignMille::CheckVolumeID(UShort_t voluid) const 
{
  ///
  Int_t k=fNModules-1;
  while (k>=0 && !(voluid==fModuleVolumeID[k]) ) k--;  
  //printf("selected element with voluid=%d : %d\n",voluid,k);
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
  // pepo da controllare...
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
  AliDebug(2,Form("\n*** Processing track with %d points ***\n",npts));
  fTrack = track;

  // checks if there are enough good points
  if (!CheckCurrentTrack()) {
    AliInfo("Track with not enough good points - will not be used...");
    return -1;
  }

  // set local starting parameters (to be substituted by ESD track parms)
  // local parms (fLocalInitParam[]) are:
  //      [0] = global x coord. of straight line intersection at y=0 plane
  //      [1] = global z coord. of straight line intersection at y=0 plane
  //      [2] = px/py  
  //      [3] = pz/py
  InitTrackParams(1);  

  for (Int_t ipt=0; ipt<npts; ipt++) {
    fTrack->GetPoint(fCluster,ipt);
    if (!CheckVolumeID(fCluster.GetVolumeID())) continue;

    // set geometry parameters for the the current module
    AliDebug(2,Form("\n--- processing point %d --- \n",ipt));    
    if (InitModuleParams()) continue;

    AliDebug(2,Form("    VolID=%d  Index=%d  InternalIdx=%d  symname=%s\n", track->GetVolumeID()[ipt], fCurrentModuleIndex ,fCurrentModuleInternalIndex, AliGeomManager::SymName(track->GetVolumeID()[ipt]) ));
    AliDebug(2,Form("    Point = ( %f , %f , %f ) \n",track->GetX()[ipt],track->GetY()[ipt],track->GetZ()[ipt]));
    
    if (SetLocalEquations()) return -1;    

  } // end loop over points
  
  return 0;
}

Int_t AliITSAlignMille::CalcIntersectionPoint(Double_t *lpar, Double_t *gpar) {
  /// calculate track intersection point in local coordinates
  /// according with a given set of parameters (local(4) and global(6))
  /// and fill fPintLoc/Glo
  ///    local are:   pgx0, pgz0, ugx0, ugz0
  ///    global are:  tx,ty,tz,psi,theta,phi (Raff's delta angles in deg.)
  /// return 0 if success
  
  AliDebug(3,Form("lpar = %g %g %g %g \ngpar= %g %g %g %g %g %g\n",lpar[0],lpar[1],lpar[2],lpar[3],gpar[0],gpar[1],gpar[2],gpar[3],gpar[4],gpar[5]));
  
  // vector of initial straight line direction in glob. coord
  // ATTENZIONE: forse va rinormalizzato tutto...
  Double_t v0g[3];
  //Double_t 
  v0g[0]=lpar[2];
  v0g[1]=1.0;
  v0g[2]=lpar[3];

  // intercept in yg=0 plane in glob coord
  Double_t p0g[3];
  p0g[0]=lpar[0];
  p0g[1]=0.0;
  p0g[2]=lpar[1];


  // prepare the TGeoHMatrix
  Double_t tr[3],ang[3];
  //Double_t rad2deg=180./TMath::Pi();
  if (fUseLocalShifts) { // just Delta matrix
    tr[0]=gpar[0]; 
    tr[1]=gpar[1]; 
    tr[2]=gpar[2];
    ang[0]=gpar[3]; // psi   (X)
    ang[1]=gpar[4]; // theta (Y)
    ang[2]=gpar[5]; // phi   (Z)
  }
  else { // total matrix with shifted parameter
    AliInfo("global shifts not implemented yet!");
    return -1;
  }

  //printf("fTempRot = 0x%x  - ang = %g %g %g \n",fTempRot,gpar[5]*rad2deg,gpar[3]*rad2deg,gpar[4]*rad2deg);

  fTempAlignObj->SetRotation(ang[0],ang[1],ang[2]);
  AliDebug(3,Form("Delta angles: psi=%f  theta=%f   phi=%f",ang[0],ang[1],ang[2]));
  TGeoHMatrix hm;
  fTempAlignObj->GetMatrix(hm);
  fTempHMat->SetRotation(hm.GetRotationMatrix());
  fTempHMat->SetTranslation(tr);
  
  // in this case the gpar[] array contains only shifts
  // and fInitModuleParam[] are set to 0
  // fCurrentModuleHMatrix is then modified as fCurrentHM*fTempHM
  if (fUseLocalShifts) 
    fTempHMat->MultiplyLeft(fCurrentModuleHMatrix);

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
  /// parlist: local  (islpar=kTRUE)  pgx0, pgz0, ugx0, ugz0
  ///          global (islpar=kFALSE) tx, ty, tz, psi, theta, phi (Raf's angles in deg)
  /// return 0 if success
  
  // copy initial parameters
  Double_t lpar[ITSMILLE_NLOCAL];
  Double_t gpar[ITSMILLE_NPARCH];
  for (Int_t i=0; i<ITSMILLE_NLOCAL; i++) lpar[i]=fLocalInitParam[i];
  for (Int_t i=0; i<ITSMILLE_NPARCH; i++) gpar[i]=fModuleInitParam[i];

  // pepopepo  
  // trial with fixed dpar...
  Double_t dpar=0.0;
  if (islpar) {
    //dpar=fLocalInitParam[paridx]*0.001;
    // set minimum dpar
    if (paridx<2) dpar=1.0e-4; // translations
    else dpar=1.0e-6; // direction
  }
  else {
    //dpar=fModuleInitParam[paridx]*0.001;
    if (paridx<3) dpar=1.0e-4; // translations
    else dpar=1.0e-2; // angles    
  }
  AliDebug(3,Form("\n+++ automatic dpar=%g\n",dpar));
  if (fDeltaPar) dpar=fDeltaPar;
  AliDebug(3,Form("+++ using dpar=%g\n\n",dpar));
  
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

Int_t AliITSAlignMille::SetLocalEquations() {
  /// Define local equation for current track and cluster in x coor.
  /// return 0 if success
  
  // store first interaction point
  CalcIntersectionPoint(fLocalInitParam, fModuleInitParam);  
  for (Int_t i=0; i<3; i++) fPintLoc0[i]=fPintLoc[i];
  
  // calculate local derivatives numerically
  Double_t dXdL[ITSMILLE_NLOCAL],dZdL[ITSMILLE_NLOCAL];
  for (Int_t i=0; i<ITSMILLE_NLOCAL; i++) {
    if (CalcDerivatives(i,kTRUE)) return -1;
    dXdL[i]=fDerivativeXLoc;
    dZdL[i]=fDerivativeZLoc;
  }

  Double_t dXdG[ITSMILLE_NPARCH],dZdG[ITSMILLE_NPARCH];
  for (Int_t i=0; i<ITSMILLE_NPARCH; i++) {
    if (CalcDerivatives(i,kFALSE)) return -1;
    dXdG[i]=fDerivativeXLoc;
    dZdG[i]=fDerivativeZLoc;
  }

  AliDebug(2,Form("\n***************\n"));
  for (Int_t i=0; i<ITSMILLE_NLOCAL; i++)
    AliDebug(2,Form("Local parameter %d - dXdpar = %g  - dZdpar = %g\n",i,dXdL[i],dZdL[i]));
  for (Int_t i=0; i<ITSMILLE_NPARCH; i++)
    AliDebug(2,Form("Global parameter %d - dXdpar = %g  - dZdpar = %g\n",i,dXdG[i],dZdG[i]));
  AliDebug(2,Form("\n***************\n"));

  
  AliDebug(2,Form("setting local equation X with fMeas=%.6f  and fSigma=%.6f",(fMeasLoc[0]-fPintLoc0[0]), fSigmaLoc[0]));
  // set equation for Xloc coordinate
  for (Int_t i=0; i<ITSMILLE_NLOCAL; i++) 
    SetLocalDerivative(i,dXdL[i]);
  for (Int_t i=0; i<ITSMILLE_NPARCH; i++)
    SetGlobalDerivative(fCurrentModuleInternalIndex*ITSMILLE_NPARCH+i,dXdG[i]);
  fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, (fMeasLoc[0]-fPintLoc0[0]), fSigmaLoc[0]);  
  

  AliDebug(2,Form("setting local equation Z with fMeas=%.6f  and fSigma=%.6f",(fMeasLoc[2]-fPintLoc0[2]), fSigmaLoc[2]));
  // set equation for Zloc coordinate
  for (Int_t i=0; i<ITSMILLE_NLOCAL; i++) 
    SetLocalDerivative(i,dZdL[i]);
  for (Int_t i=0; i<ITSMILLE_NPARCH; i++)
    SetGlobalDerivative(fCurrentModuleInternalIndex*ITSMILLE_NPARCH+i,dZdG[i]);
  fMillepede->SetLocalEquation(fGlobalDerivatives, fLocalDerivatives, (fMeasLoc[2]-fPintLoc0[2]), fSigmaLoc[2]);  
  
  return 0;
}


void AliITSAlignMille::LocalFit(Int_t iTrack, Double_t *lTrackParam, Int_t lSingleFit) {
  /// Call local fit for this track
  if (!fIsMilleInit) {
    AliInfo("Millepede has not been initialized!");
    return;
  }
  Int_t iRes = fMillepede->LocalFit(iTrack,lTrackParam,lSingleFit);
  AliDebug(2,Form("iRes = %d",iRes));
  if (iRes && !lSingleFit) {
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

