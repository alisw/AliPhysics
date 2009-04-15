/**************************************************************************
 * Copyright(c) 2007-2010, ALICE Experiment at CERN, All rights reserved. *
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

//========================================================================
//
// This class is a helper, producing ITS aligmnent objects.
// It provides also some useful functions
// See the parameters of the misalignment at the end of this script.
//
// Main author: L. Gaudichet
// Contact: andrea.dainese@lnl.infn.it
//
//========================================================================

#include <TClonesArray.h>
#include <TMath.h>
#include <TClass.h>
#include <TGeoManager.h>
#include "AliLog.h"
#include "AliAlignObjParams.h"
#include "AliITSMisAligner.h"
#include "AliITSSurveyToAlign.h"
#include "AliMathBase.h"

ClassImp(AliITSMisAligner)
  
const Int_t AliITSMisAligner::fgkNLadders[AliITSMisAligner::kNLayers] = {20,40,14,22,34,38};
const Int_t AliITSMisAligner::fgkNDetectors[AliITSMisAligner::kNLayers] = {4,4,6,8,22,25};

const Double_t kRadToDeg = 180./TMath::Pi();

//________________________________________________________________________
AliITSMisAligner::AliITSMisAligner():
    AliMisAligner(),
    fRnd(),
    fInd(0),
    fAlignObjArray(NULL),
    fStrSPD("ITS/SPD"),
    fStrSDD("ITS/SDD"),
    fStrSSD("ITS/SSD"),
    fStrStave("/Stave"),
    fStrHalfStave("/HalfStave"),
    fStrLadder("/Ladder"),
    fStrSector("/Sector"),
    fStrSensor("/Sensor"),
    fUnifSPDSector(kFALSE),
    fUnifSPDHS(kFALSE),
    fUnifSDDLadder(kFALSE),
    fUnifSSDLadder(kFALSE),
    fUnifSPDLadder(kFALSE),
    fUnifSDDModule(kFALSE),
    fUnifSSDModule(kFALSE)
{
    //
    // defaul constructor
    //
    fRnd.SetSeed(38217945);
    fAlignObjArray = new TClonesArray("AliAlignObjParams",4000);
    for(Int_t ii=0; ii<6; ii++)
    {
	fWholeITS[ii]=0.;
	fSPDSector[ii]=0.;
	fSPDHB[ii]=0.;
	fSPDBarrel[ii]=0.;
	fSPDHS[ii]=0.;
	fSPDLadder[ii]=0.;
	fSDDLayer[ii]=0.;
	fSDDBarrel[ii]=0.;
	fSDDLadder[ii]=0.;
	fSDDModule[ii]=0.;
	fSSDBarrel[ii]=0.;
	fSSDLayer[ii]=0.;
	fSSDLadder[ii]=0.;
	fSSDModule[ii]=0.;
	fSPDLadderShiftT[ii]=0.;
	fSPDLadderShiftB[ii]=0.;
	fSDDLadderShift1[ii]=0.;
	fSDDLadderShift2[ii]=0.;
	fSSDLadderShift1[ii]=0.;
	fSSDLadderShift2[ii]=0.;
    }
}

//________________________________________________________________________
AliITSMisAligner::AliITSMisAligner(const AliITSMisAligner &mAligner):
    fRnd(mAligner.fRnd),
    fInd(0),
    fAlignObjArray(mAligner.fAlignObjArray),
    fStrSPD("ITS/SPD"),
    fStrSDD("ITS/SDD"),
    fStrSSD("ITS/SSD"),
    fStrStave("/Stave"),
    fStrHalfStave("/HalfStave"),
    fStrLadder("/Ladder"),
    fStrSector("/Sector"),
    fStrSensor("/Sensor"),
    fUnifSPDSector(kFALSE),
    fUnifSPDHS(kFALSE),
    fUnifSDDLadder(kFALSE),
    fUnifSSDLadder(kFALSE),
    fUnifSPDLadder(kFALSE),
    fUnifSDDModule(kFALSE),
    fUnifSSDModule(kFALSE)
{
    //
    // copy constructor
    //
}

//________________________________________________________________________
AliITSMisAligner &AliITSMisAligner::operator= (const AliITSMisAligner &mAligner)
{
    //
    // assignment operator
    //
    fRnd = mAligner.fRnd,
    fInd = 0;
    fAlignObjArray = mAligner.fAlignObjArray;
    fStrSPD = "ITS/SPD";
    fStrSDD = "ITS/SDD";
    fStrSSD = "ITS/SSD";
    fStrStave = "/Stave";
    fStrHalfStave = "/HalfStave";
    fStrLadder = "/Ladder";
    fStrSector = "/Sector";
    fStrSensor = "/Sensor";
    fUnifSPDSector = mAligner.fUnifSPDSector;
    fUnifSPDHS = kFALSE;
    fUnifSDDLadder = kFALSE;
    fUnifSSDLadder = kFALSE;
    fUnifSPDLadder = kFALSE;
    fUnifSDDModule = kFALSE;
    fUnifSSDModule = kFALSE;
    return (*this);
}

//_______________________________________________________________________________________
TClonesArray* AliITSMisAligner::MakeAlObjsArray() {
    // Make the array of alignment objects, depending on the misalignment scenario (Zero, Residual, Full)
    //

    // Setting default values for ideal, residual or full misalignment
    SetWholeITSMisAlignment();
    SetSPDMisAlignment();
    SetSDDMisAlignment();
    SetSSDMisAlignment(); 

    // Get array of alignment objects from survey (based on which we build later the SSD objects)
    AliITSSurveyToAlign* s2a = new AliITSSurveyToAlign();
    s2a->Run();
    TClonesArray* surveyArray = dynamic_cast<TClonesArray*> (s2a->GetAlignObjsArray());
    if(!surveyArray){
	Printf("SSD survey array was not build! Probably you missed to connect to alien");
	return 0;
    }else{
	Printf("survey array contains %d entries", surveyArray->GetEntriesFast());
    }
    //(AliGeomManager::GetGeometry())->UnlockGeometry();

    AddAlignObj("ITS",fWholeITS[0],fWholeITS[1],fWholeITS[2],fWholeITS[3],fWholeITS[4],fWholeITS[5],"fixed");

    AddSectorAlignObj(1,5,fSPDSector[0],fSPDSector[1],fSPDSector[2],fSPDSector[3],fSPDSector[4],fSPDSector[5],
	    fSPDLadderShiftT[0],fSPDLadderShiftT[1],fSPDLadderShiftT[2],fSPDLadderShiftT[3],fSPDLadderShiftT[4],fSPDLadderShiftT[5],fUnifSPDSector);
    AddSectorAlignObj(6,10,fSPDSector[0],fSPDSector[1],fSPDSector[2],fSPDSector[3],fSPDSector[4],fSPDSector[5],
	    fSPDLadderShiftB[0],fSPDLadderShiftB[1],fSPDLadderShiftB[2],fSPDLadderShiftB[3],fSPDLadderShiftB[4],fSPDLadderShiftB[5],fUnifSPDSector);

    //=****************************************
    // misalignment at the level of half-staves (SPD)/ladders (SDD,SSD) :
    //=****************************************
    AddAlignObj(0,-1,fSPDHS[0],fSPDHS[1],fSPDHS[2],fSPDHS[3],fSPDHS[4],fSPDHS[5],
	    0,0,0,0,0,0,fUnifSPDHS); // all SPD1 half-staves
    AddAlignObj(1,-1,fSPDHS[0],fSPDHS[1],fSPDHS[2],fSPDHS[3],fSPDHS[4],fSPDHS[5],
	    0,0,0,0,0,0,fUnifSPDHS); // all SPD2 half-staves

    AddAlignObj(2,-1,fSDDLadder[0],fSDDLadder[1],fSDDLadder[2],fSDDLadder[3],
	    fSDDLadder[4],fSDDLadder[5],fSDDLadderShift1[0],fSDDLadderShift1[1],fSDDLadderShift1[2],fSDDLadderShift1[3],fSDDLadderShift1[4],fSDDLadderShift1[5],fUnifSDDLadder); // all SDD1 ladders
    AddAlignObj(3,-1,fSDDLadder[0],fSDDLadder[1],fSDDLadder[2],fSDDLadder[3],
	    fSDDLadder[4],fSDDLadder[5],fSDDLadderShift2[0],fSDDLadderShift2[1],fSDDLadderShift2[2],fSDDLadderShift2[3],fSDDLadderShift2[4],fSDDLadderShift2[5],fUnifSDDLadder); // all SDD2 ladders

    // all SSD ladders
    for(Int_t ii=0; ii<surveyArray->GetEntriesFast(); ii++)
    {
	AliAlignObjParams* aop = dynamic_cast<AliAlignObjParams*> (surveyArray->UncheckedAt(ii));
//	if(!aop){
//	    Printf("Unexpected missing object at %d!", ii);
//	    continue;
//	}
	TString sName(aop->GetSymName());

	// First we shift all SSD ladders by the same quantity to reproduce a barrel shift
	ShiftAlignObj(*aop,fSSDBarrel[0],fSSDBarrel[1],fSSDBarrel[2],fSSDBarrel[3],fSSDBarrel[4],fSSDBarrel[5]);
	if(sName.Contains("SSD4") && !sName.Contains("Sensor")){
	    // we correct with the factor 1.13 for the fact that, in the inner SSD layer, z lever arm is 45.135cm instead of 51cm
	    SmearAlignObj(*aop,fSSDLadder[0],fSSDLadder[1],fSSDLadder[2],fSSDLadder[3]*1.13,fSSDLadder[4]*1.13,fSSDLadder[5]);
	    new((*fAlignObjArray)[fInd++]) AliAlignObjParams(*aop);
	}else if(sName.Contains("SSD5") && !sName.Contains("Sensor")){
	    SmearAlignObj(*aop,fSSDLadder[0],fSSDLadder[1],fSSDLadder[2],fSSDLadder[3],fSSDLadder[4],fSSDLadder[5]);
	    new((*fAlignObjArray)[fInd++]) AliAlignObjParams(*aop);
	}
	aop->ApplyToGeometry();
    }

    //=****************************************
    // misalignment at the level of ladders (SPD)/modules (SDD,SSD) :
    //=****************************************
    AddAlignObj(0,fSPDLadder[0],fSPDLadder[1],fSPDLadder[2],fSPDLadder[3],fSPDLadder[4],fSPDLadder[5],fUnifSPDLadder);// all SPD1 ladders
    AddAlignObj(1,fSPDLadder[0],fSPDLadder[1],fSPDLadder[2],fSPDLadder[3],fSPDLadder[4],fSPDLadder[5],fUnifSPDLadder);// all SPD2 ladders

    AddAlignObj(2,fSDDModule[0],fSDDModule[1],fSDDModule[2],fSDDModule[3],fSDDModule[4],fSDDModule[5],fUnifSDDModule);// all SDD1 modules
    AddAlignObj(3,fSDDModule[0],fSDDModule[1],fSDDModule[2],fSDDModule[3],fSDDModule[4],fSDDModule[5],fUnifSDDModule);// all SDD2 modules

    // all SSD modules
    for(Int_t ii=0; ii<surveyArray->GetEntriesFast(); ii++)
    {
	AliAlignObjParams* aop = dynamic_cast<AliAlignObjParams*> (surveyArray->UncheckedAt(ii));
	TString sName(aop->GetSymName());
	if(sName.Contains("SSD") && sName.Contains("Sensor"))
	{
	    SmearAlignObj(*aop,fSSDModule[0],fSSDModule[1],fSSDModule[2],fSSDModule[3],fSSDModule[4],fSSDModule[5]);
	    new((*fAlignObjArray)[fInd++]) AliAlignObjParams(*aop);
	}
    }

    return fAlignObjArray;
}

//_______________________________________________________________________________________
void AliITSMisAligner::ShiftAlignObj(AliAlignObjParams &alObj, Double_t dx, Double_t dy, Double_t dz, Double_t dpsi, Double_t dtheta, Double_t dphi)
{
    // 
    // Shift the parameters of the alignment object passed as first argument by the quantities defined by the arguments
    //
    Double_t shifts[3]; Double_t angles[3];
    alObj.GetPars(shifts, angles);
    alObj.SetPars(shifts[0]+dx, shifts[1]+dy, shifts[2]+dz, angles[0]+dpsi, angles[1]+dtheta, angles[2]+dphi);
}

//_______________________________________________________________________________________
void AliITSMisAligner::SmearAlignObj(AliAlignObjParams &alObj, Double_t sx, Double_t sy, Double_t sz, Double_t spsi, Double_t stheta, Double_t sphi)
{
    // 
    // Smear the parameters of the alignment object passed as first argument in the range defined by the subsequent arguments
    //
    Double_t shifts[3]; Double_t angles[3];
    alObj.GetLocalPars(shifts, angles);
    Double_t x = AliMathBase::TruncatedGaus(shifts[0], sx, 3.*sx);
    Double_t y = AliMathBase::TruncatedGaus(shifts[1], sy, 3.*sy);
    Double_t z = AliMathBase::TruncatedGaus(shifts[2], sz, 3.*sz);
    Double_t psi = AliMathBase::TruncatedGaus(angles[0], spsi, 3.*spsi);
    Double_t theta = AliMathBase::TruncatedGaus(angles[1], stheta, 3.*stheta);
    Double_t phi = AliMathBase::TruncatedGaus(angles[2], sphi, 3.*sphi);
    alObj.SetLocalPars(x, y, z, psi, theta, phi);
}

//_______________________________________________________________________________________
void AliITSMisAligner::SetWholeITSMisAlignment()
{
    // Set misalignment for the whole ITS for the standar scenarios (zero, residual, full)
    // To make custom misalignments set directly the misalignments at each level (methods "SetS?D*Pars")
    //
    if(TString(GetMisalType())=="ideal")
    {
	// overall ITS misalignment according to survey as reported by Werner Riegler (18/07/2008)
	SetWholeITSPars(-0.12, -0.07, 0.29, 0., 0.03, 0.04);
    }else if(TString(GetMisalType())=="residual"){
	// overall ITS misalignment according to survey as reported by Werner Riegler (18/07/2008)
	// no smearing added (would clash with vertex constraint)
	SetWholeITSPars(-0.12, -0.07, 0.29, 0., 0.03, 0.04);
    }else if(TString(GetMisalType())=="full"){
	// overall ITS misalignment according to survey as reported by Werner Riegler (18/07/2008) plus smearing
	Double_t sigmatrW = 0.01;
	Double_t sigmarotW = 0.006;
	SetWholeITSPars(AliMathBase::TruncatedGaus(-0.12,sigmatrW,3.*sigmatrW),
		AliMathBase::TruncatedGaus(-0.07,sigmatrW,3.*sigmatrW),
		AliMathBase::TruncatedGaus(0.29,sigmatrW,3.*sigmatrW),
		AliMathBase::TruncatedGaus(0.,sigmarotW,3.*sigmarotW),
		AliMathBase::TruncatedGaus(0.03,sigmarotW,3.*sigmarotW),
		AliMathBase::TruncatedGaus(0.04,sigmarotW,3.*sigmarotW));
    }
}

//_______________________________________________________________________________________
void AliITSMisAligner::SetSPDMisAlignment()
{
    // Set misalignment for SPD alignable volumes for the standar scenarios (zero, residual, full)
    // To make custom misalignments set directly the misalignments at each level (methods "SetSPD*Pars")
    //
    if(TString(GetMisalType())=="ideal")
    {
	// misalignment for SPD at all levels equal to zero (identical transformations)
	SetSPDBarrelSigmas(0., 0., 0., 0., 0., 0.); 
	SetSPDHBSigmas(0., 0., 0., 0., 0., 0.); 
	SetSPDSectorSigmas(0., 0., 0., 0., 0., 0.); 
	SetSPDHSSigmas(0., 0., 0., 0., 0., 0.); 
	SetSPDLadderSigmas(0., 0., 0., 0., 0., 0.); 
    }else if(TString(GetMisalType())=="residual"){
	// misalignment at the level of SPD barrel and half-barrels left to zero
	SetSPDBarrelSigmas(0., 0., 0., 0., 0., 0.); 
	SetSPDHBSigmas(0., 0., 0., 0., 0., 0.); 
	
	// misalignment at the level of SPD sectors (source: A.Pepato)
	SetSPDSectorSigmas( 0.0050/5., //  50 micron (~tangetial, i.e. rphi)
		0.0100/5., // 100 micron (~radial)
		0.0100/5., // 100 micron
		0.0100/30.*kRadToDeg/5., // so as to have 100 micron difference at the two extremes
		0.0100/30.*kRadToDeg/5., // so as to have 100 micron difference at the two extremes
		0.0050/1.5*kRadToDeg/5.); // so as to have 50 micron difference at the two extremes
	fUnifSPDSector=kFALSE;

	// misalignment at the level of half-staves (SPD) (source: S.Moretto)
	SetSPDHSSigmas( 0.0100/4., // 100 micron
		0.0020/4., // 20 micron
		0.0020/4., // 20 micron
		0.0020/7.*kRadToDeg/4., // so as to have 20 micron difference at the two extremes
		0.0050/7.*kRadToDeg/4., // so as to have 50 micron difference at the two extremes
		0.0050/0.7*kRadToDeg/4.);// so as to have 50 micron difference at the two extremes
	fUnifSPDHS=kFALSE;
	
	// misalignment at the level of ladders (SPD) (source: R.Santoro)
	SetSPDLadderSigmas( 0.0010/5., // 10 micron
		0.0050/5., // 50 micron
		0.0010/5., // 10 micron
		0.0001*kRadToDeg/5., // 0.1 mrad
		0.0001*kRadToDeg/5., // 0.1 mrad
		0.0001*kRadToDeg/5.);// 0.1 mrad
	fUnifSPDLadder=kFALSE;

    }else if(TString(GetMisalType())=="full"){
	// misalignment at the level of SPD barrel (source: A.Pepato)
	SetSPDBarrelSigmas( 0.1000, // 1 mm (very pessimistic)
		0.1000, // 1 mm (very pessimistic)
		0.1000, // 1 mm (very pessimistic)
		0.0500/30.*kRadToDeg, // so as to have 500 micron difference at the two extremes
		0.0500/30.*kRadToDeg, // so as to have 500 micron difference at the two extremes
		0.0500/7.*kRadToDeg); // so as to have 500 micron difference at the two extremes

	// misalignment at the level of SPD half-barrels (source: A.Pepato)
	SetSPDHBSigmas( 0.0200, // 200 micron
		0.0200, // 200 micron
		0.0200, // 200 micron
		0.0100/30.*kRadToDeg, // so as to have 100 micron difference at the two extremes
		0.0100/30.*kRadToDeg, // so as to have 100 micron difference at the two extremes
		0.0100/7.*kRadToDeg); // so as to have 100 micron difference at the two extremes

	// misalignment at the level of SPD sectors (source: A.Pepato)
	SetSPDSectorSigmas( 0.0050, //  50 micron (~tangetial, i.e. rphi)
		0.0100, // 100 micron (~radial)
		0.0100, // 100 micron
		0.0100/30.*kRadToDeg, // so as to have 100 micron difference at the two extremes
		0.0100/30.*kRadToDeg, // so as to have 100 micron difference at the two extremes
		0.0050/1.5*kRadToDeg); // so as to have 50 micron difference at the two extremes
	fUnifSPDSector=kTRUE;

	// misalignment at the level of half-staves (SPD) (source: S.Moretto)
	SetSPDHSSigmas( 0.0100, // 100 micron // normal to plane
		0.0020, // 20 micron
		0.0020, // 20 micron
		0.0020/7.*kRadToDeg, // so as to have 20 micron difference at the two extremes
		0.0050/7.*kRadToDeg, // so as to have 50 micron difference at the two extremes
		0.0050/0.7*kRadToDeg); // so as to have 50 micron difference at the two extremes
	fUnifSPDHS=kTRUE;
	
	// misalignment at the level of ladders (SPD) (source: R.Santoro)
	SetSPDLadderSigmas( 0.0010, // 10 micron
		0.0030, // 50 micron
		0.0010, // 10 micron
		0.0001*kRadToDeg, // 0.1 mrad
		0.0001*kRadToDeg, // 0.1 mrad
		0.0001*kRadToDeg); // 0.1 mrad
	fUnifSPDLadder=kTRUE;
	

	// misalignment at the level of SPD barrel, half-barrels, and at the level
	// of SPD sectors
	Double_t shBtop[6], shBbot[6]; //top and bottom barrel shifts
	for(Int_t ii=0; ii<6; ii++){
	    shBtop[ii] = AliMathBase::TruncatedGaus(0.,fSPDBarrel[ii]/3,fSPDBarrel[ii]);
	    shBbot[ii] = shBtop[ii];
	}

	for(Int_t ii=0; ii<6; ii++){
	    shBtop[ii] += AliMathBase::TruncatedGaus(0.,fSPDHB[ii]/3,fSPDHB[ii]);
	    shBbot[ii] += AliMathBase::TruncatedGaus(0.,fSPDHB[ii]/3,fSPDHB[ii]);
	}
	SetSPDLadderShiftT(shBtop);
	SetSPDLadderShiftB(shBbot);
    }
}

//_______________________________________________________________________________________
void AliITSMisAligner::SetSDDMisAlignment() 
{
    // Set misalignment for SDD alignable volumes for the standar scenarios (zero, residual, full)
    // To make custom misalignments set directly the misalignments at each level (methods "SetSDD*Pars")
    //
    if(TString(GetMisalType())=="ideal")
    {
	// misalignment for SDD at all levels equal to zero
	SetSDDLayerSigmas(0., 0., 0., 0., 0., 0.);
	SetSDDBarrelSigmas(0., 0., 0., 0., 0., 0.);
	SetSDDLadderSigmas(0., 0., 0., 0., 0., 0.);
	SetSDDModuleSigmas(0., 0., 0., 0., 0., 0.);
    }else if(TString(GetMisalType())=="residual"){
	// misalignment at the level of SDD and SSD layers(source: B.Giraudo)
	SetSDDLayerSigmas(0., 0., 0., 0., 0., 0.);
	SetSDDBarrelSigmas(0., 0., 0., 0., 0., 0.);
	
	// misalignment at the level of half-staves (SPD) (source: S.Moretto)
	SetSDDLadderSigmas( 0.0005, // 5 micron
		0.0005, // 5 micron
		0.0005, // 5 micron
		0.00, //  ?
		0.00, //  ?
		0.00); //  ?
	fUnifSDDLadder=kFALSE;

	// misalignment at the level of SDD modules(source: L.Gaudichet)
	SetSDDModuleSigmas( 0.0045/5., // 45 micron
		0.0045/5., // 45 micron
		0.0105/5., // 105 micron
		0.00, // ?
		0.00, //  ?
		0.00);//  ?
	fUnifSDDModule=kFALSE;

    }else if(TString(GetMisalType())=="full"){
	// misalignment at the level of SDD layers(source: B.Giraudo)
	SetSDDBarrelSigmas( 0.0020, // 20 micron
		0.0020, // 20 micron
		0.0020, // 20 micron
		0.0020/52.*kRadToDeg,  // so as to have 20 micron difference at the two extremes
		0.0020/52.*kRadToDeg,  // so as to have 20 micron difference at the two extremes
		0.0020/20.*kRadToDeg); // so as to have 20 micron difference at the two extremes

	SetSDDLayerSigmas( 0.0010, // 10 micron
		0.0010, // 10 micron
		0.0010, // 10 micron
		0.0010/52.*kRadToDeg, // so as to have 10 micron difference at the two extremes
		0.0010/52.*kRadToDeg, // so as to have 10 micron difference at the two extremes
		0.0010/20.*kRadToDeg);  // so as to have 10 micron difference at the two extremes

	// misalignment at the level of SDD ladders
	SetSDDLadderSigmas( 0.0005, // 5 micron
		0.0005, // 5 micron
		0.0005, // 5 micron
		0.00, //  ?
		0.00, //  ?
		0.00);//  ?
	fUnifSDDLadder=kTRUE;

	// misalignment at the level of SDD modules (source: L.Gaudichet)
	SetSDDModuleSigmas( 0.0045, // 45 micron
		0.0045, // 45 micron
		0.0105, // 105 micron
		0.00, // ?
		0.00, //  ?
		0.00);//  ?
	fUnifSDDModule=kTRUE;
    }

  fSDDLadderShift1[0] = GetUnif(-fSDDBarrel[0],fSDDBarrel[0]);
  fSDDLadderShift1[1] = GetUnif(-fSDDBarrel[1],fSDDBarrel[1]);
  fSDDLadderShift1[2] = GetUnif(-fSDDBarrel[2],fSDDBarrel[2]);
  fSDDLadderShift1[3] = GetUnif(-fSDDBarrel[3],fSDDBarrel[3]);
  fSDDLadderShift1[4] = GetUnif(-fSDDBarrel[4],fSDDBarrel[4]);
  fSDDLadderShift1[5] = GetUnif(-fSDDBarrel[5],fSDDBarrel[5]);
  
  for(Int_t ii=0; ii<6; ii++)
      fSDDLadderShift2[ii] = fSDDLadderShift1[ii];

  //  layer SDD1
  fSDDLadderShift1[0] += GetUnif(-fSDDLayer[0],fSDDLayer[0]);
  fSDDLadderShift1[1] += GetUnif(-fSDDLayer[1],fSDDLayer[1]);
  fSDDLadderShift1[2] += GetUnif(-fSDDLayer[2],fSDDLayer[2]);
  fSDDLadderShift1[3] += GetUnif(-fSDDLayer[3],fSDDLayer[3]);
  fSDDLadderShift1[4] += GetUnif(-fSDDLayer[4],fSDDLayer[4]);
  fSDDLadderShift1[5] += GetUnif(-fSDDLayer[5],fSDDLayer[5]);

  //  layer SDD2
  fSDDLadderShift2[0] += GetUnif(-fSDDLayer[0],fSDDLayer[0]);
  fSDDLadderShift2[1] += GetUnif(-fSDDLayer[1],fSDDLayer[1]);
  fSDDLadderShift2[2] += GetUnif(-fSDDLayer[2],fSDDLayer[2]);
  fSDDLadderShift2[3] += GetUnif(-fSDDLayer[3],fSDDLayer[3]);
  fSDDLadderShift2[4] += GetUnif(-fSDDLayer[4],fSDDLayer[4]);
  fSDDLadderShift2[5] += GetUnif(-fSDDLayer[5],fSDDLayer[5]);

}
	
//_______________________________________________________________________________________
void AliITSMisAligner::SetSSDMisAlignment() 
{
    // Set misalignment for SSD alignable volumes for the standar scenarios (zero, residual, full)
    // To make custom misalignments set directly the misalignments at each level (methods "SetSSD*Pars")
    //
    if(TString(GetMisalType())=="ideal"){

	// zero misalignment at the level of SSD barrel
	SetSSDBarrelPars(0.,0.,0.,0.,0.,0.);
	// zero misalignment at the level of SSD ladders
	SetSSDLadderSigmas(0.,0.,0.,0.,0.,0.);
	// zero misalignment at the level of SSD modules
	SetSSDModuleSigmas(0.,0.,0.,0.,0.,0.);
	
    }else if(TString(GetMisalType())=="residual"){

	// zero misalignment at the level of SSD barrel
	SetSSDBarrelPars(0.,0.,0.,0.,0.,0.);
	// misalignment at the level of SSD ladders (source: M. Van Leeuwen)
	// values set so that overall maximum displacement (combined effect of shift and rotation
	// of the ladder) for any point of the ladder cannot exceed 10um in x, 100 um in y, 50um in z
	SetSSDLadderSigmas( 0.0005, // 5 microns
		0.0033, // 33 microns
		0.0050, // 50 microns
		0.000067*kRadToDeg, // 0.067 mrads
		0.00001*kRadToDeg, // 0.01 mrads
		0.001*kRadToDeg); // 1 mrad
	fUnifSSDLadder=kTRUE;

	// misalignment at the level of SSD modules (source: M. Van Leeuwen)
	// values set so that overall maximum displacement (combined effect of shift and rotation
	// of the ladder) for any point of the module cannot exceed 5um in x, 10 um in y, 5um in z
	SetSSDModuleSigmas( 0.00025, // 2.5 microns
		0.00034, // 3.4 microns
		0.0005, // 5 microns
		0.00017*kRadToDeg, // 0.17 mrads
		0.000125*kRadToDeg, // 0.125 mrads
		0.0001*kRadToDeg); // 0.1 mrads
	fUnifSSDModule=kTRUE;

    }else if(TString(GetMisalType())=="full"){
	// misalignment at the level of SSD layers (source: B. Giraudo)
	SetSSDBarrelPars( GetUnif(-0.0020,0.0020), // 20 micron
		GetUnif(-0.0020,0.0020), // 20 micron
		GetUnif(-0.0020,0.0020), // 20 micron
		GetUnif(-0.0020/90.*kRadToDeg,0.0020), // so as to have 20 micron difference at the two extremes
		GetUnif(-0.0020/90.*kRadToDeg,0.0020), // so as to have 20 micron difference at the two extremes
		GetUnif(-0.0020/40.*kRadToDeg,0.0020));  // so as to have 20 micron difference at the two extremes
	
	// misalignment at the level of SSD ladders (source: M. Van Leeuwen)
	// values set so that overall maximum displacement (combined effect of shift and rotation
	// of the ladder) for any point of the ladder cannot exceed 20um in x, 100 um in y, 50um in z
	SetSSDLadderSigmas( 0.0010, // 10 microns
		0.0033, // 33 microns
		0.0050, // 50 microns
		0.000067*kRadToDeg, // 0.067 mrads
		0.00002*kRadToDeg, // 0.02 mrads
		0.001*kRadToDeg); // 1 mrad
	fUnifSSDLadder=kTRUE;

	// misalignment at the level of SSD modules (source: M. Van Leeuwen)
	// values set so that overall maximum displacement (combined effect of shift and rotation
	// of the ladder) for any point of the module cannot exceed 5um in x, 10 um in y, 5um in z
	SetSSDModuleSigmas( 0.00025, // 2.5 microns
		0.00034, // 3.4 microns
		0.0005, // 5 microns
		0.00017*kRadToDeg, // 0.17 mrads
		0.000125*kRadToDeg, // 0.125 mrads
		0.0001*kRadToDeg); // 0.1 mrads
	fUnifSSDModule=kTRUE;
    }


  fSSDLadderShift1[0] = GetUnif(-fSSDBarrel[0],fSSDBarrel[0]);
  fSSDLadderShift1[1] = GetUnif(-fSSDBarrel[1],fSSDBarrel[1]);
  fSSDLadderShift1[2] = GetUnif(-fSSDBarrel[2],fSSDBarrel[2]);
  fSSDLadderShift1[3] = GetUnif(-fSSDBarrel[3],fSSDBarrel[3]);
  fSSDLadderShift1[4] = GetUnif(-fSSDBarrel[4],fSSDBarrel[4]);
  fSSDLadderShift1[5] = GetUnif(-fSSDBarrel[5],fSSDBarrel[5]);
  
  /*
  for(Int_t ii=0; ii<6; ii++)
      fSSDLadderShift2[ii] = fSSDLadderShift1[ii];

  //  layer SSD1
  fSSDLadderShift1[0] += GetUnif(-fSSDLayer[0],fSSDLayer[0]);
  fSSDLadderShift1[1] += GetUnif(-fSSDLayer[1],fSSDLayer[1]);
  fSSDLadderShift1[2] += GetUnif(-fSSDLayer[2],fSSDLayer[2]);
  fSSDLadderShift1[3] += GetUnif(-fSSDLayer[3],fSSDLayer[3]);
  fSSDLadderShift1[4] += GetUnif(-fSSDLayer[4],fSSDLayer[4]);
  fSSDLadderShift1[5] += GetUnif(-fSSDLayer[5],fSSDLayer[5]);

  //  layer SSD2
  fSSDLadderShift2[0] += GetUnif(-fSSDLayer[0],fSSDLayer[0]);
  fSSDLadderShift2[1] += GetUnif(-fSSDLayer[1],fSSDLayer[1]);
  fSSDLadderShift2[2] += GetUnif(-fSSDLayer[2],fSSDLayer[2]);
  fSSDLadderShift2[3] += GetUnif(-fSSDLayer[3],fSSDLayer[3]);
  fSSDLadderShift2[4] += GetUnif(-fSSDLayer[4],fSSDLayer[4]);
  fSSDLadderShift2[5] += GetUnif(-fSSDLayer[5],fSSDLayer[5]);
  */

}

//_______________________________________________________________________________________
AliCDBMetaData* AliITSMisAligner::GetCDBMetaData() const {
    // Returns the AliCDBMetaData to be associated with AliCDBEntry which will contain
    // the array of alignment objects for ITS misalignment
    // Presently "responsible" and "comment" are filled.
    //
  AliCDBMetaData* md = new AliCDBMetaData();
  md->SetResponsible("Andrea Dainese");

  if(TString(GetMisalType())=="ideal")
    md->SetComment("Alignment objects for ITS ideal misalignment");
  if(TString(GetMisalType())=="residual")
    md->SetComment("Alignment objects for ITS residual misalignment");
  if(TString(GetMisalType())=="full")
    md->SetComment("Alignment objects for ITS full misalignment");
  return md;
}

//________________________________________________________________________
Bool_t AliITSMisAligner::AddAlignObj(char* name,Double_t dx,Double_t dy,Double_t dz,
				Double_t dpsi,Double_t dtheta,Double_t dphi,const char* distrib) {
  //
  // misalignment by symname
  //
  Double_t vx=0.,vy=0.,vz=0.,vpsi=0.,vtheta=0.,vphi=0.;

  TString sdistrib(distrib);

  if(sdistrib==TString("gaussian")) {
    vx = AliMathBase::TruncatedGaus(0.,dx/3.,dx); // mean, sigma, max absolute value 
    vy = AliMathBase::TruncatedGaus(0.,dy/3.,dy);
    vz = AliMathBase::TruncatedGaus(0.,dz/3.,dz);
    vpsi   = AliMathBase::TruncatedGaus(0.,dpsi/3.,  dpsi );
    vtheta = AliMathBase::TruncatedGaus(0.,dtheta/3.,dtheta);
    vphi   = AliMathBase::TruncatedGaus(0.,dphi/3.,  dphi);
  }else if(sdistrib==TString("uniform")){ 
    vx = fRnd.Uniform(-dx,dx);
    vy = fRnd.Uniform(-dy,dy);
    vz = fRnd.Uniform(-dz,dz);
    vpsi = fRnd.Uniform(-dpsi,dpsi);
    vtheta = fRnd.Uniform(-dtheta,dtheta);
    vphi = fRnd.Uniform(-dphi,dphi);
  }else if(sdistrib==TString("fixed")){
    vx=dx;
    vy=dy;
    vz=dz;
    vpsi=dpsi;
    vtheta=dtheta;
    vphi=dphi;
  }else{
    AliFatal(Form("Invalid string \"%s\" specifying the misalignment type for the volume \"%s\""));
  }

  new((*fAlignObjArray)[fInd]) AliAlignObjParams(name,0,vx,vy,vz,vpsi,vtheta,vphi,kFALSE);

  AliAlignObjParams* itsalobj = (AliAlignObjParams*) fAlignObjArray->UncheckedAt(fInd);
  itsalobj->ApplyToGeometry();

  fInd++;

  return kTRUE;
}


//________________________________________________________________________
Bool_t AliITSMisAligner::AddAlignObj(Int_t lay,Double_t dx,Double_t dy,Double_t dz,
				 Double_t dpsi,Double_t dtheta,Double_t dphi,Bool_t unif) {
  //
  // misalignment at the level of sensitive alignable volumes (SPD ladders/ SDD,SSD modules)
  // done for all ladders/modules of the given layer
  //
  lay+=1; // layers are numbered from 1 to 6 in AliGeomManager

  printf("LAYER %d  MODULES %d\n",lay,AliGeomManager::LayerSize(lay));

  for (Int_t iModule = 0; iModule < AliGeomManager::LayerSize(lay); iModule++) {

    Double_t vx,vy,vz,vpsi,vtheta,vphi;
    
    if(!unif) {
      vx = AliMathBase::TruncatedGaus(0.,dx/3.,dx); // mean, sigma, max absolute value 
      vy = AliMathBase::TruncatedGaus(0.,dy/3.,dy);
      vz = AliMathBase::TruncatedGaus(0.,dz/3.,dz);
      vpsi = AliMathBase::TruncatedGaus(0.,dpsi/3.,dpsi);
      vtheta = AliMathBase::TruncatedGaus(0.,dtheta/3.,dtheta);
      vphi = AliMathBase::TruncatedGaus(0.,dphi/3.,dphi);
    } else {
      vx = fRnd.Uniform(-dx,dx);
      vy = fRnd.Uniform(-dy,dy);
      vz = fRnd.Uniform(-dz,dz);
      vpsi = fRnd.Uniform(-dpsi,dpsi);
      vtheta = fRnd.Uniform(-dtheta,dtheta);
      vphi = fRnd.Uniform(-dphi,dphi);
    }
    
    UShort_t volid = AliGeomManager::LayerToVolUID(lay,iModule);
    const char *symname = AliGeomManager::SymName(volid);
    
    new((*fAlignObjArray)[fInd]) AliAlignObjParams(symname,volid,vx,vy,vz,vpsi,vtheta,vphi,kFALSE);
    fInd++; 
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliITSMisAligner::AddAlignObj(Int_t lay,Int_t ladd,Double_t dx,Double_t dy,Double_t dz,
				       Double_t dpsi,Double_t dtheta,Double_t dphi,
				       Double_t xShift,Double_t yShift,Double_t zShift,
				       Double_t psiShift,Double_t thetaShift,Double_t phiShift,
				       Bool_t unif) {
  //
  // misalignment at the level of half-staves/ladders (ladd=-1 means that all ladders are scanned)
  //
  Double_t vx,vy,vz,vpsi,vtheta,vphi;
  Double_t tr[3],rot[3];  
  
  Int_t laddMin = ladd;
  Int_t laddMax = laddMin+1;
  if (ladd<0) {
    laddMin = 0;
    laddMax = fgkNLadders[lay];
  }

  for (Int_t iLadd=laddMin; iLadd<laddMax; iLadd++) {

    Int_t nHS = 1; 
    if (lay<2) nHS = 2;
    for (Int_t iHalfStave=0; iHalfStave<nHS; iHalfStave++) {
      
      if(!unif) {
	vx = AliMathBase::TruncatedGaus(0.,dx/3.,dx); // mean, sigma, max absolute value 
	vy = AliMathBase::TruncatedGaus(0.,dy/3.,dy);
	vz = AliMathBase::TruncatedGaus(0.,dz/3.,dz);
	vpsi = AliMathBase::TruncatedGaus(0.,dpsi/3.,dpsi);
	vtheta = AliMathBase::TruncatedGaus(0.,dtheta/3.,dtheta);
	vphi = AliMathBase::TruncatedGaus(0.,dphi/3.,dphi);
      } else {
	vx = fRnd.Uniform(-dx,dx);
	vy = fRnd.Uniform(-dy,dy);
	vz = fRnd.Uniform(-dz,dz);
	vpsi = fRnd.Uniform(-dpsi,dpsi);
	vtheta = fRnd.Uniform(-dtheta,dtheta);
	vphi = fRnd.Uniform(-dphi,dphi);
      }

      TString name(GetHalfStaveLadderSymbName(lay,iLadd,iHalfStave));

      // first apply half-stave / ladder level misalignment
      AliAlignObjParams aaop(name.Data(),0,vx,vy,vz,vpsi,vtheta,vphi,kFALSE); // set them as local
      aaop.GetPars(tr,rot); // global

      // then, apply layer-level misalignment (only for SDD and SSD)
      if(lay>1) {
	tr[0] += xShift;
	tr[1] += yShift;
	tr[2] += zShift;
	rot[0] += psiShift;
	rot[1] += thetaShift;
	rot[2] += phiShift;
      }
      new((*fAlignObjArray)[fInd]) AliAlignObjParams(name.Data(),0,tr[0],tr[1],tr[2],rot[0],rot[1],rot[2],kTRUE); // set them as global

      AliAlignObjParams* itsalobj = (AliAlignObjParams*) fAlignObjArray->UncheckedAt(fInd);
      itsalobj->ApplyToGeometry();
      fInd++;
    }
  }
  
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliITSMisAligner::AddSectorAlignObj(Int_t sectMin,Int_t sectMax,
				       Double_t dx,Double_t dy,Double_t dz,
				       Double_t dpsi,Double_t dtheta,Double_t dphi,
				       Double_t xShift,Double_t yShift,Double_t zShift,
				       Double_t psiShift,Double_t thetaShift,Double_t phiShift,Bool_t unif) {
  //
  // misalignment at the level of SPD sectors and half-barrels
  // 

  if ((sectMin<1) || (sectMax>10)) return kFALSE;
  Double_t vx,vy,vz,vpsi,vtheta,vphi;
  Double_t tr[3],rot[3];  

  for (Int_t iSect = sectMin-1; iSect<sectMax; iSect++) {

    // first, apply sector level misalignment    
    if(!unif) {
	vx = AliMathBase::TruncatedGaus(0.,dx/3.,dx); // mean, sigma, max absolute value 
	vy = AliMathBase::TruncatedGaus(0.,dy/3.,dy);
	vz = AliMathBase::TruncatedGaus(0.,dz/3.,dz);
	vpsi = AliMathBase::TruncatedGaus(0.,dpsi/3.,dpsi);
	vtheta = AliMathBase::TruncatedGaus(0.,dtheta/3.,dtheta);
	vphi = AliMathBase::TruncatedGaus(0.,dphi/3.,dphi);
    } else {
      vx = fRnd.Uniform(-dx,dx);
      vy = fRnd.Uniform(-dy,dy);
      vz = fRnd.Uniform(-dz,dz);
      vpsi = fRnd.Uniform(-dpsi,dpsi);
      vtheta = fRnd.Uniform(-dtheta,dtheta);
      vphi = fRnd.Uniform(-dphi,dphi);
    }

    TString name(GetSymbName(0));
    name += fStrSector;
    name += iSect;


    AliAlignObjParams aaop(name.Data(),0,vx,vy,vz,vpsi,vtheta,vphi,kFALSE); // set them as local
    aaop.GetPars(tr,rot); // global

    // then, apply half-barrel level misalignment
    tr[0] += xShift;
    tr[1] += yShift;
    tr[2] += zShift;
    rot[0] += psiShift;
    rot[1] += thetaShift;
    rot[2] += phiShift;

    new((*fAlignObjArray)[fInd]) AliAlignObjParams(name.Data(),0,tr[0],tr[1],tr[2],rot[0],rot[1],rot[2],kTRUE); // set them as global

    AliAlignObjParams* itsalobj = (AliAlignObjParams*) fAlignObjArray->UncheckedAt(fInd);
    itsalobj->ApplyToGeometry();
    fInd++;
  }
  return kTRUE;
}

//________________________________________________________________________
const char* AliITSMisAligner::GetSymbName(Int_t layer) const {
  //
  // be careful : SPD0 and SPD1 are not physically separated 
  //
  TString name;
  switch (layer) {
  case 0:
  case 1: name = fStrSPD; name += layer; break;
  case 2:
  case 3: name = fStrSDD; name += layer; break;
  case 4:
  case 5: name = fStrSSD; name += layer; break;
  default: AliFatal("Wrong layer index");
  }
  return name.Data();
}

//________________________________________________________________________
const char* AliITSMisAligner::GetSymbName(Int_t layer, Int_t ladder, Int_t det) const {
  //
  // symname from layer, ladder, detector
  //
  TString symname(GetHalfStaveLadderSymbName(layer,ladder,det));
  if(layer<=2){
    symname+="Ladder";
  }else if(layer<=6){
    symname+="Sensor";
  }else{
    AliError("Invalid layer!");
    return 0;
  }
  symname+=det;
  return symname.Data();
}

//________________________________________________________________________
const char* AliITSMisAligner::GetSymbName(Int_t layer,Int_t ladd) const {
  //
  // Get logical names at the level of staves / ladders
  //
  TString name(GetSymbName(layer));
  if (layer==0) { // SPD1

    int sector = ladd/2;
    name += fStrSector;
    name += sector;
    int stave = ladd-sector*2;
    name += fStrStave;
    name += stave;
  }
  else if (layer==1) { // SPD2

    int sector = ladd/4;
    name += fStrSector;
    name += sector;
    int stave = ladd-sector*4;
    name += fStrStave;
    name += stave;
  }
  else if (layer>=2 && layer<=5) { // SDD and SSD
    name += fStrLadder;
    name += ladd;
  }
  else {
    AliFatal("Wrong layer index");
  }
  return name.Data();
}

//________________________________________________________________________
const char* AliITSMisAligner::GetHalfStaveLadderSymbName(Int_t layer,Int_t ladd,Int_t halfStave) const {
  //
  // Get logical names at the level of half-staves (SPD) or ladders (SDD and SSD)
  //
  TString name(GetSymbName(layer));
  if (layer==0) { // SPD1

    int sector = ladd/2;
    name += fStrSector;
    name += sector;
    int stave = ladd-sector*2;
    name += fStrStave;
    name += stave;
    name += fStrHalfStave;
    name += halfStave;
  }
  else if (layer==1) { // SPD2

    int sector = ladd/4;
    name += fStrSector;
    name += sector;
    int stave = ladd-sector*4;
    name += fStrStave;
    name += stave;
    name += fStrHalfStave;
    name += halfStave;
  }
  else if (layer>=2 && layer<=5) { // SDD and SSD
    name += fStrLadder;
    name += ladd;
  } 
  else {
    AliFatal("Wrong layer index");
  }
  return name.Data();
}

//________________________________________________________________________
const char* AliITSMisAligner::GetParentSymName(const char* symname) {
  //
  // symnane of parent volume
  //
  TString parent(symname);
  // Give the symname of 
  if(parent.BeginsWith('/')) parent.Remove(TString::kLeading,'/');
  if(parent.EndsWith("/")) parent.Remove(TString::kTrailing,'/');
  
  if(!parent.CountChar('/')) AliErrorClass("Not a valid symbolic name");

  Int_t layer,level;
  GetLayerAndLevel(symname,layer,level);
  if(level==1) return "ITS";
  
  parent.Remove(parent.Last('/'));
  
  if((layer==0 || layer==1) && level==2){
    parent.Remove(parent.Last('/'));
    parent[7]='0';
  }
    
  return parent.Data(); 
}

//________________________________________________________________________
Bool_t AliITSMisAligner::GetLayerAndLevel(const char* symname, Int_t &layer, Int_t &level) {
  //
  // given the symbolic name set layer and level
  //
  const char* basename[6] = {"ITS/SPD0/Sector","ITS/SPD1/Sector","ITS/SDD2/Ladder","ITS/SDD3/Ladder","ITS/SSD4/Ladder","ITS/SSD5/Ladder"};
  TString strSym(symname);
  if(strSym=="ITS"){
    level=0;
    layer=-1;
    return kTRUE;
  }
  Int_t i;
  for(i=0; i<6; i++){
    if(strSym.BeginsWith(basename[i])) break;
  }

  if(i>=6){
    AliErrorClass(Form("%s is not a valid symbolic name for an ITS alignable volume",strSym.Data()));
    return kFALSE;
  }
  
  layer=i;
  //The part above could be replaced by just
  // TString seventh = strSym[7];
  // layer = seventh.Atoi();
  // if we don't need to check the validity of the symname
  
  level=1;
  switch(layer){
    case 0:
    case 1:
      if(strSym.Contains("Stave")) level=2;
      if(strSym.Contains("Ladder")) level=3;
      break;
    case 2:
    case 3:
    case 4:
    case 5:
      if(strSym.Contains("Sensor")) level=2;
  }
  
  return kTRUE;
}

//________________________________________________________________________
Int_t AliITSMisAligner::GetNSisters(const char* symname) {
  //
  // number of volumes on same level
  //
  Int_t layer,level;
  if(!GetLayerAndLevel(symname,layer,level)) return -1;
  if(level==0) return -1;
  if(level==1) return GetNLadders(layer);
  if(level==2) return GetNDetectors(layer);
  AliErrorClass(Form("Invalid layer and level"));
  return -1;
}

//________________________________________________________________________
Int_t AliITSMisAligner::GetNDaughters(const char* symname) {
  //
  // number of daughter volumes
  // 
  Int_t layer,level;
  if(!GetLayerAndLevel(symname,layer,level)) return -1;
  if(level==0) {
    Int_t nLadders = 0;
    for(Int_t lay=0; lay<6; lay++) nLadders += GetNLadders(lay);
    return nLadders;
  }
  if(level==1) return GetNDetectors(layer);
  if(level==2){
    AliWarningClass(Form("Volume %s is a sensitive volume and has no alignable dauthers",symname));
    return -1;
  }
  AliErrorClass(Form("Invalid layer and level"));
  return -1;
}

/*
//________________________________________________________________________
TString AliITSMisAligner::GetSymbName(Int_t layer,Int_t ladd,Int_t mod) const {

  // Get logical names at the level of SPD ladders / SDD and SSD modules

  Int_t halfStave = mod/2;
  TString name = GetHalfStaveLadderSymbName(layer,ladd,halfStave);

  if (layer<2) { // SPD
    name += fStrLadder;
    name += mod;
  } 
  else if (layer>=2 && layer<=5) { // SDD and SSD
    name += fStrSensor;
    name += mod;
  }
  else {
    AliFatal("Wrong layer index");
  }
  return name;
}
*/

