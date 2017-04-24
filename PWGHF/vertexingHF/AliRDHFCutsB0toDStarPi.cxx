/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliRDHFCutsB0toDStarPi.cxx $ */

/////////////////////////////////////////////////////////////
//
// Class for cuts on AOD reconstructed B0->DStarPi->D0PiPi->KPiPiPi
//
// Author: Lennart van Doremalen, l.v.r.vandoremalen@uu.nl  
//
// Based on work by A.Grelli, alessandro.grelli@uu.nl
// PID method implemented by   Y.Wang, yifei@physi.uni-heidelberg.de
//           
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsB0toDStarPi.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCutsB0toDStarPi);
/// \endcond


//--------------------------------------------------------------------------
AliRDHFCutsB0toDStarPi::AliRDHFCutsB0toDStarPi(const char* name) : 
  AliRDHFCuts(name),
  fMaxPtPid(9999.),
  fTPCflag(999.),
  fCircRadius(0.),
  fIsCutUsed(0),   
  fCutsRDD0forD0ptbin(0),
  fnPtBinsD0forD0ptbin(1),
  fnPtBinLimitsD0forD0ptbin(1),
  fPtBinLimitsD0forD0ptbin(0),
  fIsUpperCutD0forD0ptbin(0),
  fIsCutUsedD0forD0ptbin(0),
  fVarNamesD0forD0ptbin(0),
  fCutsRDD0forDStarptbin(0),
  fnPtBinsD0forDStarptbin(1),
  fnPtBinLimitsD0forDStarptbin(1),
  fPtBinLimitsD0forDStarptbin(0),
  fIsUpperCutD0forDStarptbin(0),
  fIsCutUsedD0forDStarptbin(0),
  fVarNamesD0forDStarptbin(0),
  fCutsRDDStarforDStarptbin(0),
  fnPtBinsDStarforDStarptbin(1),
  fnPtBinLimitsDStarforDStarptbin(1),
  fPtBinLimitsDStarforDStarptbin(0),
  fIsUpperCutDStarforDStarptbin(0),
  fIsCutUsedDStarforDStarptbin(0),
  fVarNamesDStarforDStarptbin(0)  
{
  //
  // Default Constructor
  //
  
  // Main cut setup as function of B0 pt bins
  const Int_t nvars=85;
  SetNVars(nvars);

  TString varNames[nvars];
  Int_t iterator = 0;

  //D0 cut variables
  varNames[iterator++]=   /*-00-*/ "inv. mass width[GeV]";    
  varNames[iterator++]=   /*-01-*/ "delta mass width  [GeV]";    //not used for D0
  varNames[iterator++]=   /*-02-*/ "pointing angle [Cos(theta)]";
  varNames[iterator++]=   /*-03-*/ "dca [cm]";                   
  varNames[iterator++]=   /*-04-*/ "Pt D0 [GeV/c]";
  varNames[iterator++]=   /*-05-*/ "Pt Kaon [GeV/c]";
  varNames[iterator++]=   /*-06-*/ "Pt Pion [GeV/c]";
  varNames[iterator++]=   /*-07-*/ "d0 D0 [cm]";
  varNames[iterator++]=   /*-08-*/ "d0 Kaon [cm]";                
  varNames[iterator++]=   /*-09-*/ "d0 Pion [cm]";                
  varNames[iterator++]=   /*-10-*/ "d0d0 [cm^2]";
  varNames[iterator++]=   /*-11-*/ "d0d0 XY [cm^2]";

  varNames[iterator++]=   /*-12-*/ "angle between both daughters"; 
  varNames[iterator++]=   /*-13-*/ "angle mother with first daughter";
  varNames[iterator++]=   /*-14-*/ "angle mother with second daughter";
  varNames[iterator++]=   /*-15-*/ "cosThetaStar";                
  varNames[iterator++]=   /*-16-*/ "vertexDistance"; 
  varNames[iterator++]=   /*-17-*/ "pseudoProperDecayTime"; 
  varNames[iterator++]=   /*-18-*/ "DecayTime"; 
  varNames[iterator++]=   /*-19-*/ "normalizedDecayTime";
  varNames[iterator++]=   /*-20-*/ "normDecayLength";
  varNames[iterator++]=   /*-21-*/ "topomatic first daughter";
  varNames[iterator++]=   /*-22-*/ "topomatic second daughter";
  varNames[iterator++]=   /*-23-*/ "topomatic max";
  varNames[iterator++]=   /*-24-*/ "topomatic min";


  varNames[iterator++]=   /*-25-*/ "pointingAngleToDStar";
  varNames[iterator++]=   /*-26-*/ "d0MotherToDStar"; 
  varNames[iterator++]=   /*-27-*/ "d0FirstDaughterToDStar"; 
  varNames[iterator++]=   /*-28-*/ "d0SecondDaughterToDStar"; 
  varNames[iterator++]=   /*-29-*/ "impactProductToDStar"; 
  varNames[iterator++]=   /*-30-*/ "impactProductXYToDStar"; 
  varNames[iterator++]=   /*-31-*/ "normDecayLengthToDStar"; 
  varNames[iterator++]=   /*-32-*/ "pseudoProperDecayTimeToDStar"; 
  varNames[iterator++]=   /*-33-*/ "DecayTimeToDStar"; 
  varNames[iterator++]=   /*-34-*/ "normalizedDecayTimeToDStar"; 


  //DStar cut variables
  varNames[iterator++]=   /*-35-*/ "inv. mass width[GeV]";    
  varNames[iterator++]=   /*-36-*/ "delta mass width  [GeV]";
  varNames[iterator++]=   /*-37-*/ "pointing angle [Cos(theta)]";
  varNames[iterator++]=   /*-38-*/ "dca [cm]";                   
  varNames[iterator++]=   /*-39-*/ "Pt DStar [GeV/c]";
  varNames[iterator++]=   /*-40-*/ "Pt D0 [GeV/c]";
  varNames[iterator++]=   /*-41-*/ "Pt Pion [GeV/c]";
  varNames[iterator++]=   /*-42-*/ "d0 DStar [cm]";
  varNames[iterator++]=   /*-43-*/ "d0 D0 [cm]";                
  varNames[iterator++]=   /*-44-*/ "d0 Pion [cm]";                
  varNames[iterator++]=   /*-45-*/ "d0d0 [cm^2]";
  varNames[iterator++]=   /*-46-*/ "d0d0 XY [cm^2]";

  varNames[iterator++]=   /*-47-*/ "angle between both daughters"; 
  varNames[iterator++]=   /*-48-*/ "angle mother with first daughter";
  varNames[iterator++]=   /*-49-*/ "angle mother with second daughter";
  varNames[iterator++]=   /*-50-*/ "cosThetaStar";                
  varNames[iterator++]=   /*-51-*/ "vertexDistance"; 
  varNames[iterator++]=   /*-52-*/ "pseudoProperDecayTime"; 
  varNames[iterator++]=   /*-53-*/ "DecayTime"; 
  varNames[iterator++]=   /*-54-*/ "normalizedDecayTime"; 
  varNames[iterator++]=   /*-55-*/ "normDecayLength";
  varNames[iterator++]=   /*-56-*/ "topomatic first daughter";
  varNames[iterator++]=   /*-57-*/ "topomatic second daughter";
  varNames[iterator++]=   /*-58-*/ "topomatic max";
  varNames[iterator++]=   /*-59-*/ "topomatic min";


  //B0 cut variables
  varNames[iterator++]=   /*-60-*/ "inv. mass width[GeV]";    
  varNames[iterator++]=   /*-61-*/ "delta mass width  [GeV]"; 
  varNames[iterator++]=   /*-62-*/ "pointing angle [Cos(theta)]";
  varNames[iterator++]=   /*-63-*/ "dca [cm]";                   
  varNames[iterator++]=   /*-64-*/ "Pt B0 [GeV/c]";
  varNames[iterator++]=   /*-65-*/ "Pt DStar [GeV/c]";
  varNames[iterator++]=   /*-66-*/ "Pt Pion [GeV/c]";
  varNames[iterator++]=   /*-67-*/ "d0 B0 [cm]";
  varNames[iterator++]=   /*-68-*/ "d0 DStar [cm]";                
  varNames[iterator++]=   /*-69-*/ "d0 Pion [cm]";                
  varNames[iterator++]=   /*-70-*/ "d0d0 [cm^2]";
  varNames[iterator++]=   /*-71-*/ "d0d0 XY [cm^2]";

  varNames[iterator++]=   /*-72-*/ "angle between both daughters"; 
  varNames[iterator++]=   /*-73-*/ "angle mother with first daughter";
  varNames[iterator++]=   /*-74-*/ "angle mother with second daughter";
  varNames[iterator++]=   /*-75-*/ "cosThetaStar";                
  varNames[iterator++]=   /*-76-*/ "vertexDistance"; 
  varNames[iterator++]=   /*-77-*/ "pseudoProperDecayTime"; 
  varNames[iterator++]=   /*-78-*/ "DecayTime"; 
  varNames[iterator++]=   /*-79-*/ "normalizedDecayTime";
  varNames[iterator++]=   /*-80-*/ "normDecayLength";
  varNames[iterator++]=   /*-81-*/ "topomatic first daughter";
  varNames[iterator++]=   /*-82-*/ "topomatic second daughter";
  varNames[iterator++]=   /*-83-*/ "topomatic max";
  varNames[iterator++]=   /*-84-*/ "topomatic min"; 

  Bool_t isUpperCut[nvars]={0};

  SetVarNames(nvars,varNames,isUpperCut);

  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);


  //
  // Initialization of D0 cuts for D0 pt bins
  //

  const Int_t nvarsD0forD0ptbin=25;
  SetNVarsD0forD0ptbin(nvarsD0forD0ptbin);

  TString varNamesD0forD0ptbin[nvarsD0forD0ptbin];
  iterator = 0;

  //D0 cut variables
  varNamesD0forD0ptbin[iterator++]=   /*-00-*/ "inv. mass width[GeV]";    
  varNamesD0forD0ptbin[iterator++]=   /*-01-*/ "delta mass width  [GeV]";    //not used for D0
  varNamesD0forD0ptbin[iterator++]=   /*-02-*/ "pointing angle [Cos(theta)]";
  varNamesD0forD0ptbin[iterator++]=   /*-03-*/ "dca [cm]";                   
  varNamesD0forD0ptbin[iterator++]=   /*-04-*/ "Pt D0 [GeV/c]";
  varNamesD0forD0ptbin[iterator++]=   /*-05-*/ "Pt Kaon [GeV/c]";
  varNamesD0forD0ptbin[iterator++]=   /*-06-*/ "Pt Pion [GeV/c]";
  varNamesD0forD0ptbin[iterator++]=   /*-07-*/ "d0 D0 [cm]";
  varNamesD0forD0ptbin[iterator++]=   /*-08-*/ "d0 Kaon [cm]";                
  varNamesD0forD0ptbin[iterator++]=   /*-09-*/ "d0 Pion [cm]";                
  varNamesD0forD0ptbin[iterator++]=   /*-10-*/ "d0d0 [cm^2]";
  varNamesD0forD0ptbin[iterator++]=   /*-11-*/ "d0d0 XY [cm^2]";

  varNamesD0forD0ptbin[iterator++]=   /*-12-*/ "angle between both daughters"; 
  varNamesD0forD0ptbin[iterator++]=   /*-13-*/ "angle mother with first daughter";
  varNamesD0forD0ptbin[iterator++]=   /*-14-*/ "angle mother with second daughter";
  varNamesD0forD0ptbin[iterator++]=   /*-15-*/ "cosThetaStar";                
  varNamesD0forD0ptbin[iterator++]=   /*-16-*/ "vertexDistance"; 
  varNamesD0forD0ptbin[iterator++]=   /*-17-*/ "pseudoProperDecayTime"; 
  varNamesD0forD0ptbin[iterator++]=   /*-18-*/ "DecayTime"; 
  varNamesD0forD0ptbin[iterator++]=   /*-19-*/ "normalizedDecayTime";
  varNamesD0forD0ptbin[iterator++]=   /*-20-*/ "normDecayLength";
  varNamesD0forD0ptbin[iterator++]=   /*-21-*/ "topomatic first daughter";
  varNamesD0forD0ptbin[iterator++]=   /*-22-*/ "topomatic second daughter";
  varNamesD0forD0ptbin[iterator++]=   /*-23-*/ "topomatic max";
  varNamesD0forD0ptbin[iterator++]=   /*-24-*/ "topomatic min";

  Bool_t isUpperCutD0forD0ptbin[nvarsD0forD0ptbin]={0};

  SetVarNamesD0forD0ptbin(nvarsD0forD0ptbin,varNamesD0forD0ptbin,isUpperCutD0forD0ptbin);

  Float_t limitsD0forD0ptbin[2]={0,999999999.};
  SetPtBinsD0forD0ptbin(2,limitsD0forD0ptbin);


  //
  // Initialization of D0 cuts for DStar pt bins
  //

  const Int_t nvarsD0forDStarptbin=35;
  SetNVarsD0forDStarptbin(nvarsD0forDStarptbin);

  TString varNamesD0forDStarptbin[nvarsD0forDStarptbin];
  iterator = 0;

  //D0 cut variables
  varNamesD0forDStarptbin[iterator++]=   /*-00-*/ "inv. mass width[GeV]";    
  varNamesD0forDStarptbin[iterator++]=   /*-01-*/ "delta mass width  [GeV]";    //not used for D0
  varNamesD0forDStarptbin[iterator++]=   /*-02-*/ "pointing angle [Cos(theta)]";
  varNamesD0forDStarptbin[iterator++]=   /*-03-*/ "dca [cm]";                   
  varNamesD0forDStarptbin[iterator++]=   /*-04-*/ "Pt D0 [GeV/c]";
  varNamesD0forDStarptbin[iterator++]=   /*-05-*/ "Pt Kaon [GeV/c]";
  varNamesD0forDStarptbin[iterator++]=   /*-06-*/ "Pt Pion [GeV/c]";
  varNamesD0forDStarptbin[iterator++]=   /*-07-*/ "d0 D0 [cm]";
  varNamesD0forDStarptbin[iterator++]=   /*-08-*/ "d0 Kaon [cm]";                
  varNamesD0forDStarptbin[iterator++]=   /*-09-*/ "d0 Pion [cm]";                
  varNamesD0forDStarptbin[iterator++]=   /*-10-*/ "d0d0 [cm^2]";
  varNamesD0forDStarptbin[iterator++]=   /*-11-*/ "d0d0 XY [cm^2]";

  varNamesD0forDStarptbin[iterator++]=   /*-12-*/ "angle between both daughters"; 
  varNamesD0forDStarptbin[iterator++]=   /*-13-*/ "angle mother with first daughter";
  varNamesD0forDStarptbin[iterator++]=   /*-14-*/ "angle mother with second daughter";
  varNamesD0forDStarptbin[iterator++]=   /*-15-*/ "cosThetaStar";                
  varNamesD0forDStarptbin[iterator++]=   /*-16-*/ "vertexDistance"; 
  varNamesD0forDStarptbin[iterator++]=   /*-17-*/ "pseudoProperDecayTime"; 
  varNamesD0forDStarptbin[iterator++]=   /*-18-*/ "DecayTime"; 
  varNamesD0forDStarptbin[iterator++]=   /*-19-*/ "normalizedDecayTime";
  varNamesD0forDStarptbin[iterator++]=   /*-20-*/ "normDecayLength";
  varNamesD0forDStarptbin[iterator++]=   /*-21-*/ "topomatic first daughter";
  varNamesD0forDStarptbin[iterator++]=   /*-22-*/ "topomatic second daughter";
  varNamesD0forDStarptbin[iterator++]=   /*-23-*/ "topomatic max";
  varNamesD0forDStarptbin[iterator++]=   /*-24-*/ "topomatic min";

  varNamesD0forDStarptbin[iterator++]=   /*-25-*/ "pointingAngleToDStar";
  varNamesD0forDStarptbin[iterator++]=   /*-26-*/ "d0MotherToDStar"; 
  varNamesD0forDStarptbin[iterator++]=   /*-27-*/ "d0FirstDaughterToDStar"; 
  varNamesD0forDStarptbin[iterator++]=   /*-28-*/ "d0SecondDaughterToDStar"; 
  varNamesD0forDStarptbin[iterator++]=   /*-29-*/ "impactProductToDStar"; 
  varNamesD0forDStarptbin[iterator++]=   /*-30-*/ "impactProductXYToDStar"; 
  varNamesD0forDStarptbin[iterator++]=   /*-31-*/ "normDecayLengthToDStar"; 
  varNamesD0forDStarptbin[iterator++]=   /*-32-*/ "pseudoProperDecayTimeToDStar"; 
  varNamesD0forDStarptbin[iterator++]=   /*-33-*/ "DecayTimeToDStar"; 
  varNamesD0forDStarptbin[iterator++]=   /*-34-*/ "normalizedDecayTimeToDStar"; 

  Bool_t isUpperCutD0forDStarptbin[nvarsD0forDStarptbin]={0};

  SetVarNamesD0forDStarptbin(nvarsD0forDStarptbin,varNamesD0forDStarptbin,isUpperCutD0forDStarptbin);

  Float_t limitsD0forDStarptbin[2]={0,999999999.};
  SetPtBinsD0forDStarptbin(2,limitsD0forDStarptbin);


  //
  // Initialization of DStar cuts for DStar pt bins
  //

  const Int_t nvarsDStarforDStarptbin=25;
  SetNVarsDStarforDStarptbin(nvarsDStarforDStarptbin);

  TString varNamesDStarforDStarptbin[nvarsDStarforDStarptbin];
  iterator = 0;

  //DStar cut variables
  varNamesDStarforDStarptbin[iterator++]=   /*-00-*/ "inv. mass width[GeV]";    
  varNamesDStarforDStarptbin[iterator++]=   /*-01-*/ "delta mass width  [GeV]";    
  varNamesDStarforDStarptbin[iterator++]=   /*-02-*/ "pointing angle [Cos(theta)]";
  varNamesDStarforDStarptbin[iterator++]=   /*-03-*/ "dca [cm]";                   
  varNamesDStarforDStarptbin[iterator++]=   /*-04-*/ "Pt DStar [GeV/c]";
  varNamesDStarforDStarptbin[iterator++]=   /*-05-*/ "Pt D0 [GeV/c]";
  varNamesDStarforDStarptbin[iterator++]=   /*-06-*/ "Pt Pion [GeV/c]";
  varNamesDStarforDStarptbin[iterator++]=   /*-07-*/ "d0 DStar [cm]";
  varNamesDStarforDStarptbin[iterator++]=   /*-08-*/ "d0 D0 [cm]";                
  varNamesDStarforDStarptbin[iterator++]=   /*-09-*/ "d0 Pion [cm]";                
  varNamesDStarforDStarptbin[iterator++]=   /*-10-*/ "d0d0 [cm^2]";
  varNamesDStarforDStarptbin[iterator++]=   /*-11-*/ "d0d0 XY [cm^2]";

  varNamesDStarforDStarptbin[iterator++]=   /*-12-*/ "angle between both daughters"; 
  varNamesDStarforDStarptbin[iterator++]=   /*-13-*/ "angle mother with first daughter";
  varNamesDStarforDStarptbin[iterator++]=   /*-14-*/ "angle mother with second daughter";
  varNamesDStarforDStarptbin[iterator++]=   /*-15-*/ "cosThetaStar";                
  varNamesDStarforDStarptbin[iterator++]=   /*-16-*/ "vertexDistance"; 
  varNamesDStarforDStarptbin[iterator++]=   /*-17-*/ "pseudoProperDecayTime"; 
  varNamesDStarforDStarptbin[iterator++]=   /*-18-*/ "DecayTime"; 
  varNamesDStarforDStarptbin[iterator++]=   /*-19-*/ "normalizedDecayTime";
  varNamesDStarforDStarptbin[iterator++]=   /*-20-*/ "normDecayLength";
  varNamesDStarforDStarptbin[iterator++]=   /*-21-*/ "topomatic first daughter";
  varNamesDStarforDStarptbin[iterator++]=   /*-22-*/ "topomatic second daughter";
  varNamesDStarforDStarptbin[iterator++]=   /*-23-*/ "topomatic max";
  varNamesDStarforDStarptbin[iterator++]=   /*-24-*/ "topomatic min";

  Bool_t isUpperCutDStarforDStarptbin[nvarsDStarforDStarptbin]={0};

  SetVarNamesDStarforDStarptbin(nvarsDStarforDStarptbin,varNamesDStarforDStarptbin,isUpperCutDStarforDStarptbin);

  Float_t limitsDStarforDStarptbin[2]={0,999999999.};
  SetPtBinsDStarforDStarptbin(2,limitsDStarforDStarptbin);






  Bool_t forOpt[16]={0}; //not yet used for B0 analysis
  SetVarsForOpt(16,forOpt);

}
//--------------------------------------------------------------------------
AliRDHFCutsB0toDStarPi::AliRDHFCutsB0toDStarPi(const AliRDHFCutsB0toDStarPi &source) :
  AliRDHFCuts(source),
  fMaxPtPid(source.fMaxPtPid),
  fTPCflag(source.fTPCflag),
  fCircRadius(source.fCircRadius),
  fIsCutUsed(source.fIsCutUsed),   
  fCutsRDD0forD0ptbin(source.fCutsRDD0forD0ptbin),
  fnPtBinsD0forD0ptbin(source.fnPtBinsD0forD0ptbin),
  fPtBinLimitsD0forD0ptbin(source.fPtBinLimitsD0forD0ptbin),
  fIsUpperCutD0forD0ptbin(source.fIsUpperCutD0forD0ptbin),
  fIsCutUsedD0forD0ptbin(source.fIsCutUsedD0forD0ptbin),
  fVarNamesD0forD0ptbin(source.fVarNamesD0forD0ptbin),
  fCutsRDD0forDStarptbin(source.fCutsRDD0forDStarptbin),
  fnPtBinsD0forDStarptbin(source.fnPtBinsD0forDStarptbin),
  fPtBinLimitsD0forDStarptbin(source.fPtBinLimitsD0forDStarptbin),
  fIsUpperCutD0forDStarptbin(source.fIsUpperCutD0forDStarptbin),
  fIsCutUsedD0forDStarptbin(source.fIsCutUsedD0forDStarptbin),
  fVarNamesD0forDStarptbin(source.fVarNamesD0forDStarptbin),
  fCutsRDDStarforDStarptbin(source.fCutsRDDStarforDStarptbin),
  fnPtBinsDStarforDStarptbin(source.fnPtBinsDStarforDStarptbin),
  fPtBinLimitsDStarforDStarptbin(source.fPtBinLimitsDStarforDStarptbin),
  fIsUpperCutDStarforDStarptbin(source.fIsUpperCutDStarforDStarptbin),
  fIsCutUsedDStarforDStarptbin(source.fIsCutUsedDStarforDStarptbin),
  fVarNamesDStarforDStarptbin(source.fVarNamesDStarforDStarptbin),  
  fnPtBinLimitsD0forD0ptbin(source.fnPtBinLimitsD0forD0ptbin),
  fnPtBinLimitsD0forDStarptbin(source.fnPtBinLimitsD0forDStarptbin),
  fnPtBinLimitsDStarforDStarptbin(source.fnPtBinLimitsDStarforDStarptbin)
{
  //
  // Copy constructor
  // 
}
//--------------------------------------------------------------------------
AliRDHFCutsB0toDStarPi::~AliRDHFCutsB0toDStarPi() {
  //  
  // Default Destructor
  //
  if(fIsCutUsed) { delete [] fIsCutUsed; fIsCutUsed=0; }
  if(fCutsRDD0forD0ptbin) { delete [] fCutsRDD0forD0ptbin; fCutsRDD0forD0ptbin=0; }
  if(fPtBinLimitsD0forD0ptbin) { delete [] fPtBinLimitsD0forD0ptbin; fPtBinLimitsD0forD0ptbin=0; }
  if(fIsUpperCutD0forD0ptbin) { delete [] fIsUpperCutD0forD0ptbin; fIsUpperCutD0forD0ptbin=0; }
  if(fIsCutUsedD0forD0ptbin) { delete [] fIsCutUsedD0forD0ptbin; fIsCutUsedD0forD0ptbin=0; }
  if(fVarNamesD0forD0ptbin) { delete [] fVarNamesD0forD0ptbin; fVarNamesD0forD0ptbin=0; }
  if(fCutsRDD0forDStarptbin) { delete [] fCutsRDD0forDStarptbin; fCutsRDD0forDStarptbin=0; }
  if(fPtBinLimitsD0forDStarptbin) { delete [] fPtBinLimitsD0forDStarptbin; fPtBinLimitsD0forDStarptbin=0; }
  if(fIsUpperCutD0forDStarptbin) { delete [] fIsUpperCutD0forDStarptbin; fIsUpperCutD0forDStarptbin=0; }
  if(fIsCutUsedD0forDStarptbin) { delete [] fIsCutUsedD0forDStarptbin; fIsCutUsedD0forDStarptbin=0; }
  if(fVarNamesD0forDStarptbin) { delete [] fVarNamesD0forDStarptbin; fVarNamesD0forDStarptbin=0; }
  if(fCutsRDDStarforDStarptbin) { delete [] fCutsRDDStarforDStarptbin; fCutsRDDStarforDStarptbin=0; }
  if(fPtBinLimitsDStarforDStarptbin) { delete [] fPtBinLimitsDStarforDStarptbin; fPtBinLimitsDStarforDStarptbin=0; }
  if(fIsUpperCutDStarforDStarptbin) { delete [] fIsUpperCutDStarforDStarptbin; fIsUpperCutDStarforDStarptbin=0; }
  if(fIsCutUsedDStarforDStarptbin) { delete [] fIsCutUsedDStarforDStarptbin; fIsCutUsedDStarforDStarptbin=0; }
  if(fVarNamesDStarforDStarptbin) { delete [] fVarNamesDStarforDStarptbin; fVarNamesDStarforDStarptbin=0; }
}
//--------------------------------------------------------------------------
AliRDHFCutsB0toDStarPi &AliRDHFCutsB0toDStarPi::operator=(const AliRDHFCutsB0toDStarPi &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliRDHFCuts::operator=(source);

  return *this;
}
//--------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  // not yet used

  return;
}
//--------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::IsSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[85]) {
  //
  // In this function we apply the selection cuts on the B0 candidate and its daughters.
  // The function returns 0 if the candidate is cut and is able to return information on which cuts the candidate passes.
  //

  fIsSelectedCuts=0;
  fIsSelectedPID=0;

  // The cuts used in this class have to be set in the maketfile.
  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoCascadeHF* candidateB0 = (AliAODRecoCascadeHF*)obj;
  if(!candidateB0){
    cout<<"candidateB0 null"<<endl;
    return 0;
  } 

  AliAODRecoCascadeHF* candidateDStar = (AliAODRecoCascadeHF*)candidateB0->GetDaughter(1);  
  if(!candidateDStar){
    cout<<"candidateDStar null"<<endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateB0->GetDaughter(0);
  if(!candidatePion){
    cout<<"candidatePion null"<<endl;
    return 0;
  }
  
  AliAODVertex * vertexB0 = candidateB0->GetSecondaryVtx();
  if(!vertexB0){
    cout<<"vertexB0 null"<<endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if(!primaryVertex){
    cout<<"primaryVertex null"<<endl;
    return 0;
  }

  Int_t returnvalue=1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField(); 

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {   

    // We check to which pt bin the candidate belongs
    Int_t ptbin = PtBin(candidateB0->Pt());
    if(ptbin==-1) return -1;    
    
    // We obtain the variable values in the section below 
    // D0Mass and B0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t mB0PDG = TDatabasePDG::Instance()->GetParticle(511)->Mass();
    
    // delta mass PDG
    Double_t deltaPDG = mD0PDG-mB0PDG;
   
    // Half width B0 mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 413;
    Double_t invMassB0 = candidateB0->InvMass(2,prongs);
    Double_t invMassDifference = TMath::Abs(mB0PDG - invMassB0);
    Double_t invMassDelta = TMath::Abs(deltaPDG-(DeltaInvMassB0Kpipipi(candidateB0)));

    Double_t pointingAngle = candidateB0->CosPointingAngle();
    Double_t dcaMother = candidateB0->GetDCA();
    Double_t ptMother = candidateB0->Pt();
    Double_t momentumMother = candidateB0->P();
    Double_t ptDStar = candidateDStar->Pt();
    Double_t ptPion = candidatePion->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateB0);
    Double_t d0z0[2],covd0z0[3],d0[2];
    motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateB0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateB0->Getd0Prong(1));  
    
    Double_t impactProduct = candidateB0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateB0->ImpParXY());  

    Double_t angleBetweenBothDaughters  = (candidateDStar->Px() * candidatePion->Px() + candidateDStar->Py() * candidatePion->Py() + candidateDStar->Pz() * candidatePion->Pz()) /(candidateDStar->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateB0->Px() * candidatePion->Px() + candidateB0->Py() * candidatePion->Py() + candidateB0->Pz() * candidatePion->Pz()) /(candidateB0->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateB0->Px() * candidateDStar->Px() + candidateB0->Py() * candidateDStar->Py() + candidateB0->Pz() * candidateDStar->Pz()) /(candidateB0->P() * candidateDStar->P());

    Double_t cosThetaStar = candidateB0->CosThetaStar(0,511,211,413);
    Double_t vertexDistance = vertexB0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateB0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(511)->Mass();
    Double_t pseudoProperDecayLength = ((vertexB0->GetX() - primaryVertex->GetX()) * candidateB0->Px() / TMath::Abs(candidateB0->Pt())) + ((vertexB0->GetY() - primaryVertex->GetY()) * candidateB0->Py() / TMath::Abs(candidateB0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = candidateB0->Phi();
    Double_t theta = candidateB0->Theta();
    Double_t covMatrix[21];
    candidateB0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTime = candidateB0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    //Topomatic
    Double_t dd0pr1=0.;
    Double_t dd0pr2=0.;
    Double_t dd0max=0.;
    Double_t dd0min=0.;
    for(Int_t ipr=0; ipr<2; ipr++) 
    {
      Double_t diffIP, errdiffIP;
      candidateB0->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
      Double_t normdd0=0.;
      if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
      if(ipr==0) dd0pr1=normdd0;
      if(ipr==1) dd0pr2=normdd0;
    }
    if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
    else {dd0max=dd0pr2; dd0min=dd0pr1;}


    // We apply the cuts 
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 60;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 61;
    cutVariableValue = invMassDelta;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 62;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 63;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt B0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 64;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt DStar [GeV/c]" -------------------------------------------------
    nCutIndex = 65;
    cutVariableValue = ptDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 66;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 B0 [cm]" -------------------------------------------------------
    nCutIndex = 67;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 DStar [cm]"-----------------------------------------------------
    nCutIndex = 68;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 69;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 70;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 71;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 72;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 73;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 74;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 75;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 76;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 77;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 78;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 79;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 80;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 81;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 82;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 83;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 84;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // select DStar that passes D0 cuts
    bPassedCut = IsDStarFromB0Selected(ptMother,candidateDStar,selectionLevel,aod,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;    
  }

  if(bPassedCut==kFALSE)
  {
    returnvalue = 0;
  } else
  {
    for (Int_t i = 60; i < 85; ++i)
    {
      if(bCutArray[i]==kTRUE){
        returnvalue = 0;
        break;
      }
    }
  }

  fIsSelectedCuts = returnvalue;
  
  return returnvalue;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::IsDStarFromB0Selected(Double_t ptB0, TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[85]) {
  //
  // Apply selection for DStar candidate
  //

  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoCascadeHF* candidateDStar = (AliAODRecoCascadeHF*)obj;
  if(!candidateDStar){
    cout<<"candidateDStar null"<<endl;
    return 0;
  }
  
  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateDStar->Get2Prong();  
  if(!candidateD0){
    cout<<"candidateD0 null"<<endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateDStar->GetDaughter(0);
  if(!candidatePion){
    cout<<"candidatePion null"<<endl;
    return 0;
  }

  AliAODVertex * vertexDStar = candidateDStar->GetSecondaryVtx();
  if(!vertexDStar){
    cout<<"vertexDStar null"<<endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if(!primaryVertex){
    cout<<"primaryVertex null"<<endl;
    return 0;
  }
  
  Int_t returnvalue=1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField(); 

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Int_t ptbin=PtBin(ptB0);    
    
    // DStarMass and D0mass
    Double_t mDSPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    // delta mass PDG
    Double_t deltaPDG = mDSPDG-mD0PDG;
      
    // Half width DStar mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 421;
    Double_t invMassDStar = candidateDStar->InvMass(2,prongs);
    Double_t invMassDifference = TMath::Abs(mDSPDG - invMassDStar);
    Double_t invMassDelta = TMath::Abs(deltaPDG-(DeltaInvMassDStarKpipi(candidateDStar)));

    Double_t pointingAngle = candidateDStar->CosPointingAngle();
    Double_t dcaMother = candidateDStar->GetDCA();
    Double_t ptMother = candidateDStar->Pt();
    Double_t momentumMother = candidateDStar->P();
    Double_t ptD0 = candidateD0->Pt();
    Double_t ptPion = candidatePion->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateDStar);
    Double_t d0z0[2],covd0z0[3],d0[2];
    motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateDStar->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateDStar->Getd0Prong(1));  
    
    Double_t impactProduct = candidateDStar->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateDStar->ImpParXY());  

    Double_t angleBetweenBothDaughters  = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) /(candidateD0->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateDStar->Px() * candidatePion->Px() + candidateDStar->Py() * candidatePion->Py() + candidateDStar->Pz() * candidatePion->Pz()) /(candidateDStar->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateDStar->Px() * candidateD0->Px() + candidateDStar->Py() * candidateD0->Py() + candidateDStar->Pz() * candidateD0->Pz()) /(candidateDStar->P() * candidateD0->P());

    Double_t cosThetaStar = candidateDStar->CosThetaStar(0,413,211,421);
    Double_t vertexDistance = vertexDStar->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateDStar->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t pseudoProperDecayLength = ((vertexDStar->GetX() - primaryVertex->GetX()) * candidateDStar->Px() / TMath::Abs(candidateDStar->Pt())) + ((vertexDStar->GetY() - primaryVertex->GetY()) * candidateDStar->Py() / TMath::Abs(candidateDStar->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = candidateDStar->Phi();
    Double_t theta = candidateDStar->Theta();
    Double_t covMatrix[21];
    candidateDStar->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTime = candidateDStar->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    //Topomatic
    Double_t dd0pr1=0.;
    Double_t dd0pr2=0.;
    Double_t dd0max=0.;
    Double_t dd0min=0.;
    for(Int_t ipr=0; ipr<2; ipr++) 
    {
      Double_t diffIP, errdiffIP;
      candidateDStar->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
      Double_t normdd0=0.;
      if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
      if(ipr==0) dd0pr1=normdd0;
      if(ipr==1) dd0pr2=normdd0;
    }
    if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
    else {dd0max=dd0pr2; dd0min=dd0pr1;}


    // We apply the cuts 
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 35;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 36;
    cutVariableValue = invMassDelta;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 37;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 38;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt DStar [GeV/c]" ----------------------------------------------------
    nCutIndex = 39;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" -------------------------------------------------
    nCutIndex = 40;
    cutVariableValue = ptD0;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 41;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 DStar [cm]" -------------------------------------------------------
    nCutIndex = 42;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]"-----------------------------------------------------
    nCutIndex = 43;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 44;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 45;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 46;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 47;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 48;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 49;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 50;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 51;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 52;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 53;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 54;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 55;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 56;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 57;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 58;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 59;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // select D0 - have to pass DStar candidate to get variables w.r.t. DStar vertex.
    bPassedCut = IsD0FromDStarSelected(ptB0,candidateDStar,selectionLevel, aod, bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;    
  }

  if(bPassedCut==kFALSE)
  {
    returnvalue = 0;
  } else
  {
    for (Int_t i = 35; i < 60; ++i)
    {
      if(bCutArray[i]==kTRUE){
        returnvalue = 0;
        break;
      }
    }
  }

  
  return returnvalue;

}
//_________________________________________________________________________________________________
Int_t AliRDHFCutsB0toDStarPi::IsD0FromDStarSelected(Double_t ptB0, TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[85]) {
  //
  // Apply selection on D0 candidate from DStar candidate. We have to pass the DStar candidate to this function to get variables w.r.t. DStar vertex.
  // 
  
  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoCascadeHF* candidateDStar = (AliAODRecoCascadeHF*)obj;
  if(!candidateDStar){
    cout<<"candidateDStar null"<<endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateDStar->Get2Prong();  
  if(!candidateD0){
    cout<<"candidateD0 null"<<endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateD0->GetDaughter(0);
  if(!candidatePion){
    cout<<"candidatePion null"<<endl;
    return 0;
  }

  AliAODTrack *candidateKaon = (AliAODTrack*)candidateD0->GetDaughter(1);
  if(!candidateKaon){
    cout<<"candidateKaon null"<<endl;
    return 0;
  }

  AliAODVertex * vertexDStar = candidateDStar->GetSecondaryVtx();
  if(!vertexDStar){
    cout<<"vertexDStar null"<<endl;
    return 0;
  }

  AliAODVertex * vertexD0 = candidateD0->GetSecondaryVtx();
  if(!vertexD0){
    cout<<"vertexD0 null"<<endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if(!primaryVertex){
    cout<<"primaryVertex null"<<endl;
    return 0;
  }

  Int_t returnvalue=1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField(); 


  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Int_t ptbin=PtBin(ptB0);    
    
    // D0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  
    // Half width DStar mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 321;
    Double_t invMassD0 = candidateD0->InvMass(2,prongs);
    Double_t invMassDifference = TMath::Abs(mD0PDG - invMassD0);

    Double_t pointingAngle = candidateD0->CosPointingAngle();
    Double_t dcaMother = candidateD0->GetDCA();
    Double_t ptMother = candidateD0->Pt();
    Double_t momentumMother = candidateD0->P();
    Double_t ptPion = candidatePion->Pt();
    Double_t ptKaon = candidateKaon->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateD0);
    Double_t d0z0[2],covd0z0[3],d0[2];
    motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateD0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateD0->Getd0Prong(1));  
    
    Double_t impactProduct = candidateD0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateD0->ImpParXY());  

    Double_t angleBetweenBothDaughters  = (candidateKaon->Px() * candidatePion->Px() + candidateKaon->Py() * candidatePion->Py() + candidateKaon->Pz() * candidatePion->Pz()) /(candidateKaon->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) /(candidateD0->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateD0->Px() * candidateKaon->Px() + candidateD0->Py() * candidateKaon->Py() + candidateD0->Pz() * candidateKaon->Pz()) /(candidateD0->P() * candidateKaon->P());

    Double_t cosThetaStar = candidateD0->CosThetaStar(0,421,211,321);
    Double_t vertexDistance = vertexD0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateD0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - primaryVertex->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - primaryVertex->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = candidateD0->Phi();
    Double_t theta = candidateD0->Theta();
    Double_t covMatrix[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTime = candidateD0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    //Topomatic
    Double_t dd0pr1=0.;
    Double_t dd0pr2=0.;
    Double_t dd0max=0.;
    Double_t dd0min=0.;
    for(Int_t ipr=0; ipr<2; ipr++) 
    {
      Double_t diffIP, errdiffIP;
      candidateD0->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
      Double_t normdd0=0.;
      if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
      if(ipr==0) dd0pr1=normdd0;
      if(ipr==1) dd0pr2=normdd0;
    }
    if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
    else {dd0max=dd0pr2; dd0min=dd0pr1;}


    // We apply the cuts 
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    // nCutIndex = 1; // not used for D0
    // cutVariableValue = invMassDelta;
    // bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    // if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = ptKaon;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]" -------------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------
    
    AliAODRecoDecay* candidateD0toDStar = (AliAODRecoDecay*)candidateD0;
    AliExternalTrackParam pionD0Track;
    AliExternalTrackParam kaonD0Track;

    Double_t d0z0DSVert[2],covd0z0DSVert[3],d0DSVert[2];

    pionD0Track.CopyFromVTrack(candidatePion);
    pionD0Track.PropagateToDCA(vertexDStar,bz,100.,d0z0DSVert,covd0z0DSVert);
    d0DSVert[0] = d0z0DSVert[0];

    kaonD0Track.CopyFromVTrack(candidateKaon);
    kaonD0Track.PropagateToDCA(vertexDStar,bz,100.,d0z0DSVert,covd0z0DSVert);
    d0DSVert[1] = d0z0DSVert[0];

    AliExternalTrackParam D0Track;
    D0Track.CopyFromVTrack(candidateD0);
    Double_t d0z0D0DSVert[2],covd0z0D0DSVert[3],d0D0DSVert;
    motherTrack.PropagateToDCA(vertexDStar,bz,100.,d0z0D0DSVert,covd0z0D0DSVert);
    d0D0DSVert = TMath::Abs(d0z0D0DSVert[0]);

    Double_t impactProductToDStar = d0DSVert[0]*d0DSVert[1];
    Double_t impactProductXYToDStar = candidateD0toDStar->ImpParXY(vertexDStar);

    Double_t pointingAngleToDStar = candidateD0toDStar->CosPointingAngle(vertexDStar);
    Double_t d0FirstDaughterToDStar = TMath::Abs(d0DSVert[0]);
    Double_t d0SecondDaughterToDStar = TMath::Abs(d0DSVert[1]);
    Double_t normDecayLengthToDStar = candidateD0toDStar->NormalizedDecayLength(vertexDStar);

    Double_t pseudoProperDecayLengthDSVert = ((vertexD0->GetX() - vertexDStar->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - vertexDStar->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTimeToDStar = pseudoProperDecayLengthDSVert * pdgMassMother/ptMother;
    Double_t DecayTimeToDStar = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phiDSVert = candidateD0->Phi();
    Double_t thetaDSVert = candidateD0->Theta();
    Double_t covMatrixDSVert[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrixDSVert);

    cp = TMath::Cos(phiDSVert);
    sp = TMath::Sin(phiDSVert);
    ct = TMath::Cos(thetaDSVert);
    st = TMath::Sin(thetaDSVert);

    errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTimeToDStar = candidateD0toDStar->NormalizedDecayLength(vertexDStar) / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    // "pointingAngleToDStar" ---------------------------------------------
    nCutIndex = 25;
    cutVariableValue = pointingAngleToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0MotherToDStar" --------------------------------------------------
    nCutIndex = 26;
    cutVariableValue = d0D0DSVert;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0FirstDaughterToDStar" -------------------------------------------
    nCutIndex = 27;
    cutVariableValue = d0FirstDaughterToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0SecondDaughterToDStar" ------------------------------------------
    nCutIndex = 28;
    cutVariableValue = d0SecondDaughterToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "impactProductToDStar" ---------------------------------------------
    nCutIndex = 29;
    cutVariableValue = impactProductToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "impactProductXYToDStar" -------------------------------------------
    nCutIndex = 30;
    cutVariableValue = impactProductXYToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normDecayLengthToDStar" -------------------------------------------
    nCutIndex = 31;
    cutVariableValue = normDecayLengthToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTimeToDStar" -------------------------------------
    nCutIndex = 32;
    cutVariableValue = pseudoProperDecayTimeToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "DecayTimeToDStar" -------------------------------------------------
    nCutIndex = 33;
    cutVariableValue = DecayTimeToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTimeToDStar" ---------------------------------------------
    nCutIndex = 34;
    cutVariableValue = normalizedDecayTimeToDStar;
    bPassedCut = ApplyCutOnVariable(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------
  }
 
  for (Int_t i = 0; i < 35; ++i)
  {
    if(bCutArray[i]==kTRUE){
      returnvalue = 0;
      break;
    }
  }

  return returnvalue;
}
//----------------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::IsD0forD0ptbinSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[25]) {
 //
  // Apply selection on D0 candidate.
  // 
  
  if(!fCutsRDD0forD0ptbin){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)obj; 
  if(!candidateD0){
    cout<<"candidateD0 null"<<endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateD0->GetDaughter(0);
  if(!candidatePion){
    cout<<"candidatePion null"<<endl;
    return 0;
  }

  AliAODTrack *candidateKaon = (AliAODTrack*)candidateD0->GetDaughter(1);
  if(!candidateKaon){
    cout<<"candidateKaon null"<<endl;
    return 0;
  }

  AliAODVertex * vertexD0 = candidateD0->GetSecondaryVtx();
  if(!vertexD0){
    cout<<"vertexD0 null"<<endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if(!primaryVertex){
    cout<<"primaryVertex null"<<endl;
    return 0;
  }

  Int_t returnvalue=1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField(); 


  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Int_t ptbin=PtBinD0forD0ptbin(candidateD0->Pt());    
    if(ptbin==-1) return -1;  

    // D0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  
    // Half width DStar mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 321;
    Double_t invMassD0 = candidateD0->InvMass(2,prongs);
    Double_t invMassDifference = TMath::Abs(mD0PDG - invMassD0);

    Double_t pointingAngle = candidateD0->CosPointingAngle();
    Double_t dcaMother = candidateD0->GetDCA();
    Double_t ptMother = candidateD0->Pt();
    Double_t momentumMother = candidateD0->P();
    Double_t ptPion = candidatePion->Pt();
    Double_t ptKaon = candidateKaon->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateD0);
    Double_t d0z0[2],covd0z0[3],d0[2];
    motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateD0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateD0->Getd0Prong(1));  
    
    Double_t impactProduct = candidateD0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateD0->ImpParXY());  

    Double_t angleBetweenBothDaughters  = (candidateKaon->Px() * candidatePion->Px() + candidateKaon->Py() * candidatePion->Py() + candidateKaon->Pz() * candidatePion->Pz()) /(candidateKaon->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) /(candidateD0->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateD0->Px() * candidateKaon->Px() + candidateD0->Py() * candidateKaon->Py() + candidateD0->Pz() * candidateKaon->Pz()) /(candidateD0->P() * candidateKaon->P());

    Double_t cosThetaStar = candidateD0->CosThetaStar(0,421,211,321);
    Double_t vertexDistance = vertexD0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateD0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - primaryVertex->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - primaryVertex->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = candidateD0->Phi();
    Double_t theta = candidateD0->Theta();
    Double_t covMatrix[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTime = candidateD0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    //Topomatic
    Double_t dd0pr1=0.;
    Double_t dd0pr2=0.;
    Double_t dd0max=0.;
    Double_t dd0min=0.;
    for(Int_t ipr=0; ipr<2; ipr++) 
    {
      Double_t diffIP, errdiffIP;
      candidateD0->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
      Double_t normdd0=0.;
      if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
      if(ipr==0) dd0pr1=normdd0;
      if(ipr==1) dd0pr2=normdd0;
    }
    if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
    else {dd0max=dd0pr2; dd0min=dd0pr1;}


    // We apply the cuts 
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    // nCutIndex = 1; // not used for D0
    // cutVariableValue = invMassDelta;
    // bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    // if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = ptKaon;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]" -------------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariableD0forD0ptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------
    
  }
 
  for (Int_t i = 0; i < 25; ++i)
  {
    if(bCutArray[i]==kTRUE){
      returnvalue = 0;
      break;
    }
  }

  return returnvalue;
}
//----------------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::IsD0forDStarptbinSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[35]) {
  //
  // Apply selection on D0 candidate from DStar candidate. We have to pass the DStar candidate to this function to get variables w.r.t. DStar vertex.
  // 
  
  if(!fCutsRDD0forDStarptbin){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoCascadeHF* candidateDStar = (AliAODRecoCascadeHF*)obj;
  if(!candidateDStar){
    cout<<"candidateDStar null"<<endl;
    return 0;
  }

  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateDStar->Get2Prong();  
  if(!candidateD0){
    cout<<"candidateD0 null"<<endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateD0->GetDaughter(0);
  if(!candidatePion){
    cout<<"candidatePion null"<<endl;
    return 0;
  }

  AliAODTrack *candidateKaon = (AliAODTrack*)candidateD0->GetDaughter(1);
  if(!candidateKaon){
    cout<<"candidateKaon null"<<endl;
    return 0;
  }

  AliAODVertex * vertexDStar = candidateDStar->GetSecondaryVtx();
  if(!vertexDStar){
    cout<<"vertexDStar null"<<endl;
    return 0;
  }

  AliAODVertex * vertexD0 = candidateD0->GetSecondaryVtx();
  if(!vertexD0){
    cout<<"vertexD0 null"<<endl;
    return 0;
  }

  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if(!primaryVertex){
    cout<<"primaryVertex null"<<endl;
    return 0;
  }

  Int_t returnvalue=1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField(); 


  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Int_t ptbin=PtBinD0forDStarptbin(candidateDStar->Pt());       
    if(ptbin==-1) return -1;    
    // D0mass
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  
    // Half width DStar mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 321;
    Double_t invMassD0 = candidateD0->InvMass(2,prongs);
    Double_t invMassDifference = TMath::Abs(mD0PDG - invMassD0);

    Double_t pointingAngle = candidateD0->CosPointingAngle();
    Double_t dcaMother = candidateD0->GetDCA();
    Double_t ptMother = candidateD0->Pt();
    Double_t momentumMother = candidateD0->P();
    Double_t ptPion = candidatePion->Pt();
    Double_t ptKaon = candidateKaon->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateD0);
    Double_t d0z0[2],covd0z0[3],d0[2];
    motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateD0->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateD0->Getd0Prong(1));  
    
    Double_t impactProduct = candidateD0->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateD0->ImpParXY());  

    Double_t angleBetweenBothDaughters  = (candidateKaon->Px() * candidatePion->Px() + candidateKaon->Py() * candidatePion->Py() + candidateKaon->Pz() * candidatePion->Pz()) /(candidateKaon->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) /(candidateD0->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateD0->Px() * candidateKaon->Px() + candidateD0->Py() * candidateKaon->Py() + candidateD0->Pz() * candidateKaon->Pz()) /(candidateD0->P() * candidateKaon->P());

    Double_t cosThetaStar = candidateD0->CosThetaStar(0,421,211,321);
    Double_t vertexDistance = vertexD0->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateD0->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();
    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - primaryVertex->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - primaryVertex->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = candidateD0->Phi();
    Double_t theta = candidateD0->Theta();
    Double_t covMatrix[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTime = candidateD0->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    //Topomatic
    Double_t dd0pr1=0.;
    Double_t dd0pr2=0.;
    Double_t dd0max=0.;
    Double_t dd0min=0.;
    for(Int_t ipr=0; ipr<2; ipr++) 
    {
      Double_t diffIP, errdiffIP;
      candidateD0->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
      Double_t normdd0=0.;
      if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
      if(ipr==0) dd0pr1=normdd0;
      if(ipr==1) dd0pr2=normdd0;
    }
    if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
    else {dd0max=dd0pr2; dd0min=dd0pr1;}


    // We apply the cuts 
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    // nCutIndex = 1; // not used for D0
    // cutVariableValue = invMassDelta;
    // bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    // if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Kaon [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = ptKaon;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]" -------------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Kaon [cm]"-----------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------
    
    AliAODRecoDecay* candidateD0toDStar = (AliAODRecoDecay*)candidateD0;
    AliExternalTrackParam pionD0Track;
    AliExternalTrackParam kaonD0Track;

    Double_t d0z0DSVert[2],covd0z0DSVert[3],d0DSVert[2];

    pionD0Track.CopyFromVTrack(candidatePion);
    pionD0Track.PropagateToDCA(vertexDStar,bz,100.,d0z0DSVert,covd0z0DSVert);
    d0DSVert[0] = d0z0DSVert[0];

    kaonD0Track.CopyFromVTrack(candidateKaon);
    kaonD0Track.PropagateToDCA(vertexDStar,bz,100.,d0z0DSVert,covd0z0DSVert);
    d0DSVert[1] = d0z0DSVert[0];

    AliExternalTrackParam D0Track;
    D0Track.CopyFromVTrack(candidateD0);
    Double_t d0z0D0DSVert[2],covd0z0D0DSVert[3],d0D0DSVert;
    motherTrack.PropagateToDCA(vertexDStar,bz,100.,d0z0D0DSVert,covd0z0D0DSVert);
    d0D0DSVert = TMath::Abs(d0z0D0DSVert[0]);

    Double_t impactProductToDStar = d0DSVert[0]*d0DSVert[1];
    Double_t impactProductXYToDStar = candidateD0toDStar->ImpParXY(vertexDStar);

    Double_t pointingAngleToDStar = candidateD0toDStar->CosPointingAngle(vertexDStar);
    Double_t d0FirstDaughterToDStar = TMath::Abs(d0DSVert[0]);
    Double_t d0SecondDaughterToDStar = TMath::Abs(d0DSVert[1]);
    Double_t normDecayLengthToDStar = candidateD0toDStar->NormalizedDecayLength(vertexDStar);

    Double_t pseudoProperDecayLengthDSVert = ((vertexD0->GetX() - vertexDStar->GetX()) * candidateD0->Px() / TMath::Abs(candidateD0->Pt())) + ((vertexD0->GetY() - vertexDStar->GetY()) * candidateD0->Py() / TMath::Abs(candidateD0->Pt()));
    Double_t pseudoProperDecayTimeToDStar = pseudoProperDecayLengthDSVert * pdgMassMother/ptMother;
    Double_t DecayTimeToDStar = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phiDSVert = candidateD0->Phi();
    Double_t thetaDSVert = candidateD0->Theta();
    Double_t covMatrixDSVert[21];
    candidateD0->GetCovarianceXYZPxPyPz(covMatrixDSVert);

    cp = TMath::Cos(phiDSVert);
    sp = TMath::Sin(phiDSVert);
    ct = TMath::Cos(thetaDSVert);
    st = TMath::Sin(thetaDSVert);

    errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTimeToDStar = candidateD0toDStar->NormalizedDecayLength(vertexDStar) / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    // "pointingAngleToDStar" ---------------------------------------------
    nCutIndex = 25;
    cutVariableValue = pointingAngleToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0MotherToDStar" --------------------------------------------------
    nCutIndex = 26;
    cutVariableValue = d0D0DSVert;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0FirstDaughterToDStar" -------------------------------------------
    nCutIndex = 27;
    cutVariableValue = d0FirstDaughterToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0SecondDaughterToDStar" ------------------------------------------
    nCutIndex = 28;
    cutVariableValue = d0SecondDaughterToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "impactProductToDStar" ---------------------------------------------
    nCutIndex = 29;
    cutVariableValue = impactProductToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "impactProductXYToDStar" -------------------------------------------
    nCutIndex = 30;
    cutVariableValue = impactProductXYToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normDecayLengthToDStar" -------------------------------------------
    nCutIndex = 31;
    cutVariableValue = normDecayLengthToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTimeToDStar" -------------------------------------
    nCutIndex = 32;
    cutVariableValue = pseudoProperDecayTimeToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "DecayTimeToDStar" -------------------------------------------------
    nCutIndex = 33;
    cutVariableValue = DecayTimeToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTimeToDStar" ---------------------------------------------
    nCutIndex = 34;
    cutVariableValue = normalizedDecayTimeToDStar;
    bPassedCut = ApplyCutOnVariableD0forDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------
  }
 
  for (Int_t i = 0; i < 35; ++i)
  {
    if(bCutArray[i]==kTRUE){
      returnvalue = 0;
      break;
    }
  }

  return returnvalue;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::IsDStarforDStarptbinSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod, Bool_t bCutArray[25]) { 
  //
  // Apply selection for DStar candidate
  //

  if(!fCutsRDDStarforDStarptbin){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  
  AliAODRecoCascadeHF* candidateDStar = (AliAODRecoCascadeHF*)obj;
  if(!candidateDStar){
    cout<<"candidateDStar null"<<endl;
    return 0;
  }
  
  AliAODRecoDecayHF2Prong* candidateD0 = (AliAODRecoDecayHF2Prong*)candidateDStar->Get2Prong();  
  if(!candidateD0){
    cout<<"candidateD0 null"<<endl;
    return 0;
  }

  AliAODTrack *candidatePion = (AliAODTrack*)candidateDStar->GetDaughter(0);
  if(!candidatePion){
    cout<<"candidatePion null"<<endl;
    return 0;
  }

  AliAODVertex * vertexDStar = candidateDStar->GetSecondaryVtx();
  if(!vertexDStar){
    cout<<"vertexDStar null"<<endl;
    return 0;
  }
  
  AliAODVertex * primaryVertex = aod->GetPrimaryVertex();
  if(!primaryVertex){
    cout<<"primaryVertex null"<<endl;
    return 0;
  }

  Int_t returnvalue=1;
  Bool_t bPassedCut = kFALSE;

  //get the magnetic field
  Double_t bz = (Double_t)aod->GetMagneticField(); 

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {
    
    Int_t ptbin=PtBinDStarforDStarptbin(candidateDStar->Pt());    
    if(ptbin==-1) return -1;  

    // DStarMass and D0mass
    Double_t mDSPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    // delta mass PDG
    Double_t deltaPDG = mDSPDG-mD0PDG;

    // Half width DStar mass
    UInt_t prongs[2];
    prongs[0] = 211; prongs[1] = 421;
    Double_t invMassDStar = candidateDStar->InvMass(2,prongs);
    Double_t invMassDifference = TMath::Abs(mDSPDG - invMassDStar);
    Double_t invMassDelta = TMath::Abs(deltaPDG-(DeltaInvMassDStarKpipi(candidateDStar)));

    Double_t pointingAngle = candidateDStar->CosPointingAngle();
    Double_t dcaMother = candidateDStar->GetDCA();
    Double_t ptMother = candidateDStar->Pt();
    Double_t momentumMother = candidateDStar->P();
    Double_t ptD0 = candidateD0->Pt();
    Double_t ptPion = candidatePion->Pt();

    AliExternalTrackParam motherTrack;
    motherTrack.CopyFromVTrack(candidateDStar);
    Double_t d0z0[2],covd0z0[3],d0[2];
    motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];
    Double_t d0Mother = TMath::Abs(d0[0]);
    Double_t d0firstTrack = TMath::Abs(candidateDStar->Getd0Prong(0));
    Double_t d0secondTrack = TMath::Abs(candidateDStar->Getd0Prong(1));  

    Double_t impactProduct = candidateDStar->Prodd0d0();
    Double_t impactProductXY = TMath::Abs(candidateDStar->ImpParXY());  

    Double_t angleBetweenBothDaughters  = (candidateD0->Px() * candidatePion->Px() + candidateD0->Py() * candidatePion->Py() + candidateD0->Pz() * candidatePion->Pz()) /(candidateD0->P() * candidatePion->P());
    Double_t angleMotherFirstDaughter = (candidateDStar->Px() * candidatePion->Px() + candidateDStar->Py() * candidatePion->Py() + candidateDStar->Pz() * candidatePion->Pz()) /(candidateDStar->P() * candidatePion->P());
    Double_t angleMotherSecondDaughter = (candidateDStar->Px() * candidateD0->Px() + candidateDStar->Py() * candidateD0->Py() + candidateDStar->Pz() * candidateD0->Pz()) /(candidateDStar->P() * candidateD0->P());

    Double_t cosThetaStar = candidateDStar->CosThetaStar(0,413,211,421);
    Double_t vertexDistance = vertexDStar->DistanceToVertex(primaryVertex);
    Double_t normDecayLength = candidateDStar->NormalizedDecayLength();
    Double_t pdgMassMother = TDatabasePDG::Instance()->GetParticle(413)->Mass();
    Double_t pseudoProperDecayLength = ((vertexDStar->GetX() - primaryVertex->GetX()) * candidateDStar->Px() / TMath::Abs(candidateDStar->Pt())) + ((vertexDStar->GetY() - primaryVertex->GetY()) * candidateDStar->Py() / TMath::Abs(candidateDStar->Pt()));
    Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = candidateDStar->Phi();
    Double_t theta = candidateDStar->Theta();
    Double_t covMatrix[21];
    candidateDStar->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normalizedDecayTime = candidateDStar->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    //Topomatic
    Double_t dd0pr1=0.;
    Double_t dd0pr2=0.;
    Double_t dd0max=0.;
    Double_t dd0min=0.;
    for(Int_t ipr=0; ipr<2; ipr++) 
    {
      Double_t diffIP, errdiffIP;
      candidateDStar->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
      Double_t normdd0=0.;
      if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
      if(ipr==0) dd0pr1=normdd0;
      if(ipr==1) dd0pr2=normdd0;
    }
    if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
    else {dd0max=dd0pr2; dd0min=dd0pr1;}


    // We apply the cuts 
    Int_t nCutIndex = 0;
    Double_t cutVariableValue = 0.0;

    // "inv. mass width [GeV]" --------------------------------------------
    nCutIndex = 0;
    cutVariableValue = invMassDifference;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "delta mass width [GeV]" -------------------------------------------
    nCutIndex = 1;
    cutVariableValue = invMassDelta;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pointing angle [Cos(theta)]" --------------------------------------
    nCutIndex = 2;
    cutVariableValue = pointingAngle;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "dca [cm]" ---------------------------------------------------------
    nCutIndex = 3;
    cutVariableValue = dcaMother;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt DStar [GeV/c]" ----------------------------------------------------
    nCutIndex = 4;
    cutVariableValue = ptMother;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt D0 [GeV/c]" -------------------------------------------------
    nCutIndex = 5;
    cutVariableValue = ptD0;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "Pt Pion [GeV/c]" --------------------------------------------------
    nCutIndex = 6;
    cutVariableValue = ptPion;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 DStar [cm]" -------------------------------------------------------
    nCutIndex = 7;
    cutVariableValue = d0Mother;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 D0 [cm]"-----------------------------------------------------
    nCutIndex = 8;
    cutVariableValue = d0firstTrack;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0 Pion [cm]" -----------------------------------------------------
    nCutIndex = 9;
    cutVariableValue = d0secondTrack;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 [cm^2]" ------------------------------------------------------
    nCutIndex = 10;
    cutVariableValue = impactProduct;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "d0d0 XY [cm^2]" ---------------------------------------------------
    nCutIndex = 11;
    cutVariableValue = impactProductXY;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle between both daughters" -------------------------------------
    nCutIndex = 12;
    cutVariableValue = angleBetweenBothDaughters;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with first daughter" ---------------------------------
    nCutIndex = 13;
    cutVariableValue = angleMotherFirstDaughter;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "angle mother with second daughter" --------------------------------
    nCutIndex = 14;
    cutVariableValue = angleMotherSecondDaughter;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "cosThetaStar" -----------------------------------------------------
    nCutIndex = 15;
    cutVariableValue = cosThetaStar;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "vertexDistance" ---------------------------------------------------
    nCutIndex = 16;
    cutVariableValue = vertexDistance;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "pseudoProperDecayTime" --------------------------------------------
    nCutIndex = 17;
    cutVariableValue = pseudoProperDecayTime;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "DecayTime" --------------------------------------------------------
    nCutIndex = 18;
    cutVariableValue = decayTime;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normalizedDecayTime" ----------------------------------------------------
    nCutIndex = 19;
    cutVariableValue = normalizedDecayTime;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "normDecayLength" --------------------------------------------------
    nCutIndex = 20;
    cutVariableValue = normDecayLength;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic first daughter" -----------------------------------------
    nCutIndex = 21;
    cutVariableValue = dd0pr1;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic second daughter" ----------------------------------------
    nCutIndex = 22;
    cutVariableValue = dd0pr2;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic max" ----------------------------------------------------
    nCutIndex = 23;
    cutVariableValue = dd0max;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------

    // "topomatic min" ----------------------------------------------------
    nCutIndex = 24;
    cutVariableValue = dd0min;
    bPassedCut = ApplyCutOnVariableDStarforDStarptbin(nCutIndex,ptbin,cutVariableValue,bCutArray);
    if(!bPassedCut && !fGetCutInfo) return 0;
    //---------------------------------------------------------------------
  }

  for (Int_t i = 0; i < 25; ++i)
  {
    if(bCutArray[i]==kTRUE){
      returnvalue = 0;
      break;
    }
  }
  
  return returnvalue;
}
//----------------------------------------------------------------------------------
Bool_t AliRDHFCutsB0toDStarPi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  // D* fiducial acceptance region // not (yet) used for the B0
  //

  // if(fMaxRapidityCand>-998.){
  //   if(TMath::Abs(y) > fMaxRapidityCand) return kFALSE;
  //   else return kTRUE;
  // }

  // if(pt > 5.) {
  //   // applying cut for pt > 5 GeV
  //   AliDebug(4,Form("pt of D* = %f (> 5), cutting at |y| < 0.8\n",pt)); 
  //   if (TMath::Abs(y) > 0.8){
  //     return kFALSE;
  //   }
  // } else {    
  //   // appliying smooth cut for pt < 5 GeV
  //   Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5; 
  //   Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;		
  //   AliDebug(2,Form("pt of D* = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY)); 
  //   if (y < minFiducialY || y > maxFiducialY){
  //     return kFALSE;
  //   }
  // }
    
  return kTRUE;
}
//_______________________________________________________________________________-
Int_t AliRDHFCutsB0toDStarPi::IsSelectedPID(AliAODRecoDecayHF* obj)
{
  //
  // PID method, n sigma approach default // not used for B0, done seperately for each daughter
  //
  
  // AliAODRecoCascadeHF* dstar = (AliAODRecoCascadeHF*)obj;
  // if(!dstar){
  //   cout<<"AliAODRecoCascadeHF null"<<endl;
  //   return 0;
  // } 
 
  // if(!fUsePID || dstar->Pt() > fMaxPtPid) return 3;
  
  // AliAODRecoDecayHF2Prong* d0 = (AliAODRecoDecayHF2Prong*)dstar->Get2Prong();  
  // if(!d0){
  //   cout<<"AliAODRecoDecayHF2Prong null"<<endl;
  //   return 0;
  // }

  // //  here the PID
  // AliAODTrack *pos = (AliAODTrack*)dstar->Get2Prong()->GetDaughter(0);
  // AliAODTrack *neg = (AliAODTrack*)dstar->Get2Prong()->GetDaughter(1);

  // if (dstar->Charge()>0){
  //   if(!SelectPID(pos,2)) return 0;//pion+
  //   if(!SelectPID(neg,3)) return 0;//kaon-
  // }else{
  //   if(!SelectPID(pos,3)) return 0;//kaon+
  //   if(!SelectPID(neg,2)) return 0;//pion-
  // }

  // if ((fPidHF->GetMatch() == 10 || fPidHF->GetMatch() == 11) && fPidHF->GetITS()) { //ITS n sigma band
  //   AliAODTrack *softPion = (AliAODTrack*) dstar->GetBachelor();

  //   if (fPidHF->CheckBands(AliPID::kPion, AliPIDResponse::kITS, softPion) == -1) {
  //     return 0;
  //   }
  // }

  return 3;
}
//_______________________________________________________________________________-
Int_t AliRDHFCutsB0toDStarPi::SelectPID(AliAODTrack *track, Int_t type)
{
  //
  //  here the PID
    
  Bool_t isParticle=kTRUE;
  Int_t match = fPidHF->GetMatch();

  if(match == 1){//n-sigma
    Bool_t TPCon=TMath::Abs(2)>1e-4?kTRUE:kFALSE;
    Bool_t TOFon=TMath::Abs(3)>1e-4?kTRUE:kFALSE;
    
    Bool_t isTPC=kTRUE;
    Bool_t isTOF=kTRUE;

    if (TPCon){//TPC
      if(fPidHF->CheckStatus(track,"TPC")){
      	if(type==2) isTPC=fPidHF->IsPionRaw(track,"TPC");
      	if(type==3) isTPC=fPidHF->IsKaonRaw(track,"TPC");
      }
    }
    if (TOFon){//TOF
      if(fPidHF->CheckStatus(track,"TOF")){
      	if(type==2) isTOF=fPidHF->IsPionRaw(track,"TOF");
      	if(type==3) isTOF=fPidHF->IsKaonRaw(track,"TOF");
      }
    }

    //--------------------------------
    // cut on high momentum in the TPC
    //--------------------------------
    Double_t pPIDcut = track->P();
    if(pPIDcut>fTPCflag) isTPC=1;
    
    isParticle = isTPC&&isTOF;
  }
  
  if(match == 2){//bayesian
    //Double_t priors[5]={0.01,0.001,0.3,0.3,0.3};
    Double_t prob[5]={1.,1.,1.,1.,1.};
    
    //fPidHF->SetPriors(priors,5);
    //    fPidHF->BayesianProbability(track,prob);
    
    Double_t max=0.;
    Int_t k=-1;
    for (Int_t i=0; i<5; i++) {
      if (prob[i]>max) {k=i; max=prob[i];}
    }
    isParticle = Bool_t(k==type);
  }

  if (match == 10 || match == 11) { //Assymetric PID using histograms
    Int_t checkTPC = fPidHF->CheckBands((AliPID::EParticleType) type, AliPIDResponse::kTPC, track);
    Int_t checkTOF = fPidHF->CheckBands((AliPID::EParticleType) type, AliPIDResponse::kTOF, track);

    isParticle = checkTPC >= 0 && checkTOF >= 0 ? kTRUE : kFALSE; //Standard: both at least compatible
    if (match == 11) { //Extra requirement: at least one identified
      isParticle = isParticle && checkTPC+checkTOF >= 1 ? kTRUE : kFALSE;
    }
  }
  
  if (match == 12) { //Circular cut
    Double_t nSigmaTOF = 0;
    Double_t nSigmaTPC = 0;

    Double_t radius = fCircRadius;

    isParticle = kTRUE;
    if (radius > 0) {
      Int_t TPCok = fPidHF->GetnSigmaTPC(track, type, nSigmaTPC);
      Int_t TOFok = fPidHF->GetnSigmaTOF(track, type, nSigmaTOF);
      if (TPCok != -1 && TOFok != -1) {
        //Both detectors gave PID information
        isParticle = TMath::Sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF) <= radius ? kTRUE : kFALSE;
      }
      else {
        //Maximum one detector gave PID information
        if (TPCok != -1) {
          isParticle = nSigmaTPC <= radius ? kTRUE : kFALSE;
        }
        if (TOFok != -1) {
          isParticle = nSigmaTOF <= radius ? kTRUE : kFALSE;
        }
      }
    }
  }

  return isParticle;
  
}
//-------------------------------------------------------------------------------------
Double_t AliRDHFCutsB0toDStarPi::DeltaInvMassDStarKpipi(AliAODRecoCascadeHF * DStar) const 
{
  ///
  /// 3 prong invariant mass of the D0 daughters and the soft pion
  ///
  Double_t e[3];
  e[0]=DStar->Get2Prong()->EProng(0,211);
  e[1]=DStar->Get2Prong()->EProng(1,321);
  e[2]=DStar->EProng(0,211);

  Double_t esum = e[0]+e[1]+e[2];
  Double_t invMassDStar = TMath::Sqrt(esum*esum-DStar->P()*DStar->P());

  Double_t invMassD0 = DStar->Get2Prong()->InvMassD0();

  return invMassDStar - invMassD0; 
}
//-------------------------------------------------------------------------------------
Double_t AliRDHFCutsB0toDStarPi::DeltaInvMassB0Kpipipi(AliAODRecoCascadeHF * B0) const 
{
  ///
  /// 4 prong invariant mass of the D0 daughters, the soft pion, and the B0 pion
  ///

  AliAODRecoCascadeHF * DStar = (AliAODRecoCascadeHF*)B0->GetDaughter(1);

  Double_t e[4];
  e[0]=DStar->EProng(0,211);
  e[1]=DStar->Get2Prong()->EProng(0,211);
  e[2]=DStar->Get2Prong()->EProng(1,321);
  e[3]=B0->EProng(0,211);
  // cout << "energy 1: " << e[0] << "energy 2: " << e[1] << "energy 3: " << e[2] << "energy 4: " << e[3] << endl;
  Double_t esum = e[0]+e[1]+e[2]+e[3];
  Double_t invMassB0 = TMath::Sqrt(esum*esum-B0->P()*B0->P());

  UInt_t prongs[2];
  prongs[0] = 211;
  prongs[1] = 421;
  // Double_t invMassDStar = DStar->InvMass(2,prongs);
  Double_t invMassD0 = DStar->Get2Prong()->InvMassD0();

  return invMassB0 - invMassD0; 
}


//---------------------------------------------------------------------------
//
//  DO for D0 pt bin functions
//
//---------------------------------------------------------------------------

void AliRDHFCutsB0toDStarPi::SetCutsD0forD0ptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDD0forD0ptbin) {
  //
  // store the cuts
  //

  if(nVars!=fnVarsD0forD0ptbin) {
    printf("Wrong number of variables: it has to be %d\n",fnVarsD0forD0ptbin);
    AliFatal("exiting");
  } 
  if(nPtBins!=fnPtBinsD0forD0ptbin) {
    printf("Wrong number of pt bins: it has to be %d\n",fnPtBinsD0forD0ptbin);
    AliFatal("exiting");
  } 

  if(!fCutsRDD0forD0ptbin)  fCutsRDD0forD0ptbin = new Float_t[fGlobalIndexD0forD0ptbin];


  for(Int_t iv=0; iv<fnVarsD0forD0ptbin; iv++) 
  {
    for(Int_t ib=0; ib<fnPtBinsD0forD0ptbin; ib++) 
    {
      //check

      if(GetGlobalIndexD0forD0ptbin(iv,ib)>=fGlobalIndexD0forD0ptbin) 
      {
        cout<<"Overflow, exit..."<<endl;
        return;
      }

      fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(iv,ib)] = cutsRDD0forD0ptbin[iv][ib];

    }
  }

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::PtBinD0forD0ptbin(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin=-1;
  if(pt<fPtBinLimitsD0forD0ptbin[0])return ptbin;
  for (Int_t i=0;i<fnPtBinsD0forD0ptbin;i++){
    if(pt<fPtBinLimitsD0forD0ptbin[i+1]) {
      ptbin=i;
      break;
    }
  }
  return ptbin;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetPtBinsD0forD0ptbin(Int_t nPtBinLimits,Float_t *ptBinLimits) {
  // Set the pt bins

  if(fPtBinLimitsD0forD0ptbin) {
    delete [] fPtBinLimitsD0forD0ptbin;
    fPtBinLimitsD0forD0ptbin = NULL;
    printf("Changing the pt bins\n");
  }

  if(nPtBinLimits != fnPtBinsD0forD0ptbin+1){
    cout<<"Warning: ptBinLimits dimension "<<nPtBinLimits<<" != nPtBins+1 ("<<fnPtBinsD0forD0ptbin+1<<")\nSetting nPtBins to "<<nPtBinLimits-1<<endl;
    SetNPtBinsD0forD0ptbin(nPtBinLimits-1);
  }

  fnPtBinLimitsD0forD0ptbin = nPtBinLimits;
  SetGlobalIndexD0forD0ptbin();
  //cout<<"Changing also Global Index -> "<<fGlobalIndex<<endl;
  fPtBinLimitsD0forD0ptbin = new Float_t[fnPtBinLimitsD0forD0ptbin];
  for(Int_t ib=0; ib<nPtBinLimits; ib++) fPtBinLimitsD0forD0ptbin[ib]=ptBinLimits[ib];
  for(Int_t ib=0; ib<nPtBinLimits; ib++) std::cout << "limit " << ib << " = " << fPtBinLimitsD0forD0ptbin[ib] << std::endl;

  return;
}



//---------------------------------------------------------------------------
//
//  DO for DStar pt bin functions
//
//---------------------------------------------------------------------------


void AliRDHFCutsB0toDStarPi::SetCutsD0forDStarptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDD0forDStarptbin) {
  //
  // store the cuts
  //
  if(nVars!=fnVarsD0forDStarptbin) {
    printf("Wrong number of variables: it has to be %d\n",fnVarsD0forDStarptbin);
    AliFatal("exiting");
  } 
  if(nPtBins!=fnPtBinsD0forDStarptbin) {
    printf("Wrong number of pt bins: it has to be %d\n",fnPtBinsD0forDStarptbin);
    AliFatal("exiting");
  } 

  if(!fCutsRDD0forDStarptbin)  fCutsRDD0forDStarptbin = new Float_t[fGlobalIndexD0forDStarptbin];
  

  for(Int_t iv=0; iv<fnVarsD0forDStarptbin; iv++) 
  {
    for(Int_t ib=0; ib<fnPtBinsD0forDStarptbin; ib++) 
    {
      //check
      if(GetGlobalIndexD0forDStarptbin(iv,ib)>=fGlobalIndexD0forDStarptbin) 
      {
        cout<<"Overflow, exit..."<<endl;
        return;
      }

      fCutsRDD0forDStarptbin[GetGlobalIndexD0forDStarptbin(iv,ib)] = cutsRDD0forDStarptbin[iv][ib];

    }
  }
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::PtBinD0forDStarptbin(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin=-1;
  if(pt<fPtBinLimitsD0forDStarptbin[0])return ptbin;
  for (Int_t i=0;i<fnPtBinsD0forDStarptbin;i++){
    if(pt<fPtBinLimitsD0forDStarptbin[i+1]) {
      ptbin=i;
      break;
    }
  }
  return ptbin;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetPtBinsD0forDStarptbin(Int_t nPtBinLimits,Float_t *ptBinLimits) {
  // Set the pt bins

  if(fPtBinLimitsD0forDStarptbin) {
    delete [] fPtBinLimitsD0forDStarptbin;
    fPtBinLimitsD0forDStarptbin = NULL;
    printf("Changing the pt bins\n");
  }

  if(nPtBinLimits != fnPtBinsD0forDStarptbin+1){
    cout<<"Warning: ptBinLimits dimension "<<nPtBinLimits<<" != nPtBins+1 ("<<fnPtBinsD0forDStarptbin+1<<")\nSetting nPtBins to "<<nPtBinLimits-1<<endl;
    SetNPtBinsD0forDStarptbin(nPtBinLimits-1);
  }

  fnPtBinLimitsD0forDStarptbin = nPtBinLimits;
  SetGlobalIndexD0forDStarptbin();
  //cout<<"Changing also Global Index -> "<<fGlobalIndex<<endl;
  fPtBinLimitsD0forDStarptbin = new Float_t[fnPtBinLimitsD0forDStarptbin];
  for(Int_t ib=0; ib<nPtBinLimits; ib++) fPtBinLimitsD0forDStarptbin[ib]=ptBinLimits[ib];

  return;
}




//---------------------------------------------------------------------------
//
//  DStar for DStar pt bin functions
//
//---------------------------------------------------------------------------



void AliRDHFCutsB0toDStarPi::SetCutsDStarforDStarptbin(Int_t nVars,Int_t nPtBins,Float_t **cutsRDDStarforDStarptbin) {
  //
  // store the cuts
  //
  if(nVars!=fnVarsDStarforDStarptbin) {
    printf("Wrong number of variables: it has to be %d\n",fnVarsDStarforDStarptbin);
    AliFatal("exiting");
  } 
  if(nPtBins!=fnPtBinsDStarforDStarptbin) {
    printf("Wrong number of pt bins: it has to be %d\n",fnPtBinsDStarforDStarptbin);
    AliFatal("exiting");
  } 

  if(!fCutsRDDStarforDStarptbin)  fCutsRDDStarforDStarptbin = new Float_t[fGlobalIndexDStarforDStarptbin];
  

  for(Int_t iv=0; iv<fnVarsDStarforDStarptbin; iv++) 
  {
    for(Int_t ib=0; ib<fnPtBinsDStarforDStarptbin; ib++) 
    {
      //check
      if(GetGlobalIndexDStarforDStarptbin(iv,ib)>=fGlobalIndexDStarforDStarptbin) 
      {
        cout<<"Overflow, exit..."<<endl;
        return;
      }

      fCutsRDDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(iv,ib)] = cutsRDDStarforDStarptbin[iv][ib];

    }
  }
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::PtBinDStarforDStarptbin(Double_t pt) const {
  //
  //give the pt bin where the pt lies.
  //
  Int_t ptbin=-1;
  if(pt<fPtBinLimitsDStarforDStarptbin[0])return ptbin;
  for (Int_t i=0;i<fnPtBinsDStarforDStarptbin;i++){
    if(pt<fPtBinLimitsDStarforDStarptbin[i+1]) {
      ptbin=i;
      break;
    }
  }
  return ptbin;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetPtBinsDStarforDStarptbin(Int_t nPtBinLimits,Float_t *ptBinLimits) {
  // Set the pt bins and initialize the cut values to 0 and use of cuts to kFALSE

  if(fPtBinLimitsDStarforDStarptbin) {
    delete [] fPtBinLimitsDStarforDStarptbin;
    fPtBinLimitsDStarforDStarptbin = NULL;
    printf("Changing the pt bins\n");
  }

  if(nPtBinLimits != fnPtBinsDStarforDStarptbin+1){
    cout<<"Warning: ptBinLimits dimension "<<nPtBinLimits<<" != nPtBins+1 ("<<fnPtBinsDStarforDStarptbin+1<<")\nSetting nPtBins to "<<nPtBinLimits-1<<endl;
    SetNPtBinsDStarforDStarptbin(nPtBinLimits-1);
  }

  fnPtBinLimitsDStarforDStarptbin = nPtBinLimits;
  SetGlobalIndexDStarforDStarptbin();
 
  fPtBinLimitsDStarforDStarptbin = new Float_t[fnPtBinLimitsDStarforDStarptbin];
  for(Int_t ib=0; ib<nPtBinLimits; ib++) fPtBinLimitsDStarforDStarptbin[ib]=ptBinLimits[ib];

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::ApplyCutOnVariable(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[85]){

  if(GetIsCutUsed(nCutIndex,ptbin)==kTRUE)
  {
    Bool_t bCut = kFALSE;
    if(GetIsUpperCut(nCutIndex)==kTRUE)
    {
      if(cutVariableValue > fCutsRD[GetGlobalIndex(nCutIndex,ptbin)]) bCut = kTRUE;
    } else
    {
      if(cutVariableValue < fCutsRD[GetGlobalIndex(nCutIndex,ptbin)]) bCut = kTRUE;
    }
    if(bCut == kTRUE) {bCutArray[nCutIndex] = 1; return 0;}
  }
  return 1; 
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::ApplyCutOnVariableD0forD0ptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[25]){

  // std::cout << "index: " << nCutIndex << ", ptbin: " << ptbin << ", used: " << GetIsCutUsedD0forD0ptbin(nCutIndex,ptbin) << ", upper: " << GetIsUpperCutD0forD0ptbin(nCutIndex) << ", value: " << cutVariableValue << ", cutvalue: " << fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)] << std::endl;
  if(GetIsCutUsedD0forD0ptbin(nCutIndex,ptbin)==kTRUE)
  {
    Bool_t bCut = kFALSE;
    if(GetIsUpperCutD0forD0ptbin(nCutIndex)==kTRUE)
    {
      if(cutVariableValue > fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)]) bCut = kTRUE;
    } else
    {
      if(cutVariableValue < fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptbin)]) bCut = kTRUE;
    }
    if(bCut == kTRUE) {bCutArray[nCutIndex] = 1; return 0;}
  }
  return 1; 
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::ApplyCutOnVariableD0forDStarptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[35]){

  if(GetIsCutUsedD0forDStarptbin(nCutIndex,ptbin)==kTRUE)
  {
    Bool_t bCut = kFALSE;
    if(GetIsUpperCutD0forDStarptbin(nCutIndex)==kTRUE)
    {
      if(cutVariableValue > fCutsRDD0forDStarptbin[GetGlobalIndexD0forDStarptbin(nCutIndex,ptbin)]) bCut = kTRUE;
    } else
    {
      if(cutVariableValue < fCutsRDD0forDStarptbin[GetGlobalIndexD0forDStarptbin(nCutIndex,ptbin)]) bCut = kTRUE;
    }
    if(bCut == kTRUE) {bCutArray[nCutIndex] = 1; return 0;}
  }
  return 1; 
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsB0toDStarPi::ApplyCutOnVariableDStarforDStarptbin(Int_t nCutIndex, Int_t ptbin, Float_t cutVariableValue, Bool_t bCutArray[25]){

  if(GetIsCutUsedDStarforDStarptbin(nCutIndex,ptbin)==kTRUE)
  {
    Bool_t bCut = kFALSE;
    if(GetIsUpperCutDStarforDStarptbin(nCutIndex)==kTRUE)
    {
      if(cutVariableValue > fCutsRDDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(nCutIndex,ptbin)]) bCut = kTRUE;
    } else
    {
      if(cutVariableValue < fCutsRDDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(nCutIndex,ptbin)]) bCut = kTRUE;
    }
    if(bCut == kTRUE) {bCutArray[nCutIndex] = 1; return 0;}
  }
  return 1; 
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetVarNamesD0forD0ptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut){
  // Set the variable names

  if(fVarNamesD0forD0ptbin) {
    delete [] fVarNamesD0forD0ptbin;
    fVarNamesD0forD0ptbin = NULL;
    //printf("Changing the variable names\n");
  }
  if(nVars!=fnVarsD0forD0ptbin){
    printf("Wrong number of variables: it has to be %d\n",fnVarsD0forD0ptbin);
    return;
  }
  //fnVars=nVars;
  fVarNamesD0forD0ptbin = new TString[nVars];
  fIsUpperCutD0forD0ptbin = new Bool_t[nVars];
  for(Int_t iv=0; iv<nVars; iv++) {
    fVarNamesD0forD0ptbin[iv] = varNames[iv];
    fIsUpperCutD0forD0ptbin[iv] = isUpperCut[iv];
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetVarNamesD0forDStarptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut){
  // Set the variable names

  if(fVarNamesD0forDStarptbin) {
    delete [] fVarNamesD0forDStarptbin;
    fVarNamesD0forDStarptbin = NULL;
    //printf("Changing the variable names\n");
  }
  if(nVars!=fnVarsD0forDStarptbin){
    printf("Wrong number of variables: it has to be %d\n",fnVarsD0forDStarptbin);
    return;
  }
  //fnVars=nVars;
  fVarNamesD0forDStarptbin = new TString[nVars];
  fIsUpperCutD0forDStarptbin = new Bool_t[nVars];
  for(Int_t iv=0; iv<nVars; iv++) {
    fVarNamesD0forDStarptbin[iv] = varNames[iv];
    fIsUpperCutD0forDStarptbin[iv] = isUpperCut[iv];
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetVarNamesDStarforDStarptbin(Int_t nVars,TString *varNames,Bool_t *isUpperCut){
  // Set the variable names

  if(fVarNamesDStarforDStarptbin) {
    delete [] fVarNamesDStarforDStarptbin;
    fVarNamesDStarforDStarptbin = NULL;
    //printf("Changing the variable names\n");
  }
  if(nVars!=fnVarsDStarforDStarptbin){
    printf("Wrong number of variables: it has to be %d\n",fnVarsDStarforDStarptbin);
    return;
  }
  //fnVars=nVars;
  fVarNamesDStarforDStarptbin = new TString[nVars];
  fIsUpperCutDStarforDStarptbin = new Bool_t[nVars];
  for(Int_t iv=0; iv<nVars; iv++) {
    fVarNamesDStarforDStarptbin[iv] = varNames[iv];
    fIsUpperCutDStarforDStarptbin[iv] = isUpperCut[iv];
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::InitializeCuts(){
  // Set the use of all cuts to kFALSE and the cut value to zero. This function has to be called after setting the pt bins.

  if(fIsCutUsed) {
    delete [] fIsCutUsed;
    fIsCutUsed = NULL;
  }

  fIsCutUsed = new Bool_t[(GetNPtBins())*(GetNVars())];
  for(Int_t iv=0; iv<(GetNPtBins())*(GetNVars()); iv++) {
    fIsCutUsed[iv] = kFALSE;
  }

  if(!fCutsRD)  fCutsRD = new Float_t[fGlobalIndex];
  

  for(Int_t iv=0; iv<fnVars; iv++) {

    for(Int_t ib=0; ib<fnPtBins; ib++) {

      //check
      if(GetGlobalIndex(iv,ib)>=fGlobalIndex) {
        cout<<"Overflow, exit..."<<endl;
        return;
      }

      fCutsRD[GetGlobalIndex(iv,ib)] = 0;

    }
  }

  // D0 for D0 pt bin

  if(fIsCutUsedD0forD0ptbin) {
    delete [] fIsCutUsedD0forD0ptbin;
    fIsCutUsedD0forD0ptbin = NULL;
  }

  fIsCutUsedD0forD0ptbin = new Bool_t[(GetNPtBinsD0forD0ptbin())*(GetNVarsD0forD0ptbin())];
  for(Int_t iv=0; iv<(GetNPtBinsD0forD0ptbin())*(GetNVarsD0forD0ptbin()); iv++) {
    fIsCutUsedD0forD0ptbin[iv] = kFALSE;

  }

  if(!fCutsRDD0forD0ptbin)  fCutsRDD0forD0ptbin = new Float_t[fGlobalIndexD0forD0ptbin];

  for(Int_t iv=0; iv<fnVarsD0forD0ptbin; iv++) 
  {
    for(Int_t ib=0; ib<fnPtBinsD0forD0ptbin; ib++) 
    {
      //check

      if(GetGlobalIndexD0forD0ptbin(iv,ib)>=fGlobalIndexD0forD0ptbin) 
      {
        cout<<"Overflow, exit..."<<endl;
        return;
      }

      fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(iv,ib)] = 0;
    }
  }

  // D0 for DStar pt bin

  if(fIsCutUsedD0forDStarptbin) {
    delete [] fIsCutUsedD0forDStarptbin;
    fIsCutUsedD0forDStarptbin = NULL;
  }

  fIsCutUsedD0forDStarptbin = new Bool_t[(GetNPtBinsD0forDStarptbin())*(GetNVarsD0forDStarptbin())];
  for(Int_t iv=0; iv<(GetNPtBinsD0forDStarptbin())*(GetNVarsD0forDStarptbin()); iv++) {
    fIsCutUsedD0forDStarptbin[iv] = kFALSE;
  }

  if(!fCutsRDD0forDStarptbin)  fCutsRDD0forDStarptbin = new Float_t[fGlobalIndexD0forDStarptbin];
  

  for(Int_t iv=0; iv<fnVarsD0forDStarptbin; iv++) 
  {
    for(Int_t ib=0; ib<fnPtBinsD0forDStarptbin; ib++) 
    {
      //check
      if(GetGlobalIndexD0forDStarptbin(iv,ib)>=fGlobalIndexD0forDStarptbin) 
      {
        cout<<"Overflow, exit..."<<endl;
        return;
      }

      fCutsRDD0forDStarptbin[GetGlobalIndexD0forDStarptbin(iv,ib)] = 0;

    }
  }

  // DStar for DStar pt bin

  if(fIsCutUsedDStarforDStarptbin) {
    delete [] fIsCutUsedDStarforDStarptbin;
    fIsCutUsedDStarforDStarptbin = NULL;
  }

  fIsCutUsedDStarforDStarptbin = new Bool_t[(GetNPtBinsDStarforDStarptbin())*(GetNVarsDStarforDStarptbin())];
  for(Int_t iv=0; iv<(GetNPtBinsDStarforDStarptbin())*(GetNVarsDStarforDStarptbin()); iv++) {
    fIsCutUsedDStarforDStarptbin[iv] = kFALSE;
  }

  if(!fCutsRDDStarforDStarptbin)  fCutsRDDStarforDStarptbin = new Float_t[fGlobalIndexDStarforDStarptbin];
  

  for(Int_t iv=0; iv<fnVarsDStarforDStarptbin; iv++) 
  {
    for(Int_t ib=0; ib<fnPtBinsDStarforDStarptbin; ib++) 
    {
      //check
      if(GetGlobalIndexDStarforDStarptbin(iv,ib)>=fGlobalIndexDStarforDStarptbin) 
      {
        cout<<"Overflow, exit..."<<endl;
        return;
      }

      fCutsRDDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(iv,ib)] = 0;

    }
  }

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetCut(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue){
  // Set the cut value and direction

  fIsCutUsed[GetGlobalIndex(nCutIndex,ptBin)] = kTRUE;
  fIsUpperCut[nCutIndex] = cutDirection;
  fCutsRD[GetGlobalIndex(nCutIndex,ptBin)] = cutValue;

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetCutD0forD0ptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue){
  // Set the cut value and direction

  fIsCutUsedD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptBin)] = kTRUE;
  fIsUpperCutD0forD0ptbin[nCutIndex] = cutDirection;
  fCutsRDD0forD0ptbin[GetGlobalIndexD0forD0ptbin(nCutIndex,ptBin)] = cutValue;

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetCutD0forDStarptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue){
  // Set the cut value and direction

  fIsCutUsedD0forDStarptbin[GetGlobalIndexD0forDStarptbin(nCutIndex,ptBin)] = kTRUE;
  fIsUpperCutD0forDStarptbin[nCutIndex] = cutDirection;
  fCutsRDD0forDStarptbin[GetGlobalIndexD0forDStarptbin(nCutIndex,ptBin)] = cutValue;

  return;
}
//---------------------------------------------------------------------------
void AliRDHFCutsB0toDStarPi::SetCutDStarforDStarptbin(Int_t nCutIndex, Int_t ptBin, AliRDHFCutsB0toDStarPi::EUpperCut cutDirection, Float_t cutValue){
  // Set the cut value and direction

  fIsCutUsedDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(nCutIndex,ptBin)] = kTRUE;
  fIsUpperCutDStarforDStarptbin[nCutIndex] = cutDirection;
  fCutsRDDStarforDStarptbin[GetGlobalIndexDStarforDStarptbin(nCutIndex,ptBin)] = cutValue;

  return;
}
//---------------------------------------------------------------------------
