// *************************************************************************
// *                                                                       *
// *               corrections to Fragmentation Functions                  *
// *                                                                       *
// *************************************************************************


/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: */

#include "TMath.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TDirectory.h"
#include "AliCFUnfolding.h"
#include "AliFragmentationFunctionCorrections.h"
#include <iostream> // OB TEST!!!
#include <fstream>

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  correction class for ouput of AliAnalysisTaskFragmentationFunction       //     
//                                                                           //
//  corrections for: reconstruction efficiency, momentum resolution,         //
//                   secondaries, UE / HI background, fluctuations           // 
//                   back-correction on jet energy on dN/dxi                 //
//                                                                           //
//  read MC ouput and write out efficiency histos, response matrices etc.    //
//  read measured distributions and apply efficiency, response matrices etc. //
//                                                                           //
//  contact: o.busch@gsi.de                                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


ClassImp(AliFragmentationFunctionCorrections)

//________________________________________________________________________
AliFragmentationFunctionCorrections::AliFragmentationFunctionCorrections()
   : TObject()
   ,fDebug(0)
   ,fNJetPtSlices(0)
   ,fJetPtSlices(0)   
   ,fNJets(0)
   ,fNJetsBgr(0)
   ,fNHistoBinsSinglePt(0)
   ,fHistoBinsSinglePt(0)
   ,fNHistoBinsPt(0)
   ,fNHistoBinsZ(0)
   ,fNHistoBinsXi(0)
   ,fHistoBinsPt(0)
   ,fHistoBinsZ(0)
   ,fHistoBinsXi(0)
   ,fNCorrectionLevels(0)
   ,fCorrFF(0)
   ,fNCorrectionLevelsBgr(0)
   ,fCorrBgr(0)
   ,fNCorrectionLevelsSinglePt(0)
   ,fCorrSinglePt(0)
   ,fh1FFXiShift(0)
   ,fh1EffSinglePt(0)
   ,fh1EffPt(0) 
   ,fh1EffZ(0) 
   ,fh1EffXi(0) 
   ,fh1EffBgrPt(0)
   ,fh1EffBgrZ(0)
   ,fh1EffBgrXi(0)
   ,fh1BbBCorrSinglePt(0) 
   ,fh1BbBPt(0)
   ,fh1BbBZ(0)
   ,fh1BbBXi(0)
   ,fh1BbBBgrPt(0)
   ,fh1BbBBgrZ(0)
   ,fh1BbBBgrXi(0)
   ,fh1FoldingCorrPt(0)
   ,fh1FoldingCorrZ(0)
   ,fh1FoldingCorrXi(0)
   ,fh1SecCorrPt(0)
   ,fh1SecCorrZ(0)
   ,fh1SecCorrXi(0)
   ,fh1SecCorrBgrPt(0)
   ,fh1SecCorrBgrZ(0)
   ,fh1SecCorrBgrXi(0)
   ,fh1FFTrackPtBackFolded(0)   
   ,fh1FFZBackFolded(0)
   ,fh1FFXiBackFolded(0)        
   ,fh1FFRatioTrackPtFolded(0)
   ,fh1FFRatioZFolded(0)
   ,fh1FFRatioXiFolded(0)
   ,fh1FFRatioTrackPtBackFolded(0)
   ,fh1FFRatioZBackFolded(0)
   ,fh1FFRatioXiBackFolded(0)
   ,fh1SingleTrackPtBackFolded(0)   
   ,fh1RatioSingleTrackPtFolded(0) 
   ,fh1RatioSingleTrackPtBackFolded(0)
   ,fhnResponseSinglePt(0)
   ,fhnResponsePt(0)
   ,fhnResponseZ(0) 
   ,fhnResponseXi(0)       
   ,fh1FFTrackPtPrior(0)
   ,fh1FFZPrior(0)
   ,fh1FFXiPrior(0)
   ,fh1SecCorrSinglePt(0)
{
   // default constructor
}

//________________________________________________________________________________________________________________________
AliFragmentationFunctionCorrections::AliFragmentationFunctionCorrections(const  AliFragmentationFunctionCorrections &copy)
  : TObject()
  ,fDebug(copy.fDebug)                        
  ,fNJetPtSlices(copy.fNJetPtSlices)                 
  ,fJetPtSlices(copy.fJetPtSlices)                  
  ,fNJets(copy.fNJets)                        
  ,fNJetsBgr(copy.fNJetsBgr)
  ,fNHistoBinsSinglePt(copy.fNHistoBinsSinglePt)  
  ,fHistoBinsSinglePt(copy.fHistoBinsSinglePt) 
  ,fNHistoBinsPt(copy.fNHistoBinsPt)                 
  ,fNHistoBinsZ(copy.fNHistoBinsZ)                  
  ,fNHistoBinsXi(copy.fNHistoBinsXi)                 
  ,fHistoBinsPt(copy.fHistoBinsPt)                  
  ,fHistoBinsZ(copy.fHistoBinsZ)                   
  ,fHistoBinsXi(copy.fHistoBinsXi)                  
  ,fNCorrectionLevels(copy.fNCorrectionLevels)            
  ,fCorrFF(copy.fCorrFF)                       
  ,fNCorrectionLevelsBgr(copy.fNCorrectionLevelsBgr)         
  ,fCorrBgr(copy.fCorrBgr)    
  ,fNCorrectionLevelsSinglePt(copy.fNCorrectionLevelsSinglePt)
  ,fCorrSinglePt(copy.fCorrSinglePt)   
  ,fh1FFXiShift(copy.fh1FFXiShift)                  
  ,fh1EffSinglePt(copy.fh1EffSinglePt)
  ,fh1EffPt(copy.fh1EffPt)                      
  ,fh1EffZ(copy.fh1EffZ)                       
  ,fh1EffXi(copy.fh1EffXi)                      
  ,fh1EffBgrPt(copy.fh1EffBgrPt)                   
  ,fh1EffBgrZ(copy.fh1EffBgrZ)                    
  ,fh1EffBgrXi(copy.fh1EffBgrXi)
  ,fh1BbBCorrSinglePt(copy.fh1BbBCorrSinglePt) 
  ,fh1BbBPt(copy.fh1BbBPt)
  ,fh1BbBZ(copy.fh1BbBZ)
  ,fh1BbBXi(copy.fh1BbBXi)
  ,fh1BbBBgrPt(copy.fh1BbBBgrPt)
  ,fh1BbBBgrZ(copy.fh1BbBBgrZ)
  ,fh1BbBBgrXi(copy.fh1BbBBgrXi)
  ,fh1FoldingCorrPt(copy.fh1FoldingCorrPt)
  ,fh1FoldingCorrZ(copy.fh1FoldingCorrZ)
  ,fh1FoldingCorrXi(copy.fh1FoldingCorrXi)
  ,fh1SecCorrPt(copy.fh1SecCorrPt)
  ,fh1SecCorrZ(copy.fh1SecCorrZ)
  ,fh1SecCorrXi(copy.fh1SecCorrXi)
  ,fh1SecCorrBgrPt(copy.fh1SecCorrBgrPt)
  ,fh1SecCorrBgrZ(copy.fh1SecCorrBgrZ) 
  ,fh1SecCorrBgrXi(copy.fh1SecCorrBgrXi)
  ,fh1FFTrackPtBackFolded(copy.fh1FFTrackPtBackFolded)        
  ,fh1FFZBackFolded(copy.fh1FFZBackFolded)              
  ,fh1FFXiBackFolded(copy.fh1FFXiBackFolded)             
  ,fh1FFRatioTrackPtFolded(copy.fh1FFRatioTrackPtFolded)       
  ,fh1FFRatioZFolded(copy.fh1FFRatioZFolded)             
  ,fh1FFRatioXiFolded(copy.fh1FFRatioXiFolded)            
  ,fh1FFRatioTrackPtBackFolded(copy.fh1FFRatioTrackPtBackFolded)   
  ,fh1FFRatioZBackFolded(copy.fh1FFRatioZBackFolded)         
  ,fh1FFRatioXiBackFolded(copy.fh1FFRatioXiBackFolded)        
  ,fh1SingleTrackPtBackFolded(copy.fh1SingleTrackPtBackFolded)     
  ,fh1RatioSingleTrackPtFolded(copy.fh1RatioSingleTrackPtFolded)    
  ,fh1RatioSingleTrackPtBackFolded(copy.fh1RatioSingleTrackPtBackFolded)
  ,fhnResponseSinglePt(copy.fhnResponseSinglePt)                 
  ,fhnResponsePt(copy.fhnResponsePt)                 
  ,fhnResponseZ(copy.fhnResponseZ)                  
  ,fhnResponseXi(copy.fhnResponseXi)                 
  ,fh1FFTrackPtPrior(copy.fh1FFTrackPtPrior)   
  ,fh1FFZPrior(copy.fh1FFZPrior)
  ,fh1FFXiPrior(copy.fh1FFXiPrior)
  ,fh1SecCorrSinglePt(copy.fh1SecCorrSinglePt)
{
  // copy constructor
  
}

// ______________________________________________________________________________________________________________________________
AliFragmentationFunctionCorrections& AliFragmentationFunctionCorrections::operator=(const AliFragmentationFunctionCorrections& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    fDebug                          = o.fDebug;                        
    fNJetPtSlices                   = o.fNJetPtSlices;                 
    fJetPtSlices                    = o.fJetPtSlices;                  
    fNJets                          = o.fNJets;                        
    fNJetsBgr                       = o.fNJetsBgr;                     
    fNHistoBinsSinglePt             = o.fNHistoBinsSinglePt;  
    fHistoBinsSinglePt              = o.fHistoBinsSinglePt; 
    fNHistoBinsPt                   = o.fNHistoBinsPt;                 
    fNHistoBinsZ                    = o.fNHistoBinsZ;                  
    fNHistoBinsXi                   = o.fNHistoBinsXi;                 
    fHistoBinsPt                    = o.fHistoBinsPt;                  
    fHistoBinsZ                     = o.fHistoBinsZ;                   
    fHistoBinsXi                    = o.fHistoBinsXi;                  
    fNCorrectionLevels              = o.fNCorrectionLevels;            
    fCorrFF                         = o.fCorrFF;                       
    fNCorrectionLevelsBgr           = o.fNCorrectionLevelsBgr;         
    fCorrBgr                        = o.fCorrBgr;                      
    fNCorrectionLevelsSinglePt      = o.fNCorrectionLevelsSinglePt;
    fCorrSinglePt                   = o.fCorrSinglePt;
    fh1FFXiShift                    = o.fh1FFXiShift;                  
    fh1EffSinglePt                  = o.fh1EffSinglePt;
    fh1EffPt                        = o.fh1EffPt;                      
    fh1EffZ                         = o.fh1EffZ;                       
    fh1EffXi                        = o.fh1EffXi;                      
    fh1EffBgrPt                     = o.fh1EffBgrPt;                   
    fh1EffBgrZ                      = o.fh1EffBgrZ;                    
    fh1EffBgrXi                     = o.fh1EffBgrXi;                   
    fh1BbBCorrSinglePt              = o.fh1BbBCorrSinglePt; 
    fh1BbBPt                        = o.fh1BbBPt;
    fh1BbBZ                         = o.fh1BbBZ;
    fh1BbBXi                        = o.fh1BbBXi;
    fh1BbBBgrPt                     = o.fh1BbBBgrPt;
    fh1BbBBgrZ                      = o.fh1BbBBgrZ;
    fh1BbBBgrXi                     = o.fh1BbBBgrXi;
    fh1FoldingCorrPt                = o.fh1FoldingCorrPt;
    fh1FoldingCorrZ                 = o.fh1FoldingCorrZ;
    fh1FoldingCorrXi                = o.fh1FoldingCorrXi;
    fh1SecCorrPt                    = o.fh1SecCorrPt;
    fh1SecCorrZ                     = o.fh1SecCorrZ;
    fh1SecCorrXi                    = o.fh1SecCorrXi;
    fh1SecCorrBgrPt                 = o.fh1SecCorrBgrPt;
    fh1SecCorrBgrZ                  = o.fh1SecCorrBgrZ;
    fh1SecCorrBgrXi                 = o.fh1SecCorrBgrXi;
    fh1FFTrackPtBackFolded          = o.fh1FFTrackPtBackFolded;        
    fh1FFZBackFolded                = o.fh1FFZBackFolded;              
    fh1FFXiBackFolded               = o.fh1FFXiBackFolded;             
    fh1FFRatioTrackPtFolded         = o.fh1FFRatioTrackPtFolded;       
    fh1FFRatioZFolded               = o.fh1FFRatioZFolded;             
    fh1FFRatioXiFolded              = o.fh1FFRatioXiFolded;            
    fh1FFRatioTrackPtBackFolded     = o.fh1FFRatioTrackPtBackFolded;   
    fh1FFRatioZBackFolded           = o.fh1FFRatioZBackFolded;         
    fh1FFRatioXiBackFolded          = o.fh1FFRatioXiBackFolded;        
    fh1SingleTrackPtBackFolded      = o.fh1SingleTrackPtBackFolded;     
    fh1RatioSingleTrackPtFolded     = o.fh1RatioSingleTrackPtFolded;    
    fh1RatioSingleTrackPtBackFolded = o.fh1RatioSingleTrackPtBackFolded;
    fhnResponseSinglePt             = o.fhnResponseSinglePt;                 
    fhnResponsePt                   = o.fhnResponsePt;                 
    fhnResponseZ                    = o.fhnResponseZ;                  
    fhnResponseXi                   = o.fhnResponseXi;                 
    fh1FFTrackPtPrior               = o.fh1FFTrackPtPrior;
    fh1FFZPrior                     = o.fh1FFZPrior;
    fh1FFXiPrior                    = o.fh1FFXiPrior;
    fh1SecCorrSinglePt              = o.fh1SecCorrSinglePt;
  }
  
  return *this;
}

//_________________________________________________________________________
AliFragmentationFunctionCorrections::~AliFragmentationFunctionCorrections()
{
  // destructor  

  if(fJetPtSlices) delete fJetPtSlices;
  if(fNJets)       delete fNJets;
  if(fNJetsBgr)    delete fNJetsBgr;

  DeleteHistoArray(fh1FFXiShift);

  DeleteHistoArray(fh1EffPt);
  DeleteHistoArray(fh1EffXi);
  DeleteHistoArray(fh1EffZ );

  DeleteHistoArray(fh1EffBgrPt);
  DeleteHistoArray(fh1EffBgrXi);
  DeleteHistoArray(fh1EffBgrZ);

  // BbB
  DeleteHistoArray(fh1BbBPt);
  DeleteHistoArray(fh1BbBXi);
  DeleteHistoArray(fh1BbBZ);

  DeleteHistoArray(fh1BbBBgrPt);
  DeleteHistoArray(fh1BbBBgrXi);
  DeleteHistoArray(fh1BbBBgrZ);

  DeleteHistoArray(fh1FoldingCorrPt);
  DeleteHistoArray(fh1FoldingCorrXi);
  DeleteHistoArray(fh1FoldingCorrZ);

  // sec corr
  DeleteHistoArray(fh1SecCorrPt);
  DeleteHistoArray(fh1SecCorrXi);
  DeleteHistoArray(fh1SecCorrZ);

  DeleteHistoArray(fh1SecCorrBgrPt);
  DeleteHistoArray(fh1SecCorrBgrXi);
  DeleteHistoArray(fh1SecCorrBgrZ);

  // unfolding

  DeleteHistoArray(fh1FFTrackPtBackFolded);
  DeleteHistoArray(fh1FFZBackFolded);
  DeleteHistoArray(fh1FFXiBackFolded);

  DeleteHistoArray(fh1FFRatioTrackPtFolded);
  DeleteHistoArray(fh1FFRatioZFolded);
  DeleteHistoArray(fh1FFRatioXiFolded);

  DeleteHistoArray(fh1FFRatioTrackPtBackFolded);
  DeleteHistoArray(fh1FFRatioZBackFolded);
  DeleteHistoArray(fh1FFRatioXiBackFolded);

  DeleteTHnSparseArray(fhnResponsePt);
  DeleteTHnSparseArray(fhnResponseZ);
  DeleteTHnSparseArray(fhnResponseXi);

  DeleteHistoArray(fh1FFTrackPtPrior);
  DeleteHistoArray(fh1FFZPrior);
  DeleteHistoArray(fh1FFXiPrior);

  // clean up corrected FF 
  
  for(Int_t c=0; c<fNCorrectionLevels; c++) delete  fCorrFF[c];
  delete[] fCorrFF;

  // clean up bgr

  for(Int_t c=0; c<fNCorrectionLevelsBgr; c++) delete fCorrBgr[c];
  delete[] fCorrBgr;

  // clean up inclusive pt 
  for(Int_t c=0; c<fNCorrectionLevelsSinglePt; c++) delete fCorrSinglePt[c];
  delete[] fCorrSinglePt;

  delete[] fNHistoBinsPt;
  delete[] fNHistoBinsZ;
  delete[] fNHistoBinsXi;

  // clean up histo bins

  if(fHistoBinsSinglePt) delete fHistoBinsSinglePt; 

  for(Int_t i=0; i<fNJetPtSlices; i++) delete fHistoBinsPt[i];
  for(Int_t i=0; i<fNJetPtSlices; i++) delete fHistoBinsZ[i];
  for(Int_t i=0; i<fNJetPtSlices; i++) delete fHistoBinsXi[i];

  delete[] fHistoBinsPt;
  delete[] fHistoBinsZ;
  delete[] fHistoBinsXi;
}

//_________________________________________________________________________________
AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::AliFragFuncCorrHistos()
  : TObject()
  ,fArraySize(0)
  ,fh1CorrFFTrackPt(0)
  ,fh1CorrFFZ(0)
  ,fh1CorrFFXi(0)
  ,fCorrLabel(0)
{
  // default constructor
  
}

//__________________________________________________________________________________________________________________
AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::AliFragFuncCorrHistos(const AliFragmentationFunctionCorrections::AliFragFuncCorrHistos& copy)
  : TObject()
  ,fArraySize(copy.fArraySize)
  ,fh1CorrFFTrackPt(0)
  ,fh1CorrFFZ(0)
  ,fh1CorrFFXi(0)
  ,fCorrLabel(copy.fCorrLabel)
{
  // copy constructor

  fh1CorrFFTrackPt = new TH1F*[copy.fArraySize];
  fh1CorrFFZ       = new TH1F*[copy.fArraySize];
  fh1CorrFFXi      = new TH1F*[copy.fArraySize];

  for(Int_t i=0; i<copy.fArraySize; i++){
    fh1CorrFFTrackPt[i] = new TH1F(*(copy.fh1CorrFFTrackPt[i]));
    fh1CorrFFZ[i]       = new TH1F(*(copy.fh1CorrFFZ[i]));
    fh1CorrFFXi[i]      = new TH1F(*(copy.fh1CorrFFXi[i]));
  }
}

//_______________________________________________________________________________________________________________________________________________________________
AliFragmentationFunctionCorrections::AliFragFuncCorrHistos& AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::operator=(const AliFragmentationFunctionCorrections::AliFragFuncCorrHistos& o)
{
  // assignment
  
  if(this!=&o){
    TObject::operator=(o);
    Int_t fArraySize_new = o.fArraySize;
    fCorrLabel       = o.fCorrLabel;
 
    TH1F** fh1CorrFFTrackPt_new = new TH1F*[fArraySize_new];
    TH1F** fh1CorrFFZ_new       = new TH1F*[fArraySize_new];
    TH1F** fh1CorrFFXi_new      = new TH1F*[fArraySize_new];

    for(Int_t i=0; i<fArraySize_new; i++){
      fh1CorrFFTrackPt_new[i] = new TH1F(*(o.fh1CorrFFTrackPt[i]));
      fh1CorrFFZ_new[i]       = new TH1F(*(o.fh1CorrFFZ[i]));
      fh1CorrFFXi_new[i]      = new TH1F(*(o.fh1CorrFFXi[i]));
    }
    
    // ---

    if(fArraySize){
      for(int i=0; i<fArraySize; i++) delete fh1CorrFFTrackPt[i];
      for(int i=0; i<fArraySize; i++) delete fh1CorrFFZ[i];
      for(int i=0; i<fArraySize; i++) delete fh1CorrFFXi[i];          
    }

    if(fh1CorrFFTrackPt) delete[] fh1CorrFFTrackPt;
    if(fh1CorrFFZ)       delete[] fh1CorrFFZ;
    if(fh1CorrFFXi)      delete[] fh1CorrFFXi;
        
    // ---

    fArraySize = fArraySize_new;

    fh1CorrFFTrackPt = fh1CorrFFTrackPt_new;
    fh1CorrFFZ       = fh1CorrFFZ_new;
    fh1CorrFFXi      = fh1CorrFFXi_new;
    
    for(Int_t i=0; i<fArraySize; i++){
      fh1CorrFFTrackPt[i] = fh1CorrFFTrackPt_new[i];
      fh1CorrFFZ[i]       = fh1CorrFFZ_new[i];
      fh1CorrFFXi[i]      = fh1CorrFFXi_new[i];
    }
  }
  
  return *this;
}

//__________________________________________________________________________________
AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::~AliFragFuncCorrHistos()
{
  // destructor 

  if(fArraySize){
    for(int i=0; i<fArraySize; i++) delete fh1CorrFFTrackPt[i];
    for(int i=0; i<fArraySize; i++) delete fh1CorrFFZ[i];
    for(int i=0; i<fArraySize; i++) delete fh1CorrFFXi[i];          
  }

  if(fh1CorrFFTrackPt) delete[] fh1CorrFFTrackPt;
  if(fh1CorrFFZ)       delete[] fh1CorrFFZ;
  if(fh1CorrFFXi)      delete[] fh1CorrFFXi;

}

//___________________________________________________________________________________________________________________
AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::AliFragFuncCorrHistos(const char* label, Int_t arraySize)
  : TObject()
  ,fArraySize(0)
  ,fh1CorrFFTrackPt(0)
  ,fh1CorrFFZ(0)
  ,fh1CorrFFXi(0)
  ,fCorrLabel(label)
{
  // constructor

  fArraySize = arraySize;
  fh1CorrFFTrackPt = new TH1F*[arraySize];
  fh1CorrFFZ       = new TH1F*[arraySize];
  fh1CorrFFXi      = new TH1F*[arraySize];

  for(int i=0; i<arraySize; i++) fh1CorrFFTrackPt[i] = new TH1F;
  for(int i=0; i<arraySize; i++) fh1CorrFFZ[i]       = new TH1F;
  for(int i=0; i<arraySize; i++) fh1CorrFFXi[i]      = new TH1F;
}

//_______________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::AddCorrHistos(Int_t slice,TH1F* histPt,TH1F* histZ,TH1F* histXi)
{
  // place histo in array, add corrLabel to histo name

  if(slice>=fArraySize){
    Printf("%s:%d -- slice > array size", (char*)__FILE__,__LINE__);
    return;
  }

  if(histPt){
    TString name = histPt->GetName();
    if(fCorrLabel.Length()>0) name += "_"+fCorrLabel;
    histPt->SetName(name);
    
    if(!(histPt->GetSumw2())) histPt->Sumw2(); 

    new(fh1CorrFFTrackPt[slice]) TH1F(*histPt);
  }
  
  if(histZ){
    TString name = histZ->GetName();
    if(fCorrLabel.Length()>0) name += "_"+fCorrLabel;
    histZ->SetName(name);

    if(!(histZ->GetSumw2())) histZ->Sumw2(); 

    new(fh1CorrFFZ[slice]) TH1F(*histZ);
  }

  if(histXi){
    TString name = histXi->GetName();
    if(fCorrLabel.Length()>0) name += "_"+fCorrLabel;
    histXi->SetName(name);

    if(!(histXi->GetSumw2())) histXi->Sumw2(); 

    new(fh1CorrFFXi[slice]) TH1F(*histXi);
  }
}

//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::ReplaceCorrHistos(Int_t slice,TH1F* histPt,TH1F* histZ,TH1F* histXi)
{
  // replace histo in array 

  if(slice>=fArraySize){
    Printf("%s:%d -- slice > array size", (char*)__FILE__,__LINE__);
    return;
  }

  if(histPt){
    if(!(histPt->GetSumw2())) histPt->Sumw2(); 

    TString name = histPt->GetName();
    histPt->SetName(name);

    delete fh1CorrFFTrackPt[slice];
    fh1CorrFFTrackPt[slice] =  new TH1F;
    new(fh1CorrFFTrackPt[slice]) TH1F(*histPt);
  }
  
  if(histZ){
    if(!(histZ->GetSumw2())) histZ->Sumw2(); 

    TString name = histZ->GetName();
    histZ->SetName(name);

    delete fh1CorrFFZ[slice];
    fh1CorrFFZ[slice] = new TH1F;
    new(fh1CorrFFZ[slice]) TH1F(*histZ);
  }

  if(histXi){
    if(!(histXi->GetSumw2())) histXi->Sumw2(); 
    
    TString name = histXi->GetName();
    histXi->SetName(name);

    delete fh1CorrFFXi[slice];
    fh1CorrFFXi[slice] = new TH1F;
    new(fh1CorrFFXi[slice]) TH1F(*histXi);
  }
}

// ___________________________________________________________________________________________
TH1F* AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::GetTrackPt(const Int_t slice)
{ 
  // return pt histo 
  
  if(slice>=fArraySize){
    Printf("%s:%d -- slice > array size", (char*)__FILE__,__LINE__);
    return 0;
  }

  return fh1CorrFFTrackPt[slice]; 

}

// ______________________________________________________________________________________
TH1F* AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::GetZ(const Int_t slice)
{ 
  // return z histo 
  
  if(slice>=fArraySize){
    Printf("%s:%d -- slice > array size", (char*)__FILE__,__LINE__);
    return 0;
  }

  return fh1CorrFFZ[slice]; 
}

// ________________________________________________________________________________________
TH1F* AliFragmentationFunctionCorrections::AliFragFuncCorrHistos::GetXi(const Int_t slice)
{ 
  // return xi histo

  if(slice>=fArraySize){
    Printf("%s:%d -- slice > array size", (char*)__FILE__,__LINE__);
    return 0;
  }

  return fh1CorrFFXi[slice]; 
}

// __________________________________________________________________________
void AliFragmentationFunctionCorrections::DeleteHistoArray(TH1F** hist) const
{
  // delete array of TH1 
 
  if(hist){
    for(Int_t i=0; i<fNJetPtSlices; i++) delete hist[i];
    delete[] hist;
  }
}

// ____________________________________________________________________________________
void AliFragmentationFunctionCorrections::DeleteTHnSparseArray(THnSparse** hist) const
{
  // delete array of THnSparse 
 
  if(hist){
    for(Int_t i=0; i<fNJetPtSlices; i++) delete hist[i];
    delete[] hist;
  }
}

// _________________________________________________________
TH1F** AliFragmentationFunctionCorrections::BookHistoArray()
{
  // book array of TH1

  if(!fNJetPtSlices){
    Printf("%s:%d -- jetPtSlices not defined", (char*)__FILE__,__LINE__);
    return 0;
  }

  TH1F** hist = new TH1F*[fNJetPtSlices];
  for(Int_t i=0; i<fNJetPtSlices; i++) hist[i] = new TH1F();
  
  return hist;
}

//__________________________________________________________________
THnSparse** AliFragmentationFunctionCorrections::BookTHnSparseArray()
{
  // book array of THnSparse

  if(!fNJetPtSlices){
    Printf("%s:%d -- jetPtSlices not defined", (char*)__FILE__,__LINE__);
    return 0;
  }

  THnSparse** hist = new THnSparse*[fNJetPtSlices];
  for(Int_t i=0; i<fNJetPtSlices; i++) hist[i] = new THnSparseF(); 
  
  return hist;
}

//_____________________________________________________________________________
void AliFragmentationFunctionCorrections::AddCorrectionLevel(const char* label)
{
  // increase corr level 

  if(!fNJetPtSlices){
    Printf("%s:%d -- jetPtSlices not defined", (char*)__FILE__,__LINE__);
    return;
  }

  if(fNCorrectionLevels >= fgMaxNCorrectionLevels){
    Printf("%s:%d -- max correction level exceeded", (char*)__FILE__,__LINE__);
    return;
  }

  fCorrFF[fNCorrectionLevels] = new AliFragFuncCorrHistos(label,fNJetPtSlices);
  fNCorrectionLevels++;
}


//________________________________________________________________________________
void AliFragmentationFunctionCorrections::AddCorrectionLevelBgr(const char* label)
{
  // increase corr level for bgr FF

  if(!fNJetPtSlices){
    if(fDebug>0)  Printf("%s:%d -- jetPtSlices not defined", (char*)__FILE__,__LINE__);
    return;
  }
  
  if(fNCorrectionLevelsBgr >= fgMaxNCorrectionLevels){
    Printf("%s:%d -- max correction level exceeded", (char*)__FILE__,__LINE__);
    return;
  }

  fCorrBgr[fNCorrectionLevelsBgr] = new AliFragFuncCorrHistos(label,fNJetPtSlices);
  fNCorrectionLevelsBgr++;
}

//_____________________________________________________________________________________
void AliFragmentationFunctionCorrections::AddCorrectionLevelSinglePt(const char* label)
{
  // increase corr level single pt spec

  Int_t nJetPtSlicesSingle = 1; // pro forma 

  if(fNCorrectionLevelsSinglePt >= fgMaxNCorrectionLevels){
    Printf("%s:%d -- max correction level exceeded", (char*)__FILE__,__LINE__);
    return;
  }

  fCorrSinglePt[fNCorrectionLevelsSinglePt] = new AliFragFuncCorrHistos(label,nJetPtSlicesSingle);
  fNCorrectionLevelsSinglePt++;
}

//_____________________________________________________________________________________________
void AliFragmentationFunctionCorrections::SetJetPtSlices(Float_t* bins, const Int_t nJetPtSlices)
{ 
  // set jet pt slices
  // once slices are known, can book all other histos 

  fNJetPtSlices = nJetPtSlices;

  const Int_t size = nJetPtSlices+1;
  fJetPtSlices  = new TArrayF(size,bins);
 
  // nJets array

  fNJets = new TArrayF(size);
  fNJets->Reset(0);

  fNJetsBgr = new TArrayF(size);
  fNJetsBgr->Reset(0);

  // histos 

  fNCorrectionLevels = 0; 
  fCorrFF = new AliFragFuncCorrHistos*[fgMaxNCorrectionLevels];
  AddCorrectionLevel(); // first 'correction' level = raw FF

  fNCorrectionLevelsBgr = 0; 
  fCorrBgr = new AliFragFuncCorrHistos*[fgMaxNCorrectionLevels];
  AddCorrectionLevelBgr(); // first 'correction' level = raw bgr dist

  fh1FFXiShift = BookHistoArray();

  // eff histos

  fh1EffPt = BookHistoArray();
  fh1EffXi = BookHistoArray();
  fh1EffZ  = BookHistoArray();

  fh1EffBgrPt = BookHistoArray();
  fh1EffBgrXi = BookHistoArray();
  fh1EffBgrZ  = BookHistoArray();


  // BbB

  fh1BbBPt = BookHistoArray();
  fh1BbBXi = BookHistoArray();
  fh1BbBZ  = BookHistoArray();

  fh1BbBBgrPt = BookHistoArray();
  fh1BbBBgrXi = BookHistoArray();
  fh1BbBBgrZ  = BookHistoArray();

  fh1FoldingCorrPt = BookHistoArray();
  fh1FoldingCorrXi = BookHistoArray();
  fh1FoldingCorrZ  = BookHistoArray();

  // SecCorr

  fh1SecCorrPt = BookHistoArray();
  fh1SecCorrXi = BookHistoArray();
  fh1SecCorrZ  = BookHistoArray();

  fh1SecCorrBgrPt = BookHistoArray();
  fh1SecCorrBgrXi = BookHistoArray();
  fh1SecCorrBgrZ  = BookHistoArray();

  // unfolding

  fh1FFTrackPtBackFolded = BookHistoArray();
  fh1FFXiBackFolded      = BookHistoArray();
  fh1FFZBackFolded       = BookHistoArray();

  fh1FFRatioTrackPtFolded = BookHistoArray();
  fh1FFRatioXiFolded      = BookHistoArray();
  fh1FFRatioZFolded       = BookHistoArray();
  
  fh1FFRatioTrackPtBackFolded = BookHistoArray();
  fh1FFRatioXiBackFolded      = BookHistoArray();
  fh1FFRatioZBackFolded       = BookHistoArray();

  //

  fhnResponsePt      = BookTHnSparseArray(); 
  fhnResponseZ       = BookTHnSparseArray();
  fhnResponseXi      = BookTHnSparseArray();

  fh1FFTrackPtPrior =  BookHistoArray();
  fh1FFZPrior       =  BookHistoArray();
  fh1FFXiPrior      =  BookHistoArray();

  // histos bins 

  fNHistoBinsPt = new Int_t[fNJetPtSlices];
  fNHistoBinsXi = new Int_t[fNJetPtSlices];
  fNHistoBinsZ  = new Int_t[fNJetPtSlices];

  for(Int_t i=0; i<fNJetPtSlices; i++) fNHistoBinsPt[i] = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) fNHistoBinsXi[i] = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) fNHistoBinsZ[i]  = 0;
  
  fHistoBinsPt = new TArrayD*[fNJetPtSlices];
  fHistoBinsXi = new TArrayD*[fNJetPtSlices];
  fHistoBinsZ  = new TArrayD*[fNJetPtSlices];

  for(Int_t i=0; i<fNJetPtSlices; i++) fHistoBinsPt[i] = new TArrayD(0);
  for(Int_t i=0; i<fNJetPtSlices; i++) fHistoBinsXi[i] = new TArrayD(0);
  for(Int_t i=0; i<fNJetPtSlices; i++) fHistoBinsZ[i]  = new TArrayD(0);
}

//_____________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::SetHistoBins(const Int_t jetPtSlice, const Int_t sizeBins, Double_t* bins, const Int_t type)
{ 
  // set histo bins for jet pt slice
  // if binning undefined for any slice, original binning will be used

  if(!fNJetPtSlices){
    Printf("%s:%d -- jetPtSlices not defined", (char*)__FILE__,__LINE__);
    return;
  }
  
  if(jetPtSlice>=fNJetPtSlices){
    Printf("%s:%d -- jetPtSlice %d exceeds max",(char*)__FILE__,__LINE__,jetPtSlice);
    return;
  }
  
  if(type == kFlagPt){
    fNHistoBinsPt[jetPtSlice] = sizeBins-1;
    fHistoBinsPt[jetPtSlice]->Set(sizeBins,bins); 
  }
  else if(type==kFlagZ){
    fNHistoBinsZ[jetPtSlice] = sizeBins-1;
    fHistoBinsZ[jetPtSlice]->Set(sizeBins,bins); 
  }

  else if(type==kFlagXi){
    fNHistoBinsXi[jetPtSlice] = sizeBins-1;
    fHistoBinsXi[jetPtSlice]->Set(sizeBins,bins); 
  }
}

//__________________________________________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::SetHistoBins(const Int_t jetPtSlice, const Int_t nBinsLimits, Double_t* binsLimits, Double_t* binsWidth, const Int_t type)
{ 
  // set histo bins for jet pt slice
  // function takes array of limits and widths (e.g. 1 GeV bins up to 10 GeV, 2 GeV width up to 20 GeV, ...)  
  // array size of binsLimits: nBinsLimits 
  // array size of binsWidth: nBinsLimits-1 
  // binsLimits have to be in increasing order
  // if binning undefined for any slice, original binning will be used

  if(!fNJetPtSlices){
    Printf("%s:%d -- jetPtSlices not defined", (char*)__FILE__,__LINE__);
    return;
  }

  if(jetPtSlice>=fNJetPtSlices){
    Printf("%s:%d -- jetPtSlice %d exceeds max",(char*)__FILE__,__LINE__,jetPtSlice);
    return;
  }


  Double_t binLimitMin = binsLimits[0];
  Double_t binLimitMax = binsLimits[nBinsLimits-1];

  Double_t binLimit = binLimitMin; // start value 
  
  Int_t sizeUpperLim = 10000; //static_cast<Int_t>(binLimitMax/binsWidth[0])+1;
  TArrayD binsArray(sizeUpperLim);
  Int_t nBins = 0; 
  binsArray.SetAt(binLimitMin,nBins++);

  while(binLimit<0.999999*binLimitMax && nBins<sizeUpperLim){ // 0.999 numerical safety factor

    Int_t currentSlice = -1;
    for(Int_t i=0; i<nBinsLimits; i++){
      if(binLimit >= 0.999999*binsLimits[i]) currentSlice = i;  
    }
    
    Double_t currentBinWidth = binsWidth[currentSlice];
    binLimit += currentBinWidth;

    // if(type == kFlagZ) std::cout<<" SetHistoBins: nBins "<<nBins<<" binLimit "<<binLimit<<" binLimitMax "<<binLimitMax
    // 				<<" currentSlice "<<currentSlice<<std::endl;

    binsArray.SetAt(binLimit,nBins++);
  }
  
  Double_t* bins = binsArray.GetArray();

  SetHistoBins(jetPtSlice,nBins,bins,type); 
}

//_____________________________________________________________________________________________________________________________________
TArrayD* AliFragmentationFunctionCorrections::GetHistoBins(const Int_t jetPtSlice,  const Int_t type)
{ 
  // set histo bins for jet pt slice
  // if binning undefined for any slice, original binning will be used

  if(!fNJetPtSlices){
    Printf("%s:%d -- jetPtSlices not defined", (char*)__FILE__,__LINE__);
    return 0;
  }
  
  if(jetPtSlice>=fNJetPtSlices){
    Printf("%s:%d -- jetPtSlice %d exceeds max",(char*)__FILE__,__LINE__,jetPtSlice);
    return 0;
  }
  
  if(type == kFlagPt){
    return fHistoBinsPt[jetPtSlice]; 
  }
  else if(type==kFlagZ){
    return fHistoBinsZ[jetPtSlice]; 
  }
  
  else if(type==kFlagXi){
    return fHistoBinsXi[jetPtSlice]; 
  }
  else{ 
    Printf("%s%d unknown type",(char*)__FILE__,__LINE__);
    return 0;
  }
}

//__________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::SetHistoBinsSinglePt(const Int_t sizeBins, Double_t* bins)
{ 
  // set histo bins for inclusive pt spectra
  // if binning undefined, original binning will be used

  fNHistoBinsSinglePt = sizeBins-1;
  
  fHistoBinsSinglePt = new TArrayD(sizeBins);
  fHistoBinsSinglePt->Set(sizeBins,bins);  
}

//__________________________________________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::SetHistoBinsSinglePt(const Int_t nBinsLimits, Double_t* binsLimits, Double_t* binsWidth)
{ 
  // set histo bins for inclusive pt spectra 
  // function takes array of limits and widths (e.g. 1 GeV bins up to 10 GeV, 2 GeV width up to 20 GeV, ...)  
  // array size of binsLimits: nBinsLimits 
  // array size of binsWidth: nBinsLimits-1 
  // binsLimits have to be in increasing order
  // if binning undefined for any slice, original binning will be used


  Double_t binLimitMin = binsLimits[0];
  Double_t binLimitMax = binsLimits[nBinsLimits-1];

  Double_t binLimit = binLimitMin; // start value 
  
  Int_t sizeUpperLim = 10000; //static_cast<Int_t>(binLimitMax/binsWidth[0])+1;
  TArrayD binsArray(sizeUpperLim);
  Int_t nBins = 0; 
  binsArray.SetAt(binLimitMin,nBins++);

  while(binLimit<binLimitMax && nBins<sizeUpperLim){

    Int_t currentSlice = -1;
    for(Int_t i=0; i<nBinsLimits; i++){
      if(binLimit >= 0.999*binsLimits[i]) currentSlice = i; // 0.999 numerical saftey factor 
    }
    
    Double_t currentBinWidth = binsWidth[currentSlice];
    binLimit += currentBinWidth;

    binsArray.SetAt(binLimit,nBins++);
  }
  
  Double_t* bins = binsArray.GetArray();
  
  SetHistoBinsSinglePt(nBins,bins); 
}

//____________________________________________________________________________________
void AliFragmentationFunctionCorrections::NormalizeTH1(TH1* hist, const Float_t nJets)
{
  // FF normalization: divide by bin width and normalize to entries/jets
  // should also work for TH2, but to be tested !!!

  if(nJets == 0){ // nothing to do
    if(fDebug>0)  Printf("%s:%d -- normalize: nJets = 0, do nothing", (char*)__FILE__,__LINE__);
    return; 
  }
  
  Int_t nBins = hist->GetNbinsX()*hist->GetNbinsY()*hist->GetNbinsZ();

  for(Int_t bin=0; bin<=nBins; bin++){ // include under-/overflow (?)

    Double_t binWidth = hist->GetBinWidth(bin);
    Double_t binCont  = hist->GetBinContent(bin);
    Double_t binErr   = hist->GetBinError(bin);
    
    binCont /= binWidth;
    binErr  /= binWidth;

    hist->SetBinContent(bin,binCont);
    hist->SetBinError(bin,binErr);
  }

  hist->Scale(1/nJets);
}

//_____________________________________________________
void AliFragmentationFunctionCorrections::NormalizeFF()
{
  // normalize FF

  if(fNCorrectionLevels>1){
    Printf("%s:%d -- FF normalization should be done BEFORE any correction !!!", (char*)__FILE__,__LINE__);
  }

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    if(fDebug>0) Printf(" normalizeFF: i %d, nJets %f",i,fNJets->At(i));

    NormalizeTH1(fCorrFF[0]->GetTrackPt(i),fNJets->At(i)); // always normalize corr level 0 = raw FF
    NormalizeTH1(fCorrFF[0]->GetZ(i),fNJets->At(i));
    NormalizeTH1(fCorrFF[0]->GetXi(i),fNJets->At(i));
  } 
}

//______________________________________________________
void AliFragmentationFunctionCorrections::NormalizeBgr()
{
  // normalize bgr/UE FF

  if(fNCorrectionLevelsBgr>1){
    Printf("%s:%d -- FF normalization should be done BEFORE any correction !!!", (char*)__FILE__,__LINE__);
  }

  for(Int_t i=0; i<fNJetPtSlices; i++){
    NormalizeTH1(fCorrBgr[0]->GetTrackPt(i), fNJetsBgr->At(i)); // always normalize corr level 0 = raw FF
    NormalizeTH1(fCorrBgr[0]->GetZ(i), fNJetsBgr->At(i));
    NormalizeTH1(fCorrBgr[0]->GetXi(i),fNJetsBgr->At(i));
  } 

}

//__________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawFF(TString strfile, TString strID, TString strFFID)
{ 
  // read raw FF - standard dir/list name
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  ReadRawFF(strfile,strdir,strlist,strFFID);
}

//____________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawFF(TString strfile, TString strdir, TString strlist, TString strFFID)
{
  // get raw FF from input file, project in jet pt slice
  // normalization done separately 

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read FF from file %s, dir %s ",(char*)__FILE__,__LINE__,strfile.Data(), strdir.Data());

  gDirectory->cd(strdir);

  TList* list = 0;
  
  if(strlist && strlist.Length()){
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }

  TString hnameJetPt(Form("fh1FFJetPt%s",strFFID.Data()));
  TString hnameTrackPt(Form("fh2FFTrackPt%s",strFFID.Data()));
  TString hnameZ(Form("fh2FFZ%s",strFFID.Data()));
  TString hnameXi(Form("fh2FFXi%s",strFFID.Data()));

  TH1F* fh1FFJetPt   = 0;
  TH2F* fh2FFTrackPt = 0;
  TH2F* fh2FFZ       = 0;
  TH2F* fh2FFXi      = 0;
  
  if(list){
    fh1FFJetPt   = (TH1F*) list->FindObject(hnameJetPt);
    fh2FFTrackPt = (TH2F*) list->FindObject(hnameTrackPt);
    fh2FFZ       = (TH2F*) list->FindObject(hnameZ);  
    fh2FFXi      = (TH2F*) list->FindObject(hnameXi); 
  }
  else{
    fh1FFJetPt   = (TH1F*) gDirectory->Get(hnameJetPt);
    fh2FFTrackPt = (TH2F*) gDirectory->Get(hnameTrackPt);
    fh2FFZ       = (TH2F*) gDirectory->Get(hnameZ);  
    fh2FFXi      = (TH2F*) gDirectory->Get(hnameXi);
  }


  if(!fh1FFJetPt)  { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameJetPt.Data());   return; }
  if(!fh2FFTrackPt){ Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameTrackPt.Data()); return; }
  if(!fh2FFZ)      { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameZ.Data());       return; }
  if(!fh2FFXi)     { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameXi.Data());      return; }

  fh1FFJetPt->SetDirectory(0);
  fh2FFTrackPt->SetDirectory(0);
  fh2FFZ->SetDirectory(0);  
  fh2FFXi->SetDirectory(0); 

  f.Close();  


  // nJets per bin

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);
    
    Int_t binLo = static_cast<Int_t>(fh1FFJetPt->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh1FFJetPt->FindBin(jetPtUpLim)) - 1;
				     
    Float_t nJetsBin = fh1FFJetPt->Integral(binLo,binUp);
    
    fNJets->SetAt(nJetsBin,i); 

    if(fDebug>0) Printf("jet pt %d to %d: nJets %f",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),fNJets->At(i));
  }
  
  // projections: FF 
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);

    Int_t binLo = static_cast<Int_t>(fh2FFTrackPt->GetXaxis()->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh2FFTrackPt->GetXaxis()->FindBin(jetPtUpLim))-1;

    if(binUp > fh2FFTrackPt->GetNbinsX()){
      Printf("%s:%d -- jet pt range %0.3f exceeds histo limits",(char*)__FILE__,__LINE__,jetPtUpLim); 
      return; 
    }
    
    TString strNameFFPt(Form("fh1FFTrackPt%s_%02d_%02d",strFFID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZ(Form("fh1FFZ%s_%02d_%02d",strFFID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXi(Form("fh1FFXi%s_%02d_%02d",strFFID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    
    // appendix 'unbinned' to avoid histos with same name after rebinning
    TH1F* projPt = (TH1F*) fh2FFTrackPt->ProjectionY(strNameFFPt+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    TH1F* projZ  = (TH1F*) fh2FFZ->ProjectionY(strNameFFZ+"_unBinned",binLo,binUp,"o");
    TH1F* projXi = (TH1F*) fh2FFXi->ProjectionY(strNameFFXi+"_unBinned",binLo,binUp,"o");
    
    if(fNHistoBinsPt[i]) projPt = (TH1F*) projPt->Rebin(fNHistoBinsPt[i],strNameFFPt,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  projZ  = (TH1F*) projZ->Rebin(fNHistoBinsZ[i],strNameFFZ,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) projXi = (TH1F*) projXi->Rebin(fNHistoBinsXi[i],strNameFFXi,fHistoBinsXi[i]->GetArray());

    projPt->SetNameTitle(strNameFFPt,"");
    projZ->SetNameTitle(strNameFFZ,"");
    projXi->SetNameTitle(strNameFFXi,"");

    // raw FF = corr level 0
    fCorrFF[0]->AddCorrHistos(i,projPt,projZ,projXi);
  }


  delete fh1FFJetPt;
  delete fh2FFTrackPt;
  delete fh2FFZ;  
  delete fh2FFXi; 
}

//_____________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawBgr(TString strfile, TString strID, TString strBgrID, TString strFFID)
{ 
  // read raw FF - standard dir/list name
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  ReadRawBgr(strfile,strdir,strlist,strBgrID,strFFID);
}

//_______________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawBgr(TString strfile, TString strdir, TString strlist, TString strBgrID, TString strFFID)
{
  // get raw FF from input file, project in jet pt slice
  // use jet dN/dpt corresponding to strFFID, bgr FF to strBgrID+strFFID
  // e.g. "fh1FFJetPtRecCuts", "fh2FFXiBgrPerpRecCuts"
  // normalization done separately 

  TString strID = strBgrID + strFFID;
  
  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read Bgr %s from file %s, dir %s ",(char*)__FILE__,__LINE__,strBgrID.Data(),strfile.Data(),strdir.Data());

  gDirectory->cd(strdir);

  TList* list = 0;
  
  if(strlist && strlist.Length()){
    if(!(list = (TList*) gDirectory->Get(strlist))){
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }

  TString hnameNJets = "fh1nRecJetsCuts"; 
  TString hnameJetPt(Form("fh1FFJetPt%s",strFFID.Data())); // not: strID.Data() !!! would not be proper normalization
  TString hnameBgrTrackPt(Form("fh2FFTrackPt%s",strID.Data()));
  TString hnameBgrZ(Form("fh2FFZ%s",strID.Data()));
  TString hnameBgrXi(Form("fh2FFXi%s",strID.Data()));

  TH1F* fh1NJets;
  TH1F* fh1FFJetPtBgr;  
  TH2F* fh2FFTrackPtBgr;
  TH2F* fh2FFZBgr;
  TH2F* fh2FFXiBgr;     
  
  if(list){
    fh1NJets        = (TH1F*) list->FindObject(hnameNJets); // needed for normalization of bgr out of 2 jets
    fh1FFJetPtBgr   = (TH1F*) list->FindObject(hnameJetPt);
    fh2FFTrackPtBgr = (TH2F*) list->FindObject(hnameBgrTrackPt);
    fh2FFZBgr       = (TH2F*) list->FindObject(hnameBgrZ);  
    fh2FFXiBgr      = (TH2F*) list->FindObject(hnameBgrXi); 
  }
  else{
    fh1NJets        = (TH1F*) gDirectory->Get(hnameNJets); // needed for normalization of bgr out of 2 jets
    fh1FFJetPtBgr   = (TH1F*) gDirectory->Get(hnameJetPt);
    fh2FFTrackPtBgr = (TH2F*) gDirectory->Get(hnameBgrTrackPt);
    fh2FFZBgr       = (TH2F*) gDirectory->Get(hnameBgrZ);  
    fh2FFXiBgr      = (TH2F*) gDirectory->Get(hnameBgrXi);  
  }


  if(!fh1FFJetPtBgr)  { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameJetPt.Data());      return; }
  if(!fh1NJets)       { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameNJets.Data());              }
  if(!fh2FFTrackPtBgr){ Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameBgrTrackPt.Data()); return; }
  if(!fh2FFZBgr)      { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameBgrZ.Data());       return; }
  if(!fh2FFXiBgr)     { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameBgrXi.Data());      return; }

  fh1FFJetPtBgr->SetDirectory(0);
  fh2FFTrackPtBgr->SetDirectory(0);
  fh2FFZBgr->SetDirectory(0);  
  fh2FFXiBgr->SetDirectory(0); 
  if(fh1NJets) fh1NJets->SetDirectory(0);

  f.Close();  

  // nJets per bin

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);
    
    Int_t binLo = static_cast<Int_t>(fh1FFJetPtBgr->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh1FFJetPtBgr->FindBin(jetPtUpLim)) - 1;
				     
    Float_t nJetsBin = fh1FFJetPtBgr->Integral(binLo,binUp);
    Double_t scaleF = 1;

    //if(strBgrID.Contains("Out2Jets")){  // scale by ratio 2 jets events / all events
    //  scaleF = fh1NJets->Integral(fh1NJets->FindBin(2),fh1NJets->GetNbinsX()) 
    //  / fh1NJets->Integral(fh1NJets->FindBin(1),fh1NJets->GetNbinsX());}


    if(strBgrID.Contains("OutAllJets") ){  // scale by ratio >3 jets events / all events

      if(fh1NJets){
	scaleF = fh1NJets->Integral(fh1NJets->FindBin(4),fh1NJets->GetNbinsX()) 
	  / fh1NJets->Integral(fh1NJets->FindBin(1),fh1NJets->GetNbinsX());
      }
      else{
	Printf("%s:%d -- use bgr OutAllJets, but histo fhn1NJets not found",(char*)__FILE__,__LINE__);    
      }
    }
    
 
    fNJetsBgr->SetAt(nJetsBin*scaleF,i); 

    if(fDebug>0) Printf("bgr jet pt %d to %d: nJets %f, scaleF %.2f",
			static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),nJetsBin,scaleF);

  }
  
  // projections: FF 
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);

    Int_t binLo = static_cast<Int_t>(fh2FFTrackPtBgr->GetXaxis()->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh2FFTrackPtBgr->GetXaxis()->FindBin(jetPtUpLim))-1;

    if(binUp > fh2FFTrackPtBgr->GetNbinsX()){
      Printf("%s:%d -- jet pt range %0.3f exceeds histo limits",(char*)__FILE__,__LINE__,jetPtUpLim); 
      return; 
    }
    
    TString strNameBgrPt(Form("fh1BgrTrackPt%s_%02d_%02d",strID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameBgrZ(Form("fh1BgrZ%s_%02d_%02d",strID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameBgrXi(Form("fh1BgrXi%s_%02d_%02d",strID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    
    // appendix 'unbinned' to avoid histos with same name after rebinning
    TH1F* projPt = (TH1F*) fh2FFTrackPtBgr->ProjectionY(strNameBgrPt+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    TH1F* projZ  = (TH1F*) fh2FFZBgr->ProjectionY(strNameBgrZ+"_unBinned",binLo,binUp,"o");
    TH1F* projXi = (TH1F*) fh2FFXiBgr->ProjectionY(strNameBgrXi+"_unBinned",binLo,binUp,"o");
    
    if(fNHistoBinsPt[i]) projPt = (TH1F*) projPt->Rebin(fNHistoBinsPt[i],strNameBgrPt,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  projZ  = (TH1F*) projZ->Rebin(fNHistoBinsZ[i],strNameBgrZ,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) projXi = (TH1F*) projXi->Rebin(fNHistoBinsXi[i],strNameBgrXi,fHistoBinsXi[i]->GetArray());

    projPt->SetNameTitle(strNameBgrPt,"");
    projZ->SetNameTitle(strNameBgrZ,"");
    projXi->SetNameTitle(strNameBgrXi,"");
    
    // raw bgr = corr level 0
    fCorrBgr[0]->AddCorrHistos(i,projPt,projZ,projXi);
  }

  delete fh1FFJetPtBgr;
  if(fh1NJets) delete fh1NJets;
  delete fh2FFTrackPtBgr;
  delete fh2FFZBgr;  
  delete fh2FFXiBgr; 
}

//_____________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawBgrEmbedding(TString strfile, TString strID, TString strFFID)
{ 
  // read raw FF - standard dir/list name
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  ReadRawBgrEmbedding(strfile,strdir,strlist,strFFID);
}

//_______________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawBgrEmbedding(TString strfile, TString strdir, TString strlist, TString strFFID)
{
  // get raw FF from input file, project in jet pt slice
  // for embedding, the bgr FF are taken from histos "fh1FFJetPtRecCuts", "fh2FFXiRecCuts"
  // normalization done separately 

  TString strBgrID = "BckgEmbed";
  TString strID = strBgrID + strFFID;

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read Bgr %s from file %s ",(char*)__FILE__,__LINE__,strFFID.Data(),strfile.Data());

  gDirectory->cd(strdir);

  TList* list = 0;
  
  if(!(list = (TList*) gDirectory->Get(strlist))){
    Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
    return;
  }

  TString hnameNJets = "fh1nRecJetsCuts"; 
  TString hnameJetPt(Form("fh1FFJetPt%s",strFFID.Data())); 
  TString hnameBgrTrackPt(Form("fh2FFTrackPt%s",strFFID.Data()));
  TString hnameBgrZ(Form("fh2FFZ%s",strFFID.Data()));
  TString hnameBgrXi(Form("fh2FFXi%s",strFFID.Data()));

  TH1F* fh1NJets        = (TH1F*) list->FindObject(hnameNJets); // needed for normalization of bgr out of 2 jets
  TH1F* fh1FFJetPtBgr   = (TH1F*) list->FindObject(hnameJetPt);
  TH2F* fh2FFTrackPtBgr = (TH2F*) list->FindObject(hnameBgrTrackPt);
  TH2F* fh2FFZBgr       = (TH2F*) list->FindObject(hnameBgrZ);  
  TH2F* fh2FFXiBgr      = (TH2F*) list->FindObject(hnameBgrXi); 

  if(!fh1FFJetPtBgr)  { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameJetPt.Data());      return; }
  if(!fh1NJets)       { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameNJets.Data());      return; }
  if(!fh2FFTrackPtBgr){ Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameBgrTrackPt.Data()); return; }
  if(!fh2FFZBgr)      { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameBgrZ.Data());       return; }
  if(!fh2FFXiBgr)     { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameBgrXi.Data());      return; }

  fh1FFJetPtBgr->SetDirectory(0);
  fh1NJets->SetDirectory(0);
  fh2FFTrackPtBgr->SetDirectory(0);
  fh2FFZBgr->SetDirectory(0);  
  fh2FFXiBgr->SetDirectory(0); 

  f.Close();  

  // nJets per bin

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);
    
    Int_t binLo = static_cast<Int_t>(fh1FFJetPtBgr->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh1FFJetPtBgr->FindBin(jetPtUpLim)) - 1;
				     
    Float_t nJetsBin = fh1FFJetPtBgr->Integral(binLo,binUp);
    Double_t scaleF = 1;

    fNJetsBgr->SetAt(nJetsBin*scaleF,i); 

    if(fDebug>0) Printf("bgr jet pt %d to %d: nJets %f, scaleF %.2f",
			static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),nJetsBin,scaleF);

  }
  
  // projections: FF 
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);

    Int_t binLo = static_cast<Int_t>(fh2FFTrackPtBgr->GetXaxis()->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh2FFTrackPtBgr->GetXaxis()->FindBin(jetPtUpLim))-1;

    if(binUp > fh2FFTrackPtBgr->GetNbinsX()){
      Printf("%s:%d -- jet pt range %0.3f exceeds histo limits",(char*)__FILE__,__LINE__,jetPtUpLim); 
      return; 
    }
    
    TString strNameBgrPt(Form("fh1BgrTrackPt%s_%02d_%02d",strID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameBgrZ(Form("fh1BgrZ%s_%02d_%02d",strID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameBgrXi(Form("fh1BgrXi%s_%02d_%02d",strID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    
    // appendix 'unbinned' to avoid histos with same name after rebinning
    TH1F* projPt = (TH1F*) fh2FFTrackPtBgr->ProjectionY(strNameBgrPt+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    TH1F* projZ  = (TH1F*) fh2FFZBgr->ProjectionY(strNameBgrZ+"_unBinned",binLo,binUp,"o");
    TH1F* projXi = (TH1F*) fh2FFXiBgr->ProjectionY(strNameBgrXi+"_unBinned",binLo,binUp,"o");
    
    if(fNHistoBinsPt[i]) projPt = (TH1F*) projPt->Rebin(fNHistoBinsPt[i],strNameBgrPt,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  projZ  = (TH1F*) projZ->Rebin(fNHistoBinsZ[i],strNameBgrZ,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) projXi = (TH1F*) projXi->Rebin(fNHistoBinsXi[i],strNameBgrXi,fHistoBinsXi[i]->GetArray());

    projPt->SetNameTitle(strNameBgrPt,"");
    projZ->SetNameTitle(strNameBgrZ,"");
    projXi->SetNameTitle(strNameBgrXi,"");
    
    // raw bgr = corr level 0
    fCorrBgr[0]->AddCorrHistos(i,projPt,projZ,projXi);
  }

  delete fh1FFJetPtBgr;
  delete fh1NJets;
  delete fh2FFTrackPtBgr;
  delete fh2FFZBgr;  
  delete fh2FFXiBgr; 
}


//__________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteOutput(TString strfile, TString strdir, Bool_t updateOutfile)
{
  // write histos to file
  // skip histos with 0 entries

  TString outfileOption = "RECREATE";
  if(updateOutfile) outfileOption = "UPDATE";

  TFile f(strfile,outfileOption); 

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- write FF to file %s ",(char*)__FILE__,__LINE__,strfile.Data());

  if(strdir && strdir.Length()){
    TDirectory* dir = f.mkdir(strdir);
    dir->cd(); 
  } 

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    for(Int_t c=0; c<fNCorrectionLevels; c++) if(fCorrFF[c]->GetTrackPt(i)->GetEntries()) fCorrFF[c]->GetTrackPt(i)->Write();
    for(Int_t c=0; c<fNCorrectionLevels; c++) if(fCorrFF[c]->GetZ(i)->GetEntries())       fCorrFF[c]->GetZ(i)->Write();
    for(Int_t c=0; c<fNCorrectionLevels; c++) if(fCorrFF[c]->GetXi(i)->GetEntries())      fCorrFF[c]->GetXi(i)->Write();

    if(fh1FFXiShift[i]->GetEntries()) fh1FFXiShift[i]->Write();

    for(Int_t c=0; c<fNCorrectionLevelsBgr; c++) if(fCorrBgr[c]->GetTrackPt(i)->GetEntries()) fCorrBgr[c]->GetTrackPt(i)->Write();
    for(Int_t c=0; c<fNCorrectionLevelsBgr; c++) if(fCorrBgr[c]->GetZ(i)->GetEntries())       fCorrBgr[c]->GetZ(i)->Write();
    for(Int_t c=0; c<fNCorrectionLevelsBgr; c++) if(fCorrBgr[c]->GetXi(i)->GetEntries())      fCorrBgr[c]->GetXi(i)->Write();

    
    if(fh1FFTrackPtBackFolded[i] && fh1FFTrackPtBackFolded[i]->GetEntries()) fh1FFTrackPtBackFolded[i]->Write();
    if(fh1FFZBackFolded[i]       && fh1FFZBackFolded[i]->GetEntries())       fh1FFZBackFolded[i]->Write();
    if(fh1FFXiBackFolded[i]      && fh1FFXiBackFolded[i]->GetEntries())      fh1FFXiBackFolded[i]->Write();


    if(fh1FFRatioTrackPtFolded[i] && fh1FFRatioTrackPtFolded[i]->GetEntries()) fh1FFRatioTrackPtFolded[i]->Write();
    if(fh1FFRatioZFolded[i]       && fh1FFRatioZFolded[i]->GetEntries())       fh1FFRatioZFolded[i]->Write();
    if(fh1FFRatioXiFolded[i]      && fh1FFRatioXiFolded[i]->GetEntries())      fh1FFRatioXiFolded[i]->Write();

    if(fh1FFRatioTrackPtBackFolded[i] && fh1FFRatioTrackPtBackFolded[i]->GetEntries()) fh1FFRatioTrackPtBackFolded[i]->Write();
    if(fh1FFRatioZBackFolded[i] && fh1FFRatioZBackFolded[i]->GetEntries())             fh1FFRatioZBackFolded[i]->Write();
    if(fh1FFRatioXiBackFolded[i] &&fh1FFRatioXiBackFolded[i]->GetEntries())            fh1FFRatioXiBackFolded[i]->Write();

    if(fh1FoldingCorrPt[i] && fh1FoldingCorrPt[i]->GetEntries()) fh1FoldingCorrPt[i]->Write(); // TEST!!!
    if(fh1FoldingCorrZ[i] && fh1FoldingCorrZ[i]->GetEntries())     fh1FoldingCorrZ[i]->Write();  // TEST!!!
    if(fh1FoldingCorrXi[i] && fh1FoldingCorrXi[i]->GetEntries()) fh1FoldingCorrXi[i]->Write(); // TEST!!!
  }
  
  // inclusive track pt

  for(Int_t c=0; c<fNCorrectionLevelsSinglePt; c++) if(fCorrSinglePt[c]->GetTrackPt(0)->GetEntries()) fCorrSinglePt[c]->GetTrackPt(0)->Write();
  if(fh1SingleTrackPtBackFolded)      fh1SingleTrackPtBackFolded->Write();  
  if(fh1RatioSingleTrackPtFolded)     fh1RatioSingleTrackPtFolded->Write();  
  if(fh1RatioSingleTrackPtBackFolded) fh1RatioSingleTrackPtBackFolded->Write();  

  


  f.Close();  
}

//____________________________________________________________________________________________________________________________________
THnSparse* AliFragmentationFunctionCorrections::TH1toSparse(const TH1F* hist, TString strName, TString strTit, const Bool_t fillConst)
{
  // copy 1-dimensional histo to THnSparse 
  // if fillConst TRUE, create THnSparse with same binning as hist but all entries = 1
  // histos with variable bin size are supported

  // note: function returns pointer to 'new' THnSparse on heap, object needs to be deleted by user

  THnSparseF* fhnHist;

  Int_t nBins       = hist->GetXaxis()->GetNbins();
  Int_t nBinsVec[1] = { nBins };

  const Double_t* binsVec = hist->GetXaxis()->GetXbins()->GetArray();
  
  if(binsVec){ // binsVec only neq NULL if histo was rebinned before

    fhnHist = new THnSparseF(strName,strTit,1,nBinsVec/*,binMinVec,binMaxVec*/);
    fhnHist->SetBinEdges(0,binsVec);
  }
  else{ // uniform bin size
    
    Double_t xMin  = hist->GetXaxis()->GetXmin();
    Double_t xMax  = hist->GetXaxis()->GetXmax();
    
    Double_t binMinVec[1] = { xMin };
    Double_t binMaxVec[1] = { xMax };
    
    fhnHist = new THnSparseF(strName,strTit,1,nBinsVec,binMinVec,binMaxVec);
  }
 

  for(Int_t bin=1; bin<=nBins; bin++){ 

    Double_t binCenter = fhnHist->GetAxis(0)->GetBinCenter(bin);

    Double_t binCoord[] = {binCenter}; 
    fhnHist->Fill(binCoord,1); // initially need to create the bin 

    Long64_t binIndex = fhnHist->GetBin(binCoord);

    Double_t cont =  hist->GetBinContent(bin);
    Double_t err  =  hist->GetBinError(bin);
    
    if(fillConst){
      cont = 1;
      err  = 0;
    }

    fhnHist->SetBinContent(binIndex,cont);
    fhnHist->SetBinError(binIndex,err);
 }

 return fhnHist;
}

//______________________________________________________________________________________________________________________________________________
TH1F* AliFragmentationFunctionCorrections::Unfold(THnSparse* hnHist, const THnSparse* hnResponse, const THnSparse* hnEff, const Int_t nIter,
						  const Bool_t useCorrelatedErrors, const THnSparse* hnPrior)
{
  // unfold input histo 

  AliCFUnfolding unfolding("unfolding","",1,hnResponse,hnEff,hnHist,hnPrior);  // arg3: nVar; hnEff required, hnPrior not (defaults to 0x0)
  unfolding.SetMaxNumberOfIterations(nIter);
  // unfolding.SetMaxChi2PerDOF(1.e-07); // OBSOLETE !!!
  // if(useSmoothing) unfolding.UseSmoothing(); 

  if(useCorrelatedErrors) unfolding.SetUseCorrelatedErrors();
  unfolding.Unfold();

  THnSparse* unfolded = unfolding.GetUnfolded();

  TString hnameUnf = hnHist->GetName();
  hnameUnf.Append("_unf");
  unfolded->SetNameTitle(hnameUnf,hnHist->GetTitle());
  
  TH1F* h1Unfolded = (TH1F*) unfolded->Projection(0); 
  h1Unfolded->SetNameTitle(hnHist->GetName(),hnHist->GetTitle());
  
  return h1Unfolded;
}

//___________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::UnfoldHistos(const Int_t nIter, const Bool_t useCorrelatedErrors, const Int_t type)
{
  // loop over jet pt slices and unfold dN/dpt spectra
  
  TString labelCorr = fCorrFF[fNCorrectionLevels-1]->GetLabel(); 
  if(!labelCorr.Contains("unfold")) AddCorrectionLevel("unfold");
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
 
    TH1F* hist = 0;
    if(type == kFlagPt)      hist = fCorrFF[fNCorrectionLevels-2]->GetTrackPt(i); // level -2: before unfolding, level -1: unfolded
    else if(type == kFlagZ)  hist = fCorrFF[fNCorrectionLevels-2]->GetZ(i);       // level -2: before unfolding, level -1: unfolded
    else if(type == kFlagXi) hist = fCorrFF[fNCorrectionLevels-2]->GetXi(i);      // level -2: before unfolding, level -1: unfolded
    else{ 
      Printf("%s%d unknown type",(char*)__FILE__,__LINE__);
      return;
    }
    
    if(!hist){
      Printf("%s%d no histo found ",(char*)__FILE__,__LINE__);
      return;
    }


    THnSparse* hnResponse = 0;
    if(type == kFlagPt)      hnResponse = fhnResponsePt[i];
    else if(type == kFlagZ)  hnResponse = fhnResponseZ[i];
    else if(type == kFlagXi) hnResponse = fhnResponseXi[i];

    TH1F* hPrior = 0;
    if(type == kFlagPt && fh1FFTrackPtPrior[i]  && ((TString(fh1FFTrackPtPrior[i]->GetName())).Length() > 0) ) hPrior = fh1FFTrackPtPrior[i];
    else if(type == kFlagZ  && fh1FFZPrior[i]   && ((TString(fh1FFZPrior[i]->GetName())).Length() > 0)       ) hPrior = fh1FFZPrior[i];
    else if(type == kFlagXi && fh1FFXiPrior[i]  && ((TString(fh1FFXiPrior[i]->GetName())).Length() > 0)      ) hPrior = fh1FFXiPrior[i];


    TString histNameTHn;
    histNameTHn = hist->GetName();
    histNameTHn.ReplaceAll("TH1","THn");
    

    TString priorNameTHn; 
    if(hPrior){
      priorNameTHn = hPrior->GetName();
      priorNameTHn.ReplaceAll("TH1","THn");
    }

    TString histNameBackFolded;
    histNameBackFolded = hist->GetName();
    histNameBackFolded.Append("_backfold");
    
    TString histNameRatioFolded;
    histNameRatioFolded = hist->GetName();
    histNameRatioFolded.ReplaceAll("fh1FF","hRatioFF");
    
    histNameRatioFolded.Append("_unfold");

    TString histNameRatioBackFolded;
    histNameRatioBackFolded = hist->GetName();
    histNameRatioBackFolded.ReplaceAll("fh1FF","hRatioFF");
    histNameRatioBackFolded.Append("_backfold");
 
    THnSparse* hnHist           = TH1toSparse(hist,histNameTHn,hist->GetTitle());
    THnSparse* hnFlatEfficiency = TH1toSparse(hist,"fhnEfficiency","eff",kTRUE); // could optionally also use real eff 
    THnSparse* hnPrior          = 0;
    if(hPrior) hnPrior = TH1toSparse(hPrior,priorNameTHn,hPrior->GetTitle());
    
    //THnSparse* hnUnfolded 
    //  = Unfold(hnHist,hnResponse,hnFlatEfficiency,nIter,useCorrelatedErrors,hnPrior);  
     
    //TH1F* hUnfolded = (TH1F*) hnUnfolded->Projection(0); 
    //hUnfolded->SetNameTitle(hist->GetName(),hist->GetTitle());
    
    TH1F* hUnfolded 
      = Unfold(hnHist,hnResponse,hnFlatEfficiency,nIter,useCorrelatedErrors,hnPrior);  
     
    //TH1F* hUnfolded = (TH1F*) hnUnfolded->Projection(0); 
    //hUnfolded->SetNameTitle(hist->GetName(),hist->GetTitle());

    // errors
    for(Int_t bin=1; bin<=hist->GetNbinsX(); bin++){
      
      Double_t contOrg = hist->GetBinContent(bin);
      Double_t errOrg  = hist->GetBinError(bin);
      
      Double_t contUnf = hUnfolded->GetBinContent(bin);
      //Double_t errUnf  = hUnfolded->GetBinError(bin);

      Double_t relErrOrg = 0;
      if(contOrg) relErrOrg = errOrg/contOrg;
      
      //       std::cout<<" cont original "<<contOrg<<" err original "<<errOrg<<std::endl;
      //       std::cout<<" cont unf "<<contUnf<<" errUnf"<<errUnf<<std::endl;
      //       std::cout<<" err/cont original "<<relErrOrg<<std::endl;
      //       if(contUnf)  std::cout<<" err/conf unf "<<errUnf/contUnf<<std::endl;
      
      hUnfolded->SetBinError(bin,relErrOrg*contUnf);
      
      //       if(hUnfolded->GetBinContent(bin)){
      // 	std::cout<<" modified err unfolded "<<hUnfolded->GetBinError(bin)/hUnfolded->GetBinContent(bin)<<std::endl; 
      //       }
    }
    

    if(type == kFlagPt) fCorrFF[fNCorrectionLevels-1]->AddCorrHistos(i,hUnfolded,0,0);
    if(type == kFlagZ)  fCorrFF[fNCorrectionLevels-1]->AddCorrHistos(i,0,hUnfolded,0);
    if(type == kFlagXi) fCorrFF[fNCorrectionLevels-1]->AddCorrHistos(i,0,0,hUnfolded);

    // backfolding: apply response matrix to unfolded spectrum
    TH1F* hBackFolded = ApplyResponse(hUnfolded,hnResponse); 
    hBackFolded->SetNameTitle(histNameBackFolded,hist->GetTitle());

    if(type == kFlagPt) fh1FFTrackPtBackFolded[i] = hBackFolded;
    if(type == kFlagZ)  fh1FFZBackFolded[i]       = hBackFolded;
    if(type == kFlagXi) fh1FFXiBackFolded[i]      = hBackFolded;
    
    // ratio unfolded to original histo 
    TH1F* hRatioUnfolded = (TH1F*) hUnfolded->Clone(histNameRatioFolded);
    hRatioUnfolded->Reset();
    hRatioUnfolded->Divide(hUnfolded,hist,1,1,"B");

    if(type == kFlagPt) fh1FFRatioTrackPtFolded[i] = hRatioUnfolded;
    if(type == kFlagZ)  fh1FFRatioZFolded[i]       = hRatioUnfolded;
    if(type == kFlagXi) fh1FFRatioXiFolded[i]      = hRatioUnfolded;


    // ratio backfolded to original histo
    TH1F* hRatioBackFolded = (TH1F*) hBackFolded->Clone(histNameRatioBackFolded);
    hRatioBackFolded->Reset();
    hRatioBackFolded->Divide(hBackFolded,hist,1,1,"B");

    if(type == kFlagPt) fh1FFRatioTrackPtBackFolded[i] = hRatioBackFolded;
    if(type == kFlagZ)  fh1FFRatioZBackFolded[i]       = hRatioBackFolded;
    if(type == kFlagXi) fh1FFRatioXiBackFolded[i]      = hRatioBackFolded;
    
    delete hnHist;
    delete hnFlatEfficiency;
 }
}

//_____________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::UnfoldPt(const Int_t nIter, const Bool_t useCorrelatedErrors)
{
  
  Int_t type = kFlagPt;
  UnfoldHistos(nIter, useCorrelatedErrors, type);
}

//_____________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::UnfoldZ(const Int_t nIter, const Bool_t useCorrelatedErrors)
{
  
  Int_t type = kFlagZ;
  UnfoldHistos(nIter, useCorrelatedErrors, type);
}

//_____________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::UnfoldXi(const Int_t nIter, const Bool_t useCorrelatedErrors)
{
  
  Int_t type = kFlagXi;
  UnfoldHistos(nIter, useCorrelatedErrors, type);
}


//______________________________________________________________________________________________
TH1F* AliFragmentationFunctionCorrections::ApplyResponse(const TH1F* hist, THnSparse* hnResponse)
{
  // apply (multiply) response matrix to hist 

  TH1F* hOut = new TH1F(*hist);
  hOut->Reset();

  const Int_t axisM = 0; 
  const Int_t axisT = 1;
 
  Int_t nBinsM = hnResponse->GetAxis(axisM)->GetNbins();
  Int_t nBinsT = hnResponse->GetAxis(axisT)->GetNbins();

  // response matrix normalization
  // do it in this function and not when reading response, since for 'proper' normalization errors are difficult to assign
  // so for unfolding proper we leave it to CORRFW ...

  Double_t normFacResponse[nBinsT+1];

  for(Int_t iT=1; iT<=nBinsT; iT++){

    Double_t sumResp = 0;
    
    for(Int_t iM=1; iM<=nBinsM; iM++){
      
      Double_t coordM = hnResponse->GetAxis(axisM)->GetBinCenter(iM);
      Double_t coordT = hnResponse->GetAxis(axisT)->GetBinCenter(iT);
      
      Double_t binCoord[] = {coordM,coordT};
      
      Long64_t binIndex = hnResponse->GetBin(binCoord);
      
      Double_t resp = hnResponse->GetBinContent(binIndex); 
      
      sumResp += resp;
    }
    
    normFacResponse[iT] = 0;
    if(sumResp) normFacResponse[iT] = 1/sumResp;
  }
  
  
  
  for(Int_t iM=1; iM<=nBinsM; iM++){
    
    Double_t contM   = 0;

    for(Int_t iT=1; iT<=nBinsT; iT++){

      Double_t contT = hist->GetBinContent(iT);
      
      Double_t coordM = hnResponse->GetAxis(axisM)->GetBinCenter(iM);
      Double_t coordT = hnResponse->GetAxis(axisT)->GetBinCenter(iT);

      Double_t binCoord[] = {coordM,coordT};
      
      Long64_t binIndex = hnResponse->GetBin(binCoord);
      
      Double_t resp = hnResponse->GetBinContent(binIndex); 
      
      contM   += resp*normFacResponse[iT]*contT; 
    }

    hOut->SetBinContent(iM,contM);
  }

  return hOut;
}

//_______________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadEfficiency(TString strfile, TString strdir, TString strlist)
{
  // read reconstruction efficiency from file
  // argument strlist optional - read from directory strdir if not specified

  // temporary histos to hold histos from file
  TH1F* hEffPt[fNJetPtSlices]; 
  TH1F* hEffZ[fNJetPtSlices];
  TH1F* hEffXi[fNJetPtSlices];
  
  for(Int_t i=0; i<fNJetPtSlices; i++) hEffPt[i] = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hEffZ[i]  = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hEffXi[i] = 0;

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read efficiencies from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
   
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Int_t jetPtLoLim = static_cast<Int_t> (fJetPtSlices->At(i));
    Int_t jetPtUpLim = static_cast<Int_t> (fJetPtSlices->At(i+1));
    
    TString strNameEffPt(Form("hEffPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameEffZ(Form("hEffZ_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameEffXi(Form("hEffXi_%02d_%02d",jetPtLoLim,jetPtUpLim));
     
   
    if(list){
      hEffPt[i] = (TH1F*) list->FindObject(strNameEffPt); 
      hEffZ[i]  = (TH1F*) list->FindObject(strNameEffZ); 
      hEffXi[i] = (TH1F*) list->FindObject(strNameEffXi); 
    }
    else{
      hEffPt[i] = (TH1F*) gDirectory->Get(strNameEffPt); 
      hEffZ[i]  = (TH1F*) gDirectory->Get(strNameEffZ); 
      hEffXi[i] = (TH1F*) gDirectory->Get(strNameEffXi); 
    }
    
    if(!hEffPt[i]){
      Printf("%s:%d -- error retrieving efficiency %s", (char*)__FILE__,__LINE__,strNameEffPt.Data());
    }
  
    if(!hEffZ[i]){
      Printf("%s:%d -- error retrieving efficiency %s", (char*)__FILE__,__LINE__,strNameEffZ.Data());
    }    

    if(!hEffXi[i]){
      Printf("%s:%d -- error retrieving efficiency %s", (char*)__FILE__,__LINE__,strNameEffXi.Data());
    }


    if(fNHistoBinsPt[i]) hEffPt[i] = (TH1F*) hEffPt[i]->Rebin(fNHistoBinsPt[i],strNameEffPt+"_rebin",fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hEffZ[i]  = (TH1F*) hEffZ[i]->Rebin(fNHistoBinsZ[i],strNameEffZ+"_rebin",fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hEffXi[i] = (TH1F*) hEffXi[i]->Rebin(fNHistoBinsXi[i],strNameEffXi+"_rebin",fHistoBinsXi[i]->GetArray());

    if(hEffPt[i]) hEffPt[i]->SetDirectory(0); 
    if(hEffZ[i])  hEffZ[i]->SetDirectory(0); 
    if(hEffXi[i]) hEffXi[i]->SetDirectory(0); 

  } // jet slices loop

  f.Close();

  for(Int_t i=0; i<fNJetPtSlices; i++){ // 2nd loop: need to close input file before placing histos
    if(hEffPt[i]) new(fh1EffPt[i]) TH1F(*hEffPt[i]);
    if(hEffZ[i])  new(fh1EffZ[i])  TH1F(*hEffZ[i]);
    if(hEffXi[i]) new(fh1EffXi[i]) TH1F(*hEffXi[i]);
  }
}

//___________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadBgrEfficiency(TString strfile, TString strdir, TString strlist)
{
  // read bgr eff from file
  // argument strlist optional - read from directory strdir if not specified
 
  TH1F* hEffPtBgr[fNJetPtSlices]; 
  TH1F* hEffZBgr [fNJetPtSlices];
  TH1F* hEffXiBgr[fNJetPtSlices];

  for(Int_t i=0; i<fNJetPtSlices; i++) hEffPtBgr[i] = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hEffZBgr[i]  = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hEffXiBgr[i] = 0;


  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read bgr efficiencies from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
       
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Int_t jetPtLoLim = static_cast<Int_t> (fJetPtSlices->At(i));
    Int_t jetPtUpLim = static_cast<Int_t> (fJetPtSlices->At(i+1));
    
    TString strNameEffPtBgr(Form("hEffPtBgr%02dto%02d",jetPtLoLim,jetPtUpLim));
    TString strNameEffZBgr(Form("hEffZBgr%02dto%02d",jetPtLoLim,jetPtUpLim));
    TString strNameEffXiBgr(Form("hEffXiBgr%02dto%02d",jetPtLoLim,jetPtUpLim));
  
   
    if(list){
      hEffPtBgr[i] = (TH1F*) list->FindObject(strNameEffPtBgr); 
      hEffZBgr[i]  = (TH1F*) list->FindObject(strNameEffZBgr); 
      hEffXiBgr[i] = (TH1F*) list->FindObject(strNameEffXiBgr); 
    }
    else{
      hEffPtBgr[i] = (TH1F*) gDirectory->Get(strNameEffPtBgr); 
      hEffZBgr[i]  = (TH1F*) gDirectory->Get(strNameEffZBgr); 
      hEffXiBgr[i] = (TH1F*) gDirectory->Get(strNameEffXiBgr); 
    }
    
    if(!hEffPtBgr[i]){
      Printf("%s:%d -- error retrieving efficiency %s", (char*)__FILE__,__LINE__,strNameEffPtBgr.Data());
    }
  
    if(!hEffZBgr[i]){
      Printf("%s:%d -- error retrieving efficiency %s", (char*)__FILE__,__LINE__,strNameEffZBgr.Data());
    }    

    if(!hEffXiBgr[i]){
      Printf("%s:%d -- error retrieving efficiency %s", (char*)__FILE__,__LINE__,strNameEffXiBgr.Data());
    }


    if(fNHistoBinsPt[i]) hEffPtBgr[i] = (TH1F*) hEffPtBgr[i]->Rebin(fNHistoBinsPt[i],strNameEffPtBgr+"_rebin",fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hEffZBgr[i]  = (TH1F*) hEffZBgr[i]->Rebin(fNHistoBinsZ[i],strNameEffZBgr+"_rebin",fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hEffXiBgr[i] = (TH1F*) hEffXiBgr[i]->Rebin(fNHistoBinsXi[i],strNameEffXiBgr+"_rebin",fHistoBinsXi[i]->GetArray());

    if(hEffPtBgr[i]) hEffPtBgr[i]->SetDirectory(0); 
    if(hEffZBgr[i])  hEffZBgr[i]->SetDirectory(0); 
    if(hEffXiBgr[i]) hEffXiBgr[i]->SetDirectory(0); 

  } // jet slices loop

  f.Close();

  for(Int_t i=0; i<fNJetPtSlices; i++){ // 2nd loop: need to close input file before placing histos
    if(hEffPtBgr[i]) new(fh1EffBgrPt[i]) TH1F(*hEffPtBgr[i]);
    if(hEffZBgr[i])  new(fh1EffBgrZ[i])  TH1F(*hEffZBgr[i]);
    if(hEffXiBgr[i]) new(fh1EffBgrXi[i]) TH1F(*hEffXiBgr[i]);
  }
}

// ________________________________________________
void AliFragmentationFunctionCorrections::EffCorr()
{
  // apply efficiency correction

  AddCorrectionLevel("eff");

  for(Int_t i=0; i<fNJetPtSlices; i++){

    TH1F* histPt = fCorrFF[fNCorrectionLevels-2]->GetTrackPt(i);
    TH1F* histZ  = fCorrFF[fNCorrectionLevels-2]->GetZ(i);
    TH1F* histXi = fCorrFF[fNCorrectionLevels-2]->GetXi(i);

    TString histNamePt = histPt->GetName();
    TString histNameZ  = histZ->GetName();
    TString histNameXi = histXi->GetName();

    
    TH1F* hFFTrackPtEffCorr = (TH1F*) histPt->Clone(histNamePt);
    hFFTrackPtEffCorr->Divide(histPt,fh1EffPt[i],1,1,"");
    
    TH1F* hFFZEffCorr = (TH1F*) histZ->Clone(histNameZ);
    hFFZEffCorr->Divide(histZ,fh1EffZ[i],1,1,"");
    
    TH1F* hFFXiEffCorr = (TH1F*) histXi->Clone(histNameXi);
    hFFXiEffCorr->Divide(histXi,fh1EffXi[i],1,1,"");
    
    fCorrFF[fNCorrectionLevels-1]->AddCorrHistos(i,hFFTrackPtEffCorr,hFFZEffCorr,hFFXiEffCorr);
  }
}

//___________________________________________________
void AliFragmentationFunctionCorrections::EffCorrBgr()
{
  // apply efficiency correction to bgr distributions

  AddCorrectionLevelBgr("eff");

  Printf("%s:%d -- apply efficiency correction, corrLevel %d",(char*)__FILE__,__LINE__,fNCorrectionLevels-1);

  for(Int_t i=0; i<fNJetPtSlices; i++){

    TH1F* histPt = fCorrBgr[fNCorrectionLevelsBgr-2]->GetTrackPt(i);
    TH1F* histZ  = fCorrBgr[fNCorrectionLevelsBgr-2]->GetZ(i);
    TH1F* histXi = fCorrBgr[fNCorrectionLevelsBgr-2]->GetXi(i);
    
    TString histNamePt = histPt->GetName();
    TString histNameZ  = histZ->GetName();
    TString histNameXi = histXi->GetName();

    
    TH1F* hFFTrackPtEffCorr = (TH1F*) histPt->Clone(histNamePt);
    hFFTrackPtEffCorr->Divide(histPt,fh1EffPt[i],1,1,"");
    
    TH1F* hFFZEffCorr = (TH1F*) histZ->Clone(histNameZ);
    hFFZEffCorr->Divide(histZ,fh1EffZ[i],1,1,"");
    
    TH1F* hFFXiEffCorr = (TH1F*) histXi->Clone(histNameXi);
    hFFXiEffCorr->Divide(histXi,fh1EffXi[i],1,1,"");
    
    fCorrBgr[fNCorrectionLevelsBgr-1]->AddCorrHistos(i,hFFTrackPtEffCorr,hFFZEffCorr,hFFXiEffCorr);
  }
}

//______________________________________________________________________
void AliFragmentationFunctionCorrections::XiShift(const Int_t corrLevel)
{
  // re-evaluate jet energy after FF corrections from dN/dpt distribution
  // apply correction (shift) to dN/dxi distribution: xi = ln (pt/E) -> xi' = ln (pt/E') = ln (pt/E x E/E') = xi + ln E/E'
  // argument corrlevel: which level of correction to be corrected/shifted to 

  if(corrLevel>=fNCorrectionLevels){ 
    Printf(" calc xi shift: corrLevel exceeded - do nothing");
    return;
  }

  Double_t* jetPtUncorr = new Double_t[fNJetPtSlices];
  Double_t* jetPtCorr   = new Double_t[fNJetPtSlices];
  Double_t* deltaXi     = new Double_t[fNJetPtSlices];

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    TH1F* histPtRaw = fCorrFF[0]->GetTrackPt(i);
    TH1F* histPt    = fCorrFF[corrLevel]->GetTrackPt(i);

    Double_t ptUncorr = 0;
    Double_t ptCorr   = 0;

    for(Int_t bin = 1; bin<=histPtRaw->GetNbinsX(); bin++){

      Double_t cont   = histPtRaw->GetBinContent(bin);
      Double_t width  = histPtRaw->GetBinWidth(bin);
      Double_t meanPt = histPtRaw->GetBinCenter(bin);

      ptUncorr += meanPt*cont*width;
    }
    
    for(Int_t bin = 1; bin<=histPt->GetNbinsX(); bin++){
      
      Double_t cont   = histPt->GetBinContent(bin);
      Double_t width  = histPt->GetBinWidth(bin);
      Double_t meanPt = histPt->GetBinCenter(bin);
      
      ptCorr += meanPt*cont*width;
    }

    jetPtUncorr[i] = ptUncorr; 
    jetPtCorr[i]   = ptCorr;   
  }

  // ---

  for(Int_t i=0; i<fNJetPtSlices; i++){

    std::cout<<" jet pt bin "<<i<<" from "<<fJetPtSlices->At(i)<<" to "<<fJetPtSlices->At(i+1)
	     <<" jetPtUncorr "<<jetPtUncorr[i]<<" corrLevel "<<corrLevel<<" jet pt corr "<<jetPtCorr[i]<<" ratio "
	     <<jetPtCorr[i]/jetPtUncorr[i]<<std::endl;
  } 

  return; // for TEST!!!

  // calc dXi from dN/dpt distribution : 
  // sum over track pt for raw and corrected FF is equivalent to raw/corrected jet pt 

  for(Int_t i=0; i<fNJetPtSlices; i++){

    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);

    Double_t meanJetPt = 0.5*(jetPtUpLim+jetPtLoLim);
    
    Double_t ptUncorr = jetPtUncorr[i]; 
    Double_t ptCorr   = jetPtCorr[i]; 

    Double_t dXi = TMath::Log(ptCorr/ptUncorr);
    
    Printf(" calc xi shift: jet pt slice %d, mean jet pt %f, ptUncorr %f, ptCorr %f, ratio corr/uncorr %f, dXi %f "
	   ,i,meanJetPt,ptUncorr,ptCorr,ptCorr/ptUncorr,dXi);
    
    deltaXi[i] = dXi;
  }
  
  // book & fill new dN/dxi histos

  for(Int_t i=0; i<fNJetPtSlices; i++){

    TH1F* histXi = fCorrFF[corrLevel]->GetXi(i);
            
    Double_t dXi = deltaXi[i];

    Int_t nBins  = histXi->GetNbinsX();
    const Double_t* binsVec = histXi->GetXaxis()->GetXbins()->GetArray();
    Float_t binsVecNew[nBins+1];
    
    TString strName = histXi->GetName(); 
    strName.Append("_shift");
    TString strTit  = histXi->GetTitle();  

    TH1F* hXiCorr; 

    // create shifted histo ...

    if(binsVec){ // binsVec only neq NULL if histo was rebinned before

      for(Int_t bin=0; bin<nBins+1; bin++) binsVecNew[bin] = binsVec[bin] + dXi;    
      hXiCorr = new TH1F(strName,strTit,nBins,binsVecNew);
    }
    else{ // uniform bin size
      
      Double_t xMin  = histXi->GetXaxis()->GetXmin();
      Double_t xMax  = histXi->GetXaxis()->GetXmax();

      xMin += dXi;
      xMax += dXi;
      
      hXiCorr = new TH1F(strName,strTit,nBins,xMin,xMax);
    }

    // ... and fill

    for(Int_t bin=1; bin<nBins+1; bin++){
      Double_t cont = histXi->GetBinContent(bin);
      Double_t err  = histXi->GetBinError(bin);
      
      hXiCorr->SetBinContent(bin,cont);
      hXiCorr->SetBinError(bin,err);
    }
    
    new(fh1FFXiShift[i]) TH1F(*hXiCorr);
    delete hXiCorr;
  }

  delete[] jetPtUncorr;
  delete[] jetPtCorr;
  delete[] deltaXi;
}

//_____________________________________________________
void AliFragmentationFunctionCorrections::SubtractBgr(Double_t sysErr)
{
  // subtract bgr distribution from FF
  // the current corr level is used for both 
  
  AddCorrectionLevel("bgrSub");

  fstream ascii_out_pt;
  fstream ascii_out_z;
  fstream ascii_out_xi;

  if(sysErr) ascii_out_pt.open("sysErrUE_pt.txt",std::ios::out);
  if(sysErr) ascii_out_z.open("sysErrUE_z.txt",std::ios::out);
  if(sysErr) ascii_out_xi.open("sysErrUE_xi.txt",std::ios::out);

  for(Int_t i=0; i<fNJetPtSlices; i++){ // jet slices

    TH1F* histPt = fCorrFF[fNCorrectionLevels-2]->GetTrackPt(i);
    TH1F* histZ  = fCorrFF[fNCorrectionLevels-2]->GetZ(i);
    TH1F* histXi = fCorrFF[fNCorrectionLevels-2]->GetXi(i);
    
    TH1F* histPtBgr = fCorrBgr[fNCorrectionLevelsBgr-1]->GetTrackPt(i);
    TH1F* histZBgr  = fCorrBgr[fNCorrectionLevelsBgr-1]->GetZ(i);
    TH1F* histXiBgr = fCorrBgr[fNCorrectionLevelsBgr-1]->GetXi(i);

    TString histNamePt = histPt->GetName();
    TString histNameZ  = histZ->GetName();
    TString histNameXi = histXi->GetName();
    

    TH1F* hFFTrackPtBgrSub = (TH1F*) histPt->Clone(histNamePt);
    hFFTrackPtBgrSub->Add(histPtBgr,-1);
    
    TH1F* hFFZBgrSub =  (TH1F*) histZ->Clone(histNameZ);
    hFFZBgrSub->Add(histZBgr,-1);
    
    TH1F* hFFXiBgrSub = (TH1F*) histXi->Clone(histNameXi);
    hFFXiBgrSub->Add(histXiBgr,-1);

    fCorrFF[fNCorrectionLevels-1]->AddCorrHistos(i,hFFTrackPtBgrSub,hFFZBgrSub,hFFXiBgrSub);

    if(sysErr){

      for(Int_t bin=1; bin<=histPt->GetNbinsX(); bin++){

	//Double_t sigPlBgr = histPt->GetBinContent(bin);
	Double_t sig      = hFFTrackPtBgrSub->GetBinContent(bin);
	Double_t bgr      = histPtBgr->GetBinContent(bin);
	
	Double_t relErr = 0; 
	if(sig>0)  relErr = sysErr*bgr/sig;

// 	std::cout<<" pt bin "<<bin<<" mean "<<histPt->GetBinCenter(bin)
// 		 <<" sigPlBgr "<<sigPlBgr<<" sig "<<sig<<" bgr "<<bgr<<" relErr "<<relErr<<std::endl;
	
	ascii_out_pt<<i<<" "<<bin<<" "<<relErr<<std::endl;
      }


      for(Int_t bin=1; bin<=histZ->GetNbinsX(); bin++){
	
	//Double_t sigPlBgr = histZ->GetBinContent(bin);
	Double_t sig      = hFFZBgrSub->GetBinContent(bin);
	Double_t bgr      = histZBgr->GetBinContent(bin);
	
	Double_t relErr = 0; 
	if(sig>0)  relErr = sysErr*bgr/sig;
	
 	std::cout<<" z bin "<<bin<<" mean "<<histZ->GetBinCenter(bin)
 		 <<" sigPlBgr "<<histZ->GetBinContent(bin)<<" sig "<<sig<<" bgr "<<bgr<<" relErr "<<relErr<<std::endl;
	
	ascii_out_z<<i<<" "<<bin<<" "<<relErr<<std::endl;
      }


      for(Int_t bin=1; bin<=histXi->GetNbinsX(); bin++){
	
	//Double_t sigPlBgr = histXi->GetBinContent(bin);
	Double_t sig      = hFFXiBgrSub->GetBinContent(bin);
	Double_t bgr      = histXiBgr->GetBinContent(bin);
	
	Double_t relErr = 0; 
	if(sig>0)  relErr = sysErr*bgr/sig;
	
 	std::cout<<" xi bin "<<bin<<" mean "<<histXi->GetBinCenter(bin)
 		 <<" sigPlBgr "<<histXi->GetBinContent(bin)<<" sig "<<sig<<" bgr "<<bgr<<" relErr "<<relErr<<std::endl;
	
	ascii_out_xi<<i<<" "<<bin<<" "<<relErr<<std::endl;
      }
    }
  }

  if(sysErr) ascii_out_pt.close(); 
  if(sysErr) ascii_out_xi.close(); 


}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteSingleTrackEff(TString strInfile, TString strID, TString strOutfile, 
							      Bool_t updateOutfile, TString strOutDir,TString strPostfix)
{ 
  // read task ouput from MC and write single track eff - standard dir/list 
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  WriteSingleTrackEff(strInfile,strdir,strlist,strOutfile,updateOutfile,strOutDir,strPostfix);
}

//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteSingleTrackEff(TString strInfile, TString strdir, TString strlist, 
							      TString strOutfile, Bool_t updateOutfile, TString strOutDir,TString strPostfix)
{
  // read task output from MC and write single track eff as function of pt, eta, phi
  // argument strlist optional - read from directory strdir if not specified
  // write eff to file stroutfile - by default only 'update' file (don't overwrite)


  TH1D* hdNdptTracksMCPrimGen;
  TH2D* hdNdetadphiTracksMCPrimGen;
  
  TH1D* hdNdptTracksMCPrimRec;
  TH2D* hdNdetadphiTracksMCPrimRec;
    

  TFile f(strInfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- writeSingleTrackEff: open task ouput file %s ",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdir && strdir.Length()){
    if(fDebug>0) Printf("%s:%d -- writeSingleTrackEff: change dir to %s",(char*)__FILE__,__LINE__,strdir.Data());
    gDirectory->cd(strdir);
  }

  TList* list = 0;

  if(strlist && strlist.Length()){

    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }


  TString hnamePtRecEffGen = "fh1TrackQAPtRecEffGen";
  if(strPostfix.Length()) hnamePtRecEffGen.Form("fh1TrackQAPtRecEffGen_%s",strPostfix.Data()); 

  TString hnamePtRecEffRec = "fh1TrackQAPtRecEffRec";
  if(strPostfix.Length()) hnamePtRecEffRec.Form("fh1TrackQAPtRecEffRec_%s",strPostfix.Data());

  TString hnameEtaPhiRecEffGen = "fh2TrackQAEtaPhiRecEffGen";
  if(strPostfix.Length()) hnameEtaPhiRecEffGen.Form("fh2TrackQAEtaPhiRecEffGen_%s",strPostfix.Data());
  
  TString hnameEtaPhiRecEffRec = "fh2TrackQAEtaPhiRecEffRec";
  if(strPostfix.Length()) hnameEtaPhiRecEffRec.Form("fh2TrackQAEtaPhiRecEffRec_%s",strPostfix.Data());
  

  if(list){
    hdNdptTracksMCPrimGen       = (TH1D*) list->FindObject(hnamePtRecEffGen);
    hdNdetadphiTracksMCPrimGen  = (TH2D*) list->FindObject(hnameEtaPhiRecEffGen);
    
    hdNdptTracksMCPrimRec       = (TH1D*) list->FindObject(hnamePtRecEffRec);
    hdNdetadphiTracksMCPrimRec  = (TH2D*) list->FindObject(hnameEtaPhiRecEffRec);
  }
  else{
    hdNdptTracksMCPrimGen       = (TH1D*) gDirectory->Get(hnamePtRecEffGen);
    hdNdetadphiTracksMCPrimGen  = (TH2D*) gDirectory->Get(hnameEtaPhiRecEffGen);

    hdNdptTracksMCPrimRec       = (TH1D*) gDirectory->Get(hnamePtRecEffRec);
    hdNdetadphiTracksMCPrimRec  = (TH2D*) gDirectory->Get(hnameEtaPhiRecEffRec);
  }

  hdNdptTracksMCPrimGen->SetDirectory(0);
  hdNdetadphiTracksMCPrimGen->SetDirectory(0);
  hdNdptTracksMCPrimRec->SetDirectory(0);
  hdNdetadphiTracksMCPrimRec->SetDirectory(0);
  
  f.Close();

  // projections: dN/deta, dN/dphi 

  TH1D* hdNdetaTracksMCPrimGen = (TH1D*) hdNdetadphiTracksMCPrimGen->ProjectionX("hdNdetaTracksMcPrimGen");
  TH1D* hdNdphiTracksMCPrimGen = (TH1D*) hdNdetadphiTracksMCPrimGen->ProjectionY("hdNdphiTracksMcPrimGen");
 
  TH1D* hdNdetaTracksMCPrimRec = (TH1D*) hdNdetadphiTracksMCPrimRec->ProjectionX("hdNdetaTracksMcPrimRec");
  TH1D* hdNdphiTracksMCPrimRec = (TH1D*) hdNdetadphiTracksMCPrimRec->ProjectionY("hdNdphiTracksMcPrimRec");

  // rebin

  TString strNamePtGen = "hTrackPtGenPrim";
  TString strNamePtRec = "hTrackPtRecPrim";

  if(fNHistoBinsSinglePt) hdNdptTracksMCPrimGen = (TH1D*) hdNdptTracksMCPrimGen->Rebin(fNHistoBinsSinglePt,strNamePtGen,fHistoBinsSinglePt->GetArray());
  if(fNHistoBinsSinglePt) hdNdptTracksMCPrimRec = (TH1D*) hdNdptTracksMCPrimRec->Rebin(fNHistoBinsSinglePt,strNamePtRec,fHistoBinsSinglePt->GetArray());
 
    // efficiency: divide 

  TString hNameTrackEffPt = "hSingleTrackEffPt";
  if(strPostfix.Length()) hNameTrackEffPt.Form("hSingleTrackEffPt_%s",strPostfix.Data());
					       
  TString hNameTrackEffEta = "hSingleTrackEffEta";
  if(strPostfix.Length()) hNameTrackEffEta.Form("hSingleTrackEffEta_%s",strPostfix.Data());

  TString hNameTrackEffPhi = "hSingleTrackEffPhi";
  if(strPostfix.Length()) hNameTrackEffPhi.Form("hSingleTrackEffPhi_%s",strPostfix.Data());


  TH1F* hSingleTrackEffPt = (TH1F*) hdNdptTracksMCPrimRec->Clone(hNameTrackEffPt);
  hSingleTrackEffPt->Divide(hdNdptTracksMCPrimRec,hdNdptTracksMCPrimGen,1,1,"B"); // binominal errors

  TH1F* hSingleTrackEffEta = (TH1F*) hdNdetaTracksMCPrimRec->Clone(hNameTrackEffEta);
  hSingleTrackEffEta->Divide(hdNdetaTracksMCPrimRec,hdNdetaTracksMCPrimGen,1,1,"B"); // binominal errors
  
  TH1F* hSingleTrackEffPhi = (TH1F*) hdNdphiTracksMCPrimRec->Clone(hNameTrackEffPhi);
  hSingleTrackEffPhi->Divide(hdNdphiTracksMCPrimRec,hdNdphiTracksMCPrimGen,1,1,"B"); // binominal errors
  
  
  TString outfileOption = "RECREATE";
  if(updateOutfile)  outfileOption = "UPDATE";

  TFile out(strOutfile,outfileOption);

  if(!out.IsOpen()){
    Printf("%s:%d -- error opening efficiency output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- write efficiency to file %s ",(char*)__FILE__,__LINE__,strOutfile.Data());

  if(strOutDir && strOutDir.Length()){
    
    TDirectory* dir;
    if((dir = ((TDirectory*) gDirectory->Get(strOutDir)))) dir->cd(); 
    else{
      dir = out.mkdir(strOutDir);
      dir->cd(); 
    } 
  }

  hSingleTrackEffPt->Write();
  hSingleTrackEffEta->Write();
  hSingleTrackEffPhi->Write();
  
  out.Close();

  delete hdNdptTracksMCPrimGen;
  delete hdNdetadphiTracksMCPrimGen;
  delete hdNdptTracksMCPrimRec;
  delete hdNdetadphiTracksMCPrimRec;
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteSingleTrackSecCorr(TString strInfile, TString strID, TString strOutfile, 
								  Bool_t updateOutfile, TString strOutDir)
{ 
  // read task ouput from MC and write single track eff - standard dir/list 
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  WriteSingleTrackSecCorr(strInfile,strdir,strlist,strOutfile,updateOutfile,strOutDir);
}

//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteSingleTrackSecCorr(TString strInfile, TString strdir, TString strlist, 
								  TString strOutfile, Bool_t updateOutfile, TString strOutDir)
{
  // read task output from MC and write single track secondaries contamination correction as function of pt, eta, phi
  // argument strlist optional - read from directory strdir if not specified
  // write corr factor to file stroutfile - by default only 'update' file (don't overwrite)

  TH1D* hdNdptTracksMCPrimRec;
  TH2D* hdNdetadphiTracksMCPrimRec;
  
  TH1D* hdNdptTracksMCSecRecNS;
  TH2D* hdNdetadphiTracksMCSecRecNS;

  TH1D* hdNdptTracksMCSecRecSsc;
  TH2D* hdNdetadphiTracksMCSecRecSsc;

  TFile f(strInfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- writeSingleTrackEff: open task ouput file %s",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){

    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }


  if(list){
    hdNdptTracksMCPrimRec        = (TH1D*) list->FindObject("fh1TrackQAPtRecEffGen");
    hdNdetadphiTracksMCPrimRec   = (TH2D*) list->FindObject("fh2TrackQAEtaPhiRecEffGen");
    
    hdNdptTracksMCSecRecNS       = (TH1D*) list->FindObject("fh1TrackQAPtSecRecNS");
    hdNdetadphiTracksMCSecRecNS  = (TH2D*) list->FindObject("fh2TrackQAEtaPhiSecRecNS");

    hdNdptTracksMCSecRecSsc      = (TH1D*) list->FindObject("fh1TrackQAPtSecRecSsc");
    hdNdetadphiTracksMCSecRecSsc = (TH2D*) list->FindObject("fh2TrackQAEtaPhiSecRecSsc");

  }
  else{
    hdNdptTracksMCPrimRec        = (TH1D*) gDirectory->Get("fh1TrackQAPtRecEffGen");
    hdNdetadphiTracksMCPrimRec   = (TH2D*) gDirectory->Get("fh2TrackQAEtaPhiRecEffGen");

    hdNdptTracksMCSecRecNS       = (TH1D*) gDirectory->Get("fh1TrackQAPtSecRecNS");
    hdNdetadphiTracksMCSecRecNS  = (TH2D*) gDirectory->Get("fh2TrackQAEtaPhiSecRecNS");

    hdNdptTracksMCSecRecSsc      = (TH1D*) gDirectory->Get("fh1TrackQAPtSecRecSsc");
    hdNdetadphiTracksMCSecRecSsc = (TH2D*) gDirectory->Get("fh2TrackQAEtaPhiSecRecSsc");
  }
  
  hdNdptTracksMCPrimRec->SetDirectory(0);
  hdNdetadphiTracksMCPrimRec->SetDirectory(0);

  hdNdptTracksMCSecRecNS->SetDirectory(0);
  hdNdetadphiTracksMCSecRecNS->SetDirectory(0);

  hdNdptTracksMCSecRecSsc->SetDirectory(0);
  hdNdetadphiTracksMCSecRecSsc->SetDirectory(0);
  
  f.Close();

  // projections: dN/deta, dN/dphi 

  TH1D* hdNdetaTracksMCPrimRec   = (TH1D*) hdNdetadphiTracksMCPrimRec->ProjectionX("hdNdetaTracksMcPrimRec");
  TH1D* hdNdphiTracksMCPrimRec   = (TH1D*) hdNdetadphiTracksMCPrimRec->ProjectionY("hdNdphiTracksMcPrimRec");
 
  TH1D* hdNdetaTracksMCSecRecNS  = (TH1D*) hdNdetadphiTracksMCSecRecNS->ProjectionX("hdNdetaTracksMcSecRecNS");
  TH1D* hdNdphiTracksMCSecRecNS  = (TH1D*) hdNdetadphiTracksMCSecRecNS->ProjectionY("hdNdphiTracksMcSecRecNS");

  TH1D* hdNdetaTracksMCSecRecSsc = (TH1D*) hdNdetadphiTracksMCSecRecSsc->ProjectionX("hdNdetaTracksMcSecRecSsc");
  TH1D* hdNdphiTracksMCSecRecSsc = (TH1D*) hdNdetadphiTracksMCSecRecSsc->ProjectionY("hdNdphiTracksMcSecRecSsc");
 
  // rebin

  TString strNamePtPrim = "hTrackPtPrim";
  TString strNamePtSec  = "hTrackPtSec";

  if(fNHistoBinsSinglePt) hdNdptTracksMCPrimRec   = (TH1D*) hdNdptTracksMCPrimRec->Rebin(fNHistoBinsSinglePt,strNamePtPrim,fHistoBinsSinglePt->GetArray());
  if(fNHistoBinsSinglePt) hdNdptTracksMCSecRecNS  = (TH1D*) hdNdptTracksMCSecRecNS->Rebin(fNHistoBinsSinglePt,strNamePtSec,fHistoBinsSinglePt->GetArray());
  if(fNHistoBinsSinglePt) hdNdptTracksMCSecRecSsc = (TH1D*) hdNdptTracksMCSecRecSsc->Rebin(fNHistoBinsSinglePt,strNamePtSec,fHistoBinsSinglePt->GetArray());


  // secondary correction factor: divide prim/(prim+sec)

  TH1F* hSingleTrackSecCorrPt = (TH1F*) hdNdptTracksMCPrimRec->Clone("hSingleTrackSecCorrPt");
  TH1F* hSumPrimSecPt = (TH1F*) hdNdptTracksMCPrimRec->Clone("hSumPrimSecPt");
  hSumPrimSecPt->Add(hdNdptTracksMCSecRecNS);
  hSumPrimSecPt->Add(hdNdptTracksMCSecRecSsc);
  hSingleTrackSecCorrPt->Divide(hdNdptTracksMCPrimRec,hSumPrimSecPt,1,1,"B"); // binominal errors


  TH1F* hSingleTrackSecCorrEta = (TH1F*) hdNdetaTracksMCPrimRec->Clone("hSingleTrackSecCorrEta");
  TH1F* hSumPrimSecEta = (TH1F*) hdNdetaTracksMCPrimRec->Clone("hSumPrimSecEta");
  hSumPrimSecEta->Add(hdNdetaTracksMCSecRecNS);
  hSumPrimSecEta->Add(hdNdetaTracksMCSecRecSsc);
  hSingleTrackSecCorrEta->Divide(hdNdetaTracksMCPrimRec,hSumPrimSecEta,1,1,"B"); // binominal errors

  TH1F* hSingleTrackSecCorrPhi = (TH1F*) hdNdphiTracksMCPrimRec->Clone("hSingleTrackSecCorrPhi");
  TH1F* hSumPrimSecPhi = (TH1F*) hdNdphiTracksMCPrimRec->Clone("hSumPrimSecPhi");
  hSumPrimSecPhi->Add(hdNdphiTracksMCSecRecNS);
  hSumPrimSecPhi->Add(hdNdphiTracksMCSecRecSsc);
  hSingleTrackSecCorrPhi->Divide(hdNdphiTracksMCPrimRec,hSumPrimSecPhi,1,1,"B"); // binominal errors

  // 
  
  TString outfileOption = "RECREATE";
  if(updateOutfile)  outfileOption = "UPDATE";

  TFile out(strOutfile,outfileOption);

  if(!out.IsOpen()){
    Printf("%s:%d -- error opening secCorr output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- write secCorr to file %s ",(char*)__FILE__,__LINE__,strOutfile.Data());

  if(strOutDir && strOutDir.Length()){  

    TDirectory* dir;
    if((dir = ((TDirectory*) gDirectory->Get(strOutDir)))) dir->cd(); 
    else{
      dir = out.mkdir(strOutDir);
      dir->cd(); 
    } 
  } 

  hdNdptTracksMCSecRecNS->Write();
  hdNdetaTracksMCSecRecNS->Write();
  hdNdphiTracksMCSecRecNS->Write();

  hdNdptTracksMCSecRecSsc->Write();
  hdNdetaTracksMCSecRecSsc->Write();
  hdNdphiTracksMCSecRecSsc->Write();

  hSingleTrackSecCorrPt->Write();
  hSingleTrackSecCorrEta->Write();
  hSingleTrackSecCorrPhi->Write();
  
  out.Close();

  delete hdNdptTracksMCPrimRec;       
  delete hdNdetadphiTracksMCPrimRec;  
  delete hdNdptTracksMCSecRecNS;
  delete hdNdetadphiTracksMCSecRecNS;
  delete hdNdptTracksMCSecRecSsc;
  delete hdNdetadphiTracksMCSecRecSsc;
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteSingleResponse(TString strInfile, TString strID, TString strOutfile, 
							      Bool_t updateOutfile, TString strOutDir)
{ 
  // read task ouput from MC and write single track eff - standard dir/list 
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  WriteSingleResponse(strInfile,strdir,strlist,strOutfile,updateOutfile,strOutDir);
}


//_____________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteSingleResponse(TString strInfile, TString strdir, TString strlist,
							      TString strOutfile, Bool_t updateOutfile, TString strOutDir)
{
  // read 2d THnSparse response matrices in pt from file
  // project TH2 
  // write to strOutfile 

  THnSparse* hnResponseSinglePt; 
  TH2F*      h2ResponseSinglePt; 
 
  TFile f(strInfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read response matrices from file %s ",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
    
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  
  
  if(list) hnResponseSinglePt = (THnSparse*) list->FindObject("fhnResponseSinglePt");
  else     hnResponseSinglePt = (THnSparse*) gDirectory->Get("fhnResponseSinglePt");
  

  if(!hnResponseSinglePt){
    Printf("%s:%d -- error retrieving response matrix fhnResponseSinglePt",(char*)__FILE__,__LINE__);
    return;
  }

  f.Close();


  // project 2d histo 
  h2ResponseSinglePt = (TH2F*) hnResponseSinglePt->Projection(1,0);// note convention: yDim,xDim
  h2ResponseSinglePt->SetNameTitle("h2ResponseSinglePt",""); 
    
  
  // write 

  TString outfileOption = "RECREATE";
  if(updateOutfile)  outfileOption = "UPDATE";
  
  TFile out(strOutfile,outfileOption);
  
  if(!out.IsOpen()){
    Printf("%s:%d -- error opening response matrix output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- write response matrices to file %s ",(char*)__FILE__,__LINE__,strOutfile.Data());
  
  if(strOutDir && strOutDir.Length()){  

    TDirectory* dir;
    if((dir = ((TDirectory*)  gDirectory->Get(strOutDir)))) dir->cd(); 
    else{
      dir = out.mkdir(strOutDir);
      dir->cd(); 
    } 
  }
  
  hnResponseSinglePt->Write();
  h2ResponseSinglePt->Write();
  
  out.Close();  
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteJetTrackEff(TString strInfile, TString strID, TString strOutfile, 
							   Bool_t updateOutfile)
{ 
  // read task ouput from MC and write single track eff - standard dir/list 
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  WriteJetTrackEff(strInfile,strdir,strlist,strOutfile,updateOutfile);
}

//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteJetTrackEff(TString strInfile, TString strdir, TString strlist, 
							   TString strOutfile, Bool_t updateOutfile)
{
  // read task output from MC and write track eff in jet pt slices as function of pt, z, xi
  // argument strlist optional - read from directory strdir if not specified
  // write eff to file strOutfile - by default only 'update' file (don't overwrite)

  TH1F* hEffPt[fNJetPtSlices];
  TH1F* hEffXi[fNJetPtSlices];
  TH1F* hEffZ[fNJetPtSlices];

  TH1F* hdNdptTracksMCPrimGen[fNJetPtSlices];
  TH1F* hdNdxiMCPrimGen[fNJetPtSlices];
  TH1F* hdNdzMCPrimGen[fNJetPtSlices];
    
  TH1F* hdNdptTracksMCPrimRec[fNJetPtSlices];
  TH1F* hdNdxiMCPrimRec[fNJetPtSlices];
  TH1F* hdNdzMCPrimRec[fNJetPtSlices];


  TH1F* fh1FFJetPtRecEffGen;

  TH2F* fh2FFTrackPtRecEffGen;
  TH2F* fh2FFZRecEffGen;
  TH2F* fh2FFXiRecEffGen;
  
  TH2F* fh2FFTrackPtRecEffRec;
  TH2F* fh2FFZRecEffRec;
  TH2F* fh2FFXiRecEffRec;
 

  TFile f(strInfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- writeJetTrackEff: open task ouput file %s",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){

    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }

  if(list){
    fh1FFJetPtRecEffGen    = (TH1F*) list->FindObject("fh1FFJetPtRecEffGen");

    fh2FFTrackPtRecEffGen  = (TH2F*) list->FindObject("fh2FFTrackPtRecEffGen");
    fh2FFZRecEffGen        = (TH2F*) list->FindObject("fh2FFZRecEffGen");
    fh2FFXiRecEffGen       = (TH2F*) list->FindObject("fh2FFXiRecEffGen");
    
    fh2FFTrackPtRecEffRec  = (TH2F*) list->FindObject("fh2FFTrackPtRecEffRec");
    fh2FFZRecEffRec        = (TH2F*) list->FindObject("fh2FFZRecEffRec");
    fh2FFXiRecEffRec       = (TH2F*) list->FindObject("fh2FFXiRecEffRec");
  }
  else{
    fh1FFJetPtRecEffGen    = (TH1F*) gDirectory->Get("fh1FFJetPtRecEffGen");

    fh2FFTrackPtRecEffGen  = (TH2F*) gDirectory->Get("fh2FFTrackPtRecEffGen");
    fh2FFZRecEffGen        = (TH2F*) gDirectory->Get("fh2FFZRecEffGen");
    fh2FFXiRecEffGen       = (TH2F*) gDirectory->Get("fh2FFXiRecEffGen");
    
    fh2FFTrackPtRecEffRec  = (TH2F*) gDirectory->Get("fh2FFTrackPtRecEffRec");
    fh2FFZRecEffRec        = (TH2F*) gDirectory->Get("fh2FFZRecEffRec");
    fh2FFXiRecEffRec       = (TH2F*) gDirectory->Get("fh2FFXiRecEffRec");
  }
  
  fh1FFJetPtRecEffGen->SetDirectory(0); 

  fh2FFTrackPtRecEffGen->SetDirectory(0);
  fh2FFZRecEffGen->SetDirectory(0);
  fh2FFXiRecEffGen->SetDirectory(0);
  
  fh2FFTrackPtRecEffRec->SetDirectory(0);
  fh2FFZRecEffRec->SetDirectory(0);
  fh2FFXiRecEffRec->SetDirectory(0);
  
  f.Close();


  // projections: FF for generated and reconstructed primaries 
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);

    Int_t binLo = static_cast<Int_t>(fh2FFTrackPtRecEffGen->GetXaxis()->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh2FFTrackPtRecEffGen->GetXaxis()->FindBin(jetPtUpLim))-1;

    if(binUp > fh2FFTrackPtRecEffGen->GetNbinsX()){
      Printf("%s:%d -- jet pt range %0.3f exceeds histo limits",(char*)__FILE__,__LINE__,jetPtUpLim); 
      return; 
    }
    
    TString strNameFFPtGen(Form("fh1FFTrackPtGenPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZGen(Form("fh1FFZGenPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXiGen(Form("fh1FFXiGenPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    
    TString strNameFFPtRec(Form("fh1FFTrackPtRecPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZRec(Form("fh1FFZRecPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXiRec(Form("fh1FFXiRecPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    
    // project 
    // appendix 'unbinned' to avoid histos with same name after rebinning

    hdNdptTracksMCPrimGen[i] = (TH1F*) fh2FFTrackPtRecEffGen->ProjectionY(strNameFFPtGen+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    hdNdzMCPrimGen[i]        = (TH1F*) fh2FFZRecEffGen->ProjectionY(strNameFFZGen+"_unBinned",binLo,binUp,"o");
    hdNdxiMCPrimGen[i]       = (TH1F*) fh2FFXiRecEffGen->ProjectionY(strNameFFXiGen+"_unBinned",binLo,binUp,"o");
    
    hdNdptTracksMCPrimRec[i] = (TH1F*) fh2FFTrackPtRecEffRec->ProjectionY(strNameFFPtRec+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    hdNdzMCPrimRec[i]        = (TH1F*) fh2FFZRecEffRec->ProjectionY(strNameFFZRec+"_unBinned",binLo,binUp,"o");
    hdNdxiMCPrimRec[i]       = (TH1F*) fh2FFXiRecEffRec->ProjectionY(strNameFFXiRec+"_unBinned",binLo,binUp,"o");
    
    // rebin

    if(fNHistoBinsPt[i]) hdNdptTracksMCPrimGen[i] = (TH1F*) hdNdptTracksMCPrimGen[i]->Rebin(fNHistoBinsPt[i],strNameFFPtGen,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hdNdzMCPrimGen[i]  = (TH1F*) hdNdzMCPrimGen[i]->Rebin(fNHistoBinsZ[i],strNameFFZGen,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hdNdxiMCPrimGen[i] = (TH1F*) hdNdxiMCPrimGen[i]->Rebin(fNHistoBinsXi[i],strNameFFXiGen,fHistoBinsXi[i]->GetArray());

    if(fNHistoBinsPt[i]) hdNdptTracksMCPrimRec[i] = (TH1F*) hdNdptTracksMCPrimRec[i]->Rebin(fNHistoBinsPt[i],strNameFFPtRec,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hdNdzMCPrimRec[i]  = (TH1F*) hdNdzMCPrimRec[i]->Rebin(fNHistoBinsZ[i],strNameFFZRec,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hdNdxiMCPrimRec[i] = (TH1F*) hdNdxiMCPrimRec[i]->Rebin(fNHistoBinsXi[i],strNameFFXiRec,fHistoBinsXi[i]->GetArray());

    hdNdptTracksMCPrimGen[i]->SetNameTitle(strNameFFPtGen,"");
    hdNdzMCPrimGen[i]->SetNameTitle(strNameFFZGen,"");
    hdNdxiMCPrimGen[i]->SetNameTitle(strNameFFXiGen,"");
    
    hdNdptTracksMCPrimRec[i]->SetNameTitle(strNameFFPtRec,"");
    hdNdzMCPrimRec[i]->SetNameTitle(strNameFFZRec,"");
    hdNdxiMCPrimRec[i]->SetNameTitle(strNameFFXiRec,"");
 
    // normalize
    
    Double_t nJetsBin = fh1FFJetPtRecEffGen->Integral(binLo,binUp);

    NormalizeTH1(hdNdptTracksMCPrimGen[i],nJetsBin); 
    NormalizeTH1(hdNdzMCPrimGen[i],nJetsBin); 
    NormalizeTH1(hdNdxiMCPrimGen[i],nJetsBin); 

    NormalizeTH1(hdNdptTracksMCPrimRec[i],nJetsBin); 
    NormalizeTH1(hdNdzMCPrimRec[i],nJetsBin); 
    NormalizeTH1(hdNdxiMCPrimRec[i],nJetsBin); 
    
    // divide rec/gen : efficiency

    TString strNameEffPt(Form("hEffPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameEffZ(Form("hEffZ_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameEffXi(Form("hEffXi_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
 
    hEffPt[i] = (TH1F*) hdNdptTracksMCPrimRec[i]->Clone(strNameEffPt);
    hEffPt[i]->Divide(hdNdptTracksMCPrimRec[i],hdNdptTracksMCPrimGen[i],1,1,"B"); // binominal errors
    
    hEffXi[i] = (TH1F*) hdNdxiMCPrimRec[i]->Clone(strNameEffXi);
    hEffXi[i]->Divide(hdNdxiMCPrimRec[i],hdNdxiMCPrimGen[i],1,1,"B"); // binominal errors
    
    hEffZ[i] = (TH1F*) hdNdzMCPrimRec[i]->Clone(strNameEffZ);
    hEffZ[i]->Divide(hdNdzMCPrimRec[i],hdNdzMCPrimGen[i],1,1,"B"); // binominal errors
  } 
  
  // write 

  TString outfileOption = "RECREATE";
  if(updateOutfile)  outfileOption = "UPDATE";

  TFile out(strOutfile,outfileOption);
  
  if(!out.IsOpen()){
    Printf("%s:%d -- error opening efficiency output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- write efficiency to file %s ",(char*)__FILE__,__LINE__,strOutfile.Data());

//   if(strdir && strdir.Length()){ 
//     TDirectory* dir = out.mkdir(strdir);
//     dir->cd(); 
//   }

  for(Int_t i=0; i<fNJetPtSlices; i++){

    hdNdptTracksMCPrimGen[i]->Write();
    hdNdxiMCPrimGen[i]->Write();
    hdNdzMCPrimGen[i]->Write();
    
    hdNdptTracksMCPrimRec[i]->Write();
    hdNdxiMCPrimRec[i]->Write();
    hdNdzMCPrimRec[i]->Write();
  
    hEffPt[i]->Write();
    hEffXi[i]->Write();
    hEffZ[i]->Write();
  }

  out.Close();

 
  delete fh1FFJetPtRecEffGen; 

  delete fh2FFTrackPtRecEffGen;
  delete fh2FFZRecEffGen;
  delete fh2FFXiRecEffGen;
  
  delete fh2FFTrackPtRecEffRec;
  delete fh2FFZRecEffRec;
  delete fh2FFXiRecEffRec;
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteJetSecCorr(TString strInfile, TString strID, TString strOutfile, 
							  Bool_t updateOutfile, TString strOutDir)
{ 
  // read task ouput from MC and write secondary correction - standard dir/list 
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  WriteJetSecCorr(strInfile,strdir,strlist,strOutfile,updateOutfile, strOutDir,kFALSE,"");
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteBgrJetSecCorr(TString strInfile, TString strBgrID, TString strID, TString strOutfile, 
							     Bool_t updateOutfile, TString strOutDir)
{ 
  // read task ouput from MC and write secondary correction - standard dir/list 
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  WriteJetSecCorr(strInfile,strdir,strlist,strOutfile,updateOutfile, strOutDir,kTRUE,strBgrID);
}


//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteJetSecCorr(TString strInfile, TString strdir, TString strlist, 
							  TString strOutfile, Bool_t updateOutfile, TString strOutDir,
							  Bool_t writeBgr, TString strBgrID)
{
  // read task output from MC and write secondary correction in jet pt slices as function of pt, z, xi
  // argument strlist optional - read from directory strdir if not specified
  // write eff to file strOutfile - by default only 'update' file (don't overwrite)

  TH1F* hSecCorrPt[fNJetPtSlices];
  TH1F* hSecCorrXi[fNJetPtSlices];
  TH1F* hSecCorrZ[fNJetPtSlices];

  // corr factors using naive Pythia strangeness - for sys uncertainties
  TH1F* hSecCorrPt_nonSc[fNJetPtSlices];
  TH1F* hSecCorrXi_nonSc[fNJetPtSlices];
  TH1F* hSecCorrZ_nonSc[fNJetPtSlices];

  TH1F* hdNdptTracksMCPrimRec[fNJetPtSlices];
  TH1F* hdNdxiMCPrimRec[fNJetPtSlices];
  TH1F* hdNdzMCPrimRec[fNJetPtSlices];
    
  TH1F* hdNdptTracksMCSecRecNS[fNJetPtSlices];
  TH1F* hdNdxiMCSecRecNS[fNJetPtSlices];
  TH1F* hdNdzMCSecRecNS[fNJetPtSlices];

  TH1F* hdNdptTracksMCSecRecS[fNJetPtSlices];
  TH1F* hdNdxiMCSecRecS[fNJetPtSlices];
  TH1F* hdNdzMCSecRecS[fNJetPtSlices];

  TH1F* hdNdptTracksMCSecRecSsc[fNJetPtSlices];
  TH1F* hdNdxiMCSecRecSsc[fNJetPtSlices];
  TH1F* hdNdzMCSecRecSsc[fNJetPtSlices];
  
  // ---

  TH1F* fh1FFJetPtRecEffRec;

  TH2F* fh2FFTrackPtRecEffRec;
  TH2F* fh2FFZRecEffRec;
  TH2F* fh2FFXiRecEffRec;
  
  TH2F* fh2FFTrackPtSecRecNS;
  TH2F* fh2FFZSecRecNS;
  TH2F* fh2FFXiSecRecNS;

  TH2F* fh2FFTrackPtSecRecS;
  TH2F* fh2FFZSecRecS;
  TH2F* fh2FFXiSecRecS;

  TH2F* fh2FFTrackPtSecRecSsc;
  TH2F* fh2FFZSecRecSsc;
  TH2F* fh2FFXiSecRecSsc;
  
  // ---

  TFile f(strInfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- writeJetTrackSecCorr: open task ouput file %s",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){

    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }

  if(list){
    fh1FFJetPtRecEffRec    = (TH1F*) list->FindObject("fh1FFJetPtRecEffRec");

    if(writeBgr){ 
      fh2FFTrackPtRecEffRec  = (TH2F*) list->FindObject(Form("fh2FFTrackPt%sRecEffRec",strBgrID.Data()));
      fh2FFZRecEffRec        = (TH2F*) list->FindObject(Form("fh2FFZ%sRecEffRec",strBgrID.Data()));
      fh2FFXiRecEffRec       = (TH2F*) list->FindObject(Form("fh2FFXi%sRecEffRec",strBgrID.Data())); 

      fh2FFTrackPtSecRecNS   = (TH2F*) list->FindObject(Form("fh2FFTrackPt%sSecRecNS",strBgrID.Data()));
      fh2FFZSecRecNS         = (TH2F*) list->FindObject(Form("fh2FFZ%sSecRecNS",strBgrID.Data()));
      fh2FFXiSecRecNS        = (TH2F*) list->FindObject(Form("fh2FFXi%sSecRecNS",strBgrID.Data()));

      fh2FFTrackPtSecRecS    = (TH2F*) list->FindObject(Form("fh2FFTrackPt%sSecRecS",strBgrID.Data()));
      fh2FFZSecRecS          = (TH2F*) list->FindObject(Form("fh2FFZ%sSecRecS",strBgrID.Data()));
      fh2FFXiSecRecS         = (TH2F*) list->FindObject(Form("fh2FFXi%sSecRecS",strBgrID.Data()));
      
      fh2FFTrackPtSecRecSsc  = (TH2F*) list->FindObject(Form("fh2FFTrackPt%sSecRecSsc",strBgrID.Data()));
      fh2FFZSecRecSsc        = (TH2F*) list->FindObject(Form("fh2FFZ%sSecRecSsc",strBgrID.Data()));
      fh2FFXiSecRecSsc       = (TH2F*) list->FindObject(Form("fh2FFXi%sSecRecSsc",strBgrID.Data()));
    }
    else{
      fh2FFTrackPtRecEffRec  = (TH2F*) list->FindObject("fh2FFTrackPtRecEffRec");
      fh2FFZRecEffRec        = (TH2F*) list->FindObject("fh2FFZRecEffRec");
      fh2FFXiRecEffRec       = (TH2F*) list->FindObject("fh2FFXiRecEffRec");
      
      fh2FFTrackPtSecRecNS   = (TH2F*) list->FindObject("fh2FFTrackPtSecRecNS");
      fh2FFZSecRecNS         = (TH2F*) list->FindObject("fh2FFZSecRecNS");
      fh2FFXiSecRecNS        = (TH2F*) list->FindObject("fh2FFXiSecRecNS");
      
      fh2FFTrackPtSecRecS    = (TH2F*) list->FindObject("fh2FFTrackPtSecRecS");
      fh2FFZSecRecS          = (TH2F*) list->FindObject("fh2FFZSecRecS");
      fh2FFXiSecRecS         = (TH2F*) list->FindObject("fh2FFXiSecRecS");
      
      fh2FFTrackPtSecRecSsc  = (TH2F*) list->FindObject("fh2FFTrackPtSecRecSsc");
      fh2FFZSecRecSsc        = (TH2F*) list->FindObject("fh2FFZSecRecSsc");
      fh2FFXiSecRecSsc       = (TH2F*) list->FindObject("fh2FFXiSecRecSsc");
    }
  }
  else{
    fh1FFJetPtRecEffRec    = (TH1F*) gDirectory->Get("fh1FFJetPtRecEffRec");

    if(writeBgr){
      fh2FFTrackPtRecEffRec  = (TH2F*) gDirectory->Get(Form("fh2FFTrackPt%sRecEffRec",strBgrID.Data()));
      fh2FFZRecEffRec        = (TH2F*) gDirectory->Get(Form("fh2FFZ%sRecEffRec",strBgrID.Data()));
      fh2FFXiRecEffRec       = (TH2F*) gDirectory->Get(Form("fh2FFXi%sRecEffRec",strBgrID.Data())); 

      fh2FFTrackPtSecRecNS   = (TH2F*) gDirectory->Get(Form("fh2FFTrackPt%sSecRecNS",strBgrID.Data()));
      fh2FFZSecRecNS         = (TH2F*) gDirectory->Get(Form("fh2FFZ%sSecRecNS",strBgrID.Data()));
      fh2FFXiSecRecNS        = (TH2F*) gDirectory->Get(Form("fh2FFXi%sSecRecNS",strBgrID.Data()));

      fh2FFTrackPtSecRecS    = (TH2F*) gDirectory->Get(Form("fh2FFTrackPt%sSecRecS",strBgrID.Data()));
      fh2FFZSecRecS          = (TH2F*) gDirectory->Get(Form("fh2FFZ%sSecRecS",strBgrID.Data()));
      fh2FFXiSecRecS         = (TH2F*) gDirectory->Get(Form("fh2FFXi%sSecRecS",strBgrID.Data()));
      
      fh2FFTrackPtSecRecSsc  = (TH2F*) gDirectory->Get(Form("fh2FFTrackPt%sSecRecSsc",strBgrID.Data()));
      fh2FFZSecRecSsc        = (TH2F*) gDirectory->Get(Form("fh2FFZ%sSecRecSsc",strBgrID.Data()));
      fh2FFXiSecRecSsc       = (TH2F*) gDirectory->Get(Form("fh2FFXi%sSecRecSsc",strBgrID.Data()));
    }
    else{
      fh2FFTrackPtRecEffRec  = (TH2F*) gDirectory->Get("fh2FFTrackPtRecEffRec");
      fh2FFZRecEffRec        = (TH2F*) gDirectory->Get("fh2FFZRecEffRec");
      fh2FFXiRecEffRec       = (TH2F*) gDirectory->Get("fh2FFXiRecEffRec");
    
      fh2FFTrackPtSecRecNS   = (TH2F*) gDirectory->Get("fh2FFTrackPtSecRecNS");
      fh2FFZSecRecNS         = (TH2F*) gDirectory->Get("fh2FFZSecRecNS");
      fh2FFXiSecRecNS        = (TH2F*) gDirectory->Get("fh2FFXiSecRecNS");
      
      fh2FFTrackPtSecRecS    = (TH2F*) gDirectory->Get("fh2FFTrackPtSecRecS");
      fh2FFZSecRecS          = (TH2F*) gDirectory->Get("fh2FFZSecRecS");
      fh2FFXiSecRecS         = (TH2F*) gDirectory->Get("fh2FFXiSecRecS");
      
      fh2FFTrackPtSecRecSsc  = (TH2F*) gDirectory->Get("fh2FFTrackPtSecRecSsc");
      fh2FFZSecRecSsc        = (TH2F*) gDirectory->Get("fh2FFZSecRecSsc");
      fh2FFXiSecRecSsc       = (TH2F*) gDirectory->Get("fh2FFXiSecRecSsc");
    }
  }
  
  fh1FFJetPtRecEffRec->SetDirectory(0); 
   
  fh2FFTrackPtRecEffRec->SetDirectory(0);
  fh2FFZRecEffRec->SetDirectory(0);
  fh2FFXiRecEffRec->SetDirectory(0);
  
  fh2FFTrackPtSecRecNS->SetDirectory(0);
  fh2FFZSecRecNS->SetDirectory(0);
  fh2FFXiSecRecNS->SetDirectory(0);
  
  fh2FFTrackPtSecRecS->SetDirectory(0);
  fh2FFZSecRecS->SetDirectory(0);
  fh2FFXiSecRecS->SetDirectory(0);
  
  fh2FFTrackPtSecRecSsc->SetDirectory(0);
  fh2FFZSecRecSsc->SetDirectory(0);
  fh2FFXiSecRecSsc->SetDirectory(0);

  f.Close();


  // projections: FF for generated and reconstructed primaries 
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);

    Int_t binLo = static_cast<Int_t>(fh2FFTrackPtRecEffRec->GetXaxis()->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh2FFTrackPtRecEffRec->GetXaxis()->FindBin(jetPtUpLim))-1;

    if(binUp > fh2FFTrackPtRecEffRec->GetNbinsX()){
      Printf("%s:%d -- jet pt range %0.3f exceeds histo limits",(char*)__FILE__,__LINE__,jetPtUpLim); 
      return; 
    }
    
    TString strNameFFPtPrimRec(Form("fh1FFTrackPtRecPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZPrimRec(Form("fh1FFZRecPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXiPrimRec(Form("fh1FFXiRecPrim_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    
    TString strNameFFPtSecRecNS(Form("fh1FFTrackPtRecSecNS_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZSecRecNS(Form("fh1FFZRecSecNS_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXiSecRecNS(Form("fh1FFXiRecSecNS_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
     
    TString strNameFFPtSecRecS(Form("fh1FFTrackPtRecSecS_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZSecRecS(Form("fh1FFZRecSecS_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXiSecRecS(Form("fh1FFXiRecSecS_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
     
    TString strNameFFPtSecRecSsc(Form("fh1FFTrackPtRecSecSsc_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZSecRecSsc(Form("fh1FFZRecSecSsc_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXiSecRecSsc(Form("fh1FFXiRecSecSsc_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));

    if(writeBgr){
      strNameFFPtPrimRec.Form("fh1BgrTrackPt%sRecPrim_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      strNameFFZPrimRec.Form("fh1BgrZ%sRecPrim_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      strNameFFXiPrimRec.Form("fh1BgrXi%sRecPrim_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      
      strNameFFPtSecRecNS.Form("fh1BgrTrackPt%sRecSecNS_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      strNameFFZSecRecNS.Form("fh1BgrZ%sRecSecNS_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      strNameFFXiSecRecNS.Form("fh1BgrXi%sRecSecNS_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      
      strNameFFPtSecRecS.Form("fh1BgrTrackPt%sRecSecS_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      strNameFFZSecRecS.Form("fh1BgrZ%sRecSecS_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      strNameFFXiSecRecS.Form("fh1BgrXi%sRecSecS_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      
      strNameFFPtSecRecSsc.Form("fh1BgrTrackPt%sRecSecSsc_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      strNameFFZSecRecSsc.Form("fh1BgrZ%sRecSecSsc_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
      strNameFFXiSecRecSsc.Form("fh1BgrXi%sRecSecSsc_%02d_%02d",strBgrID.Data(),static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim));
    }
    
    // project 
    // appendix 'unbinned' to avoid histos with same name after rebinning

    hdNdptTracksMCPrimRec[i] = (TH1F*) fh2FFTrackPtRecEffRec->ProjectionY(strNameFFPtPrimRec+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    hdNdzMCPrimRec[i]        = (TH1F*) fh2FFZRecEffRec->ProjectionY(strNameFFZPrimRec+"_unBinned",binLo,binUp,"o");
    hdNdxiMCPrimRec[i]       = (TH1F*) fh2FFXiRecEffRec->ProjectionY(strNameFFXiPrimRec+"_unBinned",binLo,binUp,"o");
    
    hdNdptTracksMCSecRecNS[i]  = (TH1F*) fh2FFTrackPtSecRecNS->ProjectionY(strNameFFPtSecRecNS+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    hdNdzMCSecRecNS[i]         = (TH1F*) fh2FFZSecRecNS->ProjectionY(strNameFFZSecRecNS+"_unBinned",binLo,binUp,"o");
    hdNdxiMCSecRecNS[i]        = (TH1F*) fh2FFXiSecRecNS->ProjectionY(strNameFFXiSecRecNS+"_unBinned",binLo,binUp,"o");
    
    hdNdptTracksMCSecRecS[i]  = (TH1F*) fh2FFTrackPtSecRecS->ProjectionY(strNameFFPtSecRecS+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    hdNdzMCSecRecS[i]         = (TH1F*) fh2FFZSecRecS->ProjectionY(strNameFFZSecRecS+"_unBinned",binLo,binUp,"o");
    hdNdxiMCSecRecS[i]        = (TH1F*) fh2FFXiSecRecS->ProjectionY(strNameFFXiSecRecS+"_unBinned",binLo,binUp,"o");

    hdNdptTracksMCSecRecSsc[i]  = (TH1F*) fh2FFTrackPtSecRecSsc->ProjectionY(strNameFFPtSecRecSsc+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    hdNdzMCSecRecSsc[i]         = (TH1F*) fh2FFZSecRecSsc->ProjectionY(strNameFFZSecRecSsc+"_unBinned",binLo,binUp,"o");
    hdNdxiMCSecRecSsc[i]        = (TH1F*) fh2FFXiSecRecSsc->ProjectionY(strNameFFXiSecRecSsc+"_unBinned",binLo,binUp,"o");


    // rebin
    if(fNHistoBinsPt[i]) hdNdptTracksMCPrimRec[i] = (TH1F*) hdNdptTracksMCPrimRec[i]->Rebin(fNHistoBinsPt[i],strNameFFPtPrimRec,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hdNdzMCPrimRec[i]        = (TH1F*) hdNdzMCPrimRec[i]->Rebin(fNHistoBinsZ[i],strNameFFZPrimRec,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hdNdxiMCPrimRec[i]       = (TH1F*) hdNdxiMCPrimRec[i]->Rebin(fNHistoBinsXi[i],strNameFFXiPrimRec,fHistoBinsXi[i]->GetArray());

    if(fNHistoBinsPt[i]) hdNdptTracksMCSecRecNS[i] = (TH1F*) hdNdptTracksMCSecRecNS[i]->Rebin(fNHistoBinsPt[i],strNameFFPtSecRecNS,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hdNdzMCSecRecNS[i]        = (TH1F*) hdNdzMCSecRecNS[i]->Rebin(fNHistoBinsZ[i],strNameFFZSecRecNS,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hdNdxiMCSecRecNS[i]       = (TH1F*) hdNdxiMCSecRecNS[i]->Rebin(fNHistoBinsXi[i],strNameFFXiSecRecNS,fHistoBinsXi[i]->GetArray());

    if(fNHistoBinsPt[i]) hdNdptTracksMCSecRecS[i] = (TH1F*) hdNdptTracksMCSecRecS[i]->Rebin(fNHistoBinsPt[i],strNameFFPtSecRecS,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hdNdzMCSecRecS[i]        = (TH1F*) hdNdzMCSecRecS[i]->Rebin(fNHistoBinsZ[i],strNameFFZSecRecS,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hdNdxiMCSecRecS[i]       = (TH1F*) hdNdxiMCSecRecS[i]->Rebin(fNHistoBinsXi[i],strNameFFXiSecRecS,fHistoBinsXi[i]->GetArray());

    if(fNHistoBinsPt[i]) hdNdptTracksMCSecRecSsc[i] = (TH1F*) hdNdptTracksMCSecRecSsc[i]->Rebin(fNHistoBinsPt[i],strNameFFPtSecRecSsc,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hdNdzMCSecRecSsc[i]        = (TH1F*) hdNdzMCSecRecSsc[i]->Rebin(fNHistoBinsZ[i],strNameFFZSecRecSsc,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hdNdxiMCSecRecSsc[i]       = (TH1F*) hdNdxiMCSecRecSsc[i]->Rebin(fNHistoBinsXi[i],strNameFFXiSecRecSsc,fHistoBinsXi[i]->GetArray());


    hdNdptTracksMCPrimRec[i]->SetNameTitle(strNameFFPtPrimRec,"");
    hdNdzMCPrimRec[i]->SetNameTitle(strNameFFZPrimRec,"");
    hdNdxiMCPrimRec[i]->SetNameTitle(strNameFFXiPrimRec,"");
    
    hdNdptTracksMCSecRecNS[i]->SetNameTitle(strNameFFPtSecRecNS,"");
    hdNdzMCSecRecNS[i]->SetNameTitle(strNameFFZSecRecNS,"");
    hdNdxiMCSecRecNS[i]->SetNameTitle(strNameFFXiSecRecNS,"");
    
    hdNdptTracksMCSecRecS[i]->SetNameTitle(strNameFFPtSecRecS,"");
    hdNdzMCSecRecS[i]->SetNameTitle(strNameFFZSecRecS,"");
    hdNdxiMCSecRecS[i]->SetNameTitle(strNameFFXiSecRecS,"");

    hdNdptTracksMCSecRecSsc[i]->SetNameTitle(strNameFFPtSecRecSsc,"");
    hdNdzMCSecRecSsc[i]->SetNameTitle(strNameFFZSecRecSsc,"");
    hdNdxiMCSecRecSsc[i]->SetNameTitle(strNameFFXiSecRecSsc,"");
    
    // normalize
    Double_t nJetsBin = fh1FFJetPtRecEffRec->Integral(binLo,binUp);

    NormalizeTH1(hdNdptTracksMCPrimRec[i],nJetsBin); 
    NormalizeTH1(hdNdzMCPrimRec[i],nJetsBin); 
    NormalizeTH1(hdNdxiMCPrimRec[i],nJetsBin); 

    NormalizeTH1(hdNdptTracksMCSecRecNS[i],nJetsBin); 
    NormalizeTH1(hdNdzMCSecRecNS[i],nJetsBin); 
    NormalizeTH1(hdNdxiMCSecRecNS[i],nJetsBin); 

    NormalizeTH1(hdNdptTracksMCSecRecS[i],nJetsBin); 
    NormalizeTH1(hdNdzMCSecRecS[i],nJetsBin); 
    NormalizeTH1(hdNdxiMCSecRecS[i],nJetsBin); 

    NormalizeTH1(hdNdptTracksMCSecRecSsc[i],nJetsBin); 
    NormalizeTH1(hdNdzMCSecRecSsc[i],nJetsBin); 
    NormalizeTH1(hdNdxiMCSecRecSsc[i],nJetsBin); 

    
    // divide prim / (prim+sec) : corr factor
    TString strNameSecCorrPt(Form("hSecCorrPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameSecCorrZ(Form("hSecCorrZ_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameSecCorrXi(Form("hSecCorrXi_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    
    if(writeBgr){
      strNameSecCorrPt.Form("hSecCorrBgrPt_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
      strNameSecCorrZ.Form("hSecCorrBgrZ_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
      strNameSecCorrXi.Form("hSecCorrBgrXi_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
    }
     
    hSecCorrPt[i] = (TH1F*) hdNdptTracksMCPrimRec[i]->Clone(strNameSecCorrPt);
    TH1F* hSumPrimSecPt = (TH1F*) hdNdptTracksMCPrimRec[i]->Clone("hSumPrimSecPt");
    hSumPrimSecPt->Add(hdNdptTracksMCSecRecNS[i]);
    hSumPrimSecPt->Add(hdNdptTracksMCSecRecSsc[i]);
    hSecCorrPt[i]->Divide(hdNdptTracksMCPrimRec[i],hSumPrimSecPt,1,1,"B"); // binominal errors

    hSecCorrXi[i] = (TH1F*) hdNdxiMCPrimRec[i]->Clone(strNameSecCorrXi);
    TH1F* hSumPrimSecXi = (TH1F*) hdNdxiMCPrimRec[i]->Clone("hSumPrimSecXi");
    hSumPrimSecXi->Add(hdNdxiMCSecRecNS[i]);
    hSumPrimSecXi->Add(hdNdxiMCSecRecSsc[i]);
    hSecCorrXi[i]->Divide(hdNdxiMCPrimRec[i],hSumPrimSecXi,1,1,"B"); // binominal errors

    hSecCorrZ[i] = (TH1F*) hdNdzMCPrimRec[i]->Clone(strNameSecCorrZ);
    TH1F* hSumPrimSecZ = (TH1F*) hdNdzMCPrimRec[i]->Clone("hSumPrimSecZ");
    hSumPrimSecZ->Add(hdNdzMCSecRecNS[i]);
    hSumPrimSecZ->Add(hdNdzMCSecRecSsc[i]);
    hSecCorrZ[i]->Divide(hdNdzMCPrimRec[i],hSumPrimSecZ,1,1,"B"); // binominal errors

    // the same using unscaled strangeness
    TString strNameSecCorrPt_nonSc(Form("hSecCorrPt_nonSc_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameSecCorrZ_nonSc(Form("hSecCorrZ_nonSc_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameSecCorrXi_nonSc(Form("hSecCorrXi_nonSc_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));

    if(writeBgr){
      strNameSecCorrPt_nonSc.Form("hSecCorrBgrPt_nonSc_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
      strNameSecCorrZ_nonSc.Form("hSecCorrBgrZ_nonSc_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
      strNameSecCorrXi_nonSc.Form("hSecCorrBgrXi_nonSc_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
    }

    hSecCorrPt_nonSc[i] = (TH1F*) hdNdptTracksMCPrimRec[i]->Clone(strNameSecCorrPt_nonSc);
    TH1F* hSumPrimSecPt_nonSc = (TH1F*) hdNdptTracksMCPrimRec[i]->Clone("hSumPrimSecPt_nonSc");
    hSumPrimSecPt_nonSc->Add(hdNdptTracksMCSecRecNS[i]);
    hSumPrimSecPt_nonSc->Add(hdNdptTracksMCSecRecS[i]); // non-scaled secondaries from strangeness
    hSecCorrPt_nonSc[i]->Divide(hdNdptTracksMCPrimRec[i],hSumPrimSecPt_nonSc,1,1,"B"); // binominal errors
    
    hSecCorrZ_nonSc[i] = (TH1F*) hdNdzMCPrimRec[i]->Clone(strNameSecCorrZ_nonSc);
    TH1F* hSumPrimSecZ_nonSc = (TH1F*) hdNdzMCPrimRec[i]->Clone("hSumPrimSecZ_nonSc");
    hSumPrimSecZ_nonSc->Add(hdNdzMCSecRecNS[i]);
    hSumPrimSecZ_nonSc->Add(hdNdzMCSecRecS[i]); // non-scaled secondaries from strangeness
    hSecCorrZ_nonSc[i]->Divide(hdNdzMCPrimRec[i],hSumPrimSecZ_nonSc,1,1,"B"); // binominal errors
    
    hSecCorrXi_nonSc[i] = (TH1F*) hdNdxiMCPrimRec[i]->Clone(strNameSecCorrXi_nonSc);
    TH1F* hSumPrimSecXi_nonSc = (TH1F*) hdNdxiMCPrimRec[i]->Clone("hSumPrimSecXi_nonSc");
    hSumPrimSecXi_nonSc->Add(hdNdxiMCSecRecNS[i]);
    hSumPrimSecXi_nonSc->Add(hdNdxiMCSecRecS[i]); // non-scaled secondaries from strangeness
    hSecCorrXi_nonSc[i]->Divide(hdNdxiMCPrimRec[i],hSumPrimSecXi_nonSc,1,1,"B"); // binominal errors
  } 
  
  // write 

  TString outfileOption = "RECREATE";
  if(updateOutfile)  outfileOption = "UPDATE";

  TFile out(strOutfile,outfileOption);
  
  if(!out.IsOpen()){
    Printf("%s:%d -- error opening sec corr output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- write jet sec corr to file %s ",(char*)__FILE__,__LINE__,strOutfile.Data());

  if(strOutDir && strOutDir.Length()){  

    TDirectory* dir;
    if((dir = ((TDirectory*)  gDirectory->Get(strOutDir)))) dir->cd(); 
    else{
      dir = out.mkdir(strOutDir);
      dir->cd(); 
    } 
  }
  

  for(Int_t i=0; i<fNJetPtSlices; i++){
  
    hSecCorrPt[i]->Write();
    hSecCorrXi[i]->Write();
    hSecCorrZ[i]->Write();
  
    hSecCorrPt_nonSc[i]->Write();
    hSecCorrXi_nonSc[i]->Write();
    hSecCorrZ_nonSc[i]->Write();
  }

  TString strSpectraDir      = "spectraSecCorr";
  if(writeBgr) strSpectraDir = "spectraBgrSecCorr";
  
  TDirectory *dOut  = gDirectory->mkdir(strSpectraDir);
  dOut->cd();
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    hdNdptTracksMCPrimRec[i]->Write(); 
    hdNdzMCPrimRec[i]->Write(); 
    hdNdxiMCPrimRec[i]->Write(); 
    
    hdNdptTracksMCSecRecNS[i]->Write(); 
    hdNdzMCSecRecNS[i]->Write(); 
    hdNdxiMCSecRecNS[i]->Write(); 
    
    hdNdptTracksMCSecRecS[i]->Write(); 
    hdNdzMCSecRecS[i]->Write(); 
    hdNdxiMCSecRecS[i]->Write(); 
    
    hdNdptTracksMCSecRecSsc[i]->Write(); 
    hdNdzMCSecRecSsc[i]->Write(); 
    hdNdxiMCSecRecSsc[i]->Write();  
  }

  out.Close();
  
  delete fh1FFJetPtRecEffRec; 

  delete fh2FFTrackPtRecEffRec;
  delete fh2FFZRecEffRec;
  delete fh2FFXiRecEffRec;
  
  delete fh2FFTrackPtSecRecNS;
  delete fh2FFZSecRecNS;
  delete fh2FFXiSecRecNS;
  
  delete fh2FFTrackPtSecRecS;
  delete fh2FFZSecRecS;
  delete fh2FFXiSecRecS;
  
  delete fh2FFTrackPtSecRecSsc;
  delete fh2FFZSecRecSsc;
  delete fh2FFXiSecRecSsc;
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteJetResponse(TString strInfile, TString strID, TString strOutfile, 
							   Bool_t updateOutfile, TString strOutDir )
{ 
  // read task ouput from MC and write single track eff - standard dir/list 
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  WriteJetResponse(strInfile,strdir,strlist,strOutfile,updateOutfile, strOutDir);
}

//_____________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteJetResponse(TString strInfile, TString strdir, TString strlist,
							   TString strOutfile, Bool_t updateOutfile, TString strOutDir)
{
  // read 3d THnSparse response matrices in pt,z,xi vs jet pt from file
  // project THnSparse + TH2 in jet pt slices 
  // write to strOutfile 

  THnSparse* hn3ResponseJetPt;
  THnSparse* hn3ResponseJetZ;
  THnSparse* hn3ResponseJetXi;

  // 2D response matrices 

  THnSparse* hnResponsePt[fNJetPtSlices];
  THnSparse* hnResponseZ[fNJetPtSlices];
  THnSparse* hnResponseXi[fNJetPtSlices];

  TH2F* h2ResponsePt[fNJetPtSlices];
  TH2F* h2ResponseZ[fNJetPtSlices];
  TH2F* h2ResponseXi[fNJetPtSlices];

  // 1D projections on gen pt / rec pt axes

  TH1F* h1FFPtRec[fNJetPtSlices]; 
  TH1F* h1FFZRec[fNJetPtSlices];
  TH1F* h1FFXiRec[fNJetPtSlices];

  TH1F* h1FFPtGen[fNJetPtSlices]; 
  TH1F* h1FFZGen[fNJetPtSlices];
  TH1F* h1FFXiGen[fNJetPtSlices];

  TH1F* h1RatioPt[fNJetPtSlices]; 
  TH1F* h1RatioZ[fNJetPtSlices];
  TH1F* h1RatioXi[fNJetPtSlices];



  TFile f(strInfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read response matrices from file %s ",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
    
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  
  
  if(list){
    hn3ResponseJetPt = (THnSparse*) list->FindObject("fhnResponseJetTrackPt");
    hn3ResponseJetZ  = (THnSparse*) list->FindObject("fhnResponseJetZ");
    hn3ResponseJetXi = (THnSparse*) list->FindObject("fhnResponseJetXi");
  }
  else{
    hn3ResponseJetPt = (THnSparse*) gDirectory->Get("fhnResponseJetTrackPt");
    hn3ResponseJetZ  = (THnSparse*) gDirectory->Get("fhnResponseJetZ");
    hn3ResponseJetXi = (THnSparse*) gDirectory->Get("fhnResponseJetXi");
  }

  
  if(!hn3ResponseJetPt){
    Printf("%s:%d -- error retrieving response matrix fhnResponseJetTrackPt",(char*)__FILE__,__LINE__);
    return;
  }

  if(!hn3ResponseJetZ){
    Printf("%s:%d -- error retrieving response matrix fhnResponseJetZ",(char*)__FILE__,__LINE__);
    return;
  }

  if(!hn3ResponseJetXi){
    Printf("%s:%d -- error retrieving response matrix fhnResponseJetXi",(char*)__FILE__,__LINE__);
    return;
  }

  f.Close();  

  // axes 

  Int_t axisJetPtTHn3 = -1;
  Int_t axisGenPtTHn3 = -1;
  Int_t axisRecPtTHn3 = -1;

  for(Int_t i=0; i<hn3ResponseJetPt->GetNdimensions(); i++){
    
    TString title = hn3ResponseJetPt->GetAxis(i)->GetTitle(); 

    if(title.Contains("jet p_{T}")) axisJetPtTHn3 = i; 
    if(title.Contains("gen p_{T}")) axisGenPtTHn3 = i; 
    if(title.Contains("rec p_{T}")) axisRecPtTHn3 = i; 
  }

  if(axisJetPtTHn3 == -1){
    Printf("%s:%d -- error axisJetPtTHn3",(char*)__FILE__,__LINE__);
    return;
  }

  if(axisGenPtTHn3 == -1){
    Printf("%s:%d -- error axisGenPtTHn3",(char*)__FILE__,__LINE__);
    return;
  }

  if(axisRecPtTHn3 == -1){
    Printf("%s:%d -- error axisRecPtTHn3",(char*)__FILE__,__LINE__);
    return;
  }

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);

    Int_t binLo = static_cast<Int_t>(hn3ResponseJetPt->GetAxis(axisJetPtTHn3)->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(hn3ResponseJetPt->GetAxis(axisJetPtTHn3)->FindBin(jetPtUpLim))-1;

    if(binUp > hn3ResponseJetPt->GetAxis(axisJetPtTHn3)->GetNbins()){
      Printf("%s:%d -- jet pt range %0.3f exceeds histo limits",(char*)__FILE__,__LINE__,jetPtUpLim); 
      return; 
    }
    
    TString strNameRespPt(Form("hnResponsePt_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameRespZ(Form("hnResponseZ_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameRespXi(Form("hnResponseXi_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));

    TString strNameTH2RespPt(Form("h2ResponsePt_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameTH2RespZ(Form("h2ResponseZ_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameTH2RespXi(Form("h2ResponseXi_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
         
    TString strNameRecPt(Form("h1FFTrackPtRecPrim_recPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameRecZ(Form("h1FFZRecPrim_recPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameRecXi(Form("h1FFXiRecPrim_recPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
 
    TString strNameGenPt(Form("h1FFTrackPtRecPrim_genPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameGenZ(Form("h1FFZRecPrim_genPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameGenXi(Form("h1FFXiRecPrim_genPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
 
    TString strNameRatioPt(Form("h1RatioTrackPtRecPrim_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameRatioZ(Form("h1RatioZRecPrim_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameRatioXi(Form("h1RatioXiRecPrim_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
 
    
    // 2D projections in jet pt range
  
    hn3ResponseJetPt->GetAxis(axisJetPtTHn3)->SetRange(binLo,binUp); 
    hn3ResponseJetZ->GetAxis(axisJetPtTHn3)->SetRange(binLo,binUp); 
    hn3ResponseJetXi->GetAxis(axisJetPtTHn3)->SetRange(binLo,binUp); 

    Int_t axesProj[2] = {axisRecPtTHn3,axisGenPtTHn3};  

    hnResponsePt[i] = hn3ResponseJetPt->Projection(2,axesProj);
    hnResponseZ[i]  = hn3ResponseJetZ->Projection(2,axesProj);
    hnResponseXi[i] = hn3ResponseJetXi->Projection(2,axesProj);
  
    hnResponsePt[i]->SetNameTitle(strNameRespPt,""); 
    hnResponseZ[i]->SetNameTitle(strNameRespZ,""); 
    hnResponseXi[i]->SetNameTitle(strNameRespXi,""); 

    h2ResponsePt[i] = (TH2F*) hnResponsePt[i]->Projection(1,0);// note convention: yDim,xDim
    h2ResponseZ[i]  = (TH2F*) hnResponseZ[i]->Projection(1,0); // note convention: yDim,xDim
    h2ResponseXi[i] = (TH2F*) hnResponseXi[i]->Projection(1,0);// note convention: yDim,xDim
 
    h2ResponsePt[i]->SetNameTitle(strNameTH2RespPt,""); 
    h2ResponseZ[i]->SetNameTitle(strNameTH2RespZ,""); 
    h2ResponseXi[i]->SetNameTitle(strNameTH2RespXi,""); 


    // 1D projections

    Int_t axisGenPtTHn2 = -1;
    Int_t axisRecPtTHn2 = -1;

    for(Int_t d=0; d<hnResponsePt[i]->GetNdimensions(); d++){
    
      TString title = hnResponsePt[i]->GetAxis(d)->GetTitle(); 

      if(title.Contains("gen p_{T}")) axisGenPtTHn2 = d; 
      if(title.Contains("rec p_{T}")) axisRecPtTHn2 = d; 
    }

    
    if(axisGenPtTHn2 == -1){
      Printf("%s:%d -- error axisGenPtTHn2",(char*)__FILE__,__LINE__);
      return;
    }
    
    if(axisRecPtTHn2 == -1){
      Printf("%s:%d -- error axisRecPtTHn2",(char*)__FILE__,__LINE__);
      return;
    }
    

    h1FFPtRec[i] = (TH1F*) hnResponsePt[i]->Projection(axisRecPtTHn2);// note convention: yDim,xDim
    h1FFZRec[i]  = (TH1F*) hnResponseZ[i]->Projection(axisRecPtTHn2);// note convention: yDim,xDim
    h1FFXiRec[i] = (TH1F*) hnResponseXi[i]->Projection(axisRecPtTHn2);// note convention: yDim,xDim
    
    h1FFPtRec[i]->SetNameTitle(strNameRecPt,""); 
    h1FFZRec[i]->SetNameTitle(strNameRecZ,""); 
    h1FFXiRec[i]->SetNameTitle(strNameRecXi,""); 

    h1FFPtGen[i] = (TH1F*) hnResponsePt[i]->Projection(axisGenPtTHn2);// note convention: yDim,xDim
    h1FFZGen[i]  = (TH1F*) hnResponseZ[i]->Projection(axisGenPtTHn2);// note convention: yDim,xDim
    h1FFXiGen[i] = (TH1F*) hnResponseXi[i]->Projection(axisGenPtTHn2);// note convention: yDim,xDim
    
    h1FFPtGen[i]->SetNameTitle(strNameGenPt,"");
    h1FFZGen[i]->SetNameTitle(strNameGenZ,"");
    h1FFXiGen[i]->SetNameTitle(strNameGenXi,"");

    // normalize 1D projections

    if(fNHistoBinsPt[i]) h1FFPtRec[i] = (TH1F*) h1FFPtRec[i]->Rebin(fNHistoBinsPt[i],"",fHistoBinsPt[i]->GetArray()); 
    if(fNHistoBinsZ[i])  h1FFZRec[i]  = (TH1F*) h1FFZRec[i]->Rebin(fNHistoBinsZ[i],"",fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) h1FFXiRec[i] = (TH1F*) h1FFXiRec[i]->Rebin(fNHistoBinsXi[i],"",fHistoBinsXi[i]->GetArray());
    
    if(fNHistoBinsPt[i]) h1FFPtGen[i] = (TH1F*) h1FFPtGen[i]->Rebin(fNHistoBinsPt[i],"",fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  h1FFZGen[i]  = (TH1F*) h1FFZGen[i]->Rebin(fNHistoBinsZ[i],"",fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) h1FFXiGen[i] = (TH1F*) h1FFXiGen[i]->Rebin(fNHistoBinsXi[i],"",fHistoBinsXi[i]->GetArray());
       
    NormalizeTH1(h1FFPtRec[i],fNJets->At(i)); 
    NormalizeTH1(h1FFZRec[i],fNJets->At(i));
    NormalizeTH1(h1FFXiRec[i],fNJets->At(i));

    NormalizeTH1(h1FFPtGen[i],fNJets->At(i)); 
    NormalizeTH1(h1FFZGen[i],fNJets->At(i));
    NormalizeTH1(h1FFXiGen[i],fNJets->At(i));

    // ratios 1D projections

    h1RatioPt[i] = (TH1F*) h1FFPtRec[i]->Clone(strNameRatioPt);
    h1RatioPt[i]->Reset();
    h1RatioPt[i]->Divide(h1FFPtRec[i],h1FFPtGen[i],1,1,"B");

    h1RatioZ[i] = (TH1F*) h1FFZRec[i]->Clone(strNameRatioZ);
    h1RatioZ[i]->Reset();
    h1RatioZ[i]->Divide(h1FFZRec[i],h1FFZGen[i],1,1,"B");
    
    h1RatioXi[i] = (TH1F*) h1FFXiRec[i]->Clone(strNameRatioXi);
    h1RatioXi[i]->Reset();
    h1RatioXi[i]->Divide(h1FFXiRec[i],h1FFXiGen[i],1,1,"B");
  }
   
  
  // write 

  TString outfileOption = "RECREATE";
  if(updateOutfile)  outfileOption = "UPDATE";
  
  TFile out(strOutfile,outfileOption);
  
  if(!out.IsOpen()){
    Printf("%s:%d -- error opening response matrix output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- write response matrices to file %s ",(char*)__FILE__,__LINE__,strOutfile.Data());

  if(strOutDir && strOutDir.Length()){  
    
    TDirectory* dir;
    if((dir = ((TDirectory*) gDirectory->Get(strOutDir)))) dir->cd(); 
    else{
      dir = out.mkdir(strOutDir);
      dir->cd(); 
    } 
  }


  for(Int_t i=0; i<fNJetPtSlices; i++){

    hnResponsePt[i]->Write();
    hnResponseXi[i]->Write();
    hnResponseZ[i]->Write();

    h2ResponsePt[i]->Write();
    h2ResponseXi[i]->Write();
    h2ResponseZ[i]->Write();

    h1FFPtRec[i]->Write(); 
    h1FFZRec[i]->Write();
    h1FFXiRec[i]->Write();
    
    h1FFPtGen[i]->Write(); 
    h1FFZGen[i]->Write();
    h1FFXiGen[i]->Write();
    
    h1RatioPt[i]->Write(); 
    h1RatioZ[i]->Write();
    h1RatioXi[i]->Write();

  }

  out.Close();  
}

//______________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadResponse(TString strfile, TString strdir, TString strlist)
{
  // read response matrices from file
  // argument strlist optional - read from directory strdir if not specified
  // note: THnSparse are not rebinned 

  THnSparse* hRespPt[fNJetPtSlices]; 
  THnSparse* hRespZ[fNJetPtSlices];
  THnSparse* hRespXi[fNJetPtSlices];
  
  for(Int_t i=0; i<fNJetPtSlices; i++) hRespPt[i] = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hRespZ[i]  = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hRespXi[i] = 0;

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read response matrices from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
   
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Int_t jetPtLoLim = static_cast<Int_t> (fJetPtSlices->At(i));
    Int_t jetPtUpLim = static_cast<Int_t> (fJetPtSlices->At(i+1));
    
    TString strNameRespPt(Form("hnResponsePt_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameRespZ(Form("hnResponseZ_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameRespXi(Form("hnResponseXi_%02d_%02d",jetPtLoLim,jetPtUpLim));
        
    if(list){
      hRespPt[i] = (THnSparse*) list->FindObject(strNameRespPt); 
      hRespZ[i]  = (THnSparse*) list->FindObject(strNameRespZ); 
      hRespXi[i] = (THnSparse*) list->FindObject(strNameRespXi); 
    }
    else{
      hRespPt[i] = (THnSparse*) gDirectory->Get(strNameRespPt); 
      hRespZ[i]  = (THnSparse*) gDirectory->Get(strNameRespZ); 
      hRespXi[i] = (THnSparse*) gDirectory->Get(strNameRespXi); 
    }
    
    if(!hRespPt[i]){
      Printf("%s:%d -- error retrieving response %s", (char*)__FILE__,__LINE__,strNameRespPt.Data());
    }
  
    if(!hRespZ[i]){
      Printf("%s:%d -- error retrieving response %s", (char*)__FILE__,__LINE__,strNameRespZ.Data());
    }    

    if(!hRespXi[i]){
      Printf("%s:%d -- error retrieving response %s", (char*)__FILE__,__LINE__,strNameRespXi.Data());
    }
   
    //    if(0){ // can't rebin THnSparse ...
    //       if(fNHistoBinsPt[i]) hRespPt[i]->SetBinEdges(0,fHistoBinsPt[i]->GetArray());
    //       if(fNHistoBinsPt[i]) hRespPt[i]->SetBinEdges(1,fHistoBinsPt[i]->GetArray());
    
    //       if(fNHistoBinsZ[i])  hRespZ[i]->SetBinEdges(0,fHistoBinsZ[i]->GetArray());
    //       if(fNHistoBinsZ[i])  hRespZ[i]->SetBinEdges(1,fHistoBinsZ[i]->GetArray());
    
    //       if(fNHistoBinsXi[i]) hRespXi[i]->SetBinEdges(0,fHistoBinsXi[i]->GetArray());
    //       if(fNHistoBinsXi[i]) hRespXi[i]->SetBinEdges(1,fHistoBinsXi[i]->GetArray());
    //     }
   
 
  } // jet slices loop

  f.Close();  // THnSparse pointers still valid even if file closed

//   for(Int_t i=0; i<fNJetPtSlices; i++){  // no copy c'tor ...
//     if(hRespPt[i]) new(fhnResponsePt[i]) THnSparseF(*hRespPt[i]);
//     if(hRespZ[i])  new(fhnResponseZ[i])  THnSparseF(*hRespZ[i]);
//     if(hRespXi[i]) new(fhnResponseXi[i]) THnSparseF(*hRespXi[i]);
//   }

  for(Int_t i=0; i<fNJetPtSlices; i++){ 
    fhnResponsePt[i] = hRespPt[i];
    fhnResponseZ[i]  = hRespZ[i];
    fhnResponseXi[i] = hRespXi[i];
  }
}

//______________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadPriors(TString strfile,const Int_t type)
{
  // read priors from file: rec primaries, gen pt dist

  if(fDebug>0) Printf("%s:%d -- read priors from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
  
  // temporary histos to store pointers from file
  TH1F* hist[fNJetPtSlices]; 

  for(Int_t i=0; i<fNJetPtSlices; i++) hist[i] = 0;

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Int_t jetPtLoLim = static_cast<Int_t> (fJetPtSlices->At(i));
    Int_t jetPtUpLim = static_cast<Int_t> (fJetPtSlices->At(i+1));
    
    TString strName;
    
    if(type == kFlagPt) strName.Form("h1FFTrackPtRecPrim_genPt_%02d_%02d",jetPtLoLim,jetPtUpLim); 
    if(type == kFlagZ)  strName.Form("h1FFZRecPrim_genPt_%02d_%02d",jetPtLoLim,jetPtUpLim);       
    if(type == kFlagXi) strName.Form("h1FFXiRecPrim_genPt_%02d_%02d",jetPtLoLim,jetPtUpLim);      
    
    hist[i] = (TH1F*) gDirectory->Get(strName); 
    
    if(!hist[i]){
      Printf("%s:%d -- error retrieving prior %s", (char*)__FILE__,__LINE__,strName.Data());
    }
    

    //if(fNHistoBinsPt[i]) hist[i] = (TH1F*) hist[i]->Rebin(fNHistoBinsPt[i],hist[i]->GetName()+"_rebin",fHistoBinsPt[i]->GetArray());

    if(hist[i]) hist[i]->SetDirectory(0); 
    
  } // jet slices loop
  
  f.Close();

  
  for(Int_t i=0; i<fNJetPtSlices; i++){ // 2nd loop: need to close input file before placing histos
    if(hist[i] && type == kFlagPt) new(fh1FFTrackPtPrior[i]) TH1F(*hist[i]); 
    if(hist[i] && type == kFlagZ)  new(fh1FFZPrior[i])       TH1F(*hist[i]); 
    if(hist[i] && type == kFlagXi) new(fh1FFXiPrior[i])      TH1F(*hist[i]); 
  }
}


//_____________________________________________________
// void AliFragmentationFunctionCorrections::RatioRecGen()
// {
//   // create ratio reconstructed over generated FF
//   // use current highest corrLevel 
  
//   Printf("%s:%d -- build ratio rec.gen, corrLevel %d",(char*)__FILE__,__LINE__,fNCorrectionLevels-1);

//   for(Int_t i=0; i<fNJetPtSlices; i++){
 
//     TH1F* histPtRec =  fCorrFF[fNCorrectionLevels-1]->GetTrackPt(i); // levels -1: latest corr level
//     TH1F* histZRec  =  fCorrFF[fNCorrectionLevels-1]->GetZ(i);       // levels -1: latest corr level
//     TH1F* histXiRec =  fCorrFF[fNCorrectionLevels-1]->GetXi(i);      // levels -1: latest corr level
    
//     TH1F* histPtGen = fh1FFTrackPtGenPrim[i]; 
//     TH1F* histZGen  = fh1FFZGenPrim[i];
//     TH1F* histXiGen = fh1FFXiGenPrim[i];

//     TString histNamePt = histPtRec->GetName();
//     TString histNameZ  = histZRec->GetName();
//     TString histNameXi = histXiRec->GetName();

//     histNamePt.ReplaceAll("fh1FF","fh1FFRatioRecGen");
//     histNameZ.ReplaceAll("fh1FF","fh1FFRatioRecGen");
//     histNameXi.ReplaceAll("fh1FF","fh1FFRatioRecGen");

//     // ratio 
//     TH1F* hRatioRecGenPt = (TH1F*) histPtRec->Clone(histNamePt);
//     hRatioRecGenPt->Reset();
//     hRatioRecGenPt->Divide(histPtRec,histPtGen,1,1,"B");

//     TH1F* hRatioRecGenZ = (TH1F*) histZRec->Clone(histNameZ);
//     hRatioRecGenZ->Reset();
//     hRatioRecGenZ->Divide(histZRec,histZGen,1,1,"B");
    
//     TH1F* hRatioRecGenXi = (TH1F*) histXiRec->Clone(histNameXi);
//     hRatioRecGenXi->Reset();
//     hRatioRecGenXi->Divide(histXiRec,histXiGen,1,1,"B");

//     new(fh1FFRatioRecGenPt[i]) TH1F(*hRatioRecGenPt);
//     new(fh1FFRatioRecGenZ[i])  TH1F(*hRatioRecGenZ);
//     new(fh1FFRatioRecGenXi[i]) TH1F(*hRatioRecGenXi);
//   }
// }

// //___________________________________________________________
// void AliFragmentationFunctionCorrections::RatioRecPrimaries()
// {
//   // create ratio reconstructed tracks over reconstructed primaries 
//   // use raw FF (corrLevel 0)
  
//   Printf("%s:%d -- build ratio rec tracks /rec primaries",(char*)__FILE__,__LINE__);

//   for(Int_t i=0; i<fNJetPtSlices; i++){
 
//     const Int_t corrLevel = 0; 

//     TH1F* histPtRec =  fCorrFF[corrLevel]->GetTrackPt(i); // levels -1: latest corr level
//     TH1F* histZRec  =  fCorrFF[corrLevel]->GetZ(i);       // levels -1: latest corr level
//     TH1F* histXiRec =  fCorrFF[corrLevel]->GetXi(i);      // levels -1: latest corr level
    
//     TH1F* histPtRecPrim = fh1FFTrackPtRecPrim[i]; 
//     TH1F* histZRecPrim  = fh1FFZRecPrim[i];
//     TH1F* histXiRecPrim = fh1FFXiRecPrim[i];

//     TString histNamePt = histPtRec->GetName();
//     TString histNameZ  = histZRec->GetName();
//     TString histNameXi = histXiRec->GetName();

//     histNamePt.ReplaceAll("fh1FF","fh1FFRatioRecPrim");
//     histNameZ.ReplaceAll("fh1FF","fh1FFRatioRecPrim");
//     histNameXi.ReplaceAll("fh1FF","fh1FFRatioRecPrim");

//     // ratio 
//     TH1F* hRatioRecPrimPt = (TH1F*) histPtRec->Clone(histNamePt);
//     hRatioRecPrimPt->Reset();
//     hRatioRecPrimPt->Divide(histPtRec,histPtRecPrim,1,1,"B");

//     TH1F* hRatioRecPrimZ = (TH1F*) histZRec->Clone(histNameZ);
//     hRatioRecPrimZ->Reset();
//     hRatioRecPrimZ->Divide(histZRec,histZRecPrim,1,1,"B");
    
//     TH1F* hRatioRecPrimXi = (TH1F*) histXiRec->Clone(histNameXi);
//     hRatioRecPrimXi->Reset();
//     hRatioRecPrimXi->Divide(histXiRec,histXiRecPrim,1,1,"B");


//     new(fh1FFRatioRecPrimPt[i]) TH1F(*hRatioRecPrimPt);
//     new(fh1FFRatioRecPrimZ[i])  TH1F(*hRatioRecPrimZ);
//     new(fh1FFRatioRecPrimXi[i]) TH1F(*hRatioRecPrimXi);
//   }
// }

// __________________________________________________________________________________
void AliFragmentationFunctionCorrections::ProjectJetResponseMatrices(TString strOutfile)
{

  // project response matrices on both axes: 
  // FF for rec primaries, in terms of generated and reconstructed momentum
  // write FF and ratios to outFile

  Printf("%s:%d -- project response matrices, write to %s",(char*)__FILE__,__LINE__,strOutfile.Data());

  TH1F* hFFPtRec[fNJetPtSlices]; 
  TH1F* hFFZRec[fNJetPtSlices];
  TH1F* hFFXiRec[fNJetPtSlices];

  TH1F* hFFPtGen[fNJetPtSlices]; 
  TH1F* hFFZGen[fNJetPtSlices];
  TH1F* hFFXiGen[fNJetPtSlices];

  TH1F* hRatioPt[fNJetPtSlices]; 
  TH1F* hRatioZ[fNJetPtSlices];
  TH1F* hRatioXi[fNJetPtSlices];


  Int_t axisGenPt = 1;
  Int_t axisRecPt = 0;

  for(Int_t i=0; i<fNJetPtSlices; i++){ 

    Int_t jetPtLoLim = static_cast<Int_t> (fJetPtSlices->At(i));
    Int_t jetPtUpLim = static_cast<Int_t> (fJetPtSlices->At(i+1));
 
    TString strNameRecPt(Form("h1FFTrackPtRecPrim_recPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameRecZ(Form("h1FFZRecPrim_recPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameRecXi(Form("h1FFXiRecPrim_recPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
 
    TString strNameGenPt(Form("h1FFTrackPtRecPrim_genPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameGenZ(Form("h1FFZRecPrim_genPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameGenXi(Form("h1FFXiRecPrim_genPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
 
    TString strNameRatioPt(Form("h1RatioTrackPtRecPrim_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameRatioZ(Form("h1RatioZRecPrim_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameRatioXi(Form("h1RatioXiRecPrim_%02d_%02d",jetPtLoLim,jetPtUpLim));
 

    hFFPtRec[i] = (TH1F*) fhnResponsePt[i]->Projection(axisRecPt);// note convention: yDim,xDim
    hFFZRec[i]  = (TH1F*) fhnResponseZ[i]->Projection(axisRecPt);// note convention: yDim,xDim
    hFFXiRec[i] = (TH1F*) fhnResponseXi[i]->Projection(axisRecPt);// note convention: yDim,xDim
    
    hFFPtRec[i]->SetNameTitle(strNameRecPt,""); 
    hFFZRec[i]->SetNameTitle(strNameRecZ,""); 
    hFFXiRec[i]->SetNameTitle(strNameRecXi,""); 


    hFFPtGen[i] = (TH1F*) fhnResponsePt[i]->Projection(axisGenPt);// note convention: yDim,xDim
    hFFZGen[i]  = (TH1F*) fhnResponseZ[i]->Projection(axisGenPt);// note convention: yDim,xDim
    hFFXiGen[i] = (TH1F*) fhnResponseXi[i]->Projection(axisGenPt);// note convention: yDim,xDim
    
    hFFPtGen[i]->SetNameTitle(strNameGenPt,""); 
    hFFZGen[i]->SetNameTitle(strNameGenZ,""); 
    hFFXiGen[i]->SetNameTitle(strNameGenXi,""); 
   

    if(fNHistoBinsPt[i]) hFFPtRec[i] = (TH1F*) hFFPtRec[i]->Rebin(fNHistoBinsPt[i],"",fHistoBinsPt[i]->GetArray()); 
    if(fNHistoBinsZ[i])  hFFZRec[i]   = (TH1F*) hFFZRec[i]->Rebin(fNHistoBinsZ[i],"",fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hFFXiRec[i]  = (TH1F*) hFFXiRec[i]->Rebin(fNHistoBinsXi[i],"",fHistoBinsXi[i]->GetArray());
    
    if(fNHistoBinsPt[i]) hFFPtGen[i] = (TH1F*) hFFPtGen[i]->Rebin(fNHistoBinsPt[i],"",fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hFFZGen[i]  = (TH1F*) hFFZGen[i]->Rebin(fNHistoBinsZ[i],"",fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hFFXiGen[i] = (TH1F*) hFFXiGen[i]->Rebin(fNHistoBinsXi[i],"",fHistoBinsXi[i]->GetArray());
    
    NormalizeTH1(hFFPtGen[i],fNJets->At(i)); 
    NormalizeTH1(hFFZGen[i],fNJets->At(i));
    NormalizeTH1(hFFXiGen[i],fNJets->At(i));
   
    NormalizeTH1(hFFPtRec[i],fNJets->At(i)); 
    NormalizeTH1(hFFZRec[i],fNJets->At(i));
    NormalizeTH1(hFFXiRec[i],fNJets->At(i));


    hRatioPt[i] = (TH1F*) hFFPtRec[i]->Clone(strNameRatioPt);
    hRatioPt[i]->Reset();
    hRatioPt[i]->Divide(hFFPtRec[i],hFFPtGen[i],1,1,"B");

    hRatioZ[i] = (TH1F*) hFFZRec[i]->Clone(strNameRatioZ);
    hRatioZ[i]->Reset();
    hRatioZ[i]->Divide(hFFZRec[i],hFFZGen[i],1,1,"B");
    
    hRatioXi[i] = (TH1F*) hFFXiRec[i]->Clone(strNameRatioXi);
    hRatioXi[i]->Reset();
    hRatioXi[i]->Divide(hFFXiRec[i],hFFXiGen[i],1,1,"B");
  }
    


  // write 
  
  TFile out(strOutfile,"RECREATE");
  
  if(!out.IsOpen()){
    Printf("%s:%d -- error opening response matrix output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    hFFPtRec[i]->Write(); 
    hFFZRec[i]->Write();
    hFFXiRec[i]->Write();
    
    hFFPtGen[i]->Write(); 
    hFFZGen[i]->Write();
    hFFXiGen[i]->Write();
    
    hRatioPt[i]->Write(); 
    hRatioZ[i]->Write();
    hRatioXi[i]->Write();
  }

  out.Close();  
}

// ____________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ProjectSingleResponseMatrix(TString strOutfile, Bool_t updateOutfile, TString strOutDir)
{ 
  // project response matrix on both axes: 
  // pt spec for rec primaries, in terms of generated and reconstructed momentum
  // write spec and ratios to outFile
  
  Printf("%s:%d -- project single pt response matrix, write to %s",(char*)__FILE__,__LINE__,strOutfile.Data());
  
  TH1F* hSpecPtRec;  
  TH1F* hSpecPtGen; 
  TH1F* hRatioPt; 
  
  Int_t axisGenPt = 1;
  Int_t axisRecPt = 0;
  
  TString strNameRecPt   = "h1SpecTrackPtRecPrim_recPt";
  TString strNameGenPt   = "h1SpecTrackPtRecPrim_genPt";
  TString strNameRatioPt = "h1RatioTrackPtRecPrim";
    
  hSpecPtRec = (TH1F*) fhnResponseSinglePt->Projection(axisRecPt);// note convention: yDim,xDim
  hSpecPtRec->SetNameTitle(strNameRecPt,""); 
  
  hSpecPtGen = (TH1F*) fhnResponseSinglePt->Projection(axisGenPt);// note convention: yDim,xDim
  hSpecPtGen->SetNameTitle(strNameGenPt,""); 
  
  hRatioPt = (TH1F*) hSpecPtRec->Clone(strNameRatioPt);
  hRatioPt->Reset();
  hRatioPt->Divide(hSpecPtRec,hSpecPtGen,1,1,"B");
  
  TString outfileOption = "RECREATE";
  if(updateOutfile) outfileOption = "UPDATE";

  TFile out(strOutfile,outfileOption);
  
  if(!out.IsOpen()){
    Printf("%s:%d -- error opening reponse matrix projections output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }
      

  if(strOutDir && strOutDir.Length()){  
    
    TDirectory* dir;
    if((dir = ((TDirectory*) gDirectory->Get(strOutDir)))) dir->cd(); 
    else{
      dir = out.mkdir(strOutDir);
      dir->cd(); 
    } 
  }
  
  hSpecPtRec->Write(); 
  hSpecPtGen->Write(); 
  hRatioPt->Write(); 
  
  out.Close();  
}


//__________________________________________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::RebinHisto(const Int_t jetPtSlice, const Int_t nBinsLimits, Double_t* binsLimits, Double_t* binsWidth, const Int_t type)
{ 
  // rebin histo, rescale bins according to new width
  // only correct for input histos with equal bin size

  // args: jetPtSlice, type, use current corr level
 
  // function takes array of limits and widths (e.g. 1 GeV bins up to 10 GeV, 2 GeV width up to 20 GeV, ...)  
  // array size of binsLimits: nBinsLimits 
  // array size of binsWidth: nBinsLimits-1 
  // binsLimits have to be in increasing order
  // if binning undefined for any slice, original binning will be kept

  if(!fNJetPtSlices){
    Printf("%s:%d -- jetPtSlices not defined", (char*)__FILE__,__LINE__);
    return;
  }

  if(jetPtSlice>=fNJetPtSlices){
    Printf("%s:%d -- jetPtSlice %d exceeds max",(char*)__FILE__,__LINE__,jetPtSlice);
    return;
  }


  Double_t binLimitMin = binsLimits[0];
  Double_t binLimitMax = binsLimits[nBinsLimits-1];

  Double_t binLimit = binLimitMin; // start value 
  
  Int_t sizeUpperLim = 1000; //static_cast<Int_t>(binLimitMax/binsWidth[0])+1; -  only if first bin has smallest width, but not case for dN/dxi ...
  TArrayD binsArray(sizeUpperLim);
  Int_t nBins = 0; 
  binsArray.SetAt(binLimitMin,nBins++);

  while(binLimit<binLimitMax && nBins<sizeUpperLim){

    Int_t currentSlice = -1;
    for(Int_t i=0; i<nBinsLimits; i++){
      if(binLimit >= 0.999*binsLimits[i]) currentSlice = i; // 0.999 numerical saftey factor 
    }
    
    Double_t currentBinWidth = binsWidth[currentSlice];
    binLimit += currentBinWidth;

    binsArray.SetAt(binLimit,nBins++);
  }
  

  TH1F* hist = 0;
  if(type == kFlagPt)            hist = fCorrFF[fNCorrectionLevels-1]->GetTrackPt(jetPtSlice); 
  else if(type == kFlagZ)        hist = fCorrFF[fNCorrectionLevels-1]->GetZ(jetPtSlice);       
  else if(type == kFlagXi)       hist = fCorrFF[fNCorrectionLevels-1]->GetXi(jetPtSlice);      
  else if(type == kFlagSinglePt) hist = fCorrSinglePt[fNCorrectionLevelsSinglePt-1]->GetTrackPt(0); 
  else{ 
    Printf("%s%d unknown type",(char*)__FILE__,__LINE__);
    return;
  }
  
  Double_t binWidthNoRebin = hist->GetBinWidth(1);

  Double_t* bins = binsArray.GetArray();

  hist = (TH1F*) hist->Rebin(nBins-1,"",bins);
 
  for(Int_t bin=0; bin <= hist->GetNbinsX(); bin++){

    Double_t binWidthRebin = hist->GetBinWidth(bin);
    Double_t scaleF = binWidthNoRebin / binWidthRebin;

    Double_t binCont  = hist->GetBinContent(bin);
    Double_t binErr   = hist->GetBinError(bin);
    
    binCont *= scaleF;
    binErr  *= scaleF;

    hist->SetBinContent(bin,binCont);
    hist->SetBinError(bin,binErr);
  }



  TH1F* temp = new TH1F(*hist);

  if(type == kFlagPt)       fCorrFF[fNCorrectionLevels-1]->ReplaceCorrHistos(jetPtSlice,temp,0,0);
  if(type == kFlagZ)        fCorrFF[fNCorrectionLevels-1]->ReplaceCorrHistos(jetPtSlice,0,temp,0);
  if(type == kFlagXi)       fCorrFF[fNCorrectionLevels-1]->ReplaceCorrHistos(jetPtSlice,0,0,temp);
  if(type == kFlagSinglePt) fCorrSinglePt[fNCorrectionLevelsSinglePt-1]->ReplaceCorrHistos(0,temp,0,0);


  delete temp;
}
//__________________________________________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteJetSpecResponse(TString strInfile, TString strdir, TString strlist/*, TString strOutfile*/)
{ 
 
 if(fDebug>0) Printf("%s:%d -- read jet spectrum response matrix from file %s ",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){ 
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  
  if(list == 0)return; // catch strlist.Lenght() == 0;
   
  THnSparse* hn6ResponseJetPt  = (THnSparse*) list->FindObject("fhnCorrelation");
 
  Int_t axis6RecJetPt = 0;
  Int_t axis6GenJetPt = 3;

  hn6ResponseJetPt->GetAxis(axis6RecJetPt)->SetTitle("rec jet p_{T} (GeV/c)"); 
  hn6ResponseJetPt->GetAxis(axis6GenJetPt)->SetTitle("gen jet p_{T} (GeV/c)"); 
  
  Int_t nBinsRecPt    = hn6ResponseJetPt->GetAxis(axis6RecJetPt)->GetNbins(); 
  Double_t loLimRecPt = hn6ResponseJetPt->GetAxis(axis6RecJetPt)->GetBinLowEdge(1);
  Double_t upLimRecPt = hn6ResponseJetPt->GetAxis(axis6RecJetPt)->GetBinUpEdge(nBinsRecPt);

  Int_t nBinsGenPt    = hn6ResponseJetPt->GetAxis(axis6GenJetPt)->GetNbins(); 
  Double_t loLimGenPt = hn6ResponseJetPt->GetAxis(axis6GenJetPt)->GetBinLowEdge(1);
  Double_t upLimGenPt = hn6ResponseJetPt->GetAxis(axis6GenJetPt)->GetBinUpEdge(nBinsGenPt);

  Int_t nBinsTrackPt = 200;
  Int_t loLimTrackPt = 0;
  Int_t upLimTrackPt = 200;
  

  Int_t    nBinsResponse[4]  = {nBinsRecPt,nBinsTrackPt,nBinsGenPt,nBinsTrackPt};
  Double_t binMinResponse[4] = {loLimRecPt,loLimTrackPt,loLimGenPt,loLimTrackPt};
  Double_t binMaxResponse[4] = {upLimRecPt,upLimTrackPt,upLimGenPt,upLimTrackPt};

  const char* labelsResponseSinglePt[4] = {"rec jet p_{T} (GeV/c)", "rec track p_{T} (GeV/c)", "gen jet p_{T} (GeV/c)", "gen track p_{T} (GeV/c)"};
  
  THnSparseD* hn4ResponseTrackPtJetPt =  new THnSparseD("hn4ResponseTrackPtJetPt","",4,nBinsResponse,binMinResponse,binMaxResponse);
  
  for(Int_t i=0; i<4; i++){
    hn4ResponseTrackPtJetPt->GetAxis(i)->SetTitle(labelsResponseSinglePt[i]); 
  }
  

  // fill 
  

}

//_____________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadSingleTrackEfficiency(TString strfile, TString strdir, TString strlist, TString strname)
{

  ReadSingleTrackCorrection(strfile,strdir,strlist,strname,kFlagEfficiency);

}

//_____________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadSingleTrackResponse(TString strfile, TString strdir, TString strlist, TString strname)
{

  ReadSingleTrackCorrection(strfile,strdir,strlist,strname,kFlagResponse);

}

//_____________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadSingleTrackSecCorr(TString strfile, TString strdir, TString strlist, TString strname)
{

  ReadSingleTrackCorrection(strfile,strdir,strlist,strname,kFlagSecondaries);

}

//______________________________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadSingleTrackCorrection(TString strfile, TString strdir, TString strlist, TString strname, const Int_t type)
{
  // read single track correction (pt) from file
  // type: efficiency / response / secondaries correction 

  if(!((type == kFlagEfficiency)  || (type == kFlagResponse) || (type == kFlagSecondaries))){
    Printf("%s:%d -- no such correction ",(char*)__FILE__,__LINE__);
    return;
  }

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0 && type==kFlagEfficiency)  Printf("%s:%d -- read single track corr from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
  if(fDebug>0 && type==kFlagResponse)    Printf("%s:%d -- read single track response from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
  if(fDebug>0 && type==kFlagSecondaries) Printf("%s:%d -- read single track secondaries corr from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
    
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  
  
  TH1F*      h1CorrHist = 0; // common TObject pointer not possible, need SetDirectory() later
  THnSparse* hnCorrHist = 0;
  
  if(type == kFlagEfficiency || type == kFlagSecondaries){

    if(list) h1CorrHist = (TH1F*) list->FindObject(strname); 
    else     h1CorrHist = (TH1F*) gDirectory->Get(strname);

    if(!h1CorrHist){
      Printf("%s:%d -- error retrieving histo %s", (char*)__FILE__,__LINE__,strname.Data());
      return;
    }

  }
  else if(type == kFlagResponse){
    
    if(list) hnCorrHist = (THnSparse*) list->FindObject(strname);
    else     hnCorrHist = (THnSparse*) gDirectory->Get(strname);

    if(!hnCorrHist){
      Printf("%s:%d -- error retrieving histo %s", (char*)__FILE__,__LINE__,strname.Data());
      return;
    }
 
  }

  if(h1CorrHist) h1CorrHist->SetDirectory(0);
  //if(hnCorrHist) hnCorrHist->SetDirectory(0);
  
  f.Close();  

  if(type == kFlagEfficiency)       fh1EffSinglePt      = h1CorrHist; 
  else if(type == kFlagResponse)    fhnResponseSinglePt = hnCorrHist;
  else if(type == kFlagSecondaries) fh1SecCorrSinglePt  = h1CorrHist; 
 
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawPtSpec(TString strInfile, TString strID)
{ 
  // read track pt spec from task ouput - standard dir/list 
     
  TString strdir  = "PWGJE_FragmentationFunction_" + strID;
  TString strlist = "fracfunc_" + strID;
    
  ReadRawPtSpec(strInfile,strdir,strlist);
}

//_______________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawPtSpec(TString strfile, TString strdir, TString strlist)
{
  // get raw pt spectra from file 

  // book histos
  fNCorrectionLevelsSinglePt = 0; 
  fCorrSinglePt = new AliFragFuncCorrHistos*[fgMaxNCorrectionLevels];
  AddCorrectionLevelSinglePt(); // first 'correction' level = raw spectrum

  // get raw pt spec from input file, normalize

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read raw spectra from file %s ",(char*)__FILE__,__LINE__,strfile.Data());

  gDirectory->cd(strdir);

  TList* list = 0;
  
  if(!(list = (TList*) gDirectory->Get(strlist))){ 
    Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
    return;
  }

  TString hnameTrackPt("fh1TrackQAPtRecCuts");
  TString hnameEvtSel("fh1EvtSelection");

  TH1F* fh1TrackPt = (TH1F*) list->FindObject(hnameTrackPt);
  TH1F* fh1EvtSel  = (TH1F*) list->FindObject(hnameEvtSel);  

  if(!fh1TrackPt){ Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameTrackPt.Data()); return; }
  if(!fh1EvtSel) { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameEvtSel.Data()); return; }
   
  Float_t nEvents = fh1EvtSel->GetBinContent(fh1EvtSel->FindBin(0));
  
  fh1TrackPt->SetDirectory(0);
  
  f.Close();  
  
  // rebin + normalize
  if(fNHistoBinsSinglePt) fh1TrackPt = (TH1F*) fh1TrackPt->Rebin(fNHistoBinsSinglePt,"",fHistoBinsSinglePt->GetArray());

  NormalizeTH1(fh1TrackPt,nEvents);

  // raw FF = corr level 0
  fCorrSinglePt[0]->AddCorrHistos(0,fh1TrackPt);
}


//_______________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadRawPtSpecQATask(TString strfile, TString strdir, TString strlist)
{
  // get raw pt spectra from file 
  // for output from Martas QA task


  // book histos
  fNCorrectionLevelsSinglePt = 0; 
  fCorrSinglePt = new AliFragFuncCorrHistos*[fgMaxNCorrectionLevels];
  AddCorrectionLevelSinglePt(); // first 'correction' level = raw spectrum

  // get raw pt spec from input file, normalize

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read raw pt spec from QA task output file %s ",(char*)__FILE__,__LINE__,strfile.Data());

  gDirectory->cd(strdir);

  TList* list = 0;
  
  if(!(list = (TList*) gDirectory->Get(strlist))){ 
    Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
    return;
  }

  TString hnameTrackPt("fPtSel");
  TString hnameEvtSel("fNEventAll");

  TH1F* fh1TrackPt = (TH1F*) list->FindObject(hnameTrackPt);
  TH1F* fh1EvtSel  = (TH1F*) list->FindObject(hnameEvtSel);  

  if(!fh1TrackPt){ Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameTrackPt.Data()); return; }
  if(!fh1EvtSel) { Printf("%s:%d -- histo %s not found",(char*)__FILE__,__LINE__,hnameEvtSel.Data()); return; }
 

  // evts after physics selection
  Float_t nEvents = fh1EvtSel->GetEntries();
 
  fh1TrackPt->SetDirectory(0);
  
  f.Close();  

  
  NormalizeTH1(fh1TrackPt,nEvents);

  // raw FF = corr level 0
  fCorrSinglePt[0]->AddCorrHistos(0,fh1TrackPt);
}

// ________________________________________________________
void AliFragmentationFunctionCorrections::EffCorrSinglePt()
{
  // apply efficiency correction to inclusive track pt spec

  AddCorrectionLevelSinglePt("eff");

  TH1F* histPt = fCorrSinglePt[fNCorrectionLevelsSinglePt-2]->GetTrackPt(0);
 
  if(histPt->GetNbinsX() != fh1EffSinglePt->GetNbinsX()) Printf("%s:%d: inconsistency pt spec and eff corr bins - rebin effCorr ...", (char*)__FILE__,__LINE__);

  TString histNamePt = histPt->GetName();
  TH1F* hTrackPtEffCorr = (TH1F*) histPt->Clone(histNamePt);
  hTrackPtEffCorr->Divide(histPt,fh1EffSinglePt,1,1,"");

  fCorrSinglePt[fNCorrectionLevelsSinglePt-1]->AddCorrHistos(0,hTrackPtEffCorr);
}

//___________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::UnfoldSinglePt(const Int_t nIter, const Bool_t useCorrelatedErrors)
{
  // unfolde inclusive dN/dpt spectra
  
  AddCorrectionLevelSinglePt("unfold");
   
  TH1F* hist = fCorrSinglePt[fNCorrectionLevelsSinglePt-2]->GetTrackPt(0); // level -2: before unfolding, level -1: unfolded
  THnSparse* hnResponse = fhnResponseSinglePt;
      
  TString histNameTHn = hist->GetName();
  if(histNameTHn.Contains("TH1")) histNameTHn.ReplaceAll("TH1","THn");
  if(histNameTHn.Contains("fPt")) histNameTHn.ReplaceAll("fPt","fhnPt");
  

  TString histNameBackFolded = hist->GetName();
  histNameBackFolded.Append("_backfold");
  
  TString histNameRatioFolded = hist->GetName();
  if(histNameRatioFolded.Contains("fh1")) histNameRatioFolded.ReplaceAll("fh1","hRatio");
  if(histNameRatioFolded.Contains("fPt")) histNameRatioFolded.ReplaceAll("fPt","hRatioPt");
  histNameRatioFolded.Append("_unfold");
  
  TString histNameRatioBackFolded = hist->GetName();
  if(histNameRatioBackFolded.Contains("fh1")) histNameRatioBackFolded.ReplaceAll("fh1","hRatio");
  if(histNameRatioBackFolded.Contains("fPt")) histNameRatioBackFolded.ReplaceAll("fPt","hRatioPt");
  histNameRatioBackFolded.Append("_backfold");

  THnSparse* hnHist           = TH1toSparse(hist,histNameTHn,hist->GetTitle());
  THnSparse* hnFlatEfficiency = TH1toSparse(hist,"fhnEfficiency","eff",kTRUE); // could optionally also use real eff 
  
  TH1F* hUnfolded 
    = Unfold(hnHist,hnResponse,hnFlatEfficiency,nIter,useCorrelatedErrors);  
  
  //TH1F* hUnfolded = (TH1F*) hnUnfolded->Projection(0); 
  //hUnfolded->SetNameTitle(hist->GetName(),hist->GetTitle());
  
  fCorrSinglePt[fNCorrectionLevelsSinglePt-1]->AddCorrHistos(0,hUnfolded);
  
  // backfolding: apply response matrix to unfolded spectrum
  TH1F* hBackFolded = ApplyResponse(hUnfolded,hnResponse); 
  hBackFolded->SetNameTitle(histNameBackFolded,hist->GetTitle());
  
  fh1SingleTrackPtBackFolded = hBackFolded;
  
  
  // ratio unfolded to original histo 
  TH1F* hRatioUnfolded = (TH1F*) hUnfolded->Clone(histNameRatioFolded);
  hRatioUnfolded->Reset();
  hRatioUnfolded->Divide(hUnfolded,hist,1,1,"B");
  
  fh1RatioSingleTrackPtFolded = hRatioUnfolded;
    
  
  // ratio backfolded to original histo
  TH1F* hRatioBackFolded = (TH1F*) hBackFolded->Clone(histNameRatioBackFolded);
  hRatioBackFolded->Reset();
  hRatioBackFolded->Divide(hBackFolded,hist,1,1,"B");
  
  fh1RatioSingleTrackPtBackFolded = hRatioBackFolded;
    
  delete hnHist;
  delete hnFlatEfficiency;
  
}

// ________________________________________________________
void AliFragmentationFunctionCorrections::SecCorrSinglePt()
{
  // apply efficiency correction to inclusive track pt spec

  AddCorrectionLevelSinglePt("secCorr");

  TH1F* histPt = fCorrSinglePt[fNCorrectionLevelsSinglePt-2]->GetTrackPt(0);
  
  if(histPt->GetNbinsX() != fh1SecCorrSinglePt->GetNbinsX())
    Printf("%s:%d: inconsistency pt spec and secondaries corr bins - rebin effCorr ...", (char*)__FILE__,__LINE__);

  TString histNamePt = histPt->GetName();
  TH1F* hTrackPtSecCorr = (TH1F*) histPt->Clone(histNamePt);
  
  hTrackPtSecCorr->Multiply(histPt,fh1SecCorrSinglePt,1,1,"");
  
  fCorrSinglePt[fNCorrectionLevelsSinglePt-1]->AddCorrHistos(0,hTrackPtSecCorr);
}

//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::dNdz2dNdxi()
{
  // transform dN/dz distribution into dN/dxi
  // for current corr level, all jet pt slices
  

  for(Int_t i=0; i<fNJetPtSlices; i++){

    TH1F* histZ = fCorrFF[fNCorrectionLevels-1]->GetZ(i);
    Int_t nBins = histZ->GetNbinsX();

    Double_t* binLims = new Double_t[nBins+1];
    
    for(Int_t bin = 0; bin<nBins; bin++){
      
      Int_t binZ = nBins-bin;

      Double_t zLo = histZ->GetXaxis()->GetBinLowEdge(binZ);
      Double_t zUp = histZ->GetXaxis()->GetBinUpEdge(binZ);
      
      Double_t xiLo = TMath::Log(1/zUp);    
      Double_t xiUp = TMath::Log(1/zLo);

      if(bin == 0) binLims[0] = xiLo;
      binLims[bin+1] = xiUp;   
    }

    // for(Int_t bin = 0; bin<=nBins; bin++) std::cout<<" bin "<<bin<<" binLims "<<binLims[bin]<<std::endl;

    TString strTitle = histZ->GetTitle();
    TString strName  = histZ->GetName();

    strName.ReplaceAll("Z","XiNew");
    strTitle.ReplaceAll("Z","XiNew");
    

    TH1F* histXiNew = new TH1F("histXiNew","",nBins,binLims);
    histXiNew->SetNameTitle(strName,strTitle);
    

    for(Int_t binZ = 1; binZ<=nBins; binZ++){
      
      Double_t meanZ  = histZ->GetBinCenter(binZ); 
      Double_t cont   = histZ->GetBinContent(binZ); 
      Double_t err    = histZ->GetBinError(binZ); 

      Double_t meanXi = TMath::Log(1/meanZ);   
      Int_t binXi = histXiNew->FindBin(meanXi);
      
      histXiNew->SetBinContent(binXi,cont);
      histXiNew->SetBinError(binXi,err);
      
      //std::cout<<" binZ "<<binZ<<" meanZ "<<meanZ<<" binXi "<<binXi<<" meanXi "<<meanXi<<" cont "<<cont<<" err "<<err<<std::endl;
    }
    
    fCorrFF[fNCorrectionLevels-1]->ReplaceCorrHistos(i,0,0,histXiNew);

    delete histXiNew;
    delete[] binLims;
  }
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteBinShiftCorr(TString strInfile, TString strIDGen,  TString strIDRec,  
							    TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim, 
							    TString strOutDir)
{ 
  TString strdirGen  = "PWGJE_FragmentationFunction_" + strIDGen;
  TString strlistGen = "fracfunc_" + strIDGen;
  
  TString strdirRec  = "PWGJE_FragmentationFunction_" + strIDRec;
  TString strlistRec = "fracfunc_" + strIDRec;
  
  WriteBinShiftCorr(strInfile,strdirGen,strlistGen,strdirRec,strlistRec,strOutfile,updateOutfile,useRecPrim, strOutDir, kFALSE, "");
}

//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteBgrBinShiftCorr(TString strInfile, TString strBgrID, TString strIDGen,  TString strIDRec,  
							       TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim, 
							       TString strOutDir)
{ 
  TString strdirGen  = "PWGJE_FragmentationFunction_" + strIDGen;
  TString strlistGen = "fracfunc_" + strIDGen;
  
  TString strdirRec  = "PWGJE_FragmentationFunction_" + strIDRec;
  TString strlistRec = "fracfunc_" + strIDRec;
  
  WriteBinShiftCorr(strInfile,strdirGen,strlistGen,strdirRec,strlistRec,strOutfile,updateOutfile,useRecPrim, strOutDir, kTRUE, strBgrID);
}

//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteBinShiftCorr(TString strInfile, TString strdirGen, TString strlistGen, 
							    TString strdirRec, TString strlistRec, 
							    TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim, 
							    TString strOutDir,  Bool_t writeBgr, TString strBgrID)
{
 
  if((writeBgr && strBgrID.Length() == 0) || (!writeBgr && strBgrID.Length()>0) ){
    Printf("%s:%d -- inconsistent arguments to WriteBinShiftCorr FF/UE", (char*)__FILE__,__LINE__);
    return; 
  }

  TH1F* hCorrPt[fNJetPtSlices];
  TH1F* hCorrXi[fNJetPtSlices];
  TH1F* hCorrZ[fNJetPtSlices];

  TH1F* hdNdptTracksMCGen[fNJetPtSlices];
  TH1F* hdNdxiMCGen[fNJetPtSlices];
  TH1F* hdNdzMCGen[fNJetPtSlices];
    
  TH1F* hdNdptTracksMCRec[fNJetPtSlices];
  TH1F* hdNdxiMCRec[fNJetPtSlices];
  TH1F* hdNdzMCRec[fNJetPtSlices];

  TH1F* fh1FFJetPtMCGen = 0;
  TH1F* fh1FFJetPtMCRec = 0;

  TH2F* fh2FFTrackPtMCGen = 0;
  TH2F* fh2FFZMCGen       = 0;
  TH2F* fh2FFXiMCGen      = 0;
  
  TH2F* fh2FFTrackPtMCRec = 0;
  TH2F* fh2FFZMCRec       = 0;
  TH2F* fh2FFXiMCRec      = 0;
 
  // gen level FF

  TFile f(strInfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- writeBinShiftCorr: open task ouput file %s",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdirGen && strdirGen.Length()) gDirectory->cd(strdirGen);

  TList* listGen = 0;

  if(strlistGen && strlistGen.Length()){

    if(!(listGen = (TList*) gDirectory->Get(strlistGen))){ 
      Printf("%s:%d -- error retrieving listGen %s from directory %s", (char*)__FILE__,__LINE__,strlistGen.Data(),strdirGen.Data());
      return;
    }
  }
 
  if(listGen){ 

    fh1FFJetPtMCGen    = (TH1F*) listGen->FindObject("fh1FFJetPtGen");
    
    if(writeBgr){
      fh2FFTrackPtMCGen  = (TH2F*) listGen->FindObject(Form("fh2FFTrackPt%sGen",strBgrID.Data()));
      fh2FFZMCGen        = (TH2F*) listGen->FindObject(Form("fh2FFZ%sGen",strBgrID.Data()));
      fh2FFXiMCGen       = (TH2F*) listGen->FindObject(Form("fh2FFXi%sGen",strBgrID.Data()));  
    }
    else{
      fh2FFTrackPtMCGen  = (TH2F*) listGen->FindObject("fh2FFTrackPtGen");
      fh2FFZMCGen        = (TH2F*) listGen->FindObject("fh2FFZGen");
      fh2FFXiMCGen       = (TH2F*) listGen->FindObject("fh2FFXiGen"); 
    }
  }
  else{
    fh1FFJetPtMCGen    = (TH1F*) gDirectory->Get("fh1FFJetPtGen");
    if(writeBgr){
      fh2FFTrackPtMCGen  = (TH2F*) gDirectory->Get(Form("fh2FFTrackPt%sGen",strBgrID.Data()));
      fh2FFZMCGen        = (TH2F*) gDirectory->Get(Form("fh2FFZ%sGen",strBgrID.Data()));
      fh2FFXiMCGen       = (TH2F*) gDirectory->Get(Form("fh2FFXi%sGen",strBgrID.Data())); 
    }
    else{
      fh2FFTrackPtMCGen  = (TH2F*) gDirectory->Get("fh2FFTrackPtGen");
      fh2FFZMCGen        = (TH2F*) gDirectory->Get("fh2FFZGen");
      fh2FFXiMCGen       = (TH2F*) gDirectory->Get("fh2FFXiGen"); 
    }
  }
  
  fh1FFJetPtMCGen->SetDirectory(0); 
  fh2FFTrackPtMCGen->SetDirectory(0);
  fh2FFZMCGen->SetDirectory(0);
  fh2FFXiMCGen->SetDirectory(0);

  f.Close();

  // rec level FF

  TFile g(strInfile,"READ");

  if(!g.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- writeBinShiftCorr: open task ouput file %s",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdirRec && strdirRec.Length()) gDirectory->cd(strdirRec);

  TList* listRec = 0;

  if(strlistRec && strlistRec.Length()){

    if(!(listRec = (TList*) gDirectory->Get(strlistRec))){ 
      Printf("%s:%d -- error retrieving listRec %s from directory %s", (char*)__FILE__,__LINE__,strlistRec.Data(),strdirRec.Data());
      return;
    }
  }
 
  if(useRecPrim){
    if(listRec){ 
      fh1FFJetPtMCRec    = (TH1F*) listRec->FindObject("fh1FFJetPtRecEffRec");

      if(writeBgr){
	fh2FFTrackPtMCRec  = (TH2F*) listRec->FindObject(Form("fh2FFTrackPt%sRecEffRec",strBgrID.Data()));
	fh2FFZMCRec        = (TH2F*) listRec->FindObject(Form("fh2FFZ%sRecEffRec",strBgrID.Data()));
	fh2FFXiMCRec       = (TH2F*) listRec->FindObject(Form("fh2FFXi%sRecEffRec",strBgrID.Data())); 
      }
      else{
	fh2FFTrackPtMCRec  = (TH2F*) listRec->FindObject("fh2FFTrackPtRecEffRec");
	fh2FFZMCRec        = (TH2F*) listRec->FindObject("fh2FFZRecEffRec");
	fh2FFXiMCRec       = (TH2F*) listRec->FindObject("fh2FFXiRecEffRec"); 
      }
    }
    else{
      fh1FFJetPtMCRec    = (TH1F*) gDirectory->Get("fh1FFJetPtRecEffRec");
      if(writeBgr){
	fh2FFTrackPtMCRec  = (TH2F*) gDirectory->Get(Form("fh2FFTrackPt%sRecEffRec",strBgrID.Data()));
	fh2FFZMCRec        = (TH2F*) gDirectory->Get(Form("fh2FFZ%sRecEffRec",strBgrID.Data()));
	fh2FFXiMCRec       = (TH2F*) gDirectory->Get(Form("fh2FFXi%sRecEffRec",strBgrID.Data())); 
      }
      else{
	fh2FFTrackPtMCRec  = (TH2F*) gDirectory->Get("fh2FFTrackPtRecEffRec");
	fh2FFZMCRec        = (TH2F*) gDirectory->Get("fh2FFZRecEffRec");
	fh2FFXiMCRec       = (TH2F*) gDirectory->Get("fh2FFXiRecEffRec"); 
      }
    }
  }
  else{
    if(listRec){ 
      fh1FFJetPtMCRec    = (TH1F*) listRec->FindObject("fh1FFJetPtRecCuts");
      if(writeBgr){
	fh2FFTrackPtMCRec  = (TH2F*) listRec->FindObject(Form("fh2FFTrackPt%sRecCuts",strBgrID.Data()));
	fh2FFZMCRec        = (TH2F*) listRec->FindObject(Form("fh2FFZ%sRecCuts",strBgrID.Data()));
	fh2FFXiMCRec       = (TH2F*) listRec->FindObject(Form("fh2FFXi%sRecCuts",strBgrID.Data())); 
      }
      else{
	fh2FFTrackPtMCRec  = (TH2F*) listRec->FindObject("fh2FFTrackPtRecCuts");
	fh2FFZMCRec        = (TH2F*) listRec->FindObject("fh2FFZRecCuts");
	fh2FFXiMCRec       = (TH2F*) listRec->FindObject("fh2FFXiRecCuts");
      }
    }
    else{
      fh1FFJetPtMCRec    = (TH1F*) gDirectory->Get("fh1FFJetPtRecCuts");
      if(writeBgr){
	fh2FFTrackPtMCRec  = (TH2F*) gDirectory->Get(Form("fh2FFTrackPt%sRecCuts",strBgrID.Data()));
	fh2FFZMCRec        = (TH2F*) gDirectory->Get(Form("fh2FFZ%sRecCuts",strBgrID.Data()));
	fh2FFXiMCRec       = (TH2F*) gDirectory->Get(Form("fh2FFXi%sRecCuts",strBgrID.Data()));  
      }
      else{
	fh2FFTrackPtMCRec  = (TH2F*) gDirectory->Get("fh2FFTrackPtRecCuts");
	fh2FFZMCRec        = (TH2F*) gDirectory->Get("fh2FFZRecCuts");
	fh2FFXiMCRec       = (TH2F*) gDirectory->Get("fh2FFXiRecCuts"); 
      }
    }
  }
  

  fh1FFJetPtMCRec->SetDirectory(0); 
  fh2FFTrackPtMCRec->SetDirectory(0);
  fh2FFZMCRec->SetDirectory(0);
  fh2FFXiMCRec->SetDirectory(0);  

  g.Close();

  // projections: FF for generated and reconstructed  
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Float_t jetPtLoLim = fJetPtSlices->At(i);
    Float_t jetPtUpLim = fJetPtSlices->At(i+1);
    
    Int_t binLo = static_cast<Int_t>(fh2FFTrackPtMCGen->GetXaxis()->FindBin(jetPtLoLim));
    Int_t binUp = static_cast<Int_t>(fh2FFTrackPtMCGen->GetXaxis()->FindBin(jetPtUpLim))-1;

    if(binUp > fh2FFTrackPtMCGen->GetNbinsX()){
      Printf("%s:%d -- jet pt range %0.3f exceeds histo limits",(char*)__FILE__,__LINE__,jetPtUpLim); 
      return; 
    }
    
    TString strNameFFPtGen(Form("fh1FFTrackPtGenBinShift_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZGen(Form("fh1FFZGenBinShift_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXiGen(Form("fh1FFXiGenBinShift_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    
    TString strNameFFPtRec(Form("fh1FFTrackPtRecBinShift_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFZRec(Form("fh1FFZRecBinShift_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));
    TString strNameFFXiRec(Form("fh1FFXiRecBinShift_%02d_%02d",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim)));

    if(writeBgr){
      strNameFFPtGen.Form("fh1TrackPtGenBinShiftBgr_%02d_%02d_%s",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim),strBgrID.Data());
      strNameFFZGen.Form("fh1ZGenBinShiftBgr_%02d_%02d_%s",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim),strBgrID.Data());
      strNameFFXiGen.Form("fh1XiGenBinShiftBgr_%02d_%02d_%s",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim),strBgrID.Data());
      
      strNameFFPtRec.Form("fh1TrackPtRecBinShiftBgr_%02d_%02d_%s",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim),strBgrID.Data());
      strNameFFZRec.Form("fh1ZRecBinShiftBgr_%02d_%02d_%s",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim),strBgrID.Data());
      strNameFFXiRec.Form("fh1XiRecBinShiftBgr_%02d_%02d_%s",static_cast<Int_t> (jetPtLoLim),static_cast<Int_t> (jetPtUpLim),strBgrID.Data());
    }
    
    // project 
    // appendix 'unbinned' to avoid histos with same name after rebinning

    hdNdptTracksMCGen[i] = (TH1F*) fh2FFTrackPtMCGen->ProjectionY(strNameFFPtGen+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    hdNdzMCGen[i]        = (TH1F*) fh2FFZMCGen->ProjectionY(strNameFFZGen+"_unBinned",binLo,binUp,"o");
    hdNdxiMCGen[i]       = (TH1F*) fh2FFXiMCGen->ProjectionY(strNameFFXiGen+"_unBinned",binLo,binUp,"o");
    
    hdNdptTracksMCRec[i] = (TH1F*) fh2FFTrackPtMCRec->ProjectionY(strNameFFPtRec+"_unBinned",binLo,binUp,"o"); // option "o": original axis range 
    hdNdzMCRec[i]        = (TH1F*) fh2FFZMCRec->ProjectionY(strNameFFZRec+"_unBinned",binLo,binUp,"o");
    hdNdxiMCRec[i]       = (TH1F*) fh2FFXiMCRec->ProjectionY(strNameFFXiRec+"_unBinned",binLo,binUp,"o");
    

    // rebin

    if(fNHistoBinsPt[i]) hdNdptTracksMCGen[i] = (TH1F*) hdNdptTracksMCGen[i]->Rebin(fNHistoBinsPt[i],strNameFFPtGen,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hdNdzMCGen[i]  = (TH1F*) hdNdzMCGen[i]->Rebin(fNHistoBinsZ[i],strNameFFZGen,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hdNdxiMCGen[i] = (TH1F*) hdNdxiMCGen[i]->Rebin(fNHistoBinsXi[i],strNameFFXiGen,fHistoBinsXi[i]->GetArray());

    if(fNHistoBinsPt[i]) hdNdptTracksMCRec[i] = (TH1F*) hdNdptTracksMCRec[i]->Rebin(fNHistoBinsPt[i],strNameFFPtRec,fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hdNdzMCRec[i]  = (TH1F*) hdNdzMCRec[i]->Rebin(fNHistoBinsZ[i],strNameFFZRec,fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hdNdxiMCRec[i] = (TH1F*) hdNdxiMCRec[i]->Rebin(fNHistoBinsXi[i],strNameFFXiRec,fHistoBinsXi[i]->GetArray());


    hdNdptTracksMCGen[i]->SetNameTitle(strNameFFPtGen,"");
    hdNdzMCGen[i]->SetNameTitle(strNameFFZGen,"");
    hdNdxiMCGen[i]->SetNameTitle(strNameFFXiGen,"");
    
    hdNdptTracksMCRec[i]->SetNameTitle(strNameFFPtRec,"");
    hdNdzMCRec[i]->SetNameTitle(strNameFFZRec,"");
    hdNdxiMCRec[i]->SetNameTitle(strNameFFXiRec,"");
    
     // normalize
    
    Double_t nJetsBinGen = fh1FFJetPtMCGen->Integral(binLo,binUp);
    Double_t nJetsBinRec = fh1FFJetPtMCRec->Integral(binLo,binUp);

    NormalizeTH1(hdNdptTracksMCGen[i],nJetsBinGen); 
    NormalizeTH1(hdNdzMCGen[i],nJetsBinGen); 
    NormalizeTH1(hdNdxiMCGen[i],nJetsBinGen); 

    NormalizeTH1(hdNdptTracksMCRec[i],nJetsBinRec); 
    NormalizeTH1(hdNdzMCRec[i],nJetsBinRec); 
    NormalizeTH1(hdNdxiMCRec[i],nJetsBinRec); 

    // divide gen/rec : corr factor

    TString strNameCorrPt(Form("hBbBCorrPt_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameCorrZ(Form("hBbBCorrZ_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));
    TString strNameCorrXi(Form("hBbBCorrXi_%02d_%02d",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim)));

    if(writeBgr){
      strNameCorrPt.Form("hBbBCorrBgrPt_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
      strNameCorrZ.Form("hBbBCorrBgrZ_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
      strNameCorrXi.Form("hBbBCorrBgrXi_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
    }

    hCorrPt[i] = (TH1F*) hdNdptTracksMCGen[i]->Clone(strNameCorrPt);
    hCorrPt[i]->Divide(hdNdptTracksMCGen[i],hdNdptTracksMCRec[i],1,1,"B"); // binominal errors
    
    hCorrXi[i] = (TH1F*) hdNdxiMCGen[i]->Clone(strNameCorrXi);
    hCorrXi[i]->Divide(hdNdxiMCGen[i],hdNdxiMCRec[i],1,1,"B"); // binominal errors
    
    hCorrZ[i] = (TH1F*) hdNdzMCGen[i]->Clone(strNameCorrZ);
    hCorrZ[i]->Divide(hdNdzMCGen[i],hdNdzMCRec[i],1,1,"B"); // binominal errors
  } 
  
  // write 

  TString outfileOption = "RECREATE";
  if(updateOutfile)  outfileOption = "UPDATE";

  TFile out(strOutfile,outfileOption);
  
  if(!out.IsOpen()){
    Printf("%s:%d -- error opening efficiency output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- write bin shift correction to file %s ",(char*)__FILE__,__LINE__,strOutfile.Data());

  if(strOutDir && strOutDir.Length()){  
    
    TDirectory* dir;
    if((dir = ((TDirectory*) gDirectory->Get(strOutDir)))) dir->cd(); 
    else{
      dir = out.mkdir(strOutDir);
      dir->cd(); 
    } 
  }

  for(Int_t i=0; i<fNJetPtSlices; i++){

    hdNdptTracksMCGen[i]->Write();
    hdNdxiMCGen[i]->Write();
    hdNdzMCGen[i]->Write();
    
    hdNdptTracksMCRec[i]->Write();
    hdNdxiMCRec[i]->Write();
    hdNdzMCRec[i]->Write();
  
    hCorrPt[i]->Write();
    hCorrXi[i]->Write();
    hCorrZ[i]->Write();
  }

  // out.Close();
 
  delete fh1FFJetPtMCGen; 
  delete fh1FFJetPtMCRec; 

  delete fh2FFTrackPtMCGen;
  delete fh2FFZMCGen;
  delete fh2FFXiMCGen;
  
  delete fh2FFTrackPtMCRec;
  delete fh2FFZMCRec;
  delete fh2FFXiMCRec;
}

//________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadBgrBinShiftCorr(TString strfile, TString strBgrID, TString strdir, TString strlist)
{

  ReadBinShiftCorr(strfile, strdir, strlist, kTRUE, strBgrID);
}

//___________________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadBinShiftCorr(TString strfile, TString strdir, TString strlist, Bool_t readBgr, TString strBgrID)
{

  if((readBgr && strBgrID.Length() == 0) || (!readBgr && strBgrID.Length()>0) ){
    Printf("%s:%d -- inconsistent arguments to ReadBinShiftCorr FF/UE", (char*)__FILE__,__LINE__);
    return; 
  }

  // temporary histos to hold histos from file
  TH1F* hCorrPt[fNJetPtSlices]; 
  TH1F* hCorrZ[fNJetPtSlices];
  TH1F* hCorrXi[fNJetPtSlices];
  
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrPt[i] = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrZ[i]  = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrXi[i] = 0;

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read FF / UE bin shift correction from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
   
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Int_t jetPtLoLim = static_cast<Int_t> (fJetPtSlices->At(i));
    Int_t jetPtUpLim = static_cast<Int_t> (fJetPtSlices->At(i+1));
    
    TString strNameCorrPt(Form("hBbBCorrPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameCorrZ(Form("hBbBCorrZ_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameCorrXi(Form("hBbBCorrXi_%02d_%02d",jetPtLoLim,jetPtUpLim));

    if(readBgr){
      strNameCorrPt.Form("hBbBCorrBgrPt_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
      strNameCorrZ.Form("hBbBCorrBgrZ_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
      strNameCorrXi.Form("hBbBCorrBgrXi_%02d_%02d_%s",static_cast<Int_t>(jetPtLoLim),static_cast<Int_t>(jetPtUpLim),strBgrID.Data());
    }

   
    if(list){
      hCorrPt[i] = (TH1F*) list->FindObject(strNameCorrPt); 
      hCorrZ[i]  = (TH1F*) list->FindObject(strNameCorrZ); 
      hCorrXi[i] = (TH1F*) list->FindObject(strNameCorrXi); 
    }
    else{
      hCorrPt[i] = (TH1F*) gDirectory->Get(strNameCorrPt); 
      hCorrZ[i]  = (TH1F*) gDirectory->Get(strNameCorrZ); 
      hCorrXi[i] = (TH1F*) gDirectory->Get(strNameCorrXi); 
    }
    
    if(!hCorrPt[i]){
      Printf("%s:%d -- error retrieving bin by bin correction %s", (char*)__FILE__,__LINE__,strNameCorrPt.Data());
    }
  
    if(!hCorrZ[i]){
      Printf("%s:%d -- error retrieving bin by bin correction %s", (char*)__FILE__,__LINE__,strNameCorrZ.Data());
    }    

    if(!hCorrXi[i]){
      Printf("%s:%d -- error retrieving bin by bin correction %s", (char*)__FILE__,__LINE__,strNameCorrXi.Data());
    }


    if(fNHistoBinsPt[i]) hCorrPt[i] = (TH1F*) hCorrPt[i]->Rebin(fNHistoBinsPt[i],strNameCorrPt+"_rebin",fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hCorrZ[i]  = (TH1F*) hCorrZ[i]->Rebin(fNHistoBinsZ[i],strNameCorrZ+"_rebin",fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hCorrXi[i] = (TH1F*) hCorrXi[i]->Rebin(fNHistoBinsXi[i],strNameCorrXi+"_rebin",fHistoBinsXi[i]->GetArray());

    if(hCorrPt[i]) hCorrPt[i]->SetDirectory(0); 
    if(hCorrZ[i])  hCorrZ[i]->SetDirectory(0); 
    if(hCorrXi[i]) hCorrXi[i]->SetDirectory(0); 

  } // jet slices loop

  f.Close();

  if(readBgr){
    for(Int_t i=0; i<fNJetPtSlices; i++){ // 2nd loop: need to close input file before placing histos
      if(hCorrPt[i]) new(fh1BbBBgrPt[i]) TH1F(*hCorrPt[i]);
      if(hCorrZ[i])  new(fh1BbBBgrZ[i])  TH1F(*hCorrZ[i]);
      if(hCorrXi[i]) new(fh1BbBBgrXi[i]) TH1F(*hCorrXi[i]);
    }
  }
  else{
    for(Int_t i=0; i<fNJetPtSlices; i++){ // 2nd loop: need to close input file before placing histos
      if(hCorrPt[i]) new(fh1BbBPt[i]) TH1F(*hCorrPt[i]);
      if(hCorrZ[i])  new(fh1BbBZ[i])  TH1F(*hCorrZ[i]);
      if(hCorrXi[i]) new(fh1BbBXi[i]) TH1F(*hCorrXi[i]);
    }
  }
}

// ________________________________________________
void AliFragmentationFunctionCorrections::BbBCorr()
{
  // apply bin-by-bin correction

  AddCorrectionLevel("BbB");

  for(Int_t i=0; i<fNJetPtSlices; i++){

    TH1F* histPt = fCorrFF[fNCorrectionLevels-2]->GetTrackPt(i);
    TH1F* histZ  = fCorrFF[fNCorrectionLevels-2]->GetZ(i);
    TH1F* histXi = fCorrFF[fNCorrectionLevels-2]->GetXi(i);

    TString histNamePt = histPt->GetName();
    TString histNameZ  = histZ->GetName();
    TString histNameXi = histXi->GetName();
    
    TH1F* hFFTrackPtBbBCorr = (TH1F*) histPt->Clone(histNamePt);
    hFFTrackPtBbBCorr->Multiply(histPt,fh1BbBPt[i],1,1,"");
    
    TH1F* hFFZBbBCorr = (TH1F*) histZ->Clone(histNameZ);
    hFFZBbBCorr->Multiply(histZ,fh1BbBZ[i],1,1,"");
    
    TH1F* hFFXiBbBCorr = (TH1F*) histXi->Clone(histNameXi);
    hFFXiBbBCorr->Multiply(histXi,fh1BbBXi[i],1,1,"");
    
    fCorrFF[fNCorrectionLevels-1]->AddCorrHistos(i,hFFTrackPtBbBCorr,hFFZBbBCorr,hFFXiBbBCorr);
  }
}

// ___________________________________________________
void AliFragmentationFunctionCorrections::BbBCorrBgr()
{
  // apply bin-by-bin correction

  AddCorrectionLevelBgr("BbB");

  for(Int_t i=0; i<fNJetPtSlices; i++){

    TH1F* histPt = fCorrBgr[fNCorrectionLevelsBgr-2]->GetTrackPt(i);
    TH1F* histZ  = fCorrBgr[fNCorrectionLevelsBgr-2]->GetZ(i);
    TH1F* histXi = fCorrBgr[fNCorrectionLevelsBgr-2]->GetXi(i);

    TString histNamePt = histPt->GetName();
    TString histNameZ  = histZ->GetName();
    TString histNameXi = histXi->GetName();
    
    TH1F* hBgrTrackPtBbBCorr = (TH1F*) histPt->Clone(histNamePt);
    hBgrTrackPtBbBCorr->Multiply(histPt,fh1BbBBgrPt[i],1,1,"");
    
    TH1F* hBgrZBbBCorr = (TH1F*) histZ->Clone(histNameZ);
    hBgrZBbBCorr->Multiply(histZ,fh1BbBBgrZ[i],1,1,"");
    
    TH1F* hBgrXiBbBCorr = (TH1F*) histXi->Clone(histNameXi);
    hBgrXiBbBCorr->Multiply(histXi,fh1BbBBgrXi[i],1,1,"");
    
    fCorrBgr[fNCorrectionLevelsBgr-1]->AddCorrHistos(i,hBgrTrackPtBbBCorr,hBgrZBbBCorr,hBgrXiBbBCorr);
  }
}

//_______________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadFoldingCorr(TString strfile, TString strdir, TString strlist)
{

  // temporary histos to hold histos from file
  TH1F* hCorrPt[fNJetPtSlices]; 
  TH1F* hCorrZ[fNJetPtSlices];
  TH1F* hCorrXi[fNJetPtSlices];
  
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrPt[i] = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrZ[i]  = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrXi[i] = 0;

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read bin shift correction from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
   
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Int_t jetPtLoLim = static_cast<Int_t> (fJetPtSlices->At(i));
    Int_t jetPtUpLim = static_cast<Int_t> (fJetPtSlices->At(i+1));
    
    TString strNameCorrPt(Form("hCorrFoldingPt_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameCorrZ(Form("hCorrFoldingZ_%02d_%02d",jetPtLoLim,jetPtUpLim));
    TString strNameCorrXi(Form("hCorrFoldingXi_%02d_%02d",jetPtLoLim,jetPtUpLim));
     
   
    if(list){
      hCorrPt[i] = (TH1F*) list->FindObject(strNameCorrPt); 
      hCorrZ[i]  = (TH1F*) list->FindObject(strNameCorrZ); 
      hCorrXi[i] = (TH1F*) list->FindObject(strNameCorrXi); 
    }
    else{
      hCorrPt[i] = (TH1F*) gDirectory->Get(strNameCorrPt); 
      hCorrZ[i]  = (TH1F*) gDirectory->Get(strNameCorrZ); 
      hCorrXi[i] = (TH1F*) gDirectory->Get(strNameCorrXi); 
    }
    
    if(!hCorrPt[i]){
      //Printf("%s:%d -- error retrieving folding correction %s", (char*)__FILE__,__LINE__,strNameCorrPt.Data());
    }
  
    if(!hCorrZ[i]){
      //Printf("%s:%d -- error retrieving folding correction %s", (char*)__FILE__,__LINE__,strNameCorrZ.Data());
    }    

    if(!hCorrXi[i]){
      Printf("%s:%d -- error retrieving folding correction %s", (char*)__FILE__,__LINE__,strNameCorrXi.Data());
    }

    if(hCorrPt[i]) hCorrPt[i]->SetDirectory(0); 
    if(hCorrZ[i])  hCorrZ[i]->SetDirectory(0); 
    if(hCorrXi[i]) hCorrXi[i]->SetDirectory(0); 

  } // jet slices loop

  f.Close();

  for(Int_t i=0; i<fNJetPtSlices; i++){ // 2nd loop: need to close input file before placing histos
    if(hCorrPt[i]) new(fh1FoldingCorrPt[i]) TH1F(*hCorrPt[i]);
    if(hCorrZ[i])  new(fh1FoldingCorrZ[i])  TH1F(*hCorrZ[i]);
    if(hCorrXi[i]) new(fh1FoldingCorrXi[i]) TH1F(*hCorrXi[i]);
  }
}

// ___________________________________________________
void AliFragmentationFunctionCorrections::FoldingCorr()
{
  // apply bin-by-bin correction

  AddCorrectionLevel("FoldingCorr");

  for(Int_t i=0; i<fNJetPtSlices; i++){

    TH1F* histPt = fCorrFF[fNCorrectionLevels-2]->GetTrackPt(i);
    TH1F* histZ  = fCorrFF[fNCorrectionLevels-2]->GetZ(i);
    TH1F* histXi = fCorrFF[fNCorrectionLevels-2]->GetXi(i);

    TString histNamePt = histPt->GetName();
    TString histNameZ  = histZ->GetName();
    TString histNameXi = histXi->GetName();
        
    std::cout<<" foldingCorr: i "<<i<<" corr pt "<<fh1FoldingCorrPt[i]<<" z "<<fh1FoldingCorrZ[i]<<" xi "
	     <<fh1FoldingCorrXi[i]<<std::endl;

    std::cout<<" foldingCorr: i "<<i<<" mean corr pt "<<fh1FoldingCorrPt[i]->GetMean()<<" z "<<fh1FoldingCorrZ[i]->GetMean()<<" xi "
	     <<fh1FoldingCorrXi[i]->GetMean()<<std::endl;

    

    TH1F* hFFTrackPtFoldingCorr = (TH1F*) histPt->Clone(histNamePt);
    if(fh1FoldingCorrPt[i] && fh1FoldingCorrPt[i]->GetMean()>0) hFFTrackPtFoldingCorr->Multiply(histPt,fh1FoldingCorrPt[i],1,1,"");
    else hFFTrackPtFoldingCorr->Reset();

    TH1F* hFFZFoldingCorr = (TH1F*) histZ->Clone(histNameZ);
    if(fh1FoldingCorrZ[i] && fh1FoldingCorrZ[i]->GetMean()>0) hFFZFoldingCorr->Multiply(histZ,fh1FoldingCorrZ[i],1,1,"");
    else hFFZFoldingCorr->Reset();
    
    TH1F* hFFXiFoldingCorr = (TH1F*) histXi->Clone(histNameXi);
    if(fh1FoldingCorrXi[i]&& fh1FoldingCorrXi[i]->GetMean()>0) hFFXiFoldingCorr->Multiply(histXi,fh1FoldingCorrXi[i],1,1,"");
    else hFFXiFoldingCorr->Reset();
    
    fCorrFF[fNCorrectionLevels-1]->AddCorrHistos(i,hFFTrackPtFoldingCorr,hFFZFoldingCorr,hFFXiFoldingCorr);
  }
}

//________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadBgrJetSecCorr(TString strfile, TString strBgrID, TString strdir, TString strlist, 
							    Bool_t useScaledStrangeness)
{

  ReadJetSecCorr(strfile, strdir, strlist, useScaledStrangeness, kTRUE, strBgrID);
}

//_______________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadJetSecCorr(TString strfile, TString strdir, TString strlist, Bool_t useScaledStrangeness, 
							 Bool_t readBgr, TString strBgrID){

  // read reconstruction efficiency from file
  // argument strlist optional - read from directory strdir if not specified

  // temporary histos to hold histos from file
  TH1F* hCorrPt[fNJetPtSlices]; 
  TH1F* hCorrZ[fNJetPtSlices];
  TH1F* hCorrXi[fNJetPtSlices];
  
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrPt[i] = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrZ[i]  = 0;
  for(Int_t i=0; i<fNJetPtSlices; i++) hCorrXi[i] = 0;

  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read secondary correction from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
   
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  

  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    Int_t jetPtLoLim = static_cast<Int_t> (fJetPtSlices->At(i));
    Int_t jetPtUpLim = static_cast<Int_t> (fJetPtSlices->At(i+1));

    TString strNameCorrPt("");
    TString strNameCorrZ("");
    TString strNameCorrXi("");
    
    if(readBgr){
      strNameCorrPt.Form("hSecCorrBgrPt_%02d_%02d_%s",jetPtLoLim,jetPtUpLim,strBgrID.Data());
      strNameCorrZ.Form("hSecCorrBgrZ_%02d_%02d_%s",jetPtLoLim,jetPtUpLim,strBgrID.Data());
      strNameCorrXi.Form("hSecCorrBgrXi_%02d_%02d_%s",jetPtLoLim,jetPtUpLim,strBgrID.Data());
      
      if(!useScaledStrangeness){
	Printf("%s:%d -- readJetSecCorr bgr: using naive Pythia strangenss ", (char*)__FILE__,__LINE__);
	strNameCorrPt.Form("hSecCorrBgrPt_nonSc_%02d_%02d_%s",jetPtLoLim,jetPtUpLim,strBgrID.Data());
	strNameCorrZ.Form("hSecCorrBgrZ_nonSc_%02d_%02d_%s",jetPtLoLim,jetPtUpLim,strBgrID.Data());
	strNameCorrXi.Form("hSecCorrBgrXi_nonSc_%02d_%02d_%s",jetPtLoLim,jetPtUpLim,strBgrID.Data());
      }
    }
    else{
      strNameCorrPt.Form("hSecCorrPt_%02d_%02d",jetPtLoLim,jetPtUpLim);
      strNameCorrZ.Form("hSecCorrZ_%02d_%02d",jetPtLoLim,jetPtUpLim);
      strNameCorrXi.Form("hSecCorrXi_%02d_%02d",jetPtLoLim,jetPtUpLim);
      
      if(!useScaledStrangeness){
	Printf("%s:%d -- readJetSecCorr: using naive Pythia strangenss ", (char*)__FILE__,__LINE__);
	strNameCorrPt.Form("hSecCorrPt_nonSc_%02d_%02d",jetPtLoLim,jetPtUpLim);
	strNameCorrZ.Form("hSecCorrZ_nonSc_%02d_%02d",jetPtLoLim,jetPtUpLim);
	strNameCorrXi.Form("hSecCorrXi_nonSc_%02d_%02d",jetPtLoLim,jetPtUpLim);
      }
    }
   
    if(list){
      hCorrPt[i] = (TH1F*) list->FindObject(strNameCorrPt); 
      hCorrZ[i]  = (TH1F*) list->FindObject(strNameCorrZ); 
      hCorrXi[i] = (TH1F*) list->FindObject(strNameCorrXi); 
    }
    else{
      hCorrPt[i] = (TH1F*) gDirectory->Get(strNameCorrPt); 
      hCorrZ[i]  = (TH1F*) gDirectory->Get(strNameCorrZ); 
      hCorrXi[i] = (TH1F*) gDirectory->Get(strNameCorrXi); 
    }
    
    if(!hCorrPt[i]){
      Printf("%s:%d -- error retrieving secondaries correction %s", (char*)__FILE__,__LINE__,strNameCorrPt.Data());
    }
  
    if(!hCorrZ[i]){
      Printf("%s:%d -- error retrieving secondaries correction %s", (char*)__FILE__,__LINE__,strNameCorrZ.Data());
    }    

    if(!hCorrXi[i]){
      Printf("%s:%d -- error retrieving secondaries correction %s", (char*)__FILE__,__LINE__,strNameCorrXi.Data());
    }

    if(fNHistoBinsPt[i]) hCorrPt[i] = (TH1F*) hCorrPt[i]->Rebin(fNHistoBinsPt[i],strNameCorrPt+"_rebin",fHistoBinsPt[i]->GetArray());
    if(fNHistoBinsZ[i])  hCorrZ[i]  = (TH1F*) hCorrZ[i]->Rebin(fNHistoBinsZ[i],strNameCorrZ+"_rebin",fHistoBinsZ[i]->GetArray());
    if(fNHistoBinsXi[i]) hCorrXi[i] = (TH1F*) hCorrXi[i]->Rebin(fNHistoBinsXi[i],strNameCorrXi+"_rebin",fHistoBinsXi[i]->GetArray());

    if(hCorrPt[i]) hCorrPt[i]->SetDirectory(0); 
    if(hCorrZ[i])  hCorrZ[i]->SetDirectory(0); 
    if(hCorrXi[i]) hCorrXi[i]->SetDirectory(0); 

  } // jet slices loop

  f.Close();

  for(Int_t i=0; i<fNJetPtSlices; i++){ // 2nd loop: need to close input file before placing histos

    if(readBgr){
      if(hCorrPt[i]) new(fh1SecCorrBgrPt[i]) TH1F(*hCorrPt[i]);
      if(hCorrZ[i])  new(fh1SecCorrBgrZ[i])  TH1F(*hCorrZ[i]);
      if(hCorrXi[i]) new(fh1SecCorrBgrXi[i]) TH1F(*hCorrXi[i]);
    }
    else{
      if(hCorrPt[i]) new(fh1SecCorrPt[i]) TH1F(*hCorrPt[i]);
      if(hCorrZ[i])  new(fh1SecCorrZ[i])  TH1F(*hCorrZ[i]);
      if(hCorrXi[i]) new(fh1SecCorrXi[i]) TH1F(*hCorrXi[i]);
    }
  }
}

// ___________________________________________________
void AliFragmentationFunctionCorrections::JetSecCorr()
{
  // apply secondaries correction

  AddCorrectionLevel("SecCorr");

  Printf("%s:%d -- apply jet secondaries correction", (char*)__FILE__,__LINE__);

  for(Int_t i=0; i<fNJetPtSlices; i++){

    TH1F* histPt = fCorrFF[fNCorrectionLevels-2]->GetTrackPt(i);
    TH1F* histZ  = fCorrFF[fNCorrectionLevels-2]->GetZ(i);
    TH1F* histXi = fCorrFF[fNCorrectionLevels-2]->GetXi(i);

    TString histNamePt = histPt->GetName();
    TString histNameZ  = histZ->GetName();
    TString histNameXi = histXi->GetName();
    
    TH1F* hFFTrackPtSecCorr = (TH1F*) histPt->Clone(histNamePt);
    hFFTrackPtSecCorr->Multiply(histPt,fh1SecCorrPt[i],1,1,"");
    
    TH1F* hFFZSecCorr = (TH1F*) histZ->Clone(histNameZ);
    hFFZSecCorr->Multiply(histZ,fh1SecCorrZ[i],1,1,"");
    
    TH1F* hFFXiSecCorr = (TH1F*) histXi->Clone(histNameXi);
    hFFXiSecCorr->Multiply(histXi,fh1SecCorrXi[i],1,1,"");
    
    fCorrFF[fNCorrectionLevels-1]->AddCorrHistos(i,hFFTrackPtSecCorr,hFFZSecCorr,hFFXiSecCorr);
  }
}



// ___________________________________________________
void AliFragmentationFunctionCorrections::JetSecCorrBgr()
{
  // apply secondaries correction to UE
  
  AddCorrectionLevelBgr("SecCorr"); 

  Printf("%s:%d -- apply jet secondaries correction", (char*)__FILE__,__LINE__);
  
  for(Int_t i=0; i<fNJetPtSlices; i++){
    
    TH1F* histPt = fCorrBgr[fNCorrectionLevelsBgr-2]->GetTrackPt(i);
    TH1F* histZ  = fCorrBgr[fNCorrectionLevelsBgr-2]->GetZ(i);
    TH1F* histXi = fCorrBgr[fNCorrectionLevelsBgr-2]->GetXi(i);

    TString histNamePt = histPt->GetName();
    TString histNameZ  = histZ->GetName();
    TString histNameXi = histXi->GetName();
    
    TH1F* hFFTrackPtSecCorr = (TH1F*) histPt->Clone(histNamePt);
    hFFTrackPtSecCorr->Multiply(histPt,fh1SecCorrBgrPt[i],1,1,"");
    
    TH1F* hFFZSecCorr = (TH1F*) histZ->Clone(histNameZ);
    hFFZSecCorr->Multiply(histZ,fh1SecCorrBgrZ[i],1,1,"");
    
    TH1F* hFFXiSecCorr = (TH1F*) histXi->Clone(histNameXi);
    hFFXiSecCorr->Multiply(histXi,fh1SecCorrBgrXi[i],1,1,"");
    
    fCorrBgr[fNCorrectionLevelsBgr-1]->AddCorrHistos(i,hFFTrackPtSecCorr,hFFZSecCorr,hFFXiSecCorr);
  }
}


//________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteBinShiftCorrSinglePt(TString strInfile, TString strIDGen,  TString strIDRec,  
								    TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim, 
								    TString strOutDir)
{ 
  TString strdirGen  = "PWGJE_FragmentationFunction_" + strIDGen;
  TString strlistGen = "fracfunc_" + strIDGen;
  
  TString strdirRec  = "PWGJE_FragmentationFunction_" + strIDRec;
  TString strlistRec = "fracfunc_" + strIDRec;
  
  WriteBinShiftCorrSinglePt(strInfile,strdirGen,strlistGen,strdirRec,strlistRec,strOutfile,updateOutfile,useRecPrim,strOutDir);
}

//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::WriteBinShiftCorrSinglePt(TString strInfile, TString strdirGen, TString strlistGen, 
								    TString strdirRec, TString strlistRec, 
								    TString strOutfile, Bool_t updateOutfile, Bool_t useRecPrim, 
								    TString strOutDir){
 
  
  TH1F* hCorrPt           = 0;
  TH1F* hdNdptTracksMCGen = 0;
  TH1F* hdNdptTracksMCRec = 0;
 
  // gen level FF

  TFile f(strInfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- writeBinShiftCorrSinglePt: open task ouput file %s",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdirGen && strdirGen.Length()) gDirectory->cd(strdirGen);

  TList* listGen = 0;

  if(strlistGen && strlistGen.Length()){

    if(!(listGen = (TList*) gDirectory->Get(strlistGen))){ 
      Printf("%s:%d -- error retrieving listGen %s from directory %s", (char*)__FILE__,__LINE__,strlistGen.Data(),strdirGen.Data());
      return;
    }
  }
 
  if(listGen){
    hdNdptTracksMCGen  = (TH1F*) listGen->FindObject("fh1TrackQAPtGen");
  }
  else{
    hdNdptTracksMCGen  = (TH1F*) gDirectory->Get("fh1TrackQAPtGen");
  }

  hdNdptTracksMCGen->SetDirectory(0);
  
  f.Close();

  // rec level FF

  TFile g(strInfile,"READ");

  if(!g.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strInfile.Data());
    return;
  }
  
  if(fDebug>0) Printf("%s:%d -- writeBinShiftCorrSinglePt: open task ouput file %s",(char*)__FILE__,__LINE__,strInfile.Data());
 
  if(strdirRec && strdirRec.Length()) gDirectory->cd(strdirRec);

  TList* listRec = 0;

  if(strlistRec && strlistRec.Length()){

    if(!(listRec = (TList*) gDirectory->Get(strlistRec))){ 
      Printf("%s:%d -- error retrieving listRec %s from directory %s", (char*)__FILE__,__LINE__,strlistRec.Data(),strdirRec.Data());
      return;
    }
  }
 

  if(useRecPrim){
    if(listRec){ 
      hdNdptTracksMCRec = (TH1F*) listRec->FindObject("fh1TrackQAPtRecEffRec");
    }
    else{
      hdNdptTracksMCRec = (TH1F*) gDirectory->Get("fh1TrackQAPtRecEffRec");
    }
  }
  else{
    if(listRec){
      hdNdptTracksMCRec  = (TH1F*) listRec->FindObject("fh1TrackQAPtRecCuts");
    } 
    else{
      hdNdptTracksMCRec  = (TH1F*) gDirectory->Get("fh1TrackQAPtRecCuts");
    }
  }
  
  hdNdptTracksMCRec->SetDirectory(0);
  
  g.Close();
    
  TString strNamePtGen = "fh1SinglePtGenBbB";
  TString strNamePtRec = "fh1SinglePtRecBbB";

  // rebin
  if(fNHistoBinsSinglePt) hdNdptTracksMCGen = (TH1F*) hdNdptTracksMCGen->Rebin(fNHistoBinsSinglePt,strNamePtGen+"_rebin",fHistoBinsSinglePt->GetArray());
  if(fNHistoBinsSinglePt) hdNdptTracksMCRec = (TH1F*) hdNdptTracksMCRec->Rebin(fNHistoBinsSinglePt,strNamePtRec+"_rebin",fHistoBinsSinglePt->GetArray());
  
  hdNdptTracksMCGen->SetNameTitle(strNamePtGen,"");
  hdNdptTracksMCRec->SetNameTitle(strNamePtRec,"");
  
  // corr fac
  TString strTitCorr        = "hBbBCorrSinglePt";
  if(useRecPrim) strTitCorr = "hBbBCorrRecPrimSinglePt";

  hCorrPt = (TH1F*) hdNdptTracksMCGen->Clone(strTitCorr);
  hCorrPt->Divide(hdNdptTracksMCGen,hdNdptTracksMCRec,1,1,"B"); // binominal errors
  
  // write 

  TString outfileOption = "RECREATE";
  if(updateOutfile)  outfileOption = "UPDATE";

  TFile out(strOutfile,outfileOption);
  
  if(!out.IsOpen()){
    Printf("%s:%d -- error opening efficiency output file %s", (char*)__FILE__,__LINE__,strOutfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- write bin shift correction to file %s ",(char*)__FILE__,__LINE__,strOutfile.Data());

  if(strOutDir && strOutDir.Length()){  
    
    TDirectory* dir;
    if((dir = ((TDirectory*) gDirectory->Get(strOutDir)))) dir->cd(); 
    else{
      dir = out.mkdir(strOutDir);
      dir->cd(); 
    } 
  }
  
  
  hdNdptTracksMCGen->Write();
  hdNdptTracksMCRec->Write();
  hCorrPt->Write();
  
  out.Close();
 
  delete hdNdptTracksMCGen; 
  delete hdNdptTracksMCRec;  
  delete hCorrPt; 
}

//___________________________________________________________________________________________________________________________________
void AliFragmentationFunctionCorrections::ReadBinShiftCorrSinglePt(TString strfile, TString strdir, TString strlist, Bool_t useRecPrim)
{
  // read reconstruction efficiency from file
  // argument strlist optional - read from directory strdir if not specified

  // temporary histos to hold histos from file
  TH1F* hBbBCorrPt = 0; 
  
  TFile f(strfile,"READ");

  if(!f.IsOpen()){
    Printf("%s:%d -- error opening raw data file %s", (char*)__FILE__,__LINE__,strfile.Data());
    return;
  }

  if(fDebug>0) Printf("%s:%d -- read BbB corr from file %s ",(char*)__FILE__,__LINE__,strfile.Data());
 
  if(strdir && strdir.Length()) gDirectory->cd(strdir);

  TList* list = 0;

  if(strlist && strlist.Length()){
   
    if(!(list = (TList*) gDirectory->Get(strlist))){ 
      Printf("%s:%d -- error retrieving list %s from directory %s", (char*)__FILE__,__LINE__,strlist.Data(),strdir.Data());
      return;
    }
  }  

    
  TString strNameBbBCorrPt = "hBbBCorrSinglePt";
  if(useRecPrim) strNameBbBCorrPt = "hBbBCorrRecPrimSinglePt";     

  if(list){
    hBbBCorrPt = (TH1F*)  list->FindObject(strNameBbBCorrPt); 
  }
  else{
    hBbBCorrPt = (TH1F*)  gDirectory->Get(strNameBbBCorrPt); 
  }
    
  if(!hBbBCorrPt){
    Printf("%s:%d -- error retrieving BbB corr single pt %s", (char*)__FILE__,__LINE__,strNameBbBCorrPt.Data());
  }
 

  if(fNHistoBinsPt) hBbBCorrPt = (TH1F*) hBbBCorrPt->Rebin(fNHistoBinsSinglePt,strNameBbBCorrPt+"_rebin",fHistoBinsSinglePt->GetArray());

  if(hBbBCorrPt) hBbBCorrPt->SetDirectory(0); 
 
  f.Close();

  fh1BbBCorrSinglePt = hBbBCorrPt;
}

// ------------------------------------------------------------------

void AliFragmentationFunctionCorrections::BbBCorrSinglePt()
{
  // apply efficiency correction to inclusive track pt spec

  AddCorrectionLevelSinglePt("BbB");

  TH1F* histPt = fCorrSinglePt[fNCorrectionLevelsSinglePt-2]->GetTrackPt(0);
  
  if(histPt->GetNbinsX() != fh1BbBCorrSinglePt->GetNbinsX()){
    Printf("%s:%d: inconsistency pt spec and BbB corr binning ", (char*)__FILE__,__LINE__);
    return; 
  }

  TString histNamePt    = histPt->GetName();
  TH1F* hTrackPtBbBCorr = (TH1F*) histPt->Clone(histNamePt);
  
  hTrackPtBbBCorr->Multiply(histPt,fh1BbBCorrSinglePt,1,1,"");
  
  fCorrSinglePt[fNCorrectionLevelsSinglePt-1]->AddCorrHistos(0,hTrackPtBbBCorr);
}

