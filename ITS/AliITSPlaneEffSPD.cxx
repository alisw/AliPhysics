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
///////////////////////////////////////////////////////////////////////////
//  Plane Efficiency class for ITS                      
//  It is used for chip by chip efficiency of the SPD,        
//  evaluated by tracks
//  (Inherits from AliITSPlaneEff)
//  Author: G.E. Bruno 
//          giuseppe.bruno@ba.infn.it
//
///////////////////////////////////////////////////////////////////////////

/* $Id$ */

#include <TMath.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include "AliITSPlaneEffSPD.h"
#include "AliLog.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
//#include "AliCDBRunRange.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSCalibrationSPD.h"

ClassImp(AliITSPlaneEffSPD)	
//______________________________________________________________________
AliITSPlaneEffSPD::AliITSPlaneEffSPD():
  AliITSPlaneEff(),
  fHisResX(0),
  fHisResZ(0),
  fHisResXZ(0),
  fHisClusterSize(0),
  fHisResXclu(0),
  fHisResZclu(0),
  fHisResXchip(0),
  fHisResZchip(0),
  fProfResXvsPhi(0),
  fProfResZvsDip(0),
  fProfResXvsPhiclu(0), 
  fProfResZvsDipclu(0),
  fHisTrackErrX(0),
  fHisTrackErrZ(0),
  fHisClusErrX(0),
  fHisClusErrZ(0),
  fHisTrackXFOtrue(0),
  fHisTrackZFOtrue(0),
  fHisTrackXFOfalse(0),
  fHisTrackZFOfalse(0),
  fHisTrackXZFOtrue(0),
  fHisTrackXZFOfalse(0){
  for (UInt_t i=0; i<kNModule*kNChip*(kNClockPhase+1); i++){
    fFound[i]=0;
    fTried[i]=0;
  }
  // default constructor
  AliDebug(1,Form("Calling default constructor"));
}
//______________________________________________________________________
AliITSPlaneEffSPD::~AliITSPlaneEffSPD(){
    // destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //     none.
    DeleteHistos();
}
//______________________________________________________________________
AliITSPlaneEffSPD::AliITSPlaneEffSPD(const AliITSPlaneEffSPD &s) : AliITSPlaneEff(s), 
//fHis(s.fHis),
fHisResX(0),
fHisResZ(0),
fHisResXZ(0),
fHisClusterSize(0),
fHisResXclu(0),
fHisResZclu(0),
fHisResXchip(0),
fHisResZchip(0),
fProfResXvsPhi(0),
fProfResZvsDip(0),
fProfResXvsPhiclu(0),
fProfResZvsDipclu(0),
fHisTrackErrX(0),
fHisTrackErrZ(0),
fHisClusErrX(0),
fHisClusErrZ(0),
fHisTrackXFOtrue(0),
fHisTrackZFOtrue(0),
fHisTrackXFOfalse(0),
fHisTrackZFOfalse(0),
fHisTrackXZFOtrue(0),
fHisTrackXZFOfalse(0)
{
    //     Copy Constructor
    // Inputs:
    //    AliITSPlaneEffSPD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:

 for (UInt_t i=0; i<kNModule*kNChip*(kNClockPhase+1); i++){
    fFound[i]=s.fFound[i];
    fTried[i]=s.fTried[i];
 }
 if(fHis) { 
   InitHistos();
   for(Int_t i=0; i<kNHisto; i++) {
      s.fHisResX[i]->Copy(*fHisResX[i]);
      s.fHisResZ[i]->Copy(*fHisResZ[i]);
      s.fHisResXZ[i]->Copy(*fHisResXZ[i]);
      s.fHisClusterSize[i]->Copy(*fHisClusterSize[i]);
      for(Int_t clu=0; clu<kNclu; clu++) {  // clu=0 --> cluster size 1
        s.fHisResXclu[i][clu]->Copy(*fHisResXclu[i][clu]);
        s.fHisResZclu[i][clu]->Copy(*fHisResZclu[i][clu]);
        s.fProfResXvsPhiclu[i][clu]->Copy(*fProfResXvsPhiclu[i][clu]);
        s.fProfResZvsDipclu[i][clu]->Copy(*fProfResZvsDipclu[i][clu]);
      }
      for(Int_t chip=0; chip<kNChip; chip++) { 
        s.fHisResXchip[i][chip]->Copy(*fHisResXchip[i][chip]);
        s.fHisResZchip[i][chip]->Copy(*fHisResZchip[i][chip]);
      }
      s.fProfResXvsPhi[i]->Copy(*fProfResXvsPhi[i]);
      s.fProfResZvsDip[i]->Copy(*fProfResZvsDip[i]);
      s.fHisTrackErrX[i]->Copy(*fHisTrackErrX[i]);
      s.fHisTrackErrZ[i]->Copy(*fHisTrackErrZ[i]);
      s.fHisClusErrX[i]->Copy(*fHisClusErrX[i]);
      s.fHisClusErrZ[i]->Copy(*fHisClusErrZ[i]);
      for(Int_t phas=0; phas<kNClockPhase;phas++){
        s.fHisTrackXFOtrue[i][phas]->Copy(*fHisTrackXFOtrue[i][phas]);
        s.fHisTrackZFOtrue[i][phas]->Copy(*fHisTrackXFOtrue[i][phas]);
        s.fHisTrackXFOfalse[i][phas]->Copy(*fHisTrackXFOtrue[i][phas]);
        s.fHisTrackZFOfalse[i][phas]->Copy(*fHisTrackXFOtrue[i][phas]);
        s.fHisTrackXZFOtrue[i][phas]->Copy(*fHisTrackXFOtrue[i][phas]);
        s.fHisTrackXZFOfalse[i][phas]->Copy(*fHisTrackXFOtrue[i][phas]);
      }
   }
 }
}
//_________________________________________________________________________
AliITSPlaneEffSPD& AliITSPlaneEffSPD::operator+=(const AliITSPlaneEffSPD &add){
    //    Add-to-me operator
    // Inputs:
    //    const AliITSPlaneEffSPD &add  simulation class to be added
    // Outputs:
    //    none.
    // Return:
    //    none
    for (UInt_t i=0; i<kNModule*kNChip*(kNClockPhase+1); i++){
      fFound[i] += add.fFound[i];
      fTried[i] += add.fTried[i];
    }
    if(fHis && add.fHis) {
      for(Int_t i=0; i<kNHisto; i++) {
        fHisResX[i]->Add(add.fHisResX[i]); 
        fHisResZ[i]->Add(add.fHisResZ[i]); 
        fHisResXZ[i]->Add(add.fHisResXZ[i]); 
        fHisClusterSize[i]->Add(add.fHisClusterSize[i]); 
        for(Int_t clu=0; clu<kNclu; clu++) {  // clu=0 --> cluster size 1
          fHisResXclu[i][clu]->Add(add.fHisResXclu[i][clu]); 
          fHisResZclu[i][clu]->Add(add.fHisResZclu[i][clu]); 
          fProfResXvsPhiclu[i][clu]->Add(add.fProfResXvsPhiclu[i][clu]);
          fProfResZvsDipclu[i][clu]->Add(add.fProfResZvsDipclu[i][clu]);
        }
        for(Int_t chip=0; chip<kNChip; chip++) {  
          fHisResXchip[i][chip]->Add(add.fHisResXchip[i][chip]); 
          fHisResZchip[i][chip]->Add(add.fHisResZchip[i][chip]); 
        }
        fProfResXvsPhi[i]->Add(add.fProfResXvsPhi[i]);
        fProfResZvsDip[i]->Add(add.fProfResZvsDip[i]);
        fHisTrackErrX[i]->Add(add.fHisTrackErrX[i]);
        fHisTrackErrZ[i]->Add(add.fHisTrackErrZ[i]);
        fHisClusErrX[i]->Add(add.fHisClusErrX[i]);
        fHisClusErrZ[i]->Add(add.fHisClusErrZ[i]);
        for(Int_t phas=0; phas<kNClockPhase;phas++){
          fHisTrackXFOtrue[i][phas]->Add(add.fHisTrackXFOtrue[i][phas]);
          fHisTrackZFOtrue[i][phas]->Add(add.fHisTrackXFOtrue[i][phas]);
          fHisTrackXFOfalse[i][phas]->Add(add.fHisTrackXFOtrue[i][phas]);
          fHisTrackZFOfalse[i][phas]->Add(add.fHisTrackXFOtrue[i][phas]);
          fHisTrackXZFOtrue[i][phas]->Add(add.fHisTrackXFOtrue[i][phas]);
          fHisTrackXZFOfalse[i][phas]->Add(add.fHisTrackXFOtrue[i][phas]);
        }
      }
    }
    return *this;
}
//______________________________________________________________________
AliITSPlaneEffSPD&  AliITSPlaneEffSPD::operator=(const
                                           AliITSPlaneEffSPD &s){
    //    Assignment operator
    // Inputs:
    //    AliITSPlaneEffSPD &s The original class for which
    //                                this class is a copy of
    // Outputs:
    //    none.
    // Return:
 
    if(this==&s) return *this;
    s.Copy(*this);
    return *this;
}
//______________________________________________________________________
void AliITSPlaneEffSPD::Copy(TObject &obj) const {
  // protected method. copy this to obj
  AliITSPlaneEff::Copy(obj);
  AliITSPlaneEffSPD& target = (AliITSPlaneEffSPD &) obj;
  for(Int_t i=0;i<kNModule*kNChip*(kNClockPhase+1);i++) {
      target.fFound[i] = fFound[i];
      target.fTried[i] = fTried[i];
  }
  CopyHistos(target);
  return;
}
//_______________________________________________________________________
void AliITSPlaneEffSPD::CopyHistos(AliITSPlaneEffSPD &target) const {
  // protected method: copy histos from this to target
  target.fHis  = fHis; // this is redundant only in some cases. Leave as it is.
  if(fHis) {
    target.fHisResX=new TH1F*[kNHisto];
    target.fHisResZ=new TH1F*[kNHisto];
    target.fHisResXZ=new TH2F*[kNHisto];
    target.fHisClusterSize=new TH2I*[kNHisto];
    target.fHisResXclu=new TH1F**[kNHisto];
    target.fHisResZclu=new TH1F**[kNHisto];
    target.fHisResXchip=new TH1F**[kNHisto];
    target.fHisResZchip=new TH1F**[kNHisto];
    target.fProfResXvsPhi=new TProfile*[kNHisto];
    target.fProfResZvsDip=new TProfile*[kNHisto];
    target.fProfResXvsPhiclu=new TProfile**[kNHisto];
    target.fProfResZvsDipclu=new TProfile**[kNHisto];
    target.fHisTrackErrX=new TH1F*[kNHisto];
    target.fHisTrackErrZ=new TH1F*[kNHisto];
    target.fHisClusErrX=new TH1F*[kNHisto];
    target.fHisClusErrZ=new TH1F*[kNHisto];
    target.fHisTrackXFOtrue=new TH1F**[kNHisto];
    target.fHisTrackZFOtrue=new TH1F**[kNHisto];
    target.fHisTrackXFOfalse=new TH1F**[kNHisto];
    target.fHisTrackZFOfalse=new TH1F**[kNHisto];
    target.fHisTrackXZFOtrue=new TH2F**[kNHisto];
    target.fHisTrackXZFOfalse=new TH2F**[kNHisto];
    for(Int_t i=0; i<kNHisto; i++) {
      target.fHisResX[i] = new TH1F(*fHisResX[i]);
      target.fHisResZ[i] = new TH1F(*fHisResZ[i]);
      target.fHisResXZ[i] = new TH2F(*fHisResXZ[i]);
      target.fHisClusterSize[i] = new TH2I(*fHisClusterSize[i]);
      target.fHisResXclu[i]=new TH1F*[kNclu];
      target.fHisResZclu[i]=new TH1F*[kNclu];
      target.fProfResXvsPhiclu[i]=new TProfile*[kNclu];
      target.fProfResZvsDipclu[i]=new TProfile*[kNclu];
      for(Int_t clu=0; clu<kNclu; clu++) {  // clu=0 --> cluster size 1
        target.fHisResXclu[i][clu] = new TH1F(*fHisResXclu[i][clu]);
        target.fHisResZclu[i][clu] = new TH1F(*fHisResZclu[i][clu]);
        target.fProfResXvsPhiclu[i][clu] = new TProfile(*fProfResXvsPhiclu[i][clu]);
        target.fProfResZvsDipclu[i][clu] = new TProfile(*fProfResZvsDipclu[i][clu]);
      }
      target.fHisResXchip[i]=new TH1F*[kNChip];
      target.fHisResZchip[i]=new TH1F*[kNChip];
      for(Int_t chip=0; chip<kNChip; chip++) {  
        target.fHisResXchip[i][chip] = new TH1F(*fHisResXchip[i][chip]);
        target.fHisResZchip[i][chip] = new TH1F(*fHisResZchip[i][chip]);
      }
      target.fProfResXvsPhi[i] = new TProfile(*fProfResXvsPhi[i]);
      target.fProfResZvsDip[i] = new TProfile(*fProfResZvsDip[i]);
      target.fHisTrackErrX[i] = new TH1F(*fHisTrackErrX[i]);
      target.fHisTrackErrZ[i] = new TH1F(*fHisTrackErrZ[i]);
      target.fHisClusErrX[i] = new TH1F(*fHisClusErrX[i]);
      target.fHisClusErrZ[i] = new TH1F(*fHisClusErrZ[i]);

      target.fHisTrackXFOtrue[i]=new TH1F*[kNClockPhase];
      target.fHisTrackZFOtrue[i]=new TH1F*[kNClockPhase];
      target.fHisTrackXFOfalse[i]=new TH1F*[kNClockPhase];
      target.fHisTrackZFOfalse[i]=new TH1F*[kNClockPhase];
      target.fHisTrackXZFOtrue[i]=new TH2F*[kNClockPhase];
      target.fHisTrackXZFOfalse[i]=new TH2F*[kNClockPhase];
      for(Int_t phas=0; phas<kNClockPhase;phas++){
      target.fHisTrackXFOtrue[i][phas]=new TH1F(*fHisTrackXFOtrue[i][phas]);
      target.fHisTrackZFOtrue[i][phas]=new TH1F(*fHisTrackZFOtrue[i][phas]);
      target.fHisTrackXFOfalse[i][phas]=new TH1F(*fHisTrackXFOfalse[i][phas]);
      target.fHisTrackZFOfalse[i][phas]=new TH1F(*fHisTrackZFOfalse[i][phas]);
      target.fHisTrackXZFOtrue[i][phas]=new TH2F(*fHisTrackXZFOtrue[i][phas]);
      target.fHisTrackXZFOfalse[i][phas]=new TH2F(*fHisTrackXZFOfalse[i][phas]);
      }
    }
  }
return;
}

//_______________________________________________________________________
Int_t AliITSPlaneEffSPD::GetMissingTracksForGivenEff(Double_t eff, Double_t RelErr,
          UInt_t im, UInt_t ic) const {
   
  //   Estimate the number of tracks still to be collected to attain a 
  //   given efficiency eff, with relative error RelErr
  //   Inputs:
  //         eff    -> Expected efficiency (e.g. those from actual estimate)
  //         RelErr -> tollerance [0,1] 
  //         im     -> module number [0,239]
  //         ic     -> chip number [0,4]
  //   Outputs: none
  //   Return: the estimated n. of tracks 
  //
if (im>=kNModule || ic>=kNChip) 
 {AliError("GetMissingTracksForGivenEff: you asked for a non existing chip");
 return -1;}
else { 
  UInt_t key=GetKey(im,ic);
  if(key<kNModule*kNChip) return GetNTracksForGivenEff(eff,RelErr)-fTried[key];
  else return -1;
}
}
//_________________________________________________________________________
Double_t  AliITSPlaneEffSPD::PlaneEff(const UInt_t im,const UInt_t ic, const Bool_t fo, const UInt_t bcm4) const {
// Compute the efficiency for a basic block, 
// Inputs:
//        im     -> module number [0,239]
//        ic     -> chip number [0,4] 
//        fo     -> boolean, true in case of Fast Or studies
//        bcm4   -> for Fast Or: bunch crossing % 4
if (im>=kNModule || ic>=kNChip) 
 {AliError("PlaneEff(Uint_t,Uint_t): you asked for a non existing chip"); return -1.;}
if(fo && bcm4>=kNClockPhase)
 {AliError("PlaneEff(Uint_t,Uint_t): you asked for Fast Or in a wrong phase"); return -1.;}
Int_t nf=-1;
Int_t nt=-1;
if(fo) {
 AliWarning("PlaneEff: you asked for FO efficiency");
 UInt_t key=GetKey(im,ic,fo,bcm4);
 if(key<kNModule*kNChip*(kNClockPhase+1)) {
   nf=fFound[key];
   nt=fTried[key];
 }
} else {
 UInt_t key=GetKey(im,ic);
 if (key<kNModule*kNChip) {
  nf=fFound[key];
  nt=fTried[key];
 }
}
return AliITSPlaneEff::PlaneEff(nf,nt);
}
//_________________________________________________________________________
Double_t  AliITSPlaneEffSPD::ErrPlaneEff(const UInt_t im,const UInt_t ic, const Bool_t fo, const UInt_t bcm4) const {
    // Compute the statistical error on efficiency for a basic block,
    // using binomial statistics 
    // Inputs:
    //        im     -> module number [0,239]
    //        ic     -> chip number [0,4] 
//        fo     -> boolean, true in case of Fast Or studies
//        bcm4   -> for Fast Or: bunch crossing % 4
if (im>=kNModule || ic>=kNChip) 
 {AliError("ErrPlaneEff(Uint_t,Uint_t): you asked for a non existing chip"); return -1.;}
if(fo && bcm4>=kNClockPhase)
 {AliError("PlaneEff(Uint_t,Uint_t): you asked for Fast Or in a wrong phase"); return -1.;}
Int_t nf=-1;
Int_t nt=-1;
if(fo) {
 AliWarning("ErrPlaneEff: you asked for FO efficiency");
 UInt_t key=GetKey(im,ic,fo,bcm4);
 if(key<kNModule*kNChip*(kNClockPhase+1)) {
   nf=fFound[key];
   nt=fTried[key];
 }
} else {
 UInt_t key=GetKey(im,ic);
 if (key<kNModule*kNChip) {
   nf=fFound[key];
   nt=fTried[key];
 }
}
return AliITSPlaneEff::ErrPlaneEff(nf,nt);
} 
//_________________________________________________________________________
Bool_t AliITSPlaneEffSPD::UpDatePlaneEff(const Bool_t Kfound,
                                         const UInt_t im, const UInt_t ic, const Bool_t fo, const UInt_t bcm4) {
  // Update efficiency for a basic block
if (im>=kNModule || ic>=kNChip) 
 {AliError("UpDatePlaneEff: you asked for a non existing chip"); return kFALSE;}
if(fo && bcm4>=kNClockPhase)
 {AliError("UpDatePlaneEff: you asked for Fast Or in a wrong phase"); return kFALSE;}
if (!fo) {
 UInt_t key=GetKey(im,ic);
 if(key<kNModule*kNChip) {
   fTried[key]++;
   if(Kfound) fFound[key]++;
   return kTRUE;
 }
}
else {
 UInt_t key=GetKey(im,ic,fo,bcm4);
 if(key<kNModule*kNChip*(kNClockPhase+1)) {
   fTried[key]++;
   if(Kfound) fFound[key]++;
   return kTRUE;
 }
}
return kFALSE;
}
//_________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetChipFromCol(const UInt_t col) const {
  // get chip given the column
if(col>=kNCol*kNChip) 
 {AliDebug(1,Form("GetChipFromCol: you asked for a non existing column %d",col)); return 10;}
return col/kNCol;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetKey(const UInt_t mod, const UInt_t chip, const Bool_t FO, const UInt_t BCm4) const {
  // get key given a basic block
UInt_t key=99999;
if(mod>=kNModule || chip>=kNChip)
  {AliDebug(1,"GetKey: you asked for a non existing block"); return 99999;}
key = mod*kNChip+chip;
if(FO) { 
  if(BCm4>= kNClockPhase) {AliDebug(1,"GetKey: you have asked Fast OR and a non exisiting BC modulo 4"); return 99999;}
  key += kNModule*kNChip*(BCm4+1);
}
return key;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSPD::SwitchChipKeyNumbering(UInt_t key) const {

// methods to switch from offline chip key numbering 
// to online Raw Stream chip numbering and viceversa. 
// Used for Fast-Or studies.
// Implemented by valerio.altini@ba.infn.it

if(key>=kNModule*kNChip*(kNClockPhase+1))
  {AliDebug(1,"SwitchChipKeyNumbering: you asked for a non existing key"); return 99999;}
UInt_t mod=9999,chip=9999,phase=9999;
GetModAndChipFromKey(key,mod,chip);
if(mod<kNModuleLy1) chip = kNChip-(chip+1);
if(IsForFO(key))phase = GetBCm4FromKey(key);

return GetKey(mod,chip,IsForFO(key),phase);

}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetModFromKey(const UInt_t key) const {
  // get mod. from key
if(key>=kNModule*kNChip*(kNClockPhase+1))
  {AliError("GetModFromKey: you asked for a non existing key"); return 9999;}
return (key%(kNModule*kNChip))/kNChip;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetChipFromKey(const UInt_t key) const {
  // retrieves chip from key
if(key>=kNModule*kNChip*(kNClockPhase+1))
  {AliError("GetChipFromKey: you asked for a non existing key"); return 999;}
return ((key%(kNModule*kNChip))%(kNModule*kNChip))%kNChip;
}
//__________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetBCm4FromKey(const UInt_t key) const {
  // retrieves the "Bunch Crossing modulo 4" (for Fast Or studies)
if(key>=kNModule*kNChip*(kNClockPhase+1))
  {AliError("GetBCm4FromKey: you asked for a non existing key"); return 999;}
if(key<kNModule*kNChip) 
  {AliDebug(1,"GetBCm4FromKey: key is below 1200, why are you asking for FO related stuff"); return 999;}

return key/(kNModule*kNChip) - 1 ;
}
//__________________________________________________________________________
Bool_t AliITSPlaneEffSPD::IsForFO(const UInt_t key) const {
if(key>=kNModule*kNChip) return kTRUE;
else return kFALSE;
}
//__________________________________________________________________________
void AliITSPlaneEffSPD::GetModAndChipFromKey(const UInt_t key,UInt_t& mod,UInt_t& chip) const {
  // get module and chip from a key
if(key>=kNModule*kNChip*(kNClockPhase+1))
  {AliError("GetModAndChipFromKey: you asked for a non existing key"); 
  mod=9999;
  chip=999;
  return;}
mod=GetModFromKey(key);
chip=GetChipFromKey(key);
return;
}
//____________________________________________________________________________
Double_t AliITSPlaneEffSPD::LivePlaneEff(UInt_t key) const {
  // returns plane efficieny after adding the fraction of sensor which is bad
if(key>=kNModule*kNChip)
  {AliError("LivePlaneEff: you asked for a non existing key");
   return -1.;}
Double_t leff=AliITSPlaneEff::LivePlaneEff(0); // this just for the Warning
leff=PlaneEff(key)+GetFracBad(key);
return leff>1?1:leff;
}
//____________________________________________________________________________
Double_t AliITSPlaneEffSPD::ErrLivePlaneEff(UInt_t key) const {
  // returns error on live plane efficiency
if(key>=kNModule*kNChip)
  {AliError("ErrLivePlaneEff: you asked for a non existing key");
   return -1.;}
Int_t nf=fFound[key];
Double_t triedInLive=GetFracLive(key)*fTried[key];
Int_t nt=TMath::Max(nf,TMath::Nint(triedInLive));
return AliITSPlaneEff::ErrPlaneEff(nf,nt); // for the time being: to be checked
}
//_____________________________________________________________________________
Double_t AliITSPlaneEffSPD::GetFracLive(const UInt_t key) const {
  // returns the fraction of the sensor which is OK
if(key>=kNModule*kNChip)
  {AliError("GetFracLive: you asked for a non existing key");
   return -1.;}
    // Compute the fraction of bad (dead+noisy) detector 
UInt_t dead=0,noisy=0;
GetDeadAndNoisyInChip(key,dead,noisy);
Double_t live=dead+noisy;
live/=(kNRow*kNCol);
return 1.-live;
}
//_____________________________________________________________________________
void AliITSPlaneEffSPD::GetDeadAndNoisyInChip(const UInt_t key,
      UInt_t& nrDeadInChip, UInt_t& nrNoisyInChip) const {
  // returns the number of dead and noisy pixels
nrDeadInChip=0;
nrNoisyInChip=0;
if(key>=kNModule*kNChip)
  {AliError("GetDeadAndNoisyInChip: you asked for a non existing key");
   return;}
    // Compute the number of bad (dead+noisy) pixel in a chip
//
if(!fInitCDBCalled) 
  {AliError("GetDeadAndNoisyInChip: CDB not inizialized: call InitCDB first");
   return;};
AliCDBManager* man = AliCDBManager::Instance();
// retrieve map of dead Pixel 
AliCDBEntry *cdbSPDDead = man->Get("ITS/Calib/SPDDead", fRunNumber);
TObjArray* spdDead;
if(cdbSPDDead) {
  spdDead = (TObjArray*)cdbSPDDead->GetObject();
  if(!spdDead) 
  {AliError("GetDeadAndNoisyInChip: SPDDead not found in CDB");
   return;}
} else {
  AliError("GetDeadAndNoisyInChip: did not find Calib/SPDDead.");
  return;
}
// retrieve map of sparse dead Pixel 
AliCDBEntry *cdbSPDSparseDead = man->Get("ITS/Calib/SPDSparseDead", fRunNumber);
TObjArray* spdSparseDead;
if(cdbSPDSparseDead) {
  spdSparseDead = (TObjArray*)cdbSPDSparseDead->GetObject();
  if(!spdSparseDead) 
  {AliError("GetDeadAndNoisyInChip: SPDSparseDead not found in CDB");
   return;}
} else {
  AliError("GetDeadAndNoisyInChip: did not find Calib/SPDSparseDead.");
  return;
}

// retrieve map of noisy Pixel 
AliCDBEntry *cdbSPDNoisy = man->Get("ITS/Calib/SPDNoisy", fRunNumber);
TObjArray* spdNoisy;
if(cdbSPDNoisy) {
  spdNoisy = (TObjArray*)cdbSPDNoisy->GetObject();
  if(!spdNoisy) 
  {AliError("GetDeadAndNoisyInChip: SPDNoisy not found in CDB");
   return;}
} else {
  AliError("GetDeadAndNoisyInChip: did not find Calib/SPDNoisy.");
  return;
}
//
UInt_t mod=GetModFromKey(key);
UInt_t chip=GetChipFromKey(key);
// count number of dead
AliITSCalibrationSPD* calibSPD=(AliITSCalibrationSPD*) spdDead->At(mod);
UInt_t nrDead = calibSPD->GetNrBad();
for (UInt_t index=0; index<nrDead; index++) {
  if(GetChipFromCol(calibSPD->GetBadColAt(index))==chip) nrDeadInChip++;
}
// add the number of sparse dead to the previous dead
calibSPD=(AliITSCalibrationSPD*) spdSparseDead->At(mod);
UInt_t nrSparseDead = calibSPD->GetNrBad();
for (UInt_t index=0; index<nrSparseDead; index++) {
  if(GetChipFromCol(calibSPD->GetBadColAt(index))==chip) nrDeadInChip++;
}
calibSPD=(AliITSCalibrationSPD*) spdNoisy->At(mod);
UInt_t nrNoisy = calibSPD->GetNrBad();
for (UInt_t index=0; index<nrNoisy; index++) {
  if(GetChipFromCol(calibSPD->GetBadColAt(index))==chip) nrNoisyInChip++;
}
return;
}
//_____________________________________________________________________________
Double_t AliITSPlaneEffSPD::GetFracBad(const UInt_t key) const {
  // returns 1-fractional live
if(key>=kNModule*kNChip)
  {AliError("GetFracBad: you asked for a non existing key");
   return -1.;}
return 1.-GetFracLive(key);
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSPD::WriteIntoCDB() const {
// write onto CDB
if(!fInitCDBCalled)
  {AliError("WriteIntoCDB: CDB not inizialized. Call InitCDB first");
   return kFALSE;}
// to be written properly: now only for debugging 
  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  //md->SetObjectClassName("AliITSPlaneEff");
  md->SetResponsible("Giuseppe Eugenio Bruno");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion("head 19/11/07"); //root version
  AliCDBId id("ITS/PlaneEff/PlaneEffSPD",0,AliCDBRunRange::Infinity()); 
  AliITSPlaneEffSPD eff; 
  eff=*this;
  Bool_t r=AliCDBManager::Instance()->GetDefaultStorage()->Put(&eff,id,md);
  delete md;
  return r;
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSPD::ReadFromCDB() {
// read from CDB
if(!fInitCDBCalled)
  {AliError("ReadFromCDB: CDB not inizialized. Call InitCDB first");
   return kFALSE;}
AliCDBEntry *cdbEntry = AliCDBManager::Instance()->Get("ITS/PlaneEff/PlaneEffSPD",fRunNumber);
if(!cdbEntry) return kFALSE;
AliITSPlaneEffSPD* eff= (AliITSPlaneEffSPD*)cdbEntry->GetObject();
if(this==eff) return kFALSE;
if(fHis) CopyHistos(*eff); // If histos already exist then copy them to eff
eff->Copy(*this);          // copy everything (statistics and histos) from eff to this
return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliITSPlaneEffSPD::AddFromCDB(AliCDBId *cdbId) {
AliCDBEntry *cdbEntry=0;
if (!cdbId) {
  if(!fInitCDBCalled)  
    {AliError("ReadFromCDB: CDB not inizialized. Call InitCDB first"); return kFALSE;}
  cdbEntry = AliCDBManager::Instance()->Get("ITS/PlaneEff/PlaneEffSPD",fRunNumber);
} else {
  cdbEntry = AliCDBManager::Instance()->Get(*cdbId);
}
if(!cdbEntry) return kFALSE;
AliITSPlaneEffSPD* eff= (AliITSPlaneEffSPD*)cdbEntry->GetObject();
*this+=*eff;
return kTRUE;
}
//_____________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetKeyFromDetLocCoord(Int_t ilay, Int_t idet, 
						Float_t, Float_t locz) const {
// method to locate a basic block from Detector Local coordinate (to be used in tracking)
UInt_t key=999999;
if(ilay<0 || ilay>1) 
  {AliError("GetKeyFromDetLocCoord: you asked for a non existing layer");
   return key;}
if(ilay==0 && (idet<0 || idet>79))
 {AliError("GetKeyFromDetLocCoord: you asked for a non existing detector");
   return key;}
if(ilay==1 && (idet<0 || idet>159))
 {AliError("GetKeyFromDetLocCoord: you asked for a non existing detector");
   return key;}

UInt_t mod=idet;
if(ilay==1) mod+=80;
key=GetKey(mod,GetChipFromCol(GetColFromLocZ(locz)));
return key;
}
//_____________________________________________________________________________
UInt_t AliITSPlaneEffSPD::GetColFromLocZ(Float_t zloc) const {
// method to retrieve column number from the local z coordinate
  UInt_t col=0;
  AliITSsegmentationSPD spd;
  Int_t ix,iz;
  if(spd.LocalToDet(0,zloc,ix,iz)) col+=iz;
  else {
    AliDebug(1,Form("cannot compute column number from local z=%f",zloc));
    col=99999;}
  return col;
/*
const Float_t kconv = 1.0E-04; // converts microns to cm.
Float_t bz[160];
for(Int_t i=000;i<160;i++) bz[i] = 425.0; // most are 425 microns except below
bz[ 31] = bz[ 32] = 625.0; // first chip boundry
bz[ 63] = bz[ 64] = 625.0; // first chip boundry
bz[ 95] = bz[ 96] = 625.0; // first chip boundry
bz[127] = bz[128] = 625.0; // first chip boundry
//
Int_t j=-1;
Float_t dz=0;
for(Int_t i=000;i<160;i++) dz+=bz[i];
dz = -0.5*kconv*dz;
if(zloc<dz || zloc>-1*dz) { // outside z range
  AliDebug(1,Form("GetColFromLocZ: cannot compute column number from local z=%f",zloc));
  return 99999;}
for(j=0;j<160;j++){
  dz += kconv*bz[j];
  if(zloc<dz) break;
} // end for j
col+=j;
//
return col;
*/
}
//________________________________________________________
Bool_t AliITSPlaneEffSPD::GetBlockBoundaries(const UInt_t key, Float_t& xmn,Float_t& xmx,
                                             Float_t& zmn,Float_t& zmx) const {
//
//  This method return the geometrical boundaries of the active volume of a given 
//  basic block, in the detector reference system.
//  Input: unique key to locate a basic block.
//  
//  Output: Ymin, Ymax, Zmin, Zmax of a basic block (chip for SPD)
//  Return: kTRUE if computation was succesfully, kFALSE otherwise
//
if(key>=kNModule*kNChip)
  {AliDebug(1,"GetBlockBoundaries: you asked for a non existing key"); return kFALSE;}
UInt_t chip=GetChipFromKey(key);
zmn=GetLocZFromCol(chip*kNCol);
zmx=GetLocZFromCol((chip+1)*kNCol);
xmn=GetLocXFromRow(0);
xmx=GetLocXFromRow(kNRow);
//
Float_t tmp=zmn;
if(zmx<zmn) {zmn=zmx; zmx=tmp;}
tmp=xmn;
if(xmx<xmn) {xmn=xmx; xmx=tmp;}
return kTRUE;
}
//________________________________________________________
Float_t AliITSPlaneEffSPD::GetLocXFromRow(const UInt_t row) const {
// 
//  This method return the local (i.e. detector reference system) lower x coordinate 
//  of the row. To get the central value of a given row, you can do 
//  1/2*[LocXFromRow(row)+LocXFromRow(row+1)].
//
//  Input: row number in the range [0,kNRow] 
//  Output: lower local X coordinate of this row.
//
if(row>kNRow)  // not >= ! allow also computation of upper limit of the last row. 
  {AliError("LocYFromRow: you asked for a non existing row"); return 9999999.;}
// Use only AliITSsegmentationSPD
AliITSsegmentationSPD spd;
Double_t dummy,x;
if(row==kNRow) spd.CellBoundries((Int_t)row-1,0,dummy,x,dummy,dummy);
else spd.CellBoundries((Int_t)row,0,x,dummy,dummy,dummy);
return (Float_t)x;

}
//________________________________________________________
Float_t AliITSPlaneEffSPD::GetLocZFromCol(const UInt_t col) const {
//
//  This method return the local (i.e. detector reference system) lower Z coordinate
//  of the column. To get the central value of a given column, you can do
//  1/2*[LocZFromCol(col)+LocZFromCol(col+1)].
//
//  Input: col number in the range [0,kNChip*kNCol]
//  Output: lower local Y coordinate of this row.
//
if(col>kNChip*kNCol) // not >= ! allow also computation of upper limit of the last column
  {AliError("LocZFromCol: you asked for a non existing column"); return 9999999.;}
// Use only AliITSsegmentationSPD
AliITSsegmentationSPD spd;
Double_t dummy,y;
if(col==kNChip*kNCol) spd.CellBoundries(0,(Int_t)col-1,dummy,dummy,dummy,y);
else spd.CellBoundries(0,(Int_t)col,dummy,dummy,y,dummy);
return (Float_t)y;

}
//__________________________________________________________
void AliITSPlaneEffSPD::InitHistos() {
  // for the moment let's create the histograms 
  // module by  module
  TString histnameResX="HistResX_mod_",aux;
  TString histnameResZ="HistResZ_mod_";
  TString histnameResXZ="HistResXZ_mod_";
  TString histnameClusterType="HistClusterType_mod_";
  TString histnameResXclu="HistResX_mod_";
  TString histnameResZclu="HistResZ_mod_";
  TString histnameResXchip="HistResX_mod_";
  TString histnameResZchip="HistResZ_mod_";
  TString profnameResXvsPhi="ProfResXvsPhi_mod_";
  TString profnameResZvsDip="ProfResZvsDip_mod_";
  TString profnameResXvsPhiclu="ProfResXvsPhi_mod_";
  TString profnameResZvsDipclu="ProfResZvsDip_mod_";
  TString histnameTrackErrX="HistTrackErrX_mod_";
  TString histnameTrackErrZ="HistTrackErrZ_mod_";
  TString histnameClusErrX="HistClusErrX_mod_";
  TString histnameClusErrZ="HistClusErrZ_mod_";
  TString histnameTrackXFOtrue="HistTrackXFOok_mod_";
  TString histnameTrackZFOtrue="HistTrackZFOok_mod_";
  TString histnameTrackXFOfalse="HistTrackXFOko_mod_";
  TString histnameTrackZFOfalse="HistTrackZFOko_mod_";
  TString histnameTrackXZFOtrue="HistTrackZvsXFOok_mod_";
  TString histnameTrackXZFOfalse="HistTrackZvsXFOko_mod_";
//

  TH1::AddDirectory(kFALSE);

  fHisResX=new TH1F*[kNHisto];
  fHisResZ=new TH1F*[kNHisto];
  fHisResXZ=new TH2F*[kNHisto];
  fHisClusterSize=new TH2I*[kNHisto];
  fHisResXclu=new TH1F**[kNHisto];
  fHisResZclu=new TH1F**[kNHisto];
  fHisResXchip=new TH1F**[kNHisto];
  fHisResZchip=new TH1F**[kNHisto];
  fProfResXvsPhi=new TProfile*[kNHisto];
  fProfResZvsDip=new TProfile*[kNHisto];
  fProfResXvsPhiclu=new TProfile**[kNHisto];
  fProfResZvsDipclu=new TProfile**[kNHisto];
  fHisTrackErrX=new TH1F*[kNHisto];
  fHisTrackErrZ=new TH1F*[kNHisto];
  fHisClusErrX=new TH1F*[kNHisto];
  fHisClusErrZ=new TH1F*[kNHisto];
  fHisTrackXFOtrue=new TH1F**[kNHisto];
  fHisTrackZFOtrue=new TH1F**[kNHisto];
  fHisTrackXFOfalse=new TH1F**[kNHisto];
  fHisTrackZFOfalse=new TH1F**[kNHisto];
  fHisTrackXZFOtrue=new TH2F**[kNHisto];
  fHisTrackXZFOfalse=new TH2F**[kNHisto];

  for (Int_t nhist=0;nhist<kNHisto;nhist++){
    aux=histnameResX;
    aux+=nhist;
    fHisResX[nhist]=new TH1F("histname","histname",1600,-0.32,0.32); // +- 3200 micron; 1 bin=4 micron
    fHisResX[nhist]->SetName(aux.Data());
    fHisResX[nhist]->SetTitle(aux.Data());

    aux=histnameResZ;
    aux+=nhist;
    fHisResZ[nhist]=new TH1F("histname","histname",1200,-0.48,0.48); // +-4800 micron; 1 bin=8 micron
    fHisResZ[nhist]->SetName(aux.Data());
    fHisResZ[nhist]->SetTitle(aux.Data());

    aux=histnameResXZ;
    aux+=nhist;
    fHisResXZ[nhist]=new TH2F("histname","histname",80,-0.16,0.16,80,-0.32,0.32); // binning:
    fHisResXZ[nhist]->SetName(aux.Data());					   // 40 micron in x;
    fHisResXZ[nhist]->SetTitle(aux.Data());					   // 80 micron in z;

    aux=histnameClusterType;
    aux+=nhist;
    fHisClusterSize[nhist]=new TH2I("histname","histname",10,0.5,10.5,10,0.5,10.5);
    fHisClusterSize[nhist]->SetName(aux.Data());
    fHisClusterSize[nhist]->SetTitle(aux.Data());

    fHisResXclu[nhist]=new TH1F*[kNclu];
    fHisResZclu[nhist]=new TH1F*[kNclu];
    fHisTrackXFOtrue[nhist]=new TH1F*[kNClockPhase];
    fHisTrackZFOtrue[nhist]=new TH1F*[kNClockPhase];
    fHisTrackXFOfalse[nhist]=new TH1F*[kNClockPhase];
    fHisTrackZFOfalse[nhist]=new TH1F*[kNClockPhase];
    fHisTrackXZFOtrue[nhist]=new TH2F*[kNClockPhase];
    fHisTrackXZFOfalse[nhist]=new TH2F*[kNClockPhase];

    for(Int_t clu=0; clu<kNclu; clu++) {  // clu=0 --> cluster size 1
      aux=histnameResXclu;
      aux+=nhist;
      aux+="_clu_";
      aux+=clu+1; // clu=0 --> cluster size 1
      fHisResXclu[nhist][clu]=new TH1F("histname","histname",1600,-0.32,0.32); // +- 3200 micron; 1 bin=4 micron
      fHisResXclu[nhist][clu]->SetName(aux.Data());
      fHisResXclu[nhist][clu]->SetTitle(aux.Data());

      aux=histnameResZclu;
      aux+=nhist;
      aux+="_clu_";
      aux+=clu+1; // clu=0 --> cluster size 1
      fHisResZclu[nhist][clu]=new TH1F("histname","histname",1200,-0.48,0.48); // +-4800 micron; 1 bin=8 micron
      fHisResZclu[nhist][clu]->SetName(aux.Data());
      fHisResZclu[nhist][clu]->SetTitle(aux.Data());
    }

    fHisResXchip[nhist]=new TH1F*[kNChip];
    fHisResZchip[nhist]=new TH1F*[kNChip];
    for(Int_t chip=0; chip<kNChip; chip++) { 
      aux=histnameResXchip;
      aux+=nhist;
      aux+="_chip_";
      aux+=chip; 
      fHisResXchip[nhist][chip]=new TH1F("histname","histname",800,-0.32,0.32); // +- 3200 micron; 1 bin=8 micron
      fHisResXchip[nhist][chip]->SetName(aux.Data());
      fHisResXchip[nhist][chip]->SetTitle(aux.Data());

      aux=histnameResZchip;
      aux+=nhist;
      aux+="_chip_";
      aux+=chip;
      fHisResZchip[nhist][chip]=new TH1F("histname","histname",300,-0.48,0.48); // +-4800 micron; 1 bin=32 micron
      fHisResZchip[nhist][chip]->SetName(aux.Data());
      fHisResZchip[nhist][chip]->SetTitle(aux.Data());
    }

    aux=histnameTrackErrX;
    aux+=nhist;
    fHisTrackErrX[nhist]=new TH1F("histname","histname",400,0.,0.32); // 0-3200 micron; 1 bin=8 micron
    fHisTrackErrX[nhist]->SetName(aux.Data());
    fHisTrackErrX[nhist]->SetTitle(aux.Data());

    aux=histnameTrackErrZ;
    aux+=nhist;
    fHisTrackErrZ[nhist]=new TH1F("histname","histname",200,0.,0.32); // 0-3200 micron; 1 bin=16 micron
    fHisTrackErrZ[nhist]->SetName(aux.Data());
    fHisTrackErrZ[nhist]->SetTitle(aux.Data());

    aux=histnameClusErrX;
    aux+=nhist;
    fHisClusErrX[nhist]=new TH1F("histname","histname",400,0.,0.08); //  0-800 micron; 1 bin=2 micron
    fHisClusErrX[nhist]->SetName(aux.Data());
    fHisClusErrX[nhist]->SetTitle(aux.Data());

    aux=histnameClusErrZ;
    aux+=nhist;
    fHisClusErrZ[nhist]=new TH1F("histname","histname",400,0.,0.32); //  0-3200 micron; 1 bin=8 micron
    fHisClusErrZ[nhist]->SetName(aux.Data());
    fHisClusErrZ[nhist]->SetTitle(aux.Data());

    aux=profnameResXvsPhi;
    aux+=nhist;
    fProfResXvsPhi[nhist]=new TProfile("histname","histname",40,-40.,40.0); // binning: range:  -40°- 40°
    fProfResXvsPhi[nhist]->SetName(aux.Data());                             //          bin width: 2°
    fProfResXvsPhi[nhist]->SetTitle(aux.Data());

    aux=profnameResZvsDip;
    aux+=nhist;
    fProfResZvsDip[nhist]=new TProfile("histname","histname",48,-72.,72.0); // binning: range:  -70°-4°
    fProfResZvsDip[nhist]->SetName(aux.Data());                             //          bin width: 3°
    fProfResZvsDip[nhist]->SetTitle(aux.Data());

    fProfResXvsPhiclu[nhist]=new TProfile*[kNclu];
    fProfResZvsDipclu[nhist]=new TProfile*[kNclu];
    for(Int_t clu=0; clu<kNclu; clu++) {  // clu=0 --> cluster size 1
      aux=profnameResXvsPhiclu;
      aux+=nhist;
      aux+="_clu_";
      aux+=clu+1; // clu=0 --> cluster size 1
      fProfResXvsPhiclu[nhist][clu]=new TProfile("histname","histname",40,-40.,40.0); // binning: range:  -40°- 40
      fProfResXvsPhiclu[nhist][clu]->SetName(aux.Data());                             //          bin width: 2°
      fProfResXvsPhiclu[nhist][clu]->SetTitle(aux.Data());

      aux=profnameResZvsDipclu;
      aux+=nhist;
      aux+="_clu_";
      aux+=clu+1; // clu=0 --> cluster size 1
      fProfResZvsDipclu[nhist][clu]= new TProfile("histname","histname",48,-72.,72.0); // binning: range:  -70°-7°
      fProfResZvsDipclu[nhist][clu]->SetName(aux.Data());                              //      bin width: 3°
      fProfResZvsDipclu[nhist][clu]->SetTitle(aux.Data());
    }

    fHisTrackXFOtrue[nhist]=new TH1F*[kNClockPhase];
    fHisTrackZFOtrue[nhist]=new TH1F*[kNClockPhase];
    fHisTrackXFOfalse[nhist]=new TH1F*[kNClockPhase];
    fHisTrackZFOfalse[nhist]=new TH1F*[kNClockPhase];
    fHisTrackXZFOtrue[nhist]=new TH2F*[kNClockPhase];
    fHisTrackXZFOfalse[nhist]=new TH2F*[kNClockPhase];
    for(Int_t phas=0; phas<kNClockPhase;phas++){
      aux=histnameTrackXFOtrue;
      aux+=nhist;
      aux+="_BCmod4_";
      aux+=phas;
      fHisTrackXFOtrue[nhist][phas]=new TH1F("histname","histname",128,-0.64,0.64); // +- 6.4 mm; 1 bin=0.1 mm
      fHisTrackXFOtrue[nhist][phas]->SetName(aux.Data());
      fHisTrackXFOtrue[nhist][phas]->SetTitle(aux.Data());
      
      aux=histnameTrackZFOtrue;
      aux+=nhist;
      aux+="_BCmod4_";
      aux+=phas;
      fHisTrackZFOtrue[nhist][phas]=new TH1F("histname","histname",350,-3.5,3.5); // +- 35. mm; 1 bin=0.2 mm
      fHisTrackZFOtrue[nhist][phas]->SetName(aux.Data());
      fHisTrackZFOtrue[nhist][phas]->SetTitle(aux.Data());
      
      aux=histnameTrackXFOfalse;
      aux+=nhist;
      aux+="_BCmod4_";
      aux+=phas;
      fHisTrackXFOfalse[nhist][phas]=new TH1F("histname","histname",128,-0.64,0.64); // +- 6.4 mm; 1 bin=0.1 mm
      fHisTrackXFOfalse[nhist][phas]->SetName(aux.Data());
      fHisTrackXFOfalse[nhist][phas]->SetTitle(aux.Data());

      aux=histnameTrackZFOfalse;
      aux+=nhist;
      aux+="_BCmod4_";
      aux+=phas;
      fHisTrackZFOfalse[nhist][phas]=new TH1F("histname","histname",350,-3.5,3.5); // +- 35. mm; 1 bin=0.2 mm
      fHisTrackZFOfalse[nhist][phas]->SetName(aux.Data());
      fHisTrackZFOfalse[nhist][phas]->SetTitle(aux.Data());
    
      aux=histnameTrackXZFOtrue;
      aux+=nhist;
      aux+="_BCmod4_";
      aux+=phas;
      fHisTrackXZFOtrue[nhist][phas]=new TH2F("histname","histname",22,-3.5,3.5,32,-0.64,0.64); //  localZ +- 35. mm; 1 bin=3.2 mm
      fHisTrackXZFOtrue[nhist][phas]->SetName(aux.Data());                                      //  localX +- 6.4 mm; 1 bin=0.4 mm
      fHisTrackXZFOtrue[nhist][phas]->SetTitle(aux.Data());

      aux=histnameTrackXZFOfalse;
      aux+=nhist;
      aux+="_BCmod4_";
      aux+=phas;
      fHisTrackXZFOfalse[nhist][phas]=new TH2F("histname","histname",22,-3.5,3.5,32,-0.64,0.64); //  localZ +- 35. mm; 1 bin=3.2 mm
      fHisTrackXZFOfalse[nhist][phas]->SetName(aux.Data());                                      //  localX +- 6.4 mm; 1 bin=0.4 mm
      fHisTrackXZFOfalse[nhist][phas]->SetTitle(aux.Data());
      } 
  } // end loop on module

  TH1::AddDirectory(kTRUE);

return;
}
//__________________________________________________________
void AliITSPlaneEffSPD::DeleteHistos() {
  if(fHisResX) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fHisResX[i];
    delete [] fHisResX; fHisResX=0;
  }
  if(fHisResZ) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fHisResZ[i];
    delete [] fHisResZ; fHisResZ=0;
  }
  if(fHisResXZ) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fHisResXZ[i];
    delete [] fHisResXZ; fHisResXZ=0;
  }
  if(fHisClusterSize) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fHisClusterSize[i];
    delete [] fHisClusterSize; fHisClusterSize=0;
  }
  if(fHisResXclu) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t clu=0; clu<kNclu; clu++) if (fHisResXclu[i][clu]) delete fHisResXclu[i][clu];
      delete [] fHisResXclu[i];
    }
    delete [] fHisResXclu;
    fHisResXclu = 0;
  }
  if(fHisResZclu) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t clu=0; clu<kNclu; clu++) if (fHisResZclu[i][clu]) delete fHisResZclu[i][clu];
      delete [] fHisResZclu[i];
    }
    delete [] fHisResZclu;
    fHisResZclu = 0;
  }
  if(fHisResXchip) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t chip=0; chip<kNChip; chip++) if (fHisResXchip[i][chip]) delete fHisResXchip[i][chip];
      delete [] fHisResXchip[i];
    }
    delete [] fHisResXchip;
    fHisResXchip = 0;
  }
  if(fHisResZchip) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t chip=0; chip<kNChip; chip++) if (fHisResZchip[i][chip]) delete fHisResZchip[i][chip];
      delete [] fHisResZchip[i];
    }
    delete [] fHisResZchip;
    fHisResZchip = 0;
  }
  if(fHisTrackErrX) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fHisTrackErrX[i];
    delete [] fHisTrackErrX; fHisTrackErrX=0;
  }
  if(fHisTrackErrZ) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fHisTrackErrZ[i];
    delete [] fHisTrackErrZ; fHisTrackErrZ=0;
  }
  if(fHisClusErrX) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fHisClusErrX[i];
    delete [] fHisClusErrX; fHisClusErrX=0;
  }
  if(fHisClusErrZ) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fHisClusErrZ[i];
    delete [] fHisClusErrZ; fHisClusErrZ=0;
  }
  if(fProfResXvsPhi) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fProfResXvsPhi[i];
    delete [] fProfResXvsPhi; fProfResXvsPhi=0;
  }
  if(fProfResZvsDip) {
    for (Int_t i=0; i<kNHisto; i++ ) delete fProfResZvsDip[i];
    delete [] fProfResZvsDip; fProfResZvsDip=0;
  }
  if(fProfResXvsPhiclu) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t clu=0; clu<kNclu; clu++) if (fProfResXvsPhiclu[i][clu]) delete fProfResXvsPhiclu[i][clu];
      delete [] fProfResXvsPhiclu[i];
    }
    delete [] fProfResXvsPhiclu;
    fProfResXvsPhiclu = 0;
  }
  if(fProfResZvsDipclu) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t clu=0; clu<kNclu; clu++) if (fProfResZvsDipclu[i][clu]) delete fProfResZvsDipclu[i][clu];
      delete [] fProfResZvsDipclu[i];
    }
    delete [] fProfResZvsDipclu;
    fProfResZvsDipclu = 0;
  }
  if(fHisTrackXFOtrue) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t phas=0; phas<kNClockPhase; phas++) if (fHisTrackXFOtrue[i][phas]) delete fHisTrackXFOtrue[i][phas];
      delete [] fHisTrackXFOtrue[i];
    }
    delete [] fHisTrackXFOtrue;
    fHisTrackXFOtrue = 0;
  }
  if(fHisTrackZFOtrue) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t phas=0; phas<kNClockPhase; phas++) if (fHisTrackZFOtrue[i][phas]) delete fHisTrackZFOtrue[i][phas];
      delete [] fHisTrackZFOtrue[i];
    }
    delete [] fHisTrackZFOtrue;
    fHisTrackZFOtrue = 0;
  }
  if(fHisTrackXFOfalse) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t phas=0; phas<kNClockPhase; phas++) if (fHisTrackXFOfalse[i][phas]) delete fHisTrackXFOfalse[i][phas];
      delete [] fHisTrackXFOfalse[i];
    }
    delete [] fHisTrackXFOfalse;
    fHisTrackXFOfalse = 0;
  }
  if(fHisTrackZFOfalse) {
    for (Int_t i=0; i<kNHisto; i++ ) {
      for (Int_t phas=0; phas<kNClockPhase; phas++) if (fHisTrackZFOfalse[i][phas]) delete fHisTrackZFOfalse[i][phas];
      delete [] fHisTrackZFOfalse[i];
    }
    delete [] fHisTrackZFOfalse;
    fHisTrackZFOfalse = 0;
  }
return;
}
//__________________________________________________________
Bool_t AliITSPlaneEffSPD::FillHistos(UInt_t key, Bool_t found,
                                     Float_t *tr, Float_t *clu, Int_t *csize, Float_t *angtrkmod) {
//
// depending on the value of key this method
// either call the standard one for clusters 
// or the one for FO studies
// if key <  1200 --> call FillHistosST
// if key >= 1200 --> call FillHistosFO
if(key>=kNModule*kNChip*(kNClockPhase+1))
  {AliError("GetChipFromKey: you asked for a non existing key"); return kFALSE;}
if(key<kNModule*kNChip) return FillHistosStd(key,found,tr,clu,csize,angtrkmod);
else return FillHistosFO(key,found,tr);
return kFALSE;
}
//__________________________________________________________
Bool_t AliITSPlaneEffSPD::FillHistosFO(UInt_t key, Bool_t found, Float_t *tr) {
// this method fill the histograms for FastOr studies
// input: - key: unique key of the basic block 
//        - found: Boolean to asses whether a FastOr bit has been associated to the track or not 
//        - tr[0],tr[1] local X and Z coordinates of the track prediction, respectively
//        - tr[2],tr[3] error on local X and Z coordinates of the track prediction, respectively
// output: kTRUE if filling was succesfull kFALSE otherwise
// side effects: updating of the histograms.
  if (!fHis) {
    AliWarning("FillHistos: histograms do not exist! Call SetCreateHistos(kTRUE) first");
    return kFALSE;
  }
  if(key>=kNModule*kNChip*(kNClockPhase+1))
    {AliWarning("FillHistos: you asked for a non existing key"); return kFALSE;}
  if(key<kNModule*kNChip)
    {AliWarning("FillHistos: you asked for a key which is not for FO studies"); return kFALSE;}
  Int_t id=GetModFromKey(key);
  Int_t BCm4=GetBCm4FromKey(key);
  if(id>=kNHisto)
    {AliWarning("FillHistos: you want to fill a non-existing histos"); return kFALSE;}
  if(found) {
    fHisTrackXFOtrue[id][BCm4]->Fill(tr[0]);
    fHisTrackZFOtrue[id][BCm4]->Fill(tr[1]);
    fHisTrackXZFOtrue[id][BCm4]->Fill(tr[1],tr[0]);
  }
  else {
    fHisTrackXFOfalse[id][BCm4]->Fill(tr[0]);
    fHisTrackZFOfalse[id][BCm4]->Fill(tr[1]);
    fHisTrackXZFOfalse[id][BCm4]->Fill(tr[1],tr[0]);
  }
return kTRUE;
}
//__________________________________________________________
Bool_t AliITSPlaneEffSPD::FillHistosStd(UInt_t key, Bool_t found, 
                                     Float_t *tr, Float_t *clu, Int_t *csize, Float_t *angtrkmod) {
// this method fill the histograms
// input: - key: unique key of the basic block 
//        - found: Boolean to asses whether a cluster has been associated to the track or not 
//        - tr[0],tr[1] local X and Z coordinates of the track prediction, respectively
//        - tr[2],tr[3] error on local X and Z coordinates of the track prediction, respectively
//        - clu[0],clu[1] local X and Z coordinates of the cluster associated to the track, respectively
//        - clu[2],clu[3] error on local X and Z coordinates of the cluster associated to the track, respectively
//        - csize[0][1] cluster size in X and Z, respectively
//        - angtrkmod[0],angtrkmod[1]  
// output: kTRUE if filling was succesfull kFALSE otherwise
// side effects: updating of the histograms. 
//
  if (!fHis) {
    AliWarning("FillHistos: histograms do not exist! Call SetCreateHistos(kTRUE) first");
    return kFALSE;
  }
  if(key>=kNModule*kNChip)
    {AliWarning("FillHistos: you asked for a non existing key"); return kFALSE;}
  Int_t id=GetModFromKey(key);
  Int_t chip=GetChipFromKey(key);
  if(id>=kNHisto) 
    {AliWarning("FillHistos: you want to fill a non-existing histos"); return kFALSE;}
  if(found) {
    Float_t resx=tr[0]-clu[0];
    Float_t resz=tr[1]-clu[1];
    fHisResX[id]->Fill(resx);
    fHisResZ[id]->Fill(resz);
    fHisResXZ[id]->Fill(resx,resz);
    fHisClusterSize[id]->Fill((Double_t)csize[0],(Double_t)csize[1]);
    if(csize[0]>0 &&  csize[0]<=kNclu) fHisResXclu[id][csize[0]-1]->Fill(resx);
    if(csize[1]>0 &&  csize[1]<=kNclu) fHisResZclu[id][csize[1]-1]->Fill(resz);
    fHisResXchip[id][chip]->Fill(resx);
    fHisResZchip[id][chip]->Fill(resz);
    fProfResXvsPhi[id]->Fill(angtrkmod[0],resx);
    fProfResZvsDip[id]->Fill(angtrkmod[1],resz);
    if(csize[0]>0 &&  csize[0]<=kNclu) fProfResXvsPhiclu[id][csize[0]-1]->Fill(angtrkmod[0],resx);
    if(csize[1]>0 &&  csize[1]<=kNclu) fProfResZvsDipclu[id][csize[1]-1]->Fill(angtrkmod[1],resz);
  }
  fHisTrackErrX[id]->Fill(tr[2]);
  fHisTrackErrZ[id]->Fill(tr[3]);
  fHisClusErrX[id]->Fill(clu[2]);
  fHisClusErrZ[id]->Fill(clu[3]);
  return kTRUE;
}
//__________________________________________________________
Bool_t AliITSPlaneEffSPD::WriteHistosToFile(TString filename, Option_t* option) {
  //
  // Saves the histograms into a tree and saves the trees into a file
  //
  if (!fHis) return kFALSE;
  if (filename.IsNull() || filename.IsWhitespace()) {
     AliWarning("WriteHistosToFile: null output filename!");
     return kFALSE;
  }
  char branchname[51];
  TFile *hFile=new TFile(filename.Data(),option,
                         "The File containing the TREEs with ITS PlaneEff Histos");
  TTree *SPDTree=new TTree("SPDTree","Tree whith Residuals and Cluster Type distributions for SPD");
  TH1F *histZ,*histX;
  TH2F *histXZ;
  TH2I *histClusterType;
  TH1F *histXclu[kNclu];
  TH1F *histZclu[kNclu];
  TH1F *histXchip[kNChip];
  TH1F *histZchip[kNChip];
  TH1F *histTrErrZ,*histTrErrX;
  TH1F *histClErrZ,*histClErrX;
  TProfile *profXvsPhi,*profZvsDip;
  TProfile *profXvsPhiclu[kNclu],*profZvsDipclu[kNclu];
  TH1F *histXtrkFOtrue[kNClockPhase];
  TH1F *histZtrkFOtrue[kNClockPhase];
  TH1F *histXtrkFOfalse[kNClockPhase];
  TH1F *histZtrkFOfalse[kNClockPhase];
  TH2F *histXZtrkFOtrue[kNClockPhase];
  TH2F *histXZtrkFOfalse[kNClockPhase];

  histZ=new TH1F();
  histX=new TH1F();
  histXZ=new TH2F();
  histClusterType=new TH2I();
  for(Int_t clu=0;clu<kNclu;clu++) {
    histXclu[clu]=new TH1F();
    histZclu[clu]=new TH1F();
  }
  for(Int_t chip=0;chip<kNChip;chip++) {
    histXchip[chip]=new TH1F();
    histZchip[chip]=new TH1F();
  }

  histTrErrX=new TH1F();
  histTrErrZ=new TH1F();
  histClErrX=new TH1F();
  histClErrZ=new TH1F();
  profXvsPhi=new TProfile();
  profZvsDip=new TProfile();
  for(Int_t clu=0;clu<kNclu;clu++) {
    profXvsPhiclu[clu]=new TProfile();
    profZvsDipclu[clu]=new TProfile();
  }

  for(Int_t phas=0; phas<kNClockPhase;phas++){
    histXtrkFOtrue[phas]=new TH1F();
    histZtrkFOtrue[phas]=new TH1F();
    histXtrkFOfalse[phas]=new TH1F();
    histZtrkFOfalse[phas]=new TH1F();
    histXZtrkFOtrue[phas]=new TH2F();
    histXZtrkFOfalse[phas]=new TH2F();
  }

  SPDTree->Branch("histX","TH1F",&histX,128000,0);
  SPDTree->Branch("histZ","TH1F",&histZ,128000,0);
  SPDTree->Branch("histXZ","TH2F",&histXZ,128000,0);
  SPDTree->Branch("histClusterType","TH2I",&histClusterType,128000,0);
  for(Int_t clu=0;clu<kNclu;clu++) {
    snprintf(branchname,50,"histXclu_%d",clu+1);
    SPDTree->Branch(branchname,"TH1F",&histXclu[clu],128000,0);
    snprintf(branchname,50,"histZclu_%d",clu+1);
    SPDTree->Branch(branchname,"TH1F",&histZclu[clu],128000,0);
  }
  for(Int_t chip=0;chip<kNChip;chip++) {
    snprintf(branchname,50,"histXchip_%d",chip);
    SPDTree->Branch(branchname,"TH1F",&histXchip[chip],128000,0);
    snprintf(branchname,50,"histZchip_%d",chip);
    SPDTree->Branch(branchname,"TH1F",&histZchip[chip],128000,0);
  }
  SPDTree->Branch("histTrErrX","TH1F",&histTrErrX,128000,0);
  SPDTree->Branch("histTrErrZ","TH1F",&histTrErrZ,128000,0);
  SPDTree->Branch("histClErrX","TH1F",&histClErrX,128000,0);
  SPDTree->Branch("histClErrZ","TH1F",&histClErrZ,128000,0);
  SPDTree->Branch("profXvsPhi","TProfile",&profXvsPhi,128000,0);
  SPDTree->Branch("profZvsDip","TProfile",&profZvsDip,128000,0);
  for(Int_t clu=0;clu<kNclu;clu++) {
    snprintf(branchname,50,"profXvsPhiclu_%d",clu+1);
    SPDTree->Branch(branchname,"TProfile",&profXvsPhiclu[clu],128000,0);
    snprintf(branchname,50,"profZvsDipclu_%d",clu+1);
    SPDTree->Branch(branchname,"TProfile",&profZvsDipclu[clu],128000,0);
  }
  for(Int_t phas=0; phas<kNClockPhase;phas++){
    snprintf(branchname,50,"histTrXFOokBCmod4_%d",phas);
    SPDTree->Branch(branchname,"TH1F",&histXtrkFOtrue[phas],128000,0);
    snprintf(branchname,50,"histTrZFOokBCmod4_%d",phas);
    SPDTree->Branch(branchname,"TH1F",&histZtrkFOtrue[phas],128000,0);
    snprintf(branchname,50,"histTrXFOkoBCmod4_%d",phas);
    SPDTree->Branch(branchname,"TH1F",&histXtrkFOfalse[phas],128000,0);
    snprintf(branchname,50,"histTrZFOkoBCmod4_%d",phas);
    SPDTree->Branch(branchname,"TH1F",&histZtrkFOfalse[phas],128000,0);
    snprintf(branchname,50,"histTrXZFOokBCmod4_%d",phas);
    SPDTree->Branch(branchname,"TH2F",&histXZtrkFOtrue[phas],128000,0);
    snprintf(branchname,50,"histTrXZFOkoBCmod4_%d",phas);
    SPDTree->Branch(branchname,"TH2F",&histXZtrkFOfalse[phas],128000,0);
  }

  for(Int_t j=0;j<kNHisto;j++){
    histX=fHisResX[j];
    histZ=fHisResZ[j];
    histXZ=fHisResXZ[j];
    histClusterType=fHisClusterSize[j];
    for(Int_t clu=0;clu<kNclu;clu++) {
      histXclu[clu]=fHisResXclu[j][clu];
      histZclu[clu]=fHisResZclu[j][clu];
    }
    for(Int_t chip=0;chip<kNChip;chip++) {
      histXchip[chip]=fHisResXchip[j][chip];
      histZchip[chip]=fHisResZchip[j][chip];
    }
    histTrErrX=fHisTrackErrX[j];
    histTrErrZ=fHisTrackErrZ[j];
    histClErrX=fHisClusErrX[j];
    histClErrZ=fHisClusErrZ[j];
    profXvsPhi=fProfResXvsPhi[j];
    profZvsDip=fProfResZvsDip[j];
    for(Int_t clu=0;clu<kNclu;clu++) {
      profXvsPhiclu[clu]=fProfResXvsPhiclu[j][clu];
      profZvsDipclu[clu]=fProfResZvsDipclu[j][clu];
    }
    for(Int_t phas=0; phas<kNClockPhase;phas++){
      histXtrkFOtrue[phas]=fHisTrackXFOtrue[j][phas];
      histZtrkFOtrue[phas]=fHisTrackZFOtrue[j][phas];
      histXtrkFOfalse[phas]=fHisTrackXFOfalse[j][phas];
      histZtrkFOfalse[phas]=fHisTrackZFOfalse[j][phas];
      histXZtrkFOtrue[phas]=fHisTrackXZFOtrue[j][phas];
      histXZtrkFOfalse[phas]=fHisTrackXZFOfalse[j][phas];
    }

    SPDTree->Fill();
  }
  hFile->Write();
  hFile->Close();
return kTRUE;
}
//__________________________________________________________
Bool_t AliITSPlaneEffSPD::ReadHistosFromFile(TString filename) {
  //
  // Read histograms from an already existing file 
  //
  if (!fHis) return kFALSE;
  if (filename.IsNull() || filename.IsWhitespace()) {
     AliWarning("ReadHistosFromFile: incorrect output filename!");
     return kFALSE;
  }
  char branchname[51];

  TH1F *h  = 0;
  TH2F *h2 = 0;
  TH2I *h2i= 0;
  TProfile *p = 0;

  TFile *file=TFile::Open(filename.Data(),"READONLY");

  if (!file || file->IsZombie()) {
    AliWarning(Form("Can't open %s !",filename.Data()));
    delete file;
    return kFALSE;
  }
  TTree *tree = (TTree*) file->Get("SPDTree");

  TBranch *histX = (TBranch*) tree->GetBranch("histX");
  TBranch *histZ = (TBranch*) tree->GetBranch("histZ");
  TBranch *histXZ = (TBranch*) tree->GetBranch("histXZ");
  TBranch *histClusterType = (TBranch*) tree->GetBranch("histClusterType");
   
  TBranch *histXclu[kNclu], *histZclu[kNclu];
  for(Int_t clu=0; clu<kNclu; clu++) {
    snprintf(branchname,50,"histXclu_%d",clu+1);
    histXclu[clu]= (TBranch*) tree->GetBranch(branchname);
    snprintf(branchname,50,"histZclu_%d",clu+1);
    histZclu[clu]= (TBranch*) tree->GetBranch(branchname);
  }

  TBranch *histXchip[kNChip], *histZchip[kNChip];
  for(Int_t chip=0; chip<kNChip; chip++) {
    snprintf(branchname,50,"histXchip_%d",chip);
    histXchip[chip]= (TBranch*) tree->GetBranch(branchname);
    snprintf(branchname,50,"histZchip_%d",chip);
    histZchip[chip]= (TBranch*) tree->GetBranch(branchname);
  }

  TBranch *histTrErrX = (TBranch*) tree->GetBranch("histTrErrX");
  TBranch *histTrErrZ = (TBranch*) tree->GetBranch("histTrErrZ");
  TBranch *histClErrX = (TBranch*) tree->GetBranch("histClErrX");
  TBranch *histClErrZ = (TBranch*) tree->GetBranch("histClErrZ");
  TBranch *profXvsPhi = (TBranch*) tree->GetBranch("profXvsPhi");
  TBranch *profZvsDip = (TBranch*) tree->GetBranch("profZvsDip");

  TBranch *profXvsPhiclu[kNclu], *profZvsDipclu[kNclu];
  for(Int_t clu=0; clu<kNclu; clu++) {
    snprintf(branchname,50,"profXvsPhiclu_%d",clu+1);
    profXvsPhiclu[clu]= (TBranch*) tree->GetBranch(branchname);
    snprintf(branchname,50,"profZvsDipclu_%d",clu+1);
    profZvsDipclu[clu]= (TBranch*) tree->GetBranch(branchname);
  }

  TBranch *histXtrkFOtrue[kNClockPhase], *histZtrkFOtrue[kNClockPhase],
          *histXtrkFOfalse[kNClockPhase], *histZtrkFOfalse[kNClockPhase],
          *histXZtrkFOtrue[kNClockPhase], *histXZtrkFOfalse[kNClockPhase];
  for(Int_t phas=0; phas<kNClockPhase;phas++){
    snprintf(branchname,50,"histTrXFOokBCmod4_%d",phas);
    histXtrkFOtrue[phas] = (TBranch*) tree->GetBranch(branchname);
    snprintf(branchname,50,"histTrZFOokBCmod4_%d",phas);
    histZtrkFOtrue[phas] = (TBranch*) tree->GetBranch(branchname);
    snprintf(branchname,50,"histTrXFOkoBCmod4_%d",phas);
    histXtrkFOfalse[phas] = (TBranch*) tree->GetBranch(branchname);
    snprintf(branchname,50,"histTrZFOkoBCmod4_%d",phas);
    histZtrkFOfalse[phas] = (TBranch*) tree->GetBranch(branchname);
    snprintf(branchname,50,"histTrXZFOokBCmod4_%d",phas);
    histXZtrkFOtrue[phas] = (TBranch*) tree->GetBranch(branchname);
    snprintf(branchname,50,"histTrXZFOkoBCmod4_%d",phas);
    histXZtrkFOfalse[phas] = (TBranch*) tree->GetBranch(branchname);
  }

  gROOT->cd();

  Int_t nevent = (Int_t)histX->GetEntries();
  if(nevent!=kNHisto) 
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  histX->SetAddress(&h);
  for(Int_t j=0;j<kNHisto;j++){
    histX->GetEntry(j);
    fHisResX[j]->Add(h);
  }

  nevent = (Int_t)histZ->GetEntries();
  if(nevent!=kNHisto) 
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  histZ->SetAddress(&h);
  for(Int_t j=0;j<kNHisto;j++){
    histZ->GetEntry(j);
    fHisResZ[j]->Add(h);
  }

  nevent = (Int_t)histXZ->GetEntries();
  if(nevent!=kNHisto) 
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  histXZ->SetAddress(&h2);
  for(Int_t j=0;j<kNHisto;j++){
    histXZ->GetEntry(j);
    fHisResXZ[j]->Add(h2);
  }

  nevent = (Int_t)histClusterType->GetEntries();
  if(nevent!=kNHisto) 
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  histClusterType->SetAddress(&h2i);
  for(Int_t j=0;j<kNHisto;j++){
    histClusterType->GetEntry(j);
    fHisClusterSize[j]->Add(h2i);
  }

  for(Int_t clu=0; clu<kNclu; clu++) {

    nevent = (Int_t)histXclu[clu]->GetEntries();
    if(nevent!=kNHisto)
      {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histXclu[clu]->SetAddress(&h);
    for(Int_t j=0;j<kNHisto;j++){
      histXclu[clu]->GetEntry(j);
      fHisResXclu[j][clu]->Add(h);
    }

   nevent = (Int_t)histZclu[clu]->GetEntries();
    if(nevent!=kNHisto)
      {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histZclu[clu]->SetAddress(&h);
    for(Int_t j=0;j<kNHisto;j++){
      histZclu[clu]->GetEntry(j);
      fHisResZclu[j][clu]->Add(h);
    }
  }


    for(Int_t chip=0; chip<kNChip; chip++) {

    nevent = (Int_t)histXchip[chip]->GetEntries();
    if(nevent!=kNHisto)
      {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histXchip[chip]->SetAddress(&h);
    for(Int_t j=0;j<kNHisto;j++){
      histXchip[chip]->GetEntry(j);
      fHisResXchip[j][chip]->Add(h);
    }

    nevent = (Int_t)histZchip[chip]->GetEntries();
    if(nevent!=kNHisto)
      {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histZchip[chip]->SetAddress(&h);
    for(Int_t j=0;j<kNHisto;j++){
      histZchip[chip]->GetEntry(j);
      fHisResZchip[j][chip]->Add(h);
    }
  }

  nevent = (Int_t)histTrErrX->GetEntries(); 
  if(nevent!=kNHisto)
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  histTrErrX->SetAddress(&h);
  for(Int_t j=0;j<kNHisto;j++){
    histTrErrX->GetEntry(j);
    fHisTrackErrX[j]->Add(h);
  }

  nevent = (Int_t)histTrErrZ->GetEntries();
  if(nevent!=kNHisto)
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  histTrErrZ->SetAddress(&h);
  for(Int_t j=0;j<kNHisto;j++){
    histTrErrZ->GetEntry(j);
    fHisTrackErrZ[j]->Add(h);
  }

  nevent = (Int_t)histClErrX->GetEntries();
  if(nevent!=kNHisto)
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  histClErrX->SetAddress(&h);
  for(Int_t j=0;j<kNHisto;j++){
    histClErrX->GetEntry(j);
    fHisClusErrX[j]->Add(h);
  }

  nevent = (Int_t)histClErrZ->GetEntries();
  if(nevent!=kNHisto)
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  histClErrZ->SetAddress(&h);
  for(Int_t j=0;j<kNHisto;j++){
    histClErrZ->GetEntry(j);
    fHisClusErrZ[j]->Add(h);
  }
  nevent = (Int_t)profXvsPhi->GetEntries();
  if(nevent!=kNHisto)
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  profXvsPhi->SetAddress(&p);
  for(Int_t j=0;j<kNHisto;j++){
    profXvsPhi->GetEntry(j);
    fProfResXvsPhi[j]->Add(p);
  }

  nevent = (Int_t)profZvsDip->GetEntries();
  if(nevent!=kNHisto)
    {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
  profZvsDip->SetAddress(&p);
  for(Int_t j=0;j<kNHisto;j++){ 
    profZvsDip->GetEntry(j);
    fProfResZvsDip[j]->Add(p);
  }

    for(Int_t clu=0; clu<kNclu; clu++) {

    nevent = (Int_t)profXvsPhiclu[clu]->GetEntries();
    if(nevent!=kNHisto)
      {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    profXvsPhiclu[clu]->SetAddress(&p);
    for(Int_t j=0;j<kNHisto;j++){
      profXvsPhiclu[clu]->GetEntry(j);
      fProfResXvsPhiclu[j][clu]->Add(p);
    }

    nevent = (Int_t)profZvsDipclu[clu]->GetEntries();
    if(nevent!=kNHisto)
      {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    profZvsDipclu[clu]->SetAddress(&p);
    for(Int_t j=0;j<kNHisto;j++){ 
      profZvsDipclu[clu]->GetEntry(j);
      fProfResZvsDipclu[j][clu]->Add(p);
    }
  }

    for(Int_t phas=0; phas<kNClockPhase;phas++){

    nevent = (Int_t)histXtrkFOtrue[phas]->GetEntries();
    if(nevent!=kNHisto)
       {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histXtrkFOtrue[phas]->SetAddress(&h);
    for(Int_t j=0;j<kNHisto;j++){
      histXtrkFOtrue[phas]->GetEntry(j);
      fHisTrackXFOtrue[j][phas]->Add(h);
    }

    nevent = (Int_t)histZtrkFOtrue[phas]->GetEntries();
    if(nevent!=kNHisto)
       {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histZtrkFOtrue[phas]->SetAddress(&h);
    for(Int_t j=0;j<kNHisto;j++){
      histZtrkFOtrue[phas]->GetEntry(j);
      fHisTrackZFOtrue[j][phas]->Add(h);
    }

    nevent = (Int_t)histXtrkFOfalse[phas]->GetEntries();
    if(nevent!=kNHisto)
       {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histXtrkFOfalse[phas]->SetAddress(&h);
    for(Int_t j=0;j<kNHisto;j++){
      histXtrkFOfalse[phas]->GetEntry(j);
      fHisTrackXFOfalse[j][phas]->Add(h);
    }

    nevent = (Int_t)histZtrkFOfalse[phas]->GetEntries();
    if(nevent!=kNHisto)
       {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histZtrkFOfalse[phas]->SetAddress(&h);
    for(Int_t j=0;j<kNHisto;j++){
      histZtrkFOfalse[phas]->GetEntry(j);
      fHisTrackZFOfalse[j][phas]->Add(h);
    }

    nevent = (Int_t)histXZtrkFOtrue[phas]->GetEntries();
    if(nevent!=kNHisto)
       {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histXZtrkFOtrue[phas]->SetAddress(&h2);
    for(Int_t j=0;j<kNHisto;j++){
      histXZtrkFOtrue[phas]->GetEntry(j);
      fHisTrackXZFOtrue[j][phas]->Add(h2);
    }

    nevent = (Int_t)histXZtrkFOfalse[phas]->GetEntries();
    if(nevent!=kNHisto)
       {AliWarning("ReadHistosFromFile: trying to read too many or too few histos!"); return kFALSE;}
    histXZtrkFOfalse[phas]->SetAddress(&h2);
    for(Int_t j=0;j<kNHisto;j++){
      histXZtrkFOfalse[phas]->GetEntry(j);
      fHisTrackXZFOfalse[j][phas]->Add(h2);
    }

   }

  delete h;   
  delete h2;  
  delete h2i; 
  delete p;   

  if (file) {
    file->Close();
    delete file;
  }
return kTRUE;
}

