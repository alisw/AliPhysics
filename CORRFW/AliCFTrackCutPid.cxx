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

//This class is intended to provide a selection on
//the PID for single charged particles as electrons, muons
//pions, kaons and protons. The user is supposed to set only one particle 
//to look at. The class at the moment uses one single ESD track.
//The identification is done in Identify(), whereas the GetID() checks 
//the track status or if the combined PID should be applied.
//The Correction Framework follows the Annalysis framework. 
//The main method of this class is the IsSelected function which returns 
//true (false) if the ESD track is (not) identified as the previously 
//setted particle.
//
//
//usage:
//
// -----ID by one detector------
//AliCFTrackCutPid *pCut = new AliCFTrackCutPid("pion_sel","pion_sel");
//Double_t priors[5]={...};
//pCut->SetPartycleType(AliPID::kPion,kFALSE)  
//pCut->SetDetectors("DET");          ^^^^^^
//pCut->SetPriors(priors);
//
// -----ID combining more detectors----
//AliCFTrackCutPid *pCut = new AliCFTrackCutPid("pion_sel","pion_sel");
//Double_t priors[5]={...};
//pCut->SetPartycleType(AliPID::kPion,kTRUE)
//pCut->SetDetectors("DET1 DET2 .."); ^^^^^
//pCut->SetPriors(priors)
//
//Comments: 
//The user can choose not to identify a track if one residual beteween 
//the identified particle probability and one among all the other 
//probabilties is smaller than a predefined value. 
//The same can be done for the detector responses. 
//This happens by setting:
// 
//pCut->SetMinDiffProb(0.005,kTRUE) //for probabilities
//
//pCut->SetMinDiffResp(0.005,kTRUE) //for det responses
//
//Annalisa.Mastroserio@ba.infn.it
//



#include "AliCFTrackCutPid.h"
#include "AliLog.h"
#include <TMath.h>
#include <TList.h>

ClassImp(AliCFTrackCutPid) 

//__________________________________________________________________________________
// CUT ON TRACK PID
//__________________________________________________________________________________
AliCFTrackCutPid::AliCFTrackCutPid() :
  AliCFCutBase(),
  fCut(0.2),
  fMinDiffResponse(0.0001),
  fMinDiffProbability(0.001),
  fgParticleType(10),
  fgIsComb(kTRUE),
  fgIsAOD(kFALSE),
  fCheckResponse(kFALSE),
  fCheckSelection(kTRUE),
  fIsPpriors(kFALSE),
  fIsDetAND(kFALSE),
  fXmin(-0.005),
  fXmax(1.005),
  fNbins(101),
  fDetRestr(-1),
  fiPartRestr(-1),
  fDetProbRestr(1),
  fProbThreshold(0.)
{ 
  //
  //Default constructor 
  //
  for(Int_t j=0; j< AliPID::kSPECIESN; j++) {
    fPriors[j]=0.2;
  }
  for(Int_t j=0; j< AliPID::kSPECIES; j++) {
    fPriorsFunc[j]=0x0;
  }
  for(Int_t jDet=0; jDet< kNdets; jDet++)  {
    fDets[jDet]=kFALSE;
    fDetsInAnd[jDet]=kFALSE;
  }
  
  InitialiseHisto();
}
//______________________________
AliCFTrackCutPid::AliCFTrackCutPid(const Char_t* name, const Char_t* title) :
  AliCFCutBase(name,title),
  fCut(0.2),
  fMinDiffResponse(0.0001),
  fMinDiffProbability(0.001),
  fgParticleType(10),
  fgIsComb(kTRUE),
  fgIsAOD(kFALSE),
  fCheckResponse(kFALSE),
  fCheckSelection(kTRUE),
  fIsPpriors(kFALSE),
  fIsDetAND(kFALSE),
  fXmin(-0.005),
  fXmax(1.005),
  fNbins(101),
  fDetRestr(-1),
  fiPartRestr(-1),
  fDetProbRestr(1),
  fProbThreshold(0.)
{
  //
  //Constructor
  // 
  for(Int_t j=0; j< AliPID::kSPECIESN; j++) {
    fPriors[j]=0.2;
  }
  for(Int_t j=0; j< AliPID::kSPECIES; j++) {
    fPriorsFunc[j]=0x0;
  }
  for(Int_t jDet=0; jDet< kNdets; jDet++)  {
    fDets[jDet]=kFALSE;
    fDetsInAnd[jDet]=kFALSE;
  }
  
  InitialiseHisto();
}
//______________________________
AliCFTrackCutPid::AliCFTrackCutPid(const AliCFTrackCutPid& c) :
  AliCFCutBase(c),
  fCut(c.fCut),
  fMinDiffResponse(c.fMinDiffResponse),
  fMinDiffProbability(c.fMinDiffProbability),
  fgParticleType(c.fgParticleType),
  fgIsComb(c.fgIsComb),
  fgIsAOD(c.fgIsAOD),
  fCheckResponse(c.fCheckResponse),
  fCheckSelection(c.fCheckSelection),
  fIsPpriors(c.fIsPpriors),
  fIsDetAND(c.fIsDetAND),
  fXmin(c.fXmin),
  fXmax(c.fXmax),
  fNbins(c.fNbins),
  fDetRestr(c.fDetRestr),
  fiPartRestr(c.fiPartRestr),
  fDetProbRestr(c.fDetProbRestr),
  fProbThreshold(c.fProbThreshold)
{
  //
  //Copy constructor
  //
  for(Int_t i=0; i< kNdets ; i++ ) { 
    fDets[i]=c.fDets[i]; 
    fDetsInAnd[i]=c.fDetsInAnd[i];
    for(Int_t iP =0; iP<AliPID::kSPECIES; iP++){
      fhResp[i][iP]=c.fhResp[i][iP];
      fhProb[i][iP]=c.fhProb[i][iP];
    }
  }
  for(Int_t j=0; j< AliPID::kSPECIESN; j++){
    fPriors[j]=c.fPriors[j];
  }
  for(Int_t j=0; j< AliPID::kSPECIES; j++){
    fPriorsFunc[j]=c.fPriorsFunc[j];
    fhCombResp[j]=c.fhCombResp[j];
    fhCombProb[j]=c.fhCombProb[j];
  }
}
//______________________________
AliCFTrackCutPid& AliCFTrackCutPid::operator=(const AliCFTrackCutPid& c)
{  
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    this->fCut=c.fCut;
    this->fMinDiffResponse=c.fMinDiffResponse;
    this->fMinDiffProbability=c.fMinDiffProbability;
    this->fgParticleType=c.fgParticleType;  
    this->fgIsComb=c.fgIsComb;
    this->fgIsAOD=c.fgIsAOD;
    this->fCheckResponse=c.fCheckResponse;
    this->fCheckSelection=c.fCheckSelection;
    this->fIsPpriors=c.fIsPpriors;
    this->fIsDetAND=c.fIsDetAND;
    this->fXmin=c.fXmin;
    this->fXmax=c.fXmax;
    this->fNbins=c.fNbins;
    this->fDetRestr=c.fDetRestr;
    this->fiPartRestr=c.fiPartRestr;
    this->fDetProbRestr=c.fDetProbRestr;
    this->fProbThreshold=c.fProbThreshold;
  
    for(Int_t i=0; i< kNdets ; i++ )   {
      this->fDets[i]=c.fDets[i];
      for(Int_t iP =0; iP<AliPID::kSPECIES; iP++){
	this->fhResp[i][iP]=c.fhResp[i][iP];
	this->fhProb[i][iP]=c.fhProb[i][iP];
      }  
    }

    for(Int_t j=0; j< AliPID::kSPECIESN; j++){
      this->fPriors[j]=c.fPriors[j];
    }
    for(Int_t j=0; j< AliPID::kSPECIES; j++){
      this->fhCombResp[j]=c.fhCombResp[j];
      this->fhCombProb[j]=c.fhCombProb[j];
      this-> fPriorsFunc[j]=c.fPriorsFunc[j];

    }
  }
  return *this ;
}
//__________________________________________________________________________________
AliCFTrackCutPid::~AliCFTrackCutPid() {
  //
  //dtor
  //

  for(Int_t i=0; i< kNdets ; i++ )   {
    for(Int_t iP =0; iP<AliPID::kSPECIES; iP++){
      if(fhResp[i][iP])delete fhResp[i][iP];
      if(fhProb[i][iP])delete fhProb[i][iP];
    }  
  }
  
  for(Int_t j=0; j< AliPID::kSPECIES; j++){
    if(fhCombResp[j])delete fhCombResp[j];
    if(fhCombProb[j])delete fhCombProb[j];
    
  }
}
//__________________________________
void AliCFTrackCutPid::SetDetectors(TString dets)
{
  //
  // The string of detectors is translated into 
  // the respective booelan data members

  if(dets.Contains("ITS")) {fDets[kITS]=kTRUE;}
  if(dets.Contains("TPC")) {fDets[kTPC]=kTRUE;}
  if(dets.Contains("TRD")) {fDets[kTRD]=kTRUE;}
  if(dets.Contains("TOF")) {fDets[kTOF]=kTRUE;}
  if(dets.Contains("HMPID")) {fDets[kHMPID]=kTRUE;}

  if(dets.Contains("ALL")) for(Int_t i=0; i< kNdets ; i++) fDets[i]=kTRUE;
}
//__________________________________
void AliCFTrackCutPid::SetPriors(Double_t r[AliPID::kSPECIES])
{
  //
  // Sets the a priori concentrations
  // 
  for(Int_t i=0; i<AliPID::kSPECIES; i++) fPriors[i]=r[i];
}
//__________________________________
void AliCFTrackCutPid::SetPriorFunctions(TF1 *func[AliPID::kSPECIES])
{
  //
  // Sets the momentu dependent a priori concentrations
  //
  
  for(Int_t i=0; i<AliPID::kSPECIES; i++) fPriorsFunc[i]=func[i];
  fIsPpriors = kTRUE;
}
//____________________________________________
void AliCFTrackCutPid::SetANDstatus(TString dets)
{
  //
  //Sets the detectors to be in AND for the combined PID
  //
  if(dets.Contains("ITS") && fDets[kITS])     {fDetsInAnd[kITS]=kTRUE;}
  if(dets.Contains("TPC") && fDets[kTPC])     {fDetsInAnd[kTPC]=kTRUE;}
  if(dets.Contains("TRD") && fDets[kTRD])     {fDetsInAnd[kTRD]=kTRUE;}
  if(dets.Contains("TOF") && fDets[kTOF])     {fDetsInAnd[kTOF]=kTRUE;}
  if(dets.Contains("HMPID") && fDets[kHMPID]) {fDetsInAnd[kHMPID]=kTRUE;}
  
  fIsDetAND = kTRUE;
}
//
void AliCFTrackCutPid::SetDetectorProbabilityRestriction(TString det, Int_t iPart, Double_t upperprob)
{
  //
  // Sets the detector, the particle and the probability
  // limit.
  //
  
  if(det.Contains("ITS")) fDetRestr = kITS;
  if(det.Contains("TPC")) fDetRestr = kTPC;
  if(det.Contains("TRD")) fDetRestr = kTRD;
  if(det.Contains("TOF")) fDetRestr = kTOF;
  if(det.Contains("HMPID")) fDetRestr = kHMPID;
  fiPartRestr = iPart; 
  fDetProbRestr = upperprob;
}
//__________________________________
void AliCFTrackCutPid::TrackInfo(const AliESDtrack *pTrk, ULong_t status[kNdets+1],Double_t pid[kNdets+1][AliPID::kSPECIES]) const
{
  //
  // Loads the responses of the XXX chosen detectors for the pTrk
  // and the corresponding trk status. The final trk status is also loaded.
  // 
  if(fDets[kITS]) {
    pTrk->GetITSpid(pid[kITS]);
    status[kITS]=AliESDtrack::kITSpid;
  }
  if(fDets[kTPC]) {
    pTrk->GetTPCpid(pid[kTPC]);
    status[kTPC]=AliESDtrack::kTPCpid;
  }
  if(fDets[kTRD]) {
    pTrk->GetTRDpid(pid[kTRD]);
    status[kTRD]=AliESDtrack::kTRDpid;
  }
  if(fDets[kTOF]) {
    pTrk->GetTOFpid(pid[kTOF]);
    status[kTOF]=AliESDtrack::kTOFpid;
  }
  if(fDets[kHMPID]) {
    pTrk->GetHMPIDpid(pid[kHMPID]);
    status[kHMPID]=AliESDtrack::kHMPIDpid;
  }
  if(fDetRestr>=0){
    if(fDetRestr == kITS) pTrk->GetITSpid(pid[kITS]);
    if(fDetRestr == kTPC) pTrk->GetITSpid(pid[kTPC]);
    if(fDetRestr == kTRD) pTrk->GetITSpid(pid[kTRD]);
    if(fDetRestr == kTOF) pTrk->GetITSpid(pid[kTOF]);
    if(fDetRestr == kHMPID) pTrk->GetITSpid(pid[kHMPID]);
  }
  
  status[kNdets]=pTrk->GetStatus();
  pTrk->GetESDpid(pid[kNdets]);
}
//__________________________________
void AliCFTrackCutPid::SetPPriors(AliESDtrack *pTrk)
{
  //
  //sets the mommentum dependent a priori concentrations
  //
  
  for(Int_t i=0; i< AliPID::kSPECIES; i++) {
    if(pTrk->P()>fPriorsFunc[i]->GetXmin() && pTrk->P() < fPriorsFunc[i]->GetXmax()) fPriors[i]=fPriorsFunc[i]->Eval(pTrk->P());
    else {AliInfo("the track momentum is not in the function range. Priors are equal") fPriors[i] = 0.2;}   
  }
}   
//______________________________________
ULong_t AliCFTrackCutPid::StatusForAND(ULong_t status[kNdets+1]) const
{
  //
  //In case of AND of more detectors the AND-detector status combination. 
  //is calculated and also returned
  //

  ULong_t andstatus=0;
  for(Int_t i=0; i< kNdets; i++) {
    if(!fDetsInAnd[i]) continue;
    andstatus = andstatus | status[i];
    AliDebug(1,Form("trk status %lu  %i AND-status  combination: %lu",status[i],i,andstatus));
  }
  return andstatus;
}
//_______________________________________
Int_t AliCFTrackCutPid::GetID(ULong_t status[kNdets+1],Double_t pid[kNdets+1][AliPID::kSPECIES]) const
{
  // Identifies the track if its probability is higher than the cut
  // value. The default value is fCut=0.2, therefore the most probable
  // one is identified by default. Here all the checks on how to identify
  // the track are performed (single detector or combined PID,..., detector
  // restriction probability
  // Returns:   integer corresponding to the identified particle
  
  Int_t iPart=-1;
  
  if(!fgIsComb){  
    Bool_t isDet=kFALSE;
    for(Int_t i=0; i<kNdets; i++){
      if(!fDets[i]) continue;
      isDet=kTRUE;
      AliDebug(1,Form("trk status %lu  %i-det-pid status %lu   -> combination: %lu",status[kNdets],i,status[i],status[kNdets]&status[i]));
      if(!(status[kNdets]&status[i])){
	iPart=-10;
	AliDebug(1,Form("detector %i -> pid trk status not ok",i));
      }
      else {
        AliDebug(1,Form("resp   : %f  %f  %f  %f  %f",pid[i][0],pid[i][1],pid[i][2],pid[i][3],pid[i][4]));
	if(fIsQAOn) iPart = IdentifyQA(pid[i],i);
	else iPart = Identify(pid[i]);
      } 
    }  
    if(!isDet){
      AliDebug(1,Form("  !!!        No detector selected, the ESD-pid response is considered"));
      iPart = Identify(pid[kNdets]);
    }
  }else{
    Double_t calcprob[5];
    CombPID(status,pid,calcprob); 
    iPart = Identify(calcprob);
    
  }
  
  
  AliDebug(1,Form("selected particle: %i",iPart));
  
  if(iPart >=0 && fiPartRestr>=0) {
    AliPID restr(pid[fDetRestr]);
    restr.SetPriors(fPriors);
    AliDebug(1,Form("setted upper limit: %f  det %i : probability %f ",fDetProbRestr,fDetRestr,restr.GetProbability((AliPID::EParticleType)fiPartRestr)));
    if(restr.GetProbability((AliPID::EParticleType)fiPartRestr) > fDetProbRestr) {
      iPart = kDetRestr;
      AliDebug(1,"\n\n the detector restrictions refused the ID \n\n");
    }
  }//cross checks with one detector probability
  
  AliDebug(1,Form("after the check the selected particle is %i",iPart));
  
  return iPart;
}
//_________________________________________________________________________________
Int_t AliCFTrackCutPid::GetAODID(AliAODTrack *aodtrack) const
{
//
// Identifies the AOD Track using the combined pid responses
//

  Double_t combpid[AliPID::kSPECIES];
  for(Int_t i=0; i< AliPID::kSPECIES; i++) {
    combpid[i]= aodtrack->PID()[i];
     if(!fhCombResp[i]) AliDebug(1,Form("\n no fhCombResp[%i], check if pidcut->Init() was called",i));
     else fhCombResp[i]->Fill(combpid[i]);
  }
  return Identify(combpid);
}
//__________________________________
Bool_t AliCFTrackCutPid::Check(const Double_t *p, Int_t iPsel, Double_t minDiff) const
{
  // Checks if there are no equal values and if the valus
  // difference between the selected particle and the others i
  // is higher than  a lower limit.
  // Returns:  kTRUE= is acceptable
  
  AliDebug(2,Form("input particle: %i",iPsel));
  Bool_t ck=kTRUE;
  
  if(iPsel<0) ck=kFALSE;
  
  else {
    for(Int_t j=0; j< AliPID::kSPECIES; j++) {
      if(j!=iPsel && TMath::Abs(p[j]-p[iPsel])<minDiff) ck=kFALSE;
    }
    if(!ck) AliDebug(1,"the values are too close ");
  }
  return ck;
}
//___________________________________________
Int_t AliCFTrackCutPid::Identify(Double_t pid[AliPID::kSPECIES]) const
{
  //
  // The identification is actually performed here with possible
  // checks on the det responses and/or probabilities
  //

  Int_t iPart = -1;
 
  AliDebug(2,Form("calc response bef: %f  %f  %f  %f  %f",pid[0],pid[1],pid[2],pid[3],pid[4]));
  AliDebug(2,Form("priors           : %f  %f  %f  %f  %f",fPriors[0],fPriors[1],fPriors[2],fPriors[3],fPriors[4]));

  AliPID getpid(pid,kTRUE);
  getpid.SetPriors(fPriors); 
  
  Double_t probability[AliPID::kSPECIES]={0.,0.,0.,0.,0.};
  for(Int_t iP=0; iP<AliPID::kSPECIES; iP++) {
    probability[iP] = getpid.GetProbability((AliPID::EParticleType)iP);
    AliDebug(2,Form("prob %i %f",iP, probability[iP]));
    if(fIsQAOn) fhCombProb[iP]->Fill(probability[iP]);
  }
  

  if (fProbThreshold > 0.) {
    if (probability[fgParticleType] >= fProbThreshold) iPart=fgParticleType;
  }
  else {
    AliPID::EParticleType sel = getpid.GetMostProbable();
    if(getpid.GetProbability(sel,fPriors)>fCut) iPart= (Int_t)sel;
    AliDebug(2,Form("probabilities   : %f  %f  %f  %f  %f",probability[0],probability[1],probability[2],probability[3],probability[4]));
  }

  if(fCheckResponse && !Check(pid,iPart, fMinDiffResponse)) iPart=kCheckResp;
  if(fCheckSelection && !Check(probability,iPart,fMinDiffProbability)) iPart=kCheckProb;

  return iPart;
}
//___________________________________________
Int_t AliCFTrackCutPid::IdentifyQA(const Double_t pid[AliPID::kSPECIES], Int_t idets) const
{
  //
  // The same as Identify, but with the QA histo filling 
  //
  //
  
  Int_t iPart = -1;
  AliDebug(1,Form("resp   : %f  %f  %f  %f  %f",pid[0],pid[1],pid[2],pid[3],pid[4]));
  
  AliPID getpid(pid,kTRUE);
  getpid.SetPriors(fPriors);
      
  AliPID::EParticleType sel = getpid.GetMostProbable();
  Double_t probability[AliPID::kSPECIES];
  for(Int_t iP=0; iP<AliPID::kSPECIES; iP++) {
    probability[iP] = getpid.GetProbability((AliPID::EParticleType)iP);
    fhProb[idets][iP]->Fill(probability[iP]);
  }
  
  AliPID toresp(pid,kTRUE); 
  Double_t qapriors[10]={0.2,0.2,0.2,0.2,0.2,0,0,0,0,0};
  toresp.SetPriors(qapriors);
  for(Int_t iPr=0; iPr<AliPID::kSPECIES; iPr++) fhResp[idets][iPr]->Fill(toresp.GetProbability((AliPID::EParticleType)iPr));
  
  if(getpid.GetProbability(sel,fPriors)>fCut) iPart= (Int_t)sel;
  AliDebug(1,Form("resp   : %f  %f  %f  %f  %f",pid[0],pid[1],pid[2],pid[3],pid[4]));
  AliDebug(1,Form("probab : %f  %f  %f  %f  %f",probability[0],probability[1],probability[2],probability[3],probability[4]));
  if(fCheckResponse && !Check(pid,iPart, fMinDiffResponse)) iPart=kCheckResp;
  if(fCheckSelection && !Check(probability,iPart,fMinDiffProbability)) iPart=kCheckProb;
  return iPart;
}
//___________________________________________
Bool_t AliCFTrackCutPid::IsSelected(TObject *track){
  //
  //  method for the pid-cut selction
  //
  Bool_t sel = kFALSE;
  
  if (!track) return kFALSE ;
  TString className(track->ClassName());
  if (className.CompareTo("AliESDtrack") == 0) {
  AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track); 
  if (!esdTrack) return kFALSE;
  ULong_t status[kNdets+1]={0,0,0,0,0,0};
  Double_t pid[kNdets+1][AliPID::kSPECIES];
  TrackInfo(esdTrack,status,pid);
  if(fIsPpriors) SetPPriors(esdTrack);
  if(GetID(status,pid)==fgParticleType) sel = kTRUE;
  }

  if (className.CompareTo("AliAODTrack") == 0) {
  AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(track);
  if (!aodtrack) return kFALSE ;
  if(GetAODID(aodtrack) == fgParticleType) sel = kTRUE;
  }

 return sel;

}
//__________________________________
void  AliCFTrackCutPid::CombPID(ULong_t status[kNdets+1],Double_t pid[kNdets+1][AliPID::kSPECIES],Double_t *combpid) const
{
  //
  // Calculates the combined PID according to the chosen detectors.
  // and provides the array of probabilities
  //
  
  Bool_t isdet=kFALSE;
  Double_t prod[AliPID::kSPECIES]={1.,1.,1.,1.,1.};
  
  ULong_t andstatus =0;
  if(fIsDetAND) {
    andstatus = StatusForAND(status);
    AliDebug(1,Form("AND combination %lu",andstatus));
  } 
  
  //Products of single detector responses
  for(Int_t j=0; j<AliPID::kSPECIES; j++){
    for(Int_t i=0; i< kNdets; i++){
      if(!fDets[i]) continue;
      if(status[kNdets]&status[i]) {
        if(fIsDetAND) {
	  ULong_t checkstatus = status[kNdets]&andstatus;
	  if(checkstatus != andstatus) continue; 
	  else {
	    prod[j]*=pid[i][j];
	    isdet = kTRUE;
	    AliDebug(1,Form("-----> trk status %lu   and status %lu -> trk-ANDdetector status combination %lu",status[kNdets],andstatus,status[kNdets]&andstatus));
	    AliDebug(1,Form("In det %i ->  particle %i response is %f",i,j,pid[i][j]));
	  }
        } else {
	  prod[j]*=pid[i][j];
	  isdet=kTRUE;
	  AliDebug(2,Form("In det %i ->  particle %i response is %f",i,j,pid[i][j]));
          
           if(fIsQAOn){
             if(!fhResp[i][j]) {AliDebug(1,Form("no pointer to the histo fhResp%i%i, check if pidcut->Init() was called",i,j));}
	     else fhResp[i][j]->Fill(pid[i][j]);    
             
             if(!fhProb[i][j]) {AliDebug(1,Form("no pointer to the histo fhProb%i%i, check if pidcut->Init() was called",i,j));}
	     else {
               AliPID detprob(pid[i],kTRUE); 
               detprob.SetPriors(fPriors);
               fhProb[i][j]->Fill(detprob.GetProbability((AliPID::EParticleType)j));    
             }       
           }
	}
      }//combined mode
    }//loop on dets
  }//loop on species

  
   //no detectors found, then go to ESD pid...
  if(!isdet) {
    AliWarning("\n !! No detector found for the combined pid --> responses are from the ESDpid !! \n");
    Double_t sumesdpid=0;
     for(Int_t nn=0; nn<AliPID::kSPECIES; nn++) sumesdpid+=pid[kNdets][nn];
     if(sumesdpid<=0) {
        AliDebug(1,"priors or ESDpid are inconsistent, please check them");
        return;
      } else {
        for(Int_t k=0; k<AliPID::kSPECIES; k++){
         combpid[k] = pid[kNdets][k]/sumesdpid;
          if(fIsQAOn) {
           if(!fhCombResp[k]) AliDebug(1,Form("\n no fhCombResp[%i], check if pidcut->Init() was called",k));
           else fhCombResp[k]->Fill(combpid[k]);
          }
         }//loop on species
       }
      return;   
    }   
    
  Double_t add = 0; for(Int_t isumm=0; isumm<5; isumm++) add+=prod[isumm];
  if(add>0) {
    for(Int_t ip =0; ip < AliPID::kSPECIES; ip++) {
      combpid[ip] =  prod[ip]/add;
     if(fIsQAOn) {
           if(!fhCombResp[ip]) AliDebug(1,Form("\n no fhCombResp[%i], check if pidcut->Init() was called",ip));
           else fhCombResp[ip]->Fill(combpid[ip]); 
           
         }
   }
   AliDebug(1,Form("calculated comb response: %f %f %f %f %f",combpid[0],combpid[1],combpid[2],combpid[3],combpid[4]));
 } else {
    AliDebug(1,"single detector responses are inconsistent, please check them....");
    return; 
   }
  AliDebug(1,Form("the ESDpid response:      %f %f %f %f %f",pid[kNdets][0],pid[kNdets][1],pid[kNdets][2],pid[kNdets][3],pid[kNdets][4]));
 
}
//__________________________________________
//
//QA part
//_________________________________________
void AliCFTrackCutPid::InitialiseHisto()
{
  //
  //QA histo initialization
  //
  for(Int_t iP=0; iP<AliPID::kSPECIES; iP++){
    fhCombResp[iP]=0x0;
    fhCombProb[iP]=0x0;
    for(Int_t iDet =0; iDet<kNdets; iDet++){
      fhResp[iDet][iP]=0x0;
      fhProb[iDet][iP]=0x0;
    }
  } 
}
//______________________________________________
void AliCFTrackCutPid::DefineHistograms()
{
  //
  //QA histo booking
  //

  if(fgIsAOD){
     const char *partic[AliPID::kSPECIES]={"electron","muon","pion","kaon","proton"}; 

     for(Int_t iPart =0; iPart < AliPID::kSPECIES; iPart++)
      {
       fhCombResp[iPart] = new TH1F(Form("%s_rCombPart%i",GetName(),iPart),Form(" %s combined response  (AODTrack)  ",partic[iPart]),fNbins,fXmin,fXmax);
       fhCombResp[iPart]->SetXTitle(Form(" %s combined response ",partic[iPart]));
       fhCombResp[iPart]->SetYTitle("entries");
       AliDebug(1,Form(  "%s is booked!!",fhCombResp[iPart]->GetName()));
       fhCombProb[iPart] = new TH1F(Form("%s_pCombPart%i",GetName(),iPart),Form("%s combined probability (AODTrack) ",partic[iPart]),fNbins,fXmin,fXmax);
       fhCombProb[iPart]->SetXTitle(Form(" %s combined probability ",partic[iPart]));
       fhCombProb[iPart]->SetYTitle("entries");
       AliDebug(1,Form(  "%s is booked!!",fhCombProb[iPart]->GetName()));
     }
   }


  else {
   const char *detect[kNdets]={"ITS","TPC","TRD","TOF","HMPID"};
   const char *partic[AliPID::kSPECIES]={"electron","muon","pion","kaon","proton"};
  
    for(Int_t iDet =0; iDet< kNdets; iDet++)
     {
       if(!fDets[iDet]) continue;
       for(Int_t iP =0; iP < AliPID::kSPECIES; iP++){
 	fhResp[iDet][iP] = new TH1F(Form("%s_rDet%iPart%i",GetName(),iDet,iP),Form("%s response for %s  ",detect[iDet],partic[iP]),fNbins,fXmin,fXmax);
        fhResp[iDet][iP]->SetXTitle(Form(" %s response ",partic[iP]));
        fhResp[iDet][iP]->SetYTitle("entries");
 	fhProb[iDet][iP] = new TH1F(Form("%s_pDet%iPart%i",GetName(),iDet,iP),Form("%s calculated probability for %s",detect[iDet],partic[iP]),fNbins,fXmin,fXmax);
        fhProb[iDet][iP]->SetXTitle(Form(" %s probability ",partic[iP]));
        fhProb[iDet][iP]->SetYTitle("entries");
       }
     }
    

  if(fgIsComb)
     { 
      for(Int_t iPart =0; iPart < AliPID::kSPECIES; iPart++)
        {
	  fhCombResp[iPart] = new TH1F(Form("%s_rCombPart%i",GetName(),iPart),Form(" %s combined response    ",partic[iPart]),fNbins,fXmin,fXmax);
          fhCombResp[iPart]->SetXTitle(Form(" %s response ",partic[iPart]));
          fhCombResp[iPart]->SetYTitle("entries");
          AliDebug(1,Form(  "%s is booked!!",fhCombResp[iPart]->GetName()));
	  fhCombProb[iPart] = new TH1F(Form("%s_pCombPart%i",GetName(),iPart),Form("%s combined probability ",partic[iPart]),fNbins,fXmin,fXmax);
          fhCombProb[iPart]->SetXTitle(Form(" %s response ",partic[iPart]));
          fhCombProb[iPart]->SetYTitle("entries");
          AliDebug(1,Form(  "%s is booked!!",fhCombProb[iPart]->GetName()));
        }
     }
   }

}
//___________________________________________________

void AliCFTrackCutPid::AddQAHistograms(TList *qalist) 
{
  //
  // adds QA histograms in a TList
  //
  if(!fIsQAOn) return;
  DefineHistograms();
  
  if(fgIsComb || fgIsAOD){
    for(Int_t iPart =0; iPart<AliPID::kSPECIES; iPart++){
      qalist->Add(fhCombResp[iPart]);
      qalist->Add(fhCombProb[iPart]);
    }
  }
  
  for(Int_t iDet=0; iDet<kNdets; iDet++){
    if(!fDets[iDet]) continue;
    for(Int_t iP =0; iP<AliPID::kSPECIES; iP++){
      qalist->Add(fhResp[iDet][iP]);
      qalist->Add(fhProb[iDet][iP]);
    }
  }  

}

