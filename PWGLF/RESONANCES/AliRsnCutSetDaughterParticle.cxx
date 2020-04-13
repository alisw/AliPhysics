//
// Author: Francesca Bellini (fbellini@cern.ch)
//
//

#include <Riostream.h>

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliRsnCutSetDaughterParticle.h"

class AliRsnCutPIDNSigma;
class AliRsnPIDRange;
class AliRsnCutPhi;

ClassImp(AliRsnCutSetDaughterParticle)

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle() :
   AliRsnCutSet("AliRsnCutSetDaughterParticle", AliRsnTarget::kDaughter),
   fPID(AliPID::kPion),
   fAppliedCutSetID(AliRsnCutSetDaughterParticle::kNDaughterCuts),
   fNsigmaTPC(1E20),
   fNsigmaTOF(1E20),
   fCutQuality(0x0),
   fAODTrkCutFilterBit(0),
   fCheckOnlyFilterBit(kTRUE),
   fUseCustomQualityCuts(kFALSE),
   fIsUse2011stdQualityCuts(kFALSE),  
   fIsUse2011stdQualityCutsHighPt(kFALSE)
{
   //
   // Default constructor
  SetPtRange(0.0, 1E20);
  SetEtaRange(1E20, 1E20);
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const char *name, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID, AliPID::EParticleType pid, Float_t nSigmaFast = -1.0, Int_t AODfilterBit = 0, Bool_t useTPCCrossedRows=kTRUE) :
   AliRsnCutSet(name, AliRsnTarget::kDaughter),
   fPID(pid),
   fAppliedCutSetID(cutSetID),
   fNsigmaTPC(nSigmaFast),
   fNsigmaTOF(nSigmaFast),
   fCutQuality(new AliRsnCutTrackQuality("CutQuality")),
   fAODTrkCutFilterBit(AODfilterBit),
   fCheckOnlyFilterBit(kTRUE),
   fUseCustomQualityCuts(kFALSE),
   fIsUse2011stdQualityCuts(kFALSE),  
   fIsUse2011stdQualityCutsHighPt(kFALSE)
{
  //
  // Constructor
  //
  //set here pt and eta range
  SetPtRange(0.15, 20.0);
  SetEtaRange(-0.8, 0.8);
  
  //if nsigma not specified, sets "no-PID" cuts
  if (nSigmaFast<=0){
    fNsigmaTPC=1e20;
    fNsigmaTOF=1e20;
    AliWarning("Requested fast n-sigma PID with negative value for n. --> Setting n = 1E20");
  }
  
  //initialize quality std and PID cuts
  InitStdQualityCuts(useTPCCrossedRows);
  Init();
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const char *name, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID, AliPID::EParticleType pid, Float_t nSigmaFastTPC = -1.0, Float_t nSigmaFastTOF = -1.0, Int_t AODfilterBit = 0, Bool_t useTPCCrossedRows=kTRUE) :
   AliRsnCutSet(name, AliRsnTarget::kDaughter),
   fPID(pid),
   fAppliedCutSetID(cutSetID),
   fNsigmaTPC(nSigmaFastTPC),
   fNsigmaTOF(nSigmaFastTOF),
   fCutQuality(new AliRsnCutTrackQuality("CutQuality")),
   fAODTrkCutFilterBit(AODfilterBit),
   fCheckOnlyFilterBit(kTRUE),
   fUseCustomQualityCuts(kFALSE),
   fIsUse2011stdQualityCuts(kFALSE),  
   fIsUse2011stdQualityCutsHighPt(kFALSE)
{
  //
  // Constructor
  //
  //set here pt and eta range
  SetPtRange(0.15, 20.0);
  SetEtaRange(-0.8, 0.8);
  
  //if nsigma not specified, sets "no-PID" cuts
  if (nSigmaFastTPC<=0){
    fNsigmaTPC=1e20;
    AliWarning("Requested fast n-sigma TPC PID with negative value for n. --> Setting n = 1E20");
  }
  if (nSigmaFastTOF<=0){
    fNsigmaTOF=1e20;
    AliWarning("Requested fast n-sigma TOF PID with negative value for n. --> Setting n = 1E20");
  }
  
  //initialize quality std and PID cuts
  InitStdQualityCuts(useTPCCrossedRows);
  Init();
}


//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const char *name, AliRsnCutTrackQuality *rsnTrackQualityCut, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID, AliPID::EParticleType pid, Float_t nSigmaFast = -1.0) :
  AliRsnCutSet(name, AliRsnTarget::kDaughter),
  fPID(pid),
  fAppliedCutSetID(cutSetID),
  fNsigmaTPC(nSigmaFast),
  fNsigmaTOF(nSigmaFast),
  fCutQuality(rsnTrackQualityCut),
  fAODTrkCutFilterBit(0),
  fCheckOnlyFilterBit(kFALSE),
  fUseCustomQualityCuts(kFALSE),
  fIsUse2011stdQualityCuts(kFALSE),  
  fIsUse2011stdQualityCutsHighPt(kFALSE)
{
  //
  // Constructor: uses externally-defined track-quality cut object
  //
  if (!rsnTrackQualityCut) {
    //if external track quality cut object not defined,
    //sets default track quality to be initialised (with cut on TPC crossed rows) +
    //sets here pt and eta cuts
    InitStdQualityCuts(kTRUE);
    SetPtRange(0.15, 20.0);
    SetEtaRange(-0.8, 0.8);
  } else {
    //checks validity of passed quality-cut object
    //if Ok, inherits all cuts including Pt and Eta cut
    fCheckOnlyFilterBit=kFALSE;
    fUseCustomQualityCuts = kTRUE;
    SetPtRange(rsnTrackQualityCut->GetPtRange(0), rsnTrackQualityCut->GetPtRange(1));
    SetEtaRange(rsnTrackQualityCut->GetEtaRange(0),rsnTrackQualityCut->GetEtaRange(1));
    AliInfo("Custom quality cuts applied");
    rsnTrackQualityCut->Print();
  } 
  
  //if nsigma not specified, sets "no-PID" cuts
  if (nSigmaFast<=0){
    fNsigmaTPC=1e20;
    fNsigmaTOF=1e20;
    AliWarning("Requested fast n-sigma PID with negative value for n. --> Setting n = 1E20");
  } 
  
  //initialize PID cuts
  Init(); 
}


//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const char *name, AliRsnCutTrackQuality *rsnTrackQualityCut, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID, AliPID::EParticleType pid, Float_t nSigmaFastTPC = -1.0, Float_t nSigmaFastTOF = -1.0) :
  AliRsnCutSet(name, AliRsnTarget::kDaughter),
  fPID(pid),
  fAppliedCutSetID(cutSetID),
  fNsigmaTPC(nSigmaFastTPC),
  fNsigmaTOF(nSigmaFastTOF),
  fCutQuality(rsnTrackQualityCut),
  fAODTrkCutFilterBit(0),
  fCheckOnlyFilterBit(kFALSE),
  fUseCustomQualityCuts(kFALSE),
  fIsUse2011stdQualityCuts(kFALSE),  
  fIsUse2011stdQualityCutsHighPt(kFALSE)
{
  //
  // Constructor: uses externally-defined track-quality cut object
  //
  if (!rsnTrackQualityCut) {
    //if external track quality cut object not defined,
    //sets default track quality to be initialised (with cut on TPC crossed rows) +
    //sets here pt and eta cuts
    InitStdQualityCuts(kTRUE);
    SetPtRange(0.15, 20.0);
    SetEtaRange(-0.8, 0.8);
  } else {
    //checks validity of passed quality-cut object
    //if Ok, inherits all cuts including Pt and Eta cut
    fCheckOnlyFilterBit=kFALSE;
    fUseCustomQualityCuts = kTRUE;
    SetPtRange(rsnTrackQualityCut->GetPtRange(0), rsnTrackQualityCut->GetPtRange(1));
    SetEtaRange(rsnTrackQualityCut->GetEtaRange(0),rsnTrackQualityCut->GetEtaRange(1));
    AliInfo("Custom quality cuts applied");
    rsnTrackQualityCut->Print();
  } 
  
  //if nsigma not specified, sets "no-PID" cuts
  if (nSigmaFastTPC<=0){
    fNsigmaTPC=1e20;
    AliWarning("Requested fast n-sigma TPC PID with negative value for n. --> Setting n = 1E20");
  }
  if (nSigmaFastTOF<=0){
    fNsigmaTOF=1e20;
    AliWarning("Requested fast n-sigma TOF PID with negative value for n. --> Setting n = 1E20");
  }
  
  //initialize PID cuts
  Init(); 
}


//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const AliRsnCutSetDaughterParticle &copy) :
  AliRsnCutSet(copy),
  fPID(copy.fPID),
  fAppliedCutSetID(copy.fAppliedCutSetID),
  fNsigmaTPC(copy.fNsigmaTPC),
  fNsigmaTOF(copy.fNsigmaTOF),
  fCutQuality(copy.fCutQuality),
  fAODTrkCutFilterBit(copy.fAODTrkCutFilterBit),
  fCheckOnlyFilterBit(copy.fCheckOnlyFilterBit),
  fUseCustomQualityCuts(copy.fUseCustomQualityCuts),
  fIsUse2011stdQualityCuts(copy.fIsUse2011stdQualityCuts),  
  fIsUse2011stdQualityCutsHighPt(copy.fIsUse2011stdQualityCutsHighPt)
{
  //
  // copy constructor
  SetPtRange(copy.fPtRange[0], copy.fPtRange[1]);
  SetEtaRange(copy.fEtaRange[0], copy.fEtaRange[1]);
}

//__________________________________________________________________________________________________
AliRsnCutSetDaughterParticle &AliRsnCutSetDaughterParticle::operator=(const AliRsnCutSetDaughterParticle &copy)
{
   //
   // operator =
   //
   AliRsnCutSet::operator=(copy);
   if (this == &copy)
      return *this;
   fPID=copy.fPID;
   fAppliedCutSetID=copy.fAppliedCutSetID;
   fNsigmaTPC=copy.fNsigmaTPC;
   fNsigmaTOF=copy.fNsigmaTOF;
   fAODTrkCutFilterBit=copy.fAODTrkCutFilterBit;
   fCheckOnlyFilterBit=copy.fCheckOnlyFilterBit;
   fUseCustomQualityCuts=copy.fUseCustomQualityCuts;
   fIsUse2011stdQualityCuts=copy.fIsUse2011stdQualityCuts;
   fIsUse2011stdQualityCutsHighPt=copy.fIsUse2011stdQualityCutsHighPt; 
   fCutQuality=copy.fCutQuality;
   return (*this);
}

//----------------------------------------------------------------------------
AliRsnCutSetDaughterParticle::~AliRsnCutSetDaughterParticle()
{
   //
   //destructor
   //
   if (fCutQuality)
      delete fCutQuality;
}

//----------------------------------------------------------------------------
void AliRsnCutSetDaughterParticle::Init()
{
  
  //define TOF match cut
  AliRsnCutTOFMatch  *iCutTOFMatch     = new AliRsnCutTOFMatch("CutTOFMatch");
  //define PID cuts
  AliRsnCutPIDNSigma *iCutTPCNSigma    = new AliRsnCutPIDNSigma("CutTPCNSigma", fPID, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
  AliRsnCutPIDNSigma *iCutTPCNSigmaElectronRejection = new AliRsnCutPIDNSigma("CutTPCNSigmaElectronRejection", AliPID::kElectron, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
  AliRsnCutPIDNSigma *iCutTPCTOFNSigma = new AliRsnCutPIDNSigma("CutTPCTOFNSigma", fPID, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
  AliRsnCutPIDNSigma *iCutTOFNSigma    = new AliRsnCutPIDNSigma("CutTOFNSigma", fPID, AliRsnCutPIDNSigma::kTOF);//, AliRsnCutPIDNSigma::kP );
  //define phi (azimuthal angle) cuts for TRD presence
  AliRsnCutPhi  *iCutPhiTRD2010        = new AliRsnCutPhi("CutPhiTRD2010","InTRD CheckTOF");
  AliRsnCutPhi  *iCutPhiNoTRD2010      = new AliRsnCutPhi("CutPhiNoTRD2010","OutTRD CheckTOF");
  
  if ((fAppliedCutSetID > AliRsnCutSetDaughterParticle::kTOFMatchPPB2011) && (fAppliedCutSetID < AliRsnCutSetDaughterParticle::kNDaughterCuts)) {  
    if (fNsigmaTPC <= 0.0) {
      AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTPCNSigma->GetName()));
      SetNsigmaForFastTPCpid(10.0);
    }
    if (fNsigmaTOF <= 0.0) {
      AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
      SetNsigmaForFastTOFpid(10.0);
    }
  }
  //
  //defines cut schemes by combining quality cuts and PID cuts
  //
  switch (fAppliedCutSetID)
    {
      case AliRsnCutSetDaughterParticle::kNoCuts :
         AliInfo("No cuts applied to daughter particle");
         break;

      case AliRsnCutSetDaughterParticle::kQualityStd2010 :
         AddCut(fCutQuality);
         SetCutScheme(fCutQuality->GetName());
         break;

      case AliRsnCutSetDaughterParticle::kQualityStd2011:
         AddCut(fCutQuality);
         SetCutScheme(fCutQuality->GetName());
         break;

    case AliRsnCutSetDaughterParticle::kQualityStd2011HighPt:
         AddCut(fCutQuality);
         SetCutScheme(fCutQuality->GetName());
         break;

  case AliRsnCutSetDaughterParticle::kQualityStd2010TRD:
     AddCut(fCutQuality);
     AddCut(iCutPhiTRD2010);
     SetCutScheme( Form("%s&%s",fCutQuality->GetName(),iCutPhiTRD2010->GetName()) );
     break;
     
   case AliRsnCutSetDaughterParticle::kQualityStd2010NoTRD:
     AddCut(fCutQuality);
     AddCut(iCutPhiNoTRD2010);
     SetCutScheme( Form("%s&%s",fCutQuality->GetName(),iCutPhiNoTRD2010->GetName()) );
     break;
   

      case AliRsnCutSetDaughterParticle::kTOFMatch :
         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
         SetCutScheme( Form("%s&(%s)",fCutQuality->GetName(), iCutTOFMatch->GetName()) );
         break;

   case AliRsnCutSetDaughterParticle::kTOFMatchTRD2010 :
     AddCut(fCutQuality);
     AddCut(iCutTOFMatch);
     AddCut(iCutPhiTRD2010);
     SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTOFMatch->GetName(), iCutPhiTRD2010->GetName()) );
     break;
	 
   case AliRsnCutSetDaughterParticle::kTOFMatchNoTRD2010 :
     AddCut(fCutQuality);
     AddCut(iCutTOFMatch);
     AddCut(iCutPhiNoTRD2010);
     SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTOFMatch->GetName(), iCutPhiNoTRD2010->GetName()) );
     break;

   case AliRsnCutSetDaughterParticle::kTOFMatchPPB2011: //pA analysis
     AddCut(fCutQuality);
     AddCut(iCutTOFMatch);
     SetCutScheme( Form("%s&(%s)",fCutQuality->GetName(), iCutTOFMatch->GetName()) );
     break;

      case    AliRsnCutSetDaughterParticle::kFastTPCpidNsigma :
         iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
         AddCut(fCutQuality);
         AddCut(iCutTPCNSigma);
         SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName()) );
         break;

      case    AliRsnCutSetDaughterParticle::kFastTOFpidNsigma :
         iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
         //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
         AddCut(fCutQuality);
         AddCut(iCutTOFNSigma);
         SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTOFNSigma->GetName()) );
         break;

      case    AliRsnCutSetDaughterParticle::kTPCTOFpidKstarPP2010 :
         /* Set TPC  PID (if no TOF)*/
         // all   below   350         MeV: 5sigma
         // all   between 350 and 500 MeV: 3sigma
         // pios above   500         MeV: 2sigma
         // kaons between 500 and 700 MeV: 2sigma
         // kaons above   700         MeV: rejected
         iCutTPCNSigma->AddPIDRange(5.0, 0.0, 0.35);
         iCutTPCNSigma->AddPIDRange(3.0, 0.35, 0.5);
         if (fPID==AliPID::kPion)
            iCutTPCNSigma->AddPIDRange(2.0, 0.5, 1E20);
         if (fPID==AliPID::kKaon)
            iCutTPCNSigma->AddPIDRange(2.0, 0.5, 0.7);

         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
         AddCut(iCutTPCNSigma);

         /* set TPC+TOF PID*/
         iCutTPCTOFNSigma->SinglePIDRange(5.0);
         iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1.5);
         iCutTOFNSigma->AddPIDRange(2.0, 1.5, 1E20);

         AddCut(iCutTPCTOFNSigma);
         AddCut(iCutTOFNSigma);

         // scheme:
         // quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
         SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
         break;

      case AliRsnCutSetDaughterParticle::kTOFpidKstarPbPb2010:
         iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
         //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
         iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal

         AddCut(fCutQuality);
         AddCut(iCutTOFNSigma);
         AddCut(iCutTPCTOFNSigma);
         SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName()) );
         break;

      case AliRsnCutSetDaughterParticle::kTOFTPCmismatchKstarPbPb2010:
         iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal
         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
         AddCut(iCutTPCTOFNSigma);
         SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTOFMatch->GetName(), iCutTPCTOFNSigma->GetName()) );
         break;
           
   case AliRsnCutSetDaughterParticle::kTOFpidKstarPbPbTRD2010:
     iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
     //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
     iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal

     AddCut(fCutQuality);
     AddCut(iCutTOFNSigma);
     AddCut(iCutTPCTOFNSigma);
     AddCut(iCutPhiTRD2010);
     SetCutScheme( Form("%s&%s&%s&%s",fCutQuality->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),iCutPhiTRD2010->GetName()) );
     break;

   case AliRsnCutSetDaughterParticle::kTOFpidKstarPbPbNoTRD2010:
     iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
     //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
     iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal

     AddCut(fCutQuality);
     AddCut(iCutTOFNSigma);
     AddCut(iCutTPCTOFNSigma);
     AddCut(iCutPhiNoTRD2010);
     SetCutScheme( Form("%s&%s&%s&%s",fCutQuality->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),iCutPhiNoTRD2010->GetName()) );
     break;

   case AliRsnCutSetDaughterParticle::kTOFMatchTPCpidNsigma :
     iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
     AddCut(fCutQuality);
     AddCut(iCutTPCNSigma);
     AddCut(iCutTOFMatch);	
     SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName()) );
     break;

   case AliRsnCutSetDaughterParticle::kTPCpidKstarPPB2011:
     iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
     AddCut(fCutQuality);
     AddCut(iCutTPCNSigma);
     SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName()) );
     break;
    
   case AliRsnCutSetDaughterParticle::kTOFpidKstarPPB2011:
     iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
     //iCutTOFNSigma->AddPIDRange(3.0, 0.35, 1E20);
     iCutTPCTOFNSigma->SinglePIDRange(5.0); //5-sigma veto on tpc signal

     AddCut(fCutQuality);
     AddCut(iCutTOFNSigma);
     AddCut(iCutTPCTOFNSigma);
     SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName()) );
     break;

   case AliRsnCutSetDaughterParticle::kTPCTOFpidKstarPPB2011:
     /* Set TPC  PID (if no TOF)*/
     // all   below   350         MeV: 5sigma
     // all   between 350 and 500 MeV: 3sigma
     // pios above   500         MeV: 2sigma
     // kaons between 500 and 700 MeV: 2sigma
     // kaons above   700         MeV: rejected
     iCutTPCNSigma->AddPIDRange(5.0, 0.0, 0.35);
     iCutTPCNSigma->AddPIDRange(3.0, 0.35, 0.5);
     if (fPID==AliPID::kPion)
       iCutTPCNSigma->AddPIDRange(2.0, 0.5, 1E20);
     if (fPID==AliPID::kKaon)
       iCutTPCNSigma->AddPIDRange(2.0, 0.5, 0.7);

     AddCut(fCutQuality);
     AddCut(iCutTOFMatch);
     AddCut(iCutTPCNSigma);

     /* set TPC+TOF PID*/
     iCutTPCTOFNSigma->SinglePIDRange(5.0);
     iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1.5);
     iCutTOFNSigma->AddPIDRange(2.0, 1.5, 1E20);

     AddCut(iCutTPCTOFNSigma);
     AddCut(iCutTOFNSigma);
	 
     // scheme:
     // quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
     SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
     break;
        
           
   case AliRsnCutSetDaughterParticle::kTPCTOFtightPidKStarPPB2011:
     iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
     
     if(fNsigmaTPC==3.0)
       iCutTOFNSigma->SinglePIDRange(5.0);
     if(fNsigmaTPC==2.0)
       iCutTOFNSigma->SinglePIDRange(3.0);
     
     AddCut(fCutQuality);
     AddCut(iCutTPCNSigma);
     AddCut(iCutTOFNSigma);
     
     // scheme:
     // quality & [ (TPCsigma & TOFsigma-2) ]
     SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(),iCutTPCNSigma->GetName(),iCutTOFNSigma->GetName()) ) ;
     break;
     
   case AliRsnCutSetDaughterParticle::kTPCpidMatchPPB2011:
     iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
     AddCut(fCutQuality);
     AddCut(iCutTPCNSigma);
     AddCut(iCutTOFMatch);
     SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName(),iCutTOFMatch->GetName()) );
     break;

    case AliRsnCutSetDaughterParticle::kTPCpidTOFveto4s:
      iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
      AddCut(fCutQuality);
      iCutTOFNSigma->SinglePIDRange(4.0);

      AddCut(iCutTPCNSigma);
      AddCut(iCutTOFMatch);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ (TPCsigma & !TOFmatch) | (TPCsigma & TOFsigma) ]
      SetCutScheme( Form("%s&((%s&(!%s))|(%s&%s))",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(),iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;
       break;

    case AliRsnCutSetDaughterParticle::kTPCpidTOFveto3s:
      iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
      AddCut(fCutQuality);
      iCutTOFNSigma->SinglePIDRange(3.0);

      AddCut(iCutTPCNSigma);
      AddCut(iCutTOFMatch);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ (TPCsigma & !TOFmatch) | (TPCsigma & TOFsigma) ]
      SetCutScheme( Form("%s&((%s&(!%s))|(%s&%s))",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(),iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;
       break;     

    case AliRsnCutSetDaughterParticle::kCombinedPidBestPtDep:
      /* Set TPC  PID (if no TOF)*/
      // all   below  500 MeV: 3sigma
      // all above 500 MeV: 2sigma
      iCutTPCNSigma->AddPIDRange(3.0, 0.0, 0.5);
      iCutTPCNSigma->AddPIDRange(2.0, 0.5, 1E20);
      
      AddCut(fCutQuality);
      AddCut(iCutTOFMatch);
      AddCut(iCutTPCNSigma);
      
      /* set TPC+TOF PID*/
      // pions if TOF match: TPC 5 sigma & TOF 3 sigma
      // kaons if TOF match: 
      //            below 1.5GeV/c: TPC 5 sigma & TOF 3 sigma
      //            above 1.5GeV/c: TPC 3 sigma & TOF 3 sigma

      if (fPID==AliPID::kPion){
	iCutTPCTOFNSigma->SinglePIDRange(5.0);
      }
      
      if (fPID==AliPID::kKaon){
	iCutTPCTOFNSigma->AddPIDRange(5.0, 0.0, 1.5);
	iCutTPCTOFNSigma->AddPIDRange(3.0, 1.5, 1E20);
      }
      iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1E20);
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
	 
     // scheme:
     // quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
     SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
     break;
      
    case AliRsnCutSetDaughterParticle::kTPCPidPtDep:
      /* Set TPC  PID (if no TOF)*/
      // all   below  500 MeV: 3sigma
      // all above 500 MeV: 2sigma
      iCutTPCNSigma->AddPIDRange(3.0, 0.0, 0.5);
      iCutTPCNSigma->AddPIDRange(2.0, 0.5, 1E20);
      
      AddCut(fCutQuality);
      AddCut(iCutTPCNSigma);
      // scheme:
     // quality & TPConly
     SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName()) ) ;
     break;

    case AliRsnCutSetDaughterParticle::kTOFPidPtDep:
      /* Set TOF  PID */
      // all   below  1500 MeV: 3sigma
      // all above 1500 MeV: 2sigma
      //TPC 5 sigma always to remove mismatch
      iCutTPCTOFNSigma->SinglePIDRange(5.0);
      iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1.5);
      iCutTOFNSigma->AddPIDRange(2.0, 1.5, 1E20);
      
      AddCut(fCutQuality);
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
      // scheme:
      // quality & TPConly
      SetCutScheme( Form("%s&%s&%s", fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName()) ) ;
      break;

    case AliRsnCutSetDaughterParticle::kTPCRejPtDepTOFNsigma:
      /* Set TOF  PID */
      // TPC mismatch rejection:
      //             below  1.500 GeV: 5sigma
      //             above  1.500 GeV: 3sigma
      // TOF nsigma PID in full pT
      iCutTPCTOFNSigma->AddPIDRange(5.0, 0.0, 1.5);
      iCutTPCTOFNSigma->AddPIDRange(3.0, 1.5, 1E20);
      iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
      
      AddCut(fCutQuality);
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
      // scheme:
      // quality & TPConly
      SetCutScheme( Form("%s&%s&%s", fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName()) ) ;
      break;
      
    case AliRsnCutSetDaughterParticle::kTPCNsigmaTOFVetoPtDep:
      /* Set TPC  PID */
      // TOF veto:
      //             below  1.500 GeV: 5sigma
      //             above  1.500 GeV: 3sigma
      // TPC nsigma PID in full pT
      iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
      iCutTOFNSigma->AddPIDRange(3.0, 0.0, 1.5);
      iCutTOFNSigma->AddPIDRange(4.0, 1.5, 1E20);
      
      AddCut(fCutQuality);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTOFMatch);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ (TPCsigma & !TOFmatch) | (TPCsigma & TOFsigma) ]
      SetCutScheme( Form("%s&((%s&(!%s))|(%s&%s))",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(),iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;
      break;
      
    case    AliRsnCutSetDaughterParticle::kTPCTOFpidLstar :      
      if (fPID==AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 1.1);
      }
      if (fPID==AliPID::kKaon) {
          iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 0.6);
      }
      
      AddCut(fCutQuality);
      AddCut(iCutTOFMatch);
      AddCut(iCutTPCNSigma);
      
      /* set TPC+TOF PID*/
      iCutTPCTOFNSigma->SinglePIDRange(5.0);
      iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 10);
      
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
      SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
      break;
      
      
      // PID cuts for Lstar analysis at 13 tev pp :
	case    AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV :      
	if (fPID==AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.15, 1.1);
	}
	if (fPID==AliPID::kKaon) {
	//	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 0.6);
	
	iCutTPCNSigma->AddPIDRange(6.,0.15,0.3);
	iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
	}
	
	AddCut(fCutQuality);
	AddCut(iCutTOFMatch);
	AddCut(iCutTPCNSigma);
	
	// set TPC+TOF PID

	iCutTPCTOFNSigma->SinglePIDRange(5.0);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 10);


	AddCut(iCutTPCTOFNSigma);
	AddCut(iCutTOFNSigma);
      
	// scheme:
	// quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
	SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
	break;



	// PID cuts for Lstar analysis at 13 tev pp :
	case    AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV_test1 :      
	if (fPID==AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.15, 1.1);
	}
	if (fPID==AliPID::kKaon) {
	//	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 0.6);
	
	iCutTPCNSigma->AddPIDRange(6.,0.15,0.3);
	iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
	}
	
	AddCut(fCutQuality);
	AddCut(iCutTOFMatch);
	AddCut(iCutTPCNSigma);
	
	// set TPC+TOF PID

	iCutTPCTOFNSigma->SinglePIDRange(3.0);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 10);


	AddCut(iCutTPCTOFNSigma);
	AddCut(iCutTOFNSigma);
      
	// scheme:
	// quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
	SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
	break;





	// PID cuts for Lstar analysis at 13 tev pp :
	case    AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV_test :      
	if (fPID==AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.15, 1.1);
	}
	if (fPID==AliPID::kKaon) {
	//	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 0.6);
	
	iCutTPCNSigma->AddPIDRange(2.5,0.15,0.4);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
	}
	
	AddCut(fCutQuality);
	AddCut(iCutTOFMatch);
	AddCut(iCutTPCNSigma);
	
	// set TPC+TOF PID

	iCutTPCTOFNSigma->SinglePIDRange(5.0);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 10);


	AddCut(iCutTPCTOFNSigma);
	AddCut(iCutTOFNSigma);
      
	// scheme:
	// quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
	SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
	break;


	


	// PID cuts for Lstar analysis at 13 tev pp  with electron rejection:
      case    AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeVERejection :      
	// Set electron rejection cut - 3sigma TPC
	iCutTPCNSigmaElectronRejection->SinglePIDRange(3.0);
	

      if (fPID==AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.15, 1.1);
	}
	if (fPID==AliPID::kKaon) {
	//	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 0.6);
	
	iCutTPCNSigma->AddPIDRange(6.,0.15,0.3);
	iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
	}
	
	AddCut(fCutQuality);
	AddCut(iCutTPCNSigmaElectronRejection);
	AddCut(iCutTOFMatch);
	AddCut(iCutTPCNSigma);
	
	// set TPC+TOF PID

	iCutTPCTOFNSigma->SinglePIDRange(5.0);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 10);


	AddCut(iCutTPCTOFNSigma);
	AddCut(iCutTOFNSigma);
      
	// scheme:
	// quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
	// SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTPCNSigmaElectronRejection->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
            
    SetCutScheme( Form("%s&(!%s)&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCNSigmaElectronRejection->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ; // requested by Priyanka
	break;


	// PID cuts for Lstar analysis at 13 tev pp  with electron rejection:
      case    AliRsnCutSetDaughterParticle::kTPCTOFpidLstar13ppTeV_test2 :      
	// Set electron rejection cut - 3sigma TPC
	iCutTPCNSigmaElectronRejection->SinglePIDRange(3.0);
	

      if (fPID==AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.15, 1.1);
	}
	if (fPID==AliPID::kKaon) {
	//	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 0.6);
	
	iCutTPCNSigma->AddPIDRange(6.,0.15,0.3);
	iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
	}
	
	AddCut(fCutQuality);
	AddCut(iCutTPCNSigmaElectronRejection);
	AddCut(iCutTOFMatch);
	AddCut(iCutTPCNSigma);
	
	// set TPC+TOF PID

	iCutTPCTOFNSigma->SinglePIDRange(3.0);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 10);


	AddCut(iCutTPCTOFNSigma);
	AddCut(iCutTOFNSigma);
      
	// scheme:
	// quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
	// SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTPCNSigmaElectronRejection->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
            
    SetCutScheme( Form("%s&(!%s)&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCNSigmaElectronRejection->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ; // requested by Priyanka
	break;






    case  AliRsnCutSetDaughterParticle::kTPCTOFpidLstarPbPb2011 :

      //Set TPC Nsigma cut
    iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
      /* set mismatch rejection cut for TOF PID*/
    iCutTPCTOFNSigma->SinglePIDRange(fNsigmaTPC);
            
   
    

      //set TOF PID
      if (fPID==AliPID::kProton) {
	iCutTOFNSigma->AddPIDRange(6.0, 0.0, 1.2);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 1.2, 1E20);
      }
      
      if (fPID==AliPID::kKaon) {
	iCutTOFNSigma->AddPIDRange(6.0, 0.0, 0.65);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.65, 1E20);
      }
      
      AddCut(fCutQuality);
      AddCut(iCutTOFMatch);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
      // scheme:
      // quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
         SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
	 
      //      SetCutScheme( Form("%s&((%s&%s)|(%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) );
      break;

    case  AliRsnCutSetDaughterParticle::kTPCTOFpidLstarPbPb2011elRej :

      // Set electron rejection cut - 3sigma TPC
      iCutTPCNSigmaElectronRejection->SinglePIDRange(3.0);

      //Set TPC Nsigma cut
      iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);      

      /* set mismatch rejection cut for TOF PID*/
      iCutTPCTOFNSigma->SinglePIDRange(fNsigmaTPC);

      //set TOF PID
      if (fPID==AliPID::kProton) {
	iCutTOFNSigma->AddPIDRange(6.0, 0.0, 1.2);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 1.2, 1E20);
      }
      
      if (fPID==AliPID::kKaon) {
	iCutTOFNSigma->AddPIDRange(6.0, 0.0, 0.65);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.65, 1E20);
      }

      AddCut(fCutQuality);
      AddCut(iCutTPCNSigmaElectronRejection);
      AddCut(iCutTOFMatch);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
      // scheme:
      // quality & (!electron) & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
        SetCutScheme( Form("%s&(!%s)&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCNSigmaElectronRejection->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;


      //      SetCutScheme( Form("%s&(!%s)&((%s&%s)|(%s))",fCutQuality->GetName(), iCutTPCNSigmaElectronRejection->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;
      break;

    case AliRsnCutSetDaughterParticle::kTPCTOFpidphipp2015:

      if(fNsigmaTOF<10.){
	iCutTPCNSigma->AddPIDRange(6.,0.0,0.3);
	iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
	AddCut(fCutQuality);
	iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);

	AddCut(iCutTPCNSigma);
	AddCut(iCutTOFMatch);
	AddCut(iCutTOFNSigma);

	// scheme:
	// quality & [ (TPCsigma & !TOFmatch) | (TPCsigma & TOFsigma) ]
	SetCutScheme( Form("%s&((%s&(!%s))|(%s&%s))",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(),iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;

      }else{
	iCutTPCNSigma->AddPIDRange(6.,0.0,0.3);
	iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
	AddCut(fCutQuality);
	iCutTOFNSigma->SinglePIDRange(fNsigmaTOF-10.);

	AddCut(iCutTPCNSigma);
	AddCut(iCutTOFMatch);
	AddCut(iCutTOFNSigma);


	iCutTPCTOFNSigma->AddPIDRange(fNsigmaTPC,0.7,1.8);
	AddCut(iCutTPCTOFNSigma);

	// scheme:
	// quality & [ (TPCTOF & TOF) | {!TPCTOF && [(TPCsigma & !TOFmatch) | (TPCsigma & TOFsigma)]} ]
	SetCutScheme( Form("%s&((%s&%s)|((!%s)&((%s&(!%s))|(%s&%s))))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(),iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;
      }

      break;

    case AliRsnCutSetDaughterParticle::kTPCTOFpidphikstarpPb2016:
      
      iCutTPCNSigma->AddPIDRange(6.,0.0,0.3);
      iCutTPCNSigma->AddPIDRange(3.,0.3,0.5);
      iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.5,1.e20);
      AddCut(fCutQuality);
      iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
      
      AddCut(iCutTPCNSigma);
      AddCut(iCutTOFMatch);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ (TPCsigma & !TOFmatch) | (TPCsigma & TOFsigma) ]
      SetCutScheme( Form("%s&((%s&(!%s))|(%s&%s))",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(),iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;
      
      break;

      
    case AliRsnCutSetDaughterParticle::kTPCpidphipp2015:

      iCutTPCNSigma->AddPIDRange(6.,0.0,0.3);
      iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
      iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
      AddCut(fCutQuality);
      AddCut(iCutTPCNSigma);
      SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName()) );
      break;

      
    case AliRsnCutSetDaughterParticle::kTPCTOFpidTunedPbPbTOFneed:
      
      /* 
	 28/04/2016: Neelima Agrawal
	 Roberto Preghenella (preghenella@bo.infn.it)
	 
	 Pb-Pb for pions, kaons and protons
	 tuned to work for both Pb-Pb 2010/2011 data
	 
	 - TPC PID enlarged at 5(7) sigma at low momentum for
	 kaons   < 0.30(0.20) GeV/c
	 protons < 0.50(0.25) GeV/c
	 
	 - TOF PID veto for
	 pions   > 0.40 GeV/c
	 kaon    > 0.45 GeV/c
	 protons > 0.80 GeV
	 
	 - TOF PID required for
	 kaons   > 0.55 GeV/c
	 protons > 1.00 GeV/c
      */
      
      /* pion cuts */
      if (fPID == AliPID::kPion) {
	iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
	//
	iCutTOFNSigma->AddPIDRange(0.00,       0.00, 0.40);  
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.40, 1.e6);  
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0);
      }
      /* kaon cuts */
      if (fPID == AliPID::kKaon) {
	iCutTPCNSigma->AddPIDRange(7.00      , 0.00, 0.20);
	iCutTPCNSigma->AddPIDRange(5.00      , 0.20, 0.30);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.30, 0.55);
	//
	iCutTOFNSigma->AddPIDRange(0.00,       0.00, 0.45);  
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.45, 1.e6);  
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0);
      }
      /* proton cuts */
      if (fPID == AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(7.00      , 0.00, 0.25);
	iCutTPCNSigma->AddPIDRange(5.00      , 0.25, 0.50);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.50, 1.00);
	//
	iCutTOFNSigma->AddPIDRange(0.00,       0.00, 0.80);  
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.80, 1.e6);
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0);
      }
      
      AddCut(fCutQuality);
      AddCut(iCutTOFMatch);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ ( TOFmatch & TOF & TPCTOF ) || ( TPConly ) ]
      SetCutScheme( Form(" %s & ( ( %s & %s & %s ) | ( %s ) )",
			 fCutQuality->GetName(),
			 iCutTOFMatch->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),
			 iCutTPCNSigma->GetName()) ) ;
      break;

    case AliRsnCutSetDaughterParticle::kTPCTOFpidTunedPbPbTOFneed_2018:
      
      /* 
	 20/04/2019: Neelima Agrawal
	 Roberto Preghenella (preghenella@bo.infn.it)
	 
	 Pb-Pb for pions, kaons and protons
	 tuned for Pb-Pb 2010 data
	 
	 - TPC PID with progressive Nsigma
	 Nsigma:       4      3      3      2
	 kaons:   0.00 - 0.20 - 0.30 - 0.40 - 0.60 GeV/c
	 protons: 0.00 - 0.25 - 0.50 - 0.70 - 1.10  GeV/c
	 
	 - TOF PID veto for
	 pions   > 0.40 GeV/c
	 kaon    > 0.45 GeV/c
	 protons > 0.80 GeV
	 
	 - TOF PID required for
	 kaons   > 0.60 GeV/c
	 protons > 1.10 GeV/c
      */
      
      /* pion cuts */
      if (fPID == AliPID::kPion) {
	iCutTPCNSigma->SinglePIDRange(3.0 * fNsigmaTPC);
	//
	iCutTOFNSigma->AddPIDRange(0.00 * fNsigmaTOF, 0.00, 0.40);  
	iCutTOFNSigma->AddPIDRange(3.00 * fNsigmaTOF, 0.40, 1.e6);  
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0 * fNsigmaTPC);
      }
      /* kaon cuts */
      if (fPID == AliPID::kKaon) {
	iCutTPCNSigma->AddPIDRange(4.00 * fNsigmaTPC, 0.00, 0.20);
	iCutTPCNSigma->AddPIDRange(3.00 * fNsigmaTPC, 0.20, 0.40);
	iCutTPCNSigma->AddPIDRange(2.00 * fNsigmaTPC, 0.40, 0.60);
	//
	iCutTOFNSigma->AddPIDRange(0.00 * fNsigmaTOF, 0.00, 0.45);  
	iCutTOFNSigma->AddPIDRange(3.00 * fNsigmaTOF, 0.45, 1.e6);  
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0 * fNsigmaTPC);
      }
      /* proton cuts */
      if (fPID == AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(4.00 * fNsigmaTPC, 0.00, 0.25);
	iCutTPCNSigma->AddPIDRange(3.00 * fNsigmaTPC, 0.25, 0.70);
	iCutTPCNSigma->AddPIDRange(2.00 * fNsigmaTPC, 0.70, 1.10);
	//
	iCutTOFNSigma->AddPIDRange(0.00 * fNsigmaTOF, 0.00, 0.80);  
	iCutTOFNSigma->AddPIDRange(3.00 * fNsigmaTOF, 0.80, 1.e6);
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0 * fNsigmaTPC);
      }
      
      AddCut(fCutQuality);
      AddCut(iCutTOFMatch);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ ( TOFmatch & TOF & TPCTOF ) || ( TPConly ) ]
      SetCutScheme( Form(" %s & ( ( %s & %s & %s ) | ( %s ) )",
			 fCutQuality->GetName(),
			 iCutTOFMatch->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),
			 iCutTPCNSigma->GetName()) ) ;
      break;
      
    case AliRsnCutSetDaughterParticle::kTPCTOFpidTunedPbPbTOFveto:
      
      /* 
	 28/04/2016: Neelima Agrawal
	 Roberto Preghenella (preghenella@bo.infn.it)
	 
	 Pb-Pb for pions, kaons and protons
	 tuned to work for both Pb-Pb 2010/2011 data
	 
	 - TPC PID enlarged at 5(7) sigma at low momentum for
	 kaons   < 0.30(0.20) GeV/c
	 protons < 0.50(0.25) GeV/c
	 
	 - TOF PID veto for
	 pions   > 0.40 GeV/c
	 kaon    > 0.45 GeV/c
	 protons > 0.80 GeV
	 
      */
      
      /* pion cuts */
      if (fPID == AliPID::kPion) {
	iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
	//
	iCutTOFNSigma->AddPIDRange(0.00,       0.00, 0.40);  
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.40, 1.e6);  
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0);
      }
      /* kaon cuts */
      if (fPID == AliPID::kKaon) {
	iCutTPCNSigma->AddPIDRange(7.00      , 0.00, 0.20);
	iCutTPCNSigma->AddPIDRange(5.00      , 0.20, 0.30);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.30, 1.e6);
	//
	iCutTOFNSigma->AddPIDRange(0.00,       0.00, 0.45);  
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.45, 1.e6);  
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0);
      }
      /* proton cuts */
      if (fPID == AliPID::kProton) {
	iCutTPCNSigma->AddPIDRange(7.00      , 0.00, 0.25);
	iCutTPCNSigma->AddPIDRange(5.00      , 0.25, 0.50);
	iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.50, 1.e6);
	//
	iCutTOFNSigma->AddPIDRange(0.00,       0.00, 0.80);  
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.80, 1.e6);
	//
	iCutTPCTOFNSigma->SinglePIDRange(5.0);
      }
      
      AddCut(fCutQuality);
      AddCut(iCutTOFMatch);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTPCTOFNSigma);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ ( TOFmatch & TOF & TPCTOF ) || ( TPConly ) ]
      SetCutScheme( Form(" %s & ( ( %s & %s & %s ) | ( %s ) )",
			 fCutQuality->GetName(),
			 iCutTOFMatch->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName(),
			 iCutTPCNSigma->GetName()) ) ;
      break;

      // PID cuts for Kstar analysis:
      case    AliRsnCutSetDaughterParticle::kTOFTPCpidKstar :      
	if (fPID==AliPID::kKaon) {
	  iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 0.6);
	}
	if (fPID==AliPID::kPion) {
	  iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.0,1.E20);
	}
	
	AddCut(fCutQuality);
	AddCut(iCutTOFMatch);
	AddCut(iCutTPCNSigma);
	
	// set TPC+TOF PID

	iCutTPCTOFNSigma->SinglePIDRange(5.0);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 1.E20);


	AddCut(iCutTPCTOFNSigma);
	AddCut(iCutTOFNSigma);
      
	// scheme:
	// quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
	SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
	break;

      // PID cuts for Delta analysis:
      case    AliRsnCutSetDaughterParticle::kTOFTPCpidDelta :      
	if (fPID==AliPID::kProton) {
	  iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 1.1);
	}
	if (fPID==AliPID::kPion) {
	  iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.0,1.E20);
	}
	
	AddCut(fCutQuality);
	AddCut(iCutTOFMatch);
	AddCut(iCutTPCNSigma);
	
	// set TPC+TOF PID

	iCutTPCTOFNSigma->SinglePIDRange(5.0);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 1.E20);


	AddCut(iCutTPCTOFNSigma);
	AddCut(iCutTOFNSigma);
      
	// scheme:
	// quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
	SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
	break;

      // PID cuts for Lstar analysis:
      case    AliRsnCutSetDaughterParticle::kTOFTPCpidLstar :      
	if (fPID==AliPID::kProton) {
	  iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 1.1);
	}
	if (fPID==AliPID::kKaon) {
	  iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.0, 0.6);
	}
	
	AddCut(fCutQuality);
	AddCut(iCutTOFMatch);
	AddCut(iCutTPCNSigma);
	
	// set TPC+TOF PID

	iCutTPCTOFNSigma->SinglePIDRange(5.0);
	iCutTOFNSigma->AddPIDRange(fNsigmaTOF, 0.0, 1.E20);


	AddCut(iCutTPCTOFNSigma);
	AddCut(iCutTOFNSigma);
      
	// scheme:
	// quality & [ (TOF & TPCTOF) || (!TOFmatch & TPConly) ]
	SetCutScheme( Form("%s&((%s&%s)|((!%s)&%s))",fCutQuality->GetName(), iCutTPCTOFNSigma->GetName(), iCutTOFNSigma->GetName(), iCutTOFMatch->GetName(), iCutTPCNSigma->GetName()) ) ;
	break;


    case AliRsnCutSetDaughterParticle::kTPCTOFvetoPhiXeXe:

      /* PID for Kaon selection for Xe-Xe 2017
	 if TOF, 
	 - Nsigma cut on TOF
	 - 5sigma mismatch rejection with TPC
	 otherwise TPC only:
	 - 6 sigma, p < 0.3 GeV/c
	 - 4 sigma, 0.3 < p < 0.4 GeV/c
	 - nsigma, p > 0.4 GeV/c
      */
      iCutTPCTOFNSigma->SinglePIDRange(5.0);
      
      /*set pt-independent TOF cut*/
      iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);

      /*set pt-dependent tpc cut*/
      iCutTPCNSigma->AddPIDRange(6.,0.0,0.3);
      iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
      iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
      
	
      AddCut(fCutQuality);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTOFMatch);
      AddCut(iCutTOFNSigma);
      AddCut(iCutTPCTOFNSigma);
	
      // scheme:
      // quality & [ ( TPConly & !TOFmatch) || ( TOF & TPCTOF ) ]
      SetCutScheme( Form("%s&((%s&(!%s))|(%s&%s))", fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName()) ) ;
      
      break;

    case AliRsnCutSetDaughterParticle::kTPCTOFPhiXeXe:

      /* PID for Kaon selection for Xe-Xe 2017
	 if TOF, use TOF & TPC
	 else use TPC only
	 
	 TOF cut
	 - Nsigma cut on TOF
	 
	 TPC only:
	 - 6 sigma, p < 0.3 GeV/c
	 - 4 sigma, 0.3 < p < 0.4 GeV/c
	 - nsigma, p > 0.4 GeV/c
      */
      
      /*set pt-independent TOF cut*/
      iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
      
      /*set pt-dependent tpc cut*/
      iCutTPCNSigma->AddPIDRange(6.,0.0,0.3);
      iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
      iCutTPCNSigma->AddPIDRange(fNsigmaTPC, 0.4, 1.e20);
	
      AddCut(fCutQuality);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTOFMatch);
      AddCut(iCutTOFNSigma);
      
      // scheme:
      // quality & [ ( TPConly & !TOFmatch) || (TOF & TPConly ) ]
      SetCutScheme( Form("%s&((%s&(!%s))|(%s&%s))", fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(), iCutTOFNSigma->GetName(), iCutTPCNSigma->GetName()) ) ;
      
      break;

    case AliRsnCutSetDaughterParticle::kTPCTOFvetoElRejPhiXeXe:

      /* PID for Kaon selection for Xe-Xe 2017
	 if TOF, 
	 - Nsigma cut on TOF
	 - 5sigma mismatch rejection with TPC
	 otherwise TPC only:
	 - 6 sigma, p < 0.3 GeV/c
	 - 4 sigma, 0.3 < p < 0.4 GeV/c
	 - nsigma, p > 0.4 GeV/c
	 Electron rejection cut
	 
      */

      // Set electron rejection cut - 3sigma TPC
      iCutTPCNSigmaElectronRejection->SinglePIDRange(3.0);
      
      // Set mismatch rejection cut - 5sigma TPC
      iCutTPCTOFNSigma->SinglePIDRange(5.0);
      
      /*set pt-independent TOF cut*/
      iCutTOFNSigma->SinglePIDRange(fNsigmaTOF);
      
      /*set pt-dependent tpc cut*/
      iCutTPCNSigma->AddPIDRange(6.,0.0,0.3);
      iCutTPCNSigma->AddPIDRange(4.,0.3,0.4);
      iCutTPCNSigma->AddPIDRange(fNsigmaTPC,0.4,1.e20);
      
	
      AddCut(fCutQuality);
      AddCut(iCutTPCNSigma);
      AddCut(iCutTPCNSigmaElectronRejection);
      AddCut(iCutTOFMatch);
      AddCut(iCutTOFNSigma);
      AddCut(iCutTPCTOFNSigma);
	
      // scheme:
      // quality & electronRejection & [ ( TPConly & !TOFmatch ) || (TOF & TPCTOF ) ]
      SetCutScheme( Form("(%s)&(!%s)&((%s&(!%s))|(%s&%s))", fCutQuality->GetName(), iCutTPCNSigmaElectronRejection->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName(), iCutTOFNSigma->GetName(), iCutTPCTOFNSigma->GetName()) ) ;
      
      break;
      
    default :
      break;
    }
  
}


//-----------------------------------------------
void AliRsnCutSetDaughterParticle::PrintTrackQualityCuts()
{
  //Prints track quality cuts
  fCutQuality->Print();
  return;
}

//-----------------------------------------------
void AliRsnCutSetDaughterParticle::InitStdQualityCuts(Bool_t useTPCCrossedRows)
{
  // initialize quality std (if not externally defined) and PID cuts
  // init cut sets by setting variable params
  // 
  // Bool_t isUse2011stdQualityCutsHighPt = ((fAppliedCutSetID==AliRsnCutSetDaughterParticle::kQualityStd2011HighPt) ||
  // 					  (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kTOFMatchPPB2011) || //pA analysis
  // 					  (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kTPCpidKstarPPB2011) ||
  // 					  (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kTOFpidKstarPPB2011) ||
  // 					  (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kTPCTOFpidKstarPPB2011) ||
  // 					  (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kTPCpidTOFvetoKStarPPB2011) ||
  // 					  (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kTPCpidMatchPPB2011));
  
  // Bool_t isUse2011stdQualityCuts = (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kQualityStd2011);

  if (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kQualityStd2011) {
    fIsUse2011stdQualityCuts = kTRUE;
    fIsUse2011stdQualityCutsHighPt = kFALSE;
  } 
  if (fAppliedCutSetID==AliRsnCutSetDaughterParticle::kQualityStd2011HighPt) {
    fIsUse2011stdQualityCuts = kFALSE;
    fIsUse2011stdQualityCutsHighPt = kTRUE;
  }
  
  if (fIsUse2011stdQualityCuts) {
    AliInfo(Form("Using 2011 std quality cuts with cut on TPC %s",(useTPCCrossedRows?"crossed rows":"N clusters")));
    fCutQuality->SetDefaults2011(useTPCCrossedRows, kFALSE);//uses filter bit 5 as default
  } else {
    if (fIsUse2011stdQualityCutsHighPt) {
      AliInfo(Form("Using 2011 std quality cuts with cut on TPC %s for high-pT", (useTPCCrossedRows?"crossed rows":"N clusters")));
      fCutQuality->SetDefaultsHighPt2011(useTPCCrossedRows, kFALSE);//uses filter bit 10 as default
    } else {
      AliInfo(Form("Using 2010 std quality cuts with cut on TPC %s", (useTPCCrossedRows?"crossed rows":"N clusters")));
      fCutQuality->SetDefaults2010(useTPCCrossedRows, kFALSE);
    }
  }
  fCutQuality->SetAODTestFilterBit(fAODTrkCutFilterBit); //changes default filter bit to the chosen filter bit
  AliInfo(Form("Applying cut on AOD filter bit %i", fAODTrkCutFilterBit));
  //apply pt and eta cuts
  fCutQuality->SetPtRange(fPtRange[0], fPtRange[1]);
  fCutQuality->SetEtaRange(fEtaRange[0], fEtaRange[1]);
  AliInfo(Form("Pt range [%3.2f,%3.2f], Eta range [%3.2f, %3.2f]", fPtRange[0], fPtRange[1], fEtaRange[0], fEtaRange[1]));
  return;
}
