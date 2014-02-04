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
AliRsnCutSetDaughterParticle::AliRsnCutSetDaughterParticle(const char *name, AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID, AliPID::EParticleType pid, Float_t nSigmaFast = -1.0, Int_t AODfilterBit = 0) :
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
  InitStdQualityCuts();
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
    //sets default track quality to be initialised +
    //sets here pt and eta cuts
    InitStdQualityCuts();
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
  AliRsnCutPIDNSigma *iCutTPCTOFNSigma = new AliRsnCutPIDNSigma("CutTPCTOFNSigma", fPID, AliRsnCutPIDNSigma::kTPC);//, AliRsnCutPIDNSigma::kTPCinnerP );
  AliRsnCutPIDNSigma *iCutTOFNSigma    = new AliRsnCutPIDNSigma("CutTOFNSigma", fPID, AliRsnCutPIDNSigma::kTOF);//, AliRsnCutPIDNSigma::kP );
  //define phi (azimuthal angle) cuts for TRD presence
  AliRsnCutPhi  *iCutPhiTRD2010        = new AliRsnCutPhi("CutPhiTRD2010","InTRD CheckTOF");
  AliRsnCutPhi  *iCutPhiNoTRD2010      = new AliRsnCutPhi("CutPhiNoTRD2010","OutTRD CheckTOF");
  
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

      case AliRsnCutSetDaughterParticle::kTOFMatch :
         AddCut(fCutQuality);
         AddCut(iCutTOFMatch);
         SetCutScheme( Form("%s&(%s)",fCutQuality->GetName(), iCutTOFMatch->GetName()) );
         break;

      case    AliRsnCutSetDaughterParticle::kFastTPCpidNsigma :
         if (fNsigmaTPC <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTPCNSigma->GetName()));
            SetNsigmaForFastTPCpid(10.0);
         }
         iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
         AddCut(fCutQuality);
         AddCut(iCutTPCNSigma);
         SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName()) );
         break;

      case    AliRsnCutSetDaughterParticle::kFastTOFpidNsigma :
         if (fNsigmaTOF <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
            SetNsigmaForFastTOFpid(10.0);
         }
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
         if (fNsigmaTOF <= 0.0) {
            AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
            SetNsigmaForFastTOFpid(10.0);
         }
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

   case AliRsnCutSetDaughterParticle::kTOFpidKstarPbPbTRD2010:
     if (fNsigmaTOF <= 0.0) {
       AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
       SetNsigmaForFastTOFpid(10.0);
     }
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
     if (fNsigmaTOF <= 0.0) {
       AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
       SetNsigmaForFastTOFpid(10.0);
     }
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
     if (fNsigmaTPC <= 0.0) {
       AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTPCNSigma->GetName()));
       SetNsigmaForFastTPCpid(10.0);
     }
     iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
     AddCut(fCutQuality);
     AddCut(iCutTPCNSigma);
     AddCut(iCutTOFMatch);	
     SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName(), iCutTOFMatch->GetName()) );
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
   
   case AliRsnCutSetDaughterParticle::kTOFMatchPPB2011: //pA analysis
     AddCut(fCutQuality);
     AddCut(iCutTOFMatch);
     SetCutScheme( Form("%s&(%s)",fCutQuality->GetName(), iCutTOFMatch->GetName()) );
     break;

   case AliRsnCutSetDaughterParticle::kTPCpidKstarPPB2011:
     if (fNsigmaTPC <= 0.0) {
       AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTPCNSigma->GetName()));
       SetNsigmaForFastTPCpid(10.0);
     }
     iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
     AddCut(fCutQuality);
     AddCut(iCutTPCNSigma);
     SetCutScheme( Form("%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName()) );
     break;
    
   case AliRsnCutSetDaughterParticle::kTOFpidKstarPPB2011:
     if (fNsigmaTOF <= 0.0) {
       AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
       SetNsigmaForFastTOFpid(10.0);
     }
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
     if (fNsigmaTPC <= 0.0) {
       AliWarning(Form("Invalid number of sigmas required for %s. Setting nSigma = 5.0",iCutTPCNSigma->GetName()));
       SetNsigmaForFastTPCpid(5.0);
     }
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
     if (fNsigmaTPC <= 0.0) {
       AliWarning(Form("Invalid number of sigmas required for %s. Setting nSigma = 5.0",iCutTPCNSigma->GetName()));
       SetNsigmaForFastTPCpid(5.0);
     }
     iCutTPCNSigma->SinglePIDRange(fNsigmaTPC);
     AddCut(fCutQuality);
     AddCut(iCutTPCNSigma);
     AddCut(iCutTOFMatch);
     SetCutScheme( Form("%s&%s&%s",fCutQuality->GetName(), iCutTPCNSigma->GetName(),iCutTOFMatch->GetName()) );
     break;

    case AliRsnCutSetDaughterParticle::kTPCpidTOFveto4s:
      if (fNsigmaTPC <= 0.0) {
	AliWarning(Form("Invalid number of sigmas required for %s. Setting nSigma = 5.0",iCutTPCNSigma->GetName()));
	SetNsigmaForFastTPCpid(5.0);
      }
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
      if (fNsigmaTPC <= 0.0) {
	AliWarning(Form("Invalid number of sigmas required for %s. Setting nSigma = 5.0",iCutTPCNSigma->GetName()));
	SetNsigmaForFastTPCpid(5.0);
      }
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
      
      if (fNsigmaTOF <= 0.0) {
	AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTOFNSigma->GetName()));
	SetNsigmaForFastTOFpid(10.0);
      }
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
      if (fNsigmaTPC <= 0.0) {
	AliWarning(Form("Invalid number of sigmas required for %s. Setting default nSigma = 10",iCutTPCNSigma->GetName()));
	SetNsigmaForFastTPCpid(10.0);
      }
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
void AliRsnCutSetDaughterParticle::InitStdQualityCuts()
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
    fCutQuality->SetDefaults2011();//uses filter bit 5 as default
  } else {
    if (fIsUse2011stdQualityCutsHighPt) {
      fCutQuality->SetDefaultsHighPt2011();//uses filter bit 5 as default
    } else {
      fCutQuality->SetDefaults2010();
      fCutQuality->SetDCARPtFormula("0.0182+0.0350/pt^1.01");
      fCutQuality->SetDCAZmax(2.0);
      fCutQuality->SetSPDminNClusters(1);
      fCutQuality->SetITSminNClusters(0);
      fCutQuality->SetITSmaxChi2(36);
      fCutQuality->SetTPCminNClusters(70);
      fCutQuality->SetTPCmaxChi2(4.0);
      fCutQuality->SetRejectKinkDaughters();
      //fCutQuality->SetITSmaxChi2(36);
      //fCutQuality->SetMaxChi2TPCConstrainedGlobal(36);
    }
  }
  fCutQuality->SetAODTestFilterBit(fAODTrkCutFilterBit); //changes default filter bit to the chosen filter bit
  
  //apply pt and eta cuts
  fCutQuality->SetPtRange(fPtRange[0], fPtRange[1]);
  fCutQuality->SetEtaRange(fEtaRange[0], fEtaRange[1]);
  AliInfo("Standard quality cuts applied");
  fCutQuality->Print();
  return;
}
