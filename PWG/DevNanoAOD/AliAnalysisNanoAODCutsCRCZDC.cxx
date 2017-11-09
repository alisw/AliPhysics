#include "AliAnalysisNanoAODCutsCRCZDC.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"
#include "AliMultSelection.h"
#include "AliVVertex.h"
#include "AliAnalysisUtils.h"
#include "TMatrixDSym.h"
#include "AliAODTracklets.h"
#include <iomanip>

ClassImp(AliAnalysisNanoAODTrackCutsCRCZDC)
ClassImp(AliAnalysisNanoAODEventCutsCRCZDC)
ClassImp(AliNanoAODSimpleSetterCRCZDC)


AliAnalysisNanoAODTrackCutsCRCZDC::AliAnalysisNanoAODTrackCutsCRCZDC():
AliAnalysisCuts(), fBitMask(1), fMinPt(0), fMaxEta(10)
{
  // default ctor 
}

Bool_t AliAnalysisNanoAODTrackCutsCRCZDC::IsSelected(TObject* obj)
{
  // Returns true if the track is good!
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(obj);
  
  
  if (!track->TestFilterBit(fBitMask))    return kFALSE;
  if (track->Pt() < fMinPt)               return kFALSE;
  if (TMath::Abs(track->Eta()) > fMaxEta) return kFALSE; 
  
  return kTRUE;  

}

AliAnalysisNanoAODEventCutsCRCZDC::AliAnalysisNanoAODEventCutsCRCZDC():
  AliAnalysisCuts(), 
  fVertexRange(-1),
  fTrackCut(0),
  fMinMultiplicity(-1),
  fMaxMultiplicity(-1)
{
  // default ctor   
}

Bool_t AliAnalysisNanoAODEventCutsCRCZDC::IsSelected(TObject* obj)
{
  // Returns true if object accepted on the event level
  
  AliAODEvent * evt = dynamic_cast<AliAODEvent*>(obj);
  
  if (fVertexRange > 0)
  {
    AliAODVertex * vertex = evt->GetPrimaryVertex();
    if (!vertex) 
      return kFALSE;
    
    if (vertex->GetNContributors() < 1) 
    {
      // SPD vertex cut
      vertex = evt->GetPrimaryVertexSPD();    
      if (!vertex || vertex->GetNContributors() < 1) 
        return kFALSE;
    }    
    
    if (TMath::Abs(vertex->GetZ()) > fVertexRange) 
      return kFALSE;
  }
  
  if (fTrackCut != 0)
  {
    Int_t trackCount = 0;
    for (Int_t j=0; j<evt->GetNumberOfTracks(); j++)
      if (fTrackCut->IsSelected(evt->GetTrack(j)))
        trackCount++;
      
    if (fMinMultiplicity > 0 && trackCount < fMinMultiplicity)
      return kFALSE;
    if (fMaxMultiplicity > 0 && trackCount > fMaxMultiplicity)
      return kFALSE;
  }
      
  return kTRUE;
}


void AliNanoAODSimpleSetterCRCZDC::SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head, TString varListHeader   ) {

  AliAODHeader * header = dynamic_cast<AliAODHeader*>(event->GetHeader());
  if (!header) AliFatal("Not a standard AOD");

  // Set custom nano aod vars
  Double_t centrV0M=-1;
  Double_t centrTRK=-1;
  Double_t centrCL1=-1;
  Double_t centrCL0=-1;

      //2015 multiplicity selection
      AliMultSelection *MultSelection = 0x0; 
      MultSelection = (AliMultSelection *) event->FindListObject("MultSelection");

      if(MultSelection){
        centrV0M = MultSelection->GetMultiplicityPercentile("V0M");
        centrCL1 = MultSelection->GetMultiplicityPercentile("CL1");
        centrCL0 = MultSelection->GetMultiplicityPercentile("CL0");
        centrTRK = MultSelection->GetMultiplicityPercentile("TRK");
     
      }else{
        //2011 
        AliCentrality * centralityObj = header->GetCentralityP();
        centrV0M = centralityObj->GetCentralityPercentile("V0M");
        centrTRK = centralityObj->GetCentralityPercentile("TRK");
        centrCL1 = centralityObj->GetCentralityPercentile("CL1");
        centrCL0 = centralityObj->GetCentralityPercentile("CL0");

      }

  Double_t magfield = header->GetMagneticField();
  Int_t runNumber = event->GetRunNumber();

    static const char * validatorString[] = {"Centr","MagField","CentrTRK","CentrCL0", "CentrCL1", "RunNumber", 0};
  TObjArray * vars = varListHeader.Tokenize(",");
  //Int_t size = vars->GetSize();
  TIter it(vars);
  TObjString *token  = 0;
  Int_t index=0;

  std::map<TString,int> cstMap = head->GetMapCstVar();

  while ((token = (TObjString*) it.Next())) {
    TString var = token->GetString().Strip(TString::kBoth, ' ');

    // Check if string is in the allowed list
    Bool_t isValid = kFALSE;
    Int_t ivalidator = 0;
    while (validatorString[ivalidator]) {
      if(var == validatorString[ivalidator++]) isValid = kTRUE;
    }

    if (!( isValid || var.BeginsWith("cst"))) AliFatal(Form("Invalid var [%s]", var.Data()));
    if     (var == "Centr"      ) head->SetCentrIndex      (index);
    else if(var == "CentrTRK"   ) head->SetCentrTRKIndex   (index);
    else if(var == "CentrCL0"   ) head->SetCentrCL0Index   (index);
    else if(var == "CentrCL1"   ) head->SetCentrCL1Index   (index);
    else if(var == "MagField"   ) head->SetMagFieldIndex   (index);
    else if(var == "RunNumber"  ) head->SetRunNumberIndex  (index);
    else {
      cstMap[var] = index;
      std::cout << "ADDING " << index << " " << cstMap[var] << " " << var.Data() << std::endl;
      
    }

    index++;
  }
  //size = index;
  if(vars) vars->Delete();
  head->SetMapCstVar(cstMap);
 
  if ((head->GetCentrIndex())!=-1)     head->SetVar(head->GetCentrIndex()    ,           centrV0M );
  if ((head->GetCentrTRKIndex())!=-1)  head->SetVar(head->GetCentrTRKIndex() ,           centrTRK );
  if ((head->GetCentrCL1Index())!=-1)  head->SetVar(head->GetCentrCL1Index() ,           centrCL1 );
  if ((head->GetCentrCL0Index())!=-1)  head->SetVar(head->GetCentrCL0Index() ,           centrCL0 );
  if ((head->GetMagFieldIndex())!=-1)  head->SetVar(head->GetMagFieldIndex() ,           magfield );
  if ((head->GetRunNumberIndex())!=-1) head->SetVar(head->GetRunNumberIndex(), Double_t(runNumber));

  //Pile-up info
  AliAODEvent *eventCopy= const_cast <AliAODEvent*> (event);
  static Int_t pileUpIndex = head ->GetVarIndex("cstPileUp");
  if(SelectPileup(eventCopy)) head->SetVar(pileUpIndex,1.);
  if(!SelectPileup(eventCopy)) head->SetVar(pileUpIndex,0.);

  Double_t SumV0=0.;
  for(Int_t i=0; i<64; i++) {
    if(std::isfinite(event->GetVZEROEqMultiplicity(i))) SumV0 += event->GetVZEROEqMultiplicity(i);
  }
  static Int_t V0Index = head->GetVarIndex("cstV0"); 
  head->SetVar(V0Index,SumV0);

  UInt_t period = event->GetPeriodNumber();
  UInt_t orbit24 = event->GetOrbitNumber();

  static Int_t periodIndex = head->GetVarIndex("cstPeriod"); 
  static Int_t orbitIndex = head->GetVarIndex("cstOrbit"); 

  head->SetVar(periodIndex, Double_t(period));
  head->SetVar(orbitIndex, Double_t(orbit24));

  AliAODTracklets *trackl = event->GetTracklets();
  Int_t nTracklets = trackl->GetNumberOfTracklets();
  static Int_t trackeletsIndex = head->GetVarIndex("cstNTrackelets"); 

  head->SetVar(trackeletsIndex, Double_t(nTracklets));

  static Int_t energyZNCIndex = head->GetVarIndex("cstEnergyZNC"); 
  Double_t energyZNC = header->GetZDCN1Energy();
  static Int_t energyZNAIndex = head->GetVarIndex("cstEnergyZNA"); 
  Double_t energyZNA = header->GetZDCN2Energy();
  static Int_t energyZPCIndex = head->GetVarIndex("cstEnergyZPC"); 
  Double_t energyZPC = header->GetZDCP1Energy();
  static Int_t energyZPAIndex = head->GetVarIndex("cstEnergyZPA"); 
  Double_t energyZPA = header->GetZDCP2Energy();

  head->SetVar(energyZNCIndex, energyZNC);
  head->SetVar(energyZNAIndex, energyZNA);
  head->SetVar(energyZPCIndex, energyZPC);
  head->SetVar(energyZPAIndex, energyZPA);
}

Bool_t AliNanoAODSimpleSetterCRCZDC::SelectPileup(AliAODEvent *aod)
{
  Bool_t BisPileup=kFALSE;

  Double_t centrV0M=300., centrCL1=300.;

  AliAnalysisUtils* utils = new AliAnalysisUtils();
  utils->SetUseMVPlpSelection(kTRUE);
  utils->SetUseOutOfBunchPileUp(kTRUE);

  AliMultSelection *MultSelection = 0x0; 
  MultSelection = (AliMultSelection *) aod->FindListObject("MultSelection");

  if(!MultSelection){

    // pileup for LHC10h and LHC11h
    centrV0M = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("V0M");
    centrCL1 = ((AliVAODHeader*)aod->GetHeader())->GetCentralityP()->GetCentralityPercentile("CL1");
    
    // check anyway pileup
    if (plpMV(aod)) {
      BisPileup=kTRUE;
    }
    
    Short_t isPileup = aod->IsPileupFromSPD(3);
    if (isPileup != 0) {
      BisPileup=kTRUE;
    }
    
    if (((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
      BisPileup=kTRUE;
    }
    
    if (aod->IsIncompleteDAQ())  {
      BisPileup=kTRUE;
    }
    
    // check vertex consistency
    const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
    const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();
    
    if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)  {
      BisPileup=kTRUE;
    }
    
    double covTrc[6], covSPD[6];
    vtTrc->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);
    
    double dz = vtTrc->GetZ() - vtSPD->GetZ();
    
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = dz/errTot;
    double nsigTrc = dz/errTrc;
    
    if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
      BisPileup=kTRUE;
    }

    
  /*  if (utils->IsPileUpEvent(aod)) {
      BisPileup=kTRUE;
    }*/
    
  }
  else {
    
    // pileup for LHC15o, using AliMultSelection
    
    if(MultSelection) {
      centrV0M = MultSelection->GetMultiplicityPercentile("V0M");
      centrCL1 = MultSelection->GetMultiplicityPercentile("CL1");
    } else {
      BisPileup=kTRUE;
    }
    
    // pile-up a la Dobrin for LHC15o
    if (plpMV(aod)) {
      BisPileup=kTRUE;
    }
    
    Short_t isPileup = aod->IsPileupFromSPD(3);
    if (isPileup != 0) {
      BisPileup=kTRUE;
    }
    
    if (((AliAODHeader*)aod->GetHeader())->GetRefMultiplicityComb08() < 0) {
      BisPileup=kTRUE;
    }
    
    if (aod->IsIncompleteDAQ())  {
      BisPileup=kTRUE;
    }
    
    if(fabs(centrV0M-centrCL1)>7.5)  {
      BisPileup=kTRUE;
    }
    
    // check vertex consistency
    const AliAODVertex* vtTrc = aod->GetPrimaryVertex();
    const AliAODVertex* vtSPD = aod->GetPrimaryVertexSPD();
    
    if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)  {
      BisPileup=kTRUE;
    }
    
    double covTrc[6], covSPD[6];
    vtTrc->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);
    
    double dz = vtTrc->GetZ() - vtSPD->GetZ();
    
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = dz/errTot;
    double nsigTrc = dz/errTrc;
    
    if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)  {
      BisPileup=kTRUE;
    }
    
    // cuts on tracks
    const Int_t nTracks = aod->GetNumberOfTracks();
    Int_t multEsd = ((AliAODHeader*)aod->GetHeader())->GetNumberOfESDTracks();
    
    //Int_t multTrk = 0;
    //Int_t multTrkBefC = 0;
    //Int_t multTrkTOFBefC = 0;
    Int_t multTPC = 0;
    
    for (Int_t it = 0; it < nTracks; it++) {
      
      AliAODTrack* aodTrk = (AliAODTrack*)aod->GetTrack(it);
      if (!aodTrk){
        delete aodTrk;
        continue;
      }
      
//      if (aodTrk->TestFilterBit(32)){
//        multTrkBefC++;
//        
//        if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
//          multTrkTOFBefC++;
//        
//        if ((TMath::Abs(aodTrk->Eta()) < 0.8) && (aodTrk->GetTPCNcls() >= 70) && (aodTrk->Pt() >= 0.2) && (aodTrk->Pt() < 20.))
//          multTrk++;
//      }
      
      if (aodTrk->TestFilterBit(128))
        multTPC++;
      
    } // end of for (Int_t it = 0; it < nTracks; it++)
    
    Double_t multTPCn = multTPC;
    Double_t multEsdn = multEsd;
    Double_t multESDTPCDif = multEsdn - multTPCn*3.38;
    
    if (multESDTPCDif > (fRejectPileUpTight?700.:15000.)) { //fRejectPileUpTight???
        BisPileup=kTRUE;
    }

    if(fRejectPileUpTight) {
      if(BisPileup==kFALSE) {
        if(!MultSelection->GetThisEventIsNotPileup()) BisPileup=kTRUE;
        if(!MultSelection->GetThisEventIsNotPileupMV()) BisPileup=kTRUE;
        if(!MultSelection->GetThisEventIsNotPileupInMultBins()) BisPileup=kTRUE;
        if(!MultSelection->GetThisEventHasNoInconsistentVertices()) BisPileup=kTRUE;
        if(!MultSelection->GetThisEventPassesTrackletVsCluster()) BisPileup=kTRUE;
        if(!MultSelection->GetThisEventIsNotIncompleteDAQ()) BisPileup=kTRUE;
        if(!MultSelection->GetThisEventHasGoodVertex2016()) BisPileup=kTRUE;
      }
    }
  }
  
  return BisPileup;
}

Bool_t AliNanoAODSimpleSetterCRCZDC::plpMV(const AliAODEvent* aod)
{
  // check for multi-vertexer pile-up
  
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2 = 5.0;
  const double kMinWDist = 15;
  
  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;
  int nPlp = 0;
  
  if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
  vtPrm = aod->GetPrimaryVertex();
  if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary
  
  //int bcPrim = vtPrm->GetBC();
  
  for (int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = (const AliVVertex*)aod->GetPileupVertexTracks(ipl);
    //
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF() > kMaxPlpChi2) continue;
    //  int bcPlp = vtPlp->GetBC();
    //  if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
    //
    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist) continue;
    //
    return kTRUE; // pile-up: well separated vertices
  }
  
  return kFALSE;
}

Double_t AliNanoAODSimpleSetterCRCZDC::GetWDist(const AliVVertex* v0, const AliVVertex* v1)
{
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    printf("One of vertices is not valid\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  double dist = -1;
  double dx = v0->GetX()-v1->GetX();
  double dy = v0->GetY()-v1->GetY();
  double dz = v0->GetZ()-v1->GetZ();
  double cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
  +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}
