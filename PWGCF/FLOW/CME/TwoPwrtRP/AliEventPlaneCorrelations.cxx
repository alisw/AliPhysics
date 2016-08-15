
#include "AliEventPlaneCorrelations.h"
#include "AliChargeOnePwrtRP.h"
#include "AliCMEVarManager.h"

ClassImp(AliEventPlaneCorrelations)



//____________________________________________________________________________
AliEventPlaneCorrelations::AliEventPlaneCorrelations(TAxis* ax,TString epA, TString epB, TString epC) :
  TObject(),
  fEPcorrelation(),
  fEventPlanes(),
  fCorrelationName()
{

  fCorrelationName = epA+"x"+epB;
  if(epC!="") fCorrelationName+="x"+epC;


  for(Int_t i=0; i<3; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nharmonics; j++) for(Int_t k=0; k<4; k++) fEPcorrelation[i][j][k]=0x0;
  for(Int_t i=0; i<3; i++) fEventPlanes[i]=0x0;


  TString ep[3] = {epA,epB,epC};

  for(Int_t ih=1; ih<=AliChargeOnePwrtRP::Nharmonics; ++ih){
    for(Int_t iep=0; iep<3; iep++){

      if(ep[iep]==""||ep[(iep+1)%3]=="") continue;

      //fEPcorrelation[iep][ih-1][0] = new TProfile(Form("XX_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),Form("XX_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),100,0,100);
      //fEPcorrelation[iep][ih-1][1] = new TProfile(Form("YY_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),Form("YY_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),100,0,100);
      //fEPcorrelation[iep][ih-1][2] = new TProfile(Form("XY_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),Form("XY_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),100,0,100);
      //fEPcorrelation[iep][ih-1][3] = new TProfile(Form("YX_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),Form("YX_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),100,0,100);
      fEPcorrelation[iep][ih-1][0] = new TProfile(Form("XX_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),Form("XX_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),ax->GetNbins(),ax->GetXbins()->GetArray());
      fEPcorrelation[iep][ih-1][1] = new TProfile(Form("YY_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),Form("YY_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),ax->GetNbins(),ax->GetXbins()->GetArray());
      fEPcorrelation[iep][ih-1][2] = new TProfile(Form("XY_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),Form("XY_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),ax->GetNbins(),ax->GetXbins()->GetArray());
      fEPcorrelation[iep][ih-1][3] = new TProfile(Form("YX_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),Form("YX_%sx%s_h%d", ep[iep].Data(), ep[(iep+1)%3].Data(), ih),ax->GetNbins(),ax->GetXbins()->GetArray());



    }
  }
}


//_______________________________________________________________________________
AliEventPlaneCorrelations::~AliEventPlaneCorrelations()
{
  //
  // De-Constructor
  //
  for(Int_t ih=1; ih<=AliChargeOnePwrtRP::Nharmonics; ++ih) for(Int_t iep=0; iep<3; iep++) for(Int_t icor=0; icor<4; icor++) if(fEPcorrelation[iep][ih-1][icor]) delete fEPcorrelation[iep][ih-1][icor];
}




//_________________________________________________________________________
void AliEventPlaneCorrelations::FillCorrelations(Float_t x){

  for(Int_t ih=1; ih<=AliChargeOnePwrtRP::Nharmonics; ++ih){
    for(Int_t iep=0; iep<3; iep++){

      if(!fEventPlanes[iep]||!fEventPlanes[(iep+1)%3]) continue;
      //cout<<fEventPlanes[iep]->N()<<"   "<<fEventPlanes[iep]->Qx(2)<<endl;

      /* Adapted to the new Qn vector Correction framework */
      if(fEventPlanes[iep]->GetN()==0||fEventPlanes[(iep+1)%3]->GetN()==0) continue;

      //if(iep==1) cout<<CorrelationName()<<"  "<<iep<<"  "<<fEventPlanes[iep]->QxNorm(ih)<<"  "<<fEventPlanes[(iep+1)%3]->QxNorm(ih)<<"   "<<fEPcorrelation[iep][ih-1][0]->GetMean(2)<<endl;
      //cout<<qvecC->N()<<endl;

      fEPcorrelation[iep][ih-1][0]->Fill(x,fEventPlanes[iep]->QxNorm(ih)*fEventPlanes[(iep+1)%3]->QxNorm(ih));
      fEPcorrelation[iep][ih-1][1]->Fill(x,fEventPlanes[iep]->QyNorm(ih)*fEventPlanes[(iep+1)%3]->QyNorm(ih));
      fEPcorrelation[iep][ih-1][2]->Fill(x,fEventPlanes[iep]->QxNorm(ih)*fEventPlanes[(iep+1)%3]->QyNorm(ih));
      fEPcorrelation[iep][ih-1][3]->Fill(x,fEventPlanes[iep]->QyNorm(ih)*fEventPlanes[(iep+1)%3]->QxNorm(ih));

    }
  }


}
