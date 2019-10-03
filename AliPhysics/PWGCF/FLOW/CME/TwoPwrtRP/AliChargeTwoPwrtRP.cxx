
#include "AliChargeTwoPwrtRP.h"
#include "AliChargeOnePwrtRP.h"
#include "AliCMEVarManager.h"

ClassImp(AliChargeTwoPwrtRP)


UChar_t AliChargeTwoPwrtRP::fCharge=0;
UShort_t AliChargeTwoPwrtRP::fBins[]={0};
Float_t AliChargeTwoPwrtRP::fPair[][2]={{0.}};

//____________________________________________________________________________
AliChargeTwoPwrtRP::AliChargeTwoPwrtRP() :
  TObject(),
  fMult(),
  f2pCorrelationShort(),
  f2pCorrelationShort2(),
  f3pCorrelationShort(),
  f3pCorrelationShort2(),
  fNbinsX(0),
  fNbinsY(0),
  fNbinsZ(0),
  fNbins(0),
  fNdim(0),
  fVarX(0),
  fVarY(0),
  fVarZ(0),
  fTrackVarX(0),
  fTrackVarY(0),
  fTrackVarZ(0),
  fTrackVarMap(),
  //fVarNames(),
  fBinLimits(),
  fEventSelected(kFALSE),
  f2pCorrelationTHn(),
  f3pCorrelationTHn(),
  fCorrelationTHnMult(),
  fEventSelection(""),
  fTracking1Quality(""),
  fTracking2Quality(""),
  fPID1(""),
  fPID2(""),
  fEventPlanes(),
  fXaxisLabel(""),
  fYaxisLabel("")
{
  //fVarNames[0]="1Pt";
  //fVarNames[1]="mPt";
  //fVarNames[2]="dPt";
  //fVarNames[3]="1Eta";
  //fVarNames[4]="mEta";
  //fVarNames[5]="dEta";

  TString charges[5] = {"aa","pp","mm","pm","mp"};
  TString PIDnames[4] = {"c","pi","K","p"};

  for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) for(Int_t j=0; j<2; j++){f2pCorrelationTHn[g][i][j]=0x0;}
  for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) for(Int_t k=0; k<2; k++) {f3pCorrelationTHn[g][i][j][k]=0x0;}
  for(Int_t g=0; g<4; g++) {fCorrelationTHnMult[g]=0x0;}
  for(Int_t g=0; g<3; g++) {fTrackVarMap[g]=-1;}

}

//_______________________________________________________________________________
AliChargeTwoPwrtRP::~AliChargeTwoPwrtRP()
{
  //
  // De-Constructor
  //
  for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) delete [] f2pCorrelationShort[g][i];
  for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) delete [] f3pCorrelationShort[g][i][j];
  for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) delete [] f2pCorrelationShort2[g][i];
  for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) delete [] f3pCorrelationShort2[g][i][j];

  for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) for(Int_t j=0; j<2; j++){delete f2pCorrelationTHn[g][i][j];}
  for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) for(Int_t k=0; k<2; k++) {delete f3pCorrelationTHn[g][i][j][k];}
  for(Int_t g=0; g<4; g++) {delete fCorrelationTHnMult[g];}
  //for(Int_t g=0; g<4; g++) for(Int_t i=0; i<7; i++) for(Int_t j=0; j<4; j++) delete [] f2pCorrelation[g][i][j];
  //for(Int_t g=0; g<4; g++) for(Int_t i=0; i<6; i++) for(Int_t j=0; j<8; j++) delete [] f3pCorrelation[g][i][j];
  //for(Int_t g=0; g<4; g++) for(Int_t i=0; i<7; i++) for(Int_t j=0; j<4; j++) delete [] f2pCorrelation2[g][i][j];
  //for(Int_t g=0; g<4; g++) for(Int_t i=0; i<6; i++) for(Int_t j=0; j<8; j++) delete [] f3pCorrelation2[g][i][j];
  for(Int_t g=0; g<4; g++) delete [] fMult[g];
}



//_______________________________________________________________________________
void AliChargeTwoPwrtRP::Init(TAxis* ax, TAxis* ay, TAxis* az) {

  SetAxes(*ax,*ay,*az);

  Int_t nx=0,ny=0,nz=0;
  Int_t minTrackVar=0;

    for(Int_t idim=0;idim<3;++idim) {
      if(idim==0&&ax) {nx=fBinLimits[idim].GetNbins();fNdim=1;}
      if(idim==1&&ay) {ny=fBinLimits[idim].GetNbins();fNdim=2;}
      if(idim==2&&az) {nz=fBinLimits[idim].GetNbins();fNdim=3;}
    }

    if(GetVarX()>minTrackVar) {fTrackVarX=GetVarX();fTrackVarMap[0]=0;}
    else if(GetVarY()>minTrackVar) {fTrackVarX=GetVarY();nx=ny;ny=0;fTrackVarMap[0]=1;if(GetVarZ()>minTrackVar) {fTrackVarY=GetVarZ();ny=nz;nz=0;fTrackVarMap[1]=2;}}
    else if(GetVarZ()>minTrackVar) {fTrackVarX=GetVarZ();nx=nz;ny=0;nz=0;fTrackVarMap[0]=2;}
    if(GetVarX()>minTrackVar){
      if(GetVarY()>minTrackVar) {fTrackVarY=GetVarY();fTrackVarMap[1]=1;if(GetVarZ()>minTrackVar) {fTrackVarZ=GetVarZ();nz;fTrackVarMap[2]=2;}}
      else if(GetVarZ()>minTrackVar) {fTrackVarY=GetVarZ();ny=nz;nz=0;fTrackVarMap[1]=2;}
    }


    fNbinsX=nx;
    fNbinsY=ny;
    fNbinsZ=nz;

    if(fNbinsX==0) fNbinsX=1;
    if(fNbinsY==0) fNbinsY=1;
    if(fNbinsZ==0) fNbinsZ=1;



  //if(nz>0) {
  //  fNbins = nx*ny*nz;

  //}
  //if(ny>0&&nz==0){
  //  fNbins = nx*ny;
  //}
  //if(ny==0){
  //  fNbins = nx;
  //}



  if(nz>0) {
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) {f2pCorrelationShort[g][i] = new Float_t [nx*ny*nz]; for(Int_t k=0; k<(nx*ny*nz); k++) f2pCorrelationShort[g][i][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) {f3pCorrelationShort[g][i][j] = new Float_t [nx*ny*nz]; for(Int_t k=0; k<(nx*ny*nz); k++) f3pCorrelationShort[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) {f2pCorrelationShort2[g][i] = new Float_t [nx*ny*nz]; for(Int_t k=0; k<(nx*ny*nz); k++) f2pCorrelationShort2[g][i][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) {f3pCorrelationShort2[g][i][j] = new Float_t [nx*ny*nz]; for(Int_t k=0; k<(nx*ny*nz); k++) f3pCorrelationShort2[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<4; g++) fMult[g]= new Int_t [nx*ny*nz];
    for(Int_t g=0; g<4; g++) for(Int_t k=0; k<(nx*ny*nz); k++) fMult[g][k] = 0;
    fNbins = nx*ny*nz;

  }
  if(ny>0&&nz==0){
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) {f2pCorrelationShort[g][i] = new Float_t [nx*ny]; for(Int_t k=0; k<(nx*ny); k++) f2pCorrelationShort[g][i][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) {f3pCorrelationShort[g][i][j] = new Float_t [nx*ny]; for(Int_t k=0; k<(nx*ny); k++) f3pCorrelationShort[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) {f2pCorrelationShort2[g][i] = new Float_t [nx*ny]; for(Int_t k=0; k<(nx*ny); k++) f2pCorrelationShort2[g][i][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) {f3pCorrelationShort2[g][i][j] = new Float_t [nx*ny]; for(Int_t k=0; k<(nx*ny); k++) f3pCorrelationShort2[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<4; g++) fMult[g]= new Int_t [nx*ny];
    for(Int_t g=0; g<4; g++) for(Int_t k=0; k<(nx*ny); k++) fMult[g][k] = 0;
    fNbins = nx*ny;
  }
  if(ny==0){
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) {f2pCorrelationShort[g][i] = new Float_t [nx]; for(Int_t k=0; k<(nx); k++) f2pCorrelationShort[g][i][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) {f3pCorrelationShort[g][i][j] = new Float_t [nx]; for(Int_t k=0; k<(nx); k++) f3pCorrelationShort[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) {f2pCorrelationShort2[g][i] = new Float_t [nx]; for(Int_t k=0; k<(nx); k++) f2pCorrelationShort2[g][i][k] = 0.0;}
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) {f3pCorrelationShort2[g][i][j] = new Float_t [nx]; for(Int_t k=0; k<(nx); k++) f3pCorrelationShort2[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<4; g++) fMult[g]= new Int_t [nx];
    for(Int_t g=0; g<4; g++) for(Int_t k=0; k<nx; k++) fMult[g][k] = 0;
    fNbins = nx;
  }


  SetAxes(*ax,*ay,*az);
    TString title;
    TString charges[5] = {"aa","pp","mm","pm","mp"};
    for(Int_t g=0; g<4; g++) {
      title = Form("c_%s_%sx%s_11_%sx%s_%s_c2_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f2pCorrelationTHn[g][0][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f2pCorrelationTHn[g][0][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title,fNdim,fBinLimits);

      title = Form("c_%s_%sx%s_22_%sx%s_%s_c2_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f2pCorrelationTHn[g][1][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f2pCorrelationTHn[g][1][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      title = Form("c_%s_%sx%s_33_%sx%s_%s_c2_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f2pCorrelationTHn[g][2][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f2pCorrelationTHn[g][2][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      title = Form("c_%s_%sx%s_44_%sx%s_%s_c2_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f2pCorrelationTHn[g][3][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f2pCorrelationTHn[g][3][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      title = Form("s_%s_%sx%s_11_%sx%s_%s_c2_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f2pCorrelationTHn[g][4][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f2pCorrelationTHn[g][4][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title,fNdim,fBinLimits);

      title = Form("s_%s_%sx%s_22_%sx%s_%s_c2_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f2pCorrelationTHn[g][5][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f2pCorrelationTHn[g][5][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      title = Form("s_%s_%sx%s_33_%sx%s_%s_c2_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f2pCorrelationTHn[g][6][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f2pCorrelationTHn[g][6][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      title = Form("s_%s_%sx%s_44_%sx%s_%s_c2_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f2pCorrelationTHn[g][7][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f2pCorrelationTHn[g][7][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      //title = Form("c_%s_%sx%s_11_%sx%s_%s_c2xy_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      //title += fXaxisLabel;
      //if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      //f2pCorrelationTHn[g][2][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      //f2pCorrelationTHn[g][2][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      //title = Form("c_%s_%sx%s_11_%sx%s_%s_c2yx_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      //title += fXaxisLabel;
      //if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      //f2pCorrelationTHn[g][3][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      //f2pCorrelationTHn[g][3][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      //title = Form("c_%s_%sx%s_11_%sx%s_%s_c2xx_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      //title += fXaxisLabel;
      //if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      //f2pCorrelationTHn[g][4][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      //f2pCorrelationTHn[g][4][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

      //title = Form("c_%s_%sx%s_11_%sx%s_%s_c2yy_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
      //title += fXaxisLabel;
      //if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      //f2pCorrelationTHn[g][5][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      //f2pCorrelationTHn[g][5][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

    for(Int_t ep=0; ep<AliChargeOnePwrtRP::Nqvectors; ep++) {
      title = Form("c_%s_%sx%s_1m1_%sx%sxc_%sa_c3%s_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data(),fEventPlanes[ep].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f3pCorrelationTHn[g][0][ep][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f3pCorrelationTHn[g][0][ep][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);
    }

    for(Int_t ep=0; ep<AliChargeOnePwrtRP::Nqvectors; ep++) {
      title = Form("c_%s_%sx%s_2m2_%sx%sxc_%sa_c3%s_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data(),fEventPlanes[ep].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f3pCorrelationTHn[g][1][ep][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f3pCorrelationTHn[g][1][ep][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);
    }

    for(Int_t ep=0; ep<AliChargeOnePwrtRP::Nqvectors; ep++) {
      title = Form("c_%s_%sx%s_13_%sx%sxc_%sa_c3%s_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data(),fEventPlanes[ep].Data());
      title += fXaxisLabel;
      if(fYaxisLabel!="") title += "x"+fYaxisLabel;
      f3pCorrelationTHn[g][2][ep][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      f3pCorrelationTHn[g][2][ep][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);
    }


   // title = Form("c_%s_%sx%s_1m1_%sx%sxc_%sa_c3EPVZEROAsin_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
   // title += fXaxisLabel;
   // if(fYaxisLabel!="") title += "x"+fYaxisLabel;
   // f3pCorrelationTHn[g][1][0][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
   // f3pCorrelationTHn[g][1][0][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

   // title = Form("c_%s_%sx%s_1m1_%sx%sxc_%sa_c3EPVZEROC_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
   // title += fXaxisLabel;
   // if(fYaxisLabel!="") title += "x"+fYaxisLabel;
   // f3pCorrelationTHn[g][1][0][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
   // f3pCorrelationTHn[g][1][0][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

   // title = Form("c_%s_%sx%s_1m1_%sx%sxc_%sa_c3EPVZEROCsin_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
   // title += fXaxisLabel;
   // if(fYaxisLabel!="") title += "x"+fYaxisLabel;
   // f3pCorrelationTHn[g][3][0][0]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
   // f3pCorrelationTHn[g][3][0][1]=AliCMEVarManager::CreateTHnF(title+"_squared",title+"_squared",fNdim,fBinLimits);

    title = Form("c_%s_%sx%s_%sx%s_%s_mult_",fEventSelection.Data(),fTracking1Quality.Data(),fTracking2Quality.Data(),fPID1.Data(),fPID2.Data(),charges[g+1].Data());
    title += fXaxisLabel;
    if(fYaxisLabel!="") title += "x"+fYaxisLabel;
    fCorrelationTHnMult[g]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);}
  }



//_______________________________________________________________________________
  void AliChargeTwoPwrtRP::SetTHn(TString eventSelection, TString tracking1Quality, TString tracking2Quality,TString p1, TString p2, TString d1, TString d2, Char_t* eventplanes[]) {

    fEventSelection = eventSelection;
    fTracking1Quality = tracking1Quality;
    fTracking2Quality = tracking2Quality;
    fPID1 = p1;
    fPID2 = p2;
    fXaxisLabel = d1;
    fYaxisLabel = d2;
    for(Int_t ep=0; ep<AliChargeOnePwrtRP::Nqvectors; ep++) {
      fEventPlanes[ep] = eventplanes[ep];
    }
  }

//_________________________________________________________________
  void AliChargeTwoPwrtRP::Clear(){
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) for(Int_t k=0; k<fNbins; k++) f2pCorrelationShort[g][i][k] = 0.0;
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) for(Int_t k=0; k<fNbins; k++) f3pCorrelationShort[g][i][j][k] = 0.0;
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N2p; i++) for(Int_t k=0; k<fNbins; k++) f2pCorrelationShort2[g][i][k] = 0.0;
    for(Int_t g=0; g<4; g++) for(Int_t i=0; i<N3p; i++) for(Int_t j=0; j<AliChargeOnePwrtRP::Nqvectors; j++) for(Int_t k=0; k<fNbins; k++) f3pCorrelationShort2[g][i][j][k] = 0.0;
    for(Int_t g=0; g<4; g++) for(Int_t k=0; k<fNbins; k++) fMult[g][k] = 0;
}
