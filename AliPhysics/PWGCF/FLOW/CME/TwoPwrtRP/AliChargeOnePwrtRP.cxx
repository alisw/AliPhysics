
#include "AliChargeOnePwrtRP.h"
#include "AliCMEVarManager.h"

ClassImp(AliChargeOnePwrtRP)


//____________________________________________________________________________
AliChargeOnePwrtRP::AliChargeOnePwrtRP() :
  TNamed(),
  fNbinsX(0),
  fNbinsY(0),
  fNbinsZ(0),
  fNbins(0),
  fVarX(-1),
  fVarY(-1),
  fVarZ(-1),
  fTrackVarX(0),
  fTrackVarY(0),
  fTrackVarZ(0),
  fTrackVarMap(),
  fNdim(0),
  fBin(0),
  fBins(),
  fBinLimits(),
  fNbin(),
  fEventSelection(0x0),
  fTracking1Quality(0x0),
  fPID1(""),
  fEventPlanes(),
  fXaxisLabel(""),
  fYaxisLabel(""),
  fMinHarmonic(1000),
  fMaxHarmonic(0),
  fMultiplicity(),
  fQvector(),
  fTrackCutsOnePwrtRP(0x0),
  fSumQ(),
  fSumQ2(),
  fSumQMult(),
  fCorrelationTHn(),
  fCorrelationTHn2(),
  fCorrelationTHnMult()
{

  for(Int_t i=0; i<3; i++) fTrackVarMap[i]=-1;
  for(Int_t i=0; i<7; i++) fBins[i]=-1;
  //for(Int_t i=0; i<3; i++) for(Int_t j=0; j<Nharmonics; j++) for(Int_t k=0; k<2; k++) fQvector[i][j][k]=0x0;
  for(Int_t i=0; i<Nqvectors; i++) for(Int_t j=0; j<3; j++) for(Int_t k=0; k<Nharmonics; k++) for(Int_t l=0; l<2; l++)  fCorrelationTHn[i][j][k][l]=0x0;
  for(Int_t i=0; i<Nqvectors; i++) for(Int_t j=0; j<3; j++) for(Int_t k=0; k<Nharmonics; k++) for(Int_t l=0; l<2; l++)  fCorrelationTHn2[i][j][k][l]=0x0;
  for(Int_t i=0; i<Nqvectors; i++) for(Int_t j=0; j<3; j++) fCorrelationTHnMult[i][j]=0x0;
  for(Int_t i=0; i<Nqvectors; i++) for(Int_t j=0; j<3; j++) for(Int_t k=0; k<Nharmonics; k++) for(Int_t l=0; l<2; l++)  fSumQ[i][j][k][l]=0x0;
  for(Int_t i=0; i<Nqvectors; i++) for(Int_t j=0; j<3; j++) for(Int_t k=0; k<Nharmonics; k++) for(Int_t l=0; l<2; l++)  fSumQ2[i][j][k][l]=0x0;
  for(Int_t i=0; i<Nqvectors; i++) for(Int_t j=0; j<3; j++) fSumQMult[i][j]=0x0;

}

//_________________________________________________________________
void AliChargeOnePwrtRP::SetTHn(TString eventSelection, TString tracking1Quality, TString p1,  TString d1, TString d2/*=""*/, Char_t* eventplanes[]) {

  fEventSelection = eventSelection;
  fTracking1Quality = tracking1Quality;
  fPID1 = p1;
  fXaxisLabel = d1;
  fYaxisLabel = d2;
  for(Int_t ep=0; ep<Nqvectors; ep++) {
    fEventPlanes[ep] = eventplanes[ep];
  }


}

//_________________________________________________________________
void AliChargeOnePwrtRP::Init(TAxis* ax, TAxis* ay, TAxis* az) {

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


  if(nz>0) {
    for(Int_t g=0; g<3; g++) for(Int_t i=0; i<Nharmonics; i++) for(Int_t j=0; j<2; j++) {fQvector[g][i][j] = new Float_t [nx*ny*nz];  for(Int_t k=0; k<(nx*ny*nz); k++) fQvector[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<3; g++) fMultiplicity[g] = new Int_t [nx*ny*nz];
    for(Int_t g=0; g<3; g++) for(Int_t k=0; k<(nx*ny*nz); k++) fMultiplicity[g][k] = 0;
    fNbins = nx*ny*nz;

  }
  if(ny>0&&nz==0){
    for(Int_t g=0; g<3; g++) for(Int_t i=0; i<Nharmonics; i++) for(Int_t j=0; j<2; j++) {fQvector[g][i][j] = new Float_t [nx*ny];  for(Int_t k=0; k<(nx*ny); k++) fQvector[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<3; g++) fMultiplicity[g] = new Int_t [nx*ny];
    for(Int_t g=0; g<3; g++) for(Int_t k=0; k<(nx*ny); k++) fMultiplicity[g][k] = 0;
    fNbins = nx*ny;
  }
  if(ny==0){
    for(Int_t g=0; g<3; g++) for(Int_t i=0; i<Nharmonics; i++) for(Int_t j=0; j<2; j++) {fQvector[g][i][j] = new Float_t [nx];  for(Int_t k=0; k<nx; k++) fQvector[g][i][j][k] = 0.0;}
    for(Int_t g=0; g<3; g++) fMultiplicity[g] = new Int_t [nx];
    for(Int_t g=0; g<3; g++) for(Int_t k=0; k<nx; k++) fMultiplicity[g][k] = 0;
    fNbins = nx;
  }


     fNbin[0] = fNbinsX;
     fNbin[1] = fNbinsY;
     fNbin[2] = fNbinsZ;

  for(Int_t i=0; i<7; i++) fBins[i]=0;




    TString title,titlex,titley;
    TString titleq,titleqx,titleqy;
    TString charges[3] = {"pa","ma","aa"};
    for(Int_t ep=0; ep<Nqvectors; ep++) {
    for(Int_t g=0; g<3; g++) {
    for(Int_t ih=fMinHarmonic; ih<=fMaxHarmonic; ih++) {
      //cout<<eventSelection<<endl;
      //cout<<tracking1Quality<<endl;
      //cout<<p1<<endl;
      //cout<<charges[g]<<endl;
       title = Form("c_%s_%s_%sxc_%s_c2EP%s_mult_",fEventSelection.Data(),fTracking1Quality.Data(),fPID1.Data(),charges[g].Data(),fEventPlanes[ep].Data());
      titlex = Form("c_%s_%s_%d%d_%sxc_%s_c2EP%sxx_"   ,fEventSelection.Data(),fTracking1Quality.Data(),ih,ih,fPID1.Data(),charges[g].Data(),fEventPlanes[ep].Data());
      titley = Form("c_%s_%s_%d%d_%sxc_%s_c2EP%syy_"   ,fEventSelection.Data(),fTracking1Quality.Data(),ih,ih,fPID1.Data(),charges[g].Data(),fEventPlanes[ep].Data());
      titleq  = Form("SumQ_%s_%s_%s_%s_%s_mult_"  ,fEventSelection.Data(),fTracking1Quality.Data(),fPID1.Data(),charges[g].Data(),fEventPlanes[ep].Data());
      titleqx = Form("SumQ_%s_%s_%d_%s_%s_%sx_"   ,fEventSelection.Data(),fTracking1Quality.Data(),ih,fPID1.Data(),charges[g].Data(),fEventPlanes[ep].Data());
      titleqy = Form("SumQ_%s_%s_%d_%s_%s_%sy_"   ,fEventSelection.Data(),fTracking1Quality.Data(),ih,fPID1.Data(),charges[g].Data(),fEventPlanes[ep].Data());
      //cout<<titlex<<"  "<<titley<<endl;
      titlex += fXaxisLabel; titley += fXaxisLabel; title+=fXaxisLabel;
      if(fYaxisLabel!="") {titlex += "x"+fYaxisLabel; titley += "x"+fYaxisLabel;title+="x"+fYaxisLabel;}
      titleqx += fXaxisLabel; titleqy += fXaxisLabel; titleq+=fXaxisLabel;
      if(fYaxisLabel!="") {titleqx += "x"+fYaxisLabel; titleqy += "x"+fYaxisLabel;titleq+="x"+fYaxisLabel;}
      fCorrelationTHn[ep][g][ih-1][0] =AliCMEVarManager::CreateTHnF(titlex,titlex,fNdim,fBinLimits);
      fCorrelationTHn[ep][g][ih-1][1] =AliCMEVarManager::CreateTHnF(titley,titley,fNdim,fBinLimits);
      fCorrelationTHn2[ep][g][ih-1][0]=AliCMEVarManager::CreateTHnF(titlex+"_squared",titlex+"_squared",fNdim,fBinLimits);
      fCorrelationTHn2[ep][g][ih-1][1]=AliCMEVarManager::CreateTHnF(titley+"_squared",titley+"_squared",fNdim,fBinLimits);
      fSumQ[ep][g][ih-1][0] =AliCMEVarManager::CreateTHnF(titleqx, titleqx,fNdim,fBinLimits);
      fSumQ[ep][g][ih-1][1] =AliCMEVarManager::CreateTHnF(titleqy, titleqy,fNdim,fBinLimits);
      fSumQ2[ep][g][ih-1][0] =AliCMEVarManager::CreateTHnF(titleqx+"_squared", titleqx+"_squared",fNdim,fBinLimits);
      fSumQ2[ep][g][ih-1][1] =AliCMEVarManager::CreateTHnF(titleqy+"_squared", titleqy+"_squared",fNdim,fBinLimits);
    }
      fCorrelationTHnMult[ep][g]=AliCMEVarManager::CreateTHnF(title,title,fNdim,fBinLimits);
      fSumQMult[ep][g] =AliCMEVarManager::CreateTHnF(titleq,titleq,fNdim,fBinLimits);
    }}

   //  for(Int_t id=0; id<nDim; id++){
   //    for(Int_t ix=0; ix<(fBinLimits[id].GetSize()-1); ix++){
   //      fBinLimits[id].SetAt((fBinLimits[id].At(ix)+fBinLimits[id].At(ix+1))/2., ix);
   //  }}


}



//_______________________________________________________________________________
AliChargeOnePwrtRP::~AliChargeOnePwrtRP()
{
  //
  // De-Constructor
  //
  for(Int_t g=0; g<3; g++) for(Int_t i=0; i<Nharmonics; i++) for(Int_t j=0; j<2; j++) delete [] fQvector[g][i][j];
  for(Int_t g=0; g<3; g++) delete [] fMultiplicity[g];

  for(Int_t ep=0; ep<Nqvectors; ep++) for(Int_t g=0; g<3; g++){
    for(Int_t ih=fMinHarmonic; ih<=fMaxHarmonic; ih++) for(Int_t ico=0; ico<2; ico++){
      delete fCorrelationTHn[ep][g][ih-1][ico];
      delete fCorrelationTHn2[ep][g][ih-1][ico];
      delete fSumQ[ep][g][ih-1][ico];
      delete fSumQ2[ep][g][ih-1][ico];
    }
    delete fCorrelationTHnMult[ep][g];
    delete fSumQMult[ep][g];
  }
}
