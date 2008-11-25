#include "AliTRDclusterResolution.h"
#include "AliTRDtrackInfo/AliTRDclusterInfo.h"

#include "AliLog.h"

#include "TObjArray.h"
#include "TAxis.h"
#include "TH2I.h"
#include "TMath.h"

ClassImp(AliTRDclusterResolution)


//_______________________________________________________
AliTRDclusterResolution::AliTRDclusterResolution()
  : AliTRDrecoTask("CalibClRes", "Calibrate Cluster Resolution for Tracking")
  ,fInfo(0x0)
{
}

//_______________________________________________________
AliTRDclusterResolution::~AliTRDclusterResolution()
{

}

//_______________________________________________________
void AliTRDclusterResolution::ConnectInputData(Option_t *)
{
  fInfo = dynamic_cast<TObjArray *>(GetInputData(0));
}

//_______________________________________________________
void AliTRDclusterResolution::CreateOutputObjects()
{
  OpenFile(0, "RECREATE");
  fContainer = Histos();
}

//_______________________________________________________
void AliTRDclusterResolution::GetRefFigure(Int_t /*ifig*/)
{

}

//_______________________________________________________
TObjArray* AliTRDclusterResolution::Histos()
{
  if(fContainer) return fContainer;
  fContainer = new TObjArray(kN+1);
  //fContainer->SetOwner(kTRUE);

  TAxis at(kNTB, -0.075, (kNTB-.5)*.15);
  TAxis ad(kND, 0., .25);
  Int_t ih = 0;
  for(Int_t id=1; id<=ad.GetNbins(); id++){
    for(Int_t it=1; it<=at.GetNbins(); it++){
      fContainer->AddAt(new TH2I(Form("h_d%02dt%02d", id, it), Form("d_{wire}(%5.2f-%5.2f)[mm] x_{drift}(%5.2f-%5.2f)[mm]", 10.*ad.GetBinCenter(id)- 5.*ad.GetBinWidth(id), 10.*ad.GetBinCenter(id)+ 5.*ad.GetBinWidth(id), 10.*at.GetBinCenter(it)- 5.*at.GetBinWidth(it), 10.*at.GetBinCenter(it)+ 5.*at.GetBinWidth(it)), 30, -.15, .15, 100, -.5, .5), ih++);
    }
  }
  fContainer->AddAt(new TH2I("h_q", "", 50, 2.2, 7.5, 100, -.5, .5), ih++);
  return fContainer;
}

//_______________________________________________________
void AliTRDclusterResolution::Exec(Option_t *)
{
  Int_t det;
  Float_t x, y, z, q, dy, dydx, dzdx, cov[3];
  TAxis at(kNTB, -0.075, (kNTB-.5)*.15); Int_t it = 0;
  TAxis ad(kND, 0., .25); Int_t id = 0;
  TH2I *h2 = 0x0;
  const AliTRDclusterInfo *cli = 0x0;
  TIterator *iter=fInfo->MakeIterator();
  while((cli=dynamic_cast<AliTRDclusterInfo*>((*iter)()))){
    dy = cli->GetResolution();
    it = at.FindBin(cli->GetDriftLength());
    if(it==0 || it == at.GetNbins()+1){
      AliWarning(Form("Drift length %f outside allowed range", cli->GetDriftLength()));
      continue;
    }
    id = ad.FindBin(cli->GetAnisochronity());
    if(id==0 || id == ad.GetNbins()+1){
      AliWarning(Form("Distance to anode %f outside allowed range", cli->GetAnisochronity()));
      continue;
    }
    if(!(h2 = (TH2I*)fContainer->At((id-1)*kNTB+it-1))){
      AliWarning(Form("Missing histo at index idx[%3d] [id[%2d] it[%2d]] xd[%f] d[%f]\n", (id-1)*kNTB+it-1, id, it, cli->GetDriftLength(), cli->GetAnisochronity()));
      continue;
    }

    cli->GetGlobalPosition(y, z, dydx, dzdx, &cov[0]);
    h2->Fill(dydx, dy);

    // resolution as a function of cluster charge
    // only for phi equal exB 
    if(TMath::Abs(dydx)<.01){
      cli->GetCluster(det, x, y, z, q);
      h2 = (TH2I*)fContainer->At(kN);
      h2->Fill(TMath::Log(q), dy);
    }
  }
  PostData(0, fContainer);
}


//_______________________________________________________
Bool_t AliTRDclusterResolution::PostProcess()
{
  return kFALSE;
}


