// MyAnalysis.C 
#ifndef __CINT__
# include <AliAnalysisManager.h>
# include <AliESDEvent.h>
# include <AliMultiplicity.h>
# include <AliESDVertex.h>
# include <TH1D.h>
# include <TH2D.h>
#else 
class TH1D;
class TH2D;
#endif 
#include <AliAnalysisTaskSE.h>
class MyAnalysis : public AliAnalysisTaskSE
{
public:
  MyAnalysis() 
    : AliAnalysisTaskSE(), fList(0), fMult(0), fVz(0)
  {}
  MyAnalysis(const char* name) 
  : AliAnalysisTaskSE(name), fList(0), fMult(0), fVz(0)
  {
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class()); // For output from Terminate
    fBranchNames = "AliMultiplicity.,SPDVertex.,PrimaryVertex.";
  }
  MyAnalysis(const MyAnalysis& o) 
  : AliAnalysisTaskSE(o), fList(o.fList), fMult(o.fMult), fVz(o.fVz)
  {}
  virtual ~MyAnalysis() {}
  MyAnalysis& operator=(const MyAnalysis&) { return *this; }
  virtual void UserCreateOutputObjects()
  {
    fList = new TList();
    fList->SetName("Sums");
    fList->SetOwner();

    fMult = new TH2D("mult", "SPD tracklets", 80, -2, 2, 10, -10, 10);
    fMult->SetXTitle("#eta");
    fMult->SetYTitle("v_{z} [cm]");
    fMult->Sumw2();
    fMult->SetDirectory(0); // Disassociate from file 
    
    fVz = new TH1D("vz", "Interaction point", 10, -10, 10);
    fVz->SetXTitle("v_{z} [cm]");
    fVz->Sumw2();
    fVz->SetDirectory(0); // Disassociate from file 

    fList->Add(fMult);
    fList->Add(fVz);

    PostData(1, fList);
  }
  virtual void UserExec(Option_t* )
  {
    AliESDEvent* event = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!event) return;
    if (event->IsPileupFromSPD(3,0.8))   return;

    const AliESDVertex* vtx = event->GetPrimaryVertexSPD();
    if (!vtx || !vtx->GetStatus()) return;
    if (vtx->IsFromVertexerZ() && 
        (vtx->GetDispersion() > 0.2 ||  vtx->GetZRes() > 1.25 * 0.2))
      return;

    const AliMultiplicity* mult = event->GetMultiplicity();
    if (!mult) return;
    
    Double_t vz = vtx->GetZ();
    fVz->Fill(vz);

    Int_t nTracklets = mult->GetNumberOfTracklets();
    for (Int_t i = 0; i < nTracklets; i++) 
      fMult->Fill(mult->GetEta(i), vz);

    PostData(1, fList);
  }
  void  Terminate(Option_t *)
  {
    TList* l = dynamic_cast<TList*>(GetOutputData(1));
    if (!l) {
      Warning("Terminate", "No out data # 1 found");
      return;
    }
 
    TH2D* mult = static_cast<TH2D*>(l->FindObject("mult"));
    TH1D* vz   = static_cast<TH1D*>(l->FindObject("vz"));
    if (!mult || !vz) {
      Warning("Terminate", "Either 'mult' (%p) or 'vz' (%p) or both not found",
              mult, vz);
      return;
    }

    TList* output = new TList;   // Needed for new output from Terminate
    output->SetName("Results");  // 1st output re-opened read-only
    output->SetOwner();

    TH2D* out = static_cast<TH2D*>(mult->Clone("dndeta"));
    out->SetTitle("dN_{ch}/d#eta from SPD tracklets per vertex bin");
    out->SetZTitle("#frac{1}{N}#frac{dN_{ch}}{d#eta}");
    out->SetDirectory(0); // Disassociate from file 
    Int_t    nVz  = mult->GetNbinsY();
    Int_t    nEta = mult->GetNbinsX();
    for (Int_t iVz = 1; iVz <= nVz; iVz++) { 
      Double_t nEv = vz->GetBinContent(iVz);
      Double_t e1  = vz->GetBinError(iVz);
      Double_t sca = (nEv == 0 ? 0 : 1. / nEv);
      for (Int_t iEta = 1; iEta <= nEta; iEta++) { 
        Double_t c  = mult->GetBinContent(iEta,iVz);
        Double_t e  = mult->GetBinError(iEta,iVz);
        Double_t ee = TMath::Sqrt(c*c * e1*e1 + nEv*nEv * e*e) * sca*sca;
        out->SetBinContent(iEta, iVz, sca * c);
        out->SetBinError(iEta, iVz, ee);
      }
    }
    Double_t etaMin = mult->GetXaxis()->GetXmin();
    Double_t etaMax = mult->GetXaxis()->GetXmax();
    out->Scale(Double_t(nEta) / (etaMax-etaMin));

    output->Add(out);
    PostData(2, output);
  }
protected:
  TList*  fList;
  TH2D*   fMult;
  TH1D*   fVz;
  ClassDef(MyAnalysis, 1);  
};
//
// EOF
//

