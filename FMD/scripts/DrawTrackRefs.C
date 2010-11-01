//____________________________________________________________________
//
// $Id: DrawHits.C 30718 2009-01-22 16:07:40Z cholm $
//
// Script that contains a class to draw hits, using the
// AliFMDInputHits class in the util library. 
//
// It draws the energy loss versus the p/(mq^2).  It can be overlayed
// with the Bethe-Bloc curve to show how the simulation behaves
// relative to the expected. 
//
// Use the script `Compile.C' to compile this class using ACLic. 
//
#include <TH2D.h>
#include <AliFMDHit.h>
#include <AliFMDInput.h>
#include <iostream>
#include <TStyle.h>
#include <TArrayF.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TLegend.h>
#include <TArrow.h>
#include <TLatex.h>
#include <TF1.h>
#include <AliStack.h>
#include <AliTrackReference.h>
#include <AliFMDStripIndex.h>
#include <AliFMDGeometry.h>
#include <AliGenEventHeader.h>
#include <AliHeader.h>
#include <THStack.h>

/** @class DrawHits
    @brief Draw hit energy loss
    @code 
    Root> .L Compile.C
    Root> Compile("DrawHits.C")
    Root> DrawHits c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class DrawTrackRefs : public AliFMDInput
{
private:
  struct Hists 
  {
    TString fPrefix;
    TString fWhere;
    TH2* fAll;
    TH2* fPrimary;
    TH2* fSecondary;
    THStack* fStack;

    //__________________________________________________________________
    Hists(const char* prefix, const char* where)
      : fPrefix(prefix), fWhere(where), 
	fAll(0), fPrimary(0), fSecondary(0), fStack(0)
    {
      fAll       = MakeHist(Form("%sAll", prefix), 
			    Form("All tracks in %s", where));
      fPrimary   = MakeHist(Form("%sPrimary", prefix), 
			    Form("All primaries in %s", where));
      fSecondary = MakeHist(Form("%sSecondary", prefix), 
			    Form("All secondaries in %s", where));
    }
     //__________________________________________________________________
    TH2* MakeHist(const char* name, const char* title) 
    {
      TH2D* ret = new TH2D(name, title, 200, -4, 6, 40, 0, 2*TMath::Pi());
      ret->SetDrawOption("colz");
      return ret;
    }
    //__________________________________________________________________
    void Draw(TVirtualPad* pad, Bool_t stack=false, Option_t* option="colz") 
    {
      if (stack) { 
	pad->cd();
	pad->Divide(1,2, 0, 0);
	pad->cd(1);
	fStack = new THStack(Form("%sDisplay", fPrefix.Data()), 
			     Form("2nd over 1st particles in %s", 
				  fWhere.Data()));
	TH1D* prim = 0;
	TH1D* sec = 0;
	fStack->Add((prim = fPrimary->ProjectionX()));
	fStack->Add((sec  = fSecondary->ProjectionX()));
	prim->SetFillColor(kGray);
	sec->SetFillColor(kRed);
	sec->SetFillStyle(3001);
	fStack->Draw();

	pad->cd(2);
	TH1* ratio = new TH1D(*sec);
	ratio->SetName(Form("%sRatio", fPrefix.Data()));
	ratio->SetTitle(Form("2nd over 1st particles in %s", fWhere.Data()));
	ratio->Divide(prim);
	ratio->SetFillColor(kRed);
	ratio->SetFillStyle(3001);
	ratio->Draw();

	return;
      }
      pad->Divide(1,3,0,0);
      TVirtualPad* p1 = pad->cd(1);
      p1->SetLogz();
      fAll->Draw(option);
      p1 = pad->cd(2);
      p1->SetLogz();
      fPrimary->Draw(option);
      p1 = pad->cd(3);
      p1->SetLogz();
      fSecondary->Draw(option);
      pad->Modified();
      pad->Update();
      pad->cd();
    }
  };
  Hists fTotal;
  Hists fFMD1I;
  Hists fFMD2I;
  Hists fFMD2O;
  Hists fFMD3I;
  Hists fFMD3O;
public:
  //__________________________________________________________________
  DrawTrackRefs() 
    : AliFMDInput("galice.root"), 
      fTotal("total", "ALICE"), 
      fFMD1I("fmd1I", "FMD1i"),
      fFMD2I("fmd2I", "FMD2i"),
      fFMD2O("fmd2O", "FMD2o"),
      fFMD3I("fmd3I", "FMD3i"),
      fFMD3O("fmd3O", "FMD3o")
  { 
    AddLoad(kKinematics);
    AddLoad(kTrackRefs);
    AddLoad(kGeometry);
    AddLoad(kHeader);
  }
  //__________________________________________________________________
  Bool_t ProcessTrackRef(AliTrackReference* trackRef, TParticle* p)
  {
    if (trackRef->DetectorId() != AliTrackReference::kFMD) return kTRUE;

    if (!p->GetPDG() || p->GetPDG()->Charge() == 0) return kTRUE;

    // Double_t eta = p->Eta();
    // Double_t phi = p->Phi();

    UShort_t d   = 0;
    Char_t   r   = '\0';
    UShort_t s   = 0;
    UShort_t t   = 0;
    AliFMDStripIndex::Unpack(trackRef->UserId(),d,r,s,t);
    
    Double_t x, y, z;
    AliFMDGeometry::Instance()->Detector2XYZ(d,r,s,t,x,y,z);

    AliGenEventHeader* genHeader = fHeader->GenEventHeader();
    TArrayF v;
    genHeader->PrimaryVertex(v);
    z -= v.fArray[2];
    
    Double_t eta, phi, theta, rr;
    AliFMDGeometry::XYZ2REtaPhiTheta(x, y, z, rr, eta, phi, theta);
    if (phi < 0)               phi += 2*TMath::Pi();
    if (phi > 2 * TMath::Pi()) phi -= 2*TMath::Pi();

    Hists* hists = 0;
    switch (d) { 
    case 1: hists = &fFMD1I; break;
    case 2: hists = (r == 'I' || r == 'i') ? &fFMD2I : &fFMD2O; break;
    case 3: hists = (r == 'I' || r == 'i') ? &fFMD3I : &fFMD3O; break;
    }
    if (!hists) return kTRUE;

    hists->fAll->Fill(eta, phi);
    
    if (fStack->IsPhysicalPrimary(trackRef->GetTrack())) 
      hists->fPrimary->Fill(eta,phi); 
    else
      hists->fSecondary->Fill(eta, phi);

    return kTRUE;
  }
  //__________________________________________________________________
  Bool_t ProcessParticle(Int_t id, TParticle* p)
  {
    if (!p) return kTRUE;

    if (!p->GetPDG() || p->GetPDG()->Charge() == 0) return kTRUE;
    // std::cout << p->GetPDG()->GetName() << std::endl;

    Double_t eta = p->Eta();
    Double_t phi = p->Phi();

    fTotal.fAll->Fill(eta, phi);
    
    if (fStack->IsPhysicalPrimary(id)) fTotal.fPrimary->Fill(eta,phi); 
    else                               fTotal.fSecondary->Fill(eta, phi);

    return kTRUE;
  }

  //__________________________________________________________________
  Bool_t Finish()
  {
    gStyle->SetPalette(1);
    // gStyle->SetOptTitle(0);
    gStyle->SetTitleBorderSize(1);
    gStyle->SetTitleFillColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderSize(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderSize(0);

    TCanvas* c2d = new TCanvas("c2d", "2D plots", 1200, 700);
    c2d->Divide(6);
    fTotal.Draw(c2d->cd(1));   c2d->cd();
    fFMD1I.Draw(c2d->cd(6));   c2d->cd();
    fFMD2I.Draw(c2d->cd(5));   c2d->cd();
    fFMD2O.Draw(c2d->cd(4));   c2d->cd();
    fFMD3I.Draw(c2d->cd(2));   c2d->cd();
    fFMD3O.Draw(c2d->cd(3));   c2d->cd();

    c2d->Modified();
    c2d->Update();
    c2d->cd();
    c2d->SaveAs("kine_etaphi.png");

    TCanvas* c1d = new TCanvas("c1d", "1D plots", 1200, 700);
    c1d->Divide(6);
    fTotal.Draw(c1d->cd(1), kTRUE);   c1d->cd();
    fFMD1I.Draw(c1d->cd(6), kTRUE);   c1d->cd();
    fFMD2I.Draw(c1d->cd(5), kTRUE);   c1d->cd();
    fFMD2O.Draw(c1d->cd(4), kTRUE);   c1d->cd();
    fFMD3I.Draw(c1d->cd(2), kTRUE);   c1d->cd();
    fFMD3O.Draw(c1d->cd(3), kTRUE);   c1d->cd();

    c1d->Modified();
    c1d->Update();
    c1d->cd();
    c1d->SaveAs("kine_eta.png");


    return kTRUE;
  }
  
  ClassDef(DrawTrackRefs,0);
};

//____________________________________________________________________
//
// EOF
//
