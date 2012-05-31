//____________________________________________________________________
//
// $Id$
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
#include "Poisson.C"
#include <TMath.h>
#include <TCanvas.h>
#include <AliFMDHit.h>
#include <AliFMDGeometry.h>

/** @class PoissonHit
    @brief Make a poisson reconstruction and compare to simulated hits
    @code 
    Root> .L Compile.C
    Root> Compile("Poisson.C")
    Root> Compile("PoissonHit.C")
    Root> PoissonHit c
    Root> c.Run();
    @endcode
    @ingroup FMD_script
 */
class PoissonHit : public Poisson
{
protected:
  TH2D*  fHits; // Histogram 
  TH2D*  fDiff; // Histogram 
public:
  /** Constructor 
      @param threshold Threshold
      @param nEta      # of @f$ \eta@f$ bins
      @param minEta    minimum @f$ \eta@f$
      @param maxEta    maximum @f$ \eta@f$
      @param nPhi      # of @f$ \eta@f$ bins
      @param minPhi    minimum @f$ \varphi@f$  
      @param maxPhi    maximum @f$ \varphi@f$ */
  PoissonHit(Double_t threshold=.3,
	     Int_t nEta=120, Float_t minEta=-6, Float_t maxEta=6, 
	     Int_t nPhi=4,   Float_t minPhi=0,  Float_t maxPhi=2*TMath::Pi())
    : Poisson(threshold, nEta, minEta, maxEta, nPhi, minPhi, maxPhi)
  { 
    AddLoad(kHits);
    AddLoad(kGeometry);
    fHits  = new TH2D(*fEmpty);
    fHits->SetName("hits");
    fHits->SetTitle("# of hits");
    fDiff  = new TH2D(*fEmpty);
    fDiff->SetName("diff");
    fDiff->SetTitle("Difference between poisson and hits");
    fHits->SetXTitle("#eta");  
    fHits->SetYTitle("#phi"); 
    fHits->SetZTitle("N");
    fDiff->SetXTitle("#eta");  
    fDiff->SetYTitle("#phi"); 
    fDiff->SetZTitle("#frac{N_{hit}-N_{poisson}}{N_{hit}}");
  }
  /** Initialize the analyser. Opens the output file. 
      @return @c true on success. */
  virtual Bool_t Init() 
  {
    if (!Poisson::Init()) return kFALSE;
    AliFMDGeometry::Instance()->Init();
    AliFMDGeometry::Instance()->InitTransformations();    
    return kTRUE;
  }
  /** Get the @f$ \eta@f$ and @f$\varphi@f$ corresponding to the
      spatial coordinates @f$ \mathbf{v} = (x,y,z)@f$
      @param x    X coordinate
      @param y    Y coordinate
      @param z    Z coordinate
      @param eta  Psuedo rapidity @f$ \eta@f$
      @param phi  Azimuthal angle @f$\varphi@f$
  */
  void PhysicalCoordinates(Double_t x, Double_t y, Double_t z, 
			   Double_t& eta, Double_t& phi)
  {
    Double_t r, theta;
    phi   =  TMath::ATan2(y, x);
    r     =  TMath::Sqrt(y * y + x * x);
    theta =  TMath::ATan2(r, z);
    eta   = -TMath::Log(TMath::Tan(theta / 2));
    if (phi < 0) phi += 2 * TMath::Pi();
  }
  /** Process one hit.  Increment bin corresponding to strip. 
      @param hit Hit. 
      @return @c true on success */
  virtual Bool_t ProcessHit(AliFMDHit* hit, TParticle*)
  {
    Double_t x, y, z;
#if 0
    AliFMDGeometry* geom = AliFMDGeometry::Instance();
    geom->Detector2XYZ(hit->Detector(),hit->Ring(),hit->Sector(),
		       hit->Strip(),x,y,z);
#else
    x = hit->X();
    y = hit->Y();
    z = hit->Z();
#endif
    Double_t eta, phi;
    PhysicalCoordinates(x, y, z, eta, phi);
    fHits->Fill(eta, phi);
    return kTRUE;
  }
  /** Begining of event
      @param event Event number
      @return @c false on error */
  virtual Bool_t Begin(Int_t event) 
  {
    if (!Poisson::Begin(event)) return kFALSE;
    fHits->Clear();
    fDiff->Clear();
    return kTRUE;
  }
  /** Let the poisson code do it's job, and then compare to the
      numberr of hits. 
      @return @c true  */
  virtual Bool_t End() 
  {
    if (!Poisson::End()) return kFALSE;
    fDiff->Add(fMult,fHits,-1.,1.);
    fDiff->Divide(fHits);
    if (!gROOT->IsBatch()) { 
      gStyle->SetPalette(1);
      TCanvas* c1 = new TCanvas("hits", "Hit multiplicity");
      c1->SetFillColor(0);
      fHits->Draw("colz");
      TCanvas* c2 = new TCanvas("diff", "Difference between Hit and poisson");
      c2->SetFillColor(0);
      fDiff->Draw("colz");
      TCanvas* c3 = new TCanvas("empty", "# of Empty strips");
      c3->SetFillColor(0);
      fEmpty->Draw("colz");
      TCanvas* c4 = new TCanvas("total", "Total # of strips");
      c4->SetFillColor(0);
      fTotal->Draw("colz");
    }
    return kTRUE;
  }
  
  ClassDef(PoissonHit,0);
};

//____________________________________________________________________
//
// EOF
//
