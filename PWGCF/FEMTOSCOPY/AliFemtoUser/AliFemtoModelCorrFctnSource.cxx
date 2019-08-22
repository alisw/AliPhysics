////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// alifemtomodelcorrfctnsource - the class for correlation function which   ///
/// uses the model framework and weight generation and saves the generated   ///
/// emission source                                                          ///
/// Authors: Adam Kisiel, kisiel@mps.ohio-state.edu                          ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////
#ifdef __ROOT__
  ClassImp(AliFemtoModelCorrFctnSource)
#endif

#include "AliFemtoModelGausLCMSFreezeOutGenerator.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoModelCorrFctnSource.h"
#include "AliFemtoKTPairCut.h"
#include "AliFemtoAnalysisReactionPlane.h"

//_______________________
AliFemtoModelCorrFctnSource::AliFemtoModelCorrFctnSource():
  AliFemtoModelCorrFctnSource("CFSource", 50, 0.0, 0.5)
{
}
//_______________________
AliFemtoModelCorrFctnSource::AliFemtoModelCorrFctnSource(const char *title, Int_t aNbins, Double_t aQinvLo, Double_t aQinvHi):
  AliFemtoModelCorrFctn(title, aNbins, aQinvLo, aQinvHi),
  fHistROut(nullptr),
  fHistRSide(nullptr),
  fHistRLong(nullptr),
  fHistRStar(nullptr),
  fHistdNdR(nullptr),
  fHistNumWS(nullptr),
  fHistDenWS(nullptr),
  fUseRPSelection(0)
{
  // basic constructor
  char *outname = Form("%sOut", title);
  char *outtitle = Form("%sOut; R_{out}", title);
  fHistROut = new TH1D(outname,outtitle,100,-50.0,50.0);

  char *sidename = Form("%sSide", title);
  char *sidetitle = Form("%sSide; R_{side}", title);
  fHistRSide = new TH1D(sidename,sidetitle,100,-50.0,50.0);

  char *longname = Form("%sLong", title);
  char *longtitle = Form("%sLong; R_{long}", title);
  fHistRLong = new TH1D(longname,longtitle,100,-50.0,50.0);

  char *invname = Form("%sInv", title);
  char *invtitle = Form("%sInv; R_{inv}", title);
  fHistRStar = new TH1D(invname,invtitle,100,-50.0,50.0);

  char *dndrname = Form("%sdNdR", title);
  char *dndrtitle = Form("%sdNdR; R_{inv}; 1/R*^2", title);
  fHistdNdR = new TH1D(dndrname,dndrtitle,100,-50.0,50.0);

  char *nwsname = Form("%sNWS", title);
  char *nwstitle = Form("Weighted Numerator; q_{inv}; weight", title);
  fHistNumWS = new TH2D(nwsname,nwstitle,50,0.0,0.5,100,0.0,2.0);

  char *dwsname = Form("%sDWS", title);
  char *dwstitle = Form("Unweighted denominator; q_{inv}; weight");
  fHistDenWS = new TH2D(dwsname,dwstitle,50,0.0,0.5,100,0.0,2.0);

  fHistROut->Sumw2();
  fHistRSide->Sumw2();
  fHistRLong->Sumw2();
  fHistRStar->Sumw2();
  fHistdNdR->Sumw2();
}
//_______________________
AliFemtoModelCorrFctnSource::AliFemtoModelCorrFctnSource(const AliFemtoModelCorrFctnSource& aCorrFctn):
  AliFemtoModelCorrFctn(aCorrFctn),
  fHistROut(0),
  fHistRSide(0),
  fHistRLong(0),
  fHistRStar(0),
  fHistdNdR(0),
  fHistNumWS(0),
  fHistDenWS(0),
  fUseRPSelection(0)
{
  // copy constructor
  fHistROut = new TH1D (*aCorrFctn.fHistROut);
  fHistRSide = new TH1D(*aCorrFctn.fHistRSide);
  fHistRLong = new TH1D(*aCorrFctn.fHistRLong);
  fHistRStar = new TH1D(*aCorrFctn.fHistRStar);
  fHistdNdR = new TH1D(*aCorrFctn.fHistdNdR);
  fHistNumWS = new TH2D(*aCorrFctn.fHistNumWS);
  fHistDenWS = new TH2D(*aCorrFctn.fHistDenWS);

  fUseRPSelection = aCorrFctn.fUseRPSelection;
}
//_______________________
AliFemtoModelCorrFctnSource::~AliFemtoModelCorrFctnSource()
{
  // destructor
  if (fHistROut) delete fHistROut;
  if (fHistRSide) delete fHistRSide;
  if (fHistRLong) delete fHistRLong;
  if (fHistRStar) delete fHistRStar;
  if (fHistdNdR) delete fHistdNdR;
  if (fHistNumWS) delete fHistNumWS;
  if (fHistDenWS) delete fHistDenWS;
  if (fNumeratorTrue) delete fNumeratorTrue;
  if (fNumeratorFake) delete fNumeratorFake;
  if (fDenominator) delete fDenominator;
}

//_______________________
AliFemtoModelCorrFctnSource& AliFemtoModelCorrFctnSource::operator=(const AliFemtoModelCorrFctnSource& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

  *fHistROut = *aCorrFctn.fHistROut;
  *fHistRSide = *aCorrFctn.fHistRSide;
  *fHistRLong = *aCorrFctn.fHistRLong;

  *fHistRStar = *aCorrFctn.fHistRStar;
  *fHistdNdR = *aCorrFctn.fHistdNdR;
  *fHistNumWS = *aCorrFctn.fHistNumWS;
  *fHistDenWS = *aCorrFctn.fHistDenWS;

  fUseRPSelection = aCorrFctn.fUseRPSelection;

  return *this;
}
//_______________________
AliFemtoString AliFemtoModelCorrFctnSource::Report()
{
  // construct report
  AliFemtoString report = "AliFemtoModelCorrFctnSource report\n";

  return report;
}

//_______________________
void AliFemtoModelCorrFctnSource::AddRealPair(AliFemtoPair* aPair)
{
  // add real (effect) pair
  if (fPairCut){
    if (fUseRPSelection) {
      auto *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) {
        std::cout << "RP aware cut requested, but not connected to the CF\n";
        if (!(fPairCut->Pass(aPair))) {
          return;
        }
      }
      else {
        auto *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *>(HbtAnalysis());
        if (!arp) {
          std::cout << "RP aware cut requested, but not connected to the CF\n";
          if (!fPairCut->Pass(aPair)) {
            return;
          }
        }
        else if (!(ktc->Pass(aPair, arp->GetCurrentReactionPlane()))) {
          return;
        }
      }
    }
    else if (!fPairCut->Pass(aPair)) {
      return;
    }
  }

  AliFemtoModelCorrFctn::AddRealPair(aPair);
}
//_______________________
void AliFemtoModelCorrFctnSource::AddMixedPair(AliFemtoPair* aPair)
{
  // add mixed (background) pair
  if (fPairCut){
    if (fUseRPSelection) {
      auto *ktc = dynamic_cast<AliFemtoKTPairCut *>(fPairCut);
      if (!ktc) {
        std::cout << "RP aware cut requested, but not connected to the CF\n";
        if (!fPairCut->Pass(aPair)) {
          return;
        }
      }
      else {
        auto *arp = dynamic_cast<AliFemtoAnalysisReactionPlane *>(HbtAnalysis());
        if (!arp) {
          std::cout << "RP aware cut requested, but not connected to the CF\n";
          if (!fPairCut->Pass(aPair)) {
             return;
          }
        }
        else if (!ktc->Pass(aPair, arp->GetCurrentReactionPlane())) {
          return;
        }
      }
    }
    else if (!fPairCut->Pass(aPair)) {
      return;
    }
  }

  AliFemtoModelCorrFctn::AddMixedPair(aPair);
  // save the generated positions
  if (aPair->KStar() < 0.2) {
    fHistROut->Fill (fManager->GetWeightGenerator()->GetRStarOut());
    fHistRSide->Fill(fManager->GetWeightGenerator()->GetRStarSide());
    fHistRLong->Fill(fManager->GetWeightGenerator()->GetRStarLong());
    fHistRStar->Fill(fManager->GetWeightGenerator()->GetRStar());
    fHistdNdR->Fill (fManager->GetWeightGenerator()->GetRStar(), 1.0/(fManager->GetWeightGenerator()->GetRStar()*fManager->GetWeightGenerator()->GetRStar()));
  }

  fHistDenWS->Fill(aPair->QInv(), 1.0);
  Double_t weight = fManager->GetWeight(aPair);
  fHistNumWS->Fill(aPair->QInv(), weight);
}
//_______________________
void AliFemtoModelCorrFctnSource::Write()
{
  // write out all the histograms
  fHistROut->Write();
  fHistRSide->Write();
  fHistRLong->Write();
  fHistRStar->Write();
  fHistdNdR->Write();
  fHistNumWS->Write();
  fHistDenWS->Write();

  AliFemtoModelCorrFctn::Write();
}
//________________________
TList* AliFemtoModelCorrFctnSource::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = AliFemtoModelCorrFctn::GetOutputList();

  tOutputList->Add(fHistROut);
  tOutputList->Add(fHistRSide);
  tOutputList->Add(fHistRLong);
  tOutputList->Add(fHistRStar);
  tOutputList->Add(fHistdNdR);
  tOutputList->Add(fHistDenWS);
  tOutputList->Add(fHistNumWS);

  return tOutputList;
}
//_______________________
void AliFemtoModelCorrFctnSource::SetUseRPSelection(unsigned short aRPSel)
{
  fUseRPSelection = aRPSel;
}
