////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnPtKstar - A correlation function that analyzes        //
// two particle mass minvariant with various mass assumptions                 //
//                                                                            //
// Authors: Małgorzata Janik majanik@cern.ch
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoCorrFctnPtKstar.h"
#include <cstdio>
#include <TMath.h>

#ifdef __ROOT__
ClassImp(AliFemtoCorrFctnPtKstar)
#endif

//____________________________
AliFemtoCorrFctnPtKstar::AliFemtoCorrFctnPtKstar(const char* title):
AliFemtoCorrFctn(),
  fPtKstar(0),
  fPtKstarDen(0),
  fPtKstar2part(0),
  fPtKstarDen2part(0),
  fPairPtKstar2part(0),
  fPairPtKstarDen2part(0),
  fKstarBetaT(0)
{
  fPtKstar = new TH2D(Form("PtvsKstar1part%s",title),"Pt vs kstar (part 1)",200,0.0,2.0, 200, 0.0, 4.0);
  fPtKstarDen = new TH2D(Form("PtvsKstarDen1part%s",title),"Pt vs kstar in mixed events (part 1)",200,0.0,2.0, 200, 0.0, 4.0);

  fPtKstar2part = new TH2D(Form("PtvsKstar2part%s",title),"Pt vs kstar (part 2)",200,0.0,2.0, 200, 0.0, 4.0);
  fPtKstarDen2part = new TH2D(Form("PtvsKstarDen2part%s",title),"Pt vs kstar in mixed events (part 2)",200,0.0,2.0, 200, 0.0, 4.0);

  fPairPtKstar2part = new TH2D(Form("PairPtvsKstar%s",title),"Pair Pt vs kstar ",100,0.0,0.5, 300, 0.0, 6.0);
  fPairPtKstarDen2part = new TH2D(Form("PairPtvsKstarDen%s",title),"Pair Pt vs kstar in mixed events ",100,0.0,0.5, 300, 0.0, 6.0);

  fKstarBetaT = new TH2D(Form("KstarvsBetaT%s",title),"kstar vs BetaT ",100,0.0,0.5, 30, 0.0, 3.0);
  }

//____________________________
AliFemtoCorrFctnPtKstar::AliFemtoCorrFctnPtKstar(const AliFemtoCorrFctnPtKstar& aCorrFctn) :
  AliFemtoCorrFctn(),
   fPtKstar(0),
  fPtKstarDen(0),
  fPtKstar2part(0),
  fPtKstarDen2part(0),
  fPairPtKstar2part(0),
  fPairPtKstarDen2part(0),
  fKstarBetaT(0)
{
  // copy constructor
  if (fPtKstar) delete fPtKstar;
  fPtKstar = new TH2D(*aCorrFctn.fPtKstar);

  if (fPtKstarDen) delete fPtKstarDen;
  fPtKstarDen = new TH2D(*aCorrFctn.fPtKstarDen);
}
//____________________________
AliFemtoCorrFctnPtKstar::~AliFemtoCorrFctnPtKstar(){
  // destructor
    delete  fPtKstar;
    delete  fPtKstarDen;
    delete  fPtKstar2part;
    delete  fPtKstarDen2part;
    delete  fPairPtKstar2part;
    delete  fPairPtKstarDen2part;
    delete  fKstarBetaT;
}
//_________________________
AliFemtoCorrFctnPtKstar& AliFemtoCorrFctnPtKstar::operator=(const AliFemtoCorrFctnPtKstar& aCorrFctn)
{
  // assignment operator
  if (this == &aCorrFctn)
    return *this;

   if (fPtKstar) delete fPtKstar;
  fPtKstar = new TH2D(*aCorrFctn.fPtKstar);
  
  if (fPtKstarDen) delete fPtKstarDen;
  fPtKstarDen = new TH2D(*aCorrFctn.fPtKstarDen);

  if (fPtKstar2part) delete fPtKstar2part;
  fPtKstar2part = new TH2D(*aCorrFctn.fPtKstar2part);
  
  if (fPtKstarDen2part) delete fPtKstarDen2part;
  fPtKstarDen2part = new TH2D(*aCorrFctn.fPtKstarDen2part);
  
  if (fPairPtKstar2part) delete fPairPtKstar2part;
  fPairPtKstar2part = new TH2D(*aCorrFctn.fPairPtKstar2part);
  
  if (fPairPtKstarDen2part) delete fPairPtKstarDen2part;
  fPairPtKstarDen2part = new TH2D(*aCorrFctn.fPairPtKstarDen2part);

  if(fKstarBetaT) delete fKstarBetaT;
  fKstarBetaT = new TH2D(*aCorrFctn.fKstarBetaT);
  
  return *this;
}
//_________________________
void AliFemtoCorrFctnPtKstar::Finish(){
  // here is where we should normalize, fit, etc...
  // we should NOT Draw() the histos (as I had done it below),
  // since we want to insulate ourselves from root at this level
  // of the code.  Do it instead at root command line with browser.
  //  mShareNumerator->Draw();
  //mShareDenominator->Draw();
  //mRatio->Draw();

}

//____________________________
AliFemtoString AliFemtoCorrFctnPtKstar::Report(){
  // create report
  string stemp = "Kstar vs Pt Monitor Report\n";
  AliFemtoString returnThis = stemp;
  return returnThis;
}
//____________________________
void AliFemtoCorrFctnPtKstar::AddRealPair( AliFemtoPair* pair){

 double tKStar = fabs(pair->KStar());
 double tPairPt = fabs(pair->KT());
  
 double px = pair->Track1()->Track()->P().x();
 double py = pair->Track1()->Track()->P().y();
 double pT = TMath::Hypot(px, py);

 
 fPtKstar->Fill(tKStar,pT);
 double px2 = pair->Track2()->Track()->P().x();
 double py2 = pair->Track2()->Track()->P().y();
 double pT2 = TMath::Hypot(px2, py2);
 fPtKstar2part->Fill(tKStar,pT2);

 fPairPtKstar2part->Fill(tKStar,tPairPt);


  //--------Pritam's addition----------
 //--------Taken from AliFemtoBetaTPairCut.cxx

 
   // Calculate transverse momentum of the pair:
  double px1 = pair->Track1()->Track()->P().x();
  double px22 = pair->Track2()->Track()->P().x();
  double py1 = pair->Track1()->Track()->P().y();
  double py22 = pair->Track2()->Track()->P().y();
  double pxpair = px1 + px22;
  double pypair = py1 + py22;
  double pTpair = TMath::Sqrt(pxpair*pxpair + pypair*pypair);
  // Calculate energies of particles:
  double pz1 = pair->Track1()->Track()->P().z();
  double pz2 = pair->Track2()->Track()->P().z();
  double pzpair = pz1 + pz2;
  double p1 = TMath::Sqrt(px1*px1 + py1*py1 + pz1*pz1);
  double p2 = TMath::Sqrt(px2*px2 + py2*py2 + pz2*pz2);
  double m1 = 0.13957;  //pion mass
  double m2 = 0.493677; //kaon mass
  double e1 = TMath::Sqrt(p1*p1 + m1*m1);
  double e2 = TMath::Sqrt(p2*p2 + m2*m2);
  // Calculate transverse mass of the pair:
  double mInvpair_2 = m1*m1 + m2*m2 + 2*(e1*e2 - px1*px2 - py1*py2 - pz1*pz2);
  double mTpair = TMath::Sqrt(mInvpair_2 + pTpair*pTpair);
  // Calculate betaT:
  double betaT = pTpair / mTpair;

  fKstarBetaT->Fill(tKStar,betaT);
  
}
//____________________________
void AliFemtoCorrFctnPtKstar::AddMixedPair( AliFemtoPair* pair){
 double tKStar = fabs(pair->KStar());
 double tPairPt = fabs(pair->KT());  
  
 double px = pair->Track1()->Track()->P().x();
 double py = pair->Track1()->Track()->P().y();
 double pT = TMath::Hypot(px, py);
 
 fPtKstarDen->Fill(tKStar,pT);
 fPairPtKstarDen2part->Fill(tKStar,tPairPt);  
 
 double px2 = pair->Track2()->Track()->P().x();
 double py2 = pair->Track2()->Track()->P().y();
 double pT2 = TMath::Hypot(px2, py2);
 fPtKstarDen2part->Fill(tKStar,pT2);


}
//____________________________

void AliFemtoCorrFctnPtKstar::WriteHistos()
{
  // Write out result histograms
  fPtKstar->Write();
  fPtKstarDen->Write();
  fPtKstar2part->Write();
  fPtKstarDen2part->Write();
  fPairPtKstar2part->Write();
  fPairPtKstarDen2part->Write();
  fKstarBetaT->Write();
}

TList* AliFemtoCorrFctnPtKstar::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fPtKstar);
  tOutputList->Add(fPtKstarDen);
  tOutputList->Add(fPtKstar2part);
  tOutputList->Add(fPtKstarDen2part);
  tOutputList->Add(fPairPtKstar2part);
  tOutputList->Add(fPairPtKstarDen2part);
  tOutputList->Add(fKstarBetaT);
  return tOutputList;
}
