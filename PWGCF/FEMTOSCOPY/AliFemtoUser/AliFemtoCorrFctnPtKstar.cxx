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
  fPtKstarDen2part(0)
{
  fPtKstar = new TH2D(Form("PtvsKstar1part%s",title),"Pt vs kstar (part 1)",200,0.0,2.0, 200, 0.0, 4.0);
  fPtKstarDen = new TH2D(Form("PtvsKstarDen1part%s",title),"Pt vs kstar in mixed events (part 1)",200,0.0,2.0, 200, 0.0, 4.0);

  fPtKstar2part = new TH2D(Form("PtvsKstar2part%s",title),"Pt vs kstar (part 2)",200,0.0,2.0, 200, 0.0, 4.0);
  fPtKstarDen2part = new TH2D(Form("PtvsKstarDen2part%s",title),"Pt vs kstar in mixed events (part 2)",200,0.0,2.0, 200, 0.0, 4.0);
}

//____________________________
AliFemtoCorrFctnPtKstar::AliFemtoCorrFctnPtKstar(const AliFemtoCorrFctnPtKstar& aCorrFctn) :
  AliFemtoCorrFctn(),
   fPtKstar(0),
  fPtKstarDen(0),
  fPtKstar2part(0),
  fPtKstarDen2part(0)
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
 double px = pair->Track1()->Track()->P().x();
 double py = pair->Track1()->Track()->P().y();
 double pT = TMath::Hypot(px, py);

 
 fPtKstar->Fill(tKStar,pT);
 double px2 = pair->Track2()->Track()->P().x();
 double py2 = pair->Track2()->Track()->P().y();
 double pT2 = TMath::Hypot(px2, py2);
 fPtKstar2part->Fill(tKStar,pT2);

}
//____________________________
void AliFemtoCorrFctnPtKstar::AddMixedPair( AliFemtoPair* pair){
 double tKStar = fabs(pair->KStar());
 double px = pair->Track1()->Track()->P().x();
 double py = pair->Track1()->Track()->P().y();
 double pT = TMath::Hypot(px, py);
 
 fPtKstarDen->Fill(tKStar,pT);
 
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
}

TList* AliFemtoCorrFctnPtKstar::GetOutputList()
{
  // Prepare the list of objects to be written to the output
  TList *tOutputList = new TList();

  tOutputList->Add(fPtKstar);
  tOutputList->Add(fPtKstarDen);
  tOutputList->Add(fPtKstar2part);
  tOutputList->Add(fPtKstarDen2part);

  return tOutputList;
}
