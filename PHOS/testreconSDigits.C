#include "AliPHOSGetter.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "AliPHOSHit.h"
#include "TFolder.h"
#include "TStopwatch.h"
#include "TObjArray.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSSDigitizer.h"

void testreconSDigits(Int_t nevent = 1, const char *config="testconfig.C")
{
  cerr<<" ___________________________________________________________________ "<<endl;
  cerr<<" "<<endl;
  cerr<<"           MESS ==> Beginning of the PHOS reconstruction. "<<endl;
  cerr<<" ___________________________________________________________________ "<<endl;
  const Float_t maxSDigits = 0 ;
  cerr<<"La valeur de maxSDigits est "<<maxSDigits<<endl;
  const Float_t widSDigits = TMath::Sqrt(maxSDigits) ;
  cerr<<"La valeur de widSDigits est "<<widSDigits<<endl;
  AliPHOSSDigitizer *sd = new AliPHOSSDigitizer("galice.root","test suite");
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  sd->ExecuteTask("deb"); 
  cerr << "# of SDigits " << gime->SDigits()->GetEntries() / gime->MaxEvent() << endl;
  Float_t nSDigits =  (Float_t) (gime->SDigitizer()->GetSDigitsInRun()) / gime->MaxEvent();
  cerr<<"La valeur de nSDigits est "<<nSDigits<<endl;
  
//   if ( gime->SDigits()->GetEntries() < maxSDigits-widSDigits ||
//        gime->SDigits()->GetEntries() > maxSDigits+widSDigits ) {
//     Erreur
//       }
}
