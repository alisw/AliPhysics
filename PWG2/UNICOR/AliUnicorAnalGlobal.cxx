/************************************************************************* 
* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. * 
*                                                                        * 
* Author: The ALICE Off-line Project.                                    * 
* Contributors are mentioned in the code where appropriate.              * 
*                                                                        * 
* Permission to use, copy, modify and distribute this software and its   * 
* documentation strictly for non-commercial purposes is hereby granted   * 
* without fee, provided that the above copyright notice appears in all   * 
* copies and that both the copyright notice and this permission notice   * 
* appear in the supporting documentation. The authors make no claims     * 
* about the suitability of this software for any purpose. It is          * 
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// event global variable analyzer
// Histogram the event variables like multiplicity, flow vector, vertex z, 
// etc. 
//=============================================================================

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include "AliUnicorEvent.h"
#include "AliUnicorAnalGlobal.h"


class TH2;
ClassImp(AliUnicorAnalGlobal)
  
//=============================================================================
AliUnicorAnalGlobal::AliUnicorAnalGlobal(Char_t *nam) : AliUnicorAnal(nam) 
{
  // constructor

  TH1D *mult = new TH1D("mult","mult",5000,-0.5,4999.5);
  mult->SetXTitle("multiplicity");
  TH1D *cent = new TH1D("cent","cent",100,0,1);
  cent->SetXTitle("centrality");
  TH2D *dire = new TH2D("dire","dire",100,-40,40,100,-40,40);
  dire->SetXTitle("Qx (GeV)");
  dire->SetYTitle("Qy (GeV)");
  TH1D *zver = new TH1D("zver","zver",120,-1.2,1.2);
  zver->SetXTitle("normalized z-vertex");
  fHistos.Add(mult);
  fHistos.Add(cent);
  fHistos.Add(dire);
  fHistos.Add(zver);
  gROOT->cd();
}
//=============================================================================
void AliUnicorAnalGlobal::Process(AliUnicorEvent *ev) const
{
  // fill event variable histograms

  TH1D *mult = (TH1D*) fHistos.At(0);
  TH1D *cent = (TH1D*) fHistos.At(1);
  TH2D *dire = (TH2D*) fHistos.At(2);
  TH1D *zver = (TH1D*) fHistos.At(3);

  mult->Fill(ev->NParticles(),1.0);
  cent->Fill(ev->Centrality(),1.0);
  Double_t qx,qy;
  ev->RP(qx,qy);
  dire->Fill(qx,qy,1.0);
  zver->Fill(ev->Zver(),1.0);
}
//=============================================================================
