#include "AliFemtoCutMonitorPairMomRes.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>
#include "AliFemtoPair.h"


AliFemtoCutMonitorPairMomRes::AliFemtoCutMonitorPairMomRes():
  fMomRes(0),
  fMomResTrueMass(0),
  fMomRes_KPpairOnly(0),
  fMomResTrueMass_KPpairOnly(0),
  fMomRes_pppairOnly(0),
  fMomResTrueMass_pppairOnly(0),
  fMassPart1(0),
  fMassPart2(0)
{
  // Default constructor
  fMomRes = new TH2D("PairMomRes","PairMomRes",100,0,2,100,0,2);
  fMomResTrueMass = new TH2D("PairMomResTrueMass","PairMomResTrueMass",100,0,2,100,0,2);
  fMomRes_KPpairOnly = new TH2D("PairMomRes_KPpairOnly","PairMomRes_KPpairOnly",100,0,2,100,0,2);
  fMomResTrueMass_KPpairOnly = new TH2D("PairMomResTrueMass_KPpairOnly","PairMomResTrueMass_KPpairOnly",100,0,2,100,0,2);
  fMomRes_pppairOnly = new TH2D("PairMomRes_pppairOnly","PairMomRes_pppairOnly",100,0,2,100,0,2);
  fMomResTrueMass_pppairOnly = new TH2D("PairMomResTrueMass_pppairOnly","PairMomResTrueMass_pppairOnly",100,0,2,100,0,2);
}

AliFemtoCutMonitorPairMomRes::AliFemtoCutMonitorPairMomRes(const char *aName, double massPart1, double massPart2, double qmin, double qmax, int nbins):
  fMomRes(0),
  fMassPart1(massPart1),
  fMassPart2(massPart2)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "PairMomRes%s", aName);
  fMomRes = new TH2D(name, "PairMomRes", nbins, qmin, qmax, nbins, qmin, qmax);

  snprintf(name, 200, "PairMomResTrueMass%s", aName);
  fMomResTrueMass = new TH2D(name, "PairMomResTrueMass", nbins, qmin, qmax, nbins, qmin, qmax);

  snprintf(name, 200, "PairMomRes_KPpairOnly%s", aName);
  fMomRes_KPpairOnly = new TH2D(name, "PairMomRes_KPpairOnly", nbins, qmin, qmax, nbins, qmin, qmax);

  snprintf(name, 200, "PairMomResTrueMass_KPpairOnly%s", aName);
  fMomResTrueMass_KPpairOnly = new TH2D(name, "PairMomResTrueMass_KPpairOnly", nbins, qmin, qmax, nbins, qmin, qmax);

  snprintf(name, 200, "PairMomResModifMass_pppairOnly%s", aName);
  fMomRes_pppairOnly = new TH2D(name, "PairMomRes_pppairOnly", nbins, qmin, qmax, nbins, qmin, qmax);

  snprintf(name, 200, "PairMomResTrueMass_pppairOnly%s", aName);
  fMomResTrueMass_pppairOnly = new TH2D(name, "PairMomResTrueMass_pppairOnly", nbins, qmin, qmax, nbins, qmin, qmax);


}

AliFemtoCutMonitorPairMomRes::AliFemtoCutMonitorPairMomRes(const AliFemtoCutMonitorPairMomRes &aCut):
  AliFemtoCutMonitor(aCut),
  fMomRes(0),
  fMomResTrueMass(0),
  fMomRes_KPpairOnly(0),
  fMomResTrueMass_KPpairOnly(0),
  fMomRes_pppairOnly(0),
  fMomResTrueMass_pppairOnly(0),
  fMassPart1(aCut.fMassPart1),
  fMassPart2(aCut.fMassPart2)
{
  // copy constructor
  fMomRes = new TH2D(*aCut.fMomRes);
  fMomResTrueMass = new TH2D(*aCut.fMomResTrueMass);
  fMomRes_KPpairOnly = new TH2D(*aCut.fMomRes_KPpairOnly);
  fMomResTrueMass_KPpairOnly = new TH2D(*aCut.fMomResTrueMass_KPpairOnly);
  fMomRes_pppairOnly = new TH2D(*aCut.fMomRes_pppairOnly);
  fMomResTrueMass_pppairOnly = new TH2D(*aCut.fMomResTrueMass_pppairOnly);
}

AliFemtoCutMonitorPairMomRes::~AliFemtoCutMonitorPairMomRes()
{
  // Destructor
  delete fMomRes;
  delete fMomResTrueMass;
  delete fMomRes_KPpairOnly;
  delete fMomResTrueMass_KPpairOnly;
  delete fMomRes_pppairOnly;
  delete fMomResTrueMass_pppairOnly;
}

AliFemtoCutMonitorPairMomRes& AliFemtoCutMonitorPairMomRes::operator=(const AliFemtoCutMonitorPairMomRes& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoCutMonitor::operator=(aCut);

  if (fMomRes) delete fMomRes;
  fMomRes = new TH2D(*aCut.fMomRes);

  if (fMomResTrueMass) delete fMomResTrueMass;
  fMomResTrueMass = new TH2D(*aCut.fMomResTrueMass);

  if (fMomRes_KPpairOnly) delete fMomRes_KPpairOnly;
  fMomRes_KPpairOnly = new TH2D(*aCut.fMomRes_KPpairOnly);

  if (fMomResTrueMass_KPpairOnly) delete fMomResTrueMass_KPpairOnly;
  fMomResTrueMass_KPpairOnly = new TH2D(*aCut.fMomResTrueMass_KPpairOnly);

  if (fMomRes_pppairOnly) delete fMomRes_pppairOnly;
  fMomRes_pppairOnly = new TH2D(*aCut.fMomRes_pppairOnly);

  if (fMomResTrueMass_pppairOnly) delete fMomResTrueMass_pppairOnly;
  fMomResTrueMass_pppairOnly = new TH2D(*aCut.fMomResTrueMass_pppairOnly);

  fMassPart1 = aCut.fMassPart1;
  fMassPart2 = aCut.fMassPart2;

  return *this;
}

AliFemtoString AliFemtoCutMonitorPairMomRes::Report(){
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorMomRes report";
  AliFemtoString returnThis = stemp;
  return returnThis;
}

void AliFemtoCutMonitorPairMomRes::Fill(const AliFemtoPair* aPair)
{
   double KStar_rec = aPair ->KStar();

  AliFemtoParticle *first = (AliFemtoParticle*)aPair->Track1();
  AliFemtoParticle *second = (AliFemtoParticle*)aPair->Track2();
  
  if(first!=NULL && second!=NULL) {

    
    AliFemtoModelHiddenInfo *tInfo1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo();
    AliFemtoModelHiddenInfo *tInfo2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();

    if(tInfo1!=NULL && tInfo2!=NULL) {

      AliFemtoLorentzVector p1(::sqrt(tInfo1->GetTrueMomentum()->Mag2() + fMassPart1*fMassPart1), *tInfo1->GetTrueMomentum());
      AliFemtoLorentzVector p2(::sqrt(tInfo2->GetTrueMomentum()->Mag2() + fMassPart2*fMassPart2), *tInfo2->GetTrueMomentum());

      AliFemtoLorentzVector p1true(::sqrt(tInfo1->GetTrueMomentum()->Mag2() + tInfo1->GetMass()*tInfo1->GetMass()), *tInfo1->GetTrueMomentum());
      AliFemtoLorentzVector p2true(::sqrt(tInfo2->GetTrueMomentum()->Mag2() + tInfo2->GetMass()*tInfo2->GetMass()), *tInfo2->GetTrueMomentum());
    

      //mass set in constructor
      const double px1 = p1.x();
      const double py1 = p1.y();
      const double pz1 = p1.z();
      const double pE1  = p1.e();
      const double mass1_sqrd = std::max({0.0, p1.m2()});

      const double px2 = p2.x();
      const double py2 = p2.y();
      const double pz2 = p2.z();
      const double pE2  = p2.e();
      const double mass2_sqrd = std::max({0.0, p2.m2()});

      const double tPx = px1 + px2;
      const double  tPy = py1 + py2;
      const double  tPz = pz1 + pz2;
      const double  tPE = pE1 + pE2;

      double tPtrans = tPx*tPx + tPy*tPy;
      double tMtrans = tPE*tPE - tPz*tPz;
      if( (tMtrans - tPtrans) > 0)
	{
	  double tPinv = ::sqrt(tMtrans - tPtrans);
	  tMtrans = ::sqrt(tMtrans);
	  tPtrans = ::sqrt(tPtrans);

	  double tQinvL = (pE1-pE2)*(pE1-pE2) - (px1-px2)*(px1-px2) -
	    (py1-py2)*(py1-py2) - (pz1-pz2)*(pz1-pz2);
	  
	  
	  double tQ = (mass1_sqrd - mass2_sqrd)/tPinv;
	  
	  if((tQ*tQ - tQinvL)>=0.000000001)
	    {
	      tQ = ::sqrt( tQ*tQ - tQinvL);
	      
	      double KStar_true = tQ/2;
	      
	      if(KStar_true!=0 || KStar_rec!=0)

   
		fMomRes->Fill(KStar_rec,KStar_true);
		

	            if((abs(tInfo1->GetPDGPid()) == 321) && (abs(tInfo2->GetPDGPid()) == 2212))
		      {
			fMomRes_KPpairOnly->Fill(KStar_rec,KStar_true);
		      }

	
	    }
	}


      //true mass from model
      const double px1true = p1true.x();
      const double py1true = p1true.y();
      const double pz1true = p1true.z();
      const double pE1true  = p1true.e();
      const double mass1_sqrdtrue = std::max({0.0, p1true.m2()});

      const double px2true = p2true.x();
      const double py2true = p2true.y();
      const double pz2true = p2true.z();
      const double pE2true  = p2true.e();
      const double mass2_sqrdtrue = std::max({0.0, p2true.m2()});

      const double tPxtrue = px1true + px2true;
      const double  tPytrue = py1true + py2true;
      const double  tPztrue = pz1true + pz2true;
      const double  tPEtrue = pE1true + pE2true;

      double tPtranstrue = tPxtrue*tPxtrue + tPytrue*tPytrue;
      double tMtranstrue = tPEtrue*tPEtrue - tPztrue*tPztrue;
      if( (tMtranstrue - tPtranstrue) > 0)
	{
	  double tPinvtrue = ::sqrt(tMtranstrue - tPtranstrue);
	  tMtranstrue = ::sqrt(tMtranstrue);
	  tPtranstrue = ::sqrt(tPtranstrue);

	  double tQinvLtrue = (pE1true-pE2true)*(pE1true-pE2true) - (px1true-px2true)*(px1true-px2true) -
	    (py1true-py2true)*(py1true-py2true) - (pz1true-pz2true)*(pz1true-pz2true);
	  
	  
	  double tQtrue = (mass1_sqrdtrue - mass2_sqrdtrue)/tPinvtrue;
	  
	  if((tQtrue*tQtrue - tQinvLtrue)>=0.000000001)
	    {
	      tQtrue = ::sqrt( tQtrue*tQtrue - tQinvLtrue);
	      
	      double KStar_truetrue = tQtrue/2;
	      
	      if(KStar_truetrue!=0 || KStar_rec!=0)

   
		fMomResTrueMass->Fill(KStar_rec,KStar_truetrue);


	            if((abs(tInfo1->GetPDGPid()) == 321) && (abs(tInfo2->GetPDGPid()) == 2212))
		      {
			fMomResTrueMass_KPpairOnly->Fill(KStar_rec,KStar_truetrue);
		      }
		    if((abs(tInfo1->GetPDGPid()) == 2212) && (abs(tInfo2->GetPDGPid()) == 2212))
		      {
			fMomResTrueMass_pppairOnly->Fill(KStar_rec,KStar_truetrue);
		      }
	
	    }
	}
	
	
      //mass from constructor deuterons stuffs, test/check, Wioleta Rzęsa wioleta.rzesa@cern.ch 
      if((abs(tInfo1->GetPDGPid()) == 2212) && (abs(tInfo2->GetPDGPid()) == 2212))
        {
          double tPtrans_deu = tPx*tPx + tPy*tPy;
          double tMtrans_deu = tPE*tPE - tPz*tPz;
          if( (tMtrans_deu - tPtrans_deu) > 0)
	    {
	      double tPinv = ::sqrt(tMtrans_deu - tPtrans_deu);
	      tMtrans_deu = ::sqrt(tMtrans_deu);
	      tPtrans_deu = ::sqrt(tPtrans_deu);

	      double tQinvL = (pE1-pE2)*(pE1-pE2) - (px1-px2)*(px1-px2) -
	      (py1-py2)*(py1-py2) - (pz1-pz2)*(pz1-pz2);
	 
	      double tQ = (mass1_sqrd - mass2_sqrd)/tPinv;
	  
	      if((tQ*tQ - tQinvL)>=0.000000001)
	        {
	          tQ = ::sqrt( tQ*tQ - tQinvL);
	          double KStar_true = tQ/2;
	        
                  const AliFemtoLorentzVector
                    &p1_re = aPair->Track1()->FourMomentum(),
                    &p2_re = aPair->Track2()->FourMomentum();
                  const double
                    px1_re = p1_re.x(),
                    py1_re = p1_re.y(),
                    pz1_re = p1_re.z(),
                    pE1_re = p1_re.e(),
                    mass1_sqrd_re = std::max({0.0, p1_re.m2()}),

                    px2_re = p2_re.x(),
                    py2_re = p2_re.y(),
                    pz2_re = p2_re.z();
                    AliFemtoLorentzVector p2_re2(::sqrt(second->Track()->P().Mag2() + fMassPart2*fMassPart2), second->Track()->P());
                    const double 
                      pE2_re = p2_re2.e(),
                      mass2_sqrd_re = std::max({0.0, p2_re2.m2()});
                      
                    const double 
                      tPx_re = px1_re + px2_re,
                      tPy_re = py1_re + py2_re,
                      tPz_re = pz1_re + pz2_re,
                      tPE_re = pE1_re + pE2_re;

                    double tPtrans_re = tPx_re*tPx_re + tPy_re*tPy_re;
                    double tMtrans_re = tPE_re*tPE_re - tPz_re*tPz_re;
                    double tPinv_re = ::sqrt(tMtrans_re - tPtrans_re);
                    tMtrans_re = ::sqrt(tMtrans_re);
                    tPtrans_re = ::sqrt(tPtrans_re);

                    double tQinvL_re = (pE1_re-pE2_re)*(pE1_re-pE2_re) - (px1_re-px2_re)*(px1_re-px2_re) -
                    (py1_re-py2_re)*(py1_re-py2_re) - (pz1_re-pz2_re)*(pz1_re-pz2_re);

                    double tQ_re = (mass1_sqrd_re - mass2_sqrd_re)/tPinv_re;
                    tQ_re = ::sqrt( tQ_re*tQ_re - tQinvL_re);

                    double KStar_rec2 = tQ_re/2;
	            if(KStar_true!=0 || KStar_rec2!=0)
			fMomRes_pppairOnly->Fill(KStar_rec2,KStar_true);		   
	
	          }
	    }
        }
        
    }
  }
  
  
}

void AliFemtoCutMonitorPairMomRes::Write()
{
  // Write out the relevant histograms
  fMomRes->Write();
  fMomResTrueMass->Write();
  fMomRes_KPpairOnly->Write();
  fMomResTrueMass_KPpairOnly->Write();
  fMomRes_pppairOnly->Write();
  fMomResTrueMass_pppairOnly->Write();
}

TList *AliFemtoCutMonitorPairMomRes::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fMomRes);
  tOutputList->Add(fMomResTrueMass);
  tOutputList->Add(fMomRes_KPpairOnly);
  tOutputList->Add(fMomResTrueMass_KPpairOnly);
  tOutputList->Add(fMomRes_pppairOnly);
  tOutputList->Add(fMomResTrueMass_pppairOnly);

  return tOutputList;
}
