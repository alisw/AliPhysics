////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorPairKT - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//  rmaselek@cern.ch                                                          //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorPairKT.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>
#include "AliFemtoPair.h"


AliFemtoCutMonitorPairKT::AliFemtoCutMonitorPairKT():
  fkT(0),
  fpTsum(0),
  fpTdiff(0),
  fQinv(0),
  fMinv(0)
{
  // Default constructor
  fkT = new TH1D("kT", "k_T ", 100, 0, 4.0);
  fpTdiff = new TH1D("pTdiff", "k_T vs Q_inv", 100, 0, 4.0);
  fpTsum = new TH1D("pTsum", "p_T", 100, 0, 8.0);
  fMinv = new TH1D("Minv", "m_inv dist", 1000, 0.0, 10.0);
  fQinv = new TH1D("Qinv", "q_inv dist", 100, 0.0, 10.0);

}

AliFemtoCutMonitorPairKT::AliFemtoCutMonitorPairKT(const char *aName, float low_limit, float high_limit):
  fkT(0),
  fpTsum(0),
  fpTdiff(0),
  fQinv(0),
  fMinv(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "kT%s", aName);
  fkT = new TH1D(name, "k_T dist", 100, 0, 4.0);
  snprintf(name, 200, "pTdiff%s", aName);
  fpTdiff = new TH1D(name, "p_t diff dist", 100, 0, 4.0);
  snprintf(name, 200, "pTsum%s", aName);
  fpTsum = new TH1D(name, "p_T sum dist", 100, 0, 8.0);
  snprintf(name, 200, "fMinv%s", aName);
  fMinv = new TH1D(name, "m_inv dist", 1000, 0.0, 10.0);
  snprintf(name, 200, "fQinv%s", aName);
  fQinv = new TH1D(name, "q_inv dist", 1000, 0.0, 10.0);


}

AliFemtoCutMonitorPairKT::AliFemtoCutMonitorPairKT(const AliFemtoCutMonitorPairKT &aCut):
  AliFemtoCutMonitor(aCut),
  fkT(0),
  fpTsum(0),
  fpTdiff(0),
  fQinv(0),
  fMinv(0)
{
  // copy constructor
  fkT = new TH1D(*aCut.fkT);
  fpTdiff = new TH1D(*aCut.fpTdiff);
  fpTsum = new TH1D(*aCut.fpTsum);
  fMinv = new TH1D(*aCut.fMinv);
  fQinv = new TH1D(*aCut.fQinv);
}

AliFemtoCutMonitorPairKT::~AliFemtoCutMonitorPairKT()
{
  // Destructor
  delete fkT;
  delete fpTdiff;
  delete fpTsum;
  delete fMinv;
  delete fQinv;
}

AliFemtoCutMonitorPairKT& AliFemtoCutMonitorPairKT::operator=(const AliFemtoCutMonitorPairKT& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoCutMonitor::operator=(aCut);

  if (fkT) delete fkT;
  fkT = new TH1D(*aCut.fkT);
  if (fpTdiff) delete fpTdiff;
  fpTdiff = new TH1D(*aCut.fpTdiff);
  if (fpTsum) delete fpTsum;
  fpTsum = new TH1D(*aCut.fpTsum);
  if (fMinv) delete fMinv;
  fMinv = new TH1D(*aCut.fMinv);
  if (fQinv) delete fQinv;
  fQinv = new TH1D(*aCut.fQinv);

  return *this;
}

AliFemtoString AliFemtoCutMonitorPairKT::Report(){
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorPairKT report";
  AliFemtoString returnThis = stemp;
  return returnThis;
}

void AliFemtoCutMonitorPairKT::Fill(const AliFemtoPair* aPair)
{
   double pt1 = aPair->Track1()->Track()->Pt();
   double pt2 = aPair->Track2()->Track()->Pt();
   double kt = aPair->KT();
   double Q_inv = aPair -> QInv();
   double M_inv = aPair -> MInv();

   double pt_sum = pt1 + pt2;
   double pt_diff = TMath::Abs(pt1 - pt2);

  fpTsum -> Fill(pt_sum);
  fpTdiff -> Fill(pt_diff);
  fkT -> Fill(kt);
  fMinv -> Fill(M_inv);
  fQinv->Fill(Q_inv);

  // const PID1 int = aPair -> Track1() -> Track() -> GetPDGPid();
  // const PID2 int = aPair -> Track2() -> Track() -> GetPDGPid();
  // fPID->Fill(PID1);
  // fPID->Fill(PID2);
}

void AliFemtoCutMonitorPairKT::Write()
{
  // Write out the relevant histograms
  fpTsum->Write();
  fpTdiff->Write();
  fkT->Write();
  fMinv->Write();
  fQinv->Write();
}

TList *AliFemtoCutMonitorPairKT::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fpTsum);
  tOutputList->Add(fpTdiff);
  tOutputList->Add(fkT);
  tOutputList->Add(fMinv);
  tOutputList->Add(fQinv);

  return tOutputList;
}
