///
/// \file AliFemtoPionLambdaCutMonitor.cxx
///


#include "AliFemtoPionLambdaCutMonitor.h"

#include "AliFemtoEvent.h"


static const double PionMass = 0.13956995;

#include <TH1F.h>
#include <TH2F.h>
#include <TObjArray.h>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoPionLambdaCutMonitor);
  /// \endcond
#endif

AliFemtoPionLambdaCutMonitor::Event::Event(const bool passing):
  _centrality(NULL),
  _multiplicity(NULL),
  _vertex_z(NULL),
  _vertex_xy(NULL)
{
  const char *title_suffix = (passing ? " (PASS)" : " (FAIL)");

  const TString pf(passing ? "_P" : "_F");

  _centrality = new TH1F(
    "centrality" + pf,
    TString::Format("Event Centrality%s", title_suffix),
    100, 0, 100.0
  );

  _multiplicity = new TH1F(
    "multiplicity" + pf,
    "Event Multiplicity",
    100, 0, 10000.0
  );

  _vertex_z = new TH1F(
    "VertexZ" + pf,
    TString::Format("Vertex Z Distribution%s;z (cm);dN/dz", title_suffix),
    128, -15.0f, 15.0f
  );

  _vertex_xy = new TH2F(
    "VertexXY" + pf,
    TString::Format("Vertex XY Distribution%s;x (cm);y (cm); dN/(dx $\\cdot$ dy)", title_suffix),
    48, -1.0f, 1.0f,
    48, -1.0f, 1.0f
  );
}

TList*
AliFemtoPionLambdaCutMonitor::Event::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(_centrality);
  output->Add(_multiplicity);
  output->Add(_vertex_z);
  output->Add(_vertex_xy);

  return olist;
}

void 
AliFemtoPionLambdaCutMonitor::Event::Fill(const AliFemtoEvent* ev)
{
  const Float_t centrality = ev->CentralityV0();
  const Int_t multiplicty = ev->NumberOfTracks();
  const AliFemtoThreeVector vertex = ev->PrimVertPos();

  _centrality->Fill(centrality);
  _multiplicity->Fill(multiplicty);
  _vertex_z->Fill(vertex.z());
  _vertex_xy->Fill(vertex.x(), vertex.y());
}


AliFemtoPionLambdaCutMonitor::Pion::Pion(const bool passing, const TString& typestr):
  _minv(NULL),
  _ypt(NULL),
  fPtPhi(NULL),
  fEtaPhi(NULL)
{
  // Build 'standard' format for histogram titles
  //  <ParticleType> <Title> <Pass/Fail>; <AxisInfo>
  const TString title_format = TString::Format("%s %%s %s; %%s",
                                               typestr.Data(),
                                               (passing ? "(PASS)" : "(FAIL)"));
  const TString pf(passing ? "_P" : "_F");

  _minv = new TH1F(
    "Mass" + pf,
    TString::Format(title_format, "M_{inv}",
                                  "M_{inv} (GeV);"
                                  "dN/dM"),
    144, 0.120, .158);

  _ypt = new TH2F(
    "eta_Pt" + pf,
    TString::Format(title_format,
                    "\\eta vs p_{T}",
             /*X*/  "$\\eta$;"
             /*Y*/  "p_{T} (GeV);"
             /*Z*/  "dN/(p_{T} $\\cdot$ \\eta)"),
     140, -1.4, 1.4,
     100, 0, 3.0);

  fPtPhi = new TH2F(
    "PtPhi" + pf,
    TString::Format(title_format,
                    "Pt vs Phi",
                    "p_{T} (GeV);"
                    "Phi (rads)"),
    144,  0.0, 3.0,
    100, -TMath::Pi(), TMath::Pi());

  fEtaPhi = new TH2F(
    "EtaPhi" + pf,
    TString::Format(title_format,
                    "\\eta vs Phi",
                    "\\eta;"
                    "Phi (rads)"),
    144, -1.4, 1.4,
    144, -TMath::Pi(), TMath::Pi());
}


TList*
AliFemtoPionLambdaCutMonitor::Pion::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(_minv);
  output->Add(_ypt);
  output->Add(fPtPhi);
  output->Add(fEtaPhi);

  return olist;
}

void AliFemtoPionLambdaCutMonitor::Pion::Fill(const AliFemtoTrack* track)
{
  const float pz = track->P().z(),
              pt = track->Pt(),
             phi = track->P().Phi();
  const double energy = ::sqrt(track->P().Mag2() + PionMass * PionMass),
                  eta = 0.5 * ::log((energy + pz) / (energy - pz));


  _minv->Fill(track->GetMass());
  _ypt->Fill(eta, pt);
  fPtPhi->Fill(pt, phi);
  fEtaPhi->Fill(eta, phi);
}

AliFemtoPionLambdaCutMonitor::Lambda::Lambda(const bool passing,
                                             const TString& typestr,
                                             const AliFemtoAnalysisPionLambda::LambdaType ltype):
  fLambdaType(ltype),
  _minv(NULL),
  _ypt(NULL),
  _dedx_pt_pro(NULL),
  _dedx_pt_pi(NULL)
{
  // Build 'standard' format for histogram titles
  //  <ParticleType> <Title> <Pass/Fail>; <AxisInfo>
  const TString title_format = TString::Format("%s %%s %s; %%s",
                                         typestr.Data(),
                                         (passing ? "(PASS)" : "(FAIL)"));
  const TString pf(passing ? "_P" : "_F");

  _minv = new TH1F(
    "V0_Minv" + pf,
    TString::Format(title_format,
                    "M_{inv}",
                    "M_{inv} (GeV); dN/dM"),
    576, 1.07, 1.118);

  _ypt = new TH2F(
    "V0_YPt" + pf,
    TString::Format(title_format,
                    "\\eta vs p_{T}",
                    "\\eta;"
                    "p_{T} (GeV);"
                    "dN/(p_{T} $\\cdot$ \\eta)"),
    140, -1.4, 1.4,
    100, 0.0, 3.0);

  _dedx_pt_pro = new TH2F(
    "dEdX_Pro" + pf,
    TString::Format(title_format,
                    "dE/dx vs p_{T} (Proton Daughter)",
                    "p_{T} (GeV);"
                    "dE/dx;"
                    "dN/(p_{T} $\\cdot$ dE/dx)"),
     128, 0, 2.2,
     128, 0, 500);

  _dedx_pt_pi = new TH2F(
    "dEdX_Pi" + pf,
    TString::Format(title_format,
                    "dE/dx vs p_{T} (Pion Daughter)",
                    "p_{T} (GeV);"
                    "dE/dx;"
                    "dN/(p_{T} $\\cdot$ dE/dx)"),
     128, 0, 2.2,
     128, 0, 500.0);
}


void
AliFemtoPionLambdaCutMonitor::Lambda::Fill(const AliFemtoV0* track)
{
  const bool type_is_lambda = (fLambdaType == AliFemtoAnalysisPionLambda::kLambda);

  const Float_t eta = track->EtaV0(),
                 pt = track->PtV0(),

               minv = (type_is_lambda) ? track->MassLambda() : track->MassAntiLambda(),

            pt_pion = (type_is_lambda) ? track->PtNeg() : track->PtPos(), 
          pt_proton = (type_is_lambda) ? track->PtPos() : track->PtNeg(), 

           eta_pion = (type_is_lambda) ? track->EtaNeg() : track->EtaPos(), 
         eta_proton = (type_is_lambda) ? track->EtaPos() : track->EtaNeg(), 

          dedx_pion = (type_is_lambda) ? track->DedxNeg() : track->DedxPos(), 
        dedx_proton = (type_is_lambda) ? track->DedxPos() : track->DedxNeg();
/*
       nsig_tpc_pro = (type_is_lambda) ? track->PosNSigmaTPCP() : track->NegNSigmaTPCP(),
       nsig_tof_pro = (type_is_lambda) ? track->PosNSigmaTOFP() : track->NegNSigmaTOFP();

      nsig_tpc_pion = (type_is_lambda) ? track->NegNSigmaTPCPi() : track->PosNSigmaTPCPi(),
      nsig_tof_pion = (type_is_lambda) ? track->NegNSigmaTOFPi() : track->PosNSigmaTOFPi();
*/

  _minv->Fill(minv);
  _ypt->Fill(eta, pt);
  _dedx_pt_pro->Fill(pt_proton, dedx_proton);
  _dedx_pt_pi->Fill(pt_pion, dedx_pion);
}


TList*
AliFemtoPionLambdaCutMonitor::Lambda::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(_minv);
  output->Add(_ypt);
  output->Add(_dedx_pt_pro);
  output->Add(_dedx_pt_pi);

  return olist;
}




AliFemtoPionLambdaCutMonitor::Pair::Pair(const bool passing, const TString& typestr):
  _minv(NULL)
{
  const TString title_format = TString::Format("%s %%s %s; %%s",
                                               typestr.Data(),
                                               (passing ? "(PASS)" : "(FAIL)"));
  const TString pf(passing ? "_P" : "_F");

  _minv = new TH1F(
    "Pair_Minv" + pf,
    TString::Format(title_format, "M_{inv}", "M_{inv} (GeV)"),
    288, 0.0, 8.0);
}

void
AliFemtoPionLambdaCutMonitor::Pair::Fill(const AliFemtoPair *pair)
{
  _minv->Fill(pair->MInv());
}



TList*
AliFemtoPionLambdaCutMonitor::Pair::GetOutputList()
{
  TList *olist = new TList();
  TCollection *output = olist;

  output->Add(_minv);

  return olist;
};

/*
void AliFemtoPionLambdaCutMonitor::_init()
{

  //
  // Event Histograms
  //
  _event_centrality_pass = new TH1F(
    "centrality_pass",
    "Event Centrality",
    100, 0, 100.0);
  _event_centrality_fail = (TH1F*)_event_centrality_pass->Clone("centrality_fail");

  _event_multiplicity_pass = new TH1F(
    "MultPass",
    "Event Multiplicty",
    100, 0, 10000.0);
  _event_multiplicity_fail = (TH1F*)_event_multiplicity_pass->Clone("MultFail");

  _event_vertex_z_pass = new TH1F(
    "VertexZPass",
    "Vertex Z Distribution;z (cm);dN/dz",
    128, -15.0f, 15.0f);
  _event_vertex_z_fail = (TH1F*)_event_vertex_z_pass->Clone("VertexZFail");

  _event_vertex_xy_pass = new TH2F(
    "VertexXYPass",
    "Vertex XY Distribution (Pass);x (cm);y (cm); dN/(dx $\\cdot$ dy)",
    48, -1.0f, 1.0f,
    48, -1.0f, 1.0f);
  _event_vertex_xy_fail = (TH2F*)_event_vertex_xy_pass->Clone("VertexXYFail");


  //
  // Pion Histograms
  //

  _pion_minv_pass = new TH1F(
    "MinvPass",
    "Pion M_{inv}; M_{inv} (GeV); dN/dM",
    144, 0.120, .158);
  _pion_minv_fail = (TH1F*)_pion_minv_pass->Clone("MinvFail");

  _pion_dedx_pt_pass = new TH2F(
    "dEdX_Pt_Pass",
    "dE/dx vs p_{T}; p_{T} (GeV); dE/dx; dN/(p_{T} $\\cdot$ dE/dx)",
     128, 0, 500,
     128, 0, 1.5);
  _pion_dedx_pt_fail = (TH2F*)_pion_dedx_pt_pass->Clone("dEdX_Pt_Fail");

  _pion_ypt_pass = new TH2F(
    "YPt_Pass",
    "Rapidity vs Pt; Y; p_{T} (GeV); dN/(p_{T} $\\cdot$ Y)",
     140, -1.4, 1.4,
     100, 0, 3.0);
  _pion_ypt_fail = (TH2F*)_pion_ypt_pass->Clone("YPt_Fail");


  _lambda_minv_pass = new TH1F(
    "V0_MinvPass",
    "Lambda M_{inv}",
    576, 1.07, 1.118);
  _lambda_minv_fail = (TH1F*)_lambda_minv_pass->Clone("V0_MinvFail");

  _lambda_ypt_pass = new TH2F(
    "V0_YPt_Pass",
    "Rapidity vs p_{T}; Y; p_{T} (GeV); dN/(p_{T} $\\cdot$ Y)",
    140, -1.4, 1.4,
    100, 0.0, 3.0);
  _lambda_ypt_fail = (TH2F*)_lambda_ypt_pass->Clone("V0_YPt_Fail");

  _lambda_dedx_pt_pro_pass = new TH2F(
    "dEdX_Pro_Pass",
    "dE/dx vs p_{T} (Proton Daughter); p_{T} (GeV); dE/dx",
     128, 0, 2.2,
     128, 0, 500);
  _lambda_dedx_pt_pro_fail = (TH2F*)_lambda_dedx_pt_pro_pass->Clone("dEdX_Pro_Fail");

  _lambda_dedx_pt_pi_pass = new TH2F(
    "dEdX_Pi_Pass",
    "dE/dx vs p_{T} (Pion Daughter); p_{T} (GeV); dE/dx",
     128, 0, 2.2,
     128, 0, 500.0);
  _lambda_dedx_pt_pi_fail = (TH2F*)_lambda_dedx_pt_pi_pass->Clone("dEdX_Pi_Fail");


  //
  // Pair Histograms
  //
  _pair_minv_pass = new TH1F(
    "Pair_Minv_Pass",
    "M_{inv} (Pass);M_{inv} (GeV)",
    288, 0.0, 8.0);
  _pair_minv_fail = (TH1F*)_pair_minv_pass->Clone("Pair_Minv_Fail");
}

AliFemtoPionLambdaCutMonitor::AliFemtoPionLambdaCutMonitor()
{
  // _init();
}

TList* AliFemtoPionLambdaCutMonitor::GetOutputList()
{
  TList *olist = new TList();
  TObjArray *event_objs = new TObjArray();
  event_objs->SetName("Event");
  event_objs->Add(_event_centrality_pass);
  event_objs->Add(_event_centrality_fail);
  event_objs->Add(_event_multiplicity_pass);
  event_objs->Add(_event_multiplicity_fail);
  event_objs->Add(_event_vertex_z_pass);
  event_objs->Add(_event_vertex_z_fail);
  event_objs->Add(_event_vertex_xy_pass);
  event_objs->Add(_event_vertex_xy_fail);
  olist->Add(event_objs);

  TObjArray *pion_objs = new TObjArray();
  pion_objs->SetName("Pion");
  pion_objs->Add(_pion_minv_pass);
  pion_objs->Add(_pion_minv_fail);
  pion_objs->Add(_pion_dedx_pt_pass);
  pion_objs->Add(_pion_dedx_pt_fail);
  pion_objs->Add(_pion_ypt_pass);
  pion_objs->Add(_pion_ypt_fail);
  olist->Add(pion_objs);

  TObjArray *lambda_objs = new TObjArray();
  lambda_objs->SetName("Lambda");
  lambda_objs->Add(_lambda_minv_pass);
  lambda_objs->Add(_lambda_minv_fail);
  lambda_objs->Add(_lambda_ypt_pass);
  lambda_objs->Add(_lambda_ypt_fail);
  lambda_objs->Add(_lambda_dedx_pt_pro_pass);
  lambda_objs->Add(_lambda_dedx_pt_pro_fail);
  lambda_objs->Add(_lambda_dedx_pt_pi_pass);
  lambda_objs->Add(_lambda_dedx_pt_pi_fail);
  olist->Add(lambda_objs);

  TObjArray *pair_objs = new TObjArray();
  pair_objs->SetName("Pair");
  pair_objs->Add(_pair_minv_pass);
  pair_objs->Add(_pair_minv_fail);
  olist->Add(pair_objs);
  return olist;
}

void AliFemtoPionLambdaCutMonitor::Fill(const AliFemtoEvent* ev)
{
  const Float_t centrality = ev->CentralityV0();
  const Int_t multiplicty = ev->NumberOfTracks();
  const AliFemtoThreeVector vertex = ev->PrimVertPos();

  _event_centrality->Fill(centrality);
  _event_multiplicity->Fill(multiplicty);
  _event_vertex_z_pass->Fill(vertex.z());
  _event_vertex_xy_pass->Fill(vertex.x(), vertex.y());

}

void AliFemtoPionLambdaCutMonitor::Fill(const AliFemtoTrack* aTrack)
{}

void AliFemtoPionLambdaCutMonitor::Fill(const AliFemtoV0* aV0)
{}

void AliFemtoPionLambdaCutMonitor::Fill(const AliFemtoPair* aPair)
{}
*/
