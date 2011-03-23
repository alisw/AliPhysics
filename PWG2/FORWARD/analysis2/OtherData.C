//____________________________________________________________________
/**
 * @defgroup pwg2_forward_otherdata  External data 
 *
 * @ingroup pwg2_forward_scripts
 */
/**
 * @file 
 * 
 * @ingroup pwg2_forward_script_otherdata
 */
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>

//____________________________________________________________________
/**
 * Values used 
 * 
 * @ingroup pwg2_forward_otherdata 
 */
enum { 
  UA5, 
  CMS, 
  ALICE, 
  INEL, 
  INELGt0, 
  NSD
};
enum { 
  /** Style used for UA5 data */
  UA5Style   = 21, 
  /** Style used for CMS data */
  CMSStyle   = 29, 
  /** Style used for ALICE published data */
  ALICEStyle = 27,
  /** Color used for UA5 data */
  UA5Color   = kBlue+1,
  /** Color used for CMS data */
  CMSColor   = kGreen+1,
  /** Color used for ALICE data */
  ALICEColor = kMagenta+1,
}; 
enum { 
  /** Marker style INEL data */
  INELStyle   = 22,
  /** Marker style INEL>0 data */
  INELGt0Style= 29,
  /** Marker style NSD data */
  NSDStyle    = 23,
  /** Color used for UA5 data */
  INELColor   = kBlue+1,
  /** Color used for CMS data */
  INELGt0Color = kGreen+1,
  /** Color used for ALICE data */
  NSDColor     = kMagenta+1
};
enum {
  /** Style offset for mirror data */
  MirrorOff  = 4
};

//____________________________________________________________________
/** 
 * Set graph attributes based on trigger type and experiment. 
 * 
 * @param g        Graph
 * @param trig     Trigger (INEL, INEL>0, NSD)
 * @param exp      Experiment 
 * @param mirror   True if mirrored data 
 * @param name     Name of graph 
 * @param title    Title of graph 
 * 
 * @ingroup pwg2_forward_otherdata
 */
void
SetGraphAttributes(TGraph* g, Int_t trig, Int_t exp, bool mirror,
		   const Char_t* name, const Char_t* title)
{
  Int_t color = 0;
  switch (exp) { 
  case UA5:   color = UA5Color;   break;
  case CMS:   color = CMSColor;   break;
  case ALICE: color = ALICEColor; break;
  }
  Int_t style = 0;
  switch (exp) { 
  case UA5:   style = UA5Style;   break;
  case CMS:   style = CMSStyle;   break;
  case ALICE: style = ALICEStyle; break;
  }
  Float_t size = g->GetMarkerSize();
  switch (style) {
  case 21: 
  case 25: size *= 0.8; break;
  case 27: size *= 1.4; break;
  }
    
  if (mirror) style += MirrorOff;

  g->SetName(name);
  g->SetTitle(title);
  g->SetMarkerStyle(style);
  g->SetMarkerSize(size);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetFillColor(0);
  g->SetFillStyle(0);
  g->GetHistogram()->SetStats(kFALSE);
  g->GetHistogram()->SetXTitle("#eta");
  g->GetHistogram()->SetYTitle("#frac{1}{N} #frac{dN_{ch}}{#eta}");
}

//____________________________________________________________________
/** 
 * Get the UA5 NSD data for pp at @f$ \sqrt{s} = 900GeV@f$
 * p7886_d1x1y4 - Z.Phys.C33:1-6,1986.
 *
 * @param mirrored Wether to produce the mirrored plot 
 * 
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* UA5Nsd(Bool_t mirrored=false) 
{
  //UA5 data NSD - p7886_d1x1y4 - Z.Phys.C33:1-6,1986.
  double x[] = { 0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 
		 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125, 4.375,
		 4.625 };
  double exm[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
		   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double exp[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
		   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double y[] = { 3.48, 3.38, 3.52, 3.68, 3.71, 3.86, 3.76, 3.66, 3.72, 
		 3.69, 3.56, 3.41, 3.15, 3.09, 2.74, 2.73, 2.32, 1.99, 1.69 };
  double eym[] = { 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 
		   0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 0.13 };
  double eyp[] = { 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 
		   0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 0.13 };
  const int np = 19;
  double xm[np]; 
  double exmm[np];
  double expm[np];
  double ym[np];
  double eymm[np];
  double eypm[np];
  for (int i = 0; i < np; i++) {
    int j = np-1-i;
    xm[i]   = -x[j];
    exmm[i] = exm[j];
    expm[i] = exp[j];
    ym[i]   = y[j];
    eymm[i] = eym[j];
    eypm[i] = eyp[j];
  }

  TGraphAsymmErrors* g  = new TGraphAsymmErrors(19,x, y, exm, exp, eym, eyp);
  TGraphAsymmErrors* gm = new TGraphAsymmErrors(19,xm,ym,exmm,expm,eymm,eypm);
  SetGraphAttributes(g,  NSD, UA5, false,"ua5_nsd",         "UA5 NSD");
  SetGraphAttributes(gm, NSD, UA5, true,"ua5_nsd_mirrored",
		     "UA5 NSD (mirrored)");

  return (mirrored ? gm : g);
}

//____________________________________________________________________
/** 
 * Get the UA5 INEL data for pp at @f$ \sqrt{s} = 900GeV@f$
 * p7886_d2x1y2 - Z.Phys.C33:1-6,1986.
 *
 * @param mirrored Wether to produce the mirrored plot 
 * 
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* UA5Inel(Bool_t mirrored=false) 
{
  //UA5 data INEL - p7886_d2x1y2 - Z.Phys.C33:1-6,1986.
  double x[] = { 0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875, 2.125, 
		 2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125, 4.375,
		 4.625 };
  double exm[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
		   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double exp[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
		   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  double y[] = { 3.14, 3.04, 3.17, 3.33, 3.33, 3.53, 3.46, 3.41, 3.45, 
		 3.39, 3.07, 3.07, 2.93, 2.93, 2.55, 2.48, 2.18, 1.91, 1.52 };
  double eym[] = { 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 
		   0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 0.13 };
  double eyp[] = { 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 
		   0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 0.13 };
  const int np = 19;
  double xm[np]; 
  double exmm[np];
  double expm[np];
  double ym[np];
  double eymm[np];
  double eypm[np];
  for (int i = 0; i < np; i++) {
    int j = np-1-i;
    xm[i]   = -x[j];
    exmm[i] = exm[j];
    expm[i] = exp[j];
    ym[i]   = y[j];
    eymm[i] = eym[j];
    eypm[i] = eyp[j];
  }
  TGraphAsymmErrors* g  = new TGraphAsymmErrors(np,x, y, exm, exp, eym, eyp);
  TGraphAsymmErrors* gm = new TGraphAsymmErrors(np,xm,ym,exmm,expm,eymm,eypm);

  SetGraphAttributes(g,  INEL, UA5, false, "ua5_inel", "UA5 INEL");
  SetGraphAttributes(gm, INEL, UA5, true, "ua5_inel_mirrored", 
		     "UA5 INEL (mirrored)");
  
  return (mirrored ? gm : g);
}

//____________________________________________________________________
/** 
 * Get the ALICE INEL data in @f$ |\eta|<1.3@f$ for pp at @f$ \sqrt{s}
 * = 900GeV@f$ 
 * p7742_d1x1y1 - Eur.Phys.J.C68:89-108,2010. 
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* AliceCentralInel900()
{  
  // INEL - p7742_d1x1y1 - Eur.Phys.J.C68:89-108,2010. 
  TGraphAsymmErrors* g = 0;
  double x[] = { -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, 
				 -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3 };
  double exm[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
		   0.1, 0.1, 0.1, 0.1, 0.1 };
  double exp[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
		   0.1, 0.1, 0.1, 0.1, 0.1 };
  double y[] = { 3.28, 3.28, 3.22, 3.12, 3.06, 3.02, 2.98, 3.02, 3.02, 
		 3.05, 3.15, 3.21, 3.26, 3.33 };
  double eym[] = { 0.06324555320336758, 0.06324555320336758, 
		   0.06324555320336758, 0.06324555320336758, 
		   0.06324555320336758, 0.05385164807134505, 
		   0.05385164807134505, 0.05385164807134505, 
		   0.05385164807134505, 0.06324555320336758, 
		   0.06324555320336758, 0.06324555320336758, 
		   0.06324555320336758, 0.06324555320336758 };
  double eyp[] = { 0.08246211251235322, 0.08246211251235322, 
		   0.08246211251235322, 0.08246211251235322, 
		   0.08246211251235322, 0.08246211251235322, 
		   0.07280109889280519, 0.08246211251235322, 
		   0.08246211251235322, 0.08246211251235322, 
		   0.08246211251235322, 0.08246211251235322, 
		   0.08246211251235322, 0.08246211251235322 };
  const int np = 14;
  for (int i = 0; i < np; i++) { 
    eym[i] += 0.02;
    eyp[i] += 0.02;
  }
  g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, INEL, ALICE, false, "alice_inel", 
		     "ALICE INEL (publ.)");

  return g;
}

//____________________________________________________________________
/** 
 * Get the ALICE INEL>0 data in @f$ |\eta|<1.3@f$ for pp at @f$
 * \sqrt{s} = 900GeV@f$ 
 *
 * p7741_d4x1y1 - Eur.Phys.J.C68:345-354,2010. 
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* AliceCentralInelGt900()
{  
  // INEL>0 - p7741_d4x1y1 - Eur.Phys.J.C68:345-354,2010. 
  double x[] = { -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 
    0.9 };
  double exm[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
    0.1 };
  double exp[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
    0.1 };
  double y[] = { 4.0, 3.87, 3.8, 3.7, 3.67, 3.73, 3.72, 3.77, 3.92, 
    4.01 };
  double eym[] = { 0.07615773105863909, 0.07615773105863909, 
		   0.07615773105863909, 0.06324555320336758, 
		   0.06324555320336758, 0.06324555320336758, 
		   0.0670820393249937, 0.07615773105863909, 
		   0.07615773105863909, 0.07615773105863909 };
  double eyp[] = { 0.08544003745317531, 0.07615773105863909, 
		   0.07615773105863909, 0.07280109889280519, 
		   0.07280109889280519, 0.07280109889280519, 
		   0.07615773105863909, 0.07615773105863909, 
		   0.08544003745317531, 0.08544003745317531 };
  const int np = 10;
  for (int i = 0; i < np; i++) { 
    double stat = (i >= 3  && i<=5) ? 0.02 : 0.03;
    eym[i] += stat;
    eyp[i] += stat;
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, INELGt0, ALICE, false, "alice_inelgt900", 
		     "ALICE INEL>0 (publ.)");
  return g;

}

//____________________________________________________________________
/** 
 * Get the ALICE INEL>0 data in @f$ |\eta|<0.9@f$ for pp at @f$
 * \sqrt{s} = 2.36TeV@f$ 
 *
 * p7741_d5x1y1 - Eur.Phys.J.C68:345-354,2010. 
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* AliceCentralInelGt2360()
{  
  // INEL>0 - p7741_d5x1y1 - Eur.Phys.J.C68:345-354,2010. 
  double x[] = { -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9 };
  double exm[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
  double exp[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
  double y[] = { 4.91, 4.76, 4.63, 4.64, 4.55, 4.55, 4.64, 4.66, 4.82, 4.88 };
  double eym[] = { 0.08544003745317531, 0.08544003745317531, 
		   0.08544003745317531, 0.08544003745317531, 
		   0.08544003745317531, 0.08544003745317531, 
		   0.08544003745317531, 0.08544003745317531, 
		   0.08544003745317531, 0.08544003745317531 };
  double eyp[] = { 0.11401754250991379, 0.11401754250991379, 
		   0.1044030650891055, 0.1044030650891055, 
		   0.1044030650891055, 0.1044030650891055, 
		   0.1044030650891055, 0.1044030650891055, 
		   0.11401754250991379, 0.11401754250991379 };
  const int np = 10;
  for (int i = 0; i < np; i++) { 
    double stat = 0.3;
    eym[i] += stat;
    eyp[i] += stat;
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, INELGt0, ALICE, false, "alice_inelgt2360", 
		     "ALICE INEL>0 (publ.)");
  return g;
}

//____________________________________________________________________
/** 
 * Get the ALICE INEL>0 data in @f$ |\eta|<0.9@f$ for pp at @f$
 * \sqrt{s} = 7TeV@f$ 
 *
 * p7741_d6x1y1 - Eur.Phys.J.C68:345-354,2010. 
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* AliceCentralInelGt7000()
{  
  // INEL>0 - p7741_d6x1y1 - Eur.Phys.J.C68:345-354,2010. 
// Plot: p7741_d6x1y1
  double x[] = { -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9 };
  double exm[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
  double exp[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
  double y[] = { 6.22, 6.07, 6.01, 5.84, 5.85, 5.85, 5.91, 6.01, 6.17, 6.26 };
  double eym[] = { 0.1216552506059644, 0.1216552506059644, 
		   0.1216552506059644, 0.11180339887498948, 
		   0.11180339887498948, 0.11180339887498948, 
		   0.11180339887498948, 0.1216552506059644, 
		   0.1216552506059644, 0.1216552506059644 };
  double eyp[] = { 0.21095023109728983, 0.21095023109728983, 
		   0.2009975124224178, 0.2009975124224178, 
		   0.2009975124224178, 0.2009975124224178, 
		   0.2009975124224178, 0.2009975124224178, 
		   0.21095023109728983, 0.21095023109728983 };
  const int np = 10;
  for (int i = 0; i < np; i++) { 
    double stat = 0.2;
    eym[i] += stat;
    eyp[i] += stat;
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, INELGt0, ALICE, false, "alice_inelgt7000", 
		     "ALICE INEL>0 (publ.)");
  return g;
}

//____________________________________________________________________
/** 
 * Get the ALICE NSD data in @f$ |\eta|<1.3@f$ for pp at @f$
 * \sqrt{s} = 900GeV@f$ 
 *
 * p7742_d2x1y1 - Eur.Phys.J.C68:89-108,2010. 
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* AliceCentralNsd900()
{
  //NSD - p7742_d2x1y1 - Eur.Phys.J.C68:89-108,2010. 
  double x[] = { -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 
		 0.5, 0.7, 0.9, 1.1, 1.3 };
  double exm[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
		   0.1, 0.1, 0.1, 0.1, 0.1 };
  double exp[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
		   0.1, 0.1, 0.1, 0.1, 0.1 };
  double y[] = { 3.9, 3.89, 3.81, 3.7, 3.64, 3.59, 3.53, 3.58, 3.59, 
		 3.61, 3.74, 3.8, 3.87, 3.95 };
  double eym[] = { 0.13341664064126335, 0.13152946437965907, 
		   0.13152946437965907, 0.1216552506059644, 
		   0.1216552506059644, 0.1216552506059644, 
		   0.1216552506059644, 0.1216552506059644, 
		   0.1216552506059644, 0.1216552506059644, 
		   0.1216552506059644, 0.13152946437965907, 
		   0.13152946437965907, 0.13341664064126335 };
  double eyp[] = { 0.13341664064126335, 0.13152946437965907, 
		   0.13152946437965907, 0.1216552506059644, 
		   0.1216552506059644, 0.1216552506059644, 
		   0.1216552506059644, 0.1216552506059644, 
		   0.1216552506059644, 0.1216552506059644, 
		   0.1216552506059644, 0.13152946437965907, 
		   0.13152946437965907, 0.13341664064126335 };
  const int np = 14;
  for (int i = 0; i < np; i++) { 
    double stat = (i == 0 || i == np-1) ? 0.03 : 0.02;
    eym[i] += stat;
    eyp[i] += stat;
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, NSD, ALICE, false, "alice_nsd", "ALICE NSD (publ.)");

  return g;
}

//____________________________________________________________________
/** 
 * Get the ALICE INEL data in @f$ |\eta|<1.3@f$ for pp at @f$
 * \sqrt{s} = 2.36TeV@f$ 
 *
 * p7742_d3x1y1 - Eur.Phys.J.C68:89-108,2010. 
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* AliceCentralInel2360()
{  
  // INEL - p7742_d3x1y1 - Eur.Phys.J.C68:89-108,2010. 
  double x[] = { -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 
		 0.5, 0.7, 0.9, 1.1, 1.3 };
  double exm[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
		   0.1, 0.1, 0.1, 0.1, 0.1 };
  double exp[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
		   0.1, 0.1, 0.1, 0.1, 0.1 };
  double y[] = { 4.08, 4.01, 4.0, 3.88, 3.77, 3.8, 3.73, 3.71, 3.79, 
		 3.82, 3.94, 4.02, 4.05, 4.16 };
  double eym[] = { 0.13341664064126335, 0.12369316876852982, 
		   0.12369316876852982, 0.1216552506059644, 
		   0.1216552506059644, 0.1216552506059644, 
		   0.1216552506059644, 0.11180339887498948, 
		   0.1216552506059644, 0.1216552506059644, 
		   0.12369316876852982, 0.12369316876852982, 
		   0.13341664064126335, 0.13341664064126335 };
  double eyp[] = { 0.2716615541441225, 0.2716615541441225, 
		   0.2716615541441225, 0.260768096208106, 
		   0.25079872407968906, 0.25079872407968906, 
		   0.25079872407968906, 0.25079872407968906, 
		   0.25079872407968906, 0.260768096208106, 
		   0.261725046566048, 0.2716615541441225, 
		   0.2716615541441225, 0.2816025568065745 };
  const int np = 14;
  for (int i = 0; i < np; i++) { 
    double stat = (i < 3  || i > np-1-4) ? 0.03 : 0.02;
    eym[i] += stat;
    eyp[i] += stat;
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, NSD, ALICE, false, "alice_inel2360", 
		     "ALICE INEL (publ.)");
  return g;
}

//____________________________________________________________________
/** 
 * Get the ALICE NSD data in @f$ |\eta|<1.3@f$ for pp at @f$
 * \sqrt{s} = 2.36TeV@f$ 
 *
 * p7742_d4x1y1 - Eur.Phys.J.C68:89-108,2010. 
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* AliceCentralNsd2360()
{  
  // NSD - p7742_d4x1y1 - Eur.Phys.J.C68:89-108,2010. 
  double x[] = { -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 
		 0.5, 0.7, 0.9, 1.1, 1.3 };
  double exm[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
		   0.1, 0.1, 0.1, 0.1, 0.1 };
  double exp[] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 
		   0.1, 0.1, 0.1, 0.1, 0.1 };
  double y[] = { 4.79, 4.72, 4.7, 4.56, 4.44, 4.47, 4.39, 4.37, 4.47, 
    4.5, 4.64, 4.73, 4.76, 4.9 };
  double eym[] = { 0.13601470508735444, 0.13341664064126335, 
		   0.13341664064126335, 0.12369316876852982, 
		   0.12369316876852982, 0.12369316876852982, 
		   0.12369316876852982, 0.12369316876852982, 
		   0.12369316876852982, 0.12369316876852982,
		   0.12369316876852982, 0.13341664064126335,
		   0.13341664064126335, 0.13341664064126335 };
  double eyp[] = { 0.18439088914585774, 0.18248287590894657,
		   0.18248287590894657, 0.1726267650163207,
		   0.1726267650163207, 0.1726267650163207, 
		   0.16278820596099708, 0.16278820596099708, 
		   0.1726267650163207, 0.1726267650163207, 
		   0.1726267650163207, 0.18248287590894657, 
		   0.18248287590894657, 0.18248287590894657 };
  const int np = 14;
  
  for (int i = 0; i < np; i++) { 
    double stat = (i < 1) ? 0.03 : 0.02;
    eym[i] += stat;
    eyp[i] += stat;
  }

  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, NSD, ALICE, false, "alice_nsd2360", 
		     "ALICE NSD (publ.)");
  return g;
}

  
//____________________________________________________________________
/** 
 * Get the CMS NSD data in @f$ |\eta|<2.25@f$ for pp at @f$
 * \sqrt{s} = 900GeV@f$ 
 *
 * p7743_d8x1y1
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* CMSNsd900()
{
  // CMS published NSD data - p7743_d8x1y1
  double x[] = { -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 
    2.25 };
  double exm[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };
  double exp[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };
  double y[] = { 3.6, 3.73, 3.62, 3.54, 3.48, 3.48, 3.54, 3.62, 3.73,  3.6 };
  double eym[] = { 0.13, 0.14, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.14,0.13 };
  double eyp[] = { 0.13, 0.14, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.14, 0.13 };
  const int np = 10;
  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, NSD, CMS, false, "cms_nsd900", "CMS NSD");

  return g;
}


//____________________________________________________________________
/** 
 * Get the CMS NSD data in @f$ |\eta|<2.25@f$ for pp at @f$
 * \sqrt{s} = 2.36GeV@f$ 
 *
 * p7743_d8x1y2
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* CMSNsd2360()
{
  // CMS NSD 2360 - p7743_d8x1y2
  double x[] = { -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75,1.25,1.75,2.25 };
  double exm[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };
  double exp[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };
  double y[] = { 4.78, 4.81, 4.66, 4.61, 4.47, 4.47, 4.61, 4.66, 4.81,  4.78 };
  double eym[] = { 0.17, 0.18, 0.17, 0.17, 0.16, 0.16, 0.17, 0.17, 0.18, 0.17 };
  double eyp[] = { 0.17, 0.18, 0.17, 0.17, 0.16, 0.16, 0.17, 0.17, 0.18, 0.17 };
  const int np = 10;
  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, NSD, CMS, false, "cms_nsd2360", "CMS NSD");
  return g;
}


//____________________________________________________________________
/** 
 * Get the CMS NSD data in @f$ |\eta|<2.25@f$ for pp at @f$
 * \sqrt{s} = 7TeV@f$ 
 *
 * p7838_d5x1y1
 *
 * @return graph of data 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TGraphAsymmErrors* CMSNsd7000()
{
  // CMS NSD 7000 - Plot: p7838_d5x1y1
  double x[] = { -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75,1.25,1.75,2.25 };
  double exm[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };
  double exp[] = { 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25 };
  double y[] = { 6.18, 6.26, 6.14, 6.01, 5.78, 5.78, 6.01, 6.14, 6.26,  6.18 };
  double eym[] = { 0.25, 0.25, 0.24, 0.24, 0.23, 0.23, 0.24, 0.24, 0.25, 0.25 };
  double eyp[] = { 0.25, 0.25, 0.24, 0.24, 0.23, 0.23, 0.24, 0.24, 0.25, 0.25 };
  const int np = 10;
  TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
  SetGraphAttributes(g, NSD, CMS, false, "cms_nsd7000", "CMS NSD");
  return g;
}

//____________________________________________________________________
/** 
 * Get a multi graph of data for a given energy and trigger type 
 * 
 * @param sys    Collision system (1: pp, 2: PbPb)
 * @param energy Energy in GeV (900, 2360, 7000)
 * @param type   Bit pattern of trigger type 
 *   - 0x1 INEL 
 *   - 0x2 INEL>0
 *   - 0x4 NSD 
 * @param centLow   Low centrality cut (only for PbPB)
 * @param centHigh  High centrality cut (only for PbPB)
 * @param aliceOnly Only return other ALICE data
 * 
 * @return A multi graph with the selected data. 
 * 
 * @ingroup pwg2_forward_otherdata
 */
TMultiGraph* 
GetData(UShort_t sys, 
	UShort_t energy,
	UShort_t type=0x1, 
	UShort_t centLow=0, 
	UShort_t centHigh=0, 
	bool     aliceOnly=false)
{
  TMultiGraph* mp = new TMultiGraph(Form("dndeta_%dGeV_%d_%03d_%03d", 
					 energy, type, centLow, centHigh),"");
  TString tn;
  TString en;
  TString sn;
  TString cn;
  if (sys == 1) { 
    sn = ", pp(p#bar{p})";
    if (energy < 1000) 
      en = Form(", #sqrt{s}=%dGeV", energy);
    else 
      en = Form(", #sqrt{s}=%f4.2TeV", float(energy)/1000);
    if (!(type & 0x7)) 
      Warning("GetData", "Unknown trigger mask 0x%x", type);

    if (TMath::Abs(energy-900) < 10) {
      if (type & 0x1) { 
	tn.Append(" INEL");
	if (!aliceOnly) mp->Add(UA5Inel(false));
	if (!aliceOnly) mp->Add(UA5Inel(true));
	mp->Add(AliceCentralInel900());
      }      
      if (type & 0x4) { 
	tn.Append(" NSD");
	if (!aliceOnly) mp->Add(UA5Nsd(false));
	if (!aliceOnly) mp->Add(UA5Nsd(true));
	mp->Add(AliceCentralNsd900());
	if (!aliceOnly) mp->Add(CMSNsd900());
      }
      if (type & 0x2) { 
	tn.Append(" INEL>0");
	mp->Add(AliceCentralInelGt900());
      }
    }
    else if (TMath::Abs(energy-2360) < 10) {
      if (type & 0x1) { 
	tn.Append(" INEL");
	mp->Add(AliceCentralInel2360());
      }
      if (type & 0x4) { 
	tn.Append(" NSD");
	mp->Add(AliceCentralNsd2360());
	if (!aliceOnly) mp->Add(CMSNsd2360());
      }
      if (type & 0x2) { 
	tn.Append(" INEL>0");
	mp->Add(AliceCentralInelGt2360());
      }
    }
    else if (TMath::Abs(energy-7000) < 10) {
      if (type & 0x1) { 
	tn.Append(" INEL");
      }
      if (type & 0x4) { 
	tn.Append(" NSD");
	if (!aliceOnly) mp->Add(CMSNsd7000());
      }
      if (type & 0x2) { 
	tn.Append(" INEL>0");
	mp->Add(AliceCentralInelGt7000());
      }
    }
    else 
      Warning("GetData", "No other results for sys=%d, energy=%d",
	      sys, energy);
  }
  else if (sys == 2) { 
    // Nothing for PbPb so far 
    cn = Form(", %d%%-%d%% central", centLow, centHigh);
    sn = ", PbPb";
    if (energy < 1000) 
      en = Form(", #sqrt{s_{NN}}=%dGeV", energy);
    else 
      en = Form(", #sqrt{s_{NN}}=%f4.2TeV", float(energy)/1000);
    Warning("GetData", "No other data for PbP b yet");
  }
  else 
    Warning("GetData", "Unknown system %d", sys);
  TString tit(Form("1/N dN_{ch}/d#eta%s%s%s%s", 
		   sn.Data(), en.Data(), tn.Data(), cn.Data()));
  mp->SetTitle(tit.Data());
  if (!mp->GetListOfGraphs() || mp->GetListOfGraphs()->GetEntries() <= 0) {
    delete mp;
    mp = 0;
  }
  return mp;
}

//____________________________________________________________________
/** 
 * Plot external data for a given selection of energy and trigger type
 * (see GetData)
 * 
 * @param sys    Collision system (1: pp, 2: PbPb)
 * @param energy Energy in GeV (900, 2360, 7000)
 * @param type   Bit pattern of trigger type 
 *   - 0x1 INEL 
 *   - 0x2 INEL>0
 *   - 0x4 NSD 
 * @param centLow   Low centrality cut (only for PbPB)
 * @param centHigh  High centrality cut (only for PbPB)
 * @param aliceOnly Only return other ALICE data
 * 
 * @ingroup pwg2_forward_otherdata
 */
void
OtherData(UShort_t sys=1, 
	      UShort_t energy=900, 
	      UShort_t type=0x1, 
	      UShort_t centLow=0, 
	      UShort_t centHigh=5, 
	      bool     aliceOnly=false)
{
  TMultiGraph* mp = GetData(sys, energy, type, centLow, centHigh, aliceOnly);
  if (!mp) return;

  gStyle->SetTitleX(0.1);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleW(0.85);
  gStyle->SetTitleH(0.05);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleTextColor(kWhite);
  gStyle->SetTitleFillColor(kBlack);
  gStyle->SetTitleFontSize(0.02);
  
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  
  TCanvas* c = new TCanvas("c", "dN/deta", 800, 600);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.05);
  

  mp->SetMinimum(0);
  mp->Draw("ap");
  if (mp->GetXaxis())
    mp->GetXaxis()->SetTitle("#eta");
  if (mp->GetYaxis())
    mp->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN_{ch}}{#eta}");

  TLegend* l = c->BuildLegend(0.3, 0.15, 0.7, 0.5);
  l->SetFillColor(0);
  l->SetBorderSize(0);

  c->cd();
}

//____________________________________________________________________
//
// EOF
//
