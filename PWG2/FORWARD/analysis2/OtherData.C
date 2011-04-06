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
  PYTHIA,
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
  /** Style used for Pythia data */
  PYTHIAStyle = 28,
  /** Color used for UA5 data */
  UA5Color   = kBlue+1,
  /** Color used for Pytia data */
  PYTHIAColor = kGray+2,
  /** Color used for CMS data */
  CMSColor   = kGreen+1,
  /** Color used for ALICE data */
  ALICEColor = kMagenta+1
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
  case UA5:    color = UA5Color;    break;
  case CMS:    color = CMSColor;    break;
  case ALICE:  color = ALICEColor;  break;
  case PYTHIA: color = PYTHIAColor; break;
  }
  Int_t style = 0;
  switch (exp) { 
  case UA5:    style = UA5Style;    break;
  case CMS:    style = CMSStyle;    break;
  case ALICE:  style = ALICEStyle;  break;
  case PYTHIA: style = PYTHIAStyle; break;
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
TGraphAsymmErrors*
Pythia900INEL()
{
   
  TGraphAsymmErrors *gre = new TGraphAsymmErrors(100);
  SetGraphAttributes(gre, INEL, PYTHIA, false, "pythia900Inel",
		     "Pythia INEL");
  gre->SetPoint(0,-3.95,1.78199);
  gre->SetPointError(0, 0, 0, 0.0145305, 0.0145305);
  gre->SetPoint(1,-3.85,1.85486);
  gre->SetPointError(1,0,0,0.0148246,0.0148246);
  gre->SetPoint(2,-3.75,1.93886);
  gre->SetPointError(2,0,0,0.0151566,0.0151566);
  gre->SetPoint(3,-3.65,1.96055);
  gre->SetPointError(3,0,0,0.0152411,0.0152411);
  gre->SetPoint(4,-3.55,1.98756);
  gre->SetPointError(4,0,0,0.0153458,0.0153458);
  gre->SetPoint(5,-3.45,2.02844);
  gre->SetPointError(5,0,0,0.0155028,0.0155028);
  gre->SetPoint(6,-3.35,2.09585);
  gre->SetPointError(6,0,0,0.0157583,0.0157583);
  gre->SetPoint(7,-3.25,2.13732);
  gre->SetPointError(7,0,0,0.0159134,0.0159134);
  gre->SetPoint(8,-3.15,2.1686);
  gre->SetPointError(8,0,0,0.0160295,0.0160295);
  gre->SetPoint(9,-3.05,2.25296);
  gre->SetPointError(9,0,0,0.0163383,0.0163383);
  gre->SetPoint(10,-2.95,2.29265);
  gre->SetPointError(10,0,0,0.0164815,0.0164815);
  gre->SetPoint(11,-2.85,2.34799);
  gre->SetPointError(11,0,0,0.0166792,0.0166792);
  gre->SetPoint(12,-2.75,2.35652);
  gre->SetPointError(12,0,0,0.0167095,0.0167095);
  gre->SetPoint(13,-2.65,2.40545);
  gre->SetPointError(13,0,0,0.0168821,0.0168821);
  gre->SetPoint(14,-2.55,2.43934);
  gre->SetPointError(14,0,0,0.0170006,0.0170006);
  gre->SetPoint(15,-2.45,2.45735);
  gre->SetPointError(15,0,0,0.0170633,0.0170633);
  gre->SetPoint(16,-2.35,2.48945);
  gre->SetPointError(16,0,0,0.0171744,0.0171744);
  gre->SetPoint(17,-2.25,2.51635);
  gre->SetPointError(17,0,0,0.0172669,0.0172669);
  gre->SetPoint(18,-2.15,2.55047);
  gre->SetPointError(18,0,0,0.0173836,0.0173836);
  gre->SetPoint(19,-2.05,2.58021);
  gre->SetPointError(19,0,0,0.0174846,0.0174846);
  gre->SetPoint(20,-1.95,2.58732);
  gre->SetPointError(20,0,0,0.0175087,0.0175087);
  gre->SetPoint(21,-1.85,2.60095);
  gre->SetPointError(21,0,0,0.0175547,0.0175547);
  gre->SetPoint(22,-1.75,2.59941);
  gre->SetPointError(22,0,0,0.0175495,0.0175495);
  gre->SetPoint(23,-1.65,2.63021);
  gre->SetPointError(23,0,0,0.0176532,0.0176532);
  gre->SetPoint(24,-1.55,2.61043);
  gre->SetPointError(24,0,0,0.0175867,0.0175867);
  gre->SetPoint(25,-1.45,2.61363);
  gre->SetPointError(25,0,0,0.0175975,0.0175975);
  gre->SetPoint(26,-1.35,2.60829);
  gre->SetPointError(26,0,0,0.0175795,0.0175795);
  gre->SetPoint(27,-1.25,2.61434);
  gre->SetPointError(27,0,0,0.0175999,0.0175999);
  gre->SetPoint(28,-1.15,2.61327);
  gre->SetPointError(28,0,0,0.0175963,0.0175963);
  gre->SetPoint(29,-1.05,2.57145);
  gre->SetPointError(29,0,0,0.0174549,0.0174549);
  gre->SetPoint(30,-0.95,2.55723);
  gre->SetPointError(30,0,0,0.0174066,0.0174066);
  gre->SetPoint(31,-0.85,2.57879);
  gre->SetPointError(31,0,0,0.0174798,0.0174798);
  gre->SetPoint(32,-0.75,2.516);
  gre->SetPointError(32,0,0,0.0172657,0.0172657);
  gre->SetPoint(33,-0.65,2.53709);
  gre->SetPointError(33,0,0,0.0173379,0.0173379);
  gre->SetPoint(34,-0.55,2.51197);
  gre->SetPointError(34,0,0,0.0172519,0.0172519);
  gre->SetPoint(35,-0.45,2.44052);
  gre->SetPointError(35,0,0,0.0170047,0.0170047);
  gre->SetPoint(36,-0.35,2.44882);
  gre->SetPointError(36,0,0,0.0170336,0.0170336);
  gre->SetPoint(37,-0.25,2.45308);
  gre->SetPointError(37,0,0,0.0170484,0.0170484);
  gre->SetPoint(38,-0.15,2.4622);
  gre->SetPointError(38,0,0,0.0170801,0.0170801);
  gre->SetPoint(39,-0.05,2.45735);
  gre->SetPointError(39,0,0,0.0170633,0.0170633);
  gre->SetPoint(40,0.05,2.49254);
  gre->SetPointError(40,0,0,0.017185,0.017185);
  gre->SetPoint(41,0.15,2.49479);
  gre->SetPointError(41,0,0,0.0171928,0.0171928);
  gre->SetPoint(42,0.25,2.49289);
  gre->SetPointError(42,0,0,0.0171862,0.0171862);
  gre->SetPoint(43,0.35,2.4628);
  gre->SetPointError(43,0,0,0.0170822,0.0170822);
  gre->SetPoint(44,0.45,2.51422);
  gre->SetPointError(44,0,0,0.0172596,0.0172596);
  gre->SetPoint(45,0.55,2.51268);
  gre->SetPointError(45,0,0,0.0172543,0.0172543);
  gre->SetPoint(46,0.65,2.51066);
  gre->SetPointError(46,0,0,0.0172474,0.0172474);
  gre->SetPoint(47,0.75,2.53661);
  gre->SetPointError(47,0,0,0.0173363,0.0173363);
  gre->SetPoint(48,0.85,2.54479);
  gre->SetPointError(48,0,0,0.0173642,0.0173642);
  gre->SetPoint(49,0.95,2.55391);
  gre->SetPointError(49,0,0,0.0173953,0.0173953);
  gre->SetPoint(50,1.05,2.5872);
  gre->SetPointError(50,0,0,0.0175083,0.0175083);
  gre->SetPoint(51,1.15,2.60344);
  gre->SetPointError(51,0,0,0.0175631,0.0175631);
  gre->SetPoint(52,1.25,2.60616);
  gre->SetPointError(52,0,0,0.0175723,0.0175723);
  gre->SetPoint(53,1.35,2.62156);
  gre->SetPointError(53,0,0,0.0176242,0.0176242);
  gre->SetPoint(54,1.45,2.61173);
  gre->SetPointError(54,0,0,0.0175911,0.0175911);
  gre->SetPoint(55,1.55,2.60415);
  gre->SetPointError(55,0,0,0.0175655,0.0175655);
  gre->SetPoint(56,1.65,2.60723);
  gre->SetPointError(56,0,0,0.0175759,0.0175759);
  gre->SetPoint(57,1.75,2.60427);
  gre->SetPointError(57,0,0,0.0175659,0.0175659);
  gre->SetPoint(58,1.85,2.56765);
  gre->SetPointError(58,0,0,0.017442,0.017442);
  gre->SetPoint(59,1.95,2.58602);
  gre->SetPointError(59,0,0,0.0175043,0.0175043);
  gre->SetPoint(60,2.05,2.55936);
  gre->SetPointError(60,0,0,0.0174138,0.0174138);
  gre->SetPoint(61,2.15,2.54858);
  gre->SetPointError(61,0,0,0.0173771,0.0173771);
  gre->SetPoint(62,2.25,2.5205);
  gre->SetPointError(62,0,0,0.0172811,0.0172811);
  gre->SetPoint(63,2.35,2.49491);
  gre->SetPointError(63,0,0,0.0171932,0.0171932);
  gre->SetPoint(64,2.45,2.42773);
  gre->SetPointError(64,0,0,0.0169601,0.0169601);
  gre->SetPoint(65,2.55,2.42879);
  gre->SetPointError(65,0,0,0.0169638,0.0169638);
  gre->SetPoint(66,2.65,2.39372);
  gre->SetPointError(66,0,0,0.0168409,0.0168409);
  gre->SetPoint(67,2.75,2.38412);
  gre->SetPointError(67,0,0,0.0168071,0.0168071);
  gre->SetPoint(68,2.85,2.31896);
  gre->SetPointError(68,0,0,0.0165758,0.0165758);
  gre->SetPoint(69,2.95,2.26209);
  gre->SetPointError(69,0,0,0.0163713,0.0163713);
  gre->SetPoint(70,3.05,2.24313);
  gre->SetPointError(70,0,0,0.0163026,0.0163026);
  gre->SetPoint(71,3.15,2.20403);
  gre->SetPointError(71,0,0,0.0161599,0.0161599);
  gre->SetPoint(72,3.25,2.12855);
  gre->SetPointError(72,0,0,0.0158808,0.0158808);
  gre->SetPoint(73,3.35,2.13104);
  gre->SetPointError(73,0,0,0.01589,0.01589);
  gre->SetPoint(74,3.45,2.06339);
  gre->SetPointError(74,0,0,0.0156358,0.0156358);
  gre->SetPoint(75,3.55,1.9846);
  gre->SetPointError(75,0,0,0.0153343,0.0153343);
  gre->SetPoint(76,3.65,1.95391);
  gre->SetPointError(76,0,0,0.0152153,0.0152153);
  gre->SetPoint(77,3.75,1.87998);
  gre->SetPointError(77,0,0,0.0149247,0.0149247);
  gre->SetPoint(78,3.85,1.86256);
  gre->SetPointError(78,0,0,0.0148554,0.0148554);
  gre->SetPoint(79,3.95,1.77239);
  gre->SetPointError(79,0,0,0.0144913,0.0144913);
  gre->SetPoint(80,4.05,1.72855);
  gre->SetPointError(80,0,0,0.014311,0.014311);
  gre->SetPoint(81,4.15,1.69479);
  gre->SetPointError(81,0,0,0.0141705,0.0141705);
  gre->SetPoint(82,4.25,1.64147);
  gre->SetPointError(82,0,0,0.0139459,0.0139459);
  gre->SetPoint(83,4.35,1.58116);
  gre->SetPointError(83,0,0,0.0136873,0.0136873);
  gre->SetPoint(84,4.45,1.55735);
  gre->SetPointError(84,0,0,0.0135838,0.0135838);
  gre->SetPoint(85,4.55,1.48815);
  gre->SetPointError(85,0,0,0.0132786,0.0132786);
  gre->SetPoint(86,4.65,1.40853);
  gre->SetPointError(86,0,0,0.0129185,0.0129185);
  gre->SetPoint(87,4.75,1.36979);
  gre->SetPointError(87,0,0,0.0127396,0.0127396);
  gre->SetPoint(88,4.85,1.32666);
  gre->SetPointError(88,0,0,0.0125374,0.0125374);
  gre->SetPoint(89,4.95,1.29763);
  gre->SetPointError(89,0,0,0.0123995,0.0123995);
  gre->SetPoint(90,5.05,1.25533);
  gre->SetPointError(90,0,0,0.0121957,0.0121957);
  gre->SetPoint(91,5.15,1.20912);
  gre->SetPointError(91,0,0,0.0119692,0.0119692);
  gre->SetPoint(92,5.25,1.18839);
  gre->SetPointError(92,0,0,0.0118661,0.0118661);
  gre->SetPoint(93,5.35,1.15948);
  gre->SetPointError(93,0,0,0.0117209,0.0117209);
  gre->SetPoint(94,5.45,1.1141);
  gre->SetPointError(94,0,0,0.0114892,0.0114892);
  gre->SetPoint(95,5.55,1.06315);
  gre->SetPointError(95,0,0,0.0112235,0.0112235);
  gre->SetPoint(96,5.65,1.05213);
  gre->SetPointError(96,0,0,0.0111651,0.0111651);
  gre->SetPoint(97,5.75,1.02476);
  gre->SetPointError(97,0,0,0.011019,0.011019);
  gre->SetPoint(98,5.85,0.984834);
  gre->SetPointError(98,0,0,0.0108022,0.0108022);
  gre->SetPoint(99,5.95,0.952844);
  gre->SetPointError(99,0,0,0.0106253,0.0106253);

  return gre;
}

//____________________________________________________________________
TGraphAsymmErrors*
Pythia900NSD()
{
   
  TGraphAsymmErrors *gre = new TGraphAsymmErrors(100);
  SetGraphAttributes(gre, NSD, PYTHIA, false, "pythia900NSD",
		     "Pythia NSD");

  gre->SetPoint(0,-3.95,2.11766);
  gre->SetPointError(0,0,0,0.0179417,0.0179417);
  gre->SetPoint(1,-3.85,2.20415);
  gre->SetPointError(1,0,0,0.0183045,0.0183045);
  gre->SetPoint(2,-3.75,2.30949);
  gre->SetPointError(2,0,0,0.0187368,0.0187368);
  gre->SetPoint(3,-3.65,2.34582);
  gre->SetPointError(3,0,0,0.0188836,0.0188836);
  gre->SetPoint(4,-3.55,2.38322);
  gre->SetPointError(4,0,0,0.0190335,0.0190335);
  gre->SetPoint(5,-3.45,2.43353);
  gre->SetPointError(5,0,0,0.0192334,0.0192334);
  gre->SetPoint(6,-3.35,2.51106);
  gre->SetPointError(6,0,0,0.0195373,0.0195373);
  gre->SetPoint(7,-3.25,2.56578);
  gre->SetPointError(7,0,0,0.0197491,0.0197491);
  gre->SetPoint(8,-3.15,2.60515);
  gre->SetPointError(8,0,0,0.0199,0.0199);
  gre->SetPoint(9,-3.05,2.7105);
  gre->SetPointError(9,0,0,0.0202984,0.0202984);
  gre->SetPoint(10,-2.95,2.77008);
  gre->SetPointError(10,0,0,0.0205203,0.0205203);
  gre->SetPoint(11,-2.85,2.83332);
  gre->SetPointError(11,0,0,0.0207532,0.0207532);
  gre->SetPoint(12,-2.75,2.84715);
  gre->SetPointError(12,0,0,0.0208038,0.0208038);
  gre->SetPoint(13,-2.65,2.91693);
  gre->SetPointError(13,0,0,0.0210571,0.0210571);
  gre->SetPoint(14,-2.55,2.95797);
  gre->SetPointError(14,0,0,0.0212048,0.0212048);
  gre->SetPoint(15,-2.45,2.97499);
  gre->SetPointError(15,0,0,0.0212657,0.0212657);
  gre->SetPoint(16,-2.35,3.01345);
  gre->SetPointError(16,0,0,0.0214027,0.0214027);
  gre->SetPoint(17,-2.25,3.04659);
  gre->SetPointError(17,0,0,0.0215201,0.0215201);
  gre->SetPoint(18,-2.15,3.09341);
  gre->SetPointError(18,0,0,0.0216848,0.0216848);
  gre->SetPoint(19,-2.05,3.13187);
  gre->SetPointError(19,0,0,0.0218192,0.0218192);
  gre->SetPoint(20,-1.95,3.13917);
  gre->SetPointError(20,0,0,0.0218446,0.0218446);
  gre->SetPoint(21,-1.85,3.16911);
  gre->SetPointError(21,0,0,0.0219485,0.0219485);
  gre->SetPoint(22,-1.75,3.15665);
  gre->SetPointError(22,0,0,0.0219053,0.0219053);
  gre->SetPoint(23,-1.65,3.19693);
  gre->SetPointError(23,0,0,0.0220446,0.0220446);
  gre->SetPoint(24,-1.55,3.17002);
  gre->SetPointError(24,0,0,0.0219517,0.0219517);
  gre->SetPoint(25,-1.45,3.18538);
  gre->SetPointError(25,0,0,0.0220048,0.0220048);
  gre->SetPoint(26,-1.35,3.18066);
  gre->SetPointError(26,0,0,0.0219885,0.0219885);
  gre->SetPoint(27,-1.25,3.19754);
  gre->SetPointError(27,0,0,0.0220467,0.0220467);
  gre->SetPoint(28,-1.15,3.18021);
  gre->SetPointError(28,0,0,0.0219869,0.0219869);
  gre->SetPoint(29,-1.05,3.13111);
  gre->SetPointError(29,0,0,0.0218165,0.0218165);
  gre->SetPoint(30,-0.95,3.12153);
  gre->SetPointError(30,0,0,0.0217831,0.0217831);
  gre->SetPoint(31,-0.85,3.14798);
  gre->SetPointError(31,0,0,0.0218752,0.0218752);
  gre->SetPoint(32,-0.75,3.07912);
  gre->SetPointError(32,0,0,0.0216347,0.0216347);
  gre->SetPoint(33,-0.65,3.10207);
  gre->SetPointError(33,0,0,0.0217151,0.0217151);
  gre->SetPoint(34,-0.55,3.06346);
  gre->SetPointError(34,0,0,0.0215796,0.0215796);
  gre->SetPoint(35,-0.45,2.97651);
  gre->SetPointError(35,0,0,0.0212711,0.0212711);
  gre->SetPoint(36,-0.35,2.98715);
  gre->SetPointError(36,0,0,0.0213091,0.0213091);
  gre->SetPoint(37,-0.25,2.98548);
  gre->SetPointError(37,0,0,0.0213032,0.0213032);
  gre->SetPoint(38,-0.15,3.00555);
  gre->SetPointError(38,0,0,0.0213746,0.0213746);
  gre->SetPoint(39,-0.05,3.01193);
  gre->SetPointError(39,0,0,0.0213973,0.0213973);
  gre->SetPoint(40,0.05,3.04385);
  gre->SetPointError(40,0,0,0.0215104,0.0215104);
  gre->SetPoint(41,0.15,3.04933);
  gre->SetPointError(41,0,0,0.0215297,0.0215297);
  gre->SetPoint(42,0.25,3.04659);
  gre->SetPointError(42,0,0,0.0215201,0.0215201);
  gre->SetPoint(43,0.35,3.00813);
  gre->SetPointError(43,0,0,0.0213838,0.0213838);
  gre->SetPoint(44,0.45,3.06666);
  gre->SetPointError(44,0,0,0.0215908,0.0215908);
  gre->SetPoint(45,0.55,3.07167);
  gre->SetPointError(45,0,0,0.0216085,0.0216085);
  gre->SetPoint(46,0.65,3.0659);
  gre->SetPointError(46,0,0,0.0215881,0.0215881);
  gre->SetPoint(47,0.75,3.09159);
  gre->SetPointError(47,0,0,0.0216784,0.0216784);
  gre->SetPoint(48,0.85,3.10846);
  gre->SetPointError(48,0,0,0.0217375,0.0217375);
  gre->SetPoint(49,0.95,3.11925);
  gre->SetPointError(49,0,0,0.0217752,0.0217752);
  gre->SetPoint(50,1.05,3.15558);
  gre->SetPointError(50,0,0,0.0219016,0.0219016);
  gre->SetPoint(51,1.15,3.16911);
  gre->SetPointError(51,0,0,0.0219485,0.0219485);
  gre->SetPoint(52,1.25,3.17246);
  gre->SetPointError(52,0,0,0.0219601,0.0219601);
  gre->SetPoint(53,1.35,3.19146);
  gre->SetPointError(53,0,0,0.0220258,0.0220258);
  gre->SetPoint(54,1.45,3.17458);
  gre->SetPointError(54,0,0,0.0219675,0.0219675);
  gre->SetPoint(55,1.55,3.16866);
  gre->SetPointError(55,0,0,0.0219469,0.0219469);
  gre->SetPoint(56,1.65,3.16592);
  gre->SetPointError(56,0,0,0.0219375,0.0219375);
  gre->SetPoint(57,1.75,3.16394);
  gre->SetPointError(57,0,0,0.0219306,0.0219306);
  gre->SetPoint(58,1.85,3.11956);
  gre->SetPointError(58,0,0,0.0217762,0.0217762);
  gre->SetPoint(59,1.95,3.14646);
  gre->SetPointError(59,0,0,0.02187,0.02187);
  gre->SetPoint(60,2.05,3.10147);
  gre->SetPointError(60,0,0,0.021713,0.021713);
  gre->SetPoint(61,2.15,3.09356);
  gre->SetPointError(61,0,0,0.0216853,0.0216853);
  gre->SetPoint(62,2.25,3.05328);
  gre->SetPointError(62,0,0,0.0215437,0.0215437);
  gre->SetPoint(63,2.35,3.01953);
  gre->SetPointError(63,0,0,0.0214243,0.0214243);
  gre->SetPoint(64,2.45,2.9373);
  gre->SetPointError(64,0,0,0.0211305,0.0211305);
  gre->SetPoint(65,2.55,2.92772);
  gre->SetPointError(65,0,0,0.0210961,0.0210961);
  gre->SetPoint(66,2.65,2.89154);
  gre->SetPointError(66,0,0,0.0209653,0.0209653);
  gre->SetPoint(67,2.75,2.87619);
  gre->SetPointError(67,0,0,0.0209096,0.0209096);
  gre->SetPoint(68,2.85,2.78924);
  gre->SetPointError(68,0,0,0.0205911,0.0205911);
  gre->SetPoint(69,2.95,2.72159);
  gre->SetPointError(69,0,0,0.0203399,0.0203399);
  gre->SetPoint(70,3.05,2.69089);
  gre->SetPointError(70,0,0,0.0202248,0.0202248);
  gre->SetPoint(71,3.15,2.64939);
  gre->SetPointError(71,0,0,0.0200682,0.0200682);
  gre->SetPoint(72,3.25,2.55545);
  gre->SetPointError(72,0,0,0.0197092,0.0197092);
  gre->SetPoint(73,3.35,2.56745);
  gre->SetPointError(73,0,0,0.0197555,0.0197555);
  gre->SetPoint(74,3.45,2.47503);
  gre->SetPointError(74,0,0,0.0193967,0.0193967);
  gre->SetPoint(75,3.55,2.36741);
  gre->SetPointError(75,0,0,0.0189703,0.0189703);
  gre->SetPoint(76,3.65,2.33412);
  gre->SetPointError(76,0,0,0.0188364,0.0188364);
  gre->SetPoint(77,3.75,2.2385);
  gre->SetPointError(77,0,0,0.0184466,0.0184466);
  gre->SetPoint(78,3.85,2.21768);
  gre->SetPointError(78,0,0,0.0183606,0.0183606);
  gre->SetPoint(79,3.95,2.1055);
  gre->SetPointError(79,0,0,0.0178901,0.0178901);
  gre->SetPoint(80,4.05,2.05047);
  gre->SetPointError(80,0,0,0.0176548,0.0176548);
  gre->SetPoint(81,4.15,2.00486);
  gre->SetPointError(81,0,0,0.0174574,0.0174574);
  gre->SetPoint(82,4.25,1.94573);
  gre->SetPointError(82,0,0,0.017198,0.017198);
  gre->SetPoint(83,4.35,1.87064);
  gre->SetPointError(83,0,0,0.0168629,0.0168629);
  gre->SetPoint(84,4.45,1.83735);
  gre->SetPointError(84,0,0,0.0167122,0.0167122);
  gre->SetPoint(85,4.55,1.75314);
  gre->SetPointError(85,0,0,0.0163247,0.0163247);
  gre->SetPoint(86,4.65,1.65828);
  gre->SetPointError(86,0,0,0.0158769,0.0158769);
  gre->SetPoint(87,4.75,1.60751);
  gre->SetPointError(87,0,0,0.015632,0.015632);
  gre->SetPoint(88,4.85,1.56312);
  gre->SetPointError(88,0,0,0.0154146,0.0154146);
  gre->SetPoint(89,4.95,1.52117);
  gre->SetPointError(89,0,0,0.0152064,0.0152064);
  gre->SetPoint(90,5.05,1.46553);
  gre->SetPointError(90,0,0,0.0149257,0.0149257);
  gre->SetPoint(91,5.15,1.42038);
  gre->SetPointError(91,0,0,0.014694,0.014694);
  gre->SetPoint(92,5.25,1.38816);
  gre->SetPointError(92,0,0,0.0145263,0.0145263);
  gre->SetPoint(93,5.35,1.35046);
  gre->SetPointError(93,0,0,0.0143277,0.0143277);
  gre->SetPoint(94,5.45,1.30075);
  gre->SetPointError(94,0,0,0.0140616,0.0140616);
  gre->SetPoint(95,5.55,1.24025);
  gre->SetPointError(95,0,0,0.0137307,0.0137307);
  gre->SetPoint(96,5.65,1.21806);
  gre->SetPointError(96,0,0,0.0136073,0.0136073);
  gre->SetPoint(97,5.75,1.19435);
  gre->SetPointError(97,0,0,0.0134742,0.0134742);
  gre->SetPoint(98,5.85,1.14175);
  gre->SetPointError(98,0,0,0.0131741,0.0131741);
  gre->SetPoint(99,5.95,1.09235);
  gre->SetPointError(99,0,0,0.012886,0.012886);

  return gre;
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
	if (!aliceOnly) mp->Add(Pythia900INEL());
	if (!aliceOnly) mp->Add(UA5Inel(false));
	if (!aliceOnly) mp->Add(UA5Inel(true));
	mp->Add(AliceCentralInel900());
      }      
      if (type & 0x4) { 
	tn.Append(" NSD");
	if (!aliceOnly) mp->Add(Pythia900NSD());
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
#if 0
    else 
      Warning("GetData", "No other results for sys=%d, energy=%d",
	      sys, energy);
#endif
  }
  else if (sys == 2) { 
    // Nothing for PbPb so far 
    cn = Form(", %d%%-%d%% central", centLow, centHigh);
    sn = ", PbPb";
    if (energy < 1000) 
      en = Form(", #sqrt{s_{NN}}=%dGeV", energy);
    else 
      en = Form(", #sqrt{s_{NN}}=%f4.2TeV", float(energy)/1000);
    // Warning("GetData", "No other data for PbPb yet");
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
