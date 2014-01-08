//____________________________________________________________________
/**
 * @defgroup pwglf_forward_otherdata  External data 
 *
 * Collection of external data points for comparisons and the like 
 *
 * @ingroup pwglf_forward_scripts
 */
/**
 * @file 
 * 
 * @ingroup pwglf_forward_script_otherdata
 */
#ifndef __CINT__
# include <TGraphAsymmErrors.h>
# include <TMultiGraph.h>
# include <TStyle.h>
# include <TMath.h>
# include <TCanvas.h>
# include <TLegend.h>
# include <TLegendEntry.h>
# include <TList.h>
# include <TAxis.h>
# include <TH1F.h>
#else
class TGraphAsymmErrors;
class TMultiGraph;
class TStyle;
class TCanvas;
class TLegend;
class TLegendEntry;
class TList;
class TAxis;
class TH1F;
class TGraph;
#endif

struct RefData
{
  //____________________________________________________________________
  /**
   * Values used 
   * 
   * @ingroup pwglf_forward_otherdata 
   */
  enum { 
    UA5, 
    CMS, 
    ALICE, 
    WIP,
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
    /** Color used for ALICE work-in-progress data */
    WIPStyle = 33,
    /** Style used for Pythia data */
    PYTHIAStyle = 28,
    /** Color used for UA5 data */
    UA5Color   = kBlue+1,
    /** Color used for Pytia data */
    PYTHIAColor = kGray+2,
    /** Color used for CMS data */
    CMSColor   = kGreen+1,
    /** Color used for ALICE data */
    ALICEColor = kMagenta+1,
    /** Color used for ALICE work-in-progress data */
    WIPColor = kMagenta+3
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
   * @param exp      Experiment 
   * @param mirror   True if mirrored data 
   * @param name     Name of graph 
   * @param title    Title of graph 
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static void SetGraphAttributes(TGraph* g, 
				 Int_t /*trig*/, 
				 Int_t exp, 
				 bool mirror,
				 const Char_t* name, 
				 const Char_t* title)
  {
    Int_t color = 0;
    switch (exp) { 
    case UA5:       color = UA5Color;       break;
    case CMS:       color = CMSColor;       break;
    case ALICE:     color = ALICEColor;     break;
    case WIP:       color = WIPColor;       break;
    case PYTHIA:    color = PYTHIAColor;    break;
    }
    Int_t style = 0;
    switch (exp) { 
    case UA5:       style = UA5Style;        break;
    case CMS:       style = CMSStyle;        break;
    case ALICE:     style = ALICEStyle;      break;
    case WIP:       style = WIPStyle;        break;
    case PYTHIA:    style = PYTHIAStyle;     break;
    }
    Float_t size = g->GetMarkerSize();
    switch (style) {
    case 21: 
    case 25: size *= 0.8; break;
    case 27: size *= 1.4; break;
    case 33: size *= 1.4; break;
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
   * Get PYTHIA 900GeV INEL data 
   * 
   * 
   * @return Data graph
   */
  static TGraphAsymmErrors* Pythia900INEL()
  {
   
    TGraphAsymmErrors *gre = new TGraphAsymmErrors(100);
    SetGraphAttributes(gre, INEL, PYTHIA, false, "pythia900Inel",
		       "Pythia");
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
  /** 
   * Get PYTHIA 900GeV NSD data 
   * 
   * 
   * @return Data graph
   */
  static TGraphAsymmErrors* Pythia900NSD()
  {
   
    TGraphAsymmErrors *gre = new TGraphAsymmErrors(100);
    SetGraphAttributes(gre, NSD, PYTHIA, false, "pythia900NSD",
		       "Pythia");

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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* UA5Nsd(Bool_t mirrored=false) 
  {
    //UA5 data NSD - p7886_d1x1y4 - Z.Phys.C33:1-6,1986.
    double x[] = { 0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875,2.125,
		   2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125,4.375,
		   4.625 };
    double exm[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double exp[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double y[] = { 3.48, 3.38, 3.52, 3.68, 3.71, 3.86, 3.76, 3.66, 3.72, 
		   3.69, 3.56, 3.41, 3.15, 3.09, 2.74, 2.73, 2.32, 1.99, 1.69 };
    double eym[] = { 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 
		     0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 
		     0.13 };
    double eyp[] = { 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 
		     0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 
		     0.13 };
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
    SetGraphAttributes(g,  NSD, UA5, false,"ua5_nsd",         "UA5");
    SetGraphAttributes(gm, NSD, UA5, true,"ua5_nsd_mirrored",
		       "UA5 (mirrored)");

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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* UA5Inel(Bool_t mirrored=false) 
  {
    //UA5 data INEL - p7886_d2x1y2 - Z.Phys.C33:1-6,1986.
    double x[] = { 0.125, 0.375, 0.625, 0.875, 1.125, 1.375, 1.625, 1.875,2.125,
		   2.375, 2.625, 2.875, 3.125, 3.375, 3.625, 3.875, 4.125,4.375,
		   4.625 };
    double exm[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double exp[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
		     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double y[] = { 3.14, 3.04, 3.17, 3.33, 3.33, 3.53, 3.46, 3.41, 3.45, 
		   3.39, 3.07, 3.07, 2.93, 2.93, 2.55, 2.48, 2.18, 1.91, 1.52 };
    double eym[] = { 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 
		     0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 
		     0.13 };
    double eyp[] = { 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 
		     0.07, 0.07, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 
		     0.13 };
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

    SetGraphAttributes(g,  INEL, UA5, false, "ua5_inel", "UA5");
    SetGraphAttributes(gm, INEL, UA5, true, "ua5_inel_mirrored", 
		       "UA5 (mirrored)");
  
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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInel900()
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
		       "#it{Eur.Phys.J}.#bf{C68}:89-108,2010"
		       /* "ALICE INEL (publ.)" */);

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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInelGt900()
  {  
#if 0
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
#else
    // These are from JFGO
    TGraphAsymmErrors *g = new TGraphAsymmErrors(15);
    g->SetPoint(0,-1.5,4.12575);
    g->SetPointError(0,0.1,0.1,0.0742967,0.0814571);
    g->SetPoint(1,-1.3,3.91209);
    g->SetPointError(1,0.1,0.1,0.0697701,0.0766199);
    g->SetPoint(2,-1.1,3.98377);
    g->SetPointError(2,0.1,0.1,0.0704503,0.0774795);
    g->SetPoint(3,-0.9,4.00035);
    g->SetPointError(3,0.1,0.1,0.0702388,0.0773433);
    g->SetPoint(4,-0.7,3.87228);
    g->SetPointError(4,0.1,0.1,0.067597,0.0745103);
    g->SetPoint(5,-0.5,3.79613);
    g->SetPointError(5,0.1,0.1,0.0659771,0.0727816);
    g->SetPoint(6,-0.3,3.70489);
    g->SetPointError(6,0.1,0.1,0.0642016,0.0708603);
    g->SetPoint(7,-0.1,3.67423);
    g->SetPointError(7,0.1,0.1,0.0635759,0.0701884);
    g->SetPoint(8,0.1,3.72765);
    g->SetPointError(8,0.1,0.1,0.0645004,0.071209);
    g->SetPoint(9,0.3,3.72171);
    g->SetPointError(9,0.1,0.1,0.064493,0.071182);
    g->SetPoint(10,0.5,3.77428);
    g->SetPointError(10,0.1,0.1,0.0655974,0.0723627);
    g->SetPoint(11,0.7,3.91704);
    g->SetPointError(11,0.1,0.1,0.0683783,0.0753716);
    g->SetPoint(12,0.9,4.00674);
    g->SetPointError(12,0.1,0.1,0.0703511,0.0774669);
    g->SetPoint(13,1.1,3.97948);
    g->SetPointError(13,0.1,0.1,0.0703744,0.077396);
    g->SetPoint(14,1.3,3.99165);
    g->SetPointError(14,0.1,0.1,0.0711888,0.078178);
#endif
    SetGraphAttributes(g, INELGt0, ALICE, false, "alice_inelgt900", 
		       // "ALICE INEL>0 (publ.)"
		       "#it{Eur.Phys.J}.#bf{C68}:345-354,2010");
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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInelGt2360()
  {  
#if 0
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
#else 
    // These are from JFGO
    TGraphAsymmErrors* g = new TGraphAsymmErrors(15);
    g->SetPoint(0,-1.5,4.79047);
    g->SetPointError(0,0.1,0.1,0.0844278,0.109947);
    g->SetPoint(1,-1.3,4.91068);
    g->SetPointError(1,0.1,0.1,0.0856751,0.112038);
    g->SetPoint(2,-1.1,4.87386);
    g->SetPointError(2,0.1,0.1,0.0842846,0.110628);
    g->SetPoint(3,-0.9,4.91365);
    g->SetPointError(3,0.1,0.1,0.084339,0.111049);
    g->SetPoint(4,-0.7,4.7601);
    g->SetPointError(4,0.1,0.1,0.0812087,0.107203);
    g->SetPoint(5,-0.5,4.63355);
    g->SetPointError(5,0.1,0.1,0.078687,0.104079);
    g->SetPoint(6,-0.3,4.63885);
    g->SetPointError(6,0.1,0.1,0.0785337,0.104014);
    g->SetPoint(7,-0.1,4.55439);
    g->SetPointError(7,0.1,0.1,0.0769842,0.10203);
    g->SetPoint(8,0.1,4.55087);
    g->SetPointError(8,0.1,0.1,0.0769246,0.101951);
    g->SetPoint(9,0.3,4.64118);
    g->SetPointError(9,0.1,0.1,0.0785732,0.104066);
    g->SetPoint(10,0.5,4.66172);
    g->SetPointError(10,0.1,0.1,0.0791652,0.104711);
    g->SetPoint(11,0.7,4.81871);
    g->SetPointError(11,0.1,0.1,0.0822086,0.108523);
    g->SetPoint(12,0.9,4.88193);
    g->SetPointError(12,0.1,0.1,0.0837944,0.110332);
    g->SetPoint(13,1.1,4.89068);
    g->SetPointError(13,0.1,0.1,0.0845754,0.111009);
    g->SetPoint(14,1.3,5.05663);
    g->SetPointError(14,0.1,0.1,0.0882216,0.115368);
#endif

    SetGraphAttributes(g, INELGt0, ALICE, false, "alice_inelgt2360", 
		       // "ALICE INEL>0 (publ.)"
		       "#it{Eur.Phys.J}.#bf{C68}:345-354,2010");
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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInelGt7000()
  {  
#if 0
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
#else 
    // These are from JFGO
    TGraphAsymmErrors* g = new TGraphAsymmErrors(15);
    g->SetPoint(0,-1.5,6.28573);
    g->SetPointError(0,0.1,0.1,0.125928,0.215392);
    g->SetPoint(1,-1.3,6.25573);
    g->SetPointError(1,0.1,0.1,0.124352,0.213795);
    g->SetPoint(2,-1.1,6.28779);
    g->SetPointError(2,0.1,0.1,0.124143,0.214399);
    g->SetPoint(3,-0.9,6.21881);
    g->SetPointError(3,0.1,0.1,0.122079,0.211642);
    g->SetPoint(4,-0.7,6.0728);
    g->SetPointError(4,0.1,0.1,0.118661,0.206355);
    g->SetPoint(5,-0.5,6.011);
    g->SetPointError(5,0.1,0.1,0.117043,0.204019);
    g->SetPoint(6,-0.3,5.84071);
    g->SetPointError(6,0.1,0.1,0.11346,0.198086);
    g->SetPoint(7,-0.1,5.8532);
    g->SetPointError(7,0.1,0.1,0.113569,0.198433);
    g->SetPoint(8,0.1,5.84811);
    g->SetPointError(8,0.1,0.1,0.11347,0.198261);
    g->SetPoint(9,0.3,5.91022);
    g->SetPointError(9,0.1,0.1,0.11481,0.200444);
    g->SetPoint(10,0.5,6.00649);
    g->SetPointError(10,0.1,0.1,0.116955,0.203866);
    g->SetPoint(11,0.7,6.17115);
    g->SetPointError(11,0.1,0.1,0.120583,0.209697);
    g->SetPoint(12,0.9,6.2645);
    g->SetPointError(12,0.1,0.1,0.122976,0.213197);
    g->SetPoint(13,1.1,6.36448);
    g->SetPointError(13,0.1,0.1,0.125657,0.217014);
    g->SetPoint(14,1.3,6.39489);
    g->SetPointError(14,0.1,0.1,0.127118,0.218551);
#endif
    SetGraphAttributes(g, INELGt0, ALICE, false, "alice_inelgt7000", 
		       // "ALICE INEL>0 (publ.)"
		       "#it{Eur.Phys.J}.#bf{C68}:345-354,2010");
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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralNsd900()
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
    SetGraphAttributes(g, NSD, ALICE, false, "alice_nsd", 
		       //"ALICE NSD (publ.)"
		       "#it{Eur.Phys.J}.#bf{C68}:89-108,2010");

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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInel2360()
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
		       // "ALICE INEL (publ.)"
		       "#it{Eur.Phys.J}.#bf{C68}:89-108,2010");
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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralNsd2360()
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
		       // "ALICE NSD (publ.)"
		       "#it{Eur.Phys.J}.#bf{C68}:89-108,2010");
    return g;
  }

  //____________________________________________________________________
  /** 
   * Get the ALICE MB (CINT5, V0AND) data in @f$ |\eta|<2@f$ for pPb at 
   * @f$ \sqrt{s} = 5.023TeV@f$ 
   *
   * arXiv:1210.3615 [nucl-ex]
   *
   * @return graph of data 
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralpPb5023()
  {  
    TGraphAsymmErrors *g = new TGraphAsymmErrors(43);
    g->SetName("");
    g->SetTitle("");
    g->SetFillColor(1);
    g->SetMarkerStyle(24);
    g->SetMarkerSize(1.3);
    g->SetPoint(0,-2.05,16.1757);
    g->SetPointError(0,0.05,0.05,0.360241,0.360241);
    g->SetPoint(1,-1.95,16.8454);
    g->SetPointError(1,0.05,0.05,0.184436,0.184436);
    g->SetPoint(2,-1.85,16.9748);
    g->SetPointError(2,0.05,0.05,0.112833,0.112833);
    g->SetPoint(3,-1.75,17.2504);
    g->SetPointError(3,0.05,0.05,0.0851792,0.0851792);
    g->SetPoint(4,-1.65,17.336);
    g->SetPointError(4,0.05,0.05,0.0674453,0.0674453);
    g->SetPoint(5,-1.55,17.4649);
    g->SetPointError(5,0.05,0.05,0.0558537,0.0558537);
    g->SetPoint(6,-1.45,17.5884);
    g->SetPointError(6,0.05,0.05,0.0481065,0.0481065);
    g->SetPoint(7,-1.35,17.7509);
    g->SetPointError(7,0.05,0.05,0.0448122,0.0448122);
    g->SetPoint(8,-1.25,17.8052);
    g->SetPointError(8,0.05,0.05,0.0403298,0.0403298);
    g->SetPoint(9,-1.15,17.8357);
    g->SetPointError(9,0.05,0.05,0.0383174,0.0383174);
    g->SetPoint(10,-1.05,17.7547);
    g->SetPointError(10,0.05,0.05,0.0356689,0.0356689);
    g->SetPoint(11,-0.95,17.6859);
    g->SetPointError(11,0.05,0.05,0.034326,0.034326);
    g->SetPoint(12,-0.85,17.6665);
    g->SetPointError(12,0.05,0.05,0.0333449,0.0333449);
    g->SetPoint(13,-0.75,17.6044);
    g->SetPointError(13,0.05,0.05,0.0325356,0.0325356);
    g->SetPoint(14,-0.65,17.4815);
    g->SetPointError(14,0.05,0.05,0.0318004,0.0318004);
    g->SetPoint(15,-0.55,17.4);
    g->SetPointError(15,0.05,0.05,0.0312675,0.0312675);
    g->SetPoint(16,-0.45,17.3425);
    g->SetPointError(16,0.05,0.05,0.0310344,0.0310344);
    g->SetPoint(17,-0.35,17.2885);
    g->SetPointError(17,0.05,0.05,0.0306043,0.0306043);
    g->SetPoint(18,-0.25,17.2646);
    g->SetPointError(18,0.05,0.05,0.0303226,0.0303226);
    g->SetPoint(19,-0.15,17.316);
    g->SetPointError(19,0.05,0.05,0.0302368,0.0302368);
    g->SetPoint(20,-0.05,17.312);
    g->SetPointError(20,0.05,0.05,0.0301444,0.0301444);
    g->SetPoint(21,0.05,17.4418);
    g->SetPointError(21,0.05,0.05,0.0301526,0.0301526);
    g->SetPoint(22,0.15,17.4944);
    g->SetPointError(22,0.05,0.05,0.0303199,0.0303199);
    g->SetPoint(23,0.25,17.642);
    g->SetPointError(23,0.05,0.05,0.0303867,0.0303867);
    g->SetPoint(24,0.35,17.8153);
    g->SetPointError(24,0.05,0.05,0.0306752,0.0306752);
    g->SetPoint(25,0.45,18.0244);
    g->SetPointError(25,0.05,0.05,0.0310274,0.0310274);
    g->SetPoint(26,0.55,18.1993);
    g->SetPointError(26,0.05,0.05,0.0314353,0.0314353);
    g->SetPoint(27,0.65,18.349);
    g->SetPointError(27,0.05,0.05,0.0316803,0.0316803);
    g->SetPoint(28,0.75,18.5976);
    g->SetPointError(28,0.05,0.05,0.0322819,0.0322819);
    g->SetPoint(29,0.85,18.8045);
    g->SetPointError(29,0.05,0.05,0.0329447,0.0329447);
    g->SetPoint(30,0.95,18.9865);
    g->SetPointError(30,0.05,0.05,0.0337513,0.0337513);
    g->SetPoint(31,1.05,19.2313);
    g->SetPointError(31,0.05,0.05,0.0354009,0.0354009);
    g->SetPoint(32,1.15,19.4055);
    g->SetPointError(32,0.05,0.05,0.0367366,0.0367366);
    g->SetPoint(33,1.25,19.5893);
    g->SetPointError(33,0.05,0.05,0.0385048,0.0385048);
    g->SetPoint(34,1.35,19.8196);
    g->SetPointError(34,0.05,0.05,0.0421699,0.0421699);
    g->SetPoint(35,1.45,19.9476);
    g->SetPointError(35,0.05,0.05,0.0451541,0.0451541);
    g->SetPoint(36,1.55,20.1012);
    g->SetPointError(36,0.05,0.05,0.0513641,0.0513641);
    g->SetPoint(37,1.65,20.1082);
    g->SetPointError(37,0.05,0.05,0.060302,0.060302);
    g->SetPoint(38,1.75,20.1732);
    g->SetPointError(38,0.05,0.05,0.0739969,0.0739969);
    g->SetPoint(39,1.85,20.1964);
    g->SetPointError(39,0.05,0.05,0.0953757,0.0953757);
    g->SetPoint(40,1.95,20.0509);
    g->SetPointError(40,0.05,0.05,0.147212,0.147212);
    g->SetPoint(41,2.05,20.3151);
    g->SetPointError(41,0.05,0.05,0.272151,0.272151);
    g->SetPoint(42,2.15,20.1319);
    g->SetPointError(42,0.05,0.05,0.802706,0.802706);


    SetGraphAttributes(g, NSD, ALICE, false, "alice_ppb50230", 
		       "arXiv:1210.3615");
    return g;
  }

  //____________________________________________________________________
  /** 
   * Get the ALICE INEL data in @f$ |\eta|<1.8@f$ for pp at @f$ \sqrt{s}
   * = 900GeV@f$ 
   * Work in progress
   *
   * @return graph of data 
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInel900Work()
  {
    TGraphAsymmErrors *g = new TGraphAsymmErrors(18);
    // g->SetPoint(0,-1.9,0);   g->SetPointError(0,0.1,0.1,0,0);
    g->SetPoint(0,-1.7,3.13935);g->SetPointError(0,0.1,0.1,0.0726186,0.0525276);
    g->SetPoint(1,-1.5,3.15634);g->SetPointError(1,0.1,0.1,0.0338547,0.0380273);
    g->SetPoint(2,-1.3,3.13683);g->SetPointError(2,0.1,0.1,0.0295176,0.0295638);
    g->SetPoint(3,-1.1,3.10618);g->SetPointError(3,0.1,0.1,0.0306925,0.0329387);
    g->SetPoint(4,-0.9,3.05921);g->SetPointError(4,0.1,0.1,0.0224684,0.025408);
    g->SetPoint(5,-0.7,3.00303);g->SetPointError(5,0.1,0.1,0.0389278,0.0238328);
    g->SetPoint(6,-0.5,2.94604);g->SetPointError(6,0.1,0.1,0.0211986,0.0322219);
    g->SetPoint(7,-0.3,2.91507);g->SetPointError(7,0.1,0.1,0.030029,0.0209573);
    g->SetPoint(8,-0.1,2.88965);g->SetPointError(8,0.1,0.1,0.0286516,0.0253694);
    g->SetPoint(9,0.1,2.89731); g->SetPointError(9,0.1,0.1,0.0334615,0.0192116);
    g->SetPoint(10,0.3,2.91188);g->SetPointError(10,0.1,0.1,0.0503868,0.024911);
    g->SetPoint(11,0.5,2.96295);g->SetPointError(11,0.1,0.1,0.030009,0.0284692);
    g->SetPoint(12,0.7,3.0089);	g->SetPointError(12,0.1,0.1,0.0189095,0.026319);
    g->SetPoint(13,0.9,3.07028);g->SetPointError(13,0.1,0.1,0.0449128,0.030738);
    g->SetPoint(14,1.1,3.10215);g->SetPointError(14,0.1,0.1,0.0288688,0.026301);
    g->SetPoint(15,1.3,3.12946);g->SetPointError(15,0.1,0.1,0.0431495,0.026355);
    g->SetPoint(16,1.5,3.14549);g->SetPointError(16,0.1,0.1,0.0322482,0.033611);
    g->SetPoint(17,1.7,3.15729);g->SetPointError(17,0.1,0.1,0.105509,0.0523796);
    // g->SetPoint(19,1.9,0);	g->SetPointError(19,0.1,0.1,0,0);

    SetGraphAttributes(g, INEL, WIP, false, "alice_pp900work", 
		       "PWG-UD/MULT - work in progress");
    return g;
  }

  //____________________________________________________________________
  /** 
   * Get the ALICE NSD data in @f$ |\eta|<1.8@f$ for pp at @f$ \sqrt{s} 
   * = 900GeV@f$  
   * Work in progress
   *
   * @return graph of data 
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralNsd900Work()
  {
    TGraphAsymmErrors *g = new TGraphAsymmErrors(18);
  
    g->SetPoint(0,-1.7,3.84726);g->SetPointError(0,0.1,0.1,0.114853,0.118974);
    g->SetPoint(1,-1.5,3.87094);g->SetPointError(1,0.1,0.1,0.10574,0.108613);
    g->SetPoint(2,-1.3,3.84769);g->SetPointError(2,0.1,0.1,0.105942,0.107644);
    g->SetPoint(3,-1.1,3.8122); g->SetPointError(3,0.1,0.1,0.100838,0.101818);
    g->SetPoint(4,-0.9,3.75388);g->SetPointError(4,0.1,0.1,0.0967073,0.0972099);
    g->SetPoint(5,-0.7,3.68733);g->SetPointError(5,0.1,0.1,0.0923424,0.0925662);
    g->SetPoint(6,-0.5,3.61874);g->SetPointError(6,0.1,0.1,0.0904027,0.090477);
    g->SetPoint(7,-0.3,3.58091);g->SetPointError(7,0.1,0.1,0.0875915,0.0875992);
    g->SetPoint(8,-0.1,3.54905);g->SetPointError(8,0.1,0.1,0.086046,0.0860293);
    g->SetPoint(9,0.1,3.55968); g->SetPointError(9,0.1,0.1,0.0884776,0.0884612);
    g->SetPoint(10,0.3,3.57729);g->SetPointError(10,0.1,0.1,0.0857614,0.085769);
    g->SetPoint(11,0.5,3.63879);g->SetPointError(11,0.1,0.1,0.0879787,0.088056);
    g->SetPoint(12,0.7,3.69422);g->SetPointError(12,0.1,0.1,0.0931736,0.093396);
    g->SetPoint(13,0.9,3.76835);g->SetPointError(13,0.1,0.1,0.0928833,0.093410);
    g->SetPoint(14,1.1,3.80647);g->SetPointError(14,0.1,0.1,0.0998124,0.100799);
    g->SetPoint(15,1.3,3.83824);g->SetPointError(15,0.1,0.1,0.102549,0.104298);
    g->SetPoint(16,1.5,3.85778);g->SetPointError(16,0.1,0.1,0.10267,0.105607);
    g->SetPoint(17,1.7,3.87075);g->SetPointError(17,0.1,0.1,0.114093,0.11829);

    SetGraphAttributes(g, NSD, WIP, false, "alice_pp900nsdwork", 
		       "PWG-UD/MULT - work in progress");
    return g;
  }

  //____________________________________________________________________
  /** 
   * Get the ALICE INEL>0 data in @f$ |\eta|<1.8@f$ for pp at @f$ \sqrt{s}
   * = 900GeV@f$ 
   * Work in progress
   *
   * @return graph of data 
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInelGt900Work()
  {
    TGraphAsymmErrors *g = new TGraphAsymmErrors(10);
     
    g->SetPoint(0,-0.9,3.90755);g->SetPointError(0,0.1,0.1,0.0377085,0.0355101);
    g->SetPoint(1,-0.7,3.8357); g->SetPointError(1,0.1,0.1,0.0477674,0.0349303);
    g->SetPoint(2,-0.5,3.76291);g->SetPointError(2,0.1,0.1,0.0277709,0.040401);
    g->SetPoint(3,-0.3,3.72336);g->SetPointError(3,0.1,0.1,0.0343553,0.0250805);
    g->SetPoint(4,-0.1,3.69098);g->SetPointError(4,0.1,0.1,0.0324842,0.0324248);
    g->SetPoint(5,0.1,3.70076); g->SetPointError(5,0.1,0.1,0.0390932,0.0246738);
    g->SetPoint(6,0.3,3.71924); g->SetPointError(6,0.1,0.1,0.0576054,0.0287106);
    g->SetPoint(7,0.5,3.7844);  g->SetPointError(7,0.1,0.1,0.0316759,0.0295124);
    g->SetPoint(8,0.7,3.84319); g->SetPointError(8,0.1,0.1,0.0293134,0.0332125);
    g->SetPoint(9,0.9,3.92163); g->SetPointError(9,0.1,0.1,0.0558339,0.0394925);

    SetGraphAttributes(g, INELGt0, WIP, false, "alice_pp900inelgtwork", 
		       "PWG-UD/MULT - work in progress");
    return g;
  }

  //____________________________________________________________________
  /**
   * Get the ALICE INEL data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 2760GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInel2760Work()
  {
    TGraphAsymmErrors*    g = new TGraphAsymmErrors(18);
    g->SetPoint(0,-1.7,4.0313);  g->SetPointError(0,0.1,0.1,0.100951,0.104678);
    g->SetPoint(1,-1.5,4.0431);  g->SetPointError(1,0.1,0.1,0.10129,0.103883);
    g->SetPoint(2,-1.3,4.01251); g->SetPointError(2,0.1,0.1,0.10847,0.110089);
    g->SetPoint(3,-1.1,3.96799); g->SetPointError(3,0.1,0.1,0.105543,0.106623);
    g->SetPoint(4,-0.9,3.89669); g->SetPointError(4,0.1,0.1,0.110974,0.111625);
    g->SetPoint(5,-0.7,3.81051); g->SetPointError(5,0.1,0.1,0.108463,0.108882);
    g->SetPoint(6,-0.5,3.76537); g->SetPointError(6,0.1,0.1,0.105488,0.105773);
    g->SetPoint(7,-0.3,3.69733); g->SetPointError(7,0.1,0.1,0.110156,0.11035);
    g->SetPoint(8,-0.1,3.68148); g->SetPointError(8,0.1,0.1,0.105564,0.105733);
    g->SetPoint(9,0.1,3.67386);  g->SetPointError(9,0.1,0.1,0.1058,0.105968);
    g->SetPoint(10,0.3,3.69873); g->SetPointError(10,0.1,0.1,0.107167,0.107367);
    g->SetPoint(11,0.5,3.76377); g->SetPointError(11,0.1,0.1,0.111177,0.111448);
    g->SetPoint(12,0.7,3.81956); g->SetPointError(12,0.1,0.1,0.107198,0.107623);
    g->SetPoint(13,0.9,3.89506); g->SetPointError(13,0.1,0.1,0.105617,0.1063);
    g->SetPoint(14,1.1,3.95888); g->SetPointError(14,0.1,0.1,0.111316,0.112336);
    g->SetPoint(15,1.3,4.00176); g->SetPointError(15,0.1,0.1,0.111751,0.113315);
    g->SetPoint(16,1.5,4.03247); g->SetPointError(16,0.1,0.1,0.114383,0.116674);
    g->SetPoint(17,1.7,4.061);   g->SetPointError(17,0.1,0.1,0.107094,0.110665);
  
    SetGraphAttributes(g, INEL, WIP, false,
		       "alice_ppInel2760Work",
		       "PWG-UD/MULT - work in progress");
    return g;
  }
  //____________________________________________________________________
  /**
   * Get the ALICE NSD data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 2760GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralNsd2760Work()
  {
    TGraphAsymmErrors*    g = new TGraphAsymmErrors(18);
    g->SetPoint(0,-1.7,4.8704);  g->SetPointError(0,0.1,0.1,0.221293,0.224755);
    g->SetPoint(1,-1.5,4.88859); g->SetPointError(1,0.1,0.1,0.221269,0.223478);
    g->SetPoint(2,-1.3,4.85326); g->SetPointError(2,0.1,0.1,0.21109,0.212455);
    g->SetPoint(3,-1.1,4.80085); g->SetPointError(3,0.1,0.1,0.21041,0.211157);
    g->SetPoint(4,-0.9,4.71513); g->SetPointError(4,0.1,0.1,0.198361,0.198749);
    g->SetPoint(5,-0.7,4.61153); g->SetPointError(5,0.1,0.1,0.194009,0.194176);
    g->SetPoint(6,-0.5,4.55715); g->SetPointError(6,0.1,0.1,0.193226,0.193281);
    g->SetPoint(7,-0.3,4.47508); g->SetPointError(7,0.1,0.1,0.182433,0.182439);
    g->SetPoint(8,-0.1,4.45709); g->SetPointError(8,0.1,0.1,0.186518,0.186506);
    g->SetPoint(9,0.1,4.44707);  g->SetPointError(9,0.1,0.1,0.185747,0.185735);
    g->SetPoint(10,0.3,4.47734); g->SetPointError(10,0.1,0.1,0.185835,0.185841);
    g->SetPoint(11,0.5,4.55477); g->SetPointError(11,0.1,0.1,0.186934,0.186991);
    g->SetPoint(12,0.7,4.62236); g->SetPointError(12,0.1,0.1,0.196631,0.196796);
    g->SetPoint(13,0.9,4.71277); g->SetPointError(13,0.1,0.1,0.204034,0.20441);
    g->SetPoint(14,1.1,4.78902); g->SetPointError(14,0.1,0.1,0.20317,0.20394);
    g->SetPoint(15,1.3,4.84008); g->SetPointError(15,0.1,0.1,0.205573,0.206967);
    g->SetPoint(16,1.5,4.87453); g->SetPointError(16,0.1,0.1,0.206314,0.208667);
    g->SetPoint(17,1.7,4.90614); g->SetPointError(17,0.1,0.1,0.218996,0.222545);
  
    SetGraphAttributes(g, NSD, WIP, false,
		       "alice_ppNsd2760Work",
		       "PWG-UD/MULT - work in progress");
    return g;
  }
  //____________________________________________________________________
  /**
   * Get the ALICE INELGt0 data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 2760GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInelGt2760Work()
  {
    TGraphAsymmErrors*    g = new TGraphAsymmErrors(10);
    g->SetPoint(0,-0.9,4.96315); g->SetPointError(0,0.1,0.1,0.0439746,0.0440108);
    g->SetPoint(1,-0.7,4.8532); g->SetPointError(1,0.1,0.1,0.0426373,0.0600727);
    g->SetPoint(2,-0.5,4.79582);g->SetPointError(2,0.1,0.1,0.0475367,0.0466255);
    g->SetPoint(3,-0.3,4.70907);g->SetPointError(3,0.1,0.1,0.0313084,0.0468084);
    g->SetPoint(4,-0.1,4.68906);g->SetPointError(4,0.1,0.1,0.0413149,0.0397909);
    g->SetPoint(5,0.1,4.67937); g->SetPointError(5,0.1,0.1,0.0346151,0.0450248);
    g->SetPoint(6,0.3,4.7109);  g->SetPointError(6,0.1,0.1,0.0408403,0.0839992);
    g->SetPoint(7,0.5,4.79359); g->SetPointError(7,0.1,0.1,0.0324516,0.0357053);
    g->SetPoint(8,0.7,4.86469); g->SetPointError(8,0.1,0.1,0.0452175,0.0477304);
    g->SetPoint(9,0.9,4.96078); g->SetPointError(9,0.1,0.1,0.0566798,0.0804077);
  
    SetGraphAttributes(g, INELGt0, WIP, false,
		       "alice_ppInelGt2760Work",
		       "PWG-UD/MULT - work in progress");
    return g;
  }

  //____________________________________________________________________
  /**
   * Get the ALICE INEL data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 7000GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInel7000Work()
  {
    TGraphAsymmErrors* g = new TGraphAsymmErrors(18);
    g->SetPoint(0,-1.7,4.97541); g->SetPointError(0,0.1,0.1,0.187526,0.162049);
    g->SetPoint(1,-1.5,4.98161); g->SetPointError(1,0.1,0.1,0.128353,0.149085);
    g->SetPoint(2,-1.3,4.94853); g->SetPointError(2,0.1,0.1,0.129841,0.144762);
    g->SetPoint(3,-1.1,4.88924); g->SetPointError(3,0.1,0.1,0.137866,0.157862);
    g->SetPoint(4,-0.9,4.79998); g->SetPointError(4,0.1,0.1,0.144492,0.158783);
    g->SetPoint(5,-0.7,4.71399); g->SetPointError(5,0.1,0.1,0.132703,0.156135);
    g->SetPoint(6,-0.5,4.63098); g->SetPointError(6,0.1,0.1,0.129938,0.147085);
    g->SetPoint(7,-0.3,4.56815); g->SetPointError(7,0.1,0.1,0.129424,0.145485);
    g->SetPoint(8,-0.1,4.52372); g->SetPointError(8,0.1,0.1,0.129049,0.145285);
    g->SetPoint(9,0.1,4.52946);  g->SetPointError(9,0.1,0.1,0.131266,0.144285);
    g->SetPoint(10,0.3,4.56411); g->SetPointError(10,0.1,0.1,0.130652,0.149019);
    g->SetPoint(11,0.5,4.63554); g->SetPointError(11,0.1,0.1,0.133415,0.144298);
    g->SetPoint(12,0.7,4.71592); g->SetPointError(12,0.1,0.1,0.136436,0.151768);
    g->SetPoint(13,0.9,4.8059);  g->SetPointError(13,0.1,0.1,0.136996,0.142551);
    g->SetPoint(14,1.1,4.88457); g->SetPointError(14,0.1,0.1,0.134237,0.142764);
    g->SetPoint(15,1.3,4.92903); g->SetPointError(15,0.1,0.1,0.131933,0.152767);
    g->SetPoint(16,1.5,4.96487); g->SetPointError(16,0.1,0.1,0.140214,0.147354);
    g->SetPoint(17,1.7,4.95502); g->SetPointError(17,0.1,0.1,0.156906,0.14759);
  
    SetGraphAttributes(g, INEL, WIP, false,
		       "alice_ppInel7000Work",
		       "PWG-UD/MULT - work in progress");
    return g;
  }
  //____________________________________________________________________
  /**
   * Get the ALICE NSD data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 7000GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralNsd7000Work()
  {
    TGraphAsymmErrors* g = new TGraphAsymmErrors(18);
    g->SetPoint(0,-1.7,6.12747);g->SetPointError(0,0.1,0.1,0.157334,0.164918);
    g->SetPoint(1,-1.5,6.1353); g->SetPointError(1,0.1,0.1,0.146834,0.152011);
    g->SetPoint(2,-1.3,6.09648);g->SetPointError(2,0.1,0.1,0.140067,0.143287);
    g->SetPoint(3,-1.1,6.02552);g->SetPointError(3,0.1,0.1,0.133435,0.135281);
    g->SetPoint(4,-0.9,5.91705);g->SetPointError(4,0.1,0.1,0.129449,0.130381);
    g->SetPoint(5,-0.7,5.81246);g->SetPointError(5,0.1,0.1,0.126477,0.126883);
    g->SetPoint(6,-0.5,5.71104);g->SetPointError(6,0.1,0.1,0.124521,0.124655);
    g->SetPoint(7,-0.3,5.63422);g->SetPointError(7,0.1,0.1,0.120116,0.12013);
    g->SetPoint(8,-0.1,5.57977);g->SetPointError(8,0.1,0.1,0.119286,0.119256);
    g->SetPoint(9,0.1,5.58662); g->SetPointError(9,0.1,0.1,0.119331,0.119301);
    g->SetPoint(10,0.3,5.6291); g->SetPointError(10,0.1,0.1,0.120683,0.120697);
    g->SetPoint(11,0.5,5.7166); g->SetPointError(11,0.1,0.1,0.122787,0.122923);
    g->SetPoint(12,0.7,5.81463);g->SetPointError(12,0.1,0.1,0.126293,0.1267);
    g->SetPoint(13,0.9,5.92404);g->SetPointError(13,0.1,0.1,0.129522,0.130456);
    g->SetPoint(14,1.1,6.01958);g->SetPointError(14,0.1,0.1,0.134505,0.136333);
    g->SetPoint(15,1.3,6.07232);g->SetPointError(15,0.1,0.1,0.140728,0.143909);
    g->SetPoint(16,1.5,6.11596);g->SetPointError(16,0.1,0.1,0.14756,0.15268);
    g->SetPoint(17,1.7,6.10155);g->SetPointError(17,0.1,0.1,0.151731,0.159518);
  
    SetGraphAttributes(g, NSD, WIP, false,
		       "alice_ppNsd7000Work",
		       "PWG-UD/MULT - work in progress");
    return g;
  }

  //____________________________________________________________________
  /**
   * Get the ALICE INELGt0 data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 7000GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInelGt7000Work()
  {
    TGraphAsymmErrors* g = new TGraphAsymmErrors(10);
    g->SetPoint(0,-0.9,6.22689);g->SetPointError(0,0.1,0.1,0.0959094,0.10395);
    g->SetPoint(1,-0.7,6.11603);g->SetPointError(1,0.1,0.1,0.0609242,0.0984269);
    g->SetPoint(2,-0.5,6.00881);g->SetPointError(2,0.1,0.1,0.0595691,0.0842045);
    g->SetPoint(3,-0.3,5.9274); g->SetPointError(3,0.1,0.1,0.0560837,0.0780806);
    g->SetPoint(4,-0.1,5.86988);g->SetPointError(4,0.1,0.1,0.0552611,0.0798584);
    g->SetPoint(5,0.1,5.8773);  g->SetPointError(5,0.1,0.1,0.062512,0.077947);
    g->SetPoint(6,0.3,5.92215); g->SetPointError(6,0.1,0.1,0.0535152,0.0863595);
    g->SetPoint(7,0.5,6.01458); g->SetPointError(7,0.1,0.1,0.0578218,0.0745799);
    g->SetPoint(8,0.7,6.1186);  g->SetPointError(8,0.1,0.1,0.0767397,0.0899574);
    g->SetPoint(9,0.9,6.23468); g->SetPointError(9,0.1,0.1,0.0786932,0.073295);
  
    SetGraphAttributes(g, INELGt0, WIP, false,
		       "alice_ppInelGt7000Work",
		       "PWG-UD/MULT - work in progress");
    return g;
  }
  //____________________________________________________________________
  /**
   * Get the ALICE INEL data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 8000GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInel8000Work()
  {
    TGraphAsymmErrors*  g = new TGraphAsymmErrors(18);
    g->SetPoint(0,-1.7,5.1333); g->SetPointError(0,0.1,0.1,0.0613865,0.0707879);
    g->SetPoint(1,-1.5,5.14341);g->SetPointError(1,0.1,0.1,0.117829,0.121428);
    g->SetPoint(2,-1.3,5.13589);g->SetPointError(2,0.1,0.1,0.110807,0.113393);
    g->SetPoint(3,-1.1,5.06167);g->SetPointError(3,0.1,0.1,0.120091,0.121635);
    g->SetPoint(4,-0.9,4.97796);g->SetPointError(4,0.1,0.1,0.116141,0.117154);
    g->SetPoint(5,-0.7,4.88431);g->SetPointError(5,0.1,0.1,0.11944,0.120064);
    g->SetPoint(6,-0.5,4.81236);g->SetPointError(6,0.1,0.1,0.11049,0.110936);
    g->SetPoint(7,-0.3,4.72239);g->SetPointError(7,0.1,0.1,0.110969,0.111284);
    g->SetPoint(8,-0.1,4.66962);g->SetPointError(8,0.1,0.1,0.125108,0.125337);
    g->SetPoint(9,0.1,4.69441); g->SetPointError(9,0.1,0.1,0.113766,0.114021);
    g->SetPoint(10,0.3,4.7335); g->SetPointError(10,0.1,0.1,0.104531,0.104866);
    g->SetPoint(11,0.5,4.79917);g->SetPointError(11,0.1,0.1,0.107076,0.107534);
    g->SetPoint(12,0.7,4.88713);g->SetPointError(12,0.1,0.1,0.106124,0.106827);
    g->SetPoint(13,0.9,4.98035);g->SetPointError(13,0.1,0.1,0.120107,0.121087);
    g->SetPoint(14,1.1,5.05366);g->SetPointError(14,0.1,0.1,0.115795,0.11739);
    g->SetPoint(15,1.3,5.11276);g->SetPointError(15,0.1,0.1,0.123574,0.125877);
    g->SetPoint(16,1.5,5.16105);g->SetPointError(16,0.1,0.1,0.0979751,0.102305);
    g->SetPoint(17,1.7,5.16477);g->SetPointError(17,0.1,0.1,0.116096,0.121392);
  
    SetGraphAttributes(g, INEL, WIP, false,
		       "alice_ppInel8000Work",
		       "PWG-UD/MULT - work in progress");
    return g;
  }

  //____________________________________________________________________
  /**
   * Get the ALICE NSD data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 8000GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralNsd8000Work()
  {
    TGraphAsymmErrors* g = new TGraphAsymmErrors(18);
    g->SetPoint(0,-1.7,6.28722);g->SetPointError(0,0.1,0.1,0.220541,0.2263);
    g->SetPoint(1,-1.5,6.29211);g->SetPointError(1,0.1,0.1,0.126653,0.132922);
    g->SetPoint(2,-1.3,6.29313);g->SetPointError(2,0.1,0.1,0.142432,0.145805);
    g->SetPoint(3,-1.1,6.1944); g->SetPointError(3,0.1,0.1,0.116871,0.119093);
    g->SetPoint(4,-0.9,6.09529);g->SetPointError(4,0.1,0.1,0.121064,0.122121);
    g->SetPoint(5,-0.7,5.97811);g->SetPointError(5,0.1,0.1,0.111525,0.112012);
    g->SetPoint(6,-0.5,5.88992);g->SetPointError(6,0.1,0.1,0.119488,0.119637);
    g->SetPoint(7,-0.3,5.78296);g->SetPointError(7,0.1,0.1,0.114947,0.114962);
    g->SetPoint(8,-0.1,5.71633);g->SetPointError(8,0.1,0.1,0.0933,0.09326);
    g->SetPoint(9,0.1,5.74663); g->SetPointError(9,0.1,0.1,0.109892,0.109857);
    g->SetPoint(10,0.3,5.79472);g->SetPointError(10,0.1,0.1,0.123704,0.123718);
    g->SetPoint(11,0.5,5.87545);g->SetPointError(11,0.1,0.1,0.122522,0.122667);
    g->SetPoint(12,0.7,5.98273);g->SetPointError(12,0.1,0.1,0.128316,0.128739);
    g->SetPoint(13,0.9,6.09037);g->SetPointError(13,0.1,0.1,0.114321,0.115437);
    g->SetPoint(14,1.1,6.18105);g->SetPointError(14,0.1,0.1,0.125412,0.127476);
    g->SetPoint(15,1.3,6.24275);g->SetPointError(15,0.1,0.1,0.118631,0.122597);
    g->SetPoint(16,1.5,6.28916);g->SetPointError(16,0.1,0.1,0.144205,0.149736);
    g->SetPoint(17,1.7,6.28878);g->SetPointError(17,0.1,0.1,0.134438,0.143695);
  
    SetGraphAttributes(g, NSD, WIP, false,
		       "alice_ppNsd8000Work",
		       "PWG-UD/MULT - work in progress");
    return g;
  }
  //____________________________________________________________________
  /**
   * Get the ALICE INELGt0 data in @f$ |\eta|<1.8@f$ for pp
   * at @f$ \sqrt{s} = 8000GeV@f$
   * Work in progress
   * 
   * @return graph of data
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* AliceCentralInelGt8000Work()
  {
    TGraphAsymmErrors*  g = new TGraphAsymmErrors(10);
    g->SetPoint(0,-0.9,6.38567);g->SetPointError(0,0.1,0.1,0.0436571,0.0436571);
    g->SetPoint(1,-0.7,6.26363);g->SetPointError(1,0.1,0.1,0.0312036,0.0312036);
    g->SetPoint(2,-0.5,6.17205);g->SetPointError(2,0.1,0.1,0.0351509,0.0351509);
    g->SetPoint(3,-0.3,6.05629);g->SetPointError(3,0.1,0.1,0.0302028,0.0302028);
    g->SetPoint(4,-0.1,5.98823);g->SetPointError(4,0.1,0.1,0.0141541,0.0141541);
    g->SetPoint(5,0.1,6.02043); g->SetPointError(5,0.1,0.1,0.0256893,0.0256893);
    g->SetPoint(6,0.3,6.07111); g->SetPointError(6,0.1,0.1,0.0380304,0.0380304);
    g->SetPoint(7,0.5,6.15492); g->SetPointError(7,0.1,0.1,0.0384435,0.0384435);
    g->SetPoint(8,0.7,6.26781); g->SetPointError(8,0.1,0.1,0.0450579,0.0450579);
    g->SetPoint(9,0.9,6.38491); g->SetPointError(9,0.1,0.1,0.0396431,0.0396431);
  
    SetGraphAttributes(g, INELGt0, WIP, false,
		       "alice_ppInelGt8000Work",
		       "PWG-UD/MULT - work in progress");
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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* CMSNsd900()
  {
    // CMS published NSD data - p7743_d8x1y1
    double x[] ={ -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75,1.25,1.75,2.25};
    double exm[] ={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double exp[] ={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double y[] = { 3.6, 3.73, 3.62, 3.54, 3.48, 3.48, 3.54, 3.62, 3.73,  3.6 };
    double eym[] ={ 0.13, 0.14, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.14,0.13 };
    double eyp[] ={ 0.13, 0.14, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.14,0.13 };
    const int np = 10;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
    SetGraphAttributes(g, NSD, CMS, false, "cms_nsd900", "CMS");

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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* CMSNsd2360()
  {
    // CMS NSD 2360 - p7743_d8x1y2
    double x[] ={ -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75,1.25,1.75,2.25};
    double exm[] ={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double exp[] ={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double y[] = { 4.78, 4.81, 4.66, 4.61, 4.47, 4.47, 4.61, 4.66, 4.81,  4.78};
    double eym[] ={ 0.17, 0.18, 0.17, 0.17, 0.16, 0.16, 0.17, 0.17, 0.18, 0.17};
    double eyp[] ={ 0.17, 0.18, 0.17, 0.17, 0.16, 0.16, 0.17, 0.17, 0.18, 0.17};
    const int np = 10;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
    SetGraphAttributes(g, NSD, CMS, false, "cms_nsd2360", "CMS");
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
   * @ingroup pwglf_forward_otherdata
   */
  static TGraphAsymmErrors* CMSNsd7000()
  {
    // CMS NSD 7000 - Plot: p7838_d5x1y1
    double x[] ={ -2.25, -1.75, -1.25, -0.75, -0.25, 0.25, 0.75,1.25,1.75,2.25};
    double exm[] ={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double exp[] ={ 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
    double y[] ={ 6.18, 6.26, 6.14, 6.01, 5.78, 5.78, 6.01, 6.14, 6.26,  6.18 };
    double eym[] ={ 0.25, 0.25, 0.24, 0.24, 0.23, 0.23, 0.24, 0.24, 0.25, 0.25};
    double eyp[] ={ 0.25, 0.25, 0.24, 0.24, 0.23, 0.23, 0.24, 0.24, 0.25, 0.25};
    const int np = 10;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(np, x, y, exm, exp, eym, eyp);
    SetGraphAttributes(g, NSD, CMS, false, "cms_nsd7000", "CMS");
    return g;
  }

  //____________________________________________________________________
  static TGraphAsymmErrors* AlicePbPb2760(UShort_t centLow, UShort_t centHigh)
  {
    const Int_t  n    = 42;
    // --- /HepData/8430/d1x1y1 --------------------------------------
    const double x[]  = { -4.875, -4.625, -4.375, -4.125, -3.875, 
			  -3.625, -3.375, -3.125, -2.875, -2.625, 
			  -2.375, -2.125, -1.875, -1.625, -1.375, 
			  -1.125, -0.875, -0.625, -0.375, -0.125, 
			   +0.125, +0.375, +0.625, +0.875, +1.125,
			  +1.375, +1.625, +1.875, +2.125, +2.375, 
			  +2.625, +2.875, +3.125, +3.375, +3.625,  
			  +3.875, +4.125, +4.375, +4.625, +4.875, 
			  +5.125, +5.375 };
    const double xe[] = { 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 
			  0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 
			  0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 
			  0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 
			  0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 
			  0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 
			  0.125, 0.125, 0.125, 0.125, 0.125, 0.125 };
    const double y00_05[]  = { 888.0,  974.0, 1046.0,  1128.0, 1209.0, 1304.0, 
			       1388.0, 1449.0, 1504.0, 1578.0, 1627.0, 1674.0, 
			       1709.0, 1736.0, 1739.0, 1699.0, 1674.0, 1658.0,
			       1627.0, 1615.0, 1615.0, 1627.0, 1658.0, 1674.0,
			       1699.0, 1739.0, 1736.0, 1709.0, 1674.0, 1627.0, 
			       1578.0, 1504.0, 1449.0, 1388.0, 1304.0, 1209.0, 
			       1128.0, 1046.0, 974.0,   888.0,  770.0,  667.0};
    const double ye00_05[] = { 39.0, 41.0, 44.0, 46.0, 50.0, 54.0, 57.0, 
			       59.0, 58.0, 61.0, 62.0, 62.0, 66.0, 64.0, 
			       65.0, 68.0, 46.0, 40.0, 39.0, 39.0, 39.0, 
			       39.0, 40.0, 46.0, 68.0, 65.0, 64.0, 66.0, 
			       62.0, 62.0, 61.0, 58.0, 59.0, 57.0, 54.0, 
				50.0, 46.0, 44.0, 41.0, 39.0, 42.0, 39.0 };
    // --- /HepData/8430/d1x1y2 --------------------------------------
    const double y05_10[] = { 744.0, 807.0, 861.0, 929.0, 1001.0, 1072.0,
			      1144.0, 1196.0, 1243.0, 1298.0, 1333.0, 1361.0,
			      1393.0, 1422.0, 1418.0, 1380.0, 1362.0, 1338.0,
			      1330.0, 1318.0, 1318.0, 1330.0, 1338.0, 1362.0,
			      1380.0, 1418.0, 1422.0, 1393.0, 1361.0, 1333.0,
			      1298.0, 1243.0, 1196.0, 1144.0, 1072.0, 1001.0, 
			      929.0, 861.0, 807.0,    744.0, 651.0, 577.0 };
    const double ye05_10[]  = { 34.0, 35.0, 37.0, 39.0, 43.0, 46.0, 49.0,
				51.0, 51.0, 52.0, 54.0, 54.0, 58.0, 56.0,
				56.0, 58.0, 37.0, 32.0, 32.0, 32.0, 32.0,
				32.0, 32.0, 37.0, 58.0, 56.0, 56.0, 58.0,
				54.0, 54.0, 52.0, 51.0, 51.0, 49.0, 46.0,
				43.0, 39.0, 37.0, 35.0, 34.0, 39.0, 35.0 };
    
    // --- /HepData/8430/d1x1y3 --------------------------------------
    const double y10_20[] = { 559.0, 615.0, 656.0, 705.0, 754.0, 812.0,
			      862.0, 899.0, 928.0, 969.0, 998.0, 1022.0,
			      1048.0, 1062.0, 1064.0, 1026.0, 1014.0, 1006.0,
			      991.0, 982.0, 982.0, 991.0, 1006.0, 1014.0,
			      1026.0, 1064.0, 1062.0, 1048.0, 1022.0, 998.0,
			      969.0, 928.0, 899.0, 862.0, 812.0, 754.0,
			      705.0, 656.0, 615.0, 559.0, 491.0, 451.0 };
    const double ye10_20[]  = { 29.0, 31.0, 32.0, 33.0, 36.0, 39.0, 42.0,
				43.0, 44.0, 45.0, 45.0, 44.0, 46.0, 47.0,
				47.0, 46.0, 28.0, 24.0, 24.0, 24.0, 24.0,
				24.0, 24.0, 28.0, 46.0, 47.0, 47.0, 46.0,
				44.0, 45.0, 45.0, 44.0, 43.0, 42.0, 39.0,
				36.0, 33.0, 32.0, 31.0, 29.0, 37.0, 29.0 };

    // --- /HepData/8430/d1x1y4 --------------------------------------
    const double y20_30[] = { 389.0, 429.0, 459.0, 491.0, 526.0, 563.0,
			      596.0, 619.0, 641.0, 669.0, 683.0, 700.0,
			      716.0, 725.0, 724.0, 701.0, 692.0, 682.0,
			      673.0, 666.0, 666.0, 673.0, 682.0, 692.0,
			      701.0, 724.0, 725.0, 716.0, 700.0, 683.0,
			      669.0, 641.0, 619.0, 596.0, 563.0, 526.0,
			      491.0, 459.0, 429.0, 389.0, 343.0, 313.0 };
    const double ye20_30[]  = { 21.0, 22.0, 24.0, 25.0, 27.0, 30.0, 31.0,
				32.0, 34.0, 35.0, 33.0, 34.0, 35.0, 36.0,
				37.0, 34.0, 19.0, 16.0, 16.0, 16.0, 16.0,
				16.0, 16.0, 19.0, 34.0, 37.0, 36.0, 35.0,
				34.0, 33.0, 35.0, 34.0, 32.0, 31.0, 30.0,
				27.0, 25.0, 24.0, 22.0, 21.0, 26.0, 22.0 };

    // --- Select ----------------------------------------------------
    const double* yp  = 0;
    const double* ype = 0;
    if      (TMath::Abs(centLow- 0) < 1 && TMath::Abs(centHigh- 5) < 1) {
      yp  = y00_05;
      ype = ye00_05;
    }
    else if (TMath::Abs(centLow- 5) < 1 && TMath::Abs(centHigh-10) < 1) {
      yp  = y05_10;
      ype = ye05_10;
    }
    else if (TMath::Abs(centLow-10) < 1 && TMath::Abs(centHigh-20) < 1) {
      yp  = y10_20;
      ype = ye10_20;
    }
    else if (TMath::Abs(centLow-20) < 1 && TMath::Abs(centHigh-30) < 1) {
      yp  = y20_30;
      ype = ye20_30;
    }
    if (!yp || !ype) return 0;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(n,x,yp,xe,xe,ype,ype);
    SetGraphAttributes(g, INEL, ALICE, false, "alice_pbpb2760", 
		       "#it{Phys.Lett.} #bf{B726}:610-622,2013");
    
    return g;
  }
  //____________________________________________________________________
  /** 
   * Get a single data graph 
   * 
   * @param which    Which type
   * @param sys      Collisition system
   * @param energy   Collision energy
   * @param type     Trigger type 
   * @param centLow  Low cut on centrality 
   * @param centHigh Up cut on centraltiy
   * 
   * @return Data graph 
   */
  static TGraphAsymmErrors* GetSingle(UShort_t which, 
				      UShort_t sys, 
				      UShort_t energy, 
				      UShort_t type=0x1, 
				      UShort_t centLow=0,
				      UShort_t centHigh=0,
				      Bool_t   verbose=false) 
  {
    TGraphAsymmErrors* ret = 0;
    if (sys == 1) { 
      if (TMath::Abs(energy-900) < 10) {
	switch (type) { 
	case 1: // INEL 
	  switch (which) { 
	  case PYTHIA:    ret = Pythia900INEL(); break;
	  case UA5:       ret = UA5Inel(false);  break;
	  case UA5+10:    ret = UA5Inel(true);   break;
	  case ALICE:     ret = AliceCentralInel900(); break;
	  case WIP:       ret = AliceCentralInel900Work(); break;
	  }      
	  break;
	case 2: // INEL>0
	  switch (which) { 
	  case ALICE: ret = AliceCentralInelGt900(); break;
	  case WIP:   ret = AliceCentralInelGt900Work(); break;
	  }
	  break;
	case 4:  // NSD 
	  switch (which) { 
	  case PYTHIA: ret = Pythia900NSD(); break;
	  case UA5:    ret = UA5Nsd(false);  break;
	  case UA5+10: ret = UA5Nsd(true);   break;
	  case ALICE:  ret = AliceCentralNsd900(); break;
	  case WIP:    ret = AliceCentralNsd900Work(); break;
	  case CMS:    ret = CMSNsd900();          break;
	  }
	  break;
	} // type 
      }
      else if (TMath::Abs(energy-2360) < 10) {
	switch (type) { 
	case 1: // INEL 
	  switch (which) { 
	  case ALICE: ret = AliceCentralInel2360(); break;
	  case WIP: ret = AliceCentralInel2760Work(); break;
	  }
	  break;
	case 2: // INEL > 0
	  switch (which) {
	  case ALICE: ret = AliceCentralInelGt2360(); break;
	  case WIP: ret = AliceCentralInelGt2760Work(); break;
	  }
	  break;
	case 4: // NSD 
	  switch (which) { 
	  case ALICE: ret = AliceCentralNsd2360(); break;
	  case CMS:   ret = CMSNsd2360(); break;
	  case WIP:   ret = AliceCentralNsd2760Work(); break;
	  }
	  break;
	}
      }
      else if (TMath::Abs(energy-2760) < 10) {
	switch (type) { 
	case 1: // INEL 
	  switch (which) { 
	  case WIP: ret = AliceCentralInel2760Work(); break;
	  }
	  break;
	case 2: // INEL > 0
	  switch (which) {
	  case WIP: ret = AliceCentralInelGt2760Work(); break;
	  }
	  break;
	case 4: // NSD 
	  switch (which) { 
	  case WIP: ret = AliceCentralNsd2760Work(); break;
	  }
	  break;
	}
      }
      else if (TMath::Abs(energy-7000) < 10) {
	switch (type) { 
	case 1: 
	  switch (which) { 
	  case WIP: ret = AliceCentralInel7000Work(); break;
	  }
	  break;
	case 2: // INEL > 0
	  switch (which) { 
	  case ALICE: ret = AliceCentralInelGt7000(); break;
	  case WIP: ret = AliceCentralInelGt7000Work(); break;
	  }
	  break;
	case 4: // NSD 
	  switch (which) { 
	  case CMS: ret = CMSNsd7000(); break;
	  case WIP: ret = AliceCentralNsd7000Work(); break;
	  }
	  break;
	}
      }
      else if (TMath::Abs(energy-8000) < 10) {
	switch (type) { 
	case 1: 
	  switch (which) { 
	  case WIP: ret = AliceCentralInel8000Work(); break;
	  }
	  break;
	case 2: // INEL > 0
	  switch (which) { 
	  case WIP: ret = AliceCentralInelGt8000Work(); break;
	  }
	  break;
	case 4: // NSD 
	  switch (which) { 
	  case WIP: ret = AliceCentralNsd8000Work(); break;
	  }
	  break;
	}
      }
    } // pp 
    else if (sys == 2) { // PbPb 
      if (TMath::Abs(energy - 2760) < 10) {
	switch (which) { 
	case ALICE: ret = AlicePbPb2760(centLow, centHigh);
	}
      }
      if (ret) {
	Int_t col = CentralityColor(centLow, centHigh);
	ret->SetMarkerColor(col);
	ret->SetLineColor(col);
	ret->SetFillColor(col);
      }
    }
    else if (sys == 3) { // pPb 
      if (TMath::Abs(energy - 5023) < 10 || 
	  TMath::Abs(energy - 8000) < 10) {
	switch (which) { 
	case ALICE: ret = AliceCentralpPb5023(); break;
	}
      }
    } // pPb
    if (!ret) {
      TString w;
      switch (which) { 
      case UA5:   w = "UA5";     break;
      case CMS:   w = "CMS";     break;
      case ALICE: w = "ALICE";   break;
      case WIP:   w = "WIP";     break;
      case PYTHIA:w = "Pyhthia"; break;
      default: w = Form("unknown(%d)", which);
      }
      TString sy;
      switch (sys) { 
      case 1:  sy = "pp"; break;
      case 2:  sy = "PbPb"; break;
      case 3:  sy = "pPb"; break;
      default: sy = Form("unknown(%d)", sys);
      }
      TString tr;
      if (sys == 1 || sys == 3) {
	switch (type) { 
	case 1:        tr = "INEL"; break;
	case 2:        tr = "INEL>0"; break;
	case 4:        tr = "NSD"; break;
	default:       tr = Form("unknown(%d)", sys);
	}
      }
      else 
	tr.Form("%2d--%2d%%", centLow, centHigh);

      if (verbose) 
	Warning("GetSingle", "Nothing to get for "
		"which=%s, sys=%s, energy=%dGeV, type=%s, "
		"centLow=%d, centHigh=%d",
		w.Data(), sy.Data(), energy, tr.Data(), centLow, centHigh);
    }
#if 0
    if (ret) {
      if (!wn.IsNull()) wn.Append(",");
      switch (which) { 
      case PYTHIA: wn.Append("Pythia"); break;
      case UA5:    wn.Append("UA5");    break;
      case CMS:    wn.Append("CMS");    break;
      case ALICE:  wn.Append("ALICE");  break;
      }      
    }
#endif
    return ret;
  }

  //____________________________________________________________________
  /** 
   * Append an item to a list 
   * 
   * @param s      List to append to 
   * @param delim  Delimiter 
   * @param what   What to append 
   * 
   * @return New string value 
   */
  static TString&  AppendItem(TString& s, char delim, const char* what)	  
  {
    if (!s.IsNull()) s.Append(Form("%c", delim));
    s.Append(what);
    return s;
  }

  static void FormatEnergy(TString& en, UShort_t sys, UShort_t energy)
  {
    en.Append(Form("#sqrt{s%s}=", sys == 1 ? "" : "_{NN}"));
    if (energy < 1000) 
      en.Append(Form("%dGeV", energy));
    else {
      if (energy % 1000 == 0) 
	en.Append(Form("%dTeV", energy/1000));
      else 
	en.Append(Form("%4.2fTeV", float(energy)/1000));
    }
  }
  static void FormatSystem(TString& sn, UShort_t sys)
  {
    switch (sys) { 
    case 1:  sn = "pp";      break;
    case 2:  sn = "PbPb";    break;
    case 3:  sn = "pPb";     break;
    default: sn = "unknown"; break;
    }
  }
  static void FormatCentrality(TString& cn, UShort_t sys, UShort_t l,UShort_t u)
  {
    if (sys == 2) cn.Form("%d%% - %d%% central", l, u);
  }
  static void FormatType(TString& tn, UShort_t sys, UShort_t type)
  {
    if (sys == 2) { 
      tn = "";
      return;
    }
    if (type == 0x2000) type = 0x4;
    if (type & 0x1) AppendItem(tn, '|', "INEL");
    if (type & 0x2) AppendItem(tn, '|', "INEL>0");
    if (type & 0x4) AppendItem(tn, '|', "NSD");
  }
  static void FormatTitle(TString& title,
			  UShort_t sys, 
			  UShort_t energy,
			  UShort_t type, 
			  UShort_t centLow, 
			  UShort_t centHigh,
			  Bool_t   seenUA5)
  {
    TString sn; FormatSystem(sn, sys);
    TString en; FormatEnergy(en, sys, energy);
    TString tn; FormatType(tn, sys, type);
    TString cn; FormatCentrality(cn, sys, centLow, centHigh);
    if (seenUA5) sn.Append("(p#bar{p})");
    if (!en.IsNull()) en.Prepend(", ");
    if (!tn.IsNull()) tn.Prepend(", ");
    if (!cn.IsNull()) cn.Prepend(", ");
    title.Form("%s%s%s%s", sn.Data(), en.Data(), tn.Data(), cn.Data());
  }
  //__________________________________________________________________
  /** 
   * Get the color for a centrality bin
   * 
   * @param bin Centrality bin 
   * 
   * @return Color 
   */
  static Int_t CentralityColor(UShort_t centLow, 
			       UShort_t centHigh,
			       UShort_t /*nBins*/=0)
  {
#if 0
    if (nBins > 0 && nBins < 6) { 
      switch (bin) { 
      case 1: return kRed+2;
      case 2: return kGreen+2;
      case 3: return kBlue+1;
      case 4: return kCyan+1;
      case 5: return kMagenta+1;
      case 6: return kYellow+2;
      }
    }
#endif
    Float_t  fc       = (centLow+double(centHigh-centLow)/2) / 100;
    Int_t    nCol     = gStyle->GetNumberOfColors();
    Int_t    icol     = TMath::Min(nCol-1,int(fc * nCol + .5));
    Int_t    col      = gStyle->GetColorPalette(icol);
    //Info("GetCentralityColor","%3d: %3d-%3d -> %3d",bin,centLow,centHigh,col);
    return col;
  }
  static Bool_t GetPbPb(TMultiGraph* mp, 
			Bool_t cms, Bool_t alice, Bool_t pythia, 
			Bool_t work, UShort_t energy, UShort_t centLow, 
			UShort_t centHigh)
  {
    // Nothing for PbPb so far 
    TGraphAsymmErrors* gCMS =(cms   ?GetSingle(CMS,   2,energy,0,
					       centLow, centHigh):0);
    TGraphAsymmErrors* gALI =(alice ?GetSingle(ALICE, 2,energy,0,
					       centLow, centHigh):0);
    TGraphAsymmErrors* gPYT =(pythia?GetSingle(PYTHIA,2,energy,0,
					       centLow, centHigh):0);
    TGraphAsymmErrors* gWRK =(work  ?GetSingle(WIP,   2,energy,0,
						 centLow, centHigh):0);
    if (gCMS) mp->Add(gCMS);
    if (gALI) mp->Add(gALI);
    if (gPYT) mp->Add(gPYT);
    if (gWRK) mp->Add(gWRK);
    return (gCMS || gALI || gPYT || gWRK);
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
   * @param which     What to get
   * 
   * @return A multi graph with the selected data. 
   * 
   * @ingroup pwglf_forward_otherdata
   */
  static TMultiGraph* GetData(UShort_t sys, 
			      UShort_t energy,
			      UShort_t type=0x1, 
			      UShort_t centLow=0, 
			      UShort_t centHigh=0, 
			      UShort_t which=0x7)
  {
    TMultiGraph* mp = new TMultiGraph(Form("dndeta_%dGeV_%d_%03d_%03d", 
					   energy, type, centLow, centHigh),"");
    bool    ua5       = (which & (1 << UA5));       // 0x1
    bool    cms       = (which & (1 << CMS));       // 0x2
    bool    alice     = (which & (1 << ALICE));     // 0x4
    bool    work      = (which & (1 << WIP)); // 0x8
    bool    pythia    = (which & (1 << PYTHIA));    // 0x10
    Bool_t  seenUA5   = false;

    if (sys == 1) { 
      if (!(type & 0x7)) 
	Warning("GetData", "Unknown trigger mask 0x%x", type);

      if (TMath::Abs(energy-2750) < 11) {
	Warning("GetData", "Using 2360GeV data for %dGeV comparison", energy);
	energy = 2360;
      }
      if (!(TMath::Abs(energy-900) < 10 || 
	    TMath::Abs(energy-2360) < 10 || 
	    TMath::Abs(energy-7000) < 10 || 
	    TMath::Abs(energy-8000) < 10)) {
	Warning("GetData", "No other results for sys=%d, energy=%d",
		sys, energy);
	return 0;
      }

      // Substitute NSD for V0-AND

      for (Int_t i = 0; i < 3; i++) { 
	UShort_t mask = (1 << i);
	if ((type & mask) == 0) continue;
	TGraphAsymmErrors* gUAp =(ua5   ?GetSingle(UA5,   sys,energy,mask):0);
	TGraphAsymmErrors* gUAn =(ua5   ?GetSingle(UA5+10,sys,energy,mask):0);
	TGraphAsymmErrors* gCMS =(cms   ?GetSingle(CMS,   sys,energy,mask):0);
	TGraphAsymmErrors* gALI =(alice ?GetSingle(ALICE, sys,energy,mask):0);
	TGraphAsymmErrors* gPYT =(pythia?GetSingle(PYTHIA,sys,energy,mask):0);
	TGraphAsymmErrors* gWRK =(work  ?GetSingle(WIP,   sys,energy,mask):0);
	if (gUAp) mp->Add(gUAp);
	if (gUAn) mp->Add(gUAn);
	if (gCMS) mp->Add(gCMS);
	if (gALI) mp->Add(gALI);
	if (gPYT) mp->Add(gPYT);
	if (gWRK) mp->Add(gWRK);
	if (gUAp || gUAn) seenUA5 = true;
      }
    }
    else if (sys == 2) { 
      // Warning("GetData", "No other data for PbPb yet");
      GetPbPb(mp,cms,alice,pythia,work,energy,centLow,centHigh);
    }
    else if (sys == 3) {
      if (!(TMath::Abs(energy-5023) < 10 ||
	    TMath::Abs(energy-8000) < 10) ) {
	Warning("GetData", "No other results for sys=%d, energy=%d",
		sys, energy);
	return 0;
      }
    
      // Info("GetData", "Getting ALICE pPb data");
      TGraphAsymmErrors* gALI =(alice ?GetSingle(ALICE, sys,energy, 0):0);
      if (gALI) mp->Add(gALI);
      // Info("GetData", "Got %p", gALI);
      // Warning("GetData", "Unknown system %d", sys);
    }

    if (!mp->GetListOfGraphs() || mp->GetListOfGraphs()->GetEntries() <= 0) {
      delete mp;
      mp = 0;
      return 0;
    }
    TString title;
    FormatTitle(title, sys, energy, type, centLow, centHigh, seenUA5);
    mp->SetTitle(title.Data());
    return mp;
  }
};
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
 * @param which     Which data to show 
 * 
 * @ingroup pwglf_forward_otherdata
 */
void
OtherData(UShort_t sys=1, 
	  UShort_t energy=900, 
	  UShort_t type=0x1, 
	  UShort_t centLow=0, 
	  UShort_t centHigh=5, 
	  UShort_t which=0xf)
{
  TMultiGraph* mp = 
    RefData::GetData(sys, energy, type, centLow, centHigh, which);
  if (!mp) {
    if (sys != 2) return;

    UShort_t  cents[] = { 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100, 0 };
    UShort_t* lC      = &(cents[0]);
    UShort_t* uC      = &(cents[1]);
    TString   title;
    RefData::FormatTitle(title, sys, energy, type, centLow, centHigh, false);
    mp   = new TMultiGraph(Form("pbpb_%04dGeV_%03d_%03d", 
				energy, centLow, centHigh), title.Data());
    Bool_t filled = false;
    while (*lC < centLow) { lC++; uC++; }
    while (*uC > 0 && *uC <= centHigh) {
      TMultiGraph* mg = RefData::GetData(sys, energy, type, *lC, *uC, which);
      if (!mg) continue;
      filled = true;
      mp->Add(mg);
      uC++;
      lC++;
    }
    if (!filled) return;
  }

  gStyle->SetTitleX(0.1);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleW(0.85);
  gStyle->SetTitleH(0.05);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleTextColor(kWhite);
  gStyle->SetTitleFillColor(kBlack);
  gStyle->SetTitleFontSize(0.02);
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  TCanvas* c = new TCanvas("c", "dN/deta", 800, 600);
  c->SetFillColor(0);
  c->SetBorderSize(0);
  c->SetBorderMode(0);
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.02);
  

  mp->SetMinimum(0);
  mp->Draw("ap");
  if (mp->GetXaxis()) {
    mp->GetXaxis()->SetTitle("#eta");
    mp->GetXaxis()->SetTitleFont(12);
    mp->GetXaxis()->SetLabelFont(132);
  }
  if (mp->GetYaxis()) {
    mp->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN_{#font[132]{ch}}}{d#eta}");
    mp->GetYaxis()->SetTitleFont(12);
    mp->GetYaxis()->SetLabelFont(132);
  }
  TLegend* l = c->BuildLegend(0.3, 0.15, 0.7, 0.5, 
			      mp->GetTitle());
  l->SetFillColor(0);
  l->SetTextFont(132);
  l->SetBorderSize(0);
  TLegendEntry* h = static_cast<TLegendEntry*>(l->GetListOfPrimitives()->At(0));
  if (h) { 
    h->SetTextFont(22);
  }
  c->cd();
}

//____________________________________________________________________
//
// EOF
//
