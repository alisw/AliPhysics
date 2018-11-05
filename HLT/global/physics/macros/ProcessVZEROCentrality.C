#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCollection.h"
#include "TFile.h"
#include "TMath.h"
#include <cstdio>
#endif
//______________________________________________________
//
Double_t GetBoundaryForPercentile( TH1F *histo, Double_t lPercentileRequested ) {
    //This function returns the boundary for a specific percentile.
    Double_t lReturnValue = 0.0;
    Double_t lPercentile = 100.0 - lPercentileRequested;
    
    const Long_t lNBins = histo->GetNbinsX();
    Double_t lCountDesired = lPercentile * histo->Integral(1, histo->GetNbinsX())/100;
    Long_t lCount = 0;
    for(Long_t ibin=1;ibin<lNBins;ibin++){
        lCount += histo->GetBinContent(ibin);
        if(lCount >= lCountDesired){
            //Found bin I am looking for!
            Double_t lWidth = histo->GetBinWidth(ibin);
            Double_t lLeftPercentile = 100.*(lCount - histo->GetBinContent(ibin))/histo->Integral(1, histo->GetNbinsX());
            Double_t lRightPercentile = 100.*lCount / histo->Integral(1, histo->GetNbinsX());
            
            Double_t lProportion = (lPercentile - lLeftPercentile)/(lRightPercentile-lLeftPercentile);
            
            lReturnValue = histo->GetBinLowEdge(ibin) + lProportion*lWidth;
            break;
        }
    }
    return lReturnValue;
}

//______________________________________________________
TString GetTriggerFromTitle( TString lTitle ){
    //This function serves to extract trigger selection string from a title
    TString lReturnValue = "N/A";
    if( lTitle.Contains(" ref: ")){
        Int_t lRemInd = lTitle.Index(" ref: ");
        lTitle.Remove(0,lRemInd);
        lTitle.ReplaceAll(" ref: ", "");
        lReturnValue = lTitle;
    }else if( lTitle.Contains(" trig: ")){
        Int_t lRemInd = lTitle.Index(" trig: ");
        lTitle.Remove(0,lRemInd);
        lTitle.ReplaceAll(" trig: ", "");
        lReturnValue = lTitle;
    }
    return lReturnValue;
}

//______________________________________________________
// Reference values obtained from run 246087, LHC15o
//______________________________________________________
const Double_t lReferenceX[] = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000, 6100, 6200, 6300, 6400, 6500, 6600, 6700, 6800, 6900, 7000, 7100, 7200, 7300, 7400, 7500, 7600, 7700, 7800, 7900, 8000, 8100, 8200, 8300, 8400, 8500, 8600, 8700, 8800, 8900, 9000, 9100, 9200, 9300, 9400, 9500, 9600, 9700, 9800, 9900, 10000, 10100, 10200, 10300, 10400, 10500, 10600, 10700, 10800, 10900, 11000, 11100, 11200, 11300, 11400, 11500, 11600, 11700, 11800, 11900, 12000, 12100, 12200, 12300, 12400, 12500, 12600, 12700, 12800, 12900, 13000, 13100, 13200, 13300, 13400, 13500, 13600, 13700, 13800, 13900, 14000, 14100, 14200, 14300, 14400, 14500, 14600, 14700, 14800, 14900, 15000, 15100, 15200, 15300, 15400, 15500, 15600, 15700, 15800, 15900, 16000, 16100, 16200, 16300, 16400, 16500, 16600, 16700, 16800, 16900, 17000, 17100, 17200, 17300, 17400, 17500, 17600, 17700, 17800, 17900, 18000, 18100, 18200, 18300, 18400, 18500, 18600, 18700, 18800, 18900, 19000, 19100, 19200, 19300, 19400, 19500, 19600, 19700, 19800, 19900, 20000, 20100, 20200, 20300, 20400, 20500, 20600, 20700, 20800, 20900, 21000, 21100, 21200, 21300, 21400, 21500, 21600, 21700, 21800, 21900, 22000, 22100, 22200, 22300, 22400, 22500, 22600, 22700, 22800, 22900, 23000, 23100, 23200, 23300, 23400, 23500, 23600, 23700, 23800, 23900, 24000, 24100, 24200, 24300, 24400, 24500, 24600, 24700, 24800, 24900, 25000, 25100, 25200, 25300, 25400, 25500, 25600, 25700, 25800, 25900, 26000, 26100, 26200, 26300, 26400, 26500, 26600, 26700, 26800, 26900, 27000, 27100, 27200, 27300, 27400, 27500, 27600, 27700, 27800, 27900, 28000, 28100, 28200, 28300, 28400, 28500, 28600, 28700, 28800, 28900, 29000, 29100, 29200, 29300, 29400, 29500, 29600, 29700, 29800, 29900, 30000, 30100, 30200, 30300, 30400, 30500, 30600, 30700, 30800, 30900, 31000, 31100, 31200, 31300, 31400, 31500, 31600, 31700, 31800, 31900, 32000, 32100, 32200, 32300, 32400, 32500, 32600, 32700, 32800, 32900, 33000, 33100, 33200, 33300, 33400, 33500, 33600, 33700, 33800, 33900, 34000, 34100, 34200, 34300, 34400, 34500, 34600, 34700, 34800, 34900, 35000, 35100, 35200, 35300, 35400};

//______________________________________________________
const Double_t lReferenceY[] = {0.0667047, 0.0543492, 0.0377883, 0.0294829, 0.0243397, 0.0207749, 0.018108, 0.0162351, 0.014617, 0.0135492, 0.0124073, 0.0116362, 0.010838, 0.0102485, 0.00964084, 0.00920953, 0.00875523, 0.00839838, 0.00803603, 0.00772431, 0.00743073, 0.00720442, 0.00696567, 0.00670288, 0.00648797, 0.00632114, 0.00619881, 0.00602249, 0.00589848, 0.0057114, 0.00561417, 0.00545662, 0.00531637, 0.00526069, 0.00514216, 0.00501287, 0.00491733, 0.00485406, 0.00476611, 0.00465264, 0.00460287, 0.00451112, 0.00446915, 0.00442106, 0.00430063, 0.00427069, 0.00420066, 0.00408994, 0.004064, 0.00400958, 0.00394104, 0.0038959, 0.00386722, 0.00375185, 0.00377125, 0.00366074, 0.00369111, 0.00356203, 0.00357216, 0.00349581, 0.00347978, 0.00346544, 0.00335808, 0.00335155, 0.00330114, 0.00331126, 0.00320813, 0.00323892, 0.00321931, 0.00316257, 0.00314296, 0.00311174, 0.00302422, 0.00305227, 0.00300523, 0.00296158, 0.00295567, 0.00291518, 0.00285317, 0.00282892, 0.00279686, 0.00281626, 0.00279517, 0.00275489, 0.00272241, 0.00270553, 0.0026585, 0.00268318, 0.00263446, 0.00258426, 0.00259544, 0.00261716, 0.00256444, 0.00250981, 0.00249146, 0.00250939, 0.00244253, 0.00245624, 0.0024322, 0.00238348, 0.00241279, 0.00239023, 0.00234446, 0.00235036, 0.00231261, 0.0023105, 0.00229595, 0.00225799, 0.002293, 0.00226642, 0.00223099, 0.00223964, 0.00221496, 0.00217974, 0.00222382, 0.00215907, 0.00214768, 0.00210297, 0.00212343, 0.00208546, 0.00205404, 0.00205404, 0.0020494, 0.00201439, 0.00206964, 0.00202367, 0.00203168, 0.00197727, 0.0019798, 0.00195175, 0.00198043, 0.00195744, 0.00192285, 0.00192981, 0.00188869, 0.00188679, 0.00191357, 0.00188742, 0.00185072, 0.00187709, 0.00184144, 0.00184081, 0.00181655, 0.00183069, 0.00176045, 0.00177205, 0.0017961, 0.001771, 0.00175518, 0.00175645, 0.001713, 0.00174274, 0.00170984, 0.00170688, 0.00166681, 0.00168495, 0.00164277, 0.00166618, 0.0016415, 0.00160923, 0.00162125, 0.0016067, 0.00160733, 0.00161851, 0.00159489, 0.00160522, 0.00160501, 0.00158962, 0.00157422, 0.00154258, 0.0015428, 0.00153141, 0.00154722, 0.00153162, 0.00150104, 0.00150757, 0.00147741, 0.00148986, 0.00148037, 0.00145949, 0.00146223, 0.00146012, 0.00144894, 0.0014173, 0.00143586, 0.00145084, 0.00141393, 0.0014211, 0.00141625, 0.00139664, 0.00138968, 0.00134243, 0.00138736, 0.00141414, 0.00134138, 0.00136352, 0.00131987, 0.00130805, 0.00134032, 0.00129856, 0.00131059, 0.00129899, 0.00129224, 0.001307, 0.00131059, 0.00128802, 0.00127579, 0.00124499, 0.00126777, 0.00125469, 0.00126651, 0.00123529, 0.00124731, 0.00121757, 0.00123867, 0.00120872, 0.00121378, 0.00117708, 0.00119395, 0.0011775, 0.00120724, 0.00118763, 0.00120555, 0.00115725, 0.00115261, 0.00117455, 0.00114481, 0.00115093, 0.00113806, 0.00112351, 0.00111465, 0.00110284, 0.00110094, 0.00110115, 0.00112182, 0.00110368, 0.00108576, 0.00109272, 0.00108238, 0.00109419, 0.00104231, 0.00109778, 0.00103619, 0.00107289, 0.00104336, 0.0010092, 0.00104126, 0.00103577, 0.00102143, 0.00104421, 0.00101531, 0.000991903, 0.00101426, 0.00101363, 0.000969968, 0.000987474, 0.000963641, 0.000990004, 0.000967437, 0.000953095, 0.000946135, 0.000985786, 0.000965328, 0.000929895, 0.000925256, 0.000914288, 0.000921037, 0.000924201, 0.000905008, 0.000923146, 0.000914921, 0.000909016, 0.000896783, 0.000887081, 0.000893619, 0.000892354, 0.000882019, 0.000864936, 0.000864092, 0.000836674, 0.000862194, 0.000854812, 0.000821489, 0.000834354, 0.000837939, 0.00083351, 0.000843212, 0.000825496, 0.000834354, 0.00080103, 0.000824863, 0.000794492, 0.000808623, 0.00080567, 0.000759271, 0.000775721, 0.00078458, 0.000752943, 0.000745561, 0.000739656, 0.000743031, 0.000724682, 0.000727002, 0.00071962, 0.000705067, 0.000699583, 0.000665627, 0.000641795, 0.000650653, 0.000626398, 0.000619649, 0.000598558, 0.000574726, 0.000548573, 0.000532122, 0.000512086, 0.000463366, 0.000449446, 0.000444595, 0.00040853, 0.000385119, 0.000349264, 0.000339352, 0.000300334, 0.000277134, 0.000255621, 0.000233265, 0.000204159, 0.000174843, 0.000165141, 0.000152065, 0.000137091, 0.000117265, 9.59634e-05, 7.97234e-05, 7.57161e-05, 6.15853e-05, 5.46253e-05, 4.49235e-05, 3.64872e-05, 3.64872e-05, 2.38327e-05, 2.44654e-05, 1.75054e-05, 1.624e-05, 1.392e-05, 9.49088e-06, 8.22543e-06, 7.3818e-06, 3.58544e-06, 3.16363e-06, 2.5309e-06, 2.5309e-06, 2.95272e-06};

const Long_t lNRefBins = sizeof (lReferenceY) / sizeof(Double_t);

//______________________________________________________
// Online Reference values obtained from run 246087, LHC15o
//______________________________________________________
const Double_t lOnlineReferenceX[] = {0, 125, 250, 375, 500, 625, 750, 875, 1000, 1125, 1250, 1375, 1500, 1625, 1750, 1875, 2000, 2125, 2250, 2375, 2500, 2625, 2750, 2875, 3000, 3125, 3250, 3375, 3500, 3625, 3750, 3875, 4000, 4125, 4250, 4375, 4500, 4625, 4750, 4875, 5000, 5125, 5250, 5375, 5500, 5625, 5750, 5875, 6000, 6125, 6250, 6375, 6500, 6625, 6750, 6875, 7000, 7125, 7250, 7375, 7500, 7625, 7750, 7875, 8000, 8125, 8250, 8375, 8500, 8625, 8750, 8875, 9000, 9125, 9250, 9375, 9500, 9625, 9750, 9875, 10000, 10125, 10250, 10375, 10500, 10625, 10750, 10875, 11000, 11125, 11250, 11375, 11500, 11625, 11750, 11875, 12000, 12125, 12250, 12375, 12500, 12625, 12750, 12875, 13000, 13125, 13250, 13375, 13500, 13625, 13750, 13875, 14000, 14125, 14250, 14375, 14500, 14625, 14750, 14875, 15000, 15125, 15250, 15375, 15500, 15625, 15750, 15875, 16000, 16125, 16250, 16375, 16500, 16625, 16750, 16875, 17000, 17125, 17250, 17375, 17500, 17625, 17750, 17875, 18000, 18125, 18250, 18375, 18500, 18625, 18750, 18875, 19000, 19125, 19250, 19375, 19500, 19625, 19750, 19875, 20000, 20125, 20250, 20375, 20500, 20625, 20750, 20875, 21000, 21125, 21250, 21375, 21500, 21625, 21750, 21875, 22000, 22125, 22250, 22375, 22500, 22625, 22750, 22875, 23000, 23125, 23250, 23375, 23500, 23625, 23750, 23875, 24000, 24125, 24250, 24375, 24500, 24625, 24750, 24875, 25000, 25125, 25250, 25375, 25500, 25625, 25750, 25875, 26000, 26125, 26250, 26375, 26500, 26625, 26750, 26875, 27000, 27125, 27250, 27375, 27500, 27625, 27750, 27875, 28000, 28125, 28250, 28375, 28500, 28625, 28750, 28875, 29000, 29125, 29250, 29375, 29500, 29625, 29750, 29875, 30000, 30125, 30250, 30375, 30500, 30625, 30750, 30875, 31000, 31125, 31250, 31375, 31500, 31625, 31750, 31875, 32000, 32125, 32250, 32375, 32500, 32625};

//______________________________________________________
const Double_t lOnlineReferenceY[] = {0.0875729, 0.0549932, 0.039097, 0.0306405, 0.0255665, 0.0215825, 0.0190856, 0.017726, 0.0155973, 0.014121, 0.0135159, 0.0127031, 0.0115658, 0.0110664, 0.0103119, 0.00984169, 0.00944803, 0.00939335, 0.00893772, 0.00829618, 0.00822328, 0.00797906, 0.00775671, 0.00726098, 0.00713705, 0.00711154, 0.00688554, 0.0065101, 0.00625494, 0.00638981, 0.00618204, 0.00622943, 0.00607269, 0.00602895, 0.00577744, 0.00536919, 0.00551499, 0.00539471, 0.00524161, 0.00505207, 0.00492814, 0.00502291, 0.00482243, 0.00468392, 0.00487346, 0.00434493, 0.00453811, 0.00450895, 0.00427931, 0.00444334, 0.00430848, 0.00435951, 0.00423557, 0.00410071, 0.00407519, 0.00398771, 0.00405697, 0.003882, 0.00376171, 0.00403874, 0.0038091, 0.00369975, 0.00364872, 0.00370339, 0.00349927, 0.0034774, 0.00361227, 0.00360498, 0.00349562, 0.00326963, 0.00334618, 0.00318579, 0.00336076, 0.00334982, 0.00311289, 0.00331337, 0.00317486, 0.00319673, 0.00296344, 0.00294522, 0.00299989, 0.00299989, 0.00295615, 0.00312747, 0.00291606, 0.0029598, 0.00287961, 0.00291241, 0.00289783, 0.00274474, 0.00263174, 0.00254791, 0.00274838, 0.00259894, 0.00264632, 0.00265361, 0.00259894, 0.002588, 0.00270829, 0.00246771, 0.00269371, 0.00256613, 0.00248594, 0.002475, 0.00240939, 0.00243855, 0.00248229, 0.00252239, 0.00247136, 0.00254426, 0.00233285, 0.00240575, 0.00228182, 0.00248958, 0.00238388, 0.00229275, 0.00225994, 0.0021433, 0.00225994, 0.00224172, 0.0021834, 0.00219069, 0.00217611, 0.00216517, 0.00225265, 0.00203395, 0.00220162, 0.00204489, 0.00228911, 0.00206676, 0.00201572, 0.00217611, 0.0018845, 0.00199385, 0.00207405, 0.00200843, 0.0019975, 0.00194647, 0.00193918, 0.00193918, 0.00189908, 0.00181525, 0.00182983, 0.00200479, 0.00172047, 0.0017387, 0.00174963, 0.00175328, 0.00179338, 0.00169131, 0.00161112, 0.0017788, 0.00176057, 0.0018517, 0.00162206, 0.00179702, 0.00153093, 0.00168767, 0.00157467, 0.00167309, 0.00162935, 0.00172047, 0.00160383, 0.00158561, 0.00145438, 0.00162935, 0.00152364, 0.00150906, 0.00157467, 0.00142522, 0.00144345, 0.00156374, 0.00165122, 0.00160748, 0.00152, 0.00147261, 0.00138148, 0.00157103, 0.00149812, 0.0015528, 0.00147625, 0.00144709, 0.00149448, 0.00130858, 0.00132316, 0.0014799, 0.00134503, 0.00145438, 0.00141793, 0.00127942, 0.00134503, 0.00135232, 0.00143251, 0.00123568, 0.00132316, 0.00123203, 0.00137419, 0.0012539, 0.00120287, 0.00127578, 0.00117007, 0.00130129, 0.00126119, 0.00130494, 0.00129036, 0.001181, 0.00126484, 0.0012539, 0.00115549, 0.00114091, 0.00115549, 0.00113726, 0.00109352, 0.0011081, 0.00111904, 0.00118829, 0.00102427, 0.00115184, 0.00106436, 0.00110081, 0.00112268, 0.00111539, 0.00105343, 0.00103885, 0.000969589, 0.000980524, 0.000918558, 0.000885753, 0.000882107, 0.000962299, 0.000852947, 0.000878462, 0.000951364, 0.000914913, 0.00076911, 0.000878462, 0.000736305, 0.000798271, 0.000714434, 0.000670693, 0.000597792, 0.000612372, 0.000546761, 0.000550406, 0.00051031, 0.000459279, 0.000298896, 0.000331702, 0.000298896, 0.00026609, 0.000222349, 0.000262445, 0.000164028, 0.000138513, 0.000116642, 8.01916e-05, 5.46761e-05, 4.37409e-05, 4.37409e-05, 3.64507e-05};

const Long_t lNOnlineRefBins = sizeof (lOnlineReferenceY) / sizeof(Double_t);



TH1F *hReferenceV0M(NULL);
TH1F *hOnlineReferenceV0M(NULL);

//______________________________________________________
Double_t pdf_PbPb_kno(Double_t *x, Double_t *par)
//Calculation of probability of a certain amplitude for Pb-Pb
//based on simple X/Y axis scaling and known distribution
{
    Double_t lMultValue = par[0]*x[0]; //par[0] -> X scaling
    
    if( lMultValue < 129.59 ){
        //This point is outside the fitting range
        //(in the original reference X!)
        TF1::RejectPoint();
        return 0;
    }
    
    return par[1]*hReferenceV0M->GetBinContent(hReferenceV0M->FindBin(lMultValue));        //par[1] -> Y scaling
}

//______________________________________________________
Double_t pdf_PbPb_onlinekno(Double_t *x, Double_t *par)
//Calculation of probability of a certain amplitude for Pb-Pb
//based on simple X/Y axis scaling and known distribution
{
    Double_t lMultValue = par[0]*x[0]; //par[0] -> X scaling
    
    if( lMultValue < 129.59 ){
        //This point is outside the fitting range
        //(in the original reference X!)
        TF1::RejectPoint();
        return 0;
    }
    
    return par[1]*hOnlineReferenceV0M->GetBinContent(hOnlineReferenceV0M->FindBin(lMultValue));        //par[1] -> Y scaling
}

Int_t process(TCollection *lIn = 0x0, TCollection *lOut = 0x0) {
    
    TString lStudiedTriggers[2] = {"V0M", "V0MOnline"};
    TFile *file = 0x0;
    
    if( !lIn ){
        //For debugging purposes: grab from actual file
        file = new TFile("TriggerExample.root", "READ") ;
        if(!file) {
            cout<<"No input specified and example file not found! Exiting!"<<endl;
            return -1;
        }
        
        lIn = (TCollection*) file -> Get("lInput");
        cout<<"Opened input from example file. Listing: "<<endl;
        lIn->ls();
    }
    
    //Save to TCollection inside file (testing)
    TFile *fFile = 0x0;
    
    if( !lOut ){
        lOut = new TObjArray();
        fFile = new TFile("testingontheflycentrality.root", "RECREATE");
        if(!fFile){
            cout<<"Output not viable, exiting!"<<endl;
            return -1;
        }
    }
    
    Int_t lAtLine = 1;
    TH2F *fThresholdTable = new TH2F("fThresholdTable", "", 3, 0, 3, 4, 0, 4);
    fThresholdTable->GetXaxis()->SetBinLabel(1, "50%");
    fThresholdTable->GetXaxis()->SetBinLabel(2, "30%");
    fThresholdTable->GetXaxis()->SetBinLabel(3, "10%");
    
    fThresholdTable->SetMarkerSize(2.5);
    fThresholdTable->SetOption("text0");
    fThresholdTable->SetStats(kFALSE);
    
    //Prepare reference histogram in case it does not exist
    if(!hReferenceV0M){
        hReferenceV0M = new TH1F("fHistReferenceV0M", "", lNRefBins, lReferenceX);
        for(Int_t i=1; i<hReferenceV0M->GetNbinsX()+1; i++)
            hReferenceV0M->SetBinContent(i, lReferenceY[i-1]);
        hReferenceV0M->SetDirectory(0);
    }
    if(!hOnlineReferenceV0M){
        hOnlineReferenceV0M = new TH1F("fHistReferenceV0MOnline", "", lNOnlineRefBins, lOnlineReferenceX);
        for(Int_t i=1; i<hOnlineReferenceV0M->GetNbinsX()+1; i++)
            hOnlineReferenceV0M->SetBinContent(i, lOnlineReferenceY[i-1]);
        hOnlineReferenceV0M->SetDirectory(0);
    }
    
    for(Int_t itrig=0; itrig<2; itrig++){
        //Get basic macros
        TH1F *histo = (TH1F*) lIn->FindObject(Form("fHistMult%s",lStudiedTriggers[itrig].Data()));
        TH1F *histoc = (TH1F*) lIn->FindObject(Form("fHistMult%s_Central",lStudiedTriggers[itrig].Data()));
        TH1F *histosc = (TH1F*) lIn->FindObject(Form("fHistMult%s_SemiCentral",lStudiedTriggers[itrig].Data()));
        TH1F *histoextra = (TH1F*) lIn->FindObject(Form("fHistMult%s_Extra",lStudiedTriggers[itrig].Data()));
        
        histo->SetLineColor(kBlack);
        histoc->SetLineColor(kRed+1);
        histosc->SetLineColor(kGreen+1);
        histoextra->SetLineColor(kBlue+1);
        
        if( !histo || !histoc || !histosc || !histoextra ) {
            cout<<"Received incomplete input for type "<<lStudiedTriggers[itrig]<<". Skipping."<<endl;
            printf("histo %p histoc %p histosc %p histoextra %p\n",histo, histoc, histosc, histoextra);
            continue;
        }
        
        Float_t lCentBounds[101];
        Float_t lCentBoundsRaw[101];
        
        cout<<"=================================================="<<endl;
        cout<<" Step 1: now doing unanchored calibration of "<<lStudiedTriggers[itrig]<<endl;
        cout<<"=================================================="<<endl;
        
        cout<<"Determining calibration on-the-fly..."<<endl;
        lCentBoundsRaw[0] = 0.0;
        cout<<"Inspect boundaries: "<<flush;
        for(Int_t i=0; i<100; i++){
            lCentBounds[i] = i;
            if (i!=0) lCentBoundsRaw[i] = GetBoundaryForPercentile( histo,  100-lCentBounds[i] );
            cout<<lCentBoundsRaw[i]<<" "<<flush;
        }
        lCentBoundsRaw[100] = 1e+5;
        cout<<lCentBoundsRaw[100]<<" (end)"<<endl;
        
        cout<<"Classifying..."<<endl;
        TH1F* fHistCentrality             = new TH1F(Form("fHistCentrality%s",lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        TH1F* fHistCentrality_Central     = new TH1F(Form("fHistCentrality%s_Central",lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        TH1F* fHistCentrality_SemiCentral = new TH1F(Form("fHistCentrality%s_SemiCentral",lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        TH1F* fHistCentrality_Extra       = new TH1F(Form("fHistCentrality%s_Extra",lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        
        fHistCentrality->SetStats(0);
        fHistCentrality_Central->SetStats(0);
        fHistCentrality_SemiCentral->SetStats(0);
        fHistCentrality_Extra->SetStats(0);
        
        //Some basic cosmetics
        fHistCentrality->SetMarkerStyle(20);
        fHistCentrality->SetMarkerSize(0.6);
        fHistCentrality->SetMarkerColor(kBlack);
        fHistCentrality->SetLineColor(kBlack);
        
        fHistCentrality_Central->SetMarkerStyle(21);
        fHistCentrality_Central->SetMarkerSize(0.7);
        fHistCentrality_Central->SetMarkerColor(kRed+1);
        fHistCentrality_Central->SetLineColor(kRed+1);
        
        fHistCentrality_SemiCentral->SetMarkerStyle(21);
        fHistCentrality_SemiCentral->SetMarkerSize(0.7);
        fHistCentrality_SemiCentral->SetMarkerColor(kGreen+1);
        fHistCentrality_SemiCentral->SetLineColor(kGreen+1);
        
        fHistCentrality_Extra->SetMarkerStyle(21);
        fHistCentrality_Extra->SetMarkerSize(0.7);
        fHistCentrality_Extra->SetMarkerColor(kBlue+1);
        fHistCentrality_Extra->SetLineColor(kBlue+1);
        
        fHistCentrality->SetTitle(Form("%s, cent. distrib. for reference sample (unanchored)", histo->GetTitle()));
        fHistCentrality_Central->SetTitle(Form("%s, cent. distrib. for central  (unanchored)", histoc->GetTitle()));
        fHistCentrality_SemiCentral->SetTitle(Form("%s, cent. distrib. for semicentral (unanchored)", histosc->GetTitle()));
        fHistCentrality_Extra->SetTitle(Form("%s, cent. distrib. for extra trigger (unanchored)", histoextra->GetTitle()));
        
        fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
        fHistCentrality_Central->GetXaxis()->SetTitle("Centrality (%)");
        fHistCentrality_SemiCentral->GetXaxis()->SetTitle("Centrality (%)");
        fHistCentrality_Extra->GetXaxis()->SetTitle("Centrality (%)");
        
        fHistCentrality->GetYaxis()->SetTitle("Event count");
        fHistCentrality_Central->GetYaxis()->SetTitle("Event count");
        fHistCentrality_SemiCentral->GetYaxis()->SetTitle("Event count");
        fHistCentrality_Extra->GetYaxis()->SetTitle("Event count");
        
        fHistCentrality->Sumw2();
        fHistCentrality_Central->Sumw2();
        fHistCentrality_SemiCentral->Sumw2();
        
        TH1F *hClassifier = new TH1F(Form("hClassifier_%s", lStudiedTriggers[itrig].Data()), "", 100, lCentBoundsRaw);
        for(Int_t i=0; i<101; i++){
            hClassifier -> SetBinContent( i+1 , 99.5 - i );
        }
        hClassifier->SetTitle(Form("%s <-> centrality map (unanchored)", lStudiedTriggers[itrig].Data()));
        hClassifier->GetYaxis()->SetTitle("Centrality (%)");
        hClassifier->GetXaxis()->SetTitle(Form("%s signal (raw)", lStudiedTriggers[itrig].Data()));
        hClassifier->SetStats(0);
        
        for(Long_t i=1; i<histo->GetNbinsX()+1; i++) {
            Float_t lCentrality = hClassifier->GetBinContent(hClassifier->FindBin(histo->GetBinCenter(i)));
            fHistCentrality->Fill( lCentrality, histo->GetBinContent(i));
            fHistCentrality_Central->Fill( lCentrality, histoc->GetBinContent(i));
            fHistCentrality_SemiCentral->Fill( lCentrality, histosc->GetBinContent(i));
            fHistCentrality_Extra->Fill( lCentrality, histoextra->GetBinContent(i));
        }
        for(Long_t i=1; i<histo->GetNbinsX()+1; i++) {
            fHistCentrality->SetBinError(i, TMath::Sqrt(fHistCentrality->GetBinContent(i)) );
            fHistCentrality_Central->SetBinError(i, TMath::Sqrt(fHistCentrality_Central->GetBinContent(i)) );
            fHistCentrality_SemiCentral->SetBinError(i, TMath::Sqrt(fHistCentrality_SemiCentral->GetBinContent(i)) );
            fHistCentrality_Extra->SetBinError(i, TMath::Sqrt(fHistCentrality_Extra->GetBinContent(i)) );
        }
        
        cout<<"=================================================="<<endl;
        cout<<" Step 2: now doing anchoring via KNO-scaling of "<<lStudiedTriggers[itrig]<<endl;
        cout<<"=================================================="<<endl;
        
        //This requires an input distribution from previous data.
        //The reference chosen for this purpose is the V0M amplitude
        //distribution from run 246087, a typical large run within LHC15o.
        
        //This run is anchored at a V0M amplitude of 129.59 in our
        //usual centrality calibration.
        
        TF1 *f1_kno = 0x0;
        
        if( lStudiedTriggers[itrig].Contains("Online") == kFALSE ){
            f1_kno = new TF1(Form("f1_kno_%s",lStudiedTriggers[itrig].Data()),pdf_PbPb_kno,0,50000,2);
        }else{
            f1_kno = new TF1(Form("f1_kno_%s",lStudiedTriggers[itrig].Data()),pdf_PbPb_onlinekno,0,50000,2);
        }
        f1_kno->SetTitle(Form("ref fit to %s signal (raw);%s signal (raw);f()", lStudiedTriggers[itrig].Data(), lStudiedTriggers[itrig].Data()));
        
        cout<<"Guessing normalization to: "<<histo->GetEntries()<<endl;
        
        //guess X-axis scaling factor
        f1_kno->SetParameter(0, 1.0);
        
        //guess normalization, 25 -> comes from rebinning done at extraction
        f1_kno->SetParameter(1, histo->GetEntries()/25.);
        
        f1_kno->SetParName(0, "X-scaling factor");
        f1_kno->SetParName(1, "Y-scaling factor");
        f1_kno->SetNpx(1000);
        
        cout<<"Fitting... "<<endl;
        histo->Fit(Form("f1_kno_%s",lStudiedTriggers[itrig].Data()),"R0W");
        
        //Create a clone for drawing purposes
        cout<<"Cloning f1... ";
        TF1 *f1_kno_drawclone = (TF1*) f1_kno->Clone( Form("f1_kno_drawclone_%s",lStudiedTriggers[itrig].Data()) );
        f1_kno_drawclone->SetLineColor(kRed+1);
        
        cout<<"Cloning histo... ";
        TH1F *Plot1 = (TH1F*) histo->Clone(Form("Plot1_%s", lStudiedTriggers[itrig].Data()));
        Plot1->GetListOfFunctions()->Add(f1_kno_drawclone);
        Plot1->SetStats(kFALSE);
        
        Double_t lAnchorPointKNO = f1_kno->GetParameter(0)*129.59;
        
        //Scale back if trigger charge, please !
        if( lStudiedTriggers[itrig].Contains("Online") ) lAnchorPointKNO *= 7147./7469.;
        
        cout<<"Determined 90% anchor point via KNO-scaling: "<<lAnchorPointKNO<<endl;
        
        cout<<"Now creating anchored histogram..."<<endl;
        
        TH1F *histo_KNOanchored = (TH1F*) histo->Clone(Form("fHistMultKNO%s",lStudiedTriggers[itrig].Data()));
        
        histo_KNOanchored->GetYaxis()->SetTitle("Event count");
        histo_KNOanchored->GetXaxis()->SetTitle(Form("%s signal (raw)",lStudiedTriggers[itrig].Data()));
        histo_KNOanchored->SetTitle(Form("Anchored ref histo for %s", histo->GetTitle()));
        histo_KNOanchored->SetLineColor(kBlack);
        histo_KNOanchored->SetStats(kFALSE);
        
        //Find the bin to anchor to (precision: +/- 1 bin)
        Int_t lAnchorBinKNO = histo_KNOanchored->FindBin(lAnchorPointKNO);
        
        //Count stuff
        Double_t lNEntries = histo_KNOanchored->Integral(1, histo_KNOanchored->GetNbinsX() );
        Double_t lNEntriesPastAnchor = histo_KNOanchored->Integral(lAnchorBinKNO, histo_KNOanchored->GetNbinsX() );
        
        //Zero out everything
        for(Int_t i=1; i<lAnchorBinKNO; i++) histo_KNOanchored->SetBinContent(i,0.0);
        
        //Calculate the expected anchored-out quantity
        Double_t lAnchoredOutKNO = (1./9.) * lNEntriesPastAnchor;
        
        //Put everything in bin 1
        histo_KNOanchored->SetBinContent(1,lAnchoredOutKNO);
        
        //Repeat procedure done before
        cout<<"Classifying after KNO-scaling anchoring..."<<endl;
        TH1F* fHistCentralityKNO             = new TH1F(Form("fHistCentralityKNO%s",lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        TH1F* fHistCentralityKNO_Central     = new TH1F(Form("fHistCentralityKNO%s_Central",lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        TH1F* fHistCentralityKNO_SemiCentral = new TH1F(Form("fHistCentralityKNO%s_SemiCentral",lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        TH1F* fHistCentralityKNO_Extra       = new TH1F(Form("fHistCentralityKNO%s_Extra",lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        
        fHistCentralityKNO->SetStats(0);
        fHistCentralityKNO_Central->SetStats(0);
        fHistCentralityKNO_SemiCentral->SetStats(0);
        fHistCentralityKNO_Extra->SetStats(0);
        
        //Some basic cosmetics
        fHistCentralityKNO->SetMarkerStyle(20);
        fHistCentralityKNO->SetMarkerSize(0.6);
        fHistCentralityKNO->SetMarkerColor(kBlack);
        fHistCentralityKNO->SetLineColor(kBlack);
        
        fHistCentralityKNO_Central->SetMarkerStyle(21);
        fHistCentralityKNO_Central->SetMarkerSize(0.7);
        fHistCentralityKNO_Central->SetMarkerColor(kRed+1);
        fHistCentralityKNO_Central->SetLineColor(kRed+1);
        
        fHistCentralityKNO_SemiCentral->SetMarkerStyle(21);
        fHistCentralityKNO_SemiCentral->SetMarkerSize(0.7);
        fHistCentralityKNO_SemiCentral->SetMarkerColor(kGreen+1);
        fHistCentralityKNO_SemiCentral->SetLineColor(kGreen+1);
        
        fHistCentralityKNO_Extra->SetMarkerStyle(21);
        fHistCentralityKNO_Extra->SetMarkerSize(0.7);
        fHistCentralityKNO_Extra->SetMarkerColor(kBlue+1);
        fHistCentralityKNO_Extra->SetLineColor(kBlue+1);
        
        fHistCentralityKNO->SetTitle(Form("%s, cent. distrib. for reference sample (anchored)", histo->GetTitle()));
        fHistCentralityKNO_Central->SetTitle(Form("%s, cent. distrib. for central  (anchored)", histoc->GetTitle()));
        fHistCentralityKNO_SemiCentral->SetTitle(Form("%s, cent. distrib. for semicentral (anchored)", histosc->GetTitle()));
        fHistCentralityKNO_Extra->SetTitle(Form("%s, cent. distrib. for extra trigger (anchored)", histoextra->GetTitle()));
        
        fHistCentralityKNO->GetXaxis()->SetTitle("Centrality (%)");
        fHistCentralityKNO_Central->GetXaxis()->SetTitle("Centrality (%)");
        fHistCentralityKNO_SemiCentral->GetXaxis()->SetTitle("Centrality (%)");
        fHistCentralityKNO_Extra->GetXaxis()->SetTitle("Centrality (%)");
        
        fHistCentralityKNO->GetYaxis()->SetTitle("Event count");
        fHistCentralityKNO_Central->GetYaxis()->SetTitle("Event count");
        fHistCentralityKNO_SemiCentral->GetYaxis()->SetTitle("Event count");
        fHistCentralityKNO_Extra->GetYaxis()->SetTitle("Event count");
        
        cout<<"Determining calibration on-the-fly..."<<endl;
        Float_t lCentBoundsRawKNO[101];
        lCentBoundsRawKNO[0] = 0.0;
        cout<<"Inspect boundaries: "<<flush;
        for(Int_t i=0; i<100; i++){
            //lCentBounds[i] = i;
            if (i!=0) lCentBoundsRawKNO[i] = GetBoundaryForPercentile( histo_KNOanchored,  100-lCentBounds[i] );
            cout<<lCentBoundsRawKNO[i]<<" "<<flush;
        }
        lCentBoundsRawKNO[100] = 1e+5;
        cout<<lCentBoundsRawKNO[100]<<" (end)"<<endl;
        
        TH1F *hClassifierKNO = new TH1F(Form("hClassifierKNO_%s",lStudiedTriggers[itrig].Data()), "", 100, lCentBoundsRawKNO);
        for(Int_t i=0; i<101; i++){
            hClassifierKNO -> SetBinContent( i+1 , 99.5 - i );
        }
        hClassifierKNO->SetTitle(Form("%s <-> centrality map (anchored)", lStudiedTriggers[itrig].Data()));
        hClassifierKNO->GetYaxis()->SetTitle("Centrality (%)");
        hClassifierKNO->GetXaxis()->SetTitle(Form("%s signal (raw)", lStudiedTriggers[itrig].Data()));
        hClassifierKNO->SetStats(kFALSE);
        
        for(Long_t i=1; i<histo->GetNbinsX()+1; i++) {
            //This is the anchored-out centrality, no extra action needed
            Float_t lCentrality = hClassifierKNO->GetBinContent(hClassifierKNO->FindBin(histo_KNOanchored->GetBinCenter(i)));
            if(lCentrality > 90) lCentrality = 101;
            fHistCentralityKNO->Fill( lCentrality, histo->GetBinContent(i));
            fHistCentralityKNO_Central->Fill( lCentrality, histoc->GetBinContent(i));
            fHistCentralityKNO_SemiCentral->Fill( lCentrality, histosc->GetBinContent(i));
            fHistCentralityKNO_Extra->Fill( lCentrality, histoextra->GetBinContent(i));
        }
        for(Long_t i=1; i<histo->GetNbinsX()+1; i++) {
            fHistCentralityKNO->SetBinError(i, TMath::Sqrt(fHistCentralityKNO->GetBinContent(i)) );
            fHistCentralityKNO_Central->SetBinError(i, TMath::Sqrt(fHistCentralityKNO_Central->GetBinContent(i)) );
            fHistCentralityKNO_SemiCentral->SetBinError(i, TMath::Sqrt(fHistCentralityKNO_SemiCentral->GetBinContent(i)) );
            fHistCentralityKNO_Extra->SetBinError(i, TMath::Sqrt(fHistCentralityKNO_Extra->GetBinContent(i)) );
        }
        
        //Final step: please copy reference histograms (watch out for similarly named histograms!)
        TH1F *histo_ref = (TH1F*) histo->Clone(Form("fHistRefMult%s",lStudiedTriggers[itrig].Data()));
        TH1F *histoc_ref = (TH1F*) histoc->Clone(Form("fHistRefMult%s_Central",lStudiedTriggers[itrig].Data()));
        TH1F *histosc_ref = (TH1F*) histosc->Clone(Form("fHistRefMult%s_SemiCentral",lStudiedTriggers[itrig].Data()));
        TH1F *histoextra_ref = (TH1F*) histoextra->Clone(Form("fHistRefMult%s_Extra",lStudiedTriggers[itrig].Data()));
        
        //histo_ref->SetTitle(Form("%s signal distrib (reference)", lStudiedTriggers[itrig].Data()));
        //histoc_ref->SetTitle(Form("%s signal distrib (central)", lStudiedTriggers[itrig].Data()));
        //histosc_ref->SetTitle(Form("%s signal distrib (semicentral)", lStudiedTriggers[itrig].Data()));
        //histoextra_ref->SetTitle(Form("%s signal distrib (extra)", lStudiedTriggers[itrig].Data()));
        
        histo_ref->GetYaxis()->SetTitle("Event counts");
        histoc_ref->GetYaxis()->SetTitle("Event counts");
        histosc_ref->GetYaxis()->SetTitle("Event counts");
        histoextra_ref->GetYaxis()->SetTitle("Event counts");
        
        histo_ref->GetXaxis()->SetTitle(Form("%s signal (raw)",lStudiedTriggers[itrig].Data()));
        histoc_ref->GetXaxis()->SetTitle(Form("%s signal (raw)",lStudiedTriggers[itrig].Data()));
        histosc_ref->GetXaxis()->SetTitle(Form("%s signal (raw)",lStudiedTriggers[itrig].Data()));
        histoextra_ref->GetXaxis()->SetTitle(Form("%s signal (raw)",lStudiedTriggers[itrig].Data()));
        
        //Bis final step: save histo with reference thresholds, please
        TH1F *fThresholds = new TH1F(Form("fThresholds%s",lStudiedTriggers[itrig].Data()),"",  3, 0, 3);
        fThresholds->GetXaxis()->SetBinLabel(1, "50%");
        fThresholds->GetXaxis()->SetBinLabel(2, "30%");
        fThresholds->GetXaxis()->SetBinLabel(3, "10%");
        fThresholds->GetXaxis()->SetLabelSize(0.06);
        fThresholds->GetYaxis()->SetTitle(Form("%s signal (raw)",lStudiedTriggers[itrig].Data()));
        fThresholds->SetMarkerSize(2.5);
        fThresholds->SetLineColor(kBlack);
        fThresholds->SetBinContent( 1, GetBoundaryForPercentile( histo,  50 ) );
        fThresholds->SetBinContent( 2, GetBoundaryForPercentile( histo,  30 ) );
        fThresholds->SetBinContent( 3, GetBoundaryForPercentile( histo,  10 ) );
        fThresholds->SetOption("text0");
        fThresholds->SetTitle("Centrality thresholds (unanchored)");
        fThresholds->SetStats(kFALSE);
        
        //Bis final step: save histo with reference thresholds, please
        TH1F *fThresholdsKNO = new TH1F(Form("fThresholdsKNO%s",lStudiedTriggers[itrig].Data()),"",  3, 0, 3);
        fThresholdsKNO->GetXaxis()->SetBinLabel(1, "50%");
        fThresholdsKNO->GetXaxis()->SetBinLabel(2, "30%");
        fThresholdsKNO->GetXaxis()->SetBinLabel(3, "10%");
        fThresholdsKNO->GetXaxis()->SetLabelSize(0.06);
        fThresholdsKNO->GetYaxis()->SetTitle(Form("%s signal (raw)",lStudiedTriggers[itrig].Data()));
        fThresholdsKNO->SetMarkerSize(2.5);
        fThresholdsKNO->SetLineColor(kBlack);
        fThresholdsKNO->SetBinContent( 1, GetBoundaryForPercentile( histo_KNOanchored,  50 ) );
        fThresholdsKNO->SetBinContent( 2, GetBoundaryForPercentile( histo_KNOanchored,  30 ) );
        fThresholdsKNO->SetBinContent( 3, GetBoundaryForPercentile( histo_KNOanchored,  10 ) );
        fThresholdsKNO->SetOption("text0");
        fThresholdsKNO->SetTitle("Centrality thresholds (anchored)");
        fThresholdsKNO->SetStats(kFALSE);
        
        //table, please
        fThresholdTable->SetBinContent( 1, 2*itrig+1,  GetBoundaryForPercentile( histo,  50 ) );
        fThresholdTable->SetBinContent( 2, 2*itrig+1,  GetBoundaryForPercentile( histo,  30 ) );
        fThresholdTable->SetBinContent( 3, 2*itrig+1,  GetBoundaryForPercentile( histo,  10 ) );
        
        fThresholdTable->SetBinContent( 1, 2*itrig+2,  GetBoundaryForPercentile( histo_KNOanchored,  50 ) );
        fThresholdTable->SetBinContent( 2, 2*itrig+2,  GetBoundaryForPercentile( histo_KNOanchored,  30 ) );
        fThresholdTable->SetBinContent( 3, 2*itrig+2,  GetBoundaryForPercentile( histo_KNOanchored,  10 ) );
        
        fThresholdTable->GetYaxis()->SetBinLabel(2*itrig+1, Form("%s unanchored", lStudiedTriggers[itrig].Data()));
        fThresholdTable->GetYaxis()->SetBinLabel(2*itrig+2, Form("%s anchored", lStudiedTriggers[itrig].Data()));
        
        //==============================================================================
        //        Create Plot2, please
        //==============================================================================
        TH1F *Plot2 = new TH1F(Form("Plot2_%s", lStudiedTriggers[itrig].Data()), "", 100, 0, 100);
        
        TH1F *Plot2_sub0 = (TH1F*) fHistCentralityKNO              ->Clone(Form("Plot2_sub0_%s", lStudiedTriggers[itrig].Data()));
        TH1F *Plot2_sub1 = (TH1F*) fHistCentralityKNO_Central      ->Clone(Form("Plot2_sub1_%s", lStudiedTriggers[itrig].Data()));
        TH1F *Plot2_sub2 = (TH1F*) fHistCentralityKNO_SemiCentral  ->Clone(Form("Plot2_sub2_%s", lStudiedTriggers[itrig].Data()));
        TH1F *Plot2_sub3 = (TH1F*) fHistCentralityKNO_Extra        ->Clone(Form("Plot2_sub3_%s", lStudiedTriggers[itrig].Data()));
        
        Plot2->SetStats(kFALSE);
        Plot2->GetYaxis()->SetTitle("Event Count");
        Plot2->GetXaxis()->SetTitle("Centrality (%)");
        
        Plot2_sub0 ->SetOption("same");
        Plot2_sub1 ->SetOption("same");
        Plot2_sub2 ->SetOption("same");
        Plot2_sub3 ->SetOption("same");
        
        //determine good ranges
        Double_t lMaxY = fHistCentralityKNO->GetMaximum();
        if ( fHistCentralityKNO_Central->GetMaximum()     > lMaxY ) lMaxY = fHistCentralityKNO_Central    ->GetMaximum();
        if ( fHistCentralityKNO_SemiCentral->GetMaximum() > lMaxY ) lMaxY = fHistCentralityKNO_SemiCentral->GetMaximum();
        if ( fHistCentralityKNO_Extra->GetMaximum()       > lMaxY ) lMaxY = fHistCentralityKNO_Extra      ->GetMaximum();
        
        //Don't add it to the outgoing directory - will be added to list of functions manually!
        Plot2_sub0 -> SetDirectory(0);
        Plot2_sub1 -> SetDirectory(0);
        Plot2_sub2 -> SetDirectory(0);
        Plot2_sub3 -> SetDirectory(0);
        
        Plot2->GetYaxis()->SetRangeUser(0, lMaxY*1.1);
        Plot2->SetFillStyle(0);
        
        Plot2->GetListOfFunctions()->Add(Plot2_sub0);
        Plot2->GetListOfFunctions()->Add(Plot2_sub1);
        Plot2->GetListOfFunctions()->Add(Plot2_sub2);
        Plot2->GetListOfFunctions()->Add(Plot2_sub3);
        
        TLine *lines[3];
        for(Int_t ilin=0; ilin<3; ilin++){
            lines[ilin] = new TLine(10+20*ilin, 0, 10+20*ilin, lMaxY*1.05);
            lines[ilin] -> SetLineStyle(7);
            lines[ilin] -> SetLineColor(kGray+1);
            Plot2->GetListOfFunctions()->Add(lines[ilin]);
        }
        
        TLegend *leg = new TLegend(.1, .91, .9, .99);
        //Process names, please
        leg->SetNColumns(4);
        leg->AddEntry(Plot2_sub0, GetTriggerFromTitle(histo->GetTitle()), "lp");
        leg->AddEntry(Plot2_sub1, GetTriggerFromTitle(histoc->GetTitle()), "lp");
        leg->AddEntry(Plot2_sub2, GetTriggerFromTitle(histosc->GetTitle()), "lp");
        leg->AddEntry(Plot2_sub3, GetTriggerFromTitle(histoextra->GetTitle()), "lp");
        
        Plot2->GetListOfFunctions()->Add(leg);
        
        cout<<"=================================================="<<endl;
        cout<<" Step 3: now doing anchoring via NBD-glauber fit"<<endl;
        cout<<"=================================================="<<endl;
        
        lOut -> Add(hClassifier);
        lOut -> Add(histo_ref);
        lOut -> Add(histoc_ref);
        lOut -> Add(histosc_ref);
        lOut -> Add(histoextra_ref);
        lOut -> Add(fHistCentrality);
        lOut -> Add(fHistCentrality_Central);
        lOut -> Add(fHistCentrality_SemiCentral);
        lOut -> Add(fHistCentrality_Extra);
        lOut -> Add(fThresholds);
        
        //KNO-scaling stuff
        lOut -> Add(f1_kno);
        
        lOut -> Add(histo_KNOanchored);
        lOut -> Add(hClassifierKNO);
        lOut -> Add(fHistCentralityKNO);
        lOut -> Add(fHistCentralityKNO_Central);
        lOut -> Add(fHistCentralityKNO_SemiCentral);
        lOut -> Add(fHistCentralityKNO_Extra);
        lOut -> Add(fThresholdsKNO);
        
        //Reference plots to draw
        lOut -> Add(Plot1);
        lOut -> Add(Plot2);
        
        //delete hReferenceV0M; hReferenceV0M=NULL;
        //delete hOnlineReferenceV0M; hOnlineReferenceV0M=NULL;
        
        
        
    }
    //Nice little table summary
    lOut -> Add(fThresholdTable);
    
    //List output being exported:
    cout<<"=================================================="<<endl;
    cout<<"Output TList at end of execution: "<<endl;
    lOut->ls();
    if( fFile ){
        lOut->Write("", TObject::kSingleKey);
        fFile->Close(); //if you close, you lose the ref histos!
                          //alternatively, do ::SetDirectory(0);
        //Cleanup test
        //lOut->Delete();
        //lIn->Delete();
        //file->Close();
        //fFile->Close();
    }
    
    cout<<" Done!"<<endl;
    return 0;
}


