#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCollection.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
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
const Double_t lReferenceY[] = {0, 0, 0.058459, 0.028266, 0.024931, 0.021396, 0.018341, 0.016266, 0.014523, 0.012981, 0.011771, 0.010827, 0.010053, 0.009338, 0.008873, 0.008492, 0.007921, 0.007565, 0.006979, 0.006767, 0.006628, 0.006383, 0.006174, 0.005948, 0.005681, 0.005546, 0.00527, 0.005257, 0.005023, 0.004839, 0.004854, 0.004709, 0.004507, 0.004409, 0.004205, 0.004165, 0.004172, 0.004067, 0.003909, 0.003934, 0.003794, 0.003817, 0.003595, 0.003665, 0.003508, 0.003512, 0.003436, 0.003363, 0.003395, 0.003314, 0.003277, 0.003105, 0.003152, 0.003137, 0.00315, 0.003049, 0.003027, 0.002941, 0.002939, 0.002862, 0.002846, 0.002764, 0.002713, 0.002698, 0.002674, 0.002595, 0.002566, 0.002586, 0.002482, 0.002511, 0.002462, 0.002524, 0.002416, 0.002422, 0.002352, 0.002419, 0.002371, 0.00229, 0.002357, 0.002354, 0.002267, 0.002211, 0.002243, 0.002201, 0.002231, 0.002189, 0.002109, 0.00206, 0.002112, 0.002001, 0.002086, 0.002142, 0.002029, 0.002041, 0.001987, 0.002003, 0.001986, 0.001937, 0.00191, 0.001908, 0.001995, 0.00187, 0.001886, 0.001881, 0.00189, 0.001865, 0.001755, 0.001837, 0.00183, 0.00182, 0.001809, 0.001717, 0.001778, 0.001769, 0.001676, 0.001648, 0.001746, 0.001705, 0.001677, 0.001675, 0.001624, 0.001576, 0.001627, 0.001546, 0.001594, 0.001589, 0.001558, 0.001615, 0.001622, 0.001594, 0.001624, 0.001472, 0.001573, 0.00152, 0.001546, 0.001509, 0.001436, 0.001503, 0.0015, 0.00146, 0.001536, 0.001497, 0.001477, 0.001462, 0.001467, 0.001443, 0.001395, 0.001423, 0.001376, 0.001434, 0.00139, 0.001381, 0.001364, 0.001383, 0.001375, 0.001352, 0.001371, 0.001389, 0.001312, 0.001308, 0.001322, 0.001338, 0.001327, 0.001322, 0.001313, 0.001373, 0.001276, 0.001304, 0.001248, 0.00127, 0.001241, 0.001299, 0.00124, 0.00123, 0.001233, 0.001207, 0.001231, 0.001205, 0.001228, 0.001156, 0.001227, 0.00117, 0.001198, 0.001196, 0.001181, 0.001174, 0.001117, 0.001188, 0.001176, 0.001115, 0.00119, 0.001174, 0.001154, 0.00114, 0.00117, 0.001157, 0.001152, 0.001111, 0.001126, 0.001165, 0.001092, 0.001098, 0.001082, 0.001134, 0.001062, 0.001133, 0.001075, 0.001089, 0.001122, 0.001071, 0.001096, 0.001046, 0.001062, 0.001074, 0.001014, 0.001009, 0.000994, 0.00108, 0.001084, 0.001044, 0.001019, 0.001033, 0.001061, 0.001049, 0.001015, 0.000997, 0.001021, 0.000994, 0.001029, 0.001006, 0.000991, 0.000981, 0.001013, 0.000963, 0.001031, 0.001002, 0.000907, 0.00094, 0.000936, 0.000985, 0.001019, 0.000911, 0.000922, 0.000964, 0.000923, 0.000953, 0.000992, 0.000894, 0.000971, 0.000886, 0.000932, 0.000912, 0.000935, 0.000884, 0.000872, 0.000909, 0.000932, 0.000882, 0.000809, 0.000929, 0.000894, 0.000912, 0.000901, 0.000912, 0.00088, 0.000894, 0.00089, 0.000887, 0.000875, 0.000887, 0.000894, 0.000865, 0.000895, 0.000869, 0.000854, 0.000841, 0.000835, 0.000854, 0.000865, 0.00083, 0.000875, 0.000802, 0.000864, 0.000847, 0.000867, 0.000815, 0.0008, 0.0009, 0.000824, 0.000861, 0.000789, 0.000825, 0.000799, 0.000785, 0.000809, 0.000828, 0.000776, 0.000808, 0.000788, 0.000796, 0.000805, 0.000826, 0.00076, 0.000787, 0.000813, 0.000802, 0.000765, 0.000862, 0.000795, 0.000797, 0.000798, 0.000758, 0.00082, 0.000751, 0.000723, 0.00078, 0.000781, 0.000755, 0.000729, 0.000749, 0.000776, 0.000711, 0.000754, 0.00071, 0.000757, 0.000737, 0.00074, 0.00073, 0.000741, 0.000748, 0.000663, 0.000723, 0.000742, 0.000729, 0.000772, 0.000708, 0.000729, 0.000747, 0.000707, 0.000756, 0.000709, 0.000706, 0.000722, 0.000758, 0.000721, 0.000711, 0.000706, 0.000681, 0.000693, 0.000732, 0.000701, 0.000714, 0.00066, 0.000731, 0.000715, 0.000685, 0.000678, 0.000678, 0.00065, 0.000643, 0.000671, 0.00065, 0.000678, 0.000645, 0.000651, 0.000642, 0.000691, 0.000701, 0.000728, 0.000667, 0.000614, 0.000658, 0.000648, 0.000623, 0.000662, 0.000685, 0.000659, 0.000653, 0.000625, 0.000672, 0.000663, 0.000657, 0.000648, 0.000636, 0.000642, 0.000619, 0.000618, 0.000608, 0.000628, 0.000623, 0.000652, 0.000658, 0.000616, 0.00063, 0.000602, 0.000646, 0.000607, 0.000652, 0.000595, 0.000656, 0.000608, 0.000626, 0.000598, 0.000654, 0.00065, 0.000583, 0.000619, 0.000588, 0.000577, 0.000591, 0.00064, 0.000605, 0.000592, 0.000597, 0.000589, 0.000601, 0.000615, 0.000589, 0.000561, 0.000581, 0.000574, 0.000573, 0.000577, 0.000569, 0.000574, 0.000579, 0.000571, 0.000576, 0.000578, 0.000614, 0.000592, 0.000562, 0.000579, 0.00061, 0.000606, 0.000577, 0.000562, 0.000566, 0.000596, 0.000588, 0.000568, 0.000561, 0.000569, 0.000549, 0.000558, 0.000564, 0.000538, 0.000558, 0.000595, 0.000528, 0.000534, 0.000531, 0.000541, 0.000557, 0.000517, 0.00051, 0.000547, 0.000556, 0.000573, 0.000578, 0.000498, 0.000568, 0.000558, 0.000507, 0.00056, 0.000529, 0.000546, 0.000547, 0.000493, 0.000524, 0.000528, 0.000529, 0.000531, 0.000482, 0.000589, 0.000517, 0.00049, 0.000523, 0.000518, 0.000517, 0.000531, 0.000519, 0.000463, 0.000527, 0.000513, 0.000549, 0.000532, 0.0005, 0.000517, 0.000541, 0.000523, 0.000524, 0.000481, 0.000541, 0.000504, 0.000511, 0.000524, 0.000552, 0.000523, 0.000498, 0.000484, 0.000493, 0.00048, 0.000508, 0.000485, 0.000538, 0.000479, 0.000494, 0.000499, 0.000509, 0.000512, 0.000489, 0.000493, 0.000476, 0.000502, 0.000494, 0.000499, 0.000474, 0.000505, 0.000489, 0.000487, 0.000473, 0.000473, 0.000466, 0.000476, 0.000484, 0.000414, 0.000455, 0.000485, 0.000495, 0.000465, 0.000448, 0.000484, 0.000434, 0.000469, 0.000498, 0.000462, 0.000425, 0.000473, 0.000454, 0.00044, 0.000479, 0.000466, 0.000484, 0.000484, 0.000459, 0.000441, 0.000468, 0.00048, 0.000433, 0.000432, 0.000445, 0.000465, 0.000452, 0.000436, 0.000437, 0.000464, 0.000443, 0.000457, 0.000471, 0.000397, 0.000458, 0.000422, 0.000457, 0.000461, 0.000456, 0.00041, 0.000442, 0.000444, 0.000414, 0.000434, 0.000422, 0.000453, 0.000432, 0.000487, 0.000449, 0.000455, 0.00042, 0.00045, 0.000464, 0.000415, 0.000419, 0.000409, 0.000431, 0.000428, 0.000447, 0.000398, 0.000445, 0.000457, 0.000409, 0.000437, 0.000442, 0.00044, 0.000421, 0.000412, 0.0004, 0.000422, 0.000396, 0.000418, 0.000439, 0.000402, 0.000432, 0.000402, 0.000416, 0.000416, 0.000385, 0.000418, 0.000439, 0.00043, 0.000418, 0.000417, 0.000413, 0.000373, 0.000403, 0.000409, 0.000391, 0.000371, 0.00041, 0.000387, 0.000402, 0.000396, 0.000381, 0.000419, 0.000388, 0.000419, 0.00039, 0.000386, 0.000412, 0.000466, 0.000419, 0.000414, 0.000381, 0.000385, 0.000412, 0.000401, 0.00039, 0.000391, 0.000363, 0.000405, 0.000408, 0.000406, 0.000405, 0.000374, 0.000386, 0.00037, 0.000403, 0.00038, 0.000373, 0.000372, 0.000403, 0.000351, 0.000397, 0.00038, 0.000369, 0.000405, 0.000378, 0.000384, 0.000395, 0.000379, 0.000371, 0.000385, 0.000379, 0.000394, 0.000385, 0.000372, 0.000374, 0.000363, 0.000376, 0.00038, 0.000403, 0.00037, 0.000388, 0.0004, 0.000343, 0.000344, 0.000363, 0.000383, 0.00037, 0.000362, 0.00037, 0.000375, 0.000357, 0.000366, 0.000371, 0.000383, 0.000373, 0.000395, 0.000362, 0.000419, 0.000391, 0.000372, 0.000353, 0.000346, 0.000367, 0.000362, 0.000346, 0.000345, 0.000356, 0.000375, 0.000351, 0.000354, 0.000401, 0.000341, 0.000385, 0.000336, 0.000351, 0.000322, 0.000304, 0.000317, 0.000309, 0.000302, 0.000297, 0.000301, 0.000306, 0.000302, 0.000275, 0.000295, 0.000253, 0.000258, 0.000258, 0.000264, 0.000237, 0.000262, 0.000241, 0.00022, 0.000212, 0.000218, 0.000186, 0.000176, 0.000164, 0.00017, 0.000151, 0.000158, 0.000149, 0.000138, 0.000133, 0.000113, 0.000102, 9.8e-05, 8.7e-05, 7.6e-05, 6.7e-05, 7.1e-05, 6e-05, 6.8e-05, 7.5e-05, 6.1e-05, 5e-05, 4.7e-05, 3.8e-05, 3.8e-05, 2.7e-05, 2.9e-05, 2.5e-05, 3.3e-05, 1.9e-05, 1.6e-05, 2.6e-05, 1.6e-05, 1e-05, 1.3e-05, 1.1e-05, 1e-05, 8e-06, 1e-06, 7e-06, 6e-06, 5e-06, 4e-06, 5e-06, 2e-06, 4e-06, 4e-06, 3e-06, 0, 0, 2e-06, 2e-06, 1e-06, 0, 1e-06, 1e-06, 0, 0, 1e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

TH1F *hReferenceV0M(NULL);

//Definition of a negative binomial distribution
//(used in mult fitter function)
TF1 *fNBD(NULL);

//______________________________________________________
Double_t pdf_PbPb_kno(Double_t *x, Double_t *par)
//Calculation of probability of a certain amplitude for Pb-Pb
//based on simple X/Y axis scaling and known distribution
{
    Double_t lSmallNumber = 1e-12;
    Double_t lMultValue = TMath::Floor(x[0]+0.5);
    Double_t lProbability = 0.0;
    
    //   hNanc->Draw();
    
    for(Long_t iNanc = 1; iNanc<999; iNanc++){
        Double_t lThisMu = ((Double_t)iNanc)*par[0];
        Double_t lThisk = ((Double_t)iNanc)*par[1];
        Double_t lpval = TMath::Power(1+lThisMu/lThisk,-1);
        fNBD->SetParameter(1,lThisk);
        fNBD->SetParameter(0,lpval);
        Double_t lMult = fNBD->Eval(lMultValue);
        lProbability += hReferenceV0M->GetBinContent(hReferenceV0M->FindBin(iNanc))*lMult;
    }
    //don't forget Y axis scaling factor, please
    return par[2]*lProbability;
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
    TH2F *fThresholdTable = new TH2F("fThresholdTable", "Threshhold table", 3, 0, 3, 4, 0, 4);
    fThresholdTable->GetXaxis()->SetBinLabel(1, "50%");
    fThresholdTable->GetXaxis()->SetBinLabel(2, "30%");
    fThresholdTable->GetXaxis()->SetBinLabel(3, "10%");
    
    fThresholdTable->SetMarkerSize(2.5);
    fThresholdTable->SetOption("text0");
    fThresholdTable->SetStats(kFALSE);
    
    //Prepare reference histogram in case it does not exist
    if(!hReferenceV0M){
        hReferenceV0M = new TH1F("fHistReferenceV0M", "", 1000, -0.5, 999.5);
        for(Int_t i=1; i<900; i++){
            //cout<<"Reference i = "<<i<<", "<<lReferenceY[i]<<endl;
            hReferenceV0M->SetBinContent(i, lReferenceY[i]);
        }
        hReferenceV0M->SetDirectory(0);
    }
    if(!fNBD){
        fNBD = new TF1("NBD","ROOT::Math::negative_binomial_pdf(x,[0],[1])",0,800);
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
        
        const Int_t lCentralityBins = 200;
        
        Float_t lCentBounds[lCentralityBins+1];
        Float_t lCentBoundsRaw[lCentralityBins+1];
        
        cout<<"=================================================="<<endl;
        cout<<" Step 1: now doing unanchored calibration of "<<lStudiedTriggers[itrig]<<endl;
        cout<<"=================================================="<<endl;
        
        cout<<"Determining calibration on-the-fly..."<<endl;
        lCentBoundsRaw[0] = 0.0;
        cout<<"Inspect boundaries: "<<flush;
        for(Int_t i=0; i<lCentralityBins; i++){
            lCentBounds[i] = i*(100./((Double_t)(lCentralityBins)));
            if (i!=0) lCentBoundsRaw[i] = GetBoundaryForPercentile( histo,  100-lCentBounds[i] );
            cout<<lCentBoundsRaw[i]<<" "<<flush;
        }
        lCentBoundsRaw[lCentralityBins] = 1e+5;
        cout<<lCentBoundsRaw[lCentralityBins]<<" (end)"<<endl;
        
        cout<<"Classifying..."<<endl;
        TH1F* fHistCentrality             = new TH1F(Form("fHistCentrality%s",lStudiedTriggers[itrig].Data()), "", lCentralityBins, 0, 100);
        TH1F* fHistCentrality_Central     = new TH1F(Form("fHistCentrality%s_Central",lStudiedTriggers[itrig].Data()), "", lCentralityBins, 0, 100);
        TH1F* fHistCentrality_SemiCentral = new TH1F(Form("fHistCentrality%s_SemiCentral",lStudiedTriggers[itrig].Data()), "", lCentralityBins, 0, 100);
        TH1F* fHistCentrality_Extra       = new TH1F(Form("fHistCentrality%s_Extra",lStudiedTriggers[itrig].Data()), "", lCentralityBins, 0, 100);
        
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
        
        TH1F *hClassifier = new TH1F(Form("hClassifier_%s", lStudiedTriggers[itrig].Data()), "", lCentralityBins, lCentBoundsRaw);
        for(Int_t i=0; i<lCentralityBins+1; i++){
            hClassifier -> SetBinContent( i+1 , 100.-100./((Double_t)(lCentralityBins))*0.5 - 100./((Double_t)(lCentralityBins))*i );
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
        f1_kno = new TF1(Form("f1_kno_%s",lStudiedTriggers[itrig].Data()),pdf_PbPb_kno,300,50000,3);
        
        
        f1_kno->SetTitle(Form("ref fit to %s signal (raw);%s signal (raw);f()", lStudiedTriggers[itrig].Data(), lStudiedTriggers[itrig].Data()));
        
        cout<<"Guessing normalization to: "<<histo->GetEntries()<<endl;
        
        f1_kno->SetParameter(0,45);
        f1_kno->FixParameter(1,1.5);
        f1_kno->SetParameter(2,1.1e+6);
        
        f1_kno->SetParName(0, "mu");
        f1_kno->SetParName(1, "k");
        f1_kno->SetParName(2, "norm");
        f1_kno->SetNpx(100);
        
        cout<<"Cloning histo... ";
        TH1F *Plot1 = (TH1F*) histo->Clone(Form("Plot1_%s", lStudiedTriggers[itrig].Data()));
        Plot1->Rebin(100);
        
        //Calculate Y axis scaling factor @ 15K
        Double_t lGuess =f1_kno->Eval(15000) / Plot1->GetBinContent( Plot1->FindBin( 15000) );
        cout<<"Y axis scale = "<<f1_kno->Eval(15000) / Plot1->GetBinContent( Plot1->FindBin( 15000) ) <<endl;
        f1_kno->SetParameter(2,1.1e+6/lGuess);
        
        cout<<"Fitting... "<<endl;
        Plot1->Fit(Form("f1_kno_%s",lStudiedTriggers[itrig].Data()),"R0");
        
        //Create a clone for drawing purposes
        cout<<"Cloning f1... ";
        TF1 *f1_kno_drawclone = (TF1*) f1_kno->Clone( Form("f1_kno_drawclone_%s",lStudiedTriggers[itrig].Data()) );
        f1_kno_drawclone->SetLineColor(kRed+1);
        
        Plot1->GetListOfFunctions()->Add(f1_kno_drawclone);
        Plot1->SetStats(kFALSE);
        
        //Reference mu from 2015 :
        // 4.16229e+01
        
        
        Double_t lAnchorPointKNO = f1_kno->GetParameter(0)*129.59/4.16229e+01;
        
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
        TH1F* fHistCentralityKNO             = new TH1F(Form("fHistCentralityKNO%s",lStudiedTriggers[itrig].Data()), "", lCentralityBins, 0, 100);
        TH1F* fHistCentralityKNO_Central     = new TH1F(Form("fHistCentralityKNO%s_Central",lStudiedTriggers[itrig].Data()), "", lCentralityBins, 0, 100);
        TH1F* fHistCentralityKNO_SemiCentral = new TH1F(Form("fHistCentralityKNO%s_SemiCentral",lStudiedTriggers[itrig].Data()), "", lCentralityBins, 0, 100);
        TH1F* fHistCentralityKNO_Extra       = new TH1F(Form("fHistCentralityKNO%s_Extra",lStudiedTriggers[itrig].Data()), "", lCentralityBins, 0, 100);
        
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
        Float_t lCentBoundsRawKNO[lCentralityBins+1];
        lCentBoundsRawKNO[0] = 0.0;
        cout<<"Inspect boundaries: "<<flush;
        for(Int_t i=0; i<lCentralityBins; i++){
            //lCentBounds[i] = i;
            if (i!=0) lCentBoundsRawKNO[i] = GetBoundaryForPercentile( histo_KNOanchored,  100-lCentBounds[i] );
            cout<<lCentBoundsRawKNO[i]<<" "<<flush;
        }
        lCentBoundsRawKNO[lCentralityBins] = 1e+5;
        cout<<lCentBoundsRawKNO[lCentralityBins]<<" (end)"<<endl;
        
        TH1F *hClassifierKNO = new TH1F(Form("hClassifierKNO_%s",lStudiedTriggers[itrig].Data()), "", lCentralityBins, lCentBoundsRawKNO);
        for(Int_t i=0; i<lCentralityBins+1; i++){
            hClassifierKNO -> SetBinContent( i+1 , 100.-100./((Double_t)(lCentralityBins))*0.5 - 100./((Double_t)(lCentralityBins))*i );
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
        Plot2_sub3 ->SetOption("same");
        Plot2_sub1 ->SetOption("same");
        Plot2_sub2 ->SetOption("same");
        
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
        
        //Size options...
        Plot2_sub0 ->SetMarkerSize(0.33);
        Plot2_sub1 ->SetMarkerSize(0.33);
        Plot2_sub2 ->SetMarkerSize(0.33);
        Plot2_sub3 ->SetMarkerSize(0.33);
        
        Plot2->GetListOfFunctions()->Add(Plot2_sub0);
        Plot2->GetListOfFunctions()->Add(Plot2_sub3);
        Plot2->GetListOfFunctions()->Add(Plot2_sub1);
        Plot2->GetListOfFunctions()->Add(Plot2_sub2);
        
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


