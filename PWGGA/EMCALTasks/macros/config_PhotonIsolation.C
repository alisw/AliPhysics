vector<Double_t> makeLinBinning(Double_t min, Double_t max, Int_t num){
	if (max<min){
		Double_t tmp=min;
		min=max;
		max=tmp;
	}
	Double_t width=(max-min)/num;
	vector<Double_t> binning;

	for (Int_t i =0; i<num+1; i++){
		binning.push_back(min+width*(Double_t)i);
	}

	return binning;
}
vector<Double_t> makeLogBinning(Double_t min, Double_t max, Int_t num){
	if (min<1e-20 || max<1e-20){
		AliError("For Log binning min and max must be > 1e-20. Using linear binning instead!");
		return makeLinBinning(num, min, max);
	}
	if (max<min){
		Double_t tmp=min;
		min=max;
		max=tmp;
	}
	Double_t expMax=TMath::Log(max/min);
	vector<Double_t> binning;

	for (Int_t i =0; i<num+1; i++){
		binning.push_back(min*TMath::Exp(expMax/num*(Double_t)i));
	}

	return binning;
}

vector<Double_t> makeArbBinning(Double_t edges[], Int_t num){

	vector<Double_t> binning;
	for(Int_t i =0; i<num; i++){
		binning.push_back(edges[i]);
	}
	return binning;
}

Double_t PID_bins[102] = {-2300.5,-2212.5,-2211.5,-2112.5,-2111.5,-2011.5,-1911.5,-1811.5,-1711.5,-1611.5,-1511.5,-1411.5,-1311.5,-1211.5,-1111.5,-1011.5,-911.5,-811.5,-711.5,-611.5,-511.5,-411.5,-321.5,-320.5,-311.5,-310.5,-309.5,-221.5,-220.5,-211.5,-210.5,-130.5,-129.5,-113.5,-111.5,-110.5,-109.5,-99.5,-49.5,-38.5,-31.5,-22.5,-21.5,-20.5,-18.5,-16.5,-14.5,-12.5,-11.5,-10.5,-0.5,0.5,10.5,11.5,12.5,14.5,16.5,18.5,20.5,21.5,22.5,31.5,38.5,49.5,99.5,109.5,110.5,111.5,113.5,129.5,130.5,210.5,211.5,220.5,221.5,309.5,310.5,311.5,320.5,321.5,411.5,511.5,611.5,711.5,811.5,911.5,1011.5,1111.5,1211.5,1311.5,1411.5,1511.5,1611.5,1711.5,1811.5,1911.5,2011.5,2111.5,2112.5,2211.5,2212.5,2300.5};

vector<Double_t> ptBin = makeLinBinning(0.,70.,70);
vector<Double_t> M02Bin = makeLinBinning(0.,2.,200);
vector<Double_t> EtisoBin = makeLinBinning(-10.,100.,110);
vector<Double_t> EtueBin = makeLinBinning(-10.,100.,110);
vector<Double_t> EtaBin = makeLinBinning(-0.7.,0.7,70);
vector<Double_t> PhiBin = makeLinBinning(1.3,3.3,50);
vector<Double_t> LabelBin = makeLinBinning(0.,1500.,150);
vector<Double_t> PDGBin = makeArbBinning(PID_bins,sizeof(PID_bins)/sizeof(PID_bins[0]));
vector<Double_t> MomPDGBin = makeArbBinning(PID_bins,sizeof(PID_bins)/sizeof(PID_bins[0]));
vector<Double_t> ClustPDGBin = makeArbBinning(PID_bins,sizeof(PID_bins)/sizeof(PID_bins[0]));
vector<Double_t> DxBin = makeLinBinning(-1.,1.,2);
vector<Double_t> DzBin = makeLinBinning(-1.,1.,2);
vector<Double_t> DecayBin = makeLinBinning(0.,10.,10);
