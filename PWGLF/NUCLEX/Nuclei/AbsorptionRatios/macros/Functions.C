
Double_t background(Double_t *x, Double_t *par) {
    Float_t xval = x[0];
    Float_t par0 = par[0];
    Float_t par1 = par[1];
    Float_t par2 = par[2];
    Float_t par3 = par[3];
    return par0 + par1*xval + par2*TMath::Exp(xval*par3);
}

Double_t backGrProton(Double_t *x, Double_t *par) {
    Float_t xval = x[0];
    Float_t par0 = par[0];
    Float_t par1 = par[1];
    Float_t par2 = par[2];
    Float_t par3 = par[3];
    Float_t par4 = par[4];
    Float_t par5 = par[5];
    Float_t par6 = par[6];
    if(xval <= (par4 + par5*par6)){
        return TMath::Exp(par0+par1*xval+par2*xval*xval) + par3*TMath::Gaus(xval, par4, par5);
    }else{
        return TMath::Exp(par0+par1*xval+par2*xval*xval) + par3*TMath::Exp(-(xval-par4-par6*par5*0.5)*par6/par5);
    }
    
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t gaussian(Double_t *x, Double_t *par) {
    
    Float_t xval = x[0];
    Float_t par0 = par[0];
    Float_t par1 = par[1];
    Float_t par2 = par[2];
    return par0*TMath::Gaus(xval,par1,par2);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t gaussianSignal(Double_t *x, Double_t *par) {
    
    Float_t xval = x[0];
    Float_t par0 = par[0];
    Float_t par1 = par[1];
    Float_t par2 = par[2];
    Float_t par3 = par[3];
    if(xval < par1+par3*par2){
        return par0*1/(TMath::Sqrt(0.5*TMath::Pi())*(par2+par2*TMath::Erf(par3/TMath::Sqrt(2)))+par2/par3*TMath::Exp(-par3*par3/2))*TMath::Exp(-(xval-par1)*(xval-par1)/2/par2/par2);
    }else{
        return par0*1/(TMath::Sqrt(0.5*TMath::Pi())*(par2+par2*TMath::Erf(par3/TMath::Sqrt(2)))+par2/par3*TMath::Exp(-par3*par3/2))*TMath::Exp(-(xval-par1-par3*par2*0.5)*par3/par2);
    }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Double_t MassFitFunction(Double_t *x, Double_t *par) {
    return gaussianSignal(x,par) + background(x,&par[4]);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t MassFitFunctionProton(Double_t *x, Double_t *par) {
    return gaussianSignal(x,par) + background(x,&par[4]);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Double_t DCAGauss(Double_t *x, Double_t *par) {
    return gaussian(x,par) + gaussian(x,&par[3]);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t DCAQuadraticFunc(Double_t *x, Double_t *par){
    return par[0]*x[0]*x[0] + par[1]*x[0] + par[2];
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Double_t DCAFitFunction(Double_t *x, Double_t *par){
    return DCAQuadraticFunc(x,par) + DCAGauss(x,&par[3]);
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
