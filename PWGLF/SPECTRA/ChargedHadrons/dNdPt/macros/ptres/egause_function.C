// ExpoGaussExpo function
// to use:

/*
.L egause_function.C++g
Double_t xmin = 0;
Double_t xmax = 1;
TF1 *egause = new TF1("egause", egause_function, xmin, xmax, 5, 1);
egause->SetParNames("Constant", "Mean", "Sigma", "ExpoLow", "ExpoHigh");
*/


Double_t egause_function(Double_t* x, Double_t *p) {
    //p[0] normalization
    //p[1] gaus mean
    //p[2] gaus sigma
    //p[3] lower exp tail
    //p[4] high exp tail
    Double_t d = (x[0]-p[1])/p[2];
    if (d  < -p[3])  return p[0]*exp(0.5*p[3]*p[3]+p[3]*d);
    if (d  >  p[4])  return p[0]*exp(0.5*p[4]*p[4]-p[4]*d);
    return p[0]*exp(-0.5*d*d);
    
}
