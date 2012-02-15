//------------------------------------------------------------------------------
// divide.C
//
// macro to divide two TGraphs using Exponential interpolation
//------------------------------------------------------------------------------

using namespace std;

Double_t ExtrapolateExp(Double_t x0, Double_t y0, Double_t x1, Double_t y1, Double_t x)
{
        // extra/interpolates using an exponential defined by (x0,y0) and (x1,y1) to x

        Double_t b = (TMath::Log(y0) - TMath::Log(y1)) / (x0 - x1);
        Double_t a = TMath::Log(y0) - b * x0;
        Double_t a2 = TMath::Log(y1) - b * x1;

        //Printf("%f %f %f", b, a, a2);

        return TMath::Exp(a + b * x);
} 


TGraphErrors* divide(TGraph* dividend, TGraph* divisor, Bool_t exp = kTRUE)
{

cout << "---------------------------------------------------------" << endl;
cout << "starting macro: << divide.C >> " << endl;
cout << "---------------------------------------------------------" << endl;

Double_t x=0;
Double_t y=0;
Double_t y1;
Double_t y2;
Double_t ex;
Double_t ey;
Double_t e1;
Double_t e2;


divisor->GetPoint(0,x,y);
Double_t minx = x;
cout << "divisor: minx = " << minx << endl;

divisor->GetPoint(divisor->GetN()-1,x,y);
Double_t maxx = x;
cout << "divisor: maxx = " << maxx << endl;

dividend->GetPoint(0,x,y);
if (x > minx) { minx = x; } 

dividend->GetPoint(dividend->GetN()-1,x,y);
if (x < maxx) { maxx = x; } 

cout << "minx = " << minx << endl;
cout << "maxx = " << maxx << endl;

Int_t bins = 0;

for (int i=0; i < dividend->GetN(); i++) {
    dividend->GetPoint(i,x,y);
    if ((x <= maxx) && (x >= minx)) {bins++;}
}

TGraphErrors* result = new TGraphErrors(bins);
TGraph* errors = new TGraph(divisor->GetN(),divisor->GetX(),divisor->GetEY());

Int_t j=0;

for (int i=0; i < dividend->GetN(); i++) {
    dividend->GetPoint(i,x,y1);
    if (!((x <= maxx) && (x >= minx))) { continue; }            
    if (exp) {
        Double_t x1temp = 0;
        Double_t y1temp = 0;
        for (int n=0; n < divisor->GetN(); n++) {
            divisor->GetPoint(n,x1temp,y1temp);
            
            if (x1temp >= x)  break;
        }
        cout << "x=" << x;
        cout << "x1temp=" << x1temp;
        if (x1temp == x) { y2 = y1temp; }
        if (x1temp > x) { 
            Double_t x0temp = 0;
            Double_t y0temp = 0;
            cout << "nfortemp=" << n << endl;
            divisor->GetPoint(n-1,x0temp,y0temp);
            y2 = ExtrapolateExp(x0temp,y0temp,x1temp,y1temp,x);
            cout << "exp extrapolation used" << endl;
        }
    } else {
        y2 = divisor->Eval(x);
    }
    
    y = y1/y2;
    e1 = dividend->GetErrorY(i);
    //e2 = errors->Eval(x);
    e2 = 0;    
    ex = 0;
    ey = y*sqrt((e1/y1)*(e1/y1) + (e2/y2)*(e2/y2));    
    result->SetPoint(j,x,y);
    result->SetPointError(j,ex,ey);
    j++;
    cout << "data point # " << j << endl;
    cout << "   x= " << x << " +- " << ex <<endl;
    cout << "   y= " << y << " +- " << ey <<endl;
}
delete errors;
errors = 0;



cout << "min:" << minx <<endl;
cout << "max:" << maxx <<endl;

cout << "---------------------------------------------------------" << endl;
cout << "finished macro: << divide.C >> " << endl;
cout << "---------------------------------------------------------" << endl;
return result;
}