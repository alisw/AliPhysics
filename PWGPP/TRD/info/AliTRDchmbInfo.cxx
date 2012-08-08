////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Chamber Info Incapsulation                                             //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////



#include "TBox.h"

#include "AliLog.h"
#include "AliTRDchmbInfo.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDchmbInfo)
//____________________________________________
AliTRDchmbInfo::AliTRDchmbInfo()
  :TNamed()
  ,fDet(-1)
  ,fStatus(0)
  ,fBox(NULL)
{
//  Constructor. Reset all fields.
  memset(fPosition, 0, 4*sizeof(Double_t));
}


//____________________________________________
AliTRDchmbInfo::AliTRDchmbInfo(Int_t det, Int_t stat, Double_t p[4])
  :TNamed(Form("D%03d", det), "")
  ,fDet(det)
  ,fStatus(stat)
  ,fBox(NULL)
{
//  Constructor. Set position fields.
  SetPosition(p);
}

//____________________________________________
AliTRDchmbInfo::~AliTRDchmbInfo()
{
//  Destructor. 
  if(fBox) delete fBox;
}

//____________________________________________
void AliTRDchmbInfo::Print(Option_t */*o*/) const
{
//   Dump formated chamber info

  printf("  DET[%03d] Status[%d]\n", fDet, fStatus);
  printf("  eta[%f %f] phi[%f %f]\n", fPosition[0], fPosition[2], fPosition[1], fPosition[3]);
}

//____________________________________________
void AliTRDchmbInfo::SetDetector(Int_t det)
{
// register detector no.
  fDet = det;
  SetName(Form("D%03d", det));
}

//____________________________________________
void AliTRDchmbInfo::Draw(Option_t *opt)
{
// set style parameters according to chamber status

  Int_t style[] = {0, 3003, 3008};
  if(!fBox) fBox = new TBox();
  fBox->SetLineColor(kBlack);fBox->SetFillColor(kWhite);fBox->SetFillStyle(style[0]);

  Float_t xmin(0.), xmax(0.);
  if(strcmp(opt, "eta")==0){
    xmin=fPosition[0]; xmax=fPosition[2];
  } else if(strcmp(opt, "pad")==0){
    Int_t stk(AliTRDgeometry::GetStack(fDet));
    xmin=-0.6+16*(4-stk)-(stk<2?4:0); xmax=xmin+(stk==2?12:16)-0.2;
  } else {
    AliWarning(Form("Unknow option \"%s\"", opt));
    return;
  }

  AliDebug(2, Form("D[%03d] 0[%+4.1f(%+4.1f) %+4.1f] 1[%+4.1f(%+4.1f) %+4.1f] opt[%d]", fDet, xmin, fPosition[0], fPosition[1], xmax, fPosition[2], fPosition[3], fStatus));
  fBox->SetX1(xmin); fBox->SetX2(xmax); fBox->SetY1(fPosition[1]); fBox->SetY2(fPosition[3]);
  if(fStatus==1){
    fBox->SetFillStyle(style[1]);fBox->SetFillColor(kBlack);
  } else {
    if(fStatus==2){
      TBox *b = fBox->DrawBox(xmin, fPosition[1], xmax, 0.5*(fPosition[3]+fPosition[1]));
      b->SetFillStyle(style[2]);b->SetFillColor(kGreen);
    } else if(fStatus==3){
      TBox *b = fBox->DrawBox(xmin, 0.5*(fPosition[3]+fPosition[1]), xmax, fPosition[3]);
      b->SetFillStyle(style[2]);b->SetFillColor(kRed);
    } else if(fStatus!=0){
      AliError(Form("Wrong chmb. status[%d] for det[%03d]", fStatus, fDet));
      return;
    }
  }
  fBox->Draw();
}
