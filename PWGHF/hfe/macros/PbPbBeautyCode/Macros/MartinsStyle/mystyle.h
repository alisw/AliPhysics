using namespace std;
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"
#include "TLine.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TGraph.h"
#include "iostream"
#include "TF1.h"
#include "TAxis.h"

void  mystyle()
{
    // to compile the lib do once
    // .L mystyle.C++
    // if lib has been compiled load the lib
    // add this lines to your macro
    // gSystem->Load("/misc/pachmay/bin/mystyle_C.so");
    // mystyle();


    cout<<"LOADING STYLE SETTING!"<<endl;

    Float_t titlesize=0.06;
    Float_t titleoffset=1.25;
    Float_t titleoffsety=1.5;
    Float_t labelsize=0.05;

    // defining canvas options ------------------------
    gStyle->SetCanvasColor(00);
    gStyle->SetPadColor(00);

    // defining size of pad ---------------------------
    gStyle->SetPadTopMargin(0.03);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadLeftMargin(0.12);

    gStyle->SetPadBorderMode(0);

    // defining outfit of axis ------------------------
    // offset of titles
    gStyle->SetTitleXOffset(titleoffset);
    gStyle->SetTitleYOffset(titleoffsety);
    // size of titles
    gStyle->SetTitleSize(titlesize,"X");
    gStyle->SetTitleSize(titlesize,"Y");
    gStyle->SetTitleSize(titlesize,"Z");
    // size of axis labels
    gStyle->SetLabelSize(labelsize,"X");
    gStyle->SetLabelSize(labelsize,"Y");
    gStyle->SetLabelSize(labelsize,"Z");

    // switch options of the histogram ----------------
    gStyle->SetOptStat(0);         // switch off stat box
    gStyle->SetErrorX(0);          // switch off error in x
    gStyle->SetTitleH(0.08);       // make title larger
    //gStyle->SetTitleW(0.8);       // make title larger
    gStyle->SetTitleStyle(0);      // make title box transparent
    gStyle->SetTitleBorderSize(0); // make border of box disapear
}

void setOPT_canvas(TCanvas *canvas)
{

    gStyle->SetLineStyleString(22,"80 18 12 18 12 12");
    gStyle->SetPadColor(0);
    /*
    canvas->SetLeftMargin(0.143005);
    canvas->SetRightMargin(0.0559586);
    canvas->SetTopMargin(0.0695266);
    canvas->SetBottomMargin(0.130178);

    canvas->Range(-173.575,-8.80888,1076.42,-1.15554);
    */
    canvas->SetFillColor(0);
    canvas->SetBorderMode(0);
    canvas->SetBorderSize(0);
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetFrameFillColor(0);
    canvas->SetFrameLineWidth(2);
    canvas->SetFrameBorderMode(0);
    canvas->SetFrameBorderSize(0);
 

}
void setOPT_hists(TH1F *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510, Int_t marker_style =20, Float_t marker_size =1.3, Int_t color=1)
{
    gStyle->SetPadRightMargin(0.05);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetXaxis()->SetTitleOffset(0.8);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
}

void setOPT_hists2D(TH2F *hist2, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510)
{
    gStyle->SetPadRightMargin(0.15);
    hist2->GetXaxis()->SetTitle(xAxisTitle);
    hist2->GetYaxis()->SetTitle(yAxisTitle);
    hist2->GetXaxis()->SetTitleSize(0.06);
    hist2->GetYaxis()->SetTitleSize(0.06);
    hist2->GetZaxis()->SetTitleSize(0.06);
    hist2->GetXaxis()->SetTitleFont(42);
    hist2->GetYaxis()->SetTitleFont(42);
    hist2->GetZaxis()->SetTitleFont(42);
    hist2->GetXaxis()->SetNdivisions(Ndevision);
    hist2->GetYaxis()->SetTitleOffset(1.2);
//    hist2->GetXaxis()->CenterTitle();
//    hist2->GetYaxis()->CenterTitle();
    hist2->GetXaxis()->SetLabelFont(42);
    hist2->GetYaxis()->SetLabelFont(42);
    hist2->GetZaxis()->SetLabelFont(42);
    hist2->GetXaxis()->SetLabelSize(0.05);
    hist2->GetYaxis()->SetLabelSize(0.05);
    hist2->GetZaxis()->SetLabelSize(0.05);
}


void setOPT_histsF(TF1 *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndevision=510, Int_t marker_style =20, Float_t marker_size =1.3, Int_t color=1)
{
    gStyle->SetPadRightMargin(0.05);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndevision);
    hist->GetYaxis()->SetTitleOffset(1.4);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
}

void setOPT_graph(TGraph *hist, TString xAxisTitle=" ", TString yAxisTitle=" ", Int_t Ndivision=510, Int_t marker_style =20, Float_t marker_size =1.3, Int_t color=1)
{
    gStyle->SetPadRightMargin(0.05);
    hist->GetXaxis()->SetTitle(xAxisTitle);
    hist->GetYaxis()->SetTitle(yAxisTitle);
    hist->GetXaxis()->SetTitleSize(0.06);
    hist->GetYaxis()->SetTitleSize(0.06);
    hist->GetXaxis()->SetTitleFont(42);
    hist->GetYaxis()->SetTitleFont(42);
    hist->GetXaxis()->SetNdivisions(Ndivision);
    hist->GetYaxis()->SetTitleOffset(1.4);
//    hist->GetXaxis()->CenterTitle();
//    hist->GetYaxis()->CenterTitle();
    hist->GetXaxis()->SetLabelFont(42);
    hist->GetYaxis()->SetLabelFont(42);
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);
    hist->SetMarkerStyle(marker_style);
    hist->SetMarkerSize(marker_size);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
}

TLegend* plotLegend(TString pos="right_top",TString Title="No Title",Float_t scaleX=0.9,Float_t scaleY=0.9,Float_t offsetX=0.0,Float_t offsetY=0.0,TString Comment="",Int_t commentcolor=1)
{
    // pos places the Legend and can be
    // right_top, mid_top, left_top, right_bottom,mid_bottom,left_bottom
    // Title gives you the Title of the legend
    // scaleX,scaleY defines the size (=1 means 1/3 pad with,1/3 pad height (according to pos))
    // offsetX , offsetY (NDC coordinates of pad -> 0.1) shift by 10% of pad with) shifts the legend by the value in x and y
    // comment is optional text below title before legend entris will start
    // commentcolor defines the text color of the comment

    Float_t left  =gPad->GetLeftMargin()*1.15;
    Float_t right =1-gPad->GetRightMargin();
    Float_t top   =1-gPad->GetTopMargin();
    Float_t bottom=gPad->GetBottomMargin()*1.15;
    Float_t mid   =gPad->GetLeftMargin() + (1-(gPad->GetRightMargin()+gPad->GetLeftMargin()))/2;
    Float_t width =(right-left)/2;
    Float_t heith =(top-bottom)/2;
    TLegend* legend;
    TLine* dummy=new TLine();
    dummy->SetLineColor(0);

    if(pos.CompareTo("left_top")==0)    legend=new TLegend(left+offsetX,top+offsetY-(scaleY*heith),left+offsetX+(scaleX*width),top+offsetY,Title.Data());                           // left top
    if(pos.CompareTo("mid_top")==0)     legend=new TLegend(mid+offsetX-((scaleX*width)/2),top+offsetY-(scaleY*heith),mid+offsetX+((scaleX*width)/2),top+offsetY,Title.Data());      // mid up
    if(pos.CompareTo("right_top")==0)   legend=new TLegend(right+offsetX-(scaleX*width),top+offsetY-(scaleY*heith),right+offsetX,top+offsetY,Title.Data());                         // right top
    if(pos.CompareTo("left_bottom")==0) legend=new TLegend(left+offsetX,bottom+offsetY,left+offsetX+(scaleX*width),bottom+offsetY+(scaleY*heith),Title.Data());                     // left bottom
    if(pos.CompareTo("mid_bottom")==0)  legend=new TLegend(mid+offsetX-((scaleX*width)/2),bottom+offsetY,mid+offsetX+((scaleX*width)/2),bottom+offsetY+(scaleY*heith),Title.Data());// mid down
    if(pos.CompareTo("right_bottom")==0)legend=new TLegend(right+offsetX-(scaleX*width),bottom+offsetY,right+offsetX,bottom+offsetY+(scaleY*heith),Title.Data());                   // right bottom

    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetMargin(0.15);



    if(Comment.CompareTo("")!=0)
    {   // comment below header
        TLegendEntry* entry=legend->AddEntry(dummy,Comment.Data());
        entry->SetTextColor(commentcolor);
        entry->SetTextFont(62);
        entry->SetOption("l");
    }
    return legend;
}
void setLegendEntry(TLegend* legend,TObject* object,char* label,Int_t col,TString opt,Float_t size=0.05)
{
    // add entry for object with label and color and option to legend
    // legend is thge pointer to the legend where you want to add the entry
    // object ist the histogram or graph which will be added
    // label is the text which will appear in the entry
    // col is th color of the text
    // opt is the option ("pL"(marker and line), "p"(marker),"l"(line))
    TLegendEntry* entry=legend->AddEntry((TObject*)object,label);
    entry->SetTextColor(col);
    entry->SetTextSize(size);
    entry->SetOption(opt.Data());
}
TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.05,Int_t color=1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

    if(x<0||y<0)
    {   // defaults
        x=gPad->GetLeftMargin()*1.15;
        y=(1-gPad->GetTopMargin())*1.02;
    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextSize(size);
    text->SetNDC();
    text->SetTextColor(color);
    text->Draw();
    return text;
}







Int_t getColor(Int_t i)
{
    // pre defined set of 10 colors
    // i is larger 9 it will start from the first color again
    Int_t colors[10]   ={2,4,6,8,38,46,14,1,30,43};
    return colors[i%10];
}
Int_t getMarker(Int_t i)
{
    // pre defined set of 10 marker styles
    // i is larger 9 it will start from the first marker style again
    Int_t style[10]   ={20,21,22,23,24,25,26,27,28,29};
    return style[i%10];
}
void setGraph(TGraph* graph,Int_t i,Int_t markerstyle=20,Int_t markercolor=1,Float_t markersize=1.,Int_t linecolor=1)
{
    // sets the graphics of graph (use full in a loop)
    // graph is the pointer of the graph
    //
    // i switches to modes:
    // (if i>0) i the number of the color of the marker and the markert style (getColor(),getMarker())
    // (if i<0) marker color and style will be  markercolor, markerstyle
    //
    // markersize, and linecolor defines the size of the marker and the line color
    // can be use in two ways:
    // values >=0 the style and color will be set via getColor(i),getMarker(i)
    // values < 0 the style and color will be set -markercolor -markerstyle

    if(i<0)
    {
        graph->SetMarkerColor(markercolor);
        graph->SetMarkerSize(markersize);
        graph->SetMarkerStyle(markerstyle);
        graph->SetLineColor(linecolor);
    }else{
        if     (markercolor>=0)graph->SetMarkerColor(getColor(i));
        else if(markercolor<0) graph->SetMarkerColor(getColor(-markercolor));
        graph->SetMarkerSize(markersize);
        if     (markerstyle>=0)graph->SetMarkerStyle(getMarker(i));
        else if(markerstyle<0) graph->SetMarkerStyle(getMarker(-markerstyle));
        graph->SetLineColor(linecolor);
    }
}
void setHist(TH1* hist,Int_t i,Int_t markerstyle=20,Int_t markercolor=1,Float_t markersize=1.,Int_t linecolor=1)
{
    // sets the graphics of histogram (use full in a loop)
    // hist is the pointer of the hiistogram
    //
    // i switches to modes:
    // (if i>0) i the number of the color of the marker and the markert style (getColor(),getMarker())
    // (if i<0) marker color and style will be  markercolor, markerstyle
    //
    // markersize, and linecolor defines the size of the marker and the line color
    // can be use in two ways:
    // values >=0 the style and color will be set via getColor(i),getMarker(i)
    // values < 0 the style and color will be set -markercolor -markerstyle
    if(i==-99)
    {
        hist->SetMarkerColor(markercolor);
        hist->SetMarkerSize(markersize);
        hist->SetMarkerStyle(markerstyle);
        hist->SetLineColor(linecolor);
    }else{
        if     (markercolor>=0)hist->SetMarkerColor(getColor(i));
        else if(markercolor<0) hist->SetMarkerColor(getMarker(-markercolor));
        hist->SetMarkerSize(markersize);
        if     (markerstyle>=0)hist->SetMarkerStyle(getMarker(i));
        else if(markerstyle<0) hist->SetMarkerStyle(getMarker(-markerstyle));
        hist->SetLineColor(linecolor);
      }
}

