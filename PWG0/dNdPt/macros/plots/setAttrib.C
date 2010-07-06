//------------------------------------------------------------------------------
// setAttrib.C
//
// helper class to set various attributes for Tgraphs/TPads
//------------------------------------------------------------------------------



void setAttrib(TGraph* g) 
{
    if (!g) return;
    if (g->GetXaxis()) {
        g->GetXaxis()->SetLabelSize(labelSize);
        g->GetXaxis()->SetTitleSize(titleSize);
        g->GetXaxis()->SetLabelFont(font);
        g->GetXaxis()->SetTitleFont(font);
        g->GetXaxis()->SetTitleOffset(3.8);
//         g->GetXaxis()->SetLabelOffset(0.02); 
    }

    if (g->GetYaxis()) {
        g->GetYaxis()->SetLabelSize(labelSize);
        g->GetYaxis()->SetTitleSize(titleSize);
        g->GetYaxis()->SetLabelFont(font);
        g->GetYaxis()->SetTitleFont(font);
        g->GetYaxis()->SetTitleOffset(2.2);       
    }
}

void setAttrib(TMultiGraph* g) 
{
    if (!g) return;
    if (g->GetXaxis()) {
        g->GetXaxis()->SetLabelSize(labelSize);
        g->GetXaxis()->SetTitleSize(titleSize);
        g->GetXaxis()->SetLabelFont(font);
        g->GetXaxis()->SetTitleFont(font);
        g->GetXaxis()->SetTitleOffset(3.8); // was 3.8
// //         g->GetXaxis()->SetLabelOffset(0.02); 
    }

    if (g->GetYaxis()) {
        g->GetYaxis()->SetLabelSize(labelSize);
        g->GetYaxis()->SetTitleSize(titleSize);
        g->GetXaxis()->SetLabelFont(font);
        g->GetYaxis()->SetTitleFont(font);
        g->GetYaxis()->SetTitleOffset(2.2); // was 2.4
//         g->GetYaxis()->SetLabelOffset(0.02); 
    }
}

void setAttrib(TPad* p) 
{
    if (!p) return;
    if (p->GetYlowNDC() > 0.01 ) {
        p->SetTopMargin(0.03);
        p->SetBottomMargin(0.0);
        p->SetLeftMargin(0.16);
        p->SetRightMargin(0.03);
        p->SetTicks(1,1);
    }
    if (p->GetYlowNDC() < 0.01 ) {
        p->SetTopMargin(0.0);
        p->SetLeftMargin(0.16);
        p->SetRightMargin(0.03);
        p->SetBottomMargin(0.28);
        p->SetTicks(1,1);
    }
}