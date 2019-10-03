## @package paint_emcal
#  paint_emcal creates the canvas and pad where the histogram output of calib_emcal can be displayed in as individual channels of the emcal and dcal in their correct real eta and phi position.

#!/usr/bin/env python

if __name__ != '__main__':
    import math, array, ROOT

## Sets up ROOT
def set_root_style():
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetNumberContours(253)
    ROOT.gStyle.SetPalette(55, ROOT.nullptr)

## Creates a canvas in order plot the emcal and dcal channels
def multipanel_pad(canvas, pad_layout_unscaled, margin = None):
    pad = []
    pad_layout = pad_layout_unscaled
    for i in range(2):
        s = sum(pad_layout_unscaled[i])
        for j in range(len(pad_layout_unscaled[i])):
            pad_layout[i][j] /= s
    pad_layout_1_reversed = list(reversed(pad_layout[1]))
    if margin == None:
        left_margin = canvas.GetLeftMargin()
        right_margin = canvas.GetRightMargin()
        top_margin = canvas.GetTopMargin()
        bottom_margin = canvas.GetBottomMargin()
    else:
        left_margin, right_margin, top_margin, bottom_margin = magin
    lr = left_margin + right_margin
    tb = top_margin + bottom_margin
    for i in range(len(pad_layout[0])):
        x0 = sum(pad_layout[0][:i]) * (1 - lr)
        x1 = lr + sum(pad_layout[0][:i+1]) * (1 - lr)
        for j in reversed(range(len(pad_layout[1]))):
            y0 = sum(pad_layout_1_reversed[:j]) * (1 - tb)
            y1 = tb + sum(pad_layout_1_reversed[:j+1]) * (1 - tb)
            name = 'pad%d' % len(pad)
            canvas.cd()
            pad.append(ROOT.TPad(name, name, x0, y0, x1, y1))
            pad[-1].SetLeftMargin(
                left_margin / (lr + pad_layout[0][i] * (1 - lr)))
            pad[-1].SetRightMargin(
                right_margin / (lr + pad_layout[0][i] * (1 - lr)))
            pad[-1].SetTopMargin(
                top_margin / (tb + pad_layout_1_reversed[j] * (1 - tb)))
            pad[-1].SetBottomMargin(
                bottom_margin / (tb + pad_layout_1_reversed[j] * (1 - tb)))
            pad[-1].SetFillStyle(0)
            pad[-1].SetFillColor(0)
            pad[-1].Draw()
            pad[-1].cd()
    return pad

## Converts from channel ID to super module number, and phi and eta coordinates in the supermodule
def to_sm_ieta_iphi(n):
    if n < 11520:
        sm = n / 1152
        n0 = sm * 1152
        nphi = 24
    elif n < 12288:
        sm = 10 + (n - 11520) / 384
        n0 = 11520 + (sm - 10) * 384
        nphi = 8
    elif n < 16896:
        sm = 12 + (n - 12288) / 768
        n0 = 12288 + (sm - 12) * 768
        nphi = 24
    else:
        sm = 18 + (n - 16896) / 384
        n0 = 16896 + (sm - 18) * 384
        nphi = 8
    n1 = n - n0
    ieta = 2 * (n1 / (2 * nphi)) + 1 - (n1 % 2);
    iphi = (n1 / 2) % nphi;
    return sm, ieta, iphi

## Converts local eta and phi, to real world eta and phi
def to_eta_phi(sm, ieta, iphi, ieta_int = None):
    if ieta_int == None:
        ieta_int = round(ieta)
    coeff_eta = [
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538761263736301, 0.0012502455282809458,
          -0.00005435850122961024, -0.013757648176291852,
          -0.000025569020644999157, 1.1116965497826313e-6 ],
        [ -0.007266662087911644, -0.000048501557965879184,
          2.108763389820745e-6, -0.013757648176291801,
          -0.000025569020645007092, 1.1116965497829259e-6 ],
        [ 0.6538955994897965, 0.0012327876984126996,
          -0.00005219165046974383, -0.013757946609494863,
          -0.000025332446808511946, 1.086380600873213e-6 ],
        [ -0.0072721088435373045, -0.00004216269841272758,
          1.131762228686001e-6, -0.01375794660949486,
          -0.00002533244680851135, 1.0863806008729037e-6 ],
        [ 0.6537747727272758, 0.0012291942148758277,
          -0.00005344322673373451, -0.01375054364989851,
          -0.00002391323211245261, 1.0397057440196336e-6 ],
        [ -0.22750791958042052, -0.0004878840193896497,
          0.000021212348669116138, -0.01375054364989848,
          -0.00002391323211246439, 1.0397057440201827e-6 ],
        [ 0.6537747727272758, 0.0012291942148758277,
          -0.00005344322673373451, -0.01375054364989851,
          -0.00002391323211245261, 1.0397057440196336e-6 ],
        [ -0.22750791958042052, -0.0004878840193896497,
          0.000021212348669116138, -0.01375054364989848,
          -0.00002391323211246439, 1.0397057440201827e-6 ],
        [ 0.6537747727272758, 0.0012291942148758277,
          -0.00005344322673373451, -0.01375054364989851,
          -0.00002391323211245261, 1.0397057440196336e-6 ],
        [ -0.22750791958042052, -0.0004878840193896497,
          0.000021212348669116138, -0.01375054364989848,
          -0.00002391323211246439, 1.0397057440201827e-6 ],
        [ 0.6538955994897965, 0.0012327876984126996,
          -0.00005219165046974383, -0.013757946609494863,
          -0.000025332446808511946, 1.086380600873213e-6 ],
        [ -0.0072721088435373045, -0.00004216269841272758,
          1.131762228686001e-6, -0.01375794660949486,
          -0.00002533244680851135, 1.0863806008729037e-6 ] ]

    coeff_phi = [
        [ 1.4153378750000025, 0.0006568888888901503,
          0.013518083333333313, -0.000057045289855166094 ],
        [ 1.4159947638888941, -0.0006568888888867641,
          0.013461038043478247, 0.00005704528985495655 ],
        [ 1.7644118472222259, 0.000646180555557041,
          0.01351747826086953, -0.00005638586956534543 ],
        [ 1.7650580277777765, -0.0006461805555546355,
          0.013461092391304308, 0.00005638586956517725 ],
        [ 2.1134724027777834, 0.0006524861111138665,
          0.013517768115941999, -0.00005654166666683834 ],
        [ 2.1141248888888984, -0.0006524861111081304,
          0.013461226449275318, 0.000056541666666427527 ],
        [ 2.462534333333343, 0.0006581805555594584,
          0.013518074275362234, -0.00005697644927563384 ],
        [ 2.463192513888897, -0.0006581805555524085,
          0.013461097826086869, 0.000056976449275180296 ],
        [ 2.8116099861111175, 0.0006471666666698442,
          0.013517443840579645, -0.00005650181159440986 ],
        [ 2.8122571527777906, -0.0006471666666621282,
          0.013460942028985465, 0.00005650181159390807 ],
        [ -3.12222569444445, 0.0006597222222129356,
          0.013424007936509517, -0.00005902777777630214 ],
        [ -3.121565972222227, -0.0006597222222309393,
          0.013364980158731608, 0.000059027777779336 ],
        [ -1.7262874166666704, 0.0008678541666655771,
          0.013520961956521787, -0.00007551086956511075 ],
        [ -1.7254195625000033, -0.000867854166667594,
          0.013445451086956564, 0.00007551086956536232 ],
        [ -1.3772168333333357, 0.000859541666666056,
          0.013520576086956563, -0.00007492391304337749 ],
        [ -1.3763572916666709, -0.000859541666667937,
          0.013445652173913077, 0.00007492391304362568 ],
        [ -1.0281504166666693, 0.0008688541666657723,
          0.013520394021739162, -0.00007519021739121045 ],
        [ -1.0272815625000018, -0.000868854166667404,
          0.013445203804347865, 0.00007519021739139448 ],
        [ -0.6787590277777789, 0.0006451388888869277,
          0.013420585317460642, -0.00005367063492026231 ],
        [ -0.67811388888889, -0.0006451388888908697,
          0.01336691468254001, 0.000053670634921005446 ] ]

    if abs(ieta - round(ieta)) >= 0.25:
        ieta_mod_2 = ieta_int % 2
    else:
        ieta_mod_2 = ieta % 2
    eta = coeff_eta[sm][0] + coeff_eta[sm][1] * iphi + \
          coeff_eta[sm][2] * iphi * iphi + \
          coeff_eta[sm][3] * ieta + \
          coeff_eta[sm][4] * ieta * iphi + \
          coeff_eta[sm][5] * ieta * iphi * iphi
    phi = coeff_phi[sm][0] + coeff_phi[sm][1] * ieta_mod_2 + \
          coeff_phi[sm][2] * iphi + \
          coeff_phi[sm][3] * iphi * ieta_mod_2
    if phi < -2:
        phi += 2 * math.pi
    return eta, phi

## Creates the emcal and dcal on the pad with real eta and phi range
def alice_emcal_canvas_pad(application_name, canvas_width = 1680):
    left_margin = 0.03125 * 1.5
    right_margin = 0.03125 * 2.1875
    top_margin = 0.03125 * 1.5
    bottom_margin = 0.03125 * 2.75
    azimuth_range = [(-1.85, -0.5), (1.3, 3.35)]
    pseudorapidity_range = (-0.8, 0.8)
    pad_layout = (map(lambda r: abs(r[1] - r[0]), azimuth_range),
                  [1.0])
    canvas_height = canvas_width * \
                    (1 - left_margin - right_margin) / \
                    sum(pad_layout[0]) * \
                    abs(pseudorapidity_range[1] -
                        pseudorapidity_range[0]) / \
                    (1 - top_margin - bottom_margin)
    canvas = ROOT.TCanvas(
        'canvas%d' % 0, application_name,
        canvas_width + 4, int(round(canvas_height)) + 28)
    canvas.SetLeftMargin(left_margin)
    canvas.SetRightMargin(right_margin)
    canvas.SetTopMargin(top_margin)
    canvas.SetBottomMargin(bottom_margin)
    return canvas, multipanel_pad(canvas, pad_layout)

## Fills the emcal and dcal canvas with either channel residue, energy/hits, chi-square, and more based on input from calib_emcal
def update(canvas, pad, root_histogram, root_histogram_index = 0,
           outline = False, text_size_sm = 0.03125 * 0.875,
           offset_sm = 3.5, text_size_ieta_iphi = 0.03125 * 0.625,
           offset_ieta_iphi = (0.5, 1),
           azimuth_range = ((-1.85, -0.5), (1.3, 3.35)),
           pseudorapidity_range = (-0.8, 0.8)):
    list_index = None
    try:
        index_list = (i for i, v in enumerate(root_histogram)
                      if isinstance(v, list)).next()
    except StopIteration:
        None
    if list_index != None:
        del root_histogram[list_index:]
    root_histogram.append([])
    scale_y = []
    for i in range(2):
        root_histogram[-1].append(ROOT.TH2D(
            'root_histogram%d_%d' %
            (len(root_histogram), len(root_histogram[-1])), '',
            1, azimuth_range[i][0], azimuth_range[i][1],
            1, pseudorapidity_range[0], pseudorapidity_range[1]))
        scale_y.append(
            (1 - pad[i].GetTopMargin() - pad[i].GetBottomMargin()) *
             pad[i].GetHNDC() * canvas.GetWh() /
             ((1 - pad[i].GetLeftMargin() - pad[i].GetRightMargin()) *
              pad[i].GetWNDC() * canvas.GetWw()))
        if i == 0:
            root_histogram[-1][-1].GetXaxis().SetLabelFont(42)
            root_histogram[-1][-1].GetYaxis().SetLabelFont(42)
            root_histogram[-1][-1].GetYaxis().SetTitle('\\eta')
            root_histogram[-1][-1].GetYaxis().SetTitleOffset(
                1.125 * scale_y[i])
        elif i == 1:
            root_histogram[-1][-1].Fill(0, 0, -1e-12)
            root_histogram[-1][-1].GetXaxis().SetLabelFont(42)
            root_histogram[-1][-1].GetXaxis().SetTitle('\\varphi')
            root_histogram[-1][-1].GetYaxis().SetLabelSize(0)
            root_histogram[-1][-1].GetZaxis().SetLabelFont(42)
        root_histogram[-1][-1].GetYaxis().SetTickLength(
            root_histogram[-1][-1].GetYaxis().GetTickLength() *
            scale_y[-1])
        root_histogram[-1][-1].GetZaxis().SetTickLength(
            root_histogram[-1][-1].GetZaxis().GetTickLength() *
            scale_y[-1])
        u0 = root_histogram[root_histogram_index].GetMinimum()
        u1 = root_histogram[root_histogram_index].GetMaximum()
        u1 = u0 + (1 + 1e-3) * (u1 - u0)
        #u0, u1 = 200, 1000
        root_histogram[-1][-1].SetMinimum(u0)
        root_histogram[-1][-1].SetMaximum(u1)
        root_histogram[-1][-1].GetZaxis().SetRangeUser(u0, u1)
    root_histogram.append([])
    ncell = 17664
    ncell_emcal = 12288
    pad[1].cd()
    root_histogram[-2][1].Draw('colz')
    pad[1].Update()
    palette = root_histogram[-2][1].GetListOfFunctions().FindObject(
        'palette')
    if palette != None:
        palette.SetX2NDC(
            palette.GetX1NDC() + scale_y[1] *
            (palette.GetX2NDC() - palette.GetX1NDC()))
    pad_index = 1
    for i in range(ncell):
        if i == ncell_emcal:
            root_histogram[-2][1].Draw('zsame')
            pad[0].cd()
            root_histogram[-2][0].Draw('colz')
            pad[0].Update()
            pad_index = 0
        content = root_histogram[root_histogram_index].\
                  GetBinContent(i + 1)
        sm, ieta, iphi = to_sm_ieta_iphi(i)
        p = [
            to_eta_phi(sm, ieta - 0.5, iphi - 0.5, ieta),
            to_eta_phi(sm, ieta - 0.5, iphi + 0.5, ieta),
            to_eta_phi(sm, ieta + 0.5, iphi + 0.5, ieta),
            to_eta_phi(sm, ieta + 0.5, iphi - 0.5, ieta),
            to_eta_phi(sm, ieta - 0.5, iphi - 0.5, ieta)
        ]
        x = array.array('d', [p[k][1] for k in range(len(p))])
        y = array.array('d', [p[k][0] for k in range(len(p))])
        root_histogram[-1].append(ROOT.TPolyLine(len(p), x, y))
        if content != 0:
            u = content
            if pad[pad_index].GetLogz() != 0 and u > 0:
                u = math.log10(u)
            u = max(root_histogram[0].GetMinimum(), min(
                root_histogram[0].GetMaximum(), u))
            color = palette.GetValueColor(u)
            root_histogram[-1][-1].SetFillColor(color)
        root_histogram[-1][-1].SetLineColor(ROOT.kGray)
        root_histogram[-1][-1].SetLineWidth(1)
        if content != 0:
            root_histogram[-1][-1].Draw('f')
        if outline:
            root_histogram[-1][-1].Draw()
        neta = sm >= 12 and sm < 18 and 32 or 48
        nphi = sm in (10, 11, 18, 19) and 8 or 24
        if ieta == 0 and iphi == 0:
            # -0.5 .. neta - 0.5
            y, x = to_eta_phi(
                sm, -(offset_sm + 0.5) + (neta + 2 * offset_sm) *
                (sm % 2),
                0.5 * (nphi - 1))
            root_histogram[-1].append(ROOT.TText(x, y, 'SM%d' % sm))
            root_histogram[-1][-1].SetTextAlign(
                sm % 2 == 0 and 21 or 23)
            root_histogram[-1][-1].SetTextFont(42)
            root_histogram[-1][-1].SetTextSize(text_size_sm)
            root_histogram[-1][-1].Draw()
            if sm % 2 == 1:
                if nphi == 8:
                    phi_list = 0, nphi - 1
                else:
                    phi_list = 0, 10, 20, nphi - 1
                for phi in phi_list:
                    y, x = to_eta_phi(
                        sm, neta + offset_ieta_iphi[0], phi)
                    root_histogram[-1].append(
                        ROOT.TText(x, y, '%d' % phi))
                    root_histogram[-1][-1].SetTextAlign(23)
                    root_histogram[-1][-1].SetTextFont(42)
                    root_histogram[-1][-1].SetTextSize(
                        text_size_ieta_iphi)
                    root_histogram[-1][-1].Draw()
            if sm in (0, 1, 12, 13):
                if neta == 32:
                    eta_list = 0, 10, 20, neta - 1
                else:
                    eta_list = 0, 10, 20, 30, 40, neta - 1
                for eta in eta_list:
                    eta_shifted = eta
                    if sm == 1 and eta == 0:
                        eta_shifted += 0.5
                    elif sm == 0 and eta == neta - 1:
                        eta_shifted -= 0.5
                    y, x = to_eta_phi(
                        sm, eta_shifted, -offset_ieta_iphi[1])
                    root_histogram[-1].append(
                        ROOT.TText(x, y, '%d' % eta))
                    root_histogram[-1][-1].SetTextAlign(32)
                    root_histogram[-1][-1].SetTextFont(42)
                    root_histogram[-1][-1].SetTextSize(
                        text_size_ieta_iphi)
                    root_histogram[-1][-1].Draw()
    canvas.Update()

#Main block
if __name__ == '__main__':
    import os, sys
    sys.path.append(os.path.join(os.environ['ROOTSYS'], 'lib'))
    import math, array, re, ROOT
    application_name = ' '.join(sys.argv)
    set_root_style()
    canvas, pad = alice_emcal_canvas_pad(application_name)
    ncell = 17664
    root_histogram = []
    log_z = False
    filename_list = []
    plot_range = None
    for f in sys.argv[1:]:
        if f in ('-l', '--log-z'):
            log_z = True
        elif f.find(':') != -1:
            plot_range = map(float, f.split(':'))
        else:
            filename_list.append(f)
    for filename in filename_list:
        root_histogram.append(ROOT.TH1D(
            'root_histogram%d' % len(root_histogram), '',
            ncell, -0.5, ncell + 0.5))
        f = open(filename, 'r')
        line = f.readline()
        while line != "":
            line_split = re.sub("+", " ", re.sub(" *$", "", re.sub("^ *", "", line[:-1]))).split(" ")
            cell_id = int(line_split[0]) - 1
            count = float(line_split[1])
            root_histogram[-1].Fill(cell_id, count)
            line = f.readline()
        if filename != filename_list[0]:
            root_histogram[0].Divide(root_histogram[0], root_histogram[-1])
        if plot_range != None:
            root_histogram[0].SetMinimum(plot_range[0])
            root_histogram[0].SetMaximum(plot_range[1])
    canvas.cd()
    for i in range(len(pad)):
        pad[i].Draw()
        if log_z:
            pad[i].SetLogz()
    update(canvas, pad, root_histogram)
    ROOT.gApplication.Run()
