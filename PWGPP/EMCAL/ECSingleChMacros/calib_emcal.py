## @package calib_emcal
#  calib_emcal requries a root file input with channel energy vs channel ID histogram. Extracts that histrogram, then creates single channel specturm, sums the good channels, and applies a power law fit to the sum to obtain the exponentatial. Then calls factorized_model to obtain single channel calibration coefficients. 

#!/usr/bin/env python

## Computes single channel calibration factor and provides before and after chi-square for each channel.
def factorized_model(hcount, hspectrum, bad_channel, gamma,hChi):
    ncell = hcount.GetNbinsX()
    content = [hcount.GetBinContent(i + 1) for i in range(ncell)]
    neta_max = 48
    nphi_max = 24
    nsm = 20
    mean_ieta = [[0 for i in range(neta_max)] for j in range(nsm)]
    count_ieta = [[0 for i in range(neta_max)] for j in range(nsm)]
    mean_iphi_2 = [[0 for i in range(2 * nphi_max)] for j in range(nsm)]
    count_iphi_2 = [[0 for i in range(2 * nphi_max)] for j in range(nsm)]
    mean_sm = [0 for j in range(nsm)]
    count_sm = [0 for j in range(nsm)]
    for i in range(ncell):
        if content[i] > 0 and i not in bad_channel:
            sm, ieta, iphi = paint_emcal.to_sm_ieta_iphi(i)
            eta, phi = paint_emcal.to_eta_phi(sm, ieta, iphi)
            mean_ieta[sm][ieta] += content[i]
            count_ieta[sm][ieta] += 1
            iphi_2 = iphi + nphi_max * (ieta % 2)
            mean_iphi_2[sm][iphi_2] += content[i]
            count_iphi_2[sm][iphi_2] += 1
            mean_sm[sm] += content[i]
            count_sm[sm] += 1
    #for loop over range(ncell) ends
    
    for sm in range(nsm):
        for ieta in range(neta_max):
            if count_ieta[sm][ieta] > 0:
                mean_ieta[sm][ieta] /= count_ieta[sm][ieta]
        for iphi_2 in range(2 * nphi_max):
            if count_iphi_2[sm][iphi_2] > 0:
                mean_iphi_2[sm][iphi_2] /= count_iphi_2[sm][iphi_2]
        if count_sm[sm] > 0:
            mean_sm[sm] /= count_sm[sm]
    #for loop over range(nsm) ends

    print_warm = True
    print_coeffs = False
    print_chis = False
    makePlots = False
    check_ietaDependence = True
    function_inner = []

    ieta_hspectrum=[]
    if check_ietaDependence:
        for i in range(ncell):
            root_histogram_spectrum.append(ROOT.TH1D(hsum.ProjectionY('e', i + 1, i + 1)))
            if not i in bad_all:
                if hsumpy == None:
                    hsumpy = ROOT.TH1D(root_histogram_spectrum[-1])
                else:
                    hsumpy.Add(root_histogram_spectrum[-1])
    
    if print_coeffs:
        file = open("lhc15o_list5_coeff.txt","w")
    if print_chis:
        file_chi = open("lhc15o_3_new_chiBeforeAfter.txt","w")
    for i in plot_cell_id + range(ncell):
        sm, ieta, iphi = paint_emcal.to_sm_ieta_iphi(i)
        eta, phi = paint_emcal.to_eta_phi(sm, ieta, iphi)
       # n(ieta, iphi) ~ m(ieta) * m(iphi) / m
        content_model = 0
        if mean_sm[sm] > 0:
            content_model = mean_ieta[sm][ieta] * mean_iphi_2[sm][iphi + nphi_max * (ieta % 2)] / mean_sm[sm]#The average energy you expect for the cell for the super module.
        # if i not in bad_channel and content[i] > 0 and content_model > 0:
        #     print (content_model / content[i])**(1.0 / 5.0), ','
        if content[i] > 0:# and content[i] < 5:
            a = 0
            if print_warm:
                if (not i in bad_channel_lhc15o_3.bad_all) and (content[i] - content_model) > 2000:
                    sys.stdout.write('%d, ' % i)
            elif content_model > 0:
                a = (content_model / content[i])**(-1.0 / gamma) #a is set as the correction
                if print_coeffs:
                    print i, a
                    file.write('%d\t%f\n' %(i, a))
                if a != 0:
                    # scale_15_17, scale_17_20, scale_20_50 = [((a * r[1])**(gamma + 1) / (gamma + 1) - (a * r[0])**(gamma + 1) / (gamma + 1)) / (r[1]**(gamma + 1) / (gamma + 1) - r[0]**(gamma + 1) / (gamma + 1)) for r in ((1.5, 17), (1.7, 2.0), (2.0, 5.0))]
                    chi_square = []
                    if len(plot_cell_id) >= 1 and i == plot_cell_id[0]:
                        histogram_model = [ROOT.TH1D(hspectrum[i]) for j in range(2)]
                        for h in histogram_model:
                            h.Reset()
                            
                    for ia1 in range(2):
                        if ia1 == 0:
                            a1 = a                            
                        else:
                            a1 = 1
                        s = 0
                        dof = 0
                        
                        
                        for j in range(hspectrum[i].GetXaxis().FindBin(1.5 / a1), hspectrum[i].GetXaxis().FindBin(5 / a1) + 1):
                            content_bin = hspectrum[i].GetBinContent(j)
                            # content_model_bin = content_model * (hspectrum[i].GetXaxis().GetBinLowEdge(j + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(j)**(gamma + 1) / (gamma + 1)) / (hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(5 / a1) + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(1.5 / a1))**(gamma + 1) / (gamma + 1))
                            # content_model_bin_orig = content_model * (hspectrum[i].GetXaxis().GetBinLowEdge(j + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(j)**(gamma + 1) / (gamma + 1)) / (hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(5) + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(1.5))**(gamma + 1) / (gamma + 1))
                            content_model_bin = content_model * (hspectrum[i].GetXaxis().GetBinLowEdge(j + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(j)**(gamma + 1) / (gamma + 1)) / (hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(5 / a1) + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(hspectrum[i].GetXaxis().FindBin(1.5 / a1))**(gamma + 1) / (gamma + 1))
                            #content_model_bin_orig = content_model * (hspectrum[i].GetXaxis().GetBinLowEdge(j + 1)**(gamma + 1) / (gamma + 1) - hspectrum[i].GetXaxis().GetBinLowEdge(j)**(gamma + 1) / (gamma + 1)) / ((5 / a1)**(gamma + 1) / (gamma + 1) - (1.5 / a1)**(gamma + 1) / (gamma + 1))

                            #if i == 1056:
                               #print i, ia1, a1, j, hspectrum[i].GetBinCenter(j), content_bin, content_model_bin
                            if len(plot_cell_id) >= 1 and i == plot_cell_id[0]:
                                if ia1 == 0:
                                    histogram_model[0].SetBinContent(j, content_model_bin)
                                else:
                                    histogram_model[1].SetBinContent(j, content_model_bin)
                            if content_bin > 0:
                               # if len(plot_cell_id) >= 1:
                                    #print i, ia1, a1, j, hspectrum[i].GetBinCenter(j), content_bin, content_model_bin#, content_model_bin_orig
                                s += (content_bin - content_model_bin)**2 /  content_bin
                                dof += 1
                        if dof > 0:
                            s /= dof
                        chi_square.append(s)

                    hChi[0].Fill(chi_square[0])
                    hChi[1].Fill(chi_square[1])
                       
                    #print i, i in bad_channel_lhc15o_3.bad_all and 1 or 0, i in bad_channel_lhc15o_mine3.hot and 1 or 0, a, chi_square[0], chi_square[1]

                    if print_chis:
                        file_chi.write("%d %d %d %f %f %f\n" %(i, i in bad_channel_lhc15o_3.bad_all and 1 or 0, i in bad_channel_lhc15o_mine3.hot and 1 or 0, a, chi_square[0], chi_square[1])) 

                    
                    if chi_square[0] / chi_square[1] < 1 and i in bad_channel_lhc15o_3.bad_all:
                        #a = chi_square[0] / chi_square[1]    
                        #if a < 1.5 and chi_square[0] < 20:
                        if chi_square[0] < 10:
                            #a = 1
                            #print i
                            if i in bad_channel_lhc15o_all_breakdown.dead_all:
                                a = 6
                            if i in bad_channel_lhc15o_all_breakdown.badC_all:
                                a = 5
                            if i in bad_channel_lhc15o_all_breakdown.warm_all:
                                a = 10
                            #if (i == 6083):
                            #print "hello", i, chi_square   
                                
                        else:
                            if i in bad_channel_lhc15o_3.bad_all:
                                a = 1
                    else:
                        a = 0
                        None
                    if len(plot_cell_id) >= 1 and i == plot_cell_id[0]:
                        canvas.cd()
                        canvas.SetLogx()
                        canvas.SetLogy()
                        hspectrum[i].SetMarkerStyle(20)
                        histogram_model[0].SetMarkerStyle(24)
                        histogram_model[1].SetMarkerStyle(25)
                        hspectrum[i].GetXaxis().SetTitle('E_T\\:(\\mathrm{GeV})')
                        hspectrum[i].GetYaxis().SetTitle('Counts')
                        hspectrum[i].GetYaxis().SetTitleFont(42)
                        hspectrum[i].GetYaxis().SetLabelFont(42)
                        hspectrum[i].SetMinimum(1)
                        hspectrum[i].Draw('e1x0')
                        histogram_model[0].Draw('e1x0same')
                        histogram_model[1].Draw('e1x0same')
                        canvas.Update()
                        ROOT.gApplication.Run()
        #print a
        #content[i] = a
        content[i] -= content_model
        #content[i] = content[i]
        content[i] = max(-2000, min(2000, content[i]))
    #for loop over  plot_cell_id + range(ncell) ends
    if print_chis:
        file_chi.close()
    if print_warm:
        sys.stdout.write('\n')
    if print_coeffs:
        file.close()
    for i in range(ncell):
        hcount.SetBinContent(i + 1, content[i])
    return hcount

## Finds gamma for a summed spectrum of all the channels in an eta for each super module. Unfinished at the moment. Do not call.
def getEtaDependence(hcount, hspectrum, bad_channel, gamma,hChi):
    ncell = hcount.GetNbinsX()
    content = [hcount.GetBinContent(i + 1) for i in range(ncell)]
    neta_max = 48
    nphi_max = 24
    nsm = 20
    ietaDict = {}
    for i in range(ncell):
        if content[i] > 0 and i not in bad_channel:
            sm, ieta, iphi = paint_emcal.to_sm_ieta_iphi(i)
            eta, phi = paint_emcal.to_eta_phi(sm, ieta, iphi)
            tmpList = [sm, ieta, iphi]
            ietaDict[i] = tmpList
            print i, ietaDict[i], '\n'
            
    #for loop over range(ncell) ends

    #sumiEta = []
    #for i in range(ncell):
    #root_histogram_spectrum.append(ROOT.TH1D(hsum.ProjectionY('e', i + 1, i + 1)))
    #    if not i in bad_all:
    #        if hsumpy == None:
    #            hsumpy = ROOT.TH1D(root_histogram_spectrum[-1])
    #        else:
    #            hsumpy.Add(root_histogram_spectrum[-1])
                                
    
#Main block
if __name__ == '__main__':
    import os, sys
    sys.path.append(os.path.join(os.environ['ROOTSYS'], 'lib'))
    import math, array, re, ROOT, paint_emcal
    from ROOT import TCanvas, TFile, TH1D, TH2D, TLine, TH1F, TColor, TStyle
    application_name = ' '.join(sys.argv)
    paint_emcal.set_root_style()
    log_z = False
    filename_list = []#Takes in all the input files
    plot_cell_id = []
    #Takes in the system arguments
    for f in sys.argv[1:]:
        if f in ('-l', '--log-z'):#Sets z axis to log
            log_z = True
        elif re.compile('^[0-9]+$').match(f):#checking for a histogram input
            plot_cell_id.append(int(f))
        else:
            filename_list.append(f)#adding in file inputs
    if len(plot_cell_id) > 0:
        canvas = canvas = ROOT.TCanvas('canvas%d' % 0, application_name,800 + 4, 600 + 28)
        pad = []
    else:
        canvas, pad = paint_emcal.alice_emcal_canvas_pad(application_name)
    ncell = 17664

    # for i in range(ncell):
    #     sm, ieta, iphi = paint_emcal.to_sm_ieta_iphi(i)
    #     if sm == 9 and ieta == 33 and iphi >= 3 and iphi <= 7:
    #         print i
    # sys.exit(0)

    #Create a 1D histogram
    root_histogram = []
    for i in range(1):
        root_histogram.append(ROOT.TH1D('root_histogram%d' % len(root_histogram), '',ncell, -0.5, ncell + 0.5))

        
    #Obtain the 2D histogram which contents cell_id on x and energy on y from the AnalysisResults.root
    hsum = None
    for filename in filename_list:
        f = ROOT.TFile.Open(filename)
        h = f.Get('AliAnalysisTaskCalibEmcal').Get('histogram').FindObject('_histogram_cell_id_amplitude')
        if hsum == None:
            hsum = ROOT.TH2D(h)
        else:
            hsum.Add(h)
        hpx = ROOT.TH1D(h.ProjectionX('hpx',h.GetYaxis().FindBin(1.5),h.GetYaxis().FindBin(5) - 1))
        print h.GetYaxis().FindBin(1.7), h.GetYaxis().FindBin(2)
        hpx_15_17 = ROOT.TH1D(h.ProjectionX('hpx',h.GetYaxis().FindBin(1.5),h.GetYaxis().FindBin(1.7) - 1))
        hpx_17_20 = ROOT.TH1D(h.ProjectionX('hpx',h.GetYaxis().FindBin(1.7),h.GetYaxis().FindBin(2) - 1))
        hpx_20_50 = ROOT.TH1D(h.ProjectionX('hpx',h.GetYaxis().FindBin(2)  ,h.GetYaxis().FindBin(5) - 1))
        #hpx.Divide(hpx_17_20, hpx_20_50)
        root_histogram[-1] = ROOT.TH1D(hpx)


    #importing the bad channel lists
    import bad_channel_lhc15o_3
    import bad_channel_lhc15o_mine3
    import bad_channel_lhc15o_all_breakdown

    #bad_channel_lhc15o_mine3.hot = []
    #print bad_channel_lhc15o_mine3.hot
    #import lhc15o_list2_coeffs

    #hot=[1160, 1275, 1682, 1683, 1685, 1686, 1961, 3544, 3849, 3858, 3878, 3884, 4245, 4566, 5836, 6111, 6800, 6801, 6802, 6804, 6805, 6806, 6807, 6809, 7104, 7874, 8811, 8812, 11088, 11101, 11142, 11144, 11146, 11871, 11950, 12063, 13253, 14396, 15316, 16340, 16346, 17070]
    hot=[1682, 3544, 7874, 8811, 8812, 11088, 11144, 11146]
    badChn = sorted(bad_channel_lhc15o_3.bad_all + hot)
    

    bad_all = sorted(bad_channel_lhc15o_3.bad_all + bad_channel_lhc15o_mine3.hot)#list of all bad and hot channels for the run set

    
    function = []


    # hpx_e = hpx.Clone()
    # cspece = TCanvas("cspece", "cspece",600, 400)
    # for i in range(1,hpx.GetXaxis().GetNbins()):
    #     if (i - 1) in bad_channel_lhc15o_3.bad_all:
    #         hpx_e.SetBinContent(i,0)
    # hpx_e.Draw()

    # hpx_d = hpx.Clone()
    # cspecd = TCanvas("cspecd", "cspecd",600, 400)
    # for i in range(1,hpx.GetXaxis().GetNbins()):
    #     if (i - 1) in  hot:
    #         hpx_d.SetBinContent(i,0)
    # hpx_d.Draw()
    
    # hpx_all = hpx.Clone()
    # cspecall = TCanvas("cspecall", "cspecall",600, 400)
    # for i in range(1,hpx.GetXaxis().GetNbins()):
    #     if (i - 1) in  badChn:
    #         hpx_all.SetBinContent(i,0)
    # hpx_all.Draw()



    #histograms to obtaina  chisqr distribution
    root_histogram_chiHists = []
    #hist at 0 is calibrated chi_sqr
    #hist at 1 is uncalibrated chi_sqr
    for i in range(2):
        root_histogram_chiHists.append(ROOT.TH1D('root_histogram_chiHists%d' % len(root_histogram), '', 10000, -0.5, 9999.5))


    #Getting the energies from the cell_amplitude histogram
    hsumpy = None
    root_histogram_spectrum = []
    for i in range(ncell):
        root_histogram_spectrum.append(ROOT.TH1D(hsum.ProjectionY('e', i + 1, i + 1)))
        if not i in bad_all:
            if hsumpy == None:
                hsumpy = ROOT.TH1D(root_histogram_spectrum[-1])
            else:
                hsumpy.Add(root_histogram_spectrum[-1])

    #The power law fitting function
    function.append(ROOT.TF1('function%d' % len(function), 'exp([0])*x^[1]', 1.5, 5))
    function[-1].SetParameter(0, 1)
    function[-1].SetParameter(1, -4)

    #Fitting the power law to the sum of spectrum in all the good channels.
    for i in range(1, hsumpy.GetNbinsX()):
        hsumpy.SetBinContent(i, hsumpy.GetBinContent(i) / hsumpy.GetXaxis().GetBinWidth(i))
        hsumpy.SetBinError(i, hsumpy.GetBinError(i) / hsumpy.GetXaxis().GetBinWidth(i))
    hsumpy.Fit(function[-1], 'r0')

    #csumpy = TCanvas("csumpy", "csumpy",600, 400)
    #hsumpy.Draw()

    #cspec = TCanvas("cspec", "cspec",600, 400)
    #root_histogram[13287].Draw()
    
    # for i in range(hsumpy.GetXaxis().FindBin(1.5),
    #                hsumpy.GetXaxis().FindBin(5)):
    #     print i, hsumpy.GetXaxis().GetBinCenter(i), hsumpy.GetBinContent(i)

    #bad_all = bad_channel_lhc15o_1.bad_all
    # if pad[0] == None: sys.exit(1)

    #cspeca = TCanvas("cspeca", "cspeca",600, 400)
    #root_histogram_spectrum[1056].Draw()
    
    
    root_histogram[-1] = factorized_model(root_histogram[-1], root_histogram_spectrum, bad_all, function[-1].GetParameter(1),root_histogram_chiHists)

    #getEtaDependence(root_histogram[-1], root_histogram_spectrum, bad_all, function[-1].GetParameter(1),root_histogram_chiHists)
                    
    
    canvas.cd()
    for i in range(len(pad)):
        pad[i].Draw()
        if log_z:
            pad[i].SetLogz()
    paint_emcal.update(canvas, pad, root_histogram)
    ROOT.gApplication.Run()

    
