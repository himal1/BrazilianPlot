import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend
import CMS_lumi, tdrStyle #himal block
import subprocess # to execute shell command
ROOT.gROOT.SetBatch(ROOT.kTRUE)
tdrStyle.setTDRStyle() #Himal block

# CMS style
CMS_lumi.cmsText = "CMS"
CMS_lumi.extraText = "Preliminary"
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = False
CMS_lumi.writeExtraText = 1
#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_13TeV = "132.9 fb^{-1}"
CMS_lumi.extraText = "Preliminary"
#CMS_lumi.lumiText = "18.3 fb^{-1}"
#








# CREATE datacards
def createMakeFileThetaB(labels,values):
    make_lines = ["void ",
                  "(int seed =5) {",
                  "using namespace RooFit;",
                  "RooRandom::randomGenerator()->SetSeed(seed);",
                  "TCanvas *c1 = new TCanvas(\"c1\",\"c1\");",
                  "RooWorkspace *w = new RooWorkspace();",
                  "w->factory(\"x[12,40]\");",
                  "reso = reso/1000.;",
                  "i=i/10.;",
                  "w->factory(\"CMS_H_width[0.12]\");",
                  "w->var(\"CMS_H_width\")->setVal(reso);",
                  "w->factory(\"CMS_H_mean[10.]\");",
                  "w->var(\"CMS_H_mean\")->setVal(i);",
                  "w->var(\"x\")->setBins(112);",
                  "w->factory(\"Exponential::background1(x,CMS_hjj_alpha[-0.0874,-1.,1.])\");//only expenential",
                  "//w->factory(\"Exponential::background2(x,CMS_hjj_alpha1[-0.0874,-1.,1.])\");//only expenential",
                  "w->factory(\"Uniform::background2(x)\");//Uniforn",
                  "w->factory(\"Voigtian:ggH_hjj(x,CMS_H_mean,CMS_H_width,CMS_H_sigma1[0.149])\");",
                  "w->factory(\"SUM::background(background1, background2)\");",
                  "w->factory(\"SUM::model_s(nB[0,8000]*background, nS[0,500]*ggH_hjj)\");",
                  "w->factory(\"SUM::model_b(nB*background)\");",
                  "RooArgSet obs(*w->var(\"x\"));",
                  "//read data",
                  "RooDataSet *data_s = new RooDataSet (\"data_s\",\"data_s\",obs,0);",
                  "data_s =  RooDataSet::read(\"BrazilianPlotData.dat\",obs);",
                  "RooDataSet *data_b = new RooDataSet (\"data_b\",\"data_b\",obs,0);",
                  "data_b =  RooDataSet::read(\"BrazilianPlotData.dat\",obs);",
                  "RooFitResult *r1 = w->pdf(\"model_s\")->fitTo(*data_b, Extended(kTRUE), RooFit::Save(true),RooFit::Minimizer(\"Minuit2\",\"Migrad\"));",  
                  "r1->Print();",
                  "//RooDataSet *data_s = w->pdf(\"model_s\")->generate(obs,Extended());",
                  "//RooDataSet *data_b = w->pdf(\"model_b\")->generate(obs,Extended());",
                  "RooPlot *frame = w->var(\"x\")->frame();",
                  "data_s->plotOn(frame);",
                  "w->pdf(\"model_s\")->plotOn(frame, LineColor(kRed));",
                  "w->pdf(\"model_s\")->plotOn(frame, Components(\"background\"));",
                  "frame->Draw();",
                  "c1->Print(\"data_s.png\");",
                  "frame = w->var(\"x\")->frame();",
                  "data_b->plotOn(frame);",
                  "w->pdf(\"model_b\")->plotOn(frame);",
                  "frame->Draw();",
                  "c1->Print(\"data_b.png\");",
                  "RooDataHist *bdata_b = new RooDataHist(\"data_obs\", \"\", obs, *data_b);",
                  "RooDataHist *bdata_s = new RooDataHist(\"data_sig\", \"\", obs, *data_s);",
                  "// ------------ Make histograms ---------------------------",
                  "TFile *allHistFile = new TFile(\"simple-shapes-TH1.root\", \"RECREATE\");",
                  "// Signal model",
                  "TH1 *signal_nominal = w->pdf(\"ggH_hjj\")->createHistogram(\"x\");",
                  "float nS1=w->var(\"nS\")->getVal();",
                  "float nB1=w->var(\"nB\")->getVal();;",
                  "float alpha=w->var(\"CMS_hjj_alpha\")->getVal();",
                  "cout<<\"alpha--->\"<<alpha<<endl;",
                  "float sigma1 = w->var(\"CMS_H_sigma1\")->getVal();",
                  "signal_nominal->SetName(\"ggH_hjj\");", 
                  "signal_nominal->Scale(nS1/signal_nominal->Integral());",
                  "sigma1 = w->var(\"CMS_H_sigma1\")->getVal();",
                  "sigma1=sigma1*1.1;",
                  "w->var(\"CMS_H_sigma1\")->setVal(sigma1);",
                  "w->var(\"CMS_H_mean\")->setVal(i);",
                  "TH1 *signal_widthUp = w->pdf(\"ggH_hjj\")->createHistogram(\"x\");",
                  "signal_widthUp->SetName(\"signal_widthUp\"); signal_widthUp->Scale(nS1/signal_widthUp->Integral());",
                  "sigma1 = w->var(\"CMS_H_sigma1\")->getVal();",
                  "sigma1=sigma1*0.85;",
                  "w->var(\"CMS_H_sigma1\")->setVal(sigma1);",
                  "TH1 *signal_widthDown = w->pdf(\"ggH_hjj\")->createHistogram(\"x\");",
                  "signal_widthDown->SetName(\"signal_widthDown\"); signal_widthDown->Scale(nS1/signal_widthDown->Integral());",
                  "w->var(\"CMS_H_sigma1\")->setVal(sigma1);",
                  "c1->Clear();",
                  "signal_widthDown->Draw(\"H\"); signal_widthDown->SetLineColor(kBlue); signal_widthDown->SetLineWidth(2);",
                  "signal_widthUp->Draw(\"H SAME\"); signal_widthUp->SetLineColor(kRed); signal_widthUp->SetLineWidth(2);",
                  "signal_nominal->Draw(\"H SAME\"); signal_nominal->SetLineColor(kBlack); signal_nominal->SetLineWidth(3);",
                  "c1->Print(\"signal_model_binned.png\");",
                  "// background model",
                  "frame = w->var(\"x\")->frame();",
                  "TH1 *background_nominal = w->pdf(\"background\")->createHistogram(\"x\");",
                  "background_nominal->SetName(\"background\"); background_nominal->Scale(nB1/background_nominal->Integral());",
                  "alpha=w->var(\"CMS_hjj_alpha\")->getVal();",
                  "float alpha1=alpha*0.90;",
                  "w->var(\"CMS_hjj_alpha\")->setVal(alpha1);",
                  "TH1 *background_alphaUp = w->pdf(\"background\")->createHistogram(\"x\");",
                  "background_alphaUp->SetName(\"background_alphaUp\"); background_alphaUp->Scale(nB1*1.15/background_alphaUp->Integral());",
                  "w->pdf(\"background\")->plotOn(frame, LineColor(kRed), LineWidth(2), Normalization(1.15));",
                  "alpha2=alpha*1.1;"
                  "w->var(\"CMS_hjj_alpha\")->setVal(alpha2);",
                  "TH1 *background_alphaDown = w->pdf(\"background\")->createHistogram(\"x\");",
                  "background_alphaDown->SetName(\"background_alphaDown\"); background_alphaDown->Scale(nB1*0.90/background_alphaDown->Integral());",
                  "w->pdf(\"background\")->plotOn(frame, LineColor(kBlue), LineWidth(2), Normalization(0.90));",
                  "w->var(\"CMS_hjj_alpha\")->setVal(alpha1);",
                  "w->pdf(\"background\")->plotOn(frame, LineColor(kBlack), LineWidth(3), Normalization(1.0));",
                  "frame->Draw(); c1->Print(\"background_model_unbinned.png\");",
                  "background_alphaDown->Draw(\"H\"); background_alphaDown->SetLineColor(kBlue); background_alphaDown->SetLineWidth(2);",
                  "background_alphaUp->Draw(\"H SAME\"); background_alphaUp->SetLineColor(kRed); background_alphaUp->SetLineWidth(2);",
                  "background_nominal->Draw(\"H SAME\"); background_nominal->SetLineColor(kBlack); background_nominal->SetLineWidth(3);",
                  "c1->Print(\"background_model_binned.png\");",
                  "// data ", 
                  "TH1 *hdata_b = bdata_b->createHistogram(\"x\"); hdata_b->SetName(\"data_obs\");",
                  "TH1 *hdata_s = bdata_s->createHistogram(\"x\"); hdata_s->SetName(\"data_sig\");",
                  "// write to file",
                  "allHistFile->WriteTObject(signal_nominal); allHistFile->WriteTObject(signal_widthUp); allHistFile->WriteTObject(signal_widthDown);",
                  "allHistFile->WriteTObject(background_nominal); allHistFile->WriteTObject(background_alphaUp); allHistFile->WriteTObject(background_alphaDown);",
                  "allHistFile->WriteTObject(hdata_b);",
                  "//allHistFile->WriteTObject(hdata_s);",
                  "// ------------ Make RooFit histograms ----------------------------------",
                  "RooWorkspace *wB = new RooWorkspace(\"w\",\"w\");",
                  "RooArgList hobs(*w->var(\"x\"));",
                  "wB->import(*bdata_b);",
                  "wB->import(*bdata_s);",
                  "RooDataHist *hsignal = new RooDataHist(\"hsignal\",\"\",hobs,signal_nominal);",
                  "RooDataHist *hsignal_widthUp = new RooDataHist(\"hsignal_widthUp\",\"\",hobs,signal_widthUp);",
                  "RooDataHist *hsignal_widthDown = new RooDataHist(\"hsignal_widthDown\",\"\",hobs,signal_widthDown);",
                  "RooDataHist *hbackground = new RooDataHist(\"hbackground\",\"\",hobs,background_nominal);",
                  "RooDataHist *hbackground_alphaUp = new RooDataHist(\"hbackground_alphaUp\",\"\",hobs,background_alphaUp);",
                  "RooDataHist *hbackground_alphaDown = new RooDataHist(\"hbackground_alphaDown\",\"\",hobs,background_alphaDown);",
                  "wB->import(*(new RooHistPdf(\"ggH_hjj\",\"\",obs,*hsignal)));",
                  "wB->import(*(new RooHistPdf(\"ggH_hjj_widthUp\",\"\",obs,*hsignal_widthUp)));",
                  "wB->import(*(new RooHistPdf(\"ggH_hjj_widthDown\",\"\",obs,*hsignal_widthDown)));",
                  "wB->import(*(new RooHistPdf(\"background\",\"\",obs,*hbackground)));",
                  "wB->import(*(new RooHistPdf(\"background_alphaUp\",\"\",obs,*hbackground_alphaUp)));",
                  "wB->import(*(new RooHistPdf(\"background_alphaDown\",\"\",obs,*hbackground_alphaDown)));",
                  "wB->writeToFile(\"simple-shapes-RooDataHist.root\");",
                  "RooWorkspace *wBP = new RooWorkspace(\"w\",\"w\");",
                  "wBP->import(*bdata_b);",
                  "wBP->import(*bdata_s);",
                  "wBP->import(*w->pdf(\"ggH_hjj\"));",
                  "wBP->import(*w->pdf(\"background\"));",
                  "wBP->writeToFile(\"simple-shapes-BinnedParam.root\");",
                  "//w->var(\"nS\")->setVal(10);",
                  "float NBkg = nB1;",
                  "RooWorkspace *wUP = new RooWorkspace(\"w \",\"w\");",
                  "wUP->var(\"x[8,40]\");",
                  "RooRealVar BkgNorm(\"background_norm\",\" \",NBkg, NBkg*0.25, NBkg*1.75);",
                  "cout <<\"Background Norm--->\" << BkgNorm.getVal()<<endl;",
                  "wUP->import(BkgNorm);",
                  "wUP->import(*data_b, Rename(\"data_obs\"));",
                  "//wUP->import(*data_s, Rename(\"data_sig\"));",
                  "wUP->import(*w->pdf(\"ggH_hjj\"));",
                  "wUP->import(*w->pdf(\"background\"));",
                  "wUP->writeToFile(\"simple-shapes-UnbinnedParam.root\");",
                  "// now we make a version in which the alpha is function of a unit gaussian,",
                  "// so that we can do normalization and parametric morphing together",
                  "//RooWorkspace *wUPN = new RooWorkspace(\"w\",\"w\");",
                  "//wUPN->var(\"x[8,40]\");",
                  "//wUPN->import(*data_b, Rename(\"data_obs\"));",
                  "//wUPN->import(*data_s, Rename(\"data_sig\"));",
                  "//wUPN->import(*w->pdf(\"ggH_hjj\"));",
                  "//RooAbsPdf *bgpdf = (RooAbsPdf *) w->pdf(\"background\")->clone(\"background_norm\");",
                  "//wUPN->factory(\"sum::param_alpha(alpha,prod(alphaNorm[0],alpha2))\");",                   
                  "// param_alpha = -0.0874 + 0.1*alphaNorm, so Gauss(-0.3,1)",
                  "//wUPN->factory(\"EDIT::background(background_norm, CMS_hjj_alpha=param_alpha)\");",
                  "//wUPN->Print(\"V\");",
                  "//wUPN->writeToFile(\"simple-shapes-UnbinnedParamNorm.root\");",
                  "gApplication -> Terminate(0);",
                  "}"
              ]

    for label, theta_B in zip(labels,values):
        make = open("makefile_"+label+".C", 'w')
        for line in make_lines:
            if line == "using namespace RooFit;":
                make.write(line+"\n")
                make.write("Float_t i = "+ label+";"+"\n")
                make.write("Float_t reso;"+"\n")
                make.write("reso = 80.31 + (2.01 * (i/10.)) + (0.1045*(i/10.)*(i/10.));"+"\n")
            elif line =="void ":
                 make.write(line)
                 make.write("makefile_"+label)
            else:
                 make.write(line+"\n")
        make.close()
        print ">>>   makefile_"+label+".C created"
    return labels

# EXECUTE datacards
def executeDataCards(labels):

    for label in labels:
        #file_name = "datacard_"+label+".txt"#Himal add
        file_name = "hjj_preapprove.txt"
        make_name = "makefile_"+label+".C"
        root_command = "root -b makefile_"+label+".C"
        combine_command = "combine -M Asymptotic -m 125 -n %s %s" % (label,file_name)
        print ""
        print ">>> " + root_command
        q = subprocess.Popen(root_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in q.stdout.readlines():
            print line.rstrip("\n")
        print ">>>   RootFile"+label+"***.root created"
        retval = q.wait()
        print ">>> " + combine_command
        p = subprocess.Popen(combine_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print line.rstrip("\n")
        print ">>>   higgsCombine"+label+".Asymptotic.mH125.root created"
        retval = p.wait()


# GET limits from root file
def getLimits(file_name):

    file = TFile(file_name)
    tree = file.Get("limit")

    limits = [ ]
    for quantile in tree:
        limits.append(tree.limit)
        print ">>>   %.2f" % limits[-1]

    return limits[:6]


# PLOT upper limits
def plotUpperLimits(labels,values):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
    
    N = len(labels)
    yellow = TGraph(2*N)    # yellow band
    green = TGraph(2*N)     # green band
    median = TGraph(N)      # median line
    observed = TGraph(N)    #observed
    up2s = [ ]
    for i in range(N):
        file_name = "higgsCombine"+labels[i]+".Asymptotic.mH125.root"
        limit = getLimits(file_name)
        up2s.append(limit[5])
        yellow.SetPoint(    i,    values[i], limit[4] ) # + 2 sigma
        green.SetPoint(     i,    values[i], limit[3] ) # + 1 sigma
        median.SetPoint(    i,    values[i], limit[2] ) # median
        green.SetPoint(  2*N-1-i, values[i], limit[1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, values[i], limit[0] ) # - 2 sigma
        observed.SetPoint(    i,    values[i], limit[5] ) # observed
    W = 800
    H  = 600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.cd()
    frame = c.DrawFrame(1.4,0.001, 4.1, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("95% upper limit on signal yield")
#    frame.GetYaxis().SetTitle("95% upper limit on #sigma #times BR / (#sigma #times BR)_{SM}")
    frame.GetXaxis().SetTitle("m_{4#mu} [GeV]")
    frame.SetMinimum(0)
    frame.SetMaximum(max(up2s)*1.20)
    frame.GetXaxis().SetLimits(min(values),max(values))

    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    yellow.Draw('F')
    
    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
    green.Draw('Fsame')

    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw('Lsame')
    
    observed.SetLineColor(1)
    observed.SetLineWidth(2)
    observed.SetLineStyle(1)
    observed.Draw('Lsame')

    CMS_lumi.CMS_lumi(c,4,11)#himal block
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')

    x1 = 0.55
    x2 = x1 + 0.24
    y2 = 0.90
    y1 = 0.65
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    legend.AddEntry(median, "Asymptotic CL_{s} Expected",'L')
    legend.AddEntry(observed, "Asymptotic CL_{s} Observed",'L')
    legend.AddEntry(green, "#pm 1 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 1 std. deviation",'f')
    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 2 std. deviation",'f')
    legend.Draw()

    print " "
    c.SaveAs("UpperLimit.png")
    c.Close()


# RANGE of floats
def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step


# MAIN
def main():

    labels = [ ]
    values = [ ]
    for theta_B in frange(12.1,39.9,0.1):
        values.append(theta_B)
        label = "%d" % (theta_B*10)
        labels.append(label)

#    createMakeFileThetaB(labels,values) #Himal block
#    executeDataCards(labels)
    plotUpperLimits(labels,values)



if __name__ == '__main__':
    main()
    
