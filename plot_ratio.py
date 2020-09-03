##-----------------------------------------
## Compare histograms of reweighted dists
## Zishuo Yang, 2020-09-02
## ----------------------------------------
import ROOT as rt

# Set output filepath
filepath = "./plots/"

# Get input
inputFile1 = rt.TFile("./HammerWeights.root")
inputFile2 = rt.TFile("./Bc2JpsiMuNu_subset.root")
tHam = inputFile1.Get("Hammer")
tTup = inputFile2.Get("DecayTree")

# Add friend to associate events
tTup.BuildIndex("runNumber","eventNumber")
tHam.AddFriend(tTup)

def draw_ratio( var, bin_range, title, filename,
		upperYmin, upperYmax,
                lowerYmin, lowerYmax):
    
    rt.gStyle.SetOptStat(0)
    c1 = rt.TCanvas("c1", "A ratio plot")

    tHam.Draw( var + ">>h1" + bin_range, "", "goff", 50000,20000)
    h1 = rt.gDirectory.Get("h1")
    h1.SetMarkerColor(rt.kBlue)
    h1.SetLineColor(rt.kBlue)
    
    tHam.Draw( var + ">>h2" + bin_range, "HammerWeight", "goff", 50000,20000)
    h2 = rt.gDirectory.Get("h2")
    h2.SetMarkerColor(rt.kRed)
    h2.SetLineColor(rt.kRed)
    
    rp = rt.TRatioPlot(h1, h2)
    rp.Draw()
    rp.GetLowerRefYaxis().SetRangeUser(lowerYmin, lowerYmax)
    rp.GetLowerRefXaxis().SetTitle(title)
    rp.GetUpperRefYaxis().SetRangeUser(upperYmin, upperYmax)
    c1.Update()

    c1.Print(filepath+filename+".png")



draw_ratio("FitVar_q2", "(50,0,12E6)", "q2 [MeV^2]", "q2", 0, 2000, 0.5, 1.2)
draw_ratio("FitVar_El", "(80,0,3.5E3)", "El* [MeV]", "El", 0, 2500, 0.5, 1.2)
 
## //TLegend *leg = new TLegend(0.55,0.70,0.90,0.90,NULL,"brNDC");
## //TLegendEntry *entry=leg->AddEntry("h1","unweighted","lpf");
## //entry=leg->AddEntry("h2","weighted","lpflpf");
## //leg->Draw();
