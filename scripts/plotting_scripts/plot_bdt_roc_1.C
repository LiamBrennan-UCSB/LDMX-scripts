#include "EffPlotTools.hh"
#include "StyleTools.hh"
#include "TFile.h"
#include "TMarker.h"
#include "TSystem.h"
#include "TTree.h"

#include <sstream>

using namespace EffPlotTools;
using namespace StyleTools;

// Macro to plot ROC curves
void plot_bdt_roc_1(const TString &bdtName1, const TString &bdtName2, const TString &flatDir1,
                    const TString &flatDir2, const TString &outDir, const TString &withSelect,
                    const double &bkgEff, const bool &isLog, const bool &isZoom) {

    // Names for each process
    vector<TString> procs = {"0.001", "0.01", "0.1", "1.0", "bkg"};

    // Labels for each process
    map<TString, TString> proclabel;
    proclabel["0.001"] = "m_{A'} = 0.001 GeV";
    proclabel["0.01"] = "m_{A'} = 0.01 GeV";
    proclabel["0.1"] = "m_{A'} = 0.1 GeV";
    proclabel["1.0"] = "m_{A'} = 1 GeV";
    proclabel["bkg"] = "Photonuclear";

    // Colors for each process
    map<TString, unsigned int> colors;
    colors["0.001"] = color_comp2;
    colors["0.01"] = color_comp3;
    colors["0.1"] = color_comp4;
    colors["1.0"] = color_comp5;
    colors["bkg"] = color_comp1;

    // Histograms for each process
    TH1D* bkghist1 = 0;
    TH1D* bkghist2 = 0;
    vector<TH1D*> sighists1;
    vector<TH1D*> sighists2;

    // Define each selection
    map<TString, TString> sel;
    sel["base"] = "1 == 1";

    // Make nice looking plots
    SetTDRStyle();
    gSystem->mkdir(outDir, true);

    // Loop over processes
    for(auto proc : procs) {

        // Read trees from input files
        TFile* file1 = new TFile(flatDir1 + "/" + proc + "_" + bdtName1 + "_eval.root");
        TTree* tree1 = (TTree*)file1->Get("EcalVeto");
        assert(tree1);
        TFile* file2 = new TFile(flatDir2 + "/" + proc + "_" + bdtName2 + "_eval.root");
        TTree* tree2 = (TTree*)file2->Get("EcalVeto");
        assert(tree2);

        // Set up histograms to be filled with the disc values
        TH1D* hist1 = new TH1D("h1_discValue_" + bdtName1 + "_" + withSelect + "_" + proc, "", 1e04, 0, 1);
        TH1D* hist2 = new TH1D("h2_discValue_" + bdtName2 + "_" + withSelect + "_" + proc, "", 1e04, 0, 1);

        // Read off disc values from trees into the histograms
        // tree1->Draw("discValue_fernand>>h1_discValue_" + bdtName1 + "_" + withSelect + "_" + proc, sel[withSelect]); // Old format
        tree1->Draw("discValue_" + bdtName1 + ">>h1_discValue_" + bdtName1 + "_" + withSelect + "_" + proc, sel[withSelect]); // New format
        // tree2->Draw("discValue_fernand>>h2_discValue_" + bdtName2 + "_" + withSelect + "_" + proc, sel[withSelect]); // Old format
        tree2->Draw("discValue_" + bdtName2 + ">>h2_discValue_" + bdtName2 + "_" + withSelect + "_" + proc, sel[withSelect]); // New format

        // Store histograms
        if(proc == "bkg") {
            bkghist1 = hist1;
            bkghist2 = hist2;
        }
        else {
            sighists1.push_back(hist1);
            sighists2.push_back(hist2);
        }
    }

    // Get background efficiencies for 0.99 cut
    float bkgeff1 = getEfficiencyForCutValue(bkghist1, 0.99)[0];
    float bkgeff2 = getEfficiencyForCutValue(bkghist2, 0.99)[0];
    cout << "Nominal efficiency for cutval = 0.99: eff(bkg) = " << bkgeff1;
    cout << "\tNew efficiency for cutval = 0.99: eff(bkg) new = " << bkgeff2;

    // Get cut values for selected background efficiency
    ostringstream ss;
    ss << bkgEff;
    string effstr = ss.str();
    float cutval1 = getCutValueForEfficiency(bkghist1, bkgEff)[0];
    float cutval2 = getCutValueForEfficiency(bkghist2, bkgEff)[0];
    cout << "\nNominal cut value for eff(bkg) = " + effstr + ": cutval =  " << cutval1;
    cout << "\tNew cut value for eff(bkg) = " + effstr + ": cutval =  " << cutval2;

    // Set up the canvas and legend
    TCanvas* c1 = MakeCanvas("c1", "", 600, 600);
    TLegend* leg = new TLegend(0.5, 0.2, 0.8, 0.5);
    SetLegendStyle(leg);
    leg->SetTextSize(0.03);

    // Loop over signals and compute ROC curve for each
    int isig = 0;
    for(auto* sighist1 : sighists1) {
        TH1D* sighist2 = sighists2[isig];
        TString hname(sighist1->GetName());
        TString procname = TString(hname(hname.Last('_') + 1, hname.Length()));

        // Make the ROC curve
        TGraph* rocgr1 = computeROCCurve(sighist1, bkghist1, "", false, false, true);
        TGraph* rocgr2 = computeROCCurve(sighist2, bkghist2, "", false, false, true);

        // Get signal efficiencies for 0.99 cut
        float sigeff1 = getEfficiencyForCutValue(sighist1, 0.99)[0];
        float sigeff2 = getEfficiencyForCutValue(sighist2, 0.99)[0];
        cout << "\nNominal efficiency for cutval = 0.99: eff(" << procname << ") = " << sigeff1;
        cout << "\tNew efficiency for cutval = 0.99: eff(" << procname << ") = " << sigeff2;

        // Get signal efficiencies for cut value above
        cout << "\nNominal efficiency for cutval = " << cutval1 << ": eff(" << procname << ") = " << getEfficiencyForCutValue(sighist1, cutval1)[0];
        cout << "\tNew efficiency for cutval = " << cutval2 << ": eff(" << procname << ") = " << getEfficiencyForCutValue(sighist2, cutval2)[0];

        // Set up colors and style
        rocgr1->SetLineColor(colors[procname]);
        rocgr1->SetMarkerColor(colors[procname]);
        rocgr1->SetLineWidth(3);
        rocgr1->SetLineStyle(2);
      
        rocgr2->SetLineColor(colors[procname]);
        rocgr2->SetMarkerColor(colors[procname]);
        rocgr2->SetLineWidth(3); 

        // Add this signal to the legend
        leg->AddEntry(rocgr1, proclabel[procname] + ", " + bdtName1, "L");
        leg->AddEntry(rocgr2, proclabel[procname] + ", " + bdtName2, "L");

        // Set axis titles
        rocgr1->GetXaxis()->SetTitle("#varepsilon(bkg)");
        rocgr1->GetYaxis()->SetTitle("#varepsilon(sig)");

        // Zoom in on the plot or set log scale if we want to
        if(isZoom) rocgr1->GetXaxis()->SetLimits(1e-04, 5e-03);
        if(isZoom) rocgr1->GetYaxis()->SetRangeUser(0.5, 1);
        if(isZoom) rocgr2->GetXaxis()->SetLimits(1e-04, 5e-03);
        if(isZoom) rocgr2->GetYaxis()->SetRangeUser(0.5, 1);
        if(isLog) c1->SetLogx();

        // If it's the first process, draw the axis as well
        c1->cd();
        if(isig == 0)
	    rocgr1->Draw("AC");
        else rocgr1->Draw("Csame");
        rocgr2->Draw("Csame");

        isig++;
    }

    // Draw the legend
    c1->cd();
    leg->Draw("same");

    // Draw LDMX text
    LDMX_lumi(c1, 0, "Simulation");

    // Save the plot
    if(isZoom and isLog) c1->SaveAs(outDir + "/" + bdtName1 + "_" + bdtName2 + "_roc_" + withSelect + "_zoom_log.pdf");
    else if(isZoom and not isLog) c1->SaveAs(outDir + "/" + bdtName1 + "_" + bdtName2 + "_roc_" + withSelect + "_zoom.pdf");
    else if(not isZoom and isLog) c1->SaveAs(outDir + "/" + bdtName1 + "_" + bdtName2 + "_roc_" + withSelect + "_log.pdf");
    else if(not isZoom and not isLog) c1->SaveAs(outDir + "/" + bdtName1 + "_" + bdtName2 + "_roc_" + withSelect + ".pdf");
    c1->Close();
}

