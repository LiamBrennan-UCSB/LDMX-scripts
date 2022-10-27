#include "EffPlotTools.hh"
#include "StyleTools.hh"
#include "TFile.h"
#include "TMarker.h"
#include "TSystem.h"
#include "TTree.h"

#include <sstream>

using namespace EffPlotTools;
using namespace StyleTools;

void plot_bdt_roc_0(const TString &bdtName, const TString &flatDir, const TString &outDir,
                    const TString &withSelect, const double &bkgEff, const int &isLog,
                    const int &isZoom) {

    // Name for each process
    vector<TString> procs = {"0.001", "0.01", "0.1", "1.0", "bkg"};

    // Label for each process
    map<TString, TString> proclabel;
    proclabel["0.001"] = "m_{A'} = 0.001 GeV";
    proclabel["0.01"] = "m_{A'} = 0.01 GeV";
    proclabel["0.1"] = "m_{A'} = 0.1 GeV";
    proclabel["1.0"] = "m_{A'} = 1 GeV";
    proclabel["bkg"] = "Photonuclear";

    // Color for each process
    map<TString, unsigned int> colors;
    colors["0.001"] = color_comp2;
    colors["0.01"] = color_comp3;
    colors["0.1"] = color_comp4;
    colors["1.0"] = color_comp5;
    colors["bkg"] = color_comp1;

    // Histogram for each process
    TH1D* bkghist = 0;
    vector<TH1D*> sighists;

    // Define each selection
    map<TString, TString> sel;
    sel["base"] = "1 == 1";

    // Make nice looking plots
    SetTDRStyle();
    gSystem->mkdir(outDir, true);
    
    // Loop over processes
    for(auto proc : procs) {

        // Read tree from input file
        TFile* file = new TFile(flatDir + "/" + proc + "_" + bdtName + "_eval.root");
        TTree* tree = (TTree*)file->Get("EcalVeto");
        assert(tree);

        // Set up histogram to be filled with the disc values
        TH1D* hist = new TH1D("h_discvalue_" + bdtName + "_" + withSelect + "_" + proc, "", 1e04, 0, 1);

        // Read off disc values from tree into the histogram
        // tree->Draw("discValue_fernand>>h_discvalue_" + bdtName + "_" + withSelect + "_" + proc, sel[withSelect]); // Old format
        tree->Draw("discValue_" + bdtName + ">>h_discvalue_" + bdtName + "_" + withSelect + "_" + proc, sel[withSelect]); // New format

        // Store histogram
        if(proc == "bkg")
            bkghist = hist;
        else
            sighists.push_back(hist);
    }

    // Get the background efficiency for 0.99 cut
    float bkgeff = getEfficiencyForCutValue(bkghist, 0.99)[0];
    cout << "Efficiency for cutval = 0.99: eff(bkg) = " << bkgeff;

    // Get the cut value for the selected background efficiency
    ostringstream ss;
    ss << bkgEff;
    string effstr = ss.str();
    float cutval = getCutValueForEfficiency(bkghist, bkgEff)[0];
    cout << "\tCut value for eff(bkg) = " + effstr + ": cutval = " << cutval;

    // Set up the canvas and legend
    TCanvas* c1 = MakeCanvas("c1", "", 600, 600);
    TLegend* leg = new TLegend(0.6, 0.2, 0.9, 0.4);
    SetLegendStyle(leg);
    leg->SetTextSize(0.03);

    // Loop over signals and compute the ROC curve for each
    int isig = 0;
    for(auto* sighist : sighists) {
        TString hname(sighist->GetName());
        TString procname = TString(hname(hname.Last('_') + 1, hname.Length()));

        // Make the ROC curve
        TGraph* rocgr = computeROCCurve(sighist, bkghist, "", false, false, true);

        //Get the signal efficiency for 0.99 cut
        float sigeff = getEfficiencyForCutValue(sighist, 0.99)[0];
        cout << "\nEfficiency for cutval = 0.99: eff(" << procname << ") = " << sigeff; 

        // Get the signal efficiency for the cut value above
        cout << "\tEfficiency for cutval = " << cutval << ": eff(" << procname << ") = " << getEfficiencyForCutValue(sighist, cutval)[0];

        // Set up colors and style
        rocgr->SetLineColor(colors[procname]);
        rocgr->SetMarkerColor(colors[procname]);
        rocgr->SetLineWidth(3);

        // Add this signal to the legend
        leg->AddEntry(rocgr, proclabel[procname], "L");

        // Set axis titles
        rocgr->GetXaxis()->SetTitle("#varepsilon(bkg)");
        rocgr->GetYaxis()->SetTitle("#varepsilon(sig)");

        // Zoom in on the plot or set log scale if we want to
        if(isZoom) rocgr->GetXaxis()->SetLimits(1e-04, 5e-03);
        if(isZoom) rocgr->GetYaxis()->SetRangeUser(0.5, 1);
        if(isLog) c1->SetLogx();

        // If it's the first process, draw the axis as well
        c1->cd();
        if(isig == 0)
            rocgr->Draw("AC");
        else rocgr->Draw("Csame");

        isig++;
    }

    // Draw the legend
    c1->cd();
    leg->Draw("same");

    // Draw LDMX text
    LDMX_lumi(c1, 0, "Simulation");

    // Save the plot
    if(isZoom and isLog) c1->SaveAs(outDir + "/" + bdtName + "_roc_" + withSelect + "_zoom_log.pdf");
    else if(isZoom and not isLog) c1->SaveAs(outDir + "/" + bdtName + "_roc_" + withSelect + "_zoom.pdf");
    else if(not isZoom and isLog) c1->SaveAs(outDir + "/" + bdtName + "_roc_" + withSelect + "_log.pdf");
    else if(not isZoom and not isLog) c1->SaveAs(outDir + "/" + bdtName + "_roc_" + withSelect + ".pdf");
    c1->Close();
}

