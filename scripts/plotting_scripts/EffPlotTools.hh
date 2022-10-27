#ifndef EFFPLOTTOOLS_HH
#define EFFPLOTTOOLS_HH

#include "TH1D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TString.h"

namespace EffPlotTools {

  TGraphAsymmErrors* computeEffGraph(TH1D* pass, TH1D* total, bool debug=false);
  TGraph* computeEffVsCutGraph(TH1D* signal, bool reversecutdir = false);
  TGraph* computeROCCurve(TH1D* signal, TH1D* background, TString title, bool reversecutdir = false, bool plotbkgrej = false, bool reverseaxes=false);
  TGraph* computeSOvSqrtBGraph(TH1D* signal, TH1D* background, TString title, bool reversecutdir = false);
  TGraph* computeSOvBGraph(TH1D* signal, TH1D* background, TString title, bool reversecutdir = false);
  double * getCutValueForEfficiency(TH1D* hist, double targeteff, bool reversecutdir=false);
  double * getEfficiencyForCutValue(TH1D* hist, double cut, bool reversecutdir=false);
  double * getCutValueForBestEffSOverSqrtEffB(TH1D* sighist, TH1D* bkghist, bool reversecutdir=false);

  TGraphAsymmErrors* computeEffGraph(TH1D* pass, TH1D* total, bool debug) {
  
    // make sure <pass> and <total> have the same binning!
  
    int npoints = total->GetNbinsX();
  
    double x[npoints], y[npoints], errx[npoints], erryl[npoints], erryh[npoints];
  
    double npass = 0.0;
    double ntotal = 0.0;
  
    for(int ibin = 1; ibin < npoints+1; ibin++) {
      x[ibin-1] = total->GetBinCenter(ibin);
      npass = pass->GetBinContent(ibin);
      ntotal = total->GetBinContent(ibin);
      y[ibin-1] = ntotal < 1.0 ? 0.0 : npass/ntotal;
      errx[ibin-1] = 0.5*total->GetBinWidth(ibin);
      if(y[ibin-1]==0.0) {
        erryl[ibin-1] = 0.0; erryh[ibin-1] = 0.0;
      } else {
        if(debug) printf("npass = %3.1f, ntotal = %3.1f, eff = %4.2f\n", npass, ntotal, y[ibin-1]);
        erryl[ibin-1] = y[ibin-1] - TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, false);
        erryh[ibin-1] = TEfficiency::ClopperPearson((unsigned int)ntotal, (unsigned int)npass, 0.683, true) - y[ibin-1];
      }
    }
  
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(npoints, x, y, errx, errx, erryl, erryh);
    gr->SetTitle("");
    gr->GetXaxis()->SetTitle(pass->GetXaxis()->GetTitle());
    gr->GetYaxis()->SetTitle("Efficiency");
  
    return gr;
  
  }
  
  TGraph* computeEffVsCutGraph(TH1D* signal, bool reversecutdir) {
  
    // make sure <pass> and <total> have the same binning!
  
    int nbins = signal->GetNbinsX();
    double xbins[nbins];
    //double binw[nbins];

    for(int ibin = 0; ibin < nbins; ibin++) {
      xbins[ibin] = signal->GetBinLowEdge(ibin+1);
      //binw[ibin] = signal->GetBinWidth(ibin+1);
    }

    //xbins[nbins] = xbins[nbins-1] + binw[nbins-1];

    double npasssig[nbins], effsig[nbins];

    double sigtotal = signal->Integral(0, nbins+1);

    for(int ibin = 0; ibin < nbins; ibin++) {
      if(reversecutdir) {
        npasssig[ibin] = signal->Integral(0, ibin);
      } else {
        npasssig[ibin] = signal->Integral(ibin, nbins+1);
      }
      effsig[ibin] = npasssig[ibin]/sigtotal;
    }

 
    TGraph *gr = new TGraph(nbins, xbins, effsig);
    gr->SetTitle("");
    gr->GetXaxis()->SetTitle(signal->GetXaxis()->GetTitle());
    gr->GetYaxis()->SetTitle("Efficiency");
  
    return gr;
  
  }
  

  TGraph* computeROCCurve(TH1D* signal, TH1D* background, TString title, bool reversecutdir, bool plotbkgrej, bool reverseaxes)
  {

    int nbins = signal->GetNbinsX();
    assert(background->GetNbinsX() == nbins);
    double xbins[nbins+1];
    double binw[nbins];

    for(int ibin = 0; ibin < nbins; ibin++) {
      xbins[ibin] = signal->GetBinLowEdge(ibin+1);
      binw[ibin] = signal->GetBinWidth(ibin+1);
    }

    xbins[nbins] = xbins[nbins-1] + binw[nbins-1];

    double npasssig[nbins], npassbkg[nbins], effsig[nbins], effbkg[nbins];

    double sigtotal = signal->Integral(0, nbins+1);
    double bkgtotal = background->Integral(0, nbins+1);

    for(int ibin = 0; ibin < nbins; ibin++) {
      if(reversecutdir) {
        npasssig[ibin] = signal->Integral(0, ibin);
        npassbkg[ibin] = background->Integral(0, ibin);
      } else {
        npasssig[ibin] = signal->Integral(ibin, nbins+1);
        npassbkg[ibin] = background->Integral(ibin, nbins+1);
      }
      effsig[ibin] = npasssig[ibin]/sigtotal;
      effbkg[ibin] = plotbkgrej ? (1.0 - (npassbkg[ibin]/bkgtotal)) : npassbkg[ibin]/bkgtotal;
    }

    TGraph* roc = 0;
    if(reverseaxes)
      roc = new TGraph(nbins, effbkg, effsig);
    else 
      roc = new TGraph(nbins, effsig, effbkg);
    roc->GetXaxis()->SetLimits(0.0, 1.0);
    roc->GetHistogram()->SetMinimum(0.0);
    roc->GetHistogram()->SetMaximum(1.0);

    roc->SetTitle(title);

    return roc;

  }

  TGraph* computeSOvSqrtBGraph(TH1D* signal, TH1D* background, TString title, bool reversecutdir)
  {

    int nbins = signal->GetNbinsX();
    assert(background->GetNbinsX() == nbins);
    double xbins[nbins];
    double binw[nbins];

    for(int ibin = 0; ibin < nbins; ibin++) {
      xbins[ibin] = signal->GetBinLowEdge(ibin+1) + 0.5*signal->GetBinWidth(ibin+1);
      binw[ibin] = signal->GetBinWidth(ibin+1);
    }

    //xbins[nbins] = xbins[nbins-1] + binw[nbins-1];

    double npasssig[nbins], npassbkg[nbins], effsig[nbins], effbkg[nbins], effsovsqrteffb[nbins];

    double sigtotal = signal->Integral(0, nbins+1);
    double bkgtotal = background->Integral(0, nbins+1);

    for(int ibin = 0; ibin < nbins; ibin++) {
      if(reversecutdir) {
        npasssig[ibin] = signal->Integral(0, ibin);
        npassbkg[ibin] = background->Integral(0, ibin);
      } else {
        npasssig[ibin] = signal->Integral(ibin, nbins+1);
        npassbkg[ibin] = background->Integral(ibin, nbins+1);
      }
      effsig[ibin] = npasssig[ibin]/sigtotal;
      effbkg[ibin] = npassbkg[ibin]/bkgtotal;
      effsovsqrteffb[ibin] = effbkg[ibin] > 0.0 ? effsig[ibin]/sqrt(effbkg[ibin]) : 0.0;
    }

    TGraph* gr = 0;
    gr = new TGraph(nbins, xbins, effsovsqrteffb);
    //gr->GetXaxis()->SetLimits(0.0, 1.0);
    gr->GetHistogram()->SetMinimum(0.0);
    gr->GetHistogram()->SetMaximum(500.0);

    gr->SetTitle(title);

    return gr;

  }

  TGraph* computeSOvBGraph(TH1D* signal, TH1D* background, TString title, bool reversecutdir)
  {

    int nbins = signal->GetNbinsX();
    assert(background->GetNbinsX() == nbins);
    double xbins[nbins];
    double binw[nbins];

    for(int ibin = 0; ibin < nbins; ibin++) {
      xbins[ibin] = signal->GetBinLowEdge(ibin+1) + 0.5*signal->GetBinWidth(ibin+1);
      binw[ibin] = signal->GetBinWidth(ibin+1);
    }

    //xbins[nbins] = xbins[nbins-1] + binw[nbins-1];

    double npasssig[nbins], npassbkg[nbins], effsig[nbins], effbkg[nbins], effsoveffb[nbins];

    double sigtotal = signal->Integral(0, nbins+1);
    double bkgtotal = background->Integral(0, nbins+1);

    for(int ibin = 0; ibin < nbins; ibin++) {
      if(reversecutdir) {
        npasssig[ibin] = signal->Integral(0, ibin);
        npassbkg[ibin] = background->Integral(0, ibin);
      } else {
        npasssig[ibin] = signal->Integral(ibin, nbins+1);
        npassbkg[ibin] = background->Integral(ibin, nbins+1);
      }
      effsig[ibin] = npasssig[ibin]/sigtotal;
      effbkg[ibin] = npassbkg[ibin]/bkgtotal;
      effsoveffb[ibin] = effbkg[ibin] > 0.0 ? effsig[ibin]/effbkg[ibin] : 0.0;
    }

    TGraph* gr = 0;
    gr = new TGraph(nbins, xbins, effsoveffb);
    //gr->GetXaxis()->SetLimits(0.0, 1.0);
    gr->GetHistogram()->SetMinimum(0.0);
    gr->GetHistogram()->SetMaximum(100.0);

    gr->SetTitle(title);

    return gr;

  }

  double * getCutValueForEfficiency(TH1D* hist, double targeteff, bool reversecutdir)
  {

    int nbins = hist->GetNbinsX();
    double xbins[nbins+1];
    double binw[nbins];

    for(int ibin = 0; ibin < nbins; ibin++) {
      xbins[ibin] = hist->GetBinLowEdge(ibin+1);
      binw[ibin] = hist->GetBinWidth(ibin+1);
    }

    xbins[nbins] = xbins[nbins-1] + binw[nbins-1];

    double npass[nbins], eff[nbins];

    double total = hist->Integral(0, nbins+1);

    double effdiff = 1.0;
    int nbin = -1;

    for(int ibin = 0; ibin < nbins; ibin++) {
      if(reversecutdir)
        npass[ibin] = hist->Integral(0, ibin);
      else 
        npass[ibin] = hist->Integral(ibin, nbins+1);
      eff[ibin] = npass[ibin]/total;
      double tmpdiff = fabs(eff[ibin] - targeteff);
      if(tmpdiff < effdiff) {
        effdiff = tmpdiff;
        nbin = ibin;
      }
    }

    static double result[2];
    result[0] = xbins[nbin];
    result[1] = eff[nbin];

    return result;

  }

  double * getEfficiencyForCutValue(TH1D* hist, double cut, bool reversecutdir)
  {

    int nbins = hist->GetNbinsX();
    double xbins[nbins+1];
    double binw[nbins];

    for(int ibin = 0; ibin < nbins; ibin++) {
      xbins[ibin] = hist->GetBinLowEdge(ibin+1);
      binw[ibin] = hist->GetBinWidth(ibin+1);
    }

    xbins[nbins] = xbins[nbins-1] + binw[nbins-1];

    double npass[nbins], eff[nbins];

    double total = hist->Integral(0, nbins+1);

    double diff = 1.0;
    int nbin = -1;

    for(int ibin = 0; ibin < nbins; ibin++) {
      if(reversecutdir)
        npass[ibin] = hist->Integral(0, ibin);
      else 
        npass[ibin] = hist->Integral(ibin, nbins+1);
      eff[ibin] = npass[ibin]/total;
      double tmpdiff = fabs(xbins[ibin] - cut);
      if(tmpdiff < diff) {
        diff = tmpdiff;
        nbin = ibin;
      }
    }

    static double result[2];

    result[0] = eff[nbin];
    result[1] = xbins[nbin];

    return result;

  }

  double * getCutValueForBestEffSOverSqrtEffB(TH1D* sighist, TH1D* bkghist, bool reversecutdir)
  {

    int nbins = sighist->GetNbinsX();
    double xbins[nbins+1];
    double binw[nbins];

    for(int ibin = 0; ibin < nbins; ibin++) {
      xbins[ibin] = sighist->GetBinLowEdge(ibin+1);
      binw[ibin] = sighist->GetBinWidth(ibin+1);
    }

    xbins[nbins] = xbins[nbins-1] + binw[nbins-1];

    double soversqrtb[nbins];

    double maxsoversqrtb = 0.0;
    int nbin = -1;

    double sigtotal = sighist->Integral(0, nbins+1);
    double bkgtotal = bkghist->Integral(0, nbins+1);

    for(int ibin = 0; ibin < nbins; ibin++) {
      if(reversecutdir)
        soversqrtb[ibin] = (sighist->Integral(0, ibin)/sigtotal)/sqrt(bkghist->Integral(0, ibin)/bkgtotal);
      else 
        soversqrtb[ibin] = (sighist->Integral(ibin, nbins+1)/sigtotal)/sqrt(bkghist->Integral(ibin, nbins+1)/bkgtotal);
      if(soversqrtb[ibin] > maxsoversqrtb) {
        maxsoversqrtb = soversqrtb[ibin];
        nbin = ibin;
      }
    }

    static double result[2];
    result[0] = xbins[nbin];
    result[1] = soversqrtb[nbin];

    return result;

  }

};

#endif
