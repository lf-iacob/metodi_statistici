// ---- ML a 2 parametri
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

using namespace std;

// dati globali
vector<double> da;
TF1 *f_global;

// ----- Definizione -logL (da minimizzare)
void logL(int &npar, double *grad, double &fval, double *par, int iflag) {
    double alfa = par[0];
    double beta = par[1];

    f_global->SetParameter(0, alfa);
    f_global->SetParameter(1, beta);

    fval = 0;
    for (size_t k = 0; k < da.size(); k++) {
        double val = f_global->Eval(da[k]);
        if (val > 0) fval -= log(val);
    }
}

void ml() {

    TFile *file = new TFile("output.root", "RECREATE");

    // ----- Def pdf
    double xmin = -0.95, xmax = 0.95;
    f_global = new TF1("f",
        "(1 + [0]*x + [1]*x*x)/"
        "(([3]-[2]) + ([0]/2)*([3]*[3]-[2]*[2]) + ([1]/3)*(pow([3],3)-pow([2],3)))",
        xmin, xmax);
    f_global->FixParameter(2, xmin);
    f_global->FixParameter(3, xmax);

    // ----- Esperimento (simulazione MC)
    double alfa_true = 0.5, beta_true = 0.5;
    f_global->SetParameter(0, alfa_true);
    f_global->SetParameter(1, beta_true);

    TTree *dd = new TTree("dd", "data");
    double e;
    dd->Branch("E", &e, "E/D");

    for (int i = 0; i < 20000; i++) {
        e = f_global->GetRandom();
        dd->Fill();
        da.push_back(e);
    }

    // ----- Minimizzazione
    TMinuit minuit(2);
    minuit.SetFCN(logL);
    minuit.DefineParameter(0, "alfa", 0.0, 0.1, -1.0, 1.0);
    minuit.DefineParameter(1, "beta", 0.0, 0.1, -1.0, 1.0);
    minuit.Migrad();
    double alfa_fit, beta_fit, err0, err1;
    minuit.GetParameter(0, alfa_fit, err0);
    minuit.GetParameter(1, beta_fit, err1);
    cout << endl << "Parametri ottimali stimati:" << endl;
    cout << "alfa= " << alfa_fit << " +- " << err0 << " (atteso: 0.5)" << endl;
    cout << "beta= " << beta_fit << " +- " << err1 << " (atteso: 0.5)" << endl << endl;

    // ----- Rappresentazione grafica
    TH1D *hdata = new TH1D("hdata", "Dati simulati", 100, xmin, xmax);
    for (auto val : da) hdata->Fill(val);
    double N = hdata->GetEntries();
    double binw = hdata->GetBinWidth(1);
    double scale_factor = N * binw;
    TF1 *f_fit = new TF1("f_fit",
        "(1 + [0]*x + [1]*x*x)/"
        "(([3]-[2]) + ([0]/2)*([3]*[3]-[2]*[2]) + ([1]/3)*(pow([3],3)-pow([2],3)))"
        "*[4]",
        xmin, xmax);
    f_fit->SetParameter(0, alfa_fit);
    f_fit->SetParameter(1, beta_fit);
    f_fit->SetParameter(2, xmin);
    f_fit->SetParameter(3, xmax);
    f_fit->SetParameter(4, scale_factor);
    f_fit->SetLineColor(kRed);
    f_fit->SetLineWidth(2);
    f_fit->SetNpx(500);
    int nbins = 200;
    vector<double> x_vals(nbins + 1), y_vals(nbins + 1), y_err(nbins + 1);
    double norm = (xmax - xmin)
                + (alfa_fit/2)*(xmax*xmax - xmin*xmin)
                + (beta_fit/3)*(pow(xmax,3) - pow(xmin,3));
    for (int i = 0; i <= nbins; i++) {
        double x = xmin + i*(xmax - xmin) / nbins;
        double y_pdf = (1 + alfa_fit*x + beta_fit*x*x) / norm;
        double dfdalpha = (x*norm - (1+alfa_fit*x+beta_fit*x*x)*(0.5*(xmax*xmax - xmin*xmin))) / (norm*norm);
        double dfdbeta  = (x*x*norm - (1+alfa_fit*x+beta_fit*x*x)*(1.0/3.0*(pow(xmax,3)-pow(xmin,3)))) / (norm*norm);
        double sigma_pdf = sqrt( pow(dfdalpha*err0,2) + pow(dfdbeta*err1,2) );
        double y = y_pdf * scale_factor;
        double sigma_y = sigma_pdf * scale_factor;
        x_vals[i] = x;
        y_vals[i] = y;
        y_err[i] = sigma_y;
    }
    TGraphErrors *g_band = new TGraphErrors(nbins + 1, &x_vals[0], &y_vals[0], 0, &y_err[0]);
    g_band->SetFillColorAlpha(kAzure-3, 0.35);
    g_band->SetLineColor(kAzure-3);
    TCanvas *c1 = new TCanvas("c1", "Fit vs dati", 800, 600);
    hdata->SetLineColor(kPink+7);
    hdata->Sumw2();
    c1->cd();
    hdata->Draw("");
    g_band->Draw("E3 SAME");
    f_fit->Draw("L SAME");

    c1->Write();


    // ----- Varianza con MC

    

    file->Write();
    //file->Close();
}

