// ---- ML a 2 parametri

#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TMinuit.h"

using namespace std;

// dati globali
vector<double> da;
TF1 *f_global;


// ----- Definizione -logL (da miminizzare)
void logL(int &npar, double *grad, double &fval, double *par, int iflag) {
    double alfa = par[0];
    double beta = par[1];

    f_global->SetParameter(0, alfa);
    f_global->SetParameter(1, beta);

    fval = 0;
    for (size_t k = 0; k < da.size(); k++) {
        double val = f_global->Eval(da[k]);
        if(val > 0) fval -= log(val);
    }
}

void ml() {
    TFile *file = new TFile("output.root", "RECREATE");

    // ----- Definizione pdf
    double xmin = -0.95, xmax = 0.95;
    f_global = new TF1("f",
           "(1+[0]*x+[1]*x*x)/(([3]-[2])+([0]/2)*([3]*[3]-[2]*[2])+([1]/3)*(pow([3],3)-pow([2],3)))",
           xmin, xmax);
    f_global->SetParameter(2, xmin);
    f_global->SetParameter(3, xmax);
    f_global->FixParameter(2, xmin);
    f_global->FixParameter(3, xmax);

    // ----- Simulazione dati
    double alfa = 0.5, beta = 0.5;
    f_global->SetParameter(0, alfa);
    f_global->SetParameter(1, beta);

    TTree *dd = new TTree("dd", "data");
    double e;
    dd->Branch("E", &e, "E/D");

    for(int i = 0; i < 20000; i++) {
        e = f_global->GetRandom();
        dd->Fill();
        da.push_back(e); // salvo anche nel vettore globale
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

    cout << endl << endl << endl << "Parametri ottimali stimati:" << endl;
    cout << "alfa= " << alfa_fit << " +- " << err0 << " (atteso: 0.5)" << endl;
    cout << "beta= " << beta_fit << " +- " << err1 << " (atteso: 0.5)" << endl << endl << endl;

    file->Write();
    file->Close();
}
