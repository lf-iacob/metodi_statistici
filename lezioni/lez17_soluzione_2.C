// Metodi stima dei parametri a confornto (soluzione 2)

#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
#include "TRandom.h"
#include "TRandom2.h"
#include "TMinuit.h"
#include "TFitResult.h"

using namespace std;

// per la definizione della pdf usata, vedere slides lezione 17
vector<double>* xVecPtr = new vector<double>();
vector<double>* alphaVec_MM = new vector<double>();
vector<double>* betaVec_MM = new vector<double>();
vector<double>* alphaVec_ML = new vector<double>();
vector<double>* betaVec_ML = new vector<double>();

vector<double>* sigma_alphaVec_MM = new vector<double>();
vector<double>* sigma_betaVec_MM = new vector<double>();
vector<double>* sigma_alphaVec_ML = new vector<double>();
vector<double>* sigma_betaVec_ML = new vector<double>();

double xMin = -1.;
double xMax = 1.;
int nbins = 100;

double func(Double_t* x, Double_t* par){ // pdf (normalizzata a 1)
  return (0.5+par[0]*x[0]);
}


void fcn(int& npar, double* deriv, double& f, double par[], int flag){
  double logLk=0;
  for (int ii=0; ii<xVecPtr->size(); ii++){ // xVec è un variabile globale definita al di fuori di tutte le funzioni e viene riempita dentro la funzione principale fitta() con le variabili x_i
    double x = xVecPtr->at(ii);
    logLk = logLk+log(func(&x,par));
  }
  f = -logLk;
}

void getData(Int_t N,Double_t alpha,int seed=0){
  int nsim = 0;
  double x, y;
  double par[1] = {alpha};
  TRandom2* r = new TRandom2(seed);
  while( nsim < N ){
    x = r->Uniform( xMin, xMax );
    y = r->Uniform( 0, 1 );
    if( y > func( &x, par ) ) continue;
    xVecPtr->push_back(x);
    nsim++;
  }

}
  
void fit_binned(int nev){

  TH1F *h = new TH1F("h","",nbins,xMin,xMax);
  h->Sumw2();

  TF1* f_fit_LSM = new TF1("f_fit_LSM","0.5+[0]*x",xMin,xMax);
  TF1* f_fit_LS = new TF1("f_fit_LS","0.5+[0]*x",xMin,xMax);
  TF1* f_fit_ML = new TF1("f_fit_ML","0.5+[0]*x",xMin,xMax);

  for (int ii=0; ii<xVecPtr->size(); ii++){
       h->Fill(xVecPtr->at(ii));
  }

  double integrale=0;
  for (int ii=1; ii<h->GetNbinsX(); ii++){
    integrale=integrale+h->GetBinContent(ii)*h->GetBinWidth(ii);
  }

  h->Scale(1./(double)integrale);  

  h->Draw();

  h->Fit(f_fit_LSM,"QS"); //fit di chi quadro modificato
  h->Fit(f_fit_LS,"QSP"); //fit di chi quadro
  h->Fit(f_fit_ML,"QSL"); //fit di Maximum likelihood (binned)

  /* FIT UNBINNED */

  int npar=1;

  double stepSize[1]={0.001};
  const char* parName[npar];
  double p[1];
  parName[0]={"alpha"};

  cout<<"FIT UNBINNED"<<endl;

  TMinuit *gMinuit = new TMinuit(npar);  //initialize TMinuit
  gMinuit->SetFCN(fcn); // gli diciamo che funzione minimizzare
  gMinuit->SetErrorDef(0.5); // altrimenti di default è 1, quello del chi2
  gMinuit->SetPrintLevel(1);
  gMinuit->SetMaxIterations(100);

  Int_t ierflg = 0;

  p[0]=0.5;
  
  for (int jj=0; jj<npar; jj++)
    {
      gMinuit->mnparm(jj,parName[jj],p[jj],stepSize[jj],0.,1.,ierflg); // setto le informazioni sui parametri, se minVale maxVal sono 0 e 0, i parametri non hanno limiti   
    }


  gMinuit->Migrad();

  double par_lk[1];
  double par_lk_err[1];

  for (int jj=0; jj<npar; jj++)
    {
	 gMinuit->GetParameter(jj,par_lk[jj],par_lk_err[jj]);
    }

  alphaVec_ML->push_back(f_fit_ML->GetParameter(0));
  sigma_alphaVec_ML->push_back(f_fit_ML->GetParError(0));
  
  cout<<"******* RESULTS *******"<<endl;
  cout<<"LSM"<<endl;
  cout<<"  alpha    "<<f_fit_LSM->GetParameter(0)<<"±"<<f_fit_LSM->GetParError(0)<<" err_fit_rel: "<< f_fit_LSM->GetParError(0)/f_fit_LSM->GetParameter(0)<<endl;
   cout<<"******************"<<endl;
  cout<<"LS"<<endl;
  cout<<"  alpha    "<<f_fit_LS->GetParameter(0)<<"±"<<f_fit_LS->GetParError(0)<<" err_fit_rel: "<< f_fit_LS->GetParError(0)/f_fit_LS->GetParameter(0)<<endl;
   cout<<"******************"<<endl;
  cout<<"ML"<<endl;
  cout<<"  alpha    "<<f_fit_ML->GetParameter(0)<<"±"<<f_fit_ML->GetParError(0)<<" err_fit_rel: "<< f_fit_ML->GetParError(0)/f_fit_ML->GetParameter(0)<<endl;
   cout<<"******************"<<endl;
  cout<<"MLU"<<endl;
  cout<<"  alpha    "<<par_lk[0]<<"±"<<par_lk_err[0]<<" err_fit_rel: "<< par_lk_err[0]/par_lk[0]<<endl;
   cout<<"******************"<<endl;

}

void metod_MM(int nev){
  //method of moments (vedi slide lezione 17)

  double sum_a1=0;

    for (int ii=0; ii<xVecPtr->size(); ii++)
      {
	sum_a1 = sum_a1 + xVecPtr->at(ii);
      }
  
  double e1 = (1./(double)xVecPtr->size())*sum_a1; // media aritmetica di a1
 
  double alpha = (3./2.)*e1;
  double sigma_alpha = sqrt((1./(double)xVecPtr->size())*(3./4.-pow(alpha,2)));

  cout<<"MM"<<endl;
  cout<<"  alpha    "<<alpha<<"±"<<sigma_alpha<<" err_fit_rel: "<< sigma_alpha/alpha<<endl;
   
  TF1* f_fit_MM = new TF1("f_fit_MM",func,xMin,xMax,1);

  f_fit_MM->SetParameter(0,alpha); // per fare il confronto corretto, devo prendere alpha e beta calcolati con gli stessi dati su cui ho fatto i fit

  double sum =0;
  TH1F *h_veri = new TH1F("h_veri","",nbins,xMin,xMax);
  h_veri->Sumw2();

  for (int xx=0; xx<h_veri->GetNbinsX(); xx++)
    sum= sum+f_fit_MM->Eval(h_veri->GetBinCenter(xx));

  f_fit_MM->SetParameter(2,(double)nev/(double)sum);
  f_fit_MM->SetLineColor(kBlue);
  f_fit_MM->Draw("same");
    
}

void DoAll(int nev, double parvalue){

  getData(nev,parvalue,2983);
  fit_binned(nev);
  metod_MM(nev);
  xVecPtr->clear();
    
}
