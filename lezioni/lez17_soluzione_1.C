// Metodi di stima dei parametri a confronto (soluzione 1)
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

// per la definizione della pdf usata, vedere slides lezione 13
vector<double>* xVecPtr = new vector<double>();
vector<double>* alphaVec_MM = new vector<double>();
vector<double>* betaVec_MM = new vector<double>();
vector<double>* alphaVec_ML = new vector<double>();
vector<double>* betaVec_ML = new vector<double>();

vector<double>* sigma_alphaVec_MM = new vector<double>();
vector<double>* sigma_betaVec_MM = new vector<double>();
vector<double>* sigma_alphaVec_ML = new vector<double>();
vector<double>* sigma_betaVec_ML = new vector<double>();

double xMin = -0.95;
double xMax = 0.95;
int nbins = 100;

double func(Double_t* x, Double_t* par){ // pdf (normalizzata a 1)
  return (1+par[0]*x[0]+par[1]*pow(x[0],2))/((xMax-xMin)+par[0]/2*(pow(xMax,2)-pow(xMin,2))+par[1]/3*(pow(xMax,3)-pow(xMin,3)));
}

double func_not_norm(Double_t* x, Double_t* par){
  return par[2]*(1+par[0]*x[0]+par[1]*pow(x[0],2));
}

void fcn(int& npar, double* deriv, double& f, double par[], int flag){
  double logLk=0;
  for (int ii=0; ii<xVecPtr->size(); ii++){ // xVec è un variabile globale definita al di fuori di tutte le funzioni e viene riempita dentro la funzione principale fitta() con le variabili x_i
    double x = xVecPtr->at(ii);
    logLk = logLk+log(func(&x,par));
  }
  f = -logLk;
}

void getData(Int_t N,Double_t alpha,Double_t beta, int seed=0){
  int nsim = 0;
  double x, y;
  double par[2] = {alpha, beta};
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

  TF1* f_fit_LSM = new TF1("f_fit_LSM",func_not_norm,xMin,xMax,3);
  TF1* f_fit_LS = new TF1("f_fit_LS",func_not_norm,xMin,xMax,3);
  TF1* f_fit_ML = new TF1("f_fit_ML",func_not_norm,xMin,xMax,3);

  for (int ii=0; ii<xVecPtr->size(); ii++){
       h->Fill(xVecPtr->at(ii));
  }
  
  h->Draw();

  f_fit_LSM->SetParLimits(0,0.1,0.9);
  f_fit_LSM->SetParLimits(1,0.1,0.9);
  f_fit_LSM->SetParLimits(2,0,9999999);
  f_fit_LS->SetParLimits(0,0.1,0.9);
  f_fit_LS->SetParLimits(1,0.1,0.9);
  f_fit_LS->SetParLimits(2,0,9999999);
  f_fit_ML->SetParLimits(0,0.1,0.9);
  f_fit_ML->SetParLimits(1,0.1,0.9);
  f_fit_ML->SetParLimits(2,0,9999999);

  h->Fit(f_fit_LSM,"QS"); //fit di chi quadro modificato
  h->Fit(f_fit_LS,"QSP"); //fit di chi quadro
  h->Fit(f_fit_ML,"QSL"); //fit di Maximum likelihood (binned)

  /* FIT UNBINNED */

  int npar=2;

  double stepSize[2]={0.001,0.001};
  double  minVal[2]={0.,0.};
  double  maxVal[2]={1.,1.};
  double p[2] ={0.,0.};
  const char* parName[npar];

  parName[0]={"alpha"};
  parName[1]={"beta"};
  cout<<"FIT UNBINNED"<<endl;
  TMinuit *gMinuit = new TMinuit(npar);  //initialize TMinuit
  gMinuit->SetFCN(fcn); // gli diciamo che funzione minimizzare
  gMinuit->SetErrorDef(0.5); // altrimenti di default è 1, quello del chi2
  gMinuit->SetPrintLevel(1);
  gMinuit->SetMaxIterations(100);

  Int_t ierflg = 0;

  for (int jj=0; jj<npar; jj++)
    {
      gMinuit->mnparm(jj,parName[jj],p[jj],stepSize[jj],0,0,ierflg); // setto le informazioni sui parametri, se minVale maxVal sono 0 e 0, i parametri non hanno limiti   
    }


  gMinuit->Migrad();

  double par_lk[2];
  double par_lk_err[2];

  for (int jj=0; jj<npar; jj++)
    {
	 gMinuit->GetParameter(jj,par_lk[jj],par_lk_err[jj]);
    }

  alphaVec_ML->push_back(f_fit_ML->GetParameter(0));
  betaVec_ML->push_back(f_fit_ML->GetParameter(1));
  sigma_alphaVec_ML->push_back(f_fit_ML->GetParError(0));
  sigma_betaVec_ML->push_back(f_fit_ML->GetParError(1));
  
  cout<<"******* RESULTS *******"<<endl;
  cout<<"LSM"<<endl;
  cout<<"  alpha    "<<f_fit_LSM->GetParameter(0)<<"±"<<f_fit_LSM->GetParError(0)<<" err_fit_rel: "<< f_fit_LSM->GetParError(0)/f_fit_LSM->GetParameter(0)<<endl;
  cout<<"  beta    "<<f_fit_LSM->GetParameter(1)<<"±"<<f_fit_LSM->GetParError(1)<<" err_fit_rel: "<< f_fit_LSM->GetParError(1)/f_fit_LSM->GetParameter(1)<<endl;
  cout<<"******************"<<endl;
  cout<<"LS"<<endl;
  cout<<"  alpha    "<<f_fit_LS->GetParameter(0)<<"±"<<f_fit_LS->GetParError(0)<<" err_fit_rel: "<< f_fit_LS->GetParError(0)/f_fit_LS->GetParameter(0)<<endl;
  cout<<"  beta    "<<f_fit_LS->GetParameter(1)<<"±"<<f_fit_LS->GetParError(1)<<" err_fit_rel: "<< f_fit_LS->GetParError(1)/f_fit_LS->GetParameter(1)<<endl;
  cout<<"******************"<<endl;
  cout<<"ML"<<endl;
  cout<<"  alpha    "<<f_fit_ML->GetParameter(0)<<"±"<<f_fit_ML->GetParError(0)<<" err_fit_rel: "<< f_fit_ML->GetParError(0)/f_fit_ML->GetParameter(0)<<endl;
  cout<<"  beta    "<<f_fit_ML->GetParameter(1)<<"±"<<f_fit_ML->GetParError(1)<<" err_fit_rel: "<< f_fit_ML->GetParError(1)/f_fit_ML->GetParameter(1)<<endl;
  cout<<"******************"<<endl;
  cout<<"MLU"<<endl;
  cout<<"  alpha    "<<par_lk[0]<<"±"<<par_lk_err[0]<<" err_fit_rel: "<< par_lk_err[0]/par_lk[0]<<endl;
  cout<<"  beta    "<<par_lk[1]<<"±"<<par_lk_err[1]<<" err_fit_rel: "<< par_lk_err[1]/par_lk[1]<<endl;
  cout<<"******************"<<endl;

}

void metod_MM(int nev){
  //method of moments (vedi slide lezione 17)

  double sum_a1=0;
  double sum_a2=0;

    for (int ii=0; ii<xVecPtr->size(); ii++)
      {
	sum_a1 = sum_a1 + xVecPtr->at(ii);
	sum_a2 = sum_a2 + pow(xVecPtr->at(ii),2);
      }
  
  double e1 = (1./(double)xVecPtr->size())*sum_a1; // media aritmetica di a1
  double e2 = (1./(double)xVecPtr->size())*sum_a2; // media artmetica di a2

  double d[5];

  for (int jj=1; jj<=5; jj++)
    d[jj-1]=(1./(double)jj)*(pow(xMax,jj)-pow(xMin,jj));

  double alpha = ((e1*d[2]-d[3])*(e2*d[0]-d[2])-(e1*d[0]-d[1])*(e2*d[2]-d[4]))/( (e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]) );
  double beta = ((e1*d[0]-d[1])*(e2*d[1]-d[3])-(e1*d[1]-d[2])*(e2*d[0]-d[2]))/( (e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]) );

  
  double d_alpha_d_e1 = (d[2]*(e2*d[0]-d[2])-d[0]*(e2*d[2]-d[4])) / ( (e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]) )-alpha*(d[1]*(e2*d[2]-d[4])-d[2]*(e2*d[1]-d[3]))/((e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]));

    double d_alpha_d_e2 = (d[0]*(e1*d[2]-d[3])-d[2]*(e1*d[0]-d[1])) / ( (e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]) )-alpha*(d[2]*(e1*d[1]-d[2])-d[1]*(e1*d[2]-d[3]))/((e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]));

    double d_beta_d_e1 = (d[0]*(e2*d[1]-d[3])-d[1]*(e2*d[0]-d[2])) / ( (e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]) )-beta*(d[1]*(e2*d[2]-d[4])-d[2]*(e2*d[1]-d[3]))/((e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]));

    double d_beta_d_e2 = (d[1]*(e1*d[0]-d[1])-d[0]*(e1*d[1]-d[2])) / ( (e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]) )-beta*(d[2]*(e1*d[1]-d[2])-d[1]*(e1*d[2]-d[3]))/((e1*d[1]-d[2])*(e2*d[2]-d[4])-(e1*d[2]-d[3])*(e2*d[1]-d[3]));
 
    double cov_e1_e1 = 0;
    double cov_e2_e2 = 0;
    double cov_e1_e2 = 0;

    for (int ii=0; ii<xVecPtr->size(); ii++)
      {
        cov_e1_e1 = cov_e1_e1 + (xVecPtr->at(ii)-e1)*(xVecPtr->at(ii)-e1);
        cov_e2_e2 = cov_e2_e2 + (xVecPtr->at(ii)*xVecPtr->at(ii)-e2*e2)*(xVecPtr->at(ii)*xVecPtr->at(ii)-e2*e2);
        cov_e1_e2 = cov_e1_e2 + (xVecPtr->at(ii)-e1)*(xVecPtr->at(ii)*xVecPtr->at(ii)-e2*e2);
      }
   
    double sigma_alpha=0;
    double sigma_beta=0;

    sigma_alpha = d_alpha_d_e1*d_alpha_d_e1*cov_e1_e1+2*d_alpha_d_e1*d_alpha_d_e2*cov_e1_e2+d_alpha_d_e2*d_alpha_d_e2*cov_e2_e2;
    sigma_beta = d_beta_d_e1*d_beta_d_e1*cov_e1_e1+2*d_beta_d_e1*d_beta_d_e2*cov_e1_e2+d_beta_d_e2*d_beta_d_e2*cov_e2_e2;

    sigma_alpha = sqrt(sigma_alpha/((nev-1)*nev));
    sigma_beta = sqrt(sigma_beta/((nev-1)*nev));

    alphaVec_MM->push_back(alpha);
    betaVec_MM->push_back(beta);
    sigma_alphaVec_MM->push_back(sigma_alpha);
    sigma_betaVec_MM->push_back(sigma_beta);
  
    cout<<"MM"<<endl;
    cout<<"  alpha    "<<alpha<<"±"<<sigma_alpha<<" err_fit_rel: "<< sigma_alpha/alpha<<endl;
    cout<<"  beta    "<<beta<<"±"<<sigma_beta<<" err_fit_rel: "<< sigma_beta/(beta)<<endl;
 
  TF1* f_fit_MM = new TF1("f_fit_MM",func_not_norm,xMin,xMax,3);

  f_fit_MM->SetParameter(0,alpha); // per fare il confronto corretto, devo prendere alpha e beta calcolati con gli stessi dati su cui ho fatto i fit
  f_fit_MM->SetParameter(1,beta);

  f_fit_MM->SetParameter(2,1);
  double sum =0;
  TH1F *h_veri = new TH1F("h_veri","",nbins,xMin,xMax);
  h_veri->Sumw2();

  for (int xx=0; xx<h_veri->GetNbinsX(); xx++)
    sum= sum+f_fit_MM->Eval(h_veri->GetBinCenter(xx));

  f_fit_MM->SetParameter(2,(double)nev/(double)sum);
  f_fit_MM->SetLineColor(kBlue);
  f_fit_MM->Draw("same");
    
}

void DoAll(int nev){
    getData(nev,0.5,0.5);
    fit_binned(nev);
    metod_MM(nev);
    xVecPtr->clear();
    
    for (int ii=0; ii<500; ii++){
      getData(nev,0.5,0.5,ii);
      fit_binned(nev);
      metod_MM(nev);
      xVecPtr->clear();
      }

  TH1F* h_alpha_MM = new TH1F("h_alpha_MM","MM;#hat{#alpha};Entries",100,0.,1.);
  h_alpha_MM->SetName("MM");
  TH1F* h_alpha_ML = new TH1F("h_alpha_ML",";#hat{#alpha};Entries",100,0.,1.);
  h_alpha_ML->SetName("ML");

  TH1F* h_beta_MM = new TH1F("h_beta_MM",";#hat{#beta};Entries",100,0.,1.);
  h_beta_MM->SetName("MM");
  TH1F* h_beta_ML = new TH1F("h_beta_ML",";#hat{#beta};Entries",100,0.,1.);
  h_beta_ML->SetName("ML");
  
  TH1F* h_sigma_alpha_MM = new TH1F("h_sigma_alpha_MM","MM;#hat{#sigma}_{#alpha};Entries",100,0.,0.1);
  h_sigma_alpha_MM->SetName("MM");
  TH1F* h_sigma_alpha_ML = new TH1F("h_sigma_alpha_ML",";#hat{#sigma}_{#alpha};Entries",100,0.,0.1);
  h_sigma_alpha_ML->SetName("ML");

  TH1F* h_sigma_beta_MM = new TH1F("h_sigma_beta_MM",";#hat{#sigma}_{#beta};Entries",100,0.,0.3);
  h_sigma_beta_MM->SetName("MM");
  TH1F* h_sigma_beta_ML = new TH1F("h_sigma_beta_ML",";#hat{#sigma}_{#beta};Entries",100,0.,0.3);
  h_sigma_beta_ML->SetName("ML");
  
  for (int jj=0; jj<alphaVec_MM->size(); jj++)
    {
      h_alpha_MM->Fill(alphaVec_MM->at(jj));
      h_alpha_ML->Fill(alphaVec_ML->at(jj));
      h_beta_MM->Fill(betaVec_MM->at(jj));
      h_beta_ML->Fill(betaVec_ML->at(jj));
      h_sigma_alpha_MM->Fill(sigma_alphaVec_MM->at(jj));
      h_sigma_alpha_ML->Fill(sigma_alphaVec_ML->at(jj));
      h_sigma_beta_MM->Fill(sigma_betaVec_MM->at(jj));
      h_sigma_beta_ML->Fill(sigma_betaVec_ML->at(jj));
    }

  TCanvas *c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  h_alpha_MM->Draw();
  h_alpha_ML->Draw("sames");
  c->cd(2);
  h_sigma_alpha_MM->Draw();
  h_sigma_alpha_ML->Draw("sames");
  c->cd(3);
  h_beta_MM->Draw();
  h_beta_ML->Draw("sames");
  c->cd(4);
  h_sigma_beta_MM->Draw();
  h_sigma_beta_ML->Draw("sames");
    
}
