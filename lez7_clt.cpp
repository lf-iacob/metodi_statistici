
// CLT: Teorema del Limite Centrale

#include <iostream>
#include <TStyle.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFile.h>


void tela(
	 int n_exp, //quante volte ripeto l'esperimento per studiare distribuzione
	 int n_mis //numero di misure effettuate durante l'esperimento (piani del rivelatori)
	 ){
  //gStyle->SetOptStat(1100); //style del box con la statistica (1100: mostra Integral + Nbins + Name + Title)

  TFile *f = new TFile(Form("out_exp%d_mis%d",n_exp,n_mis), "RECREATE");
  
  TRandom * a = new TRandom();

  TH1F* h = new TH1F("h", Form("Hist - Uniform pdf (n_exp = %d, n_mis = %d);x_{i};#", n_exp, n_mis), 500,0,3); //istogramma per i dati (uniforme -> uniforme)
  h->GetXaxis()->SetTitleFont(62);
  
  
  TH1F* h_s = new TH1F("h_s","Hist - Sum (CLT);#Sigma_{i} x_{i};#",500,0,2); //istogramma per la somma dei dati (uniforme -> gaussiana)
  h_s->GetXaxis()->SetTitleFont(62);

  double x_i; //singola misura
  for (int k=0; k<n_exp; k++){
    double sum = 0.;
    for (int N=0; N<n_mis; N++){
      x_i = a->Uniform(0,2); //misura i-esima della variabile che segue una distrubuzione uniforme, mu_i = 1 e una sigma_i = 2/√12
      sum = sum + x_i; //sum per CLT sarà guassiana, con E[] = somma di mu_i e V[] = somma di sigma_i^2
      h->Fill(x_i);
      }
    h_s->Fill(sum/n_mis);
  }

  TCanvas* c = new TCanvas();
  c->Divide(2,1);
  c->cd(1);
  h->SetLineColor(kRed+1);
  h->SetLineWidth(2);
  h->SetFillColor(kRed+1);
  h->SetFillStyle(3004);
  h->Draw();
  c->cd(2);
  h_s->SetLineColor(kRed+1);
  h_s->SetLineWidth(2);
  h_s->SetFillColor(kRed+1);
  h_s->SetFillStyle(3004);
  h_s->Draw();
}


void clt(){
  tela(10000, 1);
  tela(10000, 2);
  tela(10000, 5);
  tela(10000, 100);
}
