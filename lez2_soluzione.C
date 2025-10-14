// -- Soluzione della prof --

#include <iostream>
#include <vector>
#include <fstream>
#include <TStyle.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
using namespace std;


void get_pdf(const char* file_name) {

  double x;
  
  TTree* T = new TTree();
  T->ReadFile(file_name,"var/D");
  T->SetBranchAddress("var",&x);
  
  int tot_ev = T->GetEntries();

  TH1F* h = new TH1F("h","Data Distribution;x;Entries",100,2,8); // 100 bins con x_min=2 e x_max=8. Ogni bin è largo ∆x=(8-2)/100=0.06
  h->Sumw2();
  
  for (int ii=0; ii<tot_ev; ii++)
    {
      T->GetEntry(ii);
      h->Fill(x);
    }

  // 2 modi per trovare la PDF //

  /******* I MODO - FIT GUASSIANO E POI NORMALIZZO*/

  TF1* pdf_I = new TF1("pdf_I","[0]*exp(-0.5*pow((x-[1])/[2],2))",2,8);
  pdf_I->SetParLimits(0,0,5000);
  pdf_I->SetParLimits(1,2,8);
  pdf_I->SetParLimits(2,0,99);
  h->Fit(pdf_I);

  // Ho la gaussiana...ma per essere una pdf deve essere normalizzata a 1

  double integrale = pdf_I->Integral(2,8);
  double norm = (double)pdf_I->GetParameter(0)/(double)integrale;
  pdf_I->SetParameter(0,norm);

  /******* II MODO - NORMALIZZO E POI FIT GAUSSIANO */
  
  TH1F* h_norm = new TH1F("h_norm","",100,2,8);
  h_norm->Sumw2();

  for (int ii=0; ii<h_norm->GetNbinsX(); ii++){
    double nev_norm = (double)h->GetBinContent(ii)/((double)tot_ev*h->GetBinWidth(ii));
    h_norm->SetBinContent(ii,nev_norm);
    h_norm->SetBinError(ii,(double)h->GetBinError(ii)/((double)tot_ev*h->GetBinWidth(ii)));
  }

  TF1* pdf_II = new TF1("pdf_II","[0]*exp(-0.5*pow((x-[1])/[2],2))",2,8);
  pdf_II->SetParLimits(0,0,1);
  pdf_II->SetParLimits(1,2,8);
  pdf_II->SetParLimits(2,0,99);

  h_norm->Fit("pdf_II");
  
  TCanvas* c = new TCanvas();
  c->Divide(3,1);
  c->cd(1);
  gStyle->SetOptFit(10); // visualizza nella box della statistica i risultait del Fit

  h->Draw();
  c->cd(2);
  h_norm->Draw();

  c->cd(3);
  pdf_I->Draw();
  pdf_II->Draw("same");

  cout<<"PDF I: Mean="<<pdf_I->GetParameter(1)<<"±"<<pdf_I->GetParError(1)<<", Sigma="<<pdf_I->GetParameter(2)<<"±"<<pdf_I->GetParError(2)<<", Const="<<pdf_I->GetParameter(0)<<"±"<<pdf_I->GetParError(0)/integrale<<endl;
  cout<<"PDF II: Mean="<<pdf_II->GetParameter(1)<<"±"<<pdf_II->GetParError(1)<<", Sigma="<<pdf_II->GetParameter(2)<<"±"<<pdf_II->GetParError(2)<<", Const="<<pdf_II->GetParameter(0)<<"±"<<pdf_II->GetParError(0)<<endl;
}
