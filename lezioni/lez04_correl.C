// ---- Esercizio: De-correlazione dei dati ----

#include <iostream>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;

void decorrelate(const char* file_name){

  //---- READ DATA AND SAVE IN TREE
  double x,y;
  TTree* T = new TTree("T", "correlated_events");
  T->ReadFile(file_name,"x/D:y/D");
  T->SetBranchAddress("x",&x);
  T->SetBranchAddress("y",&y);

  int tot_ev = T->GetEntries(); //conta numero di righe di dati (coppie)
  cout<<"Total couples of data: "<<tot_ev<<endl;


  //---- HISTOGRAMS OF DATA ----
  double xmin=3, xmax=5, ymin=5.3, ymax=9.3;
  int nbinsx=100, nbinsy=100;
  double wbinx, wbiny;
  double xmean, ymean, xstd, ystd;
  
  TH1D* h1x = new TH1D("h1x", "Histogram x", nbinsx, xmin, xmax);
  h1x->Sumw2();
  h1x->GetXaxis()->SetTitle("x");
  h1x->GetYaxis()->SetTitle("counts");
  wbinx = h1x->GetBinWidth(1);
  
  TH1D* h1y = new TH1D("h1y", "Histogram y", nbinsy, ymin, ymax);
  h1y->Sumw2();
  h1y->GetXaxis()->SetTitle("y");
  h1y->GetYaxis()->SetTitle("counts");
  wbiny = h1y->GetBinWidth(1);
  
  TH2D* h2 = new TH2D("h2", "Histogram yVSx", nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  h2->GetXaxis()->SetTitle("x");
  h2->GetYaxis()->SetTitle("y");

  for (int i=0; i<tot_ev; i++){
      T->GetEntry(i);
      h1x->Fill(x);
      h1y->Fill(y);
      h2->Fill(x,y);
  }

  TCanvas* c1 = new TCanvas();
  c1->Divide(3,1);
  c1->cd(1); h1x->Draw();
  c1->cd(2); h1y->Draw();
  c1->cd(3); h2->Draw("COLZ"); //"COLZ" mostra la color bar



  //---- NOTABLE VALUES + CORRELATION
  xmean = h1x->GetMean();
  xstd = h1x->GetStdDev();
  cout<<"--- Hist 1D x ---"<<endl<<
    "Bins: "<<nbinsx<<", Min: "<<xmin<<", Max: "<<xmax<<", Bins width: "<<wbinx<<endl<<
    "Mean: "<<xmean<<", StD: "<<xstd<<endl;

  ymean = h1y->GetMean();
  ystd = h1y->GetStdDev();
  cout<<"--- Hist 1D y ---"<<endl<<
    "Bins: "<<nbinsy<<", Min: "<<ymin<<", Max: "<<ymax<<", Bins width: "<<wbiny<<endl<<
    "Mean: "<<ymean<<", StD: "<<ystd<<endl;

  double corr = h2->GetCorrelationFactor(); //corefficiente di correlazione ρ
  double cov = h2->GetCovariance(); //covarianza V_xy
  cout<<"--- Hist 2D yVSx ---"<<endl<<
    "Covariance (V_xy): "<<cov<<", Correlation factor (ρ): "<<corr<<endl;


  
  //---- DECORRELATION
  double theta = 0.5*atan((2*cov)/(xstd*xstd-ystd*ystd));
  cout<<"Theta correlation angle: "<<theta<<"rad, "<<theta*180./TMath::Pi()<<"°"<<endl;

  double z,w;
  
  TTree* Td = new TTree("Td", "decorrelated_events");
  Td->Branch("z", &z, "z/D");
  Td->Branch("w", &w, "w/D");

  double zmin=0.5, zmax=2.2, wmin=6, wmax=10.5;
  int nbinsz=100, nbinsw=100;
  double wbinz, wbinw;
  double zmean, zstd, wmean, wstd;
  
  TH1D* h1z = new TH1D("h1z", "Histogram z", nbinsz, zmin, zmax);
  h1z->Sumw2();
  h1z->GetXaxis()->SetTitle("z");
  h1z->GetYaxis()->SetTitle("counts");
  wbinz = h1z->GetBinWidth(1);
  
  TH1D* h1w = new TH1D("h1w", "Histogram w", nbinsw, wmin, wmax);
  h1w->Sumw2();
  h1w->GetXaxis()->SetTitle("w");
  h1w->GetYaxis()->SetTitle("counts");
  wbinw = h1w->GetBinWidth(1);
  
  TH2D* h2zw = new TH2D("h2zw", "Histogram wVSz", nbinsz, zmin, zmax, nbinsw, wmin, wmax);
  h2zw->GetXaxis()->SetTitle("z");
  h2zw->GetYaxis()->SetTitle("w");
  
  for (int i=0; i<tot_ev; i++){
      T->GetEntry(i);
      z = x*cos(theta)+y*sin(theta);
      w = -x*sin(theta)+y*cos(theta);
      h1z->Fill(z);
      h1w->Fill(w);
      h2zw->Fill(z,w);
      Td->Fill();
  }

  TCanvas* c2 = new TCanvas();
  c2->Divide(3,1);
  c2->cd(1); h1z->Draw();
  c2->cd(2); h1w->Draw();
  c2->cd(3); h2zw->Draw("COLZ");

  TCanvas* c3 = new TCanvas();
  c3->Divide(2,1);
  c3->cd(1); h2->Draw("COLZ");
  c3->cd(2); h2zw->Draw("COLZ");


  //---- NOTABLE VALUES -  DECORRELATED
  zmean = h1z->GetMean();
  zstd = h1z->GetStdDev();
  cout<<"--- Hist 1D z ---"<<endl<<
    "Bins: "<<nbinsz<<", Min: "<<zmin<<", Max: "<<zmax<<", Bins width: "<<wbinz<<endl<<
    "Mean: "<<zmean<<", StD: "<<zstd<<endl;

  wmean = h1w->GetMean();
  wstd = h1w->GetStdDev();
  cout<<"--- Hist 1D w ---"<<endl<<
    "Bins: "<<nbinsw<<", Min: "<<wmin<<", Max: "<<wmax<<", Bins width: "<<wbinw<<endl<<
    "Mean: "<<wmean<<", StD: "<<wstd<<endl;

  double corr_zw = h2zw->GetCorrelationFactor();
  double cov_zw = h2zw->GetCovariance();
  cout<<"--- Hist 2D wVSz ---"<<endl<<
    "Covariance (V_zw): "<<cov_zw<<", Correlation factor (ρ'~0): "<<corr_zw<<endl;


  
  //---- TFILE
  TFile* hist_file = new TFile("hist_file.root", "RECREATE");
  //"RECREATE" sovrascrive file precedente se c'era

  //Write() serve per salvare i plot nel file ROOT
  T->Write();
  h1x->Write();
  h1y->Write();
  //h1xy->Write();
  h2->Write();
  c1->Write();
  Td->Write();
  h1z->Write();
  h1w->Write();
  h2zw->Write();
  c2->Write();
  c3->Write();

  hist_file->Close();

  /* //---- PDF (alternativa al TFile)
  c1->GetPad(1)->SetBottomMargin(0.15);
  c1->GetPad(2)->SetBottomMargin(0.15);
  c1->GetPad(3)->SetBottomMargin(0.15);
  c2->GetPad(1)->SetBottomMargin(0.15);
  c2->GetPad(2)->SetBottomMargin(0.15);
  c2->GetPad(3)->SetBottomMargin(0.15);
  c3->GetPad(1)->SetBottomMargin(0.15);
  c3->GetPad(2)->SetBottomMargin(0.15);
  c1->Print("file_pdf.pdf[");
  c1->Print("file_pdf.pdf");
  c2->Print("file_pdf.pdf");
  c3->Print("file_pdf.pdf");
  c3->Print("file_pdf.pdf]");
  */
  
}


/* HOW TO OPEN FROM TERMINAL
Create file emacs with the script with ".C" extention: lez4_correl.C.
Open root with command "root".
Write command ".L lez4_correl.C" to load macro file.
Use the function with command "decorrelate("data_correlated.txt")".
Write "new TBrowser" to open the tree structure saved thanks to the TFile. Click on the plots of the data.
*/
