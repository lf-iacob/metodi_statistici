// --- Weighted events: Monte Carlo method ---

#include <iostream>
#include <TRandom3.h>
#include <TH1F.h>
#include <vector>
#include <TF1.h>
#include <TCanvas.h>
using namespace std;

void w_mc(){
  cout<<endl<<"--- Weighted events: Monte Carlo method ---"<<endl;

  //parameters of distribution
  int n_bin, N;
  cout<<endl<<"Numero di eventi totali: N=";
  cin>>N;
  cout<<"Numero di bin: n_bin=";
  cin>>n_bin;
  double x_max, x_min;
  cout<<"Estremi del dominio: x_min,x_max=";
  cin>>x_min>>x_max;
  double delta=(x_max-x_min)/n_bin;
  cout<<"Intervallo ogni bin: ∆="<<delta<<endl;

  //uniform distribution - normalised
  double ev_i;
  TRandom3 *r = new TRandom3();
  TH1F *h_u = new TH1F("h_u", "Uniform distribution;x;counts", n_bin, x_min, x_max);
  for(int i=0; i<N; i++){
    ev_i=r->Uniform(x_min, x_max);
    h_u->Fill(ev_i);
  }
  for(int i=1; i<=n_bin; i++) //normalizzazione
    h_u->SetBinContent(i,h_u->GetBinContent(i)/(delta*N));

  //wanted pdf: gaussiana standard
  double norm=1/sqrt(2*M_PI);
  TF1 *f = new TF1("f", "([0]*exp(-0.5*((x-[1])/[2])^2))", x_min, x_max);
  f->SetParameters(norm, 0, 1); //norm,mean,std dev
  f->SetTitle("Function - pdf;x;f(x)");
    
  //MC weighted events method
  vector<double>
    x_i(n_bin),  //centri dei bin
    w_i(n_bin),  //pesi
    n_i(n_bin);  //counts del nuovo istogramma con pdf voluta
  double N_pdf=0;   //total events in new hist
  for(int i=1; i<=n_bin; i++){
    x_i[i-1]=x_min+(i-0.5)*(delta);
    w_i[i-1]=f->Eval(x_i[i-1])*(x_max-x_min);
    n_i[i-1]=h_u->GetBinContent(i)*w_i[i-1];
    N_pdf+=n_i[i-1];
  }

  //pdf distribution - normalised
  TH1F *h_pdf = new TH1F("h_pdf", "Wanted distribution;x;counts", n_bin, x_min, x_max);
  for(int i=1; i<=n_bin; i++){
    //h_pdf->SetBinContent(i,n_i[i]); //not normalised
    h_pdf->SetBinContent(i,n_i[i-1]/(delta*N_pdf)); //normalised
  }
  
  //graphic representations
  TCanvas *c = new TCanvas();
  c->Divide(3,1);
  c->cd(1);
  h_u->Draw();
  c->cd(2);
  f->Draw();
  c->cd(3);
  h_pdf->Draw();
  f->Draw("same");

}
