// --- MC: Hit or Miss ---

#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <TGraph.h>
#include <TCanvas.h>
using namespace std;
const int N=10000; //#dati generati

void hit_miss(){
  cout<<endl<<"---------- MC: Hit or Miss method ----------"<<endl;

  // definisco estremi dominio
  cout<<endl<<"Definisci estremi del dominio: ";
  double x_min, x_max, x_i, y_i;
  cin>>x_min>>x_max;
  cout<<" ----> intervallo x:("<<x_min<<","<<x_max<<")"<<endl;
  
  TRandom3 *r = new TRandom3();

  // genero il dominio
  vector<double> x;
  for(int i=0; i<N; i++){
    x_i = r->Uniform(x_min, x_max);
    x.push_back(x_i);
  }

  // valuto la funzione
  vector<double> f;
  double y_max=0;
  for(int i=0; i<N; i++){
    y_i=3*x[i]+5;    // espressione analitica della pdf generanda
    f.push_back(y_i);
    if(y_i>y_max)
      y_max=y_i;
  }
  cout<<" ----> intervallo y:(0,"<<y_max<<")"<<endl;

  // genero il codominio
  vector<double> y;
  for(int i=0; i<N; i++){
    y_i = r->Uniform(0, y_max);
    y.push_back(y_i);
  }

  // seleziono solo i dati che fanno la pdf voluta
  vector<double> x_pdf, y_pdf;
  int n=0;
  for(int i=0; i<N; i++){
    if(f[i]>y[i]){
      x_pdf.push_back(x[i]);
      y_pdf.push_back(y[i]);
      n++;
    }
  }
  cout<<endl<<"Dati generati: "<<N<<"; Dati selezionati: "<<n<<endl<<" ----> Dati persi: "<<N-n<<endl<<endl;


  // plots per rappresentazione grafica
  TGraph *p_rand = new TGraph(N, &x[0], &y[0]);
  p_rand->SetTitle("Random generation;x;y");
  p_rand->SetMarkerSize(0.2);
  p_rand->SetMarkerColor(kGray);
  p_rand->SetMarkerStyle(20);

  TGraph *p_pdf = new TGraph(n, &x_pdf[0], &y_pdf[0]);
  p_pdf->SetTitle("PDF;x;y");
  p_pdf->SetMarkerSize(0.2);
  p_pdf->SetMarkerStyle(20);
  
  TGraph *p_fun = new TGraph(N, &x[0], &f[0]);
  p_fun->SetTitle("Function;x;f(x)");
  p_fun->SetLineColor(kPink+8);

  TCanvas *c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  p_fun->Draw("AL");
  c->cd(2);
  p_rand->Draw("AP");
  c->cd(3);
  p_pdf->Draw("AP");
  c->cd(4);
  p_rand->SetTitle("Results;x;counts"); 
  p_rand->Draw("AP");
  p_pdf->Draw("P SAME");
  p_fun->Draw("L SAME");
  
}
