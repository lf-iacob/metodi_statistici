
// --- Stima di \pi con Buffon Monte Carlo: studio dei parametri ---

#include <iostream>
#include <TRandom3.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <vector>
using namespace std;

vector<double> experiment( 
			  double N, //numero di lanci per ciascun esperimento
			  double l, //distanza barre verticali
			  double d  //lunghezza di ciascun ago
			  ){
  TRandom3 *r = new TRandom3();
  double n_int=0; //condizione di interazione verificata
  double theta_i;
  double x_i;
  for(int i=0; i<N; i++){
    theta_i = r->Uniform(0,M_PI/2);
    x_i = r->Uniform(0,l/2);
    if(x_i<(d*cos(theta_i))/2){
      n_int+=1;
    }
  }

  double pi = (2*d*N)/(l*n_int); //pi stimato con MC
  
  vector<double> v = {N, l, d, n_int, pi};
  return v;
}


void stat_needles_MC(){
  cout<<endl<<"----- STIMA \u03C0 con MONTE CARLO -----"<<endl;
  double N=10000; // numero di aghi totali lanciati: N
  double l_max=100; // massima distanza tra sbarre verticali simulata
  cout<<endl<<"Lanci per ogni esperimento: "<<N<<endl<<"Distanza massima righe griglia: "<<l_max<<endl<<endl;

  vector<double> l_i;
  vector<double> d_i;
  vector<double> n_int_i;
  vector<double> pi_i;

  for(int i=1; i<l_max; i++){        //i->l
    for(int j=1; j<i; j++){          //j->d
      vector<double> v = experiment(N, i, j);
      l_i.push_back(v[1]);
      d_i.push_back(v[2]);
      n_int_i.push_back(v[3]);
      pi_i.push_back(v[4]);
    }
  }

  int dots = pi_i.size();
  //cout<<"Size: "<<dots<<endl;
  vector<double> frac(dots);
  for(int i=0; i<dots; i++)
    frac[i]=l_i[i]/d_i[i];
  
  TGraph *g1 = new TGraph(dots, &frac[0], &pi_i[0]); //plot piVSl/d
  g1->SetTitle("Plot #piVS(l/d);l/d;#pi");
  g1->SetMarkerSize(2);
  TGraph *g2 = new TGraph(dots, &frac[0], &n_int_i[0]); //plot n_intVSl/d
  g2->SetTitle("Plot n_intVS(l/d);l/d;n_int");
  g2->SetMarkerSize(2);
  TH1F *h_pi = new TH1F("h_pi", "Hist #pi;#pi;counts", 50, 2.8, 3.5);
  for(int i=0; i<dots; i++)
    h_pi->Fill(pi_i[i]);
  
  TCanvas *c = new TCanvas();
  c->Divide(3,1);
  c->cd(1);
  g1->Draw("AP"); //A: mette assi; P: mette punti
  c->cd(2);
  g2->Draw("AP");
  c->cd(3);
  h_pi->Draw();
}
