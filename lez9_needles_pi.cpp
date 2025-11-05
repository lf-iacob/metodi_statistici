// --- Stima matematica di pi: esperimento di Buffon ---

#include <iostream>
#include <TRandom3.h>
#include <TH1F.h>
#include <TCanvas.h>
using namespace std;

void needles_MC(){
  cout<<endl<<"----- STIMA \u03C0 con MONTE CARLO -----"<<endl;
  
  int N=100000; // numero di aghi totali lanciati: N
  double l=7; // distanza barre verticali: l
  double d=2.5; // lunghezza degli aghi: d(<l)
  int n_int=0; // condizione di “interazione verificata”: n_int

  cout<<endl<<"Lanci totali: "<<N<<endl<<"Distanza barre: "<<l<<endl<<"Lunghezza ago: "<<d<<endl;

  /* ESPERIMENTO: lancio degli aghi
     Genero randomicamente (pdf uniforme) l'angolo theta di incli-
     nazione di ago rispetto all'orizzontale e la distanza x tra
     CM di ago e linea della griglia più vicina.
  */

  cout<<endl<<"... Esperimento in corso ..."<<endl;
  /*
  // Istogrammi per theta e x
  TH1F *h_theta = new TH1F("h_theta", "Hist #theta; #theta; counts", N, 0, M_PI/2);
  TH1F *h_x = new TH1F("h_x", "Hist x; x; counts", N, 0, l/2);
  */
  
  TRandom *r = new TRandom();
  double theta_i;
  double x_i;
  for(int i=0; i<N; i++){
    theta_i = r->Uniform(0,M_PI/2);
    //h_theta->Fill(theta_i); //fill di istogramma in theta
    x_i = r->Uniform(0, l/2);
    //h_x->Fill(x_i); //fill di istogramma in x
    if(x_i<(d*cos(theta_i))/2){
      n_int+=1;
    }
  }

  /*
  // Plot di istogrammi
  TCanvas *c = new TCanvas();
  c->Divide(2,1);
  c->cd(1);
  h_theta->Draw();
  c->cd(2);
  h_x->Draw();
  */

  cout<<endl<<"Numero di interazioni avvenute: "<<n_int<<" su "<<N<<" lanci totali."<<endl;

  // stima di pi con metodo Monte Carlo - hit or miss
  double pi = (2*d*N)/(l*n_int);

  cout<<"\u03C0 ---> Valore stimato: "<<pi<<endl<<"  ---> Valor 'vero': "<<M_PI<<endl<<endl;

}
