
#include <iostream>
#include <vector>
#include <fstream>
#include <TApplication.h>
#include <TH1D.h>
#include <TCanvas.h>
using namespace std;


int main(int argc, char** argv) {
  cout<<"Lezione 2 - VARIABILI ALEATORIE"<<endl;
  
  fstream f;
  f.open("gaus_data.txt", ios::in);
  if(f.fail()==true)
    cerr<<"File can't be opened!"<<endl;
  else{
    vector<double>vx;
    while(f.eof()==false){
      double x;
      f>>x>>ws;
      vx.push_back(x);
    } 
    f.close();
    cout<<"Operation successfully completed"<<endl;

    
    int dim=vx.size();
    cout<<"The vector has size "<<dim<<endl;

    
    cout<<"Random variables - pdf"<<endl;
    TApplication app("app", &argc, argv);
    TCanvas* c1 = new TCanvas("c1", "Canvas", 800, 600);
    TH1D* h1 = new TH1D("h1", "Random variables - pdf", 100, 0, 9);
    h1->Sumw2();
    for (double x : vx)
        h1->Fill(x);
    //pdf hypothesis: gaussian distribution
    h1->Fit("gaus");
    h1->GetXaxis()->SetTitle("x");
    h1->GetYaxis()->SetTitle("Counts");
    h1->Draw();
    app.Run();
  }
  
}







