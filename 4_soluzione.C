// SOLUZIONE proposta dalla prof

#include "TH1F.h"
#include "TH2F.h"
#include "TVector.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
using namespace std;

void readData(){

  Float_t x;
  Float_t y;
  Float_t mean_x =0.;
  Float_t mean_y =0.;
  Float_t mean_xy=0.;
  Float_t rms_x =0.;
  Float_t rms_y =0.;
  Float_t rms_xy=0.;
  
  TTree* T = new TTree();
  T->ReadFile("data_correlated.txt","b_x/F:b_y/F");
  T->SetBranchAddress("b_x",&x);
  T->SetBranchAddress("b_y",&y);

  TH2F* h2 = new TH2F("h2",";x;y",100,3,5,100,5,10);
  TH2F* h2_rot = new TH2F("h2_rot",";x;y",100,6,10,100,-4,1);
  TH1F* hxy = new TH1F("hxy",";xy",200,3,60);
  
  for (int ii=0; ii<T->GetEntries(); ii++){
    T->GetEntry(ii);
    h2->Fill(x,y);
    hxy->Fill(x*y);
    mean_x = mean_x+x;
    mean_y = mean_y+y;
    mean_xy = mean_xy+x*y;
  }
  mean_x = (Float_t)mean_x/(Float_t)T->GetEntries();
  mean_y = (Float_t)mean_y/(Float_t)T->GetEntries();
  mean_xy = (Float_t)mean_xy/(Float_t)T->GetEntries();
  
  for (int ii=0; ii<T->GetEntries(); ii++){
    T->GetEntry(ii);
    rms_x = rms_x+pow((x-mean_x),2);
    rms_y = rms_y+pow((y-mean_y),2);
    rms_xy = rms_xy+pow((x*y-mean_xy),2);

  }

    TCanvas *c = new TCanvas();
    h2->Draw("colz"); // colz è per mostrare le entries in ogni bin con diversi colori

    TCanvas *c2 = new TCanvas();
    h2->Draw("lego"); 
    
    TCanvas *c3 = new TCanvas();

    Float_t sigma_x = sqrt((Float_t)rms_x/(Float_t)T->GetEntries());
    Float_t sigma_y = sqrt((Float_t)rms_y/(Float_t)T->GetEntries());
    Float_t sigma_xy = sqrt((Float_t)rms_xy/(Float_t)T->GetEntries());
    
    cout<<"MEAN x: "<<(Float_t)mean_x<<"       from HISTO: "<<h2->GetMean(1)<<endl;
    cout<<"RMS x: "<<sigma_x<<"       from HISTO: "<<h2->GetRMS(1)<<endl;

    cout<<"MEAN y: "<<(Float_t)mean_y<<"       from HISTO: "<<h2->GetMean(2)<<endl;
    cout<<"RMS y: "<<sigma_y<<"       from HISTO: "<<h2->GetRMS(2)<<endl;

    cout<<"MEAN xy: "<<(Float_t)mean_xy<<"       from HISTO: "<<hxy->GetMean()<<endl;
    cout<<"RMS xy: "<<sigma_xy<<"       from HISTO: "<<hxy->GetRMS()<<endl;

    //rho = (E[xy]-E[x]E[y]) / σx*σy
    
    Float_t corr_factor = (mean_xy-mean_x*mean_y)/(sqrt((Float_t)rms_x/(Float_t)T->GetEntries())*sqrt((Float_t)rms_y/(Float_t)T->GetEntries()));
    cout<<" correlation factor: "<<corr_factor<<"     from HISTOS: "<<h2->GetCorrelationFactor()<<endl;


    //λ+ = 1/2(σx^2+σy^2+√( (σx^2+σy^2)^2 -4 σx^2 σy^2(1-rho^2)  )

    Float_t lambda_plus = 0.5* (pow(sigma_x,2)+pow(sigma_y,2)+sqrt(pow((pow(sigma_x,2)+pow(sigma_y,2)),2)-4*pow(sigma_x,2)*pow(sigma_y,2)*(1-pow(corr_factor,2) )) );
	
   
    // theta = 0.5*atctan-(σx^2-λ+)/Vxy = -(σx^2-λ+) / rho*σx*σy   perchè  Vxy = rho*σx*σy

    Float_t tan_theta = -(pow(sigma_x,2)-lambda_plus)/(corr_factor*sigma_x*sigma_y);
    Float_t theta = atan(tan_theta);
    cout<<"theta (deg): "<<theta*180./TMath::Pi()<<endl;
    
    for (int ii=0; ii<T->GetEntries(); ii++){
      T->GetEntry(ii);
      Float_t x_rot = x*cos(theta)+y*sin(theta);
      Float_t y_rot = -x*sin(theta)+y*cos(theta);
      h2_rot->Fill(x_rot,y_rot);
    }

    h2_rot->Draw("colz");
    
}

