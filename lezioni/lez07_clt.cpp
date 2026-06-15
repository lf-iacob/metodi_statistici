// CLT: Teorema del limite centrale

void tela(
	  int n_mis,       //numero di misure effettuate durante l'esperimento (piani del rivelatori)
	  int n_exp=10000  //quante volte ripeto l'esperimento per studiare distribuzione
	  ){
  //gStyle->SetOptStat(1100); //style del box con la statistica (1100: mostra Integral + Nbins + Name + Title)

  TFile *f = new TFile(Form("out_exp%d_mis%d.root",n_exp,n_mis), "RECREATE");
  
  TRandom * a = new TRandom();

  TH1F* h = new TH1F("h", Form("Hist - Original pdf (n_exp = %d, n_mis = %d);x_{i};#", n_exp, n_mis),500,-5,30); //hist per i dati (originale -> originale)
  h->GetXaxis()->SetTitleFont(62);
  
  
  TH1F* h_s = new TH1F("h_s","Hist - Sum (CLT);#Sigma_{i} x_{i};#",500,-5,40); //hist per la somma (originale -> gaussiana)
  h_s->GetXaxis()->SetTitleFont(62);

  TF1 *fVav = new TF1("fVav","TMath::Vavilov(x,[0],[1])",-5,30);
  fVav->SetParameters(0.1,0.2);
 
  double x_i; //singola misura
  for (int k=0; k<n_exp; k++){
    double sum = 0.;
    for (int N=0; N<n_mis; N++){
      //x_i: misura i-esima della variabile che segue una distrubuzione originale, mu_i e sigma_i scelte
      //x_i = a->Uniform(0,2); //Uniforme (hist limits = 0,2 et 0,3)
      //x_i = a->Landau(0,1);  //Landau (hist limits = -5,30 et -5,40)
      x_i = fVav->GetRandom();  //Vavilov (precedentemente definita)
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

  f->Write();
}


void clt(){
  tela(1);
  tela(2);
  tela(5);
  tela(100);
}
