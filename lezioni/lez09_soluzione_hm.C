{//metodo più semplice per fare una macro in root, ma sconsigliato. Si lancia con root hit_and_miss.C

  TF1* f = new TF1("f","(5-0.8*x+1.5*exp(-pow(x-2.5,2)/0.05))",0,5);
  
  int N =1000000;
  double x, y; // queste sono le variabili random che voglio generare

  TH1F* h = new TH1F("h","Generated Data; x; y",100,0,5); //Inizializzo l'istogramma. "TH1F", istogramma ("H") unidimensionale ("1"), che riempirò con variabili Float ("F"). L'asse x va da 0 a 5 e sarà suddiviso in 100 bins.
  h->Sumw2();// in automatico tratta le entries come variabili Poissoniane. Vuol dire che se nel bin i-esimo ci sono N eventi, associa al bin i-esimo come errore √N

  for (int i=0; i<N; i++ ){ // genero i miei dati secondo la pdf f, utilizzando il metodo dell'"hit & miss"
    x = gRandom->Uniform( 0, 5); // variabile random x distribuita uniformemente nel range dell'asse x: 0,5.
    y = gRandom->Uniform( 0, f->GetMaximum());// variabile random y distribuita uniformemente nel range dell'asse y.
    if( y > f->Eval(x) ) continue; // chiedo se y è maggiore del valore assunto dalla pdf in x. Se sì, la escludo.
    h->Fill(x); // riempio l'istogramma con i dati generati
  }

  // abbellisco lo stile
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlue+1);
  h->SetLineColor(kBlue+1);
  f->SetLineWidth(2);
  f->SetLineColor(kRed+1);
  f->SetTitle("y = 5-0.8x+1.5e^{#frac{-(x-2.5)^2}{0.05}}");
  f->GetXaxis()->SetTitle("x");
  f->GetYaxis()->SetTitle("y");
  
  TCanvas *c = new TCanvas("c","c",700,1200); // Inizializzo la Canvas (la "tela") al cui interno metterò i plots
  c->Divide(1,2); //  Divido la canvas in 1 colonna e 2 righe
  c->cd(1); // seleziono la prima riga
  h->Draw(); // disegno l'istogramma nella prima riga della Canvas
  c->cd(2); // seleziono la seconda riga
  f->Draw(); // disegno la pdf nella seconda riga della Canvas
  

}
