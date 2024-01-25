void getEfficiency(string fileInName="", string fileOutName="")
{
  TFile *fIn = new TFile(fileInName.c_str(), "read");
  TFile *fOut = new TFile(fileOutName.c_str(), "recreate");

  TH2D *h2_recYPt_proton = (TH2D*) fIn->Get("h2_recYPt_proton");
  TH2D *h2_recYPt_pionP  = (TH2D*) fIn->Get("h2_recYPt_pionP");
  TH2D *h2_recYPt_pionM  = (TH2D*) fIn->Get("h2_recYPt_pionM");
  TH2D *h2_recYPt_pions  = (TH2D*) fIn->Get("h2_recYPt_pions");
  TH2D *h2_recYPt_kaonP  = (TH2D*) fIn->Get("h2_recYPt_kaonP");
  TH2D *h2_recYPt_kaonM  = (TH2D*) fIn->Get("h2_recYPt_kaonM");
  TH2D *h2_recYPt_kaons  = (TH2D*) fIn->Get("h2_recYPt_kaons");
  TH2D *h2_simYPt_proton = (TH2D*) fIn->Get("h2_simYPt_proton");
  TH2D *h2_simYPt_pionP  = (TH2D*) fIn->Get("h2_simYPt_pionP");
  TH2D *h2_simYPt_pionM  = (TH2D*) fIn->Get("h2_simYPt_pionM");
  TH2D *h2_simYPt_pions  = (TH2D*) fIn->Get("h2_simYPt_pions");
  TH2D *h2_simYPt_kaonP  = (TH2D*) fIn->Get("h2_simYPt_kaonP");
  TH2D *h2_simYPt_kaonM  = (TH2D*) fIn->Get("h2_simYPt_kaonM");
  TH2D *h2_simYPt_kaons  = (TH2D*) fIn->Get("h2_simYPt_kaons");

  if(!h2_recYPt_proton) return;
  if(!h2_recYPt_pionP) return;
  if(!h2_recYPt_pionM) return;
  if(!h2_recYPt_pions) return;
  if(!h2_recYPt_kaonP) return;
  if(!h2_recYPt_kaonM) return;
  if(!h2_recYPt_kaons) return;
  if(!h2_simYPt_proton) return;
  if(!h2_simYPt_pionP) return;
  if(!h2_simYPt_pionM) return;
  if(!h2_simYPt_pions) return;
  if(!h2_simYPt_kaonP) return;
  if(!h2_simYPt_kaonM) return;
  if(!h2_simYPt_kaons) return;
  
  TH2D *h2_effYPt_proton = (TH2D*) h2_recYPt_proton->Clone();
  h2_effYPt_proton->SetName("h2_effYPt_proton");
  h2_effYPt_proton->SetTitle("Efficiency (Y-pT) of primary protons;y;p_{T} (GeV/c)");
  h2_effYPt_proton->Divide(h2_simYPt_proton);
  TH2D *h2_effYPt_pionP = (TH2D*) h2_recYPt_pionP->Clone();
  h2_effYPt_pionP->SetName("h2_effYPt_pionP");
  h2_effYPt_pionP->SetTitle("Efficiency (Y-pT) of primary #pi^{+};y;p_{T} (GeV/c)");
  h2_effYPt_pionP->Divide(h2_simYPt_pionP);
  TH2D *h2_effYPt_pionM = (TH2D*) h2_recYPt_pionM->Clone();
  h2_effYPt_pionM->SetName("h2_effYPt_pionM");
  h2_effYPt_pionM->SetTitle("Efficiency (Y-pT) of primary #pi^{-};y;p_{T} (GeV/c)");
  h2_effYPt_pionM->Divide(h2_simYPt_pionM);
  TH2D *h2_effYPt_pions = (TH2D*) h2_recYPt_pions->Clone();
  h2_effYPt_pions->SetName("h2_effYPt_pions");
  h2_effYPt_pions->SetTitle("Efficiency (Y-pT) of primary #pi^{#pm};y;p_{T} (GeV/c)");
  h2_effYPt_pions->Divide(h2_simYPt_pions);
  TH2D *h2_effYPt_kaonP = (TH2D*) h2_recYPt_kaonP->Clone();
  h2_effYPt_kaonP->SetName("h2_effYPt_kaonP");
  h2_effYPt_kaonP->SetTitle("Efficiency (Y-pT) of primary K^{+};y;p_{T} (GeV/c)");
  h2_effYPt_kaonP->Divide(h2_simYPt_kaonP);
  TH2D *h2_effYPt_kaonM = (TH2D*) h2_recYPt_kaonM->Clone();
  h2_effYPt_kaonM->SetName("h2_effYPt_kaonM");
  h2_effYPt_kaonM->SetTitle("Efficiency (Y-pT) of primary K^{-};y;p_{T} (GeV/c)");
  h2_effYPt_kaonM->Divide(h2_simYPt_kaonM);
  TH2D *h2_effYPt_kaons = (TH2D*) h2_recYPt_kaons->Clone();
  h2_effYPt_kaons->SetName("h2_effYPt_kaons");
  h2_effYPt_kaons->SetTitle("Efficiency (Y-pT) of primary K^{#pm};y;p_{T} (GeV/c)");
  h2_effYPt_kaons->Divide(h2_simYPt_kaons);

  fOut->cd();
  h2_effYPt_proton->Write();
  h2_effYPt_pionP->Write();
  h2_effYPt_pionM->Write();
  h2_effYPt_pions->Write();
  h2_effYPt_kaonP->Write();
  h2_effYPt_kaonM->Write();
  h2_effYPt_kaons->Write();
  fOut->Close();
}
