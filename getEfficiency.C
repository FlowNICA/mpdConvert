void getEfficiency(string fileInName="", string fileOutName="")
{
  const int niter = 5;
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

  TH2D *h2_recDcaXPt_proton = (TH2D*) fIn->Get("h2_recDcaXPt_proton");
  TH2D *h2_recDcaXPt_pionP = (TH2D*) fIn->Get("h2_recDcaXPt_pionP");
  TH2D *h2_recDcaXPt_pionM = (TH2D*) fIn->Get("h2_recDcaXPt_pionM");
  TH2D *h2_recDcaXPt_pions = (TH2D*) fIn->Get("h2_recDcaXPt_pions");
  TH2D *h2_recDcaXPt_kaonP = (TH2D*) fIn->Get("h2_recDcaXPt_kaonP");
  TH2D *h2_recDcaXPt_kaonM = (TH2D*) fIn->Get("h2_recDcaXPt_kaonM");
  TH2D *h2_recDcaXPt_kaons = (TH2D*) fIn->Get("h2_recDcaXPt_kaons");
  TH2D *h2_recDcaYPt_proton = (TH2D*) fIn->Get("h2_recDcaYPt_proton");
  TH2D *h2_recDcaYPt_pionP = (TH2D*) fIn->Get("h2_recDcaYPt_pionP");
  TH2D *h2_recDcaYPt_pionM = (TH2D*) fIn->Get("h2_recDcaYPt_pionM");
  TH2D *h2_recDcaYPt_pions = (TH2D*) fIn->Get("h2_recDcaYPt_pions");
  TH2D *h2_recDcaYPt_kaonP = (TH2D*) fIn->Get("h2_recDcaYPt_kaonP");
  TH2D *h2_recDcaYPt_kaonM = (TH2D*) fIn->Get("h2_recDcaYPt_kaonM");
  TH2D *h2_recDcaYPt_kaons = (TH2D*) fIn->Get("h2_recDcaYPt_kaons");
  TH2D *h2_recDcaZPt_proton = (TH2D*) fIn->Get("h2_recDcaZPt_proton");
  TH2D *h2_recDcaZPt_pionP = (TH2D*) fIn->Get("h2_recDcaZPt_pionP");
  TH2D *h2_recDcaZPt_pionM = (TH2D*) fIn->Get("h2_recDcaZPt_pionM");
  TH2D *h2_recDcaZPt_pions = (TH2D*) fIn->Get("h2_recDcaZPt_pions");
  TH2D *h2_recDcaZPt_kaonP = (TH2D*) fIn->Get("h2_recDcaZPt_kaonP");
  TH2D *h2_recDcaZPt_kaonM = (TH2D*) fIn->Get("h2_recDcaZPt_kaonM");
  TH2D *h2_recDcaZPt_kaons = (TH2D*) fIn->Get("h2_recDcaZPt_kaons");

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

  if (!h2_recDcaXPt_proton) return;
  if (!h2_recDcaXPt_pionP) return;
  if (!h2_recDcaXPt_pionM) return;
  if (!h2_recDcaXPt_pions) return;
  if (!h2_recDcaXPt_kaonP) return;
  if (!h2_recDcaXPt_kaonM) return;
  if (!h2_recDcaXPt_kaons) return;
  if (!h2_recDcaYPt_proton) return;
  if (!h2_recDcaYPt_pionP) return;
  if (!h2_recDcaYPt_pionM) return;
  if (!h2_recDcaYPt_pions) return;
  if (!h2_recDcaYPt_kaonP) return;
  if (!h2_recDcaYPt_kaonM) return;
  if (!h2_recDcaYPt_kaons) return;
  if (!h2_recDcaZPt_proton) return;
  if (!h2_recDcaZPt_pionP) return;
  if (!h2_recDcaZPt_pionM) return;
  if (!h2_recDcaZPt_pions) return;
  if (!h2_recDcaZPt_kaonP) return;
  if (!h2_recDcaZPt_kaonM) return;
  if (!h2_recDcaZPt_kaons) return;

  // Prepare efficiency (Ycm-Pt) maps
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

  
  // Calculate n-sigma for DCA plots
  int firstbin;
  std::vector<TH1D*> vh1_dcax_proton, vh1_dcay_proton, vh1_dcaz_proton; // vector of DCA (x,y,z) in pT bins for protons
  std::vector<TH1D*> vh1_dcax_pionP, vh1_dcay_pionP, vh1_dcaz_pionP; // vector of DCA (x,y,z) in pT bins for pionP
  std::vector<TH1D*> vh1_dcax_pionM, vh1_dcay_pionM, vh1_dcaz_pionM; // vector of DCA (x,y,z) in pT bins for pionM
  std::vector<TH1D*> vh1_dcax_pions, vh1_dcay_pions, vh1_dcaz_pions; // vector of DCA (x,y,z) in pT bins for pions
  std::vector<TH1D*> vh1_dcax_kaonP, vh1_dcay_kaonP, vh1_dcaz_kaonP; // vector of DCA (x,y,z) in pT bins for kaonP
  std::vector<TH1D*> vh1_dcax_kaonM, vh1_dcay_kaonM, vh1_dcaz_kaonM; // vector of DCA (x,y,z) in pT bins for kaonM
  std::vector<TH1D*> vh1_dcax_kaons, vh1_dcay_kaons, vh1_dcaz_kaons; // vector of DCA (x,y,z) in pT bins for kaons
  std::vector<TF1*> vf1_dcax_proton, vf1_dcay_proton, vf1_dcaz_proton;
  std::vector<TF1*> vf1_dcax_pionP, vf1_dcay_pionP, vf1_dcaz_pionP;
  std::vector<TF1*> vf1_dcax_pionM, vf1_dcay_pionM, vf1_dcaz_pionM;
  std::vector<TF1*> vf1_dcax_pions, vf1_dcay_pions, vf1_dcaz_pions;
  std::vector<TF1*> vf1_dcax_kaonP, vf1_dcay_kaonP, vf1_dcaz_kaonP;
  std::vector<TF1*> vf1_dcax_kaonM, vf1_dcay_kaonM, vf1_dcaz_kaonM;
  std::vector<TF1*> vf1_dcax_kaons, vf1_dcay_kaons, vf1_dcaz_kaons;
  std::vector<double> vsigm_dcax_proton, vsigm_dcay_proton, vsigm_dcaz_proton;
  std::vector<double> vsigm_dcax_pionP, vsigm_dcay_pionP, vsigm_dcaz_pionP;
  std::vector<double> vsigm_dcax_pionM, vsigm_dcay_pionM, vsigm_dcaz_pionM;
  std::vector<double> vsigm_dcax_pions, vsigm_dcay_pions, vsigm_dcaz_pions;
  std::vector<double> vsigm_dcax_kaonP, vsigm_dcay_kaonP, vsigm_dcaz_kaonP;
  std::vector<double> vsigm_dcax_kaonM, vsigm_dcay_kaonM, vsigm_dcaz_kaonM;
  std::vector<double> vsigm_dcax_kaons, vsigm_dcay_kaons, vsigm_dcaz_kaons;
  std::vector<double> vsigm_err_dcax_proton, vsigm_err_dcay_proton, vsigm_err_dcaz_proton;
  std::vector<double> vsigm_err_dcax_pionP, vsigm_err_dcay_pionP, vsigm_err_dcaz_pionP;
  std::vector<double> vsigm_err_dcax_pionM, vsigm_err_dcay_pionM, vsigm_err_dcaz_pionM;
  std::vector<double> vsigm_err_dcax_pions, vsigm_err_dcay_pions, vsigm_err_dcaz_pions;
  std::vector<double> vsigm_err_dcax_kaonP, vsigm_err_dcay_kaonP, vsigm_err_dcaz_kaonP;
  std::vector<double> vsigm_err_dcax_kaonM, vsigm_err_dcay_kaonM, vsigm_err_dcaz_kaonM;
  std::vector<double> vsigm_err_dcax_kaons, vsigm_err_dcay_kaons, vsigm_err_dcaz_kaons;
  std::vector<double> vpt_x_proton, vpt_err_x_proton;
  std::vector<double> vpt_x_pion, vpt_err_x_pion;

  auto h2_dcax_proton = (TH2D*)h2_recDcaXPt_proton->Clone();
  h2_dcax_proton->RebinY(10);
  firstbin = h2_dcax_proton->GetYaxis()->FindBin(0.4);
  for (int i=firstbin; i<h2_dcax_proton->GetNbinsY(); ++i) {
    vh1_dcax_proton.push_back( (TH1D*)h2_dcax_proton->ProjectionX(Form("h1_%s_dcax_proton_bin%i", h2_dcax_proton->GetName(), i-firstbin), i, i) );
    vh1_dcax_proton.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcax_proton->GetTitle(),
                                          h2_dcax_proton->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcax_proton->GetXaxis()->GetBinUpEdge(i)) );

    vpt_x_proton.push_back(h2_dcax_proton->GetYaxis()->GetBinCenter(i));
    vpt_err_x_proton.push_back(0.);
    vf1_dcax_proton.push_back( new TF1(Form("f1_dcax_proton_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcax_proton.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcay_proton = (TH2D*)h2_recDcaYPt_proton->Clone();
  h2_dcay_proton->RebinY(10);
  firstbin = h2_dcay_proton->GetYaxis()->FindBin(0.4);
  for (int i=firstbin; i<h2_dcay_proton->GetNbinsY(); ++i) {
    vh1_dcay_proton.push_back( (TH1D*)h2_dcay_proton->ProjectionX(Form("h1_%s_dcay_proton_bin%i", h2_dcay_proton->GetName(), i-firstbin), i, i) );
    vh1_dcay_proton.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcay_proton->GetTitle(),
                                          h2_dcay_proton->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcay_proton->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcay_proton.push_back( new TF1(Form("f1_dcay_proton_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcay_proton.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcaz_proton = (TH2D*)h2_recDcaZPt_proton->Clone();
  h2_dcaz_proton->RebinY(10);
  firstbin = h2_dcaz_proton->GetYaxis()->FindBin(0.4);
  for (int i=firstbin; i<h2_dcaz_proton->GetNbinsY(); ++i) {
    vh1_dcaz_proton.push_back( (TH1D*)h2_dcaz_proton->ProjectionX(Form("h1_%s_dcaz_proton_bin%i", h2_dcaz_proton->GetName(), i-firstbin), i, i) );
    vh1_dcaz_proton.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcaz_proton->GetTitle(),
                                          h2_dcaz_proton->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcaz_proton->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcaz_proton.push_back( new TF1(Form("f1_dcaz_proton_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcaz_proton.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcax_pionP = (TH2D*)h2_recDcaXPt_pionP->Clone();
  h2_dcax_pionP->RebinY(10);
  firstbin = h2_dcax_pionP->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcax_pionP->GetNbinsY(); ++i) {
    vh1_dcax_pionP.push_back( (TH1D*)h2_dcax_pionP->ProjectionX(Form("h1_%s_dcax_pionP_bin%i", h2_dcax_pionP->GetName(), i-firstbin), i, i) );
    vh1_dcax_pionP.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcax_pionP->GetTitle(),
                                          h2_dcax_pionP->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcax_pionP->GetXaxis()->GetBinUpEdge(i)) );

    vpt_x_pion.push_back(h2_dcax_pionP->GetYaxis()->GetBinCenter(i));
    vpt_err_x_pion.push_back(0.);
    vf1_dcax_pionP.push_back( new TF1(Form("f1_dcax_pionP_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcax_pionP.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcay_pionP = (TH2D*)h2_recDcaYPt_pionP->Clone();
  h2_dcay_pionP->RebinY(10);
  firstbin = h2_dcay_pionP->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcay_pionP->GetNbinsY(); ++i) {
    vh1_dcay_pionP.push_back( (TH1D*)h2_dcay_pionP->ProjectionX(Form("h1_%s_dcay_pionP_bin%i", h2_dcay_pionP->GetName(), i-firstbin), i, i) );
    vh1_dcay_pionP.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcay_pionP->GetTitle(),
                                          h2_dcay_pionP->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcay_pionP->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcay_pionP.push_back( new TF1(Form("f1_dcay_pionP_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcay_pionP.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcaz_pionP = (TH2D*)h2_recDcaZPt_pionP->Clone();
  h2_dcaz_pionP->RebinY(10);
  firstbin = h2_dcaz_pionP->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcaz_pionP->GetNbinsY(); ++i) {
    vh1_dcaz_pionP.push_back( (TH1D*)h2_dcaz_pionP->ProjectionX(Form("h1_%s_dcaz_pionP_bin%i", h2_dcaz_pionP->GetName(), i-firstbin), i, i) );
    vh1_dcaz_pionP.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcaz_pionP->GetTitle(),
                                          h2_dcaz_pionP->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcaz_pionP->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcaz_pionP.push_back( new TF1(Form("f1_dcaz_pionP_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcaz_pionP.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcax_pionM = (TH2D*)h2_recDcaXPt_pionM->Clone();
  h2_dcax_pionM->RebinY(10);
  firstbin = h2_dcax_pionM->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcax_pionM->GetNbinsY(); ++i) {
    vh1_dcax_pionM.push_back( (TH1D*)h2_dcax_pionM->ProjectionX(Form("h1_%s_dcax_pionM_bin%i", h2_dcax_pionM->GetName(), i-firstbin), i, i) );
    vh1_dcax_pionM.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcax_pionM->GetTitle(),
                                          h2_dcax_pionM->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcax_pionM->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcax_pionM.push_back( new TF1(Form("f1_dcax_pionM_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcax_pionM.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcay_pionM = (TH2D*)h2_recDcaYPt_pionM->Clone();
  h2_dcay_pionM->RebinY(10);
  firstbin = h2_dcay_pionM->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcay_pionM->GetNbinsY(); ++i) {
    vh1_dcay_pionM.push_back( (TH1D*)h2_dcay_pionM->ProjectionX(Form("h1_%s_dcay_pionM_bin%i", h2_dcay_pionM->GetName(), i-firstbin), i, i) );
    vh1_dcay_pionM.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcay_pionM->GetTitle(),
                                          h2_dcay_pionM->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcay_pionM->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcay_pionM.push_back( new TF1(Form("f1_dcay_pionM_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcay_pionM.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcaz_pionM = (TH2D*)h2_recDcaZPt_pionM->Clone();
  h2_dcaz_pionM->RebinY(10);
  firstbin = h2_dcaz_pionM->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcaz_pionM->GetNbinsY(); ++i) {
    vh1_dcaz_pionM.push_back( (TH1D*)h2_dcaz_pionM->ProjectionX(Form("h1_%s_dcaz_pionM_bin%i", h2_dcaz_pionM->GetName(), i-firstbin), i, i) );
    vh1_dcaz_pionM.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcaz_pionM->GetTitle(),
                                          h2_dcaz_pionM->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcaz_pionM->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcaz_pionM.push_back( new TF1(Form("f1_dcaz_pionM_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcaz_pionM.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcax_pions = (TH2D*)h2_recDcaXPt_pions->Clone();
  h2_dcax_pions->RebinY(10);
  firstbin = h2_dcax_pions->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcax_pions->GetNbinsY(); ++i) {
    vh1_dcax_pions.push_back( (TH1D*)h2_dcax_pions->ProjectionX(Form("h1_%s_dcax_pions_bin%i", h2_dcax_pions->GetName(), i-firstbin), i, i) );
    vh1_dcax_pions.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcax_pions->GetTitle(),
                                          h2_dcax_pions->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcax_pions->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcax_pions.push_back( new TF1(Form("f1_dcax_pions_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcax_pions.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcay_pions = (TH2D*)h2_recDcaYPt_pions->Clone();
  h2_dcay_pions->RebinY(10);
  firstbin = h2_dcay_pions->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcay_pions->GetNbinsY(); ++i) {
    vh1_dcay_pions.push_back( (TH1D*)h2_dcay_pions->ProjectionX(Form("h1_%s_dcay_pions_bin%i", h2_dcay_pions->GetName(), i-firstbin), i, i) );
    vh1_dcay_pions.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcay_pions->GetTitle(),
                                          h2_dcay_pions->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcay_pions->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcay_pions.push_back( new TF1(Form("f1_dcay_pions_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcay_pions.at(i-firstbin)->FixParameter(1, 0.);
  }
  auto h2_dcaz_pions = (TH2D*)h2_recDcaZPt_pions->Clone();
  h2_dcaz_pions->RebinY(10);
  firstbin = h2_dcaz_pions->GetYaxis()->FindBin(0.1);
  for (int i=firstbin; i<h2_dcaz_pions->GetNbinsY(); ++i) {
    vh1_dcaz_pions.push_back( (TH1D*)h2_dcaz_pions->ProjectionX(Form("h1_%s_dcaz_pions_bin%i", h2_dcaz_pions->GetName(), i-firstbin), i, i) );
    vh1_dcaz_pions.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p_{T} < %2.2f GeV/c",
                                          h2_dcaz_pions->GetTitle(),
                                          h2_dcaz_pions->GetXaxis()->GetBinLowEdge(i),
                                          h2_dcaz_pions->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dcaz_pions.push_back( new TF1(Form("f1_dcaz_pions_ptbin%i", i-firstbin), "gaus", -0.1, 0.1) );
    vf1_dcaz_pions.at(i-firstbin)->FixParameter(1, 0.);
  }

  // Fitting dca distributions
  std::cout << "Fitting dca distributions" << std::endl;
  for (int i=0; i<vh1_dcax_proton.size(); ++i){
    if (!vh1_dcax_proton.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcax_proton.at(i)->Fit(vf1_dcax_proton.at(i), "RNQ");
    }
    vsigm_dcax_proton.push_back(vf1_dcax_proton.at(i)->GetParameter(2));
    vsigm_err_dcax_proton.push_back(vf1_dcax_proton.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcay_proton.size(); ++i){
    if (!vh1_dcay_proton.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcay_proton.at(i)->Fit(vf1_dcay_proton.at(i), "RNQ");
    }
    vsigm_dcay_proton.push_back(vf1_dcay_proton.at(i)->GetParameter(2));
    vsigm_err_dcay_proton.push_back(vf1_dcay_proton.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcaz_proton.size(); ++i){
    if (!vh1_dcaz_proton.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcaz_proton.at(i)->Fit(vf1_dcaz_proton.at(i), "RNQ");
    }
    vsigm_dcaz_proton.push_back(vf1_dcaz_proton.at(i)->GetParameter(2));
    vsigm_err_dcaz_proton.push_back(vf1_dcaz_proton.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcax_pionP.size(); ++i){
    if (!vh1_dcax_pionP.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcax_pionP.at(i)->Fit(vf1_dcax_pionP.at(i), "RNQ");
    }
    vsigm_dcax_pionP.push_back(vf1_dcax_pionP.at(i)->GetParameter(2));
    vsigm_err_dcax_pionP.push_back(vf1_dcax_pionP.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcay_pionP.size(); ++i){
    if (!vh1_dcay_pionP.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcay_pionP.at(i)->Fit(vf1_dcay_pionP.at(i), "RNQ");
    }
    vsigm_dcay_pionP.push_back(vf1_dcay_pionP.at(i)->GetParameter(2));
    vsigm_err_dcay_pionP.push_back(vf1_dcay_pionP.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcaz_pionP.size(); ++i){
    if (!vh1_dcaz_pionP.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcaz_pionP.at(i)->Fit(vf1_dcaz_pionP.at(i), "RNQ");
    }
    vsigm_dcaz_pionP.push_back(vf1_dcaz_pionP.at(i)->GetParameter(2));
    vsigm_err_dcaz_pionP.push_back(vf1_dcaz_pionP.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcax_pionM.size(); ++i){
    if (!vh1_dcax_pionM.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcax_pionM.at(i)->Fit(vf1_dcax_pionM.at(i), "RNQ");
    }
    vsigm_dcax_pionM.push_back(vf1_dcax_pionM.at(i)->GetParameter(2));
    vsigm_err_dcax_pionM.push_back(vf1_dcax_pionM.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcay_pionM.size(); ++i){
    if (!vh1_dcay_pionM.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcay_pionM.at(i)->Fit(vf1_dcay_pionM.at(i), "RNQ");
    }
    vsigm_dcay_pionM.push_back(vf1_dcay_pionM.at(i)->GetParameter(2));
    vsigm_err_dcay_pionM.push_back(vf1_dcay_pionM.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcaz_pionM.size(); ++i){
    if (!vh1_dcaz_pionM.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcaz_pionM.at(i)->Fit(vf1_dcaz_pionM.at(i), "RNQ");
    }
    vsigm_dcaz_pionM.push_back(vf1_dcaz_pionM.at(i)->GetParameter(2));
    vsigm_err_dcaz_pionM.push_back(vf1_dcaz_pionM.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcax_pions.size(); ++i){
    if (!vh1_dcax_pions.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcax_pions.at(i)->Fit(vf1_dcax_pions.at(i), "RNQ");
    }
    vsigm_dcax_pions.push_back(vf1_dcax_pions.at(i)->GetParameter(2));
    vsigm_err_dcax_pions.push_back(vf1_dcax_pions.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcay_pions.size(); ++i){
    if (!vh1_dcay_pions.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcay_pions.at(i)->Fit(vf1_dcay_pions.at(i), "RNQ");
    }
    vsigm_dcay_pions.push_back(vf1_dcay_pions.at(i)->GetParameter(2));
    vsigm_err_dcay_pions.push_back(vf1_dcay_pions.at(i)->GetParError(2));
  }
  for (int i=0; i<vh1_dcaz_pions.size(); ++i){
    if (!vh1_dcaz_pions.at(i)) continue;
    for (int k=0; k<niter; ++k){
      vh1_dcaz_pions.at(i)->Fit(vf1_dcaz_pions.at(i), "RNQ");
    }
    vsigm_dcaz_pions.push_back(vf1_dcaz_pions.at(i)->GetParameter(2));
    vsigm_err_dcaz_pions.push_back(vf1_dcaz_pions.at(i)->GetParError(2));
  }

  // Make graphs and functions with DCA n-sigma
  auto gr_dcax_sigm_proton = new TGraphErrors(vpt_x_proton.size()-1,
                                              vpt_x_proton.data(), vsigm_dcax_proton.data(),
                                              vpt_err_x_proton.data(), vsigm_err_dcax_proton.data());
  gr_dcax_sigm_proton->SetName("gr_dcax_sigm_proton");
  gr_dcax_sigm_proton->SetTitle("#sigma_{DCA_{x}} of protons vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{x}}, cm");
  auto f1_dcax_sigm_proton = new TF1("f1_dcax_sigm_proton", "pol6", 0.4, 3.);
  gr_dcax_sigm_proton->Fit(f1_dcax_sigm_proton, "RNQ");
  auto gr_dcay_sigm_proton = new TGraphErrors(vpt_x_proton.size()-1,
                                              vpt_x_proton.data(), vsigm_dcay_proton.data(),
                                              vpt_err_x_proton.data(), vsigm_err_dcay_proton.data());
  gr_dcay_sigm_proton->SetName("gr_dcay_sigm_proton");
  gr_dcay_sigm_proton->SetTitle("#sigma_{DCA_{y}} of protons vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{y}}, cm");
  auto f1_dcay_sigm_proton = new TF1("f1_dcay_sigm_proton", "pol6", 0.4, 3.);
  gr_dcay_sigm_proton->Fit(f1_dcay_sigm_proton, "RNQ");
  auto gr_dcaz_sigm_proton = new TGraphErrors(vpt_x_proton.size()-1,
                                              vpt_x_proton.data(), vsigm_dcaz_proton.data(),
                                              vpt_err_x_proton.data(), vsigm_err_dcaz_proton.data());
  gr_dcaz_sigm_proton->SetName("gr_dcaz_sigm_proton");
  gr_dcaz_sigm_proton->SetTitle("#sigma_{DCA_{z}} of protons vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{z}}, cm");
  auto f1_dcaz_sigm_proton = new TF1("f1_dcaz_sigm_proton", "pol6", 0.4, 3.);
  gr_dcaz_sigm_proton->Fit(f1_dcaz_sigm_proton, "RNQ");
  auto gr_dcax_sigm_pionP = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcax_pionP.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcax_pionP.data());
  gr_dcax_sigm_pionP->SetName("gr_dcax_sigm_pionP");
  gr_dcax_sigm_pionP->SetTitle("#sigma_{DCA_{x}} of pionPs vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{x}}, cm");
  auto f1_dcax_sigm_pionP = new TF1("f1_dcax_sigm_pionP", "pol6", 0.1, 1.5);
  gr_dcax_sigm_pionP->Fit(f1_dcax_sigm_pionP, "RNQ");
  auto gr_dcay_sigm_pionP = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcay_pionP.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcay_pionP.data());
  gr_dcay_sigm_pionP->SetName("gr_dcay_sigm_pionP");
  gr_dcay_sigm_pionP->SetTitle("#sigma_{DCA_{y}} of pionPs vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{y}}, cm");
  auto f1_dcay_sigm_pionP = new TF1("f1_dcay_sigm_pionP", "pol6", 0.1, 1.5);
  gr_dcay_sigm_pionP->Fit(f1_dcay_sigm_pionP, "RNQ");
  auto gr_dcaz_sigm_pionP = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcaz_pionP.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcaz_pionP.data());
  gr_dcaz_sigm_pionP->SetName("gr_dcaz_sigm_pionP");
  gr_dcaz_sigm_pionP->SetTitle("#sigma_{DCA_{z}} of pionPs vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{z}}, cm");
  auto f1_dcaz_sigm_pionP = new TF1("f1_dcaz_sigm_pionP", "pol6", 0.1, 1.5);
  gr_dcaz_sigm_pionP->Fit(f1_dcaz_sigm_pionP, "RNQ");
  auto gr_dcax_sigm_pionM = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcax_pionM.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcax_pionM.data());
  gr_dcax_sigm_pionM->SetName("gr_dcax_sigm_pionM");
  gr_dcax_sigm_pionM->SetTitle("#sigma_{DCA_{x}} of pionMs vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{x}}, cm");
  auto f1_dcax_sigm_pionM = new TF1("f1_dcax_sigm_pionM", "pol6", 0.1, 1.5);
  gr_dcax_sigm_pionM->Fit(f1_dcax_sigm_pionM, "RNQ");
  auto gr_dcay_sigm_pionM = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcay_pionM.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcay_pionM.data());
  gr_dcay_sigm_pionM->SetName("gr_dcay_sigm_pionM");
  gr_dcay_sigm_pionM->SetTitle("#sigma_{DCA_{y}} of pionMs vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{y}}, cm");
  auto f1_dcay_sigm_pionM = new TF1("f1_dcay_sigm_pionM", "pol6", 0.1, 1.5);
  gr_dcay_sigm_pionM->Fit(f1_dcay_sigm_pionM, "RNQ");
  auto gr_dcaz_sigm_pionM = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcaz_pionM.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcaz_pionM.data());
  gr_dcaz_sigm_pionM->SetName("gr_dcaz_sigm_pionM");
  gr_dcaz_sigm_pionM->SetTitle("#sigma_{DCA_{z}} of pionMs vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{z}}, cm");
  auto f1_dcaz_sigm_pionM = new TF1("f1_dcaz_sigm_pionM", "pol6", 0.1, 1.5);
  gr_dcaz_sigm_pionM->Fit(f1_dcaz_sigm_pionM, "RNQ");
  auto gr_dcax_sigm_pions = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcax_pions.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcax_pions.data());
  gr_dcax_sigm_pions->SetName("gr_dcax_sigm_pions");
  gr_dcax_sigm_pions->SetTitle("#sigma_{DCA_{x}} of pionss vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{x}}, cm");
  auto f1_dcax_sigm_pions = new TF1("f1_dcax_sigm_pions", "pol6", 0.1, 1.5);
  gr_dcax_sigm_pions->Fit(f1_dcax_sigm_pions, "RNQ");
  auto gr_dcay_sigm_pions = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcay_pions.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcay_pions.data());
  gr_dcay_sigm_pions->SetName("gr_dcay_sigm_pions");
  gr_dcay_sigm_pions->SetTitle("#sigma_{DCA_{y}} of pionss vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{y}}, cm");
  auto f1_dcay_sigm_pions = new TF1("f1_dcay_sigm_pions", "pol6", 0.1, 1.5);
  gr_dcay_sigm_pions->Fit(f1_dcay_sigm_pions, "RNQ");
  auto gr_dcaz_sigm_pions = new TGraphErrors(vpt_x_pion.size()-1,
                                              vpt_x_pion.data(), vsigm_dcaz_pions.data(),
                                              vpt_err_x_pion.data(), vsigm_err_dcaz_pions.data());
  gr_dcaz_sigm_pions->SetName("gr_dcaz_sigm_pions");
  gr_dcaz_sigm_pions->SetTitle("#sigma_{DCA_{z}} of pionss vs p_{T};p_{T}, GeV/c;#sigma_{DCA_{z}}, cm");
  auto f1_dcaz_sigm_pions = new TF1("f1_dcaz_sigm_pions", "pol6", 0.1, 1.5);
  gr_dcaz_sigm_pions->Fit(f1_dcaz_sigm_pions, "RNQ");

  fOut->cd();

  h2_effYPt_proton->Write();
  h2_effYPt_pionP->Write();
  h2_effYPt_pionM->Write();
  h2_effYPt_pions->Write();
  h2_effYPt_kaonP->Write();
  h2_effYPt_kaonM->Write();
  h2_effYPt_kaons->Write();

  h2_dcax_proton->Write();
  gr_dcax_sigm_proton->Write();
  f1_dcax_sigm_proton->Write();
  h2_dcay_proton->Write();
  gr_dcay_sigm_proton->Write();
  f1_dcay_sigm_proton->Write();
  h2_dcaz_proton->Write();
  gr_dcaz_sigm_proton->Write();
  f1_dcaz_sigm_proton->Write();
  h2_dcax_pionP->Write();
  gr_dcax_sigm_pionP->Write();
  f1_dcax_sigm_pionP->Write();
  h2_dcay_pionP->Write();
  gr_dcay_sigm_pionP->Write();
  f1_dcay_sigm_pionP->Write();
  h2_dcaz_pionP->Write();
  gr_dcaz_sigm_pionP->Write();
  f1_dcaz_sigm_pionP->Write();
  h2_dcax_pionM->Write();
  gr_dcax_sigm_pionM->Write();
  f1_dcax_sigm_pionM->Write();
  h2_dcay_pionM->Write();
  gr_dcay_sigm_pionM->Write();
  f1_dcay_sigm_pionM->Write();
  h2_dcaz_pionM->Write();
  gr_dcaz_sigm_pionM->Write();
  f1_dcaz_sigm_pionM->Write();
  h2_dcax_pions->Write();
  gr_dcax_sigm_pions->Write();
  f1_dcax_sigm_pions->Write();
  h2_dcay_pions->Write();
  gr_dcay_sigm_pions->Write();
  f1_dcay_sigm_pions->Write();
  h2_dcaz_pions->Write();
  gr_dcaz_sigm_pions->Write();
  f1_dcaz_sigm_pions->Write();

  fOut->mkdir("dcax_proton_slices");
  fOut->cd("dcax_proton_slices");
  for (int i=0; i<vh1_dcax_proton.size(); i++) {
    vh1_dcax_proton.at(i)->Write();
    vf1_dcax_proton.at(i)->Write();
  }
  fOut->mkdir("dcay_proton_slices");
  fOut->cd("dcay_proton_slices");
  for (int i=0; i<vh1_dcay_proton.size(); i++) {
    vh1_dcay_proton.at(i)->Write();
    vf1_dcay_proton.at(i)->Write();
  }
  fOut->mkdir("dcaz_proton_slices");
  fOut->cd("dcaz_proton_slices");
  for (int i=0; i<vh1_dcaz_proton.size(); i++) {
    vh1_dcaz_proton.at(i)->Write();
    vf1_dcaz_proton.at(i)->Write();
  }
  fOut->mkdir("dcax_pionP_slices");
  fOut->cd("dcax_pionP_slices");
  for (int i=0; i<vh1_dcax_pionP.size(); i++) {
    vh1_dcax_pionP.at(i)->Write();
    vf1_dcax_pionP.at(i)->Write();
  }
  fOut->mkdir("dcay_pionP_slices");
  fOut->cd("dcay_pionP_slices");
  for (int i=0; i<vh1_dcay_pionP.size(); i++) {
    vh1_dcay_pionP.at(i)->Write();
    vf1_dcay_pionP.at(i)->Write();
  }
  fOut->mkdir("dcaz_pionP_slices");
  fOut->cd("dcaz_pionP_slices");
  for (int i=0; i<vh1_dcaz_pionP.size(); i++) {
    vh1_dcaz_pionP.at(i)->Write();
    vf1_dcaz_pionP.at(i)->Write();
  }
  fOut->mkdir("dcax_pionM_slices");
  fOut->cd("dcax_pionM_slices");
  for (int i=0; i<vh1_dcax_pionM.size(); i++) {
    vh1_dcax_pionM.at(i)->Write();
    vf1_dcax_pionM.at(i)->Write();
  }
  fOut->mkdir("dcay_pionM_slices");
  fOut->cd("dcay_pionM_slices");
  for (int i=0; i<vh1_dcay_pionM.size(); i++) {
    vh1_dcay_pionM.at(i)->Write();
    vf1_dcay_pionM.at(i)->Write();
  }
  fOut->mkdir("dcaz_pionM_slices");
  fOut->cd("dcaz_pionM_slices");
  for (int i=0; i<vh1_dcaz_pionM.size(); i++) {
    vh1_dcaz_pionM.at(i)->Write();
    vf1_dcaz_pionM.at(i)->Write();
  }
  fOut->mkdir("dcax_pions_slices");
  fOut->cd("dcax_pions_slices");
  for (int i=0; i<vh1_dcax_pions.size(); i++) {
    vh1_dcax_pions.at(i)->Write();
    vf1_dcax_pions.at(i)->Write();
  }
  fOut->mkdir("dcay_pions_slices");
  fOut->cd("dcay_pions_slices");
  for (int i=0; i<vh1_dcay_pions.size(); i++) {
    vh1_dcay_pions.at(i)->Write();
    vf1_dcay_pions.at(i)->Write();
  }
  fOut->mkdir("dcaz_pions_slices");
  fOut->cd("dcaz_pions_slices");
  for (int i=0; i<vh1_dcaz_pions.size(); i++) {
    vh1_dcaz_pions.at(i)->Write();
    vf1_dcaz_pions.at(i)->Write();
  }

  fOut->Close();
}
