void getDca(std::string inFile_qa_step2, std::string outFile) {
  TStopwatch timer;
  timer.Start();

  const int niter = 5;
  const int minentries = 10;

  std::cout << "Starting DCA fitting procedure" << std::endl;

  std::cout << "Reading input" << std::endl;
  TFile *fi = new TFile(inFile_qa_step2.c_str(), "read");
  if (!fi) {
    std::cerr << "Error: cannot open input file: " << inFile_qa_step2.c_str()
              << std::endl;
    return;
  }

  std::vector<std::unique_ptr<TH2D>> vh2_recDcaPt;
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaPt_proton"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaXPt_proton"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaYPt_proton"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaZPt_proton"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaPt_pionP"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaXPt_pionP"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaYPt_pionP"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaZPt_pionP"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaPt_pionM"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaXPt_pionM"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaYPt_pionM"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaZPt_pionM"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaPt_pions"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaXPt_pions"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaYPt_pions"))));
  vh2_recDcaPt.emplace_back(
      std::make_unique<TH2D>(*((TH2D *)fi->Get("h2_recDcaZPt_pions"))));

  for (auto &h : vh2_recDcaPt) {
    h->RebinY(10);
  }

  std::vector<std::unique_ptr<TH1D>>
      vh1_dca_prt, vh1_dcax_prt, vh1_dcay_prt, vh1_dcaz_prt; // vector of TH1D for DCA fit of protons
  std::vector<std::unique_ptr<TH1D>>
      vh1_dca_pip, vh1_dcax_pip, vh1_dcay_pip, vh1_dcaz_pip; // vector of TH1D for DCA fit of pi+
  std::vector<std::unique_ptr<TH1D>>
      vh1_dca_pim, vh1_dcax_pim, vh1_dcay_pim, vh1_dcaz_pim; // vector of TH1D for DCA fit of pi-
  std::vector<std::unique_ptr<TF1>>
      vf1_dca_prt, vf1_dcax_prt, vf1_dcay_prt, vf1_dcaz_prt; // vector of TF1 1gaus for DCA fit of protons
  std::vector<std::unique_ptr<TF1>>
      vf1_dca_pip, vf1_dcax_pip, vf1_dcay_pip, vf1_dcaz_pip; // vector of TF1 1gaus for DCA fit of pi+
  std::vector<std::unique_ptr<TF1>>
      vf1_dca_pim, vf1_dcax_pim, vf1_dcay_pim, vf1_dcaz_pim; // vector of TF1 1gaus for DCA fit of pi-

  std::vector<double> v_dcafit_x_prt, v_dcafit_ex_prt, v_dcafit_mean_prt,
      v_dcafit_sigm_prt, v_dcafit_emean_prt, v_dcafit_esigm_prt;
  std::vector<double> v_dcaxfit_x_prt, v_dcaxfit_ex_prt, v_dcaxfit_mean_prt,
      v_dcaxfit_sigm_prt, v_dcaxfit_emean_prt, v_dcaxfit_esigm_prt;
  std::vector<double> v_dcayfit_x_prt, v_dcayfit_ex_prt, v_dcayfit_mean_prt,
      v_dcayfit_sigm_prt, v_dcayfit_emean_prt, v_dcayfit_esigm_prt;
  std::vector<double> v_dcazfit_x_prt, v_dcazfit_ex_prt, v_dcazfit_mean_prt,
      v_dcazfit_sigm_prt, v_dcazfit_emean_prt, v_dcazfit_esigm_prt;
  std::vector<double> v_dcafit_x_pip, v_dcafit_ex_pip, v_dcafit_mean_pip,
      v_dcafit_sigm_pip, v_dcafit_emean_pip, v_dcafit_esigm_pip;
  std::vector<double> v_dcaxfit_x_pip, v_dcaxfit_ex_pip, v_dcaxfit_mean_pip,
      v_dcaxfit_sigm_pip, v_dcaxfit_emean_pip, v_dcaxfit_esigm_pip;
  std::vector<double> v_dcayfit_x_pip, v_dcayfit_ex_pip, v_dcayfit_mean_pip,
      v_dcayfit_sigm_pip, v_dcayfit_emean_pip, v_dcayfit_esigm_pip;
  std::vector<double> v_dcazfit_x_pip, v_dcazfit_ex_pip, v_dcazfit_mean_pip,
      v_dcazfit_sigm_pip, v_dcazfit_emean_pip, v_dcazfit_esigm_pip;
  std::vector<double> v_dcafit_x_pim, v_dcafit_ex_pim, v_dcafit_mean_pim,
      v_dcafit_sigm_pim, v_dcafit_emean_pim, v_dcafit_esigm_pim;
  std::vector<double> v_dcaxfit_x_pim, v_dcaxfit_ex_pim, v_dcaxfit_mean_pim,
      v_dcaxfit_sigm_pim, v_dcaxfit_emean_pim, v_dcaxfit_esigm_pim;
  std::vector<double> v_dcayfit_x_pim, v_dcayfit_ex_pim, v_dcayfit_mean_pim,
      v_dcayfit_sigm_pim, v_dcayfit_emean_pim, v_dcayfit_esigm_pim;
  std::vector<double> v_dcazfit_x_pim, v_dcazfit_ex_pim, v_dcazfit_mean_pim,
      v_dcazfit_sigm_pim, v_dcazfit_emean_pim, v_dcazfit_esigm_pim;
  
  int ptbins = vh2_recDcaPt.at(0)->GetNbinsX();
  for (int i = 1; i < ptbins + 1; i++) {
    vh1_dca_prt.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(0)->ProjectionX(Form("h1_dca_prt_ptbin%i", i), i, i))));
    vh1_dcax_prt.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(1)->ProjectionX(Form("h1_dcax_prt_ptbin%i", i), i, i))));
    vh1_dcay_prt.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(2)->ProjectionX(Form("h1_dcay_prt_ptbin%i", i), i, i))));
    vh1_dcaz_prt.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(3)->ProjectionX(Form("h1_dcaz_prt_ptbin%i", i), i, i))));
    vh1_dca_pip.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(4)->ProjectionX(Form("h1_dca_pip_ptbin%i", i), i, i))));
    vh1_dcax_pip.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(5)->ProjectionX(Form("h1_dcax_pip_ptbin%i", i), i, i))));
    vh1_dcay_pip.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(6)->ProjectionX(Form("h1_dcay_pip_ptbin%i", i), i, i))));
    vh1_dcaz_pip.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(7)->ProjectionX(Form("h1_dcaz_pip_ptbin%i", i), i, i))));
    vh1_dca_pim.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(8)->ProjectionX(Form("h1_dca_pim_ptbin%i", i), i, i))));
    vh1_dcax_pim.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(9)->ProjectionX(Form("h1_dcax_pim_ptbin%i", i), i, i))));
    vh1_dcay_pim.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(10)->ProjectionX(Form("h1_dcay_pim_ptbin%i", i), i, i))));
    vh1_dcaz_pim.emplace_back(std::make_unique<TH1D>(*(
        vh2_recDcaPt.at(11)->ProjectionX(Form("h1_dcaz_pim_ptbin%i", i), i, i))));
    
    vf1_dca_prt.emplace_back(
        std::make_unique<TF1>(Form("f1_dca_prt_ptbin%i", i), "gaus", 0., 1.));
    vf1_dcax_prt.emplace_back(
        std::make_unique<TF1>(Form("f1_dcax_prt_ptbin%i", i), "gaus", -0.2, 0.2));
    vf1_dcay_prt.emplace_back(
        std::make_unique<TF1>(Form("f1_dcay_prt_ptbin%i", i), "gaus", -0.2, 0.2));
    vf1_dcaz_prt.emplace_back(
        std::make_unique<TF1>(Form("f1_dcaz_prt_ptbin%i", i), "gaus", -0.2, 0.2));
    vf1_dca_pip.emplace_back(
        std::make_unique<TF1>(Form("f1_dca_pip_ptbin%i", i), "gaus", 0., 1.));
    vf1_dcax_pip.emplace_back(
        std::make_unique<TF1>(Form("f1_dcax_pip_ptbin%i", i), "gaus", -0.2, 0.2));
    vf1_dcay_pip.emplace_back(
        std::make_unique<TF1>(Form("f1_dcay_pip_ptbin%i", i), "gaus", -0.2, 0.2));
    vf1_dcaz_pip.emplace_back(
        std::make_unique<TF1>(Form("f1_dcaz_pip_ptbin%i", i), "gaus", -0.2, 0.2));
    vf1_dca_pim.emplace_back(
        std::make_unique<TF1>(Form("f1_dca_pim_ptbin%i", i), "gaus", 0., 1.));
    vf1_dcax_pim.emplace_back(
        std::make_unique<TF1>(Form("f1_dcax_pim_ptbin%i", i), "gaus", -0.2, 0.2));
    vf1_dcay_pim.emplace_back(
        std::make_unique<TF1>(Form("f1_dcay_pim_ptbin%i", i), "gaus", -0.2, 0.2));
    vf1_dcaz_pim.emplace_back(
        std::make_unique<TF1>(Form("f1_dcaz_pim_ptbin%i", i), "gaus", -0.2, 0.2));
    
    if (vh1_dca_prt.back()->GetEntries() > minentries)
      v_dcafit_x_prt.emplace_back(
          vh2_recDcaPt.at(0)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcax_prt.back()->GetEntries() > minentries)
      v_dcaxfit_x_prt.emplace_back(
          vh2_recDcaPt.at(1)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcay_prt.back()->GetEntries() > minentries)
      v_dcayfit_x_prt.emplace_back(
          vh2_recDcaPt.at(2)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcaz_prt.back()->GetEntries() > minentries)
      v_dcazfit_x_prt.emplace_back(
          vh2_recDcaPt.at(3)->GetYaxis()->GetBinCenter(i));
    if (vh1_dca_pip.back()->GetEntries() > minentries)
      v_dcafit_x_pip.emplace_back(
          vh2_recDcaPt.at(4)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcax_pip.back()->GetEntries() > minentries)
      v_dcaxfit_x_pip.emplace_back(
          vh2_recDcaPt.at(5)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcay_pip.back()->GetEntries() > minentries)
      v_dcayfit_x_pip.emplace_back(
          vh2_recDcaPt.at(6)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcaz_pip.back()->GetEntries() > minentries)
      v_dcazfit_x_pip.emplace_back(
          vh2_recDcaPt.at(7)->GetYaxis()->GetBinCenter(i));
    if (vh1_dca_pim.back()->GetEntries() > minentries)
      v_dcafit_x_pim.emplace_back(
          vh2_recDcaPt.at(8)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcax_pim.back()->GetEntries() > minentries)
      v_dcaxfit_x_pim.emplace_back(
          vh2_recDcaPt.at(9)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcay_pim.back()->GetEntries() > minentries)
      v_dcayfit_x_pim.emplace_back(
          vh2_recDcaPt.at(10)->GetYaxis()->GetBinCenter(i));
    if (vh1_dcaz_pim.back()->GetEntries() > minentries)
      v_dcazfit_x_pim.emplace_back(
          vh2_recDcaPt.at(11)->GetYaxis()->GetBinCenter(i));
    if (vh1_dca_prt.back()->GetEntries() > minentries)
      v_dcafit_ex_prt.emplace_back(
          vh2_recDcaPt.at(0)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcax_prt.back()->GetEntries() > minentries)
      v_dcaxfit_ex_prt.emplace_back(
          vh2_recDcaPt.at(1)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcay_prt.back()->GetEntries() > minentries)
      v_dcayfit_ex_prt.emplace_back(
          vh2_recDcaPt.at(2)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcaz_prt.back()->GetEntries() > minentries)
      v_dcazfit_ex_prt.emplace_back(
          vh2_recDcaPt.at(3)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dca_pip.back()->GetEntries() > minentries)
      v_dcafit_ex_pip.emplace_back(
          vh2_recDcaPt.at(4)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcax_pip.back()->GetEntries() > minentries)
      v_dcaxfit_ex_pip.emplace_back(
          vh2_recDcaPt.at(5)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcay_pip.back()->GetEntries() > minentries)
      v_dcayfit_ex_pip.emplace_back(
          vh2_recDcaPt.at(6)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcaz_pip.back()->GetEntries() > minentries)
      v_dcazfit_ex_pip.emplace_back(
          vh2_recDcaPt.at(7)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dca_pim.back()->GetEntries() > minentries)
      v_dcafit_ex_pim.emplace_back(
          vh2_recDcaPt.at(8)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcax_pim.back()->GetEntries() > minentries)
      v_dcaxfit_ex_pim.emplace_back(
          vh2_recDcaPt.at(9)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcay_pim.back()->GetEntries() > minentries)
      v_dcayfit_ex_pim.emplace_back(
          vh2_recDcaPt.at(10)->GetYaxis()->GetBinWidth(i) * 0.5);
    if (vh1_dcaz_pim.back()->GetEntries() > minentries)
      v_dcazfit_ex_pim.emplace_back(
          vh2_recDcaPt.at(11)->GetYaxis()->GetBinWidth(i) * 0.5);

  }

  // Fitting DCA distributions
  std::cout << "Fitting DCA distributions" << std::endl;
  for (int i = 0; i < vh1_dca_prt.size(); i++) {
    if (vh1_dca_prt.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dca_prt.at(i)->Fit(&*vf1_dca_prt.at(i), "RNQ");
    }
    v_dcafit_mean_prt.emplace_back(vf1_dca_prt.at(i)->GetParameter(1));
    v_dcafit_emean_prt.emplace_back(vf1_dca_prt.at(i)->GetParError(1));
    v_dcafit_sigm_prt.emplace_back(vf1_dca_prt.at(i)->GetParameter(2));
    v_dcafit_esigm_prt.emplace_back(vf1_dca_prt.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcax_prt.size(); i++) {
    if (vh1_dcax_prt.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcax_prt.at(i)->Fit(&*vf1_dcax_prt.at(i), "RNQ");
    }
    v_dcaxfit_mean_prt.emplace_back(vf1_dcax_prt.at(i)->GetParameter(1));
    v_dcaxfit_emean_prt.emplace_back(vf1_dcax_prt.at(i)->GetParError(1));
    v_dcaxfit_sigm_prt.emplace_back(vf1_dcax_prt.at(i)->GetParameter(2));
    v_dcaxfit_esigm_prt.emplace_back(vf1_dcax_prt.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcay_prt.size(); i++) {
    if (vh1_dcay_prt.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcay_prt.at(i)->Fit(&*vf1_dcay_prt.at(i), "RNQ");
    }
    v_dcayfit_mean_prt.emplace_back(vf1_dcay_prt.at(i)->GetParameter(1));
    v_dcayfit_emean_prt.emplace_back(vf1_dcay_prt.at(i)->GetParError(1));
    v_dcayfit_sigm_prt.emplace_back(vf1_dcay_prt.at(i)->GetParameter(2));
    v_dcayfit_esigm_prt.emplace_back(vf1_dcay_prt.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcaz_prt.size(); i++) {
    if (vh1_dcaz_prt.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcaz_prt.at(i)->Fit(&*vf1_dcaz_prt.at(i), "RNQ");
    }
    v_dcazfit_mean_prt.emplace_back(vf1_dcaz_prt.at(i)->GetParameter(1));
    v_dcazfit_emean_prt.emplace_back(vf1_dcaz_prt.at(i)->GetParError(1));
    v_dcazfit_sigm_prt.emplace_back(vf1_dcaz_prt.at(i)->GetParameter(2));
    v_dcazfit_esigm_prt.emplace_back(vf1_dcaz_prt.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dca_pip.size(); i++) {
    if (vh1_dca_pip.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dca_pip.at(i)->Fit(&*vf1_dca_pip.at(i), "RNQ");
    }
    v_dcafit_mean_pip.emplace_back(vf1_dca_pip.at(i)->GetParameter(1));
    v_dcafit_emean_pip.emplace_back(vf1_dca_pip.at(i)->GetParError(1));
    v_dcafit_sigm_pip.emplace_back(vf1_dca_pip.at(i)->GetParameter(2));
    v_dcafit_esigm_pip.emplace_back(vf1_dca_pip.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcax_pip.size(); i++) {
    if (vh1_dcax_pip.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcax_pip.at(i)->Fit(&*vf1_dcax_pip.at(i), "RNQ");
    }
    v_dcaxfit_mean_pip.emplace_back(vf1_dcax_pip.at(i)->GetParameter(1));
    v_dcaxfit_emean_pip.emplace_back(vf1_dcax_pip.at(i)->GetParError(1));
    v_dcaxfit_sigm_pip.emplace_back(vf1_dcax_pip.at(i)->GetParameter(2));
    v_dcaxfit_esigm_pip.emplace_back(vf1_dcax_pip.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcay_pip.size(); i++) {
    if (vh1_dcay_pip.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcay_pip.at(i)->Fit(&*vf1_dcay_pip.at(i), "RNQ");
    }
    v_dcayfit_mean_pip.emplace_back(vf1_dcay_pip.at(i)->GetParameter(1));
    v_dcayfit_emean_pip.emplace_back(vf1_dcay_pip.at(i)->GetParError(1));
    v_dcayfit_sigm_pip.emplace_back(vf1_dcay_pip.at(i)->GetParameter(2));
    v_dcayfit_esigm_pip.emplace_back(vf1_dcay_pip.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcaz_pip.size(); i++) {
    if (vh1_dcaz_pip.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcaz_pip.at(i)->Fit(&*vf1_dcaz_pip.at(i), "RNQ");
    }
    v_dcazfit_mean_pip.emplace_back(vf1_dcaz_pip.at(i)->GetParameter(1));
    v_dcazfit_emean_pip.emplace_back(vf1_dcaz_pip.at(i)->GetParError(1));
    v_dcazfit_sigm_pip.emplace_back(vf1_dcaz_pip.at(i)->GetParameter(2));
    v_dcazfit_esigm_pip.emplace_back(vf1_dcaz_pip.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dca_pim.size(); i++) {
    if (vh1_dca_pim.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dca_pim.at(i)->Fit(&*vf1_dca_pim.at(i), "RNQ");
    }
    v_dcafit_mean_pim.emplace_back(vf1_dca_pim.at(i)->GetParameter(1));
    v_dcafit_emean_pim.emplace_back(vf1_dca_pim.at(i)->GetParError(1));
    v_dcafit_sigm_pim.emplace_back(vf1_dca_pim.at(i)->GetParameter(2));
    v_dcafit_esigm_pim.emplace_back(vf1_dca_pim.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcax_pim.size(); i++) {
    if (vh1_dcax_pim.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcax_pim.at(i)->Fit(&*vf1_dcax_pim.at(i), "RNQ");
    }
    v_dcaxfit_mean_pim.emplace_back(vf1_dcax_pim.at(i)->GetParameter(1));
    v_dcaxfit_emean_pim.emplace_back(vf1_dcax_pim.at(i)->GetParError(1));
    v_dcaxfit_sigm_pim.emplace_back(vf1_dcax_pim.at(i)->GetParameter(2));
    v_dcaxfit_esigm_pim.emplace_back(vf1_dcax_pim.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcay_pim.size(); i++) {
    if (vh1_dcay_pim.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcay_pim.at(i)->Fit(&*vf1_dcay_pim.at(i), "RNQ");
    }
    v_dcayfit_mean_pim.emplace_back(vf1_dcay_pim.at(i)->GetParameter(1));
    v_dcayfit_emean_pim.emplace_back(vf1_dcay_pim.at(i)->GetParError(1));
    v_dcayfit_sigm_pim.emplace_back(vf1_dcay_pim.at(i)->GetParameter(2));
    v_dcayfit_esigm_pim.emplace_back(vf1_dcay_pim.at(i)->GetParError(2));
  }
  for (int i = 0; i < vh1_dcaz_pim.size(); i++) {
    if (vh1_dcaz_pim.at(i)->GetEntries() <= minentries)
      continue;
    for (int k = 0; k < niter; ++k) {
      vh1_dcaz_pim.at(i)->Fit(&*vf1_dcaz_pim.at(i), "RNQ");
    }
    v_dcazfit_mean_pim.emplace_back(vf1_dcaz_pim.at(i)->GetParameter(1));
    v_dcazfit_emean_pim.emplace_back(vf1_dcaz_pim.at(i)->GetParError(1));
    v_dcazfit_sigm_pim.emplace_back(vf1_dcaz_pim.at(i)->GetParameter(2));
    v_dcazfit_esigm_pim.emplace_back(vf1_dcaz_pim.at(i)->GetParError(2));
  }

  auto gr_dcafit_mean_prt = std::make_unique<TGraphErrors>(
      v_dcafit_x_prt.size(), v_dcafit_x_prt.data(), v_dcafit_mean_prt.data(),
      v_dcafit_ex_prt.data(), v_dcafit_emean_prt.data());
  gr_dcafit_mean_prt->SetName("gr_dcafit_mean_prt");
  gr_dcafit_mean_prt->SetTitle("mean of DCA peak in p_{T} bins for protons");

  auto gr_dcaxfit_mean_prt = std::make_unique<TGraphErrors>(
      v_dcaxfit_x_prt.size(), v_dcaxfit_x_prt.data(), v_dcaxfit_mean_prt.data(),
      v_dcaxfit_ex_prt.data(), v_dcaxfit_emean_prt.data());
  gr_dcaxfit_mean_prt->SetName("gr_dcaxfit_mean_prt");
  gr_dcaxfit_mean_prt->SetTitle("mean of DCA_{x} peak in p_{T} bins for protons");

  auto gr_dcayfit_mean_prt = std::make_unique<TGraphErrors>(
      v_dcayfit_x_prt.size(), v_dcayfit_x_prt.data(), v_dcayfit_mean_prt.data(),
      v_dcayfit_ex_prt.data(), v_dcayfit_emean_prt.data());
  gr_dcayfit_mean_prt->SetName("gr_dcayfit_mean_prt");
  gr_dcayfit_mean_prt->SetTitle("mean of DCA_{y} peak in p_{T} bins for protons");

  auto gr_dcazfit_mean_prt = std::make_unique<TGraphErrors>(
      v_dcazfit_x_prt.size(), v_dcazfit_x_prt.data(), v_dcazfit_mean_prt.data(),
      v_dcazfit_ex_prt.data(), v_dcazfit_emean_prt.data());
  gr_dcazfit_mean_prt->SetName("gr_dcazfit_mean_prt");
  gr_dcazfit_mean_prt->SetTitle("mean of DCA_{z} peak in p_{T} bins for protons");

  auto gr_dcafit_mean_pip = std::make_unique<TGraphErrors>(
      v_dcafit_x_pip.size(), v_dcafit_x_pip.data(), v_dcafit_mean_pip.data(),
      v_dcafit_ex_pip.data(), v_dcafit_emean_pip.data());
  gr_dcafit_mean_pip->SetName("gr_dcafit_mean_pip");
  gr_dcafit_mean_pip->SetTitle(
      "mean of DCA peak in p_{T} bins for pions (#pi^{+})");

  auto gr_dcaxfit_mean_pip = std::make_unique<TGraphErrors>(
      v_dcaxfit_x_pip.size(), v_dcaxfit_x_pip.data(), v_dcaxfit_mean_pip.data(),
      v_dcaxfit_ex_pip.data(), v_dcaxfit_emean_pip.data());
  gr_dcaxfit_mean_pip->SetName("gr_dcaxfit_mean_pip");
  gr_dcaxfit_mean_pip->SetTitle("mean of DCA_{x} peak in p_{T} bins for pions (#pi^{+})");

  auto gr_dcayfit_mean_pip = std::make_unique<TGraphErrors>(
      v_dcayfit_x_pip.size(), v_dcayfit_x_pip.data(), v_dcayfit_mean_pip.data(),
      v_dcayfit_ex_pip.data(), v_dcayfit_emean_pip.data());
  gr_dcayfit_mean_pip->SetName("gr_dcayfit_mean_pip");
  gr_dcayfit_mean_pip->SetTitle("mean of DCA_{y} peak in p_{T} bins for pions (#pi^{+})");

  auto gr_dcazfit_mean_pip = std::make_unique<TGraphErrors>(
      v_dcazfit_x_pip.size(), v_dcazfit_x_pip.data(), v_dcazfit_mean_pip.data(),
      v_dcazfit_ex_pip.data(), v_dcazfit_emean_pip.data());
  gr_dcazfit_mean_pip->SetName("gr_dcazfit_mean_pip");
  gr_dcazfit_mean_pip->SetTitle("mean of DCA_{z} peak in p_{T} bins for pions (#pi^{+})");

  auto gr_dcafit_mean_pim = std::make_unique<TGraphErrors>(
      v_dcafit_x_pim.size(), v_dcafit_x_pim.data(), v_dcafit_mean_pim.data(),
      v_dcafit_ex_pim.data(), v_dcafit_emean_pim.data());
  gr_dcafit_mean_pim->SetName("gr_dcafit_mean_pim");
  gr_dcafit_mean_pim->SetTitle(
      "mean of DCA peak in p_{T} bins for pions (#pi^{-})");

  auto gr_dcaxfit_mean_pim = std::make_unique<TGraphErrors>(
      v_dcaxfit_x_pim.size(), v_dcaxfit_x_pim.data(), v_dcaxfit_mean_pim.data(),
      v_dcaxfit_ex_pim.data(), v_dcaxfit_emean_pim.data());
  gr_dcaxfit_mean_pim->SetName("gr_dcaxfit_mean_pim");
  gr_dcaxfit_mean_pim->SetTitle("mean of DCA_{x} peak in p_{T} bins for pions (#pi^{-})");

  auto gr_dcayfit_mean_pim = std::make_unique<TGraphErrors>(
      v_dcayfit_x_pim.size(), v_dcayfit_x_pim.data(), v_dcayfit_mean_pim.data(),
      v_dcayfit_ex_pim.data(), v_dcayfit_emean_pim.data());
  gr_dcayfit_mean_pim->SetName("gr_dcayfit_mean_pim");
  gr_dcayfit_mean_pim->SetTitle("mean of DCA_{y} peak in p_{T} bins for pions (#pi^{-})");

  auto gr_dcazfit_mean_pim = std::make_unique<TGraphErrors>(
      v_dcazfit_x_pim.size(), v_dcazfit_x_pim.data(), v_dcazfit_mean_pim.data(),
      v_dcazfit_ex_pim.data(), v_dcazfit_emean_pim.data());
  gr_dcazfit_mean_pim->SetName("gr_dcazfit_mean_pim");
  gr_dcazfit_mean_pim->SetTitle("mean of DCA_{z} peak in p_{T} bins for pions (#pi^{-})");

  auto gr_dcafit_sigm_prt = std::make_unique<TGraphErrors>(
      v_dcafit_x_prt.size(), v_dcafit_x_prt.data(), v_dcafit_sigm_prt.data(),
      v_dcafit_ex_prt.data(), v_dcafit_esigm_prt.data());
  gr_dcafit_sigm_prt->SetName("gr_dcafit_sigm_prt");
  gr_dcafit_sigm_prt->SetTitle("sigma of DCA peak in p_{T} bins for protons");

  auto gr_dcaxfit_sigm_prt = std::make_unique<TGraphErrors>(
      v_dcaxfit_x_prt.size(), v_dcaxfit_x_prt.data(), v_dcaxfit_sigm_prt.data(),
      v_dcaxfit_ex_prt.data(), v_dcaxfit_esigm_prt.data());
  gr_dcaxfit_sigm_prt->SetName("gr_dcaxfit_sigm_prt");
  gr_dcaxfit_sigm_prt->SetTitle("sigma of DCA_{x} peak in p_{T} bins for protons");

  auto gr_dcayfit_sigm_prt = std::make_unique<TGraphErrors>(
      v_dcayfit_x_prt.size(), v_dcayfit_x_prt.data(), v_dcayfit_sigm_prt.data(),
      v_dcayfit_ex_prt.data(), v_dcayfit_esigm_prt.data());
  gr_dcayfit_sigm_prt->SetName("gr_dcayfit_sigm_prt");
  gr_dcayfit_sigm_prt->SetTitle("sigma of DCA_{y} peak in p_{T} bins for protons");

  auto gr_dcazfit_sigm_prt = std::make_unique<TGraphErrors>(
      v_dcazfit_x_prt.size(), v_dcazfit_x_prt.data(), v_dcazfit_sigm_prt.data(),
      v_dcazfit_ex_prt.data(), v_dcazfit_esigm_prt.data());
  gr_dcazfit_sigm_prt->SetName("gr_dcazfit_sigm_prt");
  gr_dcazfit_sigm_prt->SetTitle("sigma of DCA_{z} peak in p_{T} bins for protons");

  auto gr_dcafit_sigm_pip = std::make_unique<TGraphErrors>(
      v_dcafit_x_pip.size(), v_dcafit_x_pip.data(), v_dcafit_sigm_pip.data(),
      v_dcafit_ex_pip.data(), v_dcafit_esigm_pip.data());
  gr_dcafit_sigm_pip->SetName("gr_dcafit_sigm_pip");
  gr_dcafit_sigm_pip->SetTitle(
      "sigma of DCA peak in p_{T} bins for pions (#pi^{+})");

  auto gr_dcaxfit_sigm_pip = std::make_unique<TGraphErrors>(
      v_dcaxfit_x_pip.size(), v_dcaxfit_x_pip.data(), v_dcaxfit_sigm_pip.data(),
      v_dcaxfit_ex_pip.data(), v_dcaxfit_esigm_pip.data());
  gr_dcaxfit_sigm_pip->SetName("gr_dcaxfit_sigm_pip");
  gr_dcaxfit_sigm_pip->SetTitle("sigma of DCA_{x} peak in p_{T} bins for pions (#pi^{+})");

  auto gr_dcayfit_sigm_pip = std::make_unique<TGraphErrors>(
      v_dcayfit_x_pip.size(), v_dcayfit_x_pip.data(), v_dcayfit_sigm_pip.data(),
      v_dcayfit_ex_pip.data(), v_dcayfit_esigm_pip.data());
  gr_dcayfit_sigm_pip->SetName("gr_dcayfit_sigm_pip");
  gr_dcayfit_sigm_pip->SetTitle("sigma of DCA_{y} peak in p_{T} bins for pions (#pi^{+})");

  auto gr_dcazfit_sigm_pip = std::make_unique<TGraphErrors>(
      v_dcazfit_x_pip.size(), v_dcazfit_x_pip.data(), v_dcazfit_sigm_pip.data(),
      v_dcazfit_ex_pip.data(), v_dcazfit_esigm_pip.data());
  gr_dcazfit_sigm_pip->SetName("gr_dcazfit_sigm_pip");
  gr_dcazfit_sigm_pip->SetTitle("sigma of DCA_{z} peak in p_{T} bins for pions (#pi^{+})");

  auto gr_dcafit_sigm_pim = std::make_unique<TGraphErrors>(
      v_dcafit_x_pim.size(), v_dcafit_x_pim.data(), v_dcafit_sigm_pim.data(),
      v_dcafit_ex_pim.data(), v_dcafit_esigm_pim.data());
  gr_dcafit_sigm_pim->SetName("gr_dcafit_sigm_pim");
  gr_dcafit_sigm_pim->SetTitle(
      "sigma of DCA peak in p_{T} bins for pions (#pi^{-})");

  auto gr_dcaxfit_sigm_pim = std::make_unique<TGraphErrors>(
      v_dcaxfit_x_pim.size(), v_dcaxfit_x_pim.data(), v_dcaxfit_sigm_pim.data(),
      v_dcaxfit_ex_pim.data(), v_dcaxfit_esigm_pim.data());
  gr_dcaxfit_sigm_pim->SetName("gr_dcaxfit_sigm_pim");
  gr_dcaxfit_sigm_pim->SetTitle("sigma of DCA_{x} peak in p_{T} bins for pions (#pi^{-})");

  auto gr_dcayfit_sigm_pim = std::make_unique<TGraphErrors>(
      v_dcayfit_x_pim.size(), v_dcayfit_x_pim.data(), v_dcayfit_sigm_pim.data(),
      v_dcayfit_ex_pim.data(), v_dcayfit_esigm_pim.data());
  gr_dcayfit_sigm_pim->SetName("gr_dcayfit_sigm_pim");
  gr_dcayfit_sigm_pim->SetTitle("sigma of DCA_{y} peak in p_{T} bins for pions (#pi^{-})");

  auto gr_dcazfit_sigm_pim = std::make_unique<TGraphErrors>(
      v_dcazfit_x_pim.size(), v_dcazfit_x_pim.data(), v_dcazfit_sigm_pim.data(),
      v_dcazfit_ex_pim.data(), v_dcazfit_esigm_pim.data());
  gr_dcazfit_sigm_pim->SetName("gr_dcazfit_sigm_pim");
  gr_dcazfit_sigm_pim->SetTitle("sigma of DCA_{z} peak in p_{T} bins for pions (#pi^{-})");
  
  // Writing output
  std::cout << "Writing output" << std::endl;
  auto fo = new TFile(outFile.c_str(), "recreate");
  fo->cd();

  gr_dcafit_mean_prt->Write();
  gr_dcafit_mean_pip->Write();
  gr_dcafit_mean_pim->Write();
  gr_dcafit_sigm_prt->Write();
  gr_dcafit_sigm_pip->Write();
  gr_dcafit_sigm_pim->Write();

  gr_dcaxfit_mean_prt->Write();
  gr_dcaxfit_mean_pip->Write();
  gr_dcaxfit_mean_pim->Write();
  gr_dcaxfit_sigm_prt->Write();
  gr_dcaxfit_sigm_pip->Write();
  gr_dcaxfit_sigm_pim->Write();

  gr_dcayfit_mean_prt->Write();
  gr_dcayfit_mean_pip->Write();
  gr_dcayfit_mean_pim->Write();
  gr_dcayfit_sigm_prt->Write();
  gr_dcayfit_sigm_pip->Write();
  gr_dcayfit_sigm_pim->Write();

  gr_dcazfit_mean_prt->Write();
  gr_dcazfit_mean_pip->Write();
  gr_dcazfit_mean_pim->Write();
  gr_dcazfit_sigm_prt->Write();
  gr_dcazfit_sigm_pip->Write();
  gr_dcazfit_sigm_pim->Write();

  fo->mkdir("PtSlicesPrt");
  fo->cd("PtSlicesPrt");
  for (int i = 0; i < vh1_dca_prt.size(); i++) {
    vh1_dca_prt.at(i)->Write();
    vf1_dca_prt.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcax_prt.size(); i++) {
    vh1_dcax_prt.at(i)->Write();
    vf1_dcax_prt.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcay_prt.size(); i++) {
    vh1_dcay_prt.at(i)->Write();
    vf1_dcay_prt.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcaz_prt.size(); i++) {
    vh1_dcaz_prt.at(i)->Write();
    vf1_dcaz_prt.at(i)->Write();
  }

  fo->mkdir("PtSlicesPip");
  fo->cd("PtSlicesPip");
  for (int i = 0; i < vh1_dca_pip.size(); i++) {
    vh1_dca_pip.at(i)->Write();
    vf1_dca_pip.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcax_pip.size(); i++) {
    vh1_dcax_pip.at(i)->Write();
    vf1_dcax_pip.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcay_pip.size(); i++) {
    vh1_dcay_pip.at(i)->Write();
    vf1_dcay_pip.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcaz_pip.size(); i++) {
    vh1_dcaz_pip.at(i)->Write();
    vf1_dcaz_pip.at(i)->Write();
  }

  fo->mkdir("PtSlicesPim");
  fo->cd("PtSlicesPim");
  for (int i = 0; i < vh1_dca_pim.size(); i++) {
    vh1_dca_pim.at(i)->Write();
    vf1_dca_pim.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcax_pim.size(); i++) {
    vh1_dcax_pim.at(i)->Write();
    vf1_dcax_pim.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcay_pim.size(); i++) {
    vh1_dcay_pim.at(i)->Write();
    vf1_dcay_pim.at(i)->Write();
  }
  for (int i = 0; i < vh1_dcaz_pim.size(); i++) {
    vh1_dcaz_pim.at(i)->Write();
    vf1_dcaz_pim.at(i)->Write();
  }

  fo->Close();

  timer.Stop();
  timer.Print();
}
