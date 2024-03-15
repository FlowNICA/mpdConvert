void getPID(std::string inFile_qa_step1, std::string outFile_pid)
{
  std::cout << "Starting PID" << std::endl;
  const int niter = 5;

  std::vector<double> m2_pdg = {
    pow(TDatabasePDG::Instance()->GetParticle(211)->Mass(), 2), // pions
    pow(TDatabasePDG::Instance()->GetParticle(321)->Mass(), 2), // kaons
    pow(TDatabasePDG::Instance()->GetParticle(2212)->Mass(), 2), // protons
  };

  std::cout << "Reading input" << std::endl;
  TFile *fi = new TFile(inFile_qa_step1.c_str(),"read");
  if (!fi){
    std::cerr << "Error: cannot open input file: " << inFile_qa_step1.c_str() << std::endl;
    return;
  }

  auto h2_recDedxPq = (TH2D*)fi->Get("h2_recDedxPq");
  auto h2_recM2Pq = (TH2D*)fi->Get("h2_recM2Pq");
  auto h3_recDedxM2Pq = (TH3D*)fi->Get("h3_recDedxM2Pq");

  if (!h2_recDedxPq || !h2_recM2Pq || !h3_recDedxM2Pq){
    std::cerr << "Error: cannot open needed histograms!" << std::endl;
    return;
  }

  //Prepare dE/dx profiles for PID fit
  h3_recDedxM2Pq->GetXaxis()->SetRangeUser(0.6, 1.2);
  h3_recDedxM2Pq->GetZaxis()->SetRangeUser(0.3, 5.);
  auto h2_recDedxPq_proton = (TH2D*)h3_recDedxM2Pq->Project3D("yz")->Clone();
  auto p_recDedxPq_proton = (TProfile*)h2_recDedxPq_proton->ProfileX(Form("p_%s_proton", h2_recDedxPq_proton->GetName()));

  //Prepare parametrization for the Bethe-Bloch distribution (from ALEPH, ALICE experiments)
  auto f1_dedx_proton = new TF1("f1_dedx_proton", "[0]/pow(x*x/(0.938*0.938+x*x),[3]*0.5)*( [1] + pow(x*x/(0.938*0.938+x*x),[3]*0.5)-TMath::Log([2]+pow(0.938/x,[4])))", 0.3, 5.);

  //Fit dE/dx vs P/q distribution
  for (int i=0; i<niter; i++)
    h2_recDedxPq_proton->Fit(f1_dedx_proton,"RNQ");

  auto h2_recDedxNormPq = std::make_unique<TH2D>("h2_recDedxNormPq", Form("Normalized %s", h2_recDedxPq->GetTitle()),
                                                 h2_recDedxPq->GetNbinsX(), h2_recDedxPq->GetXaxis()->GetBinLowEdge(1),
                                                 h2_recDedxPq->GetXaxis()->GetBinLowEdge(h2_recDedxPq->GetNbinsX()+1),
                                                 200, -2., 2.);

  auto h2_recDedxNormPq_proton = std::make_unique<TH2D>("h2_recDedxNormPq_proton", Form("Normalized %s for protons", h2_recDedxPq->GetTitle()),
                                                        h2_recDedxPq_proton->GetNbinsX(), h2_recDedxPq_proton->GetXaxis()->GetBinLowEdge(1),
                                                        h2_recDedxPq_proton->GetXaxis()->GetBinLowEdge(h2_recDedxPq_proton->GetNbinsX()+1),
                                                        200, -2., 2.);

  //Resampling of dE/dx vs P/q to get ((dE/dx)meas - (dE/dx)calc)/(dE/dx)calc vs P/q distribution
  std::cout << "Resampling dE/dx vs P/q distributions" << std::endl;
  long nentries = (long) h2_recDedxPq->GetEntries();
  double x_samp, y_samp;
  for (long i=0; i<nentries; ++i){
    h2_recDedxPq->GetRandom2(x_samp, y_samp);
    h2_recDedxNormPq->Fill(x_samp, (y_samp - f1_dedx_proton->Eval(x_samp))/f1_dedx_proton->Eval(x_samp));
  }
  nentries = (long) h2_recDedxPq_proton->GetEntries();
  for (long i=0; i<nentries; ++i){ 
    h2_recDedxPq_proton->GetRandom2(x_samp, y_samp);
    h2_recDedxNormPq_proton->Fill(x_samp, (y_samp - f1_dedx_proton->Eval(x_samp))/f1_dedx_proton->Eval(x_samp));
  }

  //Prepare for normalized dE/dx fit
  std::vector<std::unique_ptr<TH1D>> vh1_dedx_p; // vector of TH1D with dedx distribution for protons
  std::vector<std::unique_ptr<TF1>> vf1_dedx_1gaus_p; // vector of TF1 1gaus for proton dedx peaks
  std::vector<std::unique_ptr<TF1>> vf1_dedx_1gaus_pip; // vector of TF1 1gaus for pi+ dedx peaks
  std::vector<std::unique_ptr<TF1>> vf1_dedx_1gaus_kap; // vector of TF1 1gaus for K+ dedx peaks
  std::vector<double> vsigm_dedx_p_x, vsigm_dedx_p_y, vsigm_dedx_p_ex, vsigm_dedx_p_ey;

  auto h2_dedxpq_proton = std::make_unique<TH2D>(*((TH2D*)h2_recDedxNormPq_proton->Clone()));
  h2_dedxpq_proton->RebinX(5);

  int firstbin = h2_dedxpq_proton->GetXaxis()->FindBin(0.31);
  for (int i=firstbin; i<h2_dedxpq_proton->GetNbinsX(); ++i) {
    vh1_dedx_p.push_back( std::make_unique<TH1D>( *(h2_dedxpq_proton->ProjectionY(Form("h1_%s_p_pqbin%i", h2_dedxpq_proton->GetName(), i-firstbin), i, i)) ) );
    vh1_dedx_p.at(i-firstbin)->SetTitle( Form("%s for protons for %2.2f < p/q < %2.2f GeV/c",
                                          h2_dedxpq_proton->GetTitle(),
                                          h2_dedxpq_proton->GetXaxis()->GetBinLowEdge(i),
                                          h2_dedxpq_proton->GetXaxis()->GetBinUpEdge(i)) );

    vf1_dedx_1gaus_p.push_back( std::make_unique<TF1>(Form("f1_dedx_1gaus_p_pqbin%i", i-firstbin), "gaus", -0.5, 0.5) );
    vf1_dedx_1gaus_pip.push_back( std::make_unique<TF1>(Form("f1_dedx_1gaus_pip_pqbin%i", i-firstbin), "gaus", -0.5, 0.5) );
    vf1_dedx_1gaus_kap.push_back( std::make_unique<TF1>(Form("f1_dedx_1gaus_kap_pqbin%i", i-firstbin), "gaus", -0.5, 0.5) );
    vf1_dedx_1gaus_p.at(i-firstbin)->SetParameter(1, 0.);
    vf1_dedx_1gaus_pip.at(i-firstbin)->SetParameter(1, 0.);
    vf1_dedx_1gaus_kap.at(i-firstbin)->SetParameter(1, 0.);

    vsigm_dedx_p_x.push_back(h2_dedxpq_proton->GetXaxis()->GetBinCenter(i));
    vsigm_dedx_p_ex.push_back(h2_dedxpq_proton->GetXaxis()->GetBinWidth(i)*0.5);
  }

  // Fitting dE/dx distributions
  std::cout << "Fitting dE/dx distributions" << std::endl;
  for (int i=0; i<vh1_dedx_p.size(); ++i){
    std::cout << "\tfitting bin " << i << " out of " << vh1_dedx_p.size() << std::endl;
    for (int k=0; k<niter; ++k){
      vh1_dedx_p.at(i)->Fit(&*vf1_dedx_1gaus_p.at(i), "RNQ");
    }

    vsigm_dedx_p_y.push_back(vf1_dedx_1gaus_p.at(i)->GetParameter(2));
    vsigm_dedx_p_ey.push_back(vf1_dedx_1gaus_p.at(i)->GetParError(2));
  }

  auto gr_dedx_sigm_p = std::make_unique<TGraphErrors>(vsigm_dedx_p_x.size(),
                                                       &vsigm_dedx_p_x[0], &vsigm_dedx_p_y[0],
                                                       &vsigm_dedx_p_ex[0], &vsigm_dedx_p_ey[0]);

  gr_dedx_sigm_p->SetName("gr_dedx_sigm_p");
  gr_dedx_sigm_p->SetTitle("#sigma_{dE/dx} of protons vs P/q;p/q, GeV/c;#sigma_{dE/dx}, a.u.");

  auto f1_dedx_sigm_p = std::make_unique<TF1>("f1_dedx_sigm_p", "pol0", 0.3, 2.);
  gr_dedx_sigm_p->Fit(&*f1_dedx_sigm_p, "RNQ");

  //Prepare for M^2 fit
  std::vector<std::unique_ptr<TH1D>> vh1_m2; // vector of TH1D with m2 distribution
  std::vector<std::unique_ptr<TF1>> vf1_m2_1gaus_p; // vector of TF1 1gaus for proton m2 peaks
  std::vector<std::unique_ptr<TF1>> vf1_m2_1gaus_pip; // vector of TF1 1gaus for pi+ m2 peaks
  std::vector<std::unique_ptr<TF1>> vf1_m2_1gaus_kap; // vector of TF1 1gaus for K+ m2 peaks
  std::vector<std::unique_ptr<TF1>> vf1_m2_2gaus_p; // vector of TF1 2gaus for proton m2 peaks
  std::vector<std::unique_ptr<TF1>> vf1_m2_2gaus_pip; // vector of TF1 2gaus for pi+ m2 peaks
  std::vector<std::unique_ptr<TF1>> vf1_m2_2gaus_kap; // vector of TF1 2gaus for K+ m2 peaks
  std::vector<std::unique_ptr<TF1>> vf1_m2_2gaus_all; // vector of TF1 2gaus for all 3 m2 peaks
  std::vector<double> vsigm_m2_p_x, vsigm_m2_p_y, vsigm_m2_p_ex, vsigm_m2_p_ey;
  
  auto h2_m2pq = (TH2D*)h2_recM2Pq->Clone();
  h2_m2pq->RebinX(10);

  firstbin = h2_m2pq->GetXaxis()->FindBin(0.3);
  for (int i=firstbin; i<h2_m2pq->GetNbinsX(); ++i) {
    vh1_m2.push_back( std::make_unique<TH1D>( *(h2_m2pq->ProjectionY(Form("h1_%s_pqbin%i", h2_m2pq->GetName(), i-firstbin), i, i)) ) );
    vh1_m2.at(i-firstbin)->SetTitle( Form("%s for %2.2f < p/q < %2.2f GeV/c",
                                          h2_m2pq->GetTitle(),
                                          h2_m2pq->GetXaxis()->GetBinLowEdge(i),
                                          h2_m2pq->GetXaxis()->GetBinUpEdge(i)) );

    vf1_m2_1gaus_p.push_back( std::make_unique<TF1>(Form("f1_m2_1gaus_p_pqbin%i", i-firstbin), "gaus", m2_pdg.at(2)-0.1, m2_pdg.at(2)+0.1) );
    vf1_m2_1gaus_pip.push_back( std::make_unique<TF1>(Form("f1_m2_1gaus_pip_pqbin%i", i-firstbin), "gaus", m2_pdg.at(0)-0.1, m2_pdg.at(0)+0.1) );
    vf1_m2_1gaus_kap.push_back( std::make_unique<TF1>(Form("f1_m2_1gaus_kap_pqbin%i", i-firstbin), "gaus", m2_pdg.at(1)-0.1, m2_pdg.at(1)+0.1) );
    vf1_m2_1gaus_p.at(i-firstbin)->SetParameter(1, m2_pdg.at(2));
    vf1_m2_1gaus_pip.at(i-firstbin)->SetParameter(1, m2_pdg.at(1));
    vf1_m2_1gaus_kap.at(i-firstbin)->SetParameter(1, m2_pdg.at(0));

    vf1_m2_2gaus_p.push_back( std::make_unique<TF1>(Form("f1_m2_2gaus_p_pqbin%i", i-firstbin), "gaus(0)+gaus(3)", m2_pdg.at(2)-0.3, m2_pdg.at(2)+0.3) );
    vf1_m2_2gaus_pip.push_back( std::make_unique<TF1>(Form("f1_m2_2gaus_pip_pqbin%i", i-firstbin), "gaus(0)+gaus(3)", m2_pdg.at(0)-0.2, m2_pdg.at(0)+0.2) );
    vf1_m2_2gaus_kap.push_back( std::make_unique<TF1>(Form("f1_m2_2gaus_kap_pqbin%i", i-firstbin), "gaus(0)+gaus(3)", m2_pdg.at(1)-0.1, m2_pdg.at(1)+0.1) );
    vf1_m2_2gaus_p.at(i-firstbin)->SetParameter(1, m2_pdg.at(2));
    vf1_m2_2gaus_p.at(i-firstbin)->SetParameter(4, m2_pdg.at(2));
    vf1_m2_2gaus_pip.at(i-firstbin)->SetParameter(1, m2_pdg.at(0));
    vf1_m2_2gaus_pip.at(i-firstbin)->SetParameter(4, m2_pdg.at(0));
    vf1_m2_2gaus_kap.at(i-firstbin)->SetParameter(1, m2_pdg.at(1));
    vf1_m2_2gaus_kap.at(i-firstbin)->SetParameter(4, m2_pdg.at(1));

    vf1_m2_2gaus_all.push_back( std::make_unique<TF1>(Form("f1_m2_2gaus_all_pqbin%i", i-firstbin),
                                                      "gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",
                                                      m2_pdg.at(0)-0.2, m2_pdg.at(2)*1.2) );
    vf1_m2_2gaus_all.at(i-firstbin)->SetParameter(1, m2_pdg.at(0));
    vf1_m2_2gaus_all.at(i-firstbin)->SetParameter(4, m2_pdg.at(0));
    vf1_m2_2gaus_all.at(i-firstbin)->SetParameter(7, m2_pdg.at(1));
    vf1_m2_2gaus_all.at(i-firstbin)->SetParameter(10, m2_pdg.at(1));
    vf1_m2_2gaus_all.at(i-firstbin)->SetParameter(13, m2_pdg.at(2));
    vf1_m2_2gaus_all.at(i-firstbin)->SetParameter(16, m2_pdg.at(2));

    vsigm_m2_p_x.push_back(h2_m2pq->GetXaxis()->GetBinCenter(i));
    vsigm_m2_p_ex.push_back(h2_m2pq->GetXaxis()->GetBinWidth(i)*0.5);
  }
  
  // Fitting m2 distributions
  std::cout << "Fitting m2 distributions" << std::endl;
  for (int i=0; i<vh1_m2.size(); ++i){
    std::cout << "\tfitting bin " << i << " out of " << vh1_m2.size() << std::endl;
    for (int k=0; k<niter; ++k){
      vh1_m2.at(i)->Fit(&*vf1_m2_1gaus_p.at(i), "RNQ");
      vh1_m2.at(i)->Fit(&*vf1_m2_1gaus_pip.at(i), "RNQ");
      vh1_m2.at(i)->Fit(&*vf1_m2_1gaus_kap.at(i), "RNQ");
    }
    vf1_m2_2gaus_p.at(i)->SetParameter(0, vf1_m2_1gaus_p.at(i)->GetParameter(0));
    vf1_m2_2gaus_p.at(i)->SetParameter(3, vf1_m2_1gaus_p.at(i)->GetParameter(0));
    vf1_m2_2gaus_p.at(i)->SetParameter(2, vf1_m2_1gaus_p.at(i)->GetParameter(2));
    vf1_m2_2gaus_p.at(i)->SetParameter(5, 5.*vf1_m2_1gaus_p.at(i)->GetParameter(2));

    vf1_m2_2gaus_pip.at(i)->SetParameter(0, vf1_m2_1gaus_pip.at(i)->GetParameter(0));
    vf1_m2_2gaus_pip.at(i)->SetParameter(3, vf1_m2_1gaus_pip.at(i)->GetParameter(0));
    vf1_m2_2gaus_pip.at(i)->SetParameter(2, vf1_m2_1gaus_pip.at(i)->GetParameter(2));
    vf1_m2_2gaus_pip.at(i)->SetParameter(5, 5.*vf1_m2_1gaus_pip.at(i)->GetParameter(2));

    vf1_m2_2gaus_kap.at(i)->SetParameter(0, vf1_m2_1gaus_kap.at(i)->GetParameter(0));
    vf1_m2_2gaus_kap.at(i)->SetParameter(3, vf1_m2_1gaus_kap.at(i)->GetParameter(0));
    vf1_m2_2gaus_kap.at(i)->SetParameter(2, vf1_m2_1gaus_kap.at(i)->GetParameter(2));
    vf1_m2_2gaus_kap.at(i)->SetParameter(5, 5.*vf1_m2_1gaus_kap.at(i)->GetParameter(2));

    for (int k=0; k<niter; ++k){
      vh1_m2.at(i)->Fit(&*vf1_m2_2gaus_p.at(i), "RNQ");
      vh1_m2.at(i)->Fit(&*vf1_m2_2gaus_pip.at(i), "RNQ");
      vh1_m2.at(i)->Fit(&*vf1_m2_2gaus_kap.at(i), "RNQ");
    }

    vf1_m2_2gaus_all.at(i)->SetParameter(0, vf1_m2_2gaus_pip.at(i)->GetParameter(0));
    vf1_m2_2gaus_all.at(i)->SetParameter(2, vf1_m2_2gaus_pip.at(i)->GetParameter(2));
    vf1_m2_2gaus_all.at(i)->SetParameter(3, vf1_m2_2gaus_pip.at(i)->GetParameter(3));
    vf1_m2_2gaus_all.at(i)->SetParameter(5, vf1_m2_2gaus_pip.at(i)->GetParameter(5));
    vf1_m2_2gaus_all.at(i)->SetParameter(6, vf1_m2_2gaus_kap.at(i)->GetParameter(0));
    vf1_m2_2gaus_all.at(i)->SetParameter(8, vf1_m2_2gaus_kap.at(i)->GetParameter(2));
    vf1_m2_2gaus_all.at(i)->SetParameter(9, vf1_m2_2gaus_kap.at(i)->GetParameter(3));
    vf1_m2_2gaus_all.at(i)->SetParameter(11, vf1_m2_2gaus_kap.at(i)->GetParameter(5));
    vf1_m2_2gaus_all.at(i)->SetParameter(12, vf1_m2_2gaus_p.at(i)->GetParameter(0));
    vf1_m2_2gaus_all.at(i)->SetParameter(14, vf1_m2_2gaus_p.at(i)->GetParameter(2));
    vf1_m2_2gaus_all.at(i)->SetParameter(15, vf1_m2_2gaus_p.at(i)->GetParameter(3));
    vf1_m2_2gaus_all.at(i)->SetParameter(17, vf1_m2_2gaus_p.at(i)->GetParameter(5));

    for (int k=0; k<niter; ++k){
      vh1_m2.at(i)->Fit(&*vf1_m2_2gaus_all.at(i), "RNQ");
    }

    vsigm_m2_p_y.push_back(vf1_m2_1gaus_p.at(i)->GetParameter(2));
    vsigm_m2_p_ey.push_back(vf1_m2_1gaus_p.at(i)->GetParError(2));
    // vsigm_m2_p_y.push_back(vf1_m2_2gaus_p.at(i)->GetParameter(2));
    // vsigm_m2_p_ey.push_back(vf1_m2_2gaus_p.at(i)->GetParError(2));
    // vsigm_m2_p_y.push_back(vf1_m2_2gaus_all.at(i)->GetParameter(14));
    // vsigm_m2_p_ey.push_back(vf1_m2_2gaus_all.at(i)->GetParError(14));
  }

  auto gr_m2_sigm_p = std::make_unique<TGraphErrors>(vsigm_m2_p_x.size(),
                                                     &vsigm_m2_p_x[0], &vsigm_m2_p_y[0],
                                                     &vsigm_m2_p_ex[0], &vsigm_m2_p_ey[0]);

  gr_m2_sigm_p->SetName("gr_m2_sigm_p");
  gr_m2_sigm_p->SetTitle("#sigma_{m^{2}} of protons vs P/q;p/q, GeV/c;#sigma_{m^{2}}, (GeV/c^{2})^{2}");

  auto f1_m2_sigm_p = std::make_unique<TF1>("f1_m2_sigm_p", "pol2", 0.3, 2.);
  gr_m2_sigm_p->Fit(&*f1_m2_sigm_p, "RNQ");

  // Writing output
  std::cout << "Writing output" << std::endl;
  auto fo = new TFile(outFile_pid.c_str(), "recreate");
  fo->cd();


  h2_recDedxPq->Write();
  h2_recDedxNormPq->Write();
  h2_recDedxPq_proton->Write();
  h2_recDedxNormPq_proton->Write();
  f1_dedx_proton->Write();
  gr_dedx_sigm_p->Write();
  f1_dedx_sigm_p->Write();

  fo->mkdir("dedx_slices");
  fo->cd("dedx_slices");
  for (int i=0; i<vh1_dedx_p.size(); ++i){
    vh1_dedx_p.at(i)->Write();
    vf1_dedx_1gaus_p.at(i)->Write();
  }

  fo->cd();
  h2_recM2Pq->Write();
  gr_m2_sigm_p->Write();
  f1_m2_sigm_p->Write();
  
  fo->mkdir("m2_slices");
  fo->cd("m2_slices");
  for (int i=0; i<vh1_m2.size(); ++i){
    vh1_m2.at(i)->Write();
    vf1_m2_1gaus_p.at(i)->Write();
    vf1_m2_1gaus_pip.at(i)->Write();
    vf1_m2_1gaus_kap.at(i)->Write();
    vf1_m2_2gaus_p.at(i)->Write();
    vf1_m2_2gaus_pip.at(i)->Write();
    vf1_m2_2gaus_kap.at(i)->Write();
    vf1_m2_2gaus_all.at(i)->Write();
  }
  fo->Close();
}
