1;115;0cfloat EnCorr(float e);

void readEmcFromDst(std::string inFileName, std::string outFileName)
{
  TStopwatch timer;
    timer.Start();

    int A1 = 124, Z1 = 54; // Xe
    int A2 = 184, Z2 = 74; // W

    std::cout << "A1 = " << A1 << ", Z1 = " << Z1 << ", N1 = " << (A1-Z1) << std::endl;
    std::cout << "A2 = " << A2 << ", Z2 = " << Z2 << ", N2 = " << (A2-Z2) << std::endl;

    TChain *dstTree = new TChain("mpdsim");
    std::ifstream infile(inFileName.c_str());
    std::string line;
    while (std::getline(infile,line)){
      dstTree->Add(line.c_str());
    }

    // Activate branches
    MpdEvent *event = nullptr;
    dstTree->SetBranchAddress("MPDEvent.", &event);
    FairMCEventHeader *MCHeader = nullptr;
    dstTree->SetBranchAddress("MCEventHeader.", &MCHeader);
    TClonesArray *MCTracks = nullptr;
    dstTree->SetBranchAddress("MCTrack", &MCTracks);
    TObjArray * mpdEMCClusters = new TObjArray() ;
    MpdEmcClusterKI* EMCCluster;
    dstTree->SetBranchAddress("EmcCluster",&mpdEMCClusters);
    TClonesArray *FHCalHits = nullptr;
    dstTree->SetBranchAddress("ZdcDigi", &FHCalHits);
    TClonesArray *vertexes = nullptr;
    dstTree->SetBranchAddress("Vertex", &vertexes);

    int N_b_bins = 200;
    int N_Efhcal_bins = 500;
    int N_Eemc_bins = 500;

    std::pair<double,double> b_bins = {0.,20.};
    std::pair<double,double> Efhcal_bins = {0.,10.};
    std::pair<double,double> Eemc_bins = {0., 100.};

    TH1D *h_primN_before = new TH1D("h_primN_before", "Number of primary neutrons before empty event cut", 500, 0., 500.);
    TH1D *h_primP_before = new TH1D("h_primP_before", "Number of primary protons before empty event cut", 500, 0., 500.);
    TH1D *h_primPart_before = new TH1D("h_primPart_before", "Number of primary particles before empty event cut", 1000, 0., 1000.);
    TH1D *h_primN_after = new TH1D("h_primN_after", "Number of primary neutrons after empty event cut", 500, 0., 500.);
    TH1D *h_primP_after = new TH1D("h_primP_after", "Number of primary protons after empty event cut", 500, 0., 500.);
    TH1D *h_primPart_after = new TH1D("h_primPart_after", "Number of primary particles after empty event cut", 1000, 0., 1000.);

    TH1D *h_bimp = new TH1D("h_bimp","impact parameter", N_b_bins, b_bins.first, b_bins.second);
    TH1D *h_Efhcal = new TH1D("h_Efhcal", "E_{FHCal}", N_Efhcal_bins, Efhcal_bins.first, Efhcal_bins.second);
    TH1D *h_Eemc = new TH1D("h_Eemc", "E_{EMC}", N_Eemc_bins, Eemc_bins.first, Eemc_bins.second);

    TH2D *h2_b_Efhcal = new TH2D("h2_b_Efhcal", "b vs E_{FHCal};b, fm; E_{FHCal}, GeV", N_b_bins, b_bins.first, b_bins.second, N_Efhcal_bins, Efhcal_bins.first, Efhcal_bins.second);
    TH2D *h2_b_Eemc = new TH2D("h2_b_Eemc", "b vs E_{EMC};b, fm;E_{EMC}, GeV", 200, b_bins.first, b_bins.second, N_Eemc_bins, Eemc_bins.first, Eemc_bins.second);
    TH2D *h2_Efhcal_Eemc = new TH2D("h2_Efhcal_Eemc", "E_{FHCal} vs E_{EMC};E_{FHCal}, GeV;E_{EMC}, GeV", N_Efhcal_bins, Efhcal_bins.first, Efhcal_bins.second, N_Eemc_bins, Eemc_bins.first, Eemc_bins.second);
    
    TH3D *h3_b_Efhcal_Eemc = new TH3D("h3_b_Efhcal_Eemc", "b vs E_{FHCal} vs E_{EMC};b, fm;E_{FHCal}, GeV;E_{EMC}, GeV", N_b_bins, b_bins.first, b_bins.second, N_Efhcal_bins, Efhcal_bins.first, Efhcal_bins.second, N_Eemc_bins, Eemc_bins.first, Eemc_bins.second);

    Int_t Num_Of_Modules = 90;

    Int_t events = dstTree->GetEntries();
    cout << " Number of events in DST file = " << events << endl;
    TVector3 primaryVertex;
    for (Int_t i = 0; i < events; i++) {
        dstTree->GetEntry(i);
        if (i&100) std::cout << "Event [" << i << "/" << events << "]" << std::endl;

        // Reading Reco Event
        MpdVertex *vertex = (MpdVertex *)vertexes->First();
        vertex->Position(primaryVertex);
        // Select events with VtxZ around target (Z= -85 cm)
        if (primaryVertex.Z() > -80.) continue;
        if (primaryVertex.Z() < -90.) continue;

        Int_t Ntracks = event->GetGlobalTracks()->GetEntriesFast();
        Int_t NMCtracks = MCTracks->GetEntriesFast();

        Int_t RefMult = 0, NmcNeutrons = 0, NmcProtons = 0, Nparticles = 0;

        // Mc tracks loop
        for (Int_t iMcTrack = 0; iMcTrack < NMCtracks; iMcTrack++) {
          MpdMCTrack *mctrack = (MpdMCTrack*) MCTracks->UncheckedAt(iMcTrack);
          if (!mctrack) continue;
          if (mctrack->GetMotherId() == -1 && mctrack->GetPdgCode() == 2112) NmcNeutrons++;
          if (mctrack->GetMotherId() == -1 && mctrack->GetPdgCode() == 2212) NmcProtons++;
          if (mctrack->GetMotherId() == -1) Nparticles++;
        }
        h_primN_before->Fill(NmcNeutrons);
        h_primP_before->Fill(NmcProtons);
        h_primPart_before->Fill(Nparticles);
        //std::cout << "\tNeutrons = " << NmcNeutrons << std::endl;
        if (NmcNeutrons == (A1+A2-(Z1+Z2))) continue;
        //if (NmcProtons == (Z1+Z2)) continue;
        //if (Nparticles <= (A1+A2)) continue;
        h_primN_after->Fill(NmcNeutrons);
        h_primP_after->Fill(NmcProtons);
        h_primPart_after->Fill(Nparticles);

        // Reco tracks loop
        for (Int_t iTrack = 0; iTrack < Ntracks; iTrack++) {
            MpdTrack* track = (MpdTrack*) event->GetGlobalTracks()->UncheckedAt(iTrack);
            
        } // track loop

        // FHCal loop
        Int_t number_of_FHCal_hits = FHCalHits->GetEntriesFast();
        float total_fhcal_energy = 0.;
        for(int ihit=0; ihit<number_of_FHCal_hits; ihit++) {
          MpdZdcDigi *FHCalHit = (MpdZdcDigi*) FHCalHits->UncheckedAt(ihit);
          Int_t DetId = FHCalHit->GetDetectorID();
          Int_t ModId = FHCalHit->GetModuleID()-1;
          Int_t ModNumber = ModId + (Num_Of_Modules/2) * (DetId-1);
          if (ModNumber == 22) continue;
          if (ModNumber == 67) continue;
          if (ModNumber > 45) continue;
          total_fhcal_energy += FHCalHit->GetELoss();
        }

        // EMC loop
        float total_emc_energy = 0.;
        for (long int j=0; j<mpdEMCClusters->GetEntries(); j++){
          EMCCluster = (MpdEmcClusterKI*)mpdEMCClusters->At(j);

          if (EMCCluster->GetMultiplicity() <= 1) continue;

          float conv = 1.0/0.3065; // only ~ 30% energy is collected
          float E_true = EnCorr( EMCCluster->GetE() * conv );
          total_emc_energy += E_true;
        }

        float bimp = MCHeader->GetB();
        h_bimp->Fill(bimp);
        h_Efhcal->Fill(total_fhcal_energy);
        h_Eemc->Fill(total_emc_energy);
        h2_b_Efhcal->Fill(bimp, total_fhcal_energy);
        h2_b_Eemc->Fill(bimp, total_emc_energy);
        h2_Efhcal_Eemc->Fill(total_fhcal_energy, total_emc_energy);
        h3_b_Efhcal_Eemc->Fill(bimp, total_fhcal_energy, total_emc_energy);
    } // event loop

    TFile *fo = new TFile(outFileName.c_str(),"recreate");
    fo->cd();

    h_bimp->Write();
    h_Efhcal->Write();
    h_Eemc->Write();

    h2_b_Efhcal->Write();
    h2_b_Eemc->Write();
    h2_Efhcal_Eemc->Write();

    h3_b_Efhcal_Eemc->Write();

    h_primN_before->Write();
    h_primP_before->Write();
    h_primPart_before->Write();
    h_primN_after->Write();
    h_primP_after->Write();
    h_primPart_after->Write();

    fo->Close();

    timer.Stop();
    timer.Print();
}

// (EMC) Function EnCorr() compensates energy for non-linearity (3% effect)
float EnCorr(float e) {
  float x = e;
  if (x > 2.5) x = 2.5; // parameterization is available up to 2.5 GeV only
  float corr[6] = { -1.321760e-002, 8.475623e-002, -9.007282e-002, 4.742396e-002, -1.279561e-002, 1.397394e-003 };
  float ecorr = corr[0] + corr[1]*x + corr[2]*x*x + corr[3]*x*x*x + corr[4]*x*x*x*x + corr[5]*x*x*x*x*x;
  return e/(1.0 + ecorr);
}
