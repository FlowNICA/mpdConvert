using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

//MpdFieldMap* magField{nullptr};
MpdConstField *magField{nullptr};

TChain* makeChain(string& filename, const char* treename) {
  cout << "Adding files to chain:" << endl;
  TChain *chain = new TChain(treename);
  if (filename.rfind(".root") < filename.size())
    chain->Add(filename.data());
  else {
    TFileCollection fc("fc", "", filename.c_str());
    chain->AddFileInfoList((TCollection*)fc.GetList());
  }
  chain->ls();
  return chain;
}

RVec<float> getPt(vector<fourVector> _p)
{
  vector <float> pt_;
  for (auto& mom:_p)
    pt_.push_back(mom.Pt());
  return pt_;
}

RVec<int> isGoodTrack(vector<fourVector> _p, RVec<int> _nhits)
{
  vector<int> vec_tracks;
  for(int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto nhits = _nhits.at(i);
    if (nhits > 16)
      vec_tracks.push_back(1);
    else
      vec_tracks.push_back(0);
  }
  return vec_tracks;
}

RVec<float> trGoodPid(RVec<float> _tracks, RVec<int> _isGoodTrack, RVec<int> _isPrimary, RVec<int> _isPid){
  std::vector<float> vec_tracks;
  for (int i=0; i<_tracks.size(); ++i){
    auto track = _tracks.at(i);
    auto goodtrack = _isGoodTrack.at(i);
    auto isPrimary = _isPrimary.at(i);
    auto isPid = _isPid.at(i);
    if (goodtrack && isPrimary && isPid)
      vec_tracks.push_back(track);
  }
  return vec_tracks;
}

RVec<float> trPrimPid(RVec<float> _tracks, RVec<int> _isPrimary, RVec<int> _isPid){
  std::vector<float> vec_tracks;
  for (int i=0; i<_tracks.size(); ++i){
    auto track = _tracks.at(i);
    auto isPrimary = _isPrimary.at(i);
    auto isPid = _isPid.at(i);
    if (isPrimary && isPid)
      vec_tracks.push_back(track);
  }
  return vec_tracks;
}

void runQaMpd(string fileIn="", string fileOut="", std::string cm_energy="2.5", std::string str_nucleus_mass="209")
{
  TStopwatch timer;
  timer.Start();

  const double sNN = std::stod( cm_energy ); // in GeV
  const double M = 0.938; // in GeV/c^2
  const double T = sNN*sNN/(2.*M) - 2.*M;
  const double E = T + M;
  const double P = sqrt( E*E - M*M );
  const double Y_BEAM = 0.25 * log( (E + P) / (E - P) );
  const double nucleus_mass = std::stod(str_nucleus_mass);
  const double NUCLEUS_RADIUS = 1.25 * pow( nucleus_mass, 1.0 / 3.0 );

  std::cout << "sqrtSnn = " << sNN << " GeV; T = " << T << "A GeV; Y_BEAM = " << Y_BEAM << std::endl;
  std::cout << "A = " << nucleus_mass << "; R = " << NUCLEUS_RADIUS << std::endl;

  ROOT::RDataFrame d("t", fileIn.c_str());
  TFile fOut(fileOut.c_str(),"recreate");

  vector <RResultPtr<::TH1D >> hists;
  vector <RResultPtr<::TH2D >> hists2d;

  auto dd = d
    .Filter("mcB<16.")
    .Define("isGoodTrack", isGoodTrack, {"recoGlobalMom", "recoGlobalNhits"})
    .Define("isPrimary", "recoGlobalSimMotherId==-1")
    .Define("isProton", "recoGlobalSimPdg==2212")
    .Define("isPionP", "recoGlobalSimPdg==211")
    .Define("isPionM", "recoGlobalSimPdg==-211")
    .Define("isKaonP", "recoGlobalSimPdg==321")
    .Define("isKaonM", "recoGlobalSimPdg==-321")
    .Define("recPt", getPt, {"recoGlobalMom"})
    .Define("recY", [Y_BEAM]( const RVec<int> vec_pdg, vector<fourVector> vec_momentum ){
      RVec<float> vec_y;
      vec_y.reserve( vec_pdg.size() );
      for( int i=0; i<vec_pdg.size(); ++i ){
        auto pdg = vec_pdg.at(i);
        if( pdg == 0 ){
          vec_y.push_back( -999. );
          continue;
        }
        auto pz = vec_momentum.at(i).Pz();
        auto p = vec_momentum.at(i).P();
        TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg);
        if (!particle) {
          vec_y.push_back( -999. );
          continue;
        }
        auto m = particle->Mass();
        auto E = sqrt( p*p + m*m );
        auto y = 0.5 * log( (E+pz)/(E-pz) ) - Y_BEAM;
        vec_y.push_back(y);
      }
      return vec_y;
    },{"recoGlobalSimPdg", "recoGlobalMom"})
    .Define("recPtGoodProton", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isProton"})
    .Define("recPtGoodPionP", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isPionP"})
    .Define("recPtGoodPionM", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isPionM"})
    .Define("recPtGoodKaonP", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isKaonP"})
    .Define("recPtGoodKaonM", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isKaonM"})
    .Define("recYGoodProton", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isProton"})
    .Define("recYGoodPionP", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isPionP"})
    .Define("recYGoodPionM", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isPionM"})
    .Define("recYGoodKaonP", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isKaonP"})
    .Define("recYGoodKaonM", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isKaonM"})
    .Define("isSimPrimary", "simMotherId==-1")
    .Define("isSimProton", "simPdg==2212")
    .Define("isSimPionP", "simPdg==211")
    .Define("isSimPionM", "simPdg==-211")
    .Define("isSimKaonP", "simPdg==321")
    .Define("isSimKaonM", "simPdg==-321")
    .Define("simPt", getPt, {"simMom"})
    .Define("simY", [Y_BEAM]( const RVec<int> vec_pdg, vector<fourVector> vec_momentum ){
      RVec<float> vec_y;
      vec_y.reserve( vec_pdg.size() );
      for( int i=0; i<vec_pdg.size(); ++i ){
        auto pdg = vec_pdg.at(i);
        if( pdg == 0 ){
          vec_y.push_back( -999. );
          continue;
        }
        auto pz = vec_momentum.at(i).Pz();
        auto p = vec_momentum.at(i).P();
        TParticlePDG *particle = TDatabasePDG::Instance()->GetParticle(pdg);
        if (!particle) {
          vec_y.push_back( -999. );
          continue;
        }
        auto m = particle->Mass();
        auto E = sqrt( p*p + m*m );
        auto y = 0.5 * log( (E+pz)/(E-pz) ) - Y_BEAM;
        vec_y.push_back(y);
      }
      return vec_y;
    },{"simPdg", "simMom"})
    .Define("simPtGoodProton", trPrimPid, {"simPt", "isSimPrimary", "isSimProton"})
    .Define("simPtGoodPionP",  trPrimPid, {"simPt", "isSimPrimary", "isSimPionP"})
    .Define("simPtGoodPionM",  trPrimPid, {"simPt", "isSimPrimary", "isSimPionM"})
    .Define("simPtGoodKaonP",  trPrimPid, {"simPt", "isSimPrimary", "isSimKaonP"})
    .Define("simPtGoodKaonM",  trPrimPid, {"simPt", "isSimPrimary", "isSimKaonM"})
    .Define("simYGoodProton",  trPrimPid, {"simY",  "isSimPrimary", "isSimProton"})
    .Define("simYGoodPionP",   trPrimPid, {"simY",  "isSimPrimary", "isSimPionP"})
    .Define("simYGoodPionM",   trPrimPid, {"simY",  "isSimPrimary", "isSimPionM"})
    .Define("simYGoodKaonP",   trPrimPid, {"simY",  "isSimPrimary", "isSimKaonP"})
    .Define("simYGoodKaonM",   trPrimPid, {"simY",  "isSimPrimary", "isSimKaonM"})
  ;

  dd.Foreach([](ULong64_t evtId){if (evtId % 100 == 0) cout << "\r" << evtId;}, {"rdfentry_"});

  // Make lists of histograms for QA
  hists2d.push_back(dd.Histo2D({"h2_recVtx_XY","Reconstructed vertex XY;x (cm);y (cm)",500,-1,1,500,-1,1}, "recoPrimVtxX", "recoPrimVtxY")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_X","Reconstructed vertex X;x (cm)",500,-1,1}, "recoPrimVtxX")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_Y","Reconstructed vertex Y;y (cm)",500,-1,1}, "recoPrimVtxY")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_Z","Reconstructed vertex Z;z (cm)",500,-1,1}, "recoPrimVtxZ"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_proton", "Reconstructed protons Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_pionP", "Reconstructed pions (#pi^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_pionM", "Reconstructed pions (#pi^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_kaonP", "Reconstructed kaons (K^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_kaonM", "Reconstructed kaons (K^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_proton", "Simulated protons Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodProton", "simPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_pionP", "Simulated pions (#pi^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodPionP", "simPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_pionM", "Simulated pions (#pi^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodPionM", "simPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_kaonP", "Simulated kaons (K^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodKaonP", "simPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_kaonM", "Simulated kaons (K^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodKaonM", "simPtGoodKaonM"));

  // Write QA histograms to the output file
  fOut.cd();
  for (auto& hist:hists)
    hist->Write();
  for (auto& hist:hists2d)
    hist->Write();
  fOut.Close();

  std::cout << std::endl;
  timer.Stop();
  timer.Print();
}
