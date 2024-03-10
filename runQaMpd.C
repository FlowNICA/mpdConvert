using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

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

RVec<float> getEta(vector<fourVector> _p)
{
  vector <float> pt_;
  for (auto& mom:_p)
    pt_.push_back(mom.Eta());
  return pt_;
}

RVec<float> getPhi(vector<fourVector> _p)
{
  vector <float> pt_;
  for (auto& mom:_p)
    pt_.push_back(mom.Phi());
  return pt_;
}

RVec<float> getDcaMag(vector<XYZVector> _dca)
{
  vector <float> dca_;
  for (auto& dca:_dca)
    dca_.push_back(sqrt(dca.Mag2()));
  return dca_;
}

RVec<int> isGoodTrack(vector<fourVector> _p, RVec<int> _nhits, RVec<float> _dca)
{
  vector<int> vec_tracks;
  for(int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto nhits = std::round(_nhits.at(i));
    auto dca = _dca.at(i);
    if (nhits > 16 && dca < 100.)
      vec_tracks.push_back(1);
    else
      vec_tracks.push_back(0);
  }
  return vec_tracks;
}

RVec<int> isGoodTrack4Protons(vector<fourVector> _p, RVec<int> _nhits, RVec<float> _dca)
{
  vector<int> vec_tracks;
  for(int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto nhits = std::round(_nhits.at(i));
    auto dca = _dca.at(i);
    if (nhits > 27 && dca < 1.)
      vec_tracks.push_back(1);
    else
      vec_tracks.push_back(0);
  }
  return vec_tracks;
}

RVec<int> isGoodTrack4Pions(vector<fourVector> _p, RVec<int> _nhits, RVec<float> _dca)
{
  vector<int> vec_tracks;
  for(int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto nhits = std::round(_nhits.at(i));
    auto dca = _dca.at(i);
    if (nhits > 22 && dca < 3.5)
      vec_tracks.push_back(1);
    else
      vec_tracks.push_back(0);
  }
  return vec_tracks;
}

RVec<int> isGoodPosEta(vector<fourVector> _p, RVec<int> _isGoodTrack)
{
  vector<int> vec_tracks;
  for(int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto isgood = _isGoodTrack.at(i);
    if (isgood && mom.Eta() > 0)
      vec_tracks.push_back(1);
    else
      vec_tracks.push_back(0);
  }
  return vec_tracks;
}

RVec<int> isGoodNegEta(vector<fourVector> _p, RVec<int> _isGoodTrack)
{
  vector<int> vec_tracks;
  for(int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto isgood = _isGoodTrack.at(i);
    if (isgood && mom.Eta() < 0)
      vec_tracks.push_back(1);
    else
      vec_tracks.push_back(0);
  }
  return vec_tracks;
}

RVec<float> trGoodPid(RVec<float> _tracks, RVec<int> _isGoodTrack, RVec<int> _isPrimary, RVec<int> _isPid){
  std::vector<float> vec_tracks(_tracks.size(), -999.);
  for (int i=0; i<_tracks.size(); ++i){
    auto track = _tracks.at(i);
    auto goodtrack = std::round(_isGoodTrack.at(i));
    auto isPrimary = std::round(_isPrimary.at(i));
    auto isPid = std::round(_isPid.at(i));
    if (goodtrack && isPrimary && isPid)
      vec_tracks.at(i) = track;
  }
  return vec_tracks;
}

RVec<float> trPrimPid(RVec<float> _tracks, RVec<int> _isPrimary, RVec<int> _isPid){
  std::vector<float> vec_tracks(_tracks.size(), -999.);
  for (int i=0; i<_tracks.size(); ++i){
    auto track = _tracks.at(i);
    auto isPrimary = std::round(_isPrimary.at(i));
    auto isPid = std::round(_isPid.at(i));
    if (isPrimary && isPid)
      vec_tracks.at(i) = track;
  }
  return vec_tracks;
}

RVec<float> trPrimFHCalPid(RVec<float> _tracks, RVec<int> _isPrimary, RVec<int> _isPid, RVec<bool> _hasHit){
  std::vector<float> vec_tracks(_tracks.size(), -999.);
  for (int i=0; i<_tracks.size(); ++i){
    auto track = _tracks.at(i);
    auto isPrimary = std::round(_isPrimary.at(i));
    auto isPid = std::round(_isPid.at(i));
    auto hit = _hasHit.at(i);
    if (isPrimary && isPid && hit)
      vec_tracks.at(i) = track;
  }
  return vec_tracks;
}

RVec<float> getGoodTrackResolutionPid(RVec<float> _reco_tracks, RVec<float> _sim_tracks, RVec<int> _sim_ids, RVec<int> _isGoodTrack, RVec<int> _isPrimary, RVec<int> _isPid)
{
  std::vector<float> vec_delta(_reco_tracks.size(), -999.);
  for (int i = 0; i < _sim_ids.size(); ++i){
    auto rec_track = _reco_tracks.at(i);
    auto sim_id = _sim_ids.at(i);
    if (sim_id >= _sim_tracks.size()){
      continue;
    }
    if (sim_id < 0){
      continue;
    }
    auto sim_track = _sim_tracks.at(sim_id);
    auto delta = (sim_track <= std::numeric_limits<float>::min() || sim_track == -999. || rec_track == -999.) ? -999. : abs(rec_track - sim_track)/sim_track;
    auto goodtrack = std::round(_isGoodTrack.at(i));
    auto isPrimary = std::round(_isPrimary.at(i));
    auto isPid = std::round(_isPid.at(i));
    if (goodtrack && isPrimary && isPid)
      vec_delta.at(i) = delta; 
  }
  return vec_delta;
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
  vector <RResultPtr<::TProfile >> profs;
  vector <RResultPtr<::TProfile2D >> profs2d;

  auto dd = d
    .Filter("mcB<16.")
    .Define("recDca", getDcaMag, {"recoGlobalDca"})
    .Define("isGoodTrack", isGoodTrack, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodTrackProton", isGoodTrack4Protons, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodTrackPionP", isGoodTrack4Pions, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodTrackPionM", isGoodTrack4Pions, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodTrackPions", isGoodTrack4Pions, {"recoGlobalMom", "recoGlobalNhits", "recDca"})
    .Define("isGoodPosEta", isGoodPosEta, {"recoGlobalMom", "isGoodTrack"})
    .Define("isGoodNegEta", isGoodNegEta, {"recoGlobalMom", "isGoodTrack"})
    .Define("isPrimary", "recoGlobalSimMotherId==-1")
    .Define("isProton", "recoGlobalSimPdg==2212")
    .Define("isPionP", "recoGlobalSimPdg==211")
    .Define("isPionM", "recoGlobalSimPdg==-211")
    .Define("isPions", "abs(recoGlobalSimPdg)==211")
    .Define("isKaonP", "recoGlobalSimPdg==321")
    .Define("isKaonM", "recoGlobalSimPdg==-321")
    .Define("isKaons", "abs(recoGlobalSimPdg)==321")
    .Define("recPt", getPt, {"recoGlobalMom"})
    .Define("recEta", getEta, {"recoGlobalMom"})
    .Define("recPhi", getPhi, {"recoGlobalMom"})
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
    .Define("recPtGoodProton", trGoodPid, {"recPt", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recPtGoodPionP", trGoodPid, {"recPt", "isGoodTrackPionP", "isPrimary", "isPionP"})
    .Define("recPtGoodPionM", trGoodPid, {"recPt", "isGoodTrackPionM", "isPrimary", "isPionM"})
    .Define("recPtGoodPions", trGoodPid, {"recPt", "isGoodTrackPions", "isPrimary", "isPions"})
    .Define("recPtGoodKaonP", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isKaonP"})
    .Define("recPtGoodKaonM", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isKaonM"})
    .Define("recPtGoodKaons", trGoodPid, {"recPt", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recYGoodProton", trGoodPid, {"recY", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recYGoodPionP", trGoodPid, {"recY", "isGoodTrackPionP", "isPrimary", "isPionP"})
    .Define("recYGoodPionM", trGoodPid, {"recY", "isGoodTrackPionM", "isPrimary", "isPionM"})
    .Define("recYGoodPions", trGoodPid, {"recY", "isGoodTrackPions", "isPrimary", "isPions"})
    .Define("recYGoodKaonP", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isKaonP"})
    .Define("recYGoodKaonM", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isKaonM"})
    .Define("recYGoodKaons", trGoodPid, {"recY", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recEtaGoodProton", trGoodPid,{"recEta", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recEtaGoodPionP", trGoodPid, {"recEta", "isGoodTrackPionP", "isPrimary", "isPionP"})
    .Define("recEtaGoodPionM", trGoodPid, {"recEta", "isGoodTrackPionM", "isPrimary", "isPionM"})
    .Define("recEtaGoodPions", trGoodPid, {"recEta", "isGoodTrackPions", "isPrimary", "isPions"})
    .Define("recEtaGoodKaonP", trGoodPid, {"recEta", "isGoodTrack", "isPrimary", "isKaonP"})
    .Define("recEtaGoodKaonM", trGoodPid, {"recEta", "isGoodTrack", "isPrimary", "isKaonM"})
    .Define("recEtaGoodKaons", trGoodPid, {"recEta", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("isSimPrimary", "simMotherId==-1")
    .Define("isSimProton", "simPdg==2212")
    .Define("isSimPionP", "simPdg==211")
    .Define("isSimPionM", "simPdg==-211")
    .Define("isSimPions", "abs(simPdg)==211")
    .Define("isSimKaonP", "simPdg==321")
    .Define("isSimKaonM", "simPdg==-321")
    .Define("isSimKaons", "abs(simPdg)==321")
    .Define("simPt", getPt, {"simMom"})
    .Define("simEta", getEta, {"simMom"})
    .Define("simPhi", getPhi, {"simMom"})
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
    .Define("simPtGoodPions",  trPrimPid, {"simPt", "isSimPrimary", "isSimPions"})
    .Define("simPtGoodKaonP",  trPrimPid, {"simPt", "isSimPrimary", "isSimKaonP"})
    .Define("simPtGoodKaonM",  trPrimPid, {"simPt", "isSimPrimary", "isSimKaonM"})
    .Define("simPtGoodKaons",  trPrimPid, {"simPt", "isSimPrimary", "isSimKaons"})
    .Define("simYGoodProton",  trPrimPid, {"simY",  "isSimPrimary", "isSimProton"})
    .Define("simYGoodPionP",   trPrimPid, {"simY",  "isSimPrimary", "isSimPionP"})
    .Define("simYGoodPionM",   trPrimPid, {"simY",  "isSimPrimary", "isSimPionM"})
    .Define("simYGoodPions",   trPrimPid, {"simY",  "isSimPrimary", "isSimPions"})
    .Define("simYGoodKaonP",   trPrimPid, {"simY",  "isSimPrimary", "isSimKaonP"})
    .Define("simYGoodKaonM",   trPrimPid, {"simY",  "isSimPrimary", "isSimKaonM"})
    .Define("simYGoodKaons",   trPrimPid, {"simY",  "isSimPrimary", "isSimKaons"})
    .Define("simEtaGoodProton", trPrimPid, {"simEta", "isSimPrimary", "isSimProton"})
    .Define("simEtaGoodPionP",  trPrimPid, {"simEta", "isSimPrimary", "isSimPionP"})
    .Define("simEtaGoodPionM",  trPrimPid, {"simEta", "isSimPrimary", "isSimPionM"})
    .Define("simEtaGoodPions",  trPrimPid, {"simEta", "isSimPrimary", "isSimPions"})
    .Define("simEtaGoodKaonP",  trPrimPid, {"simEta", "isSimPrimary", "isSimKaonP"})
    .Define("simEtaGoodKaonM",  trPrimPid, {"simEta", "isSimPrimary", "isSimKaonM"})
    .Define("simEtaGoodKaons",  trPrimPid, {"simEta", "isSimPrimary", "isSimKaons"})
    .Define("simPtGoodFHCalProton",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimProton", "simHasHitFHCal"})
    .Define("simPtGoodFHCalPionP",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simPtGoodFHCalPionM",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simPtGoodFHCalKaonP",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simPtGoodFHCalKaonM",  trPrimFHCalPid, {"simPt", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simYGoodFHCalProton",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimProton", "simHasHitFHCal"})
    .Define("simYGoodFHCalPionP",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simYGoodFHCalPionM",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simYGoodFHCalKaonP",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simYGoodFHCalKaonM",  trPrimFHCalPid, {"simY", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalProton",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimProton", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalPionP",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalPionM",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimPionP", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalKaonP",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("simEtaGoodFHCalKaonM",  trPrimFHCalPid, {"simEta", "isSimPrimary", "isSimKaonP", "simHasHitFHCal"})
    .Define("recDPtGoodProton", getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recDPtGoodPionP",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrackPionP",  "isPrimary", "isPionP"}) 
    .Define("recDPtGoodPionM",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrackPionM",  "isPrimary", "isPionM"}) 
    .Define("recDPtGoodPions",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrackPions",  "isPrimary", "isPions"}) 
    .Define("recDPtGoodKaonP",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonP"}) 
    .Define("recDPtGoodKaonM",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonM"}) 
    .Define("recDPtGoodKaons",  getGoodTrackResolutionPid, {"recPt", "simPt", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recDYGoodProton", getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recDYGoodPionP",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrackPionP",  "isPrimary", "isPionP"}) 
    .Define("recDYGoodPionM",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrackPionM",  "isPrimary", "isPionM"}) 
    .Define("recDYGoodPions",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrackPions",  "isPrimary", "isPions"}) 
    .Define("recDYGoodKaonP",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonP"}) 
    .Define("recDYGoodKaonM",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonM"}) 
    .Define("recDYGoodKaons",  getGoodTrackResolutionPid, {"recY", "simY", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recDEtaGoodProton", getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recDEtaGoodPionP",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrackPionP",  "isPrimary", "isPionP"}) 
    .Define("recDEtaGoodPionM",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrackPionM",  "isPrimary", "isPionM"}) 
    .Define("recDEtaGoodPions",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrackPions",  "isPrimary", "isPions"}) 
    .Define("recDEtaGoodKaonP",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonP"}) 
    .Define("recDEtaGoodKaonM",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonM"}) 
    .Define("recDEtaGoodKaons",  getGoodTrackResolutionPid, {"recEta", "simEta", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaons"})
    .Define("recDPhiGoodProton", getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrackProton", "isPrimary", "isProton"})
    .Define("recDPhiGoodPionP",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrackPionP",  "isPrimary", "isPionP"}) 
    .Define("recDPhiGoodPionM",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrackPionM",  "isPrimary", "isPionM"}) 
    .Define("recDPhiGoodPions",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrackPions",  "isPrimary", "isPions"}) 
    .Define("recDPhiGoodKaonP",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonP"}) 
    .Define("recDPhiGoodKaonM",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaonM"}) 
    .Define("recDPhiGoodKaons",  getGoodTrackResolutionPid, {"recPhi", "simPhi", "recoGlobalSimIndex", "isGoodTrack", "isPrimary", "isKaons"})
  ;

  dd.Foreach([](ULong64_t evtId){if (evtId % 100 == 0) cout << "\r" << evtId;}, {"rdfentry_"});

  // Make lists of histograms for QA
  hists.push_back(dd.Histo1D({"h1_recVtx_X","Reconstructed vertex X;x (cm)",500,-1,1}, "recoPrimVtxX")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_Y","Reconstructed vertex Y;y (cm)",500,-1,1}, "recoPrimVtxY")); 
  hists.push_back(dd.Histo1D({"h1_recVtx_Z","Reconstructed vertex Z;z (cm)",500,-1,1}, "recoPrimVtxZ"));
  profs.push_back(dd.Profile1D({"p1_DPt_recNhits_proton", "Pt-resolution for reconstructed protons vs. Nhits;N_{hits};#deltap_{T}", 50, 0., 50., -998., 1000.}, "recoGlobalNhits", "recDPtGoodProton"));
  profs.push_back(dd.Profile1D({"p1_DPt_recDCA_proton", "Pt-resolution for reconstructed protons vs. DCA;DCA (cm);#deltap_{T}", 50, 0., 5., -998., 1000.}, "recDca", "recDPtGoodProton"));
  profs.push_back(dd.Profile1D({"p1_DPt_recChi2_proton", "Pt-resolution for reconstructed protons vs. Chi2;#chi_{2}/ndf;#deltap_{T}", 5000, 0., 5000., -998., 1000.}, "recoGlobalChi2", "recDPtGoodProton"));
  profs.push_back(dd.Profile1D({"p1_DPt_recNhits_pionP", "Pt-resolution for reconstructed pions (#pi^{+}) vs. Nhits;N_{hits};#deltap_{T}", 50, 0., 50., -998., 1000.}, "recoGlobalNhits", "recDPtGoodPionP"));
  profs.push_back(dd.Profile1D({"p1_DPt_recDCA_pionP", "Pt-resolution for reconstructed pions (#pi^{+}) vs. DCA;DCA (cm);#deltap_{T}", 50, 0., 5., -998., 1000.}, "recDca", "recDPtGoodPionP"));
  profs.push_back(dd.Profile1D({"p1_DPt_recChi2_pionP", "Pt-resolution for reconstructed pions (#pi^{+}) vs. Chi2;#chi_{2}/ndf;#deltap_{T}", 5000, 0., 5000., -998., 1000.}, "recoGlobalChi2", "recDPtGoodPionP"));
  profs.push_back(dd.Profile1D({"p1_DPt_recNhits_pionM", "Pt-resolution for reconstructed pions (#pi^{-}) vs. Nhits;N_{hits};#deltap_{T}", 50, 0., 50., -998., 1000.}, "recoGlobalNhits", "recDPtGoodPionM"));
  profs.push_back(dd.Profile1D({"p1_DPt_recDCA_pionM", "Pt-resolution for reconstructed pions (#pi^{-}) vs. DCA;DCA (cm);#deltap_{T}", 50, 0., 5., -998., 1000.}, "recDca", "recDPtGoodPionM"));
  profs.push_back(dd.Profile1D({"p1_DPt_recChi2_pionM", "Pt-resolution for reconstructed pions (#pi^{-}) vs. Chi2;#chi_{2}/ndf;#deltap_{T}", 5000, 0., 5000., -998., 1000.}, "recoGlobalChi2", "recDPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recVtx_XY","Reconstructed vertex XY;x (cm);y (cm)",500,-1,1,500,-1,1}, "recoPrimVtxX", "recoPrimVtxY"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_proton", "Simulated protons Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalProton", "simPtGoodFHCalProton"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_pionP", "Simulated pions (#pi^{+}) Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalPionP", "simPtGoodFHCalPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_pionM", "Simulated pions (#pi^{+}) Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalPionM", "simPtGoodFHCalPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_kaonP", "Simulated Kaons (#pi^{+}) Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalKaonP", "simPtGoodFHCalKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_inFHCal_kaonM", "Simulated Kaons (#pi^{+}) Ycm-pT in FHCal;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodFHCalKaonM", "simPtGoodFHCalKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_proton", "Simulated protons #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalProton", "simPtGoodFHCalProton"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_pionP", "Simulated pions (#pi^{+}) #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalPionP", "simPtGoodFHCalPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_pionM", "Simulated pions (#pi^{+}) #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalPionM", "simPtGoodFHCalPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_kaonP", "Simulated Kaons (#pi^{+}) #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalKaonP", "simPtGoodFHCalKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_inFHCal_kaonM", "Simulated Kaons (#pi^{+}) #eta-pT in FHCal;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodFHCalKaonM", "simPtGoodFHCalKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_proton", "Reconstructed protons Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_pionP", "Reconstructed pions (#pi^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_pionM", "Reconstructed pions (#pi^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_pions", "Reconstructed pions (#pi^{#pm}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_kaonP", "Reconstructed kaons (K^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_kaonM", "Reconstructed kaons (K^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_recYPt_kaons", "Reconstructed kaons (K^{#pm}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_proton", "Reconstructed protons #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodProton", "recPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_pionP", "Reconstructed pions (#pi^{+}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodPionP", "recPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_pionM", "Reconstructed pions (#pi^{-}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodPionM", "recPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_pions", "Reconstructed pions (#pi^{#pm}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodPions", "recPtGoodPions"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_kaonP", "Reconstructed kaons (K^{+}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodKaonP", "recPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_kaonM", "Reconstructed kaons (K^{-}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodKaonM", "recPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_recEtaPt_kaons", "Reconstructed kaons (K^{#pm}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "recEtaGoodKaons", "recPtGoodKaons"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_proton", "Simulated protons Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodProton", "simPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_pionP", "Simulated pions (#pi^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodPionP", "simPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_pionM", "Simulated pions (#pi^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodPionM", "simPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_pions", "Simulated pions (#pi^{#pm}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodPions", "simPtGoodPions"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_kaonP", "Simulated kaons (K^{+}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodKaonP", "simPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_kaonM", "Simulated kaons (K^{-}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodKaonM", "simPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_simYPt_kaons", "Simulated kaons (K^{#pm}) Ycm-pT;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "simYGoodKaons", "simPtGoodKaons"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_proton", "Simulated protons #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodProton", "simPtGoodProton"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_pionP", "Simulated pions (#pi^{+}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodPionP", "simPtGoodPionP"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_pionM", "Simulated pions (#pi^{-}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodPionM", "simPtGoodPionM"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_pions", "Simulated pions (#pi^{#pm}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodPions", "simPtGoodPions"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_kaonP", "Simulated kaons (K^{+}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodKaonP", "simPtGoodKaonP"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_kaonM", "Simulated kaons (K^{-}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodKaonM", "simPtGoodKaonM"));
  hists2d.push_back(dd.Histo2D({"h2_simEtaPt_kaons", "Simulated kaons (K^{#pm}) #eta-pT;#eta; p_{T} (GeV/c)", 340, -0.2, 3.2, 300, 0., 3.}, "simEtaGoodKaons", "simPtGoodKaons"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_proton", "Pt-resolution for reconstructed protons in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton", "recDPtGoodProton"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_pionP", "Pt-resolution for reconstructed pions (#pi^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP", "recDPtGoodPionP"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_pionM", "Pt-resolution for reconstructed pions (#pi^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM", "recDPtGoodPionM"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_pions", "Pt-resolution for reconstructed pions (#pi^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions", "recDPtGoodPions"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_kaonP", "Pt-resolution for reconstructed kaons (K^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP", "recDPtGoodKaonP"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_kaonM", "Pt-resolution for reconstructed kaons (K^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM", "recDPtGoodKaonM"));
  profs2d.push_back(dd.Profile2D({"p2_DPt_recYPt_kaons", "Pt-resolution for reconstructed kaons (K^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons", "recDPtGoodKaons"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_proton", "Phi-resolution for reconstructed protons in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton", "recDPhiGoodProton"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_pionP", "Phi-resolution for reconstructed pions (#pi^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP", "recDPhiGoodPionP"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_pionM", "Phi-resolution for reconstructed pions (#pi^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM", "recDPhiGoodPionM"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_pions", "Phi-resolution for reconstructed pions (#pi^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions", "recDPhiGoodPions"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_kaonP", "Phi-resolution for reconstructed kaons (K^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP", "recDPhiGoodKaonP"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_kaonM", "Phi-resolution for reconstructed kaons (K^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM", "recDPhiGoodKaonM"));
  profs2d.push_back(dd.Profile2D({"p2_DPhi_recYPt_kaons", "Phi-resolution for reconstructed kaons (K^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons", "recDPhiGoodKaons"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_proton", "Ycm-resolution for reconstructed protons in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton", "recDYGoodProton"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_pionP", "Ycm-resolution for reconstructed pions (#pi^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP", "recDYGoodPionP"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_pionM", "Ycm-resolution for reconstructed pions (#pi^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM", "recDYGoodPionM"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_pions", "Ycm-resolution for reconstructed pions (#pi^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions", "recDYGoodPions"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_kaonP", "Ycm-resolution for reconstructed kaons (K^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP", "recDYGoodKaonP"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_kaonM", "Ycm-resolution for reconstructed kaons (K^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM", "recDYGoodKaonM"));
  profs2d.push_back(dd.Profile2D({"p2_DY_recYPt_kaons", "Ycm-resolution for reconstructed kaons (K^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons", "recDYGoodKaons"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_proton", "Eta-resolution for reconstructed protons in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodProton", "recPtGoodProton", "recDEtaGoodProton"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_pionP", "Eta-resolution for reconstructed pions (#pi^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionP", "recPtGoodPionP", "recDEtaGoodPionP"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_pionM", "Eta-resolution for reconstructed pions (#pi^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPionM", "recPtGoodPionM", "recDEtaGoodPionM"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_pions", "Eta-resolution for reconstructed pions (#pi^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodPions", "recPtGoodPions", "recDEtaGoodPions"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_kaonP", "Eta-resolution for reconstructed kaons (K^{+}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonP", "recPtGoodKaonP", "recDEtaGoodKaonP"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_kaonM", "Eta-resolution for reconstructed kaons (K^{-}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaonM", "recPtGoodKaonM", "recDEtaGoodKaonM"));
  profs2d.push_back(dd.Profile2D({"p2_DEta_recYPt_kaons", "Eta-resolution for reconstructed kaons (K^{#pm}) in Ycm-pT plane;y_{CM}; p_{T} (GeV/c)", 300, -1.5, 1.5, 300, 0., 3.}, "recYGoodKaons", "recPtGoodKaons", "recDEtaGoodKaons"));

  // Write QA histograms to the output file
  fOut.cd();
  for (auto& hist:hists)
    hist->Write();
  for (auto& hist:hists2d)
    hist->Write();
  for (auto& hist:profs)
    hist->Write();
  for (auto& hist:profs2d)
    hist->Write();
  fOut.Close();

  std::cout << std::endl;
  timer.Stop();
  timer.Print();
}
