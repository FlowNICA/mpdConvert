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

int getRefMultFxt(vector<fourVector> _p, RVec<int> _nhits)
{
  int Mult = 0;
  for (int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto nhits = _nhits.at(i);
    if (nhits<=16) continue;
    //if (abs(mom.Eta())<0.5) // Collider Mode
    if (mom.Eta()>0. && mom.Eta()<2.) // FXT Mode
      Mult++;
  }
  return Mult;
}

int getRefMultColl(vector<fourVector> _p, RVec<int> _nhits)
{
  int Mult = 0;
  for (int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto nhits = _nhits.at(i);
    if (nhits<=16) continue;
    if (abs(mom.Eta())<0.5) // Collider Mode
    //if (mom.Eta()>0. && mom.Eta()<2.) // FXT Mode
      Mult++;
  }
  return Mult;
}


RVec<float> getPt(vector<fourVector> _p)
{
  vector <float> pt_;
  for (auto& mom:_p)
    pt_.push_back(mom.Pt());
  return pt_;
}

RVec<float> getP(vector<fourVector> _p)
{
  vector <float> p_;
  for (auto& mom:_p)
    p_.push_back(mom.P());
  return p_;
}

RVec<float> getRigidity(vector<fourVector> _p,RVec<short> _q)
{
  vector <float> pq_;
  for (int i=0; i<_p.size(); ++i){
    auto mom = _p.at(i);
    auto charge = _q.at(i);
    if (charge != 0)
      pq_.push_back(mom.P()/(float)(charge));
    else
      pq_.push_back(-9999.);
  }
  return pq_;
}

RVec<float> getRigidityReverse(vector<fourVector> _p,RVec<short> _q)
{
  vector <float> pq_;
  for (int i=0; i<_p.size(); ++i){
    auto mom = _p.at(i);
    auto charge = _q.at(i);
    if (charge != 0)
      pq_.push_back(-1.*mom.P()/(float)(charge));
    else
      pq_.push_back(-9999.);
  }
  return pq_;
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

RVec<int> isRecPrim(RVec<float> _dca)
{
  vector<int> vec_prim;
  for (auto& dca:_dca){
    int prim = (dca < 1.) ? 1 : 0;
    vec_prim.push_back(prim);
  }
  return vec_prim;
}

RVec<int> isGoodTrack(vector<fourVector> _p, RVec<int> _nhits, RVec<float> _dca, RVec<bool> _tof)
{
  vector<int> vec_tracks;
  for(int i=0; i<_p.size(); ++i) {
    auto mom = _p.at(i);
    auto nhits = std::round(_nhits.at(i));
    auto dca = _dca.at(i);
    auto tof = _tof.at(i);
    if (nhits > 16 && dca < 100. && tof)
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
