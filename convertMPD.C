using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using ROOT::VecOps::Map;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

MpdFieldMap* magField{nullptr};

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

bool isGoodSimTrack(const MpdMCTrack &track)
{
  return true;
}

double getMass(int pdg)
{
  if(pdg<1000000) 
    return TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  else 
    return 0.931*(pdg/10%1000);
}

int getCharge(int pdg)
{
  if(pdg<1000000)
    return TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
  else 
    return pdg/10000%1000;
}

vector<fourVector> simMomentum(RVec<MpdMCTrack> tracks)
{
  vector<fourVector> momenta;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    TVector3 mom;
    track.GetMomentum(mom);
    double mass=getMass(track.GetPdgCode());
    momenta.push_back({mom.Pt(),mom.Eta(),mom.Phi(),sqrt(mass*mass+mom.Mag2())});
  }
  return momenta;
}

vector<XYZTVector> simPosStart(const RVec<MpdMCTrack> tracks)
{
  vector<XYZTVector> pos;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    pos.push_back({track.GetStartX(),track.GetStartY(),track.GetStartZ(),track.GetStartT()});
  }
  return pos;
}

RVec<int> simMotherId(const RVec<MpdMCTrack> tracks)
{
  vector<int> mothId;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    mothId.push_back(track.GetMotherId());
  }
  return mothId;
}

RVec<int> simPdg(const RVec<MpdMCTrack> tracks)
{
  vector<int> pdg;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    pdg.push_back(track.GetPdgCode());
  }
  return pdg;
}

RVec<short> simCharge (const RVec<int> pdg)
{
  vector<short> ch;
  for (auto &p:pdg)
    ch.push_back(getCharge(p));
  return ch;
}

vector<float> trackP(const RVec<float> &pt, const RVec<float> &theta)
try {
  vector<float> momenta;
  for (int i=0; i<pt.size(); i++) {
    float tt = TMath::Tan(theta.at(i));
    float pz = tt ? TMath::Abs(pt.at(i)) / tt : 9999.;
    momenta.push_back(sqrt(pt.at(i)*pt.at(i) + pz*pz));
  }
  return momenta;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<fourVector> trackMomentum(const RVec<float> &pt, const RVec<float> &theta, const RVec<float> &phi)
try {
  vector<fourVector> momenta;
  for (int i=0; i<pt.size(); i++) {
    momenta.push_back({abs(pt.at(i)), -TMath::Log(TMath::Tan(0.5 * theta.at(i))), phi.at(0), 0});
  }
  return momenta;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<short> recCharge(const RVec<float> &pt)
try {
  vector<short> charge;
  for (auto trPt:pt) {
    int q = (trPt > 0) ? -1. : 1.;
    charge.push_back(q);
  }
  return charge;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> recDca(const RVec<float> &dcaX, const RVec<float> &dcaY, const RVec<float> &dcaZ)
try {
  vector<XYZVector> dca;
  for (int i=0; i<dcaX.size(); i++) {
    dca.push_back({dcaX.at(i),dcaY.at(i),dcaZ.at(i)});
  }
  return dca;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<vector<float>> trGlobalParam(const RVec<float> &pointX, const RVec<float> &pointY, const RVec<float> &pointZ)
try {
  vector<vector<float>> parameters;
  for (int i=0; i<pointX.size(); i++) {
    parameters.emplace_back();
    parameters.back().push_back( pointX.at(i) );
    parameters.back().push_back( pointY.at(i) );
    parameters.back().push_back( pointZ.at(i) );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> recSimIndex(const RVec<int> &recId, const RVec<MpdMCTrack> simTracks)
try {
  vector<int> newIndex;
  int shift=0;
  int nSimTracks = simTracks.size();
  for (int i=0;i<nSimTracks;i++) {
    if (!isGoodSimTrack(simTracks.at(i)))
    {
      shift++;
      newIndex.push_back(-1);
    }
    else
      newIndex.push_back(i-shift);
  }
  vector<int> simIndex;
  for (auto trId:recId) {
    int oldIndex=trId;
    if (oldIndex<0 || oldIndex>=nSimTracks)
      simIndex.push_back(-1);
    else
      simIndex.push_back(newIndex.at(oldIndex));
  }
  return simIndex;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> trSimPdg(const RVec<int> sim_index, const RVec<int> sim_pdg)
try {
  std::vector<int> pdg;
  for( auto idx : sim_index ) {
    if( idx < 0 ) {
      pdg.push_back(-1);
      continue;
    }
    if( idx > sim_pdg.size() ) {
      pdg.push_back(-1);
      continue;
    }
    pdg.push_back(sim_pdg.at(idx));
  }
  return pdg;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> trSimMotherId(const RVec<int> sim_index, const RVec<int> sim_motherId)
try {
  std::vector<int> motherId;
  for( auto idx : sim_index ) {
    if( idx < 0 ) {
      motherId.push_back(-1);
      continue;
    }
    if( idx > sim_motherId.size() ) {
      motherId.push_back(-1);
      continue;
    }
    motherId.push_back(sim_motherId.at(idx));
  }
  return motherId;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

XYZVector GetFHCalPos(int iModule)
{
  const int Nmodules = 45;
  int xAxisSwitch = (iModule < Nmodules) ? 1 : -1;
  int module = (iModule < Nmodules) ? iModule : iModule - Nmodules;
  float x, y, z;
  z = (iModule < Nmodules) ? 320. : -320.;
  if (module >= 0 && module <= 4)
  {
    y = 45.;
    x = -1. * (module - 2) * 15.;
  }
  else if ((module >= 5) && (module <= 39))
  {
    y = (3 - (module + 2) / 7) * 15.;
    x = (3 - (module + 2) % 7) * 15.;
  }
  else if ((module >= 40) && (module <= 44))
  {
    y = -45.;
    x = -1. * (module - 42) * 15.;
  }
  XYZVector vModPos(x * xAxisSwitch, y, z);

  return vModPos;
}

vector<XYZVector> modulePos ()
try {
  bool verbose = false;
  const int nModules = 90;
  vector<XYZVector> modulePosVector(nModules,{0.,0.,0.});
  for (int i=0; i<nModules; i++) {
    modulePosVector.at(i) = GetFHCalPos(i);
  }
  if (verbose)
  {
    printf("%d module positions:\n", nModules);
    for(int i=0;i<nModules;i++)
      printf("%d: (%f, %f, %f)\n", i, modulePosVector.at(i).x(), modulePosVector.at(i).y(), modulePosVector.at(i).z());
  }
  return modulePosVector;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> fhcalModE(const RVec<MpdZdcDigi> fhcalHits)
try {
  const int nModules = 90;
  vector<float> fhcalModEnergy(nModules, 0.);
  for (auto hit:fhcalHits) {
    int det_id = hit.GetDetectorID();
    int mod_id = hit.GetModuleID()-1;
    int i_module = mod_id + (nModules/2) * (det_id-1);
    fhcalModEnergy.at(i_module) += hit.GetELoss();
  }
  return fhcalModEnergy;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> moduleId (const vector<XYZVector> modulePos)
try {
  vector <int> moduleIds;
  for (int i=0;i<modulePos.size();i++)
    moduleIds.push_back(i+1);
  return moduleIds;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<bool> hasHitFhcal (const RVec<MpdMCTrack> particles)
{
  vector<bool> hasHit;
  for(auto &part:particles)
    hasHit.push_back(part.GetNPoints(kZDC)>0);
  return hasHit;
}

void convertMPD(string inDst="", string fileOut="")
{
  TStopwatch timer;
  timer.Start();
  
  TChain *chainRec=makeChain(inDst, "mpdsim");
  ROOT::RDataFrame d(*chainRec);

  int nEvents = chainRec->GetEntries();

  auto fhcalModPos = modulePos();

  auto dd = d
    .Define("evtId","MCEventHeader.fEventId")
    .Define("recoPrimVtxX","MPDEvent.PrimaryVerticesX")
    .Define("recoPrimVtxY","MPDEvent.PrimaryVerticesY")
    .Define("recoPrimVtxZ","MPDEvent.PrimaryVerticesZ")
    .Define("recoPrimVtxChi2", "MPDEvent.PrimaryVerticesChi2")
    .Define("recoVtxX","Vertex.fX")
    .Define("recoVtxY","Vertex.fY")
    .Define("recoVtxZ","Vertex.fZ")
    .Define("recoVtxChi2", "Vertex.fChi2")
    .Define("recoVtxNDF", "Vertex.fNDF")
    .Define("recoVtxNtracks", "Vertex.fNTracks")
    .Define("mcVtxX","MCEventHeader.fX")
    .Define("mcVtxY","MCEventHeader.fY")
    .Define("mcVtxZ","MCEventHeader.fZ")
    .Define("mcB", "MCEventHeader.fB")
    .Define("mcRP", "MCEventHeader.fRotZ")
    .Define("recoGlobalMom", trackMomentum, {"MPDEvent.fGlobalTracks.fPt", "MPDEvent.fGlobalTracks.fTheta", "MPDEvent.fGlobalTracks.fPhi"})
    .Define("recoGlobalNhits", "MPDEvent.fGlobalTracks.fNofHits")
    .Define("recoGlobalNhitsPoss", "MPDEvent.fGlobalTracks.fNofHitsPossTpc")
    .Define("recoGlobalNhitsFit", "MPDEvent.fGlobalTracks.fNofHitsFitTpc")
    .Define("recoGlobalChi2", "MPDEvent.fGlobalTracks.fChi2")
    .Define("recoGlobalP", trackP, {"MPDEvent.fGlobalTracks.fPt", "MPDEvent.fGlobalTracks.fTheta"})
    .Define("recoGlobalTofFlag", "MPDEvent.fGlobalTracks.fTofFlag")
    .Define("recoGlobalCharge",recCharge,{"MPDEvent.fGlobalTracks.fPt"})
    .Define("recoGlobalDca", recDca, {"MPDEvent.fGlobalTracks.fDCAX", "MPDEvent.fGlobalTracks.fDCAY", "MPDEvent.fGlobalTracks.fDCAZ"})
    .Define("recoGlobalTofMass2", "MPDEvent.fGlobalTracks.fTofMass2")
    //.Define("recoGlobalParamFirst", trGlobalParam, {"MPDEvent.fGlobalTracks.fFirstPointX", "MPDEvent.fGlobalTracks.fFirstPointY", "MPDEvent.fGlobalTracks.fFirstPointZ"})
    .Define("simMom", simMomentum, {"MCTrack"})
    .Define("simPosStart", simPosStart, {"MCTrack"})
    .Define("simMotherId", simMotherId, {"MCTrack"})
    .Define("simPdg", simPdg, {"MCTrack"})
    .Define("simCharge", simCharge, {"simPdg"})
    .Define("simHasHitFHCal", hasHitFhcal, {"MCTrack"})
    .Define("recoGlobalSimIndex", recSimIndex, {"MPDEvent.fGlobalTracks.fID", "MCTrack"})
    .Define("recoGlobalSimPdg", trSimPdg, {"recoGlobalSimIndex", "simPdg"})
    .Define("recoGlobalSimMotherId", trSimMotherId, {"recoGlobalSimIndex", "simMotherId"})
    .Define("fhcalModPos", [fhcalModPos](){return fhcalModPos; })
    .Define("fhcalModId", moduleId, {"fhcalModPos"})
    .Define("fhcalModE",fhcalModE,{"ZdcDigi"})
  ;
  dd.Foreach([](uint evtId){if (evtId % 100 == 0) cout << "\n" << evtId;}, {"evtId"}); // progress display 
  cout << endl;

  vector<string> definedNames;
  vector<string> toExclude={""};
  for (auto& definedName:dd.GetDefinedColumnNames())
  {
    bool exclude=false;
    for (auto &nameToExclude:toExclude)
      if (definedName==nameToExclude)
        exclude=true;
    if (!exclude)
      definedNames.push_back(definedName);
  }
  dd.Snapshot("t", fileOut, definedNames);
  
  timer.Stop();
  timer.Print();
}
