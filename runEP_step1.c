// Reads the flat tree "t" produced by convertMPD.C / convertMPD_legacy.C
// with a manual TChain loop into Event objects. Only the branches connected
// in EventReader are read.
#include <ROOT/RVec.hxx>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

using namespace ROOT;
using namespace ROOT::Math;
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

template <typename T>
void connectBranch(TChain *chain, const char *name, T *address)
{
  // Activates the branch (with all its sub-branches) and sets its address
  auto branch = chain->GetBranch(name);
  if (!branch)
    throw runtime_error(Form("ERROR: No branch %s found in the input tree!", name));
  chain->SetBranchStatus(name, 1);
  if (branch->GetListOfBranches()->GetEntries() > 0)
    chain->SetBranchStatus(Form("%s.*", name), 1);
  chain->SetBranchAddress(name, address);
}

struct Track
{
  fourVector mom;
  XYZVector dca;
  int nhits;
  short charge;
  bool tofFlag;
  float tofMass2;
  float tpcDedx;
  int simPdg;       // 0 if no associated sim track
  int simMotherId;  // -1 for primary
};

struct FhcalModule
{
  int id;
  XYZVector pos;
  float energy;
};

struct TofHit
{
  XYZVector pos;
  float time;
  float dx;  // tofDX/DY/DZ: meaning not documented in this repo's converters
  float dy;
  float dz;
};

struct EmcCluster
{
  int mult;      // number of towers in the cluster; recommended cut: mult > 1
  float energy;  // non-linearity corrected energy; recommended cut: > 0.05 GeV
  XYZVector pos; // global coordinates
  float chi2;    // chi2/ndf; recommended cut for photons: chi2 < 4
  float time;    // smeared simulated time of the shower
  float dPhi;    // distance to the closest TPC track in phi
  float dZ;      // distance to the closest TPC track in z
};

class Event
{
public:
  int evId;
  XYZVector recoVtx;
  XYZVector mcVtx;
  double mcB;
  double mcRP;
  vector<Track> tracks;
  vector<FhcalModule> fhcalModules;
  vector<TofHit> tofHits;
  vector<EmcCluster> emcClusters;

  void Clear()
  {
    tracks.clear();
    fhcalModules.clear();
    tofHits.clear();
    emcClusters.clear();
  }
};

class EventReader
{
public:
  EventReader(string fileIn)
  {
    fChain = makeChain(fileIn, "t");

    // The tree is huge (~50x compression), so read only what is needed
    fChain->SetBranchStatus("*", 0);
    connectBranch(fChain, "evId", &evId);
    connectBranch(fChain, "recoPrimVtxX", &recoPrimVtxX);
    connectBranch(fChain, "recoPrimVtxY", &recoPrimVtxY);
    connectBranch(fChain, "recoPrimVtxZ", &recoPrimVtxZ);
    connectBranch(fChain, "mcVtxX", &mcVtxX);
    connectBranch(fChain, "mcVtxY", &mcVtxY);
    connectBranch(fChain, "mcVtxZ", &mcVtxZ);
    connectBranch(fChain, "mcB", &mcB);
    connectBranch(fChain, "mcRP", &mcRP);
    connectBranch(fChain, "recoGlobalMom", &recoGlobalMom);
    connectBranch(fChain, "recoGlobalDca", &recoGlobalDca);
    connectBranch(fChain, "recoGlobalNhits", &recoGlobalNhits);
    connectBranch(fChain, "recoGlobalCharge", &recoGlobalCharge);
    connectBranch(fChain, "recoGlobalTofFlag", &recoGlobalTofFlag);
    connectBranch(fChain, "recoGlobalTofMass2", &recoGlobalTofMass2);
    connectBranch(fChain, "recoGlobalTpcDedx", &recoGlobalTpcDedx);
    connectBranch(fChain, "recoGlobalSimPdg", &recoGlobalSimPdg);
    connectBranch(fChain, "recoGlobalSimMotherId", &recoGlobalSimMotherId);
    connectBranch(fChain, "fhcalModPos", &fhcalModPos);
    connectBranch(fChain, "fhcalModId", &fhcalModId);
    connectBranch(fChain, "fhcalModE", &fhcalModE);
    connectBranch(fChain, "tofPos", &tofPos);
    connectBranch(fChain, "tofTime", &tofTime);
    connectBranch(fChain, "tofDX", &tofDX);
    connectBranch(fChain, "tofDY", &tofDY);
    connectBranch(fChain, "tofDZ", &tofDZ);
    connectBranch(fChain, "emcMult", &emcMult);
    connectBranch(fChain, "emcEnergy", &emcEnergy);
    connectBranch(fChain, "emcPos", &emcPos);
    connectBranch(fChain, "emcChi2", &emcChi2);
    connectBranch(fChain, "emcTime", &emcTime);
    connectBranch(fChain, "emcDPhi", &emcDPhi);
    connectBranch(fChain, "emcDZ", &emcDZ);
  }

  ~EventReader() { delete fChain; }

  Long64_t GetEntries() { return fChain->GetEntries(); }

  bool ReadEvent(Long64_t entry, Event &event)
  {
    event.Clear();
    if (fChain->GetEntry(entry) <= 0)
      return false;

    event.evId = evId;
    event.recoVtx.SetXYZ(recoPrimVtxX, recoPrimVtxY, recoPrimVtxZ);
    event.mcVtx.SetXYZ(mcVtxX, mcVtxY, mcVtxZ);
    event.mcB = mcB;
    event.mcRP = mcRP;

    const auto nTracks = recoGlobalMom->size();
    event.tracks.reserve(nTracks);
    for (size_t i=0; i<nTracks; ++i)
      event.tracks.push_back({recoGlobalMom->at(i), recoGlobalDca->at(i),
                              recoGlobalNhits->at(i), recoGlobalCharge->at(i),
                              recoGlobalTofFlag->at(i), recoGlobalTofMass2->at(i),
                              recoGlobalTpcDedx->at(i), recoGlobalSimPdg->at(i),
                              recoGlobalSimMotherId->at(i)});

    const auto nModules = fhcalModE->size();
    event.fhcalModules.reserve(nModules);
    for (size_t i=0; i<nModules; ++i)
      event.fhcalModules.push_back({fhcalModId->at(i), fhcalModPos->at(i), fhcalModE->at(i)});

    const auto nTofHits = tofPos->size();
    event.tofHits.reserve(nTofHits);
    for (size_t i=0; i<nTofHits; ++i)
      event.tofHits.push_back({tofPos->at(i), tofTime->at(i),
                               tofDX->at(i), tofDY->at(i), tofDZ->at(i)});

    const auto nEmcClusters = emcPos->size();
    event.emcClusters.reserve(nEmcClusters);
    for (size_t i=0; i<nEmcClusters; ++i)
      event.emcClusters.push_back({emcMult->at(i), emcEnergy->at(i), emcPos->at(i),
                                   emcChi2->at(i), emcTime->at(i),
                                   emcDPhi->at(i), emcDZ->at(i)});

    return true;
  }

private:
  TChain *fChain{nullptr};

  // Branch buffers
  int evId;
  float recoPrimVtxX, recoPrimVtxY, recoPrimVtxZ;
  double mcVtxX, mcVtxY, mcVtxZ;
  double mcB, mcRP;
  vector<fourVector> *recoGlobalMom{nullptr};
  vector<XYZVector> *recoGlobalDca{nullptr};
  RVec<int> *recoGlobalNhits{nullptr};
  RVec<short> *recoGlobalCharge{nullptr};
  RVec<bool> *recoGlobalTofFlag{nullptr};
  RVec<float> *recoGlobalTofMass2{nullptr};
  RVec<float> *recoGlobalTpcDedx{nullptr};
  RVec<int> *recoGlobalSimPdg{nullptr};
  RVec<int> *recoGlobalSimMotherId{nullptr};
  vector<XYZVector> *fhcalModPos{nullptr};
  RVec<int> *fhcalModId{nullptr};
  RVec<float> *fhcalModE{nullptr};
  vector<XYZVector> *tofPos{nullptr};
  RVec<float> *tofTime{nullptr};
  RVec<float> *tofDX{nullptr};
  RVec<float> *tofDY{nullptr};
  RVec<float> *tofDZ{nullptr};
  RVec<int> *emcMult{nullptr};
  RVec<float> *emcEnergy{nullptr};
  vector<XYZVector> *emcPos{nullptr};
  RVec<float> *emcChi2{nullptr};
  RVec<float> *emcTime{nullptr};
  RVec<float> *emcDPhi{nullptr};
  RVec<float> *emcDZ{nullptr};
};

void runEP_step1(string fileIn="", string fileOut="", Long64_t nMaxEvents=-1)
{
  TStopwatch timer;
  timer.Start();

  EventReader reader(fileIn);
  Long64_t nEvents = reader.GetEntries();
  if (nMaxEvents > 0 && nMaxEvents < nEvents)
    nEvents = nMaxEvents;
  std::cout << "Number of events to read: " << nEvents << std::endl;
  if (nEvents <= 0) {
    std::cerr << "No events found in the input." << std::endl;
    return;
  }

  TFile fOut(fileOut.c_str(), "recreate");

  Event event;
  for (Long64_t iEv=0; iEv<nEvents; ++iEv) {
    if (iEv % 100 == 0) cout << "\r" << iEv << flush;
    if (!reader.ReadEvent(iEv, event))
      continue;

    for (const auto &track : event.tracks) {
      // TODO: track loop
    }

    for (const auto &module : event.fhcalModules) {
      // TODO: FHCal module loop
    }

    for (const auto &hit : event.tofHits) {
      // TODO: TOF hit loop
    }

    for (const auto &cluster : event.emcClusters) {
      // TODO: EMC cluster loop
    }
  }
  cout << endl;

  fOut.Close();

  timer.Stop();
  timer.Print();
}
