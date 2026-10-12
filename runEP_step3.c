// Reads the flat tree "t" produced by convertMPD.C / convertMPD_legacy.C
// with a manual TChain loop into Event objects. Only the branches connected
// in EventReader are read.
#include <ROOT/RVec.hxx>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

using namespace ROOT;
using namespace ROOT::Math;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

std::vector<int> x1_modules = { 0, 1, 5, 6, 7, 12, 13, 14 };
std::vector<int> x2_modules = { 26, 27, 28, 33, 34, 35, 40, 41 };
std::vector<int> x3_modules = { 30, 31, 32, 37, 38, 39, 43, 44 };
std::vector<int> x4_modules = { 3, 4, 9, 10, 11, 16, 17, 18 };

const std::pair<double, double> tof_pos_cut = {-200., -100.};
const std::pair<double, double> emc_pos_cut = {-200., -100.};


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

namespace ep{
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

  bool ReadEvent(Long64_t entry, ep::Event &event)
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

Float_t GetCentMult(Int_t mult)
// Returns the centrality (middle of the class, in percent) for a given multiplicity, or -1
{
  Float_t minCentPercent [20] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95 };
  Float_t maxCentPercent [20] = { 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100 };
  Int_t minMult [20] = { 263, 233, 210, 190, 172, 157, 144, 132, 122, 113, 105, 98, 91, 86, 81, 76, 71, 66, 59, 1 };
  Int_t maxMult [20] = { 526, 263, 233, 210, 190, 172, 157, 144, 132, 122, 113, 105, 98, 91, 86, 81, 76, 71, 66, 59 };
        Int_t centBin = -1;
        for (Int_t i = 0; i < 20; i++)
        {
                if (mult >= minMult[i] && mult < maxMult[i])
                        centBin = i;
        }
        if (centBin == -1) return -1.;
        return (maxCentPercent[centBin] - minCentPercent[centBin]) / 2. + minCentPercent[centBin];
}

  }

void runEP_step3(string fileIn="", string fileOut="", string fileInStep1="", string fileInStep2="", Long64_t nMaxEvents=-1)
{
  TStopwatch timer;
  timer.Start();

  ep::EventReader reader(fileIn);
  Long64_t nEvents = reader.GetEntries();
  if (nMaxEvents > 0 && nMaxEvents < nEvents)
    nEvents = nMaxEvents;
  std::cout << "Number of events to read: " << nEvents << std::endl;
  if (nEvents <= 0) {
    std::cerr << "No events found in the input." << std::endl;
    return;
  }
  TFile *fi_step1 = new TFile(fileInStep1.c_str(), "read");
  if (!fi_step1) {
    std::cerr << "Cannot open file from step1!" << std::endl;
    return;
  }
  TFile *fi_step2 = new TFile(fileInStep2.c_str(), "read");
  if (!fi_step2) {
    std::cerr << "Cannot open file from step2!" << std::endl;
    return;
  }

  TFile fOut(fileOut.c_str(), "recreate");

  TH1D *h_cent = new TH1D("h_cent", "Centrality diestribution", 20, 0., 100.);

  TH2D *h2_fhcal_xy_x1 = new TH2D("h2_fhcal_xy_x1","(XY) positions of FHCal (X1 sub-event)", 7, -52.5, 52.5, 7, -52.5, 52.5);
  TH2D *h2_fhcal_xy_x2 = new TH2D("h2_fhcal_xy_x2","(XY) positions of FHCal (X2 sub-event)", 7, -52.5, 52.5, 7, -52.5, 52.5);
  TH2D *h2_fhcal_xy_x3 = new TH2D("h2_fhcal_xy_x3","(XY) positions of FHCal (X3 sub-event)", 7, -52.5, 52.5, 7, -52.5, 52.5);
  TH2D *h2_fhcal_xy_x4 = new TH2D("h2_fhcal_xy_x4","(XY) positions of FHCal (X4 sub-event)", 7, -52.5, 52.5, 7, -52.5, 52.5);

  TH1D *h1_tof_posz = new TH1D("h1_tof_posz", "Z-position of TOF hits", 600, -300., 300.);
  TH1D *h1_tof_posz_after_cut = new TH1D("h1_tof_posz_after_cut", "Z-position of TOF hits after Z-pos cut", 600, -300., 300.);
  TH1D *h1_emc_posz = new TH1D("h1_emc_posz", "Z-position of EMC hits", 600, -300., 300.);
  TH1D *h1_emc_posz_after_cut = new TH1D("h1_emc_posz_after_cut", "Z-position of EMC hits after Z-pos cut", 600, -300., 300.);

  // Profiles for recentering
  TProfile *p_avg_qx_fhcal_x1 = (TProfile*)fi_step1->Get("p_avg_qx_fhcal_x1");
  TProfile *p_avg_qx_fhcal_x2 = (TProfile*)fi_step1->Get("p_avg_qx_fhcal_x2");
  TProfile *p_avg_qx_fhcal_x3 = (TProfile*)fi_step1->Get("p_avg_qx_fhcal_x3");
  TProfile *p_avg_qx_fhcal_x4 = (TProfile*)fi_step1->Get("p_avg_qx_fhcal_x4");
  TProfile *p_avg_qx_tof = (TProfile*)fi_step1->Get("p_avg_qx_tof");
  TProfile *p_avg_qx_emc = (TProfile*)fi_step1->Get("p_avg_qx_emc");
  TProfile *p_avg_qy_fhcal_x1 = (TProfile*)fi_step1->Get("p_avg_qy_fhcal_x1");
  TProfile *p_avg_qy_fhcal_x2 = (TProfile*)fi_step1->Get("p_avg_qy_fhcal_x2");
  TProfile *p_avg_qy_fhcal_x3 = (TProfile*)fi_step1->Get("p_avg_qy_fhcal_x3");
  TProfile *p_avg_qy_fhcal_x4 = (TProfile*)fi_step1->Get("p_avg_qy_fhcal_x4");
  TProfile *p_avg_qy_tof = (TProfile*)fi_step1->Get("p_avg_qy_tof");
  TProfile *p_avg_qy_emc = (TProfile*)fi_step1->Get("p_avg_qy_emc");

  // Profiles for shift correction (flattening)
  TProfile *p_avg_cos1_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_cos1_fhcal_x1");
  TProfile *p_avg_cos2_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_cos2_fhcal_x1");
  TProfile *p_avg_cos3_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_cos3_fhcal_x1");
  TProfile *p_avg_cos4_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_cos4_fhcal_x1");
  TProfile *p_avg_cos5_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_cos5_fhcal_x1");
  TProfile *p_avg_cos6_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_cos6_fhcal_x1");
  TProfile *p_avg_cos7_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_cos7_fhcal_x1");
  TProfile *p_avg_cos8_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_cos8_fhcal_x1");
  TProfile *p_avg_sin1_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_sin1_fhcal_x1");
  TProfile *p_avg_sin2_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_sin2_fhcal_x1");
  TProfile *p_avg_sin3_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_sin3_fhcal_x1");
  TProfile *p_avg_sin4_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_sin4_fhcal_x1");
  TProfile *p_avg_sin5_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_sin5_fhcal_x1");
  TProfile *p_avg_sin6_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_sin6_fhcal_x1");
  TProfile *p_avg_sin7_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_sin7_fhcal_x1");
  TProfile *p_avg_sin8_fhcal_x1 = (TProfile*)fi_step2->Get("p_avg_sin8_fhcal_x1");

  TProfile *p_avg_cos1_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_cos1_fhcal_x2");
  TProfile *p_avg_cos2_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_cos2_fhcal_x2");
  TProfile *p_avg_cos3_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_cos3_fhcal_x2");
  TProfile *p_avg_cos4_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_cos4_fhcal_x2");
  TProfile *p_avg_cos5_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_cos5_fhcal_x2");
  TProfile *p_avg_cos6_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_cos6_fhcal_x2");
  TProfile *p_avg_cos7_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_cos7_fhcal_x2");
  TProfile *p_avg_cos8_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_cos8_fhcal_x2");
  TProfile *p_avg_sin1_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_sin1_fhcal_x2");
  TProfile *p_avg_sin2_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_sin2_fhcal_x2");
  TProfile *p_avg_sin3_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_sin3_fhcal_x2");
  TProfile *p_avg_sin4_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_sin4_fhcal_x2");
  TProfile *p_avg_sin5_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_sin5_fhcal_x2");
  TProfile *p_avg_sin6_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_sin6_fhcal_x2");
  TProfile *p_avg_sin7_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_sin7_fhcal_x2");
  TProfile *p_avg_sin8_fhcal_x2 = (TProfile*)fi_step2->Get("p_avg_sin8_fhcal_x2");

  TProfile *p_avg_cos1_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_cos1_fhcal_x3");
  TProfile *p_avg_cos2_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_cos2_fhcal_x3");
  TProfile *p_avg_cos3_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_cos3_fhcal_x3");
  TProfile *p_avg_cos4_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_cos4_fhcal_x3");
  TProfile *p_avg_cos5_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_cos5_fhcal_x3");
  TProfile *p_avg_cos6_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_cos6_fhcal_x3");
  TProfile *p_avg_cos7_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_cos7_fhcal_x3");
  TProfile *p_avg_cos8_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_cos8_fhcal_x3");
  TProfile *p_avg_sin1_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_sin1_fhcal_x3");
  TProfile *p_avg_sin2_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_sin2_fhcal_x3");
  TProfile *p_avg_sin3_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_sin3_fhcal_x3");
  TProfile *p_avg_sin4_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_sin4_fhcal_x3");
  TProfile *p_avg_sin5_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_sin5_fhcal_x3");
  TProfile *p_avg_sin6_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_sin6_fhcal_x3");
  TProfile *p_avg_sin7_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_sin7_fhcal_x3");
  TProfile *p_avg_sin8_fhcal_x3 = (TProfile*)fi_step2->Get("p_avg_sin8_fhcal_x3");

  TProfile *p_avg_cos1_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_cos1_fhcal_x4");
  TProfile *p_avg_cos2_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_cos2_fhcal_x4");
  TProfile *p_avg_cos3_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_cos3_fhcal_x4");
  TProfile *p_avg_cos4_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_cos4_fhcal_x4");
  TProfile *p_avg_cos5_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_cos5_fhcal_x4");
  TProfile *p_avg_cos6_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_cos6_fhcal_x4");
  TProfile *p_avg_cos7_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_cos7_fhcal_x4");
  TProfile *p_avg_cos8_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_cos8_fhcal_x4");
  TProfile *p_avg_sin1_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_sin1_fhcal_x4");
  TProfile *p_avg_sin2_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_sin2_fhcal_x4");
  TProfile *p_avg_sin3_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_sin3_fhcal_x4");
  TProfile *p_avg_sin4_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_sin4_fhcal_x4");
  TProfile *p_avg_sin5_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_sin5_fhcal_x4");
  TProfile *p_avg_sin6_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_sin6_fhcal_x4");
  TProfile *p_avg_sin7_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_sin7_fhcal_x4");
  TProfile *p_avg_sin8_fhcal_x4 = (TProfile*)fi_step2->Get("p_avg_sin8_fhcal_x4");

  TProfile *p_avg_cos1_tof = (TProfile*)fi_step2->Get("p_avg_cos1_tof");
  TProfile *p_avg_cos2_tof = (TProfile*)fi_step2->Get("p_avg_cos2_tof");
  TProfile *p_avg_cos3_tof = (TProfile*)fi_step2->Get("p_avg_cos3_tof");
  TProfile *p_avg_cos4_tof = (TProfile*)fi_step2->Get("p_avg_cos4_tof");
  TProfile *p_avg_cos5_tof = (TProfile*)fi_step2->Get("p_avg_cos5_tof");
  TProfile *p_avg_cos6_tof = (TProfile*)fi_step2->Get("p_avg_cos6_tof");
  TProfile *p_avg_cos7_tof = (TProfile*)fi_step2->Get("p_avg_cos7_tof");
  TProfile *p_avg_cos8_tof = (TProfile*)fi_step2->Get("p_avg_cos8_tof");
  TProfile *p_avg_sin1_tof = (TProfile*)fi_step2->Get("p_avg_sin1_tof");
  TProfile *p_avg_sin2_tof = (TProfile*)fi_step2->Get("p_avg_sin2_tof");
  TProfile *p_avg_sin3_tof = (TProfile*)fi_step2->Get("p_avg_sin3_tof");
  TProfile *p_avg_sin4_tof = (TProfile*)fi_step2->Get("p_avg_sin4_tof");
  TProfile *p_avg_sin5_tof = (TProfile*)fi_step2->Get("p_avg_sin5_tof");
  TProfile *p_avg_sin6_tof = (TProfile*)fi_step2->Get("p_avg_sin6_tof");
  TProfile *p_avg_sin7_tof = (TProfile*)fi_step2->Get("p_avg_sin7_tof");
  TProfile *p_avg_sin8_tof = (TProfile*)fi_step2->Get("p_avg_sin8_tof");

  TProfile *p_avg_cos1_emc = (TProfile*)fi_step2->Get("p_avg_cos1_emc");
  TProfile *p_avg_cos2_emc = (TProfile*)fi_step2->Get("p_avg_cos2_emc");
  TProfile *p_avg_cos3_emc = (TProfile*)fi_step2->Get("p_avg_cos3_emc");
  TProfile *p_avg_cos4_emc = (TProfile*)fi_step2->Get("p_avg_cos4_emc");
  TProfile *p_avg_cos5_emc = (TProfile*)fi_step2->Get("p_avg_cos5_emc");
  TProfile *p_avg_cos6_emc = (TProfile*)fi_step2->Get("p_avg_cos6_emc");
  TProfile *p_avg_cos7_emc = (TProfile*)fi_step2->Get("p_avg_cos7_emc");
  TProfile *p_avg_cos8_emc = (TProfile*)fi_step2->Get("p_avg_cos8_emc");
  TProfile *p_avg_sin1_emc = (TProfile*)fi_step2->Get("p_avg_sin1_emc");
  TProfile *p_avg_sin2_emc = (TProfile*)fi_step2->Get("p_avg_sin2_emc");
  TProfile *p_avg_sin3_emc = (TProfile*)fi_step2->Get("p_avg_sin3_emc");
  TProfile *p_avg_sin4_emc = (TProfile*)fi_step2->Get("p_avg_sin4_emc");
  TProfile *p_avg_sin5_emc = (TProfile*)fi_step2->Get("p_avg_sin5_emc");
  TProfile *p_avg_sin6_emc = (TProfile*)fi_step2->Get("p_avg_sin6_emc");
  TProfile *p_avg_sin7_emc = (TProfile*)fi_step2->Get("p_avg_sin7_emc");
  TProfile *p_avg_sin8_emc = (TProfile*)fi_step2->Get("p_avg_sin8_emc");

  // Profiles for <QQ> and <cos(1-2)>, <sin(1-2)> correlations
  TProfile *p_cos_fhcal_x1_fhcal_x2 = new TProfile("p_cos_fhcal_x1_fhcal_x2", "<cos(PsiA-PsiB)> for FHCal X1 & FHCal X2 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x1_fhcal_x3 = new TProfile("p_cos_fhcal_x1_fhcal_x3", "<cos(PsiA-PsiB)> for FHCal X1 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x1_fhcal_x4 = new TProfile("p_cos_fhcal_x1_fhcal_x4", "<cos(PsiA-PsiB)> for FHCal X1 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x1_tof = new TProfile("p_cos_fhcal_x1_tof", "<cos(PsiA-PsiB)> for FHCal X1 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x1_emc = new TProfile("p_cos_fhcal_x1_emc", "<cos(PsiA-PsiB)> for FHCal X1 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x2_fhcal_x3 = new TProfile("p_cos_fhcal_x2_fhcal_x3", "<cos(PsiA-PsiB)> for FHCal X2 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x2_fhcal_x4 = new TProfile("p_cos_fhcal_x2_fhcal_x4", "<cos(PsiA-PsiB)> for FHCal X2 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x2_tof = new TProfile("p_cos_fhcal_x2_tof", "<cos(PsiA-PsiB)> for FHCal X2 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x2_emc = new TProfile("p_cos_fhcal_x2_emc", "<cos(PsiA-PsiB)> for FHCal X2 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x3_fhcal_x4 = new TProfile("p_cos_fhcal_x3_fhcal_x4", "<cos(PsiA-PsiB)> for FHCal X3 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x3_tof = new TProfile("p_cos_fhcal_x3_tof", "<cos(PsiA-PsiB)> for FHCal X3 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x3_emc = new TProfile("p_cos_fhcal_x3_emc", "<cos(PsiA-PsiB)> for FHCal X3 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x4_tof = new TProfile("p_cos_fhcal_x4_tof", "<cos(PsiA-PsiB)> for FHCal X4 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_cos_fhcal_x4_emc = new TProfile("p_cos_fhcal_x4_emc", "<cos(PsiA-PsiB)> for FHCal X4 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x1_fhcal_x2 = new TProfile("p_sin_fhcal_x1_fhcal_x2", "<sin(PsiA-PsiB)> for FHCal X1 & FHCal X2 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x1_fhcal_x3 = new TProfile("p_sin_fhcal_x1_fhcal_x3", "<sin(PsiA-PsiB)> for FHCal X1 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x1_fhcal_x4 = new TProfile("p_sin_fhcal_x1_fhcal_x4", "<sin(PsiA-PsiB)> for FHCal X1 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x1_tof = new TProfile("p_sin_fhcal_x1_tof", "<sin(PsiA-PsiB)> for FHCal X1 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x1_emc = new TProfile("p_sin_fhcal_x1_emc", "<sin(PsiA-PsiB)> for FHCal X1 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x2_fhcal_x3 = new TProfile("p_sin_fhcal_x2_fhcal_x3", "<sin(PsiA-PsiB)> for FHCal X2 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x2_fhcal_x4 = new TProfile("p_sin_fhcal_x2_fhcal_x4", "<sin(PsiA-PsiB)> for FHCal X2 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x2_tof = new TProfile("p_sin_fhcal_x2_tof", "<sin(PsiA-PsiB)> for FHCal X2 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x2_emc = new TProfile("p_sin_fhcal_x2_emc", "<sin(PsiA-PsiB)> for FHCal X2 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x3_fhcal_x4 = new TProfile("p_sin_fhcal_x3_fhcal_x4", "<sin(PsiA-PsiB)> for FHCal X3 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x3_tof = new TProfile("p_sin_fhcal_x3_tof", "<sin(PsiA-PsiB)> for FHCal X3 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x3_emc = new TProfile("p_sin_fhcal_x3_emc", "<sin(PsiA-PsiB)> for FHCal X3 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x4_tof = new TProfile("p_sin_fhcal_x4_tof", "<sin(PsiA-PsiB)> for FHCal X4 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_sin_fhcal_x4_emc = new TProfile("p_sin_fhcal_x4_emc", "<sin(PsiA-PsiB)> for FHCal X4 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x1_fhcal_x2 = new TProfile("p_xx_fhcal_x1_fhcal_x2", "<cos(PsiA)cos(PsiB)> for FHCal X1 & FHCal X2 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x1_fhcal_x3 = new TProfile("p_xx_fhcal_x1_fhcal_x3", "<cos(PsiA)cos(PsiB)> for FHCal X1 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x1_fhcal_x4 = new TProfile("p_xx_fhcal_x1_fhcal_x4", "<cos(PsiA)cos(PsiB)> for FHCal X1 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x1_tof = new TProfile("p_xx_fhcal_x1_tof", "<cos(PsiA)cos(PsiB)> for FHCal X1 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x1_emc = new TProfile("p_xx_fhcal_x1_emc", "<cos(PsiA)cos(PsiB)> for FHCal X1 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x2_fhcal_x3 = new TProfile("p_xx_fhcal_x2_fhcal_x3", "<cos(PsiA)cos(PsiB)> for FHCal X2 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x2_fhcal_x4 = new TProfile("p_xx_fhcal_x2_fhcal_x4", "<cos(PsiA)cos(PsiB)> for FHCal X2 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x2_tof = new TProfile("p_xx_fhcal_x2_tof", "<cos(PsiA)cos(PsiB)> for FHCal X2 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x2_emc = new TProfile("p_xx_fhcal_x2_emc", "<cos(PsiA)cos(PsiB)> for FHCal X2 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x3_fhcal_x4 = new TProfile("p_xx_fhcal_x3_fhcal_x4", "<cos(PsiA)cos(PsiB)> for FHCal X3 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x3_tof = new TProfile("p_xx_fhcal_x3_tof", "<cos(PsiA)cos(PsiB)> for FHCal X3 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x3_emc = new TProfile("p_xx_fhcal_x3_emc", "<cos(PsiA)cos(PsiB)> for FHCal X3 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x4_tof = new TProfile("p_xx_fhcal_x4_tof", "<cos(PsiA)cos(PsiB)> for FHCal X4 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xx_fhcal_x4_emc = new TProfile("p_xx_fhcal_x4_emc", "<cos(PsiA)cos(PsiB)> for FHCal X4 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x1_fhcal_x2 = new TProfile("p_yy_fhcal_x1_fhcal_x2", "<sin(PsiA)sin(PsiB)> for FHCal X1 & FHCal X2 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x1_fhcal_x3 = new TProfile("p_yy_fhcal_x1_fhcal_x3", "<sin(PsiA)sin(PsiB)> for FHCal X1 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x1_fhcal_x4 = new TProfile("p_yy_fhcal_x1_fhcal_x4", "<sin(PsiA)sin(PsiB)> for FHCal X1 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x1_tof = new TProfile("p_yy_fhcal_x1_tof", "<sin(PsiA)sin(PsiB)> for FHCal X1 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x1_emc = new TProfile("p_yy_fhcal_x1_emc", "<sin(PsiA)sin(PsiB)> for FHCal X1 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x2_fhcal_x3 = new TProfile("p_yy_fhcal_x2_fhcal_x3", "<sin(PsiA)sin(PsiB)> for FHCal X2 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x2_fhcal_x4 = new TProfile("p_yy_fhcal_x2_fhcal_x4", "<sin(PsiA)sin(PsiB)> for FHCal X2 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x2_tof = new TProfile("p_yy_fhcal_x2_tof", "<sin(PsiA)sin(PsiB)> for FHCal X2 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x2_emc = new TProfile("p_yy_fhcal_x2_emc", "<sin(PsiA)sin(PsiB)> for FHCal X2 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x3_fhcal_x4 = new TProfile("p_yy_fhcal_x3_fhcal_x4", "<sin(PsiA)sin(PsiB)> for FHCal X3 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x3_tof = new TProfile("p_yy_fhcal_x3_tof", "<sin(PsiA)sin(PsiB)> for FHCal X3 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x3_emc = new TProfile("p_yy_fhcal_x3_emc", "<sin(PsiA)sin(PsiB)> for FHCal X3 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x4_tof = new TProfile("p_yy_fhcal_x4_tof", "<sin(PsiA)sin(PsiB)> for FHCal X4 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yy_fhcal_x4_emc = new TProfile("p_yy_fhcal_x4_emc", "<sin(PsiA)sin(PsiB)> for FHCal X4 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x1_fhcal_x2 = new TProfile("p_xy_fhcal_x1_fhcal_x2", "<cos(PsiA)sin(PsiB)> for FHCal X1 & FHCal X2 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x1_fhcal_x3 = new TProfile("p_xy_fhcal_x1_fhcal_x3", "<cos(PsiA)sin(PsiB)> for FHCal X1 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x1_fhcal_x4 = new TProfile("p_xy_fhcal_x1_fhcal_x4", "<cos(PsiA)sin(PsiB)> for FHCal X1 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x1_tof = new TProfile("p_xy_fhcal_x1_tof", "<cos(PsiA)sin(PsiB)> for FHCal X1 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x1_emc = new TProfile("p_xy_fhcal_x1_emc", "<cos(PsiA)sin(PsiB)> for FHCal X1 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x2_fhcal_x3 = new TProfile("p_xy_fhcal_x2_fhcal_x3", "<cos(PsiA)sin(PsiB)> for FHCal X2 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x2_fhcal_x4 = new TProfile("p_xy_fhcal_x2_fhcal_x4", "<cos(PsiA)sin(PsiB)> for FHCal X2 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x2_tof = new TProfile("p_xy_fhcal_x2_tof", "<cos(PsiA)sin(PsiB)> for FHCal X2 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x2_emc = new TProfile("p_xy_fhcal_x2_emc", "<cos(PsiA)sin(PsiB)> for FHCal X2 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x3_fhcal_x4 = new TProfile("p_xy_fhcal_x3_fhcal_x4", "<cos(PsiA)sin(PsiB)> for FHCal X3 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x3_tof = new TProfile("p_xy_fhcal_x3_tof", "<cos(PsiA)sin(PsiB)> for FHCal X3 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x3_emc = new TProfile("p_xy_fhcal_x3_emc", "<cos(PsiA)sin(PsiB)> for FHCal X3 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x4_tof = new TProfile("p_xy_fhcal_x4_tof", "<cos(PsiA)sin(PsiB)> for FHCal X4 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_xy_fhcal_x4_emc = new TProfile("p_xy_fhcal_x4_emc", "<cos(PsiA)sin(PsiB)> for FHCal X4 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x1_fhcal_x2 = new TProfile("p_yx_fhcal_x1_fhcal_x2", "<sin(PsiA)cos(PsiB)> for FHCal X1 & FHCal X2 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x1_fhcal_x3 = new TProfile("p_yx_fhcal_x1_fhcal_x3", "<sin(PsiA)cos(PsiB)> for FHCal X1 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x1_fhcal_x4 = new TProfile("p_yx_fhcal_x1_fhcal_x4", "<sin(PsiA)cos(PsiB)> for FHCal X1 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x1_tof = new TProfile("p_yx_fhcal_x1_tof", "<sin(PsiA)cos(PsiB)> for FHCal X1 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x1_emc = new TProfile("p_yx_fhcal_x1_emc", "<sin(PsiA)cos(PsiB)> for FHCal X1 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x2_fhcal_x3 = new TProfile("p_yx_fhcal_x2_fhcal_x3", "<sin(PsiA)cos(PsiB)> for FHCal X2 & FHCal X3 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x2_fhcal_x4 = new TProfile("p_yx_fhcal_x2_fhcal_x4", "<sin(PsiA)cos(PsiB)> for FHCal X2 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x2_tof = new TProfile("p_yx_fhcal_x2_tof", "<sin(PsiA)cos(PsiB)> for FHCal X2 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x2_emc = new TProfile("p_yx_fhcal_x2_emc", "<sin(PsiA)cos(PsiB)> for FHCal X2 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x3_fhcal_x4 = new TProfile("p_yx_fhcal_x3_fhcal_x4", "<sin(PsiA)cos(PsiB)> for FHCal X3 & FHCal X4 vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x3_tof = new TProfile("p_yx_fhcal_x3_tof", "<sin(PsiA)cos(PsiB)> for FHCal X3 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x3_emc = new TProfile("p_yx_fhcal_x3_emc", "<sin(PsiA)cos(PsiB)> for FHCal X3 & EMC vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x4_tof = new TProfile("p_yx_fhcal_x4_tof", "<sin(PsiA)cos(PsiB)> for FHCal X4 & TOF vs centrality from TOF", 20, 0., 100.);
  TProfile *p_yx_fhcal_x4_emc = new TProfile("p_yx_fhcal_x4_emc", "<sin(PsiA)cos(PsiB)> for FHCal X4 & EMC vs centrality from TOF", 20, 0., 100.);

  ep::Event event;
  double qx_fhcal_x1, qx_fhcal_x2, qx_fhcal_x3, qx_fhcal_x4;
  double qx_tof, qx_emc;
  double qy_fhcal_x1, qy_fhcal_x2, qy_fhcal_x3, qy_fhcal_x4;
  double qy_tof, qy_emc;
  double sumw_fhcal_x1, sumw_fhcal_x2, sumw_fhcal_x3, sumw_fhcal_x4;
  double sumw_tof, sumw_emc;
  double psi_fhcal_x1, psi_fhcal_x2, psi_fhcal_x3, psi_fhcal_x4;
  double psi_tof, psi_emc;

  for (Long64_t iEv=0; iEv<nEvents; ++iEv) {
    if (iEv % 100 == 0) cout << "\r" << iEv << flush;
    if (!reader.ReadEvent(iEv, event))
      continue;

    int tofmult = event.tofHits.size();
    float cent = ep::GetCentMult(tofmult);
    if (cent < 0) continue;

    h_cent->Fill(cent);

    for (const auto &track : event.tracks) {
      // TODO: track loop
    }

    // for (const auto &module : event.fhcalModules) {
    //   // TODO: FHCal module loop
    // }

    qx_fhcal_x1 = 0.; qy_fhcal_x1 = 0.; sumw_fhcal_x1 = 0.; psi_fhcal_x1 = 0.;
    for (const auto &modid:x1_modules) {
      auto fhcalmod = event.fhcalModules.at(modid);
      h2_fhcal_xy_x1->Fill(fhcalmod.pos.x(), fhcalmod.pos.y(), fhcalmod.energy);
      qx_fhcal_x1 += fhcalmod.energy * cos( 1. * fhcalmod.pos.phi() );
      qy_fhcal_x1 += fhcalmod.energy * sin( 1. * fhcalmod.pos.phi() );
      sumw_fhcal_x1 += fhcalmod.energy;
    }

    qx_fhcal_x2 = 0.; qy_fhcal_x2 = 0.; sumw_fhcal_x2 = 0.; psi_fhcal_x2 = 0.;
    for (const auto &modid:x2_modules) {
      auto fhcalmod = event.fhcalModules.at(modid);
      h2_fhcal_xy_x2->Fill(fhcalmod.pos.x(), fhcalmod.pos.y(), fhcalmod.energy);
      qx_fhcal_x2 += fhcalmod.energy * cos( 1. * fhcalmod.pos.phi() );
      qy_fhcal_x2 += fhcalmod.energy * sin( 1. * fhcalmod.pos.phi() );
      sumw_fhcal_x2 += fhcalmod.energy;
    }

    qx_fhcal_x3 = 0.; qy_fhcal_x3 = 0.; sumw_fhcal_x3 = 0.; psi_fhcal_x3 = 0.;
    for (const auto &modid:x3_modules) {
      auto fhcalmod = event.fhcalModules.at(modid);
      h2_fhcal_xy_x3->Fill(fhcalmod.pos.x(), fhcalmod.pos.y(), fhcalmod.energy);
      qx_fhcal_x3 += fhcalmod.energy * cos( 1. * fhcalmod.pos.phi() );
      qy_fhcal_x3 += fhcalmod.energy * sin( 1. * fhcalmod.pos.phi() );
      sumw_fhcal_x3 += fhcalmod.energy;
    }

    qx_fhcal_x4 = 0.; qy_fhcal_x4 = 0.; sumw_fhcal_x4 = 0.; psi_fhcal_x4 = 0.;
    for (const auto &modid:x4_modules) {
      auto fhcalmod = event.fhcalModules.at(modid);
      h2_fhcal_xy_x4->Fill(fhcalmod.pos.x(), fhcalmod.pos.y(), fhcalmod.energy);
      qx_fhcal_x4 += fhcalmod.energy * cos( 1. * fhcalmod.pos.phi() );
      qy_fhcal_x4 += fhcalmod.energy * sin( 1. * fhcalmod.pos.phi() );
      sumw_fhcal_x4 += fhcalmod.energy;
    }

    qx_tof = 0.; qy_tof = 0.; sumw_tof = 0.; psi_tof = 0.;
    for (const auto &hit : event.tofHits) {
      h1_tof_posz->Fill(hit.pos.z());
      if (hit.pos.z() < tof_pos_cut.first) continue;
      if (hit.pos.z() > tof_pos_cut.second) continue;
      h1_tof_posz_after_cut->Fill(hit.pos.z());
      qx_tof += cos( 1. * hit.pos.phi() );
      qy_tof += sin( 1. * hit.pos.phi() );
      sumw_tof++;
    }

  qx_emc = 0.; qy_emc = 0.; sumw_emc = 0.; psi_emc = 0.;
    for (const auto &cluster : event.emcClusters) {
      h1_emc_posz->Fill(cluster.pos.z());
      if (cluster.pos.z() < emc_pos_cut.first) continue;
      if (cluster.pos.z() > emc_pos_cut.second) continue;
      if (cluster.energy < 0.05) continue;
      if (cluster.mult <= 1) continue;
      h1_emc_posz_after_cut->Fill(cluster.pos.z());
      qx_emc += cluster.energy * cos( 1. * cluster.pos.phi() );
      qy_emc += cluster.energy * sin( 1. * cluster.pos.phi() );
      sumw_emc += cluster.energy;
    }

    // Applying recentering
    qx_fhcal_x1 -= p_avg_qx_fhcal_x1->GetBinContent(p_avg_qx_fhcal_x1->FindBin(cent));
    qx_fhcal_x2 -= p_avg_qx_fhcal_x2->GetBinContent(p_avg_qx_fhcal_x2->FindBin(cent));
    qx_fhcal_x3 -= p_avg_qx_fhcal_x3->GetBinContent(p_avg_qx_fhcal_x3->FindBin(cent));
    qx_fhcal_x4 -= p_avg_qx_fhcal_x4->GetBinContent(p_avg_qx_fhcal_x4->FindBin(cent));
    qx_tof -= p_avg_qx_tof->GetBinContent(p_avg_qx_tof->FindBin(cent));
    qx_emc -= p_avg_qx_emc->GetBinContent(p_avg_qx_emc->FindBin(cent));

    qy_fhcal_x1 -= p_avg_qy_fhcal_x1->GetBinContent(p_avg_qy_fhcal_x1->FindBin(cent));
    qy_fhcal_x2 -= p_avg_qy_fhcal_x2->GetBinContent(p_avg_qy_fhcal_x2->FindBin(cent));
    qy_fhcal_x3 -= p_avg_qy_fhcal_x3->GetBinContent(p_avg_qy_fhcal_x3->FindBin(cent));
    qy_fhcal_x4 -= p_avg_qy_fhcal_x4->GetBinContent(p_avg_qy_fhcal_x4->FindBin(cent));
    qy_tof -= p_avg_qy_tof->GetBinContent(p_avg_qy_tof->FindBin(cent));
    qy_emc -= p_avg_qy_emc->GetBinContent(p_avg_qy_emc->FindBin(cent));

    psi_fhcal_x1 = 1./1. * atan2(qy_fhcal_x1, qx_fhcal_x1);
    psi_fhcal_x2 = 1./1. * atan2(qy_fhcal_x2, qx_fhcal_x2);
    psi_fhcal_x3 = 1./1. * atan2(qy_fhcal_x3, qx_fhcal_x3);
    psi_fhcal_x4 = 1./1. * atan2(qy_fhcal_x4, qx_fhcal_x4);
    psi_tof = 1./1. * atan2(qy_tof, qx_tof);
    psi_emc = 1./1. * atan2(qy_emc, qx_emc);

    // Applying shift (flattening)
    psi_fhcal_x1 +=  1./1. * (
      2./1. * (p_avg_cos1_fhcal_x1->GetBinContent(p_avg_cos1_fhcal_x1->FindBin(cent))*sin(1.*psi_fhcal_x1) - p_avg_sin1_fhcal_x1->GetBinContent(p_avg_sin1_fhcal_x1->FindBin(cent))*cos(1.*psi_fhcal_x1) ) +
      2./2. * (p_avg_cos2_fhcal_x1->GetBinContent(p_avg_cos2_fhcal_x1->FindBin(cent))*sin(2.*psi_fhcal_x1) - p_avg_sin2_fhcal_x1->GetBinContent(p_avg_sin2_fhcal_x1->FindBin(cent))*cos(2.*psi_fhcal_x1) ) +
      2./3. * (p_avg_cos3_fhcal_x1->GetBinContent(p_avg_cos3_fhcal_x1->FindBin(cent))*sin(3.*psi_fhcal_x1) - p_avg_sin3_fhcal_x1->GetBinContent(p_avg_sin3_fhcal_x1->FindBin(cent))*cos(3.*psi_fhcal_x1) ) +
      2./4. * (p_avg_cos4_fhcal_x1->GetBinContent(p_avg_cos4_fhcal_x1->FindBin(cent))*sin(4.*psi_fhcal_x1) - p_avg_sin4_fhcal_x1->GetBinContent(p_avg_sin4_fhcal_x1->FindBin(cent))*cos(4.*psi_fhcal_x1) ) +
      2./5. * (p_avg_cos5_fhcal_x1->GetBinContent(p_avg_cos5_fhcal_x1->FindBin(cent))*sin(5.*psi_fhcal_x1) - p_avg_sin5_fhcal_x1->GetBinContent(p_avg_sin5_fhcal_x1->FindBin(cent))*cos(5.*psi_fhcal_x1) ) +
      2./6. * (p_avg_cos6_fhcal_x1->GetBinContent(p_avg_cos6_fhcal_x1->FindBin(cent))*sin(6.*psi_fhcal_x1) - p_avg_sin6_fhcal_x1->GetBinContent(p_avg_sin6_fhcal_x1->FindBin(cent))*cos(6.*psi_fhcal_x1) ) +
      2./7. * (p_avg_cos7_fhcal_x1->GetBinContent(p_avg_cos7_fhcal_x1->FindBin(cent))*sin(7.*psi_fhcal_x1) - p_avg_sin7_fhcal_x1->GetBinContent(p_avg_sin7_fhcal_x1->FindBin(cent))*cos(7.*psi_fhcal_x1) ) +
      2./8. * (p_avg_cos8_fhcal_x1->GetBinContent(p_avg_cos8_fhcal_x1->FindBin(cent))*sin(8.*psi_fhcal_x1) - p_avg_sin8_fhcal_x1->GetBinContent(p_avg_sin8_fhcal_x1->FindBin(cent))*cos(8.*psi_fhcal_x1) )
    );

    psi_fhcal_x2 +=  1./1. * (
      2./1. * (p_avg_cos1_fhcal_x2->GetBinContent(p_avg_cos1_fhcal_x2->FindBin(cent))*sin(1.*psi_fhcal_x2) - p_avg_sin1_fhcal_x2->GetBinContent(p_avg_sin1_fhcal_x2->FindBin(cent))*cos(1.*psi_fhcal_x2) ) +
      2./2. * (p_avg_cos2_fhcal_x2->GetBinContent(p_avg_cos2_fhcal_x2->FindBin(cent))*sin(2.*psi_fhcal_x2) - p_avg_sin2_fhcal_x2->GetBinContent(p_avg_sin2_fhcal_x2->FindBin(cent))*cos(2.*psi_fhcal_x2) ) +
      2./3. * (p_avg_cos3_fhcal_x2->GetBinContent(p_avg_cos3_fhcal_x2->FindBin(cent))*sin(3.*psi_fhcal_x2) - p_avg_sin3_fhcal_x2->GetBinContent(p_avg_sin3_fhcal_x2->FindBin(cent))*cos(3.*psi_fhcal_x2) ) +
      2./4. * (p_avg_cos4_fhcal_x2->GetBinContent(p_avg_cos4_fhcal_x2->FindBin(cent))*sin(4.*psi_fhcal_x2) - p_avg_sin4_fhcal_x2->GetBinContent(p_avg_sin4_fhcal_x2->FindBin(cent))*cos(4.*psi_fhcal_x2) ) +
      2./5. * (p_avg_cos5_fhcal_x2->GetBinContent(p_avg_cos5_fhcal_x2->FindBin(cent))*sin(5.*psi_fhcal_x2) - p_avg_sin5_fhcal_x2->GetBinContent(p_avg_sin5_fhcal_x2->FindBin(cent))*cos(5.*psi_fhcal_x2) ) +
      2./6. * (p_avg_cos6_fhcal_x2->GetBinContent(p_avg_cos6_fhcal_x2->FindBin(cent))*sin(6.*psi_fhcal_x2) - p_avg_sin6_fhcal_x2->GetBinContent(p_avg_sin6_fhcal_x2->FindBin(cent))*cos(6.*psi_fhcal_x2) ) +
      2./7. * (p_avg_cos7_fhcal_x2->GetBinContent(p_avg_cos7_fhcal_x2->FindBin(cent))*sin(7.*psi_fhcal_x2) - p_avg_sin7_fhcal_x2->GetBinContent(p_avg_sin7_fhcal_x2->FindBin(cent))*cos(7.*psi_fhcal_x2) ) +
      2./8. * (p_avg_cos8_fhcal_x2->GetBinContent(p_avg_cos8_fhcal_x2->FindBin(cent))*sin(8.*psi_fhcal_x2) - p_avg_sin8_fhcal_x2->GetBinContent(p_avg_sin8_fhcal_x2->FindBin(cent))*cos(8.*psi_fhcal_x2) )
    );

    psi_fhcal_x3 +=  1./1. * (
      2./1. * (p_avg_cos1_fhcal_x3->GetBinContent(p_avg_cos1_fhcal_x3->FindBin(cent))*sin(1.*psi_fhcal_x3) - p_avg_sin1_fhcal_x3->GetBinContent(p_avg_sin1_fhcal_x3->FindBin(cent))*cos(1.*psi_fhcal_x3) ) +
      2./2. * (p_avg_cos2_fhcal_x3->GetBinContent(p_avg_cos2_fhcal_x3->FindBin(cent))*sin(2.*psi_fhcal_x3) - p_avg_sin2_fhcal_x3->GetBinContent(p_avg_sin2_fhcal_x3->FindBin(cent))*cos(2.*psi_fhcal_x3) ) +
      2./3. * (p_avg_cos3_fhcal_x3->GetBinContent(p_avg_cos3_fhcal_x3->FindBin(cent))*sin(3.*psi_fhcal_x3) - p_avg_sin3_fhcal_x3->GetBinContent(p_avg_sin3_fhcal_x3->FindBin(cent))*cos(3.*psi_fhcal_x3) ) +
      2./4. * (p_avg_cos4_fhcal_x3->GetBinContent(p_avg_cos4_fhcal_x3->FindBin(cent))*sin(4.*psi_fhcal_x3) - p_avg_sin4_fhcal_x3->GetBinContent(p_avg_sin4_fhcal_x3->FindBin(cent))*cos(4.*psi_fhcal_x3) ) +
      2./5. * (p_avg_cos5_fhcal_x3->GetBinContent(p_avg_cos5_fhcal_x3->FindBin(cent))*sin(5.*psi_fhcal_x3) - p_avg_sin5_fhcal_x3->GetBinContent(p_avg_sin5_fhcal_x3->FindBin(cent))*cos(5.*psi_fhcal_x3) ) +
      2./6. * (p_avg_cos6_fhcal_x3->GetBinContent(p_avg_cos6_fhcal_x3->FindBin(cent))*sin(6.*psi_fhcal_x3) - p_avg_sin6_fhcal_x3->GetBinContent(p_avg_sin6_fhcal_x3->FindBin(cent))*cos(6.*psi_fhcal_x3) ) +
      2./7. * (p_avg_cos7_fhcal_x3->GetBinContent(p_avg_cos7_fhcal_x3->FindBin(cent))*sin(7.*psi_fhcal_x3) - p_avg_sin7_fhcal_x3->GetBinContent(p_avg_sin7_fhcal_x3->FindBin(cent))*cos(7.*psi_fhcal_x3) ) +
      2./8. * (p_avg_cos8_fhcal_x3->GetBinContent(p_avg_cos8_fhcal_x3->FindBin(cent))*sin(8.*psi_fhcal_x3) - p_avg_sin8_fhcal_x3->GetBinContent(p_avg_sin8_fhcal_x3->FindBin(cent))*cos(8.*psi_fhcal_x3) )
    );

    psi_fhcal_x4 +=  1./1. * (
      2./1. * (p_avg_cos1_fhcal_x4->GetBinContent(p_avg_cos1_fhcal_x4->FindBin(cent))*sin(1.*psi_fhcal_x4) - p_avg_sin1_fhcal_x4->GetBinContent(p_avg_sin1_fhcal_x4->FindBin(cent))*cos(1.*psi_fhcal_x4) ) +
      2./2. * (p_avg_cos2_fhcal_x4->GetBinContent(p_avg_cos2_fhcal_x4->FindBin(cent))*sin(2.*psi_fhcal_x4) - p_avg_sin2_fhcal_x4->GetBinContent(p_avg_sin2_fhcal_x4->FindBin(cent))*cos(2.*psi_fhcal_x4) ) +
      2./3. * (p_avg_cos3_fhcal_x4->GetBinContent(p_avg_cos3_fhcal_x4->FindBin(cent))*sin(3.*psi_fhcal_x4) - p_avg_sin3_fhcal_x4->GetBinContent(p_avg_sin3_fhcal_x4->FindBin(cent))*cos(3.*psi_fhcal_x4) ) +
      2./4. * (p_avg_cos4_fhcal_x4->GetBinContent(p_avg_cos4_fhcal_x4->FindBin(cent))*sin(4.*psi_fhcal_x4) - p_avg_sin4_fhcal_x4->GetBinContent(p_avg_sin4_fhcal_x4->FindBin(cent))*cos(4.*psi_fhcal_x4) ) +
      2./5. * (p_avg_cos5_fhcal_x4->GetBinContent(p_avg_cos5_fhcal_x4->FindBin(cent))*sin(5.*psi_fhcal_x4) - p_avg_sin5_fhcal_x4->GetBinContent(p_avg_sin5_fhcal_x4->FindBin(cent))*cos(5.*psi_fhcal_x4) ) +
      2./6. * (p_avg_cos6_fhcal_x4->GetBinContent(p_avg_cos6_fhcal_x4->FindBin(cent))*sin(6.*psi_fhcal_x4) - p_avg_sin6_fhcal_x4->GetBinContent(p_avg_sin6_fhcal_x4->FindBin(cent))*cos(6.*psi_fhcal_x4) ) +
      2./7. * (p_avg_cos7_fhcal_x4->GetBinContent(p_avg_cos7_fhcal_x4->FindBin(cent))*sin(7.*psi_fhcal_x4) - p_avg_sin7_fhcal_x4->GetBinContent(p_avg_sin7_fhcal_x4->FindBin(cent))*cos(7.*psi_fhcal_x4) ) +
      2./8. * (p_avg_cos8_fhcal_x4->GetBinContent(p_avg_cos8_fhcal_x4->FindBin(cent))*sin(8.*psi_fhcal_x4) - p_avg_sin8_fhcal_x4->GetBinContent(p_avg_sin8_fhcal_x4->FindBin(cent))*cos(8.*psi_fhcal_x4) )
    );

    psi_tof +=  1./1. * (
      2./1. * (p_avg_cos1_tof->GetBinContent(p_avg_cos1_tof->FindBin(cent))*sin(1.*psi_tof) - p_avg_sin1_tof->GetBinContent(p_avg_sin1_tof->FindBin(cent))*cos(1.*psi_tof) ) +
      2./2. * (p_avg_cos2_tof->GetBinContent(p_avg_cos2_tof->FindBin(cent))*sin(2.*psi_tof) - p_avg_sin2_tof->GetBinContent(p_avg_sin2_tof->FindBin(cent))*cos(2.*psi_tof) ) +
      2./3. * (p_avg_cos3_tof->GetBinContent(p_avg_cos3_tof->FindBin(cent))*sin(3.*psi_tof) - p_avg_sin3_tof->GetBinContent(p_avg_sin3_tof->FindBin(cent))*cos(3.*psi_tof) ) +
      2./4. * (p_avg_cos4_tof->GetBinContent(p_avg_cos4_tof->FindBin(cent))*sin(4.*psi_tof) - p_avg_sin4_tof->GetBinContent(p_avg_sin4_tof->FindBin(cent))*cos(4.*psi_tof) ) +
      2./5. * (p_avg_cos5_tof->GetBinContent(p_avg_cos5_tof->FindBin(cent))*sin(5.*psi_tof) - p_avg_sin5_tof->GetBinContent(p_avg_sin5_tof->FindBin(cent))*cos(5.*psi_tof) ) +
      2./6. * (p_avg_cos6_tof->GetBinContent(p_avg_cos6_tof->FindBin(cent))*sin(6.*psi_tof) - p_avg_sin6_tof->GetBinContent(p_avg_sin6_tof->FindBin(cent))*cos(6.*psi_tof) ) +
      2./7. * (p_avg_cos7_tof->GetBinContent(p_avg_cos7_tof->FindBin(cent))*sin(7.*psi_tof) - p_avg_sin7_tof->GetBinContent(p_avg_sin7_tof->FindBin(cent))*cos(7.*psi_tof) ) +
      2./8. * (p_avg_cos8_tof->GetBinContent(p_avg_cos8_tof->FindBin(cent))*sin(8.*psi_tof) - p_avg_sin8_tof->GetBinContent(p_avg_sin8_tof->FindBin(cent))*cos(8.*psi_tof) )
    );

    psi_emc +=  1./1. * (
      2./1. * (p_avg_cos1_emc->GetBinContent(p_avg_cos1_emc->FindBin(cent))*sin(1.*psi_emc) - p_avg_sin1_emc->GetBinContent(p_avg_sin1_emc->FindBin(cent))*cos(1.*psi_emc) ) +
      2./2. * (p_avg_cos2_emc->GetBinContent(p_avg_cos2_emc->FindBin(cent))*sin(2.*psi_emc) - p_avg_sin2_emc->GetBinContent(p_avg_sin2_emc->FindBin(cent))*cos(2.*psi_emc) ) +
      2./3. * (p_avg_cos3_emc->GetBinContent(p_avg_cos3_emc->FindBin(cent))*sin(3.*psi_emc) - p_avg_sin3_emc->GetBinContent(p_avg_sin3_emc->FindBin(cent))*cos(3.*psi_emc) ) +
      2./4. * (p_avg_cos4_emc->GetBinContent(p_avg_cos4_emc->FindBin(cent))*sin(4.*psi_emc) - p_avg_sin4_emc->GetBinContent(p_avg_sin4_emc->FindBin(cent))*cos(4.*psi_emc) ) +
      2./5. * (p_avg_cos5_emc->GetBinContent(p_avg_cos5_emc->FindBin(cent))*sin(5.*psi_emc) - p_avg_sin5_emc->GetBinContent(p_avg_sin5_emc->FindBin(cent))*cos(5.*psi_emc) ) +
      2./6. * (p_avg_cos6_emc->GetBinContent(p_avg_cos6_emc->FindBin(cent))*sin(6.*psi_emc) - p_avg_sin6_emc->GetBinContent(p_avg_sin6_emc->FindBin(cent))*cos(6.*psi_emc) ) +
      2./7. * (p_avg_cos7_emc->GetBinContent(p_avg_cos7_emc->FindBin(cent))*sin(7.*psi_emc) - p_avg_sin7_emc->GetBinContent(p_avg_sin7_emc->FindBin(cent))*cos(7.*psi_emc) ) +
      2./8. * (p_avg_cos8_emc->GetBinContent(p_avg_cos8_emc->FindBin(cent))*sin(8.*psi_emc) - p_avg_sin8_emc->GetBinContent(p_avg_sin8_emc->FindBin(cent))*cos(8.*psi_emc) )
    );

    p_cos_fhcal_x1_fhcal_x2->Fill(cent, cos(1.*( psi_fhcal_x1 - psi_fhcal_x2 )));
    p_cos_fhcal_x1_fhcal_x3->Fill(cent, cos(1.*( psi_fhcal_x1 - psi_fhcal_x3 )));
    p_cos_fhcal_x1_fhcal_x4->Fill(cent, cos(1.*( psi_fhcal_x1 - psi_fhcal_x4 )));
    p_cos_fhcal_x1_tof->Fill(cent, cos(1.*( psi_fhcal_x1 - psi_tof )));
    p_cos_fhcal_x1_emc->Fill(cent, cos(1.*( psi_fhcal_x1 - psi_emc )));
    p_cos_fhcal_x2_fhcal_x3->Fill(cent, cos(1.*( psi_fhcal_x2 - psi_fhcal_x3 )));
    p_cos_fhcal_x2_fhcal_x4->Fill(cent, cos(1.*( psi_fhcal_x2 - psi_fhcal_x4 )));
    p_cos_fhcal_x2_tof->Fill(cent, cos(1.*( psi_fhcal_x2 - psi_tof )));
    p_cos_fhcal_x2_emc->Fill(cent, cos(1.*( psi_fhcal_x2 - psi_emc )));
    p_cos_fhcal_x3_fhcal_x4->Fill(cent, cos(1.*( psi_fhcal_x3 - psi_fhcal_x4 )));
    p_cos_fhcal_x3_tof->Fill(cent, cos(1.*( psi_fhcal_x3 - psi_tof )));
    p_cos_fhcal_x3_emc->Fill(cent, cos(1.*( psi_fhcal_x3 - psi_emc )));
    p_cos_fhcal_x4_tof->Fill(cent, cos(1.*( psi_fhcal_x4 - psi_tof )));
    p_cos_fhcal_x4_emc->Fill(cent, cos(1.*( psi_fhcal_x4 - psi_emc )));

    p_sin_fhcal_x1_fhcal_x2->Fill(cent, sin(1.*( psi_fhcal_x1 - psi_fhcal_x2 )));
    p_sin_fhcal_x1_fhcal_x3->Fill(cent, sin(1.*( psi_fhcal_x1 - psi_fhcal_x3 )));
    p_sin_fhcal_x1_fhcal_x4->Fill(cent, sin(1.*( psi_fhcal_x1 - psi_fhcal_x4 )));
    p_sin_fhcal_x1_tof->Fill(cent, sin(1.*( psi_fhcal_x1 - psi_tof )));
    p_sin_fhcal_x1_emc->Fill(cent, sin(1.*( psi_fhcal_x1 - psi_emc )));
    p_sin_fhcal_x2_fhcal_x3->Fill(cent, sin(1.*( psi_fhcal_x2 - psi_fhcal_x3 )));
    p_sin_fhcal_x2_fhcal_x4->Fill(cent, sin(1.*( psi_fhcal_x2 - psi_fhcal_x4 )));
    p_sin_fhcal_x2_tof->Fill(cent, sin(1.*( psi_fhcal_x2 - psi_tof )));
    p_sin_fhcal_x2_emc->Fill(cent, sin(1.*( psi_fhcal_x2 - psi_emc )));
    p_sin_fhcal_x3_fhcal_x4->Fill(cent, sin(1.*( psi_fhcal_x3 - psi_fhcal_x4 )));
    p_sin_fhcal_x3_tof->Fill(cent, sin(1.*( psi_fhcal_x3 - psi_tof )));
    p_sin_fhcal_x3_emc->Fill(cent, sin(1.*( psi_fhcal_x3 - psi_emc )));
    p_sin_fhcal_x4_tof->Fill(cent, sin(1.*( psi_fhcal_x4 - psi_tof )));
    p_sin_fhcal_x4_emc->Fill(cent, sin(1.*( psi_fhcal_x4 - psi_emc )));

    p_xx_fhcal_x1_fhcal_x2->Fill(cent, cos(1.*(psi_fhcal_x1)) * cos(1.*(psi_fhcal_x2)));
    p_xx_fhcal_x1_fhcal_x3->Fill(cent, cos(1.*(psi_fhcal_x1)) * cos(1.*(psi_fhcal_x3)));
    p_xx_fhcal_x1_fhcal_x4->Fill(cent, cos(1.*(psi_fhcal_x1)) * cos(1.*(psi_fhcal_x4)));
    p_xx_fhcal_x1_tof->Fill(cent,      cos(1.*(psi_fhcal_x1)) * cos(1.*(psi_tof)));
    p_xx_fhcal_x1_emc->Fill(cent,      cos(1.*(psi_fhcal_x1)) * cos(1.*(psi_emc)));
    p_xx_fhcal_x2_fhcal_x3->Fill(cent, cos(1.*(psi_fhcal_x2)) * cos(1.*(psi_fhcal_x3)));
    p_xx_fhcal_x2_fhcal_x4->Fill(cent, cos(1.*(psi_fhcal_x2)) * cos(1.*(psi_fhcal_x4)));
    p_xx_fhcal_x2_tof->Fill(cent,      cos(1.*(psi_fhcal_x2)) * cos(1.*(psi_tof)));
    p_xx_fhcal_x2_emc->Fill(cent,      cos(1.*(psi_fhcal_x2)) * cos(1.*(psi_emc)));
    p_xx_fhcal_x3_fhcal_x4->Fill(cent, cos(1.*(psi_fhcal_x3)) * cos(1.*(psi_fhcal_x4)));
    p_xx_fhcal_x3_tof->Fill(cent,      cos(1.*(psi_fhcal_x3)) * cos(1.*(psi_tof)));
    p_xx_fhcal_x3_emc->Fill(cent,      cos(1.*(psi_fhcal_x3)) * cos(1.*(psi_tof)));
    p_xx_fhcal_x4_tof->Fill(cent,      cos(1.*(psi_fhcal_x4)) * cos(1.*(psi_tof)));
    p_xx_fhcal_x4_emc->Fill(cent,      cos(1.*(psi_fhcal_x4)) * cos(1.*(psi_tof)));

    p_yy_fhcal_x1_fhcal_x2->Fill(cent, sin(1.*(psi_fhcal_x1)) * sin(1.*(psi_fhcal_x2)));
    p_yy_fhcal_x1_fhcal_x3->Fill(cent, sin(1.*(psi_fhcal_x1)) * sin(1.*(psi_fhcal_x3)));
    p_yy_fhcal_x1_fhcal_x4->Fill(cent, sin(1.*(psi_fhcal_x1)) * sin(1.*(psi_fhcal_x4)));
    p_yy_fhcal_x1_tof->Fill(cent,      sin(1.*(psi_fhcal_x1)) * sin(1.*(psi_tof)));
    p_yy_fhcal_x1_emc->Fill(cent,      sin(1.*(psi_fhcal_x1)) * sin(1.*(psi_emc)));
    p_yy_fhcal_x2_fhcal_x3->Fill(cent, sin(1.*(psi_fhcal_x2)) * sin(1.*(psi_fhcal_x3)));
    p_yy_fhcal_x2_fhcal_x4->Fill(cent, sin(1.*(psi_fhcal_x2)) * sin(1.*(psi_fhcal_x4)));
    p_yy_fhcal_x2_tof->Fill(cent,      sin(1.*(psi_fhcal_x2)) * sin(1.*(psi_tof)));
    p_yy_fhcal_x2_emc->Fill(cent,      sin(1.*(psi_fhcal_x2)) * sin(1.*(psi_emc)));
    p_yy_fhcal_x3_fhcal_x4->Fill(cent, sin(1.*(psi_fhcal_x3)) * sin(1.*(psi_fhcal_x4)));
    p_yy_fhcal_x3_tof->Fill(cent,      sin(1.*(psi_fhcal_x3)) * sin(1.*(psi_tof)));
    p_yy_fhcal_x3_emc->Fill(cent,      sin(1.*(psi_fhcal_x3)) * sin(1.*(psi_tof)));
    p_yy_fhcal_x4_tof->Fill(cent,      sin(1.*(psi_fhcal_x4)) * sin(1.*(psi_tof)));
    p_yy_fhcal_x4_emc->Fill(cent,      sin(1.*(psi_fhcal_x4)) * sin(1.*(psi_tof)));

    p_xy_fhcal_x1_fhcal_x2->Fill(cent, cos(1.*(psi_fhcal_x1)) * sin(1.*(psi_fhcal_x2)));
    p_xy_fhcal_x1_fhcal_x3->Fill(cent, cos(1.*(psi_fhcal_x1)) * sin(1.*(psi_fhcal_x3)));
    p_xy_fhcal_x1_fhcal_x4->Fill(cent, cos(1.*(psi_fhcal_x1)) * sin(1.*(psi_fhcal_x4)));
    p_xy_fhcal_x1_tof->Fill(cent,      cos(1.*(psi_fhcal_x1)) * sin(1.*(psi_tof)));
    p_xy_fhcal_x1_emc->Fill(cent,      cos(1.*(psi_fhcal_x1)) * sin(1.*(psi_emc)));
    p_xy_fhcal_x2_fhcal_x3->Fill(cent, cos(1.*(psi_fhcal_x2)) * sin(1.*(psi_fhcal_x3)));
    p_xy_fhcal_x2_fhcal_x4->Fill(cent, cos(1.*(psi_fhcal_x2)) * sin(1.*(psi_fhcal_x4)));
    p_xy_fhcal_x2_tof->Fill(cent,      cos(1.*(psi_fhcal_x2)) * sin(1.*(psi_tof)));
    p_xy_fhcal_x2_emc->Fill(cent,      cos(1.*(psi_fhcal_x2)) * sin(1.*(psi_emc)));
    p_xy_fhcal_x3_fhcal_x4->Fill(cent, cos(1.*(psi_fhcal_x3)) * sin(1.*(psi_fhcal_x4)));
    p_xy_fhcal_x3_tof->Fill(cent,      cos(1.*(psi_fhcal_x3)) * sin(1.*(psi_tof)));
    p_xy_fhcal_x3_emc->Fill(cent,      cos(1.*(psi_fhcal_x3)) * sin(1.*(psi_tof)));
    p_xy_fhcal_x4_tof->Fill(cent,      cos(1.*(psi_fhcal_x4)) * sin(1.*(psi_tof)));
    p_xy_fhcal_x4_emc->Fill(cent,      cos(1.*(psi_fhcal_x4)) * sin(1.*(psi_tof)));

    p_yx_fhcal_x1_fhcal_x2->Fill(cent, sin(1.*(psi_fhcal_x1)) * cos(1.*(psi_fhcal_x2)));
    p_yx_fhcal_x1_fhcal_x3->Fill(cent, sin(1.*(psi_fhcal_x1)) * cos(1.*(psi_fhcal_x3)));
    p_yx_fhcal_x1_fhcal_x4->Fill(cent, sin(1.*(psi_fhcal_x1)) * cos(1.*(psi_fhcal_x4)));
    p_yx_fhcal_x1_tof->Fill(cent,      sin(1.*(psi_fhcal_x1)) * cos(1.*(psi_tof)));
    p_yx_fhcal_x1_emc->Fill(cent,      sin(1.*(psi_fhcal_x1)) * cos(1.*(psi_emc)));
    p_yx_fhcal_x2_fhcal_x3->Fill(cent, sin(1.*(psi_fhcal_x2)) * cos(1.*(psi_fhcal_x3)));
    p_yx_fhcal_x2_fhcal_x4->Fill(cent, sin(1.*(psi_fhcal_x2)) * cos(1.*(psi_fhcal_x4)));
    p_yx_fhcal_x2_tof->Fill(cent,      sin(1.*(psi_fhcal_x2)) * cos(1.*(psi_tof)));
    p_yx_fhcal_x2_emc->Fill(cent,      sin(1.*(psi_fhcal_x2)) * cos(1.*(psi_emc)));
    p_yx_fhcal_x3_fhcal_x4->Fill(cent, sin(1.*(psi_fhcal_x3)) * cos(1.*(psi_fhcal_x4)));
    p_yx_fhcal_x3_tof->Fill(cent,      sin(1.*(psi_fhcal_x3)) * cos(1.*(psi_tof)));
    p_yx_fhcal_x3_emc->Fill(cent,      sin(1.*(psi_fhcal_x3)) * cos(1.*(psi_tof)));
    p_yx_fhcal_x4_tof->Fill(cent,      sin(1.*(psi_fhcal_x4)) * cos(1.*(psi_tof)));
    p_yx_fhcal_x4_emc->Fill(cent,      sin(1.*(psi_fhcal_x4)) * cos(1.*(psi_tof)));
  }

  // Write output
  h_cent->Write();

  h2_fhcal_xy_x1->Write();
  h2_fhcal_xy_x2->Write();
  h2_fhcal_xy_x3->Write();
  h2_fhcal_xy_x4->Write();

  h1_tof_posz->Write();
  h1_tof_posz_after_cut->Write();
  h1_emc_posz->Write();
  h1_emc_posz_after_cut->Write();

  p_cos_fhcal_x1_fhcal_x2->Write();
  p_cos_fhcal_x1_fhcal_x3->Write();
  p_cos_fhcal_x1_fhcal_x4->Write();
  p_cos_fhcal_x1_tof->Write();
  p_cos_fhcal_x1_emc->Write();
  p_cos_fhcal_x2_fhcal_x3->Write();
  p_cos_fhcal_x2_fhcal_x4->Write();
  p_cos_fhcal_x2_tof->Write();
  p_cos_fhcal_x2_emc->Write();
  p_cos_fhcal_x3_fhcal_x4->Write();
  p_cos_fhcal_x3_tof->Write();
  p_cos_fhcal_x3_emc->Write();
  p_cos_fhcal_x4_tof->Write();
  p_cos_fhcal_x4_emc->Write();

  p_sin_fhcal_x1_fhcal_x2->Write();
  p_sin_fhcal_x1_fhcal_x3->Write();
  p_sin_fhcal_x1_fhcal_x4->Write();
  p_sin_fhcal_x1_tof->Write();
  p_sin_fhcal_x1_emc->Write();
  p_sin_fhcal_x2_fhcal_x3->Write();
  p_sin_fhcal_x2_fhcal_x4->Write();
  p_sin_fhcal_x2_tof->Write();
  p_sin_fhcal_x2_emc->Write();
  p_sin_fhcal_x3_fhcal_x4->Write();
  p_sin_fhcal_x3_tof->Write();
  p_sin_fhcal_x3_emc->Write();
  p_sin_fhcal_x4_tof->Write();
  p_sin_fhcal_x4_emc->Write();

  p_xx_fhcal_x1_fhcal_x2->Write();
  p_xx_fhcal_x1_fhcal_x3->Write();
  p_xx_fhcal_x1_fhcal_x4->Write();
  p_xx_fhcal_x1_tof->Write();
  p_xx_fhcal_x1_emc->Write();
  p_xx_fhcal_x2_fhcal_x3->Write();
  p_xx_fhcal_x2_fhcal_x4->Write();
  p_xx_fhcal_x2_tof->Write();
  p_xx_fhcal_x2_emc->Write();
  p_xx_fhcal_x3_fhcal_x4->Write();
  p_xx_fhcal_x3_tof->Write();
  p_xx_fhcal_x3_emc->Write();
  p_xx_fhcal_x4_tof->Write();
  p_xx_fhcal_x4_emc->Write();

  p_yy_fhcal_x1_fhcal_x2->Write();
  p_yy_fhcal_x1_fhcal_x3->Write();
  p_yy_fhcal_x1_fhcal_x4->Write();
  p_yy_fhcal_x1_tof->Write();
  p_yy_fhcal_x1_emc->Write();
  p_yy_fhcal_x2_fhcal_x3->Write();
  p_yy_fhcal_x2_fhcal_x4->Write();
  p_yy_fhcal_x2_tof->Write();
  p_yy_fhcal_x2_emc->Write();
  p_yy_fhcal_x3_fhcal_x4->Write();
  p_yy_fhcal_x3_tof->Write();
  p_yy_fhcal_x3_emc->Write();
  p_yy_fhcal_x4_tof->Write();
  p_yy_fhcal_x4_emc->Write();

  p_xy_fhcal_x1_fhcal_x2->Write();
  p_xy_fhcal_x1_fhcal_x3->Write();
  p_xy_fhcal_x1_fhcal_x4->Write();
  p_xy_fhcal_x1_tof->Write();
  p_xy_fhcal_x1_emc->Write();
  p_xy_fhcal_x2_fhcal_x3->Write();
  p_xy_fhcal_x2_fhcal_x4->Write();
  p_xy_fhcal_x2_tof->Write();
  p_xy_fhcal_x2_emc->Write();
  p_xy_fhcal_x3_fhcal_x4->Write();
  p_xy_fhcal_x3_tof->Write();
  p_xy_fhcal_x3_emc->Write();
  p_xy_fhcal_x4_tof->Write();
  p_xy_fhcal_x4_emc->Write();

  p_yx_fhcal_x1_fhcal_x2->Write();
  p_yx_fhcal_x1_fhcal_x3->Write();
  p_yx_fhcal_x1_fhcal_x4->Write();
  p_yx_fhcal_x1_tof->Write();
  p_yx_fhcal_x1_emc->Write();
  p_yx_fhcal_x2_fhcal_x3->Write();
  p_yx_fhcal_x2_fhcal_x4->Write();
  p_yx_fhcal_x2_tof->Write();
  p_yx_fhcal_x2_emc->Write();
  p_yx_fhcal_x3_fhcal_x4->Write();
  p_yx_fhcal_x3_tof->Write();
  p_yx_fhcal_x3_emc->Write();
  p_yx_fhcal_x4_tof->Write();
  p_yx_fhcal_x4_emc->Write();

  fOut.Close();

  timer.Stop();
  timer.Print();
}
