#include <iostream>
#include <vector>
#include <map>

#include <TMath.h>
#include <TChain.h>
#include <TDatabasePDG.h>
#include <TCollection.h>
#include <TFileCollection.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiE4D.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TMatrixDSym.h>
#include <TMatrixDUtilsfwd.h>

#include "MpdConstField.h"
#include "MpdFieldMap.h"
#include "MpdEvent.h"
#include "MpdVertex.h"
#include "MpdTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdMCTrack.h"
#include "MpdZdcDigi.h"

using namespace ROOT;
using namespace ROOT::Math;
using namespace ROOT::RDF;
using ROOT::VecOps::Map;
using fourVector=LorentzVector<PtEtaPhiE4D<double>>;

#ifdef __CINT__

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ class vector<vector<float>>+;

#endif

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

vector<fourVector> simMomentum(RVec<MpdMCTrack> &tracks)
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

vector<XYZTVector> simPosStart(const RVec<MpdMCTrack> &tracks)
{
  vector<XYZTVector> pos;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    pos.push_back({track.GetStartX(),track.GetStartY(),track.GetStartZ(),track.GetStartT()});
  }
  return pos;
}

RVec<int> simMotherId(const RVec<MpdMCTrack> &tracks)
{
  vector<int> mothId;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track))
      continue;
    mothId.push_back(track.GetMotherId());
  }
  return mothId;
}

RVec<int> simPdg(const RVec<MpdMCTrack> &tracks)
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

RVec<MpdTrack> getGlobalTracks(MpdEvent &event)
{
  vector<MpdTrack> tracks;
  auto trArray = (TClonesArray*) event.GetGlobalTracks();
  int Ntracks = trArray->GetEntriesFast();
  for (int i=0; i<Ntracks; i++) {
    auto track = (MpdTrack*) trArray->At(i);
    tracks.push_back(*track);
  }
  return tracks;
}

RVec<int> recNhits(const RVec<MpdTrack> &tracks)
try {
  vector<int> nhits;
  for (auto &track:tracks) {
    int nhit = track.GetNofHits();
    nhits.push_back(nhit);
  }
  return nhits;
} catch(const std::exception& e ){
  std::cout << __func__ << std::endl;
}

RVec<int> recNhitsFit(const RVec<MpdTrack> &tracks)
try {
  vector<int> nhits;
  for (auto &track:tracks) {
    int nhit = track.GetNofHitsFitTpc();
    nhits.push_back(nhit);
  }
  return nhits;
} catch(const std::exception& e ){
  std::cout << __func__ << std::endl;
}

RVec<int> recNhitsPoss(const RVec<MpdTrack> &tracks)
try {
  vector<int> nhits;
  for (auto &track:tracks) {
    int nhit = track.GetNofHitsPossTpc();
    nhits.push_back(nhit);
  }
  return nhits;
} catch(const std::exception& e ){
  std::cout << __func__ << std::endl;
}

vector<float> trackP(const RVec<MpdTrack> &tracks)
try {
  vector<float> momenta;
  for (auto &track:tracks) {
    momenta.push_back(sqrt(track.GetPt()*track.GetPt()+track.GetPz()*track.GetPz()));
  }
  return momenta;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<fourVector> trackMomentum(const RVec<MpdTrack> &tracks)
try {
  vector<fourVector> momenta;
  for (auto &track:tracks) {
    momenta.push_back({track.GetPt(), track.GetEta(), track.GetPhi(), 0});
  }
  return momenta;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<short> recCharge(const RVec<MpdTrack> &tracks)
try {
  vector<short> charge;
  for (auto track:tracks) {
    int q = track.GetCharge();
    charge.push_back(q);
  }
  return charge;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> recDca(const RVec<MpdTrack> &tracks)
try {
  vector<XYZVector> dca;
  for (auto &track:tracks) {
    dca.push_back({track.GetDCAX(), track.GetDCAY(), track.GetDCAZ()});
  }
  return dca;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<vector<float>> trGlobalFirstParam(const RVec<MpdTrack> &tracks)
try {
  vector<vector<float>> parameters;
  for (auto &track:tracks) {
    parameters.emplace_back();
    parameters.back().push_back( track.GetFirstPointX() );
    parameters.back().push_back( track.GetFirstPointY() );
    parameters.back().push_back( track.GetFirstPointY() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<vector<float>> trGlobalLastParam(const RVec<MpdTrack> &tracks)
try {
  vector<vector<float>> parameters;
  for (auto &track:tracks) {
    parameters.emplace_back();
    parameters.back().push_back( track.GetLastPointX() );
    parameters.back().push_back( track.GetLastPointY() );
    parameters.back().push_back( track.GetLastPointY() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<vector<float> > covMatrix(const RVec<MpdTpcKalmanTrack> &kalman_tracks)
{
  vector<vector<float>> covariance_matrix;
  for (auto &kalman_track:kalman_tracks) {
    auto *cov = kalman_track.GetCovariance();
    int ncols = cov->GetNcols();
    // int nrows = cov->GetNrows();
    // cout << "\t\tNcols = " << ncols << ", Nrows = " << nrows << endl;
    covariance_matrix.emplace_back();
    for (int i=0; i<ncols; ++i) {
      for (int j=0; j<i; ++j) {
        covariance_matrix.back().push_back( (float)TMatrixDRow(*cov, i)(j) );
      }
    }
  }
  return covariance_matrix;
}

RVec<int> recSimIndex(const RVec<MpdTrack> &recoTracks, const RVec<MpdMCTrack> simTracks)
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
  for (auto track:recoTracks) {
    int oldIndex=track.GetID();
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

RVec<int> trSimPdg(const RVec<int> &sim_index, const RVec<int> &sim_pdg)
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

RVec<int> trSimMotherId(const RVec<int> &sim_index, const RVec<int> &sim_motherId)
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

RVec<float> fhcalModE(const RVec<MpdZdcDigi> &fhcalHits)
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

RVec<int> moduleId (const vector<XYZVector> &modulePos)
try {
  vector <int> moduleIds;
  for (int i=0;i<modulePos.size();i++)
    moduleIds.push_back(i+1);
  return moduleIds;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<bool> hasHitFhcal (const RVec<MpdMCTrack> &particles)
{
  vector<bool> hasHit;
  for(auto &part:particles)
    hasHit.push_back(part.GetNPoints(kZDC)>0);
  return hasHit;
}

XYZVector vtxPos(const MpdVertex &vertex)
{
  return {vertex.GetX(), vertex.GetY(), vertex.GetZ()};
}

double vtxChi2(const MpdVertex &vertex)
{
  return vertex.GetChi2();
}

int vtxNdf(const MpdVertex &vertex)
{
  return vertex.GetNDF();
}

double vtxChi2Ndf(const MpdVertex &vertex)
{
  return (double)(vertex.GetChi2())/(double)(vertex.GetNDF());
}

int vtxNtracks(const MpdVertex &vertex)
{
  return vertex.GetNTracks();
}

void convertMPD(string inDst="", string fileOut="")
{
  TStopwatch timer;
  timer.Start();
  
  TChain *chainRec=makeChain(inDst, "mpdsim");
  ROOT::RDataFrame d(*chainRec);

  int nEvents = chainRec->GetEntries();

  auto fhcalModPos = modulePos();

  //magField = new MpdFieldMap("B-field_v2", "A");
  magField = new MpdConstField();
  magField->SetField(0., 0., 5.); // values are in kG:  1T = 10kG
  magField->SetFieldRegion(-230, 230, -230, 230, -375, 375); // values in cm
  cout << "FIELD at (0., 0., 0.) = (" << magField->GetBx(0., 0., 0.)
       << "; " << magField->GetBy(0., 0., 0.)
       << "; " << magField->GetBz(0., 0., 0.) << ")" << endl;


  auto dd = d
    .Define("evtId","MCEventHeader.fEventId")
    .Define("recoPrimVtxX","MPDEvent.PrimaryVerticesX")
    .Define("recoPrimVtxY","MPDEvent.PrimaryVerticesY")
    .Define("recoPrimVtxZ","MPDEvent.PrimaryVerticesZ")
    .Define("recoPrimVtxChi2", "MPDEvent.PrimaryVerticesChi2")
    //.Define("recoVtxPos", vtxPos, {"Vertex"})
    //.Define("recoVtxChi2", vtxChi2, {"Vertex"})
    //.Define("recoVtxNDF", vtxNdf, {"Vertex"})
    //.Define("recoVtxChi2NDF", vtxChi2Ndf, {"Vertex"})
    //.Define("recoVtxNtracks", vtxNtracks, {"Vertex"})
    .Define("mcVtxX","MCEventHeader.fX")
    .Define("mcVtxY","MCEventHeader.fY")
    .Define("mcVtxZ","MCEventHeader.fZ")
    .Define("mcB", "MCEventHeader.fB")
    .Define("mcRP", "MCEventHeader.fRotZ")
    .Define("recoGlobalTracks", getGlobalTracks, {"MPDEvent."})
    .Define("recoGlobalMom", trackMomentum, {"recoGlobalTracks"})
    .Define("recoGlobalNhits", recNhits, {"recoGlobalTracks"})
    .Define("recoGlobalNhitsFit", recNhitsFit, {"recoGlobalTracks"})
    .Define("recoGlobalNhitsPoss", recNhitsPoss, {"recoGlobalTracks"})
    .Define("recoGlobalChi2", "MPDEvent.fGlobalTracks.fChi2")
    .Define("recoGlobalP", trackP, {"recoGlobalTracks"})
    .Define("recoGlobalTofFlag", "MPDEvent.fGlobalTracks.fTofFlag")
    .Define("recoGlobalCharge",recCharge,{"recoGlobalTracks"})
    .Define("recoGlobalDca", recDca, {"recoGlobalTracks"})
    .Define("recoGlobalTofMass2", "MPDEvent.fGlobalTracks.fTofMass2")
    .Define("recoGlobalParamFirst", trGlobalFirstParam, {"recoGlobalTracks"})
    .Define("recoGlobalParamLast", trGlobalLastParam, {"recoGlobalTracks"})
    .Define("recoKalmanCovMatrix", covMatrix, {"TpcKalmanTrack"})
    .Define("recoKalmanNofWrong", "TpcKalmanTrack.fNofWrong")
    .Define("recoKalmanLength", "TpcKalmanTrack.fLength")
    .Define("recoKalmanCh2Ndf", "TpcKalmanTrack.fChi2")
    .Define("simMom", simMomentum, {"MCTrack"})
    .Define("simPosStart", simPosStart, {"MCTrack"})
    .Define("simMotherId", simMotherId, {"MCTrack"})
    .Define("simPdg", simPdg, {"MCTrack"})
    .Define("simCharge", simCharge, {"simPdg"})
    .Define("simHasHitFHCal", hasHitFhcal, {"MCTrack"})
    .Define("recoGlobalSimIndex", recSimIndex, {"recoGlobalTracks", "MCTrack"})
    .Define("recoGlobalSimPdg", trSimPdg, {"recoGlobalSimIndex", "simPdg"})
    .Define("recoGlobalSimMotherId", trSimMotherId, {"recoGlobalSimIndex", "simMotherId"})
    .Define("fhcalModPos", [fhcalModPos](){return fhcalModPos; })
    .Define("fhcalModId", moduleId, {"fhcalModPos"})
    .Define("fhcalModE",fhcalModE,{"ZdcDigi"})
  ;
  dd.Foreach([](uint evtId){if (evtId % 100 == 0) cout << "\n" << evtId;}, {"evtId"}); // progress display 
  cout << endl;

  vector<string> definedNames;
  vector<string> toExclude={"recoGlobalTracks"};
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
