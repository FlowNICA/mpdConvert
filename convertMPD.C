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

RVec<MpdMCTrack> getAssocSimTracks(const RVec<MpdTrack> &mpdtracks, const RVec<MpdMCTrack> &mctracks)
{
  // Returns vector of mc tracks that have assotiation with mpdtrack
  vector<MpdMCTrack> goodtracks;
  for (auto &track:mpdtracks) {
    auto mctrack = mctracks.at(track.GetID());
    goodtracks.push_back(mctrack);
  }
  return goodtracks;
}

bool isGoodSimTrack(const MpdMCTrack &track, const RVec<MpdMCTrack> &assocMcTracks, double mcvtxX, double mcvtxY, double mcvtxZ)
{
  //bool isInReco = (std::find(assocMcTracks.begin(), assocMcTracks.end(), track) == assocMcTracks.end());
  //if(isInReco && track.GetMotherId() != -1) return false;
  //if (track.GetNPoints(kTPC)<=0) return false; // only mc tracks that have at least 1 point in TPC!
  auto startX = track.GetStartX();
  auto startY = track.GetStartY();
  auto startZ = track.GetStartZ();
  auto startR = sqrt(startX*startX+startY*startY);
  auto dist = sqrt(pow(mcvtxX-startX,2)+pow(mcvtxY-startY,2)+pow(mcvtxZ-startZ,2));
  //if ( startR>=140. && abs(startZ)>=140. ) return false; // reject mc tracks outside of the TPC region
  if (dist>=5. && track.GetNPoints(kTPC)<=0) return false; // reject mc tracks with dist>5 cm and 0 mc points in TPC
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

vector<fourVector> simMomentum(RVec<MpdMCTrack> &tracks, const RVec<MpdMCTrack> &assocMcTracks, double mcvtxX, double mcvtxY, double mcvtxZ)
{
  vector<fourVector> momenta;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track, assocMcTracks, mcvtxX, mcvtxY, mcvtxZ))
      continue;
    TVector3 mom;
    track.GetMomentum(mom);
    double mass=getMass(track.GetPdgCode());
    momenta.push_back({mom.Pt(),mom.Eta(),mom.Phi(),sqrt(mass*mass+mom.Mag2())});
  }
  return momenta;
}

vector<XYZTVector> simPosStart(const RVec<MpdMCTrack> &tracks, const RVec<MpdMCTrack> &assocMcTracks, double mcvtxX, double mcvtxY, double mcvtxZ)
{
  vector<XYZTVector> pos;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track, assocMcTracks, mcvtxX, mcvtxY, mcvtxZ))
      continue;
    pos.push_back({track.GetStartX(),track.GetStartY(),track.GetStartZ(),track.GetStartT()});
  }
  return pos;
}

RVec<int> simMotherId(const RVec<MpdMCTrack> &tracks, const RVec<MpdMCTrack> &assocMcTracks, double mcvtxX, double mcvtxY, double mcvtxZ)
{
  vector<int> mothId;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track, assocMcTracks, mcvtxX, mcvtxY, mcvtxZ))
      continue;
    mothId.push_back(track.GetMotherId());
  }
  return mothId;
}

RVec<int> simPdg(const RVec<MpdMCTrack> &tracks, const RVec<MpdMCTrack> &assocMcTracks, double mcvtxX, double mcvtxY, double mcvtxZ)
{
  vector<int> pdg;
  for (auto& track:tracks) {
    if (!isGoodSimTrack(track, assocMcTracks, mcvtxX, mcvtxY, mcvtxZ))
      continue;
    pdg.push_back(track.GetPdgCode());
  }
  return pdg;
}

RVec<bool> hasHitFhcal (const RVec<MpdMCTrack> &tracks, const RVec<MpdMCTrack> &assocMcTracks, double mcvtxX, double mcvtxY, double mcvtxZ)
{
  vector<bool> hasHit;
  for(auto &track:tracks) {
    if (!isGoodSimTrack(track, assocMcTracks, mcvtxX, mcvtxY, mcvtxZ))
      continue;
    hasHit.push_back(track.GetNPoints(kZDC)>0);
  }
  return hasHit;
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
  throw e;
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
  throw e;
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
  throw e;
}

RVec<bool> recHasTofHit(const RVec<MpdTrack> &tracks)
try {
  vector<bool> tofhits;
  for (auto &track:tracks) {
    bool hit = (track.GetTofFlag() == 2 || track.GetTofFlag() == 6);
    tofhits.push_back(hit);
  }
  return tofhits;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> recTofMass2(const RVec<MpdTrack> &tracks)
try {
  vector<float> mass2;
  for (auto &track:tracks) {
    mass2.push_back(track.GetTofMass2());
  }
  return mass2;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> recTpcDedx(const RVec<MpdTrack> &tracks)
try{
  vector<float> dedx;
  for (auto& track:tracks) {
    dedx.push_back(track.GetdEdXTPC());
  }
  return dedx;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
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

RVec<float> recChi2(const RVec<MpdTrack> &tracks)
try {
  vector<float> chi2;
  for (auto &track:tracks) {
    chi2.push_back(track.GetChi2());
  }
  return chi2;
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

vector< vector<float> > trGlobalFirstParam(const RVec<MpdTrack> &tracks)
try {
  vector<vector<float>> parameters;
  for (auto &track:tracks) {
    parameters.emplace_back();
    parameters.back().push_back( track.GetFirstPointX() );
    parameters.back().push_back( track.GetFirstPointY() );
    parameters.back().push_back( track.GetFirstPointZ() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector< vector<float> > trGlobalLastParam(const RVec<MpdTrack> &tracks)
try {
  vector<vector<float>> parameters;
  for (auto &track:tracks) {
    parameters.emplace_back();
    parameters.back().push_back( track.GetLastPointX() );
    parameters.back().push_back( track.GetLastPointY() );
    parameters.back().push_back( track.GetLastPointZ() );
  }
  return parameters;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

bool propagateKalman2vertex(MpdTpcKalmanTrack in, MpdTpcKalmanTrack &out)
{
  // Propagate kalman track to the dca to vertex
  MpdTpcKalmanTrack tmp(in);
  bool verbose = false;

  // Resetting params and cov. matrix to those ones at dca2beamline ...
  tmp.SetParam(*in.GetParam());
  tmp.SetCovariance(*in.GetCovariance());

  if (verbose) {
    cout << "Track params. at 2D-DCA to beamline (vertex):" << endl;
    in.GetParam()->Print();
    in.GetCovariance()->Print();
    cout << "Propagated track params. at 2D-DCA to beamline (vertex):" << endl;
    tmp.GetParam()->Print();
    tmp.GetCovariance()->Print();
  }

  out = tmp;
  return true;
}

bool propagateKalman2first(MpdTpcKalmanTrack in, MpdTpcKalmanTrack &out)
{
  // Propagate kalman track to the inner TPC shell
  MpdTpcKalmanTrack tmp(in);
  bool verbose = false;
  tmp.SetDirection(MpdKalmanTrack::kOutward);

  /// Propagate TPC track to 27 cm inward (fPos == 27)
  MpdKalmanHit hitTmp;
  hitTmp.SetType(MpdKalmanHit::kFixedR);
  hitTmp.SetPos(tmp.GetPos()); // fPos ~= 27 cm
  
  tmp.SetParamNew(*tmp.GetParam());
  tmp.SetPos(tmp.GetPosNew());
  tmp.ReSetWeight();
  TMatrixDSym w = *tmp.GetWeight(); // save current weight matrix

  Bool_t ok = MpdKalmanFilter::Instance()->PropagateToHit(&tmp, &hitTmp, false, false, -1.0);
  if (!ok) {
    cerr << "WARNING: Could not propagate kalman track " << tmp.GetTrackID() << ". Skipping this one." << endl;
    out = tmp;
    return ok;
  }
  
  // After that GetCovariance() and GetParamNew() will give the values at the GetPosNew() point near the inner shell (R ~= 27 cm)
  TMatrixDSym cov(*tmp.Weight2Cov());
  TMatrixD param(*tmp.GetParamNew());
  
  // Resetting params and cov. matrix by the propagated ones ...
  tmp.SetParam(param);
  tmp.SetCovariance(cov);

  if (verbose) {
    cout << "Track params. at 2D-DCA to beamline (vertex):" << endl;
    in.GetParam()->Print();
    in.GetCovariance()->Print();
    cout << "Propagated track params. at the inner TPC shell (first):" << endl;
    tmp.GetParam()->Print();
    tmp.GetCovariance()->Print();
  }

  out = tmp;
  return true;
}

bool propagateKalman2last(MpdTpcKalmanTrack in, MpdTpcKalmanTrack &out)
{
  // Propagate kalman track to the last point (point of the last hit)
  MpdTpcKalmanTrack tmp(in);
  bool verbose = false;

  // Resetting params and cov. matrix to those ones  at the last hit ...
  tmp.SetParam(*in.GetParamAtHit());

  TMatrixDSym* weightAtLastHit = in.GetWeightAtHit();
  TMatrixDSym& covarAtLastHit = weightAtLastHit->InvertFast();

  tmp.SetCovariance(covarAtLastHit);

  if (verbose) {
    cout << "Track params. at 2D-DCA to beamline (vertex):" << endl;
    in.GetParam()->Print();
    in.GetCovariance()->Print();
    cout << "Propagated track params. at the last hit (last):" << endl;
    tmp.GetParam()->Print();
    tmp.GetCovariance()->Print();
  }

  out = tmp;
  return true;
}

RVec<MpdTpcKalmanTrack> getKalmanFirst(const RVec<MpdTpcKalmanTrack> &kalman_tracks)
try {
  RVec<MpdTpcKalmanTrack> tracks;
  for (auto &kalman_track:kalman_tracks) {
    MpdTpcKalmanTrack track;
    propagateKalman2first(kalman_track, track);
    tracks.push_back(track);
  }
  return tracks;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<MpdTpcKalmanTrack> getKalmanLast(const RVec<MpdTpcKalmanTrack> &kalman_tracks)
try {
  RVec<MpdTpcKalmanTrack> tracks;
  for (auto &kalman_track:kalman_tracks) {
    MpdTpcKalmanTrack track;
    propagateKalman2last(kalman_track, track);
    tracks.push_back(track);
  }
  return tracks;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector< vector<float> > covMatrix(const RVec<MpdTpcKalmanTrack> &kalman_tracks)
try {
  // Get covariance matrix in standard 6x6 format (x, y, z, px, py, pz)
  //     state vector in MPD: (r*phi, z, phi, lambda=pi/2-theta, -q/pT)
  //     Cout = J Cin J^T - A->B transformation
  vector<vector<float>> covariance_matrix;
  bool verbose = false;

  // Getting Jacobian ...
  float J[6][6];
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      J[i][j] = 0;

  for (auto &kalman_track:kalman_tracks) {
    auto *cov = kalman_track.GetCovariance();
    auto *par = kalman_track.GetParam();
    // Getting corresponding track parameters ...
    int q = kalman_track.Charge();
    double pt = kalman_track.Pt();
    double phi = (*par)(2, 0);
    double lambda = (*par)(3, 0);

    // Coordinate transformations ...
    J[0][0] = 1.; // dx / dx
    J[1][1] = 1.; // dy / dy
    J[2][2] = 1.; // dz / dz
    // Momentum transformations ...
    J[3][3] = -pt * sin(phi); // dPx / d\Phi
    J[3][5] = pt * pt * cos(phi) / q; // dPx / d(-q / Pt)
    J[4][3] = pt * cos(phi); // dPy / d\Phi
    J[4][5] = pt * pt * sin(phi) / q; // dPy / d(-q / Pt)
    J[5][4] = pt / (cos(lambda) * cos(lambda)); // dPz / dLambda
    J[5][5] = tan(lambda) * pt * pt / q; // dPz / d(-q / pt)

    // Extending track covariance matrix by one row and one column ...
    float CovIn[6][6]; // triangular -> symmetric matrix
    CovIn[0][0] = 1e-4; // dx. start from nowhere

    for (int i = 1; i < 6; i++) {
      CovIn[i][0] = 0;
      CovIn[0][i] = 0;
    }

    for (int i = 1; i < 6; i++) {
      for (int j = 1; j <= i; j++) {
        CovIn[i][j] = (*cov)(i - 1, j - 1);
        CovIn[j][i] = (*cov)(i - 1, j - 1);
      }
    }

    float CovInJt[6][6]; // CovInJt = CovIn * J^t
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) {
        CovInJt[i][j] = 0;
        for (int k = 0; k < 6; k++)
          CovInJt[i][j] += CovIn[i][k] * J[j][k];
      }
    
    float CovOut[6][6]; // CovOut = J * CovInJt
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++) {
        CovOut[i][j] = 0;
        for (int k = 0; k < 6; k++)
          CovOut[i][j] += J[i][k] * CovInJt[k][j];
      }

    // Lower triangle of the symmetric 6x6 covariance matrix (21 elements)
    // C[x, y, z, px, py, pz]
    // { c00, c1[0..1], c2[0..2], ... c4[0..5] }
    covariance_matrix.emplace_back();
    for (Int_t i = 0; i < 6; i++)
      for (Int_t j = 0; j <= i; j++)
        covariance_matrix.back().push_back(CovOut[i][j]);
  }

  if (verbose) {
    cout << "------------------------" << endl;
    cout << "Track cov. matrix (Cout): " << endl;
    for (auto &v:covariance_matrix) {
      for (auto &x:v) {
        printf("%f ", x);
      }
      cout << endl;
    }
    cout << "------------------------" << endl;
  }
  
  return covariance_matrix;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector< vector<float> > getKalmanParams(const RVec<MpdTpcKalmanTrack> &kalman_tracks)
try {
  // Get a set of parameters from the kalman track: (x,y,z,px,py,pz)
  vector<vector<float>> params;
  for (auto &kalman_track:kalman_tracks) {
    auto phi = kalman_track.GetParam(0) / kalman_track.GetPosNew();
    auto x = kalman_track.GetPosNew() * cos(phi);
    auto y = kalman_track.GetPosNew() * sin(phi);
    auto z = kalman_track.GetZ();
    auto px = kalman_track.Momentum3().X();
    auto py = kalman_track.Momentum3().Y();
    auto pz = kalman_track.Momentum3().Z();
    auto q = kalman_track.Charge();
    
    params.emplace_back();
    params.back().push_back( x );
    params.back().push_back( y );
    params.back().push_back( z );
    params.back().push_back( px );
    params.back().push_back( py );
    params.back().push_back( pz );
    params.back().push_back( q );
  }
  return params;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> kalmanNwrong(const RVec<MpdTpcKalmanTrack> &kalman_tracks)
try {
  vector<int> Nwrong;
  for (auto &track:kalman_tracks) {
    Nwrong.push_back(track.GetNofWrong());
  }
  return Nwrong;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> kalmanLength(const RVec<MpdTpcKalmanTrack> &kalman_tracks)
try {
  vector<float> length;
  for (auto &track:kalman_tracks) {
    length.push_back(track.GetLength());
  }
  return length;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> kalmanChi2(const RVec<MpdTpcKalmanTrack> &kalman_tracks)
try {
  vector<float> chi2;
  for (auto &track:kalman_tracks) {
    chi2.push_back(track.GetChi2());
  }
  return chi2;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

float determinant3x3( const std::array<std::array<float, 3>, 3>& matrix ) try {
  auto x_0 = matrix[0][0] * ( matrix[1][1]*matrix[2][2] - matrix[1][2]*matrix[2][1]  );
  auto x_1 = matrix[0][1] * ( matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0]  );
  auto x_2 = matrix[0][2] * ( matrix[1][0]*matrix[2][1] - matrix[1][1]*matrix[2][0]  );

  return x_0 - x_1 + x_2;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

std::array<float, 3> cramerFieldSolver3x3( std::array<float, 3> field, std::array<float, 3> coordinate ) try {
  // Solving the system of equation to extract parameters of quadratic extrapolation of the magnetic field
  // Ax = B
  // xi = detAi / detA
  std::array<std::array<float, 3>, 3> A;
  A[0] = {1.0f, 1.0f, 1.0f };
  A[1] = { coordinate[0], coordinate[1], coordinate[2] };
  A[2] = { coordinate[0]*coordinate[0], coordinate[1]*coordinate[1], coordinate[2]*coordinate[2] };

  auto A0 = A;
  A0[0] = field;
  auto A1 = A;
  A1[1] = field;
  auto A2 = A;
  A2[2] = field;

  auto detA = determinant3x3( A );
  auto detA0 = determinant3x3( A0 );
  auto detA1 = determinant3x3( A1 );
  auto detA2 = determinant3x3( A2 );

  auto p0 = detA0 / detA;
  auto p1 = detA1 / detA;
  auto p2 = detA2 / detA;

  return {p0, p1, p2};
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector< vector<float> > magneticField(const RVec<MpdTpcKalmanTrack> kalman_tracks)
try{
  vector<vector<float>> magnetic_field;
  for (auto &kalman_track:kalman_tracks) {
    std::array<float, 3> hit_z;
    std::array<float, 3> hit_bx;
    std::array<float, 3> hit_by;
    std::array<float, 3> hit_bz;

    for( int i=0; i<3; ++i ){
      // It seems size of the hitmap cannot be less than 4, but just to be safe
      if( i > kalman_track.GetHits()->GetEntriesFast() )
        magnetic_field.push_back( std::vector<float>(10, 0.0f) );

      auto hit = (MpdKalmanHit *)kalman_track.GetHits()->At(i);
      auto r = hit->GetDist();
      auto phi = hit->GetPhi();
      auto x = r*cos(phi);
      auto y = r*sin(phi);
      auto z = hit->GetMeas(1);

      hit_z.at(i) = z;
      hit_bx.at(i) = magField->GetBx( x, y, z ); // kGs
      hit_by.at(i) = magField->GetBy( x, y, z ); // kGs
      hit_bz.at(i) = magField->GetBz( x, y, z ); // kGs
    }

    auto parameters_bx = cramerFieldSolver3x3( hit_bx, hit_z );
    auto parameters_by = cramerFieldSolver3x3( hit_by, hit_z );
    auto parameters_bz = cramerFieldSolver3x3( hit_bz, hit_z );

    magnetic_field.emplace_back();
    for( const auto& c : parameters_bx )
      magnetic_field.back().push_back( c );
    for( const auto& c : parameters_by )
      magnetic_field.back().push_back( c );
    for( const auto& c : parameters_bz )
      magnetic_field.back().push_back( c );
    magnetic_field.back().push_back( 0.0 ); // z0
  }
  return magnetic_field;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> recSimIndex(const RVec<MpdTrack> &recoTracks, const RVec<MpdMCTrack> simTracks, const RVec<MpdMCTrack> &assocMcTracks, double mcvtxX, double mcvtxY, double mcvtxZ)
try {
  vector<int> newIndex;
  int shift=0;
  int nSimTracks = simTracks.size();
  for (int i=0;i<nSimTracks;i++) {
    if (!isGoodSimTrack(simTracks.at(i), assocMcTracks, mcvtxX, mcvtxY, mcvtxZ))
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
    if( idx > (int)sim_pdg.size() ) {
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
    if( idx > (int)sim_motherId.size() ) {
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

RVec<int> moduleId (const vector<XYZVector> &modulePos)
try {
  bool verbose = false;
  vector <int> moduleIds;
  for (int i=0;i<(int)modulePos.size();i++)
    moduleIds.push_back(i+1);
  if (verbose) {
    for(int i=0;i<(int)moduleIds.size();i++)
      printf("\t%d: (id = %d)\n", i, moduleIds.at(i));
  }
  return moduleIds;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> modulePos (const char *geoFile, const char *detectorTag)
{
  bool verbose = true;
  map <int,XYZVector> modulePosMap;
  printf("Reading %s geometry from geometry file\n", detectorTag);
  TFile *fiGeo = new TFile(geoFile, "read");
  TGeoVolumeAssembly *tgva = (TGeoVolumeAssembly *) fiGeo->Get("TOP");
  if (!tgva)
    throw runtime_error(Form("ERROR: No TGeoVolumeAssembly in file %s", geoFile));
  TGeoNode *topNode = tgva->GetNode(0);
  if (!topNode)
    throw runtime_error(Form("ERROR: No top node found in file %s", geoFile));
  // Declare nodes for left and right parts of the detector
  TGeoNode *det1Node = nullptr;
  TGeoNode *det2Node = nullptr;
  TString nodeName1, nodeName2;
  bool node1Found=false, node2Found=false;

  int nDetectors = topNode->GetNdaughters();
  if (nDetectors < 1)
    throw runtime_error(Form("ERROR: No detector nodes found in file %s", geoFile));
  if (nDetectors < 2)
    std::cerr << Form("WARNING: 2 detectors expected but only 1 was found in the top node!") << std::endl;
  for (int i = 0; i < nDetectors; i++) {
    det1Node = topNode->GetDaughter(i);
    nodeName1 = det1Node->GetName();
    nodeName1.ToLower();
    if (nodeName1.Contains(detectorTag)) {
      node1Found = true;
      for (int j = i+1; j < nDetectors; j++) {
        det2Node = topNode->GetDaughter(j);
        nodeName2 = det2Node->GetName();
        nodeName2.ToLower();
        if (nodeName2.Contains(detectorTag)) {
          node2Found = true;
          break;
        }
      }
      break;
    }
  }

  TVector3 frontFaceLocal, frontFaceGlobal;
  int nModules, nModules1, nModules2;
  nModules1 = 45; // default value
  nModules2 = 45; // default value
  if (node1Found) {
    auto geoMatrix = det1Node->GetMatrix();
    auto geoBox = (TGeoBBox*) det1Node->GetVolume()->GetShape();
    frontFaceLocal.SetXYZ(0, 0, -geoBox->GetDZ());
    geoMatrix->LocalToMaster(&frontFaceLocal[0], &frontFaceGlobal[0]);
    printf("%s node name: %s\n", detectorTag, nodeName1.Data());

    nModules1 = det1Node->GetNdaughters();
    for (int i = 0; i < nModules1; ++i) {
      auto *daughter = det1Node->GetDaughter(i);
      auto geoMatrix = daughter->GetMatrix();
      TVector3 translation(geoMatrix->GetTranslation());
      int modId = daughter->GetNumber();
      double x  = translation.X();
      double y  = translation.Y();
      translation.SetZ(frontFaceGlobal.Z());
      double z  = translation.Z();
      if (verbose) printf("\tModule %i, position: (%2.1f, %2.1f, %2.1f)\n", modId, x, y, z);
      modulePosMap.insert({modId, {x,y,z}});
    }
  }
  if (node2Found) {
    auto geoMatrix = det2Node->GetMatrix();
    auto geoBox = (TGeoBBox*) det2Node->GetVolume()->GetShape();
    frontFaceLocal.SetXYZ(0, 0, -geoBox->GetDZ());
    geoMatrix->LocalToMaster(&frontFaceLocal[0], &frontFaceGlobal[0]);
    printf("%s node name: %s\n", detectorTag, nodeName2.Data());

    nModules2 = det2Node->GetNdaughters();
    for (int i = 0; i < nModules2; ++i) {
      auto *daughter = det2Node->GetDaughter(i);
      auto geoMatrix = daughter->GetMatrix();
      TVector3 translation(geoMatrix->GetTranslation());
      int modId = daughter->GetNumber();
      double x  = translation.X();
      double y  = translation.Y();
      translation.SetZ(frontFaceGlobal.Z());
      double z  = translation.Z();
      if (verbose) printf("\tModule %i, position: (%3.1f, %3.1f, %3.1f)\n", modId + nModules1 + 1, x, y, z);
      modulePosMap.insert({modId + nModules1 + 1, {x,y,z}});
    }
  }

  fiGeo->Close();
  nModules = modulePosMap.rbegin()->first;
    vector <XYZVector> modulePosVector(nModules,{0.,0.,0.});
  for(auto &modulePos:modulePosMap)
    modulePosVector.at(modulePos.first-1)=modulePos.second;
  if (verbose)
  {
    printf("%d module positions:\n", nModules);
    for(int i=0;i<nModules;i++)
      printf("%d: (%3.1f, %3.1f, %3.1f)\n", i, modulePosVector.at(i).x(), modulePosVector.at(i).y(), modulePosVector.at(i).z());
  }
  return modulePosVector;
}

RVec<float> fhcalModE(const RVec<MpdZdcDigi> &fhcalHits, RVec<int> moduleIds)
try {
  bool verbose = false;
  const int nModules = moduleIds.size();
  vector<float> fhcalModEnergy(nModules, 0.);
  for (auto hit:fhcalHits) {
    int det_id = hit.GetDetectorID();
    int mod_id = hit.GetModuleID()-1;
    int i_module = mod_id + (nModules/2) * (det_id-1);
    fhcalModEnergy.at(i_module) += hit.GetELoss();
  }
  if (verbose) {
    for(int i=0;i<nModules;i++)
      printf("\t%d: (E = %f)\n", i, fhcalModEnergy.at(i));
  }
  return fhcalModEnergy;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<int> emcClusterMult(const TObjArray &emcClusters)
//  number of towers in the cluster.
// Recommended cut: mult > 1
try{
  vector<int> mults;
  for (int i=0; i<emcClusters.GetEntries(); i++) {
    auto cluster = (MpdEmcClusterKI*)emcClusters.At(i);
    mults.push_back(cluster->GetMultiplicity());
  }
  return mults;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> emcClusterEnergy(const TObjArray &emcClusters)
//  reconstructed energy of the cluster
// Recommended cut: e_corrected > 0.05 GeV
try{
  vector<float> energy;
  float conv = 1.0/0.3065; // only ~30% energy is collected in EMC
  for (int i=0; i<emcClusters.GetEntries(); i++) {
    auto cluster = (MpdEmcClusterKI*)emcClusters.At(i);
    float e_raw = cluster->GetE();
    float e1 = e_raw*conv;
    float x = e1;
    if (x > 2.5) x = 2.5; // parametrization is available for up to 2.5 GeV
    float corr[6] = { -1.321760e-002, 8.475623e-002, -9.007282e-002, 4.742396e-002, -1.279561e-002, 1.397394e-003 };
    float corr_coeff = corr[0] + corr[1]*x + corr[2]*x*x + corr[3]*x*x*x + corr[4]*x*x*x*x + corr[5]*x*x*x*x*x;
    float e_corrected = e1/(1. + corr_coeff); // corrected on non-linearity
    energy.push_back(e_corrected);
  }
  return energy;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

vector<XYZVector> emcClusterPos(const TObjArray &emcClusters)
//  XYZ position of the cluster in EMC in the global coordinate system.
// Recommended: use XYZ coordinates to calculate pz,py,pz
try{
  vector<XYZVector> pos;
  for (int i=0; i<emcClusters.GetEntries(); i++) {
    auto cluster = (MpdEmcClusterKI*)emcClusters.At(i);
    pos.push_back({cluster->GetX(), cluster->GetY(), cluster->GetZ()});
  }
  return pos;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> emcClusterChi2(const TObjArray &emcClusters)
//  Chi2/ndf of the cluster in EMC.
// Recommended cuts: chi2<4 for photons
try{
  vector<float> chi;
  for (int i=0; i<emcClusters.GetEntries(); i++) {
    auto cluster = (MpdEmcClusterKI*)emcClusters.At(i);
    chi.push_back(cluster->GetChi2());
  }
  return chi;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> emcClusterTime(const TObjArray &emcClusters)
//  simulated time of the shower.
// Recommended cuts: chi2<4 for photons
try{
  vector<float> time;
  TRandom RND;
  for (int i=0; i<emcClusters.GetEntries(); i++) {
    auto cluster = (MpdEmcClusterKI*)emcClusters.At(i);
    // smear time with 0.5 ns gaus
    time.push_back(RND.Gaus(cluster->GetTime()));
  }
  return time;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> emcClusterDPhi(const TObjArray &emcClusters)
//  distance from the cluster to the closets TPC track in dPhi.
try{
  vector<float> dphi;
  for (int i=0; i<emcClusters.GetEntries(); i++) {
    auto cluster = (MpdEmcClusterKI*)emcClusters.At(i);
    dphi.push_back(cluster->GetDPhi());
  }
  return dphi;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

RVec<float> emcClusterDZ(const TObjArray &emcClusters)
//  distance from the cluster to the closets TPC track in dZ.
try{
  vector<float> dz;
  for (int i=0; i<emcClusters.GetEntries(); i++) {
    auto cluster = (MpdEmcClusterKI*)emcClusters.At(i);
    dz.push_back(cluster->GetDZ());
  }
  return dz;
} catch( const std::exception& e ){
  std::cout << __func__ << std::endl;
  throw e;
}

void convertMPD(string inDst="", string fileOut="", string inGeo="")
{
  TStopwatch timer;
  timer.Start();

  cout << "Creating additional dictionary for IO..." << endl;
  gInterpreter->GenerateDictionary("vector<vector<float>>");
  
  TChain *chainRec=makeChain(inDst, "mpdsim");
  ROOT::RDataFrame d(*chainRec);

  int nEvents = chainRec->GetEntries();

  auto fhcalModPos = modulePos(inGeo.c_str(), "zdc");

  //magField = new MpdFieldMap("B-field_v2", "A");
  magField = new MpdConstField();
  magField->SetField(0., 0., 5.); // values are in kG:  1T = 10kG
  magField->SetFieldRegion(-230, 230, -230, 230, -375, 375); // values in cm
  magField->Init();
  magField->Print();

  auto dd = d
    .Define("evId", [](unsigned int id){ return (int)id; }, {"MCEventHeader.fEventId"})
    .Define("recoPrimVtxX", {"MPDEvent.PrimaryVerticesX"})
    .Define("recoPrimVtxY", {"MPDEvent.PrimaryVerticesY"})
    .Define("recoPrimVtxZ", {"MPDEvent.PrimaryVerticesZ"})
    .Define("recoPrimVtxChi2", {"MPDEvent.PrimaryVerticesChi2"})
    .Define("recoVtxNtracks", "Vertex.fNTracks")
    .Define("mcVtxX", {"MCEventHeader.fX"})
    .Define("mcVtxY", {"MCEventHeader.fY"})
    .Define("mcVtxZ", {"MCEventHeader.fZ"})
    .Define("mcB", { "MCEventHeader.fB"})
    .Define("mcRP", { "MCEventHeader.fRotZ"})
    .Define("recoGlobalTracks", getGlobalTracks, {"MPDEvent."})
    .Define("recoGlobalMom", trackMomentum, {"recoGlobalTracks"})
    .Define("recoGlobalNhits", recNhits, {"recoGlobalTracks"})
    .Define("recoGlobalNhitsFit", recNhitsFit, {"recoGlobalTracks"})
    .Define("recoGlobalNhitsPoss", recNhitsPoss, {"recoGlobalTracks"})
    .Define("recoGlobalChi2",  recChi2, {"recoGlobalTracks"})
    //.Define("recoGlobalP", trackP, {"recoGlobalTracks"})
    .Define("recoGlobalTofFlag", recHasTofHit, {"recoGlobalTracks"})
    .Define("recoGlobalCharge",recCharge,{"recoGlobalTracks"})
    .Define("recoGlobalDca", recDca, {"recoGlobalTracks"})
    .Define("recoGlobalTofMass2", recTofMass2, {"recoGlobalTracks"})
    .Define("recoGlobalTpcDedx", recTpcDedx, {"recoGlobalTracks"})
    //.Define("recoGlobalParamFirst", trGlobalFirstParam, {"recoGlobalTracks"})
    //.Define("recoGlobalParamLast", trGlobalLastParam, {"recoGlobalTracks"})
    //.Define("recoKalmanFirst", getKalmanFirst, {"TpcKalmanTrack"})
    //.Define("recoKalmanLast", getKalmanLast, {"TpcKalmanTrack"})
    //.Define("recoKalmanParamFirst", getKalmanParams, {"recoKalmanFirst"})
    //.Define("recoKalmanParamLast", getKalmanParams, {"recoKalmanLast"})
    .Define("recoKalmanParamVertex", getKalmanParams, {"TpcKalmanTrack"})
    .Define("recoKalmanCovMtxVertex", covMatrix, {"TpcKalmanTrack"})
    .Define("recoKalmanNofWrong", kalmanNwrong, {"TpcKalmanTrack"})
    .Define("recoKalmanLength", kalmanLength, {"TpcKalmanTrack"})
    .Define("recoKalmanCh2Ndf", kalmanChi2, {"TpcKalmanTrack"})
    .Define("recoKalmanMagField", magneticField, {"TpcKalmanTrack"})
    .Define("simAssocTracks", getAssocSimTracks, {"recoGlobalTracks", "MCTrack"})
    .Define("simMom", simMomentum, {"MCTrack", "simAssocTracks", "MCEventHeader.fX", "MCEventHeader.fY", "MCEventHeader.fZ"})
    .Define("simPosStart", simPosStart, {"MCTrack", "simAssocTracks", "MCEventHeader.fX", "MCEventHeader.fY", "MCEventHeader.fZ"})
    .Define("simMotherId", simMotherId, {"MCTrack", "simAssocTracks", "MCEventHeader.fX", "MCEventHeader.fY", "MCEventHeader.fZ"})
    .Define("simPdg", simPdg, {"MCTrack", "simAssocTracks", "MCEventHeader.fX", "MCEventHeader.fY", "MCEventHeader.fZ"})
    .Define("simCharge", simCharge, {"simPdg"})
    .Define("simHasHitFHCal", hasHitFhcal, {"MCTrack", "simAssocTracks", "MCEventHeader.fX", "MCEventHeader.fY", "MCEventHeader.fZ"})
    .Define("recoGlobalSimIndex", recSimIndex, {"recoGlobalTracks", "MCTrack", "simAssocTracks", "MCEventHeader.fX", "MCEventHeader.fY", "MCEventHeader.fZ"})
    .Define("recoGlobalSimPdg", trSimPdg, {"recoGlobalSimIndex", "simPdg"})
    .Define("recoGlobalSimMotherId", trSimMotherId, {"recoGlobalSimIndex", "simMotherId"})
    .Define("fhcalModPos", [fhcalModPos](){return fhcalModPos; })
    .Define("fhcalModId", moduleId, {"fhcalModPos"})
    .Define("fhcalModE",fhcalModE,{"ZdcDigi", "fhcalModId"})
    .Define("emcMult", emcClusterMult,    {"EmcCluster"})
    .Define("emcEnergy", emcClusterEnergy,{"EmcCluster"})
    .Define("emcPos", emcClusterPos,      {"EmcCluster"})
    .Define("emcChi2", emcClusterChi2,    {"EmcCluster"})
    .Define("emcTime", emcClusterTime,    {"EmcCluster"})
    .Define("emcDPhi", emcClusterDPhi,    {"EmcCluster"})
    .Define("emcDZ", emcClusterDZ,        {"EmcCluster"})
  ;
  dd.Foreach([](ULong64_t evtId){if (evtId % 100 == 0) cout << "\r" << evtId;}, {"rdfentry_"}); // progress display 
  cout << endl;

  vector<string> definedNames;
  vector<string> toExclude={"recoGlobalTracks", "simAssocTracks", "recoKalmanFirst", "recoKalmanLast"};
  for (auto& definedName:dd.GetDefinedColumnNames())
  {
    bool exclude=false;
    for (auto &nameToExclude:toExclude)
      if (definedName==nameToExclude)
        exclude=true;
    if (!exclude)
      definedNames.push_back(definedName);
  }
  auto n_events_filtered = *(dd.Count());
  std::cout << "Number of filtered events: " << n_events_filtered << std::endl;
  if( n_events_filtered > 0 )
    dd.Snapshot("t", fileOut, definedNames);

  cout << "Removing generated additional dictionaries..." << endl;
  gSystem->Exec("rm AutoDict_*");
  
  timer.Stop();
  timer.Print();
}
