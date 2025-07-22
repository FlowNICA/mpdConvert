#include <Rtypes.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TStopwatch.h>

#include "FairMCEventHeader.h"

R__ADD_INCLUDE_PATH($VMCWORKDIR)

void MpdDst2UniGen(std::string inFileName, std::string outFileName)
{
  if (inFileName == "") {
    std::cerr << "Input file path is empty! Abort..." << std::endl;
    return;
  }
  if (outFileName == "") {
    std::cerr << "Output file path is empty! Abort..." << std::endl;
    return;
  }

  TStopwatch timer;
  timer.Start();

  TChain *dstTree = new TChain("mpdsim");
  dstTree->Add(inFileName.c_str());
  
  // Activate branches
  MpdEvent *event{nullptr};
  dstTree->SetBranchAddress("MPDEvent.", &event);
  FairMCEventHeader *MCHeader{nullptr};
  dstTree->SetBranchAddress("MCEventHeader.", &MCHeader);
  TClonesArray *fMCTracks{nullptr};
  dstTree->SetBranchAddress("MCTrack", &fMCTracks);

  TFile *fo = new TFile(outFileName.c_str(), "recreate");
  URun uheader;
  UEvent *uevent = new UEvent();
  TTree *tree = new TTree("events", generator);
  tree->Branch("event", "UEvent", event);
  
  long events = dstTree->GetEntries();
  std::cout << " Number of events in the DST file = " << events << std::endl;

  for (long i=0; i<events; i++){
    dstTree->GetEntry(i);

    uevent->Clear();
    // fill event info:     evtID,                  impact par            psiRP        # of event steps  ev step number   ev step time
    uevent->SetParameters(MCHeader->GetEventID(), MCHeader->GetB(), MCHeader->GetRotZ(), 0               ,    0          , 0 );

    int mctracks = fMCTracks->GetEntries();
    for (int itr=0; itr<mctracks; itr++){
      int tr_counter = 0;
      MpdMCTrack *mctrack = (MpdMCTrack*)fMCTracks->UncheckedAt(itr);
      if (!mctrack) continue;
      if (mctrack->GetMotherId() != -1) continue;

      int status = 0, parent = 0, parentDecay = 0, mate = 0, decay = 0, child[2] = {0, 0};
      double weight = 1.;
      TLorentzVector position(mctrack->GetStartX(), mctrack->GetStartY(), mctrack->GetStartZ(), mctrack->GetStartT());
      TLorentzVector momentum(mctrack->GetPx(), mctrack->GetPy(), mctrack->GetPz(), mctrack->GetEnergy());

      UParticle uparticle(tr_counter, mctrack->GetPdgCode(), status, parent, parentDecay, mate, decay, child, momentum, position, weight);
      
      uevent->AddParticle(uparticle);
      tree->Fill();
      tr_counter++;
    }
  }
  uheader.Write();
  tree->Write();

}
