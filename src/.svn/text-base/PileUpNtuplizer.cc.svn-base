#include "../interface/PileUpNtuplizer.h"

//===================================================================================================================
PileUpNtuplizer::PileUpNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , pileUpLabel_( labels[0] )
{

}

//===================================================================================================================
PileUpNtuplizer::~PileUpNtuplizer( void )
{

}

//===================================================================================================================
void PileUpNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByLabel(pileUpLabel_, pileUpInfo_);  
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  
  for( PVI = pileUpInfo_->begin(); PVI != pileUpInfo_->end(); ++PVI ) {
     nBranches_->nPuVtxTrue.push_back(PVI->getTrueNumInteractions());
     nBranches_->nPuVtx.push_back(PVI->getPU_NumInteractions());
     nBranches_->bX.push_back(PVI->getBunchCrossing());
     nBranches_->PuVtxZ.push_back(PVI->getPU_zpositions());
     nBranches_->sumpTlow.push_back(PVI->getPU_sumpT_lowpT());
     nBranches_->sumpThigh.push_back(PVI->getPU_sumpT_highpT());
     nBranches_->ntrksLowpT.push_back(PVI->getPU_ntrks_lowpT());
     nBranches_->ntrksHighpT.push_back(PVI->getPU_ntrks_highpT());
  }

}
