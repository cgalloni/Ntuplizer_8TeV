#include "../interface/TriggersNtuplizer.h"

//===================================================================================================================        
TriggersNtuplizer::TriggersNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , HLTtriggersLabel_( labels[0] )
{
   
}

//===================================================================================================================
TriggersNtuplizer::~TriggersNtuplizer( void )
{

}

//===================================================================================================================
void TriggersNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

   /* here we want to save triggers info*/ 
  event.getByLabel(HLTtriggersLabel_, HLTtriggers_);
  
  edm::TriggerResults tr = *HLTtriggers_;        
  for( unsigned int t = 0; t < HLTtriggers_->size(); ++t){   	   
    nBranches_->HLTtrig.push_back(tr.accept(t)); 
  }
  bool isFired_HLT_HT650 = false;
  bool isFired_HLT_PFJet320 = false;
  bool isFired_HLT = false;  
  
  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
  unsigned int TrggIndex_PFHT650_v5( trigNames.triggerIndex("HLT_PFHT650_v5"));
  unsigned int TrggIndex_PFHT650_v6( trigNames.triggerIndex("HLT_PFHT650_v6"));
  unsigned int TrggIndex_PFHT650_v7( trigNames.triggerIndex("HLT_PFHT650_v7"));
  unsigned int TrggIndex_PFHT650_v8( trigNames.triggerIndex("HLT_PFHT650_v8"));
  unsigned int TrggIndex_PFHT650_v9( trigNames.triggerIndex("HLT_PFHT650_v9"));
  unsigned int TrggIndex_PFNoPUHT650_v1( trigNames.triggerIndex("HLT_PFNoPUHT650_v1"));
  unsigned int TrggIndex_PFNoPUHT650_v3( trigNames.triggerIndex("HLT_PFNoPUHT650_v3"));
  unsigned int TrggIndex_PFNoPUHT650_v4( trigNames.triggerIndex("HLT_PFNoPUHT650_v4"));
  unsigned int TrggIndex_PFJet320_v3( trigNames.triggerIndex("HLT_PFJet320_v3"));
  unsigned int TrggIndex_PFJet320_v4( trigNames.triggerIndex("HLT_PFJet320_v4"));
  unsigned int TrggIndex_PFJet320_v5( trigNames.triggerIndex("HLT_PFJet320_v5"));
  unsigned int TrggIndex_PFJet320_v6( trigNames.triggerIndex("HLT_PFJet320_v6"));
  unsigned int TrggIndex_PFJet320_v8( trigNames.triggerIndex("HLT_PFJet320_v8"));
  unsigned int TrggIndex_PFJet320_v9( trigNames.triggerIndex("HLT_PFJet320_v9"));
  if(TrggIndex_PFHT650_v5 < HLTtriggers_->size()) isFired_HLT_HT650 = HLTtriggers_->accept(TrggIndex_PFHT650_v5);
  if(TrggIndex_PFHT650_v6 < HLTtriggers_->size()) isFired_HLT_HT650 = HLTtriggers_->accept(TrggIndex_PFHT650_v6);
  if(TrggIndex_PFHT650_v7 < HLTtriggers_->size()) isFired_HLT_HT650 = HLTtriggers_->accept(TrggIndex_PFHT650_v7);
  if(TrggIndex_PFHT650_v8 < HLTtriggers_->size()) isFired_HLT_HT650 = HLTtriggers_->accept(TrggIndex_PFHT650_v8);
  if(TrggIndex_PFHT650_v9 < HLTtriggers_->size()) isFired_HLT_HT650 = HLTtriggers_->accept(TrggIndex_PFHT650_v9);
  if(TrggIndex_PFNoPUHT650_v1 < HLTtriggers_->size()) isFired_HLT_HT650 = HLTtriggers_->accept(TrggIndex_PFNoPUHT650_v1);
  if(TrggIndex_PFNoPUHT650_v3 < HLTtriggers_->size()) isFired_HLT_HT650 = HLTtriggers_->accept(TrggIndex_PFNoPUHT650_v3);
  if(TrggIndex_PFNoPUHT650_v4 < HLTtriggers_->size()) isFired_HLT_HT650 = HLTtriggers_->accept(TrggIndex_PFNoPUHT650_v4);
  if(TrggIndex_PFJet320_v3 < HLTtriggers_->size()) isFired_HLT_PFJet320 = HLTtriggers_->accept(TrggIndex_PFJet320_v3);
  if(TrggIndex_PFJet320_v4 < HLTtriggers_->size()) isFired_HLT_PFJet320 = HLTtriggers_->accept(TrggIndex_PFJet320_v4);
  if(TrggIndex_PFJet320_v5 < HLTtriggers_->size()) isFired_HLT_PFJet320 = HLTtriggers_->accept(TrggIndex_PFJet320_v5);
  if(TrggIndex_PFJet320_v6 < HLTtriggers_->size()) isFired_HLT_PFJet320 = HLTtriggers_->accept(TrggIndex_PFJet320_v6);
  if(TrggIndex_PFJet320_v8 < HLTtriggers_->size()) isFired_HLT_PFJet320 = HLTtriggers_->accept(TrggIndex_PFJet320_v8);
  if(TrggIndex_PFJet320_v9 < HLTtriggers_->size()) isFired_HLT_PFJet320 = HLTtriggers_->accept(TrggIndex_PFJet320_v9);
  if(isFired_HLT_PFJet320 || isFired_HLT_HT650) isFired_HLT=true;
   
  nBranches_->HLT_Jet320=isFired_HLT_PFJet320;
  nBranches_->HLT_Jet650= isFired_HLT_HT650;
  nBranches_->HLT_JetTrig= isFired_HLT;
}
