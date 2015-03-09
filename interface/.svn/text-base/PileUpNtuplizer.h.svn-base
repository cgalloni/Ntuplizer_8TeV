#ifndef PileUpNtuplizer_H
#define PileUpNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class PileUpNtuplizer : public CandidateNtuplizer {

public:
  PileUpNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches );
  ~PileUpNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::InputTag pileUpLabel_; 
     
   edm::Handle<std::vector< PileupSummaryInfo > >  pileUpInfo_;
      
};

#endif // PileUpNtuplizer_H
