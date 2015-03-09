#ifndef GenJetsNtuplizer_H
#define GenJetsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class GenJetsNtuplizer : public CandidateNtuplizer {

public:
   GenJetsNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches );
   ~GenJetsNtuplizer( void );
   
   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::InputTag genJetsAK5Label_    ;
   edm::InputTag genJetsAK5nonuLabel_;
   edm::InputTag genJetsCA8Label_    ;
   edm::InputTag genJetsCA8nonuLabel_;
  
   edm::Handle< std::vector<reco::GenJet> > genJetsAK5_    ;
   edm::Handle< std::vector<reco::GenJet> > genJetsAK5nonu_;
   edm::Handle< std::vector<reco::GenJet> > genJetsCA8_    ;
   edm::Handle< std::vector<reco::GenJet> > genJetsCA8nonu_; 
       
};

#endif // GenJetsNtuplizer_H
