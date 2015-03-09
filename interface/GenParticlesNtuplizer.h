#ifndef GenParticlesNtuplizer_H
#define GenParticlesNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class GenParticlesNtuplizer : public CandidateNtuplizer {

public:
  GenParticlesNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches );
  ~GenParticlesNtuplizer( void ); 

  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

private:
   edm::InputTag genParticlesLabel_;
   edm::Handle< std::vector<reco::GenParticle> > genParticles_;
      
};

#endif // GenParticlesNtuplizer_H
