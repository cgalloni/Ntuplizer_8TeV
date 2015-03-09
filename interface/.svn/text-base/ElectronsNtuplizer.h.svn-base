#ifndef ElectronsNtuplizer_H
#define ElectronsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class ElectronsNtuplizer : public CandidateNtuplizer {

public:
  ElectronsNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches );
  ~ElectronsNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::InputTag eleLabel_        ;
   edm::InputTag verticesLabel_   ;
   edm::InputTag rhoLabel_        ;
   edm::InputTag beamspotLabel_   ;
   edm::InputTag conversionsLabel_;
   edm::InputTag eleNotBoostedLabel_;
   
   edm::Handle<std::vector<pat::Electron> >  electrons_  ;
   edm::Handle< std::vector<reco::Vertex> >  vertices_   ;
   edm::Handle< reco::ConversionCollection > conversions_; 
   edm::Handle< reco::BeamSpot >             beamspot_   ; 
   edm::Handle< double >                     rho_        ;
   edm::Handle<std::vector<pat::Electron> >  electronsNotBoosted_;

      
};

#endif // ElectronsNtuplizer_H
