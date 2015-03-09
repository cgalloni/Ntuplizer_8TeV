#include "../interface/VerticesNtuplizer.h"

#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"

#include <cmath>

//===================================================================================================================
VerticesNtuplizer::VerticesNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , verticesLabel_( labels[0] )
{

}

//===================================================================================================================
VerticesNtuplizer::~VerticesNtuplizer( void )
{

}

//===================================================================================================================
void VerticesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByLabel(verticesLabel_, vertices_);
    
  nBranches_->nPVs = vertices_->size();
  
  for( unsigned int v = 0; v < vertices_->size(); ++v ){
  
     reco::Vertex vtx = (*vertices_)[v];
     
     nBranches_->PV_x 	 .push_back(vtx.x());
     nBranches_->PV_y 	 .push_back(vtx.y());
     nBranches_->PV_z 	 .push_back(vtx.z());
     nBranches_->PV_xerr .push_back(vtx.xError());
     nBranches_->PV_yerr .push_back(vtx.yError());
     nBranches_->PV_zerr .push_back(vtx.zError());
     nBranches_->PV_chi2 .push_back(vtx.chi2());
     nBranches_->PV_ndof .push_back(vtx.ndof());
     
  }
  
}
