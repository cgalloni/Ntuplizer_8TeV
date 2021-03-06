#include "../interface/TausNtuplizer.h"

#include <cmath>

#include "TVector3.h"
#include "TLorentzVector.h"

//===================================================================================================================        
TausNtuplizer::TausNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , tausLabel_( labels[0] )
   , eleTausLabel_( labels[1] )
   , muTausLabel_( labels[2] )
   , rhoLabel_( labels[3] )
   , verticesLabel_( labels[4] )
  
{
   
}

//===================================================================================================================
TausNtuplizer::~TausNtuplizer( void )
{

}

//===================================================================================================================
void TausNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
 
    event.getByLabel(tausLabel_, taus_); 
    event.getByLabel(eleTausLabel_, eleTaus_); 
    event.getByLabel(muTausLabel_, muTaus_);  
    event.getByLabel(rhoLabel_, rho_);
    event.getByLabel(verticesLabel_, vertices_);

    for( size_t t = 0; t < taus_->size(); ++t ){

      pat::Tau tau = (*taus_)[t];

    
      nBranches_->lep_type   .push_back(tau.pdgId());
      nBranches_->lep_charge .push_back(tau.charge());
      nBranches_->lep_e      .push_back(tau.energy());
      nBranches_->lep_eta    .push_back(tau.eta());
      nBranches_->lep_etaTrack    .push_back(tau.eta());
      nBranches_->lep_mass   .push_back(tau.mass());
      nBranches_->lep_pt     .push_back(tau.pt());
      nBranches_->lep_phi    .push_back(tau.phi());
      nBranches_->lep_TauType.push_back(1);
     
      //    TLorentzVector tau4V;
      //       tau4V.SetPtEtaPhiE(tau.pt(),tau.eta(),tau.phi(),tau.energy());
      //       float p= pow((tau4V.Px()* tau4V.Px()+ 
      // 		    tau4V.Py() * tau4V.Py() + tau4V.Pz() * tau4V.Pz()),0.5);
      //       float d= p* 2.9 * pow(10,-13) / (1.78) *3*pow(10,11);
      //       std::cout<<" d " << d<< " p "<< p <<std::endl;
      
      //Single top
      double rho = *(rho_.product());
      float  deltaR = 0.3;
      double energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->lep_pfRhoCorrRelIso03.push_back((tau.chargedHadronIso() + std::max(0.,tau.neutralHadronIso() + tau.photonIso() - energy))/tau.pt());      
      nBranches_->lep_pfRhoCorrRelIso03Boost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - energy))/tau.pt()); 
      deltaR = 0.4;
      energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->lep_pfRhoCorrRelIso04.push_back((tau.chargedHadronIso() + std::max(0.,tau.neutralHadronIso() + tau.photonIso() - energy))/tau.pt());
      nBranches_->lep_pfRhoCorrRelIso04Boost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - energy))/tau.pt()); 
            
      nBranches_->lep_pfDeltaCorrRelIso     .push_back((tau.chargedHadronIso() + std::max(0., tau.neutralHadronIso() + tau.photonIso() - 0.5*tau.puChargedHadronIso()))/tau.pt());
      nBranches_->lep_pfRelIso              .push_back((tau.chargedHadronIso() + tau.neutralHadronIso()+ tau.photonIso())/tau.pt()) ; 
      nBranches_->lep_photonIso             .push_back(tau.photonIso());
      nBranches_->lep_neutralHadIso         .push_back(tau.neutralHadronIso());
      nBranches_->lep_chargedHadIso         .push_back(tau.chargedHadronIso());
      nBranches_->lep_trackIso              .push_back(tau.trackIso());      
      nBranches_->lep_pfDeltaCorrRelIsoBoost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0., tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - 0.5*tau.userIsolation(pat::PfPUChargedHadronIso)))/tau.pt());
      nBranches_->lep_pfRelIsoBoost         .push_back((tau.userIsolation(pat::PfChargedHadronIso) + tau.userIsolation(pat::PfNeutralHadronIso)+ tau.userIsolation(pat::PfGammaIso))/tau.pt()) ; 
      nBranches_->lep_photonIsoBoost        .push_back(tau.userIsolation(pat::PfGammaIso));
      nBranches_->lep_neutralHadIsoBoost    .push_back(tau.userIsolation(pat::PfNeutralHadronIso));
      nBranches_->lep_chargedHadIsoBoost    .push_back(tau.userIsolation(pat::PfChargedHadronIso));

      nBranches_->lep_normChi2         .push_back(-99);
   
      nBranches_->lep_trackerHits      .push_back(-99);
      nBranches_->lep_matchedStations  .push_back(-99);
      nBranches_->lep_pixelHits        .push_back(-99);
      nBranches_->lep_globalHits       .push_back(-99);
      nBranches_->lep_isGlobalMuon     .push_back(-99);
      nBranches_->lep_isLooseMuon      .push_back(-99);
      nBranches_->lep_isSoftMuon       .push_back(-99);
      nBranches_->lep_isTightMuon      .push_back(-99);
      nBranches_->lep_isHighPtMuon     .push_back(-99);
      nBranches_->lep_dxy	       .push_back(-99);
      nBranches_->lep_dz	       .push_back(-99);
      nBranches_->lep_idMVAtrig        .push_back(-99);
      nBranches_->lep_looseId	       .push_back(-99);
      nBranches_->lep_vetoId	       .push_back(-99);
      nBranches_->lep_tightId	       .push_back(-99);
      nBranches_->lep_trigTightId      .push_back(-99);
      nBranches_->lep_isHEEP           .push_back(-99);
      
          
      nBranches_->decayModeFindingNewDMs                     .push_back(tau.tauID("decayModeFindingNewDMs"  		    ));
      nBranches_->decayModeFindingOldDMs	             .push_back(tau.tauID("decayModeFindingOldDMs"  		    ));
      nBranches_->decayModeFinding	                     .push_back(tau.tauID("decayModeFinding"			    ));
      nBranches_->byLooseIsolation                           .push_back(tau.tauID("byLooseIsolation"			    ));
      nBranches_->byVLooseCombinedIsolationDeltaBetaCorr     .push_back(tau.tauID("byVLooseCombinedIsolationDeltaBetaCorr"  ));
      nBranches_->byLooseCombinedIsolationDeltaBetaCorr      .push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr"   ));
      nBranches_->byMediumCombinedIsolationDeltaBetaCorr     .push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr"  ));
      nBranches_->byTightCombinedIsolationDeltaBetaCorr      .push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr"   ));
      nBranches_->byCombinedIsolationDeltaBetaCorrRaw        .push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw"	    ));
      nBranches_->byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits" ));
      nBranches_->byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"    ));
      nBranches_->chargedIsoPtSum                            .push_back(tau.tauID("chargedIsoPtSum" 			    ));
      nBranches_->neutralIsoPtSum                            .push_back(tau.tauID("neutralIsoPtSum" 			    ));
      nBranches_->puCorrPtSum                                .push_back(tau.tauID("puCorrPtSum"				    ));
      nBranches_->byIsolationMVA3oldDMwoLTraw                .push_back(tau.tauID("byIsolationMVA3oldDMwoLTraw"		    )); 
      nBranches_->byVLooseIsolationMVA3oldDMwoLT             .push_back(tau.tauID("byVLooseIsolationMVA3oldDMwoLT"  	    ));
      nBranches_->byLooseIsolationMVA3oldDMwoLT              .push_back(tau.tauID("byLooseIsolationMVA3oldDMwoLT"		    ));
      nBranches_->byMediumIsolationMVA3oldDMwoLT             .push_back(tau.tauID("byMediumIsolationMVA3oldDMwoLT"  	    ));
      nBranches_->byTightIsolationMVA3oldDMwoLT              .push_back(tau.tauID("byTightIsolationMVA3oldDMwoLT"		    ));
      nBranches_->byVTightIsolationMVA3oldDMwoLT             .push_back(tau.tauID("byVTightIsolationMVA3oldDMwoLT"  	    ));
      nBranches_->byVVTightIsolationMVA3oldDMwoLT            .push_back(tau.tauID("byVVTightIsolationMVA3oldDMwoLT" 	    ));
      nBranches_->byIsolationMVA3oldDMwLTraw                 .push_back(tau.tauID("byIsolationMVA3oldDMwLTraw"		    ));
      nBranches_->byVLooseIsolationMVA3oldDMwLT              .push_back(tau.tauID("byVLooseIsolationMVA3oldDMwLT"		    ));
      nBranches_->byLooseIsolationMVA3oldDMwLT               .push_back(tau.tauID("byLooseIsolationMVA3oldDMwLT"		    ));
      nBranches_->byMediumIsolationMVA3oldDMwLT              .push_back(tau.tauID("byMediumIsolationMVA3oldDMwLT"		    ));
      nBranches_->byTightIsolationMVA3oldDMwLT               .push_back(tau.tauID("byTightIsolationMVA3oldDMwLT"		    ));
      nBranches_->byVTightIsolationMVA3oldDMwLT              .push_back(tau.tauID("byVTightIsolationMVA3oldDMwLT"		    ));
      nBranches_->byVVTightIsolationMVA3oldDMwLT             .push_back(tau.tauID("byVVTightIsolationMVA3oldDMwLT"  	    ));
      nBranches_->byIsolationMVA3newDMwoLTraw                .push_back(tau.tauID("byIsolationMVA3newDMwoLTraw"		    ));
      nBranches_->byVLooseIsolationMVA3newDMwoLT             .push_back(tau.tauID("byVLooseIsolationMVA3newDMwoLT"  	    ));
      nBranches_->byLooseIsolationMVA3newDMwoLT              .push_back(tau.tauID("byLooseIsolationMVA3newDMwoLT"		    ));
      nBranches_->byMediumIsolationMVA3newDMwoLT             .push_back(tau.tauID("byMediumIsolationMVA3newDMwoLT"  	    ));
      nBranches_->byTightIsolationMVA3newDMwoLT              .push_back(tau.tauID("byTightIsolationMVA3newDMwoLT"		    ));
      nBranches_->byVTightIsolationMVA3newDMwoLT             .push_back(tau.tauID("byVTightIsolationMVA3newDMwoLT"  	    ));
      nBranches_->byVVTightIsolationMVA3newDMwoLT            .push_back(tau.tauID("byVVTightIsolationMVA3newDMwoLT" 	    ));
      nBranches_->byIsolationMVA3newDMwLTraw                 .push_back(tau.tauID("byIsolationMVA3newDMwLTraw"		    ));
      nBranches_->byVLooseIsolationMVA3newDMwLT              .push_back(tau.tauID("byVLooseIsolationMVA3newDMwLT"		    ));
      nBranches_->byLooseIsolationMVA3newDMwLT               .push_back(tau.tauID("byLooseIsolationMVA3newDMwLT"		    ));
      nBranches_->byMediumIsolationMVA3newDMwLT              .push_back(tau.tauID("byMediumIsolationMVA3newDMwLT"		    ));
      nBranches_->byTightIsolationMVA3newDMwLT               .push_back(tau.tauID("byTightIsolationMVA3newDMwLT"		    ));
      nBranches_->byVTightIsolationMVA3newDMwLT              .push_back(tau.tauID("byVTightIsolationMVA3newDMwLT"		    ));
      nBranches_->byVVTightIsolationMVA3newDMwLT             .push_back(tau.tauID("byVVTightIsolationMVA3newDMwLT"  	    ));
      nBranches_->againstElectronLoose                       .push_back(tau.tauID("againstElectronLoose"			    ));
      nBranches_->againstElectronMedium                      .push_back(tau.tauID("againstElectronMedium"			    ));
      nBranches_->againstElectronTight                       .push_back(tau.tauID("againstElectronTight"			    ));
      nBranches_->againstElectronMVA5raw                     .push_back(tau.tauID("againstElectronMVA5raw"  		    ));
      nBranches_->againstElectronMVA5category                .push_back(tau.tauID("againstElectronMVA5category"		    ));
      nBranches_->againstElectronVLooseMVA5                  .push_back(tau.tauID("againstElectronVLooseMVA5"		    ));
      nBranches_->againstElectronLooseMVA5                   .push_back(tau.tauID("againstElectronLooseMVA5"		    ));
      nBranches_->againstElectronMediumMVA5                  .push_back(tau.tauID("againstElectronMediumMVA5"		    ));
      nBranches_->againstElectronTightMVA5                   .push_back(tau.tauID("againstElectronTightMVA5"		    ));
      nBranches_->againstElectronVTightMVA5                  .push_back(tau.tauID("againstElectronVTightMVA5"		    ));
      nBranches_->againstElectronDeadECAL                    .push_back(tau.tauID("againstElectronDeadECAL" 		    ));
      nBranches_->againstMuonLoose                           .push_back(tau.tauID("againstMuonLoose"			    ));
      nBranches_->againstMuonMedium                          .push_back(tau.tauID("againstMuonMedium"			    ));
      nBranches_->againstMuonTight                           .push_back(tau.tauID("againstMuonTight"			    ));
      nBranches_->againstMuonLoose2                          .push_back(tau.tauID("againstMuonLoose2"			    ));
      nBranches_->againstMuonMedium2                         .push_back(tau.tauID("againstMuonMedium2"			    ));
      nBranches_->againstMuonTight2                          .push_back(tau.tauID("againstMuonTight2"			    ));
      nBranches_->againstMuonLoose3                          .push_back(tau.tauID("againstMuonLoose3"			    ));
      nBranches_->againstMuonTight3                          .push_back(tau.tauID("againstMuonTight3"			    ));
      nBranches_->againstMuonMVAraw                          .push_back(tau.tauID("againstMuonMVAraw"			    ));
      nBranches_->againstMuonLooseMVA                        .push_back(tau.tauID("againstMuonLooseMVA"			    ));
      nBranches_->againstMuonMediumMVA                       .push_back(tau.tauID("againstMuonMediumMVA"			    ));
      nBranches_->againstMuonTightMVA                        .push_back(tau.tauID("againstMuonTightMVA"			    ));  
      int                  nCCand          =-99   ;
      int                  nNCand          =-99   ;
      int                  nGCand          =-99   ;
    
      
      float  vSubtrack_innerposX   =-99  ;
      float  vSubtrack_innerposY   =-99  ;
      float  vSubtrack_innerposZ   =-99  ;
      
  

      nCCand=tau.signalPFChargedHadrCands().size();
      nNCand=tau.signalPFNeutrHadrCands().size();
      nGCand=tau.signalPFGammaCands().size();


      //       if ( tau.signalPFChargedHadrCands().size()!=0){
      // 	std::cout <<"(*tau.signalPFChargedHadrCands())[0].pt"<<(tau.signalPFChargedHadrCands())[0]->pt()<< "vx "<<(tau.signalPFChargedHadrCands())[0]->vx()<<
      //  	  " dxy "<<( - (tau.signalPFChargedHadrCands())[0]->vx() * (tau.signalPFChargedHadrCands())[0]->py() + 
      // (tau.signalPFChargedHadrCands())[0]->vy() * (tau.signalPFChargedHadrCands())[0]->px() )/ (tau.signalPFChargedHadrCands())[0]->pt() <<
      // 	  std::endl;
      //       }
      
      // float dxy = pow(tau.flightLength().X() * tau.flightLength().X() + tau.flightLength().Y() * tau.flightLength().Y(),0.5 );
      
      //  std::cout <<"nlep " << nBranches_->nlep<< "dxy " <<dxy<< " " << pow(tau.secondaryVertexPos().X() * tau.secondaryVertexPos().X() +
      // tau.secondaryVertexPos().Y() * tau.secondaryVertexPos().Y(),0.5 )<< " " <<tau.dxy() <<std:: endl;
      

      //       std::cout <<"tau.secondaryVertexPos().X() "<<tau.secondaryVertexPos().X()<<" vertices_->at(0).position().X() "<< vertices_->at(0).position().X()
      //  <<"  tau.secondaryVertexPos().X() - vertices_->at(0).position().X()" << (tau.secondaryVertexPos().X() - vertices_->at(0).position().X() ) 
      //  <<" tau flight length X "<<tau.flightLength().X()<<" tau.hasSecondaryVertex() "<<tau.hasSecondaryVertex()<<std::endl;
     

      if(tau.hasSecondaryVertex()){
	vSubtrack_innerposX    =(tau.secondaryVertexPos().X() // - vertices_->at(0).position().X()
			       );     
	vSubtrack_innerposY    =(tau.secondaryVertexPos().Y()//  - vertices_->at(0).position().Y()
				 );
	vSubtrack_innerposZ    =(tau.secondaryVertexPos().Z()//  - vertices_->at(0).position().Z()
				 );
      }

      nBranches_->nCharCand                 .push_back(nCCand        );
      nBranches_->nNeuCand                  .push_back(nNCand        );
      nBranches_->nGamCand                  .push_back(nGCand        );

      nBranches_->leadTrack_dxy             .push_back(tau.dxy()     );
      nBranches_->leadTrack_dxySig          .push_back(tau.dxy_Sig() );

     
      nBranches_->secVtxX                   .push_back(vSubtrack_innerposX);
      nBranches_->secVtxY                   .push_back(vSubtrack_innerposY);
      nBranches_->secVtxZ                   .push_back(vSubtrack_innerposZ);
      
      nBranches_->decayMode                 .push_back(tau.decayMode()    );
     
	
      nBranches_->flightLengthX             .push_back(tau.flightLength().X() );
      nBranches_->flightLengthY             .push_back(tau.flightLength().Y() );
      nBranches_->flightLengthZ             .push_back(tau.flightLength().Z() );
    

      
      nBranches_->flightLength              .push_back(pow( tau.flightLength().X() *  tau.flightLength().X() 
							    + tau.flightLength().Y() * tau.flightLength().Y()
							    + tau.flightLength().Z() * tau.flightLength().Z(),0.5));					
      nBranches_->flightLengthSig	    .push_back(tau.flightLengthSig());
    
      //       std::cout << "tau.dxy "    <<   tau.dxy()<<
      // 	" tau.dxy_Sig "          <<  tau.dxy_Sig() <<
      // 	" If Tau has Sec Vtx = " <<  tau.hasSecondaryVertex()<< 
      // 	" flightLength XY "      << pow( tau.flightLength().X() *  tau.flightLength().X()+ tau.flightLength().Y() * tau.flightLength().Y() ,0.5) << 
      // 	" flight length XYZ "    <<  pow(tau.flightLength().X() * tau.flightLength().X() + tau.flightLength().Y() * tau.flightLength().Y()
      //  + tau.flightLength().Z() * tau.flightLength().Z(),0.5 )<< 
      // 	"  pow(tau.vx() * tau.vx()+ tau.vy() * tau.vy() ,0.5) " << pow(tau.vx() * tau.vx()+ tau.vy() * tau.vy() ,0.5) <<
      // 	" pow(tau.vx() * tau.vx()+ tau.vy() * tau.vy() + tau.vz() * tau.vz() ,0.5) " << pow(tau.vx() * tau.vx()
      //   + tau.vy() * tau.vy() + tau.vz() * tau.vz() ,0.5) << 
      // 	"tau Vtx less PV XYZ " <<  pow(pow(tau.vx() - vertices_->at(0).position().X(),2)+ pow( tau.vy()- vertices_->at(0).position().Y(),2)
      // + pow(tau.vz()- vertices_->at(0).position().Z(),2)  ,0.5) << std::endl;

      bool  isATauSelection=false;  
     
      
      if(// (tau.pt()>20) &&
	 (fabs(tau.eta())<2.4) &&
	 (tau.tauID("decayModeFindingNewDMs")>0.5) &&
	 (tau.tauID("againstMuonLoose")>0.5) &&
	 (tau.tauID("againstElectronLoose")>0.5) &&
	 (tau.tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) ){
	isATauSelection=true;}
      
      nBranches_->lep_isATauSelection.push_back(isATauSelection);
      nBranches_->lep_isAEleTauSelection.push_back(false);
      nBranches_->lep_isAMuTauSelection.push_back(false);
      nBranches_->lep_isASelection.push_back(isATauSelection);  
      nBranches_->lep_isCSelection.push_back(isATauSelection);
     
      nBranches_->lep_SemileptonicPFIso        .push_back(-99);
      nBranches_->lep_SemileptonicCorrPFIso    .push_back(-99);
      nBranches_->lep_DilPFIso                 .push_back(-99);
      nBranches_->lep_DilDETIso                .push_back(-99);

      nBranches_->lep_isAGlobal   .push_back(-99);
      nBranches_->lep_isATracker  .push_back(-99);
      nBranches_->lep_cktTrack_pt.push_back (-99); 
 }
  
  for( size_t t = 0; t < eleTaus_->size(); ++t ){
 
      pat::Tau eleTau = (*eleTaus_)[t];

              
      nBranches_->lep_type   .push_back(eleTau.pdgId());
      nBranches_->lep_charge .push_back(eleTau.charge());
      nBranches_->lep_e      .push_back(eleTau.energy());
      nBranches_->lep_eta    .push_back(eleTau.eta()); 
      nBranches_->lep_etaTrack    .push_back(eleTau.eta());
      nBranches_->lep_mass   .push_back(eleTau.mass());
      nBranches_->lep_pt     .push_back(eleTau.pt());
      nBranches_->lep_phi    .push_back(eleTau.phi());
      nBranches_->lep_TauType.push_back(2);

      //Single top
      double rho = *(rho_.product());
      float  deltaR = 0.3;
      double energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->lep_pfRhoCorrRelIso03.push_back((eleTau.chargedHadronIso() + std::max(0.,eleTau.neutralHadronIso() + eleTau.photonIso() - energy))/eleTau.pt());      
      nBranches_->lep_pfRhoCorrRelIso03Boost.push_back((eleTau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,eleTau.userIsolation(pat::PfNeutralHadronIso) + eleTau.userIsolation(pat::PfGammaIso) - energy))/eleTau.pt()); 
      deltaR = 0.4;
      energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->lep_pfRhoCorrRelIso04.push_back((eleTau.chargedHadronIso() + std::max(0.,eleTau.neutralHadronIso() + eleTau.photonIso() - energy))/eleTau.pt());
      nBranches_->lep_pfRhoCorrRelIso04Boost.push_back((eleTau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,eleTau.userIsolation(pat::PfNeutralHadronIso) + eleTau.userIsolation(pat::PfGammaIso) - energy))/eleTau.pt()); 
            
      nBranches_->lep_pfDeltaCorrRelIso     .push_back((eleTau.chargedHadronIso() + std::max(0., eleTau.neutralHadronIso() + eleTau.photonIso() - 0.5*eleTau.puChargedHadronIso()))/eleTau.pt());
      nBranches_->lep_pfRelIso              .push_back((eleTau.chargedHadronIso() + eleTau.neutralHadronIso()+ eleTau.photonIso())/eleTau.pt()) ; 
      nBranches_->lep_photonIso             .push_back(eleTau.photonIso());
      nBranches_->lep_neutralHadIso         .push_back(eleTau.neutralHadronIso());
      nBranches_->lep_chargedHadIso         .push_back(eleTau.chargedHadronIso());
      nBranches_->lep_trackIso              .push_back(eleTau.trackIso());            
      nBranches_->lep_pfDeltaCorrRelIsoBoost.push_back((eleTau.userIsolation(pat::PfChargedHadronIso) + std::max(0., eleTau.userIsolation(pat::PfNeutralHadronIso) + eleTau.userIsolation(pat::PfGammaIso) - 0.5*eleTau.userIsolation(pat::PfPUChargedHadronIso)))/eleTau.pt());
      nBranches_->lep_pfRelIsoBoost         .push_back((eleTau.userIsolation(pat::PfChargedHadronIso) + eleTau.userIsolation(pat::PfNeutralHadronIso)+ eleTau.userIsolation(pat::PfGammaIso))/eleTau.pt()) ; 
      nBranches_->lep_photonIsoBoost        .push_back(eleTau.userIsolation(pat::PfGammaIso));
      nBranches_->lep_neutralHadIsoBoost    .push_back(eleTau.userIsolation(pat::PfNeutralHadronIso));
      nBranches_->lep_chargedHadIsoBoost    .push_back(eleTau.userIsolation(pat::PfChargedHadronIso));
      
      nBranches_->lep_normChi2       .push_back(-99);
      nBranches_->lep_isGlobalMuon   .push_back(-99);
      nBranches_->lep_trackerHits    .push_back(-99);
      nBranches_->lep_matchedStations.push_back(-99);
      nBranches_->lep_pixelHits      .push_back(-99);
      nBranches_->lep_globalHits     .push_back(-99);
      nBranches_->lep_isLooseMuon    .push_back(-99);
      nBranches_->lep_isSoftMuon     .push_back(-99);
      nBranches_->lep_isTightMuon    .push_back(-99);
      nBranches_->lep_isHighPtMuon   .push_back(-99);
      nBranches_->lep_dxy	     .push_back(-99);
     
      nBranches_->lep_dz	     .push_back(-99);
      nBranches_->lep_idMVAtrig      .push_back(-99);
      nBranches_->lep_looseId	     .push_back(-99);
      nBranches_->lep_vetoId	     .push_back(-99);
      nBranches_->lep_tightId	     .push_back(-99);
      nBranches_->lep_trigTightId    .push_back(-99);
      nBranches_->lep_isHEEP         .push_back(-99);

      nBranches_->decayModeFindingNewDMs                     .push_back(eleTau.tauID("decayModeFindingNewDMs"		       ));
      nBranches_->decayModeFindingOldDMs	             .push_back(eleTau.tauID("decayModeFindingOldDMs"		       ));
      nBranches_->decayModeFinding	                     .push_back(eleTau.tauID("decayModeFinding"			       ));
      nBranches_->byLooseIsolation                           .push_back(eleTau.tauID("byLooseIsolation"			       ));
      nBranches_->byVLooseCombinedIsolationDeltaBetaCorr     .push_back(eleTau.tauID("byVLooseCombinedIsolationDeltaBetaCorr"      ));
      nBranches_->byLooseCombinedIsolationDeltaBetaCorr      .push_back(eleTau.tauID("byLooseCombinedIsolationDeltaBetaCorr"       ));
      nBranches_->byMediumCombinedIsolationDeltaBetaCorr     .push_back(eleTau.tauID("byMediumCombinedIsolationDeltaBetaCorr"      ));
      nBranches_->byTightCombinedIsolationDeltaBetaCorr      .push_back(eleTau.tauID("byTightCombinedIsolationDeltaBetaCorr"       ));
      nBranches_->byCombinedIsolationDeltaBetaCorrRaw        .push_back(eleTau.tauID("byCombinedIsolationDeltaBetaCorrRaw"         ));
      nBranches_->byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(eleTau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(eleTau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits" ));
      nBranches_->byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(eleTau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(eleTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"    ));
      nBranches_->chargedIsoPtSum                            .push_back(eleTau.tauID("chargedIsoPtSum"			       ));
      nBranches_->neutralIsoPtSum                            .push_back(eleTau.tauID("neutralIsoPtSum"			       ));
      nBranches_->puCorrPtSum                                .push_back(eleTau.tauID("puCorrPtSum"  			       ));
      nBranches_->byIsolationMVA3oldDMwoLTraw                .push_back(eleTau.tauID("byIsolationMVA3oldDMwoLTraw"  	       )); 
      nBranches_->byVLooseIsolationMVA3oldDMwoLT             .push_back(eleTau.tauID("byVLooseIsolationMVA3oldDMwoLT"	       ));
      nBranches_->byLooseIsolationMVA3oldDMwoLT              .push_back(eleTau.tauID("byLooseIsolationMVA3oldDMwoLT"	       ));
      nBranches_->byMediumIsolationMVA3oldDMwoLT             .push_back(eleTau.tauID("byMediumIsolationMVA3oldDMwoLT"	       ));
      nBranches_->byTightIsolationMVA3oldDMwoLT              .push_back(eleTau.tauID("byTightIsolationMVA3oldDMwoLT"	       ));
      nBranches_->byVTightIsolationMVA3oldDMwoLT             .push_back(eleTau.tauID("byVTightIsolationMVA3oldDMwoLT"	       ));
      nBranches_->byVVTightIsolationMVA3oldDMwoLT            .push_back(eleTau.tauID("byVVTightIsolationMVA3oldDMwoLT"	       ));
      nBranches_->byIsolationMVA3oldDMwLTraw                 .push_back(eleTau.tauID("byIsolationMVA3oldDMwLTraw"		       ));
      nBranches_->byVLooseIsolationMVA3oldDMwLT              .push_back(eleTau.tauID("byVLooseIsolationMVA3oldDMwLT"	       ));
      nBranches_->byLooseIsolationMVA3oldDMwLT               .push_back(eleTau.tauID("byLooseIsolationMVA3oldDMwLT" 	       ));
      nBranches_->byMediumIsolationMVA3oldDMwLT              .push_back(eleTau.tauID("byMediumIsolationMVA3oldDMwLT"	       ));
      nBranches_->byTightIsolationMVA3oldDMwLT               .push_back(eleTau.tauID("byTightIsolationMVA3oldDMwLT" 	       ));
      nBranches_->byVTightIsolationMVA3oldDMwLT              .push_back(eleTau.tauID("byVTightIsolationMVA3oldDMwLT"	       ));
      nBranches_->byVVTightIsolationMVA3oldDMwLT             .push_back(eleTau.tauID("byVVTightIsolationMVA3oldDMwLT"	       ));
      nBranches_->byIsolationMVA3newDMwoLTraw                .push_back(eleTau.tauID("byIsolationMVA3newDMwoLTraw"  	       ));
      nBranches_->byVLooseIsolationMVA3newDMwoLT             .push_back(eleTau.tauID("byVLooseIsolationMVA3newDMwoLT"	       ));
      nBranches_->byLooseIsolationMVA3newDMwoLT              .push_back(eleTau.tauID("byLooseIsolationMVA3newDMwoLT"	       ));
      nBranches_->byMediumIsolationMVA3newDMwoLT             .push_back(eleTau.tauID("byMediumIsolationMVA3newDMwoLT"	       ));
      nBranches_->byTightIsolationMVA3newDMwoLT              .push_back(eleTau.tauID("byTightIsolationMVA3newDMwoLT"	       ));
      nBranches_->byVTightIsolationMVA3newDMwoLT             .push_back(eleTau.tauID("byVTightIsolationMVA3newDMwoLT"	       ));
      nBranches_->byVVTightIsolationMVA3newDMwoLT            .push_back(eleTau.tauID("byVVTightIsolationMVA3newDMwoLT"	       ));
      nBranches_->byIsolationMVA3newDMwLTraw                 .push_back(eleTau.tauID("byIsolationMVA3newDMwLTraw"		       ));
      nBranches_->byVLooseIsolationMVA3newDMwLT              .push_back(eleTau.tauID("byVLooseIsolationMVA3newDMwLT"	       ));
      nBranches_->byLooseIsolationMVA3newDMwLT               .push_back(eleTau.tauID("byLooseIsolationMVA3newDMwLT" 	       ));
      nBranches_->byMediumIsolationMVA3newDMwLT              .push_back(eleTau.tauID("byMediumIsolationMVA3newDMwLT"	       ));
      nBranches_->byTightIsolationMVA3newDMwLT               .push_back(eleTau.tauID("byTightIsolationMVA3newDMwLT" 	       ));
      nBranches_->byVTightIsolationMVA3newDMwLT              .push_back(eleTau.tauID("byVTightIsolationMVA3newDMwLT"	       ));
      nBranches_->byVVTightIsolationMVA3newDMwLT             .push_back(eleTau.tauID("byVVTightIsolationMVA3newDMwLT"	       ));
      nBranches_->againstElectronLoose                       .push_back(eleTau.tauID("againstElectronLoose" 		       ));
      nBranches_->againstElectronMedium                      .push_back(eleTau.tauID("againstElectronMedium"		       ));
      nBranches_->againstElectronTight                       .push_back(eleTau.tauID("againstElectronTight" 		       ));
      nBranches_->againstElectronMVA5raw                     .push_back(eleTau.tauID("againstElectronMVA5raw"		       ));
      nBranches_->againstElectronMVA5category                .push_back(eleTau.tauID("againstElectronMVA5category"  	       ));
      nBranches_->againstElectronVLooseMVA5                  .push_back(eleTau.tauID("againstElectronVLooseMVA5"		       ));
      nBranches_->againstElectronLooseMVA5                   .push_back(eleTau.tauID("againstElectronLooseMVA5"		       ));
      nBranches_->againstElectronMediumMVA5                  .push_back(eleTau.tauID("againstElectronMediumMVA5"		       ));
      nBranches_->againstElectronTightMVA5                   .push_back(eleTau.tauID("againstElectronTightMVA5"		       ));
      nBranches_->againstElectronVTightMVA5                  .push_back(eleTau.tauID("againstElectronVTightMVA5"		       ));
      nBranches_->againstElectronDeadECAL                    .push_back(eleTau.tauID("againstElectronDeadECAL"		       ));
      nBranches_->againstMuonLoose                           .push_back(eleTau.tauID("againstMuonLoose"			       ));
      nBranches_->againstMuonMedium                          .push_back(eleTau.tauID("againstMuonMedium"			       ));
      nBranches_->againstMuonTight                           .push_back(eleTau.tauID("againstMuonTight"			       ));
      nBranches_->againstMuonLoose2                          .push_back(eleTau.tauID("againstMuonLoose2"			       ));
      nBranches_->againstMuonMedium2                         .push_back(eleTau.tauID("againstMuonMedium2"			       ));
      nBranches_->againstMuonTight2                          .push_back(eleTau.tauID("againstMuonTight2"			       ));
      nBranches_->againstMuonLoose3                          .push_back(eleTau.tauID("againstMuonLoose3"			       ));
      nBranches_->againstMuonTight3                          .push_back(eleTau.tauID("againstMuonTight3"			       ));
      nBranches_->againstMuonMVAraw                          .push_back(eleTau.tauID("againstMuonMVAraw"			       ));
      nBranches_->againstMuonLooseMVA                        .push_back(eleTau.tauID("againstMuonLooseMVA"  		       ));
      nBranches_->againstMuonMediumMVA                       .push_back(eleTau.tauID("againstMuonMediumMVA" 		       ));
      nBranches_->againstMuonTightMVA                        .push_back(eleTau.tauID("againstMuonTightMVA"  		       ));   
        
      int                  nCCand          =-99   ;
      int                  nNCand          =-99   ;
      int                  nGCand          =-99   ;
    
      
      
      float  vSubtrack_innerposX   =-99  ;
      float  vSubtrack_innerposY   =-99  ;
      float  vSubtrack_innerposZ   =-99  ;
     


      nCCand=eleTau.signalPFChargedHadrCands().size();
      nNCand=eleTau.signalPFNeutrHadrCands().size();
      nGCand=eleTau.signalPFGammaCands().size();



      if(eleTau.hasSecondaryVertex()){
	vSubtrack_innerposX    =(eleTau.secondaryVertexPos().X()//  - vertices_->at(0).position().X()
				 );     
	vSubtrack_innerposY    =(eleTau.secondaryVertexPos().Y() // - vertices_->at(0).position().Y()
				 );
	vSubtrack_innerposZ    =(eleTau.secondaryVertexPos().Z() // - vertices_->at(0).position().Z()
				 );
      }
      
      nBranches_->leadTrack_dxy             .push_back(eleTau.dxy());
      nBranches_->leadTrack_dxySig          .push_back(eleTau.dxy_Sig());


      nBranches_->nCharCand                 .push_back(nCCand        );
      nBranches_->nNeuCand                  .push_back(nNCand        );
      nBranches_->nGamCand                  .push_back(nGCand        );
						      
      nBranches_->secVtxX                   .push_back(vSubtrack_innerposX);
      nBranches_->secVtxY                   .push_back(vSubtrack_innerposY);
      nBranches_->secVtxZ                   .push_back(vSubtrack_innerposZ);
      
      nBranches_->decayMode                 .push_back(eleTau.decayMode()    );
      
      nBranches_->flightLengthX             .push_back(eleTau.flightLength().X() );
      nBranches_->flightLengthY             .push_back(eleTau.flightLength().Y() );
      nBranches_->flightLengthZ             .push_back(eleTau.flightLength().Z() );
      nBranches_->flightLength              .push_back(pow( eleTau.flightLength().X() * eleTau.flightLength().X()+ eleTau.flightLength().Y() * eleTau.flightLength().Y() + eleTau.flightLength().Z() * eleTau.flightLength().Z(),0.5));
      nBranches_->flightLengthSig	    .push_back(eleTau.flightLengthSig());
    

      bool  isAEleTauSelection=false;  
     

      if(// (eleTau.pt()>20) &&
	 (fabs(eleTau.eta())<2.4) &&
	 (eleTau.tauID("decayModeFindingNewDMs")>0.5) &&
	 (eleTau.tauID("againstMuonLoose")>0.5) &&
	 (eleTau.tauID("againstElectronLoose")>0.5) &&
	 (eleTau.tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) ){
	isAEleTauSelection=true;}
      
      nBranches_->lep_isATauSelection.push_back(false);
      nBranches_->lep_isAEleTauSelection.push_back(isAEleTauSelection);
      nBranches_->lep_isAMuTauSelection.push_back(false);
      nBranches_->lep_isASelection.push_back(isAEleTauSelection);  
      nBranches_->lep_isCSelection.push_back(isAEleTauSelection);

      nBranches_->lep_SemileptonicPFIso        .push_back(-99);
      nBranches_->lep_SemileptonicCorrPFIso    .push_back(-99);
      nBranches_->lep_DilPFIso                 .push_back(-99);
      nBranches_->lep_DilDETIso                .push_back(-99);
    
      nBranches_->lep_isAGlobal   .push_back(-99);
      nBranches_->lep_isATracker  .push_back(-99);
      nBranches_->lep_cktTrack_pt.push_back (-99);
  }
  
  for( size_t t = 0; t < muTaus_->size(); ++t ){

      pat::Tau muTau = (*muTaus_)[t];

      /*std::cout << "mutau " << t << "/" << muTaus->size();
      std::cout << ": pt = " << muTau.pt() << " - eta = " << muTau.eta() << " - phi = " << muTau.phi() << " - e = " << muTau.energy() << std::endl;*/
       
      nBranches_->lep_type   .push_back(muTau.pdgId());
      nBranches_->lep_charge .push_back(muTau.charge());
      nBranches_->lep_e      .push_back(muTau.energy());
      nBranches_->lep_eta    .push_back(muTau.eta()); 
      nBranches_->lep_etaTrack    .push_back(muTau.eta());
      nBranches_->lep_mass   .push_back(muTau.mass());
      nBranches_->lep_pt     .push_back(muTau.pt());
      nBranches_->lep_phi    .push_back(muTau.phi());
      nBranches_->lep_TauType.push_back(3);
     
      //Single top
      double rho = *(rho_.product());
      float  deltaR = 0.3;
      double energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->lep_pfRhoCorrRelIso03.push_back((muTau.chargedHadronIso() + std::max(0.,muTau.neutralHadronIso() + muTau.photonIso() - energy))/muTau.pt());      
      nBranches_->lep_pfRhoCorrRelIso03Boost.push_back((muTau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,muTau.userIsolation(pat::PfNeutralHadronIso) + muTau.userIsolation(pat::PfGammaIso) - energy))/muTau.pt()); 
      deltaR = 0.4;
      energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->lep_pfRhoCorrRelIso04.push_back((muTau.chargedHadronIso() + std::max(0.,muTau.neutralHadronIso() + muTau.photonIso() - energy))/muTau.pt());
      nBranches_->lep_pfRhoCorrRelIso04Boost.push_back((muTau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,muTau.userIsolation(pat::PfNeutralHadronIso) + muTau.userIsolation(pat::PfGammaIso) - energy))/muTau.pt()); 
            
      nBranches_->lep_pfDeltaCorrRelIso     .push_back((muTau.chargedHadronIso() + std::max(0., muTau.neutralHadronIso() + muTau.photonIso() - 0.5*muTau.puChargedHadronIso()))/muTau.pt());
      nBranches_->lep_pfRelIso              .push_back((muTau.chargedHadronIso() + muTau.neutralHadronIso()+ muTau.photonIso())/muTau.pt()) ; 
      nBranches_->lep_photonIso             .push_back(muTau.photonIso());
      nBranches_->lep_neutralHadIso         .push_back(muTau.neutralHadronIso());
      nBranches_->lep_chargedHadIso         .push_back(muTau.chargedHadronIso());
      nBranches_->lep_trackIso              .push_back(muTau.trackIso());            
      nBranches_->lep_pfDeltaCorrRelIsoBoost.push_back((muTau.userIsolation(pat::PfChargedHadronIso) + std::max(0., muTau.userIsolation(pat::PfNeutralHadronIso) + muTau.userIsolation(pat::PfGammaIso) - 0.5*muTau.userIsolation(pat::PfPUChargedHadronIso)))/muTau.pt());
      nBranches_->lep_pfRelIsoBoost         .push_back((muTau.userIsolation(pat::PfChargedHadronIso) + muTau.userIsolation(pat::PfNeutralHadronIso)+ muTau.userIsolation(pat::PfGammaIso))/muTau.pt()) ; 
      nBranches_->lep_photonIsoBoost        .push_back(muTau.userIsolation(pat::PfGammaIso));
      nBranches_->lep_neutralHadIsoBoost    .push_back(muTau.userIsolation(pat::PfNeutralHadronIso));
      nBranches_->lep_chargedHadIsoBoost    .push_back(muTau.userIsolation(pat::PfChargedHadronIso));
      
      nBranches_->lep_normChi2       .push_back(-99);
      nBranches_->lep_isGlobalMuon   .push_back(-99);
      nBranches_->lep_trackerHits    .push_back(-99);
      nBranches_->lep_matchedStations.push_back(-99);
      nBranches_->lep_pixelHits      .push_back(-99);
      nBranches_->lep_globalHits     .push_back(-99);
      nBranches_->lep_isLooseMuon    .push_back(-99);
      nBranches_->lep_isSoftMuon     .push_back(-99);
      nBranches_->lep_isTightMuon    .push_back(-99);
      nBranches_->lep_isHighPtMuon   .push_back(-99);
      nBranches_->lep_dxy	     .push_back(-99);
      nBranches_->lep_dz	     .push_back(-99);
      nBranches_->lep_idMVAtrig      .push_back(-99);
      nBranches_->lep_looseId	     .push_back(-99);
      nBranches_->lep_vetoId	     .push_back(-99);
      nBranches_->lep_tightId	     .push_back(-99);
      nBranches_->lep_trigTightId    .push_back(-99);
      nBranches_->lep_isHEEP         .push_back(-99);
      
      nBranches_->decayModeFindingNewDMs                     .push_back(muTau.tauID("decayModeFindingNewDMs"		      ));
      nBranches_->decayModeFindingOldDMs			         .push_back(muTau.tauID("decayModeFindingOldDMs"		      ));
      nBranches_->decayModeFinding	                         .push_back(muTau.tauID("decayModeFinding"			      ));
      nBranches_->byLooseIsolation                           .push_back(muTau.tauID("byLooseIsolation"			      ));
      nBranches_->byVLooseCombinedIsolationDeltaBetaCorr     .push_back(muTau.tauID("byVLooseCombinedIsolationDeltaBetaCorr"      ));
      nBranches_->byLooseCombinedIsolationDeltaBetaCorr      .push_back(muTau.tauID("byLooseCombinedIsolationDeltaBetaCorr"       ));
      nBranches_->byMediumCombinedIsolationDeltaBetaCorr     .push_back(muTau.tauID("byMediumCombinedIsolationDeltaBetaCorr"      ));
      nBranches_->byTightCombinedIsolationDeltaBetaCorr      .push_back(muTau.tauID("byTightCombinedIsolationDeltaBetaCorr"       ));
      nBranches_->byCombinedIsolationDeltaBetaCorrRaw        .push_back(muTau.tauID("byCombinedIsolationDeltaBetaCorrRaw"	      ));
      nBranches_->byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(muTau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(muTau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits" ));
      nBranches_->byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(muTau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(muTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"    ));
      nBranches_->chargedIsoPtSum                            .push_back(muTau.tauID("chargedIsoPtSum"			      ));
      nBranches_->neutralIsoPtSum                            .push_back(muTau.tauID("neutralIsoPtSum"			      ));
      nBranches_->puCorrPtSum                                .push_back(muTau.tauID("puCorrPtSum"				      ));
      nBranches_->byIsolationMVA3oldDMwoLTraw                .push_back(muTau.tauID("byIsolationMVA3oldDMwoLTraw"		      )); 
      nBranches_->byVLooseIsolationMVA3oldDMwoLT             .push_back(muTau.tauID("byVLooseIsolationMVA3oldDMwoLT"	      ));
      nBranches_->byLooseIsolationMVA3oldDMwoLT              .push_back(muTau.tauID("byLooseIsolationMVA3oldDMwoLT" 	      ));
      nBranches_->byMediumIsolationMVA3oldDMwoLT             .push_back(muTau.tauID("byMediumIsolationMVA3oldDMwoLT"	      ));
      nBranches_->byTightIsolationMVA3oldDMwoLT              .push_back(muTau.tauID("byTightIsolationMVA3oldDMwoLT" 	      ));
      nBranches_->byVTightIsolationMVA3oldDMwoLT             .push_back(muTau.tauID("byVTightIsolationMVA3oldDMwoLT"	      ));
      nBranches_->byVVTightIsolationMVA3oldDMwoLT            .push_back(muTau.tauID("byVVTightIsolationMVA3oldDMwoLT"	      ));
      nBranches_->byIsolationMVA3oldDMwLTraw                 .push_back(muTau.tauID("byIsolationMVA3oldDMwLTraw"		      ));
      nBranches_->byVLooseIsolationMVA3oldDMwLT              .push_back(muTau.tauID("byVLooseIsolationMVA3oldDMwLT" 	      ));
      nBranches_->byLooseIsolationMVA3oldDMwLT               .push_back(muTau.tauID("byLooseIsolationMVA3oldDMwLT"  	      ));
      nBranches_->byMediumIsolationMVA3oldDMwLT              .push_back(muTau.tauID("byMediumIsolationMVA3oldDMwLT" 	      ));
      nBranches_->byTightIsolationMVA3oldDMwLT               .push_back(muTau.tauID("byTightIsolationMVA3oldDMwLT"  	      ));
      nBranches_->byVTightIsolationMVA3oldDMwLT              .push_back(muTau.tauID("byVTightIsolationMVA3oldDMwLT" 	      ));
      nBranches_->byVVTightIsolationMVA3oldDMwLT             .push_back(muTau.tauID("byVVTightIsolationMVA3oldDMwLT"	      ));
      nBranches_->byIsolationMVA3newDMwoLTraw                .push_back(muTau.tauID("byIsolationMVA3newDMwoLTraw"		      ));
      nBranches_->byVLooseIsolationMVA3newDMwoLT             .push_back(muTau.tauID("byVLooseIsolationMVA3newDMwoLT"	      ));
      nBranches_->byLooseIsolationMVA3newDMwoLT              .push_back(muTau.tauID("byLooseIsolationMVA3newDMwoLT" 	      ));
      nBranches_->byMediumIsolationMVA3newDMwoLT             .push_back(muTau.tauID("byMediumIsolationMVA3newDMwoLT"	      ));
      nBranches_->byTightIsolationMVA3newDMwoLT              .push_back(muTau.tauID("byTightIsolationMVA3newDMwoLT" 	      ));
      nBranches_->byVTightIsolationMVA3newDMwoLT             .push_back(muTau.tauID("byVTightIsolationMVA3newDMwoLT"	      ));
      nBranches_->byVVTightIsolationMVA3newDMwoLT            .push_back(muTau.tauID("byVVTightIsolationMVA3newDMwoLT"	      ));
      nBranches_->byIsolationMVA3newDMwLTraw                 .push_back(muTau.tauID("byIsolationMVA3newDMwLTraw"		      ));
      nBranches_->byVLooseIsolationMVA3newDMwLT              .push_back(muTau.tauID("byVLooseIsolationMVA3newDMwLT" 	      ));
      nBranches_->byLooseIsolationMVA3newDMwLT               .push_back(muTau.tauID("byLooseIsolationMVA3newDMwLT"  	      ));
      nBranches_->byMediumIsolationMVA3newDMwLT              .push_back(muTau.tauID("byMediumIsolationMVA3newDMwLT" 	      ));
      nBranches_->byTightIsolationMVA3newDMwLT               .push_back(muTau.tauID("byTightIsolationMVA3newDMwLT"  	      ));
      nBranches_->byVTightIsolationMVA3newDMwLT              .push_back(muTau.tauID("byVTightIsolationMVA3newDMwLT" 	      ));
      nBranches_->byVVTightIsolationMVA3newDMwLT             .push_back(muTau.tauID("byVVTightIsolationMVA3newDMwLT"	      ));
      nBranches_->againstElectronLoose                       .push_back(muTau.tauID("againstElectronLoose"  		      ));
      nBranches_->againstElectronMedium                      .push_back(muTau.tauID("againstElectronMedium" 		      ));
      nBranches_->againstElectronTight                       .push_back(muTau.tauID("againstElectronTight"  		      ));
      nBranches_->againstElectronMVA5raw                     .push_back(muTau.tauID("againstElectronMVA5raw"		      ));
      nBranches_->againstElectronMVA5category                .push_back(muTau.tauID("againstElectronMVA5category"		      ));
      nBranches_->againstElectronVLooseMVA5                  .push_back(muTau.tauID("againstElectronVLooseMVA5"		      ));
      nBranches_->againstElectronLooseMVA5                   .push_back(muTau.tauID("againstElectronLooseMVA5"		      ));
      nBranches_->againstElectronMediumMVA5                  .push_back(muTau.tauID("againstElectronMediumMVA5"		      ));
      nBranches_->againstElectronTightMVA5                   .push_back(muTau.tauID("againstElectronTightMVA5"		      ));
      nBranches_->againstElectronVTightMVA5                  .push_back(muTau.tauID("againstElectronVTightMVA5"		      ));
      nBranches_->againstElectronDeadECAL                    .push_back(muTau.tauID("againstElectronDeadECAL"		      ));
      nBranches_->againstMuonLoose                           .push_back(muTau.tauID("againstMuonLoose"			      ));
      nBranches_->againstMuonMedium                          .push_back(muTau.tauID("againstMuonMedium"			      ));
      nBranches_->againstMuonTight                           .push_back(muTau.tauID("againstMuonTight"			      ));
      nBranches_->againstMuonLoose2                          .push_back(muTau.tauID("againstMuonLoose2"			      ));
      nBranches_->againstMuonMedium2                         .push_back(muTau.tauID("againstMuonMedium2"			      ));
      nBranches_->againstMuonTight2                          .push_back(muTau.tauID("againstMuonTight2"			      ));
      nBranches_->againstMuonLoose3                          .push_back(muTau.tauID("againstMuonLoose3"			      ));
      nBranches_->againstMuonTight3                          .push_back(muTau.tauID("againstMuonTight3"			      ));
      nBranches_->againstMuonMVAraw                          .push_back(muTau.tauID("againstMuonMVAraw"			      ));
      nBranches_->againstMuonLooseMVA                        .push_back(muTau.tauID("againstMuonLooseMVA"			      ));
      nBranches_->againstMuonMediumMVA                       .push_back(muTau.tauID("againstMuonMediumMVA"  		      ));
      nBranches_->againstMuonTightMVA                        .push_back(muTau.tauID("againstMuonTightMVA"			      ));         
   
      int    nCCand	           = -99  ;
      int    nNCand	           = -99  ;
      int    nGCand	           = -99  ;	       
      float  vSubtrack_innerposX   = -99  ;
      float  vSubtrack_innerposY   = -99  ;
      float  vSubtrack_innerposZ   = -99  ;
  
      nCCand=muTau.signalPFChargedHadrCands().size();
      nNCand=muTau.signalPFNeutrHadrCands().size();
      nGCand=muTau.signalPFGammaCands().size();

      // float dxy = pow(muTau.flightLength().X() * muTau.flightLength().X() + muTau.flightLength().Y() * muTau.flightLength().Y(),0.5 );
     
     //  std::cout <<"nlep " << nBranches_->nlep<< "dxy " <<dxy<< " " << pow(muTau.secondaryVertexPos().X() * muTau.secondaryVertexPos().X() + muTau.secondaryVertexPos().Y() * muTau.secondaryVertexPos().Y(),0.5 )<< " " <<muTau.dxy() <<std:: endl;


      if(muTau.hasSecondaryVertex()){
	vSubtrack_innerposX    =(muTau.secondaryVertexPos().X() // - vertices_->at(0).position().X()
				 );
	vSubtrack_innerposY    =(muTau.secondaryVertexPos().Y() // - vertices_->at(0).position().Y()
				 );
	vSubtrack_innerposZ    =(muTau.secondaryVertexPos().Z() // - vertices_->at(0).position().Z()
				 );
      } 

      nBranches_->leadTrack_dxy             .push_back(muTau.dxy());
      nBranches_->leadTrack_dxySig          .push_back(muTau.dxy_Sig());
      
      nBranches_->nCharCand                 .push_back(nCCand        );
      nBranches_->nNeuCand                  .push_back(nNCand        );
      nBranches_->nGamCand                  .push_back(nGCand        );
     			      
      nBranches_->secVtxX                   .push_back(vSubtrack_innerposX);
      nBranches_->secVtxY                   .push_back(vSubtrack_innerposY);
      nBranches_->secVtxZ                   .push_back(vSubtrack_innerposZ);
      nBranches_->decayMode                 .push_back(muTau.decayMode()    );
      nBranches_->flightLengthX             .push_back(muTau.flightLength().X() );
      nBranches_->flightLengthY             .push_back(muTau.flightLength().Y() );
      nBranches_->flightLengthZ             .push_back(muTau.flightLength().Z() );
      nBranches_->flightLength              .push_back(pow( muTau.flightLength().X() *  muTau.flightLength().X()+ muTau.flightLength().Y() * muTau.flightLength().Y() + muTau.flightLength().Z() * muTau.flightLength().Z(),0.5));
      nBranches_->flightLengthSig	    .push_back(muTau.flightLengthSig());
  

      //TAU SELECTION - MUTAU
      bool  isAMuTauSelection=false;

      if(// (muTau.pt()>20) &&
	 (fabs(muTau.eta())<2.4) &&
	 (muTau.tauID("decayModeFindingNewDMs")>0.5) &&
	 (muTau.tauID("againstMuonLoose")>0.5) &&
	 (muTau.tauID("againstElectronLoose")>0.5) &&
	 (muTau.tauID("byVLooseIsolationMVA3newDMwoLT")>0.5) ){
	isAMuTauSelection=true;}
      
      //   if(foundJet){
      //       // if(ROOT::Math::VectorUtil::DeltaR(patTau->p4(),SelectedJet->p4())>0.8){
      // // 	foundTauMuTau=true;
      // // 	if(patTau->pt()>ptTauMuTau){
      // // 	  SelectedTauMuTau=patTau;
      // // 	  ptTauMuTau=patTau->pt();
      // // 	}
      // //       }
      //     }
	
      nBranches_->lep_isATauSelection.push_back(false);
      nBranches_->lep_isAEleTauSelection.push_back(false);
      nBranches_->lep_isAMuTauSelection.push_back(isAMuTauSelection);
      nBranches_->lep_isASelection.push_back(isAMuTauSelection);
      nBranches_->lep_isCSelection.push_back(isAMuTauSelection);
   
      nBranches_->lep_SemileptonicPFIso        .push_back(-99);
      nBranches_->lep_SemileptonicCorrPFIso    .push_back(-99);
      nBranches_->lep_DilPFIso                 .push_back(-99);
      nBranches_->lep_DilDETIso                .push_back(-99);
      nBranches_->lep_isAGlobal   .push_back(-99);
      nBranches_->lep_isATracker  .push_back(-99);
      nBranches_->lep_cktTrack_pt.push_back (-99);
  }
    
    
    nBranches_->nlep +=  taus_->size() +eleTaus_->size()+ muTaus_->size();
}  
