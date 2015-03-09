#include "../interface/MuonsNtuplizer.h"

#include <cmath>

//===================================================================================================================        
MuonsNtuplizer::MuonsNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , muonsLabel_( labels[0] )
   , verticesLabel_( labels[1] )
   , rhoLabel_( labels[2] )
   , muonsNotBoostedLabel_( labels[3] )
{
   
}

//===================================================================================================================
MuonsNtuplizer::~MuonsNtuplizer( void )
{

}

//===================================================================================================================
float MuonPFIso(pat::Muon muon, bool highpt){
  float sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon.pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt + std::max( 0. ,sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt))// / muon.pt()
    ;
 //  if(highpt){
//     reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 17., 40., 0.25)).first;
//     iso = (sumChargedHadronPt+ std::max ( 0. , sumNeutralHadronEt + sumPhotonEt - 0.5* sumPUPt ) )/cktTrack->pt();
//   }
  return iso;
}

float MuonCorrPFIso(pat::Muon muon, bool highpt){
  float sumChargedHadronPt = muon.userIsolation(pat::PfChargedHadronIso);//pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon.userIsolation(pat::PfNeutralHadronIso);//pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon.userIsolation(pat::PfGammaIso);//pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon.userIsolation(pat::User2Iso);//pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt+ std::max( 0., sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt ))// /muon.pt()
    ;
  // if(highpt){ 
//     reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 17., 40., 0.25)).first;
//     iso = (sumChargedHadronPt+ std::max(0., sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt))/cktTrack->pt();
//   }
  return iso;
}

void SelectTrackerGlobalID(pat::Muon muon ,  reco::Vertex primaryVertex,  bool & hasAtLeastOneHighPtMuon){
  hasAtLeastOneHighPtMuon = false;
  if((muon.isGlobalMuon())) {
    reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 17., 40., 0.25)).first;
    if(// (cktTrack->pt()>=10)
//        &&
       (fabs(cktTrack->eta())<=2.4)  
       &&(fabs(cktTrack->phi())<=3.2)  
       &&((cktTrack->ptError()/cktTrack->pt())<=0.3)  ////////
       &&(muon.globalTrack()->hitPattern().numberOfValidMuonHits()>0)  
       &&(muon.numberOfMatches()>1)  
       &&(fabs(cktTrack->dxy(primaryVertex.position()))<0.2)  
       &&(fabs(cktTrack->dz(primaryVertex.position()))<0.5)  
       &&(muon.innerTrack()->hitPattern().numberOfValidPixelHits()>0)  
       &&(muon.innerTrack()->hitPattern().trackerLayersWithMeasurement()>5) )
      hasAtLeastOneHighPtMuon=true;
  }
}

// for(pat::MuonCollection::const_iterator muon = muoH->begin(); muon != muoH->end(); ++muon) {
//     if(!(muon->isGlobalMuon())) continue;
//     reco::TrackRef cktTrack = (muon::tevOptimized(*muon, 200, 17., 40., 0.25)).first;
//     if(cktTrack->pt()<10) continue;
//     if(abs(cktTrack->eta())>2.4) continue;
//     if(abs(cktTrack->phi())>3.2) continue;
//     if((cktTrack->ptError()/cktTrack->pt())>0.3) continue;
//     if(muon->globalTrack()->hitPattern().numberOfValidMuonHits()<=0) continue;
//     if(muon->numberOfMatches()<=1) continue;
//     if(fabs(cktTrack->dxy(primaryVertex.position()))>=0.2) continue;
//     if(fabs(cktTrack->dz( primaryVertex.position()))>=0.5) continue;
//     if(muon->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
//     if(muon->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
//     if(ROOT::Math::VectorUtil::DeltaR(muon->p4(),SelectedJet->p4())<0.8) continue;
// }






void SelectTrackerMuon(pat::Muon muon , reco::Vertex primaryVertex,  bool & hasAtTrackerMuon){
  hasAtTrackerMuon = false;
  if(// (muon.pt()>=10)  &&
     (fabs(muon.eta())<=2.4)  &&
     (fabs(muon.phi())<=3.2)  &&
     ((muon.isTrackerMuon()))  &&
     (muon.numberOfMatches()>1)  &&
     (fabs(muon.muonBestTrack()->dz(primaryVertex.position()))<0.5)  &&
     (fabs(muon.dB())<0.2 )  &&
     (muon.innerTrack()->hitPattern().numberOfValidPixelHits()>0)  &&
     (muon.innerTrack()->hitPattern().trackerLayersWithMeasurement()>5)  &&
     ((muon.muonBestTrack()->ptError()/muon.muonBestTrack()->pt())<0.3))//  &&(MuonDETIso<0.2)
    hasAtTrackerMuon = true;
}
  
  
  
//===================================================================================================================
void MuonsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

    event.getByLabel(muonsLabel_, muons_); 
    event.getByLabel(verticesLabel_, vertices_);
    event.getByLabel(rhoLabel_, rho_);
    event.getByLabel(muonsNotBoostedLabel_, muonsNotBoosted_); 

    //std::cout << "//////////////////////" << muons->size() << std::endl;
    for( size_t m = 0; m < muons_->size(); ++m ){
     
      pat::Muon muon = (*muons_)[m];
      pat::Muon muonNotBoosted = (*muonsNotBoosted_)[m];
      //std::cout << "muon " << m << "/" << muons->size();
      //std::cout << ": pt = " << muon.pt() << " - eta = " << muon.eta() << " - phi = " << muon.phi() << " - e = " << muon.energy() << std::endl;
            
      nBranches_->lep_type   .push_back(muon.pdgId() );   	 
      nBranches_->lep_charge .push_back(muon.charge());   	 
      nBranches_->lep_e      .push_back(muon.energy());   	 
      nBranches_->lep_eta    .push_back(muon.eta()   );
      nBranches_->lep_etaTrack    .push_back(muon.eta()   );   	 
      nBranches_->lep_mass   .push_back(muon.mass()  );   	 
      nBranches_->lep_pt     .push_back(muon.pt()    );   	 
      nBranches_->lep_phi    .push_back(muon.phi()   );
      nBranches_->lep_TauType.push_back(0);      
      	  	 
      float dxy = fabs(muon.muonBestTrack()->dxy(vertices_->at(0).position()));
      nBranches_->lep_dxy.push_back(dxy);
      float dz= fabs(muon.muonBestTrack()->dz(vertices_->at(0).position()));
      nBranches_->lep_dz.push_back(dz);

      nBranches_->lep_isGlobalMuon.push_back(muon.isGlobalMuon());      
      nBranches_->lep_isLooseMuon .push_back(muon.isLooseMuon()); 
      nBranches_->lep_isSoftMuon  .push_back(muon.isSoftMuon(vertices_->at(0)));  
      nBranches_->lep_isTightMuon .push_back(muon.isTightMuon(vertices_->at(0))); 
      nBranches_->lep_isHighPtMuon.push_back(muon.isHighPtMuon(vertices_->at(0),reco::improvedTuneP));

   
//       std::cout<< "muon.pt() " << muon.pt() << std::endl;
//       std::cout<<" muon.eta() "<< muon.eta()<<std::endl;
//       std::cout<< "muon.phi() " << muon.phi() << std::endl;


   /*===== ISO ====*/
      //Muon POG
      //const reco::MuonPFIsolation& pfIsolationR03 = muon.pfIsolationR03();
      //const reco::MuonPFIsolation& pfIsolationR04 = muon.pfIsolationR04();
      
      //VHbb
      //pfRhoCorrRelIso = (muon.chargedHadronIso() + std::max(0., muon.neutralHadronIso() + muon.photonIso() - 0.5*muon.puChargedHadronIso()))/muon.pt();
      
      //Single Top
      double rho = *(rho_.product()); 
      float  deltaR = 0.3;
      double energy = TMath::Pi()*deltaR*deltaR*rho;
      nBranches_->lep_pfRhoCorrRelIso03.push_back((muon.chargedHadronIso() + std::max(0., muon.neutralHadronIso() + muon.photonIso() - energy))/muon.pt());
      nBranches_->lep_pfRhoCorrRelIso03Boost.push_back((muon.userIsolation(pat::PfChargedHadronIso) + std::max(0., muon.userIsolation(pat::PfNeutralHadronIso) + muon.userIsolation(pat::PfGammaIso) - energy))/muon.pt());
      deltaR = 0.4;
      energy = TMath::Pi()*deltaR*deltaR*rho;
      nBranches_->lep_pfRhoCorrRelIso04.push_back((muon.chargedHadronIso() + std::max(0., muon.neutralHadronIso() + muon.photonIso() - energy))/muon.pt());
      nBranches_->lep_pfRhoCorrRelIso04Boost.push_back((muon.userIsolation(pat::PfChargedHadronIso) + std::max(0., muon.userIsolation(pat::PfNeutralHadronIso) + muon.userIsolation(pat::PfGammaIso) - energy))/muon.pt());
      /*===== ISO ====*/
      
      nBranches_->lep_pfDeltaCorrRelIso     .push_back((muon.chargedHadronIso() + std::max(0., muon.neutralHadronIso() + muon.photonIso() - 0.5*muon.puChargedHadronIso()))/muon.pt());
      nBranches_->lep_pfRelIso              .push_back((muon.chargedHadronIso() + muon.neutralHadronIso()+ muon.photonIso())/muon.pt()) ; 
      nBranches_->lep_photonIso             .push_back(muon.photonIso());
      nBranches_->lep_neutralHadIso         .push_back(muon.neutralHadronIso());
      nBranches_->lep_chargedHadIso         .push_back(muon.chargedHadronIso());
      nBranches_->lep_trackIso              .push_back(muon.trackIso());
   

      nBranches_->lep_pfDeltaCorrRelIsoBoost.push_back((muon.userIsolation(pat::PfChargedHadronIso) + std::max(0., muon.userIsolation(pat::PfNeutralHadronIso) + muon.userIsolation(pat::PfGammaIso) - 0.5*muon.userIsolation(pat::PfPUChargedHadronIso)))/muon.pt());
      nBranches_->lep_pfRelIsoBoost         .push_back((muon.userIsolation(pat::PfChargedHadronIso) + muon.userIsolation(pat::PfNeutralHadronIso)+ muon.userIsolation(pat::PfGammaIso))/muon.pt()) ; 
      nBranches_->lep_photonIsoBoost        .push_back(muon.userIsolation(pat::PfGammaIso));
      nBranches_->lep_neutralHadIsoBoost    .push_back(muon.userIsolation(pat::PfNeutralHadronIso));
      nBranches_->lep_chargedHadIsoBoost    .push_back(muon.userIsolation(pat::PfChargedHadronIso));
     
      //////////// DilMuonPFIso



      /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 Apr9: with Aniello we decide to use the collection "patMuonsWithTrigger"  ,
	 that differs from the "patMuonsWithTriggerBoosted" just for the DIL isolations.
	 I am changing those and save the new one in the ntuple. 
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
	 nBranches_->lep_DilPFIso              .push_back(MuonPFIso(muon,false));

      */ 
     
 
      nBranches_->lep_DilPFIso              .push_back(MuonPFIso(muonNotBoosted,false));
      
      ////DilMuonDETIso Filled later
      nBranches_->lep_SemileptonicPFIso            .push_back(MuonPFIso(muon,true));


      nBranches_->lep_SemileptonicCorrPFIso             .push_back(MuonCorrPFIso(muon,true));
      
      // /////// cout Aniello's & ours Isolation
      //       std::cout << " ********************* Muon **************************" <<std::endl;
      //       std::cout << " pat::PfGammaIso "<<pat::PfGammaIso<<" pat::PfNeutralHadronIso "<< pat::PfNeutralHadronIso<<" pat::PfChargedHadronIso "<< pat::PfChargedHadronIso <<" pat::User1Iso "<<pat::User1Iso<<" pat::User2Iso " <<pat::User2Iso<< std::endl;
      
      //       std::cout << " lep_pfRhoCorrRelIso03 "<< (muon.chargedHadronIso() + std::max(0., muon.neutralHadronIso() + muon.photonIso() - TMath::Pi()*0.3*0.3*rho))/muon.pt() <<std::endl;
      //       std::cout << " lep_pfRhoCorrRelIso03Boost "<<(muon.userIsolation(pat::PfChargedHadronIso) + std::max(0., muon.userIsolation(pat::PfNeutralHadronIso) + muon.userIsolation(pat::PfGammaIso) - TMath::Pi()*0.3*0.3*rho))/muon.pt()<< std::endl;
      //       std::cout << " lep_pfRhoCorrRelIso04 "<< (muon.chargedHadronIso() + std::max(0., muon.neutralHadronIso() + muon.photonIso() - TMath::Pi()*0.4*0.4*rho))/muon.pt() <<std::endl;
      //       std::cout << " lep_pfRhoCorrRelIso04Boost "<<(muon.userIsolation(pat::PfChargedHadronIso) + std::max(0., muon.userIsolation(pat::PfNeutralHadronIso) + muon.userIsolation(pat::PfGammaIso) - TMath::Pi()*0.4*0.4*rho))/muon.pt()<< std::endl;
      
      //       std::cout << " muon.chargedHadronIso() " << muon.chargedHadronIso() <<" muon.pfIsolationR04().sumChargedHadronPt "<< muon.pfIsolationR04().sumChargedHadronPt<<" muon.userIsolation(pat::PfChargedHadronIso) " <<muon.userIsolation(pat::PfChargedHadronIso)<<std::endl;
      //       std::cout << " muon.neutralHadronIso() " <<  muon.neutralHadronIso()<<"  muon.pfIsolationR04().sumNeutralHadronEt "<< muon.pfIsolationR04().sumNeutralHadronEt<< " muon.userIsolation(pat::PfNeutralHadronIso) "<<  muon.userIsolation(pat::PfNeutralHadronIso)<<std::endl;
      //       std::cout << " muon.photonIso() "<<   muon.photonIso()<< " muon.pfIsolationR04().sumPhotonEt "<< muon.pfIsolationR04().sumPhotonEt <<" muon.userIsolation(pat::PfGammaIso) " <<muon.userIsolation(pat::PfGammaIso)<<std::endl;
      
      //       std::cout << " muon.pfIsolationR04().sumPUPt "<<  muon.pfIsolationR04().sumPUPt<< " muon.userIsolation(pat::User2Iso) "<<muon.userIsolation(pat::User2Iso) <<std::endl;
      //       std::cout << " muon.userIsolation( pat::User1Iso) "<<  muon.userIsolation( pat::User1Iso)<<std::endl;
      //       std::cout << " MuonCorrPFIso(muon,true) " << MuonCorrPFIso(muon,true) <<" MuonPFIso(muon,true) "<< MuonPFIso(muon,true)<<std::endl;
      //       std::cout << " MuonCorrPFIso(muon,false) " << MuonCorrPFIso(muon,false) <<" MuonPFIso(muon,false) "<< MuonPFIso(muon,false)<<std::endl;
      //       std::cout << " ********************* Muon End **************************" <<std::endl;


      double normChi2        = -99;
      int    trackerHits     = -99;
      int    pixelHits       = -99;
      int    globalMuonHits  = -99;
      
      if( muon.isGlobalMuon() ) 
	normChi2=muon.normChi2();
              
      if( !muon.track().isNull() )
        trackerHits = (muon.track())->hitPattern().trackerLayersWithMeasurement();
      
      if( !muon.innerTrack().isNull() )
        pixelHits = (muon.innerTrack())->hitPattern().numberOfValidPixelHits();
      
      if( !muon.globalTrack().isNull() )
        globalMuonHits = (muon.globalTrack())->hitPattern().numberOfValidMuonHits();
                  
      nBranches_->lep_normChi2       .push_back(normChi2);
      nBranches_->lep_trackerHits    .push_back(trackerHits);
      nBranches_->lep_matchedStations.push_back(muon.numberOfMatchedStations());
      nBranches_->lep_pixelHits      .push_back(pixelHits);
      nBranches_->lep_globalHits     .push_back(globalMuonHits);

      nBranches_->lep_idMVAtrig  .push_back(-99);
      nBranches_->lep_looseId	 .push_back(-99);
      nBranches_->lep_vetoId	 .push_back(-99);
      nBranches_->lep_tightId	 .push_back(-99);
      nBranches_->lep_trigTightId.push_back(-99);
      nBranches_->lep_isHEEP     .push_back(-99);
      
      /////////Semileptonic Selection
      bool isASelection=false;       
      bool isAGlobal= false ;
     
      SelectTrackerGlobalID(muon,vertices_->at(0),isAGlobal);
     
      nBranches_->lep_isASelection    .push_back(isAGlobal  );
      nBranches_->lep_isAGlobal      .push_back(isAGlobal   );

  //     std::cout << "isASelection "<<isASelection<< std::endl;
      if (isAGlobal==1){  

	reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 17., 40., 0.25)).first;

	nBranches_->lep_cktTrack_pt.push_back (cktTrack->pt());}
      else 	nBranches_->lep_cktTrack_pt.push_back (-99);

      std::vector<pat::MuonCollection::const_iterator> SelectedHighptMuo;
      for(pat::MuonCollection::const_iterator muon_tmp = muons_->begin(); muon_tmp != muons_->end(); ++muon_tmp) {
	if(muon_tmp->pt()<10) continue;
	if(fabs(muon_tmp->eta())>2.4) continue;
	if(fabs(muon_tmp->phi())>3.2) continue;
	if(!(muon_tmp->isTrackerMuon())) continue;
	if(muon_tmp->numberOfMatches()<=1) continue;
	if(fabs(muon_tmp->muonBestTrack()->dz(vertices_->at(0).position()))>=0.5) continue;
	if(fabs(muon_tmp->dB())>=0.2 ) continue;
	if(muon_tmp->innerTrack()->hitPattern().numberOfValidPixelHits()<=0) continue;
	if(muon_tmp->innerTrack()->hitPattern().trackerLayersWithMeasurement()<=5) continue;
	if((muon_tmp->muonBestTrack()->ptError()/muon_tmp->muonBestTrack()->pt())>0.3) continue;
	SelectedHighptMuo.push_back(muon_tmp);
      }
      
      /////////Dilepton Selection

      bool isADilSelection=false;
      /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 Apr9: with Aniello we decide to use the collection "patMuonsWithTrigger"  ,
	 that differs from the "patMuonsWithTriggerBoosted" just for the DIL isolations.
	 I am changing those and save the new one in the ntuple. 
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 
	 float MuonDETIso = muon.trackIso()
	 
      */ 
     
      float MuonDETIso = muonNotBoosted.trackIso()// /muon.pt()
	;
      for(unsigned int j = 0; j< SelectedHighptMuo.size();++j){
	if((muon.pt()==SelectedHighptMuo[j]->pt()) && (muon.eta()==SelectedHighptMuo[j]->eta()) && (muon.phi()==SelectedHighptMuo[j]->phi())) continue;
	double dR = ROOT::Math::VectorUtil::DeltaR(muon.p4(),SelectedHighptMuo[j]->p4());
	if(dR < 0.3) MuonDETIso = MuonDETIso - ((SelectedHighptMuo[j]->track()->pt())// /muon.pt()
						);
      }
//       std ::cout    <<" MuonDETIso " <<  MuonDETIso <<std::endl;

      nBranches_->lep_DilDETIso    .push_back(MuonDETIso);
      
    
      bool hasAtLeastTrackerMuon;
      SelectTrackerMuon(muon , vertices_->at(0), hasAtLeastTrackerMuon);
      
      isADilSelection=hasAtLeastTrackerMuon;
   //    nBranches_->lep_isADilSelection.push_back(isADilSelection);

      nBranches_->lep_isATracker.push_back(hasAtLeastTrackerMuon);
      
 

      //Cesar's selection 
      bool isCSelection=false;



      isCSelection=isASelection ;
      nBranches_->lep_isCSelection.push_back(isCSelection);

      nBranches_->decayModeFindingNewDMs                     .push_back(-1);
   nBranches_->decayModeFindingOldDMs	             .push_back(-1);
   nBranches_->decayModeFinding	                     .push_back(-1);
   nBranches_->byLooseIsolation                           .push_back(-1);
   nBranches_->byVLooseCombinedIsolationDeltaBetaCorr     .push_back(-1);
   nBranches_->byLooseCombinedIsolationDeltaBetaCorr      .push_back(-1);
   nBranches_->byMediumCombinedIsolationDeltaBetaCorr     .push_back(-1);
   nBranches_->byTightCombinedIsolationDeltaBetaCorr      .push_back(-1);
   nBranches_->byCombinedIsolationDeltaBetaCorrRaw        .push_back(-1);
   nBranches_->byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(-1);
   nBranches_->byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(-1);
   nBranches_->byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(-1);
   nBranches_->byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(-1);
   nBranches_->chargedIsoPtSum                            .push_back(-1);
   nBranches_->neutralIsoPtSum                            .push_back(-1);
   nBranches_->puCorrPtSum                                .push_back(-1);
   nBranches_->byIsolationMVA3oldDMwoLTraw                .push_back(-1); 
   nBranches_->byVLooseIsolationMVA3oldDMwoLT             .push_back(-1);
   nBranches_->byLooseIsolationMVA3oldDMwoLT              .push_back(-1);
   nBranches_->byMediumIsolationMVA3oldDMwoLT             .push_back(-1);
   nBranches_->byTightIsolationMVA3oldDMwoLT              .push_back(-1);
   nBranches_->byVTightIsolationMVA3oldDMwoLT             .push_back(-1);
   nBranches_->byVVTightIsolationMVA3oldDMwoLT            .push_back(-1);
   nBranches_->byIsolationMVA3oldDMwLTraw                 .push_back(-1);
   nBranches_->byVLooseIsolationMVA3oldDMwLT              .push_back(-1);
   nBranches_->byLooseIsolationMVA3oldDMwLT               .push_back(-1);
   nBranches_->byMediumIsolationMVA3oldDMwLT              .push_back(-1);
   nBranches_->byTightIsolationMVA3oldDMwLT               .push_back(-1);
   nBranches_->byVTightIsolationMVA3oldDMwLT              .push_back(-1);
   nBranches_->byVVTightIsolationMVA3oldDMwLT             .push_back(-1);
   nBranches_->byIsolationMVA3newDMwoLTraw                .push_back(-1);
   nBranches_->byVLooseIsolationMVA3newDMwoLT             .push_back(-1);
   nBranches_->byLooseIsolationMVA3newDMwoLT              .push_back(-1);
   nBranches_->byMediumIsolationMVA3newDMwoLT             .push_back(-1);
   nBranches_->byTightIsolationMVA3newDMwoLT              .push_back(-1);
   nBranches_->byVTightIsolationMVA3newDMwoLT             .push_back(-1);
   nBranches_->byVVTightIsolationMVA3newDMwoLT            .push_back(-1);
   nBranches_->byIsolationMVA3newDMwLTraw                 .push_back(-1);
   nBranches_->byVLooseIsolationMVA3newDMwLT              .push_back(-1);
   nBranches_->byLooseIsolationMVA3newDMwLT               .push_back(-1);
   nBranches_->byMediumIsolationMVA3newDMwLT              .push_back(-1);
   nBranches_->byTightIsolationMVA3newDMwLT               .push_back(-1);
   nBranches_->byVTightIsolationMVA3newDMwLT              .push_back(-1);
   nBranches_->byVVTightIsolationMVA3newDMwLT             .push_back(-1);
   nBranches_->againstElectronLoose                       .push_back(-1);
   nBranches_->againstElectronMedium                      .push_back(-1);
   nBranches_->againstElectronTight                       .push_back(-1);
   nBranches_->againstElectronMVA5raw                     .push_back(-1);
   nBranches_->againstElectronMVA5category                .push_back(-1);
   nBranches_->againstElectronVLooseMVA5                  .push_back(-1);
   nBranches_->againstElectronLooseMVA5                   .push_back(-1);
   nBranches_->againstElectronMediumMVA5                  .push_back(-1);
   nBranches_->againstElectronTightMVA5                   .push_back(-1);
   nBranches_->againstElectronVTightMVA5                  .push_back(-1);
   nBranches_->againstElectronDeadECAL                    .push_back(-1);
   nBranches_->againstMuonLoose                           .push_back(-1);
   nBranches_->againstMuonMedium                          .push_back(-1);
   nBranches_->againstMuonTight                           .push_back(-1);
   nBranches_->againstMuonLoose2                          .push_back(-1);
   nBranches_->againstMuonMedium2                         .push_back(-1);
   nBranches_->againstMuonTight2                          .push_back(-1);
   nBranches_->againstMuonLoose3                          .push_back(-1);
   nBranches_->againstMuonTight3                          .push_back(-1);
   nBranches_->againstMuonMVAraw                          .push_back(-1);
   nBranches_->againstMuonLooseMVA                        .push_back(-1);
   nBranches_->againstMuonMediumMVA                       .push_back(-1);
   nBranches_->againstMuonTightMVA                        .push_back(-1);  

    } 



    
    nBranches_->nlep +=  muons_->size();
    
}

