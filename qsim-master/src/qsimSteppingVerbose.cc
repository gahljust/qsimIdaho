//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: qsimSteppingVerbose.cc,v 1.4 2006-06-29 17:54:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "qsimSteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4UnitsTable.hh" 

#include <TFile.h>
#include <TTree.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimSteppingVerbose::qsimSteppingVerbose()
{
       outfile = new TFile("/home/fisibula/install_aps/QSIM/build/qsim_trackEvent.root","RECREATE");
       tree = new TTree("tree", "qsim traking file");
       tree->Branch("EventID",&EventID,"EventID/I");
       tree->Branch("PDGEncoding",&PDGEncoding,"PDGEncoding/I");
       tree->Branch("ParticleName",&ParticleName);
       tree->Branch("StepNumber",&StepNumber,"StepNumber/I");
       tree->Branch("TrackID",&TrackID,"TrackID/I");
       tree->Branch("ParentID",&ParentID,"ParentID/I");
       tree->Branch("Position_x",&Position_x,"Position_x/D");
       tree->Branch("Position_y",&Position_y,"Position_y/D");
       tree->Branch("Position_z",&Position_z,"Position_z/D");
       tree->Branch("Momentum_x",&Momentum_x,"Momentum_x/D");
       tree->Branch("Momentum_y",&Momentum_y,"Momentum_y/D");
       tree->Branch("Momentum_z",&Momentum_z,"Momentum_z/D");
       tree->Branch("KineticEnergy",&KineticEnergy,"KineticEnergy/D");
       tree->Branch("VolumeName",&VolumeName);
       tree->Branch("ProcessName",&ProcessName);

       //tree->Write();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimSteppingVerbose::~qsimSteppingVerbose()
{
 //outfile->cd();
 //tree->Write();
 outfile->Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void qsimSteppingVerbose::StepInfo()
{
  CopyState();
  
  G4int prec = G4cout.precision(3);

  //if(fTrack->GetVolume()->GetName()=="Quartz" && fTrack->GetCurrentStepNumber()==1 && fTrack->GetKineticEnergy()>0)
    //{
      EventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
      PDGEncoding = fTrack->GetDefinition()->GetPDGEncoding();
      ParticleName = fTrack->GetDefinition()->GetParticleName();
      StepNumber = fTrack->GetCurrentStepNumber();
      TrackID = fTrack->GetTrackID();
      ParentID = fStep->GetTrack()->GetParentID();
      Position_x = fTrack->GetPosition().x();
      Position_y = fTrack->GetPosition().y();
      Position_z = fTrack->GetPosition().z();
      Momentum_x = fTrack->GetMomentum().x();
      Momentum_y = fTrack->GetMomentum().y();
      Momentum_z = fTrack->GetMomentum().z();
      KineticEnergy = fTrack->GetKineticEnergy();
      VolumeName = fTrack->GetVolume()->GetName();
      ProcessName = fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

      tree->Fill();      
    //}

  if( verboseLevel >= 1 ){
    //if( verboseLevel >= 4 ) VerboseTrack();
    /*if( verboseLevel >= 3 ){
      G4cout << G4endl;
      G4cout << std::setw( 5) << "#EventID#"     << " "
             << std::setw( 5) << "#Step#"     << " "
	     << std::setw( 6) << "X"          << "    "
	     << std::setw( 6) << "Y"          << "    "  
	     << std::setw( 6) << "Z"          << "    "
             << std::setw( 6) << "Px"          << "    "
             << std::setw( 6) << "Py"          << "    "
             << std::setw( 6) << "Pz"          << "    "
	     << std::setw( 9) << "KineE"      << " "
	     << std::setw( 9) << "dEStep"     << " "  
	     << std::setw(10) << "StepLeng"     
	     << std::setw(10) << "TrakLeng" 
	     << std::setw(10) << "Volume"    << "  "
	     << std::setw(10) << "Process"   << G4endl;    	          
    }*/

    /*G4cout << std::setw(5) << G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID() << " "
           << std::setw(5) << fTrack->GetCurrentStepNumber() << " "
           << std::setw(6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
           << std::setw(6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
           << std::setw(6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
           << std::setw(6) << G4BestUnit(fTrack->GetMomentum().x(),"Energy")
           << std::setw(6) << G4BestUnit(fTrack->GetMomentum().y(),"Energy")
           << std::setw(6) << G4BestUnit(fTrack->GetMomentum().z(),"Energy")
           << std::setw(6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
           << std::setw(6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
           << std::setw(6) << G4BestUnit(fStep->GetStepLength(),"Length")
           << std::setw(6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
           << "   ";*/

    // if( fStepStatus != fWorldBoundary){ 
    if( fTrack->GetNextVolume() != 0 ) { 
      //G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
    } else {
      //G4cout << std::setw(10) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
      //G4cout << "  "
      //       << std::setw(10)
	//     << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	  //                               ->GetProcessName();
    } else {
      //G4cout << "   UserLimit";
    }

    //G4cout << G4endl;

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	/*G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot 
	       << "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << G4endl;*/

	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot; 
                        lp1<(*fSecondary).size(); lp1++){
	  /*G4cout << "    : "
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << std::setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  G4cout << G4endl;*/
	}
              
	/*G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << G4endl;*/
      }
    }
    
  }
  G4cout.precision(prec);
  outfile->Write("",TObject::kOverwrite);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void qsimSteppingVerbose::TrackingStarted()
{

  CopyState();
  G4int prec = G4cout.precision(3);   
      EventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
      PDGEncoding = fTrack->GetDefinition()->GetPDGEncoding();
      ParticleName = fTrack->GetDefinition()->GetParticleName();
      StepNumber = fTrack->GetCurrentStepNumber();
      TrackID = fTrack->GetTrackID();
      ParentID = fStep->GetTrack()->GetParentID();
      Position_x = fTrack->GetPosition().x();
      Position_y = fTrack->GetPosition().y();
      Position_z = fTrack->GetPosition().z();
      Momentum_x = fTrack->GetMomentum().x();
      Momentum_y = fTrack->GetMomentum().y();
      Momentum_z = fTrack->GetMomentum().z();
      KineticEnergy = fTrack->GetKineticEnergy();
      VolumeName = fTrack->GetVolume()->GetName();
      ProcessName = "initStep";

      tree->Fill();
      
  if( verboseLevel > 0 ){

    /*G4cout << std::setw( 5) << "EventID#"      << " "
           << std::setw( 5) << "Step#"      << " "
           << std::setw( 6) << "X"          << "    "
	   << std::setw( 6) << "Y"          << "    "  
	   << std::setw( 6) << "Z"          << "    "
           << std::setw( 6) << "Px"          << "    "
           << std::setw( 6) << "Py"          << "    "
           << std::setw( 6) << "Pz"          << "    "
	   << std::setw( 9) << "KineE"      << " "
	   << std::setw( 9) << "dEStep"     << " "  
	   << std::setw(10) << "StepLeng"  
	   << std::setw(10) << "TrakLeng"
	   << std::setw(10) << "Volume"     << "  "
	   << std::setw(10) << "Process"    << G4endl;	     

   G4cout << std::setw( 5) << G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID() << " "
          << std::setw( 5) << fTrack->GetCurrentStepNumber() << " "
          << std::setw( 6) << G4BestUnit(fTrack->GetPosition().x(),"Length")
          << std::setw( 6) << G4BestUnit(fTrack->GetPosition().y(),"Length")
          << std::setw( 6) << G4BestUnit(fTrack->GetPosition().z(),"Length")
          << std::setw(6) << G4BestUnit(fTrack->GetMomentum().x(),"Energy")
          << std::setw(6) << G4BestUnit(fTrack->GetMomentum().y(),"Energy")
          << std::setw(6) << G4BestUnit(fTrack->GetMomentum().z(),"Energy")
          << std::setw( 6) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
          << std::setw( 6) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
          << std::setw( 6) << G4BestUnit(fStep->GetStepLength(),"Length")
          << std::setw( 6) << G4BestUnit(fTrack->GetTrackLength(),"Length")
          << "  ";*/

    if(fTrack->GetNextVolume()){
      //G4cout << std::setw(10) << fTrack->GetVolume()->GetName();
    } else {
      //G4cout << "OutOfWorld";
    }
    //G4cout << "    initStep" << G4endl;
  }
  G4cout.precision(prec);
  outfile->Write("",TObject::kOverwrite);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

