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
// $Id: qsimSteppingVerbose.hh,v 1.2 2006-06-29 17:54:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class qsimSteppingVerbose;

#ifndef qsimSteppingVerbose_h
#define qsimSteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

#include "TROOT.h"
#include "TObject.h"

class TFile;
class TTree;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class qsimSteppingVerbose : public G4SteppingVerbose
{
 public:   

   qsimSteppingVerbose();
  ~qsimSteppingVerbose();

   void StepInfo();
   void TrackingStarted();

 private:
   TFile *outfile;
   TTree *tree;

   int EventID;
   int PDGEncoding;
   std::string ParticleName;
   int StepNumber;
   int TrackID;
   int ParentID;
   double Position_x;
   double Position_y;
   double Position_z;
   double Momentum_x;
   double Momentum_y;
   double Momentum_z;
   double KineticEnergy;
   std::string VolumeName;
   std::string ProcessName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

