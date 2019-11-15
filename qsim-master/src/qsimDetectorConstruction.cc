#include "qsimDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "qsimDetector.hh"
#include "qsimScintDetector.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4GenericTrap.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Line 494 Quarts and Tungsten



void qsimDetectorConstruction::DetModeSet(G4int detMode = 3) {

    fDetMode = detMode;
    // 0 is PREX-I design
    // 1 is PREX-II prototype (so-called "design 3")
        // 2 SAM
        // 3 showerMax
        // 4 tandem mount

}

void qsimDetectorConstruction::QModeSet(G4int qMode = 2) {

    fQMode = qMode;
    // 0 is PREX-I design
    // 0 is PREX-II prototype (so-called "design 3")
        // 1 SAM
        // 2 showerMax

}

void qsimDetectorConstruction::StandModeSet(G4int standMode = 0) {

    fStandMode = standMode;
    // 1 cosmic setup (detector, lead, scintillators)
    // 0 beam setup (detector only)

}

void qsimDetectorConstruction::ConfModeSet(G4int confMode = 0) {

    fConfMode = confMode;
    // 0 is 1Q
    // 1  is 1 stack
    // 2 is 2 stack
    // 3  is 3 stack
    // 4 is 4 stack

}

qsimDetectorConstruction::qsimDetectorConstruction() {

    DetModeSet();
    QModeSet();
    StandModeSet();
    ConfModeSet();

    fQuartzPolish = 0.98;
    fDetAngle = 0.0*deg;
    fDetPosX = 0.0*cm;
    fDetPosY = 0.0*cm;

    // fNewStand = false; // Default setting is for the setup to reflect to old cosmic stand. True will go to the new design. Messenger has commands to switch between these at command line or at batch mode run time as well.
    // fAccBeamStand = false; // Only affects stand components: true deletes the lead block.

    // Quartz (defined for all modes) Quartz2 is only for tandem mount, otherwise it is not used.
    if(fQMode == 0){
    quartz_x = 1.75*cm; // CSC measures in SolidWorks 0.689 x 2.952 x 0.197 cm
    quartz2_x = 1.75*cm;
    quartz_y = 7.5*cm;  //2.5
    quartz2_y = 7.5*cm;
    //Half cm
    //quartz_z = 0.3*cm;
    //One cm
    quartz_z = 0.3*cm;//downstream
    quartz2_z = 0.5*cm;//upstream
    }

    if(fQMode == 1){
    quartz_x = 19.5*mm/2;
    quartz2_x = 1.75*cm;
    quartz_y = 1.85*cm/2;//1.65*cm;
    quartz2_y = 7.5*cm;
    //Change quartz thickness here.
    quartz_z = 3*mm;//0.65*cm;
    quartz2_z = 0.5*cm;
    }
// ===================================================================REPLACE THIS STUFF===============================================================================
    if (fQMode == 2) {
        //quartz_x = 246*mm/2; // replaces lines 90-93
	quartz_x = 25.40*mm/2; //benchmarck
	quartz2_x = 1.75*cm;
        quartz_y = 117.5*mm/2; // 117.5 (open) ---> 152.5 (close)
	//Change quartz thickness here.
	quartz2_y = 7.5*cm;
        quartz_z = 10.0*mm/2; // REPLACE END
	quartz2_z = 0.5*cm;
    }
// ====================================================================================================================================================================
    quartz_zPos = 0.0*cm; //-.9*cm;//-1.1*cm; //-.9*cm; //-.6*cm;

    cone_rmin1 = 2.1*cm;
    cone_rmax1 = cone_rmin1+.05*cm;
    cone_rmin2 = 2.5*cm;  // normally 2.5*cm;
    cone_rmax2 = cone_rmin2+.05*cm;
    cone_z = quartz_y+.5*cm;    //3
    cone_sphi = 0.;
    cone_fphi = 2*3.1415;

    rin = cone_rmin2;  // normally 2.5*cm;
    rout = rin+.05*cm;
    lngth = 1.6*cm;  // PMT dist. = 2*lngth +1cm  (10.4 == 4.5, 6.8 == 2.9)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qsimDetectorConstruction::~qsimDetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* qsimDetectorConstruction::Construct() {

    // Define materials

    G4double a, z, density;
    G4int nelements;

    // Air
    G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
    G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

    G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);//density=1.29*mg/cm3
    Air->AddElement(N, 79.*perCent);
    Air->AddElement(O, 21.*perCent);

    // Quartz
    G4Element* Si = new G4Element("Silicon", "Si", z=14 , a=28*g/mole);

    G4Material* Quartz = new G4Material("Quartz", density= 2.203*g/cm3, nelements=2);
    Quartz->AddElement(Si, 1);
    Quartz->AddElement(O, 2);

    // Aluminum for mirror and stand (need separate materials so that mirror can be reflective)
    G4Element* Al = new G4Element("Aluminum", "Al", z=13 , a=27*g/mole);
    G4Material* Alu_Mat = new G4Material("Alu_Mat", 2.7*g/cm3, nelements=1);
    Alu_Mat->AddElement(Al, 1);
    G4Material* Mirror = new G4Material("Mirror", density= 2.7*g/cm3, nelements=1);
    Mirror->AddElement(Al, 1);

    // Lead
    G4Element* Pb = new G4Element("Lead", "Pb", z=82 , a=207.2*g/mole);
    G4Material* Pb_Mat = new G4Material("Pb_Mat", 11.34*g/cm3, nelements=1);
    Pb_Mat->AddElement(Pb, 1);

    // Let us make cathode from a special metal (reflectivity 0, efficiency of photoelectrons 25%)
    G4Material* CATH = new G4Material("CATH", density= 2.7*g/cm3, nelements=1);
    CATH->AddElement(Al, 1);

    //Tungsten
    G4Element* W = new G4Element("Tungsten", "W", z=74 , a=183.84*g/mole);
    G4Material* Tungsten = new G4Material("Tungsten", density= 19.25*g/cm3, nelements=1);
    Tungsten->AddElement(W, 1);

    //Iron
    G4Element* Fe = new G4Element("Iron", "Fe", z=26 , a=55.846*g/mole);

    //Chromium
    G4Element* Cr = new G4Element("Chromium", "Cr", z=24 , a=51.996*g/mole);

    //Nickel
    G4Element* Ni = new G4Element("Nickel", "Ni", z=28 , a=58.693*g/mole);

    //Carbon
    G4Element* C = new G4Element("Carbon", "C", z=6 , a=12.011*g/mole);

    //Steel
    G4Material* Steel = new G4Material("Steel", density= 8.06*g/cm3, nelements=4);
    G4double fractionmass;
    Steel->AddElement(Fe, fractionmass=0.712);
    Steel->AddElement(Cr, fractionmass=0.18);
    Steel->AddElement(Ni, fractionmass=0.09);
    Steel->AddElement(C, fractionmass=0.018);




    // Define optical property tables

    const G4int nEntries = 49;//49 R7723Q -actualized //47 9305QKFL//205 all - old//101 280nm cut off//75 320nm cut off

    // Array of photon energies for R7723Q and R375
    /*G4double PhotonEnergy[nEntries] =
    {   2.4,2.42,2.44,2.46,2.48,2.5,2.52,2.54,2.56,2.58,//515nm
        2.6,2.62,2.64,2.66,2.68,2.7,2.72,2.74,2.76,2.78,
        2.8,2.82,2.84,2.86,2.88,2.9,2.92,2.94,2.96,2.98,
        3.0,3.02,3.04,3.06,3.08,3.1,3.12,3.14,3.16,3.18,
        3.2,3.22,3.24,3.26,3.28,3.3,3.32,3.34,3.36,3.38,//Cut off -> 3.2eV ~ 380nm, 41 entries
        3.4,3.42,3.44,3.46,3.48,3.5,3.52,3.54,3.56,3.58,
        3.6,3.62,3.64,3.66,3.68,3.7,3.72,3.74,3.76,3.78,
        3.8,3.82,3.84,3.86,3.88,3.9,3.92,3.94,3.96,3.98,//Cut off -> 3.88eV ~ 320nm, 75 entries
        4.0,4.02,4.04,4.06,4.08,4.1,4.12,4.14,4.16,4.18,         //Glass cuts off above 4.135eV, 87 entries
        4.2,4.22,4.24,4.26,4.28,4.3,4.32,4.34,4.36,4.38,
        4.4,4.42,4.44,4.46,4.48,4.5,4.52,4.54,4.56,4.58,         //Cut off -> 4.4eV ~ 280nm, 101 entries
        4.6,4.62,4.64,4.66,4.68,4.7,4.72,4.74,4.76,4.78,
        4.8,4.82,4.84,4.86,4.88,4.9,4.92,4.94,4.96,4.98,        // Cut off -> 4.96eV ~ 250nm
        5,5.02,5.04,5.06,5.08,5.1,5.12,5.14,5.16,5.18,    // 5.04eV = 246 nm is the 30% cutoff, 133 entries
        5.2,5.22,5.24,5.26,5.28,5.3,5.32,5.34,5.36,5.38,
        5.4,5.42,5.44,5.46,5.48,5.5,5.52,5.54,5.56,5.58,
        5.6,5.62,5.64,5.66,5.68,5.7,5.72,5.74,5.76,5.78,
        5.8,5.82,5.84,5.86,5.88,5.9,5.92,5.94,5.96,5.98,
        6,6.02,6.04,6.06,6.08,6.1,6.12,6.14,6.16,6.18,
        6.21,6.29,6.38,6.48,6.57,6.67,6.78,6.87,6.98,
        7.08,7.20,7.32,7.44,7.56,7.69//161nm
    };*/

    // Array of photon energies for 9305QKFL (97 entries)
    //G4double PhotonEnergy[nEntries] = {1.88, 1.91, 1.94, 1.97, 2.00, 2.03, 2.07, 2.10, 2.14, 2.18, 2.21, 2.25, 2.30, 2.34, 2.38, 2.43, 2.48, 2.53, 2.58, 2.64, 2.70, 2.76, 2.82, 2.88, 2.95, 3.02, 3.10, 3.18, 3.26, 3.35, 3.44, 3.54, 3.65, 3.76, 3.87, 4.00, 4.13, 4.28, 4.43, 4.59, 4.77, 4.96, 5.17, 5.39, 5.64, 5.90, 6.20};

    //Array of photon energies for R7723Q - actualized (49 entries)
    G4double PhotonEnergy[nEntries] = {1.8233, 1.85051, 1.87855, 1.90745, 1.93725, 1.968, 1.99975, 2.03253, 2.0664, 2.10143, 2.13766, 2.17516, 2.214, 2.25426, 2.296, 2.33932, 2.38431, 2.43106, 2.47968, 2.53029, 2.583, 2.63796, 2.69531, 2.7552, 2.81782, 2.88335, 2.952, 3.024, 3.09961, 3.17908, 3.26274, 3.35092, 3.44401, 3.54241, 3.64659, 3.7571, 3.87451, 3.99949, 4.13281, 4.27532, 4.42801, 4.59201, 4.76862, 4.95937, 5.16601, 5.39062, 5.63565, 5.90401, 6.19921};

    // Cathode quantum efficiency
    // Response obtained from the plot of the quantum efficiency as a function of wavelength and then changed to eV for the Bialkali photocathode (synthetic silica)


    //R7723Q
    /*G4double EfficiencyArrayPercent[nEntries] =
    {   11.0,12.0,12.5,13.1,13.5,14.5,15.2,16.0,16.5,17.0, // percentages here
        17.5,18.0,18.5,19.0,19.2,19.7,20.1,20.7,20.9,21.1,
        21.6,22.0,22.5,22.7,23.0,23.5,23.7,24.0,24.0,24.2,
        24.2,24.5,25.0,25.0,25.3,25.5,25.5,25.5,25.5,25.5,
        25.5,25.5,25.5,25.7,26.1,26.1,26.1,26.1,25.6,25.6, //41 entries 3.2eV
        25.6,25.6,25.6,25.6,25.6,25.6,25.6,25.6,25.6,25.6,
        25.6,25.6,25.6,26.1,26.1,26.1,26.1,26.1,26.1,26.1,
        25.6,25.0,25.0,25.0,25.0,24.5,24.5,24.5,24.5,24.3, //25.0 75 entries 3.88eV
        24.0,24.0,24.0,24.0,24.0,24.0,23.5,23.5,23.5,23.5,
        23.5,23.5,23.5,23.3,23.1,22.8,22.6,22.6,22.6,22.6, // 4.38 eV
        22.6,22.6,22.3,22.1,22.1,22.1,22.0,21.8,21.7,21.3, //100 entries 4.58 eV//22.6 101 entries 4.4eV
        21.2,21.0,20.8,20.8,20.8,20.8,20.8,20.8,20.8,20.8, // 4.78 eV
        20.4,20.4,20.4,20.4,20.4,20.4,20.4,20.4,20.2,20.0, // 4.98 eV
        20.0,20.0,20.0,20.0,20.0,20.0,20.0,19.5,19.5,19.5, // 5.18 eV
        19.5,19.5,19.5,19.5,19.1,19.1,19.1,19.1,19.1,19.1, // 5.38 eV
        19.1,19.1,19.1,19.0,18.8,18.8,18.8,18.8,18.8,18.8, // 5.58 eV
        18.8,18.4,18.4,18.4,18.4,18.4,18.4,18.4,18.4,18.4, // 5.78 eV
        18.4,18.4,18.4,18.4,18.4,18.4,18.4,18.4,18.4,18.4, // 5.98eV
        18.4,18.2,18.0,18.0,18.0,18.0,18.0,18.0,18.0,18.0,
        18.0,17.6,17.6,17.6,17.6,17.2,16.5,16.2,15.9,15.2,
        14.9,14.3,12.1,10.2,9.6};*/                         // 6.18 eV

    //R375 quantum efficiency.
    /*G4double EfficiencyArrayPercent[nEntries] =
    {24.9629, 24.7014, 24.7014, 24.7014, 24.4454,
     24.1893, 24.4454, 24.7014, 24.7014, 24.7014,
      24.4454, 24.1893, 24.1893, 24.1893, 23.9386,
     23.9386, 23.6879, 23.6879, 23.6879, 23.4423,//20
     23.4423, 23.1968, 23.1968, 23.1968, 22.9564,
     22.9564, 22.7159, 22.4805, 22.4805, 22.245,
     22.245,  22.0144, 22.0144, 21.7839, 21.7839,
     21.7839, 21.5581, 21.5581, 21.3323, 21.3323,//40
     21.3323, 21.1112, 21.1112, 20.89,   20.89,
     20.89,   20.6735, 20.6735, 20.457,  20.457,
     20.457,  20.457,  20.457,  20.2449, 20.2449,
      20.0329, 20.0329, 20.0329, 20.0329, 19.6176, //80
     19.6176, 19.6176, 19.6176, 19.6176, 19.4143,
     19.4143, 19.2109, 19.2109, 19.2109, 19.2109,
     19.0118, 19.0118, 19.0118, 18.8127, 18.8127,
      18.8127, 18.8127, 18.8127, 18.4227, 18.4227,//100
      18.4227, 18.4227, 18.4227, 18.2317, 18.2317,
     18.0408, 18.0408, 18.0408, 18.0408, 18.0408,
     18.0408, 17.8538, 17.8538, 17.8538, 17.6668,
     17.6668, 17.6668, 17.4836, 17.4836, 17.4836,//120
     17.3005, 17.3005, 17.3005, 16.9419, 16.9419,
     16.9419, 16.7663, 16.7663, 16.7663, 16.7663,
     16.4187, 16.4187, 16.4187, 16.2467, 16.2467,
     16.2467, 16.2467, 15.9099, 15.9099, 15.9099,
     15.745,  15.745,  15.745,  15.745,  15.2571,
     15.2571, 15.2571, 15.2571, 14.4794, 14.4794,
     14.4794, 14.4794, 13.5975, 13.5975, 13.5975,
     13.5975, 13.0396, 13.0396, 13.0396, 13.0396};*/

        //9305QKFL
        //G4double EfficiencyArrayPercent[nEntries] = {0.0107, 0.0777, 0.172, 0.354, 0.653, 1.10, 1.70, 2.39, 3.23, 4.08, 5.09, 6.10, 7.47, 9.71, 12.7, 15.4, 16.4, 17.4, 18.4, 19.8, 21.4, 23.1, 24.7, 25.9, 27.1, 27.6, 28.3, 29.1, 29.1, 28.8, 30.6, 30.9, 31.2, 30.8, 31.3, 30.4, 30.9, 29.8, 28.7, 27.2, 25.7, 25.4, 25.8, 27.3, 29.0, 32.0, 38.8};

    //R7723Q - actualized
    G4double EfficiencyArrayPercent[nEntries] = {0, 0.1, 0.2, 0.3, 0.6, 0.9, 1.4, 2, 2.7, 3.4, 4.2, 5, 5.8, 6.7, 7.7, 9.3, 12, 14.6, 15.9, 16.5, 17.4, 18.3, 19.7, 20.9, 21.8, 22.5, 23.1, 23.6, 24, 24.2, 24.4, 24.2, 24.1, 24.3, 24.4, 24.3, 24.1, 23.7, 22.9, 21.7, 20.1, 18.3, 17.2, 16.7, 15.8, 15, 15.1, 15.4, 16.2};

    //KCsBs Cathode reflectivity
    G4double CathodeReflectivity[nEntries] = {19.9068, 19.9606, 20.1601, 20.2915, 20.5281, 20.8154, 21.0835, 21.4784, 21.886, 22.4168, 23.043, 23.7001, 24.448, 25.1821, 25.7158, 25.1727, 22.9509, 20.9488, 20.5015, 20.8426, 21.4012, 21.7948, 21.2349, 19.832, 19.0381, 18.6857, 18.4112, 17.8851, 17.1713, 16.3226, 15.532, 14.964, 14.5585, 14.1401, 13.9015, 13.44, 12.526, 10.4672, 7.63344, 5.23522, 4.54473, 4.57344, 4.65438, 4.7869, 5.10269, 5.78094, 7.08912, 9.42469, 13.3784};

    G4double EfficiencyArray[nEntries];

    G4double RefractiveIndex1[nEntries];    // quartz refractive index
    G4double Absorption1[nEntries];         // quartz absorption length
    G4double RefractiveIndex2[nEntries];    // air refractive index
    G4double Reflectivity1[nEntries]; //= {0.7016, 0.7034, 0.7053, 0.7053, 0.7053, 0.7053, 0.7062, 0.7071, 0.7071, 0.7071, 0.7071, 0.7071, 0.709, 0.7107, 0.7107, 0.7126, 0.7126, 0.7127, 0.7127, 0.7127, 0.7127, 0.7127, 0.7127, 0.7044, 0.7035, 0.7053, 0.7017, 0.6971, 0.6916, 0.6861, 0.6815, 0.6751, 0.66955, 0.665, 0.6613, 0.6576, 0.65305, 0.64745, 0.6439, 0.6375, 0.6375, 0.6356, 0.6338, 0.632, 0.6282, 0.6237, 0.6182, 0.6127, 0.60895, 0.60355, 0.59805, 0.59255, 0.58705, 0.57785, 0.57235, 0.55585, 0.53565, 0.5246, 0.4898, 0.45125, 0.4146, 0.3522, 0.3173, 0.2329, 0.2201, 0.1623, 0.1522, 0.14485, 0.1687, 0.17605, 0.1852, 0.19625, 0.1972, 0.1944, 0.1926, 0.1889, 0.18705, 0.1843, 0.1825, 0.1788, 0.177, 0.1751, 0.1724, 0.1733, 0.1751, 0.1779, 0.1807, 0.1825, 0.1843, 0.1871, 0.1889, 0.19075, 0.1935, 0.1944, 0.19625, 0.199, 0.2009, 0.2027, 0.2027, 0.2045, 0.2045, 0.20545, 0.2064, 0.2073, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2082, 0.2064, 0.2045, 0.2027, 0.2018, 0.2009, 0.19995, 0.199, 0.1972, 0.1954, 0.19445, 0.1935, 0.1935, 0.1917, 0.1908, 0.18895, 0.188, 0.1862, 0.1862, 0.1853, 0.1844, 0.1844, 0.18345, 0.1825, 0.1825, 0.1816, 0.1807, 0.1798, 0.1789, 0.1789, 0.17795, 0.177, 0.177, 0.1761, 0.1752, 0.1761, 0.18165, 0.1844, 0.18715, 0.1945, 0.19815, 0.2018, 0.2055, 0.20915, 0.2128, 0.21465, 0.21835, 0.2257, 0.23485, 0.2422, 0.25135, 0.25135, 0.26055, 0.2679, 0.27705, 0.28615, 0.2954, 0.30455, 0.31375, 0.3211, 0.32845, 0.3376, 0.34675, 0.35595, 0.35595, 0.3651, 0.3725, 0.38165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165, 0.4165};       // mirror reflectivity (Miro 27 - 30 deg)
//G4double Reflectivity1[nEntries];// = {0.54295, 0.5359, 0.5417, 0.5481, 0.5565, 0.5642, 0.5707, 0.579, 0.58545, 0.59315, 0.5996, 0.60735, 0.608, 0.6016, 0.5964, 0.5906, 0.58545, 0.5796, 0.57385, 0.5694, 0.5642, 0.56615, 0.57385, 0.5822, 0.5887, 0.59705, 0.6034, 0.6118, 0.61895, 0.62665, 0.62145, 0.615, 0.60865, 0.6028, 0.5964, 0.5899, 0.58415, 0.5777, 0.58415, 0.5912, 0.5983, 0.6047, 0.61245, 0.6189, 0.62535, 0.63175, 0.6395, 0.63435, 0.62795, 0.62145, 0.6137, 0.60605, 0.59955, 0.5931, 0.58795, 0.5931, 0.6002, 0.606, 0.6137, 0.6188, 0.62535, 0.63175, 0.6369, 0.64205, 0.63685, 0.6291, 0.62535, 0.6181, 0.6131, 0.6078, 0.6028, 0.59565, 0.59955, 0.60465, 0.60985, 0.6143, 0.6188, 0.6253, 0.63045, 0.6356, 0.63045, 0.6253, 0.6188, 0.61245, 0.6053, 0.59825, 0.59185, 0.5873, 0.5853, 0.59055, 0.59435, 0.59825, 0.6028, 0.60595, 0.6105, 0.6085, 0.60215, 0.59435, 0.58795, 0.57895, 0.57245, 0.56345, 0.55705, 0.55185, 0.5499, 0.5544, 0.55705, 0.5628, 0.56735, 0.5699, 0.5744, 0.57245, 0.56735, 0.56085, 0.55575, 0.54925, 0.5441, 0.53515, 0.52995, 0.52615, 0.52865, 0.53125, 0.53385, 0.53645, 0.53645, 0.53255, 0.5274, 0.52355, 0.51705, 0.5106, 0.50805, 0.50165, 0.4958, 0.49005, 0.4836, 0.47715, 0.4694, 0.46425, 0.45015, 0.4372, 0.42435, 0.41795, 0.405, 0.39345, 0.38055, 0.37415, 0.36765, 0.35995, 0.35735, 0.35095, 0.34835, 0.34065, 0.3355, 0.33165, 0.32515, 0.32265, 0.31485, 0.30845, 0.30585, 0.29945, 0.29945, 0.30715, 0.31105, 0.31745, 0.32135, 0.32515, 0.33165, 0.3342, 0.34195, 0.34575, 0.3529, 0.35605, 0.35995, 0.36765, 0.3715, 0.37415, 0.38185, 0.38575, 0.38955, 0.39605, 0.39985, 0.4037, 0.41145, 0.4147, 0.41785, 0.42435, 0.4282, 0.43205, 0.4353, 0.44235, 0.44625, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425, 0.46425}; // mirror reflectivity (UVS - 30 deg)
G4double Reflectivity_laterals[nEntries];// = {0.7612, 0.7621, 0.764, 0.764, 0.7658, 0.76675, 0.7677, 0.7694, 0.7694, 0.7713, 0.7713, 0.7732, 0.7732, 0.775, 0.775, 0.7769, 0.7769, 0.7787, 0.7787, 0.775, 0.7732, 0.7694, 0.7668, 0.764, 0.7613, 0.7585, 0.7549, 0.7448, 0.73375, 0.7228, 0.71545, 0.7044, 0.6943, 0.6852, 0.6751, 0.66765, 0.6567, 0.64745, 0.6402, 0.62915, 0.6219, 0.6127, 0.598, 0.5815, 0.565, 0.543, 0.52645, 0.50995, 0.49895, 0.48245, 0.4659, 0.44755, 0.42545, 0.4017, 0.3889, 0.36685, 0.343, 0.33015, 0.30815, 0.31365, 0.321, 0.3329, 0.33935, 0.35035, 0.3595, 0.37055, 0.3779, 0.38515, 0.38515, 0.377, 0.3687, 0.3302, 0.288, 0.2696, 0.26045, 0.2476, 0.23845, 0.2293, 0.2201, 0.21275, 0.20175, 0.1935, 0.18525, 0.1825, 0.1834, 0.1843, 0.1862, 0.188, 0.19075, 0.1926, 0.1935, 0.1953, 0.1981, 0.199, 0.2008, 0.20175, 0.2036, 0.2064, 0.2073, 0.2082, 0.2064, 0.20545, 0.2045, 0.2045, 0.2027, 0.2027, 0.2018, 0.2009, 0.19995, 0.199, 0.199, 0.1972, 0.1972, 0.1954, 0.1954, 0.1935, 0.1935, 0.1926, 0.1917, 0.1908, 0.1899, 0.1899, 0.18895, 0.188, 0.1862, 0.1862, 0.1844, 0.18345, 0.1825, 0.1798, 0.1789, 0.1789, 0.177, 0.1761, 0.1752, 0.1743, 0.17245, 0.1715, 0.1706, 0.1697, 0.1679, 0.16695, 0.166, 0.166, 0.16325, 0.1623, 0.1614, 0.1605, 0.1605, 0.1614, 0.1623, 0.1623, 0.1623, 0.16325, 0.1642, 0.1642, 0.1642, 0.1651, 0.166, 0.166, 0.1688, 0.1743, 0.1853, 0.1908, 0.1963, 0.2018, 0.20915, 0.2165, 0.22385, 0.22935, 0.2367, 0.244, 0.25135, 0.2624, 0.2624, 0.2807, 0.30635, 0.3073, 0.3101, 0.31375, 0.3165, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925, 0.31925}; //reflectivity for lateral mirror's walls, Miro 27 60 deg
    G4double Reflectivity2[nEntries];       // cathode reflectivity
    G4double QWReflectivity[nEntries];

    for (int i = 0; i < nEntries; i++) {

        PhotonEnergy[i] = PhotonEnergy[i]*eV;
    EfficiencyArray[i] = 0.01*1.0438*EfficiencyArrayPercent[i];
        RefractiveIndex1[i]= 1.438 + (.01197*PhotonEnergy[i]/eV) - (.001955*PhotonEnergy[i]*PhotonEnergy[i]/eV/eV) + (.0004793*PhotonEnergy[i]*PhotonEnergy[i]*PhotonEnergy[i]/eV/eV/eV);

        // *** need to update this
        //Quartz
        Absorption1[i] = (exp(4.325)*exp(1.191*PhotonEnergy[i]/eV)*exp(-.213*PhotonEnergy[i]*PhotonEnergy[i]/eV/eV)*exp(-.04086*PhotonEnergy[i]*PhotonEnergy[i]*PhotonEnergy[i]/eV/eV/eV))*m;
        if (Absorption1[i] > 25*m) {
            Absorption1[i] = 25*m;
        }

        // *** need to update this
        //if (PhotonEnergy[i] < 4.135*eV) {
	Reflectivity1[i] = 0.8;//0.9825;
        if (PhotonEnergy[i] < 3.53*eV) {
        //if (PhotonEnergy[i] < 3.93*eV) {
	  //Reflectivity1[i] = 0.9825;// 0.625;
            Reflectivity_laterals[i] = 0.85;//0.775;//0.725;//0.85;// 0.825;
        //} else if (PhotonEnergy[i] >= 4.135*eV && PhotonEnergy[i] < 6.203*eV) {
        //} else if (PhotonEnergy[i] >= 3.1*eV && PhotonEnergy[i] < 4.135*eV) {
        //} else if (PhotonEnergy[i] >= 3.93*eV && PhotonEnergy[i] < 5.93*eV) {
            //Reflectivity1[i] = 0.9825;//0.97;//0.5;  // .7
            //Reflectivity_laterals[i] = 0.4;//0.725;//0.75;//0.8;
        } else {
            //Reflectivity1[i] = 0.9825;//0.4;  // .6
            Reflectivity_laterals[i] = 0.2;//0.2;//0.8;
        }

        RefractiveIndex2[i] = 1.;
        Reflectivity2[i] = 0.125;//0.01*CathodeReflectivity[i]; //cathode reflectivity
        QWReflectivity[i] = 0.0;
    }


    // Quartz material property table

    G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
    myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
    myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
    Quartz->SetMaterialPropertiesTable(myMPT1);

    // Air material properties table
    // *** need to add absorption length?
    G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
    myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
    Air->SetMaterialPropertiesTable(myMPT2);


    G4SDManager* SDman = G4SDManager::GetSDMpointer();


    ////////////////////////////////////////////
    // Defining geometry
    ////////////////////////////////////////////

    // Volumes

    // World and detector box

    double world_x, world_y, world_z;
    double det_x, det_y, det_z;

    world_x = world_y = world_z = 5*m;
    if (fDetMode == 0) {
    det_x = 15*cm;
    det_y = 6*cm;
    det_z = 6*cm;
    }
    if (fDetMode == 1) {
    det_x = 15*cm;
    det_y = 6*cm;
    det_z = 6*cm;
    }
    if (fDetMode == 2) {
    det_x = 90*cm;
    det_y = 6*cm;
    det_z = 6*cm;
    }

    if (fDetMode == 3) {
    det_x = 20*cm;
    det_y = 20*cm;
    det_z = 10*cm;
    }

    if (fDetMode == 4) {
    det_x = 20*cm;
    det_y = 20*cm;
    det_z = 5*cm;
    }

    G4Box* world_box = new G4Box("World",world_x,world_y,world_z);

    G4LogicalVolume* world_log
    = new G4LogicalVolume(world_box,Air,"World",0,0,0);

    world_log->SetVisAttributes(G4VisAttributes::Invisible);

    G4VPhysicalVolume* world_phys
    = new G4PVPlacement(0,G4ThreeVector(),world_log,"World",0,false,0);

    G4Box* det_box = new G4Box("Detector",det_x,det_y,det_z);
    G4Box* det_box_c1 = new G4Box("Detector",det_x,det_y,det_z);

    G4LogicalVolume* det_log
    = new G4LogicalVolume(det_box,Air,"Detector_log",0,0,0);

    G4LogicalVolume* det_log_c1
    = new G4LogicalVolume(det_box_c1,Air,"Detector_log_c1",0,0,0);

    G4LogicalVolume* det_log_c2
    = new G4LogicalVolume(det_box_c1,Air,"Detector_log_c2",0,0,0);

    det_log->SetVisAttributes(G4VisAttributes::Invisible);
    det_log_c1->SetVisAttributes(G4VisAttributes::Invisible);
    det_log_c2->SetVisAttributes(G4VisAttributes::Invisible);

    //G4double maxStep = 0.1*mm, minEkin = 0.0*eV, minRange = 1*m;
    //G4UserLimits *MyLimits = new G4UserLimits();
    //MyLimits->SetMaxAllowedStep(maxStep);
    //MyLimits->SetUserMinEkine(minEkin);
    //MyLimits->SetUserMinRange(minRange);


    // First, create solids and logical volumes
 // QUARTZ DIMENSIONS ======================================================================================================================================================
    G4double q_yLB = quartz_y - (quartz_z);
    G4double q_yLB2 = quartz2_y - (quartz2_z);

    G4Trap* quartz_box = new G4Trap("Quartz", 2*quartz_x, 2*quartz_z, 2*quartz_y, 2*q_yLB);

    G4LogicalVolume* quartz_log
    = new G4LogicalVolume(quartz_box,Quartz,"Quartz",0,0,0);

    //quartz_log->SetUserLimits(MyLimits);

    G4Trap* quartz_box2 = new G4Trap("Quartz2", 2*quartz2_x, 2*quartz2_z, 2*quartz2_y, 2*q_yLB2);

    G4LogicalVolume* quartz_log2
    = new G4LogicalVolume(quartz_box2,Quartz,"Quartz2",0,0,0);

    //qsimScintDetector* quartzSD = new qsimScintDetector("QuartzSD", 10);

    //SDman->AddNewDetector(quartzSD);
    //quartz_log->SetSensitiveDetector(quartzSD);
 //=====================================================================================================================================================================
    G4RotationMatrix* rotQ = new G4RotationMatrix;

    rotQ->rotateX(M_PI/2.*rad);
    rotQ->rotateZ(0.*deg);
    if(fDetMode == 0) {
        rotQ->rotateZ(M_PI*rad);
    }

    G4VPhysicalVolume* quartz_phys
    = new G4PVPlacement(rotQ,G4ThreeVector(0,0,quartz_zPos),quartz_log,"Quartz", det_log,false,0);

    G4RotationMatrix* rotQ2 = new G4RotationMatrix;

    rotQ2->rotateX(M_PI/2.*rad);
    rotQ2->rotateZ(0.*deg);

    G4VPhysicalVolume* quartz2_phys
    = new G4PVPlacement(rotQ2,G4ThreeVector(0,0,quartz_zPos),quartz_log2,"Quartz2",det_log_c1,false,0);


 // TUNGSTEN DIMENSIONS ===================================================================================================================================================
         Wthickness = 8*mm/2-0.01*mm/2; // REPLACES LINE 374 (this is new)
//=================================================================================
    G4Box* tungsten_box = new G4Box("tungsten_box",quartz_y-quartz_z,quartz_x,Wthickness);

    G4LogicalVolume* tungsten_box_log = new G4LogicalVolume(tungsten_box,Tungsten,"tungsten_box_log",0,0,0);

    G4VisAttributes *tungsten_boxx = new G4VisAttributes();
    tungsten_boxx->SetColour(0.2, 0.2, 0.2);
    tungsten_box_log->SetVisAttributes(tungsten_boxx);

    if (fDetMode == 3) {

      if (fConfMode == 0){ // 1 quartz

       }

       if (fConfMode == 1){// 1 stack
         G4RotationMatrix* rotTung = new G4RotationMatrix;

         rotTung->rotateX(0.*rad);

         // QUARTZ AND TUNGSTEN PLACEMENT ===================================================================================================================================================
              G4VPhysicalVolume* tungsten_box_phys_1 // Replaces lines 390-400
              = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-quartz_z-8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);

         G4RotationMatrix* rotQ_c1 = new G4RotationMatrix;
         rotQ_c1->rotateX(M_PI/2.*rad);
         rotQ_c1->rotateZ(0.*deg);
         rotQ_c1->rotateX(180.*deg);

        }

        if (fConfMode == 2){ //2 stack
          G4RotationMatrix* rotTung = new G4RotationMatrix;

          rotTung->rotateX(0.*rad);

          // QUARTZ AND TUNGSTEN PLACEMENT ===================================================================================================================================================
               G4VPhysicalVolume* tungsten_box_phys_1 // Replaces lines 390-400
               = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-quartz_z-8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);
               G4VPhysicalVolume* tungsten_box_phys_2
               = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-3*quartz_z-3*8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);

          G4RotationMatrix* rotQ_c1 = new G4RotationMatrix;
          rotQ_c1->rotateX(M_PI/2.*rad);
          rotQ_c1->rotateZ(0.*deg);
          rotQ_c1->rotateX(180.*deg);

           G4VPhysicalVolume* quartz_phys_c1 // REPLACES lines 408-425 (THIS IS NEW)
           = new G4PVPlacement(rotQ_c1,G4ThreeVector(0,0,quartz_zPos-2*quartz_z-8*mm),quartz_log,"Quartz", det_log,false,0);



         }

         if (fConfMode == 3){ //3 stack
           G4RotationMatrix* rotTung = new G4RotationMatrix;

           rotTung->rotateX(0.*rad);

           // QUARTZ AND TUNGSTEN PLACEMENT ===================================================================================================================================================
                G4VPhysicalVolume* tungsten_box_phys_1 // Replaces lines 390-400
                = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-quartz_z-8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);
                G4VPhysicalVolume* tungsten_box_phys_2
                = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-3*quartz_z-3*8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);
                G4VPhysicalVolume* tungsten_box_phys_3
                = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-5*quartz_z-5*8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);


            G4RotationMatrix* rotQ_c1 = new G4RotationMatrix;
            rotQ_c1->rotateX(M_PI/2.*rad);
            rotQ_c1->rotateZ(0.*deg);
            rotQ_c1->rotateX(180.*deg);

            G4VPhysicalVolume* quartz_phys_c1 // REPLACES lines 408-425 (THIS IS NEW)
            = new G4PVPlacement(rotQ_c1,G4ThreeVector(0,0,quartz_zPos-2*quartz_z-8*mm),quartz_log,"Quartz", det_log,false,0);

            G4RotationMatrix* rotQ_c2 = new G4RotationMatrix;
            rotQ_c2->rotateX(M_PI/2.*rad);
            rotQ_c2->rotateZ(0.*deg);
            //rotQ_c1->rotateX(180.*deg);

            G4VPhysicalVolume* quartz_phys_c2
            = new G4PVPlacement(rotQ_c2,G4ThreeVector(0,0,quartz_zPos-4*quartz_z-2*8*mm),quartz_log,"Quartz", det_log,false,0);




          }

          if (fConfMode == 4){ //4 stack
            G4RotationMatrix* rotTung = new G4RotationMatrix;

            rotTung->rotateX(0.*rad);

            // QUARTZ AND TUNGSTEN PLACEMENT ===================================================================================================================================================
                 G4VPhysicalVolume* tungsten_box_phys_1 // Replaces lines 390-400
                 = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-quartz_z-8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);
                 G4VPhysicalVolume* tungsten_box_phys_2
                 = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-3*quartz_z-3*8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);
                 G4VPhysicalVolume* tungsten_box_phys_3
                 = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-5*quartz_z-5*8*mm/2),tungsten_box_log,"tungsten", det_log,false,0);
                 G4VPhysicalVolume* tungsten_box_phys_4
                 = new G4PVPlacement(rotTung,G4ThreeVector(-12.5*mm/4,0,quartz_zPos-7*quartz_z-7*8*mm/2-0.5*mm),tungsten_box_log,"tungsten", det_log,false,0); // REPLACE END


             G4RotationMatrix* rotQ_c1 = new G4RotationMatrix;
             rotQ_c1->rotateX(M_PI/2.*rad);
             rotQ_c1->rotateZ(0.*deg);
             rotQ_c1->rotateX(180.*deg);

             G4VPhysicalVolume* quartz_phys_c1 // REPLACES lines 408-425 (THIS IS NEW)
             = new G4PVPlacement(rotQ_c1,G4ThreeVector(0,0,quartz_zPos-2*quartz_z-8*mm),quartz_log,"Quartz", det_log,false,0);

             G4RotationMatrix* rotQ_c2 = new G4RotationMatrix;
             rotQ_c2->rotateX(M_PI/2.*rad);
             rotQ_c2->rotateZ(0.*deg);
             //rotQ_c1->rotateX(180.*deg);

             G4VPhysicalVolume* quartz_phys_c2
             = new G4PVPlacement(rotQ_c2,G4ThreeVector(0,0,quartz_zPos-4*quartz_z-2*8*mm),quartz_log,"Quartz", det_log,false,0);


             G4RotationMatrix* rotQ_c3 = new G4RotationMatrix;
             rotQ_c3->rotateX(M_PI/2.*rad);
             rotQ_c3->rotateZ(0.*deg);
             rotQ_c3->rotateX(180.*deg);

             G4VPhysicalVolume* quartz_phys_c3
             = new G4PVPlacement(rotQ_c3,G4ThreeVector(0,0,quartz_zPos-6*quartz_z-3*8*mm),quartz_log,"Quartz", det_log,false,0);

           }



     //===============================






 // REPLACE END
//======================================================================================================================================================================

    }

    // Light guide and tube mirror (only for PREX-I design)

    G4Trap *lightguide_big = new G4Trap("lightguide_big",11.21*cm,2.97456*deg,90.0*deg,3.15*cm+0.05*cm,3.25*cm+0.05*cm,3.25*cm+0.05*cm,0.0*deg,0.225*cm+0.05*cm,2.0*cm+0.05*cm,2.0*cm+0.05*cm,0.0*deg);
    G4Trap *lightguide_small = new G4Trap("lightguide_small",11.22*cm,2.97456*deg,90.0*deg,3.15*cm,3.25*cm,3.25*cm,0.0*deg,0.225*cm,2.0*cm,2.0*cm,0.0*deg);

    G4VSolid *lightguide_virt = new G4SubtractionSolid("lightguide_virt", lightguide_big, lightguide_small);

    G4LogicalVolume *lightguide_log = new G4LogicalVolume(lightguide_virt, Mirror, "LG_log",0,0,0);

    G4Cons* mirror_tube = new G4Cons("TMirror",cone_rmin2,cone_rmax2,2.5*cm,2.55*cm,lngth,cone_sphi,cone_fphi);

    G4LogicalVolume* tmirror_log = new G4LogicalVolume(mirror_tube,Mirror,"TMirror",0,0,0);

    //SAM Light guide

    //Front
    G4Box* frontPlate_1 = new G4Box("frontPlate_1", 0.5*mm/2, 21.06*mm/2, 7.06*mm/2);

    G4LogicalVolume* frontPlate_1_log = new G4LogicalVolume(frontPlate_1, Mirror, "FrontPlate_1_log",0,0,0);

    G4Box* frontPlate_2 = new G4Box("frontPlate_2", 0.5*mm/2, 20.0*mm/2+0.1*mm, 7.0*mm/2+0.1*mm);

    G4LogicalVolume* frontPlate_2_log = new G4LogicalVolume(frontPlate_2, Mirror, "FrontPlate_2_log",0,0,0);

    //Top
   G4Box* topPlate_1 = new G4Box("topPlate_1", 12.5*mm/2, 21.06*mm/2, 0.5*mm/2);

   G4LogicalVolume* topPlate_1_log = new G4LogicalVolume(topPlate_1, Mirror, "TopPlate_1_log",0,0,0);

   G4Box* topPlate_2 = new G4Box("topPlate_2", 45.27*mm/2 , 21.06*mm/2, 0.5*mm/2);

   G4LogicalVolume* topPlate_2_log = new G4LogicalVolume(topPlate_2, Mirror, "TopPlate_2_log",0,0,0);
//---------------------------------------------------
//Coverting topPlate_2 into a detector

    //qsimDetector* quartzSD_topPlate_2 = new qsimDetector("QuartzSD_topPlate_2", 4);

    //SDman->AddNewDetector(quartzSD_topPlate_2);
    //topPlate_2_log->SetSensitiveDetector(quartzSD_topPlate_2);
//---------------------------------------------------


   G4Box* topPlate_3 = new G4Box("topPlate_3", 300.78*mm/2, 21.06*mm/2, 0.5*mm/2);

   G4LogicalVolume* topPlate_3_log = new G4LogicalVolume(topPlate_3, Mirror, "TopPlate_3_log",0,0,0);

   //Bottom
   //G4Box* botPlate_1 = new G4Box("botPlate_1", 384.05*mm/2+0.05*cm, 20.0*mm/2+0.05*cm+0.1*mm, 0.05*cm/2);
   G4Box* botPlate_1 = new G4Box("botPlate_1", 19.56*mm/2, 21.06*mm/2, 0.5*mm/2);
   G4LogicalVolume* botPlate_1_log = new G4LogicalVolume(botPlate_1, Mirror, "botPlate_1_log",0,0,0);

    G4Box* botPlate_2 = new G4Box("botPlate_2", 355.10*mm/2-19.56*mm/2, 21.06*mm/2, 0.5*mm/2);
    G4LogicalVolume* botPlate_2_log = new G4LogicalVolume(botPlate_2, Mirror, "botPlate_2_log",0,0,0);

   //right

    G4Trap* RPlate_1 = new G4Trap("RPlate_1", 0.5*mm, 7.06*mm, 19.56*mm, 12.5*mm);

    G4LogicalVolume* RPlate_1_log = new G4LogicalVolume(RPlate_1, Mirror, "RPlate_1_log",0,0,0);

    G4Trap* RPlate_2 = new G4Trap("RPlate_2", 0.5*mm, 24.38*mm,  355.10*mm-19.56*mm, 300.78*mm);

    G4LogicalVolume* RPlate_2_log = new G4LogicalVolume(RPlate_2, Mirror, "RPlate_2_log",0,0,0);

    G4int nCVtx = 8;
    std::vector<G4TwoVector> cvtx(nCVtx);
    cvtx[0] = G4TwoVector(0.0*cm, 0.0*mm);
    cvtx[1] = G4TwoVector(-7.06*mm, 7.06*mm);
    cvtx[2] = G4TwoVector(34.76*mm, 24.38*mm);
    cvtx[3] = G4TwoVector(34.76*mm, 24.38*mm);
    cvtx[4] = G4TwoVector(0.0*cm, 0.0*mm);
    cvtx[5] = G4TwoVector(-7.06*mm,  7.06*mm);
    cvtx[6] = G4TwoVector(34.76*mm, 24.38*mm);
    cvtx[7] = G4TwoVector(34.76*mm, 24.38*mm);

    G4GenericTrap* RPlate_3 = new G4GenericTrap("RPlate_3",0.5*mm/2,cvtx);

    G4LogicalVolume* RPlate_3_log = new G4LogicalVolume(RPlate_3, Mirror, "RPlate_3_log",0,0,0);

    G4int nCVtx_2 = 8;
    std::vector<G4TwoVector> cvtx_2(nCVtx_2);
    cvtx_2[0] = G4TwoVector(  0.0*cm,   0.0*mm);
    cvtx_2[1] = G4TwoVector(0.0*mm-0.1*mm,  7*mm+0.1*mm);
    cvtx_2[2] = G4TwoVector( 7*mm+0.1*mm,   0*mm+0.1*mm);
    cvtx_2[3] = G4TwoVector( 7*mm+0.1*mm,   0*mm+0.1*mm);
    cvtx_2[4] = G4TwoVector(  0.0*cm,   0.0*mm);
    cvtx_2[5] = G4TwoVector(0.0*mm-0.1*mm,  7.0*mm+0.1*mm);
    cvtx_2[6] = G4TwoVector( 7*mm+0.1*mm,   0*mm+0.1*mm);
    cvtx_2[7] = G4TwoVector( 7*mm+0.1*mm,   0*mm+0.1*mm);

    G4GenericTrap* RPlate_4 = new G4GenericTrap("RPlate_4",0.5*mm,cvtx_2);

    G4LogicalVolume* RPlate_4_log = new G4LogicalVolume(RPlate_4, Mirror, "RPlate_4_log",0,0,0);


    //left

    G4LogicalVolume* LPlate_1_log = new G4LogicalVolume(RPlate_1, Mirror, "LPlate_log",0,0,0);

    G4LogicalVolume* LPlate_2_log = new G4LogicalVolume(RPlate_2, Mirror, "LPlate_2_log",0,0,0);

    G4GenericTrap* LPlate_3 = new G4GenericTrap("LPlate_3",0.5*mm/2,cvtx);

    G4LogicalVolume* LPlate_3_log = new G4LogicalVolume(LPlate_3, Mirror, "LPlate_3_log",0,0,0);

    G4GenericTrap* LPlate_4 = new G4GenericTrap("LPlate_4",0.5*mm,cvtx_2);

    G4LogicalVolume* LPlate_4_log = new G4LogicalVolume(LPlate_4, Mirror, "LPlate_4_log",0,0,0);


    //Tungsten

    G4Box* W_box = new G4Box("W_blox",20.0*mm/2, 20.0*mm/2, 10.0*mm/2);

    G4LogicalVolume* W_log = new G4LogicalVolume(W_box,Tungsten,"Tungsten",0,0,0);

    G4VisAttributes *W_boxx = new G4VisAttributes();
    W_boxx->SetColour(0.2, 0.2, 0.2);
    W_log->SetVisAttributes(W_boxx);

    //G4double maxStep = 0.00001*mm, minEkin = 1*eV;
    //G4UserLimits *MyLimits = new G4UserLimits();
    //MyLimits->SetMaxAllowedStep(maxStep);
    //MyLimits->SetUserMinEkine(minEkin);
    //W_log->SetUserLimits(MyLimits);

    //G4double maxStep = 100.0*mm, maxLength = 0.1*mm, maxTime = 0.1*ns, minEkin = 10*eV;
    //W_log->SetUserLimits(new G4UserLimits(maxStep, maxLength, maxTime, minEkin));

    //
    //
    //
    //
    //
    //shower-max mirrors
    //
    //
    //
    //
    //


//
//
// REPLACE FROM HERE
//
//

        G4Trap* mirror_box_1 = new G4Trap("mirror_box_1",130.85*mm/2,74*mm/2,0.5*mm/2,0.5*mm/2,102.42*mm/2); // Replaces lines 540-570

        G4LogicalVolume* mirror_box_1_log = new G4LogicalVolume(mirror_box_1,Mirror,"mirror_box_1_log",0,0,0);

        G4VisAttributes *mirror_boxx_1 = new G4VisAttributes();
        mirror_boxx_1->SetForceWireframe(true);
        mirror_box_1_log->SetVisAttributes(mirror_boxx_1);

        G4Trap* mirror_box_2 = new G4Trap("mirror_box_2",70*mm/2,130.85*mm/2,0.5*mm/2,0.5*mm/2,172.34*mm/2);

        G4LogicalVolume* mirror_box_2_log = new G4LogicalVolume(mirror_box_2,Mirror,"mirror_box_1_1_log",0,0,0);

        G4VisAttributes *mirror_boxx_2 = new G4VisAttributes();
        mirror_boxx_2->SetForceWireframe(true);
        mirror_box_2_log->SetVisAttributes(mirror_boxx_2);

        G4Trap* mirror_box_4 = new G4Trap("mirror_box_4",246*mm/2,246*mm/2,0.5*mm/2,0.5*mm/2,106.29*mm/2);

         G4LogicalVolume* mirror_box_4_log = new G4LogicalVolume(mirror_box_4,Mirror,"mirror_box_4_log",0,0,0);

        G4VisAttributes *mirror_boxx_4 = new G4VisAttributes();
        mirror_boxx_4->SetForceWireframe(true);
        mirror_box_4_log->SetVisAttributes(mirror_boxx_4);

        G4Trap* mirror_box_6 = new G4Trap("mirror_box_6",68*mm/2,246*mm/2,0.5*mm/2,0.5*mm/2,150.68*mm/2);

        G4LogicalVolume* mirror_box_6_log = new G4LogicalVolume(mirror_box_6,Mirror,"mirror_box_6_log",0,0,0);

        G4VisAttributes *mirror_boxx_6 = new G4VisAttributes();
        mirror_boxx_6->SetForceWireframe(true);
        mirror_box_6_log->SetVisAttributes(mirror_boxx_6); // REPLACE END

//
//
// TO HERE
//
//

    G4Trap* mirror_box_3 = new G4Trap("mirror_box_3",94.42*mm/2,70.15/2*mm,0.5*mm/2,0.5*mm/2,59.24*mm/2);

    G4LogicalVolume* mirror_box_3_log
    = new G4LogicalVolume(mirror_box_3,Mirror,"mirror_box_3_log",0,0,0);

    G4VisAttributes *mirror_boxx_3 = new G4VisAttributes();
    mirror_boxx_3->SetForceWireframe(true);
    mirror_box_3_log->SetVisAttributes(mirror_boxx_3);

    G4Trap* mirror_box_5 = new G4Trap("mirror_box_5",131.31*mm/2,168.84*mm/2,0.5*mm/2,0.5*mm/2,67.48*mm/2);

    G4LogicalVolume* mirror_box_5_log
    = new G4LogicalVolume(mirror_box_5,Mirror,"mirror_box_4_log",0,0,0);

    G4VisAttributes *mirror_boxx_5 = new G4VisAttributes();
    mirror_boxx_5->SetForceWireframe(true);
    mirror_box_5_log->SetVisAttributes(mirror_boxx_5);
// SUITCASE MIRROR DIMENSIONS =====================================================================================================================================================
        G4Box* mirror_box_7 = new G4Box("mirror_box_5",117.5*mm/2,0.5*mm/2,74*mm/2+1.0*mm); // Replaces lines 596-614

        G4LogicalVolume* mirror_box_7_log
        = new G4LogicalVolume(mirror_box_7,Mirror,"mirror_box_7_log",0,0,0);

        G4VisAttributes *mirror_boxx_7 = new G4VisAttributes();
        mirror_boxx_7->SetForceWireframe(true);
        mirror_box_7_log->SetVisAttributes(mirror_boxx_7);

        G4Box* mirror_box_8 = new G4Box("mirror_box_8",117.5*mm/2,23*mm/2,0.5*mm/2);

        G4LogicalVolume* mirror_box_8_log
        = new G4LogicalVolume(mirror_box_8,Mirror,"mirror_box_8_log",0,0,0);

        G4VisAttributes *mirror_boxx_8 = new G4VisAttributes();
        mirror_boxx_8->SetForceWireframe(true);
        mirror_box_8_log->SetVisAttributes(mirror_boxx_8);

        G4Box* mirror_box_9 = new G4Box("mirror_box_6",0.5*mm/2,23*mm/2,74*mm/2); // REPLACE END
//======================================
    G4LogicalVolume* mirror_box_9_log
    = new G4LogicalVolume(mirror_box_9,Mirror,"mirror_box_9_log",0,0,0);

    G4VisAttributes *mirror_boxx_9 = new G4VisAttributes();
    mirror_boxx_9->SetForceWireframe(true);
    mirror_box_9_log->SetVisAttributes(mirror_boxx_9);
//========================================================================================================================================================
    // PMT and cathode

    G4double anini = 0.*deg;
    G4double anspan = 360.*deg;

    G4double prin = 0;
    G4double prout = 2.6*cm;
    if (fDetMode == 3){
    prout = 3.81*cm;
    }
    G4double plngth = 3.0*mm/2;

    G4Tubs* pmt = new G4Tubs("PMT",prin,prout,plngth,anini,anspan);

    G4LogicalVolume* pmt_log
    = new G4LogicalVolume(pmt,Quartz,"PMT",0,0,0);

    G4String DetSDname = "tracker1";

    G4double cin = 0;
    G4double cout = 2.6*cm;
    if (fDetMode == 3){
    cout = 3.81*cm;
    }
    G4double clngth = 2.0E-5*mm/2;

    G4Tubs* cath = new G4Tubs("CATH",cin,cout,clngth,anini,anspan);

    G4LogicalVolume* cath_log
    = new G4LogicalVolume(cath,CATH,"CATH",0,0,0);

    G4LogicalVolume* cath_log_c1
    = new G4LogicalVolume(cath,CATH,"CATH",0,0,0);

    qsimDetector* cathSD = new qsimDetector("cath", 2);

    SDman->AddNewDetector(cathSD);

    cath_log->SetSensitiveDetector(cathSD);
    //cath_log_c1->SetSensitiveDetector(cathSD);//6mm tandem quartz

    G4VisAttributes *cathatt = new G4VisAttributes();
    cathatt->SetColour(1.0, 1.0, 0.2);
    cath_log->SetVisAttributes(cathatt);
    cath_log_c1->SetVisAttributes(cathatt);

    // Scintillators
    // Coincidence volumes **** NOTE: Upper scint is below the quartz (First coincidence w/ e-)


    G4Box* upperScint = new G4Box("upperScint",4.5*cm,4.5*cm,0.75*cm);
    G4LogicalVolume* uScint_log = new G4LogicalVolume(upperScint,Air,"upperScint",0,0,0);

    // Make sensitive
    DetSDname = "tracker2";

    qsimScintDetector* upScintSD = new qsimScintDetector(DetSDname, 1);

    SDman->AddNewDetector(upScintSD);
    uScint_log->SetSensitiveDetector(upScintSD);

    G4double upScint_pos;

    upScint_pos = quartz_z-50*cm; //45*cm; // changed to 45 cm from 50 cm as a rough estimate based on CAD measurements of the PMT model 1 + quartz design on the new stand design.




    G4Box* lowScint = new G4Box("lowScint",4.5*cm,4.5*cm,0.75*cm);
    G4LogicalVolume* lScint_log = new G4LogicalVolume(lowScint,Air,"lowScint",0,0,0);

    // Make sensitive
    DetSDname = "tracker3";

    qsimScintDetector* loScintSD = new qsimScintDetector(DetSDname, 2);

    SDman->AddNewDetector(loScintSD);
    lScint_log->SetSensitiveDetector(loScintSD);

    G4double loScint_pos;
    loScint_pos = upScint_pos+1.006*m; //new setup is 1.02874*m; // measured to be 1.02874 m between the two scintillators in the CAD drawings. Previously was just 1.0 m



    // LEAD BLOCK
    ///////////////////////////////////////////////////////////////////////////////////////

    G4Box* Pb_blox = new G4Box("Pb_blox", 10.16*cm,7.62*cm, 10.16*cm);
    //   expanded to ensure nothing
    //   can hit the scint. w/o the lead.
    G4LogicalVolume* Pb_log = new G4LogicalVolume(Pb_blox,Pb_Mat,"Lead",0,0,0);

    G4double Pb_pos;
    Pb_pos = loScint_pos-15.35*cm; // new setup is = loScint_pos-18.554*cm; //(-1*quartz_z)+(30.0*cm-(quartz_y*sin(scintAngle)))*sin(scintAngle);
    // If fAccBeamStand == true then remove the lead bricks, else leave them

    // Detector

    // Mylar (PREX-II)

    G4Trap* mylar_box_1 = new G4Trap("mylar_box_1", 0.05*mm, 2*quartz_z+2*0.01*mm+2*0.05*mm, 2*quartz_y, 2*q_yLB);

    G4LogicalVolume* mylar_box_1_log
    = new G4LogicalVolume(mylar_box_1,Mirror,"mylar_box_1_log",0,0,0);

    G4Box* mylar_box_2 = new G4Box("mylar_box_2", q_yLB + quartz_z, quartz_x + 0.01*mm, 0.05*mm/2);

    G4LogicalVolume* mylar_box_2_log
    = new G4LogicalVolume(mylar_box_2,Mirror,"mylar_box_2_log",0,0,0);

    G4Box* mylar_box_3 = new G4Box("mylar_box_3", q_yLB, quartz_x + 0.01*mm, 0.05*mm/2);

    G4LogicalVolume* mylar_box_3_log
    = new G4LogicalVolume(mylar_box_3,Mirror,"mylar_box_3_log",0,0,0);

    G4Box* mylar_box_4 = new G4Box("mylar_box_4", quartz_z + 0.01*mm + 0.05*mm, quartz_x + 0.01*mm + 0.05*mm, 0.05*mm/2);

    G4LogicalVolume* mylar_box_4_log
    = new G4LogicalVolume(mylar_box_4,Mirror,"mylar_box_4_log",0,0,0);

//------------------

    G4Trap* mylar2_box_1 = new G4Trap("mylar2_box_1", 0.05*mm, 2*quartz2_z+2*0.01*mm+2*0.05*mm, 2*quartz2_y, 2*q_yLB2);

    G4LogicalVolume* mylar2_box_1_log
    = new G4LogicalVolume(mylar2_box_1,Mirror,"mylar2_box_1_log",0,0,0);

    G4Box* mylar2_box_2 = new G4Box("mylar2_box_2", q_yLB2 + quartz2_z, quartz2_x + 0.01*mm, 0.05*mm/2);

    G4LogicalVolume* mylar2_box_2_log
    = new G4LogicalVolume(mylar2_box_2,Mirror,"mylar2_box_2_log",0,0,0);

    G4Box* mylar2_box_3 = new G4Box("mylar2_box_3", q_yLB2, quartz2_x + 0.01*mm, 0.05*mm/2);

    G4LogicalVolume* mylar2_box_3_log
    = new G4LogicalVolume(mylar2_box_3,Mirror,"mylar2_box_3_log",0,0,0);

    G4Box* mylar2_box_4 = new G4Box("mylar2_box_4", quartz2_z + 0.01*mm + 0.05*mm, quartz2_x + 0.01*mm + 0.05*mm, 0.05*mm/2);

    G4LogicalVolume* mylar2_box_4_log
    = new G4LogicalVolume(mylar2_box_4,Mirror,"mylar2_box_4_log",0,0,0);

    //SLAC Test beampipe exit window
    G4Tubs* beampipe_window = new G4Tubs("beampipe_window",0*mm,10*cm,0.127*mm,0*deg,360*deg);

    G4LogicalVolume* beampipe_window_log
    = new G4LogicalVolume(beampipe_window,Steel,"beampipe_window_log",0,0,0);

    qsimScintDetector* WindowSD = new qsimScintDetector("window", 15);

    SDman->AddNewDetector(WindowSD);
    beampipe_window_log->SetSensitiveDetector(WindowSD);

    // Place physical volumes

    G4RotationMatrix* detrot = new G4RotationMatrix;
    G4RotationMatrix* rot_pmt = new G4RotationMatrix;
    G4RotationMatrix* rot_pmt_c1 = new G4RotationMatrix;

    if (fDetMode == 0) {

        detrot->rotateY(45.*deg);
        G4RotationMatrix* rotlg = new G4RotationMatrix;
        rotlg->rotateY(M_PI/2.*rad);
        rotlg->rotateZ(-M_PI/2.*rad);

        rot_pmt->rotateY(M_PI/2.*rad);
        G4VPhysicalVolume* tmirror_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(7.25*cm+lngth+2.*cm,0.,.9*cm),tmirror_log,"TMirror",det_log,false,0);
        G4VPhysicalVolume* lightguide_phys = new G4PVPlacement(rotlg,G4ThreeVector(0.*cm,0,-0.375*cm+.9*cm),lightguide_log,"lightguide_phys", det_log,false,0);
        G4VPhysicalVolume* pmt_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(7.25*cm+plngth+7.*cm,0.,.9*cm),pmt_log,"PMT",det_log,false,0);
        G4VPhysicalVolume* cath_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(7.25*cm+2.*plngth+7.*cm,0.,.9*cm),cath_log,"CATH",det_log,false,0);
    }

    if (fDetMode == 1) {

        detrot->rotateY(fDetAngle);
        rot_pmt -> rotateY(M_PI/4.*rad);

        G4VPhysicalVolume* pmt_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(7.5*cm + 0.01*mm,0.*cm,0.*mm),pmt_log,"PMT",det_log,false,0);
        G4VPhysicalVolume* cath_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(7.5*cm+plngth*cos(M_PI/4*rad)*mm + 0.01*mm+0.15*mm,0.*cm,-plngth*cos(M_PI/4.*rad)*mm),cath_log,"CATH",det_log,false,0);

        G4RotationMatrix* rotMylar_1 = new G4RotationMatrix;
        rotMylar_1->rotateY(0*deg);

        G4VPhysicalVolume* topMylar_1_phys
        = new G4PVPlacement(rotMylar_1,G4ThreeVector(-quartz_z/2,0*mm,-quartz_z-0.05*mm/2-0.01*mm),mylar_box_3_log,"topMylar_1_phys",det_log,false,0);

        G4RotationMatrix* rotMylar_2 = new G4RotationMatrix;
        rotMylar_2->rotateY(0*deg);

        //G4VPhysicalVolume* botMylar_1_phys
        //= new G4PVPlacement(rotMylar_2,G4ThreeVector(quartz_z/2,0*mm,quartz_z+0.05*mm/2+0.01*mm),mylar_box_2_log,"botMylar_1_phys",det_log,false,0);

        G4RotationMatrix* rotMylar_3 = new G4RotationMatrix;
        rotMylar_3->rotateX(90*deg);

        //G4VPhysicalVolume* leftMylar_1_phys
        //= new G4PVPlacement(rotMylar_3,G4ThreeVector(0*mm,-quartz_x - 0.05*mm/2- 0.01*mm,0*mm),mylar_box_1_log,"leftMylar_1_phys",det_log,false,0);

        G4RotationMatrix* rotMylar_4 = new G4RotationMatrix;
        rotMylar_4->rotateX(90*deg);

        //G4VPhysicalVolume* rightMylar_1_phys
        //= new G4PVPlacement(rotMylar_4,G4ThreeVector(0*mm,quartz_x + 0.05*mm/2 +  0.01*mm,0*mm),mylar_box_1_log,"rightMylar_1_phys",det_log,false,0);

        G4RotationMatrix* rotMylar_5 = new G4RotationMatrix;
        rotMylar_5->rotateY(90*deg);

        //G4VPhysicalVolume* frontMylar_1_phys
        //= new G4PVPlacement(rotMylar_5,G4ThreeVector(-quartz_y + quartz_z/2-0.01*mm-0.05*mm/2,0*mm,0*mm),mylar_box_4_log,"rightMylar_1_phys",det_log,false,0);

    }

    if (fDetMode == 2){

        G4RotationMatrix* rotF_1 = new G4RotationMatrix;
    rotF_1->rotateY(0*deg);

        G4VPhysicalVolume* frontPlate_1_phys
        = new G4PVPlacement(rotF_1,G4ThreeVector(-8.03*mm,0.0*mm,0.0*mm),frontPlate_1_log,"FrontPlate_1_phys",det_log,false,0);

        G4RotationMatrix* rotF_2 = new G4RotationMatrix;
    rotF_2->rotateY(0*deg);

        //G4VPhysicalVolume* frontPlate_2_phys
        //= new G4PVPlacement(rotF_2,G4ThreeVector(6.75*mm,0.0*mm,-3.25*mm-0.25*mm),frontPlate_2_log,"FrontPlate_2_phys",det_log,false,0);

        G4RotationMatrix* rotT_1 = new G4RotationMatrix;
        rotT_1->rotateY(0*deg);

    G4VPhysicalVolume* topPlate_1_phys
        = new G4PVPlacement(rotT_1,G4ThreeVector(-2.03*mm,0,-3.53*mm+0.5*mm/2),topPlate_1_log,"TopPlate_1_phys",det_log,false,0);

    G4RotationMatrix* rotT_2 = new G4RotationMatrix;
    rotT_2->rotateY(-22.5*deg);

    G4VPhysicalVolume* topPlate_2_phys
    = new G4PVPlacement(rotT_2,G4ThreeVector(25.13*mm+0.1*mm,0.0*mm,-12.19*mm+0.5*mm/2),topPlate_2_log,"TopPlate_2_phys",det_log,false,0);

    G4RotationMatrix* rotT_3 = new G4RotationMatrix;
    rotT_3->rotateY(0*deg);

    G4VPhysicalVolume* topPlate_3_phys
    = new G4PVPlacement(rotT_3,G4ThreeVector(196.68*mm-0.5*mm,0.0*mm,-20.85*mm+0.5*mm/2),topPlate_3_log,"TopPlate_3_phys",det_log,false,0);

    G4RotationMatrix* rotB_1 = new G4RotationMatrix;
    rotB_1->rotateY(0*deg);
        rotB_1->rotateY(0*deg);

    //G4VPhysicalVolume* botPlate_1_phys
    //= new G4PVPlacement(rotB_1,G4ThreeVector(178.77*mm, 0.0*mm,6.5*mm+0.1*mm+0.05*cm/2),botPlate_1_log,"botPlate_1_phys",det_log,false,0);

    G4VPhysicalVolume* botPlate_1_phys
    = new G4PVPlacement(rotB_1,G4ThreeVector(1.5*mm, 0.0*mm,3.53*mm-0.5*mm/2),botPlate_1_log,"botPlate_1_phys",det_log,false,0);

    G4RotationMatrix* rotB_2 = new G4RotationMatrix;
    rotB_2->rotateY(0*deg);
        rotB_2->rotateY(0*deg);

    G4VPhysicalVolume* botPlate_2_phys
    = new G4PVPlacement(rotB_2,G4ThreeVector(179.05*mm, 0.0*mm, 3.53*mm-0.5*mm/2),botPlate_2_log,"botPlate_2_phys",det_log,false,0);

    G4RotationMatrix* rotR_1 = new G4RotationMatrix;
    rotR_1->rotateY(0*deg);
    rotR_1->rotateX(90*deg);

    G4VPhysicalVolume* RPlate_1_phys
    = new G4PVPlacement(rotR_1,G4ThreeVector(-0.26*mm,-10.28*mm, 0.0*mm),RPlate_1_log,"RPlate_1_phys",det_log,false,0);

    G4RotationMatrix* rotR_2 = new G4RotationMatrix;
        rotR_2->rotateX(90*deg);
    rotR_2->rotateY(180*deg);
    rotR_2->rotateZ(0*deg);

    G4VPhysicalVolume* RPlate_2_phys
    = new G4PVPlacement(rotR_2,G4ThreeVector(187.99*mm-0.25*mm, -10.28*mm, -8.66*mm),RPlate_2_log,"RPlate_2_phys",det_log,false,0);

        G4RotationMatrix* rotR_3 = new G4RotationMatrix;
        rotR_3->rotateX(90*deg);
    rotR_3->rotateY(0*deg);
    rotR_3->rotateZ(0.0*deg);

    G4VPhysicalVolume* RPlate_3_phys
    = new G4PVPlacement(rotR_3,G4ThreeVector(11.28*mm, -10.28*mm, 3.53*mm),RPlate_3_log,"RPlate_3_phys",det_log,false,0);

        G4RotationMatrix* rotR_4 = new G4RotationMatrix;
        rotR_4->rotateX(90*deg);
    rotR_4->rotateY(0*deg);
    rotR_4->rotateZ(0.0*deg);

    //G4VPhysicalVolume* RPlate_4_phys
    //= new G4PVPlacement(rotR_4,G4ThreeVector(8.5*mm-1.75*mm,-10.0*mm-0.05*cm/2-0.1*mm/2,-3*mm+3.5*mm),RPlate_4_log,"RPlate_4_phys",det_log,false,0);

    G4RotationMatrix* rotL_1 = new G4RotationMatrix;
    rotL_1->rotateY(0*deg);
    rotL_1->rotateX(90*deg);

    G4VPhysicalVolume* LPlate_1_phys
    = new G4PVPlacement(rotL_1,G4ThreeVector(-0.26*mm, 10.28*mm, 0.0*mm),LPlate_1_log,"RPlate_phys",det_log,false,0);

    G4RotationMatrix* rotL_2 = new G4RotationMatrix;
        rotL_2->rotateX(90*deg);
    rotL_2->rotateY(180*deg);
    rotL_2->rotateZ(0*deg);

    G4VPhysicalVolume* LPlate_2_phys
    = new G4PVPlacement(rotL_2,G4ThreeVector(187.99*mm-0.25*mm, 10.28*mm, -8.66*mm),LPlate_2_log,"LPlate_2_phys",det_log,false,0);

        G4RotationMatrix* rotL_3 = new G4RotationMatrix;
        rotL_3->rotateX(90*deg);
    rotL_3->rotateY(0*deg);
    rotL_3->rotateZ(0.0*deg);

    G4VPhysicalVolume* LPlate_3_phys
    = new G4PVPlacement(rotL_3,G4ThreeVector(11.28*mm, 10.28*mm, 3.53*mm),LPlate_3_log,"LPlate_3_phys",det_log,false,0);

        G4RotationMatrix* rotL_4 = new G4RotationMatrix;
        rotL_4->rotateX(90*deg);
    rotL_4->rotateY(0*deg);
    rotL_4->rotateZ(0.0*deg);

    //G4VPhysicalVolume* LPlate_4_phys
    //= new G4PVPlacement(rotL_4,G4ThreeVector(8.5*mm-1.75*mm,10.0*mm+0.05*cm/2+0.1*mm/2,-3*mm+3.5*mm),LPlate_4_log,"LPlate_3_phys",det_log,false,0);


    rot_pmt -> rotateY(M_PI/2.*rad);

    G4VPhysicalVolume* pmt_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(347.32*mm+2.5*mm+plngth,0.0*mm,-8.66*mm),pmt_log,"PMT",det_log,false,0);


    G4VPhysicalVolume* cath_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(347.32*mm+2.5*mm+2.*plngth,0.0*mm,-8.66*mm),cath_log,"CATH",det_log,false,0);

    G4RotationMatrix* WRoll = new G4RotationMatrix;
    WRoll->rotateY(0.*deg);

    //G4PVPlacement* W_phys;
        //W_phys = new G4PVPlacement(WRoll,G4ThreeVector(-0.05*cm-0.1*mm-(53/4)*mm+10*mm,0.0*mm,-(13/2)*mm-0.1*mm-0.05*cm-0.05*cm-10*mm/2),W_log,"W",det_log,false,0);



        }



    //
    //
    //
    //
    //
    //                   SHOWERMAX TRANSLATIONS
    //
    //
    //
    //
    //
        if (fDetMode == 3){

//
//
//REPLACE FROM HERE
//
//

        G4RotationMatrix* rotM_1 = new G4RotationMatrix; // Replaces lines 860-938
        rotM_1->rotateY(90*deg);
        rotM_1->rotateX(0*deg);

        G4ThreeVector zTrans_1(113.09*mm,123*mm,-30.75*mm);

        //G4VPhysicalVolume* mirror_box_1_phys
        //= new G4PVPlacement(rotM_1,zTrans_1,mirror_box_1_log,"mirror_box_1_phys",det_log,false,0);



        G4RotationMatrix* rotM_1_1 = new G4RotationMatrix;
        rotM_1_1->rotateY(90*deg);
        rotM_1_1->rotateX(0*deg);

        G4ThreeVector zTrans_1_1(113.09*mm,-123*mm,-30.75*mm);

        //G4VPhysicalVolume* mirror_box_1_1_phys
        //= new G4PVPlacement(rotM_1_1,zTrans_1_1,mirror_box_1_log,"mirror_box_1_1_phys",det_log,false,0);


        G4RotationMatrix* rotM_2 = new G4RotationMatrix;
        rotM_2->rotateY(90*deg);
        rotM_2->rotateX(31.09*deg);

        G4ThreeVector zTrans_2(238.09*mm,78.5*mm,-30.75*mm);

        //G4VPhysicalVolume* mirror_box_2_phys
        //= new G4PVPlacement(rotM_2,zTrans_2,mirror_box_2_log,"mirror_box_2_phys",det_log,false,0);


        G4RotationMatrix* rotM_2_1 = new G4RotationMatrix;
        rotM_2_1->rotateY(90*deg);
        rotM_2_1->rotateX(-31.09*deg);

        G4ThreeVector zTrans_2_1(238.09*mm,-78.5*mm,-30.75*mm);

        //G4VPhysicalVolume* mirror_box_2_1_phys
        //= new G4PVPlacement(rotM_2_1,zTrans_2_1,mirror_box_2_log,"mirror_box_2_1_phys",det_log,false,0);


        G4RotationMatrix* rotM_4 = new G4RotationMatrix;
        rotM_4->rotateZ(90*deg);
        rotM_4->rotateX(-105.51*deg);

        G4ThreeVector zTrans_4(113.09*mm,0*mm,20.46*mm);

        //G4VPhysicalVolume* mirror_box_4_phys
        //= new G4PVPlacement(rotM_4,zTrans_4,mirror_box_4_log,"mirror_box_4_phys",det_log,false,0);


        G4RotationMatrix* rotM_4_1 = new G4RotationMatrix;
        rotM_4_1->rotateZ(90*deg);
        rotM_4_1->rotateX(-74.49*deg);

        G4ThreeVector zTrans_4_1(113.09*mm,0*mm,-81.96*mm);

        //G4VPhysicalVolume* mirror_box_4_1_phys
        //= new G4PVPlacement(rotM_4_1,zTrans_4_1,mirror_box_4_log,"mirror_box_4_1_phys",det_log,false,0);


        G4RotationMatrix* rotM_6 = new G4RotationMatrix;
        rotM_6->rotateZ(90*deg);
        rotM_6->rotateX(-78.35*deg);

        G4ThreeVector zTrans_6(238.09*mm,0*mm,19.46*mm);

        //G4VPhysicalVolume* mirror_box_6_phys
        //= new G4PVPlacement(rotM_6,zTrans_6,mirror_box_6_log,"mirror_box_6_phys",det_log,false,0);


        G4RotationMatrix* rotM_6_1 = new G4RotationMatrix;
        rotM_6_1->rotateZ(90*deg);
        rotM_6_1->rotateX(-101.65*deg);

        G4ThreeVector zTrans_6_1(238.09*mm,0*mm,-80.96*mm);

        //G4VPhysicalVolume* mirror_box_6_1_phys
        //= new G4PVPlacement(rotM_6_1,zTrans_6_1,mirror_box_6_log,"mirror_box_6_1_phys",det_log,false,0); // REPLACE END


//
//
// RPLACE TO HERE
//
//


        G4ThreeVector zTrans_3((289.77/*+17.5*/)*mm,49.79*mm,-27.75*mm);

        G4RotationMatrix* rotM_3 = new G4RotationMatrix;
        rotM_3->rotateY(90*deg);
        rotM_3->rotateX(32.38*deg+180*deg);

        //G4VPhysicalVolume* mirror_box_3_phys
        //= new G4PVPlacement(rotM_3,zTrans_3,mirror_box_3_log,"mirror_box_3_phys",
          //              det_log,false,0);

        G4ThreeVector zTrans_3_1((289.77/*+17.5*/)*mm,-49.79*mm,-27.75*mm);

        G4RotationMatrix* rotM_3_1 = new G4RotationMatrix;
        rotM_3_1->rotateY(90*deg);
        rotM_3_1->rotateX(-32.38*deg+180*deg);

        //G4VPhysicalVolume* mirror_box_3_1_phys
        //= new G4PVPlacement(rotM_3_1,zTrans_3_1,mirror_box_3_log,"mirror_box_3_1_phys",
          //              det_log,false,0);



        G4ThreeVector zTrans_5((231.96/*+17.5*/)*mm,0.0*mm,-82.91*mm);

        G4RotationMatrix* rotM_5 = new G4RotationMatrix;
        rotM_5->rotateZ(90*deg);
        rotM_5->rotateX(-90*deg-13.63*deg);

        //G4VPhysicalVolume* mirror_box_5_phys
        //= new G4PVPlacement(rotM_5,zTrans_5,mirror_box_5_log,"mirror_box_5_phys",
          //              det_log,false,0);

        G4ThreeVector zTrans_5_1((231.96/*+17.5*/)*mm,0.0*mm,27.41*mm);

        G4RotationMatrix* rotM_5_1 = new G4RotationMatrix;
        rotM_5_1->rotateZ(90*deg);
        rotM_5_1->rotateX(-90*deg+13.63*deg);

        //G4VPhysicalVolume* mirror_box_5_1_phys
        //= new G4PVPlacement(rotM_5_1,zTrans_5_1,mirror_box_5_log,"mirror_box_5_1_phys",
          //              det_log,false,0);


// SUITCASE MIRROR TRANSLATIONS ====================================================================================================================================================
//=================================================================================================================================================================================
        G4ThreeVector zTrans_7(6.25*mm/2,25.4*mm/2+0.5*mm/2+0.01*mm,-30.75*mm); // Replaces lines 993-1034

         G4RotationMatrix* rotM_7 = new G4RotationMatrix;
         rotM_7->rotateY(0.0*deg);
         rotM_7->rotateZ(0.0*deg);

         //G4VPhysicalVolume* mirror_box_7_phys
         //= new G4PVPlacement(rotM_7,zTrans_7,mirror_box_7_log,"mirror_box_7_phys",                          det_log,false,0);

        G4ThreeVector zTrans_7_1(6.25*mm/2,-25.4*mm/2-0.5*mm/2-0.01*mm,-30.75*mm);

         G4RotationMatrix* rotM_7_1 = new G4RotationMatrix;
         rotM_7_1->rotateY(0.0*deg);
         rotM_7_1->rotateZ(0.0*deg);

         //G4VPhysicalVolume* mirror_box_7_1_phys
         //= new G4PVPlacement(rotM_7_1,zTrans_7_1,mirror_box_7_log,"mirror_box_7_1_phys",                        det_log,false,0);

        G4ThreeVector zTrans_8(6.25*mm/2,0.0*mm,12.5*mm/2+0.5*mm/2+0.01*mm);

         G4RotationMatrix* rotM_8 = new G4RotationMatrix;
         rotM_8->rotateY(0.0*deg);
         rotM_8->rotateZ(0.0*deg);

         //G4VPhysicalVolume* mirror_box_8_phys
         //= new G4PVPlacement(rotM_8,zTrans_8,mirror_box_8_log,"mirror_box_8_phys",                          det_log,false,0);

        G4ThreeVector zTrans_8_1(6.25*mm/2,0.0*mm,-7*12.5*mm/2-3*8*mm-0.5*mm/2-0.01*mm);

         G4RotationMatrix* rotM_8_1 = new G4RotationMatrix;
         rotM_8_1->rotateY(0.0*deg);
         rotM_8_1->rotateZ(0.0*deg);

        //G4VPhysicalVolume* mirror_box_8_1_phys
        //= new G4PVPlacement(rotM_8_1,zTrans_8_1,mirror_box_8_log,"mirror_box_8_1_phys",                         det_log,false,0);

        G4ThreeVector zTrans_9(-55.725*mm-1*mm,0.0*mm,-30.75*mm); // END REPLACE
//=======================================
        G4RotationMatrix* rotM_9 = new G4RotationMatrix;
        rotM_9->rotateY(0.0*deg);
        rotM_9->rotateZ(0.0*deg);

        //G4VPhysicalVolume* mirror_box_9_phys
        //= new G4PVPlacement(rotM_9,zTrans_9,mirror_box_9_log,"mirror_box_9_phys", det_log,false,0);
//===============================================================================================================================================================================
    rot_pmt -> rotateY(M_PI/2.*rad);
// PMT
 //=============================================================================================================================================================================
        G4VPhysicalVolume* pmt_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(62.88*mm+plngth,0.0*mm,-30.75*mm),pmt_log,"PMT",det_log,false,0); // Replaces lines 1046-1047
        G4VPhysicalVolume* cath_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(62.88*mm+2*plngth+clngth+0.0*mm,0.0*mm,-30.75*mm),cath_log,"CATH",det_log,false,0); // REPLACE END// PMT

        G4OpticalSurface* CTHOpSurface = new G4OpticalSurface("CathodeOpSurface");
        CTHOpSurface -> SetType(dielectric_metal);
        CTHOpSurface -> SetFinish(polished);
        CTHOpSurface -> SetModel(glisur);
        G4MaterialPropertiesTable* COpSurfaceProperty = new G4MaterialPropertiesTable();
        COpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity2,nEntries);
        COpSurfaceProperty -> AddProperty("EFFICIENCY",PhotonEnergy,EfficiencyArray,nEntries);
        CTHOpSurface -> SetMaterialPropertiesTable(COpSurfaceProperty);
        G4LogicalBorderSurface* CathodeSurface =
        new G4LogicalBorderSurface("CathodeSurface",pmt_phys,cath_phys,CTHOpSurface);
//=============================================================================================================================================================================
        }

    double rotation = 0.0*deg;
    double rotation_c1 = 0.0*deg;

    if (fDetMode == 4) {

        detrot->rotateY(fDetAngle);
        rot_pmt -> rotateY(M_PI/4.*rad);
        rot_pmt_c1 -> rotateY(M_PI/4.*rad);

        G4VPhysicalVolume* pmt_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(7.5*cm+3.6*mm,0.*cm,0.0*mm),pmt_log,"PMT",det_log,false,0);//-159.6*mm
        G4VPhysicalVolume* cath_phys = new G4PVPlacement(rot_pmt,G4ThreeVector(7.5*cm+(plngth+clngth+0.0*mm)*cos(M_PI/4*rad)+3.6*mm,0.*cm,-(plngth+clngth+0.0*mm)*sin(M_PI/4.*rad)-0.0*mm),cath_log,"CATH",det_log,false,0);


        G4OpticalSurface* CTHOpSurface = new G4OpticalSurface("CathodeOpSurface");
        CTHOpSurface -> SetType(dielectric_metal);
        CTHOpSurface -> SetFinish(polished);
        CTHOpSurface -> SetModel(glisur);
        G4MaterialPropertiesTable* COpSurfaceProperty = new G4MaterialPropertiesTable();
        COpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity2,nEntries);
        COpSurfaceProperty -> AddProperty("EFFICIENCY",PhotonEnergy,EfficiencyArray,nEntries);
        CTHOpSurface -> SetMaterialPropertiesTable(COpSurfaceProperty);
        G4LogicalBorderSurface* CathodeSurface =
        new G4LogicalBorderSurface("CathodeSurface",pmt_phys,cath_phys,CTHOpSurface);


        G4VPhysicalVolume* pmt_phys_c1 = new G4PVPlacement(rot_pmt_c1,G4ThreeVector(7.5*cm+3.6*mm,0.*cm,0.0*mm),pmt_log,"PMT",det_log_c1,false,0);//-159.6*mm
        G4VPhysicalVolume* cath_phys_c1 = new G4PVPlacement(rot_pmt_c1,G4ThreeVector(7.5*cm+(plngth+clngth+0.0*mm)*cos(M_PI/4*rad)+3.6*mm,0.*cm,-(plngth+clngth+0.0*mm)*sin(M_PI/4.*rad) - 0.0*mm),cath_log_c1,"CATH",det_log_c1,false,0);

        G4LogicalBorderSurface* CathodeSurface_c1 =
        new G4LogicalBorderSurface("CathodeSurface_c1",pmt_phys_c1,cath_phys_c1,CTHOpSurface);

        G4RotationMatrix* rotMylar_1 = new G4RotationMatrix;
        rotMylar_1->rotateY(0*deg);

        //G4VPhysicalVolume* topMylar_1_phys
        //= new G4PVPlacement(rotMylar_1,G4ThreeVector(-quartz_z/2-0.01*mm+0.05*mm,0*mm,-quartz_z-0.05*mm/2-0.01*mm),mylar_box_3_log,"topMylar_1_phys",det_log,false,0);

        G4RotationMatrix* rotMylar_2 = new G4RotationMatrix;
        rotMylar_2->rotateY(0*deg);

        //G4VPhysicalVolume* botMylar_1_phys
        //= new G4PVPlacement(rotMylar_2,G4ThreeVector(quartz_z/2-0.01*mm,0*mm,quartz_z+0.05*mm/2+0.01*mm),mylar_box_2_log,"botMylar_1_phys",det_log,false,0);

        G4RotationMatrix* rotMylar_3 = new G4RotationMatrix;
        rotMylar_3->rotateX(90*deg);

        //G4VPhysicalVolume* leftMylar_1_phys
        //= new G4PVPlacement(rotMylar_3,G4ThreeVector(0*mm,-quartz_x - 0.05*mm/2- 0.01*mm,0*mm),mylar_box_1_log,"leftMylar_1_phys",det_log,false,0);

        G4RotationMatrix* rotMylar_4 = new G4RotationMatrix;
        rotMylar_4->rotateX(90*deg);

        //G4VPhysicalVolume* rightMylar_1_phys
        //= new G4PVPlacement(rotMylar_4,G4ThreeVector(0*mm,quartz_x + 0.05*mm/2 +  0.01*mm,0*mm),mylar_box_1_log,"rightMylar_1_phys",det_log,false,0);

        G4RotationMatrix* rotMylar_5 = new G4RotationMatrix;
        rotMylar_5->rotateY(90.0*deg);

        //G4VPhysicalVolume* frontMylar_1_phys
        //= new G4PVPlacement(rotMylar_5,G4ThreeVector(-quartz_y + quartz_z/2-0.01*mm-0.05*mm/2,0*mm,0*mm),mylar_box_4_log,"rightMylar_1_phys",det_log,false,0);

//-----------

        G4RotationMatrix* rotMylar2_1 = new G4RotationMatrix;
        rotMylar2_1->rotateY(0*deg);

        //G4VPhysicalVolume* topMylar2_1_phys
        //= new G4PVPlacement(rotMylar2_1,G4ThreeVector(-quartz2_z/2-0.01*mm,0*mm,-quartz2_z-0.05*mm/2-0.01*mm-quartz_zPos),mylar2_box_3_log,"topMylar2_1_phys",det_log_c1,false,0);

        G4RotationMatrix* rotMylar2_2 = new G4RotationMatrix;
        rotMylar2_2->rotateY(0*deg);

        //G4VPhysicalVolume* botMylar2_1_phys
        //= new G4PVPlacement(rotMylar2_2,G4ThreeVector(quartz2_z/2-0.01*mm,0*mm,quartz2_z+0.05*mm/2+0.01*mm+quartz_zPos),mylar2_box_2_log,"botMylar2_1_phys",det_log_c1,false,0);

        G4RotationMatrix* rotMylar2_3 = new G4RotationMatrix;
        rotMylar2_3->rotateX(90*deg);

        //G4VPhysicalVolume* leftMylar2_1_phys
        //= new G4PVPlacement(rotMylar2_3,G4ThreeVector(0*mm,-quartz_x - 0.05*mm/2 - 0.01*mm,quartz_zPos),mylar2_box_1_log,"leftMylar2_1_phys",det_log_c1,false,0);

        G4RotationMatrix* rotMylar2_4 = new G4RotationMatrix;
        rotMylar2_4->rotateX(90*deg);

        //G4VPhysicalVolume* rightMylar2_1_phys
        //= new G4PVPlacement(rotMylar2_4,G4ThreeVector(0*mm,quartz_x + 0.05*mm/2 + 0.01*mm,quartz_zPos),mylar2_box_1_log,"rightMylar2_1_phys",det_log_c1,false,0);

        G4RotationMatrix* rotMylar2_5 = new G4RotationMatrix;
        rotMylar2_5->rotateY(90*deg);

        //G4VPhysicalVolume* frontMylar2_1_phys
        //= new G4PVPlacement(rotMylar2_5,G4ThreeVector(-quartz_y + quartz2_z/2-0.01*mm-0.05*mm/2,0*mm,0*mm-quartz_zPos),mylar2_box_4_log,"rightMylar2_1_phys",det_log_c1,false,0);

    }

    G4RotationMatrix* detrot_c1 = new G4RotationMatrix;
    detrot_c1->rotateY(fDetAngle+rotation);

    G4RotationMatrix* detrot_c3 = new G4RotationMatrix;
    detrot_c3->rotateY(fDetAngle+rotation_c1);

    G4RotationMatrix* rot_window = new G4RotationMatrix;
    rot_window -> rotateY(0*rad);

    G4VPhysicalVolume* beampipe_window_phys = new G4PVPlacement(rot_window,G4ThreeVector(0*mm,0.0*mm,-3.9624*m),beampipe_window_log,"BEAMPIPE_WINDOW",world_log,false,0); // SLAC beampipe window

    G4VPhysicalVolume* det_phys
    = new G4PVPlacement(detrot_c1,G4ThreeVector(fDetPosX,fDetPosY,0.*cm),det_log,"detector_phys",world_log,false,0);


    if (fStandMode == 1) {
        G4PVPlacement* uScint_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,upScint_pos-1.*cm),uScint_log,"upperScint",world_log,false,0);
        G4PVPlacement* lScint_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0,loScint_pos),lScint_log,"lowerScint",world_log,false,0);
        G4PVPlacement* Pb_phys = new G4PVPlacement(0,G4ThreeVector(0.0,0.0*cm,Pb_pos),Pb_log,"Pb",world_log,false,0);
    }

    // Surfaces

    // quartz
    G4OpticalSurface* OpQuartzSurface = new G4OpticalSurface("QuartzSurface");
    OpQuartzSurface->SetType(dielectric_dielectric);
    OpQuartzSurface->SetFinish(ground);
    //OpQuartzSurface->SetFinish(polished);
    //OpQuartzSurface->SetModel(glisur);
    //OpQuartzSurface->SetPolish(fQuartzPolish);
    OpQuartzSurface->SetModel(unified);
    OpQuartzSurface->SetSigmaAlpha(0.05);

    //G4LogicalBorderSurface* QuartzSurface =
    //new G4LogicalBorderSurface("QuartzSurface",quartz_phys,det_phys,OpQuartzSurface);

    //Changed from border surface to skin surface
    G4LogicalSkinSurface* QuartzSurface =
    new G4LogicalSkinSurface("QuartzSurface",quartz_log,OpQuartzSurface);

    if (fDetMode == 4) {
    G4VPhysicalVolume* det_phys_c1
    = new G4PVPlacement(detrot_c3,G4ThreeVector(fDetPosX,fDetPosY,-159.6*mm),det_log_c1,"detector_phys_c1",world_log,false,0);
    // quartz optical properties for quart2_phys
    G4OpticalSurface* OpQuartzSurface_c1 = new G4OpticalSurface("QuartzSurface_c1");
    OpQuartzSurface_c1->SetType(dielectric_dielectric);
    OpQuartzSurface_c1->SetFinish(ground);
    //OpQuartzSurface_c1->SetFinish(polished);
    OpQuartzSurface_c1->SetModel(glisur);
    OpQuartzSurface_c1->SetPolish(fQuartzPolish);

    G4LogicalBorderSurface* QuartzSurface2 =
    new G4LogicalBorderSurface("QuartzSurface2",quartz2_phys,det_phys_c1,OpQuartzSurface_c1);
    }

    // mirrors and cathode
    G4OpticalSurface* MOpSurface = new G4OpticalSurface("MirrorOpSurface");
    G4OpticalSurface* MOpSurface_laterals = new G4OpticalSurface("MirrorOpSurface");
    //G4OpticalSurface* CTHOpSurface = new G4OpticalSurface("CathodeOpSurface");
    // quartz window
    G4OpticalSurface* QWOpSurface = new G4OpticalSurface("QuartzWindowOpSurface");

    MOpSurface -> SetType(dielectric_metal);
    MOpSurface -> SetFinish(ground);
    MOpSurface -> SetModel(glisur);
    //MOpSurface -> SetFinish(polished);
    MOpSurface -> SetPolish(0.2);

    MOpSurface_laterals -> SetType(dielectric_metal);
    MOpSurface_laterals -> SetFinish(polished);
    MOpSurface_laterals -> SetModel(glisur);
    MOpSurface_laterals -> SetPolish(1);

    //CTHOpSurface -> SetType(dielectric_metal);
    //CTHOpSurface -> SetFinish(polishedlumirrorair);
    //CTHOpSurface -> SetModel(unified);

    QWOpSurface -> SetType(dielectric_dielectric);
    QWOpSurface -> SetFinish(ground);
    QWOpSurface -> SetModel(glisur);
    QWOpSurface -> SetPolish(0.98);

    const G4int num = 2;
    G4double Ephoton[num] = {2.038*eV, 4.144*eV};

    G4MaterialPropertiesTable* MOpSurfaceProperty = new G4MaterialPropertiesTable();
    G4MaterialPropertiesTable* MOpSurfaceProperty_laterals = new G4MaterialPropertiesTable();
    //G4MaterialPropertiesTable* COpSurfaceProperty = new G4MaterialPropertiesTable();
    G4MaterialPropertiesTable* TubeSurfaceProperty = new G4MaterialPropertiesTable();

    G4MaterialPropertiesTable* QWOpSurfaceProperty = new G4MaterialPropertiesTable();

    MOpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity1,nEntries);

    MOpSurface -> SetMaterialPropertiesTable(MOpSurfaceProperty);

    MOpSurfaceProperty_laterals -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity_laterals,nEntries);

    MOpSurface_laterals -> SetMaterialPropertiesTable(MOpSurfaceProperty_laterals);

    //COpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity2,nEntries);
    //COpSurfaceProperty -> AddProperty("EFFICIENCY",PhotonEnergy,EfficiencyArray,nEntries);


    //CTHOpSurface -> SetMaterialPropertiesTable(COpSurfaceProperty);

    //QWOpSurfaceProperty -> AddProperty("REFLECTIVITY",PhotonEnergy,QWReflectivity,nEntries);
    QWOpSurfaceProperty -> AddProperty("RINDEX",PhotonEnergy,RefractiveIndex1,nEntries);
    QWOpSurfaceProperty -> AddProperty("ABSLENGTH",PhotonEnergy,Absorption1,nEntries);

    if (fDetMode == 0) {
        G4LogicalSkinSurface* TubeSurface_1 = new
        G4LogicalSkinSurface("TubeMirrorOpS_1",tmirror_log,MOpSurface);
        G4LogicalSkinSurface* lightguideSurface = new
        G4LogicalSkinSurface("lightguideOps",lightguide_log,MOpSurface);
    }

    if (fDetMode == 1) {
        G4LogicalSkinSurface* MylarSurface_1 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar_box_1_log,MOpSurface);

        G4LogicalSkinSurface* MylarSurface_2 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar_box_2_log,MOpSurface);

        G4LogicalSkinSurface* MylarSurface_3 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar_box_3_log,MOpSurface);

        G4LogicalSkinSurface* MylarSurface_4 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar_box_4_log,MOpSurface);

    }

    if (fDetMode == 2) {

    G4LogicalSkinSurface* FrontSurface_1 = new
    G4LogicalSkinSurface("FrontMirrorOpS_1",frontPlate_1_log,MOpSurface);

    G4LogicalSkinSurface* TopSurface_1 = new
    G4LogicalSkinSurface("TopMirrorOpS_1",topPlate_1_log,MOpSurface);

    G4LogicalSkinSurface* TopSurface_2 = new
    G4LogicalSkinSurface("TopMirrorOpS_2",topPlate_2_log,MOpSurface_laterals);

    G4LogicalSkinSurface* TopSurface_3 = new
    G4LogicalSkinSurface("TopMirrorOpS_3",topPlate_3_log,MOpSurface_laterals);

    G4LogicalSkinSurface* BotSurface_1 = new
    G4LogicalSkinSurface("BotMirrorOpS_1",botPlate_1_log,MOpSurface);

    G4LogicalSkinSurface* BotSurface_2 = new
    G4LogicalSkinSurface("BotMirrorOpS_2",botPlate_2_log,MOpSurface_laterals);

    G4LogicalSkinSurface* LSurface_1 = new
    G4LogicalSkinSurface("LMirrorOpS_1",LPlate_1_log,MOpSurface);

        G4LogicalSkinSurface* LSurface_2 = new
    G4LogicalSkinSurface("LMirrorOpS_2",LPlate_2_log,MOpSurface_laterals);

        G4LogicalSkinSurface* LSurface_3 = new
    G4LogicalSkinSurface("LMirrorOpS_3",LPlate_3_log,MOpSurface_laterals);

    G4LogicalSkinSurface* RSurface_1 = new
    G4LogicalSkinSurface("RMirrorOpS_1",RPlate_1_log,MOpSurface);

        G4LogicalSkinSurface* RSurface_2 = new
    G4LogicalSkinSurface("RMirrorOpS_2",RPlate_2_log,MOpSurface_laterals);

        G4LogicalSkinSurface* RSurface_3 = new
    G4LogicalSkinSurface("RMirrorOpS_3",RPlate_3_log,MOpSurface_laterals);

    }

    if (fDetMode == 3) {

        G4LogicalSkinSurface* MSurface_1 = new
    G4LogicalSkinSurface("MirrorOpS_1",mirror_box_1_log,MOpSurface);

        G4LogicalSkinSurface* MSurface_2 = new
    G4LogicalSkinSurface("MirrorOpS_2",mirror_box_2_log,MOpSurface);

        G4LogicalSkinSurface* MSurface_3 = new
    G4LogicalSkinSurface("MirrorOpS_3",mirror_box_3_log,MOpSurface);

        G4LogicalSkinSurface* MSurface_4 = new
    G4LogicalSkinSurface("MirrorOpS_4",mirror_box_4_log,MOpSurface);

        G4LogicalSkinSurface* MSurface_5 = new
    G4LogicalSkinSurface("MirrorOpS_5",mirror_box_5_log,MOpSurface);

        G4LogicalSkinSurface* MSurface_6 = new
    G4LogicalSkinSurface("MirrorOpS_6",mirror_box_6_log,MOpSurface);

        G4LogicalSkinSurface* MSurface_7 = new
    G4LogicalSkinSurface("MirrorOpS_7",mirror_box_7_log,MOpSurface);

        G4LogicalSkinSurface* MSurface_8 = new
    G4LogicalSkinSurface("MirrorOpS_8",mirror_box_8_log,MOpSurface);

        G4LogicalSkinSurface* MSurface_9 = new
    G4LogicalSkinSurface("MirrorOpS_9",mirror_box_9_log,MOpSurface);

    }

    if (fDetMode == 4) {
        G4LogicalSkinSurface* MylarSurface_1 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar_box_1_log,MOpSurface);

        G4LogicalSkinSurface* MylarSurface_2 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar_box_2_log,MOpSurface);

        G4LogicalSkinSurface* MylarSurface_3 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar_box_3_log,MOpSurface);

        G4LogicalSkinSurface* MylarSurface_4 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar_box_4_log,MOpSurface);
//-------
        G4LogicalSkinSurface* Mylar2Surface_1 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar2_box_1_log,MOpSurface);

        G4LogicalSkinSurface* Mylar2Surface_2 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar2_box_2_log,MOpSurface);

        G4LogicalSkinSurface* Mylar2Surface_3 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar2_box_3_log,MOpSurface);

        G4LogicalSkinSurface* Mylar2Surface_4 = new
        G4LogicalSkinSurface("MylarOpS_1",mylar2_box_4_log,MOpSurface);

    }

    //G4LogicalSkinSurface* QuartzWindowSurface = new
    //G4LogicalSkinSurface("QuartzWindowOpS",pmt_log,QWOpSurface);

    //G4LogicalSkinSurface* CathSurface = new
    //G4LogicalSkinSurface("CathOpS1", cath_log,CTHOpSurface);

    //G4LogicalSkinSurface* CathSurface_c1 = new
    //G4LogicalSkinSurface("CathOpS2", cath_log_c1,CTHOpSurface);

    // Generate & Add Material Properties Table attached to the optical surfaces

    //OpticalQuartzSurface
    G4double RefractiveIndex[num] = {1.46, 1.46};
    //G4double Reflectivity[num]    = {1, 1};
    G4double SpecularLobe[num]    = {0.3, 0.3};//{0.3, 0.3}
    G4double SpecularSpike[num]   = {0.2, 0.2};//{0.2, 0.2}
    G4double Backscatter[num]     = {0.2, 0.2};//{0.2, 0.2}

    G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

    myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
    //myST1->AddProperty("REFLECTIVITY",          PhotonEnergy, Reflectivity,    num);
    myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
    myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
    myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);

    OpQuartzSurface->SetMaterialPropertiesTable(myST1);

    //always return the physical World
    return world_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
