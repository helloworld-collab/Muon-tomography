/* Tanguy Pierog and Ralf Ulrich  18.12.2020
This is a simple template to read a CORSIKA DAT File and print part of the input
using COAST libraries. It is meant as a simple example that user should adapt
to their use. Most of the variables have been indicated but more can be found
in COAST doxygen documentation (in particular for the missing (here) longitudinal
profiles, Cherenkov or muon information)*/

#include <crsRead/MCorsikaReader.h>

#include <crs/TSubBlock.h>
#include <crs/MRunHeader.h>
#include <crs/MEventHeader.h>
#include <crs/MEventEnd.h>
#include <crs/MParticleBlock.h>
#include <crs/MLongitudinalBlock.h>
#include <crs/MParticle.h>

#include <TFile.h>
#include <TH1D.h>
#include <TNtupleD.h>

#include <iostream>
#include <sstream>
#include <map>
using namespace std;


// to hold data of one observation level
struct ObsLevel {
  TNtupleD* data;
  double x;
  double y;
  double x2;
  double y2;
  double w;
};
          



int
main (int argc, char **argv) 
{  
  if (argc<2) {
    cout << "  please specify Corsika file " << endl;
    return 1;
  }
    
  string fname(argv[1]);
  crsRead::MCorsikaReader cr(fname, 3);
  
  string inputFile;
  if (fname.rfind('/')==string::npos) {
    inputFile = fname;
  } else {
    inputFile = fname.substr(fname.rfind('/')+1, 
                             fname.size()-fname.rfind('/'));
  }
  
  ostringstream oFileName;
  oFileName << inputFile.c_str() << ".root";

  cout << " Writing summary to output file: " << oFileName.str() << endl;

  TFile oFile (oFileName.str ().c_str (), "RECREATE");
      
  // create TTree for global event properties

  double Eprim,Xmax,Nmax,X0,Nph,Nel,Nmu,Nha,Ntot,Zfirst,zenith,azimuth;
  static TTree * Event;
  Event = new TTree("Event","Event Properties");
  Event->Branch("Eprim",&Eprim,"Eprim/D");    // Theta (zenith) angle
  Event->Branch("Theta",&zenith,"Theta/D");    // Theta (zenith) angle
  Event->Branch("Phi",&azimuth,"Phi/D");    // Phi (azimuth) angle
  Event->Branch("Xmax",&Xmax,"Xmax/D");    // position of maximum
  Event->Branch("Nmax",&Nmax,"Nmax/D");    // height of maximum
  Event->Branch("X0",&X0,"X0/D");    // position of first interaction from GH fit
  Event->Branch("Zfirst",&Zfirst,"Zfirst/D");    // real position of first interaction
  Event->Branch("Nph",&Nph,"Nph/D");    // Number of photon in all observation levels
  Event->Branch("Nel",&Nel,"Nel/D");    // Number of electrons+positrons in all observation levels
  Event->Branch("Nmu",&Nmu,"Nmu/D");    // Number of muons in all observation levels
  Event->Branch("Nha",&Nha,"Nha/D");    // Number of hadrons in all observation levels
  Event->Branch("Ntot",&Ntot,"Ntot/D");    // Number of particles in all observation levels
  
  int ShowerCounter = 0;
  
  crs::MRunHeader Run;
  while (cr.GetRun (Run)) {

    const int nRun = Run.GetRunID();
      cout << "---------------------------------\n"
           << " Run info:\n"
           << "  run number = " << nRun << "\n";
      
    
                  /* Full list of variables available : 
CREAL 	GetRunID () const 
CINT    GetDateStart () const
CREAL 	GetVersion () const 
CINT 	GetNObservationLevels () const 
CREAL 	GetObservationHeight (int index) const 
CREAL 	GetSpectralSlope () const 
CREAL 	GetEMin () const 
CREAL 	GetEMax () const 
CREAL 	GetFlagEGS4 () const 
CREAL 	GetFlagNKG () const 
CREAL 	GetCutoffHadrons () const 
CREAL 	GetCutoffMuons () const 
CREAL 	GetCutoffElectrons () const 
CREAL 	GetCutoffPhotons () const 
CREAL 	GetSamplingPlanePointX () const 
CREAL 	GetSamplingPlanePointY () const 
CREAL 	GetSamplingPlanePointZ () const 
CREAL 	GetSamplingPlaneTheta () const 
CREAL 	GetSamplingPlanePhi () const 
double 	GetSamplingPlaneNormalX () const 
double 	GetSamplingPlaneNormalY () const 
double 	GetSamplingPlaneNormalZ () const 
CREAL 	GetRotationAngle () const 
CINT 	GetNumberOfShowers () const 
CREAL 	GetAtmosphereLayerBoundary (int index) const 
CREAL 	GetAtmosphereA (int index) const 
CREAL 	GetAtmosphereB (int index) const 
CREAL 	GetAtmosphereC (int index) const 
CREAL 	GetVerticalDepth (const CREAL heightAboveSeaLevel) const
CREAL   GetConstNFLAIN () const 
CREAL   GetConstNFLDIF () const 
CREAL   GetConstNFLPI0 () const 
CREAL   GetConstNFLPIF () const 
CREAL   GetConstNFLCHE () const 
CREAL   GetConstNFRAGM () const 
 */
   
    crs::MEventHeader Shower;
    while (cr.GetShower(Shower)) {
      
      ++ShowerCounter;
      /*      
      ostringstream oFileName;
      oFileName << inputFile.c_str() << "_"
                << Shower.GetEventNumber () << ".root";

      cout << " Writing summary to output file: " << oFileName.str() << endl;

      TFile oFile (oFileName.str ().c_str (), "RECREATE");
      */
      const int nObsLevel = Shower.GetNObservationLevels();
      map<int, ObsLevel> obsLevel;
      
      for (int iObsLevel=1; iObsLevel<=nObsLevel; ++iObsLevel) { 
        
        double height = Shower.GetObservationHeight(iObsLevel-1);            
        cout << " init obs-level " << iObsLevel << " at h=" << height << "cm" << endl;
        ObsLevel emptyLevel;
        ostringstream tTitle, tName;
        tTitle << "Data at level " << iObsLevel;
        tName << "data_" << iObsLevel;
        emptyLevel.data = new TNtupleD(tName.str().c_str(), tTitle.str().c_str(), "id:ek:px:py:pz:x:y:t:w:theta");
        
          emptyLevel.x  = 0;
          emptyLevel.y  = 0;
          emptyLevel.w  = 0;
          emptyLevel.x2 = 0;
          emptyLevel.y2 = 0;
        
        obsLevel[iObsLevel] = emptyLevel;

      } // end loop observation levels

      Eprim = Shower.GetEnergy();
      zenith = Shower.GetTheta();
      azimuth = Shower.GetPhi();
      Zfirst = Shower.GetZFirst();
      cout << "---------------------------------\n"
           << " Shower with:\n"
           << "  Zfirst = " << Zfirst
           << "  Theta = " << zenith
           << "  Phi = " << azimuth << "\n";
      
      
                  /* Full list of variables available : 
CINT 	GetEventNumber () const 
CREAL 	GetParticleId () const 
CREAL 	GetEnergy () const 
CREAL   GetStartingAltitude () const
CREAL 	GetFirstTarget () const 
CREAL   GetZFirst () const
CREAL   GetPx () const   
CREAL   GetPy () const   
CREAL   GetPz () const   
CREAL   GetTheta () const
CREAL   GetPhi () const
CINT 	GetNRandomSequences () const 
CINT 	GetSeed (int index) const 
CINT 	GetInitialCallsMod (int index) const 
CINT 	GetInitialCallsDiv (int index) const 
CREAL 	GetRunNumber () const 
CINT    GetDateStart () const
CREAL 	GetVersion () const 
CINT 	GetNObservationLevels () const 
CREAL 	GetObservationHeight (int index) const 
CREAL 	GetSpectralSlope () const 
CREAL   GetEMin () const
CREAL   GetEMax () const
CREAL   GetCutoffHadrons () const
CREAL   GetCutoffMuons () const
CREAL   GetCutoffElectrons () const
CREAL   GetCutoffPhotons () const
CREAL 	GetNFLAIN () const 
CREAL 	GetNFLDIF () const 
CREAL 	GetNFLPI0 () const 
CREAL 	GetNFLPIF () const 
CREAL 	GetNFLCHE () const 
CREAL 	GetNFRAGM () const 
CREAL   GetBx () const
CREAL   GetBz () const
CREAL 	GetFlagEGS4 () const 
CREAL 	GetFlagNKG () const 
CREAL 	GetHadronicLowEModell () const 
CREAL 	GetHadronicHighEModell () const 
CREAL 	GetFlagCherenkov () const 
CREAL 	GetFlagNeutrino () const 
CREAL 	GetFlagCurved () const 
CINT 	GetFlagComputer () const 	1: IBM, 2: Transputer, 3: DEC/UNIX, 4: Mac, 5: VAX/VMS, 6: GNU/Linux
CREAL   GetThetaMin () const    
CREAL   GetThetaMax () const    
CREAL   GetPhiMin () const        
CREAL   GetPhiMax () const        
CREAL 	GetCherenkovBunch () const 
CREAL 	GetCherenkovNumberX () const 
CREAL 	GetCherenkovNumberY () const 
CREAL 	GetCherenkovGridX () const 
CREAL   GetCherenkovGridY () const
CREAL 	GetCherenkovDetectorX () const 
CREAL   GetCherenkovDetectorY () const
CREAL   GetCherenkovOutputFlag () const
CREAL 	GetArrayRotation () const 
CREAL 	GetFlagExtraMuonInformation () const 
CREAL 	GetMultipleScatteringStep () const 
CREAL 	GetCherenkovBandwidthMin () const 
CREAL 	GetCherenkovBandwidthMax () const 
CINT 	GetNUsesOfEvent () const 
CREAL 	GetCherenkovCoreX (int index) const 
CREAL 	GetCherenkovCoreY (int index) const 
CREAL 	GetFlagSIBYLL () const 
CREAL 	GetFlagSIBYLLCross () const 
CREAL 	GetFlagQGSJET () const 
CREAL 	GetFlagQGSJETCross () const 
CREAL 	GetFlagDPMJET () const 
CREAL 	GetFlagDPMJETCross () const 
CREAL 	GetFlagVENUSCross () const 
CREAL 	GetFlagMuonMultiple () const 
CREAL 	GetNKGRadialRange () const 
CREAL 	GetEFractionThinningH () const 
CREAL 	GetEFractionThinningEM () const 
CREAL 	GetWMaxHadronic () const 
CREAL 	GetWMaxEM () const 
CREAL 	GetRMaxThinning () const 
CREAL 	GetInnerAngle () const 
CREAL 	GetOuterAngle () const 
CREAL 	GetTransitionEnergy () const 
CREAL 	GetSkimmingIncidence () const 
CREAL 	GetSkimmingAltitude () const 
CREAL 	GetStartingHeight () const
CREAL   GetVersionFLUKA () const                   */
      
      crs::TSubBlock Data;
      while (cr.GetData (Data)) {
        
        switch (Data.GetBlockType ()) {
          
            case crs::TSubBlock::ePARTDATA:
            {
              const crs::MParticleBlock& ParticleData = Data;
              crs::MParticleBlock::ParticleListConstIterator iEntry;
              for (iEntry = ParticleData.FirstParticle();
                   iEntry != ParticleData.LastParticle();
                   ++iEntry) {

		// DUMP
		//iEntry->Dump();
                
                if (iEntry->IsParticle()) {
                  
                  crs::MParticle iPart(*iEntry);
                  
                  const int id    = iPart.GetParticleID();
                  const int level = iPart.GetObservationLevel();
                  const double w  = iPart.GetWeight();
                  const double ek  = iPart.GetKinEnergy();
                  const double px = iPart.GetPx();
                  const double py = iPart.GetPy();
                  const double pz = iPart.GetPz();
                  const double x  = iPart.GetX();
                  const double y  = iPart.GetY();
                  const double t  = iPart.GetTime();
                  const double theta = iPart.GetTheta();

                  /* Full list of variables available : 
CREAL 	GetWeight () const 
virtual std::string GetParticleName () const 
int 	GetParticleID () const 
int 	GetParticleId () const 
bool 	IsParticle () const 
bool 	IsNucleus () const 
bool 	IsCherenkov () const 
bool 	IsMuonProductionInfo () const 
bool 	IsEmpty () const 
ParticleType 	GetType () const 
virtual int 	GetObservationLevel () const 
virtual int 	GetHadronicGeneration () const 
virtual CREAL 	GetPx () const 
virtual CREAL 	GetPy () const 
virtual CREAL 	GetPz () const 
virtual CREAL 	GetX () const 
virtual CREAL 	GetY () const 
virtual CREAL 	GetTime () const 
double 	GetMass () const 	mass in GeV
int 	GetPDGCode () const 
double 	GetKinEnergy () const 	kin. energy in GeV
double 	GetTheta () const 	zenith angle in rad
double 	GetPSquared () const 	squared of momentum in GeV
double 	GetE () const 
                  */
                
		  if (obsLevel.count(level)==0) {
		   cout << " detected new obs-level " << level 
			 << ". Possibly INCLIN-Option " << endl;
		    ObsLevel emptyLevel;
		    ostringstream tTitle, tName;
		    tTitle << "Data at level " << level;
		    tName << "data_" << level;
		    obsLevel[level] = emptyLevel;
		    obsLevel[level].data = new TNtupleD(tName.str().c_str(), 
							tTitle.str().c_str(),
							"id:ek:px:py:pz:x:y:t:w:theta");
		  }
		  
                  obsLevel[level].data->Fill(id,ek,px,py,pz,x,y,t,w,theta);
                  obsLevel[level].x  += x*w;
                  obsLevel[level].y  += y*w;
                  obsLevel[level].x2 += x*x*w;
                  obsLevel[level].y2 += y*y*w;
                  obsLevel[level].w  += w;
                  
                }
                
              } // end particle loop
              
              break;
            }
            
            case crs::TSubBlock::eLONG:
              break;
              
            default:
              break;
        } // end data block
        
      } // loop data
      
      crs::MEventEnd ShowerSummary;
      cr.GetShowerSummary(ShowerSummary);
      Xmax = ShowerSummary.GetXmax();
      Nmax = ShowerSummary.GetNmax();
      X0 = ShowerSummary.GetX0();
      Nph = ShowerSummary.GetPhotons();
      Nel = ShowerSummary.GetElectrons();
      Nmu = ShowerSummary.GetMuons();
      Nha = ShowerSummary.GetHadrons();
      Ntot = ShowerSummary.GetParticles();
      Event->Fill();
      
                  /* Full list of variable available : 

CREAL 	GetEventNumber () const 
CREAL 	GetPhotons () const 
CREAL 	GetElectrons () const 
CREAL 	GetHadrons () const 
CREAL 	GetMuons () const 
CREAL 	GetParticles () const 
CREAL 	GetNmax () const 	Longitudinal distribution (see manual p. 64)
CREAL 	GetX0 () const 
CREAL 	GetXmax () const 
CREAL 	GetLongi_a () const 
CREAL 	GetLongi_b () const 
CREAL 	GetLongi_c () const 
CREAL 	GetLongi_Chi2 () const 
CREAL 	GetWeightedPhotons () const 	Added according to the CORSIKA manual.
CREAL 	GetWeightedElectrons () const 
CREAL 	GetWeightedHadrons () const 
CREAL 	GetWeightedMuons () const 
CREAL 	GetNPreshower () const                  */
      
      /*cout << "---------------------------------\n"
           << " Shower info:\n"
           << "  Xmax = " << Xmax << "\n";*/
      
      oFile.cd();
      for (map<int, ObsLevel>::iterator iLevel = obsLevel.begin();
           iLevel != obsLevel.end();
           ++iLevel) {
        
        
        double npart=iLevel->second.w;
        if(npart>0){
        cout << "   observation level: " << iLevel->first 
	     << " with " << npart << " particles "
	     << "    <x>= " << iLevel->second.x / npart
             << " +/- " << sqrt(iLevel->second.x2/npart-pow(iLevel->second.x/npart,2))/max(1.,npart-1.)
             << " cm <y>= " << iLevel->second.y /npart
             << " +/- " << sqrt(iLevel->second.y2/npart-pow(iLevel->second.y/npart,2))/max(1.,npart-1.)
	     << " cm \n";
     	iLevel->second.data->Write();
        }                  
        
         
      } // loop observation levels
      
      //      oFile.Close();
      
    } // loop shower
    
  } // loop runs (usually just 1)    

  oFile.cd();
  Event->Write();
  oFile.Close();
  cout << " Read " << ShowerCounter << " showers from file " << endl;
  
  return 0;
}


  
