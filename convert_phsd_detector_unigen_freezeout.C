//************************************************************************************************************************************************
/** author: Susanne Glaessel (Universitaet Frankfurt)
 ** Macro for creating the detector input from PHSD-PHQMD output with freeze-out
 ** coordinates without clusters
 **
 ** This macro converts the PHQMD output-files for the final timestep to the
 ** UniGen format including the FREEZEOUT-TIME & -POSITION. The UniGen-output
 ** can be selected with the flag: WriteUnigen = kTRUE.
 ** (Note: Momentum in UniGen-ouput is the final momentum, not the momentum at 
 ** freeze-out.)
 ** 
 ** To write out the FREEZEOUT-MOMENTUM in addition another output format needs to  
 ** be used. It can be selected with the flag: WriteEventFreeze = kTRUE. 
 **
 ** The output is:
 ** - Unigen: [dataset].phsd.root
 ** - root-file including freeze-out momentum: [dataset].phsd_freeze.root
 **
 ** inputfiles are: 
 ** inputPHSD           : input information for PHQMD simulation
 ** phsd.dat            : hadrons (at last timestep)
 **
 ** ProcessIds for channels with deuterons are changed to make them positive 3-digits:
 ** PHQMD processId -> UniGen fDecay: 1101 -> 701; -1101 -> 801; 1301 -> 703; -1301 -> 803.
 **
 **/             

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>

struct PHadron {
  Int_t fPdgId;
  TVector3 fP;
  Float_t fEnergy;
  Int_t fProcessId;    
  Int_t fParentId; 
  Int_t fBaryonId;
  TLorentzVector fXTFreeze;
  TLorentzVector fPEFreeze;
  PHadron() : fPdgId(0), fEnergy(0.), fProcessId(-1), fParentId(-1), fBaryonId(-1) { fP.SetXYZ(0.,0.,0.); fXTFreeze.SetXYZT(0.,0.,0.,0.); fPEFreeze.SetXYZT(0.,0.,0.,0.);  };
  PHadron(Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Int_t processId, Int_t parentId, Int_t baryonId, Float_t xposfo, Float_t yposfo, Float_t zposfo , Float_t timefo, Float_t xpfo, Float_t ypfo, Float_t zpfo, Float_t energyfo) : fPdgId(pdgId), fEnergy(energy), fProcessId(processId), fParentId(parentId), fBaryonId(baryonId) { fP.SetXYZ(Px,Py,Pz); fXTFreeze.SetXYZT(xposfo,yposfo,zposfo,timefo); fPEFreeze.SetXYZT(xpfo,ypfo,zpfo,energyfo);  };
};

class PEvent : public TObject  {
public:
  Int_t fEventId;
  Int_t fNParticipants;
  Float_t fB;
  Float_t fPhi;
  vector<PHadron> fhadrons;
  ClassDef(PEvent, 1);
};

class ParticleFreeze : public TObject {
public:
  Int_t fIndex;
  Int_t fPdgId;
  Int_t fParent;
  Int_t fDecay;
  TVector3 fP;
  Float_t fEnergy;
  Float_t fTimeFreeze;
  TVector3 fXFreeze;
  TVector3 fPFreeze;
  Float_t  fEnergyFreeze; 
  Int_t fOrigin;
  Int_t fWeight;
  ParticleFreeze() : fIndex(-1), fPdgId(0), fParent(-1), fDecay(-1), fEnergy(0.), fTimeFreeze(0.), fOrigin(0), fWeight(0) {fP.SetXYZ(0.,0.,0.); fXFreeze.SetXYZ(0.,0.,0.); fPFreeze.SetXYZ(0.,0.,0.); };
  ParticleFreeze(Int_t index, Int_t pdgId, Int_t parent, Int_t decay, TVector3 P, Float_t energy, TLorentzVector XTFreeze, TLorentzVector PEFreeze, Int_t origin, Int_t weight) : fIndex(index), fPdgId(pdgId), fParent(parent), fDecay(decay), fP(P), fEnergy(energy), fOrigin(origin), fWeight(weight) {fTimeFreeze = XTFreeze.T(); fXFreeze.SetXYZ(XTFreeze.X(), XTFreeze.Y(), XTFreeze.Z()); fPFreeze.SetXYZ(PEFreeze.X(), PEFreeze.Y(), PEFreeze.Z()); fEnergyFreeze = PEFreeze.T();};
  ParticleFreeze(Int_t index, Int_t pdgId, Int_t parent, Int_t decay, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Float_t TimeFreeze, TVector3 XFreeze, TVector3 PFreeze, Float_t energyFreeze, Int_t origin, Int_t weight) : fIndex(index), fPdgId(pdgId), fParent(parent), fDecay(decay), fEnergy(energy), fTimeFreeze(TimeFreeze), fXFreeze(XFreeze), fPFreeze(PFreeze), fEnergyFreeze(energyFreeze), fOrigin(origin), fWeight(weight) {fP.SetXYZ(Px, Py, Pz); };
  ClassDef(ParticleFreeze, 1);
};

class EventFreeze : public TObject  {
public:
  Int_t fEventId;
  Float_t fB;
  Int_t fNParticipants;
  Float_t fTime;
  Float_t fPhi;
  Int_t fNpa;
  vector<ParticleFreeze> fParticles;
  void SetParameters(Int_t eventId, Float_t b, Int_t nParticipants, Float_t time, Float_t phi) { fEventId = eventId; fB = b; fTime = time; fPhi = phi; fNParticipants = nParticipants;};
  void AddParticle(Int_t index, Int_t pdgId, Int_t parent, Int_t decay, TVector3 P, Float_t energy, TLorentzVector XTFreeze, TLorentzVector PEFreeze, Int_t origin, Int_t weight) { fParticles.push_back(ParticleFreeze(index, pdgId, parent, decay, P, energy, XTFreeze, PEFreeze, origin, weight)); fNpa += 1; }; 
  void AddParticle(Int_t index, Int_t pdgId, Int_t parent, Int_t decay, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Float_t timefo, TVector3 posfo, TVector3 pfo, Float_t energyfo, Int_t origin, Int_t weight) { fParticles.push_back(ParticleFreeze(index, pdgId, parent, decay, Px, Py, Pz, energy, timefo, posfo, pfo, energyfo, origin, weight)); fNpa += 1; }; 
  void Clear() { fParticles.clear(); fNpa = 0; };
  ClassDef(EventFreeze, 1);
};

Int_t ChangeProcessId3digits(Int_t processId) {

  /** Change processIds for channels with deuterons to make them positive 3-digits **/
  
  Int_t processId_3digits = processId;
  if (TMath::Abs(processId) > 999) {
    if ( processId == 1101 ) processId_3digits = 701; 
    if ( processId == -1101) processId_3digits = 801;
    if ( processId == 1301 ) processId_3digits = 703;
    if ( processId == -1301) processId_3digits = 803; 
  }
  return processId_3digits;
}

void convert_phsd_detector_unigen_freezeout(TString indir = "",
					     TString dataset = "",
					     Int_t firstevent = 0,
		      			     Bool_t WriteUnigen = kTRUE,
	      				     Bool_t WriteEventFreeze = kTRUE)
{

  cout << "********************************************************************" <<endl;
  cout << "Start run " << dataset << endl;
  cout << "********************************************************************" <<endl;
  
  // -----   In- and output file names   ------------------------------------
 
  TString inputFileInfo = Form("%s/inputPHSD",indir.Data());
  TString inputFileBulk = Form("%s/%s/phsd.dat",indir.Data(),dataset.Data());

  TString outdir = Form("%s/unigen",indir.Data());	 

  TString rootFileTmp = Form("%s/%s.phsd_freeze_tmp.root",outdir.Data(),dataset.Data());
  TString rootFileDet = Form("%s/%s.phsd.root",outdir.Data(),dataset.Data());
  TString rootFileFreeze;
  if (WriteEventFreeze == kTRUE) {
    rootFileFreeze = Form("%s/%s.phsd_freeze.root",outdir.Data(),dataset.Data());
  }

  cout << endl;
  cout << "Input files are: " << endl;
  cout << inputFileInfo << endl;
  cout << inputFileBulk << endl;

  // ------------------------------------------------------------------------

  Int_t NUM, ISUBS, NTIME;
  Int_t aProj, zProj, aTarg, zTarg;  
  Double_t eLab, bMin, bMax;
    
  FILE *inputInfo = fopen(inputFileInfo, "r");
  
  fscanf(inputInfo, "%i %*[^\n]%*c", &aTarg);
  fscanf(inputInfo, "%i %*[^\n]%*c", &zTarg);
  fscanf(inputInfo, "%i %*[^\n]%*c", &aProj);
  fscanf(inputInfo, "%i %*[^\n]%*c", &zProj);
  fscanf(inputInfo, "%lf %*[^\n]%*c", &eLab);
  fscanf(inputInfo, "%lf %*[^\n]%*c", &bMin);
  fscanf(inputInfo, "%lf %*[^\n]%*c", &bMax);
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%i %*[^\n]%*c", &NUM);
  fscanf(inputInfo, "%i %*[^\n]%*c", &ISUBS);
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%*[^\n]%*c");
  fscanf(inputInfo, "%i %*[^\n]%*c", &NTIME);
  
  Int_t nEvents=NUM*ISUBS;
  cout << endl;
  cout << "Conversion of " << nEvents << " events with " << NTIME << " timesteps" <<endl;
  cout << endl; 

  const Float_t Ebin_max = 0.0;
  const Float_t kProtonMass = 0.938272321;
  Float_t pProj = TMath::Sqrt(eLab*kProtonMass/2);
  Float_t pTarg = -pProj;
  Float_t energyCM = TMath::Sqrt((eLab + (2 * kProtonMass)) * (2 * kProtonMass));
  Float_t pBeam = 2 * pProj / TMath::Sqrt(1 - pProj*pProj / (pProj*pProj + kProtonMass*kProtonMass));

  cout << "-----------Collision System-------------------------------- "                           << endl;
  cout << endl;
  cout << "Mass of target (GeV/c*c)         : " << aTarg                           << endl;
  cout << "Number of protons in target      : " << zTarg                           << endl;
  cout << "Mass of projectile (GeV/c*c)     : " << aProj                           << endl;
  cout << "Number of protons in projectile  : " << zProj                           << endl;
  cout << "Centre of mass energy (GeV)      : " << setprecision(3) << energyCM     << endl;
  cout << "Laboratory energy (AGeV)         : " << setprecision(3) << eLab         << endl;
  cout << "Beam momentum (AGeV/c)           : " << setprecision(3) << pBeam        << endl;
  cout << "Projectile momentum (AGeV/c)     : " << setprecision(3) << pProj        << endl;
  cout << "Target momentum (AGeV/c)         : " << setprecision(3) << pTarg        << endl;
  cout << "Impact parameter b (fm)          : " << setprecision(2) << bMin << "-"  << setprecision(2) << bMax << endl; 
  cout << "----------------------------------------------------------- "                           << endl;
  cout << endl;

  TFile *outputTmp = new TFile (rootFileTmp, "recreate");
  TTree *treeTmp = new TTree ("events", "events");
  PEvent *event = new PEvent();
  treeTmp->Branch ("event", "PEvent", &event, 12800000);

  FILE *BulkFile = fopen(inputFileBulk, "r");
  
  for (int isub = 0; isub < ISUBS; isub++) {  // loop over all subsequent runs
    for (int irun = 0; irun < NUM; irun ++) {  // loop over all parallel runs
      outputTmp->cd();
      event->fEventId = firstevent + isub*NUM + irun;
	  
      Int_t nPart = 0;
      Int_t pdgId, charge, baryonId, processId, parentId;
      Float_t Px, Py, Pz, energy, xposfo, yposfo, zposfo, timefo, xpfo, ypfo, zpfo, energyfo;

      //Get Hadrons from phsd.dat
      if(fscanf(BulkFile, "%i %*i %*i %f %*i %*i %*i %*f %*f\n", &nPart, &event->fB)==EOF)
	throw runtime_error("Unexpected end of file phsd.dat at run " + to_string(irun));
      if(fscanf(BulkFile, "%i %f %*[^\n]%*c", &event->fNParticipants, &event->fPhi)==EOF)   
	throw runtime_error("Unexpected end of file phsd.dat at run " + to_string(irun));
	  
      outputTmp->cd();
      event->fhadrons.clear();
	  
      for (int i = 0; i < nPart; i++) {

	if(fscanf(BulkFile, "%i %i %f %f %f %f %i %i %i %*f %*f %*f %f %f %f %f %f %f %f %f %*f %*f\n", &pdgId, &charge, &Px, &Py, &Pz, &energy, &processId, &parentId, &baryonId, &xposfo, &yposfo, &zposfo, &timefo, &xpfo, &ypfo, &zpfo, &energyfo)==EOF) {	
	  throw runtime_error("Unexpected end of file phsd.dat at run " + to_string(irun) + " particle " + to_string(i));
	}
	if((pdgId<1000 && pdgId>-1000) || TMath::Abs(pdgId) == 100121) baryonId = -1;
	if(TMath::Abs(pdgId) == 100121) pdgId = 1000010020*charge; // correct pdg-code for kinetic deuterons

	if (TMath::Abs(processId) > 999)
	  processId = ChangeProcessId3digits(processId); // change processIds for channels with deuterons to make them positive 3-digits
	 	    
	outputTmp->cd();
	event->fhadrons.push_back(PHadron(pdgId, Px, Py, Pz, energy, processId, parentId, baryonId, xposfo, yposfo, zposfo ,timefo, xpfo,ypfo,zpfo, energyfo));
      }
      outputTmp->cd();
      treeTmp->Fill();
    }  // end loop parallel runs
  }  // end loop subsequent runs

  int check_eof;
  if(fscanf(BulkFile, "%i %*[^\n]%*c", &check_eof) != EOF) 
    throw runtime_error("\n  Error when reading " + inputFileBulk + ": File not read until the end. Check input format.\n ");
  
  outputTmp->cd();
  treeTmp->Write();

  fclose(BulkFile);

  //********Fill final trees *****************
  
  URun *header; UEvent *uevent; TFile *output; TTree *tree;
 
  if (WriteUnigen == kTRUE) {
    header = new URun ("phqmd", "without clusters", aProj, zProj, pProj, aTarg, zTarg, pTarg, bMin, bMax, 0, 0, 0, 0, nEvents);
    uevent = new UEvent;
    output = new TFile (rootFileDet, "recreate");
    tree = new TTree ("events", "signal");
    header->Write();
    tree->Branch ("event", "UEvent", uevent);
  }
  
  TFile *outputFreeze; TTree *treeFreeze; EventFreeze *eventFreeze;  
  
  if (WriteEventFreeze == kTRUE) {
    outputFreeze = new TFile (rootFileFreeze, "recreate");
    treeFreeze = new TTree ("events", "events");
    eventFreeze = new EventFreeze();
    treeFreeze->Branch ("event", "EventFreeze", &eventFreeze, 12800000);
  }

  for (int ievent = 0 ; ievent < treeTmp->GetEntries(); ievent++) {

    outputTmp->cd();
    treeTmp->GetEntry(ievent);

    Int_t index = 0 ;
    Int_t parentId = -1;
    Int_t child[2] = {0,0};
      
    if (WriteUnigen == kTRUE) {
      output->cd();
      uevent->Clear();
      uevent->SetParameters(event->fEventId, event->fB, event->fPhi, event->fNParticipants, 1, 0.0, 0);
    }
      
    if (WriteEventFreeze == kTRUE) {
      outputFreeze->cd();
      eventFreeze->Clear();
      eventFreeze->SetParameters(event->fEventId, event->fB, event->fNParticipants, 0.0, event->fPhi);
    }
      
    outputTmp->cd(); // Loop over hadrons
    for (auto hadron : event->fhadrons) { 	

      if (WriteUnigen == kTRUE) {
	output->cd();
	uevent->AddParticle (index, hadron.fPdgId, 0, hadron.fParentId, -1, -1, hadron.fProcessId, child, hadron.fP.X(), hadron.fP.Y(),hadron.fP.Z(), hadron.fEnergy, hadron.fXTFreeze.X(), hadron.fXTFreeze.Y(), hadron.fXTFreeze.Z(), hadron.fXTFreeze.T(), 1);
      }
	  
      if (WriteEventFreeze == kTRUE) {
	outputFreeze->cd();
	eventFreeze->AddParticle(index, hadron.fPdgId, hadron.fParentId, hadron.fProcessId, hadron.fP, hadron.fEnergy, hadron.fXTFreeze, hadron.fPEFreeze, 0, 1);
      }
      index++;    		  
    } // end loop hadrons
      
    if (WriteUnigen == kTRUE) { 
      output->cd();
      tree->Fill();
    }
    if (WriteEventFreeze == kTRUE) {
      outputFreeze->cd();
      treeFreeze->Fill();
    }	
  } 

  outputTmp->Close();
  gSystem->Unlink(rootFileTmp);
  
  if (WriteUnigen == kTRUE) {
    output->cd();
    tree->Write();
    output->Close();
  }
  if (WriteEventFreeze == kTRUE) {
    outputFreeze->cd();
    treeFreeze->Write();
    outputFreeze->Close();
  }
  
  cout << endl;
  cout << "Macro finished successfully." << endl;
  if (WriteUnigen == kTRUE) cout << "Output file is: " <<rootFileDet<< endl;
  if (WriteEventFreeze == kTRUE) cout << "Output file including freezeout momentum is: " <<rootFileFreeze<< endl;
}
  
  
