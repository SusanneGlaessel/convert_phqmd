//************************************************************************************************************************************************
/** author: Susanne Glaessel (Universitaet Frankfurt)
 ** Macro for creating the detector input from PHSD-PHQMD output with freeze-out
 ** coordinates  
 **
 ** This macro converts the PHQMD output-files after stabilization for the final 
 ** timestep to the UniGen format including the FREEZEOUT-TIME & -POSITION. The 
 ** UniGen-output can be selected with the flag: WriteUnigen = kTRUE.
 ** (Note: Momentum in UniGen-ouput is the final momentum, not the momentum at 
 ** freeze-out.)
 ** 
 ** To write out the FREEZEOUT-MOMENTUM in addition another output format needs to  
 ** be used. It can be selected with the flag: WriteEventFreeze = kTRUE. 
 **
 ** In PHQMD clusters / anticlusters are recognised independently of their physical 
 ** existence. The cluster-baryons are listed separately in the outputfile. This 
 ** routine builds clusters from the single baryons based on their cluster-ID and 
 ** identifies physical clusters according to the cluster_table.root. Baryons from 
 ** unphysical clusters (eg. p-p) are counted as single baryons.
 ** With option "ConvertAllClusters": Clusters with A > 7 are counted independent of 
 ** their physical existince.
 **
 ** PHQMD writes baryons and anti-baryons into two separate files. The conversion 
 ** of anti-baryons is optional and can be switched off with option 
 ** ConvertAnti = kFALSE.
 **
 ** The output is: 
 ** - for clusters & anticlusters: [dataset].phqmd.root
 ** - for clusters only: [dataset].phqmd_noanti.root
 **
 ** inputfiles are: 
 ** inputPHSD           : input information for PHQMD simulation
 ** fort.891            : baryons that are entering the cluster routine  
 **                       (for every timestep)
 ** fort.881 (optional) : antibaryons that are entering the cluster routine 
 ** phsd.dat            : hadrons (at last timestep)
 ** cluster_table.dat   : contains a list with nuclei & hypernuclei with up to 
 **                       9 baryons with their baryon content and branching ratio for  
 **                       the identification of physical clusters
 **
 ** Deuterons are not only produced with the MST-algorithm, but also through kinetic
 ** interactions. The information of the origin can be taken from the UniGen-output:
 ** UParticle* particle -> GetStatus(): = 0 kinetic deuteron, = 1 MST deuteron
 **
 ** ProcessIds and ParentIdsfor clusters: The UniGen variables fDecay and fParent allow  
 ** to store 3 Ids with3 digits each. 3 cluster-baryons to be stored are selected according  
 ** to the following rule: First all Ids for Sigma0s are stored, then for Lambdas, protons  
 ** and neutrons.
 **
 ** ProcessIds for channels with deuterons are changed to make them positive 3-digits:
 ** PHQMD processId -> UniGen fDecay: 1101 -> 701; -1101 -> 801; 1301 -> 703; -1301 -> 803.
 **
 ** ParentIds are changed into 3-digits by keeping only the last 3-digits (The digits 
 ** specific to the event are removed.).
 **
 ** Inputfiles are from PHQMD MST-mode, for SACA-mode, input-files need to be 
 ** replaced by fort.893 & fort.883
 **/             

#include "TROOT.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>

std::map<int,int[30]> eventId2Entry;
std::vector<std::map<int,int>> baryonId2pos;
std::vector<std::map<int,int>> baryons2hadrons;
std::vector<std::map<int,vector<int>>> clusterId2baryonIds;
map<int, vector<float>> eventId2time;

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

struct PBaryon {
  Int_t fPdgId;
  Int_t fnBary;
  TVector3 fP;
  TVector3 fX;
  Float_t fMass;
  Float_t fEnergy;
  Int_t fBaryonId;
  Int_t fClusterId;
  Float_t fEbin;
  Float_t fTimeFreeze;
  PBaryon() : fPdgId(0),fnBary(0),fMass(0.),fEnergy(0.),fBaryonId(0),fClusterId(0),fEbin(0), fTimeFreeze(0.) { fP.SetXYZ(0.,0.,0.); fX.SetXYZ(0.,0.,0.); };
  PBaryon(Int_t PdgId, Int_t nBary, Float_t Px, Float_t Py, Float_t Pz, Float_t X, Float_t Y, Float_t Z, Float_t Mass, Int_t baryonId, Int_t clusterId, Float_t TimeFreeze, Float_t Ebin) : fPdgId(PdgId), fnBary(nBary), fMass(Mass), fBaryonId(baryonId), fClusterId(clusterId), fTimeFreeze(TimeFreeze), fEbin(Ebin) { fP.SetXYZ(Px,Py,Pz); fX.SetXYZ(X,Y,Z); fEnergy = TMath::Sqrt(Mass*Mass +Px*Px + Py*Py + Pz*Pz); };
};

struct PBaryon_cluster {
  Int_t fPdgId;
  TVector3 fP;
  Float_t fMass;
  Float_t fEnergy;
  Float_t fEbin;
  Float_t fTimeFreeze;
  TLorentzVector fXTFreeze;
  TLorentzVector fPEFreeze;
  Int_t fBaryonId;
  Int_t fProcessId;
  Int_t fParentId;
  PBaryon_cluster(Int_t baryonId, Int_t PdgId, Int_t processId, Int_t parentId, TVector3 P, Float_t mass, Float_t energy, TLorentzVector XTFreeze, TLorentzVector PEFreeze, Float_t Ebin) : fBaryonId(baryonId), fPdgId(PdgId), fProcessId(processId), fParentId(parentId), fP(P), fMass(mass), fEnergy(energy), fXTFreeze(XTFreeze), fPEFreeze(PEFreeze), fEbin(Ebin) {};
};

class PEvent : public TObject  {
public:
  Int_t fEventId;
  Int_t fStepNr;
  Int_t fNParticipants;
  Float_t fB;
  Float_t fTime;
  Float_t fPhi;
  vector<PHadron> fhadrons;
  vector<PBaryon> fbaryons;
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

struct ClusterEntry {
  Int_t fPdgId;
  Int_t fNProt;
  Int_t fNBary0;
  Int_t fNLamb;
  Int_t fNSigm;
  Double_t fBR;
  ClusterEntry() : fPdgId(-1), fNProt(-1), fNBary0(-1), fNLamb(-1), fNSigm(-1), fBR(-1) {};
  ClusterEntry(Int_t pdgId, Int_t nProt, Int_t nBary0, Int_t nLamb, Int_t nSigm, Double_t br) : fPdgId(pdgId), fNProt(nProt), fNBary0(nBary0), fNLamb(nLamb), fNSigm(nSigm), fBR(br) {};
};

void GetPdgIdBaryon(Int_t charge, Int_t &pdgId, Bool_t IsAnti)
{
  if(TMath::Abs(charge) == 1) pdgId = 2212;
  if(charge == 0) pdgId = 2112;
  if(TMath::Abs(charge) == 17) pdgId = 3122;
  if (TMath::Abs(charge) == 18) pdgId = 3212;

  if (IsAnti == kTRUE) pdgId *= -1;
}

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

void GetClusterList(TString clustertablename, std::vector<ClusterEntry> &clusterlist)
{
  cout << "Clustertable used: " << clustertablename << endl;
  
  char clustername [80];
  Int_t pdgId, nProt, nBary0, nLamb, nSigm;
  Double_t br;

  std::map<Int_t, std::vector<Double_t>> content2br;
  content2br.clear();

  clusterlist.clear();

  FILE *fClusterTable = fopen(clustertablename, "r");

  printf("------------------------------------\n");
  printf("name             pdg            br\n");
  printf("------------------------------------\n");
  while(1) {
    if (fscanf(fClusterTable, "%s %i %i %i %i %i %lf\n", clustername, &pdgId, &nProt, &nBary0, &nLamb, &nSigm, &br)==EOF) break;
    clusterlist.push_back(ClusterEntry(pdgId, nProt, nBary0, nLamb, nSigm, br));
    printf("%-16s %-10i %8.2f\n", clustername, pdgId, br);
    Int_t content = 1000000000 + (nBary0 + nProt) * 10 + nProt * 10000 + nLamb * 10000000 + nSigm * 100000000;
    if (br < 1.0) content2br[content].push_back(br);
  }
  printf("------------------------------------\n");

  for (auto content : content2br) {
    Double_t br_total = 0;
    for (int i = 0; i < content.second.size(); i++) {
      br_total += content.second[i];
    }
    if (br_total != 1) {
      Int_t pdgA   = 10; Int_t pdgZ = 10000; Int_t pdgL = 10000000; 
      throw runtime_error("\n Clustertable: Sum of branching ratios is " + to_string(br_total) + " != 1 for: A = " + to_string((content.first % pdgZ) / pdgA) + ", Z = " + to_string((content.first % pdgL) / pdgZ) + ", L = " + to_string((content.first % (10 * pdgL)) / pdgL) + ", S = " + to_string((content.first % (100 * pdgL)) / (pdgL*10)) + "!");
    }
  }
  
  fclose(fClusterTable);
}

void GetClusterPdg(std::vector<ClusterEntry> clusterlist, std::vector<PBaryon_cluster> baryons_cluster, Int_t clusterId, Bool_t ConvertAllClusters, Int_t &pdgIdCl)
{
  Int_t nProtCl = 0; Int_t nBary0Cl = 0; Int_t nLambCl = 0; Int_t nSigmCl = 0;
  Int_t nbary = baryons_cluster.size();
  for (int ibary = 0; ibary < nbary; ibary++) {
    Int_t pdgId = baryons_cluster.at(ibary).fPdgId;
    if (TMath::Abs(pdgId) == 2212) nProtCl ++;
    if (TMath::Abs(pdgId) == 2112 || TMath::Abs(pdgId) == 3122 || TMath::Abs(pdgId) == 3212) nBary0Cl ++;
    if (TMath::Abs(pdgId) == 3122) nLambCl ++;
    if (TMath::Abs(pdgId) == 3212) nSigmCl ++;
  }
   
  std::map<Int_t, Double_t> pdg2br;
  pdg2br.clear();
  Double_t br_total = 0;
    
  for (auto cluster : clusterlist) {
    if (nProtCl == cluster.fNProt && nBary0Cl == cluster.fNBary0 && nLambCl == cluster.fNLamb && nSigmCl == cluster.fNSigm) {
      if (cluster.fBR == 1) {
	br_total = cluster.fBR;
	pdgIdCl = TMath::Sign(cluster.fPdgId, clusterId);
	break;
      }
      else {
	pdg2br [cluster.fPdgId] = cluster.fBR;
	br_total += cluster.fBR;
	if (br_total == 1) {
	  Double_t br_integral = 0;
	  Double_t rndm = gRandom->Rndm();
	  for (auto pdg : pdg2br) {
	    br_integral += pdg.second;
	    if (rndm <= br_integral) {
	      pdgIdCl = TMath::Sign(pdg.first, clusterId);
	      break;
	    }
	  }
	  break;
	}
      }	
    }
    else {
      pdgIdCl = 99999;
    }
  }

  if (ConvertAllClusters == kTRUE && pdgIdCl == 99999 && nbary > 7) {                    
    pdgIdCl = 1000000000 + nbary * 10 + nProtCl * 10000 + (nLambCl+nSigmCl) * 10000000;  
    pdgIdCl = TMath::Sign(pdgIdCl, clusterId);                                                                                  
  } 
}

Float_t CalculateClusterBindingEnergy(std::vector<PBaryon_cluster> baryons_cluster)
{
  Float_t Ebin = 0.;
  for (int ibary = 0; ibary < baryons_cluster.size(); ibary++) {
    Ebin += baryons_cluster.at(ibary).fEbin;
  }
  return Ebin; 
}

void CalculateClusterKin(std::vector<PBaryon_cluster> baryons_cluster, Float_t &Px, Float_t &Py, Float_t &Pz, Float_t &energy)
{
  Px = 0.; Py = 0.; Pz = 0.; energy = 0;
  Float_t Mass = 0.;
  for (int ibary = 0; ibary < baryons_cluster.size(); ibary++) {
    Int_t pdgId = baryons_cluster.at(ibary).fPdgId;
    Px += baryons_cluster.at(ibary).fP.X();
    Py += baryons_cluster.at(ibary).fP.Y();
    Pz += baryons_cluster.at(ibary).fP.Z();
    Mass += baryons_cluster.at(ibary).fMass;
  }
  energy = TMath::Sqrt(Mass*Mass+Px*Px+Py*Py+Pz*Pz);
}

void GetClusterProcessParentId(std::vector<PBaryon_cluster> baryons_cluster, Int_t &processId, Int_t &parentId)
{  
  const Int_t npos = 3;
  Int_t nbary = baryons_cluster.size();

  std::vector<PBaryon_cluster> baryons_cluster_sorted =  baryons_cluster;

  for (int ibary = 1; ibary < nbary ; ibary++) {  
    PBaryon_cluster baryon_tmp = baryons_cluster_sorted.at(ibary);
    Int_t pdgId_tmp = baryons_cluster_sorted.at(ibary).fPdgId;
    Int_t ibary_tmp = ibary;
    for (int jbary = ibary - 1; jbary >= 0; jbary --) {
      if (pdgId_tmp > baryons_cluster_sorted.at(jbary).fPdgId) {
	baryons_cluster_sorted.at(ibary_tmp) = baryons_cluster_sorted.at(jbary);
	baryons_cluster_sorted.at(jbary) = baryon_tmp;	
	ibary_tmp --; 
      }
    }
  }
   
  std::vector<Int_t> poswrite;
  poswrite.resize(nbary);
    
  if (nbary <= 3) {
    for (int ibary = 0; ibary < nbary ; ibary++)
      poswrite.at(ibary) = 1;
  }
  else {
    std::array<Int_t, 4> nspecies = {0, 0, 0, 0}; // number of sigma0, lambda, proton, neutron
    for (int ibary = 0; ibary < nbary; ibary++) {
      Int_t pdgId = baryons_cluster.at(ibary).fPdgId;
      if (TMath::Abs(pdgId) == 3212) nspecies.at(0) ++;
      if (TMath::Abs(pdgId) == 3122) nspecies.at(1) ++;
      if (TMath::Abs(pdgId) == 2212) nspecies.at(2) ++;
      if (TMath::Abs(pdgId) == 2112) nspecies.at(3) ++;
    }
    Int_t nfree = npos; Int_t ipos_last = 0;
    for (int ispec = 0; ispec < 4; ispec ++) {
      if (nspecies.at(ispec) <= nfree) {
	nfree -= nspecies.at(ispec);
	for (int ipos = ipos_last; ipos < ipos_last + nspecies.at(ispec); ipos++) {
	  poswrite.at(ipos) = 1;
	}
	ipos_last += nspecies.at(ispec);
	if (nfree == 0) {
	  for (int ipos = ipos_last; ipos < nbary; ipos++)
	    poswrite.at(ipos) = 0;
	  break;
	}
      }
      else {
	for (int ipos = ipos_last; ipos < ipos_last + nspecies.at(ispec); ipos++) 
	  poswrite.at(ipos) = 0;
	ipos_last += nspecies.at(ispec);
      }

    }
  }

  Int_t ipos = 0;
  processId = 1e9; parentId = 1e9;
  for (int ibary = 0; ibary < nbary; ibary++) {
    if (poswrite.at(ibary) == 1) {
      processId += TMath::Abs(baryons_cluster_sorted.at(ibary).fProcessId) * 1e6 / TMath::Power(10,3*ipos);
      Int_t parentId_long = baryons_cluster_sorted.at(ibary).fParentId;
      Int_t parentId_short = (parentId_long % 1000); //use only last 3 digits to store parentIds of 3 baryons
      parentId  += parentId_short * 1e6 / TMath::Power(10,3*ipos);
      ipos ++;
    }
  }

  if (ipos < npos) {
  Int_t ipos_last = ipos;
  for (int ipos = ipos_last; ipos < npos; ipos++)
    processId += 999 * 1e6 / TMath::Power(10,3*ipos);
  }
}

void CalculateClusterProductionTime(Int_t eventId, Int_t final_step, Int_t clusterId, Int_t nbary, Float_t &TimeProductionCluster)
{
  /** Calculates formation time of a cluster. **/
 
  auto it_ievent = eventId2Entry.find(eventId);
  Int_t ievent_final = it_ievent->second[final_step]; 
  auto it_cluster = clusterId2baryonIds[ievent_final].find(clusterId);
  
  std::vector<int> cluster_baryons_last;
  cluster_baryons_last.clear();
  for (auto baryId : it_cluster->second) 
    cluster_baryons_last.push_back(baryId);

  Float_t time = eventId2time.find(eventId)->second.at(final_step);
    
  for (int istep = final_step; istep >= 0; istep --) {
    Int_t ievent_ts = it_ievent->second[istep];
    
    if (clusterId2baryonIds[ievent_ts].find(clusterId) != clusterId2baryonIds[ievent_ts].end()) {
   
      auto it = clusterId2baryonIds[ievent_ts].find(clusterId);
      Int_t cluster_size_ts = it->second.size();

      if (cluster_size_ts != nbary) break;

      Bool_t foundId = kFALSE;
      for (int baryIdlast : cluster_baryons_last) {
	foundId = kFALSE;
	for (int baryIdTs : it->second) {
	  if (baryIdlast == baryIdTs) foundId = kTRUE;
	}
	if (foundId == kFALSE) break;
      }
      if (foundId == kFALSE) break;
      time = eventId2time.find(eventId)->second.at(istep);
    }
  }
  TimeProductionCluster = time;
}

void CalculateClusterFreezeOutTime(Int_t eventId, std::vector<PBaryon_cluster> baryons_cluster, Int_t nbary, Float_t TimeProductionCluster, Int_t &TsFreeze, Float_t &TimeFreezeCluster, Float_t &deltaT)
{  
  std::vector<float> timefo_bary;
  for (int ibary = 0; ibary < nbary; ibary ++)
    timefo_bary.push_back(baryons_cluster.at(ibary).fXTFreeze.T());

  timefo_bary.push_back(TimeProductionCluster);
  TimeFreezeCluster = *max_element(timefo_bary.begin(), timefo_bary.end());

  auto it_TsFreeze_up = std::lower_bound(eventId2time.find(eventId)->second.begin(), eventId2time.find(eventId)->second.end(), TimeFreezeCluster);
  TsFreeze = std::distance(eventId2time.find(eventId)->second.begin(), it_TsFreeze_up); // first timestep after freezeout-time
  if (it_TsFreeze_up == eventId2time.find(eventId)->second.end()) TsFreeze --;
  deltaT = TimeFreezeCluster - eventId2time.find(eventId)->second.at(TsFreeze);
}

void convert_phqmd_detector_unigen_freezeout(TString indir = "",
					     TString dataset = "",
					     Int_t firstevent = 0,
					     Bool_t ConvertAllClusters = kTRUE,
				      	     Bool_t ConvertAnti = kTRUE,
		      			     Bool_t WriteUnigen = kTRUE,
	      				     Bool_t WriteEventFreeze = kTRUE)
{

  cout << "********************************************************************" <<endl;
  cout << "Start run " << dataset << endl;
  cout << "********************************************************************" <<endl;
  
  // -----   In- and output file names   ------------------------------------


  TString inputFileInfo = Form("%s/inputPHSD",indir.Data());
  TString inputFileBulk = Form("%s/%s/phsd.dat",indir.Data(),dataset.Data());
  TString inputFileBaryonFriga = Form("%s/%s/fort.891",indir.Data(),dataset.Data());
  TString inputFileBaryonFrigaAnti = Form("%s/%s/fort.881",indir.Data(),dataset.Data());

  TString outdir;
  if (ConvertAllClusters == kTRUE)
    outdir = Form("%s/unigen/allclusters",indir.Data());	
  else
    outdir = Form("%s/unigen/smallclusters",indir.Data()); 	  

  TString rootFileTmp =  Form("%s/%s.phqmd_freeze_tmp.root",outdir.Data(),dataset.Data());
  TString rootFileDet;
  if (WriteUnigen == kTRUE) {
    if (ConvertAnti == kTRUE) rootFileDet = Form("%s/%s.phqmd.root",outdir.Data(),dataset.Data());
    if (ConvertAnti == kFALSE) rootFileDet = Form("%s/%s.phqmd_noanti.root",outdir.Data(),dataset.Data());
  }
  TString rootFileFreeze;
  if (WriteEventFreeze == kTRUE) {
    if (ConvertAnti == kTRUE) rootFileFreeze = Form("%s/%s.phqmd_freeze.root",outdir.Data(),dataset.Data());
    if (ConvertAnti == kFALSE) rootFileFreeze = Form("%s/%s.phqmd_freeze_noanti.root",outdir.Data(),dataset.Data());
  }

  TString clustertable = Form("%s/cluster_table.dat",outdir.Data());

  cout << endl;
  cout << "Input files are: " << endl;
  cout << inputFileInfo << endl;
  cout << inputFileBulk << endl;
  cout << inputFileBaryonFriga << endl;
  if (ConvertAnti == kTRUE)  cout << inputFileBaryonFrigaAnti << endl;

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
  FILE *BaryonFrigaFile = fopen(inputFileBaryonFriga, "r");
  FILE *BaryonFrigaFileAnti;
  if (ConvertAnti == kTRUE) BaryonFrigaFileAnti = fopen(inputFileBaryonFrigaAnti, "r");
  
  for (int isub = 0; isub < ISUBS; isub++) {  // loop over all subsequent runs
    for (int it = 0; it < NTIME+1; it++) { // loop over all timesteps
      for (int irun = 0; irun < NUM; irun ++) {  // loop over all parallel runs
	outputTmp->cd();
	event->fEventId = firstevent + isub*NUM + irun;
	event->fStepNr = it;
	  
	Int_t nPart = 0;
	Int_t pdgId, charge, baryonId, processId, parentId;
	Float_t Px, Py, Pz, energy, xposfo, yposfo, zposfo, timefo, xpfo, ypfo, zpfo, energyfo;

	//Get Hadrons from phsd.dat

	if (it == NTIME) {
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
	} 

	// Get Baryons Friga from fort.891
	
	Int_t nBaryEntries, nBary, clusterId;
	Float_t time, Xpos, Ypos, Zpos, Mass, TimeFreeze, Ebin;

	Bool_t IsAnti = kFALSE;

	if(fscanf(BaryonFrigaFile, "%*i %*i %*f %f %*[^\n]%*c", &time)==EOF)
	  throw runtime_error("Unexpected end of file fort.891 at run " + to_string(irun));
 	if(fscanf(BaryonFrigaFile, "%*[^\n]%*c")==EOF)
 	  throw runtime_error("Unexpected end of file fort.891 at run " + to_string(irun));
 	if(fscanf(BaryonFrigaFile, "%i %*[^\n]%*c", &nBaryEntries)==EOF)
	  throw runtime_error("Unexpected end of file fort.891 at run " + to_string(irun));

	outputTmp->cd();
	event->fTime = time;
	event->fbaryons.clear();
	
	for (int i = 0; i < nBaryEntries; i++) {
	  if(fscanf(BaryonFrigaFile, "%*i %i %f %f %f %f %f %f %f %i %i %i %*i %*i %*i %*i %f %f\n", &charge, &Px, &Py, &Pz, &Xpos, &Ypos, &Zpos, &Mass, &clusterId, &nBary, &baryonId, &TimeFreeze, &Ebin)==EOF) {
	    throw runtime_error("Unexpected end of file fort.891 at run " + to_string(irun) + "timestep" + to_string(it) + " particle " + to_string(i));
	  }

	  GetPdgIdBaryon(charge, pdgId, IsAnti);
	  outputTmp->cd();
	  event->fbaryons.push_back(PBaryon(pdgId, nBary, Px, Py, Pz, Xpos, Ypos, Zpos, Mass, baryonId, clusterId, TimeFreeze, Ebin));
	} 

	// Get Anti-Baryons Friga from fort.881
	
	if (ConvertAnti == kTRUE) {

	  IsAnti = kTRUE;

	  if(fscanf(BaryonFrigaFileAnti, "%*[^\n]%*c")==EOF)
	    throw runtime_error("Unexpected end of file fort.881 at run " + to_string(irun));
	  if(fscanf(BaryonFrigaFileAnti, "%*[^\n]%*c")==EOF)
	    throw runtime_error("Unexpected end of file fort.881 at run " + to_string(irun));
	  if(fscanf(BaryonFrigaFileAnti, "%i %*[^\n]%*c", &nBaryEntries)==EOF)
	    throw runtime_error("Unexpected end of file fort.881 at run " + to_string(irun));

	  for (int i = 0; i < nBaryEntries; i++) {
	    if(fscanf(BaryonFrigaFileAnti, "%*i %i %f %f %f %f %f %f %f %i %i %i %*i %*i %*i %*i %f %f\n", &charge, &Px, &Py, &Pz, &Xpos, &Ypos, &Zpos, &Mass, &clusterId, &nBary, &baryonId, &TimeFreeze, &Ebin)==EOF) {
	      throw runtime_error("Unexpected end of file fort.881 at run " + to_string(irun) + "timestep" + to_string(it) + " particle " + to_string(i));
	    }

	    GetPdgIdBaryon(charge, pdgId, IsAnti);
	    outputTmp->cd();
	    event->fbaryons.push_back(PBaryon(pdgId, nBary, Px, Py, Pz, Xpos, Ypos, Zpos, Mass, baryonId, -clusterId, TimeFreeze, Ebin));
	  } 
	} 
	treeTmp->Fill();
      }  // end loop parallel runs
    }  // end loop timesteps
  }  // end loop subsequent runs

  int check_eof;
  if(fscanf(BulkFile, "%i %*[^\n]%*c", &check_eof) != EOF) 
    throw runtime_error("\n  Error when reading " + inputFileBulk + ": File not read until the end. Check input format.\n ");
  
  outputTmp->cd();
  treeTmp->Write();

  fclose(BulkFile);
  fclose(BaryonFrigaFile);
  if (ConvertAnti == kTRUE) fclose(BaryonFrigaFileAnti);

  //********Make Maps to find baryons across timesteps & in phsd.dat and to assign baryons to clusters *****************

  if (NTIME > 29) throw runtime_error("Size of eventmap is too small for " + to_string(NTIME) + " timesteps in fort.891");
  baryonId2pos.resize(treeTmp->GetEntries());
  baryons2hadrons.resize(treeTmp->GetEntries());
  clusterId2baryonIds.resize(treeTmp->GetEntries());
  eventId2time.clear();
 
  outputTmp->cd();
  for (int ievent=0;ievent<treeTmp->GetEntries();ievent++) {

    outputTmp->cd();
    treeTmp->GetEntry(ievent);
 
    eventId2Entry[event->fEventId][event->fStepNr] = ievent;

    auto it_time = eventId2time.find(event->fEventId);
    if (it_time != eventId2time.end())
      eventId2time[event->fEventId].push_back(event->fTime);
    else 
      eventId2time[event->fEventId] = {event->fTime};
    
    for (int ibaryon = 0 ; ibaryon < event->fbaryons.size() ; ibaryon++) {
      PBaryon	baryon = event->fbaryons[ibaryon];
      baryonId2pos[ievent][baryon.fBaryonId] = ibaryon;
      Int_t clusterId = baryon.fClusterId;
      Int_t baryonId = baryon.fBaryonId;
      auto it = clusterId2baryonIds[ievent].find(clusterId);
      if (it != clusterId2baryonIds[ievent].end()) 
	it->second.emplace_back(baryonId);
      else 
	clusterId2baryonIds[ievent][clusterId] = {static_cast<int>(baryonId)};

      for (int ihadron = 0 ; ihadron < event->fhadrons.size() ; ihadron++) {
	PHadron hadron = event->fhadrons[ihadron];
	if (hadron.fBaryonId == baryon.fBaryonId) {
	  baryons2hadrons[ievent][baryon.fBaryonId] = ihadron;
	  break;
	}
      }
    }
  }

  //********Fill final trees *****************
  
  URun *header; UEvent *uevent; TFile *output; TTree *tree;
 
  if (WriteUnigen == kTRUE) {
    if (ConvertAllClusters == kTRUE)
      header = new URun ("phqmd", "all clusters", aProj, zProj, pProj, aTarg, zTarg, pTarg, bMin, bMax, 0, 0, 0, 0, nEvents);
    else
      header = new URun ("phqmd", "physical clusters A < 10", aProj, zProj, pProj, aTarg, zTarg, pTarg, bMin, bMax, 0, 0, 0, 0, nEvents);
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
  
  std::vector<ClusterEntry> clusterlist;
  GetClusterList(clustertable, clusterlist);
  gRandom->SetSeed(0);
  
  if (ConvertAllClusters == kTRUE)
    cout << "Conversion of physical clusters according to clustertable.\nConversion of clusters A > 7 independent of their physical existence." << endl;
  else
    cout << "Conversion of physical clusters according to clustertable.\nBaryons from all other clusters are counted as single baryons." << endl;
  outputTmp->cd();
  for (int ievent = 0 ; ievent < treeTmp->GetEntries(); ievent++) {

    outputTmp->cd();
    treeTmp->GetEntry(ievent);
   
    if(event->fStepNr == NTIME) {

      Int_t index = 0 ;
      Int_t parentId = -1;
      Int_t child[2] = {0,0};
      
      if (WriteUnigen == kTRUE) {
	output->cd();
	uevent->Clear();
	uevent->SetParameters(event->fEventId, event->fB, event->fPhi, event->fNParticipants, 1, event->fTime, 0);
      }
      
      if (WriteEventFreeze == kTRUE) {
	outputFreeze->cd();
	eventFreeze->Clear();
	eventFreeze->SetParameters(event->fEventId, event->fB, event->fNParticipants, event->fTime, event->fPhi);
      }
      
      outputTmp->cd(); // Loop over hadrons
      for (auto hadron : event->fhadrons) { 	
	auto it_bar2had = baryons2hadrons[ievent].find(hadron.fBaryonId);	
	if (it_bar2had != baryons2hadrons[ievent].end()) continue; // baryon is part of a cluster

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
      
      outputTmp->cd(); // Loop over cluster baryons & check if clusters are physical
      for (auto it_clId : clusterId2baryonIds[ievent]) {

	Int_t clusterId = it_clId.first;
	Int_t nbary = it_clId.second.size();

	std::vector<PBaryon_cluster> baryons_cluster;
	baryons_cluster.clear();
	
	for (Int_t baryonId : it_clId.second) {
	  auto it_bary = baryonId2pos[ievent].find(baryonId);
	  PBaryon baryon = event->fbaryons[it_bary->second];

	  auto it_had = baryons2hadrons[ievent].find(baryon.fBaryonId);
	  Int_t hadronId = it_had->second;
	  PHadron hadron = event->fhadrons[hadronId];
	  
	  baryons_cluster.push_back(PBaryon_cluster(baryonId, baryon.fPdgId, hadron.fProcessId, hadron.fParentId, baryon.fP, baryon.fMass, baryon.fEnergy, hadron.fXTFreeze, hadron.fPEFreeze, baryon.fEbin));	    
	}
	
	if (nbary == 1) {
	  Int_t ibary  = 0;
	  if (WriteUnigen == kTRUE) {
	    output->cd();
	    uevent->AddParticle (index, baryons_cluster.at(ibary).fPdgId, 1, baryons_cluster.at(ibary).fParentId, -1, -1, baryons_cluster.at(ibary).fProcessId, child, baryons_cluster.at(ibary).fP.X(), baryons_cluster.at(ibary).fP.Y(), baryons_cluster.at(ibary).fP.Z(), baryons_cluster.at(ibary).fEnergy, baryons_cluster.at(ibary).fXTFreeze.X(), baryons_cluster.at(ibary).fXTFreeze.Y(), baryons_cluster.at(ibary).fXTFreeze.Z(), baryons_cluster.at(ibary).fXTFreeze.T(), 1);
	  }	  
	  if (WriteEventFreeze == kTRUE) {
	    outputFreeze->cd();
	    eventFreeze->AddParticle(index, baryons_cluster.at(ibary).fPdgId, baryons_cluster.at(ibary).fParentId, baryons_cluster.at(ibary).fProcessId, baryons_cluster.at(ibary).fP, baryons_cluster.at(ibary).fEnergy, baryons_cluster.at(ibary).fXTFreeze, baryons_cluster.at(ibary).fPEFreeze, 1, 1);
	  }	  
	  index++;
	}
	
	if (nbary > 1) {
	  Int_t pdgId = 99999;
	  GetClusterPdg(clusterlist, baryons_cluster, clusterId, ConvertAllClusters, pdgId);
	  Float_t Ebin = CalculateClusterBindingEnergy(baryons_cluster);
	  
	  if (pdgId != 99999 && Ebin/nbary < Ebin_max) {
	    Float_t Px = 0; Float_t Py = 0; Float_t Pz = 0; Float_t energy = 0;
	    CalculateClusterKin(baryons_cluster, Px, Py, Pz, energy);

	    Int_t clusterProcessId = -1; Int_t clusterParentId = -1;
	    GetClusterProcessParentId(baryons_cluster, clusterProcessId, clusterParentId);
	    
	    //Calculate Cluster freeze-out coordinates
	    Float_t TimeProductionCluster = 0.0;
	    CalculateClusterProductionTime(event->fEventId, NTIME, clusterId, nbary, TimeProductionCluster);
	    Float_t TimeFreezeCluster = 0.0 ;Float_t deltaT = 0.0; Int_t TsFreeze = 0;
	    CalculateClusterFreezeOutTime(event->fEventId, baryons_cluster, nbary, TimeProductionCluster, TsFreeze, TimeFreezeCluster, deltaT);
 
	    auto it_eventnr = eventId2Entry.find(event->fEventId);
	    Int_t eventNrFreeze = it_eventnr -> second[TsFreeze];
	    treeTmp->GetEntry(eventNrFreeze);

	    TVector3 posfo_cluster = {0.0, 0.0, 0.0}; TVector3 pfo_cluster = {0.0, 0.0, 0.0}; Float_t mass_cluster = 0;
	    for (int ibary = 0; ibary < nbary; ibary++) {	      
		auto it_bary = baryonId2pos[eventNrFreeze].find(baryons_cluster.at(ibary).fBaryonId);
		PBaryon baryon_fo = event->fbaryons[it_bary->second];
		for (int i = 0; i < 3; i++) {
		  if (baryons_cluster.at(ibary).fXTFreeze.T() < TimeFreezeCluster)
		    posfo_cluster(i) += baryon_fo.fX(i) + deltaT * baryon_fo.fP(i) / baryon_fo.fEnergy;
		  else
		    posfo_cluster(i) +=  baryons_cluster.at(ibary).fXTFreeze(i);
		  pfo_cluster(i) +=  baryon_fo.fP(i);
		}
		mass_cluster += baryon_fo.fMass;
	    }
	    for (int i = 0; i < 3; i++)
	      posfo_cluster(i) /= nbary;

	    Float_t energyfo_cluster = TMath::Sqrt(mass_cluster*mass_cluster + pfo_cluster(0)* pfo_cluster(0) + pfo_cluster(1)* pfo_cluster(1) + pfo_cluster(2)* pfo_cluster(2));

	    treeTmp->GetEntry(ievent);

	    Int_t weight;
	    
	    if (WriteUnigen == kTRUE) {
	      output->cd();
	      weight = 1;
	      uevent->AddParticle (index, pdgId, 1, clusterParentId, -1, -1, clusterProcessId, child, Px, Py, Pz, energy, posfo_cluster.X(), posfo_cluster.Y(), posfo_cluster.Z(), TimeFreezeCluster, weight);
	      weight = 0;
	      for (int ibary=0;ibary<nbary;ibary++) // store single baryons with weight = 0
	      	uevent->AddParticle (index, baryons_cluster.at(ibary).fPdgId, 1, baryons_cluster.at(ibary).fParentId, -1, -1, baryons_cluster.at(ibary).fProcessId, child, baryons_cluster.at(ibary).fP.X(), baryons_cluster.at(ibary).fP.Y(), baryons_cluster.at(ibary).fP.Z(), baryons_cluster.at(ibary).fEnergy, baryons_cluster.at(ibary).fXTFreeze.X(), baryons_cluster.at(ibary).fXTFreeze.Y(), baryons_cluster.at(ibary).fXTFreeze.Z(), baryons_cluster.at(ibary).fXTFreeze.T(), weight);
 
	    }	    
	    if (WriteEventFreeze == kTRUE) {
	      outputFreeze->cd();
	      weight = 1;
	      eventFreeze->AddParticle(index, pdgId, clusterParentId, clusterProcessId, Px, Py, Pz, energy, TimeFreezeCluster, posfo_cluster, pfo_cluster, energyfo_cluster, 1, weight);
	      weight = 0;
	      for (int ibary=0;ibary<nbary;ibary++) // store single baryons with weight = 0
		eventFreeze->AddParticle(index, baryons_cluster.at(ibary).fPdgId, baryons_cluster.at(ibary).fParentId, baryons_cluster.at(ibary).fProcessId, baryons_cluster.at(ibary).fP, baryons_cluster.at(ibary).fEnergy, baryons_cluster.at(ibary).fXTFreeze, baryons_cluster.at(ibary).fPEFreeze, 1, weight);
	      
	    }	    
	    index++;
	  }
	  else {
	    for (int ibary=0;ibary<nbary;ibary++) {
	      if (WriteUnigen == kTRUE) {
		output->cd();
		uevent->AddParticle (index, baryons_cluster.at(ibary).fPdgId, 1, baryons_cluster.at(ibary).fParentId, -1, -1, baryons_cluster.at(ibary).fProcessId, child, baryons_cluster.at(ibary).fP.X(), baryons_cluster.at(ibary).fP.Y(), baryons_cluster.at(ibary).fP.Z(), baryons_cluster.at(ibary).fEnergy, baryons_cluster.at(ibary).fXTFreeze.X(), baryons_cluster.at(ibary).fXTFreeze.Y(), baryons_cluster.at(ibary).fXTFreeze.Z(), baryons_cluster.at(ibary).fXTFreeze.T(), 1);
	      }	      
	      if (WriteEventFreeze == kTRUE) {
		outputFreeze->cd();
		eventFreeze->AddParticle(index, baryons_cluster.at(ibary).fPdgId, baryons_cluster.at(ibary).fParentId, baryons_cluster.at(ibary).fProcessId, baryons_cluster.at(ibary).fP, baryons_cluster.at(ibary).fEnergy, baryons_cluster.at(ibary).fXTFreeze, baryons_cluster.at(ibary).fPEFreeze, 1, 1);
	      }	      
	      index++;
	    }
	  }
	}
      } 
      if (WriteUnigen == kTRUE) { 
        output->cd();
        tree->Fill();
      }
      if (WriteEventFreeze == kTRUE) {
        outputFreeze->cd();
        treeFreeze->Fill();
      }	
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
  
  
