//************************************************************************************************************************************************
/** author: Susanne Glaessel (Universitaet Frankfurt)
 ** Macro for creating the detector input from PHQMD output (PHQMD versions 5.2 
 ** & 5.2 Winn) 
 **
 ** This macro converts the PHQMD output-files after stabilization for the final 
 ** timestep to the UniGen format.
 ** 
 ** Options:
 **
 ** 1) <ConvertAllClusters>
 ** In PHQMD clusters / anticlusters are recognised independently of their physical 
 ** existence. The cluster-baryons are listed separately in the outputfile. This 
 ** routine builds clusters from the single baryons based on their cluster-ID and 
 ** identifies physical clusters according to the cluster_table.root. Baryons from 
 ** unphysical clusters (eg. p-p) are counted as single baryons.
 ** With option "ConvertAllClusters": Clusters with A > 7 are counted independent of 
 ** their physical existince.
 **
 ** 2) <ConvertAnti>
 ** PHQMD writes baryons and anti-baryons into two separate files. The conversion 
 ** of anti-baryons is optional and can be switched off with option 
 ** ConvertAnti = kFALSE.
 **
 ** 3) <FreezeCoords>
 ** In the PHQMD-code there is the option to write out freeze-out coordinates 
 ** (IFreezeOut=1). If the PHQMD-output contains freeze-out coordinates, option  
 ** FreezeCoords = kTRUE needs to be set. 
 ** (Note: Momentum in UniGen-ouput is the final momentum, not the momentum at 
 ** freeze-out.)
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
 ** UParticle* particle -> GetStatus(): = 0 kinetic deuteron, = 1 MST deuteron.
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
  Int_t fBaryonId;
  TLorentzVector fXTFreeze;
  TVector3 fPFreeze;
  PHadron() : fPdgId(0), fEnergy(0.), fBaryonId(0) { fP.SetXYZ(0.,0.,0.); fXTFreeze.SetXYZT(0.,0.,0.,0.); fPFreeze.SetXYZ(0.,0.,0.);  };
  PHadron(Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Int_t baryonId) : fPdgId(pdgId), fEnergy(energy), fBaryonId(baryonId) { fP.SetXYZ(0.0, 0.0, 0.0); fXTFreeze.SetXYZT(0.0, 0.0, 0.0, 0.0); fPFreeze.SetXYZ(0.0, 0.0, 0.0);  };
  PHadron(Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Int_t baryonId, Float_t xposfo, Float_t yposfo, Float_t zposfo , Float_t timefo, Float_t xpfo, Float_t ypfo, Float_t zpfo) : fPdgId(pdgId), fEnergy(energy), fBaryonId(baryonId) { fP.SetXYZ(Px,Py,Pz); fXTFreeze.SetXYZT(xposfo,yposfo,zposfo,timefo); fPFreeze.SetXYZ(xpfo,ypfo,zpfo);  };
};

struct PBaryon {
  Int_t fPdgId;
  Int_t fnBary;
  TVector3 fP;
  TVector3 fX;
  Float_t fEnergy;
  Int_t fBaryonId;
  Int_t fClusterId;
  Float_t fEbin;
  Float_t fTimeFreeze;
  PBaryon() : fPdgId(0),fnBary(0),fEnergy(0.),fBaryonId(0),fClusterId(0),fEbin(0), fTimeFreeze(0.) { fP.SetXYZ(0.,0.,0.); fX.SetXYZ(0.,0.,0.); };
  PBaryon(Int_t PdgId, Int_t nBary, Float_t Px, Float_t Py, Float_t Pz, Float_t Mass, Float_t X, Float_t Y, Float_t Z, Int_t baryonId, Int_t clusterId, Float_t TimeFreeze, Float_t Ebin) : fPdgId(PdgId), fnBary(nBary), fBaryonId(baryonId), fClusterId(clusterId), fTimeFreeze(TimeFreeze), fEbin(Ebin) { fP.SetXYZ(Px,Py,Pz); fX.SetXYZ(X,Y,Z); fEnergy = TMath::Sqrt(Mass*Mass +Px*Px + Py*Py + Pz*Pz); };
};

struct PBaryon_cluster {
  Int_t fPdgId;
  TVector3 fP;
  Float_t fEnergy;
  Float_t fEbin;
  Float_t fTimeFreeze;
  TLorentzVector fXTFreeze;
  TVector3 fPFreeze;
  Int_t fBaryonId;
  Int_t fClusterId;
  PBaryon_cluster(Int_t baryonId, Int_t PdgId, TVector3 P, Float_t energy, TLorentzVector XTFreeze, TVector3 PFreeze, Float_t Ebin) : fBaryonId(baryonId), fPdgId(PdgId), fP(P), fEnergy(energy), fXTFreeze(XTFreeze), fPFreeze(PFreeze), fEbin(Ebin) {};
};

class PEvent : public TObject  {
public:
  Int_t fEventId;
  Int_t fStepNr;
  Int_t fnParticipants;
  Float_t fb;
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
  TVector3 fP;
  Float_t fEnergy;
  Float_t fTimeFreeze;
  TVector3 fXFreeze;
  TVector3 fPFreeze;
  Int_t fOrigin;
  ParticleFreeze() : fIndex(0), fPdgId(0), fEnergy(0.), fTimeFreeze(0.), fOrigin(0) {fP.SetXYZ(0.,0.,0.); fXFreeze.SetXYZ(0.,0.,0.); fPFreeze.SetXYZ(0.,0.,0.); };
  ParticleFreeze(Int_t index, Int_t pdgId, TVector3 P, Float_t energy, TLorentzVector XFreeze, TVector3 PFreeze, Int_t Origin) : fIndex(index), fPdgId(pdgId), fP(P), fEnergy(energy), fPFreeze(PFreeze), fOrigin(Origin) {fTimeFreeze = XFreeze.T(); fXFreeze.SetXYZ(XFreeze.X(), XFreeze.Y(), XFreeze.Z()); };
  ParticleFreeze(Int_t index, Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Float_t TimeFreeze, TVector3 XFreeze, TVector3 PFreeze, Int_t Origin) : fIndex(index), fPdgId(pdgId), fEnergy(energy), fTimeFreeze(TimeFreeze), fXFreeze(XFreeze), fPFreeze(PFreeze), fOrigin(Origin) {fP.SetXYZ(Px, Py, Pz); };
  ParticleFreeze(Int_t index, Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, TLorentzVector XFreeze, TVector3 PFreeze, Int_t Origin) : fIndex(index), fPdgId(pdgId), fEnergy(energy), fPFreeze(PFreeze), fOrigin(Origin) {fP.SetXYZ(Px, Py, Pz); fTimeFreeze = XFreeze.T(); fXFreeze.SetXYZ(XFreeze.X(), XFreeze.Y(), XFreeze.Z()); };
  ClassDef(ParticleFreeze, 1);
};

class EventFreeze : public TObject  {
public:
  Int_t fEventId;
  Float_t fb;
  Float_t fTime;
  Float_t fPhi;
  Int_t fnParticipants;
  Int_t fNpa;
  vector<ParticleFreeze> fparticles;
  void SetParameters(Int_t eventId, Float_t b, Float_t time, Float_t phi, Int_t nParticipants) { fEventId = eventId; fb = b; fTime = time; fPhi = phi; fnParticipants = nParticipants;};
  void AddParticle(Int_t index, Int_t pdgId, TVector3 P, Float_t energy, TLorentzVector XTFreeze, TVector3 PFreeze, Int_t Origin) { fparticles.push_back(ParticleFreeze(index, pdgId, P, energy, XTFreeze, PFreeze, Origin)); fNpa += 1; }; 
  void AddParticle(Int_t index, Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, Float_t TimeFreezeCluster, TVector3 posfo_cluster, TVector3 pfo_cluster, Int_t Origin) { fparticles.push_back(ParticleFreeze(index, pdgId, Px, Py, Pz, energy, TimeFreezeCluster, posfo_cluster, pfo_cluster, Origin)); fNpa += 1; }; 
  void AddParticle(Int_t index, Int_t pdgId, Float_t Px, Float_t Py, Float_t Pz, Float_t energy, TLorentzVector XTFreeze, TVector3 PFreeze, Int_t Origin) {
    fparticles.push_back(ParticleFreeze(index, pdgId, Px, Py, Pz, energy, XTFreeze, PFreeze, Origin)); fNpa += 1; };
  void Clear() { fparticles.clear(); fNpa = 0; };
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
    if (TMath::Abs(pdgId) == 2212 || TMath::Abs(pdgId) == 2112) Mass += 0.938;
    if (TMath::Abs(pdgId) == 3122) Mass += 1.116;
    if (TMath::Abs(pdgId) == 3212) Mass += 1.193;
  }
  energy = TMath::Sqrt(Mass*Mass+Px*Px+Py*Py+Pz*Pz);
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

void convert_phqmd_detector_unigen(TString indir = "",
					     TString dataset = "00001",
					     Int_t firstevent = 0,
					     Bool_t FreezeCoords = kTRUE,
					     Bool_t ConvertAllClusters = kTRUE,
					     Bool_t ConvertAnti = kTRUE)
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

  TString rootFileTmp =  Form("%s/%s.phqmd_tmp.root",outdir.Data(),dataset.Data());
  TString rootFileDet;
  if (ConvertAnti == kTRUE) rootFileDet = Form("%s/%s.phqmd.root",outdir.Data(),dataset.Data());
  if (ConvertAnti == kFALSE) rootFileDet = Form("%s/%s.phqmd_noanti.root",outdir.Data(),dataset.Data());
  
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
	  
	Int_t nPart, pdgId, charge, baryonId;
	Float_t Px, Py, Pz, energy, xposfo, yposfo, zposfo, timefo, xpfo, ypfo, zpfo;

	//Get Hadrons from phsd.dat

	if (it == NTIME) {

	  if(fscanf(BulkFile, "%i %*i %*i %f %*i %*i %*i %*f %*f\n", &nPart, &event->fb)==EOF)
	    throw runtime_error("Unexpected end of file phsd.dat at run " + to_string(irun));	
	  if(fscanf(BulkFile, "%i %f %*[^\n]%*c", &event->fnParticipants, &event->fPhi)==EOF)   
	    throw runtime_error("Unexpected end of file phsd.dat at run " + to_string(irun));

	  outputTmp->cd();
	  event->fhadrons.clear();
	  
	  for (int i = 0; i < nPart; i++) {

	    if (FreezeCoords == kFALSE)
	      if(fscanf(BulkFile, "%i %i %f %f %f %f %*i %*i %i\n", &pdgId, &charge, &Px, &Py, &Pz, &energy, &baryonId)==EOF) {
		throw runtime_error("Unexpected end of file phsd.dat at run " + to_string(irun) + " particle " + to_string(i));
	      }

	    if (FreezeCoords == kTRUE)
	      if(fscanf(BulkFile, "%i %i %f %f %f %f %*i %*i %i %*f %*f %*f %f %f %f %f %f %f %f %*f %*f %*f\n", &pdgId, &charge, &Px, &Py, &Pz, &energy, &baryonId, &xposfo, &yposfo, &zposfo, &timefo, &xpfo, &ypfo, &zpfo)==EOF) {	
		throw runtime_error("Unexpected end of file phsd.dat at run " + to_string(irun) + " particle " + to_string(i));
	      }
	      
	    if((pdgId<1000 && pdgId>-1000) || TMath::Abs(pdgId) == 100121) baryonId = -1;
	    if(TMath::Abs(pdgId) == 100121) pdgId = 1000010020*charge; // correct pdg-code for kinetic deuterons    
	    outputTmp->cd();
	    if (FreezeCoords == kFALSE)
	      event->fhadrons.push_back(PHadron(pdgId, Px, Py, Pz, energy, baryonId));
	    if (FreezeCoords == kTRUE)
	      event->fhadrons.push_back(PHadron(pdgId, Px, Py, Pz, energy, baryonId, xposfo, yposfo, zposfo ,timefo, xpfo,ypfo,zpfo));
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
	  event->fbaryons.push_back(PBaryon(pdgId, nBary, Px, Py, Pz, Mass, Xpos, Ypos, Zpos, baryonId, clusterId, TimeFreeze, Ebin));
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
	    event->fbaryons.push_back(PBaryon(pdgId, nBary, Px, Py, Pz, Mass, Xpos, Ypos, Zpos, baryonId, -clusterId, TimeFreeze, Ebin));
	  } 
	} 
	treeTmp->Fill();
      }  // end loop parallel runs
    }  // end loop timesteps
  }  // end loop subsequent runs
  
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
 
  if (ConvertAllClusters == kTRUE)
    header = new URun ("phqmd", "all clusters", aProj, zProj, pProj, aTarg, zTarg, pTarg, bMin, bMax, 0, 0, 0, 0, nEvents);
  else
    header = new URun ("phqmd", "physical clusters A < 10", aProj, zProj, pProj, aTarg, zTarg, pTarg, bMin, bMax, 0, 0, 0, 0, nEvents);
  uevent = new UEvent;
  output = new TFile (rootFileDet, "recreate");
  tree = new TTree ("events", "signal");
  header->Write();
  tree->Branch ("event", "UEvent", uevent);
  
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

      Int_t index=0;
      Int_t child[2] = {0,0};
      
      output->cd();
      uevent->Clear();
      uevent->SetParameters(event->fEventId, event->fb, event->fPhi, event->fnParticipants, 1, event->fTime, 0);
      
      outputTmp->cd(); // Loop over hadrons
      for (auto hadron : event->fhadrons) { 	
	auto it_bar2had = baryons2hadrons[ievent].find(hadron.fBaryonId);	
	if (it_bar2had != baryons2hadrons[ievent].end()) continue; // baryon is part of a cluster

	output->cd();
	uevent->AddParticle (index, hadron.fPdgId, 0, -1, -1, -1, -1, child, hadron.fP.X(), hadron.fP.Y(),hadron.fP.Z(), hadron.fEnergy, hadron.fXTFreeze.X(), hadron.fXTFreeze.Y(), hadron.fXTFreeze.Z(), hadron.fXTFreeze.T(), 1);

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

	  baryons_cluster.push_back(PBaryon_cluster(baryonId, baryon.fPdgId, baryon.fP, baryon.fEnergy, hadron.fXTFreeze, hadron.fPFreeze, baryon.fEbin));	    
	}
	
	if (nbary == 1) {
	  Int_t ibary  = 0;
	  output->cd();
	  uevent->AddParticle (index, baryons_cluster.at(ibary).fPdgId, 1, -1, -1, -1, -1, child, baryons_cluster.at(ibary).fP.X(), baryons_cluster.at(ibary).fP.Y(), baryons_cluster.at(ibary).fP.Z(), baryons_cluster.at(ibary).fEnergy, baryons_cluster.at(ibary).fXTFreeze.X(), baryons_cluster.at(ibary).fXTFreeze.Y(), baryons_cluster.at(ibary).fXTFreeze.Z(), baryons_cluster.at(ibary).fXTFreeze.T(), 1); 	  
	  index++;
	}
	
	if (nbary > 1) {
	  Int_t pdgId = 99999;
	  GetClusterPdg(clusterlist, baryons_cluster, clusterId, ConvertAllClusters, pdgId);
	  Float_t Ebin = CalculateClusterBindingEnergy(baryons_cluster);
	  
	  if (pdgId != 99999 && Ebin/nbary < Ebin_max) {	   
	    Float_t Px = 0; Float_t Py = 0; Float_t Pz = 0; Float_t energy = 0; 
	    CalculateClusterKin(baryons_cluster, Px, Py, Pz, energy);
	    
	    //Calculate Cluster freeze-out coordinates
	    
	    Float_t TimeFreezeCluster = 0.0 ;
	    TVector3 posfo_cluster = {0.0, 0.0, 0.0}; TVector3 pfo_cluster = {0.0, 0.0, 0.0};
	    
	    if (FreezeCoords == kTRUE) {
	      Float_t TimeProductionCluster = 0.0;
	      CalculateClusterProductionTime(event->fEventId, NTIME, clusterId, nbary, TimeProductionCluster);
	      Float_t deltaT = 0.0; Int_t TsFreeze = 0;
	      CalculateClusterFreezeOutTime(event->fEventId, baryons_cluster, nbary, TimeProductionCluster, TsFreeze, TimeFreezeCluster, deltaT);
	    
	      auto it_eventnr = eventId2Entry.find(event->fEventId);
	      Int_t eventNrFreeze = it_eventnr -> second[TsFreeze];
	      treeTmp->GetEntry(eventNrFreeze);
	    
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
	      }
	      for (int i = 0; i < 3; i++)
		posfo_cluster(i) /= nbary;

	      treeTmp->GetEntry(ievent);
	    }
	    
	    output->cd(); 
	    uevent->AddParticle (index, pdgId, 1, -1, -1, -1, -1, child, Px, Py, Pz, energy, posfo_cluster.X(), posfo_cluster.Y(), posfo_cluster.Z(), TimeFreezeCluster, 1);
	   
	    index++;
	  }
	  else {
	    for (int ibary=0;ibary<nbary;ibary++) {
	      output->cd();
	      uevent->AddParticle (index, baryons_cluster.at(ibary).fPdgId, 1, -1, -1, -1, -1, child, baryons_cluster.at(ibary).fP.X(), baryons_cluster.at(ibary).fP.Y(), baryons_cluster.at(ibary).fP.Z(), baryons_cluster.at(ibary).fEnergy, baryons_cluster.at(ibary).fXTFreeze.X(), baryons_cluster.at(ibary).fXTFreeze.Y(), baryons_cluster.at(ibary).fXTFreeze.Z(), baryons_cluster.at(ibary).fXTFreeze.T(), 1);    
	      index++;
	    }
	  }
	}
      } 
 
      output->cd();
      tree->Fill();
   
    } 
  }
  outputTmp->Close();
  gSystem->Unlink(rootFileTmp);

  output->cd();
  tree->Write();
  output->Close();
  
  cout << endl;
  cout << "Macro finished successfully." << endl;
  cout << "Output file is: " <<rootFileDet<< endl;
}
  
  
