///////////////////////
//// Weights demo ////
///////////////////////
#include "Hammer/Tools/HammerRoot.hh"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TBranch.h"
#include "Hammer/Hammer.hh"
#include "Hammer/Process.hh"
#include "Hammer/Particle.hh"
#include "Hammer/Math/FourMomentum.hh"
#include <TLorentzVector.h>

using namespace std;
using namespace Hammer;

int main() {

    ////////////////////////////////////////////////////////
    //// Initialize the sample and store tensor weights ////
    ////////////////////////////////////////////////////////

    Double_t Bc_TRUEP_X, Bc_TRUEP_Y, Bc_TRUEP_Z, Bc_TRUEP_E;
    Double_t Jpsi_TRUEP_X, Jpsi_TRUEP_Y, Jpsi_TRUEP_Z, Jpsi_TRUEP_E;
    Double_t BachMu_TRUEP_X, BachMu_TRUEP_Y, BachMu_TRUEP_Z, BachMu_TRUEP_E;
    Float_t Bc_Nu_TRUEP_X, Bc_Nu_TRUEP_Y, Bc_Nu_TRUEP_Z, Bc_Nu_TRUEP_E;
    Double_t FitVar_El, FitVar_q2, FitVar_Mmiss2;
    Int_t Bc_TRUEID, Jpsi_TRUEID, BachMu_TRUEID, Bc_Nu_TRUEID;
    UInt_t runNumber;
    ULong64_t eventNumber;

    std::cout << "Reading input file" << std::endl;
    TFile* inputfile{new TFile{"Bc2JpsiMuNu_subset.root", "READ"}};
    std::cout << "Input file read" << std::endl;
    std::cout << "Reading input tree" << std::endl;
    TTree* inputtree = (TTree*) inputfile->Get("DecayTree");
    std::cout << "Tree read" << std::endl;

    Int_t inputEvents = inputtree->GetEntries();
    inputtree->SetBranchAddress("Bc_TRUEP_X",&Bc_TRUEP_X);
    inputtree->SetBranchAddress("Bc_TRUEP_Y",&Bc_TRUEP_Y);
    inputtree->SetBranchAddress("Bc_TRUEP_Z",&Bc_TRUEP_Z);
    inputtree->SetBranchAddress("Bc_TRUEP_E",&Bc_TRUEP_E);
    inputtree->SetBranchAddress("Bc_TRUEID",&Bc_TRUEID);
    inputtree->SetBranchAddress("Jpsi_TRUEP_X",&Jpsi_TRUEP_X);
    inputtree->SetBranchAddress("Jpsi_TRUEP_Y",&Jpsi_TRUEP_Y);
    inputtree->SetBranchAddress("Jpsi_TRUEP_Z",&Jpsi_TRUEP_Z);
    inputtree->SetBranchAddress("Jpsi_TRUEP_E",&Jpsi_TRUEP_E);
    inputtree->SetBranchAddress("Jpsi_TRUEID",&Jpsi_TRUEID);
    inputtree->SetBranchAddress("BachMu_TRUEP_X",&BachMu_TRUEP_X);
    inputtree->SetBranchAddress("BachMu_TRUEP_Y",&BachMu_TRUEP_Y);
    inputtree->SetBranchAddress("BachMu_TRUEP_Z",&BachMu_TRUEP_Z);
    inputtree->SetBranchAddress("BachMu_TRUEP_E",&BachMu_TRUEP_E);
    inputtree->SetBranchAddress("BachMu_TRUEID",&BachMu_TRUEID);
    inputtree->SetBranchAddress("Bc_Nu_TRUEP_X",&Bc_Nu_TRUEP_X);
    inputtree->SetBranchAddress("Bc_Nu_TRUEP_Y",&Bc_Nu_TRUEP_Y);
    inputtree->SetBranchAddress("Bc_Nu_TRUEP_Z",&Bc_Nu_TRUEP_Z);
    inputtree->SetBranchAddress("Bc_Nu_TRUEP_E",&Bc_Nu_TRUEP_E);
    inputtree->SetBranchAddress("Bc_Nu_TRUEID",&Bc_Nu_TRUEID);

    inputtree->SetBranchAddress("FitVar_El", &FitVar_El);
    inputtree->SetBranchAddress("FitVar_Mmiss2", &FitVar_Mmiss2);
    inputtree->SetBranchAddress("FitVar_q2", &FitVar_q2);

    inputtree->SetBranchAddress("runNumber", &runNumber);
    inputtree->SetBranchAddress("eventNumber", &eventNumber);

    //Outfile
    TFile* outFile{new TFile{"./DemoWeights.root", "RECREATE"}};
    // prepares the tree
    TTree* tree{new TTree{"Hammer", "Hammer output Tree"}};
    tree->SetDirectory(outFile);
    //gInterpreter->GenerateDictionary("std::set<size_t>","set");
    Int_t evNumber = 0;
    //set<size_t> processIds;
    Hammer::RootIOBuffer treeBuffer;
    Int_t processId = 0;
    treeBuffer.maxLength = 256*1024*1024;
    treeBuffer.length = 0;
    treeBuffer.start = new UChar_t[static_cast<size_t>(treeBuffer.maxLength)];
    tree->Branch("ev_number", &evNumber, "ev_number/I");
    tree->Branch("runNumber", &runNumber, "runNumber/I");
    tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    tree->Branch("proc_ids", &processId);
    tree->Branch("type", &treeBuffer.kind, "type/B");
    tree->Branch("length", &treeBuffer.length, "length/I");
    tree->Branch("record", treeBuffer.start, "record[length]/b");

    Hammer::Hammer ham{};
    Hammer::IOBuffer fbBuffer;
    //Declare included processes
    ham.includeDecay(string("BcJpsiMuNu"));
    ham.addFFScheme("Scheme1", {{"BcJpsi", "BGL"}});
    //ham.setOptions("BctoJpsiBGL: {Chim: 0.007168, Chip: 0.012539}"); //Set options for FF parameters
    //ham.setOptions("BctoJpsiBGL: {avec: [0.003,-0.022,0.150] }"); //Set options for FF parameters
    ham.setOptions("BctoJpsiBGL: {dvec: [0.,0.,0.] }"); //Set options for FF parameters
    //Can add more schemes eg: ham.addFFScheme("Scheme2", {{"BD", "CLN"}, {"BD*", "BGL"}});
    //Declare the input FF scheme
    ham.setFFInputScheme({{"BcJpsi", "EFG"}});

    ham.addHistogram("q2_vs_El_vs_Mmiss", {4, 20, 18}, true,
                     {{0., 11800000.}, {0., 2600.}, {-2000000., 10000000.}});
    ham.addHistogram("q2",{4}, true, {{0., 11800000.}});
    ham.addHistogram("q2_true",{4}, true, {{0., 11800000.}});
    ham.addHistogram("El", {20}, true, {{0., 2600.}});
    ham.addHistogram("El_toFit", {20}, true, {{0., 2600.}});
    ham.addHistogram("Mmiss", {18}, true, {{-2000000., 10000000.}});

    ham.setOptions("Histos: {KeepErrors: false}"); //This line will enable the errors in the bin histograms
    ham.setUnits("MeV");
    ham.initRun();
    //ham.setOptions("BtoDCLN_num: {RhoSq: 2.}");
    //ham.setOptions("BtoDCLN: {RhoSq: 1.185}");

    //ham.addPurePSVertices({"BDMuNu"},false); //false meaning denominator, true num, default=true



    //Saves the FF scheme definitions & other run information
    fbBuffer = ham.saveRunHeader();
    treeBuffer = fbBuffer; // necessary to make sure buffer starts always at same location (alternatives?)
    tree->Fill();
    //Loop over events. Can be more than one process per event if there are two taus, or zero!
    size_t savedEvents = 0ul;

    Double_t Mcorr;

    /*
    TBranch* B0_PX_branch = (TBranch*) inputtree->GetBranch();
    TBranch* B0_PY_branch = (TBranch*) inputtree->GetBranch();

    B0_PX_branch->SetBranchAddress();
    B0_PY_branch->SetBranchAddress();
    */

    TLorentzVector vmu, vnu;
    for(Int_t i =0; i < inputEvents; ++i) {
        //if (io->rdstate() == 0 && io->fill_next_event(&ge)) {

            inputtree->GetEntry(i);

            evNumber = static_cast<Int_t>(i);
            processId = 0;
            if (i % 1000 == 0) {
                cout<< "processing event " << i << endl;
            }
            Particle Bc({Bc_TRUEP_E,Bc_TRUEP_X,Bc_TRUEP_Y,Bc_TRUEP_Z},Bc_TRUEID);
            Particle Jpsi({Jpsi_TRUEP_E,Jpsi_TRUEP_X,Jpsi_TRUEP_Y,Jpsi_TRUEP_Z},Jpsi_TRUEID);
            Particle mu({BachMu_TRUEP_E,BachMu_TRUEP_X,BachMu_TRUEP_Y,BachMu_TRUEP_Z},BachMu_TRUEID);
            Particle neutrino({Bc_Nu_TRUEP_E,Bc_Nu_TRUEP_X,Bc_Nu_TRUEP_Y,Bc_Nu_TRUEP_Z},Bc_Nu_TRUEID);

            vmu.SetPxPyPzE(BachMu_TRUEP_X,BachMu_TRUEP_Y,BachMu_TRUEP_Z,BachMu_TRUEP_E);
            vnu.SetPxPyPzE(Bc_Nu_TRUEP_X,Bc_Nu_TRUEP_Y,Bc_Nu_TRUEP_Z,Bc_Nu_TRUEP_E);

            Process proc;

            auto Bc_idx = proc.addParticle(Bc);
            auto Jpsi_idx = proc.addParticle(Jpsi);
            auto mu_idx = proc.addParticle(mu);
            auto neutrino_idx = proc.addParticle(neutrino);

            proc.addVertex(Bc_idx,{Jpsi_idx,mu_idx,neutrino_idx});

            ham.initEvent(); //Here goes the weight
            auto procId = ham.addProcess(proc);
            // look for decay processes starting with B mesons
            // can also use Hammer PDG constants for increased readability
            // auto processes = parseGenEvent(ge, {PID::BPLUS, PID::BMINUS, PID::BZERO, -PID::BZERO});
                                             // that can be used to retrieve process-specific information later
            if (procId != 0){ // Make sure process isn't forbidden and add it to the procIds list
               processId = procId;

               //Compute the bin
               //size_t Mcorr_bin = binNumber(Mcorr);

               //ham.setEventHistogramBin("Mcorr", {Mcorr_bin});

               ham.fillEventHistogram("q2_vs_El_vs_Mmiss", {FitVar_q2, FitVar_El, FitVar_Mmiss2});
               ham.fillEventHistogram("q2", {FitVar_q2});
               ham.fillEventHistogram("q2_true", {(vmu+vnu).M2()});
               ham.fillEventHistogram("El", {FitVar_El});
               ham.fillEventHistogram("El_toFit", {FitVar_El});
               ham.fillEventHistogram("Mmiss", {FitVar_Mmiss2});

               //std::cout << FitVar_q2 <<" " << q2_bin <<  std::endl;

               ham.processEvent();
               fbBuffer = ham.saveEventWeights();
               treeBuffer = fbBuffer;
               tree->Fill();
               savedEvents +=1;
	       //std::cout << runNumber << "," << eventNumber << "," << evNumber <<std::endl;
            }
    }
    // save the total rates (optional)
    processId = 0;
    evNumber = 0;
    fbBuffer = ham.saveRates();
    treeBuffer = fbBuffer;
    tree->Fill();
//    // save the sum of weights of processed events (the ones for which processEvent() has been called)
//    fbBuffer = ham.saveHistogram("Total Sum of Weights");
//    treeBuffer = fbBuffer;
//    tree->Fill();
//
//    //fbBuffer = ham.saveHistogram("Mcorr");
//    fbBuffer = ham.saveHistogram("q2_vs_El_vs_Mmiss");
//    treeBuffer = fbBuffer;
//    tree->Fill();
//
//    fbBuffer = ham.saveHistogram("Mmiss");
//    treeBuffer = fbBuffer;
//    tree->Fill();
//
//    fbBuffer = ham.saveHistogram("El");
//    treeBuffer = fbBuffer;
//    tree->Fill();
//
//    fbBuffer = ham.saveHistogram("El_toFit");
//    treeBuffer = fbBuffer;
//    tree->Fill();
//
//    fbBuffer = ham.saveHistogram("q2");
//    treeBuffer = fbBuffer;
//    tree->Fill();
//
//    fbBuffer = ham.saveHistogram("q2_true");
//    treeBuffer = fbBuffer;
//    tree->Fill();

    // save the histograms
//    for(auto& elem: ham.saveHistogram("Total Sum of Weights")) {
//        treeBuffer = elem;
//        tree->Fill();
//    }
    for(auto& elem: ham.saveHistogram("q2_true")) {
        treeBuffer = elem;
        tree->Fill();
    }


    // write and close file
    tree->Write();
    outFile->Write();
    outFile->Close();

    cout << "Events Stored: " << savedEvents << endl;

    // perform cleanup
    delete outFile;
    delete[] treeBuffer.start;

} // int main()
