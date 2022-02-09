///////////////////////
//// Weights demo ////
///////////////////////
#include "Hammer/Tools/HammerRoot.hh"

#include "TFile.h"
#include "TTree.h"
#include "Hammer/Hammer.hh"

#include "TSystem.h"
// #include "TInterpreter.h"

using namespace std;

int main() {

    /////////////////////////////////////////////////////
    //// Now reread event weights and play with them ////
    /////////////////////////////////////////////////////

    Hammer::Hammer ham{};

    // gInterpreter->GenerateDictionary("std::set<size_t>","set");
    TFile* inFile{new TFile{"./DemoWeights.root", "READ"}};
    // prepares the tree
    if (!inFile) { return EXIT_FAILURE; }
    TTree* tree;
    inFile->GetObject("Hammer",tree);
    tree->Print();
    TBranch* brecord = nullptr;
    TBranch* btype = nullptr;

    // Make a clone and save weights there
    TFile* outFile{new TFile{"./HammerWeights.root", "RECREATE"}};
    TTree *out_tree = tree->CloneTree(); // Can change to (0) to not copy data
    Double_t weight;
    auto weight_branch = out_tree->Branch("HammerWeight", &weight, "HammerWeight/D");

    Hammer::IOBuffer buf{Hammer::RecordType::UNDEFINED, 0ul, new uint8_t[16*1024*1024]};
    tree->SetBranchAddress("record",buf.start, &brecord);
    tree->SetBranchAddress("type",&buf.kind, &btype);
    Long64_t nevents = tree->GetEntries() - 3;
    for (size_t idW = 0; idW <= 0; ++idW) {
        auto rwgtstart = std::chrono::system_clock::now();
        // Begin run
        double val = (static_cast<double>(idW))*0.2;
        ham.setUnits("MeV");
        ham.initRun();
        auto entry = tree->LoadTree(0);
        brecord->GetEntry(entry);
        btype->GetEntry(entry);
        if(!ham.loadRunHeader(buf)) {
            inFile->Close();
            delete inFile;
            return EXIT_FAILURE;
        }
        //Set the WCs
        ////ham.setWilsonCoefficients("BtoCMuNu", {{"S_qLlL", val*1i}, {"T_qLlL", val/4.}});
        ham.setWilsonCoefficients("BtoCMuNu", {{"SM", 1.}});
        //ham.setWilsonCoefficients("BtoCMuNu", {1., val*1i, 0., 0., 0., val/4., 0., 0., 0., 0., 0.});
        //Container for evtwgts
        vector<double> evtwgts;
        evtwgts.reserve(static_cast<size_t>(nevents));
        for(Int_t i = 1; i<= nevents; ++i) {
            entry = tree->LoadTree(i);
            brecord->GetEntry(entry);
            btype->GetEntry(entry);
            if (i % 1000 == 0) {
                cout << "." << flush;
            }
            ham.initEvent();
            ham.loadEventWeights(buf);
            double evtwgt = ham.getWeight("Scheme1");

	    //We could have also done instead, without using evtIds:
            //double evtwgt = 1.;
            //auto wgtmap = ham.getWeights("Scheme1"); //get all process weights in the event map<HashId, double>
            //for(auto& elem : wgtmap){
            //    evtwgt *= elem.second;
            //}
            //Or we could have requested the weight restricted to specific subset of processes in the event
            //by passing the list of their IDs E.g. for two processes (proc1, proc2) one would write:
            //double evtwgt = ham.getWeight("Scheme1", {proc1, proc2});
            //Store the computed weight in a vector<double>, or do whatever you want with it!
            evtwgts.push_back(evtwgt);
	    if (i<10) {cout <<"i = "<< i << "," <<"evtwgt = "<< evtwgt << endl;}
        }
        cout << endl << "Reweighted " << evtwgts.size() << " events to S_qLlL = " << val << "i, T_qLlL = " << val/4. << ": ";
        for(auto i =0; i < 10; i++) {
            cout << evtwgts[i] << ", ";
        }
        cout << " ...." << endl;
	// Update weight vector with correct order to fill tree
	// because the tree has n+3 entries: [header, [weights], ?, ?]
	evtwgts.insert( evtwgts.begin(), 0. );
	evtwgts.insert( evtwgts.end(), 2, 0. );
	//cout << "n entries: " << evtwgts.size() << " vs " << out_tree->GetEntries() << endl;
        // Save weights into the new tree
	for(int i =0; i < out_tree->GetEntries(); i++) {
	    weight = evtwgts[i];
	    weight_branch->Fill();
        }
        auto rwgtend = std::chrono::system_clock::now();
        auto durrwgt = rwgtend - rwgtstart;
        typedef std::chrono::duration<float> float_seconds;
        auto secsrwgt = std::chrono::duration_cast<float_seconds>(durrwgt);
        cout << "Reweight Time: " << secsrwgt.count() << endl;
    }
    tree->ResetBranchAddresses();
    inFile->Close();
    out_tree->Write("", TObject::kOverwrite);// overwrites existing tree so only 1 tree header saved
    outFile->Close();

    delete inFile;
    delete[] buf.start;
} // int main()
