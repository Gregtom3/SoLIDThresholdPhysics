import ROOT
import uproot
import glob
import os 
import numpy as np

def get_lumi_from_integrated(integrated_luminosity,TARGET_TYPE,days):
    luminosity = integrated_luminosity / ((0.5 if TARGET_TYPE=="d" else 1) * 1.0e-24 * 1.0e-9 * 3600.0 * 24 * days)
    formatted_luminosity = f"{luminosity:.2e}".replace("+","")
    return formatted_luminosity
    
def get_beamE_lumi_events_days(PROJECT_NAME,PROJECT_DIR,version="bh"):
    
    if version not in ["bh","electroproduction","photoproduction"]:
        raise ValueError("'version' for 'get_beamE_lumi_events_days()' must be either 'bh' or 'electroproduction' or 'photoproduction'")
        
    sample_file = os.listdir(f"{PROJECT_DIR}/{PROJECT_NAME}/{version}/data/")[0]
    sample_tfile = ROOT.TFile(f"{PROJECT_DIR}/{PROJECT_NAME}/{version}/data/{sample_file}")
    beam_energy = np.round(sample_tfile.Get("beam_energy").GetBinContent(1),2)
    lumi = sample_tfile.Get("integrated_luminosity").GetBinContent(1)
    events = sample_tfile.Get("num_events").GetBinContent(1) * len(glob.glob(f"{PROJECT_DIR}/{PROJECT_NAME}/{version}/data/*.root"))
    days = sample_tfile.Get("days").GetBinContent(1)
    
    return beam_energy, lumi, events, days

def get_tchains(PROJECT_NAME, PROJECT_DIR, version="acc"):
    
    if version not in ["acc", "all"]:
        raise ValueError("'version' for 'get_tchains()' must be either 'acc' or 'all'")
    
    bh_chain = ROOT.TChain("tree"+("Acc" if version=="acc" else ""))
    photo_chain = ROOT.TChain("tree"+("Acc" if version=="acc" else ""))
    electro_chain = ROOT.TChain("tree"+("Acc" if version=="acc" else ""))
    
    bh_chain.Add(f"{PROJECT_DIR}/{PROJECT_NAME}/bh/data/*.root")
    photo_chain.Add(f"{PROJECT_DIR}/{PROJECT_NAME}/photoproduction/data/*.root")
    electro_chain.Add(f"{PROJECT_DIR}/{PROJECT_NAME}/electroproduction/data/*.root")
    
    return bh_chain, photo_chain, electro_chain

def get_uproots(PROJECT_NAME, PROJECT_DIR, version="acc"):
    
    if version not in ["acc", "all"]:
        raise ValueError("'version' for 'get_uproot()' must be either 'acc' or 'all'")
    
    tree_name = "tree" + ("Acc" if version == "acc" else "")

    bh_files = glob.glob(f"{PROJECT_DIR}/{PROJECT_NAME}/bh/data/*.root")
    photo_files = glob.glob(f"{PROJECT_DIR}/{PROJECT_NAME}/photoproduction/data/*.root")
    electro_files = glob.glob(f"{PROJECT_DIR}/{PROJECT_NAME}/electroproduction/data/*.root")

    bh_trees = [uproot.open(f)[tree_name] for f in bh_files]
    photo_trees = [uproot.open(f)[tree_name] for f in photo_files]
    electro_trees = [uproot.open(f)[tree_name] for f in electro_files]

    return bh_trees, photo_trees, electro_trees