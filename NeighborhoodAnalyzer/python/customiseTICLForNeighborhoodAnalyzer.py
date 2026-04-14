from PlayGround.NeighborhoodAnalyzer.NeighborhoodAnalyzer_cfi import *
#from PlayGround.PCaloHitAnalyzer.PCaloHitAnalyzer_cfi import *

def customiseTICLForNeighborhoodAnalyzer(process):
    process.NeighborhoodAnalyzer = NeighborhoodAnalyzer.clone()
    return process