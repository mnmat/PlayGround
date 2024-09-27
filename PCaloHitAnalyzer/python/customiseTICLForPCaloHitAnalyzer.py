from PlayGround.PCaloHitAnalyzer.PCaloHitAnalyzer_cfi import *

def customiseTICLForPCaloHitAnalyzer(process):
    process.pcalohitAnalyzer = pcalohitAnalyzer.clone()
    return process
