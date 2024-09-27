from PlayGround.IsolationRequirementAnalyzer.IsolationRequirementAnalyzer_cfi import *

def customiseTICLForKalmanFilterIsolationRequirements(process):
    process.kfIsolationRequirements = isolationRequirementAnalyzer.clone()
    return process
