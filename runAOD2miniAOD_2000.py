# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --filein file:B2G-RunIISpring15DR74-00001_step2.root --fileout file:B2G-RunIISpring15DR74-00001.root --mc --eventcontent MINIAODSIM --runUnscheduled --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1,Configuration/DataProcessing/Utils.addMonitoring --datatier MINIAODSIM --conditions MCRUN2_74_V9 --step PAT --python_filename /afs/cern.ch/cms/PPD/PdmV/work/McM/submit/B2G-RunIISpring15DR74-00001/B2G-RunIISpring15DR74-00001_3_cfg.py --no_exec -n 82
import FWCore.ParameterSet.Config as cms

process = cms.Process('PAT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        "/store/mc/RunIISpring15DR74/RSGravToZZToLLQQ_kMpl01_M-2000_TuneCUETP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/02450391-E205-E511-8314-AC853D9DACF7.root", 
        "/store/mc/RunIISpring15DR74/RSGravToZZToLLQQ_kMpl01_M-2000_TuneCUETP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/063C73E8-CA05-E511-A39D-525400E50CD9.root",
        "/store/mc/RunIISpring15DR74/RSGravToZZToLLQQ_kMpl01_M-2000_TuneCUETP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/380584A1-E205-E511-9DE9-AC853D9DAC25.root",
        "/store/mc/RunIISpring15DR74/RSGravToZZToLLQQ_kMpl01_M-2000_TuneCUETP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/68C197FF-CA05-E511-A52D-F45214C748B6.root",
        "/store/mc/RunIISpring15DR74/RSGravToZZToLLQQ_kMpl01_M-2000_TuneCUETP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/7EAED283-F405-E511-B2B0-000F530E4644.root",
        "/store/mc/RunIISpring15DR74/RSGravToZZToLLQQ_kMpl01_M-2000_TuneCUETP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/A643BDA7-E205-E511-8BB8-000F532734B0.root",
        "/store/mc/RunIISpring15DR74/RSGravToZZToLLQQ_kMpl01_M-2000_TuneCUETP8M1_13TeV-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v1/70000/EAAE079A-E205-E511-8BAF-AC853D9DAC41.root",
        #"/store/relval/CMSSW_7_4_0/RelValProdTTbar/AODSIM/MCRUN1_74_V4-v1/00000/3C61F496-DEDA-E411-A468-0025905A6134.root",
        #           "/store/relval/CMSSW_7_4_0/RelValProdTTbar/AODSIM/MCRUN1_74_V4-v1/00000/FAB2339D-DEDA-E411-80F7-002618943962.root"
        ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:82'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:ZZ_MuTauEleTau_2000.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
process.load('Configuration.StandardSequences.PATMC_cff')

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions
