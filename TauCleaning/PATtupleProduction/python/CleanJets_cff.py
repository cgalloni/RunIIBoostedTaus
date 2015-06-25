import FWCore.ParameterSet.Config as cms

##MU-TAUH CHANNEL
ak4PFJetsNoMu = cms.EDProducer(
    "CleanJetsMuTauProducer"
    )
CleanJetsMuTauSequence = cms.Sequence(ak4PFJetsNoMu)


##ELE-TAUH CHANNEL
ak4PFJetsNoEle = cms.EDProducer(
    "CleanJetsETauProducer"
    )
CleanJetsETauSequence = cms.Sequence(ak4PFJetsNoEle)
