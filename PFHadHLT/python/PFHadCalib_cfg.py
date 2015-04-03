import FWCore.ParameterSet.Config as cms

process = cms.Process("PFHadCalib")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.load("PFHadCalib.PFHadHLT.hlt_testJetsMET_Method3_keepAllCollectionsForPFHadronCalibration")

from PFHadCalib.PFHadHLT.hlt_testJetsMET_Method3_keepAllCollectionsForPFHadronCalibration import *
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring(
#        'file:/d3/scratch/swlee/works/JERC/CMSSW_7_4_0_pre9/src/PFHadCalib/PFHadHLT/python/outputA.root'
#    )
#)

process.pfhadhlt = cms.EDAnalyzer('PFHadHLT',
                   HLTPFCandidates  = cms.InputTag("hltParticleFlow","","TEST"),
                   PFSimParticles   = cms.InputTag("particleFlowSimParticle","","TEST"),

                   ptMin = cms.double(1.),                     # Minimum pt
                   pMin = cms.double(3.),                      # Minimum p
                   nPixMin = cms.int32(2),                     # Nb of pixel hits
                   #nHitMin = cms.vint32(14,17,20,17),          # Nb of track hits
                   nHitMin = cms.vint32(9,9,9,9),          # Nb of track hits
                   nEtaMin = cms.vdouble(1.4, 1.6, 2.0, 2.4),  # in these eta ranges
                   hcalMin = cms.double(1.),                   # Minimum hcal energy
                   ecalMax = cms.double(0.2),                  # Maximum ecal energy 
                   rootOutputFile = cms.string('PFHadCalibration.root')
)


#process.p = cms.Path(process.pfhadhlt)

process.AOutput = cms.EndPath( process.hltPreAOutput + process.pfhadhlt )
